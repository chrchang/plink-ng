// This library is part of PLINK 2.00, copyright (C) 2005-2023 Shaun Purcell,
// Christopher Chang.
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


#include "pgenlib_write.h"

#include <unistd.h>  // unlink()

#ifdef __cplusplus
namespace plink2 {
#endif

static inline PgenWriterCommon* GetPwcp(STPgenWriter* spgwp) {
  return &GET_PRIVATE(*spgwp, pwc);
}

static inline FILE** GetPgenOutfilep(STPgenWriter* spgwp) {
  return &GET_PRIVATE(*spgwp, pgen_outfile);
}

static inline FILE** GetPgiOrFinalPgenOutfilep(STPgenWriter* spgwp) {
  return &GET_PRIVATE(*spgwp, pgi_or_final_pgen_outfile);
}

static inline char** GetFnameBufp(STPgenWriter* spgwp) {
  return &GET_PRIVATE(*spgwp, fname_buf);
}

void GenovecInvertCopyUnsafe(const uintptr_t* __restrict genovec, uint32_t sample_ct, uintptr_t* __restrict genovec_inverted_copy) {
  // flips 0 to 2 and vice versa.
  // "unsafe" because trailing bits are not zeroed out.
  const uint32_t vec_ct = NypCtToVecCt(sample_ct);
  assert(VecIsAligned(genovec));
  const VecW not_m1 = VCONST_W(kMaskAAAA);
  const VecW* vin_ptr = R_CAST(const VecW*, genovec);
  VecW* vout_ptr = R_CAST(VecW*, genovec_inverted_copy);
  for (uint32_t vidx = 0; vidx != vec_ct; ++vidx) {
    const VecW cur_vec = vin_ptr[vidx];
    // flip high bit iff low bit is unset
    vout_ptr[vidx] = cur_vec ^ vecw_and_notfirst(vecw_slli(cur_vec, 1), not_m1);
  }
}

void PreinitSpgw(STPgenWriter* spgwp) {
  *GetPgenOutfilep(spgwp) = nullptr;
  *GetPgiOrFinalPgenOutfilep(spgwp) = nullptr;
  *GetFnameBufp(spgwp) = nullptr;
}

void PreinitMpgw(MTPgenWriter* mpgwp) {
  mpgwp->pgen_outfile = nullptr;
  mpgwp->pgi_or_final_pgen_outfile = nullptr;
  mpgwp->fname_buf = nullptr;
}

// OpenFail and WriteFail errors can now refer to either the .pgen, .pgen.tmp,
// or .pgen.pgi file.
PglErr PwcInitPhase1(const char* __restrict fname, uintptr_t* explicit_nonref_flags, uint32_t variant_ct_limit, uint32_t sample_ct, PgenWriteMode write_mode, PgenGlobalFlags phase_dosage_gflags, uint32_t nonref_flags_storage, uintptr_t vrec_len_byte_ct, PgenWriterCommon* pwcp, FILE** pgen_outfile_ptr, FILE** pgi_or_final_pgen_outfile_ptr, char** fname_buf_ptr) {
  pwcp->explicit_nonref_flags = nullptr;
  if (nonref_flags_storage == 3) {
    if (unlikely(!explicit_nonref_flags)) {
      return kPglRetImproperFunctionCall;
    }
    pwcp->explicit_nonref_flags = explicit_nonref_flags;
  }
  pwcp->variant_ct_limit = variant_ct_limit;
  pwcp->sample_ct = sample_ct;
  pwcp->phase_dosage_gflags = phase_dosage_gflags;
  pwcp->nonref_flags_storage = nonref_flags_storage;
  pwcp->vrec_len_byte_ct = vrec_len_byte_ct;
#ifndef NDEBUG
  pwcp->vblock_fpos = nullptr;
  pwcp->vrec_len_buf = nullptr;
  pwcp->vrtype_buf = nullptr;
  pwcp->fwrite_buf = nullptr;
  pwcp->fwrite_bufp = nullptr;
  pwcp->genovec_hets_buf = nullptr;
  pwcp->genovec_invert_buf = nullptr;
  pwcp->ldbase_genovec = nullptr;
  pwcp->ldbase_raregeno = nullptr;
  pwcp->ldbase_difflist_sample_ids = nullptr;
#endif
  pwcp->vidx = 0;

  *pgen_outfile_ptr = nullptr;
  *pgi_or_final_pgen_outfile_ptr = nullptr;
  *fname_buf_ptr = nullptr;
  if (write_mode != kPgenWriteBackwardSeek) {
    const uint32_t fname_slen = strlen(fname);
    if (fname_slen > kPglFnamesize - 5) {
      return kPglRetMalformedInput;
    }
    pwcp->vblock_fpos_offset = 3;
    if (write_mode == kPgenWriteAndCopy) {
      if (unlikely(pgl_malloc(fname_slen + 5, fname_buf_ptr))) {
        return kPglRetNomem;
      }
      char* fname_iter = memcpya(*fname_buf_ptr, fname, fname_slen);
      // In the kPgenWriteAndCopy case, pgen_outfile is the *temporary* .pgen
      // output.
      strcpy_k(fname_iter, ".tmp");
      *pgen_outfile_ptr = fopen(*fname_buf_ptr, FOPEN_WB);
      if (unlikely(!(*pgen_outfile_ptr))) {
        return kPglRetOpenFail;
      }
      if (unlikely(fwrite_checked("l\x1b\x20", 3, *pgen_outfile_ptr))) {
        return kPglRetWriteFail;
      }
    } else {
      // kPgenWriteSeparateIndex
      char fname_buf[kPglFnamesize];
      char* fname_iter = memcpya(fname_buf, fname, fname_slen);
      strcpy_k(fname_iter, ".pgi");
      *pgi_or_final_pgen_outfile_ptr = fopen(fname_buf, FOPEN_WB);
      if (unlikely(!(*pgi_or_final_pgen_outfile_ptr))) {
        return kPglRetOpenFail;
      }
      if (unlikely(fwrite_checked("l\x1b\x30", 3, *pgi_or_final_pgen_outfile_ptr))) {
        return kPglRetWriteFail;
      }
    }
  }
  FILE* header_ff = fopen(fname, FOPEN_WB);
  if (unlikely(!header_ff)) {
    return kPglRetOpenFail;
  }
  if (write_mode != kPgenWriteAndCopy) {
    *pgen_outfile_ptr = header_ff;
  } else {
    *pgi_or_final_pgen_outfile_ptr = header_ff;
  }
  fwrite_unlocked("l\x1b", 2, 1, header_ff);
  const int32_t third_byte = (write_mode == kPgenWriteSeparateIndex)? 0x20 : 0x10;
  if (unlikely(putc_checked(third_byte, header_ff))) {
    return kPglRetWriteFail;
  }
  if (write_mode != kPgenWriteBackwardSeek) {
    return kPglRetSuccess;
  }

  const uint32_t variant_ct = variant_ct_limit;
  const uint32_t vblock_ct = DivUp(variant_ct, kPglVblockSize);
  uintptr_t header_bytes_left = 9 + vblock_ct * sizeof(int64_t) + variant_ct * vrec_len_byte_ct;
  if (phase_dosage_gflags) {
    // 8-bit vrtypes
    header_bytes_left += variant_ct;
  } else {
    // 4-bit vrtypes
    header_bytes_left += DivUp(variant_ct, 2);
  }
  if (nonref_flags_storage == 3) {
    header_bytes_left += DivUp(variant_ct, CHAR_BIT);
  }

  // this should be the position of the first variant
  pwcp->vblock_fpos_offset = 3 + header_bytes_left;

  uintptr_t zeroed_cachelines_needed = DivUp(header_bytes_left, kCacheline);
  if (zeroed_cachelines_needed > (kPglFwriteBlockSize / kCacheline)) {
    zeroed_cachelines_needed = kPglFwriteBlockSize / kCacheline;
  }
  // could wait until fwrite_buf is allocated, and make sure it's aligned?
  unsigned char zerobuf[kPglFwriteBlockSize];
  memset(zerobuf, 0, zeroed_cachelines_needed * kCacheline);
  while (header_bytes_left > kPglFwriteBlockSize) {
    fwrite_unlocked(zerobuf, kPglFwriteBlockSize, 1, header_ff);
    header_bytes_left -= kPglFwriteBlockSize;
  }
  if (unlikely(fwrite_checked(zerobuf, header_bytes_left, header_ff))) {
    return kPglRetWriteFail;
  }
  return kPglRetSuccess;
}

uintptr_t CountSpgwAllocCachelinesRequired(uint32_t variant_ct_limit, uint32_t sample_ct, PgenGlobalFlags phase_dosage_gflags, uint32_t max_vrec_len) {
  // vblock_fpos
  const uint32_t vblock_ct_limit = DivUp(variant_ct_limit, kPglVblockSize);
  uintptr_t cachelines_required = Int64CtToCachelineCt(vblock_ct_limit);

  // vrec_len_buf
  const uintptr_t vrec_len_byte_ct = BytesToRepresentNzU32(max_vrec_len);
  cachelines_required += DivUp(variant_ct_limit * vrec_len_byte_ct, kCacheline);

  // vrtype_buf
  if (phase_dosage_gflags) {
    cachelines_required += DivUp(variant_ct_limit, kCacheline);
  } else {
    cachelines_required += DivUp(variant_ct_limit, kCacheline * 2);
  }

  // genovec_hets_buf, genovec_invert_buf, ldbase_genovec
  cachelines_required += 3 * NypCtToCachelineCt(sample_ct);

  const uint32_t max_difflist_len = 2 * (sample_ct / kPglMaxDifflistLenDivisor);
  // ldbase_raregeno
  cachelines_required += NypCtToCachelineCt(max_difflist_len);

  // ldbase_difflist_sample_ids
  cachelines_required += 1 + (max_difflist_len / kInt32PerCacheline);

  // fwrite_buf
  // + (5 + sizeof(AlleleCode)) * kPglDifflistGroupSize to avoid buffer
  // overflow in middle of difflist writing
  cachelines_required += DivUp(max_vrec_len + kPglFwriteBlockSize + (5 + sizeof(AlleleCode)) * kPglDifflistGroupSize, kCacheline);
  // possible todo: dosage (doesn't currently need an allocation, but that's
  // unlikely to remain true--e.g. get_ref_nonref_genotype_counts_and_dosages
  // tends to use workspace_vec when a function it calls doesn't use it...)

  // probable todo: multiallelic (phased) dosage
  return cachelines_required;
}

static_assert(kPglMaxAlleleCt == 255, "Need to update SpgwInitPhase1().");
PglErr SpgwInitPhase1(const char* __restrict fname, const uintptr_t* __restrict allele_idx_offsets, uintptr_t* __restrict explicit_nonref_flags, uint32_t variant_ct_limit, uint32_t sample_ct, uint32_t allele_ct_upper_bound, PgenWriteMode write_mode, PgenGlobalFlags phase_dosage_gflags, uint32_t nonref_flags_storage, STPgenWriter* spgwp, uintptr_t* alloc_cacheline_ct_ptr, uint32_t* max_vrec_len_ptr) {
  assert(variant_ct_limit);
  assert(sample_ct);

  // separate from MpgwInitPhase1's version of this computation since the
  // latter wants a better bound on the compressed size of an entire vblock
  // than max_vrec_len * kPglVblockSize...
  uint64_t max_vrec_len = NypCtToByteCt(sample_ct);
  uintptr_t max_alt_ct_p1 = 2;
  if (allele_idx_offsets) {
    const uint32_t variant_ct = variant_ct_limit;
    if (allele_idx_offsets[variant_ct] != 2 * variant_ct) {
      assert(allele_idx_offsets[0] == 0);
      assert(allele_idx_offsets[variant_ct] > 2 * variant_ct);
      max_alt_ct_p1 = 3;
      uintptr_t prev_offset = 0;
      for (uint32_t vidx = 1; vidx <= variant_ct; ++vidx) {
        const uintptr_t cur_offset = allele_idx_offsets[vidx];
        if (cur_offset - prev_offset > max_alt_ct_p1) {
          max_alt_ct_p1 = cur_offset - prev_offset;
        }
        prev_offset = cur_offset;
      }
    }
  } else if (allele_ct_upper_bound) {
    max_alt_ct_p1 = allele_ct_upper_bound;
  }
  if (max_alt_ct_p1 > 2) {
    // see comments in middle of MpgwInitPhase1()
    max_vrec_len += 2 + sizeof(AlleleCode) + GetAux1bAlleleEntryByteCt(max_alt_ct_p1, sample_ct - 1);
    // try to permit uncompressed records to be larger than this, only error
    // out when trying to write a larger compressed record.
  }
  if (phase_dosage_gflags & kfPgenGlobalHardcallPhasePresent) {
    // phasepresent, phaseinfo
    max_vrec_len += 2 * DivUp(sample_ct, CHAR_BIT);
  }
  if (phase_dosage_gflags & kfPgenGlobalDosagePresent) {
    const uint32_t dphase_gflag = (phase_dosage_gflags / kfPgenGlobalDosagePhasePresent) & 1;
    // aux3, aux7
    max_vrec_len += (1 + dphase_gflag) * DivUp(sample_ct, 8);
    // aux4, aux8
    max_vrec_len += (2 + 2 * dphase_gflag) * S_CAST(uint64_t, sample_ct);
    // todo: multiallelic dosage
  }
#ifdef __LP64__
  if (max_vrec_len >= kPglMaxBytesPerVariant) {
    max_vrec_len = kPglMaxBytesPerVariant;
  }
#else
  if (unlikely(max_vrec_len > kMaxBytesPerIO + 1 - kPglFwriteBlockSize)) {
    return kPglRetNomem;
  }
#endif
  *max_vrec_len_ptr = max_vrec_len;
  const uintptr_t vrec_len_byte_ct = BytesToRepresentNzU32(max_vrec_len);

  PgenWriterCommon* pwcp = GetPwcp(spgwp);
  FILE** pgen_outfilep = GetPgenOutfilep(spgwp);
  FILE** pgi_or_final_pgen_outfilep = GetPgiOrFinalPgenOutfilep(spgwp);
  char** fname_bufp = GetFnameBufp(spgwp);
  PglErr reterr = PwcInitPhase1(fname, explicit_nonref_flags, variant_ct_limit, sample_ct, write_mode, phase_dosage_gflags, nonref_flags_storage, vrec_len_byte_ct, pwcp, pgen_outfilep, pgi_or_final_pgen_outfilep, fname_bufp);
  if (likely(!reterr)) {
    *alloc_cacheline_ct_ptr = CountSpgwAllocCachelinesRequired(variant_ct_limit, sample_ct, phase_dosage_gflags, max_vrec_len);
  }
  return reterr;
}

static_assert(kPglMaxAlleleCt == 255, "Need to update MpgwInitPhase1().");
void MpgwInitPhase1(const uintptr_t* __restrict allele_idx_offsets, uint32_t variant_ct, uint32_t sample_ct, PgenGlobalFlags phase_dosage_gflags, uintptr_t* alloc_base_cacheline_ct_ptr, uint64_t* alloc_per_thread_cacheline_ct_ptr, uint32_t* vrec_len_byte_ct_ptr, uint64_t* vblock_cacheline_ct_ptr) {
  // variant_ct must be exact due to how MpgwFlush() works.
  assert(variant_ct);
  assert(sample_ct);
  // vblock_fpos
  const uint32_t vblock_ct = DivUp(variant_ct, kPglVblockSize);
  uint32_t alloc_base_cacheline_ct = Int64CtToCachelineCt(vblock_ct);

  // vrtype_buf
  if (phase_dosage_gflags) {
    alloc_base_cacheline_ct += DivUp(variant_ct, kCacheline);
  } else {
    alloc_base_cacheline_ct += DivUp(variant_ct, kCacheline * 2);
  }

  // pwcs
  uint64_t alloc_per_thread_cacheline_ct = DivUp(sizeof(PgenWriterCommon), kCacheline);

  // genovec_hets_buf, genovec_invert_buf, ldbase_genovec
  // (could make genovec_hets_buf allocation conditional, but not worth the
  // additional load this puts on import functions)
  alloc_per_thread_cacheline_ct += 3 * NypCtToCachelineCt(sample_ct);

  const uint32_t max_difflist_len = 2 * (sample_ct / kPglMaxDifflistLenDivisor);
  // ldbase_raregeno
  alloc_per_thread_cacheline_ct += NypCtToCachelineCt(max_difflist_len);

  // ldbase_difflist_sample_ids
  alloc_per_thread_cacheline_ct += 1 + (max_difflist_len / kInt32PerCacheline);

  uint64_t max_vrec_len = NypCtToByteCt(sample_ct);
  if (phase_dosage_gflags & kfPgenGlobalHardcallPhasePresent) {
    max_vrec_len += 2 * DivUp(sample_ct, CHAR_BIT);
  }
  const uint32_t dosage_gflag = (phase_dosage_gflags / kfPgenGlobalDosagePresent) & 1;
  const uint32_t dosage_phase_gflag = (phase_dosage_gflags / kfPgenGlobalDosagePhasePresent) & 1;
  if (dosage_gflag) {
    max_vrec_len += ((1 + dosage_phase_gflag) * DivUp(sample_ct, CHAR_BIT)) + (2 + 2 * dosage_phase_gflag) * S_CAST(uint64_t, sample_ct);
  }
  const uint32_t max_vblock_size = MINV(variant_ct, kPglVblockSize);
  // max_vrec_len is a uint64_t
  uint64_t max_vblock_byte_ct = max_vrec_len * max_vblock_size;
  if (allele_idx_offsets && (allele_idx_offsets[variant_ct] != 2 * variant_ct)) {
    assert(allele_idx_offsets[0] == 0);
    assert(allele_idx_offsets[variant_ct] > 2 * variant_ct);
    // When multiallelic variants are present, larger write buffers are needed.
    // we compute the largest possible size here.
    //
    // One way to hit the maximum size is to have exactly 1 ref/altx genotype,
    // all other genotypes are altx/alty, and all genotypes include a rarealt.
    // Then, we need 2 + sizeof(AlleleCode) bytes for the format byte and
    // aux1a, (sample_ct + 6) / 8 bytes for the bitarray at the front of aux1b,
    // and
    //   alt ct  aux1_last_quarter bytes required
    //   ------  --------------------------------
    //        2           (sample_ct - 1 + 7) / 8
    //      3-4           (sample_ct - 1 + 1) / 2
    //     5-16                   (sample_ct - 1)
    //   17-256               2 * (sample_ct - 1)
    //
    // (todo: also take multiallelic dosage into account)
    assert(!dosage_gflag);
    uintptr_t prev_offset = 0;
    uint32_t vidx = 0;
    const uint32_t extra_bytes_base = 2 + sizeof(AlleleCode) + (sample_ct + 6) / 8;
    const uint64_t extra_bytes_max = kPglMaxBytesPerVariant - max_vrec_len;
    uint64_t extra_byte_cts[4];
    uint32_t extra_alt_ceil = kPglMaxAlleleCt;

    // alt_ct == 2
    uint64_t cur_extra_byte_ct = extra_bytes_base + (sample_ct + 6) / 8;
    extra_byte_cts[0] = cur_extra_byte_ct;
    extra_byte_cts[1] = 0;  // force initialization
    extra_byte_cts[2] = 0;
    extra_byte_cts[3] = 0;
    if (cur_extra_byte_ct >= extra_bytes_max) {
      extra_alt_ceil = 2;
    } else {
      // alt_ct in [3, 4]
      cur_extra_byte_ct = extra_bytes_base + (sample_ct / 2);
      extra_byte_cts[1] = cur_extra_byte_ct;
      if (cur_extra_byte_ct >= extra_bytes_max) {
        extra_alt_ceil = 3;
      } else {
        // alt_ct in [5, 16]
        cur_extra_byte_ct = extra_bytes_base + sample_ct - 1;
        extra_byte_cts[2] = cur_extra_byte_ct;
        if (cur_extra_byte_ct >= extra_bytes_max) {
          extra_alt_ceil = 5;
        } else {
          // alt_ct in [17, 256]
          cur_extra_byte_ct = extra_bytes_base + 2 * (sample_ct - 1);
          extra_byte_cts[3] = cur_extra_byte_ct;
          if (cur_extra_byte_ct >= extra_bytes_max) {
            extra_alt_ceil = 16;
          }
        }
      }
    }
    uint32_t extra_alt_ceil_ct = 0;
    const uint64_t uncompressed_biallelic_vrec_len = max_vrec_len;
    uint32_t altx_seen_mask = 0;
    uint32_t max_alt_ct_p1 = 3;
    for (uint32_t vblock_end = vidx + kPglVblockSize; ; vblock_end += kPglVblockSize) {
      if (vblock_end > variant_ct) {
        if (vidx == variant_ct) {
          break;
        }
        vblock_end = variant_ct;
      }
      uint64_t cur_vblock_byte_ct = uncompressed_biallelic_vrec_len * (vblock_end - vidx);
      uint32_t altx_seen[4];
      ZeroU32Arr(4, altx_seen);
      for (; vidx < vblock_end; ) {
        const uintptr_t cur_offset = allele_idx_offsets[++vidx];
        const uint32_t alt_ct_p1 = cur_offset - prev_offset;
        if (alt_ct_p1 > 2) {
          if (alt_ct_p1 >= extra_alt_ceil) {
            ++extra_alt_ceil_ct;
          } else {
            // don't need to track this when we hit the ceiling
            if (alt_ct_p1 > max_alt_ct_p1) {
              max_alt_ct_p1 = alt_ct_p1;
            }

            if (alt_ct_p1 < 6) {
              altx_seen[(alt_ct_p1 != 3)] += 1;
            } else {
              altx_seen[2 + (alt_ct_p1 >= 17)] += 1;
            }
          }
        }
        prev_offset = cur_offset;
      }
      cur_vblock_byte_ct += extra_alt_ceil_ct * extra_bytes_max;
      for (uint32_t uii = 0; uii != 4; ++uii) {
        if (altx_seen[uii]) {
          const uint32_t cur_seen_ct = altx_seen[uii];
          altx_seen_mask |= 1 << uii;
          cur_vblock_byte_ct += cur_seen_ct * extra_byte_cts[uii];
        }
      }
      if (cur_vblock_byte_ct > max_vblock_byte_ct) {
        max_vblock_byte_ct = cur_vblock_byte_ct;
      }
    }
    if (extra_alt_ceil_ct) {
      max_vrec_len = kPglMaxBytesPerVariant;
    } else {
      max_vrec_len = uncompressed_biallelic_vrec_len + extra_byte_cts[bsru32(altx_seen_mask)];
    }
  }
  if (max_vrec_len >= kPglMaxBytesPerVariant) {
    max_vrec_len = kPglMaxBytesPerVariant;
    max_vblock_byte_ct = kPglMaxBytesPerVariant * S_CAST(uint64_t, max_vblock_size);
  }

  // vrec_len_buf
  // previously used overlapping uint32_t writes-to-memory, but that was
  // incompatible with multithreaded compression
  *vrec_len_byte_ct_ptr = BytesToRepresentNzU32(max_vrec_len);
  *alloc_base_cacheline_ct_ptr = alloc_base_cacheline_ct + DivUp(S_CAST(uintptr_t, variant_ct) * (*vrec_len_byte_ct_ptr), kCacheline);

  // main write buffer
  *vblock_cacheline_ct_ptr = DivUpU64(max_vblock_byte_ct, kCacheline);
  *alloc_per_thread_cacheline_ct_ptr = alloc_per_thread_cacheline_ct + (*vblock_cacheline_ct_ptr);
}


void PwcInitPhase2(uintptr_t fwrite_cacheline_ct, uint32_t thread_ct, PgenWriterCommon** pwcs, unsigned char* pwc_alloc) {
  const uint32_t variant_ct_limit = pwcs[0]->variant_ct_limit;

  unsigned char* alloc_iter = pwc_alloc;
  const uint32_t vblock_ct = DivUp(variant_ct_limit, kPglVblockSize);
  const PgenGlobalFlags phase_dosage_gflags = pwcs[0]->phase_dosage_gflags;
  uint32_t vrtype_buf_byte_ct;
  if (phase_dosage_gflags) {
    vrtype_buf_byte_ct = RoundUpPow2(variant_ct_limit, kCacheline);
  } else {
    vrtype_buf_byte_ct = DivUp(variant_ct_limit, kCacheline * 2) * kCacheline;
  }
  pwcs[0]->vblock_fpos = S_CAST(uint64_t*, arena_alloc_raw_rd(vblock_ct * sizeof(int64_t), &alloc_iter));

  pwcs[0]->vrec_len_buf = S_CAST(unsigned char*, arena_alloc_raw_rd(variant_ct_limit * pwcs[0]->vrec_len_byte_ct, &alloc_iter));

  pwcs[0]->vrtype_buf = S_CAST(uintptr_t*, arena_alloc_raw(vrtype_buf_byte_ct, &alloc_iter));
  // the PwcAppend... functions assume these bytes are zeroed out
  memset(pwcs[0]->vrtype_buf, 0, vrtype_buf_byte_ct);

  const uint32_t sample_ct = pwcs[0]->sample_ct;
  const uint32_t genovec_byte_alloc = NypCtToCachelineCt(sample_ct) * kCacheline;
  const uint32_t max_difflist_len = 2 * (sample_ct / kPglMaxDifflistLenDivisor);
  for (uint32_t tidx = 0; tidx != thread_ct; ++tidx) {
    if (tidx) {
      pwcs[tidx]->vblock_fpos = pwcs[0]->vblock_fpos;
      pwcs[tidx]->vrec_len_buf = pwcs[0]->vrec_len_buf;
      pwcs[tidx]->vrtype_buf = pwcs[0]->vrtype_buf;
    }
    pwcs[tidx]->genovec_hets_buf = S_CAST(uintptr_t*, arena_alloc_raw(genovec_byte_alloc, &alloc_iter));
    pwcs[tidx]->genovec_invert_buf = S_CAST(uintptr_t*, arena_alloc_raw(genovec_byte_alloc, &alloc_iter));
    pwcs[tidx]->ldbase_genovec = S_CAST(uintptr_t*, arena_alloc_raw(genovec_byte_alloc, &alloc_iter));

    pwcs[tidx]->ldbase_raregeno = S_CAST(uintptr_t*, arena_alloc_raw(NypCtToCachelineCt(max_difflist_len) * kCacheline, &alloc_iter));
    pwcs[tidx]->ldbase_difflist_sample_ids = S_CAST(uint32_t*, arena_alloc_raw_rd((max_difflist_len + 1) * sizeof(int32_t), &alloc_iter));

    pwcs[tidx]->fwrite_buf = alloc_iter;
    pwcs[tidx]->fwrite_bufp = alloc_iter;
    alloc_iter = &(alloc_iter[fwrite_cacheline_ct * kCacheline]);
  }
}

void SpgwInitPhase2(uint32_t max_vrec_len, STPgenWriter* spgwp, unsigned char* spgw_alloc) {
  const uintptr_t fwrite_cacheline_ct = DivUp(max_vrec_len + kPglFwriteBlockSize + (5 + sizeof(AlleleCode)) * kPglDifflistGroupSize, kCacheline);
  PgenWriterCommon* pwcp = GetPwcp(spgwp);
  PwcInitPhase2(fwrite_cacheline_ct, 1, &pwcp, spgw_alloc);
}

PglErr MpgwInitPhase2(const char* __restrict fname, uintptr_t* __restrict explicit_nonref_flags, uint32_t variant_ct, uint32_t sample_ct, PgenWriteMode write_mode, PgenGlobalFlags phase_dosage_gflags, uint32_t nonref_flags_storage, uint32_t vrec_len_byte_ct, uintptr_t vblock_cacheline_ct, uint32_t thread_ct, unsigned char* mpgw_alloc, MTPgenWriter* mpgwp) {
  assert(thread_ct);
  const uintptr_t pwc_byte_ct = RoundUpPow2(sizeof(PgenWriterCommon), kCacheline);
  for (uint32_t tidx = 0; tidx != thread_ct; ++tidx) {
    mpgwp->pwcs[tidx] = R_CAST(PgenWriterCommon*, &(mpgw_alloc[tidx * pwc_byte_ct]));
  }
  PglErr reterr = PwcInitPhase1(fname, explicit_nonref_flags, variant_ct, sample_ct, write_mode, phase_dosage_gflags, nonref_flags_storage, vrec_len_byte_ct, mpgwp->pwcs[0], &(mpgwp->pgen_outfile), &(mpgwp->pgi_or_final_pgen_outfile), &(mpgwp->fname_buf));
  if (unlikely(reterr)) {
    return reterr;
  }
  mpgwp->thread_ct = thread_ct;
  for (uint32_t tidx = 1; tidx != thread_ct; ++tidx) {
    *(mpgwp->pwcs[tidx]) = *(mpgwp->pwcs[0]);
    mpgwp->pwcs[tidx]->vidx = tidx * kPglVblockSize;
  }
  PwcInitPhase2(vblock_cacheline_ct, thread_ct, mpgwp->pwcs, &(mpgw_alloc[thread_ct * pwc_byte_ct]));
  return kPglRetSuccess;
}


void CountLdAndInvertedLdDiffs(const uintptr_t* __restrict ldbase_genovec, const uintptr_t* __restrict genovec, uint32_t sample_ct, uint32_t* ld_diff_ctp, uint32_t* ld_inv_diff_ctp) {
  // Requires trailing bits to be zeroed out.
  const uint32_t word_ct = NypCtToWordCt(sample_ct);
  const uintptr_t* genovec_end = &(genovec[word_ct]);
  uint32_t ld_diff_ct = 0;
  uint32_t ld_inv_diff_ct = 0;
  // construct the words we want to popcount_nypvec_01 on the fly
  const VecW m1 = VCONST_W(kMask5555);
  const VecW m2 = VCONST_W(kMask3333);
  const VecW m4 = VCONST_W(kMask0F0F);
  const VecW* ldbase_vvec_iter = R_CAST(const VecW*, ldbase_genovec);
  const VecW* geno_vvec_iter = R_CAST(const VecW*, genovec);
  VecW acc_ld = vecw_setzero();
  VecW acc_ld_inv = vecw_setzero();
  uintptr_t cur_incr = 60;
  for (uintptr_t full_vecs_left = 3 * (word_ct / (3 * kWordsPerVec)); ; full_vecs_left -= cur_incr) {
    if (full_vecs_left < 60) {
      if (!full_vecs_left) {
        ld_diff_ct = HsumW(acc_ld);
        ld_inv_diff_ct = HsumW(acc_ld_inv);
        break;
      }
      cur_incr = full_vecs_left;
    }
    VecW inner_acc_ld = vecw_setzero();
    VecW inner_acc_ld_inv = vecw_setzero();
    const VecW* geno_vvec_stop = &(geno_vvec_iter[cur_incr]);
    do {
      VecW loader_ldbase1 = *ldbase_vvec_iter++;
      VecW loader_geno1 = *geno_vvec_iter++;
      VecW loader_ldbase2 = *ldbase_vvec_iter++;
      VecW loader_geno2 = *geno_vvec_iter++;
      VecW xor1 = loader_ldbase1 ^ loader_geno1;
      VecW xor2 = loader_ldbase2 ^ loader_geno2;
      VecW xor_shifted1 = vecw_srli(xor1, 1);
      VecW xor_shifted2 = vecw_srli(xor2, 1);
      // xor(_low)  xor_shifted  loader_geno   result
      //         1                                  1
      //         0            0            0        1
      //         0            0            1        0
      //         0            1            0        0
      //         0            1            1        1
      // gah, don't see a way to avoid throwing in an extra xor for
      // loader_geno...
      VecW count_ld_inv = (xor1 | (xor_shifted1 ^ loader_geno1 ^ m1)) & m1;
      loader_ldbase1 = *ldbase_vvec_iter++;
      VecW count_ld = (xor1 | xor_shifted1) & m1;
      loader_geno1 = *geno_vvec_iter++;
      count_ld_inv = count_ld_inv + ((xor2 | (xor_shifted2 ^ loader_geno2 ^ m1)) & m1);
      xor1 = loader_ldbase1 ^ loader_geno1;
      count_ld = count_ld + ((xor2 | xor_shifted2) & m1);
      xor_shifted1 = vecw_srli(xor1, 1);
      count_ld_inv = count_ld_inv + ((xor1 | (xor_shifted1 ^ loader_geno1 ^ m1)) & m1);
      count_ld = count_ld + ((xor1 | xor_shifted1) & m1);
      // now count_ld and count_ld_inv each have 64 2-bit values from 0-3

      count_ld_inv = (count_ld_inv & m2) + (vecw_srli(count_ld_inv, 2) & m2);
      count_ld = (count_ld & m2) + (vecw_srli(count_ld, 2) & m2);
      // now they have 32 4-bit values from 0-6

      inner_acc_ld_inv = inner_acc_ld_inv + ((count_ld_inv + vecw_srli(count_ld_inv, 4)) & m4);
      inner_acc_ld = inner_acc_ld + ((count_ld + vecw_srli(count_ld, 4)) & m4);
    } while (geno_vvec_iter < geno_vvec_stop);
    const VecW m0 = vecw_setzero();
    acc_ld_inv = acc_ld_inv + vecw_bytesum(inner_acc_ld_inv, m0);
    acc_ld = acc_ld + vecw_bytesum(inner_acc_ld, m0);
  }
  const uintptr_t* ldbase_iter = R_CAST(const uintptr_t*, ldbase_vvec_iter);
  const uintptr_t* genovec_iter = R_CAST(const uintptr_t*, geno_vvec_iter);
  while (genovec_iter < genovec_end) {
    uintptr_t ldbase_word = *ldbase_iter++;
    uintptr_t geno_word = *genovec_iter++;
    uintptr_t xor_result = ldbase_word ^ geno_word;
    uintptr_t xor_result_shifted = xor_result >> 1;
    ld_diff_ct += Popcount01Word((xor_result | xor_result_shifted) & kMask5555);
    ld_inv_diff_ct += Popcount01Word((xor_result | (xor_result_shifted ^ (~geno_word))) & kMask5555);
  }
  *ld_diff_ctp = ld_diff_ct;
  // trailing entries in last word are always "different"
  *ld_inv_diff_ctp = ld_inv_diff_ct - ((-sample_ct) & (kBitsPerWordD2 - 1));
}

void CountLdAndInvertedLdDiffsList(const uintptr_t* __restrict ldbase_raregeno, const uint32_t* __restrict ldbase_difflist_sample_ids, const uintptr_t* __restrict raregeno, const uint32_t* __restrict difflist_sample_ids, uint32_t ldbase_difflist_len, uint32_t difflist_len, uint32_t* ld_diff_ctp, uint32_t* ld_inv_diff_ctp) {
  // assumes ldbase_difflist_sample_ids[ldbase_difflist_len] == sample_ct

  // only the count(s) with aligned common_geno values are valid.  e.g. if
  // ldbase_common_geno and difflist_common_geno are both zero, the ld_inv_diff
  // return value can be anything, while if they're both three, ld_diff and
  // ld_inv_diff are both accurate.

  // some similarities to ParseLdAndMergeDifflistSubset(), but much simpler.
  // noticeably slower than CountLdAndInvertedLdDiffs() when the lists aren't
  // tiny.
  // possible todo: take threshold into account?

  uint32_t collision_ct = 0;
  uint32_t ld_diff_ct = 0;
  uint32_t ld_inv_diff_ct = 0;
  uint32_t ldbase_sample_idx = ldbase_difflist_sample_ids[0];
  uint32_t ldbase_difflist_idx = 1;
  // this loop is a bit slow.  attempt to bail halfway through?
  for (uint32_t difflist_idx = 0; difflist_idx != difflist_len; ++difflist_idx) {
    const uint32_t raw_sample_idx = difflist_sample_ids[difflist_idx];
    while (ldbase_sample_idx < raw_sample_idx) {
      ldbase_sample_idx = ldbase_difflist_sample_ids[ldbase_difflist_idx++];
    }
    if (ldbase_sample_idx > raw_sample_idx) {
      continue;
    }
    const uint32_t cur_raregeno = GetNyparrEntry(raregeno, difflist_idx);
    const uint32_t cur_ldbase_raregeno = GetNyparrEntry(ldbase_raregeno, ldbase_difflist_idx - 1);
    const uint32_t cur_inv_raregeno = (6 - cur_raregeno) & 3;
    ld_diff_ct += (cur_ldbase_raregeno != cur_raregeno);
    ldbase_sample_idx = ldbase_difflist_sample_ids[ldbase_difflist_idx++];
    ++collision_ct;
    ld_inv_diff_ct += (cur_ldbase_raregeno != cur_inv_raregeno);
  }
  // no more collisions, don't actually need to look at rest of
  // ldbase_difflist
  const uint32_t base_diff_ct = ldbase_difflist_len + difflist_len - 2 * collision_ct;
  *ld_diff_ctp = base_diff_ct + ld_diff_ct;
  *ld_inv_diff_ctp = base_diff_ct + ld_inv_diff_ct;
}

uint32_t SaveLdDifflist(const uintptr_t* __restrict genovec, const uintptr_t* __restrict ldbase_genovec, uintptr_t common_geno, uint32_t difflist_len, PgenWriterCommon* pwcp) {
  unsigned char* fwrite_bufp = pwcp->fwrite_bufp;
  if (!difflist_len) {
    *fwrite_bufp = 0;
    pwcp->fwrite_bufp = &(fwrite_bufp[1]);
    return 1;
  }
  unsigned char* fwrite_bufp_start = fwrite_bufp;
  fwrite_bufp = Vint32Append(difflist_len, fwrite_bufp);
  const uint32_t sample_id_byte_ct = BytesToRepresentNzU32(pwcp->sample_ct);
  const uintptr_t common_geno_word = common_geno * kMask5555;
  const uint32_t group_ct = DivUp(difflist_len, kPglDifflistGroupSize);
  unsigned char* group_first_sample_ids_iter = fwrite_bufp;
  unsigned char* extra_byte_cts_iter = &(fwrite_bufp[group_ct * sample_id_byte_ct]);
#ifdef NO_UNALIGNED
#  error "Unaligned accesses in SaveLdDifflist()."
#endif
  uintptr_t* raregeno_iter = R_CAST(uintptr_t*, &(extra_byte_cts_iter[group_ct - 1]));
  fwrite_bufp = &(extra_byte_cts_iter[group_ct + (difflist_len - 1) / 4]);
  unsigned char* last_group_vint_start = fwrite_bufp;
  uintptr_t raregeno_word = 0;
  uint32_t last_sample_idx = 0;
  uint32_t difflist_idx = 0;
  for (uint32_t widx = 0; ; ++widx) {
    const uintptr_t cur_geno_word = genovec[widx];
    uintptr_t xor_word = ldbase_genovec? ldbase_genovec[widx] : common_geno_word;
    xor_word ^= cur_geno_word;
    if (xor_word) {
      const uint32_t sample_idx_base = widx * kBitsPerWordD2;
      do {
        const uint32_t sample_idx_lowbits = ctzw(xor_word) / 2;
        const uint32_t new_sample_idx = sample_idx_base + sample_idx_lowbits;
        raregeno_word |= ((cur_geno_word >> (2 * sample_idx_lowbits)) & 3) << (2 * (difflist_idx % kBitsPerWordD2));
        if (!(difflist_idx % kPglDifflistGroupSize)) {
          group_first_sample_ids_iter = memcpyua(group_first_sample_ids_iter, &new_sample_idx, sample_id_byte_ct);
          if (difflist_idx) {
            *extra_byte_cts_iter++ = S_CAST(uintptr_t, fwrite_bufp - last_group_vint_start) - (kPglDifflistGroupSize - 1);
          }
          last_group_vint_start = fwrite_bufp;
        } else {
          assert(new_sample_idx >= last_sample_idx + 1);
          fwrite_bufp = Vint32Append(new_sample_idx - last_sample_idx, fwrite_bufp);
        }
        ++difflist_idx;
        last_sample_idx = new_sample_idx;
        if (difflist_idx == difflist_len) {
          SubwordStore(raregeno_word, 1 + (((difflist_len - 1) / 4) % sizeof(intptr_t)), raregeno_iter);
          pwcp->fwrite_bufp = fwrite_bufp;
          return fwrite_bufp - fwrite_bufp_start;
        }
        if (!(difflist_idx % kBitsPerWordD2)) {
          *raregeno_iter++ = raregeno_word;
          raregeno_word = 0;
        }
        xor_word &= (~(3 * k1LU)) << (2 * sample_idx_lowbits);
      } while (xor_word);
    }
  }
}

void OnebitPreprocessBuf(const uintptr_t* __restrict genovec, uint32_t sample_ct, uint32_t common2_code, uintptr_t* __restrict genovec_buf) {
  assert(sample_ct);
  const uint32_t vec_ct = NypCtToVecCt(sample_ct);
  // todo: look for better ways to perform some of these operations
  const VecW* geno_vvec_iter = R_CAST(const VecW*, genovec);
  const VecW* geno_vvec_end = &(geno_vvec_iter[vec_ct]);
  VecW* write_iter = R_CAST(VecW*, genovec_buf);
  const VecW m1 = VCONST_W(kMask5555);
  if (common2_code < 5) {
    if (common2_code == 1) {
      // 11 -> 10, everything else unchanged
      // todo: check if these loops are actually faster as simple while loops
      // todo: check if it's better to unroll these loops to process 2 __m128is
      //       at a time
      do {
        const VecW cur_geno = *geno_vvec_iter++;
        *write_iter++ = vecw_and_notfirst(m1 & vecw_srli(cur_geno, 1), cur_geno);
      } while (geno_vvec_iter < geno_vvec_end);
    } else if (common2_code == 3) {
      // 00 -> 00, 01 -> 10, 10 -> 10, 11 -> 01
      do {
        const VecW cur_geno = *geno_vvec_iter++;
        const VecW cur_geno_rshift = vecw_srli(cur_geno, 1);
        const VecW cur_geno_xor_masked = (cur_geno ^ cur_geno_rshift) & m1;
        const VecW cur_geno_or_masked = (cur_geno | cur_geno_rshift) & m1;
        *write_iter++ = cur_geno_xor_masked + cur_geno_or_masked;
      } while (geno_vvec_iter < geno_vvec_end);
    } else {
      assert(common2_code == 2);
      // 00 -> 00, 01 -> 10, 10 -> 01, 11 -> 10
      do {
        const VecW cur_geno = *geno_vvec_iter++;
        const VecW cur_geno_or_masked = (cur_geno | vecw_srli(cur_geno, 1)) & m1;
        const VecW cur_geno_lowbits = cur_geno & m1;
        *write_iter++ = cur_geno_lowbits + cur_geno_or_masked;
      } while (geno_vvec_iter < geno_vvec_end);
    }
  } else {
    if (common2_code == 5) {
      // 00 -> 10, 01 -> 00, 10 -> 01, 11 -> 10
      do {
        const VecW cur_geno = *geno_vvec_iter++;
        const VecW cur_geno_rshift = vecw_srli(cur_geno, 1);
        const VecW cur_geno_not_xor_masked = vecw_and_notfirst(cur_geno ^ cur_geno_rshift, m1);
        const VecW cur_geno_rshift_masked = cur_geno_rshift & m1;
        *write_iter++ = cur_geno_not_xor_masked + (cur_geno_not_xor_masked | cur_geno_rshift_masked);
      } while (geno_vvec_iter < geno_vvec_end);
    } else if (common2_code == 9) {
      // 00 -> 10, 01 -> 10, 10 -> 00, 11 -> 01
      const VecW not_m1 = VCONST_W(kMaskAAAA);
      do {
        const VecW cur_geno = *geno_vvec_iter++;
        *write_iter++ = (cur_geno ^ not_m1) - vecw_and_notfirst(not_m1, vecw_and_notfirst(vecw_srli(cur_geno, 1), cur_geno));
      } while (geno_vvec_iter < geno_vvec_end);
    } else {
      assert(common2_code == 6);
      // 00 -> 10, 01 -> 00, 10 -> 10, 11 -> 01
      do {
        const VecW cur_geno = *geno_vvec_iter++;
        const VecW cur_geno_not_lowbits = vecw_and_notfirst(cur_geno, m1);
        const VecW cur_geno_rshift_masked = vecw_srli(cur_geno, 1) & m1;
        *write_iter++ = cur_geno_not_lowbits + (cur_geno_not_lowbits | cur_geno_rshift_masked);
      } while (geno_vvec_iter < geno_vvec_end);
    }
  }
}

uint32_t SaveOnebit(const uintptr_t* __restrict genovec, uint32_t common2_code, uint32_t onebit_difflist_len, PgenWriterCommon* pwcp) {
  // Uses ldbase_genovec as a temporary buffer.

  // common2_code is expected to have the difference between the common
  // genotype values in bits 0-1, and the smaller common genotype value in bits
  // 2-3.
  unsigned char* fwrite_bufp_start = pwcp->fwrite_bufp;
  *fwrite_bufp_start = common2_code;
  const uint32_t sample_ct = pwcp->sample_ct;
  uintptr_t* __restrict genovec_buf = pwcp->ldbase_genovec;
  // There's a 4-byte-interleaved format which is slightly more efficient for
  // unsubsetted handling (~10 fewer compression/decompression operations per
  // 32 genotypes), but that's only a 1-2% speedup, which probably isn't worth
  // the more annoying subsetting.
  //
  // Any 10s and 11s are saved as 00 in this part.
  // Similar logic is used to handle the other five possibilities (00/10,
  // 00/11, 01/10, 01/11, 10/11); all of them should be expected to actually
  // happen.  (E.g. 01/11 can happen at a high MAF variant when there's lots of
  // missing data.)  To reduce branching, we preprocess genovec_buf to have
  // even bit set iff the corresponding genotype is equal to the high common
  // genotype value, and odd bit set iff the corresponding genotype is one of
  // the two uncommon values.  (There may be a better way to do this, analogous
  // to the simpler decompression algorithm.)
  OnebitPreprocessBuf(genovec, sample_ct, common2_code, genovec_buf);
  ZeroTrailingNyps(sample_ct, genovec_buf);
  const uint32_t word_read_ct = NypCtToWordCt(sample_ct);
#ifdef NO_UNALIGNED
#  error "Unaligned accesses in SaveOnebit()."
#endif
  Halfword* fwrite_bufp_alias_halfword = R_CAST(Halfword*, &(fwrite_bufp_start[1]));
  PackWordsToHalfwordsMask(genovec_buf, word_read_ct, fwrite_bufp_alias_halfword);
  const uint32_t onebit_block_len = (sample_ct + 15) / CHAR_BIT;
  unsigned char* fwrite_bufp = Vint32Append(onebit_difflist_len, &(fwrite_bufp_start[onebit_block_len]));
  // the rest is almost identical to SaveLdDifflist()
  if (!onebit_difflist_len) {
    pwcp->fwrite_bufp = fwrite_bufp;
    return (onebit_block_len + 1);
  }
  const uint32_t sample_id_byte_ct = BytesToRepresentNzU32(pwcp->sample_ct);
  const uint32_t group_ct = DivUp(onebit_difflist_len, kPglDifflistGroupSize);
  unsigned char* group_first_sample_ids_iter = fwrite_bufp;
  unsigned char* extra_byte_cts_iter = &(fwrite_bufp[group_ct * sample_id_byte_ct]);
  uintptr_t* raregeno_iter = R_CAST(uintptr_t*, &(extra_byte_cts_iter[group_ct - 1]));
  fwrite_bufp = &(extra_byte_cts_iter[group_ct + (onebit_difflist_len - 1) / 4]);
  unsigned char* last_group_vint_start = fwrite_bufp;
  uintptr_t raregeno_word = 0;
  uint32_t last_sample_idx = 0;
  uint32_t difflist_idx = 0;
  for (uint32_t widx = 0; ; ++widx) {
    uintptr_t xor_word = genovec_buf[widx] & kMaskAAAA;
    if (xor_word) {
      const uintptr_t cur_geno_word = genovec[widx];
      const uint32_t sample_idx_base = widx * kBitsPerWordD2;

      do {
        const uint32_t sample_idx_lowbits = ctzw(xor_word) / 2;
        const uint32_t new_sample_idx = sample_idx_base + sample_idx_lowbits;
        if (!(difflist_idx % kBitsPerWordD2)) {
          if (!(difflist_idx % kPglDifflistGroupSize)) {
            SubU32StoreMov(new_sample_idx, sample_id_byte_ct, &group_first_sample_ids_iter);
            if (difflist_idx) {
              *extra_byte_cts_iter++ = S_CAST(uintptr_t, fwrite_bufp - last_group_vint_start) - (kPglDifflistGroupSize - 1);
              *raregeno_iter++ = raregeno_word;
              raregeno_word = 0;
            }
            last_group_vint_start = fwrite_bufp;
            goto SaveOnebit_skip_delta_write;
          }
          *raregeno_iter++ = raregeno_word;
          raregeno_word = 0;
        }
        assert(new_sample_idx >= last_sample_idx + 1);
        fwrite_bufp = Vint32Append(new_sample_idx - last_sample_idx, fwrite_bufp);
      SaveOnebit_skip_delta_write:
        raregeno_word |= ((cur_geno_word >> (2 * sample_idx_lowbits)) & 3) << (2 * (difflist_idx % kBitsPerWordD2));
        ++difflist_idx;
        last_sample_idx = new_sample_idx;
        xor_word &= xor_word - 1;
      } while (xor_word);
      // trailing bits of genovec_buf guaranteed to be zeroed out
      if (difflist_idx == onebit_difflist_len) {
        SubwordStore(raregeno_word, 1 + (((onebit_difflist_len - 1) / 4) % sizeof(intptr_t)), raregeno_iter);
        pwcp->fwrite_bufp = fwrite_bufp;
        return fwrite_bufp - fwrite_bufp_start;
      }
    }
  }
}

// returns vrec_len
uint32_t PwcAppendBiallelicGenovecMain(const uintptr_t* __restrict genovec, uint32_t vidx, PgenWriterCommon* pwcp, uint32_t* het_ct_ptr, uint32_t* altxy_ct_ptr, unsigned char* vrtype_ptr) {
  const uint32_t sample_ct = pwcp->sample_ct;
  assert((!(sample_ct % kBitsPerWordD2)) || (!(genovec[sample_ct / kBitsPerWordD2] >> (2 * (sample_ct % kBitsPerWordD2)))));
  STD_ARRAY_DECL(uint32_t, 4, genocounts);
  GenoarrCountFreqsUnsafe(genovec, sample_ct, genocounts);
  if (het_ct_ptr) {
    *het_ct_ptr = genocounts[1];
    if (altxy_ct_ptr) {
      *altxy_ct_ptr = genocounts[2];
    }
  }
  uint32_t most_common_geno = (genocounts[1] > genocounts[0]);
  uint32_t second_most_common_geno = 1 - most_common_geno;
  uint32_t largest_geno_ct = genocounts[most_common_geno];
  uint32_t second_largest_geno_ct = genocounts[second_most_common_geno];
  for (uint32_t cur_geno = 2; cur_geno != 4; ++cur_geno) {
    const uint32_t cur_geno_ct = genocounts[cur_geno];
    if (cur_geno_ct > second_largest_geno_ct) {
      if (cur_geno_ct > largest_geno_ct) {
        second_largest_geno_ct = largest_geno_ct;
        second_most_common_geno = most_common_geno;
        largest_geno_ct = cur_geno_ct;
        most_common_geno = cur_geno;
      } else {
        second_largest_geno_ct = cur_geno_ct;
        second_most_common_geno = cur_geno;
      }
    }
  }
  const uint32_t difflist_len = sample_ct - largest_geno_ct;
  const uint32_t rare_2_geno_ct_sum = difflist_len - second_largest_geno_ct;
  // average of 10-11 bits per difflist entry
  const uint32_t sample_ctd8 = sample_ct / 8;
  const uint32_t sample_ctd64 = sample_ct / 64;
  uint32_t max_difflist_len = sample_ctd8 - 2 * sample_ctd64 + rare_2_geno_ct_sum;
  if (max_difflist_len > sample_ctd8) {
    max_difflist_len = sample_ctd8;
  }
  const uint32_t difflist_viable = (most_common_geno != 1) && (difflist_len <= max_difflist_len);

  uintptr_t* ldbase_genovec = pwcp->ldbase_genovec;
  STD_ARRAY_REF(uint32_t, 4) ldbase_genocounts = pwcp->ldbase_genocounts;
  if (!(vidx % kPglVblockSize)) {
    // beginning of a variant block.  save raw fpos in header; LD compression
    // prohibited.

    // er, need to use a relative offset in the multithreaded case, absolute
    // position isn't known
    pwcp->vblock_fpos[vidx / kPglVblockSize] = pwcp->vblock_fpos_offset + S_CAST(uintptr_t, pwcp->fwrite_bufp - pwcp->fwrite_buf);
  } else if (difflist_len > sample_ctd64) {
    // do not use LD compression if there are at least this many differences.
    // tune this threshold in the future.
    const uint32_t ld_diff_threshold = difflist_viable? (difflist_len - sample_ctd64) : max_difflist_len;
    // number of changes between current genovec and LD reference is bounded
    // below by sum(genocounts[x] - ldbase_genocounts[x]) / 2
    const int32_t count02_limit = 2 * ld_diff_threshold - abs_i32(genocounts[1] - ldbase_genocounts[1]) + abs_i32(genocounts[3] - ldbase_genocounts[3]);
    if ((S_CAST(int32_t, abs_i32(genocounts[0] - ldbase_genocounts[0]) + abs_i32(genocounts[2] - ldbase_genocounts[2])) < count02_limit) || (S_CAST(int32_t, abs_i32(genocounts[0] - ldbase_genocounts[2]) + abs_i32(genocounts[2] - ldbase_genocounts[0])) < count02_limit)) {
      uint32_t ld_diff_ct;
      uint32_t ld_inv_diff_ct;
      // okay, perform a brute-force diff
      // (could check LD vs. inverted LD separately?)
      if (pwcp->ldbase_common_geno < 4) {
        // unpack to ldbase_genovec
        PgrDifflistToGenovecUnsafe(pwcp->ldbase_raregeno, pwcp->ldbase_difflist_sample_ids, pwcp->ldbase_common_geno, sample_ct, pwcp->ldbase_difflist_len, ldbase_genovec);
        ZeroTrailingNyps(sample_ct, ldbase_genovec);
        pwcp->ldbase_common_geno = UINT32_MAX;
      }
      CountLdAndInvertedLdDiffs(ldbase_genovec, genovec, sample_ct, &ld_diff_ct, &ld_inv_diff_ct);
      if ((ld_diff_ct < ld_diff_threshold) || (ld_inv_diff_ct < ld_diff_threshold)) {
        const uintptr_t invert_before_compressing = (ld_inv_diff_ct < ld_diff_ct);
        *vrtype_ptr = 2 + invert_before_compressing;
        if (invert_before_compressing) {
          GenovecInvertCopyUnsafe(genovec, sample_ct, pwcp->genovec_invert_buf);
          ld_diff_ct = ld_inv_diff_ct;
        }
        return SaveLdDifflist(invert_before_compressing? pwcp->genovec_invert_buf : genovec, ldbase_genovec, 0, ld_diff_ct, pwcp);
      }
    }
  }
  const uint32_t genovec_word_ct = NypCtToWordCt(sample_ct);
  STD_ARRAY_COPY(genocounts, 4, ldbase_genocounts);
  pwcp->ldbase_common_geno = UINT32_MAX;
  if ((!difflist_viable) && (rare_2_geno_ct_sum < sample_ct / (2 * kPglMaxDifflistLenDivisor))) {
    *vrtype_ptr = 1;
    uint32_t larger_common_geno = second_most_common_geno;
    uint32_t smaller_common_geno = most_common_geno;
    if (most_common_geno > second_most_common_geno) {
      larger_common_geno = most_common_geno;
      smaller_common_geno = second_most_common_geno;
    }
    const uint32_t vrec_len = SaveOnebit(genovec, larger_common_geno + (smaller_common_geno * 3), rare_2_geno_ct_sum, pwcp);
    memcpy(ldbase_genovec, genovec, genovec_word_ct * sizeof(intptr_t));
    return vrec_len;
  }
  memcpy(ldbase_genovec, genovec, genovec_word_ct * sizeof(intptr_t));
  if (difflist_viable) {
#ifdef FUTURE_ENCODER
    if ((!difflist_len) && (!most_common_geno)) {
      *vrtype_ptr = 5;
      return 0;
    }
#endif
    *vrtype_ptr = 4 + most_common_geno;
    return SaveLdDifflist(genovec, nullptr, most_common_geno, difflist_len, pwcp);
  }
  *vrtype_ptr = 0;
  const uint32_t vrec_len = NypCtToByteCt(sample_ct);
  pwcp->fwrite_bufp = memcpyua(pwcp->fwrite_bufp, genovec, vrec_len);
  return vrec_len;
}

void PwcAppendBiallelicGenovec(const uintptr_t* __restrict genovec, PgenWriterCommon* pwcp) {
  const uint32_t vidx = pwcp->vidx;
  unsigned char vrtype;
  const uint32_t vrec_len = PwcAppendBiallelicGenovecMain(genovec, vidx, pwcp, nullptr, nullptr, &vrtype);
  const uintptr_t vrec_len_byte_ct = pwcp->vrec_len_byte_ct;
  pwcp->vidx += 1;
  SubU32Store(vrec_len, vrec_len_byte_ct, &(pwcp->vrec_len_buf[vidx * vrec_len_byte_ct]));
  // could have a single expression which branchlessly handles both cases, but
  // doubt that's worthwhile
  if (!pwcp->phase_dosage_gflags) {
    pwcp->vrtype_buf[vidx / kBitsPerWordD4] |= S_CAST(uintptr_t, vrtype) << (4 * (vidx % kBitsPerWordD4));
  } else {
    R_CAST(unsigned char*, pwcp->vrtype_buf)[vidx] = vrtype;
  }
}

BoolErr SpgwFlush(STPgenWriter* spgwp) {
  PgenWriterCommon* pwcp = GetPwcp(spgwp);
  if (pwcp->fwrite_bufp >= &(pwcp->fwrite_buf[kPglFwriteBlockSize])) {
    const uintptr_t cur_byte_ct = pwcp->fwrite_bufp - pwcp->fwrite_buf;
    FILE** pgen_outfilep = GetPgenOutfilep(spgwp);
    if (unlikely(fwrite_checked(pwcp->fwrite_buf, cur_byte_ct, *pgen_outfilep))) {
      return 1;
    }
    pwcp->vblock_fpos_offset += cur_byte_ct;
    pwcp->fwrite_bufp = pwcp->fwrite_buf;
  }
  return 0;
}


uint32_t SaveLdTwoListDelta(const uintptr_t* __restrict difflist_raregeno, const uint32_t* __restrict difflist_sample_ids, uint32_t ld_diff_ct, PgenWriterCommon* pwcp) {
  // assumes ldbase_difflist_sample_ids[ldbase_difflist_len] == sample_ct, and
  // difflist_sample_ids[ldbase_difflist_len] == sample_ct.
  // assumes biallelic data.

  // similar to SaveLdDifflist() and, to a lesser degree,
  // ParseLdAndMergeDifflistSubset()
  unsigned char* fwrite_bufp = pwcp->fwrite_bufp;
  if (!ld_diff_ct) {
    *fwrite_bufp = 0;
    pwcp->fwrite_bufp = &(fwrite_bufp[1]);
    return 1;
  }
  unsigned char* fwrite_bufp_start = fwrite_bufp;
  fwrite_bufp = Vint32Append(ld_diff_ct, fwrite_bufp);
  const uint32_t sample_id_byte_ct = BytesToRepresentNzU32(pwcp->sample_ct);
  const uint32_t group_ct = DivUp(ld_diff_ct, kPglDifflistGroupSize);
  const uint32_t ldbase_common_geno = pwcp->ldbase_common_geno;
  assert(ldbase_common_geno < 4);
  const uintptr_t* __restrict ldbase_raregeno = pwcp->ldbase_raregeno;
  const uint32_t* __restrict ldbase_sample_ids = pwcp->ldbase_difflist_sample_ids;
  unsigned char* group_first_sample_ids_iter = fwrite_bufp;
  unsigned char* extra_byte_cts_iter = &(fwrite_bufp[group_ct * sample_id_byte_ct]);
#ifdef NO_UNALIGNED
#  error "Unaligned accesses in SaveLdTwoListDelta()."
#endif
  uintptr_t* raregeno_write_iter = R_CAST(uintptr_t*, &(extra_byte_cts_iter[group_ct - 1]));
  fwrite_bufp = &(extra_byte_cts_iter[group_ct + (ld_diff_ct - 1) / 4]);
  unsigned char* last_group_vint_start = fwrite_bufp;
  uintptr_t ldbase_raregeno_word = 0;
  uintptr_t difflist_raregeno_word = 0;
  uintptr_t raregeno_write_word = 0;
  uint32_t last_sample_idx = 0;

  uint32_t next_ldbase_sample_idx = ldbase_sample_ids[0];
  uint32_t next_difflist_sample_idx = difflist_sample_ids[0];
  uint32_t ldbase_idx = 0;
  uint32_t difflist_idx = 0;
  uint32_t diff_written_ct = 0;
  while (diff_written_ct < ld_diff_ct) {
    uintptr_t cur_geno;
    uint32_t new_sample_idx;
    if (next_ldbase_sample_idx <= next_difflist_sample_idx) {
      ldbase_raregeno_word >>= 2;
      if (!(ldbase_idx % kBitsPerWordD2)) {
        ldbase_raregeno_word = ldbase_raregeno[ldbase_idx / kBitsPerWordD2];
      }
      ++ldbase_idx;
    }
    if (next_difflist_sample_idx <= next_ldbase_sample_idx) {
      difflist_raregeno_word >>= 2;
      if (!(difflist_idx % kBitsPerWordD2)) {
        difflist_raregeno_word = difflist_raregeno[difflist_idx / kBitsPerWordD2];
      }
      new_sample_idx = next_difflist_sample_idx;
      ++difflist_idx;
      cur_geno = difflist_raregeno_word & 3;
      next_difflist_sample_idx = difflist_sample_ids[difflist_idx];
      if (next_ldbase_sample_idx == new_sample_idx) {
        next_ldbase_sample_idx = ldbase_sample_ids[ldbase_idx];
        if (cur_geno == (ldbase_raregeno_word & 3)) {
          continue;
        }
      }
    } else {
      cur_geno = ldbase_common_geno;
      new_sample_idx = next_ldbase_sample_idx;
      next_ldbase_sample_idx = ldbase_sample_ids[ldbase_idx];
    }
    raregeno_write_word |= cur_geno << (2 * (diff_written_ct % kBitsPerWordD2));
    if (!(diff_written_ct % kPglDifflistGroupSize)) {
      SubU32StoreMov(new_sample_idx, sample_id_byte_ct, &group_first_sample_ids_iter);
      if (diff_written_ct) {
        *extra_byte_cts_iter++ = S_CAST(uintptr_t, fwrite_bufp - last_group_vint_start) - (kPglDifflistGroupSize - 1);
      }
      last_group_vint_start = fwrite_bufp;
    } else {
      fwrite_bufp = Vint32Append(new_sample_idx - last_sample_idx, fwrite_bufp);
    }
    last_sample_idx = new_sample_idx;
    ++diff_written_ct;
    if (!(diff_written_ct % kBitsPerWordD2)) {
      *raregeno_write_iter++ = raregeno_write_word;
      raregeno_write_word = 0;
    }
  }
  if (diff_written_ct % kBitsPerWordD2) {
    SubwordStore(raregeno_write_word, 1 + (((ld_diff_ct - 1) / 4) % kBytesPerWord), raregeno_write_iter);
  }
  pwcp->fwrite_bufp = fwrite_bufp;
  return fwrite_bufp - fwrite_bufp_start;
}

uint32_t SaveLdInputList(PgenWriterCommon* pwcp) {
  // simply "copies" ldbase_{raregeno,difflist_sample_ids,difflist_len} to the
  // write buffer.
  unsigned char* fwrite_bufp = pwcp->fwrite_bufp;
  const uint32_t difflist_len = pwcp->ldbase_difflist_len;
  if (!difflist_len) {
    *fwrite_bufp = 0;
    pwcp->fwrite_bufp = &(fwrite_bufp[1]);
    return 1;
  }
  unsigned char* fwrite_bufp_start = fwrite_bufp;
  fwrite_bufp = Vint32Append(difflist_len, fwrite_bufp);
  const uint32_t sample_id_byte_ct = BytesToRepresentNzU32(pwcp->sample_ct);
  const uint32_t group_ct = DivUp(difflist_len, kPglDifflistGroupSize);
  const uint32_t* __restrict difflist_sample_ids = pwcp->ldbase_difflist_sample_ids;
  unsigned char* group_first_sample_ids_iter = fwrite_bufp;
  unsigned char* extra_byte_cts_iter = &(fwrite_bufp[group_ct * sample_id_byte_ct]);
  fwrite_bufp = memcpyua(&(extra_byte_cts_iter[group_ct - 1]), pwcp->ldbase_raregeno, NypCtToByteCt(difflist_len));
  unsigned char* last_group_vint_start = nullptr;
  uint32_t last_sample_idx = 0;
  for (uint32_t difflist_idx = 0; difflist_idx != difflist_len; ++difflist_idx) {
    const uint32_t new_sample_idx = difflist_sample_ids[difflist_idx];
    if (!(difflist_idx % kPglDifflistGroupSize)) {
      SubU32StoreMov(new_sample_idx, sample_id_byte_ct, &group_first_sample_ids_iter);
      if (difflist_idx) {
        *extra_byte_cts_iter++ = S_CAST(uintptr_t, fwrite_bufp - last_group_vint_start) - (kPglDifflistGroupSize - 1);
      }
      last_group_vint_start = fwrite_bufp;
    } else {
      // assert(new_sample_idx >= last_sample_idx + 1);
      fwrite_bufp = Vint32Append(new_sample_idx - last_sample_idx, fwrite_bufp);
    }
    last_sample_idx = new_sample_idx;
  }
  pwcp->fwrite_bufp = fwrite_bufp;
  return fwrite_bufp - fwrite_bufp_start;
}

uint32_t PwcAppendBiallelicDifflistLimitedMain(const uintptr_t* __restrict raregeno, const uint32_t* __restrict difflist_sample_ids, uint32_t vidx, uint32_t difflist_common_geno, uint32_t difflist_len, PgenWriterCommon* pwcp, unsigned char* vrtype_ptr) {
  const uint32_t sample_ct = pwcp->sample_ct;
  // caller's responsibility not to exceed this limit
  assert(difflist_len <= 2 * (sample_ct / kPglMaxDifflistLenDivisor));

  // trailing bits of raregeno must be zeroed out

  assert(difflist_common_geno < 4);
  assert((!(difflist_len % kBitsPerWordD2)) || (!(raregeno[difflist_len / kBitsPerWordD2] >> (2 * (difflist_len % kBitsPerWordD2)))));
  assert(difflist_sample_ids[difflist_len] == sample_ct);
  STD_ARRAY_DECL(uint32_t, 4, genocounts);
  GenoarrCountFreqsUnsafe(raregeno, difflist_len, genocounts);
  assert(!genocounts[difflist_common_geno]);
  genocounts[difflist_common_geno] = sample_ct - difflist_len;
  uint32_t second_most_common_geno = difflist_common_geno? 0 : 1;
  uint32_t second_largest_geno_ct = genocounts[second_most_common_geno];
  for (uint32_t cur_geno = second_most_common_geno + 1; cur_geno != 4; ++cur_geno) {
    if (cur_geno == difflist_common_geno) {
      continue;
    }
    const uint32_t cur_geno_ct = genocounts[cur_geno];
    if (cur_geno_ct > second_largest_geno_ct) {
      second_most_common_geno = cur_geno;
      second_largest_geno_ct = cur_geno_ct;
    }
  }
  const uint32_t rare_2_geno_ct_sum = difflist_len - second_largest_geno_ct;
  const uint32_t sample_ctd8 = sample_ct / 8;
  const uint32_t sample_ctd64 = sample_ct / 64;
  uint32_t max_difflist_len = sample_ctd8 - 2 * sample_ctd64 + rare_2_geno_ct_sum;
  if (max_difflist_len > sample_ctd8) {
    max_difflist_len = sample_ctd8;
  }
  const uint32_t difflist_viable = (difflist_common_geno != 1) && (difflist_len <= max_difflist_len);
  STD_ARRAY_REF(uint32_t, 4) ldbase_genocounts = pwcp->ldbase_genocounts;
  if (!(vidx % kPglVblockSize)) {
    pwcp->vblock_fpos[vidx / kPglVblockSize] = pwcp->vblock_fpos_offset + S_CAST(uintptr_t, pwcp->fwrite_bufp - pwcp->fwrite_buf);
  } else if (difflist_len > sample_ctd64) {
    const uint32_t ld_diff_threshold = difflist_viable? (difflist_len - sample_ctd64) : max_difflist_len;
    // number of changes between current genovec and LD reference is bounded
    // below by sum(genocounts[x] - ldbase_genocounts[x]) / 2
    const int32_t count02_limit = 2 * ld_diff_threshold - abs_i32(genocounts[1] - ldbase_genocounts[1]) + abs_i32(genocounts[3] - ldbase_genocounts[3]);
    if ((S_CAST(int32_t, abs_i32(genocounts[0] - ldbase_genocounts[0]) + abs_i32(genocounts[2] - ldbase_genocounts[2])) < count02_limit) || (S_CAST(int32_t, abs_i32(genocounts[0] - ldbase_genocounts[2]) + abs_i32(genocounts[2] - ldbase_genocounts[0])) < count02_limit)) {
      uint32_t ld_diff_ct;
      uint32_t ld_inv_diff_ct;
      if (pwcp->ldbase_common_geno < 4) {
        pwcp->ldbase_difflist_sample_ids[pwcp->ldbase_difflist_len] = sample_ct;
        CountLdAndInvertedLdDiffsList(pwcp->ldbase_raregeno, pwcp->ldbase_difflist_sample_ids, raregeno, difflist_sample_ids, pwcp->ldbase_difflist_len, difflist_len, &ld_diff_ct, &ld_inv_diff_ct);
        const uint32_t difflist_common_geno_inv = (6 - difflist_common_geno) & 3;
        if (pwcp->ldbase_common_geno != difflist_common_geno) {
          ld_diff_ct = ld_diff_threshold;
        }
        if (pwcp->ldbase_common_geno != difflist_common_geno_inv) {
          ld_inv_diff_ct = ld_diff_threshold;
        }
        if ((ld_diff_ct < ld_diff_threshold) || (ld_inv_diff_ct < ld_diff_threshold)) {
          const uintptr_t invert_before_compressing = (ld_inv_diff_ct < ld_diff_ct);
          *vrtype_ptr = 2 + invert_before_compressing;
          if (invert_before_compressing) {
            GenovecInvertCopyUnsafe(raregeno, difflist_len, pwcp->genovec_invert_buf);
            // difflist_common_geno = difflist_common_geno_inv;
            ld_diff_ct = ld_inv_diff_ct;
          }
          return SaveLdTwoListDelta(invert_before_compressing? pwcp->genovec_invert_buf : raregeno, difflist_sample_ids, ld_diff_ct, pwcp);
        }
      } else {
        uintptr_t* __restrict genobuf = pwcp->genovec_invert_buf;
        PgrDifflistToGenovecUnsafe(raregeno, difflist_sample_ids, difflist_common_geno, sample_ct, difflist_len, genobuf);
        ZeroTrailingNyps(sample_ct, genobuf);
        CountLdAndInvertedLdDiffs(pwcp->ldbase_genovec, genobuf, sample_ct, &ld_diff_ct, &ld_inv_diff_ct);
        if ((ld_diff_ct < ld_diff_threshold) || (ld_inv_diff_ct < ld_diff_threshold)) {
          const uintptr_t invert_before_compressing = (ld_inv_diff_ct < ld_diff_ct);
          *vrtype_ptr = 2 + invert_before_compressing;
          if (invert_before_compressing) {
            GenovecInvertUnsafe(sample_ct, genobuf);
            ld_diff_ct = ld_inv_diff_ct;
          }
          return SaveLdDifflist(genobuf, pwcp->ldbase_genovec, 0, ld_diff_ct, pwcp);
        }
      }
    }
  }
  STD_ARRAY_COPY(genocounts, 4, ldbase_genocounts);
  if (difflist_viable) {
    *vrtype_ptr = 4 + difflist_common_geno;
    memcpy(pwcp->ldbase_raregeno, raregeno, NypCtToByteCt(difflist_len));
    memcpy(pwcp->ldbase_difflist_sample_ids, difflist_sample_ids, difflist_len * sizeof(int32_t));
    // memcpy(pwcp->ldbase_difflist_sample_ids, difflist_sample_ids, (difflist_len + 1) * sizeof(int32_t));
    pwcp->ldbase_common_geno = difflist_common_geno;
    pwcp->ldbase_difflist_len = difflist_len;
    return SaveLdInputList(pwcp);
  }
  pwcp->ldbase_common_geno = UINT32_MAX;
  const uint32_t use_onebit = (rare_2_geno_ct_sum < sample_ct / (2 * kPglMaxDifflistLenDivisor));
  uintptr_t* genobuf = use_onebit? pwcp->genovec_invert_buf : pwcp->ldbase_genovec;
  PgrDifflistToGenovecUnsafe(raregeno, difflist_sample_ids, difflist_common_geno, sample_ct, difflist_len, genobuf);
  ZeroTrailingNyps(sample_ct, genobuf);
  *vrtype_ptr = use_onebit;
  if (use_onebit) {
    uint32_t larger_common_geno = second_most_common_geno;
    uint32_t smaller_common_geno = difflist_common_geno;
    if (difflist_common_geno > second_most_common_geno) {
      larger_common_geno = difflist_common_geno;
      smaller_common_geno = second_most_common_geno;
    }
    const uint32_t vrec_len = SaveOnebit(genobuf, larger_common_geno + (smaller_common_geno * 3), rare_2_geno_ct_sum, pwcp);
    memcpy(pwcp->ldbase_genovec, genobuf, NypCtToWordCt(sample_ct) * sizeof(uintptr_t));
    return vrec_len;
  }
  const uint32_t vrec_len = NypCtToByteCt(sample_ct);
  pwcp->fwrite_bufp = memcpyua(pwcp->fwrite_bufp, genobuf, vrec_len);
  return vrec_len;
}

void PwcAppendBiallelicDifflistLimited(const uintptr_t* __restrict raregeno, const uint32_t* __restrict difflist_sample_ids, uint32_t difflist_common_geno, uint32_t difflist_len, PgenWriterCommon* pwcp) {
  const uint32_t vidx = pwcp->vidx;
  unsigned char vrtype;
  const uint32_t vrec_len = PwcAppendBiallelicDifflistLimitedMain(raregeno, difflist_sample_ids, vidx, difflist_common_geno, difflist_len, pwcp, &vrtype);
  const uintptr_t vrec_len_byte_ct = pwcp->vrec_len_byte_ct;
  pwcp->vidx += 1;
  SubU32Store(vrec_len, vrec_len_byte_ct, &(pwcp->vrec_len_buf[vidx * vrec_len_byte_ct]));
  if (!pwcp->phase_dosage_gflags) {
    pwcp->vrtype_buf[vidx / kBitsPerWordD4] |= S_CAST(uintptr_t, vrtype) << (4 * (vidx % kBitsPerWordD4));
  } else {
    R_CAST(unsigned char*, pwcp->vrtype_buf)[vidx] = vrtype;
  }
}


static inline BoolErr CheckedVrecLenIncr(uintptr_t incr, uint32_t* vrec_len_ptr) {
  // maybe track vrec_left instead of vrec_len...
#ifdef __LP64__
  const uintptr_t new_vrec_len = *vrec_len_ptr + incr;
  if (unlikely(new_vrec_len > kPglMaxBytesPerVariant)) {
    return 1;
  }
  *vrec_len_ptr = new_vrec_len;
  return 0;
#else
  *vrec_len_ptr += incr;
  return 0;
#endif
}

// ok if delta_bitarr trailing bits set
BoolErr PwcAppendDeltalist(const uintptr_t* delta_bitarr, uint32_t deltalist_len, PgenWriterCommon* pwcp, uint32_t* vrec_len_ptr) {
  assert(deltalist_len);
  unsigned char* fwrite_bufp = pwcp->fwrite_bufp;
  unsigned char* fwrite_bufp_start = fwrite_bufp;
  fwrite_bufp = Vint32Append(deltalist_len, fwrite_bufp);
  const uint32_t sample_id_byte_ct = BytesToRepresentNzU32(pwcp->sample_ct);
  const uint32_t group_ct = DivUp(deltalist_len, kPglDifflistGroupSize);
  unsigned char* group_first_sample_ids_iter = fwrite_bufp;
  if (unlikely(CheckedVrecLenIncr(S_CAST(uintptr_t, fwrite_bufp - fwrite_bufp_start) + group_ct * sample_id_byte_ct + (group_ct - 1), vrec_len_ptr))) {
    return 1;
  }
  unsigned char* extra_byte_cts_iter = &(fwrite_bufp[group_ct * sample_id_byte_ct]);
  fwrite_bufp_start = &(extra_byte_cts_iter[group_ct - 1]);
  fwrite_bufp = fwrite_bufp_start;
#ifdef __LP64__
  const intptr_t vrec_left = kPglMaxBytesPerVariant - (*vrec_len_ptr);
#endif
  unsigned char* last_group_vint_start = nullptr;
  uintptr_t last_sample_idx = 0;
  uintptr_t new_sample_idx_base = 0;
  uintptr_t delta_bits = delta_bitarr[0];
  for (uint32_t deltalist_idx = 0; deltalist_idx != deltalist_len; ++deltalist_idx) {
    const uintptr_t new_sample_idx = BitIter1(delta_bitarr, &new_sample_idx_base, &delta_bits);
    if (!(deltalist_idx % kPglDifflistGroupSize)) {
      SubU32StoreMov(new_sample_idx, sample_id_byte_ct, &group_first_sample_ids_iter);
      if (deltalist_idx) {
        *extra_byte_cts_iter++ = S_CAST(uintptr_t, fwrite_bufp - last_group_vint_start) - (kPglDifflistGroupSize - 1);
      }
#ifdef __LP64__
      // Note that this inequality goes in the opposite direction from the
      // usual PtrCheck.
      if (unlikely(S_CAST(intptr_t, fwrite_bufp - fwrite_bufp_start) > vrec_left)) {
        return 1;
      }
#endif
      last_group_vint_start = fwrite_bufp;
    } else {
      assert(new_sample_idx >= last_sample_idx + 1);
      fwrite_bufp = Vint32Append(new_sample_idx - last_sample_idx, fwrite_bufp);
    }
    last_sample_idx = new_sample_idx;
  }
  pwcp->fwrite_bufp = fwrite_bufp;
  return CheckedVrecLenIncr(fwrite_bufp - fwrite_bufp_start, vrec_len_ptr);
}

static_assert(sizeof(AlleleCode) == 1, "PwcAppendMultiallelicMain() needs to be updated.");
BoolErr PwcAppendMultiallelicMain(const uintptr_t* __restrict genovec, const uintptr_t* __restrict patch_01_set, const AlleleCode* __restrict patch_01_vals, const uintptr_t* __restrict patch_10_set, const AlleleCode* __restrict patch_10_vals, uint32_t allele_ct, uint32_t patch_01_ct, uint32_t patch_10_ct, uint32_t vidx, PgenWriterCommon* pwcp, const uintptr_t** genovec_hets_ptr, uint32_t* het_ct_ptr, unsigned char* vrtype_ptr, uint32_t* vrec_len_ptr) {
  uint32_t genovec_het_ct;
  uint32_t genovec_altxy_ct;
  uint32_t vrec_len = PwcAppendBiallelicGenovecMain(genovec, vidx, pwcp, &genovec_het_ct, &genovec_altxy_ct, vrtype_ptr);
  if ((!patch_01_ct) && (!patch_10_ct)) {
    if (genovec_hets_ptr) {
      *genovec_hets_ptr = genovec;
      *het_ct_ptr = genovec_het_ct;
    }
    *vrec_len_ptr = vrec_len;
    return 0;
  }
  assert(allele_ct > 2);
  unsigned char* format_byte_ptr = pwcp->fwrite_bufp;
  ++vrec_len;
  pwcp->fwrite_bufp += 1;
  unsigned char format_byte = 0;
  if (!patch_01_ct) {
    format_byte = 15;
  } else {
    // We could use a different bitarray vs. deltalist threshold, since the
    // usual threshold was chosen under the assumption that
    //   (maximum allowed deltalist index + 1) == bitarray length.
    // However, the subsetting-speed advantage of the deltalist representation
    // can be expected to make up for its larger size, so we may actually want
    // to increase instead of decrease the current 1/9 fraction.  For now, I'll
    // punt on this decision.
    const uint32_t max_deltalist_entry_ct = genovec_het_ct / kPglMaxDeltalistLenDivisor;
    if (patch_01_ct < max_deltalist_entry_ct) {
      // impossible for this to run out of space
      PwcAppendDeltalist(patch_01_set, patch_01_ct, pwcp, &vrec_len);
      format_byte = 1;
    } else {
      const uint32_t genovec_het_ctb = DivUp(genovec_het_ct, CHAR_BIT);
#ifdef NO_UNALIGNED
#  error "Unaligned accesses in PwcAppendMultiallelicMain()."
#endif
      CopyGenomatchSubset(patch_01_set, genovec, kMask5555, 0, genovec_het_ct, R_CAST(uintptr_t*, pwcp->fwrite_bufp));
      pwcp->fwrite_bufp = &(pwcp->fwrite_bufp[genovec_het_ctb]);
      vrec_len += genovec_het_ctb;
    }
    if (allele_ct > 3) {
      if (allele_ct <= 18) {
        const AlleleCode* patch_01_vals_iter = patch_01_vals;

        // er, no, we want this to be a different pointer type depending on the
        // allele count and build type.  but make those optimizations later,
        // after a simple unified version of this code has been tested.
        uintptr_t* write_alias = R_CAST(uintptr_t*, pwcp->fwrite_bufp);
        uint32_t written_byte_ct;
        if (allele_ct == 4) {
          // 1 bit per entry, clear = ref/alt2, set = ref/alt3
#ifdef __LP64__
          const uint32_t fullvec_ct = patch_01_ct / kBytesPerVec;
          Vec8thUint* write_aliasv = R_CAST(Vec8thUint*, write_alias);
          if (fullvec_ct) {
            // parallel-equality-to-3 check, movemask
            // (+125 or <<7 followed by movemask also works)
            for (uint32_t vec_idx = 0; vec_idx != fullvec_ct; ++vec_idx) {
              VecUc cur_vec = vecuc_loadu(patch_01_vals_iter);
              patch_01_vals_iter = &(patch_01_vals_iter[kBytesPerVec]);
              write_aliasv[vec_idx] = vecuc_movemask(R_CAST(VecUc, vecw_slli(R_CAST(VecW, cur_vec), 7)));
            }
          }
          const uint32_t remainder = patch_01_ct % kBytesPerVec;
          if (remainder) {
            const uintptr_t* cur_read_iter = R_CAST(const uintptr_t*, patch_01_vals_iter);
            unsigned char* write_uc = R_CAST(unsigned char*, &(write_aliasv[fullvec_ct]));
            const uint32_t word_ct_m1 = (remainder - 1) / kBytesPerWord;
            for (uint32_t widx = 0; ; ++widx) {
              uintptr_t cur_word;
              if (widx >= word_ct_m1) {
                if (widx > word_ct_m1) {
                  break;
                }
                cur_word = SubwordLoad(cur_read_iter, ModNz(remainder, kBytesPerWord));
              } else {
                cur_word = *cur_read_iter++;
              }
#  ifdef USE_AVX2
              cur_word = _pext_u64(cur_word, kMask0101);
#  else
              // extract the low bit from each byte
              cur_word = cur_word & kMask0101;
              cur_word = (cur_word * 0x102040810204080LLU) >> 56;
#  endif
              write_uc[widx] = cur_word;
            }
          }
#else  // !__LP64__
          const uint32_t word_ct_m1 = (patch_01_ct - 1) / kBitsPerWord;
          uint32_t loop_len = kBitsPerWord;
          for (uint32_t widx = 0; ; ++widx) {
            if (widx >= word_ct_m1) {
              if (widx > word_ct_m1) {
                break;
              }
              loop_len = ModNz(patch_01_ct, kBitsPerWord);
            }
            uintptr_t cur_word = 0;
            for (uint32_t uii = loop_len - 1; ; --uii) {
              cur_word |= patch_01_vals_iter[uii] & 1;
              if (!uii) {
                break;
              }
              cur_word <<= 1;
            }
            write_alias[widx] = cur_word;
            patch_01_vals_iter = &(patch_01_vals_iter[loop_len]);
          }
#endif  // !__LP64__
          written_byte_ct = DivUp(patch_01_ct, 8);
        } else if (allele_ct <= 6) {
          // 2 bits per entry
          const uint32_t word_ct_m1 = (patch_01_ct - 1) / kBitsPerWordD2;
          uint32_t loop_len = kBitsPerWordD2;
          for (uint32_t widx = 0; ; ++widx) {
            if (widx >= word_ct_m1) {
              if (widx > word_ct_m1) {
                break;
              }
              loop_len = ModNz(patch_01_ct, kBitsPerWordD2);
            }
            uintptr_t cur_word = 0;
            for (uint32_t uii = loop_len - 1; ; --uii) {
              // todo:
              // 1. parallel-subtract (now all bytes are in 0..3)
              // 2. bitwise or with right-shift-6
              // 3. bitwise or with right-shift-12
              // 4. _mm{256}_shuffle_epi8() gather from bytes 0, 4, 8, ...
              // SSE4.2:
              // 5. store result of _mm_extract_epi32(, 0)
              // AVX2:
              // 5. store results of _mm256_extract_epi32(, 0) and (, 4)
              cur_word |= patch_01_vals_iter[uii] - 2;
              if (!uii) {
                break;
              }
              cur_word <<= 2;
            }
            write_alias[widx] = cur_word;
            patch_01_vals_iter = &(patch_01_vals_iter[loop_len]);
          }
          written_byte_ct = DivUp(patch_01_ct, 4);
        } else {
          // 4 bits per entry
          const uint32_t word_ct_m1 = (patch_01_ct - 1) / kBitsPerWordD4;
          uint32_t loop_len = kBitsPerWordD4;
          for (uint32_t widx = 0; ; ++widx) {
            if (widx >= word_ct_m1) {
              if (widx > word_ct_m1) {
                break;
              }
              loop_len = ModNz(patch_01_ct, kBitsPerWordD4);
            }
            uintptr_t cur_word = 0;
            for (uint32_t uii = loop_len - 1; ; --uii) {
              // todo: replace this with parallel-subtract, followed by bitwise
              // or with right-shift-4, followed by _mm{256}_shuffle_epi8()
              // gather from even bytes, followed by _mm{256}_extract_epi64()
              // (or _mm256_permute4x64_epi64() followed by
              // _mm256_extracti128_si256()?  test this.)
              // Or in non-SSE4.2 case, use end of PackWordToHalfword logic.
              cur_word |= patch_01_vals_iter[uii] - 2;
              if (!uii) {
                break;
              }
              cur_word <<= 4;
            }
            write_alias[widx] = cur_word;
            patch_01_vals_iter = &(patch_01_vals_iter[loop_len]);
          }
          written_byte_ct = DivUp(patch_01_ct, 2);
        }
        pwcp->fwrite_bufp = &(pwcp->fwrite_bufp[written_byte_ct]);
        vrec_len += written_byte_ct;
      } else {
        // 1 byte per entry
        // need 2-byte and 3-byte cases here for larger AlleleCode, those also
        // need to check for vrec_len overflow
        unsigned char* payload = pwcp->fwrite_bufp;
        for (uint32_t patch_idx = 0; patch_idx != patch_01_ct; ++patch_idx) {
          // todo: check if the compiler vectorizes this
          payload[patch_idx] = patch_01_vals[patch_idx] - 2;
        }
        pwcp->fwrite_bufp = &(pwcp->fwrite_bufp[patch_01_ct]);
        vrec_len += patch_01_ct;
      }
    }
  }
  if (!patch_10_ct) {
    format_byte |= 240;
    if (genovec_hets_ptr) {
      *genovec_hets_ptr = genovec;
      *het_ct_ptr = genovec_het_ct;
    }
  } else {
    const uint32_t max_deltalist_entry_ct = genovec_altxy_ct / kPglMaxDeltalistLenDivisor;
    if (patch_10_ct < max_deltalist_entry_ct) {
      // can't actually fail for now, but cannot ignore overflow check if
      // sizeof(AlleleCode) > 1
      if (unlikely(PwcAppendDeltalist(patch_10_set, patch_10_ct, pwcp, &vrec_len))) {
        return 1;
      }
      format_byte |= 16;
    } else {
      const uintptr_t genovec_altxy_ctb = DivUp(genovec_altxy_ct, CHAR_BIT);
      if (unlikely(CheckedVrecLenIncr(genovec_altxy_ctb, &vrec_len))) {
        return 1;
      }
      CopyGenomatchSubset(patch_10_set, genovec, kMaskAAAA, 0, genovec_altxy_ct, R_CAST(uintptr_t*, pwcp->fwrite_bufp));
      pwcp->fwrite_bufp = &(pwcp->fwrite_bufp[genovec_altxy_ctb]);
    }
    if (allele_ct <= 17) {
      const AlleleCode* patch_10_vals_iter = patch_10_vals;
      uintptr_t* write_alias = R_CAST(uintptr_t*, pwcp->fwrite_bufp);
      uint32_t bytes_to_write;
      if (allele_ct == 3) {
        // 1 bit per entry, clear = alt1/alt2, set = alt2/alt2
        bytes_to_write = DivUp(patch_10_ct, 8);
        if (unlikely(CheckedVrecLenIncr(bytes_to_write, &vrec_len))) {
          return 1;
        }
#ifdef __LP64__
        const uint32_t fullvec_ct = patch_10_ct / (kBytesPerVec / 2);
        Vec16thUint* write_aliasv = R_CAST(Vec16thUint*, write_alias);
        if (fullvec_ct) {
#  if defined(USE_SSE42) && !defined(USE_AVX2)
          // SSE4.2: _mm_shuffle_epi8() to gather even bytes, parallel
          //         equality-to-2 check, movemask
          // (+126 or <<6 followed by movemask also works)
          const VecUc gather_even = vecuc_setr8(0, 2, 4, 6, 8, 10, 12, 14,
                                                -1, -1, -1, -1, -1, -1, -1, -1);
          for (uint32_t vec_idx = 0; vec_idx != fullvec_ct; ++vec_idx) {
            const VecUc cur_vec = vecuc_loadu(patch_10_vals_iter);
            patch_10_vals_iter = &(patch_10_vals_iter[kBytesPerVec]);
            const VecUc even_vec = vecuc_shuffle8(cur_vec, gather_even);
            const uint32_t cur_u32 = vecuc_movemask(R_CAST(VecUc, vecw_slli(R_CAST(VecW, even_vec), 6)));
            write_aliasv[vec_idx] = cur_u32;
          }
#  else
          // SSE2/AVX2: parallel equality-to-{2, 2} check, movemask,
          //            PackWordToHalfword
          // (uint16_t +126 followed by movemask also works, and avoids an
          // extra masking step in SSE2 case)
          // (<<6 doesn't work here)
          const VecU16 all126 = vecu16_set1(126);
          for (uint32_t vec_idx = 0; vec_idx != fullvec_ct; ++vec_idx) {
            VecU16 cur_vec = vecu16_loadu(patch_10_vals_iter);
            patch_10_vals_iter = &(patch_10_vals_iter[kBytesPerVec]);
            const uint32_t cur_u32 = vecu16_movemask(cur_vec + all126);
            write_aliasv[vec_idx] = PackVec8thUintTo16th(cur_u32);
          }
#  endif
        }
        const uint32_t remainder = patch_10_ct % (kBytesPerVec / 2);
        if (remainder) {
          uint32_t cur_u32 = 0;
          for (uint32_t uii = remainder - 1; ; --uii) {
            cur_u32 += patch_10_vals_iter[2 * uii];
            if (!uii) {
              break;
            }
            cur_u32 <<= 1;
          }
          // this effectively subtracts 1 from each entry
          cur_u32 += 1 - (1U << remainder);
          write_aliasv[fullvec_ct] = cur_u32;
          patch_10_vals_iter = &(patch_10_vals_iter[remainder * 2]);
        }
#else  // !__LP64__
        const uint32_t word_ct_m1 = (patch_10_ct - 1) / kBitsPerWord;
        uint32_t loop_len = kBitsPerWord;
        for (uint32_t widx = 0; ; ++widx) {
          if (widx >= word_ct_m1) {
            if (widx > word_ct_m1) {
              break;
            }
            loop_len = ModNz(patch_10_ct, kBitsPerWord);
          }
          uintptr_t cur_word = 0;
          for (uint32_t uii = loop_len - 1; ; --uii) {
            cur_word |= patch_10_vals_iter[2 * uii] - 1;
            if (!uii) {
              break;
            }
            cur_word <<= 1;
          }
          write_alias[widx] = cur_word;
          patch_10_vals_iter = &(patch_10_vals_iter[loop_len * 2]);
        }
#endif  // !__LP64__
      } else if (allele_ct <= 5) {
        // 2 bits per half-entry
        bytes_to_write = DivUp(patch_10_ct, 2);
        if (unlikely(CheckedVrecLenIncr(bytes_to_write, &vrec_len))) {
          return 1;
        }
        // this may be worse than simple loop, check this later
        const uint32_t word_ct_m1 = (patch_10_ct - 1) / kBitsPerWordD4;
        uint32_t loop_len = kBitsPerWordD2;
        for (uint32_t widx = 0; ; ++widx) {
          if (widx >= word_ct_m1) {
            if (widx > word_ct_m1) {
              break;
            }
            loop_len = 2 * ModNz(patch_10_ct, kBitsPerWordD4);
          }
          uintptr_t cur_word = 0;
          for (uint32_t uii = loop_len - 1; ; --uii) {
            // todo: same as 2-bit patch_01 case, except starting with
            // parallel-subtract 1 instead of 2
            cur_word |= patch_10_vals_iter[uii] - 1;
            if (!uii) {
              break;
            }
            cur_word <<= 2;
          }
          write_alias[widx] = cur_word;
          patch_10_vals_iter = &(patch_10_vals_iter[loop_len]);
        }
      } else {
        // 4 bits per half-entry
        bytes_to_write = patch_10_ct;
        if (unlikely(CheckedVrecLenIncr(bytes_to_write, &vrec_len))) {
          return 1;
        }
        const uint32_t word_ct_m1 = (patch_10_ct - 1) / kBytesPerWord;
        uint32_t loop_len = kBitsPerWordD4;
        for (uint32_t widx = 0; ; ++widx) {
          if (widx >= word_ct_m1) {
            if (widx > word_ct_m1) {
              break;
            }
            loop_len = 2 * ModNz(patch_10_ct, kBytesPerWord);
          }
          uintptr_t cur_word = 0;
          for (uint32_t uii = loop_len - 1; ; --uii) {
            cur_word |= patch_10_vals_iter[uii] - 1;
            if (!uii) {
              break;
            }
            cur_word <<= 4;
          }
          write_alias[widx] = cur_word;
          patch_10_vals_iter = &(patch_10_vals_iter[loop_len]);
        }
      }
      pwcp->fwrite_bufp = &(pwcp->fwrite_bufp[bytes_to_write]);
    } else {
      // 1 byte per half-entry
      // need 2- and 3-byte cases for larger AlleleCode
      unsigned char* payload = pwcp->fwrite_bufp;
      if (unlikely(CheckedVrecLenIncr(patch_10_ct * (2 * k1LU), &vrec_len))) {
        return 1;
      }
      const uint32_t patch_10_ct_x2 = patch_10_ct * 2;
      // hopefully the compiler automatically vectorizes this when
      // sizeof(AlleleCode) == 1
      for (uint32_t uii = 0; uii != patch_10_ct_x2; ++uii) {
        payload[uii] = patch_10_vals[uii] - 1;
      }
      pwcp->fwrite_bufp = &(pwcp->fwrite_bufp[patch_10_ct * 2]);
    }
    if (genovec_hets_ptr) {
      // This is somewhat redundant with the earlier steps, but let's get
      // everything working before doing more optimization.
      const uint32_t sample_ct = pwcp->sample_ct;
      const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
      const Halfword* patch_10_set_alias = R_CAST(const Halfword*, patch_10_set);
      const AlleleCode* patch_10_vals_iter = patch_10_vals;
      const AlleleCode* patch_10_vals_end = &(patch_10_vals[patch_10_ct * 2]);
      uintptr_t* genovec_hets = nullptr;
      uint32_t het_ct = genovec_het_ct;
      // todo: try inverting this loop, at least in 64-bit case.  perform
      // movemask-based scan for heterozygous patch_10_vals entries, then
      // advance widx/detect_10 as needed.
      for (uint32_t widx = 0; widx != sample_ctl2; ++widx) {
        uint32_t patch_10_hw = patch_10_set_alias[widx];
        if (patch_10_hw) {
          do {
            const AlleleCode val1 = *patch_10_vals_iter++;
            const AlleleCode val2 = *patch_10_vals_iter++;
            const uintptr_t lowbit = patch_10_hw & (-patch_10_hw);
            if (val1 != val2) {
              if (!genovec_hets) {
                genovec_hets = pwcp->genovec_hets_buf;
                memcpy(genovec_hets, genovec, sample_ctl2 * sizeof(intptr_t));
              }
              // clear high bit, set low bit
              genovec_hets[widx] ^= lowbit * lowbit * 3;
              ++het_ct;
            }
            patch_10_hw ^= lowbit;
          } while (patch_10_hw);
          if (patch_10_vals_iter == patch_10_vals_end) {
            break;
          }
        }
      }
      *het_ct_ptr = het_ct;
      if (genovec_hets) {
        *genovec_hets_ptr = genovec_hets;
      } else {
        *genovec_hets_ptr = genovec;
      }
    }
  }
  *format_byte_ptr = format_byte;
  *vrtype_ptr |= 8;
  *vrec_len_ptr = vrec_len;
  return 0;
}

BoolErr PwcAppendMultiallelicSparse(const uintptr_t* __restrict genovec, const uintptr_t* __restrict patch_01_set, const AlleleCode* __restrict patch_01_vals, const uintptr_t* __restrict patch_10_set, const AlleleCode* __restrict patch_10_vals, uint32_t allele_ct, uint32_t patch_01_ct, uint32_t patch_10_ct, PgenWriterCommon* pwcp) {
  const uint32_t vidx = pwcp->vidx;
  unsigned char vrtype;
  uint32_t vrec_len;
  if (unlikely(PwcAppendMultiallelicMain(genovec, patch_01_set, patch_01_vals, patch_10_set, patch_10_vals, allele_ct, patch_01_ct, patch_10_ct, vidx, pwcp, nullptr, nullptr, &vrtype, &vrec_len))) {
    return 1;
  }
  pwcp->vidx += 1;
  const uintptr_t vrec_len_byte_ct = pwcp->vrec_len_byte_ct;
  SubU32Store(vrec_len, vrec_len_byte_ct, &(pwcp->vrec_len_buf[vidx * vrec_len_byte_ct]));
  if (!pwcp->phase_dosage_gflags) {
    pwcp->vrtype_buf[vidx / kBitsPerWordD4] |= S_CAST(uintptr_t, vrtype) << (4 * (vidx % kBitsPerWordD4));
  } else {
    R_CAST(unsigned char*, pwcp->vrtype_buf)[vidx] = vrtype;
  }
  return 0;
}

void PglMultiallelicDenseToSparse(const AlleleCode* __restrict wide_codes, uint32_t sample_ct, uintptr_t* __restrict genoarr, uintptr_t* __restrict patch_01_set, AlleleCode* __restrict patch_01_vals, uintptr_t* __restrict patch_10_set, AlleleCode* __restrict patch_10_vals, uint32_t* __restrict patch_01_ct_ptr, uint32_t* __restrict patch_10_ct_ptr) {
  const uint32_t word_ct_m1 = (sample_ct - 1) / kBitsPerWordD2;
  const AlleleCode* wide_codes_iter = wide_codes;
  Halfword* patch_01_set_alias = R_CAST(Halfword*, patch_01_set);
  Halfword* patch_10_set_alias = R_CAST(Halfword*, patch_10_set);
  AlleleCode* patch_01_vals_iter = patch_01_vals;
  AlleleCode* patch_10_vals_iter = patch_10_vals;
  uint32_t loop_len = kBitsPerWordD2;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= word_ct_m1) {
      if (widx > word_ct_m1) {
        *patch_01_ct_ptr = patch_01_vals_iter - patch_01_vals;
        *patch_10_ct_ptr = S_CAST(uintptr_t, patch_10_vals_iter - patch_10_vals) / 2;
        return;
      }
      loop_len = ModNz(sample_ct, kBitsPerWordD2);
    }
    uintptr_t geno_word = 0;
    uint32_t patch_01_hw = 0;
    uint32_t patch_10_hw = 0;
    // ascending loop may be worse than descending here
    for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits) {
      // todo: try making this movemask-based in 64-bit case
      // parallel checks for equality-to-255, equality-to-0, inequality-to-1
      // rarealts are inequality-to-1 and (not 0 or 255); can use ctzu32() and
      // &= x - 1 to iterate
      // hom-ref iff double-equality to 0
      // het-ref iff low equality-to-0 bit set, high clear
      // missing iff equality-to-255 (better be double)
      // altx-alty otherwise
      const AlleleCode first_code = *wide_codes_iter++;
      const AlleleCode second_code = *wide_codes_iter++;
      uintptr_t cur_geno;
      if (first_code) {
        if (first_code == kMissingAlleleCode) {
          cur_geno = 3;
        } else {
          cur_geno = 2;
          if (second_code > 1) {
            patch_10_hw |= 1U << sample_idx_lowbits;
            *patch_10_vals_iter++ = first_code;
            *patch_10_vals_iter++ = second_code;
          }
        }
      } else if (!second_code) {
        continue;
      } else {
        cur_geno = 1;
        if (second_code > 1) {
          patch_01_hw |= 1U << sample_idx_lowbits;
          *patch_01_vals_iter++ = second_code;
        }
      }
      geno_word |= cur_geno << (2 * sample_idx_lowbits);
    }
    genoarr[widx] = geno_word;
    patch_01_set_alias[widx] = patch_01_hw;
    patch_10_set_alias[widx] = patch_10_hw;
  }
}

static_assert(sizeof(AlleleCode) == 1, "PglMultiallelicSparseToDense() needs to be updated.");
void PglMultiallelicSparseToDense(const uintptr_t* __restrict genoarr, const uintptr_t* __restrict patch_01_set, const AlleleCode* __restrict patch_01_vals, const uintptr_t* __restrict patch_10_set, const AlleleCode* __restrict patch_10_vals, const AlleleCode* __restrict remap, uint32_t sample_ct, uint32_t patch_01_ct, uint32_t patch_10_ct, uintptr_t* __restrict flipped, AlleleCode* __restrict wide_codes) {
  if (flipped) {
    ZeroWArr(BitCtToWordCt(sample_ct), flipped);
  }
  if ((!remap) || ((!remap[0]) && (remap[1] == 1))) {
    GenoarrLookup256x2bx4(genoarr, kHcToAlleleCodes, sample_ct, wide_codes);
    if (!remap) {
      if (patch_01_ct) {
        uintptr_t sample_idx_base = 0;
        uintptr_t cur_bits = patch_01_set[0];
        AlleleCode* wide_codes1 = &(wide_codes[1]);
        for (uint32_t uii = 0; uii != patch_01_ct; ++uii) {
          const uintptr_t sample_idx = BitIter1(patch_01_set, &sample_idx_base, &cur_bits);
          wide_codes1[2 * sample_idx] = patch_01_vals[uii];
        }
      }
      if (patch_10_ct) {
        uintptr_t sample_idx_base = 0;
        uintptr_t cur_bits = patch_10_set[0];
        const DoubleAlleleCode* patch_10_vals_alias = R_CAST(const DoubleAlleleCode*, patch_10_vals);
        DoubleAlleleCode* wide_codes_alias = R_CAST(DoubleAlleleCode*, wide_codes);
        for (uint32_t uii = 0; uii != patch_10_ct; ++uii) {
          const uintptr_t sample_idx = BitIter1(patch_10_set, &sample_idx_base, &cur_bits);
          wide_codes_alias[sample_idx] = patch_10_vals_alias[uii];
        }
      }
      return;
    }
  } else {
    // 0 -> (remap[0], remap[0])
    // 1 -> (remap[0], remap[1])
    // 2 -> (remap[1], remap[1])
    // 3 -> (255, 255)
    // 256x4 lookup table usually too expensive to reinitialize for every
    // single variant.  16x2 should be a good compromise?
    // todo: try Expand1bitTo8 strategy in SSE4.2+ case, see if that's fast
    // enough to remove 0/1 special case
    const uint32_t remap0 = remap[0];
    const uint32_t remap1 = remap[1];
    uint32_t table4[4];
    table4[0] = remap0 * 257;
    table4[1] = remap0 + remap1 * 256;
    table4[2] = remap1 * 257;
    table4[3] = 0xffff;
    if (remap1 < remap0) {
      table4[1] = remap1 + remap0 * 256;
      if (flipped) {
        // See GenoarrToMissingnessUnsafe().
        const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
        Halfword* flipped_alias = R_CAST(Halfword*, flipped);
        for (uint32_t widx = 0; widx != sample_ctl2; ++widx) {
          const uintptr_t cur_geno_word = genoarr[widx];
          flipped_alias[widx] = PackWordToHalfwordMask5555(cur_geno_word & (~(cur_geno_word >> 1)));
        }
      }
    }
    // could have an InitLookup16x2bx2 function for this
    uint32_t table16[16];
    for (uint32_t uii = 0; uii != 4; ++uii) {
      const uint32_t high_bits = table4[uii] << 16;
      uint32_t* table_row = &(table16[uii * 4]);
      for (uint32_t ujj = 0; ujj != 4; ++ujj) {
        table_row[ujj] = high_bits | table4[ujj];
      }
    }
    const uint32_t sample_ctd4 = sample_ct / 4;
    const unsigned char* genoarr_alias = R_CAST(const unsigned char*, genoarr);
    uint32_t* wide_codes_alias_iter = R_CAST(uint32_t*, wide_codes);
    for (uint32_t byte_idx = 0; byte_idx != sample_ctd4; ++byte_idx) {
      const uint32_t geno_byte = genoarr_alias[byte_idx];
      *wide_codes_alias_iter++ = table16[geno_byte & 15];
      *wide_codes_alias_iter++ = table16[geno_byte >> 4];
    }
    const uint32_t remainder = sample_ct % 4;
    if (remainder) {
      uint32_t geno_byte = genoarr_alias[sample_ctd4];
      if (remainder & 2) {
        *wide_codes_alias_iter++ = table16[geno_byte & 15];
        geno_byte = geno_byte >> 4;
      }
      if (remainder & 1) {
        uint16_t* last_code_pair_ptr = R_CAST(uint16_t*, wide_codes_alias_iter);
        // guess it's easy enough to defend against dirty trailing bits for now
        *last_code_pair_ptr = table4[geno_byte & 3];
      }
    }
  }
  if (patch_01_vals) {
    uintptr_t sample_idx_base = 0;
    uintptr_t cur_bits = patch_01_set[0];
    const AlleleCode remap0 = remap[0];
    const AlleleCode remap1 = remap[1];
    if ((!remap0) || ((!flipped) && (remap0 == 1) && (!remap1))) {
      // no flips possible
      AlleleCode* wide_codes1 = &(wide_codes[1]);
      for (uint32_t uii = 0; uii != patch_01_ct; ++uii) {
        const uintptr_t sample_idx = BitIter1(patch_01_set, &sample_idx_base, &cur_bits);
        wide_codes1[2 * sample_idx] = remap[patch_01_vals[uii]];
      }
    } else if (!flipped) {
      for (uint32_t uii = 0; uii != patch_01_ct; ++uii) {
        const uintptr_t sample_idx = BitIter1(patch_01_set, &sample_idx_base, &cur_bits);
        const AlleleCode ac = remap[patch_01_vals[uii]];
        if (ac > remap0) {
          wide_codes[2 * sample_idx + 1] = ac;
        } else {
          wide_codes[2 * sample_idx] = ac;
          wide_codes[2 * sample_idx + 1] = remap0;
        }
      }
    } else if (remap1 >= remap0) {
      for (uint32_t uii = 0; uii != patch_01_ct; ++uii) {
        const uintptr_t sample_idx = BitIter1(patch_01_set, &sample_idx_base, &cur_bits);
        const AlleleCode ac = remap[patch_01_vals[uii]];
        if (ac > remap0) {
          wide_codes[2 * sample_idx + 1] = ac;
        } else {
          wide_codes[2 * sample_idx] = ac;
          wide_codes[2 * sample_idx + 1] = remap0;
          SetBit(sample_idx, flipped);
        }
      }
    } else {
      for (uint32_t uii = 0; uii != patch_01_ct; ++uii) {
        const uintptr_t sample_idx = BitIter1(patch_01_set, &sample_idx_base, &cur_bits);
        const AlleleCode ac = remap[patch_01_vals[uii]];
        if (ac > remap0) {
          wide_codes[2 * sample_idx + 1] = ac;
          ClearBit(sample_idx, flipped);
        } else {
          wide_codes[2 * sample_idx] = ac;
          wide_codes[2 * sample_idx + 1] = remap0;
        }
      }
    }
  }
  if (patch_10_ct) {
    uintptr_t sample_idx_base = 0;
    uintptr_t cur_bits = patch_10_set[0];
    if (!flipped) {
      for (uint32_t uii = 0; uii != patch_10_ct; ++uii) {
        const uintptr_t sample_idx = BitIter1(patch_10_set, &sample_idx_base, &cur_bits);
        const AlleleCode ac0 = remap[patch_10_vals[2 * uii]];
        const AlleleCode ac1 = remap[patch_10_vals[2 * uii + 1]];
        if (ac0 <= ac1) {
          wide_codes[2 * sample_idx] = ac0;
          wide_codes[2 * sample_idx + 1] = ac1;
        } else {
          wide_codes[2 * sample_idx] = ac1;
          wide_codes[2 * sample_idx + 1] = ac0;
        }
      }
    } else {
      for (uint32_t uii = 0; uii != patch_10_ct; ++uii) {
        const uintptr_t sample_idx = BitIter1(patch_10_set, &sample_idx_base, &cur_bits);
        const AlleleCode ac0 = remap[patch_10_vals[2 * uii]];
        const AlleleCode ac1 = remap[patch_10_vals[2 * uii + 1]];
        if (ac0 <= ac1) {
          wide_codes[2 * sample_idx] = ac0;
          wide_codes[2 * sample_idx + 1] = ac1;
        } else {
          wide_codes[2 * sample_idx] = ac1;
          wide_codes[2 * sample_idx + 1] = ac0;
          SetBit(sample_idx, flipped);
        }
      }
    }
  }
}

// tolerates extraneous phaseinfo bits
BoolErr AppendHphase(const uintptr_t* __restrict genoarr_hets, const uintptr_t* __restrict phasepresent, const uintptr_t* __restrict phaseinfo, uint32_t het_ct, uint32_t phasepresent_ct, PgenWriterCommon* pwcp, unsigned char* vrtype_ptr, uint32_t* vrec_len_ptr) {
  assert(phasepresent_ct);
  const uint32_t sample_ct = pwcp->sample_ct;
  *vrtype_ptr += 16;
  const uint32_t het_ctp1_8 = 1 + (het_ct / CHAR_BIT);
#ifdef NO_UNALIGNED
#  error "Unaligned accesses in AppendHphase()."
#endif
  uintptr_t* fwrite_bufp_alias = R_CAST(uintptr_t*, pwcp->fwrite_bufp);
  uintptr_t phaseinfo_write_word = 0;
  uint32_t phaseinfo_write_idx_lowbits;
  unsigned char* fwrite_bufp_final;
  if (het_ct == phasepresent_ct) {
    // no need to write phasepresent; just write phaseinfo directly to output
    // buffer
    if (unlikely(CheckedVrecLenIncr(het_ctp1_8, vrec_len_ptr))) {
      return 1;
    }
    CopyGenomatchSubset(phaseinfo, genoarr_hets, kMask5555, 1, het_ct, fwrite_bufp_alias);
    fwrite_bufp_final = &(pwcp->fwrite_bufp[het_ctp1_8]);
  } else {
    // this is a minor variant of ExpandThenSubsetBytearr()
    if (unlikely(CheckedVrecLenIncr(het_ctp1_8 + DivUp(phasepresent_ct, 8), vrec_len_ptr))) {
      return 1;
    }
    const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
    uintptr_t* phaseinfo_tmp = pwcp->genovec_invert_buf;
    uintptr_t* phaseinfo_tmp_iter = phaseinfo_tmp;
    uint32_t phasepresent_write_idx_lowbits = 1;
    phaseinfo_write_idx_lowbits = 0;
    uintptr_t phasepresent_write_word = 1;  // first bit set
    for (uint32_t widx = 0; widx != sample_ctl2; ++widx) {
      const uintptr_t geno_word = genoarr_hets[widx];
      uintptr_t geno_hets = Word01(geno_word);
      if (geno_hets) {
        const uint32_t phasepresent_halfword = R_CAST(const Halfword*, phasepresent)[widx];
        if (phasepresent_halfword) {
          const uint32_t phaseinfo_halfword = R_CAST(const Halfword*, phaseinfo)[widx];
          do {
            const uint32_t sample_idx_lowbits = ctzw(geno_hets) / 2;
            if ((phasepresent_halfword >> sample_idx_lowbits) & 1) {
              phasepresent_write_word |= k1LU << phasepresent_write_idx_lowbits;
              phaseinfo_write_word |= S_CAST(uintptr_t, (phaseinfo_halfword >> sample_idx_lowbits) & k1LU) << phaseinfo_write_idx_lowbits;
              if (++phaseinfo_write_idx_lowbits == kBitsPerWord) {
                *phaseinfo_tmp_iter++ = phaseinfo_write_word;
                phaseinfo_write_word = 0;
                phaseinfo_write_idx_lowbits = 0;
              }
            }
            if (++phasepresent_write_idx_lowbits == kBitsPerWord) {
              *fwrite_bufp_alias++ = phasepresent_write_word;
              phasepresent_write_word = 0;
              phasepresent_write_idx_lowbits = 0;
            }
            geno_hets &= geno_hets - 1;
          } while (geno_hets);
        } else {
          phasepresent_write_idx_lowbits += PopcountWord(geno_hets);
          if (phasepresent_write_idx_lowbits >= kBitsPerWord) {
            *fwrite_bufp_alias++ = phasepresent_write_word;
            phasepresent_write_word = 0;
            phasepresent_write_idx_lowbits -= kBitsPerWord;
          }
        }
      }
    }
    fwrite_bufp_final = R_CAST(unsigned char*, fwrite_bufp_alias);
    if (phasepresent_write_idx_lowbits) {
      const uint32_t cur_byte_ct = DivUp(phasepresent_write_idx_lowbits, CHAR_BIT);
      // er, safe to write the entire word...
      SubwordStoreMov(phasepresent_write_word, cur_byte_ct, &fwrite_bufp_final);
    }
    fwrite_bufp_final = memcpyua(fwrite_bufp_final, phaseinfo_tmp, sizeof(intptr_t) * (phaseinfo_tmp_iter - phaseinfo_tmp));
    if (phaseinfo_write_idx_lowbits) {
      const uint32_t cur_byte_ct = DivUp(phaseinfo_write_idx_lowbits, CHAR_BIT);
      SubwordStoreMov(phaseinfo_write_word, cur_byte_ct, &fwrite_bufp_final);
    }
    assert(S_CAST(uintptr_t, fwrite_bufp_final - pwcp->fwrite_bufp) == het_ctp1_8 + DivUp(phasepresent_ct, 8));
  }
  pwcp->fwrite_bufp = fwrite_bufp_final;
  return 0;
}

void PwcAppendBiallelicGenovecHphase(const uintptr_t* __restrict genovec, const uintptr_t* __restrict phasepresent, const uintptr_t* __restrict phaseinfo, PgenWriterCommon* pwcp) {
  // assumes phase_dosage_gflags is nonzero
  const uint32_t vidx = pwcp->vidx;
  unsigned char* vrtype_dest = &(R_CAST(unsigned char*, pwcp->vrtype_buf)[vidx]);
  uint32_t het_ct;
  uint32_t vrec_len = PwcAppendBiallelicGenovecMain(genovec, vidx, pwcp, &het_ct, nullptr, vrtype_dest);
  const uintptr_t vrec_len_byte_ct = pwcp->vrec_len_byte_ct;
  const uint32_t sample_ct = pwcp->sample_ct;
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  pwcp->vidx += 1;
  unsigned char* vrec_len_dest = &(pwcp->vrec_len_buf[vidx * vrec_len_byte_ct]);
  const uint32_t phasepresent_ct = phasepresent? PopcountWords(phasepresent, sample_ctl) : het_ct;
  if (phasepresent_ct) {
    AppendHphase(genovec, phasepresent, phaseinfo, het_ct, phasepresent_ct, pwcp, vrtype_dest, &vrec_len);
  }
  SubU32Store(vrec_len, vrec_len_byte_ct, vrec_len_dest);
}

BoolErr PwcAppendMultiallelicGenovecHphase(const uintptr_t* __restrict genovec, const uintptr_t* __restrict patch_01_set, const AlleleCode* __restrict patch_01_vals, const uintptr_t* __restrict patch_10_set, const AlleleCode* __restrict patch_10_vals, const uintptr_t* __restrict phasepresent, const uintptr_t* __restrict phaseinfo, uint32_t allele_ct, uint32_t patch_01_ct, uint32_t patch_10_ct, PgenWriterCommon* pwcp) {
  const uint32_t vidx = pwcp->vidx;
  const uintptr_t* genovec_hets;
  unsigned char vrtype;
  uint32_t het_ct;
  uint32_t vrec_len;
  if (unlikely(PwcAppendMultiallelicMain(genovec, patch_01_set, patch_01_vals, patch_10_set, patch_10_vals, allele_ct, patch_01_ct, patch_10_ct, vidx, pwcp, &genovec_hets, &het_ct, &vrtype, &vrec_len))) {
    return 1;
  }
  const uint32_t sample_ct = pwcp->sample_ct;
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  const uint32_t phasepresent_ct = phasepresent? PopcountWords(phasepresent, sample_ctl) : het_ct;
  unsigned char* vrtype_dest = &(R_CAST(unsigned char*, pwcp->vrtype_buf)[vidx]);
  *vrtype_dest = vrtype;
  if (phasepresent_ct) {
    if (unlikely(AppendHphase(genovec_hets, phasepresent, phaseinfo, het_ct, phasepresent_ct, pwcp, vrtype_dest, &vrec_len))) {
      return 1;
    }
  }
  pwcp->vidx += 1;
  const uintptr_t vrec_len_byte_ct = pwcp->vrec_len_byte_ct;
  SubU32Store(vrec_len, vrec_len_byte_ct, &(pwcp->vrec_len_buf[vidx * vrec_len_byte_ct]));
  return 0;
}

// ok if dosage_present trailing bits set
BoolErr AppendDosage16(const uintptr_t* __restrict dosage_present, const uint16_t* dosage_main, uint32_t dosage_ct, uint32_t dphase_ct, PgenWriterCommon* pwcp, unsigned char* vrtype_ptr, uint32_t* vrec_len_ptr) {
  const uint32_t sample_ct = pwcp->sample_ct;
  const uint32_t max_deltalist_entry_ct = sample_ct / kPglMaxDeltalistLenDivisor;
  if (dosage_ct < max_deltalist_entry_ct) {
    // case 1: store dosage IDs as deltalist.
    if (unlikely(PwcAppendDeltalist(dosage_present, dosage_ct, pwcp, vrec_len_ptr))) {
      return 1;
    }
    *vrtype_ptr += 0x20;
  } else if ((dosage_ct == sample_ct) && ((!dphase_ct) || (dphase_ct == sample_ct))) {
    // case 2: fixed-width, no need to store dosage IDs at all.
    // dosage_main permitted to have 65535 = missing
    // bugfix (9 Dec 2018): since this forces fixed-width dosage-phase storage,
    // also require dphase_ct == 0 or dphase_ct == sample_ct.
    *vrtype_ptr += 0x40;
  } else {
    // case 3: save dosage_present bitarray directly.
    const uint32_t sample_ctb = DivUp(sample_ct, CHAR_BIT);
    if (unlikely(CheckedVrecLenIncr(sample_ctb, vrec_len_ptr))) {
      return 1;
    }
    pwcp->fwrite_bufp = memcpyua(pwcp->fwrite_bufp, dosage_present, sample_ctb);
    *vrtype_ptr += 0x60;
  }
  if (unlikely(CheckedVrecLenIncr(dosage_ct * sizeof(int16_t), vrec_len_ptr))) {
    return 1;
  }
  pwcp->fwrite_bufp = memcpyua(pwcp->fwrite_bufp, dosage_main, dosage_ct * sizeof(int16_t));
  return 0;
}

// ok if dosage_present trailing bits set
BoolErr PwcAppendBiallelicGenovecDosage16(const uintptr_t* __restrict genovec, const uintptr_t* __restrict dosage_present, const uint16_t* dosage_main, uint32_t dosage_ct, PgenWriterCommon* pwcp) {
  // safe to call this even when entire file has no phase/dosage info
  const uint32_t vidx = pwcp->vidx;
  unsigned char vrtype;
  uint32_t vrec_len = PwcAppendBiallelicGenovecMain(genovec, vidx, pwcp, nullptr, nullptr, &vrtype);
  const uintptr_t vrec_len_byte_ct = pwcp->vrec_len_byte_ct;
  pwcp->vidx += 1;
  unsigned char* vrec_len_dest = &(pwcp->vrec_len_buf[vidx * vrec_len_byte_ct]);
  if (dosage_ct) {
    if (unlikely(AppendDosage16(dosage_present, dosage_main, dosage_ct, 0, pwcp, &vrtype, &vrec_len))) {
      return 1;
    }
  }
  SubU32Store(vrec_len, vrec_len_byte_ct, vrec_len_dest);
  if (!pwcp->phase_dosage_gflags) {
    pwcp->vrtype_buf[vidx / kBitsPerWordD4] |= S_CAST(uintptr_t, vrtype) << (4 * (vidx % kBitsPerWordD4));
  } else {
    R_CAST(unsigned char*, pwcp->vrtype_buf)[vidx] = vrtype;
  }
  return 0;
}

// NOT ok if phasepresent trailing bits set, but ok for dosage_present
BoolErr PwcAppendBiallelicGenovecHphaseDosage16(const uintptr_t* __restrict genovec, const uintptr_t* __restrict phasepresent, const uintptr_t* __restrict phaseinfo, const uintptr_t* __restrict dosage_present, const uint16_t* dosage_main, uint32_t dosage_ct, PgenWriterCommon* pwcp) {
  // assumes there is phase and/or dosage data in output file, otherwise
  // vrtype_dest needs to be replaced

  // this mostly overlaps with PwcAppendBiallelicGenovecHphase(); probably
  // get rid of the latter
  const uint32_t vidx = pwcp->vidx;
  unsigned char* vrtype_dest = &(R_CAST(unsigned char*, pwcp->vrtype_buf)[vidx]);
  uint32_t het_ct;
  uint32_t vrec_len = PwcAppendBiallelicGenovecMain(genovec, vidx, pwcp, &het_ct, nullptr, vrtype_dest);
  const uintptr_t vrec_len_byte_ct = pwcp->vrec_len_byte_ct;
  const uint32_t sample_ct = pwcp->sample_ct;
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  pwcp->vidx += 1;
  unsigned char* vrec_len_dest = &(pwcp->vrec_len_buf[vidx * vrec_len_byte_ct]);
  const uint32_t phasepresent_ct = phasepresent? PopcountWords(phasepresent, sample_ctl) : het_ct;
  if (phasepresent_ct) {
    AppendHphase(genovec, phasepresent, phaseinfo, het_ct, phasepresent_ct, pwcp, vrtype_dest, &vrec_len);
  }
  if (dosage_ct) {
    if (unlikely(AppendDosage16(dosage_present, dosage_main, dosage_ct, 0, pwcp, vrtype_dest, &vrec_len))) {
      return 1;
    }
  }
  SubU32Store(vrec_len, vrec_len_byte_ct, vrec_len_dest);
  return 0;
}

BoolErr AppendDphase16(const uintptr_t* __restrict dosage_present, const uintptr_t* __restrict dphase_present, const int16_t* dphase_delta, uint32_t dosage_ct, uint32_t dphase_ct, PgenWriterCommon* pwcp, unsigned char* vrtype_ptr, uint32_t* vrec_len_ptr) {
  assert(dphase_ct);
  *vrtype_ptr += 0x80;
  if (dphase_ct != pwcp->sample_ct) {
    // bugfix (9 Dec 2018): forgot to separate fixed-width phased-dosage case
    // bugfix (12 Feb 2019): double-added dphase_ct * sizeof(int16_t)
    const uint32_t dphase_present_byte_ct = DivUp(dosage_ct, CHAR_BIT);
    if (unlikely(CheckedVrecLenIncr(dphase_present_byte_ct, vrec_len_ptr))) {
      return 1;
    }
    CopyBitarrSubset(dphase_present, dosage_present, dosage_ct, R_CAST(uintptr_t*, pwcp->fwrite_bufp));
    pwcp->fwrite_bufp = &(pwcp->fwrite_bufp[dphase_present_byte_ct]);
  }
  // bugfix (2 Jan 2019): forgot to update vrec_len here
  if (unlikely(CheckedVrecLenIncr(dphase_ct * sizeof(int16_t), vrec_len_ptr))) {
    return 1;
  }
  pwcp->fwrite_bufp = memcpyua(pwcp->fwrite_bufp, dphase_delta, dphase_ct * sizeof(int16_t));
  return 0;
}

BoolErr PwcAppendBiallelicGenovecDphase16(const uintptr_t* __restrict genovec, const uintptr_t* __restrict phasepresent, const uintptr_t* __restrict phaseinfo, const uintptr_t* __restrict dosage_present, const uintptr_t* __restrict dphase_present, const uint16_t* dosage_main, const int16_t* dphase_delta, uint32_t dosage_ct, uint32_t dphase_ct, PgenWriterCommon* pwcp) {
  // assumes there is phase and/or dosage data in output file, otherwise
  // vrtype_dest needs to be replaced
  const uint32_t vidx = pwcp->vidx;
  unsigned char* vrtype_dest = &(R_CAST(unsigned char*, pwcp->vrtype_buf)[vidx]);
  uint32_t het_ct;
  uint32_t vrec_len = PwcAppendBiallelicGenovecMain(genovec, vidx, pwcp, &het_ct, nullptr, vrtype_dest);
  const uintptr_t vrec_len_byte_ct = pwcp->vrec_len_byte_ct;
  const uint32_t sample_ct = pwcp->sample_ct;
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  pwcp->vidx += 1;
  unsigned char* vrec_len_dest = &(pwcp->vrec_len_buf[vidx * vrec_len_byte_ct]);
  const uint32_t phasepresent_ct = phasepresent? PopcountWords(phasepresent, sample_ctl) : het_ct;
  if (phasepresent_ct) {
    AppendHphase(genovec, phasepresent, phaseinfo, het_ct, phasepresent_ct, pwcp, vrtype_dest, &vrec_len);
  }
  if (dosage_ct) {
    if (unlikely(AppendDosage16(dosage_present, dosage_main, dosage_ct, dphase_ct, pwcp, vrtype_dest, &vrec_len))) {
      return 1;
    }
    if (dphase_ct) {
      if (unlikely(AppendDphase16(dosage_present, dphase_present, dphase_delta, dosage_ct, dphase_ct, pwcp, vrtype_dest, &vrec_len))) {
        return 1;
      }
    }
  }
  SubU32Store(vrec_len, vrec_len_byte_ct, vrec_len_dest);
  return 0;
}

uint32_t PhaseOrDosagePresent(const uintptr_t* vrtype_buf, uint32_t variant_ct) {
  // Assumes trailing bits of vrtype_buf are zeroed out.  (This happens during
  // vrtype_buf initialization.)
  const uint32_t word_ct = DivUp(variant_ct, kBytesPerWord);
  for (uint32_t widx = 0; widx != word_ct; ++widx) {
    if (vrtype_buf[widx] & (kMask0101 * 0xf0)) {
      return 1;
    }
  }
  return 0;
}

PglErr PwcFinish(PgenWriterCommon* pwcp, FILE** pgen_outfile_ptr, FILE** pgi_or_final_pgen_outfile_ptr, char** fname_buf_ptr) {
  const uint32_t variant_ct = pwcp->vidx;
  FILE** header_ff_ptr;
  if (*pgi_or_final_pgen_outfile_ptr) {
    // kPgenWriteSeparateIndex, kPgenWriteAndCopy
    assert(variant_ct);
    assert(pwcp->variant_ct_limit >= variant_ct);
    if (pwcp->phase_dosage_gflags) {
      // Check if there is at least one variant with a phase or dosage data
      // track.
      if (!PhaseOrDosagePresent(pwcp->vrtype_buf, variant_ct)) {
        // If not, shrink the index.
        Reduce8to4bitInplaceUnsafe(variant_ct, pwcp->vrtype_buf);
        pwcp->phase_dosage_gflags = kfPgenGlobal0;
      }
    }
    if (unlikely(fclose_null(pgen_outfile_ptr))) {
      return kPglRetWriteFail;
    }
    header_ff_ptr = pgi_or_final_pgen_outfile_ptr;
    if (*fname_buf_ptr) {
      *pgen_outfile_ptr = fopen(*fname_buf_ptr, FOPEN_RB);
      if (unlikely(!(*pgen_outfile_ptr))) {
        return kPglRetOpenFail;
      }
    }
  } else {
    // kPgenWriteBackwardSeek
    assert(pwcp->variant_ct_limit == variant_ct);
    header_ff_ptr = pgen_outfile_ptr;
    if (unlikely(fseeko(*header_ff_ptr, 3, SEEK_SET))) {
      return kPglRetWriteFail;
    }
  }
  FILE* header_ff = *header_ff_ptr;
  const PgenGlobalFlags phase_dosage_gflags = pwcp->phase_dosage_gflags;
  fwrite_unlocked(&variant_ct, sizeof(int32_t), 1, header_ff);
  fwrite_unlocked(&(pwcp->sample_ct), sizeof(int32_t), 1, header_ff);

  const uint32_t vrec_len_byte_ct = pwcp->vrec_len_byte_ct;
  const uint32_t control_byte = (vrec_len_byte_ct - 1) + (4 * (phase_dosage_gflags != 0)) + (pwcp->nonref_flags_storage << 6);
  putc_unlocked(control_byte, header_ff);

  const uint32_t vblock_ct = DivUp(variant_ct, kPglVblockSize);
  if (*fname_buf_ptr) {
    // Current vblock_fpos values are with no index.  Determine the size of the
    // index and add it to all array entries.

    // variant_ct, sample_ct, control_byte, vblock_fpos itself
    uintptr_t index_size = 9 + vblock_ct * sizeof(int64_t);
    // vrtypes
    if (phase_dosage_gflags != 0) {
      index_size += variant_ct;
    } else {
      index_size += DivUp(variant_ct, 2);
    }
    // vrec_lens
    index_size += vrec_len_byte_ct * S_CAST(uintptr_t, variant_ct);
    // nonref_flags
    if (pwcp->explicit_nonref_flags) {
      index_size += DivUp(variant_ct, CHAR_BIT);
    }
    for (uint32_t vblock_idx = 0; vblock_idx != vblock_ct; ++vblock_idx) {
      pwcp->vblock_fpos[vblock_idx] += index_size;
    }
  }
  fwrite_unlocked(pwcp->vblock_fpos, vblock_ct * sizeof(int64_t), 1, header_ff);
  const unsigned char* vrtype_buf_iter = R_CAST(unsigned char*, pwcp->vrtype_buf);
  const unsigned char* vrec_len_buf_iter = pwcp->vrec_len_buf;
  uint32_t vrec_iter_incr = kPglVblockSize * vrec_len_byte_ct;
  uint32_t vrtype_buf_iter_incr = phase_dosage_gflags? kPglVblockSize : (kPglVblockSize / 2);
  uint32_t nonref_flags_write_byte_ct = kPglVblockSize / CHAR_BIT;
  const unsigned char* vrec_len_buf_last = &(vrec_len_buf_iter[S_CAST(uintptr_t, vblock_ct - 1) * vrec_iter_incr]);
  uintptr_t* explicit_nonref_flags = pwcp->explicit_nonref_flags;
  uintptr_t* explicit_nonref_flags_iter = explicit_nonref_flags;
  for (; ; vrec_len_buf_iter = &(vrec_len_buf_iter[vrec_iter_incr])) {
    if (vrec_len_buf_iter >= vrec_len_buf_last) {
      if (vrec_len_buf_iter > vrec_len_buf_last) {
        break;
      }
      const uint32_t vblock_size = ModNz(variant_ct, kPglVblockSize);
      vrtype_buf_iter_incr = phase_dosage_gflags? vblock_size : DivUp(vblock_size, 2);
      vrec_iter_incr = vblock_size * vrec_len_byte_ct;
      nonref_flags_write_byte_ct = DivUp(vblock_size, CHAR_BIT);
    }
    // 4b(i): array of 4-bit or 1-byte vrtypes
    fwrite_unlocked(vrtype_buf_iter, vrtype_buf_iter_incr, 1, header_ff);
    vrtype_buf_iter = &(vrtype_buf_iter[vrtype_buf_iter_incr]);

    // 4b(ii): array of variant record lengths
    if (unlikely(fwrite_checked(vrec_len_buf_iter, vrec_iter_incr, header_ff))) {
      return kPglRetWriteFail;
    }

    // If we were writing allele_idx_offsets, that would happen here.

    // 4b(iii): explicit nonref flags
    if (explicit_nonref_flags) {
      if (unlikely(fwrite_checked(explicit_nonref_flags_iter, nonref_flags_write_byte_ct, header_ff))) {
        return kPglRetWriteFail;
      }
      explicit_nonref_flags_iter = &(explicit_nonref_flags_iter[kPglVblockSize / kBitsPerWord]);
    }
  }
  if (!(*fname_buf_ptr)) {
    return fclose_null(header_ff_ptr)? kPglRetWriteFail : kPglRetSuccess;
  }
  // Append temporary-.pgen body (excluding the 3 bytes containing the leading
  // magic number) to final output file.
  unsigned char copybuf[kPglFwriteBlockSize + 3];
  uintptr_t nbyte = fread(copybuf, 1, kPglFwriteBlockSize + 3, *pgen_outfile_ptr);
  if (unlikely((nbyte <= 3) || (!memequal_k(copybuf, "l\x1b\x20", 3)))) {
    return kPglRetRewindFail;
  }
  nbyte -= 3;
  if (unlikely(!fwrite_unlocked(&(copybuf[3]), nbyte, 1, header_ff))) {
    return kPglRetWriteFail;
  }
  while (1) {
    nbyte = fread(copybuf, 1, kPglFwriteBlockSize, *pgen_outfile_ptr);
    if (!nbyte) {
      break;
    }
    if (unlikely(!fwrite_unlocked(copybuf, nbyte, 1, header_ff))) {
      return kPglRetWriteFail;
    }
  }
  if (unlikely(fclose_null(pgen_outfile_ptr) ||
               unlink(*fname_buf_ptr))) {
    return kPglRetWriteFail;
  }
  free(*fname_buf_ptr);
  *fname_buf_ptr = nullptr;
  return kPglRetSuccess;
}

PglErr SpgwFinish(STPgenWriter* spgwp) {
  PgenWriterCommon* pwcp = GetPwcp(spgwp);
  FILE** pgen_outfilep = GetPgenOutfilep(spgwp);
  if (unlikely(fwrite_checked(pwcp->fwrite_buf, pwcp->fwrite_bufp - pwcp->fwrite_buf, *pgen_outfilep))) {
    return kPglRetWriteFail;
  }
  FILE** pgi_or_final_pgen_outfilep = GetPgiOrFinalPgenOutfilep(spgwp);
  char** fname_bufp = GetFnameBufp(spgwp);
  return PwcFinish(pwcp, pgen_outfilep, pgi_or_final_pgen_outfilep, fname_bufp);
}

PglErr MpgwFlush(MTPgenWriter* mpgwp) {
  PgenWriterCommon* pwcp = mpgwp->pwcs[0];
  uint32_t vidx = RoundDownPow2(pwcp->vidx - 1, kPglVblockSize);
  uint32_t thread_ct = mpgwp->thread_ct;
  const uint32_t variant_ct = pwcp->variant_ct_limit;
  const uint32_t is_last_flush = ((vidx + thread_ct * kPglVblockSize) >= variant_ct);
  if (is_last_flush) {
    thread_ct = DivUp(variant_ct - vidx, kPglVblockSize);
  }
  uint64_t* vblock_fpos = pwcp->vblock_fpos;
  FILE* pgen_outfile = mpgwp->pgen_outfile;
  const uint32_t vidx_incr = (thread_ct - 1) * kPglVblockSize;
  uint64_t cur_vblock_fpos = ftello(pgen_outfile);
  for (uint32_t tidx = 0; tidx != thread_ct; ++tidx) {
    vblock_fpos[(vidx / kPglVblockSize) + tidx] = cur_vblock_fpos;
    PgenWriterCommon* cur_pwcp = mpgwp->pwcs[tidx];
    uintptr_t cur_vblock_byte_ct = cur_pwcp->fwrite_bufp - cur_pwcp->fwrite_buf;
    if (unlikely(fwrite_checked(cur_pwcp->fwrite_buf, cur_vblock_byte_ct, pgen_outfile))) {
      return kPglRetWriteFail;
    }
    cur_pwcp->vidx += vidx_incr;
    cur_pwcp->fwrite_bufp = cur_pwcp->fwrite_buf;
    cur_vblock_fpos += cur_vblock_byte_ct;
  }
  if (!is_last_flush) {
    return kPglRetSuccess;
  }
  pwcp->vidx = variant_ct;
  return PwcFinish(pwcp, &(mpgwp->pgen_outfile), &(mpgwp->pgi_or_final_pgen_outfile), &(mpgwp->fname_buf));
}

BoolErr CleanupSpgw(STPgenWriter* spgwp, PglErr* reterrp) {
  // assume file is open if spgw.pgen_outfile is not null
  // memory is the responsibility of the caller for now
  FILE** pgen_outfilep = GetPgenOutfilep(spgwp);
  FILE** pgi_or_final_pgen_outfilep = GetPgiOrFinalPgenOutfilep(spgwp);
  char** fname_bufp = GetFnameBufp(spgwp);
  BoolErr fclose_err = 0;
  if (*pgi_or_final_pgen_outfilep) {
    fclose_err = fclose_null(pgi_or_final_pgen_outfilep);
  }
  if (*pgen_outfilep) {
    if (fclose_null(pgen_outfilep)) {
      fclose_err = 1;
    }
  }
  if (*fname_bufp) {
    free(*fname_bufp);
    *fname_bufp = nullptr;
  }
  if (unlikely(fclose_err && (*reterrp == kPglRetSuccess))) {
    *reterrp = kPglRetWriteFail;
  }
  return fclose_err;
}

BoolErr CleanupMpgw(MTPgenWriter* mpgwp, PglErr* reterrp) {
  if (!mpgwp) {
    return 0;
  }
  FILE** pgen_outfilep = &(mpgwp->pgen_outfile);
  FILE** pgi_or_final_pgen_outfilep = &(mpgwp->pgi_or_final_pgen_outfile);
  char** fname_bufp = &(mpgwp->fname_buf);
  BoolErr fclose_err = 0;
  if (*pgi_or_final_pgen_outfilep) {
    fclose_err = fclose_null(pgi_or_final_pgen_outfilep);
  }
  if (*pgen_outfilep) {
    if (fclose_null(pgen_outfilep)) {
      fclose_err = 1;
    }
  }
  if (*fname_bufp) {
    free(*fname_bufp);
    *fname_bufp = nullptr;
  }
  if (unlikely(fclose_err && (*reterrp == kPglRetSuccess))) {
    *reterrp = kPglRetWriteFail;
  }
  return fclose_err;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
