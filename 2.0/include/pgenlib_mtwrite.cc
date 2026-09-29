// This library is part of PLINK 2.0, copyright (C) 2005-2026 Shaun Purcell,
// Christopher Chang, Benjamin Demaille.
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

#include "pgenlib_mtwrite.h"

#include <assert.h>
#include <limits.h>
#include <stdlib.h>

#include "plink2_bits.h"

#ifdef __cplusplus
namespace plink2 {
#endif

void PreinitMpgw(MTPgenWriter* mpgwp) {
  mpgwp->pgen_outfile = nullptr;
  mpgwp->pgi_or_final_pgen_outfile = nullptr;
  mpgwp->fname_buf = nullptr;
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

PglErr MpgwInitPhase2Ex(const char* __restrict fname, uintptr_t* __restrict explicit_nonref_flags, PgenExtensionLl* header_exts, PgenExtensionLl* footer_exts, uint32_t variant_ct, uint32_t sample_ct, PgenWriteMode write_mode, PgenGlobalFlags phase_dosage_gflags, uint32_t nonref_flags_storage, uint32_t vrec_len_byte_ct, uintptr_t vblock_cacheline_ct, uint32_t thread_ct, unsigned char* mpgw_alloc, MTPgenWriter* mpgwp) {
  assert(thread_ct);
  const uintptr_t pwc_byte_ct = RoundUpPow2(sizeof(PgenWriterCommon), kCacheline);
  for (uint32_t tidx = 0; tidx != thread_ct; ++tidx) {
    mpgwp->pwcs[tidx] = R_CAST(PgenWriterCommon*, &(mpgw_alloc[tidx * pwc_byte_ct]));
  }
  PglErr reterr = PwcInitPhase1(fname, explicit_nonref_flags, header_exts, footer_exts, variant_ct, sample_ct, write_mode, phase_dosage_gflags, nonref_flags_storage, vrec_len_byte_ct, mpgwp->pwcs[0], &(mpgwp->pgen_outfile), &(mpgwp->pgi_or_final_pgen_outfile), &(mpgwp->fname_buf));
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
