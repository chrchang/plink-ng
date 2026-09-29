#ifndef __PGENLIB_MTWRITE_H__
#define __PGENLIB_MTWRITE_H__

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


// pgenlib_mtwrite contains multithreaded-writer-specific code.
// It has been separated out from pgenlib_write due to CRAN coding standards:
// pgenlib_write is now used by pgenlibr, and CRAN disallows
// MTPgenWriterStruct'st flexible array member.

#include "pgenlib_misc.h"
#include "pgenlib_write.h"
#include "plink2_base.h"

#ifdef __cplusplus
namespace plink2 {
#endif

typedef struct MTPgenWriterStruct {
  MOVABLE_BUT_NONCOPYABLE(MTPgenWriterStruct);
  FILE* pgen_outfile;
  FILE* pgi_or_final_pgen_outfile;
  char* fname_buf;
  uint32_t thread_ct;
  PgenWriterCommon* pwcs[];
} MTPgenWriter;

void PreinitMpgw(MTPgenWriter* mpgwp);

// moderately likely that there isn't enough memory to use the maximum number
// of threads, so this returns per-thread memory requirements before forcing
// the caller to specify thread count
// (eventually should write code which falls back on STPgenWriter
// when there isn't enough memory for even a single 64k variant block, at least
// for the most commonly used plink 2.0 functions)
void MpgwInitPhase1(const uintptr_t* __restrict allele_idx_offsets, uint32_t variant_ct, uint32_t sample_ct, PgenGlobalFlags phase_dosage_gflags, uintptr_t* alloc_base_cacheline_ct_ptr, uint64_t* alloc_per_thread_cacheline_ct_ptr, uint32_t* vrec_len_byte_ct_ptr, uint64_t* vblock_cacheline_ct_ptr);

// Caller is responsible for printing open-fail error message.
PglErr MpgwInitPhase2Ex(const char* __restrict fname, uintptr_t* __restrict explicit_nonref_flags, PgenExtensionLl* header_exts, PgenExtensionLl* footer_exts, uint32_t variant_ct, uint32_t sample_ct, PgenWriteMode write_mode, PgenGlobalFlags phase_dosage_gflags, uint32_t nonref_flags_storage, uint32_t vrec_len_byte_ct, uintptr_t vblock_cacheline_ct, uint32_t thread_ct, unsigned char* mpgw_alloc, MTPgenWriter* mpgwp);

HEADER_INLINE PglErr MpgwInitPhase2(const char* __restrict fname, uintptr_t* __restrict explicit_nonref_flags, uint32_t variant_ct, uint32_t sample_ct, PgenWriteMode write_mode, PgenGlobalFlags phase_dosage_gflags, uint32_t nonref_flags_storage, uint32_t vrec_len_byte_ct, uintptr_t vblock_cacheline_ct, uint32_t thread_ct, unsigned char* mpgw_alloc, MTPgenWriter* mpgwp) {
  return MpgwInitPhase2Ex(fname, explicit_nonref_flags, nullptr, nullptr, variant_ct, sample_ct, write_mode, phase_dosage_gflags, nonref_flags_storage, vrec_len_byte_ct, vblock_cacheline_ct, thread_ct, mpgw_alloc, mpgwp);
}

// Last flush automatically writes footer if present, backfills header, and
// closes the file.
// (caller should set mpgwp = nullptr after that)
PglErr MpgwFlush(MTPgenWriter* mpgwp);


// this closes the file if open, but does not free any memory
// handles mpgwp == nullptr, since it shouldn't be allocated on the stack
// error-return iff reterr was success and was changed to kPglRetWriteFail
// (i.e. an error message should be printed), though this is not relevant for
// plink2
BoolErr CleanupMpgw(MTPgenWriter* mpgwp, PglErr* reterrp);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PGENLIB_MTWRITE_H__
