// This file is part of PLINK 2.0, copyright (C) 2005-2025 Shaun Purcell,
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

#include "plink2_import.h"

#include <assert.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "include/pgenlib_misc.h"
#include "include/pgenlib_write.h"
#include "include/plink2_bgzf.h"
#include "include/plink2_bits.h"
#include "include/plink2_htable.h"
#include "include/plink2_memory.h"
#include "include/plink2_string.h"
#include "include/plink2_text.h"
#include "include/plink2_thread.h"
#include "plink2_cmdline.h"
#include "plink2_compress_stream.h"
#include "plink2_decompress.h"
#include "plink2_import_legacy.h"
#include "plink2_psam.h"
#include "plink2_pvar.h"
#include "plink2_random.h"

#ifdef __cplusplus
namespace plink2 {
#endif

void InitPlink1Dosage(Plink1DosageInfo* plink1_dosage_info_ptr) {
  plink1_dosage_info_ptr->flags = kfPlink1Dosage0;
  STD_ARRAY_FILL0(plink1_dosage_info_ptr->skips);
  plink1_dosage_info_ptr->chr_col_idx = UINT32_MAX;
  plink1_dosage_info_ptr->pos_col_idx = UINT32_MAX;
  plink1_dosage_info_ptr->id_delim = '\0';
}

void InitGenDummy(GenDummyInfo* gendummy_info_ptr) {
  gendummy_info_ptr->flags = kfGenDummy0;
  gendummy_info_ptr->pheno_ct = 1;
  gendummy_info_ptr->geno_mfreq_ct = 0;
  gendummy_info_ptr->geno_mfreqs = nullptr; // defensive
  gendummy_info_ptr->pheno_mfreq = 0.0;
  gendummy_info_ptr->phase_freq = 0.0;
  gendummy_info_ptr->dosage_freq = 0.0;
}

void CleanupGenDummy(GenDummyInfo* gendummy_info_ptr) {
  free_cond(gendummy_info_ptr->geno_mfreqs);
}


// Basic workflow for fully-powered import functions (VCF/BCF, BGEN-1.3):
// * Perform an initial scan of the file to count the number of variants we're
//   actually importing, save allele counts of each of those variants, and
//   determine maximum input and output record sizes (checking whether phase or
//   dosage information is present in the process).
// * Rewind, shrink file-read buffer, set worker thread count, size and
//   allocate file-parsing double-buffer, initialize .pgen writer.
// * Iterate through file again.  On the first main loop iteration, the I/O
//   thread just appends raw records to the first half of the double-buffer,
//   spacing them so that the parsed results can be written back in-place.  On
//   the second iteration, the I/O thread appends raw records to the second
//   half of the double-buffer while the worker threads replace the raw bytes
//   in the first half with a parsed form.  On subsequent iterations, the I/O
//   thread uses the relevant SpgwAppend... functions to flush parsed data in
//   its side of the double-buffer, then overwrites it with new raw records.
// This does not currently make sense for --make-pgen, since (unless a
// compute-heavy operation like --indiv-sort is involved) the worker threads
// usually wouldn't have enough to do.

// allele count determined from global allele_idx_offsets[] array
FLAGSET_DEF_START()
  kfGparse0,
  kfGparseNull = (1 << 0),  // just create an all-missing record
  kfGparseHphase = (1 << 1),
  kfGparseDosage = (1 << 2),
  kfGparseDphase = (1 << 3)
FLAGSET_DEF_END(GparseFlags);

typedef struct GparseReadBgenMetadataStruct {
  uint32_t input_byte_ct;
  uint32_t uncompressed_byte_ct;
} GparseReadBgenMetadata;

typedef struct GparseReadVcfMetadataStruct {
  STD_ARRAY_DECL(uint32_t, 2, qual_field_idxs);
  uint16_t gt_exists;
  uint16_t qual_exists;
  uint32_t dosage_field_idx;
  uint32_t hds_field_idx;
  uintptr_t line_idx;  // for error reporting
} GparseReadVcfMetadata;

typedef struct GparseReadBcfMetadataStruct {
  // [0] = GQ, [1] = DP
  STD_ARRAY_DECL(uint32_t, 2, qual_vec_offsets);
  uint32_t gt_vec_offset;
  uint32_t dosage_vec_offset;
  uint32_t hds_vec_offset;
  uintptr_t rec_idx;  // for error reporting
} GparseReadBcfMetadata;

// The "parsed form" is as follows (all arrays vector-aligned, vector counts in
// parentheses):
//   genovec (sample_ctv2)
//   iff allele_ct > 2 set:
//     patch_01_set (sample_ctv)
//     patch_01_vals (DivUp(sample_ct, kDosagePerVec))
//     patch_10_set (sample_ctv)
//     patch_10_vals (DivUp(sample_ct * 2, kDosagePerVec))
//   iff kfGparseHphase set:
//     phasepresent (sample_ctv)
//     phaseinfo (sample_ctv)
//   iff kfGparseDosage set:
//     dosage_present (sample_ctv)
//     dosage_main (DivUp(sample_ct, kDosagePerVec))
//     iff allele_ct > 2:
//       todo; probably want to exploit alternate size limit
//     iff kfGparseDphase set:
//       dphase_present (sample_ctv)
//       dphase_delta (DivUp(sample_ct, kDosagePerVec))
//       iff allele_ct > 2:
//         todo; probably want to exploit alternate size limit

uint64_t GparseWriteByteCt(uint32_t sample_ct, uint32_t allele_ct, GparseFlags flags) {
  const uint32_t is_multiallelic = (allele_ct > 2);
  uint64_t vec_ct = NypCtToVecCt(sample_ct);
  if ((flags & (kfGparseHphase | kfGparseDosage | kfGparseDphase)) || is_multiallelic) {
    const uint32_t sample_ctv = BitCtToVecCt(sample_ct);
    if (is_multiallelic) {
      // patch_01_set, patch_10_set: sample_ctv each
      // patch_01_vals: sample_ct AlleleCodes
      // patch_10_vals: 2 * sample_ct AlleleCodes
      vec_ct += 2 * sample_ctv + AlleleCodeCtToVecCt(sample_ct) + AlleleCodeCtToVecCt(sample_ct * 2);
    }
    if (flags & kfGparseHphase) {
      vec_ct += 2 * sample_ctv;
    }
    if (flags & kfGparseDosage) {
      const uint32_t dosage_ctv = DosageCtToVecCt(sample_ct);
      uintptr_t vec_incr = sample_ctv + dosage_ctv;
      if (is_multiallelic) {
        // todo
        // note that in VCF/BCF cases, we can compute an alternate upper bound
        // based on input record size, and this can buy us a lot
        logerrprintf("Error: GparseWriteByteCt multiallelic dosage size request.\n");
        exit(S_CAST(int32_t, kPglRetInternalError));
      }
      vec_ct += vec_incr * (1 + ((flags / kfGparseDphase) & 1));
    }
  }
  return vec_ct * kBytesPerVec;
}

typedef struct GparseWriteMetadataStruct {
  uint32_t patch_01_ct;
  uint32_t patch_10_ct;
  uint32_t phasepresent_exists;  // may want count later
  uint32_t dosage_ct;
  uint32_t multiallelic_dosage_ct;
  uint32_t dphase_ct;
  uint32_t multiallelic_dphase_ct;
} GparseWriteMetadata;

union GparseMetadata {
  GparseReadVcfMetadata read_vcf;
  GparseReadBcfMetadata read_bcf;
  GparseReadBgenMetadata read_bgen;
  GparseWriteMetadata write;
};

typedef struct GparseRecordStruct {
  unsigned char* record_start;  // vector-aligned
  GparseFlags flags;
  union GparseMetadata metadata;
} GparseRecord;

// returns genovec
uintptr_t* GparseGetPointers(unsigned char* record_start, uint32_t sample_ct, uint32_t allele_ct, GparseFlags flags, uintptr_t** patch_01_set_ptr, AlleleCode** patch_01_vals_ptr, uintptr_t** patch_10_set_ptr, AlleleCode** patch_10_vals_ptr, uintptr_t** phasepresent_ptr, uintptr_t** phaseinfo_ptr, uintptr_t** dosage_present_ptr, Dosage** dosage_main_ptr, uintptr_t** dphase_present_ptr, SDosage** dphase_delta_ptr) {
  const uint32_t is_multiallelic = (allele_ct > 2);
  uintptr_t* genovec = R_CAST(uintptr_t*, record_start);
  if ((flags & (kfGparseHphase | kfGparseDosage | kfGparseDphase)) || is_multiallelic) {
    uintptr_t* record_iter = &(genovec[NypCtToAlignedWordCt(sample_ct)]);
    const uint32_t sample_ctaw = BitCtToAlignedWordCt(sample_ct);
    if (is_multiallelic) {
      *patch_01_set_ptr = record_iter;
      record_iter = &(record_iter[sample_ctaw]);
      *patch_01_vals_ptr = R_CAST(AlleleCode*, record_iter);
      record_iter = &(record_iter[AlleleCodeCtToAlignedWordCt(sample_ct)]);
      *patch_10_set_ptr = record_iter;
      record_iter = &(record_iter[sample_ctaw]);
      *patch_10_vals_ptr = R_CAST(AlleleCode*, record_iter);
      record_iter = &(record_iter[kWordsPerVec * DivUp(2 * sample_ct, kAlleleCodesPerVec)]);
    }
    if (flags & kfGparseHphase) {
      *phasepresent_ptr = record_iter;
      record_iter = &(record_iter[sample_ctaw]);
      *phaseinfo_ptr = record_iter;
      record_iter = &(record_iter[sample_ctaw]);
    }
    if (flags & kfGparseDosage) {
      const uint32_t dosage_ctaw = DosageCtToAlignedWordCt(sample_ct);
      *dosage_present_ptr = record_iter;
      record_iter = &(record_iter[sample_ctaw]);
      *dosage_main_ptr = R_CAST(Dosage*, record_iter);
      record_iter = &(record_iter[dosage_ctaw]);
      if (is_multiallelic) {
        // todo
      }
      if (flags & kfGparseDphase) {
        *dphase_present_ptr = record_iter;
        record_iter = &(record_iter[sample_ctaw]);
        *dphase_delta_ptr = R_CAST(SDosage*, record_iter);
        // record_iter = &(record_iter[dosage_ctaw]);
        if (is_multiallelic) {
          // todo
        }
      }
    }
  }
  return genovec;
}

// Might want to change allele_idx_offsets parameter to
// block_allele_idx_offsets, to facilitate single-pass .bgen parsing.
PglErr GparseFlush(const GparseRecord* grp, const uintptr_t* allele_idx_offsets, uint32_t write_block_size, STPgenWriter* spgwp) {
  PglErr reterr = kPglRetSuccess;
  {
    if (allele_idx_offsets) {
      allele_idx_offsets = &(allele_idx_offsets[SpgwGetVidx(spgwp)]);
    }
    const uint32_t sample_ct = SpgwGetSampleCt(spgwp);
    uint32_t allele_ct = 2;
    uintptr_t* patch_01_set = nullptr;
    AlleleCode* patch_01_vals = nullptr;
    uintptr_t* patch_10_set = nullptr;
    AlleleCode* patch_10_vals = nullptr;
    uintptr_t* phasepresent = nullptr;
    uintptr_t* phaseinfo = nullptr;
    uintptr_t* dosage_present = nullptr;
    Dosage* dosage_main = nullptr;
    uintptr_t* dphase_present = nullptr;
    SDosage* dphase_delta = nullptr;
    // todo: multiallelic dosage buffers
    for (uintptr_t write_block_vidx = 0; write_block_vidx != write_block_size; ++write_block_vidx) {
      const GparseRecord* cur_gparse_rec = &(grp[write_block_vidx]);
      const GparseFlags flags = cur_gparse_rec->flags;
      if (allele_idx_offsets) {
        allele_ct = allele_idx_offsets[write_block_vidx + 1] - allele_idx_offsets[write_block_vidx];
      }
      uintptr_t* genovec = GparseGetPointers(cur_gparse_rec->record_start, sample_ct, allele_ct, flags, &patch_01_set, &patch_01_vals, &patch_10_set, &patch_10_vals, &phasepresent, &phaseinfo, &dosage_present, &dosage_main, &dphase_present, &dphase_delta);
      const GparseWriteMetadata* cur_gwmp = &(cur_gparse_rec->metadata.write);
      if (!(flags & (kfGparseHphase | kfGparseDosage))) {
        if (allele_ct == 2) {
          if (unlikely(SpgwAppendBiallelicGenovec(genovec, spgwp))) {
            goto GparseFlush_ret_WRITE_FAIL;
          }
        } else {
          reterr = SpgwAppendMultiallelicSparse(genovec, patch_01_set, patch_01_vals, patch_10_set, patch_10_vals, allele_ct, cur_gwmp->patch_01_ct, cur_gwmp->patch_10_ct, spgwp);
          if (unlikely(reterr)) {
            goto GparseFlush_ret_1;
          }
        }
      } else {
        const uint32_t cur_dosage_ct = cur_gwmp->dosage_ct;
        const uint32_t cur_dphase_ct = cur_gwmp->dphase_ct;
        if (!cur_dphase_ct) {
          if (!cur_dosage_ct) {
            if (allele_ct == 2) {
              if (!cur_gwmp->phasepresent_exists) {
                if (unlikely(SpgwAppendBiallelicGenovec(genovec, spgwp))) {
                  goto GparseFlush_ret_WRITE_FAIL;
                }
              } else {
                if (unlikely(SpgwAppendBiallelicGenovecHphase(genovec, phasepresent, phaseinfo, spgwp))) {
                  goto GparseFlush_ret_WRITE_FAIL;
                }
              }
            } else {
              if (!cur_gwmp->phasepresent_exists) {
                reterr = SpgwAppendMultiallelicSparse(genovec, patch_01_set, patch_01_vals, patch_10_set, patch_10_vals, allele_ct, cur_gwmp->patch_01_ct, cur_gwmp->patch_10_ct, spgwp);
              } else {
                reterr = SpgwAppendMultiallelicGenovecHphase(genovec, patch_01_set, patch_01_vals, patch_10_set, patch_10_vals, phasepresent, phaseinfo, allele_ct, cur_gwmp->patch_01_ct, cur_gwmp->patch_10_ct, spgwp);
              }
              if (unlikely(reterr)) {
                goto GparseFlush_ret_1;
              }
            }
          } else {
            if (!cur_gwmp->phasepresent_exists) {
              reterr = SpgwAppendBiallelicGenovecDosage16(genovec, dosage_present, dosage_main, cur_dosage_ct, spgwp);
            } else {
              reterr = SpgwAppendBiallelicGenovecHphaseDosage16(genovec, phasepresent, phaseinfo, dosage_present, dosage_main, cur_dosage_ct, spgwp);
            }
            if (unlikely(reterr)) {
              goto GparseFlush_ret_1;
            }
          }
        } else {
          reterr = SpgwAppendBiallelicGenovecDphase16(genovec, phasepresent, phaseinfo, dosage_present, dphase_present, dosage_main, dphase_delta, cur_dosage_ct, cur_dphase_ct, spgwp);
          if (unlikely(reterr)) {
            goto GparseFlush_ret_1;
          }
        }
      }
    }
  }
  while (0) {
  GparseFlush_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
 GparseFlush_ret_1:
  return reterr;
}


// AlwaysOmit: always present, always zero, could be either 2 or 3 input fields
// OmitWhenPresent: --iid-sid, FID always 0 when present
// CopyOr0: --iid-sid, present and sometimes nonzero when 3 input fields,
//          write "0" when only 2 input fields
// AlwaysCopy: always present, sometimes nonzero, could be either 2 or 3 inputs
ENUM_U31_DEF_START()
  kImportFidDelimModeAlwaysOmit,
  kImportFidDelimModeOmitWhenPresent,
  kImportFidDelimModeCopyOr0,
  kImportFidDelimModeAlwaysCopy
ENUM_U31_DEF_END(ImportFidDelimMode);

// NonexistOrOmit: don't generate SID column
//                 (caller now responsible for confirming during scanning pass
//                 that no IDs have 4+ parts)
// CopyOr0: present and sometimes nonzero when 3 input fields, write "0" when
//          only 2 input fields
// AlwaysCopy: always present, sometimes nonzero
ENUM_U31_DEF_START()
  kImportSidDelimModeNonexistOrOmit,
  kImportSidDelimModeCopyOr0,
  kImportSidDelimModeAlwaysCopy
ENUM_U31_DEF_END(ImportSidDelimMode);

typedef struct ImportSampleIdContextStruct {
  const char* const_fid;
  uint32_t const_fid_slen;
  uint32_t double_id;

  ImportFidDelimMode fid_delim_mode;
  ImportSidDelimMode sid_delim_mode;
  char id_delim;
} ImportSampleIdContext;

void InitImportSampleIdContext(const char* const_fid, ImportFlags import_flags, char id_delim, ImportSampleIdContext* isicp) {
  isicp->const_fid = nullptr;
  isicp->const_fid_slen = 0;
  if (const_fid && (!strequal_k_unsafe(const_fid, "0"))) {
    isicp->const_fid_slen = strlen(const_fid);
    isicp->const_fid = const_fid;
  }
  isicp->double_id = (import_flags / kfImportDoubleId) & 1;
  isicp->fid_delim_mode = kImportFidDelimModeAlwaysOmit;
  isicp->sid_delim_mode = kImportSidDelimModeNonexistOrOmit;
  isicp->id_delim = id_delim;
}

PglErr ImportSampleId(const char* input_id_iter, const char* input_id_end, const ImportSampleIdContext* isicp, char** write_iterp) {
  PglErr reterr = kPglRetSuccess;
  {
    const char id_delim = isicp->id_delim;
    char* write_iter = *write_iterp;
    if (id_delim) {
      const char* first_delim = S_CAST(const char*, memchr(input_id_iter, ctou32(id_delim), input_id_end - input_id_iter));
      assert(first_delim);  // previously validated
      if (unlikely(*input_id_iter == id_delim)) {
        snprintf(g_logbuf, kLogbufSize, "Error: '%c' at beginning of sample ID.\n", id_delim);
        goto ImportSampleId_ret_INCONSISTENT_INPUT_2;
      }
      if (unlikely(input_id_end[-1] == id_delim)) {
        snprintf(g_logbuf, kLogbufSize, "Error: '%c' at end of sample ID.\n", id_delim);
        goto ImportSampleId_ret_INCONSISTENT_INPUT_2;
      }
      const uint32_t first_part_slen = first_delim - input_id_iter;
      if (unlikely(first_part_slen > kMaxIdSlen)) {
        // strictly speaking, you could have e.g. a 20k char ID which
        // splits into a valid pair with the right delimiter, but you're
        // not supposed to have sample IDs anywhere near that length so
        // I'll classify this as MalformedInput.
        goto ImportSampleId_ret_MALFORMED_INPUT_LONG_ID;
      }
      const char* second_part_start = &(first_delim[1]);
      const char* second_part_end = AdvToDelimOrEnd(second_part_start, input_id_end, id_delim);
      const uint32_t second_part_slen = second_part_end - second_part_start;
      if (unlikely(second_part_slen > kMaxIdSlen)) {
        goto ImportSampleId_ret_MALFORMED_INPUT_LONG_ID;
      }
      const char* iid_start = second_part_start;
      uint32_t iid_slen = second_part_slen;
      if (second_part_end == input_id_end) {
        if (isicp->fid_delim_mode == kImportFidDelimModeOmitWhenPresent) {
          iid_start = input_id_iter;
          iid_slen = first_part_slen;
        } else if (isicp->fid_delim_mode == kImportFidDelimModeCopyOr0) {
          iid_start = input_id_iter;
          iid_slen = first_part_slen;
          write_iter = strcpya_k(write_iter, "0\t");
        } else if (isicp->fid_delim_mode == kImportFidDelimModeAlwaysCopy) {
          write_iter = memcpyax(write_iter, input_id_iter, first_part_slen, '\t');
        }
        write_iter = memcpya(write_iter, iid_start, iid_slen);
        if (isicp->sid_delim_mode == kImportSidDelimModeCopyOr0) {
          write_iter = strcpya_k(write_iter, "\t0");
        } else if (isicp->sid_delim_mode == kImportSidDelimModeAlwaysCopy) {
          *write_iter++ = '\t';
          write_iter = memcpya(write_iter, second_part_start, second_part_slen);
        }
      } else {
        if (unlikely(second_part_slen == 0)) {
          snprintf(g_logbuf, kLogbufSize, "Error: Consecutive instances of '%c' in sample ID.\n", id_delim);
          goto ImportSampleId_ret_INCONSISTENT_INPUT_2;
        }
        if (isicp->fid_delim_mode >= kImportFidDelimModeCopyOr0) {
          write_iter = memcpyax(write_iter, input_id_iter, first_part_slen, '\t');
        }
        write_iter = memcpya(write_iter, iid_start, iid_slen);
        if (isicp->sid_delim_mode >= kImportSidDelimModeCopyOr0) {
          const char* sid_start = &(second_part_end[1]);
          const uint32_t sid_slen = input_id_end - sid_start;
          if (unlikely(sid_slen > kMaxIdSlen)) {
            goto ImportSampleId_ret_MALFORMED_INPUT_LONG_ID;
          }
          *write_iter++ = '\t';
          write_iter = memcpya(write_iter, sid_start, sid_slen);
        }
      }
      if (unlikely((iid_slen == 1) && (iid_start[0] == '0'))) {
        logerrputs("Error: Sample ID induces an invalid IID of '0'.\n");
        goto ImportSampleId_ret_INCONSISTENT_INPUT;
      }
    } else {
      const uint32_t token_slen = input_id_end - input_id_iter;
      if (unlikely((*input_id_iter == '0') && (token_slen == 1))) {
        logerrputs("Error: Sample ID cannot be '0'.\n");
        goto ImportSampleId_ret_MALFORMED_INPUT;
      }
      if (unlikely(token_slen > kMaxIdSlen)) {
        goto ImportSampleId_ret_MALFORMED_INPUT_LONG_ID;
      }
      if (isicp->double_id) {
        write_iter = memcpyax(write_iter, input_id_iter, token_slen, '\t');
      } else if (isicp->const_fid) {
        write_iter = memcpyax(write_iter, isicp->const_fid, isicp->const_fid_slen, '\t');
      }
      write_iter = memcpya(write_iter, input_id_iter, token_slen);
    }
    *write_iterp = write_iter;
  }
  while (0) {
  ImportSampleId_ret_MALFORMED_INPUT_LONG_ID:
    logerrputs("Error: FIDs, IIDs, and SIDs are limited to " MAX_ID_SLEN_STR " characters.\n");
  ImportSampleId_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  ImportSampleId_ret_INCONSISTENT_INPUT_2:
    logerrputsb();
  ImportSampleId_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
  return reterr;
}

PglErr ImportIidFromSampleId(const char* input_id_iter, const char* input_id_end, const ImportSampleIdContext* isicp, const char** iid_start_ptr, uint32_t* iid_slen_ptr) {
  PglErr reterr = kPglRetSuccess;
  {
    const char id_delim = isicp->id_delim;
    const char* iid_start = input_id_iter;
    uint32_t iid_slen;
    if (id_delim) {
      const char* first_delim = S_CAST(const char*, memchr(input_id_iter, ctou32(id_delim), input_id_end - input_id_iter));
      if (unlikely(!first_delim)) {
        snprintf(g_logbuf, kLogbufSize, "Error: No '%c' in sample ID.\n", id_delim);
        goto ImportIidFromSampleId_ret_INCONSISTENT_INPUT_2;
      }
      const char* second_part_start = &(first_delim[1]);
      const char* second_part_end = AdvToDelimOrEnd(second_part_start, input_id_end, id_delim);
      uint32_t iid_second;
      if ((isicp->fid_delim_mode == kImportFidDelimModeAlwaysCopy) ||
          (isicp->fid_delim_mode == kImportFidDelimModeAlwaysOmit)) {
        iid_second = 1;
      } else {
        iid_second = (second_part_end != input_id_end);
      }
      const char* iid_end = first_delim;
      if (iid_second) {
        iid_start = second_part_start;
        iid_end = second_part_end;
      }
      iid_slen = iid_end - iid_start;
    } else {
      iid_slen = input_id_end - iid_start;
    }
    if (unlikely((iid_slen == 1) && (*iid_start == '0'))) {
      logerrputs("Error: IID cannot be '0'.\n");
      goto ImportIidFromSampleId_ret_MALFORMED_INPUT;
    }
    if (unlikely(iid_slen > kMaxIdSlen)) {
      goto ImportIidFromSampleId_ret_MALFORMED_INPUT_LONG_ID;
    }
    *iid_start_ptr = iid_start;
    *iid_slen_ptr = iid_slen;
  }
  while (0) {
  ImportIidFromSampleId_ret_MALFORMED_INPUT_LONG_ID:
    logerrputs("Error: FIDs, IIDs, and SIDs are limited to " MAX_ID_SLEN_STR " characters.\n");
  ImportIidFromSampleId_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  ImportIidFromSampleId_ret_INCONSISTENT_INPUT_2:
    logerrputsb();
    reterr = kPglRetInconsistentInput;
    break;
  }
  return reterr;
}

PglErr VcfSampleLine(const char* preexisting_psamname, const char* const_fid, MiscFlags misc_flags, ImportFlags import_flags, FamCol fam_cols, char id_delim, char idspace_to, char flag_char, char* sample_line_first_id, char* outname, char* outname_end, uint32_t* sample_ct_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* psamfile = nullptr;
  uintptr_t line_idx = 0;
  PglErr reterr = kPglRetSuccess;
  TextStream psam_txs;
  PreinitTextStream(&psam_txs);
  {
    ImportSampleIdContext isic;
    InitImportSampleIdContext(const_fid, import_flags, id_delim, &isic);
    uint32_t write_fid = 0;
    uint32_t write_sid = 0;
    // Alpha 2 changes:
    // * --id-delim can no longer be used with --const-fid/--double-id.  All
    //   samples must be multipart.
    // * New default is --const-fid 0, which causes no FID column to be written
    //   at all.
    if (id_delim) {
      if (id_delim != ' ') {
        char* sample_line_iter = strchrnul_n(sample_line_first_id, ' ');
        if (*sample_line_iter == ' ') {
          if (unlikely(!idspace_to)) {
            logerrputs("Error: VCF/BCF2 sample ID contains space(s).  Use --idspace-to to convert them\nto another character, or \"--id-delim ' '\" to interpret the spaces as FID/IID\nor IID/SID delimiters.\n");
            goto VcfSampleLine_ret_INCONSISTENT_INPUT;
          }
          do {
            *sample_line_iter = idspace_to;
            sample_line_iter = strchrnul_n(&(sample_line_iter[1]), ' ');
          } while (*sample_line_iter == ' ');
        }
      }
      const char* sample_line_iter = sample_line_first_id;
      const uint32_t iid_sid = (misc_flags / kfMiscIidSid) & 1;
      while (ctou32(sample_line_iter[0]) >= ' ') {
        const char* token_end = strchrnul_n(sample_line_iter, '\t');
        const char* first_delim = S_CAST(const char*, memchr(sample_line_iter, ctou32(id_delim), token_end - sample_line_iter));
        if (unlikely(!first_delim)) {
          snprintf(g_logbuf, kLogbufSize, "Error: No '%c' in sample ID.\n", id_delim);
          goto VcfSampleLine_ret_INCONSISTENT_INPUT_2;
        }
        const uint32_t first_slen = first_delim - sample_line_iter;
        const char* second_part_start = &(first_delim[1]);
        const char* maybe_second_part_end = S_CAST(const char*, memchr(second_part_start, ctou32(id_delim), token_end - second_part_start));
        if (maybe_second_part_end == nullptr) {
          if (!iid_sid) {
            write_fid |= (first_slen != 1) || (sample_line_iter[0] != '0');
          } else {
            const uint32_t sid_slen = token_end - second_part_start;
            write_sid |= (sid_slen != 1) || (second_part_start[0] != '0');
          }
        } else {
          write_fid |= (first_slen != 1) || (sample_line_iter[0] != '0');
          const char* sid_start = &(maybe_second_part_end[1]);
          const uint32_t sid_slen = token_end - sid_start;
          if (unlikely(memchr(sid_start, ctou32(id_delim), sid_slen))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Too many instances of --id-delim argument '%c' in sample ID.\n", id_delim);
            goto VcfSampleLine_ret_INCONSISTENT_INPUT_2;
          }
          write_sid |= (sid_slen != 1) || (sid_start[0] != '0');
        }
        if (*token_end != '\t') {
          break;
        }
        sample_line_iter = &(token_end[1]);
      }
      if (!iid_sid) {
        isic.fid_delim_mode = write_fid? kImportFidDelimModeAlwaysCopy : kImportFidDelimModeAlwaysOmit;
        isic.sid_delim_mode = write_sid? kImportSidDelimModeCopyOr0 : kImportSidDelimModeNonexistOrOmit;
      } else {
        isic.fid_delim_mode = write_fid? kImportFidDelimModeCopyOr0 : kImportFidDelimModeOmitWhenPresent;
        isic.sid_delim_mode = write_sid? kImportSidDelimModeAlwaysCopy : kImportSidDelimModeNonexistOrOmit;
      }
    } else if (isic.const_fid || isic.double_id) {
      write_fid = 1;
    }
    const char* sample_line_iter = sample_line_first_id;
    uint32_t sample_ct = 0;
    if (!preexisting_psamname) {
      snprintf(outname_end, kMaxOutfnameExtBlen, ".psam");
      if (unlikely(fopen_checked(outname, FOPEN_WB, &psamfile))) {
        goto VcfSampleLine_ret_OPEN_FAIL;
      }
      char* write_iter = g_textbuf;
      char* textbuf_flush = &(write_iter[kMaxMediumLine]);
      *write_iter++ = '#';
      if (write_fid) {
        write_iter = strcpya_k(write_iter, "FID\t");
      }
      write_iter = strcpya_k(write_iter, "IID");
      if (write_sid) {
        write_iter = strcpya_k(write_iter, "\tSID");
      }
      write_iter = strcpya_k(write_iter, "\tSEX");
      AppendBinaryEoln(&write_iter);
      while (ctou32(sample_line_iter[0]) >= ' ') {
        ++sample_ct;
        const char* token_end = NextPrespace(sample_line_iter);
        reterr = ImportSampleId(sample_line_iter, token_end, &isic, &write_iter);
        if (unlikely(reterr)) {
          goto VcfSampleLine_ret_1;
        }
        // PAT/MAT/PHENO1 not required in .psam file
        // SEX now included, so that --vcf + --out has the same effect as --vcf
        // + --make-pgen + --out
        write_iter = strcpya_k(write_iter, "\tNA");
        AppendBinaryEoln(&write_iter);
        if (unlikely(fwrite_ck(textbuf_flush, psamfile, &write_iter))) {
          goto VcfSampleLine_ret_WRITE_FAIL;
        }
        if (*token_end != '\t') {
          break;
        }
        sample_line_iter = &(token_end[1]);
      }
      if (unlikely(fclose_flush_null(textbuf_flush, write_iter, &psamfile))) {
        goto VcfSampleLine_ret_WRITE_FAIL;
      }
    } else {
      // check consistency of IIDs between VCF and .psam file.
      reterr = SizeAndInitTextStream(preexisting_psamname, bigstack_left(), 1, &psam_txs);
      if (unlikely(reterr)) {
        goto VcfSampleLine_ret_TSTREAM_FAIL;
      }
      uint32_t sample_line_eoln = IsEoln(sample_line_iter[0]);
      char* psam_line_start;
      do {
        ++line_idx;
        psam_line_start = TextGet(&psam_txs);
        if (!psam_line_start) {
          if (TextStreamErrcode2(&psam_txs, &reterr)) {
            goto VcfSampleLine_ret_TSTREAM_FAIL;
          }
          if (unlikely(!sample_line_eoln)) {
            snprintf(g_logbuf, kLogbufSize, "Error: --%ccf file contains more sample IDs than %s.\n", flag_char, preexisting_psamname);
            goto VcfSampleLine_ret_INCONSISTENT_INPUT_WW;
          }
          *sample_ct_ptr = 0;
          goto VcfSampleLine_ret_1;
        }
      } while ((psam_line_start[0] == '#') && (!tokequal_k(&(psam_line_start[1]), "FID")) && (!tokequal_k(&(psam_line_start[1]), "IID")));
      uint32_t fid_exists;
      if (psam_line_start[0] == '#') {
        // only check for matching IIDs for now.
        fid_exists = (psam_line_start[1] == 'F');
        // bugfix (12 Apr 2018): forgot to skip this line
        ++line_idx;
        psam_line_start = TextGet(&psam_txs);
      } else {
        fid_exists = fam_cols & kfFamCol1;
      }
      for (; psam_line_start; ++line_idx, psam_line_start = TextGet(&psam_txs)) {
        if (unlikely(psam_line_start[0] == '#')) {
          snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of %s starts with a '#'. (This is only permitted before the first nonheader line, and a #FID/IID header line is present it must denote the end of the header block.)\n", line_idx, preexisting_psamname);
          goto VcfSampleLine_ret_MALFORMED_INPUT_WW;
        }
        const char* psam_iid_start = psam_line_start;
        if (fid_exists) {
          psam_iid_start = FirstNonTspace(CurTokenEnd(psam_iid_start));
          if (unlikely(IsEolnKns(*psam_iid_start))) {
            goto VcfSampleLine_ret_MISSING_TOKENS;
          }
        }
        if (unlikely(sample_line_eoln)) {
          snprintf(g_logbuf, kLogbufSize, "Error: --%ccf file contains fewer sample IDs than %s.\n", flag_char, preexisting_psamname);
          goto VcfSampleLine_ret_INCONSISTENT_INPUT_WW;
        }
        ++sample_ct;
        const char* sample_line_token_end = NextPrespace(sample_line_iter);
        const char* sample_line_iid_start;
        uint32_t sample_line_iid_slen;
        reterr = ImportIidFromSampleId(sample_line_iter, sample_line_token_end, &isic, &sample_line_iid_start, &sample_line_iid_slen);
        if (unlikely(reterr)) {
          goto VcfSampleLine_ret_1;
        }
        if (unlikely(!memequal(sample_line_iid_start, psam_iid_start, sample_line_iid_slen) || (ctou32(psam_iid_start[sample_line_iid_slen]) > 32))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Mismatched IDs between --%ccf file and %s.\n", flag_char, preexisting_psamname);
          goto VcfSampleLine_ret_INCONSISTENT_INPUT_WW;
        }
        sample_line_eoln = (*sample_line_token_end != '\t');
        sample_line_iter = &(sample_line_token_end[1]);
      }
      if (unlikely(TextStreamErrcode2(&psam_txs, &reterr))) {
        goto VcfSampleLine_ret_TSTREAM_FAIL;
      }
      if (unlikely(!sample_line_eoln)) {
        snprintf(g_logbuf, kLogbufSize, "Error: --%ccf file contains more sample IDs than %s.\n", flag_char, preexisting_psamname);
        goto VcfSampleLine_ret_INCONSISTENT_INPUT_WW;
      }
    }
    *sample_ct_ptr = sample_ct;
  }
  while (0) {
  VcfSampleLine_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  VcfSampleLine_ret_TSTREAM_FAIL:
    TextStreamErrPrint(preexisting_psamname, &psam_txs);
    break;
  VcfSampleLine_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  VcfSampleLine_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetMalformedInput;
    break;
  VcfSampleLine_ret_MISSING_TOKENS:
    logerrprintfww("Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, preexisting_psamname);
    reterr = kPglRetMalformedInput;
    break;
  VcfSampleLine_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
  VcfSampleLine_ret_INCONSISTENT_INPUT_2:
    logerrputsb();
  VcfSampleLine_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 VcfSampleLine_ret_1:
  CleanupTextStream2(preexisting_psamname, &psam_txs, &reterr);
  fclose_cond(psamfile);
  BigstackReset(bigstack_mark);
  return reterr;
}

uint32_t VcfIsHetShort(const char* first_gchar_ptr, VcfHalfCall vcf_half_call) {
  // '.' == ascii 46, '0' == ascii 48
  // if kVcfHalfCallReference, .|0 is not phased, but .|1 is
  const uint32_t first_gchar = ctou32(first_gchar_ptr[0]);
  const uint32_t second_gchar = ctou32(first_gchar_ptr[2]);
  return (first_gchar != second_gchar) && (((first_gchar != 46) && (second_gchar != 46)) || ((vcf_half_call == kVcfHalfCallReference) && ((first_gchar > 48) || (second_gchar > 48))));
}

uint32_t VcfIsHetLong(const char* first_allele_idx_str_iter, const char* second_allele_idx_str_iter, VcfHalfCall vcf_half_call) {
  uint32_t first_allele_idx = ctou32(*first_allele_idx_str_iter++);
  uint32_t second_allele_idx = ctou32(*second_allele_idx_str_iter++);
  if (first_allele_idx == 46) {
    return (vcf_half_call == kVcfHalfCallReference) && (second_allele_idx > 48);
  }
  if (second_allele_idx == 46) {
    return (vcf_half_call == kVcfHalfCallReference) && (first_allele_idx > 48);
  }
  first_allele_idx -= 48;
  for (; ; ++first_allele_idx_str_iter) {
    const uint32_t next_digit = ctou32(*first_allele_idx_str_iter) - 48;
    if (next_digit >= 10) {
      break;
    }
    first_allele_idx = first_allele_idx * 10 + next_digit;
  }
  second_allele_idx -= 48;
  for (; ; ++second_allele_idx_str_iter) {
    const uint32_t next_digit = ctou32(*second_allele_idx_str_iter) - 48;
    if (next_digit >= 10) {
      break;
    }
    second_allele_idx = second_allele_idx * 10 + next_digit;
  }
  return (first_allele_idx != second_allele_idx);
}

uint32_t GetVcfFormatPosition(const char* __restrict needle, const char* format_start, const char* format_end, uint32_t needle_slen) {
  // no longer assumes first field is GT
  const char* token_start = format_start;
  for (uint32_t field_idx = 0; ; ++field_idx) {
    const char* token_end = AdvToDelimOrEnd(token_start, format_end, ':');
    if ((S_CAST(uintptr_t, token_end - token_start) == needle_slen) && memequal(token_start, needle, needle_slen)) {
      return field_idx;
    }
    if (token_end == format_end) {
      return UINT32_MAX;
    }
    token_start = &(token_end[1]);
  }
}

uint32_t VcfQualScanInit1(const char* format_start, const char* format_end, int32_t vcf_min_gq, int32_t vcf_min_dp, int32_t vcf_max_dp, STD_ARRAY_REF(uint32_t, 2) qual_field_idxs) {
  uint32_t qual_exists = 0;
  qual_field_idxs[0] = UINT32_MAX;
  qual_field_idxs[1] = UINT32_MAX;
  if (vcf_min_gq >= 0) {
    qual_field_idxs[0] = GetVcfFormatPosition("GQ", format_start, format_end, 2);
    qual_exists = (qual_field_idxs[0] != UINT32_MAX);
  }
  if ((vcf_min_dp >= 0) || (vcf_max_dp != 0x7fffffff)) {
    qual_field_idxs[1] = GetVcfFormatPosition("DP", format_start, format_end, 2);
    if (qual_field_idxs[1] != UINT32_MAX) {
      qual_exists = 1;
    }
  }
  return qual_exists;
}

uint32_t VcfQualScanInit2(STD_ARRAY_KREF(uint32_t, 2) qual_field_idxs, STD_ARRAY_KREF(int32_t, 2) qual_mins, STD_ARRAY_KREF(int32_t, 2) qual_maxs, STD_ARRAY_REF(uint32_t, 2) qual_field_skips, STD_ARRAY_REF(int32_t, 2) qual_line_mins, STD_ARRAY_REF(int32_t, 2) qual_line_maxs) {
  // handcoded for now, but can be switched to a std::sort call if necessary
  const uint32_t gq_field_idx = qual_field_idxs[0];
  uint32_t dp_field_idx = qual_field_idxs[1];
  uint32_t qual_field_ct = 0;
  if (dp_field_idx < gq_field_idx) {
    qual_field_skips[0] = dp_field_idx;
    qual_line_mins[0] = qual_mins[1];
    qual_line_maxs[0] = qual_maxs[1];
    qual_field_ct = 1;
    dp_field_idx = UINT32_MAX;
  }
  if (gq_field_idx != UINT32_MAX) {
    qual_field_skips[qual_field_ct] = gq_field_idx;
    qual_line_mins[qual_field_ct] = qual_mins[0];
    qual_line_maxs[qual_field_ct] = 0x7fffffff;
    ++qual_field_ct;
    if (dp_field_idx != UINT32_MAX) {
      qual_field_skips[qual_field_ct] = dp_field_idx;
      qual_line_mins[qual_field_ct] = qual_mins[1];
      qual_line_maxs[qual_field_ct] = qual_maxs[1];
      ++qual_field_ct;
    }
  }
  if (qual_field_ct == 2) {
    qual_field_skips[1] -= qual_field_skips[0];
  }
  return qual_field_ct;
}

// These values are constant over a single VCF line, and relevant to parsing
// inner loop(s).  VcfImportBaseContext covers what's needed for dosageless
// import, while VcfImportContext covers the general case.
typedef struct VcfImportContextBaseStruct {
  // constant across file
  uint32_t sample_ct;
  VcfHalfCall halfcall_mode;
  uint32_t error_on_polyploid;

  // line-specific
  uint32_t gt_exists;
  STD_ARRAY_DECL(uint32_t, 2, qual_field_skips);
  STD_ARRAY_DECL(int32_t, 2, qual_line_mins);
  STD_ARRAY_DECL(int32_t, 2, qual_line_maxs);
  uint32_t qual_field_ct;  // must be set to zero if no qual fields
} VcfImportBaseContext;

typedef struct VcfImportContext {
  VcfImportBaseContext vibc;

  // constant across file
  uint32_t dosage_is_gp;
  uint32_t dosage_erase_halfdist;
  double import_dosage_certainty;

  // line-specific
  uint32_t dosage_field_idx;
  uint32_t hds_field_idx;
} VcfImportContext;

// returns 1 if a quality check failed
// assumes either 1 or 2 qual fields, otherwise change this to a loop
uint32_t VcfCheckQuals(STD_ARRAY_KREF(uint32_t, 2) qual_field_skips, STD_ARRAY_KREF(int32_t, 2) qual_line_mins, STD_ARRAY_KREF(int32_t, 2) qual_line_maxs, const char* gtext_iter, const char* gtext_end, uint32_t qual_field_ct) {
  const uint32_t skip0 = qual_field_skips[0];  // this can now be zero
  if (skip0) {
    gtext_iter = AdvToNthDelimChecked(gtext_iter, gtext_end, skip0, ':');
    if (!gtext_iter) {
      return 0;
    }
    ++gtext_iter;
  }
  int32_t ii;
  if ((!ScanInt32(gtext_iter, &ii)) && ((ii < qual_line_mins[0]) || (ii > qual_line_maxs[0]))) {
    return 1;
  }
  if (qual_field_ct == 1) {
    return 0;
  }
  gtext_iter = AdvToNthDelimChecked(gtext_iter, gtext_end, qual_field_skips[1], ':');
  if (!gtext_iter) {
    return 0;
  }
  ++gtext_iter;
  return (!ScanInt32(gtext_iter, &ii)) && ((ii < qual_line_mins[1]) || (ii > qual_line_maxs[1]));
}

// kDosageParseForceMissing = --import-dosage-certainty filter applied
ENUM_U31_DEF_START()
  kDosageParseOk,
  kDosageParsePolyploid,
  kDosageParseMissing,
  kDosageParseForceMissing
ENUM_U31_DEF_END(DosageParseResult);

BoolErr ParseVcfBiallelicGp(const char* gp_iter, uint32_t is_haploid, double import_dosage_certainty, DosageParseResult* dpr_ptr, double* alt_dosage_ptr) {
  // P(0/0), P(0/1), P(1/1), etc.
  // assumes dpr initialized to kDosageParseOk
  double prob_0alt;
  gp_iter = ScanadvDouble(gp_iter, &prob_0alt);
  if (unlikely((!gp_iter) || (prob_0alt < 0.0) || (prob_0alt > 1.0) || (*gp_iter != ','))) {
    return 1;
  }
  double prob_1alt;
  gp_iter = ScanadvDouble(&(gp_iter[1]), &prob_1alt);
  if (unlikely((!gp_iter) || (prob_1alt < 0.0) || (prob_1alt > 1.0))) {
    return 1;
  }
  if (is_haploid) {
    const double denom = prob_0alt + prob_1alt;
    if (denom <= 2 * import_dosage_certainty) {
      if ((prob_0alt <= import_dosage_certainty) && (prob_1alt <= import_dosage_certainty)) {
        *dpr_ptr = kDosageParseForceMissing;
        return 1;
      }
    }
    *alt_dosage_ptr = 2 * prob_1alt / denom;
    return 0;
  }
  double prob_2alt;
  if (unlikely((*gp_iter != ',') || (!ScanadvDouble(&(gp_iter[1]), &prob_2alt)) || (prob_2alt < 0.0) || (prob_2alt > 1.0))) {
    return 1;
  }
  const double denom = prob_0alt + prob_1alt + prob_2alt;
  if (denom <= 3 * import_dosage_certainty) {
    if ((prob_0alt <= import_dosage_certainty) && (prob_1alt <= import_dosage_certainty) && (prob_2alt <= import_dosage_certainty)) {
      // force-missing
      // ok to use <= since we multiplied by (1 - epsilon)
      // during command-line parsing.  this lets us avoid
      // special-casing denom=0.
      *dpr_ptr = kDosageParseForceMissing;
      return 1;  // not really an error
    }
  }
  *alt_dosage_ptr = (prob_1alt + 2 * prob_2alt) / denom;
  return 0;
}

BoolErr ParseVcfBiallelicDosage(const char* gtext_iter, const char* gtext_end, uint32_t dosage_field_idx, uint32_t is_haploid, uint32_t dosage_is_gp, double import_dosage_certainty, DosageParseResult* dpr_ptr, uint32_t* dosage_int_ptr) {
  // assumes dosage_field_idx != UINT32_MAX
  // assumes dpr initialized to kDosageParseOk
  // returns 1 if missing OR parsing error.  error: dpr still kDosageParseOk.
  if (dosage_field_idx) {
    gtext_iter = AdvToNthDelimChecked(gtext_iter, gtext_end, dosage_field_idx, ':');
    if (!gtext_iter) {
      *dpr_ptr = kDosageParseMissing;
      return 1;
    }
    ++gtext_iter;
  }
  if ((gtext_iter[0] == '?') || ((gtext_iter[0] == '.') && (ctou32(gtext_iter[1]) - 48 >= 10))) {
    // missing field (dot/'?' followed by non-digit)
    // could enforce gtext_iter[1] == colon, comma, etc.?
    *dpr_ptr = kDosageParseMissing;
    return 1;
  }
  double alt_dosage;
  if (dosage_is_gp) {
    if (ParseVcfBiallelicGp(gtext_iter, is_haploid, import_dosage_certainty, dpr_ptr, &alt_dosage)) {
      return 1;
    }
  } else {
    if (unlikely((!ScanadvDouble(gtext_iter, &alt_dosage)) || (alt_dosage < 0.0))) {
      return 1;
    }
    if (is_haploid) {
      // possible todo: allow this to be suppressed (maybe upstream of this
      // function); 1000 Genomes phase 1 haploid dosages are still on 0..2
      // scale
      // right now the best approach for importing those files is commenting
      // out this line and recompiling...
      if (import_dosage_certainty != 0.0) {
        // quasi-bugfix (19 Feb 2019): dosage=DS import should respect
        // --import-dosage-certainty
        if (((1.0 - alt_dosage) <= import_dosage_certainty) && (alt_dosage <= import_dosage_certainty)) {
          *dpr_ptr = kDosageParseForceMissing;
          return 1;
        }
      }
      alt_dosage *= 2;
    } else {
      if (import_dosage_certainty != 0.0) {
        const double dist_from_1 = fabs(1.0 - alt_dosage);
        if ((1.0 - dist_from_1 <= import_dosage_certainty) && (dist_from_1 <= import_dosage_certainty)) {
          *dpr_ptr = kDosageParseForceMissing;
          return 1;
        }
      }
    }
    if (unlikely(alt_dosage > 2.0)) {
      return 1;
    }
  }
  *dosage_int_ptr = S_CAST(int32_t, alt_dosage * kDosageMid + 0.5);
  return 0;
}

BoolErr ParseVcfBiallelicHds(const char* gtext_iter, const char* gtext_end, uint32_t dosage_field_idx, uint32_t hds_field_idx, uint32_t is_haploid, uint32_t dosage_is_gp, double import_dosage_certainty, DosageParseResult* dpr_ptr, uint32_t* dosage_int_ptr, int32_t* cur_dphase_delta_ptr, uint32_t* hds_valid_ptr) {
  // assumes dpr initialized to kDosageParseOk
  // assumes cur_dphase_delta initialized to 0
  // assumes hds_valid initialized to 0
  // assumes dosage_field_idx != UINT32_MAX and/or hds_field_idx != UINT32_MAX
  // returns 1 if missing OR parsing error.  error: dpr still kDosageParseOk,
  // unless kDosageParsePolyploid.

  if (hds_field_idx != UINT32_MAX) {
    // search for HDS first, then DS
    const char* hds_gtext_iter = gtext_iter;
    if (hds_field_idx) {
      hds_gtext_iter = AdvToNthDelimChecked(gtext_iter, gtext_end, hds_field_idx, ':');
      if (!hds_gtext_iter) {
        goto ParseVcfBiallelicHdsSkip;
      }
      ++hds_gtext_iter;
    }
    if ((hds_gtext_iter[0] != '?') && ((hds_gtext_iter[0] != '.') || (ctou32(hds_gtext_iter[1]) - 48 < 10))) {
      double dosage1;
      // tried implementing fast path for '0', '1', '0,0', '1,1'; only helped
      // ~10% in HDS-force case
      hds_gtext_iter = ScanadvDouble(hds_gtext_iter, &dosage1);
      if (unlikely((!hds_gtext_iter) || (dosage1 < 0.0) || (dosage1 > 1.0))) {
        return 1;
      }
      // if hds_valid and (cur_dphase_delta == 0), caller should override
      // hardcall-phase
      *hds_valid_ptr = 1;
      if (*hds_gtext_iter != ',') {
        if (import_dosage_certainty != 0.0) {
          if ((1.0 - dosage1 <= import_dosage_certainty) && (dosage1 <= import_dosage_certainty)) {
            *dpr_ptr = kDosageParseForceMissing;
            return 1;
          }
        }
        *dosage_int_ptr = S_CAST(int32_t, dosage1 * kDosageMax + 0.5);
        return 0;
      }
      ++hds_gtext_iter;
      // don't permit HDS half-calls
      double dosage2;
      hds_gtext_iter = ScanadvDouble(hds_gtext_iter, &dosage2);
      if (unlikely((!hds_gtext_iter) || (dosage2 < 0.0) || (dosage2 > 1.0))) {
        return 1;
      }
      if (unlikely(*hds_gtext_iter == ',')) {
        *dpr_ptr = kDosageParsePolyploid;
        return 1;
      }
      const double dosage_sum = dosage1 + dosage2;
      if (import_dosage_certainty != 0.0) {
        // Assume maximal het probability.
        const double dist_from_1 = fabs(1.0 - dosage_sum);
        if ((1.0 - dist_from_1 <= import_dosage_certainty) && (dist_from_1 <= import_dosage_certainty)) {
          *dpr_ptr = kDosageParseForceMissing;
          return 1;
        }
      }

      // force this to be nonnegative, since static_cast<int32_t> rounds
      // negative numbers toward zero
      const double dosage_diffp1 = 1.0 + dosage1 - dosage2;

      *dosage_int_ptr = S_CAST(int32_t, dosage_sum * kDosageMid + 0.5);
      *cur_dphase_delta_ptr = S_CAST(int32_t, dosage_diffp1 * kDosageMid + 0.5) - kDosageMid;
      return 0;
    }
  ParseVcfBiallelicHdsSkip:
    if (dosage_field_idx == UINT32_MAX) {
      *dpr_ptr = kDosageParseMissing;
      return 1;
    }
  }
  return ParseVcfBiallelicDosage(gtext_iter, gtext_end, dosage_field_idx, is_haploid, dosage_is_gp, import_dosage_certainty, dpr_ptr, dosage_int_ptr);
}

ENUM_U31_DEF_START()
  kVcfParseOk,
  kVcfParseMissingTokens,
  kVcfParseInvalidGt,
  kVcfParseHalfCallError,
  kVcfParseInvalidDosage,
  kVcfParsePolyploidError
ENUM_U31_DEF_END(VcfParseErr);

VcfParseErr VcfScanBiallelicHdsLine(const VcfImportContext* vicp, const char* format_end, uint32_t* phase_or_dosage_found_ptr, char** line_iter_ptr) {
  // Either just DS, or DS+HDS.
  // only need to find phase *or* dosage
  const uint32_t gt_exists = vicp->vibc.gt_exists;
  const uint32_t sample_ct = vicp->vibc.sample_ct;
  const uint32_t hds_field_idx = vicp->hds_field_idx;
  // special case: if there is no GT or HDS field, there can't be any phase
  // information; we're only here to check for a non-integer dosage.  So if
  // there either (i) aren't any '.' characters, or (ii) every second character
  // is a '\t', we can early-exit.
  if ((!gt_exists) && (hds_field_idx == UINT32_MAX)) {
    const char* dot_or_lf_ptr = S_CAST(const char*, rawmemchr2(format_end, '.', '\n'));
    if (*dot_or_lf_ptr == '\n') {
      *line_iter_ptr = K_CAST(char*, dot_or_lf_ptr);
      return kVcfParseOk;
    }
    const char* lf_ptr = AdvToDelim(dot_or_lf_ptr, '\n');
    if (lf_ptr[-1] == '\r') {
      --lf_ptr;
    }
    const uintptr_t body_len = lf_ptr - format_end;
    if (body_len <= 2 * sample_ct) {
      if (unlikely(body_len < 2 * sample_ct)) {
        return kVcfParseMissingTokens;
      }
      // TODO: check for alternate tabs
    }
  }
  const VcfHalfCall halfcall_mode = vicp->vibc.halfcall_mode;
  const uint32_t error_on_polyploid = vicp->vibc.error_on_polyploid;
  STD_ARRAY_KREF(uint32_t, 2) qual_field_skips = vicp->vibc.qual_field_skips;
  STD_ARRAY_KREF(int32_t, 2) qual_line_mins = vicp->vibc.qual_line_mins;
  STD_ARRAY_KREF(int32_t, 2) qual_line_maxs = vicp->vibc.qual_line_maxs;
  const uint32_t qual_field_ct = vicp->vibc.qual_field_ct;

  const uint32_t dosage_erase_halfdist = vicp->dosage_erase_halfdist;
  const double import_dosage_certainty = vicp->import_dosage_certainty;
  const uint32_t dosage_field_idx = vicp->dosage_field_idx;

  const char* dosagescan_iter = format_end;
  uint32_t cur_gt_phased = 0;
  uint32_t is_haploid = 0;
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    const char* cur_gtext_start = ++dosagescan_iter;
    const char* cur_gtext_end = FirstPrespace(dosagescan_iter);
    if (unlikely((*cur_gtext_end != '\t') && (sample_idx + 1 != sample_ct))) {
      return kVcfParseMissingTokens;
    }
    dosagescan_iter = cur_gtext_end;
    if (qual_field_ct && VcfCheckQuals(qual_field_skips, qual_line_mins, qual_line_maxs, cur_gtext_start, cur_gtext_end, qual_field_ct)) {
      continue;
    }
    if (gt_exists) {
      cur_gt_phased = (cur_gtext_start[1] == '|');
      is_haploid = (cur_gtext_start[1] != '/') && (!cur_gt_phased);
    }
    DosageParseResult dpr = kDosageParseOk;
    int32_t cur_dphase_delta = 0;
    uint32_t hds_valid = 0;
    uint32_t dosage_int;
    if (ParseVcfBiallelicHds(cur_gtext_start, cur_gtext_end, dosage_field_idx, hds_field_idx, is_haploid, 0, import_dosage_certainty, &dpr, &dosage_int, &cur_dphase_delta, &hds_valid)) {
      if (unlikely(!dpr)) {
        return kVcfParseInvalidDosage;
      }
      if (dpr == kDosageParsePolyploid) {
        if (unlikely(error_on_polyploid)) {
          return kVcfParsePolyploidError;
        }
        dpr = kDosageParseForceMissing;
      } else if ((dpr != kDosageParseForceMissing) && cur_gt_phased && VcfIsHetShort(cur_gtext_start, halfcall_mode)) {
        goto VcfScanBiallelicHdsLine_found;
      }
    } else if (cur_dphase_delta) {
      const uint32_t dphase_side1 = dosage_int + cur_dphase_delta;
      const uint32_t dphase_side2 = dosage_int - cur_dphase_delta;
      const uint32_t dphase_halfdist1 = DphaseHalfdist(dphase_side1);
      const uint32_t dphase_halfdist2 = DphaseHalfdist(dphase_side2);
      const uint32_t dphase_erase_halfdist = dosage_erase_halfdist + kDosage4th;
      if ((dphase_halfdist1 < dphase_erase_halfdist) ||
          (dphase_halfdist2 < dphase_erase_halfdist) ||
          (((dphase_side1 + kDosageMid) ^ (dphase_side2 + kDosageMid)) & kDosageMax)) {
        goto VcfScanBiallelicHdsLine_found;
      }
    } else {
      const uint32_t cur_halfdist = BiallelicDosageHalfdist(dosage_int);
      if ((cur_halfdist < dosage_erase_halfdist) ||
          ((!hds_valid) && cur_gt_phased && VcfIsHetShort(cur_gtext_start, halfcall_mode))) {
        goto VcfScanBiallelicHdsLine_found;
      }
    }
  }
  while (0) {
  VcfScanBiallelicHdsLine_found:
    *phase_or_dosage_found_ptr = 1;
  }
  *line_iter_ptr = K_CAST(char*, dosagescan_iter);
  return kVcfParseOk;
}

// either HDS, or both dosage and phase flags set
VcfParseErr VcfConvertPhasedBiallelicDosageLine(const VcfImportContext* vicp, const char* linebuf_iter, uintptr_t* genovec, uintptr_t* phasepresent, uintptr_t* phaseinfo, uintptr_t* dosage_present, uintptr_t* dphase_present, Dosage** dosage_main_iter_ptr, SDosage** dphase_delta_iter_ptr) {
  const uint32_t sample_ct = vicp->vibc.sample_ct;
  const uint32_t sample_ctl2_m1 = (sample_ct - 1) / kBitsPerWordD2;
  const VcfHalfCall halfcall_mode = vicp->vibc.halfcall_mode;
  const uint32_t error_on_polyploid = vicp->vibc.error_on_polyploid;
  const uint32_t gt_exists = vicp->vibc.gt_exists;
  STD_ARRAY_KREF(uint32_t, 2) qual_field_skips = vicp->vibc.qual_field_skips;
  STD_ARRAY_KREF(int32_t, 2) qual_line_mins = vicp->vibc.qual_line_mins;
  STD_ARRAY_KREF(int32_t, 2) qual_line_maxs = vicp->vibc.qual_line_maxs;
  const uint32_t qual_field_ct = vicp->vibc.qual_field_ct;

  Halfword* phasepresent_alias = R_CAST(Halfword*, phasepresent);
  Halfword* phaseinfo_alias = R_CAST(Halfword*, phaseinfo);
  const uint32_t dosage_is_gp = vicp->dosage_is_gp;
  const uint32_t dosage_erase_halfdist = vicp->dosage_erase_halfdist;
  const uint32_t dphase_erase_halfdist = dosage_erase_halfdist + kDosage4th;
  const double import_dosage_certainty = vicp->import_dosage_certainty;
  const uint32_t dosage_field_idx = vicp->dosage_field_idx;
  const uint32_t hds_field_idx = vicp->hds_field_idx;
  Halfword* dosage_present_alias = R_CAST(Halfword*, dosage_present);
  Halfword* dphase_present_alias = R_CAST(Halfword*, dphase_present);
  Dosage* dosage_main_iter = *dosage_main_iter_ptr;
  SDosage* dphase_delta_iter = *dphase_delta_iter_ptr;
  uint32_t inner_loop_last = kBitsPerWordD2 - 1;
  uint32_t cur_gt_phased = 0;
  uint32_t is_haploid = 0;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= sample_ctl2_m1) {
      if (widx > sample_ctl2_m1) {
        if (widx % 2) {
          phasepresent_alias[widx] = 0;
          // phaseinfo doesn't matter
          dosage_present_alias[widx] = 0;
          dphase_present_alias[widx] = 0;
        }
        break;
      }
      inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
    }
    uintptr_t genovec_word = 0;
    uint32_t phasepresent_hw = 0;
    uint32_t phaseinfo_hw = 0;
    uint32_t dosage_present_hw = 0;
    uint32_t dphase_present_hw = 0;
    for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
      const char* cur_gtext_end = FirstPrespace(linebuf_iter);
      if (unlikely((*cur_gtext_end != '\t') && ((sample_idx_lowbits != inner_loop_last) || (widx != sample_ctl2_m1)))) {
        return kVcfParseMissingTokens;
      }
      uintptr_t cur_geno = 3;
      if (!(qual_field_ct && VcfCheckQuals(qual_field_skips, qual_line_mins, qual_line_maxs, linebuf_iter, cur_gtext_end, qual_field_ct))) {
        // We now parse dosage first.  Only care about the hardcall if
        // (i) there's no dosage, or
        // (ii) there's no phased dosage, and it's a phased het.
        if (gt_exists) {
          cur_gt_phased = (linebuf_iter[1] == '|');
          is_haploid = (linebuf_iter[1] != '/') && (!cur_gt_phased);
        }
        const uint32_t shifted_bit = 1U << sample_idx_lowbits;
        DosageParseResult dpr = kDosageParseOk;
        int32_t cur_dphase_delta = 0;
        uint32_t hds_valid = 0;
        uint32_t dosage_int;
        if (!ParseVcfBiallelicHds(linebuf_iter, cur_gtext_end, dosage_field_idx, hds_field_idx, is_haploid, dosage_is_gp, import_dosage_certainty, &dpr, &dosage_int, &cur_dphase_delta, &hds_valid)) {
          if (hds_valid) {
            const uint32_t dphase_halfdist1 = DphaseHalfdist(dosage_int + cur_dphase_delta);
            const uint32_t dphase_halfdist2 = DphaseHalfdist(dosage_int - cur_dphase_delta);
            if ((dphase_halfdist1 < dphase_erase_halfdist) || (dphase_halfdist2 < dphase_erase_halfdist)) {
              // No need to fill cur_geno here, since it'll get corrected by
              // --hard-call-threshold.
              dosage_present_hw |= shifted_bit;
              *dosage_main_iter++ = dosage_int;
              if (cur_dphase_delta) {
                dphase_present_hw |= shifted_bit;
                *dphase_delta_iter++ = cur_dphase_delta;
              }
            } else {
              // Not saving dosage, since it's too close to an integer
              // (--dosage-erase-threshold).  Just directly synthesize the
              // hardcall we need.
              cur_geno = (dosage_int + kDosage4th) / kDosageMid;
              if (cur_geno == 1) {
                // Since dphase_erase_halfdist >= 8193, dphase_halfdist1 and
                // dphase_halfdist2 are both in [0, 8191] or [24577, 32768], so
                // it's always appropriate to save hardcall-phase.
                phasepresent_hw |= shifted_bit;
                if (cur_dphase_delta > 0) {
                  phaseinfo_hw |= shifted_bit;
                }
              }
            }
            goto VcfConvertPhasedBiallelicDosageLine_geno_done;
          }
          // defer handling of unphased dosage
        } else if (unlikely(!dpr)) {
          return kVcfParseInvalidDosage;
        } else if (dpr == kDosageParseForceMissing) {
          // bugfix (20 Feb 2019): if forced missing due to
          // --import-dosage-certainty, gt_exists must be ignored
          goto VcfConvertPhasedBiallelicDosageLine_geno_done;
        } else if (dpr == kDosageParsePolyploid) {
          if (unlikely(error_on_polyploid)) {
            return kVcfParsePolyploidError;
          }
          goto VcfConvertPhasedBiallelicDosageLine_geno_done;
        }
        if (gt_exists) {
          cur_geno = ctow(*linebuf_iter) - 48;
          // this nasty code should be spelled out in as few places as possible
          if (cur_geno <= 1) {
            if (is_haploid) {
              cur_geno *= 2;
            } else {
              const char polyploid_char = linebuf_iter[3];
              // Update (26 Sep 2023): "x/y/." is no longer treated as "x/y".
              // If someone's depending on that behavior, add another
              // --polyploid mode for it.
              if ((polyploid_char == '/') || (polyploid_char == '|')) {
                if (unlikely(error_on_polyploid)) {
                  return kVcfParsePolyploidError;
                }
                cur_geno = 3;
              } else {
                const uintptr_t second_allele_idx = ctow(linebuf_iter[2]) - 48;
                if (second_allele_idx <= 1) {
                  cur_geno += second_allele_idx;
                  if (cur_gt_phased && (cur_geno == 1)) {
                    phasepresent_hw |= shifted_bit;
#ifdef USE_AVX2
                    // See comment in VcfConvertPhasedBiallelicLine().
                    if (!second_allele_idx) {
                      // 1|0
                      phaseinfo_hw |= shifted_bit;
                    }
#else
                    phaseinfo_hw |= shifted_bit & (second_allele_idx - 1);
#endif
                  }
                } else if (unlikely(second_allele_idx != (~k0LU) * 2)) {
                  // not '.'
                  return kVcfParseInvalidGt;
                } else if (halfcall_mode == kVcfHalfCallMissing) {
                  cur_geno = 3;
                } else if (unlikely(halfcall_mode == kVcfHalfCallError)) {
                  return kVcfParseHalfCallError;
                } else {
                  // kVcfHalfCallHaploid, kVcfHalfCallReference
                  cur_geno <<= halfcall_mode;
                }
              }
            }
          } else if (unlikely(cur_geno != (~k0LU) * 2)) {
            // not '.'
            return kVcfParseInvalidGt;
          } else if (halfcall_mode != kVcfHalfCallMissing) {
            const char second_allele_char = linebuf_iter[2];
            if ((second_allele_char != '.') && ((linebuf_iter[1] == '/') || (linebuf_iter[1] == '|'))) {
              cur_geno = ctow(second_allele_char) - 48;
              if (unlikely(cur_geno > 1)) {
                return kVcfParseInvalidGt;
              }
              if (unlikely(halfcall_mode == kVcfHalfCallError)) {
                return kVcfParseHalfCallError;
              }
              // kVcfHalfCallHaploid, kVcfHalfCallReference
              cur_geno <<= halfcall_mode;
            } else {
              cur_geno = 3;
            }
          } else {
            cur_geno = 3;
          }
        }
        if (!dpr) {
          // now actually handle the unphased dosage
          const uint32_t cur_halfdist = BiallelicDosageHalfdist(dosage_int);
          if (cur_halfdist < dosage_erase_halfdist) {
            // ok for cur_geno to be 'wrong' for now, since it'll get corrected
            // by --hard-call-threshold
            dosage_present_hw |= shifted_bit;
            *dosage_main_iter++ = dosage_int;
          } else {
            // Not saving dosage, since it's too close to an integer, except
            // possibly in the implicit-phased-dosage edge case.
            // If that integer actually conflicts with the hardcall, we must
            // override the hardcall.
            cur_geno = (dosage_int + kDosage4th) / kDosageMid;
            if (phasepresent_hw & shifted_bit) {
              if (cur_geno != 1) {
                // Hardcall-phase no longer applies.
                phasepresent_hw ^= shifted_bit;
              } else if (cur_halfdist * 2 < dphase_erase_halfdist) {
                // Implicit phased-dosage, e.g. 0|0.99.  More stringent
                // dosage_erase_halfdist applies.
                dosage_present_hw |= shifted_bit;
                *dosage_main_iter++ = dosage_int;
              }
            }
          }
        }
      }
    VcfConvertPhasedBiallelicDosageLine_geno_done:
      genovec_word |= cur_geno << (2 * sample_idx_lowbits);
      linebuf_iter = &(cur_gtext_end[1]);
    }
    genovec[widx] = genovec_word;
    phasepresent_alias[widx] = phasepresent_hw;
    phaseinfo_alias[widx] = phaseinfo_hw;
    dosage_present_alias[widx] = dosage_present_hw;
    dphase_present_alias[widx] = dphase_present_hw;
  }
  *dosage_main_iter_ptr = dosage_main_iter;
  *dphase_delta_iter_ptr = dphase_delta_iter;
  return kVcfParseOk;
}

// alt_ct < 10
uintptr_t VcfScanShortallelicLine(const VcfImportBaseContext* vibcp, const char* format_end, char** line_iter_ptr) {
  // Just check for a phased het.
  const VcfHalfCall halfcall_mode = vibcp->halfcall_mode;
  STD_ARRAY_KREF(uint32_t, 2) qual_field_skips = vibcp->qual_field_skips;
  STD_ARRAY_KREF(int32_t, 2) qual_line_mins = vibcp->qual_line_mins;
  STD_ARRAY_KREF(int32_t, 2) qual_line_maxs = vibcp->qual_line_maxs;
  const uint32_t qual_field_ct = vibcp->qual_field_ct;
  const char* phasescan_iter = format_end;
  uintptr_t retval = 0;
  while (1) {
    // this should quickly fail if there are no phased hardcalls at all.
    if (incr_strchrnul_n_mov('|', &phasescan_iter)) {
      goto VcfScanShortallelicLine_no_phased_het;
    }
    if (phasescan_iter[-2] != '\t') {
      // at least one other FORMAT field uses the '|' character.
      // switch to iterating over tabs.
      break;
    }
    if (VcfIsHetShort(&(phasescan_iter[-1]), halfcall_mode)) {
      if ((!qual_field_ct) || (!VcfCheckQuals(qual_field_skips, qual_line_mins, qual_line_maxs, phasescan_iter, FirstPrespace(&(phasescan_iter[2])), qual_field_ct))) {
        goto VcfScanShortallelicLine_phased_het_found;
      }
    }
    if (incr_strchrnul_n_mov('\t', &phasescan_iter)) {
      goto VcfScanShortallelicLine_no_phased_het;
    }
  }
  while (!incr_strchrnul_n_mov('\t', &phasescan_iter)) {
    if ((phasescan_iter[2] != '|') || (!VcfIsHetShort(&(phasescan_iter[1]), halfcall_mode))) {
      continue;
    }
    if ((!qual_field_ct) || (!VcfCheckQuals(qual_field_skips, qual_line_mins, qual_line_maxs, phasescan_iter, FirstPrespace(&(phasescan_iter[4])), qual_field_ct))) {
      goto VcfScanShortallelicLine_phased_het_found;
    }
  }
  while (0) {
  VcfScanShortallelicLine_phased_het_found:
    retval = 1;
  }
 VcfScanShortallelicLine_no_phased_het:
  *line_iter_ptr = K_CAST(char*, phasescan_iter);
  return retval;
}

VcfParseErr VcfConvertUnphasedBiallelicLine(const VcfImportBaseContext* vibcp, const char* linebuf_iter, uintptr_t* genovec) {
  const uint32_t sample_ct = vibcp->sample_ct;
  const uint32_t sample_ctl2_m1 = (sample_ct - 1) / kBitsPerWordD2;
  const VcfHalfCall halfcall_mode = vibcp->halfcall_mode;
  const uint32_t error_on_polyploid = vibcp->error_on_polyploid;
  STD_ARRAY_KREF(uint32_t, 2) qual_field_skips = vibcp->qual_field_skips;
  STD_ARRAY_KREF(int32_t, 2) qual_line_mins = vibcp->qual_line_mins;
  STD_ARRAY_KREF(int32_t, 2) qual_line_maxs = vibcp->qual_line_maxs;
  const uint32_t qual_field_ct = vibcp->qual_field_ct;

  uint32_t inner_loop_last = kBitsPerWordD2 - 1;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= sample_ctl2_m1) {
      if (widx > sample_ctl2_m1) {
        break;
      }
      inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
    }
    uintptr_t genovec_word = 0;
    for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
      const char* cur_gtext_end = FirstPrespace(linebuf_iter);
      if (unlikely((*cur_gtext_end != '\t') && ((sample_idx_lowbits != inner_loop_last) || (widx != sample_ctl2_m1)))) {
        return kVcfParseMissingTokens;
      }
      uintptr_t cur_geno;
      if (qual_field_ct && VcfCheckQuals(qual_field_skips, qual_line_mins, qual_line_maxs, linebuf_iter, cur_gtext_end, qual_field_ct)) {
        // skipping polyploid check for now
        cur_geno = 3;
      } else {
        // Still must check for '|', since phasing_flags bit is unset when all
        // entries are e.g. 0|0.  We just don't bother distinguishing it from
        // '/'.
        const uint32_t is_haploid = (linebuf_iter[1] != '/') && (linebuf_iter[1] != '|');
        cur_geno = ctow(*linebuf_iter) - 48;
        if (cur_geno <= 1) {
          if (is_haploid) {
            cur_geno *= 2;
          } else {
            const char polyploid_char = linebuf_iter[3];
            if ((polyploid_char == '/') || (polyploid_char == '|')) {
              if (unlikely(error_on_polyploid)) {
                return kVcfParsePolyploidError;
              }
              cur_geno = 3;
            } else {
              // bugfix (26 Jul 2018): this needs to be ctow, not ctou32
              const uintptr_t second_allele_idx = ctow(linebuf_iter[2]) - 48;
              if (second_allele_idx <= 1) {
                cur_geno += second_allele_idx;
              } else if (unlikely(second_allele_idx != (~k0LU) * 2)) {
                // not '.'
                return kVcfParseInvalidGt;
              } else if (halfcall_mode == kVcfHalfCallMissing) {
                cur_geno = 3;
              } else if (unlikely(halfcall_mode == kVcfHalfCallError)) {
                return kVcfParseHalfCallError;
              } else {
                // kVcfHalfCallHaploid, kVcfHalfCallReference
                cur_geno <<= halfcall_mode;
              }
            }
          }
        } else {
          if (unlikely(cur_geno != (~k0LU) * 2)) {
            // not '.'
            return kVcfParseInvalidGt;
          }
          cur_geno = 3;
          if (!is_haploid) {
            // ./x
            const char polyploid_char = linebuf_iter[3];
            if ((polyploid_char == '/') || (polyploid_char == '|')) {
              if (unlikely(error_on_polyploid)) {
                return kVcfParsePolyploidError;
              }
              // could perform InvalidGt check?
            } else {
              const char second_allele_char = linebuf_iter[2];
              if ((second_allele_char != '.') && (halfcall_mode != kVcfHalfCallMissing)) {
                cur_geno = ctow(second_allele_char) - 48;
                if (unlikely(cur_geno > 1)) {
                  return kVcfParseInvalidGt;
                }
                if (unlikely(halfcall_mode == kVcfHalfCallError)) {
                  return kVcfParseHalfCallError;
                }
                // kVcfHalfCallHaploid, kVcfHalfCallReference
                cur_geno <<= halfcall_mode;
              }
            }
          }
        }
      }
      genovec_word |= cur_geno << (2 * sample_idx_lowbits);
      linebuf_iter = &(cur_gtext_end[1]);
    }
    genovec[widx] = genovec_word;
  }
  return kVcfParseOk;
}

VcfParseErr VcfConvertUnphasedMultiallelicLine(const VcfImportBaseContext* vibcp, const char* linebuf_iter, uintptr_t allele_ct, uint32_t* __restrict patch_01_ctp, uint32_t* __restrict patch_10_ctp, uintptr_t* __restrict genovec, uintptr_t* __restrict patch_01_set, AlleleCode* __restrict patch_01_vals, uintptr_t* __restrict patch_10_set, AlleleCode* __restrict patch_10_vals) {
  const uint32_t sample_ct = vibcp->sample_ct;
  const uint32_t sample_ctl2_m1 = (sample_ct - 1) / kBitsPerWordD2;
  const VcfHalfCall halfcall_mode = vibcp->halfcall_mode;
  const uint32_t error_on_polyploid = vibcp->error_on_polyploid;
  STD_ARRAY_KREF(uint32_t, 2) qual_field_skips = vibcp->qual_field_skips;
  STD_ARRAY_KREF(int32_t, 2) qual_line_mins = vibcp->qual_line_mins;
  STD_ARRAY_KREF(int32_t, 2) qual_line_maxs = vibcp->qual_line_maxs;
  const uint32_t qual_field_ct = vibcp->qual_field_ct;
  Halfword* patch_01_set_alias = R_CAST(Halfword*, patch_01_set);
  Halfword* patch_10_set_alias = R_CAST(Halfword*, patch_10_set);
  AlleleCode* patch_01_vals_iter = patch_01_vals;
  AlleleCode* patch_10_vals_iter = patch_10_vals;

  uint32_t inner_loop_last = kBitsPerWordD2 - 1;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= sample_ctl2_m1) {
      if (widx > sample_ctl2_m1) {
        if (widx % 2) {
          // might not need this, but lets play it safe
          patch_01_set_alias[widx] = 0;
          patch_10_set_alias[widx] = 0;
        }
        break;
      }
      inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
    }
    uintptr_t genovec_word = 0;
    uint32_t patch_01_hw = 0;
    uint32_t patch_10_hw = 0;
    if (allele_ct <= 10) {
      for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
        const char* cur_gtext_end = FirstPrespace(linebuf_iter);
        if (unlikely((*cur_gtext_end != '\t') && ((sample_idx_lowbits != inner_loop_last) || (widx != sample_ctl2_m1)))) {
          return kVcfParseMissingTokens;
        }
        uintptr_t cur_geno;
        if (qual_field_ct && VcfCheckQuals(qual_field_skips, qual_line_mins, qual_line_maxs, linebuf_iter, cur_gtext_end, qual_field_ct)) {
          // skipping polyploid check for now
          cur_geno = 3;
        } else {
          const uint32_t is_haploid = (linebuf_iter[1] != '/') && (linebuf_iter[1] != '|');
          cur_geno = ctow(*linebuf_iter) - 48;
          if (cur_geno < allele_ct) {
            if (is_haploid) {
              if (cur_geno <= 1) {
                cur_geno *= 2;
              } else {
                const uint32_t shifted_bit = 1U << sample_idx_lowbits;
                patch_10_hw |= shifted_bit;
                *patch_10_vals_iter++ = cur_geno;
                *patch_10_vals_iter++ = cur_geno;
                cur_geno = 2;
              }
            } else {
              const char polyploid_char = linebuf_iter[3];
              if ((polyploid_char == '/') || (polyploid_char == '|')) {
                if (unlikely(error_on_polyploid)) {
                  return kVcfParsePolyploidError;
                }
                cur_geno = 3;
              } else {
                uintptr_t second_allele_idx = ctow(linebuf_iter[2]) - 48;
                if (second_allele_idx < allele_ct) {
                  if (cur_geno <= 1) {
                    if (second_allele_idx <= 1) {
                      cur_geno += second_allele_idx;
                    } else {
                      const uint32_t shifted_bit = 1U << sample_idx_lowbits;
                      if (!cur_geno) {
                        patch_01_hw |= shifted_bit;
                        *patch_01_vals_iter++ = second_allele_idx;
                      } else {
                        patch_10_hw |= shifted_bit;
                        *patch_10_vals_iter++ = 1;
                        *patch_10_vals_iter++ = second_allele_idx;
                      }
                      ++cur_geno;
                    }
                  } else {
                    const uint32_t shifted_bit = 1U << sample_idx_lowbits;
                    if (!second_allele_idx) {
                      patch_01_hw |= shifted_bit;
                      *patch_01_vals_iter++ = cur_geno;
                      cur_geno = 1;
                    } else {
                      patch_10_hw |= shifted_bit;
                      if (second_allele_idx < cur_geno) {
                        // this is not explicitly disallowed, but should almost
                        // never happen
                        const uintptr_t ulii = cur_geno;
                        cur_geno = second_allele_idx;
                        second_allele_idx = ulii;
                      }
                      *patch_10_vals_iter++ = cur_geno;
                      *patch_10_vals_iter++ = second_allele_idx;
                      cur_geno = 2;
                    }
                  }
                } else {
                  // second_allele_idx >= allele_ct
                  if (unlikely(second_allele_idx != (~k0LU) * 2)) {
                    // not '.'
                    return kVcfParseInvalidGt;
                  }
                  if (halfcall_mode == kVcfHalfCallMissing) {
                    cur_geno = 3;
                  } else if (unlikely(halfcall_mode == kVcfHalfCallError)) {
                    return kVcfParseHalfCallError;
                  } else {
                    // kVcfHalfCallHaploid, kVcfHalfCallReference
                    if (cur_geno <= 1) {
                      cur_geno <<= halfcall_mode;
                    } else {
                      const uint32_t shifted_bit = 1U << sample_idx_lowbits;
                      if (halfcall_mode == kVcfHalfCallReference) {
                        patch_01_hw |= shifted_bit;
                        *patch_01_vals_iter++ = cur_geno;
                        cur_geno = 1;
                      } else {
                        patch_10_hw |= shifted_bit;
                        *patch_10_vals_iter++ = cur_geno;
                        *patch_10_vals_iter++ = cur_geno;
                        cur_geno = 2;
                      }
                    }
                  }
                }
              }
            }
          } else if (unlikely(cur_geno != (~k0LU) * 2)) {
            // not '.'
            return kVcfParseInvalidGt;
          } else {
            cur_geno = 3;
            if (!is_haploid) {
              // ./x
              const char polyploid_char = linebuf_iter[3];
              if ((polyploid_char == '/') || (polyploid_char == '|')) {
                if (unlikely(error_on_polyploid)) {
                  return kVcfParsePolyploidError;
                }
                // could perform InvalidGt check?
              } else {
                const char second_allele_char = linebuf_iter[2];
                if ((second_allele_char != '.') && (halfcall_mode != kVcfHalfCallMissing)) {
                  cur_geno = ctow(second_allele_char) - 48;
                  if (unlikely(cur_geno >= allele_ct)) {
                    return kVcfParseInvalidGt;
                  }
                  if (unlikely(halfcall_mode == kVcfHalfCallError)) {
                    return kVcfParseHalfCallError;
                  }
                  // kVcfHalfCallHaploid, kVcfHalfCallReference
                  if (cur_geno <= 1) {
                    cur_geno <<= halfcall_mode;
                  } else {
                    const uint32_t shifted_bit = 1U << sample_idx_lowbits;
                    if (halfcall_mode == kVcfHalfCallReference) {
                      patch_01_hw |= shifted_bit;
                      *patch_01_vals_iter++ = cur_geno;
                      cur_geno = 1;
                    } else {
                      patch_10_hw |= shifted_bit;
                      *patch_10_vals_iter++ = cur_geno;
                      *patch_10_vals_iter++ = cur_geno;
                      cur_geno = 2;
                    }
                  }
                }
              }
            }
          }
        }
        genovec_word |= cur_geno << (2 * sample_idx_lowbits);
        linebuf_iter = &(cur_gtext_end[1]);
      }
    } else {
      // allele_ct > 10
      for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
        const char* cur_gtext_end = FirstPrespace(linebuf_iter);
        if (unlikely((*cur_gtext_end != '\t') && ((sample_idx_lowbits != inner_loop_last) || (widx != sample_ctl2_m1)))) {
          return kVcfParseMissingTokens;
        }
        uintptr_t cur_geno;
        if (qual_field_ct && VcfCheckQuals(qual_field_skips, qual_line_mins, qual_line_maxs, linebuf_iter, cur_gtext_end, qual_field_ct)) {
          // skipping polyploid check for now
          cur_geno = 3;
        } else {
          cur_geno = ctow(*linebuf_iter) - 48;
          if (cur_geno < 10) {
            const char* gt_first_nondigit = &(linebuf_iter[1]);
            uintptr_t next_digit;
            for (; ; ++gt_first_nondigit) {
              next_digit = ctow(*gt_first_nondigit) - 48;
              if (next_digit >= 10) {
                break;
              }
              cur_geno = cur_geno * 10 + next_digit;
              if (cur_geno >= allele_ct) {
                return kVcfParseInvalidGt;
              }
            }
            // '/' = ascii 47, '|' = ascii 124 (i.e. 48 + 76)
            const uint32_t is_haploid = (next_digit != (~k0LU)) && (next_digit != 76);
            if (is_haploid) {
              if (cur_geno <= 1) {
                cur_geno *= 2;
              } else {
                const uint32_t shifted_bit = 1U << sample_idx_lowbits;
                patch_10_hw |= shifted_bit;
                *patch_10_vals_iter++ = cur_geno;
                *patch_10_vals_iter++ = cur_geno;
                cur_geno = 2;
              }
            } else {
              const char* second_allele_iter = &(gt_first_nondigit[1]);
              uintptr_t second_allele_idx = ctow(*second_allele_iter++) - 48;
              if (second_allele_idx < 10) {
                for (; ; ++second_allele_iter) {
                  next_digit = ctow(*second_allele_iter) - 48;
                  if (next_digit >= 10) {
                    break;
                  }
                  second_allele_idx = second_allele_idx * 10 + next_digit;
                  if (second_allele_idx >= allele_ct) {
                    return kVcfParseInvalidGt;
                  }
                }
                if ((next_digit == ~k0LU) || (next_digit == 76)) {
                  if (unlikely(error_on_polyploid)) {
                    return kVcfParsePolyploidError;
                  }
                  cur_geno = 3;
                } else if (cur_geno <= 1) {
                  if (second_allele_idx <= 1) {
                    cur_geno += second_allele_idx;
                  } else {
                    const uint32_t shifted_bit = 1U << sample_idx_lowbits;
                    if (!cur_geno) {
                      patch_01_hw |= shifted_bit;
                      *patch_01_vals_iter++ = second_allele_idx;
                    } else {
                      patch_10_hw |= shifted_bit;
                      *patch_10_vals_iter++ = 1;
                      *patch_10_vals_iter++ = second_allele_idx;
                    }
                    ++cur_geno;
                  }
                } else {
                  // first_allele_idx == cur_geno >= 2
                  const uint32_t shifted_bit = 1U << sample_idx_lowbits;
                  if (!second_allele_idx) {
                    patch_01_hw |= shifted_bit;
                    *patch_01_vals_iter++ = cur_geno;
                    cur_geno = 1;
                  } else {
                    patch_10_hw |= shifted_bit;
                    if (second_allele_idx < cur_geno) {
                      const uintptr_t ulii = cur_geno;
                      cur_geno = second_allele_idx;
                      second_allele_idx = ulii;
                    }
                    *patch_10_vals_iter++ = cur_geno;
                    *patch_10_vals_iter++ = second_allele_idx;
                    cur_geno = 2;
                  }
                }
              } else if (unlikely(second_allele_idx != (~k0LU) * 2)) {
                // not '.'
                return kVcfParseInvalidGt;
              } else {
                const char polyploid_char = *second_allele_iter;
                if ((polyploid_char == '/') || (polyploid_char == '|')) {
                  if (unlikely(error_on_polyploid)) {
                    return kVcfParsePolyploidError;
                  }
                  cur_geno = 3;
                } else if (halfcall_mode == kVcfHalfCallMissing) {
                  cur_geno = 3;
                } else if (unlikely(halfcall_mode == kVcfHalfCallError)) {
                  return kVcfParseHalfCallError;
                } else {
                  // kVcfHalfCallHaploid, kVcfHalfCallReference
                  if (cur_geno <= 1) {
                    cur_geno <<= halfcall_mode;
                  } else {
                    const uint32_t shifted_bit = 1U << sample_idx_lowbits;
                    if (halfcall_mode == kVcfHalfCallReference) {
                      patch_01_hw |= shifted_bit;
                      *patch_01_vals_iter++ = cur_geno;
                      cur_geno = 1;
                    } else {
                      patch_10_hw |= shifted_bit;
                      *patch_10_vals_iter++ = cur_geno;
                      *patch_10_vals_iter++ = cur_geno;
                      cur_geno = 2;
                    }
                  }
                }
              }
            }
          } else if (unlikely(cur_geno != (~k0LU) * 2)) {
            // not '.'
            return kVcfParseInvalidGt;
          } else {
            cur_geno = 3;
            const uint32_t is_haploid = (linebuf_iter[1] != '/') && (linebuf_iter[1] != '|');
            if (!is_haploid) {
              const char second_allele_first_char = linebuf_iter[2];
              if (second_allele_first_char == '.') {
                if (error_on_polyploid) {
                  const char polyploid_char = linebuf_iter[3];
                  if (unlikely((polyploid_char == '/') || (polyploid_char == '|'))) {
                    return kVcfParsePolyploidError;
                  }
                }
              } else {
                // don't check halfcall_mode == kVcfHalfCallMissing yet because
                // we want error_on_polyploid to work
                cur_geno = ctow(second_allele_first_char) - 48;
                if (unlikely(cur_geno >= 10)) {
                  return kVcfParseInvalidGt;
                }
                uintptr_t next_digit;
                for (const char* second_allele_iter = &(linebuf_iter[3]); ; ++second_allele_iter) {
                  next_digit = ctow(*second_allele_iter) - 48;
                  if (next_digit >= 10) {
                    break;
                  }
                  cur_geno = cur_geno * 10 + next_digit;
                  if (unlikely(cur_geno >= allele_ct)) {
                    return kVcfParseInvalidGt;
                  }
                }
                if ((next_digit == ~k0LU) || (next_digit == 76)) {
                  if (unlikely(error_on_polyploid)) {
                    return kVcfParsePolyploidError;
                  }
                  cur_geno = 3;
                } else if (halfcall_mode != kVcfHalfCallMissing) {
                  if (unlikely(halfcall_mode == kVcfHalfCallError)) {
                    return kVcfParseHalfCallError;
                  }
                  // kVcfHalfCallHaploid, kVcfHalfCallReference
                  if (cur_geno <= 1) {
                    cur_geno <<= halfcall_mode;
                  } else {
                    const uint32_t shifted_bit = 1U << sample_idx_lowbits;
                    if (halfcall_mode == kVcfHalfCallReference) {
                      patch_01_hw |= shifted_bit;
                      *patch_01_vals_iter++ = cur_geno;
                      cur_geno = 1;
                    } else {
                      patch_10_hw |= shifted_bit;
                      *patch_10_vals_iter++ = cur_geno;
                      *patch_10_vals_iter++ = cur_geno;
                      cur_geno = 2;
                    }
                  }
                } else {
                  // bugfix (3 Feb 2024)
                  cur_geno = 3;
                }
              }
            }
          }
        }
        genovec_word |= cur_geno << (2 * sample_idx_lowbits);
        linebuf_iter = &(cur_gtext_end[1]);
      }
    }
    genovec[widx] = genovec_word;
    patch_01_set_alias[widx] = patch_01_hw;
    patch_10_set_alias[widx] = patch_10_hw;
  }
  *patch_01_ctp = patch_01_vals_iter - patch_01_vals;
  *patch_10_ctp = S_CAST(uintptr_t, patch_10_vals_iter - patch_10_vals) / 2;
  return kVcfParseOk;
}

VcfParseErr VcfConvertPhasedBiallelicLine(const VcfImportBaseContext* vibcp, const char* linebuf_iter, uintptr_t* genovec, uintptr_t* phasepresent, uintptr_t* phaseinfo) {
  const uint32_t sample_ct = vibcp->sample_ct;
  const uint32_t sample_ctl2_m1 = (sample_ct - 1) / kBitsPerWordD2;
  const VcfHalfCall halfcall_mode = vibcp->halfcall_mode;
  const uint32_t error_on_polyploid = vibcp->error_on_polyploid;
  STD_ARRAY_KREF(uint32_t, 2) qual_field_skips = vibcp->qual_field_skips;
  STD_ARRAY_KREF(int32_t, 2) qual_line_mins = vibcp->qual_line_mins;
  STD_ARRAY_KREF(int32_t, 2) qual_line_maxs = vibcp->qual_line_maxs;
  const uint32_t qual_field_ct = vibcp->qual_field_ct;
  Halfword* phasepresent_alias = R_CAST(Halfword*, phasepresent);
  Halfword* phaseinfo_alias = R_CAST(Halfword*, phaseinfo);

  uint32_t inner_loop_last = kBitsPerWordD2 - 1;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= sample_ctl2_m1) {
      if (widx > sample_ctl2_m1) {
        if (widx % 2) {
          // might not need this, but lets play it safe
          phasepresent_alias[widx] = 0;
        }
        break;
      }
      inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
    }
    uintptr_t genovec_word = 0;
    uint32_t phasepresent_hw = 0;
    uint32_t phaseinfo_hw = 0;
    for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
      const char* cur_gtext_end = FirstPrespace(linebuf_iter);
      if (unlikely((*cur_gtext_end != '\t') && ((sample_idx_lowbits != inner_loop_last) || (widx != sample_ctl2_m1)))) {
        return kVcfParseMissingTokens;
      }
      uintptr_t cur_geno;
      if (qual_field_ct && VcfCheckQuals(qual_field_skips, qual_line_mins, qual_line_maxs, linebuf_iter, cur_gtext_end, qual_field_ct)) {
        // skipping polyploid check for now
        cur_geno = 3;
      } else {
        const uint32_t is_phased = (linebuf_iter[1] == '|');
        const uint32_t is_haploid = (!is_phased) && (linebuf_iter[1] != '/');
        cur_geno = ctow(*linebuf_iter) - 48;
        if (cur_geno <= 1) {
          if (is_haploid) {
            cur_geno *= 2;
          } else {
            const char polyploid_char = linebuf_iter[3];
            if ((polyploid_char == '/') || (polyploid_char == '|')) {
              if (unlikely(error_on_polyploid)) {
                return kVcfParsePolyploidError;
              }
              cur_geno = 3;
            } else {
              const uintptr_t second_allele_idx = ctow(linebuf_iter[2]) - 48;
              if (second_allele_idx <= 1) {
                cur_geno += second_allele_idx;
                if (is_phased && (cur_geno == 1)) {
                  const uint32_t shifted_bit = 1U << sample_idx_lowbits;
                  phasepresent_hw |= shifted_bit;
#ifdef USE_AVX2
                  // compiler seems to automatically handle this well iff it
                  // has access to BMI/BMI2 instructions, at least on my main
                  // dev Mac
                  if (!second_allele_idx) {
                    // 1|0
                    phaseinfo_hw |= shifted_bit;
                  }
#else
                  phaseinfo_hw |= shifted_bit & (second_allele_idx - 1);
#endif
                }
              } else if (unlikely(second_allele_idx != (~k0LU) * 2)) {
                // not '.'
                return kVcfParseInvalidGt;
              } else if (halfcall_mode == kVcfHalfCallMissing) {
                cur_geno = 3;
              } else if (unlikely(halfcall_mode == kVcfHalfCallError)) {
                return kVcfParseHalfCallError;
              } else {
                // kVcfHalfCallHaploid, kVcfHalfCallReference
                if (is_phased && (halfcall_mode == kVcfHalfCallReference) && (cur_geno == 1)) {
                  const uint32_t shifted_bit = 1U << sample_idx_lowbits;
                  phasepresent_hw |= shifted_bit;
                  phaseinfo_hw |= shifted_bit;
                } else {
                  cur_geno <<= halfcall_mode;
                }
              }
            }
          }
        } else if (unlikely(cur_geno != (~k0LU) * 2)) {
          // not '.'
          return kVcfParseInvalidGt;
        } else {
          cur_geno = 3;
          if (!is_haploid) {
            // ./x
            const char second_allele_char = linebuf_iter[2];
            const char polyploid_char = linebuf_iter[3];
            if ((polyploid_char == '/') || (polyploid_char == '|')) {
              if (unlikely(error_on_polyploid)) {
                return kVcfParsePolyploidError;
              }
              // could perform InvalidGt check?
            } else if ((second_allele_char != '.') && (halfcall_mode != kVcfHalfCallMissing)) {
              cur_geno = ctow(second_allele_char) - 48;
              if (unlikely(cur_geno > 1)) {
                return kVcfParseInvalidGt;
              }
              if (unlikely(halfcall_mode == kVcfHalfCallError)) {
                return kVcfParseHalfCallError;
              }
              // kVcfHalfCallHaploid, kVcfHalfCallReference
              if (is_phased && (halfcall_mode == kVcfHalfCallReference) && (cur_geno == 1)) {
                const uint32_t shifted_bit = 1U << sample_idx_lowbits;
                phasepresent_hw |= shifted_bit;
              } else {
                cur_geno <<= halfcall_mode;
              }
            }
          }
        }
      }
      genovec_word |= cur_geno << (2 * sample_idx_lowbits);
      linebuf_iter = &(cur_gtext_end[1]);
    }
    genovec[widx] = genovec_word;
    phasepresent_alias[widx] = phasepresent_hw;
    phaseinfo_alias[widx] = phaseinfo_hw;
  }
  return kVcfParseOk;
}

// todo: try merging this with VcfConvertUnphasedMultiallelicLine, most of the
// code is duplicated
VcfParseErr VcfConvertPhasedMultiallelicLine(const VcfImportBaseContext* vibcp, const char* linebuf_iter, uintptr_t allele_ct, uint32_t* __restrict patch_01_ctp, uint32_t* __restrict patch_10_ctp, uintptr_t* __restrict genovec, uintptr_t* __restrict patch_01_set, AlleleCode* __restrict patch_01_vals, uintptr_t* __restrict patch_10_set, AlleleCode* __restrict patch_10_vals, uintptr_t* __restrict phasepresent, uintptr_t* __restrict phaseinfo) {
  const uint32_t sample_ct = vibcp->sample_ct;
  const uint32_t sample_ctl2_m1 = (sample_ct - 1) / kBitsPerWordD2;
  const VcfHalfCall halfcall_mode = vibcp->halfcall_mode;
  const uint32_t error_on_polyploid = vibcp->error_on_polyploid;
  STD_ARRAY_KREF(uint32_t, 2) qual_field_skips = vibcp->qual_field_skips;
  STD_ARRAY_KREF(int32_t, 2) qual_line_mins = vibcp->qual_line_mins;
  STD_ARRAY_KREF(int32_t, 2) qual_line_maxs = vibcp->qual_line_maxs;
  const uint32_t qual_field_ct = vibcp->qual_field_ct;
  Halfword* patch_01_set_alias = R_CAST(Halfword*, patch_01_set);
  Halfword* patch_10_set_alias = R_CAST(Halfword*, patch_10_set);
  AlleleCode* patch_01_vals_iter = patch_01_vals;
  AlleleCode* patch_10_vals_iter = patch_10_vals;
  Halfword* phasepresent_alias = R_CAST(Halfword*, phasepresent);
  Halfword* phaseinfo_alias = R_CAST(Halfword*, phaseinfo);

  uint32_t inner_loop_last = kBitsPerWordD2 - 1;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= sample_ctl2_m1) {
      if (widx > sample_ctl2_m1) {
        if (widx % 2) {
          // might not need this, but lets play it safe
          patch_01_set_alias[widx] = 0;
          patch_10_set_alias[widx] = 0;
          phasepresent_alias[widx] = 0;
          // phaseinfo doesn't matter
        }
        break;
      }
      inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
    }
    uintptr_t genovec_word = 0;
    uint32_t patch_01_hw = 0;
    uint32_t patch_10_hw = 0;
    uint32_t phasepresent_hw = 0;
    uint32_t phaseinfo_hw = 0;
    uint32_t shifted_bit = 1;
    if (allele_ct <= 10) {
      for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
        const char* cur_gtext_end = FirstPrespace(linebuf_iter);
        if (unlikely((*cur_gtext_end != '\t') && ((sample_idx_lowbits != inner_loop_last) || (widx != sample_ctl2_m1)))) {
          return kVcfParseMissingTokens;
        }
        uintptr_t cur_geno;
        if (qual_field_ct && VcfCheckQuals(qual_field_skips, qual_line_mins, qual_line_maxs, linebuf_iter, cur_gtext_end, qual_field_ct)) {
          // skipping polyploid check for now
          cur_geno = 3;
        } else {
          const uint32_t is_phased = (linebuf_iter[1] == '|');
          const uint32_t is_haploid = (!is_phased) && (linebuf_iter[1] != '/');
          cur_geno = ctow(*linebuf_iter) - 48;
          if (cur_geno < allele_ct) {
            if (is_haploid) {
              if (cur_geno <= 1) {
                cur_geno *= 2;
              } else {
                patch_10_hw |= shifted_bit;
                *patch_10_vals_iter++ = cur_geno;
                *patch_10_vals_iter++ = cur_geno;
                cur_geno = 2;
              }
            } else {
              const char polyploid_char = linebuf_iter[3];
              if ((polyploid_char == '/') || (polyploid_char == '|')) {
                if (unlikely(error_on_polyploid)) {
                  return kVcfParsePolyploidError;
                }
                cur_geno = 3;
              } else {
                uintptr_t second_allele_idx = ctow(linebuf_iter[2]) - 48;
                if (second_allele_idx < allele_ct) {
                  if (cur_geno <= 1) {
                    if (is_phased && (cur_geno != second_allele_idx)) {
                      phasepresent_hw |= shifted_bit;
                      if (!second_allele_idx) {
                        // has to be 1|0 here, since cur_geno <= 1
                        phaseinfo_hw |= shifted_bit;
                      }
                    }
                    if (second_allele_idx <= 1) {
                      cur_geno += second_allele_idx;
                    } else {
                      if (!cur_geno) {
                        patch_01_hw |= shifted_bit;
                        *patch_01_vals_iter++ = second_allele_idx;
                      } else {
                        patch_10_hw |= shifted_bit;
                        *patch_10_vals_iter++ = 1;
                        *patch_10_vals_iter++ = second_allele_idx;
                      }
                      ++cur_geno;
                    }
                  } else {
                    // first_allele_idx == cur_geno >= 2
                    if (is_phased && (cur_geno != second_allele_idx)) {
                      phasepresent_hw |= shifted_bit;
                      if (second_allele_idx < cur_geno) {
                        phaseinfo_hw |= shifted_bit;
                      }
                    }
                    if (!second_allele_idx) {
                      patch_01_hw |= shifted_bit;
                      *patch_01_vals_iter++ = cur_geno;
                      cur_geno = 1;
                    } else {
                      patch_10_hw |= shifted_bit;
                      if (second_allele_idx < cur_geno) {
                        const uintptr_t ulii = cur_geno;
                        cur_geno = second_allele_idx;
                        second_allele_idx = ulii;
                      }
                      *patch_10_vals_iter++ = cur_geno;
                      *patch_10_vals_iter++ = second_allele_idx;
                      cur_geno = 2;
                    }
                  }
                } else {
                  // second_allele_idx >= allele_ct
                  if (unlikely(second_allele_idx != (~k0LU) * 2)) {
                    // not '.'
                    return kVcfParseInvalidGt;
                  }
                  if (halfcall_mode == kVcfHalfCallMissing) {
                    cur_geno = 3;
                  } else if (unlikely(halfcall_mode == kVcfHalfCallError)) {
                    return kVcfParseHalfCallError;
                  } else {
                    // kVcfHalfCallHaploid, kVcfHalfCallReference
                    if (cur_geno <= 1) {
                      if (is_phased && (halfcall_mode == kVcfHalfCallReference) && cur_geno) {
                        phasepresent_hw |= shifted_bit;
                        phaseinfo_hw |= shifted_bit;
                      } else {
                        cur_geno <<= halfcall_mode;
                      }
                    } else if (halfcall_mode == kVcfHalfCallReference) {
                      if (is_phased) {
                        phasepresent_hw |= shifted_bit;
                        phaseinfo_hw |= shifted_bit;
                      }
                      patch_01_hw |= shifted_bit;
                      *patch_01_vals_iter++ = cur_geno;
                      cur_geno = 1;
                    } else {
                      patch_10_hw |= shifted_bit;
                      *patch_10_vals_iter++ = cur_geno;
                      *patch_10_vals_iter++ = cur_geno;
                      cur_geno = 2;
                    }
                  }
                }
              }
            }
          } else if (unlikely(cur_geno != (~k0LU) * 2)) {
            // not '.'
            return kVcfParseInvalidGt;
          } else {
            cur_geno = 3;
            if (!is_haploid) {
              const char polyploid_char = linebuf_iter[3];
              if ((polyploid_char == '/') || (polyploid_char == '|')) {
                if (unlikely(error_on_polyploid)) {
                  return kVcfParsePolyploidError;
                }
                // could perform InvalidGt check?
              } else {
                const char second_allele_char = linebuf_iter[2];
                if ((second_allele_char != '.') && (halfcall_mode != kVcfHalfCallMissing)) {
                  cur_geno = ctow(second_allele_char) - 48;
                  if (unlikely(cur_geno >= allele_ct)) {
                    return kVcfParseInvalidGt;
                  }
                  if (unlikely(halfcall_mode == kVcfHalfCallError)) {
                    return kVcfParseHalfCallError;
                  }
                  if (is_phased && (halfcall_mode == kVcfHalfCallReference) && cur_geno) {
                    phasepresent_hw |= shifted_bit;
                  }
                  // kVcfHalfCallHaploid, kVcfHalfCallReference
                  if (cur_geno <= 1) {
                    cur_geno <<= halfcall_mode;
                  } else {
                    if (halfcall_mode == kVcfHalfCallReference) {
                      patch_01_hw |= shifted_bit;
                      *patch_01_vals_iter++ = cur_geno;
                      cur_geno = 1;
                    } else {
                      patch_10_hw |= shifted_bit;
                      *patch_10_vals_iter++ = cur_geno;
                      *patch_10_vals_iter++ = cur_geno;
                      cur_geno = 2;
                    }
                  }
                }
              }
            }
          }
        }
        genovec_word |= cur_geno << (2 * sample_idx_lowbits);
        shifted_bit *= 2;
        linebuf_iter = &(cur_gtext_end[1]);
      }
    } else {
      // allele_ct > 10
      for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
        const char* cur_gtext_end = FirstPrespace(linebuf_iter);
        if (unlikely((*cur_gtext_end != '\t') && ((sample_idx_lowbits != inner_loop_last) || (widx != sample_ctl2_m1)))) {
          return kVcfParseMissingTokens;
        }
        uintptr_t cur_geno;
        if (qual_field_ct && VcfCheckQuals(qual_field_skips, qual_line_mins, qual_line_maxs, linebuf_iter, cur_gtext_end, qual_field_ct)) {
          // skipping polyploid check for now
          cur_geno = 3;
        } else {
          cur_geno = ctow(*linebuf_iter) - 48;
          if (cur_geno < 10) {
            const char* gt_first_nondigit = &(linebuf_iter[1]);
            uintptr_t next_digit;
            for (; ; ++gt_first_nondigit) {
              next_digit = ctow(*gt_first_nondigit) - 48;
              if (next_digit >= 10) {
                break;
              }
              cur_geno = cur_geno * 10 + next_digit;
              if (cur_geno >= allele_ct) {
                return kVcfParseInvalidGt;
              }
            }
            // '/' = ascii 47, '|' = ascii 124 (i.e. 48 + 76)
            const uint32_t is_phased = (next_digit == 76);
            const uint32_t is_haploid = (!is_phased) && (next_digit != (~k0LU));
            if (is_haploid) {
              if (cur_geno <= 1) {
                cur_geno *= 2;
              } else {
                patch_10_hw |= shifted_bit;
                *patch_10_vals_iter++ = cur_geno;
                *patch_10_vals_iter++ = cur_geno;
                cur_geno = 2;
              }
            } else {
              const char* second_allele_iter = &(gt_first_nondigit[1]);
              uintptr_t second_allele_idx = ctow(*second_allele_iter++) - 48;
              if (second_allele_idx < 10) {
                for (; ; ++second_allele_iter) {
                  next_digit = ctow(*second_allele_iter) - 48;
                  if (next_digit >= 10) {
                    break;
                  }
                  second_allele_idx = second_allele_idx * 10 + next_digit;
                  if (second_allele_idx >= allele_ct) {
                    return kVcfParseInvalidGt;
                  }
                }
                if ((next_digit == ~k0LU) || (next_digit == 76)) {
                  if (unlikely(error_on_polyploid)) {
                    return kVcfParsePolyploidError;
                  }
                  cur_geno = 3;
                } else if (cur_geno <= 1) {
                  if (is_phased && (cur_geno != second_allele_idx)) {
                    phasepresent_hw |= shifted_bit;
                    if (!second_allele_idx) {
                      // has to be 1|0 here, since cur_geno <= 1
                      phaseinfo_hw |= shifted_bit;
                    }
                  }
                  if (second_allele_idx <= 1) {
                    cur_geno += second_allele_idx;
                  } else {
                    if (!cur_geno) {
                      patch_01_hw |= shifted_bit;
                      *patch_01_vals_iter++ = second_allele_idx;
                    } else {
                      patch_10_hw |= shifted_bit;
                      *patch_10_vals_iter++ = 1;
                      *patch_10_vals_iter++ = second_allele_idx;
                    }
                    ++cur_geno;
                  }
                } else {
                  // first_allele_idx == cur_geno >= 2
                  if (is_phased && (cur_geno != second_allele_idx)) {
                    phasepresent_hw |= shifted_bit;
                    if (second_allele_idx < cur_geno) {
                      phaseinfo_hw |= shifted_bit;
                    }
                  }
                  if (!second_allele_idx) {
                    patch_01_hw |= shifted_bit;
                    *patch_01_vals_iter++ = cur_geno;
                    cur_geno = 1;
                  } else {
                    patch_10_hw |= shifted_bit;
                    if (second_allele_idx < cur_geno) {
                      const uintptr_t ulii = cur_geno;
                      cur_geno = second_allele_idx;
                      second_allele_idx = ulii;
                    }
                    *patch_10_vals_iter++ = cur_geno;
                    *patch_10_vals_iter++ = second_allele_idx;
                    cur_geno = 2;
                  }
                }
              } else if (unlikely(second_allele_idx != (~k0LU) * 2)) {
                // not '.'
                return kVcfParseInvalidGt;
              } else {
                const char polyploid_char = *second_allele_iter;
                if ((polyploid_char == '/') || (polyploid_char == '|')) {
                  if (unlikely(error_on_polyploid)) {
                    return kVcfParsePolyploidError;
                  }
                  cur_geno = 3;
                } else if (halfcall_mode == kVcfHalfCallMissing) {
                  cur_geno = 3;
                } else if (unlikely(halfcall_mode == kVcfHalfCallError)) {
                  return kVcfParseHalfCallError;
                } else {
                  // kVcfHalfCallHaploid, kVcfHalfCallReference
                  if (cur_geno <= 1) {
                    if (is_phased && (halfcall_mode == kVcfHalfCallReference) && cur_geno) {
                      phasepresent_hw |= shifted_bit;
                      phaseinfo_hw |= shifted_bit;
                    } else {
                      cur_geno <<= halfcall_mode;
                    }
                  } else if (halfcall_mode == kVcfHalfCallReference) {
                    if (is_phased) {
                      phasepresent_hw |= shifted_bit;
                      phaseinfo_hw |= shifted_bit;
                    }
                    patch_01_hw |= shifted_bit;
                    *patch_01_vals_iter++ = cur_geno;
                    cur_geno = 1;
                  } else {
                    patch_10_hw |= shifted_bit;
                    *patch_10_vals_iter++ = cur_geno;
                    *patch_10_vals_iter++ = cur_geno;
                    cur_geno = 2;
                  }
                }
              }
            }
          } else if (unlikely(cur_geno != (~k0LU) * 2)) {
            // not '.'
            return kVcfParseInvalidGt;
          } else {
            cur_geno = 3;
            const uint32_t is_haploid = (linebuf_iter[1] != '/') && (linebuf_iter[1] != '|');
            if (!is_haploid) {
              const char second_allele_first_char = linebuf_iter[2];
              if (second_allele_first_char == '.') {
                if (error_on_polyploid) {
                  const char polyploid_char = linebuf_iter[3];
                  if (unlikely((polyploid_char == '/') || (polyploid_char == '|'))) {
                    return kVcfParsePolyploidError;
                  }
                }
              } else {
                // don't check halfcall_mode == kVcfHalfCallMissing yet because
                // we want error_on_polyploid to work
                cur_geno = ctow(second_allele_first_char) - 48;
                if (unlikely(cur_geno >= 10)) {
                  return kVcfParseInvalidGt;
                }
                uintptr_t next_digit;
                for (const char* second_allele_iter = &(linebuf_iter[3]); ; ++second_allele_iter) {
                  next_digit = ctow(*second_allele_iter) - 48;
                  if (next_digit >= 10) {
                    break;
                  }
                  cur_geno = cur_geno * 10 + next_digit;
                  if (cur_geno >= allele_ct) {
                    return kVcfParseInvalidGt;
                  }
                }
                if ((next_digit == ~k0LU) || (next_digit == 76)) {
                  if (unlikely(error_on_polyploid)) {
                    return kVcfParsePolyploidError;
                  }
                  cur_geno = 3;
                } else if (halfcall_mode != kVcfHalfCallMissing) {
                  if (unlikely(halfcall_mode == kVcfHalfCallError)) {
                    return kVcfParseHalfCallError;
                  }
                  if ((linebuf_iter[1] == '|') && (halfcall_mode == kVcfHalfCallReference) && cur_geno) {
                    phasepresent_hw |= shifted_bit;
                  }
                  // kVcfHalfCallHaploid, kVcfHalfCallReference
                  if (cur_geno <= 1) {
                    cur_geno <<= halfcall_mode;
                  } else {
                    if (halfcall_mode == kVcfHalfCallReference) {
                      patch_01_hw |= shifted_bit;
                      *patch_01_vals_iter++ = cur_geno;
                      cur_geno = 1;
                    } else {
                      patch_10_hw |= shifted_bit;
                      *patch_10_vals_iter++ = cur_geno;
                      *patch_10_vals_iter++ = cur_geno;
                      cur_geno = 2;
                    }
                  }
                } else {
                  // bugfix (3 Feb 2024)
                  cur_geno = 3;
                }
              }
            }
          }
        }
        genovec_word |= cur_geno << (2 * sample_idx_lowbits);
        shifted_bit *= 2;
        linebuf_iter = &(cur_gtext_end[1]);
      }
    }
    genovec[widx] = genovec_word;
    patch_01_set_alias[widx] = patch_01_hw;
    patch_10_set_alias[widx] = patch_10_hw;
    phasepresent_alias[widx] = phasepresent_hw;
    phaseinfo_alias[widx] = phaseinfo_hw;
  }
  *patch_01_ctp = patch_01_vals_iter - patch_01_vals;
  *patch_10_ctp = S_CAST(uintptr_t, patch_10_vals_iter - patch_10_vals) / 2;
  return kVcfParseOk;
}

uintptr_t VcfScanLongallelicLine(const VcfImportBaseContext* vibcp, const char* format_end, char** line_iter_ptr) {
  const VcfHalfCall halfcall_mode = vibcp->halfcall_mode;
  STD_ARRAY_KREF(uint32_t, 2) qual_field_skips = vibcp->qual_field_skips;
  STD_ARRAY_KREF(int32_t, 2) qual_line_mins = vibcp->qual_line_mins;
  STD_ARRAY_KREF(int32_t, 2) qual_line_maxs = vibcp->qual_line_maxs;
  const uint32_t qual_field_ct = vibcp->qual_field_ct;
  const char* phasescan_iter = format_end;
  uintptr_t retval = 0;
  while (1) {
    const char* next_vbar = phasescan_iter;
    // this should quickly fail if there are no phased hardcalls at all.
    if (strchrnul_n_mov('|', &next_vbar)) {
      goto VcfScanLongallelicLine_no_phased_het;
    }
    phasescan_iter = &(next_vbar[-1]);
    unsigned char ucc = *phasescan_iter;
    if (ucc == '.') {
      --phasescan_iter;
    } else {
      while (IsDigit(ucc)) {
        ucc = *(--phasescan_iter);
      }
    }
    if (*phasescan_iter != '\t') {
      // at least one other FORMAT field uses the '|' character.
      // switch to iterating over tabs.
      break;
    }
    if (VcfIsHetLong(&(phasescan_iter[1]), &(next_vbar[1]), halfcall_mode)) {
      if ((!qual_field_ct) || (!VcfCheckQuals(qual_field_skips, qual_line_mins, qual_line_maxs, phasescan_iter, FirstPrespace(&(next_vbar[2])), qual_field_ct))) {
        goto VcfScanLongallelicLine_phased_het_found;
      }
    }
    if (incr_strchrnul_n_mov('\t', &phasescan_iter)) {
      goto VcfScanLongallelicLine_no_phased_het;
    }
  }
  while (!incr_strchrnul_n_mov('\t', &phasescan_iter)) {
    const char* next_vbar = &(phasescan_iter[1]);
    while (IsDigit(*next_vbar)) {
      ++next_vbar;
    }
    if ((*next_vbar != '|') || (!VcfIsHetLong(&(phasescan_iter[1]), &(next_vbar[1]), halfcall_mode))) {
      continue;
    }
    if ((!qual_field_ct) || (!VcfCheckQuals(qual_field_skips, qual_line_mins, qual_line_maxs, phasescan_iter, FirstPrespace(&(next_vbar[2])), qual_field_ct))) {
      goto VcfScanLongallelicLine_phased_het_found;
    }
  }
  while (0) {
  VcfScanLongallelicLine_phased_het_found:
    retval = 1;
  }
 VcfScanLongallelicLine_no_phased_het:
  *line_iter_ptr = K_CAST(char*, phasescan_iter);
  return retval;
}

typedef struct VcfGenoToPgenCtxStruct {
  uint32_t sample_ct;
  uint32_t hard_call_halfdist;
  uint32_t dosage_erase_halfdist;
  double import_dosage_certainty;
  VcfHalfCall halfcall_mode;
  uint32_t error_on_polyploid;
  uint32_t dosage_is_gp;
  STD_ARRAY_DECL(int32_t, 2, qual_mins);
  STD_ARRAY_DECL(int32_t, 2, qual_maxs);

  unsigned char** thread_wkspaces;

  uint32_t* thread_bidxs[2];
  GparseRecord* gparse[2];
  const uintptr_t* block_allele_idx_offsets[2];

  // PglErr set by main thread
  VcfParseErr* vcf_parse_errs;
  uintptr_t* err_line_idxs;
  uint32_t parse_failed;
} VcfGenoToPgenCtx;

THREAD_FUNC_DECL VcfGenoToPgenThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  VcfGenoToPgenCtx* ctx = S_CAST(VcfGenoToPgenCtx*, arg->sharedp->context);

  VcfImportContext vic;
  const uint32_t sample_ct = ctx->sample_ct;
  const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  vic.vibc.sample_ct = sample_ct;
  vic.vibc.halfcall_mode = ctx->halfcall_mode;
  vic.vibc.error_on_polyploid = ctx->error_on_polyploid;
  vic.dosage_is_gp = ctx->dosage_is_gp;
  vic.dosage_erase_halfdist = ctx->dosage_erase_halfdist;
  vic.import_dosage_certainty = ctx->import_dosage_certainty;
  const uint32_t hard_call_halfdist = ctx->hard_call_halfdist;
  STD_ARRAY_DECL(int32_t, 2, qual_mins);
  STD_ARRAY_DECL(int32_t, 2, qual_maxs);
  STD_ARRAY_COPY(ctx->qual_mins, 2, qual_mins);
  STD_ARRAY_COPY(ctx->qual_maxs, 2, qual_maxs);
  unsigned char* thread_wkspace = ctx->thread_wkspaces[tidx];
  uintptr_t* patch_01_set = nullptr;
  AlleleCode* patch_01_vals = nullptr;
  uintptr_t* patch_10_set = nullptr;
  AlleleCode* patch_10_vals = nullptr;
  uintptr_t* phasepresent = nullptr;
  uintptr_t* phaseinfo = nullptr;
  uintptr_t* dosage_present = nullptr;
  Dosage* dosage_main = nullptr;
  uintptr_t* dphase_present = nullptr;
  SDosage* dphase_delta = nullptr;
  SDosage* tmp_dphase_delta = R_CAST(SDosage*, thread_wkspace);
  thread_wkspace = &(thread_wkspace[RoundUpPow2(sample_ct * sizeof(SDosage), kBytesPerVec)]);
  uintptr_t* write_patch_01_set = nullptr;
  AlleleCode* write_patch_01_vals = nullptr;
  uintptr_t* write_patch_10_set = nullptr;
  AlleleCode* write_patch_10_vals = nullptr;
  uintptr_t* write_phasepresent = nullptr;
  uintptr_t* write_phaseinfo = nullptr;
  uintptr_t* write_dosage_present = nullptr;
  Dosage* write_dosage_main = nullptr;
  uintptr_t* write_dphase_present = nullptr;
  SDosage* write_dphase_delta = nullptr;
  uint32_t cur_allele_ct = 2;
  uint32_t parity = 0;
  VcfParseErr vcf_parse_err = kVcfParseOk;
  uintptr_t line_idx = 0;
  do {
    const uintptr_t* block_allele_idx_offsets = ctx->block_allele_idx_offsets[parity];
    const uint32_t bidx_end = ctx->thread_bidxs[parity][tidx + 1];
    GparseRecord* cur_gparse = ctx->gparse[parity];

    for (uint32_t bidx = ctx->thread_bidxs[parity][tidx]; bidx != bidx_end; ++bidx) {
      GparseRecord* grp = &(cur_gparse[bidx]);
      uint32_t patch_01_ct = 0;
      uint32_t patch_10_ct = 0;
      uint32_t cur_phasepresent_exists = 0;
      uint32_t dosage_ct = 0;
      uint32_t dphase_ct = 0;
      GparseFlags gparse_flags = grp->flags;
      if (gparse_flags == kfGparseNull) {
        SetAllBits(2 * sample_ct, R_CAST(uintptr_t*, grp->record_start));
      } else {
        if (block_allele_idx_offsets) {
          cur_allele_ct = block_allele_idx_offsets[bidx + 1] - block_allele_idx_offsets[bidx];
        }
        uintptr_t* genovec = GparseGetPointers(thread_wkspace, sample_ct, cur_allele_ct, gparse_flags, &patch_01_set, &patch_01_vals, &patch_10_set, &patch_10_vals, &phasepresent, &phaseinfo, &dosage_present, &dosage_main, &dphase_present, &dphase_delta);
        uintptr_t* write_genovec = GparseGetPointers(grp->record_start, sample_ct, cur_allele_ct, gparse_flags, &write_patch_01_set, &write_patch_01_vals, &write_patch_10_set, &write_patch_10_vals, &write_phasepresent, &write_phaseinfo, &write_dosage_present, &write_dosage_main, &write_dphase_present, &write_dphase_delta);
        char* genotext_start = R_CAST(char*, grp->record_start);
        ++genotext_start;
        uint32_t qual_field_ct = grp->metadata.read_vcf.qual_exists;
        if (qual_field_ct) {
          qual_field_ct = VcfQualScanInit2(grp->metadata.read_vcf.qual_field_idxs, qual_mins, qual_maxs, vic.vibc.qual_field_skips, vic.vibc.qual_line_mins, vic.vibc.qual_line_maxs);
        }
        vic.vibc.gt_exists = grp->metadata.read_vcf.gt_exists;
        vic.vibc.qual_field_ct = qual_field_ct;
        vic.dosage_field_idx = grp->metadata.read_vcf.dosage_field_idx;
        vic.hds_field_idx = grp->metadata.read_vcf.hds_field_idx;
        if ((vic.hds_field_idx == UINT32_MAX) && (vic.dosage_field_idx == UINT32_MAX)) {
          if (!(gparse_flags & kfGparseHphase)) {
            if (cur_allele_ct == 2) {
              vcf_parse_err = VcfConvertUnphasedBiallelicLine(&(vic.vibc), genotext_start, genovec);
            } else {
              vcf_parse_err = VcfConvertUnphasedMultiallelicLine(&(vic.vibc), genotext_start, cur_allele_ct, &patch_01_ct, &patch_10_ct, genovec, patch_01_set, patch_01_vals, patch_10_set, patch_10_vals);
            }
          } else {
            if (cur_allele_ct == 2) {
              vcf_parse_err = VcfConvertPhasedBiallelicLine(&(vic.vibc), genotext_start, genovec, phasepresent, phaseinfo);
            } else {
              vcf_parse_err = VcfConvertPhasedMultiallelicLine(&(vic.vibc), genotext_start, cur_allele_ct, &patch_01_ct, &patch_10_ct, genovec, patch_01_set, patch_01_vals, patch_10_set, patch_10_vals, phasepresent, phaseinfo);
            }
            cur_phasepresent_exists = !AllWordsAreZero(phasepresent, sample_ctl);
          }
          if (unlikely(vcf_parse_err)) {
            line_idx = grp->metadata.read_vcf.line_idx;
            goto VcfGenoToPgenThread_malformed;
          }
        } else {
          Dosage* dosage_main_iter = dosage_main;
          SDosage* dphase_delta_iter = dphase_delta;
          if (cur_allele_ct == 2) {
            vcf_parse_err = VcfConvertPhasedBiallelicDosageLine(&vic, genotext_start, genovec, phasepresent, phaseinfo, dosage_present, dphase_present, &dosage_main_iter, &dphase_delta_iter);
          } else {
            // multiallelic dosage: shouldn't be possible to get here yet
            exit(S_CAST(int32_t, kPglRetInternalError));
          }
          if (unlikely(vcf_parse_err)) {
            line_idx = grp->metadata.read_vcf.line_idx;
            goto VcfGenoToPgenThread_malformed;
          }
          dosage_ct = dosage_main_iter - dosage_main;
          if (dosage_ct) {
            dphase_ct = ApplyHardCallThreshPhased(dosage_present, dosage_main, dosage_ct, hard_call_halfdist, genovec, phasepresent, phaseinfo, dphase_present, dphase_delta, tmp_dphase_delta);
            memcpy(write_dosage_present, dosage_present, sample_ctl * sizeof(intptr_t));
            memcpy(write_dosage_main, dosage_main, dosage_ct * sizeof(Dosage));
            if (dphase_ct) {
              memcpy(write_dphase_present, dphase_present, sample_ctl * sizeof(intptr_t));
              memcpy(write_dphase_delta, dphase_delta, dphase_ct * sizeof(SDosage));
            }
          }
          cur_phasepresent_exists = !AllWordsAreZero(phasepresent, sample_ctl);
        }
        memcpy(write_genovec, genovec, sample_ctl2 * sizeof(intptr_t));
        if (patch_01_ct) {
          memcpy(write_patch_01_set, patch_01_set, sample_ctl * sizeof(intptr_t));
          memcpy(write_patch_01_vals, patch_01_vals, patch_01_ct * sizeof(AlleleCode));
        }
        if (patch_10_ct) {
          memcpy(write_patch_10_set, patch_10_set, sample_ctl * sizeof(intptr_t));
          memcpy(write_patch_10_vals, patch_10_vals, patch_10_ct * sizeof(AlleleCode) * 2);
        }
        if (cur_phasepresent_exists || dphase_ct) {
          memcpy(write_phasepresent, phasepresent, sample_ctl * sizeof(intptr_t));
          memcpy(write_phaseinfo, phaseinfo, sample_ctl * sizeof(intptr_t));
        }
      }

      grp->metadata.write.patch_01_ct = patch_01_ct;
      grp->metadata.write.patch_10_ct = patch_10_ct;
      grp->metadata.write.phasepresent_exists = cur_phasepresent_exists;
      grp->metadata.write.dosage_ct = dosage_ct;
      grp->metadata.write.multiallelic_dosage_ct = 0;
      grp->metadata.write.dphase_ct = dphase_ct;
      grp->metadata.write.multiallelic_dphase_ct = 0;
    }
    while (0) {
    VcfGenoToPgenThread_malformed:
      ctx->vcf_parse_errs[tidx] = vcf_parse_err;
      ctx->err_line_idxs[tidx] = line_idx;
      ctx->parse_failed = 1;
      break;
    }
    parity = 1 - parity;
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

static const char kGpText[] = "GP";

// pgen_generated and psam_generated assumed to be initialized to 1.
static_assert(!kVcfHalfCallReference, "VcfToPgen() assumes kVcfHalfCallReference == 0.");
static_assert(kVcfHalfCallHaploid == 1, "VcfToPgen() assumes kVcfHalfCallHaploid == 1.");
PglErr VcfToPgen(const char* vcfname, const char* preexisting_psamname, const char* const_fid, const char* dosage_import_field, MiscFlags misc_flags, ImportFlags import_flags, LoadFilterLogFlags load_filter_log_import_flags, uint32_t no_samples_ok, uint32_t is_update_or_impute_sex, uint32_t is_splitpar, uint32_t is_sortvars, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, double import_dosage_certainty, char id_delim, char idspace_to, int32_t vcf_min_gq, int32_t vcf_min_dp, int32_t vcf_max_dp, VcfHalfCall halfcall_mode, FamCol fam_cols, uint32_t import_max_allele_ct, uint32_t max_thread_ct, char* outname, char* outname_end, ChrInfo* cip, uint32_t* pgen_generated_ptr, uint32_t* psam_generated_ptr) {
  // Performs a 2-pass load.  Probably staying that way after sequential writer
  // is implemented since header lines are a pain.
  //
  // preexisting_psamname should be nullptr if no such file was specified.
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  char* pvar_cswritep = nullptr;
  uintptr_t line_idx = 1;
  const uint32_t half_call_explicit_error = (halfcall_mode == kVcfHalfCallError);
  // == 2 when searched for and found in header
  // then becomes a boolean
  uint32_t format_hds_search = 0;
  PglErr reterr = kPglRetSuccess;
  VcfParseErr vcf_parse_err = kVcfParseOk;
  CompressStreamState pvar_css;
  PreinitCstream(&pvar_css);
  ThreadGroup tg;
  PreinitThreads(&tg);
  VcfGenoToPgenCtx ctx;
  TextStream vcf_txs;
  STPgenWriter spgw;
  PreinitTextStream(&vcf_txs);
  PreinitSpgw(&spgw);
  {
    uint32_t max_line_blen;
    if (StandardizeMaxLineBlen(bigstack_left() / 4, &max_line_blen)) {
      goto VcfToPgen_ret_NOMEM;
    }
    // probable todo: if chromosome filter specified, take advantage of an
    // index file with the standard extension if it's present.  (The index
    // reader can be very simple since *only* chromosome filters are supported
    // by import functions.)  Do the same for .bgen and .bcf files.

    reterr = ForceNonFifo(vcfname);
    if (unlikely(reterr)) {
      if (reterr == kPglRetOpenFail) {
        const uint32_t slen = strlen(vcfname);
        if ((!StrEndsWith(vcfname, ".vcf", slen)) &&
            (!StrEndsWith(vcfname, ".vcf.gz", slen))) {
          logerrprintfww("Error: Failed to open %s : %s. (--vcf expects a complete filename; did you forget '.vcf' at the end?)\n", vcfname, strerror(errno));
        } else {
          logerrprintfww(kErrprintfFopen, vcfname, strerror(errno));
        }
      } else {
        logerrprintfww(kErrprintfRewind, vcfname);
      }
      goto VcfToPgen_ret_1;
    }
    reterr = InitTextStreamEx(vcfname, 1, kMaxLongLine, max_line_blen, ClipU32(max_thread_ct - 1, 1, 4), &vcf_txs);
    if (unlikely(reterr)) {
      goto VcfToPgen_ret_TSTREAM_FAIL;
    }
    const uint32_t prohibit_extra_chr = (misc_flags / kfMiscProhibitExtraChr) & 1;
    uint32_t dosage_import_field_slen = 0;

    uint32_t unforced_gp = 0;
    if (dosage_import_field) {
      dosage_import_field_slen = strlen(dosage_import_field);
      ctx.dosage_is_gp = 0;
      if (strequal_k(dosage_import_field, "HDS", dosage_import_field_slen)) {
        format_hds_search = 1;
        // special case: search for DS and HDS
        dosage_import_field = &(dosage_import_field[1]);
        dosage_import_field_slen = 2;
      } else if (strequal_k(dosage_import_field, "GP-force", dosage_import_field_slen)) {
        dosage_import_field = kGpText;
        dosage_import_field_slen = 2;
        ctx.dosage_is_gp = 1;
      } else if (strequal_k(dosage_import_field, "GP", dosage_import_field_slen)) {
        unforced_gp = (import_dosage_certainty == 0.0);
        // bugfix (20 Aug 2018): forgot to initialize dosage_is_gp
        ctx.dosage_is_gp = 1;
      }
    }

    ctx.qual_mins[0] = vcf_min_gq;
    ctx.qual_mins[1] = vcf_min_dp;
    ctx.qual_maxs[0] = 0x7fffffff;
    ctx.qual_maxs[1] = vcf_max_dp;
    VcfImportContext vic;
    // sample_ct set later
    if (halfcall_mode == kVcfHalfCallDefault) {
      halfcall_mode = kVcfHalfCallError;
    }
    ctx.halfcall_mode = halfcall_mode;  // bugfix (28 Nov 2019)
    ctx.error_on_polyploid = !(import_flags & kfImportPolyploidMissing);
    vic.vibc.halfcall_mode = halfcall_mode;
    vic.vibc.error_on_polyploid = ctx.error_on_polyploid;

    vic.dosage_is_gp = ctx.dosage_is_gp;

    // always positive
    vic.dosage_erase_halfdist = kDosage4th - dosage_erase_thresh;

    vic.import_dosage_certainty = import_dosage_certainty;

    vic.dosage_field_idx = UINT32_MAX;
    vic.hds_field_idx = UINT32_MAX;

    max_line_blen = 1;  // Now means "max_observed_line_blen".
    uint32_t format_gt_exists = 0;
    uint32_t format_gq_relevant = 0;
    uint32_t format_dp_relevant = 0;
    uint32_t format_dosage_relevant = 0;
    uint32_t info_pr_exists = 0;
    uint32_t info_pr_nonflag_exists = 0;
    uint32_t info_nonpr_exists = 0;
    uint32_t chrset_exists = 0;
    char* line_iter = TextLineEnd(&vcf_txs);
    char* prev_line_start;
    for (; ; ++line_idx) {
      // don't tolerate leading spaces
      reterr = TextNextLineUnsafe(&vcf_txs, &line_iter);
      if (unlikely(reterr)) {
        if (reterr == kPglRetEof) {
          logerrputs("Error: No #CHROM header line or variant records in --vcf file.\n");
          goto VcfToPgen_ret_MALFORMED_INPUT;
        }
        goto VcfToPgen_ret_TSTREAM_FAIL;
      }
      prev_line_start = line_iter;
      if (unlikely(*line_iter != '#')) {
        if ((line_idx == 1) && memequal_sk(line_iter, "BCF")) {
          // this is more informative than "missing header line"...
          if (line_iter[3] == 2) {
            snprintf(g_logbuf, kLogbufSize, "Error: %s appears to be a BCF2 file. Try --bcf instead of --vcf.\n", vcfname);
            goto VcfToPgen_ret_MALFORMED_INPUT_WW;
          }
          if (line_iter[3] == 4) {
            snprintf(g_logbuf, kLogbufSize, "Error: %s appears to be a BCF1 file. Use 'bcftools view' to convert it to a PLINK-readable VCF.\n", vcfname);
            goto VcfToPgen_ret_MALFORMED_INPUT_WW;
          }
        }
        logerrputs("Error: No #CHROM header line in --vcf file.\n");
        goto VcfToPgen_ret_MALFORMED_INPUT;
      }
      if (line_iter[1] != '#') {
        break;
      }
      // Recognized header lines:
      // ##fileformat: discard (regenerate; todo: conditionally error out)
      // ##fileDate: discard (regenerate)
      // ##source: discard (regenerate)
      // ##contig: conditionally keep
      // ##INFO: note presence of INFO/PR, note presence of at least one non-PR
      //         field, keep data (though, if INFO/PR is the *only* field,
      //         omit it from the .pvar for consistency with --make-pgen
      //         default)
      //         update (8 Sep 2017): nonflag INFO/PR is noted, and not treated
      //         specially unless provisional-reference INFO/PR output would
      //         conflict with it
      // ##FORMAT: note presence of FORMAT/GT and FORMAT/GP, discard
      //           (regenerate)
      // ##chrSet: if recognized, perform consistency check and/or update
      //           chr_info
      //
      // Everything else (##FILTER, ##reference, etc.) is passed through
      // unchanged.  FILTER values in the VCF body do not have to be mentioned
      // in the header (since only BCF, not VCF, spec requires that).
      //
      // Because of how ##contig is handled (we only keep the lines which
      // correspond to chromosomes/contigs actually present in the VCF, and not
      // filtered out), we wait until the second pass to write the .pvar.
      if (StrStartsWithUnsafe(&(line_iter[2]), "chrSet=<")) {
        if (unlikely(chrset_exists)) {
          logerrputs("Error: Multiple ##chrSet header lines in --vcf file.\n");
          goto VcfToPgen_ret_MALFORMED_INPUT;
        }
        chrset_exists = 1;
        // .pvar loader will print a warning if necessary
        reterr = ReadChrsetHeaderLine(&(line_iter[2 + strlen("chrSet=<")]), "--vcf file", misc_flags, line_idx, cip);
        if (unlikely(reterr)) {
          goto VcfToPgen_ret_1;
        }
      } else if (StrStartsWithUnsafe(&(line_iter[2]), "FORMAT=<")) {
        line_iter = &(line_iter[2 + strlen("FORMAT=<")]);
        char* idval;
        uint32_t id_slen;
        if (unlikely(HkvlineId(&line_iter, &idval, &id_slen))) {
          goto VcfToPgen_ret_MALFORMED_HEADER_LINE;
        }
        char* numstr;
        char* typestr;
        uint32_t num_slen;
        uint32_t type_slen;
        if (unlikely(HkvlineNumType(line_iter, &numstr, &num_slen, &typestr, &type_slen))) {
          goto VcfToPgen_ret_MALFORMED_HEADER_LINE;
        }
        if (strequal_k(idval, "GT", id_slen)) {
          if (unlikely(format_gt_exists)) {
            logerrputs("Error: Duplicate FORMAT/GT header line in --vcf file.\n");
            goto VcfToPgen_ret_MALFORMED_INPUT;
          }
          if (unlikely((!strequal_k(numstr, "1", num_slen)) || (!strequal_k(typestr, "String", type_slen)))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Header line %" PRIuPTR " of --vcf file does not have expected FORMAT/GT format.\n", line_idx);
            goto VcfToPgen_ret_MALFORMED_INPUT_WW;
          }
          // bugfix (21 Nov 2023)
          format_gt_exists = 1;
        } else if ((vcf_min_gq != -1) && strequal_k(idval, "GQ", id_slen)) {
          if (unlikely(format_gq_relevant)) {
            logerrputs("Error: Duplicate FORMAT/GQ header line in --vcf file.\n");
            goto VcfToPgen_ret_MALFORMED_INPUT;
          }
          if (unlikely(!strequal_k(numstr, "1", num_slen))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Header line %" PRIuPTR " of --vcf file does not have expected FORMAT/GQ format.\n", line_idx);
            goto VcfToPgen_ret_MALFORMED_INPUT_WW;
          }
          format_gq_relevant = 1;
        } else if (((vcf_min_dp != -1) || (vcf_max_dp != 0x7fffffff)) && strequal_k(idval, "DP", id_slen)) {
          if (unlikely(format_dp_relevant)) {
            logerrputs("Error: Duplicate FORMAT/DP header line in --vcf file.\n");
            goto VcfToPgen_ret_MALFORMED_INPUT;
          }
          if (unlikely(!strequal_k(numstr, "1", num_slen))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Header line %" PRIuPTR " of --vcf file does not have expected FORMAT/DP format.\n", line_idx);
          }
          format_dp_relevant = 1;
        } else if (dosage_import_field) {
          if ((id_slen == dosage_import_field_slen) && memequal(idval, dosage_import_field, id_slen)) {
            if (unlikely(format_dosage_relevant)) {
              snprintf(g_logbuf, kLogbufSize, "Error: Duplicate FORMAT/%s header line in --vcf file.\n", dosage_import_field);
              goto VcfToPgen_ret_MALFORMED_INPUT_WW;
            }
            format_dosage_relevant = 1;
          } else if (format_hds_search && strequal_k(idval, "HDS", id_slen)) {
            if (unlikely(format_hds_search == 2)) {
              logerrputs("Error: Duplicate FORMAT/HDS header line in --vcf file.\n");
              goto VcfToPgen_ret_MALFORMED_INPUT;
            }
            format_hds_search = 2;
          } else if (unforced_gp && strequal_k(idval, "DS", id_slen)) {
            logerrputs("Error: --vcf dosage=GP specified, but --import-dosage-certainty was not and\nFORMAT/DS header line is present.\nSince " PROG_NAME_STR " collapses genotype probabilities down to dosages (even when\nperforming simple operations like \"" PROG_NAME_STR " --vcf ... --export bgen-1.2 ...\"),\n'dosage=GP' almost never makes sense in this situation.  Either change it to\n'dosage=DS' (if dosages are good enough for your analysis), or use another\nprogram to work with the genotype probabilities.\nThere is one notable exception to the preceding recommendation: you are writing\na script to process VCF files that are guaranteed to have FORMAT/GP, but may or\nmay not have FORMAT/DS.  You can use 'dosage=GP-force' to suppress this error\nin that situation.\n");
            goto VcfToPgen_ret_INCONSISTENT_INPUT;
          }
        }
      } else if (StrStartsWithUnsafe(&(line_iter[2]), "INFO=<")) {
        line_iter = &(line_iter[2 + strlen("INFO=<")]);
        char* idval;
        uint32_t id_slen;
        if (unlikely(HkvlineId(&line_iter, &idval, &id_slen))) {
          goto VcfToPgen_ret_MALFORMED_HEADER_LINE;
        }
        if (strequal_k(idval, "PR", id_slen)) {
          if (unlikely(info_pr_exists || info_pr_nonflag_exists)) {
            logerrputs("Error: Duplicate INFO/PR header line in --vcf file.\n");
            goto VcfToPgen_ret_MALFORMED_INPUT;
          }
          char* typestr;
          uint32_t type_slen;
          if (unlikely(HkvlineFind(line_iter, "Type", &typestr, &type_slen))) {
            goto VcfToPgen_ret_MALFORMED_HEADER_LINE;
          }
          info_pr_nonflag_exists = !strequal_k(typestr, "Flag", type_slen);
          info_pr_exists = 1 - info_pr_nonflag_exists;
          if (info_pr_nonflag_exists) {
            logerrprintfww("Warning: Header line %" PRIuPTR " of --vcf file has an unexpected definition of INFO/PR. This interferes with a few merge and liftover operations.\n", line_idx);
          }
        } else {
          info_nonpr_exists = 1;
        }
      }
      line_iter = AdvPastDelim(line_iter, '\n');
      const uint32_t prev_line_blen = line_iter - prev_line_start;
      if (prev_line_blen > max_line_blen) {
        max_line_blen = prev_line_blen;
      }
    }
    const uint32_t ref_n_missing = (import_flags / kfImportVcfRefNMissing) & 1;
    if (unlikely(ref_n_missing && (!info_pr_exists))) {
      logerrputs("Error: --vcf-ref-n-missing was specified, but the VCF does not have the\nINFO/PR header line that should be present in any .ped-derived VCF.\n");
      goto VcfToPgen_ret_INCONSISTENT_INPUT;
    }
    const uint32_t require_gt = (load_filter_log_import_flags / kfLoadFilterLogVcfRequireGt) & 1;
    if (unlikely((!format_gt_exists) && require_gt)) {
      logerrputs("Error: No FORMAT/GT key found in --vcf file header, when --vcf-require-gt was\nspecified.\n(If this header line is actually present, but with extra spaces or unusual\nfield ordering, standardize the header with e.g. bcftools.)\n");
      goto VcfToPgen_ret_INCONSISTENT_INPUT;
    }
    if ((!format_gq_relevant) && (vcf_min_gq != -1)) {
      logerrputs("Warning: No FORMAT/GQ key found in --vcf file header.  --vcf-min-gq ignored.\n(If this header line is actually present, but with extra spaces or unusual\nfield ordering, standardize the header with e.g. bcftools.)\n");
      vcf_min_gq = -1;
    }
    if ((!format_dp_relevant) && ((vcf_max_dp != 0x7fffffff) || (vcf_min_dp != -1))) {
      logerrputs("Warning: No FORMAT/DP key found in --vcf file header.  --vcf-{max,min}-dp\nignored.\n(If this header line is actually present, but with extra spaces or unusual\nfield ordering, standardize the header with e.g. bcftools.)\n");
      vcf_max_dp = 0x7fffffff;
      vcf_min_dp = -1;
    }
    const uint32_t format_gq_or_dp_relevant = format_gq_relevant || format_dp_relevant;
    if (format_hds_search) {
      --format_hds_search;
      if (!format_hds_search) {
        if (!format_dosage_relevant) {
          logerrputs("Warning: No FORMAT/DS or /HDS key found in --vcf file header.  Dosages will not\nbe imported.\n(If these header line(s) are actually present, but with extra spaces or unusual\nfield ordering, standardize the header with e.g. bcftools.)\n");
        } else {
          // since we did find FORMAT/DS, trailing parenthetical is very
          // unlikely to be relevant
          logerrputs("Warning: No FORMAT/HDS key found in --vcf file header.  Dosages will be imported\n(from FORMAT/DS), but phase information will be limited or absent.\n");
        }
      }
    } else if ((!format_dosage_relevant) && dosage_import_field) {
      logerrprintfww("Warning: No FORMAT/%s key found in --vcf file header. Dosages will not be imported. (If this header line is actually present, but with extra spaces or unusual field ordering, standardize the header with e.g. bcftools.)\n", dosage_import_field);
    }
    FinalizeChrset(load_filter_log_import_flags, cip);
    // don't call FinalizeChrInfo here, since this may be followed by --pmerge,
    // etc.

    if (unlikely(!StrStartsWithUnsafe(line_iter, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"))) {
      snprintf(g_logbuf, kLogbufSize, "Error: Header line %" PRIuPTR " of --vcf file does not have expected field sequence after #CHROM.\n", line_idx);
      goto VcfToPgen_ret_MALFORMED_INPUT_WW;
    }
    char* linebuf_iter = &(line_iter[strlen("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")]);
    uint32_t sample_ct = 0;
    if (StrStartsWithUnsafe(linebuf_iter, "\tFORMAT\t")) {
      reterr = VcfSampleLine(preexisting_psamname, const_fid, misc_flags, import_flags, fam_cols, id_delim, idspace_to, 'v', &(linebuf_iter[strlen("\tFORMAT\t")]), outname, outname_end, &sample_ct);
      if (unlikely(reterr)) {
        goto VcfToPgen_ret_1;
      }
    }
    if (unlikely((!sample_ct) && (!no_samples_ok))) {
      logerrputs("Error: No samples in --vcf file.  (This is only permitted when you haven't\nspecified another operation which requires genotype or sample information.)\n");
      goto VcfToPgen_ret_DEGENERATE_DATA;
    }
    vic.vibc.sample_ct = sample_ct;
    // bugfix (5 Jun 2018): must initialize qual_field_ct to zero
    vic.vibc.qual_field_ct = 0;

    uint32_t variant_ct = 0;
    uintptr_t max_variant_ct = bigstack_left() / sizeof(intptr_t);
    if (info_pr_exists) {
      // nonref_flags
      max_variant_ct -= BitCtToAlignedWordCt(max_variant_ct) * kWordsPerVec;
    }
#ifdef __LP64__
    if (max_variant_ct > kPglMaxVariantCt) {
      max_variant_ct = kPglMaxVariantCt;
    }
#endif
    uintptr_t base_chr_present[kChrExcludeWords];
    ZeroWArr(kChrExcludeWords, base_chr_present);

    const uintptr_t header_line_ct = line_idx;
    const uint32_t max_variant_ctaw = BitCtToAlignedWordCt(max_variant_ct);
    // don't need dosage_flags or dphase_flags; dosage overrides GT so slow
    // parse needed
    uintptr_t* nonref_flags = nullptr;
    if (info_pr_exists) {
      nonref_flags = S_CAST(uintptr_t*, bigstack_alloc_raw_rd(max_variant_ctaw * sizeof(intptr_t)));
    }
    uintptr_t* nonref_flags_iter = nonref_flags;
    uintptr_t* allele_idx_offsets = R_CAST(uintptr_t*, g_bigstack_base);
    uintptr_t max_postformat_blen = 1;  // starting from tab at end of FORMAT
    uintptr_t variant_skip_ct = 0;
    uintptr_t nonref_word = 0;
    uintptr_t allele_idx_end = 0;
    uint32_t max_alt_ct = 1;
    uint32_t max_allele_slen = 1;
    uint32_t max_qualfilterinfo_slen = 6;
    uint32_t phase_or_dosage_found = 0;
    uint32_t not_single_sample_no_nonvar = (sample_ct != 1) || format_dosage_relevant || format_hds_search || (import_flags & kfImportVcfAllowNoNonvar);
    uint32_t nonvar_nonmissing_ct = 0;

    while (1) {
      ++line_idx;
      line_iter = AdvPastDelim(line_iter, '\n');
      const uint32_t prev_line_blen = line_iter - prev_line_start;
      if (prev_line_blen > max_line_blen) {
        max_line_blen = prev_line_blen;
      }
      reterr = TextNextLineUnsafe(&vcf_txs, &line_iter);
      if (reterr) {
        if (likely(reterr == kPglRetEof)) {
          // reterr = kPglRetSuccess;
          break;
        }
        goto VcfToPgen_ret_TSTREAM_FAIL;
      }
      prev_line_start = line_iter;
      // we were previously tolerating trailing newlines here, but there wasn't
      // a good reason for doing so.
      if (unlikely(ctou32(*line_iter) <= 32)) {
        if ((*line_iter == ' ') || (*line_iter == '\t')) {
          snprintf(g_logbuf, kLogbufSize, "Error: Leading space or tab on line %" PRIuPTR " of --vcf file.\n", line_idx);
          goto VcfToPgen_ret_MALFORMED_INPUT_2N;
        }
        goto VcfToPgen_ret_MISSING_TOKENS;
      }
      linebuf_iter = line_iter;
      char* chr_code_end = NextPrespace(line_iter);
      if (unlikely(*chr_code_end != '\t')) {
        goto VcfToPgen_ret_MISSING_TOKENS;
      }
      // QUAL/FILTER enforcement is now postponed till .pvar loading.  only
      // other things we do during the scanning pass are (i) count alt alleles,
      // and (ii) check whether any phased genotype calls are present.

      char* pos_end = NextPrespace(chr_code_end);
      if (unlikely(*pos_end != '\t')) {
        goto VcfToPgen_ret_MISSING_TOKENS;
      }

      // may as well check ID length here
      // postpone POS validation till second pass so we only have to parse it
      // once
      char* id_end = NextPrespace(pos_end);
      if (unlikely(*id_end != '\t')) {
        goto VcfToPgen_ret_MISSING_TOKENS;
      }
      if (unlikely(S_CAST(uintptr_t, id_end - pos_end) > kMaxIdBlen)) {
        snprintf(g_logbuf, kLogbufSize, "Error: Invalid ID on line %" PRIuPTR " of --vcf file (max " MAX_ID_SLEN_STR " chars).\n", line_idx);
        goto VcfToPgen_ret_MALFORMED_INPUT_WWN;
      }

      // note REF length
      char* ref_allele_start = &(id_end[1]);
      linebuf_iter = FirstPrespace(ref_allele_start);
      if (unlikely(*linebuf_iter != '\t')) {
        goto VcfToPgen_ret_MISSING_TOKENS;
      }
      uint32_t cur_max_allele_slen = linebuf_iter - ref_allele_start;
      if (unlikely(memchr(ref_allele_start, ',', cur_max_allele_slen) != nullptr)) {
        snprintf(g_logbuf, kLogbufSize, "Error: Invalid REF allele on line %" PRIuPTR " of --vcf file.\n", line_idx);
        goto VcfToPgen_ret_MALFORMED_INPUT_WWN;
      }

      uint32_t alt_ct = 1;
      unsigned char ucc;
      // treat ALT=. as if it were an actual allele for now
      for (; ; ++alt_ct) {
        char* cur_allele_start = ++linebuf_iter;
        ucc = *linebuf_iter;
        if (unlikely((ucc <= ',') && (ucc != '*'))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Invalid alternate allele on line %" PRIuPTR " of --vcf file.\n", line_idx);
          goto VcfToPgen_ret_MALFORMED_INPUT_2N;
        }
        do {
          ucc = *(++linebuf_iter);
          // allow GATK 3.4 <*:DEL> symbolic allele
        } while ((ucc > ',') || (ucc == '*'));
        const uint32_t cur_allele_slen = linebuf_iter - cur_allele_start;
        if (cur_allele_slen > cur_max_allele_slen) {
          cur_max_allele_slen = cur_allele_slen;
        }
        if (ucc != ',') {
          break;
        }
      }

      if (unlikely(ucc != '\t')) {
        snprintf(g_logbuf, kLogbufSize, "Error: Malformed ALT field on line %" PRIuPTR " of --vcf file.\n", line_idx);
        goto VcfToPgen_ret_MALFORMED_INPUT_2N;
      }
      if (alt_ct > max_alt_ct) {
        if (alt_ct >= import_max_allele_ct) {
          ++variant_skip_ct;
          line_iter = linebuf_iter;
          continue;
        }
        max_alt_ct = alt_ct;
      }

      // skip QUAL, FILTER
      char* qual_start_m1 = linebuf_iter;
      for (uint32_t uii = 0; uii != 2; ++uii) {
        linebuf_iter = NextPrespace(linebuf_iter);
        if (unlikely(*linebuf_iter != '\t')) {
          goto VcfToPgen_ret_MISSING_TOKENS;
        }
      }

      // --vcf-require-gt
      char* info_start = &(linebuf_iter[1]);
      char* info_end = FirstPrespace(info_start);
      if (sample_ct) {
        if (unlikely(*info_end != '\t')) {
          goto VcfToPgen_ret_MISSING_TOKENS;
        }
        linebuf_iter = &(info_end[1]);
        vic.vibc.gt_exists = memequal_sk(linebuf_iter, "GT") && ((linebuf_iter[2] == ':') || (linebuf_iter[2] == '\t'));
        if (require_gt && (!vic.vibc.gt_exists)) {
          ++variant_skip_ct;
          line_iter = linebuf_iter;
          continue;
        }
      }
      const uint32_t cur_qualfilterinfo_slen = info_end - qual_start_m1;

      // all converters *do* respect chromosome filters
      // wait till this point to apply it, since we don't want to
      // add a contig name to the hash table unless at least one variant on
      // that contig wasn't filtered out for other reasons.
      uint32_t cur_chr_code;
      reterr = GetOrAddChrCodeDestructive("--vcf file", line_idx, prohibit_extra_chr, line_iter, chr_code_end, cip, &cur_chr_code);
      if (unlikely(reterr)) {
        goto VcfToPgen_ret_1;
      }
      if (!IsSet(cip->chr_mask, cur_chr_code)) {
        ++variant_skip_ct;
        line_iter = info_end;
        continue;
      }
      if (cur_max_allele_slen > max_allele_slen) {
        max_allele_slen = cur_max_allele_slen;
      }
      if (cur_qualfilterinfo_slen > max_qualfilterinfo_slen) {
        max_qualfilterinfo_slen = cur_qualfilterinfo_slen;
      }
      if (cur_chr_code <= cip->max_code) {
        SetBit(cur_chr_code, base_chr_present);
      }

      allele_idx_offsets[variant_ct] = allele_idx_end;
      allele_idx_end += alt_ct + 1;
      const uint32_t variant_idx_lowbits = variant_ct % kBitsPerWord;
      if (info_pr_exists) {
        if (PrInInfo(info_end - info_start, info_start)) {
          nonref_word |= k1LU << variant_idx_lowbits;
        }
        if (variant_idx_lowbits == (kBitsPerWord - 1)) {
          *nonref_flags_iter++ = nonref_word;
          nonref_word = 0;
        }
      }
      if (sample_ct) {
        // linebuf_iter currently points to beginning of FORMAT field
        const char* format_end = FirstPrespace(linebuf_iter);
        if (unlikely(*format_end != '\t')) {
          goto VcfToPgen_ret_MISSING_TOKENS;
        }
        if (phase_or_dosage_found && not_single_sample_no_nonvar) {
          goto VcfToPgen_linescan_done;
        }

        if (format_dosage_relevant) {
          vic.dosage_field_idx = GetVcfFormatPosition(dosage_import_field, linebuf_iter, format_end, dosage_import_field_slen);
        }
        if (format_hds_search) {
          // theoretically possible for HDS to be in VCF header without
          // accompanying DS
          vic.hds_field_idx = GetVcfFormatPosition("HDS", linebuf_iter, format_end, 3);
        }
        if (format_gq_or_dp_relevant) {
          STD_ARRAY_DECL(uint32_t, 2, qual_field_idxs);
          uint32_t qual_field_ct = VcfQualScanInit1(linebuf_iter, format_end, vcf_min_gq, vcf_min_dp, vcf_max_dp, qual_field_idxs);
          // bugfix (5 Jun 2018): must initialize qual_field_ct to zero
          vic.vibc.qual_field_ct = 0;
          if (qual_field_ct) {
            vic.vibc.qual_field_ct = VcfQualScanInit2(qual_field_idxs, ctx.qual_mins, ctx.qual_maxs, vic.vibc.qual_field_skips, vic.vibc.qual_line_mins, vic.vibc.qual_line_maxs);
          }
        }

        // Check if there's at least one phased het call, and/or at least one
        // relevant dosage.
        // If there's only one sample, maybe also check whether the VCF has any
        // hom-REF calls at all.
        // Don't bother multithreading this since it's trivial.
        if ((vic.hds_field_idx != UINT32_MAX) || (vic.dosage_field_idx != UINT32_MAX)) {
          if (alt_ct == 1) {
            vcf_parse_err = VcfScanBiallelicHdsLine(&vic, format_end, &phase_or_dosage_found, &line_iter);
          } else {
            putc_unlocked('\n', stdout);
            logerrputs("Error: --vcf multiallelic dosage import is under development.\n");
            reterr = kPglRetNotYetSupported;
            goto VcfToPgen_ret_1;
          }
          if (unlikely(vcf_parse_err)) {
            goto VcfToPgen_ret_PARSE;
          }
        } else {
          if (!vic.vibc.gt_exists) {
            goto VcfToPgen_linescan_done;
          }
          if (!not_single_sample_no_nonvar) {
            if (format_end[1] == '0') {
              const char cc = format_end[2];
              not_single_sample_no_nonvar = (((cc == '/') || (cc == '|')) && (format_end[3] == '0')) || (cc == ':') || (cc == '\t');
            }
            // This can miss ./0 and the like, but that should practically
            // never matter.
            nonvar_nonmissing_ct += (format_end[1] != '.');
          }
          if (!phase_or_dosage_found) {
            if (alt_ct < 10) {
              phase_or_dosage_found = VcfScanShortallelicLine(&(vic.vibc), format_end, &line_iter);
            } else {
              phase_or_dosage_found = VcfScanLongallelicLine(&(vic.vibc), format_end, &line_iter);
            }
          }
        }
      VcfToPgen_linescan_done:
        line_iter = AdvToDelim(line_iter, '\n');
        const uint32_t cur_postformat_slen = line_iter - format_end;
        if (cur_postformat_slen >= max_postformat_blen) {
          max_postformat_blen = cur_postformat_slen + 1;
        }
      }
      if (unlikely(variant_ct++ == max_variant_ct)) {
#ifdef __LP64__
        if (variant_ct == kPglMaxVariantCt) {
          putc_unlocked('\n', stdout);
          logerrputs("Error: " PROG_NAME_STR " does not support more than 2^31 - 3 variants.  We recommend using\nother software for very deep studies of small numbers of genomes.\n");
          goto VcfToPgen_ret_MALFORMED_INPUT;
        }
#endif
        goto VcfToPgen_ret_NOMEM;
      }
      if (!(variant_ct % 1000)) {
        printf("\r--vcf: %uk variants scanned.", variant_ct / 1000);
        fflush(stdout);
      }
    }
    if (variant_ct % kBitsPerWord) {
      if (nonref_flags_iter) {
        *nonref_flags_iter = nonref_word;
      }
    } else if (unlikely(!variant_ct)) {
      if (!variant_skip_ct) {
        logerrputs("Error: No variants in --vcf file.\n");
        goto VcfToPgen_ret_DEGENERATE_DATA;
      }
      char* write_iter = strcpya_k(g_logbuf, "Error: All ");
      write_iter = wtoa(variant_skip_ct, write_iter);
      write_iter = strcpya_k(write_iter, " variant");
      if (variant_skip_ct != 1) {
        *write_iter++ = 's';
      }
      write_iter = strcpya_k(write_iter, " in --vcf file excluded by ");
      AppendLoadFilterFlagnames(load_filter_log_import_flags, &write_iter);
      strcpy_k(write_iter, ".\n");
      goto VcfToPgen_ret_INCONSISTENT_INPUT_WW;
    }

    putc_unlocked('\r', stdout);
    {
      char* write_iter = strcpya_k(g_logbuf, "--vcf: ");
      write_iter = wtoa(variant_ct + variant_skip_ct, write_iter);
      write_iter = strcpya_k(write_iter, " variant");
      if (variant_ct + variant_skip_ct != 1) {
        *write_iter++ = 's';
      }
      write_iter = strcpya_k(write_iter, " scanned");
      if (variant_skip_ct) {
        write_iter = strcpya_k(write_iter, "; ");
        write_iter = wtoa(variant_skip_ct, write_iter);
        write_iter = strcpya_k(write_iter, " excluded by ");
        AppendLoadFilterFlagnames(load_filter_log_import_flags, &write_iter);
        write_iter = strcpya_k(write_iter, ", ");
        write_iter = u32toa(variant_ct, write_iter);
        write_iter = strcpya_k(write_iter, " remaining");
      } else if (load_filter_log_import_flags) {
        write_iter = strcpya_k(write_iter, " (");
        AppendLoadFilterFlagnames(load_filter_log_import_flags, &write_iter);
        write_iter = strcpya_k(write_iter, " had no effect)");
      }
      strcpy_k(write_iter, ".\n");
      WordWrapB(0);
      logputsb();
    }
    if ((!not_single_sample_no_nonvar) && (nonvar_nonmissing_ct >= 1000)) {
      logerrputs("Error: All genotypes in single-sample VCF contain at least one ALT allele; this\nimplies the VCF was incorrectly generated.  You probably need to backtrack and\ne.g. rerun GATK GenotypeGVCFs with the --include-non-variant-sites flag added.\n");
      goto VcfToPgen_ret_INCONSISTENT_INPUT;
    }

    if (allele_idx_end > 2 * variant_ct) {
      allele_idx_offsets[variant_ct] = allele_idx_end;
      BigstackFinalizeW(allele_idx_offsets, variant_ct + 1);
    } else {
      allele_idx_offsets = nullptr;
    }

    // Close file, then reopen with a smaller line-load buffer and (if bgzf)
    // reduce decompression thread count.  2 is good in the simplest cases
    // (no GQ/DP filter, no dosage), otherwise limit to 1.
    uint32_t decompress_thread_ct = 1;
    uint32_t calc_thread_ct;
    {
      if ((vcf_min_gq != -1) || (vcf_min_dp != -1) || (phase_or_dosage_found && (format_dosage_relevant || format_hds_search))) {
        // "are lines expensive to parse?"  will add a multiallelic condition
        // to the disjunction soon
        // this is based on a bunch of DS-force measurements
        calc_thread_ct = 1 + (sample_ct > 5) + (sample_ct > 12) + (sample_ct > 32) + (sample_ct > 512);
      } else {
        if (TextIsMt(&vcf_txs) && (max_thread_ct > 1)) {
          decompress_thread_ct = 2;
        }
        // this seems to saturate around 3 threads.
        calc_thread_ct = 1 + (sample_ct > 40) + (sample_ct > 320);
      }
      if (unlikely(CleanupTextStream2(vcfname, &vcf_txs, &reterr))) {
        goto VcfToPgen_ret_1;
      }
      BigstackEndReset(bigstack_end_mark);
      reterr = InitTextStreamEx(vcfname, 1, kMaxLongLine, MAXV(max_line_blen, kTextStreamBlenFast), decompress_thread_ct, &vcf_txs);
      if (unlikely(reterr)) {
        goto VcfToPgen_ret_TSTREAM_FAIL;
      }
      if (calc_thread_ct + decompress_thread_ct > max_thread_ct) {
        calc_thread_ct = MAXV(1, max_thread_ct - decompress_thread_ct);
      }
    }

    const uint32_t variant_ctl = BitCtToWordCt(variant_ct);
    PgenGlobalFlags phase_dosage_gflags = kfPgenGlobal0;
    GparseFlags gparse_flags = kfGparse0;  // yeah, this is a bit redundant
    // bugfix (22 Jun 2018): if dosage= was specified, we need to reserve
    // phasepresent/dosage_present buffers for each variant even when we know
    // they'll be irrelevant, since they're needed during conversion.
    if (phase_or_dosage_found || format_hds_search || format_dosage_relevant) {
      if ((!format_hds_search) && (!format_dosage_relevant)) {
        // TODO: if phase_or_dosage_found is false, why do we still have to set
        // phase_dosage_gflags here?  with this behavior, what was the point of
        // scanning the lines in the first pass?
        phase_dosage_gflags = kfPgenGlobalHardcallPhasePresent;
        gparse_flags = kfGparseHphase;
      } else {
        // thanks to automatic --hard-call-threshold, we may need to save
        // dosage-phase even when there's no HDS field
        phase_dosage_gflags = kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent;
        gparse_flags = kfGparseHphase | kfGparseDosage | kfGparseDphase;
      }
    }
    uint32_t nonref_flags_storage = 1;
    if (nonref_flags) {
      const uint32_t variant_ctl_m1 = variant_ctl - 1;
      const uintptr_t last_nonref_flags_word = nonref_flags[variant_ctl_m1];
      if (!last_nonref_flags_word) {
        for (uint32_t widx = 0; widx != variant_ctl_m1; ++widx) {
          if (nonref_flags[widx]) {
            nonref_flags_storage = 3;
            break;
          }
        }
      } else if (!((~last_nonref_flags_word) << ((-variant_ct) & (kBitsPerWord - 1)))) {
        nonref_flags_storage = 2;
        for (uint32_t widx = 0; widx != variant_ctl_m1; ++widx) {
          if (~nonref_flags[widx]) {
            nonref_flags_storage = 3;
            break;
          }
        }
      } else {
        nonref_flags_storage = 3;
      }
      if (nonref_flags_storage != 3) {
        // yeah, we may now have a temporary "memory leak" here (if
        // multiallelic variants are present, and thus allele_idx_offsets[]
        // must be kept around), but this array is typically only 1/64 the size
        // of allele_idx_offsets[].
        if (!allele_idx_offsets) {
          BigstackReset(nonref_flags);
        }
        nonref_flags = nullptr;
      }
    }
    snprintf(outname_end, kMaxOutfnameExtBlen, ".pvar");
    const uint32_t output_zst = ((import_flags & (kfImportKeepAutoconv | kfImportKeepAutoconvVzs)) != kfImportKeepAutoconv);
    if (output_zst) {
      snprintf(&(outname_end[5]), kMaxOutfnameExtBlen - 5, ".zst");
    }
    const uintptr_t overflow_buf_size = kCompressStreamBlock + MAXV(2 * max_allele_slen + max_qualfilterinfo_slen + kMaxIdSlen + 32, kCompressStreamBlock);
    reterr = InitCstreamAlloc(outname, 0, output_zst, sample_ct? 1 : MAXV(1, max_thread_ct - decompress_thread_ct), overflow_buf_size, &pvar_css, &pvar_cswritep);
    if (unlikely(reterr)) {
      goto VcfToPgen_ret_1;
    }
    for (line_idx = 1, line_iter = TextLineEnd(&vcf_txs); ; ++line_idx, line_iter = AdvPastDelim(line_iter, '\n')) {
      reterr = TextNextLineUnsafe(&vcf_txs, &line_iter);
      if (unlikely(reterr)) {
        goto VcfToPgen_ret_TSTREAM_FAIL;
      }
      if (line_idx == header_line_ct) {
        break;
      }
      // chrSet skipped here since we call AppendChrsetLine after this loop
      if (StrStartsWithUnsafe(line_iter, "##fileformat=") || StrStartsWithUnsafe(line_iter, "##fileDate=") || StrStartsWithUnsafe(line_iter, "##source=") || StrStartsWithUnsafe(line_iter, "##FORMAT=") || StrStartsWithUnsafe(line_iter, "##chrSet=")) {
        continue;
      }
      if (StrStartsWithUnsafe(line_iter, "##contig=<ID=")) {
        char* contig_name_start = &(line_iter[strlen("##contig=<ID=")]);
        char* contig_name_end = strchrnul_n(contig_name_start, ',');
        if (*contig_name_end != ',') {
          contig_name_end = Memrchr(contig_name_start, '>', contig_name_end - contig_name_start);
          if (unlikely(!contig_name_end)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Header line %" PRIuPTR " of --vcf file does not have expected ##contig format.\n", line_idx);
            goto VcfToPgen_ret_MALFORMED_INPUT_WW;
          }
        }
        const uint32_t cur_chr_code = GetChrCodeCounted(cip, contig_name_end - contig_name_start, contig_name_start);
        if (IsI32Neg(cur_chr_code)) {
          continue;
        }
        if (cur_chr_code <= cip->max_code) {
          if (!IsSet(base_chr_present, cur_chr_code)) {
            continue;
          }
        } else {
          if (!IsSet(cip->chr_mask, cur_chr_code)) {
            continue;
          }
        }
        // Note that, when --output-chr is specified, we don't update the
        // ##contig header line chromosome code in the .pvar file, since
        // ##contig is not an explicit part of the .pvar specification, it's
        // just another blob of text as far as the main body of plink2 is
        // concerned.  However, the codes are brought in sync during VCF/BCF
        // export.
      }
      // force OS-appropriate eoln
      char* line_last = AdvToDelim(line_iter, '\n');
#ifdef _WIN32
      if (line_last[-1] == '\r') {
        --line_last;
      }
      // NOT safe to use AppendBinaryEoln here.
      if (unlikely(CsputsStd(line_iter, line_last - line_iter, &pvar_css, &pvar_cswritep))) {
        goto VcfToPgen_ret_WRITE_FAIL;
      }
      pvar_cswritep = strcpya_k(pvar_cswritep, "\r\n");
#else
      char* line_write_end;
      if (line_last[-1] == '\r') {
        line_write_end = line_last;
        line_last[-1] = '\n';
      } else {
        line_write_end = &(line_last[1]);
      }
      if (unlikely(CsputsStd(line_iter, line_write_end - line_iter, &pvar_css, &pvar_cswritep))) {
        goto VcfToPgen_ret_WRITE_FAIL;
      }
#endif
      line_iter = line_last;
    }
    if (cip->chrset_source) {
      AppendChrsetLine(cip, &pvar_cswritep);
    }
    pvar_cswritep = strcpya_k(pvar_cswritep, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER");
    if (info_nonpr_exists) {
      pvar_cswritep = strcpya_k(pvar_cswritep, "\tINFO");
    }
    AppendBinaryEoln(&pvar_cswritep);
    if (unlikely(Cswrite(&pvar_css, &pvar_cswritep))) {
      goto VcfToPgen_ret_WRITE_FAIL;
    }
    if (max_alt_ct > kPglMaxAltAlleleCt) {
      logerrprintfww("Error: VCF file has a variant with %u ALT alleles; this build of " PROG_NAME_STR " is limited to " PGL_MAX_ALT_ALLELE_CT_STR ". (You can use \"--import-max-alleles " PGL_MAX_ALLELE_CT_STR "\" to filter out such variants.)\n", max_alt_ct);
      reterr = kPglRetNotYetSupported;
      goto VcfToPgen_ret_1;
    }
    const uint32_t max_allele_ct = max_alt_ct + 1;
    // may as well have a functional progress meter in no-samples case
    uint32_t main_block_size = MINV(65536, variant_ct);
    uint32_t per_thread_block_limit = main_block_size;
    uint32_t cur_thread_block_vidx_limit = 1;
    uintptr_t per_thread_byte_limit = 0;
    unsigned char* geno_bufs[2];
    // defensive
    geno_bufs[0] = nullptr;
    geno_bufs[1] = nullptr;
    if (hard_call_thresh == UINT32_MAX) {
      hard_call_thresh = kDosageMid / 10;
    }
    const uint32_t hard_call_halfdist = kDosage4th - hard_call_thresh;
    if (sample_ct) {
      snprintf(outname_end, kMaxOutfnameExtBlen, ".pgen");
      uintptr_t spgw_alloc_cacheline_ct;
      uint32_t max_vrec_len;
      reterr = SpgwInitPhase1(outname, allele_idx_offsets, nonref_flags, variant_ct, sample_ct, max_allele_ct, kPgenWriteBackwardSeek, phase_dosage_gflags, nonref_flags_storage, &spgw, &spgw_alloc_cacheline_ct, &max_vrec_len);
      if (unlikely(reterr)) {
        if (reterr == kPglRetOpenFail) {
          logerrprintfww(kErrprintfFopen, outname, strerror(errno));
        }
        goto VcfToPgen_ret_1;
      }
      unsigned char* spgw_alloc;
      if (unlikely(bigstack_alloc_uc(spgw_alloc_cacheline_ct * kCacheline, &spgw_alloc))) {
        goto VcfToPgen_ret_NOMEM;
      }
      SpgwInitPhase2(max_vrec_len, &spgw, spgw_alloc);

      if (unlikely(bigstack_alloc_ucp(calc_thread_ct, &ctx.thread_wkspaces) ||
                   bigstack_alloc_u32(calc_thread_ct + 1, &(ctx.thread_bidxs[0])) ||
                   bigstack_alloc_u32(calc_thread_ct + 1, &(ctx.thread_bidxs[1])) ||
                   bigstack_calloc_w(calc_thread_ct, &ctx.err_line_idxs))) {
        goto VcfToPgen_ret_NOMEM;
      }
      ctx.vcf_parse_errs = S_CAST(VcfParseErr*, bigstack_alloc_raw_rd(calc_thread_ct * sizeof(VcfParseErr)));
      if (unlikely(!ctx.vcf_parse_errs)) {
        goto VcfToPgen_ret_NOMEM;
      }
      ctx.sample_ct = sample_ct;
      ctx.hard_call_halfdist = hard_call_halfdist;
      ctx.dosage_erase_halfdist = vic.dosage_erase_halfdist;
      ctx.import_dosage_certainty = import_dosage_certainty;
      ctx.parse_failed = 0;
      // defensive
      ctx.gparse[0] = nullptr;
      ctx.gparse[1] = nullptr;
      ctx.block_allele_idx_offsets[0] = nullptr;
      ctx.block_allele_idx_offsets[1] = nullptr;
      // Finished with all other memory allocations, so all remaining workspace
      // can be spent on multithreaded parsing.  Spend up to 1/6 on
      // g_thread_wkspaces (tune this fraction later).
      // Probable todo: factor out common parts with bgen-1.3 initialization
      // into separate function(s).
      uint64_t max_write_byte_ct = GparseWriteByteCt(sample_ct, max_allele_ct, gparse_flags);
      // always allocate tmp_dphase_delta for now
      uint64_t thread_wkspace_cl_ct = DivUp(max_write_byte_ct + sample_ct * sizeof(SDosage), kCacheline);
      uintptr_t cachelines_avail = bigstack_left() / (6 * kCacheline);
      if (calc_thread_ct * thread_wkspace_cl_ct > cachelines_avail) {
        if (unlikely(thread_wkspace_cl_ct > cachelines_avail)) {
          goto VcfToPgen_ret_NOMEM;
        }
        calc_thread_ct = cachelines_avail / thread_wkspace_cl_ct;
      }
      for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
        ctx.thread_wkspaces[tidx] = S_CAST(unsigned char*, bigstack_alloc_raw(thread_wkspace_cl_ct * kCacheline));
        ctx.vcf_parse_errs[tidx] = kVcfParseOk;
      }

      // be pessimistic re: rounding
      cachelines_avail = (bigstack_left() / kCacheline) - 4;
      const uint64_t max_bytes_req_per_variant = sizeof(GparseRecord) + MAXV(max_postformat_blen, max_write_byte_ct) + calc_thread_ct;
      if (unlikely(cachelines_avail * kCacheline < 2 * max_bytes_req_per_variant)) {
        goto VcfToPgen_ret_NOMEM;
      }
      // use worst-case gparse_flags since lines will usually be similar
      uintptr_t min_bytes_req_per_variant = sizeof(GparseRecord) + GparseWriteByteCt(sample_ct, 2, gparse_flags);
      main_block_size = (cachelines_avail * kCacheline) / (min_bytes_req_per_variant * 2);
      // this is arbitrary, there's no connection to kPglVblockSize
      if (main_block_size > 65536) {
        main_block_size = 65536;
      }
      // divide by 2 for better parallelism in small-variant-count case
      // round up per_thread_block_limit so we only have two blocks
      if (main_block_size > DivUp(variant_ct, 2)) {
        main_block_size = DivUp(variant_ct, 2) + calc_thread_ct - 1;
      }
      // may as well guarantee divisibility
      per_thread_block_limit = main_block_size / calc_thread_ct;
      main_block_size = per_thread_block_limit * calc_thread_ct;
      if (unlikely(SetThreadCt(calc_thread_ct, &tg))) {
        goto VcfToPgen_ret_NOMEM;
      }
      ctx.gparse[0] = S_CAST(GparseRecord*, bigstack_alloc_raw_rd(main_block_size * sizeof(GparseRecord)));
      ctx.gparse[1] = S_CAST(GparseRecord*, bigstack_alloc_raw_rd(main_block_size * sizeof(GparseRecord)));
      SetThreadFuncAndData(VcfGenoToPgenThread, &ctx, &tg);
      cachelines_avail = bigstack_left() / (kCacheline * 2);
      geno_bufs[0] = S_CAST(unsigned char*, bigstack_alloc_raw(cachelines_avail * kCacheline));
      geno_bufs[1] = S_CAST(unsigned char*, bigstack_alloc_raw(cachelines_avail * kCacheline));
      // This is only used for comparison purposes, so it is unnecessary to
      // round it down to a multiple of kBytesPerVec even though every actual
      // record will be vector-aligned.
      per_thread_byte_limit = (cachelines_avail * kCacheline) / calc_thread_ct;
    }

    const uint32_t can_fail_on_ds_only = strequal_k(dosage_import_field, "DS", dosage_import_field_slen) && (!require_gt);
    uint32_t fail_on_ds_only = 0;

    // If --lax-chrx-import wasn't specified, we now error out when all of the
    // following conditions hold:
    // 1. At least one sample is present.
    // 2a. Either no .fam/.psam or --update-sex file was provided, --impute-sex
    //     was not specified, and there's at least one chrX variant, or
    // 2b. The basic chromosome set is default/human, --split-par wasn't
    //     specified, and there's a variant in the intersection of the b37 and
    //     b38 PARs.
    //
    // Without this error, it's way too easy for users to process chrX
    // incorrectly without realizing it, especially since VCF files directly
    // encode ploidy information and it's not immediately obvious why plink2
    // fails to read that.
    //
    // If --lax-chrx-import wasn't specified, and conditions (1) and (2a) hold,
    // we can error out immediately.
    // If conditions (1) holds, (2a) doesn't hold, and (2b) might hold, we
    // initialize par_warn_code, and in the main loop we error out if we
    // encounter a chrX POS in the intersection of the b37 and b38 PARs.
    // (Update, 11 Apr 2025: we reduce the latter error to a warning when
    // --sort-vars is part of the command line.  Otherwise, --lax-chrx-import
    // would become necessary when dealing with an unsorted VCF.)

    // Could also require the .fam/.psam to contain non-NA sex information?
    // But this should already be enough to address the main footgun.
    const uint32_t sex_info_avail = preexisting_psamname || is_update_or_impute_sex;
    uint32_t par_warn_code = UINT32_MAX;
    uint32_t print_splitpar_warning = 0;
    if ((!(import_flags & kfImportLaxChrX)) && sample_ct) {
      const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
      if ((!IsI32Neg(x_code)) && IsSet(base_chr_present, x_code)) {
        if (unlikely(!sex_info_avail)) {
          logerrputs("Error: chrX is present in the input file, but no sex information was provided;\nrerun this import with --psam, --update-sex, or --impute-sex.  --split-par may\nalso be appropriate.\n");
          goto VcfToPgen_ret_INCONSISTENT_INPUT;
        } else if (IsHumanChrset(cip) && (!is_splitpar)) {
          par_warn_code = x_code;
        }
      }
    }

    // Main workflow:
    // 1. Set n=0, load genotype data for first main_block_size variants
    //    while writing .pvar
    //
    // 2. Spawn threads processing batch n genotype data
    // 3. If n>0, write results for block (n-1)
    // 4. Increment n by 1
    // 5. Load/write-.pvar for batch (n+1) unless eof
    // 6. Join threads
    // 7. Goto step 2 unless eof
    //
    // 8. Write results for last block
    uint32_t prev_block_write_ct = 0;
    uint32_t genotext_byte_ct = 0;
    uintptr_t record_byte_ct = 0;
    uint32_t allele_ct = 0;
    uint32_t* thread_bidxs = nullptr;
    GparseRecord* cur_gparse = nullptr;
    unsigned char* geno_buf_iter = nullptr;
    unsigned char* cur_thread_byte_stop = nullptr;
    uint32_t parity = 0;
    for (uint32_t vidx_start = 0; ; ) {
      uint32_t cur_block_write_ct = 0;
      if (!IsLastBlock(&tg)) {
        const uint32_t block_vidx_limit = variant_ct - vidx_start;
        cur_thread_block_vidx_limit = MINV(block_vidx_limit, per_thread_block_limit);
        uint32_t cur_thread_fill_idx = 0;
        if (sample_ct) {
          thread_bidxs = ctx.thread_bidxs[parity];
          cur_gparse = ctx.gparse[parity];
          if (allele_idx_offsets) {
            ctx.block_allele_idx_offsets[parity] = &(allele_idx_offsets[vidx_start]);
          }
          geno_buf_iter = geno_bufs[parity];
          cur_thread_byte_stop = &(geno_buf_iter[per_thread_byte_limit]);
          thread_bidxs[0] = 0;
        }
        uint32_t block_vidx = 0;
        GparseRecord* grp;
        if (!genotext_byte_ct) {
          goto VcfToPgen_load_start;
        }
        // we may stop before main_block_size due to insufficient space in
        // geno_bufs[parity].  if so, we copy over the post-FORMAT part of the
        // current line before proceeding.
        while (1) {
          grp = &(cur_gparse[block_vidx]);
          grp->record_start = geno_buf_iter;
          grp->flags = gparse_flags;
          STD_ARRAY_COPY(vic.vibc.qual_field_skips, 2, grp->metadata.read_vcf.qual_field_idxs);
          grp->metadata.read_vcf.gt_exists = vic.vibc.gt_exists;
          grp->metadata.read_vcf.qual_exists = vic.vibc.qual_field_ct;
          grp->metadata.read_vcf.dosage_field_idx = vic.dosage_field_idx;
          grp->metadata.read_vcf.hds_field_idx = vic.hds_field_idx;
          grp->metadata.read_vcf.line_idx = line_idx;
          if (gparse_flags != kfGparseNull) {
            memcpy(geno_buf_iter, linebuf_iter, genotext_byte_ct);
          }
          geno_buf_iter = &(geno_buf_iter[record_byte_ct]);
          ++block_vidx;

          // true iff this is the last variant we're keeping in the entire
          // file
          if (block_vidx == block_vidx_limit) {
            for (; cur_thread_fill_idx != calc_thread_ct; ) {
              // save endpoint for current thread, and tell any leftover
              // threads to do nothing
              thread_bidxs[++cur_thread_fill_idx] = block_vidx;
            }
            break;
          }
        VcfToPgen_load_start:
          ++line_idx;
          line_iter = AdvPastDelim(line_iter, '\n');
          // In principle, it shouldn't be necessary to check the exact value
          // of reterr, but this may be useful for bug investigation.
          reterr = TextNextLineUnsafe(&vcf_txs, &line_iter);
          if (unlikely(reterr)) {
            goto VcfToPgen_ret_TSTREAM_FAIL;
          }

          // 1. check if we skip this variant.  chromosome filter,
          //    import_max_allele_ct, and require_gt can cause this.
          char* chr_code_end = AdvToDelim(line_iter, '\t');
          uint32_t chr_code_base = GetChrCodeRaw(line_iter);
          if (chr_code_base == UINT32_MAX) {
            // skip hash table lookup if we know we aren't skipping the variant
            if (variant_skip_ct) {
              *chr_code_end = '\0';
              // can't overread, nonstd_names not in main workspace
              const uint32_t chr_code = IdHtableFind(line_iter, TO_CONSTCPCONSTP(cip->nonstd_names), cip->nonstd_id_htable, chr_code_end - line_iter, kChrHtableSize);
              if ((chr_code == UINT32_MAX) || (!IsSet(cip->chr_mask, chr_code))) {
                line_iter = chr_code_end;
                goto VcfToPgen_load_start;
              }
              *chr_code_end = '\t';
            }
            if (can_fail_on_ds_only) {
              fail_on_ds_only = cip->haploid_mask[0] & 1;
            }
          } else {
            if (chr_code_base >= kMaxContigs) {
              chr_code_base = cip->xymt_codes[chr_code_base - kMaxContigs];
            }
            if (IsI32Neg(chr_code_base) || (!IsSet(base_chr_present, chr_code_base))) {
              assert(variant_skip_ct);
              line_iter = chr_code_end;
              goto VcfToPgen_load_start;
            }
            if (can_fail_on_ds_only) {
              fail_on_ds_only = IsSet(cip->haploid_mask, chr_code_base);
            }
          }
          // chr_code_base is now a proper numeric chromosome index for
          // non-contigs, and UINT32_MAX if it's a contig name
          char* pos_str = &(chr_code_end[1]);
          char* pos_str_end = AdvToDelim(pos_str, '\t');
          // copy ID, REF verbatim...
          linebuf_iter = AdvToNthDelim(&(pos_str_end[1]), 2, '\t');
          if (ref_n_missing && memequal_sk(&(linebuf_iter[-2]), "\tN")) {
            // ...unless --vcf-ref-n-missing applies.
            linebuf_iter[-1] = '.';
          }

          // ALT, QUAL, FILTER, INFO
          char* alt_start = &(linebuf_iter[1]);
          char* alt_end = AdvToDelim(alt_start, '\t');
          if ((import_max_allele_ct < 0x7ffffffe) && variant_skip_ct) {
            const uint32_t alt_ct = 1 + CountByte(alt_start, ',', alt_end - alt_start);
            if (alt_ct >= import_max_allele_ct) {
              line_iter = alt_end;
              goto VcfToPgen_load_start;
            }
          }
          char* filter_end = AdvToNthDelim(&(alt_end[1]), 2, '\t');
          char* format_start = nullptr;
          char* info_end;
          if (sample_ct) {
            info_end = AdvToDelim(&(filter_end[1]), '\t');
            format_start = &(info_end[1]);
            vic.vibc.gt_exists = memequal_sk(format_start, "GT") && ((format_start[2] == ':') || (format_start[2] == '\t'));
            if (require_gt && (!vic.vibc.gt_exists)) {
              line_iter = format_start;
              goto VcfToPgen_load_start;
            }
          } else {
            info_end = NextPrespace(filter_end);
          }
          // No variant-skipping past this point; safe to start writing to
          // pvar_cswritep.

          // make sure POS starts with an integer, apply --output-chr setting
          uint32_t cur_bp;
          if (unlikely(ScanUintDefcap(pos_str, &cur_bp))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid POS on line %" PRIuPTR " of --vcf file.\n", line_idx);
            goto VcfToPgen_ret_MALFORMED_INPUT_2N;
          }

          if (chr_code_base == UINT32_MAX) {
            pvar_cswritep = memcpya(pvar_cswritep, line_iter, chr_code_end - line_iter);
          } else {
            if (par_warn_code == chr_code_base) {
              if ((cur_bp <= kPAR1IntersectionLast) || (cur_bp >= kPAR2IntersectionFirst)) {
                if (unlikely(!is_sortvars)) {
                  putc_unlocked('\n', stdout);
                  logerrputs("Error: Human chrX pseudoautosomal variant(s) appear to be present in the input\nVCF, but --split-par was not specified.\n");
                  goto VcfToPgen_ret_INCONSISTENT_INPUT;
                }
                print_splitpar_warning = 1;
                par_warn_code = UINT32_MAX;
              }
            }
            pvar_cswritep = chrtoa(cip, chr_code_base, pvar_cswritep);
          }
          *pvar_cswritep++ = '\t';
          pvar_cswritep = u32toa(cur_bp, pvar_cswritep);

          // first copy includes both REF and ALT1
          char* copy_start = pos_str_end;
          uint32_t alt_ct;
          for (alt_ct = 1; ; ++alt_ct) {
            ++linebuf_iter;
            unsigned char ucc;
            do {
              ucc = *(++linebuf_iter);
              // allow GATK 3.4 <*:DEL> symbolic allele
            } while ((ucc > ',') || (ucc == '*'));
            pvar_cswritep = memcpya(pvar_cswritep, copy_start, linebuf_iter - copy_start);
            if (unlikely(Cswrite(&pvar_css, &pvar_cswritep))) {
              goto VcfToPgen_ret_WRITE_FAIL;
            }
            if (ucc != ',') {
              break;
            }
            copy_start = linebuf_iter;
          }

          if (info_nonpr_exists) {
            // VCF specification permits whitespace in INFO field, while PVAR
            // does not.  Check for whitespace and error out if necessary.
            if (unlikely(memchr(filter_end, ' ', info_end - filter_end))) {
              snprintf(g_logbuf, kLogbufSize, "Error: INFO field on line %" PRIuPTR " of --vcf file contains a space; this cannot be imported by " PROG_NAME_STR ". Remove or reformat the field before reattempting import.\n", line_idx);
              goto VcfToPgen_ret_MALFORMED_INPUT_WWN;
            }
            pvar_cswritep = memcpya(pvar_cswritep, linebuf_iter, info_end - linebuf_iter);
          } else {
            pvar_cswritep = memcpya(pvar_cswritep, linebuf_iter, filter_end - linebuf_iter);
          }
          AppendBinaryEoln(&pvar_cswritep);
          if (!sample_ct) {
            if (++block_vidx == cur_thread_block_vidx_limit) {
              break;
            }
            line_iter = info_end;
            goto VcfToPgen_load_start;
          }
          if ((!vic.vibc.gt_exists) && (!format_dosage_relevant) && (!format_hds_search)) {
            line_iter = AdvToDelim(format_start, '\n');
            if (unlikely(CountByte(format_start, '\t', line_iter - format_start) != sample_ct)) {
              goto VcfToPgen_ret_MISSING_TOKENS;
            }
            gparse_flags = kfGparseNull;
            genotext_byte_ct = 1;
          } else {
            linebuf_iter = AdvToDelim(format_start, '\t');
            if (format_gq_or_dp_relevant) {
              vic.vibc.qual_field_ct = VcfQualScanInit1(format_start, linebuf_iter, vcf_min_gq, vcf_min_dp, vcf_max_dp, vic.vibc.qual_field_skips);
            }
            if ((!phase_or_dosage_found) && (!format_dosage_relevant) && (!format_hds_search)) {
              gparse_flags = kfGparse0;
            } else {
              if (format_dosage_relevant) {
                vic.dosage_field_idx = GetVcfFormatPosition(dosage_import_field, format_start, linebuf_iter, dosage_import_field_slen);
              }
              if (format_hds_search) {
                vic.hds_field_idx = GetVcfFormatPosition("HDS", format_start, linebuf_iter, 3);
              }
              if (unlikely(fail_on_ds_only && (!vic.vibc.gt_exists) && (vic.dosage_field_idx != UINT32_MAX) && (vic.hds_field_idx == UINT32_MAX))) {
                // could allow the chrX all-female case, but let's keep this
                // simple for now
                snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of --vcf file is for a chrX, chrM, or fully-haploid variant, and has a DS field without a companion GT field to clarify whether each DS value is on a 0..1 or 0..2 scale. This cannot be imported by " PROG_NAME_STR "; please e.g. regenerate the file with GT present.\n", line_idx);
                goto VcfToPgen_ret_MALFORMED_INPUT_WWN;
              }
              gparse_flags = ((vic.dosage_field_idx != UINT32_MAX) || (vic.hds_field_idx != UINT32_MAX))? (kfGparseHphase | kfGparseDosage | kfGparseDphase) : kfGparseHphase;
            }
            line_iter = AdvToDelim(linebuf_iter, '\n');
            genotext_byte_ct = 1 + S_CAST(uintptr_t, line_iter - linebuf_iter);
          }
          allele_ct = alt_ct + 1;
          const uintptr_t write_byte_ct_limit = GparseWriteByteCt(sample_ct, allele_ct, gparse_flags);
          record_byte_ct = MAXV(RoundUpPow2(genotext_byte_ct, kBytesPerVec), write_byte_ct_limit);

          if ((block_vidx == cur_thread_block_vidx_limit) || (S_CAST(uintptr_t, cur_thread_byte_stop - geno_buf_iter) < record_byte_ct)) {
            thread_bidxs[++cur_thread_fill_idx] = block_vidx;
            if (cur_thread_fill_idx == calc_thread_ct) {
              break;
            }
            cur_thread_byte_stop = &(cur_thread_byte_stop[per_thread_byte_limit]);
            cur_thread_block_vidx_limit = MINV(cur_thread_block_vidx_limit + per_thread_block_limit, block_vidx_limit);
          }
        }
        cur_block_write_ct = block_vidx;
      }
      if (sample_ct) {
        if (vidx_start) {
          JoinThreads(&tg);
          if (unlikely(ctx.parse_failed)) {
            goto VcfToPgen_ret_THREAD_PARSE;
          }
        }
        if (!IsLastBlock(&tg)) {
          if (vidx_start + cur_block_write_ct == variant_ct) {
            DeclareLastThreadBlock(&tg);
          }
          if (unlikely(SpawnThreads(&tg))) {
            goto VcfToPgen_ret_THREAD_CREATE_FAIL;
          }
        }
        parity = 1 - parity;
        if (vidx_start) {
          // write *previous* block results
          reterr = GparseFlush(ctx.gparse[parity], allele_idx_offsets, prev_block_write_ct, &spgw);
          if (unlikely(reterr)) {
            goto VcfToPgen_ret_1;
          }
        }
      } else if (vidx_start + cur_block_write_ct == variant_ct) {
        break;
      }
      if (vidx_start == variant_ct) {
        break;
      }
      if (vidx_start) {
        printf("\r--vcf: %uk variants converted.", vidx_start / 1000);
        if (vidx_start <= main_block_size) {
          fputs("    \b\b\b\b", stdout);
        }
        fflush(stdout);
      }
      vidx_start += cur_block_write_ct;
      prev_block_write_ct = cur_block_write_ct;
    }
    if (unlikely(CswriteCloseNull(&pvar_css, pvar_cswritep))) {
      goto VcfToPgen_ret_WRITE_FAIL;
    }
    if (sample_ct) {
      reterr = SpgwFinish(&spgw);
      if (unlikely(reterr)) {
        goto VcfToPgen_ret_1;
      }
    }
    putc_unlocked('\r', stdout);
    char* write_iter = strcpya_k(g_logbuf, "--vcf: ");
    const uint32_t outname_base_slen = outname_end - outname;
    if (sample_ct) {
      write_iter = memcpya(write_iter, outname, outname_base_slen + 5);
      write_iter = strcpya_k(write_iter, " + ");
    } else {
      *pgen_generated_ptr = 0;
    }
    write_iter = memcpya(write_iter, outname, outname_base_slen);
    write_iter = strcpya_k(write_iter, ".pvar");
    if (output_zst) {
      write_iter = strcpya_k(write_iter, ".zst");
    }
    if (sample_ct && (!preexisting_psamname)) {
      write_iter = strcpya_k(write_iter, " + ");
      write_iter = memcpya(write_iter, outname, outname_base_slen);
      write_iter = strcpya_k(write_iter, ".psam");
    } else {
      *psam_generated_ptr = 0;
    }
    write_iter = strcpya_k(write_iter, " written");
    if (!sample_ct) {
      write_iter = strcpya_k(write_iter, " (no samples present)");
    }
    strcpy_k(write_iter, ".\n");
    WordWrapB(0);
    logputsb();
    if (print_splitpar_warning) {
      logerrputs("Warning: Human chrX pseudoautosomal variant(s) appear to be present in the\ninput VCF.  You probably want to include --split-par in your next command.\n");
    }
  }
  while (0) {
  VcfToPgen_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  VcfToPgen_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  VcfToPgen_ret_TSTREAM_FAIL:
    putc_unlocked('\n', stdout);
    TextStreamErrPrint("--vcf file", &vcf_txs);
    break;
  VcfToPgen_ret_THREAD_PARSE:
    {
      // Doesn't really cost us anything to report the first error if multiple
      // converter threads error out in the same block, though we can't
      // generally guarantee that we'll report the first error in the file.
      for (uint32_t tidx = 0; ; ++tidx) {
        vcf_parse_err = ctx.vcf_parse_errs[tidx];
        if (vcf_parse_err) {
          line_idx = ctx.err_line_idxs[tidx];
          break;
        }
      }
    }
  VcfToPgen_ret_PARSE:
    if (vcf_parse_err == kVcfParseInvalidGt) {
      putc_unlocked('\n', stdout);
      logerrprintf("Error: Line %" PRIuPTR " of --vcf file has an invalid GT field.\n", line_idx);
      reterr = kPglRetMalformedInput;
      break;
    } else if (vcf_parse_err == kVcfParseHalfCallError) {
      putc_unlocked('\n', stdout);
      logerrprintf("Error: Line %" PRIuPTR " of --vcf file has a GT half-call.\n", line_idx);
      if (!half_call_explicit_error) {
        logerrputs("Use --vcf-half-call to specify how these should be processed.\n");
      }
      reterr = kPglRetMalformedInput;
      break;
    } else if (vcf_parse_err == kVcfParseInvalidDosage) {
      // probable todo: distinguish HDS errors (right now, it just prints
      // "invalid DS field" on all HDS errors).
      putc_unlocked('\n', stdout);
      logerrprintfww("Error: Line %" PRIuPTR " of --vcf file has an invalid %s%s field.\n", line_idx, dosage_import_field, format_hds_search? " or HDS" : "");
      reterr = kPglRetInconsistentInput;
      break;
    } else if (vcf_parse_err == kVcfParsePolyploidError) {
      putc_unlocked('\n', stdout);
      logerrprintfww("Error: Line %" PRIuPTR " of --vcf file has a polyploid genotype.%s\n", line_idx, (import_flags & kfImportPolyploidExplicitError)? "" : " (Use '--polyploid-mode missing' to treat these as missing values.)");
      reterr = kPglRetInconsistentInput;
      break;
    }
  VcfToPgen_ret_MISSING_TOKENS:
    putc_unlocked('\n', stdout);
    logerrprintf("Error: Line %" PRIuPTR " of --vcf file has fewer tokens than expected.\n", line_idx);
    reterr = kPglRetMalformedInput;
    break;
  VcfToPgen_ret_MALFORMED_INPUT_2N:
    putc_unlocked('\n', stdout);
    logerrputsb();
  VcfToPgen_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  VcfToPgen_ret_MALFORMED_HEADER_LINE:
    logerrprintf("Error: Header line %" PRIuPTR " of --vcf file is malformed.\n", line_idx);
    reterr = kPglRetMalformedInput;
    break;
  VcfToPgen_ret_MALFORMED_INPUT_WWN:
    putc_unlocked('\n', stdout);
  VcfToPgen_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetMalformedInput;
    break;
  VcfToPgen_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  VcfToPgen_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  VcfToPgen_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  VcfToPgen_ret_DEGENERATE_DATA:
    reterr = kPglRetDegenerateData;
    break;
  }
 VcfToPgen_ret_1:
  CleanupSpgw(&spgw, &reterr);
  CleanupThreads(&tg);
  CleanupTextStream2("--vcf file", &vcf_txs, &reterr);
  CswriteCloseCond(&pvar_css, pvar_cswritep);
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

PglErr BcfHeaderLineIdxCheck(const char* line_iter, uint32_t header_line_idx) {
  while (1) {
    const char* tag_start = &(line_iter[1]);
    if (unlikely(memequal_sk(tag_start, "IDX="))) {
      // Although this is text, it's not supposed to be human-edited, so we
      // deliberately discourage spec-conforming behavior that's different from
      // bcftools's straightforward approach of always putting IDX= at the end.
      logerrprintfww("Error: Line %u in BCF text header block has IDX= in the center instead of the end of the line; this is not currently supported by " PROG_NAME_STR ". Contact us if you need this to work.\n", header_line_idx);
      return kPglRetNotYetSupported;
    }
    line_iter = AdvPastDelim(tag_start, '=');
    if (*line_iter != '"') {
      line_iter = strchrnul_n(line_iter, ',');
      if (*line_iter == ',') {
        continue;
      }
      return kPglRetSuccess;
    }
    // Need to worry about backslash-escaped characters.
    ++line_iter;
    while (1) {
      line_iter = strchrnul2_n(line_iter, '\\', '"');
      const char cc = *line_iter;
      if (cc == '"') {
        break;
      }
      if (unlikely((cc == '\n') || (line_iter[1] == '\n'))) {
        goto BcfHeaderLineIdxCheck_FAIL;
      }
      line_iter = &(line_iter[2]);
    }
    ++line_iter;
    const char cc = *line_iter;
    if (cc == ',') {
      continue;
    }
    // bugfix (19 Feb 2020)
    if (likely(cc == '>')) {
      return kPglRetSuccess;
    }
    break;
  }
 BcfHeaderLineIdxCheck_FAIL:
  logerrprintf("Error: Line %u in BCF text header block is malformed.\n", header_line_idx);
  return kPglRetMalformedInput;
}

// Caller's responsibility to check for overread and nonnegativity.
BoolErr ScanBcfTypedInt(const unsigned char** vrec_iterp, uint32_t* uint_ptr) {
  const unsigned char* vrec_iter = *vrec_iterp;
  const uint32_t type_descriptor_byte = *vrec_iter++;
  int32_t ii;
  if (type_descriptor_byte == 0x11) {
    ii = *R_CAST(const int8_t*, vrec_iter);
    ++vrec_iter;
  } else if (type_descriptor_byte == 0x12) {
    ii = *R_CAST(const int16_t*, vrec_iter);
    vrec_iter = &(vrec_iter[2]);
  } else if (likely(type_descriptor_byte == 0x13)) {
    ii = *R_CAST(const int32_t*, vrec_iter);
    vrec_iter = &(vrec_iter[4]);
  } else {
    return 1;
  }
  *vrec_iterp = vrec_iter;
  *uint_ptr = ii;
  return 0;
}

// value_type guaranteed to be in 0..15, but otherwise unvalidated.  value_ct
// not validated.
static inline BoolErr ScanBcfType(const unsigned char** vrec_iterp, uint32_t* value_type_ptr, uint32_t* value_ct_ptr) {
  const uint32_t type_descriptor_byte = **vrec_iterp;
  *vrec_iterp += 1;
  *value_type_ptr = type_descriptor_byte & 0xf;
  *value_ct_ptr = type_descriptor_byte >> 4;
  if (*value_ct_ptr != 15) {
    return 0;
  }
  return ScanBcfTypedInt(vrec_iterp, value_ct_ptr);
}

static inline BoolErr ScanBcfTypeAligned(const unsigned char** vrec_iterp, uint32_t* value_type_ptr, uint32_t* value_ct_ptr) {
  const uint32_t type_descriptor_byte = **vrec_iterp;
  *value_type_ptr = type_descriptor_byte & 0xf;
  *value_ct_ptr = type_descriptor_byte >> 4;
  if (*value_ct_ptr != 15) {
    *vrec_iterp += kBytesPerVec;
    return 0;
  }
#ifdef __LP64__
  const unsigned char* vrec_iter_tmp = *vrec_iterp;
  ++vrec_iter_tmp;
  *vrec_iterp += kBytesPerVec;
  // bugfix (5 Apr 2023): this was not scanning the correct value.  Previously
  // unnoticed since it would only come up for very high ploidy.
  return ScanBcfTypedInt(&vrec_iter_tmp, value_ct_ptr);
#else
  *vrec_iterp += 1;
  BoolErr ret_boolerr = ScanBcfTypedInt(vrec_iterp, value_ct_ptr);
  AlignKUcToVec(vrec_iterp);
  return ret_boolerr;
#endif
}

// Check overread within this function, since it can be long enough for integer
// overflow, etc.
BoolErr ScanBcfTypedString(const unsigned char* vrec_end, const unsigned char** vrec_iterp, const char** string_startp, uint32_t* slen_ptr) {
  const unsigned char* vrec_iter = *vrec_iterp;
  const uint32_t type_descriptor_byte = *vrec_iter++;
  if (unlikely((type_descriptor_byte & 0xf) != 7)) {
    return 1;
  }
  uint32_t slen = type_descriptor_byte >> 4;
  if (slen == 15) {
    if (unlikely(ScanBcfTypedInt(&vrec_iter, &slen) || (S_CAST(int32_t, slen) < 15))) {
      return 1;
    }
  }
  if (unlikely(S_CAST(intptr_t, vrec_end - vrec_iter) < S_CAST(intptr_t, slen))) {
    return 1;
  }
  *string_startp = R_CAST(const char*, vrec_iter);
  *vrec_iterp = &(vrec_iter[slen]);
  *slen_ptr = slen;
  return 0;
}

// Unlike the VCF case, everything in BcfImportContext is constant across the
// file.  Variation in FORMAT field presence and order is dealt with before the
// handoff to BcfGenoToPgenThread, and UINT32_MAX offset values are used to
// indicate missing fields.
typedef struct BcfImportBaseContextStruct {
  uint32_t sample_ct;
  VcfHalfCall halfcall_mode;
  uint32_t error_on_polyploid;
  // [0] = GQ, [1] = DP
  STD_ARRAY_DECL(int32_t, 2, qual_mins);
  STD_ARRAY_DECL(int32_t, 2, qual_maxs);
#ifdef USE_AVX2
  unsigned char biallelic_gt_lookup[32];
#else
  unsigned char biallelic_gt_lookup[16];
#endif
} BcfImportBaseContext;

typedef struct BcfImportContextStruct {
  BcfImportBaseContext bibc;
  uint32_t dosage_is_gp;
  uint32_t dosage_erase_halfdist;
  double import_dosage_certainty;
} BcfImportContext;

ENUM_U31_DEF_START()
  kBcfParseOk,
  kBcfParseMalformedGeneric,
  kBcfParseHalfCallError,
  kBcfParseInvalidDosage,
  kBcfParsePolyploidError,
  kBcfParseWideGt,
  kBcfParseFloatDp,
  kBcfParseNonfloatDosage
ENUM_U31_DEF_END(BcfParseErr);

// Sets genovec bits to 0b11 whenever GQ/DP doesn't pass.  (Caller should pass
// in a separate array if they need to distinguish between ordinary genovec
// missingness vs. missingness due to this filter.)
BcfParseErr BcfParseGqDpMain(const unsigned char* qual_main, uint32_t sample_ct, uint32_t qual_min, uint32_t qual_max, uint32_t qual_value_type, uintptr_t* genovec) {
  // need to support int8, int16, int32, and float.
#ifndef __LP64__
  const uint32_t word_ct_m1 = (sample_ct - 1) / kBitsPerWordD2;
  uint32_t loop_len = kBitsPerWordD2;
#endif
  if (qual_value_type == 1) {
    // int8
#ifdef __LP64__
    const uint32_t fullvec_ct = sample_ct / kBytesPerVec;
    Vec4thUint* genovec_alias = R_CAST(Vec4thUint*, genovec);
#endif
    if (qual_max >= 0x7f) {
      // Usual case: safe to skip qual_max comparisons.
      // Subtract 1 with wraparound, and then check for < -1, since we need
      // 0x80 (missing) to always pass.
      const int8_t effective_qual_min_m1 = (qual_min >= 0x80)? 0x7f : (S_CAST(int8_t, qual_min) - 1);
#ifdef __LP64__
      const VecUc* qual_alias = R_CAST(const VecUc*, qual_main);
      const VecUc all1 = vecuc_set1(0xff);
      const VecI8 min_m1_vec = veci8_set1(effective_qual_min_m1);
      for (uint32_t vidx = 0; vidx != fullvec_ct; ++vidx) {
        // don't start with VecI8, since we start with a wraparound
        // subtraction.
        const VecUc vv = vecuc_loadu(&(qual_alias[vidx]));
        const VecI8 vv_m1 = R_CAST(VecI8, vv + all1);
        VecI8 fail_vec = (min_m1_vec > vv_m1);
        fail_vec = veci8_permute0xd8_if_avx2(fail_vec);
        const VecI8 fail_vec_lo = veci8_unpacklo8(fail_vec, fail_vec);
        const VecI8 fail_vec_hi = veci8_unpackhi8(fail_vec, fail_vec);
        // todo: better ARM implementation
        const Vec4thUint fail_bits_lo = veci8_movemask(fail_vec_lo);
        const Vec4thUint fail_bits_hi = veci8_movemask(fail_vec_hi);
        genovec_alias[vidx] |= fail_bits_lo | (fail_bits_hi << kBytesPerVec);
      }
      const uint32_t remainder = sample_ct % kBytesPerVec;
      if (remainder) {
        const unsigned char* trailing_start = &(qual_main[fullvec_ct * kBytesPerVec]);
        Vec4thUint fail_bits = 0;
        for (uint32_t uii = 0; uii != remainder; ++uii) {
          const int8_t cur_qual_m1 = S_CAST(int8_t, trailing_start[uii] - 1);
          if (effective_qual_min_m1 > cur_qual_m1) {
            fail_bits |= (3 * k1LU) << (2 * uii);
          }
        }
        genovec_alias[fullvec_ct] |= fail_bits;
      }
#else
      const unsigned char* qual_iter = qual_main;
      for (uint32_t widx = 0; ; ++widx) {
        if (widx >= word_ct_m1) {
          if (widx > word_ct_m1) {
            break;
          }
          loop_len = ModNz(sample_ct, kBitsPerWordD2);
        }
        uintptr_t geno_word = genovec[widx];
        for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits) {
          const int8_t cur_qual_m1 = S_CAST(int8_t, qual_iter[sample_idx_lowbits] - 1);
          if (effective_qual_min_m1 > cur_qual_m1) {
            geno_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
          }
        }
        qual_iter = &(qual_iter[kBitsPerWordD2]);
        genovec[widx] = geno_word;
      }
#endif
    } else {
      // Two comparisons required, and we need to special-case the missing
      // value.
      // (There is currently no case where we enforce a maximum but not a
      // minimum.  Even when --vcf-max-dp is specified without --vcf-min-dp, a
      // minimum depth of zero is always enforced, and we're scanning a vector
      // of signed values.)
#ifdef __LP64__
      const VecI8* qual_alias = R_CAST(const VecI8*, qual_main);
      const VecI8 min_vec = veci8_set1(qual_min);
      const VecI8 max_vec = veci8_set1(qual_max);
      const VecI8 missing_vec = veci8_set1(0x80);
      for (uint32_t vidx = 0; vidx != fullvec_ct; ++vidx) {
        const VecI8 vv = veci8_loadu(&(qual_alias[vidx]));
        const VecI8 fail_or_missing_vec = (min_vec > vv) | (max_vec < vv);
        const VecI8 cur_missing_vec = (vv == missing_vec);
        VecI8 fail_vec = veci8_and_notfirst(cur_missing_vec, fail_or_missing_vec);
        fail_vec = veci8_permute0xd8_if_avx2(fail_vec);
        const VecI8 fail_vec_lo = veci8_unpacklo8(fail_vec, fail_vec);
        const VecI8 fail_vec_hi = veci8_unpackhi8(fail_vec, fail_vec);
        const Vec4thUint fail_bits_lo = veci8_movemask(fail_vec_lo);
        const Vec4thUint fail_bits_hi = veci8_movemask(fail_vec_hi);
        genovec_alias[vidx] |= fail_bits_lo | (fail_bits_hi << kBytesPerVec);
      }
      const uint32_t remainder = sample_ct % kBytesPerVec;
      if (remainder) {
        const int8_t* trailing_start = R_CAST(const int8_t*, &(qual_main[fullvec_ct * kBytesPerVec]));
        const int8_t qual_min_i8 = S_CAST(int8_t, qual_min);
        const int8_t qual_max_i8 = S_CAST(int8_t, qual_max);
        Vec4thUint fail_bits = 0;
        for (uint32_t uii = 0; uii != remainder; ++uii) {
          const int8_t cur_qual = trailing_start[uii];
          if (((qual_min_i8 > cur_qual) || (qual_max_i8 < cur_qual)) && (cur_qual != -128)) {
            fail_bits |= (3 * k1LU) << (2 * uii);
          }
        }
        genovec_alias[fullvec_ct] |= fail_bits;
      }
#else
      const int8_t* qual_iter = R_CAST(const int8_t*, qual_main);
      const int8_t qual_min_i8 = S_CAST(int8_t, qual_min);
      const int8_t qual_max_i8 = S_CAST(int8_t, qual_max);
      for (uint32_t widx = 0; ; ++widx) {
        if (widx >= word_ct_m1) {
          if (widx > word_ct_m1) {
            break;
          }
          loop_len = ModNz(sample_ct, kBitsPerWordD2);
        }
        uintptr_t geno_word = genovec[widx];
        for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits) {
          const int8_t cur_qual = qual_iter[sample_idx_lowbits];
          if (((qual_min_i8 > cur_qual) || (qual_max_i8 < cur_qual)) && (cur_qual != -128)) {
            geno_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
          }
        }
        qual_iter = &(qual_iter[kBitsPerWordD2]);
        genovec[widx] = geno_word;
      }
#endif
    }
  } else if (qual_value_type == 2) {
    // int16
#ifdef __LP64__
    const uint32_t fullvec_ct = sample_ct / kInt16PerVec;
    Vec8thUint* genovec_alias = R_CAST(Vec8thUint*, genovec);
#endif
    if (qual_max >= 0x7fff) {
      const int16_t effective_qual_min_m1 = (qual_min >= 0x8000)? 0x7fff : (S_CAST(int16_t, qual_min) - 1);
#ifdef __LP64__
      const VecU16* qual_alias = R_CAST(const VecU16*, qual_main);
      const VecU16 all1 = vecu16_set1(0xffff);
      const VecI16 min_m1_vec = veci16_set1(effective_qual_min_m1);
      for (uint32_t vidx = 0; vidx != fullvec_ct; ++vidx) {
        const VecU16 vv = vecu16_loadu(&(qual_alias[vidx]));
        const VecI16 vv_m1 = R_CAST(VecI16, vv + all1);
        const VecI16 fail_vec = (min_m1_vec > vv_m1);
        const Vec8thUint fail_bits = veci16_movemask(fail_vec);
        genovec_alias[vidx] |= fail_bits;
      }
      const uint32_t remainder = sample_ct % kInt16PerVec;
      if (remainder) {
        const uint16_t* trailing_start = R_CAST(const uint16_t*, &(qual_main[fullvec_ct * kInt16PerVec]));
        Vec8thUint fail_bits = 0;
        for (uint32_t uii = 0; uii != remainder; ++uii) {
          const int16_t cur_qual_m1 = S_CAST(int16_t, trailing_start[uii]);
          if (effective_qual_min_m1 > cur_qual_m1) {
            fail_bits |= 3U << (2 * uii);
          }
        }
        genovec_alias[fullvec_ct] |= fail_bits;
      }
#else
      const uint16_t* qual_iter = R_CAST(const uint16_t*, qual_main);
      for (uint32_t widx = 0; ; ++widx) {
        if (widx >= word_ct_m1) {
          if (widx > word_ct_m1) {
            break;
          }
          loop_len = ModNz(sample_ct, kBitsPerWordD2);
        }
        uintptr_t geno_word = genovec[widx];
        for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits) {
          const int16_t cur_qual_m1 = S_CAST(int16_t, qual_iter[sample_idx_lowbits] - 1);
          if (effective_qual_min_m1 > cur_qual_m1) {
            geno_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
          }
        }
        qual_iter = &(qual_iter[kBitsPerWordD2]);
        genovec[widx] = geno_word;
      }
#endif
    } else {
      // Two comparisons required, and we need to special-case the missing
      // value.
#ifdef __LP64__
      const VecI16* qual_alias = R_CAST(const VecI16*, qual_main);
      const VecI16 min_vec = veci16_set1(qual_min);
      const VecI16 max_vec = veci16_set1(qual_max);
      const VecI16 missing_vec = veci16_set1(0x8000);
      for (uint32_t vidx = 0; vidx != fullvec_ct; ++vidx) {
        const VecI16 vv = veci16_loadu(&(qual_alias[vidx]));
        const VecI16 fail_or_missing_vec = (min_vec > vv) | (max_vec < vv);
        const VecI16 cur_missing_vec = (vv == missing_vec);
        const VecI16 fail_vec = veci16_and_notfirst(cur_missing_vec, fail_or_missing_vec);
        const Vec8thUint fail_bits = veci16_movemask(fail_vec);
        genovec_alias[vidx] |= fail_bits;
      }
      const uint32_t remainder = sample_ct % kInt16PerVec;
      if (remainder) {
        const int16_t* trailing_start = R_CAST(const int16_t*, &(qual_main[fullvec_ct * kInt16PerVec]));
        const int16_t qual_min_i16 = S_CAST(int16_t, qual_min);
        const int16_t qual_max_i16 = S_CAST(int16_t, qual_max);
        Vec8thUint fail_bits = 0;
        for (uint32_t uii = 0; uii != remainder; ++uii) {
          const int16_t cur_qual = trailing_start[uii];
          if (((qual_min_i16 > cur_qual) || (qual_max_i16 < cur_qual)) && (cur_qual != -32768)) {
            fail_bits |= 3U << (2 * uii);
          }
        }
        genovec_alias[fullvec_ct] |= fail_bits;
      }
#else
      const int16_t* qual_iter = R_CAST(const int16_t*, qual_main);
      const int16_t qual_min_i16 = S_CAST(int16_t, qual_min);
      const int16_t qual_max_i16 = S_CAST(int16_t, qual_max);
      for (uint32_t widx = 0; ; ++widx) {
        if (widx >= word_ct_m1) {
          if (widx > word_ct_m1) {
            break;
          }
          loop_len = ModNz(sample_ct, kBitsPerWordD2);
        }
        uintptr_t geno_word = genovec[widx];
        for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits) {
          const int16_t cur_qual = qual_iter[sample_idx_lowbits];
          if (((qual_min_i16 > cur_qual) || (qual_max_i16 < cur_qual)) && (cur_qual != -32768)) {
            geno_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
          }
        }
        qual_iter = &(qual_iter[kBitsPerWordD2]);
        genovec[widx] = geno_word;
      }
#endif
    }
  } else if (unlikely((qual_value_type == 5) && (qual_max != 0x7fffffff))) {
    // Leave out the two-comparison float subcase for now, it's messier and
    // floating-point DP is unlikely.
    return kBcfParseFloatDp;
  } else if ((qual_value_type == 3) || ((qual_value_type == 5) && (!qual_min))) {
    // int32, or float qual_min == 0.
#ifdef __LP64__
    const uint32_t fullvec_ct = sample_ct / kInt32PerVec;
    Vec16thUint* genovec_alias = R_CAST(Vec16thUint*, genovec);
#endif
    if (qual_max == 0x7fffffff) {
      const int32_t qual_min_m1 = S_CAST(int32_t, qual_min) - 1;
#ifdef __LP64__
      const VecU32* qual_alias = R_CAST(const VecU32*, qual_main);
      const VecU32 all1 = vecu32_set1(UINT32_MAX);
      const VecI32 min_m1_vec = veci32_set1(qual_min_m1);
      for (uint32_t vidx = 0; vidx != fullvec_ct; ++vidx) {
        const VecU32 vv = vecu32_loadu(&(qual_alias[vidx]));
        const VecI32 vv_m1 = R_CAST(VecI32, vv + all1);
        const VecI32 fail_vec = (min_m1_vec > vv_m1);
#  ifdef USE_AVX2
        const Vec16thUint fail_bits = _pext_u32(veci32_movemask(fail_vec), 0x33333333);
#  else
        const __m128i vec_packed = _mm_packs_epi16(R_CAST(__m128i, fail_vec), R_CAST(__m128i, fail_vec));
        // unwanted high bits get truncated here.
        const Vec16thUint fail_bits = _mm_movemask_epi8(vec_packed);
#  endif
        genovec_alias[vidx] |= fail_bits;
      }
      const uint32_t remainder = sample_ct % kInt32PerVec;
      if (remainder) {
        const uint32_t* trailing_start = R_CAST(const uint32_t*, &(qual_main[fullvec_ct * kInt32PerVec]));
        Vec16thUint fail_bits = 0;
        for (uint32_t uii = 0; uii != remainder; ++uii) {
          const int32_t cur_qual_m1 = S_CAST(int32_t, trailing_start[uii] - 1);
          if (qual_min_m1 > cur_qual_m1) {
            fail_bits |= 3U << (2 * uii);
          }
        }
        genovec_alias[fullvec_ct] |= fail_bits;
      }
#else
      const uint32_t* qual_iter = R_CAST(const uint32_t*, qual_main);
      for (uint32_t widx = 0; ; ++widx) {
        if (widx >= word_ct_m1) {
          if (widx > word_ct_m1) {
            break;
          }
          loop_len = ModNz(sample_ct, kBitsPerWordD2);
        }
        uintptr_t geno_word = genovec[widx];
        for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits) {
          const int32_t cur_qual_m1 = S_CAST(int32_t, qual_iter[sample_idx_lowbits] - 1);
          if (qual_min_m1 > cur_qual_m1) {
            geno_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
          }
        }
        qual_iter = &(qual_iter[kBitsPerWordD2]);
        genovec[widx] = geno_word;
      }
#endif
    } else {
      // Two comparisons required, and we need to special-case the missing
      // value.
#ifdef __LP64__
      const VecI32* qual_alias = R_CAST(const VecI32*, qual_main);
      const VecI32 min_vec = veci32_set1(qual_min);
      const VecI32 max_vec = veci32_set1(qual_max);
      const VecI32 missing_vec = veci32_set1(0x80000000);
      for (uint32_t vidx = 0; vidx != fullvec_ct; ++vidx) {
        const VecI32 vv = veci32_loadu(&(qual_alias[vidx]));
        const VecI32 fail_or_missing_vec = (min_vec > vv) | (max_vec < vv);
        const VecI32 cur_missing_vec = (vv == missing_vec);
        const VecI32 fail_vec = veci32_and_notfirst(cur_missing_vec, fail_or_missing_vec);
#  ifdef USE_AVX2
        const Vec16thUint fail_bits = _pext_u32(veci32_movemask(fail_vec), 0x33333333);
#  else
        const __m128i vec_packed = _mm_packs_epi16(R_CAST(__m128i, fail_vec), R_CAST(__m128i, fail_vec));
        const Vec16thUint fail_bits = _mm_movemask_epi8(vec_packed);
#  endif
        genovec_alias[vidx] |= fail_bits;
      }
      const uint32_t remainder = sample_ct % kInt32PerVec;
      if (remainder) {
        const int32_t* trailing_start = R_CAST(const int32_t*, &(qual_main[fullvec_ct * kInt32PerVec]));
        const int32_t qual_min_i32 = S_CAST(int32_t, qual_min);
        const int32_t qual_max_i32 = S_CAST(int32_t, qual_max);
        Vec16thUint fail_bits = 0;
        for (uint32_t uii = 0; uii != remainder; ++uii) {
          const int32_t cur_qual = trailing_start[uii];
          // see https://stackoverflow.com/questions/9941261/warning-this-decimal-constant-is-unsigned-only-in-iso-c90
          if (((qual_min_i32 > cur_qual) || (qual_max_i32 < cur_qual)) && (cur_qual != (-2147483647 - 1))) {
            fail_bits |= 3U << (2 * uii);
          }
        }
        genovec_alias[fullvec_ct] |= fail_bits;
      }
#else
      const int32_t* qual_iter = R_CAST(const int32_t*, qual_main);
      const int32_t qual_min_i32 = S_CAST(int32_t, qual_min);
      const int32_t qual_max_i32 = S_CAST(int32_t, qual_max);
      for (uint32_t widx = 0; ; ++widx) {
        if (widx >= word_ct_m1) {
          if (widx > word_ct_m1) {
            break;
          }
          loop_len = ModNz(sample_ct, kBitsPerWordD2);
        }
        uintptr_t geno_word = genovec[widx];
        for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits) {
          const int32_t cur_qual = qual_iter[sample_idx_lowbits];
          if (((qual_min_i32 > cur_qual) || (qual_max_i32 < cur_qual)) && (cur_qual != (-2147483647 - 1))) {
            geno_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
          }
        }
        qual_iter = &(qual_iter[kBitsPerWordD2]);
        genovec[widx] = geno_word;
      }
#endif
    }
  } else if (likely(qual_value_type == 5)) {
    // float, qual_min > 0.
    // 0x80000000 does not pass in this case.
    if (qual_value_type == 5) {
      // Still more efficient to perform the comparison in integer-space, since
      // we treat all NaN values (0x7f800000...0x7fffffff) as passing.

      // Convert to the appropriate floating-point bit pattern, rounding up.
      const uint32_t orig_qual_min = qual_min;
      const float fxx = S_CAST(float, qual_min);
      memcpy(&qual_min, &fxx, 4);
      if (orig_qual_min > 0x1000000) {
        // We need to add 1 to qual_min iff high_bit_idx > low_bit_idx + 23.
        const uint32_t high_bit_idx = bsru32(orig_qual_min);
        const uint32_t low_bit_idx = ctzu32(orig_qual_min);
        qual_min += (high_bit_idx > (low_bit_idx + 23));
      }
    }
#ifdef __LP64__
    const uint32_t fullvec_ct = sample_ct / kInt32PerVec;
    Vec16thUint* genovec_alias = R_CAST(Vec16thUint*, genovec);
    const VecI32* qual_alias = R_CAST(const VecI32*, qual_main);
    const VecI32 min_vec = veci32_set1(qual_min);
    for (uint32_t vidx = 0; vidx != fullvec_ct; ++vidx) {
      const VecI32 vv = veci32_loadu(&(qual_alias[vidx]));
      const VecI32 fail_vec = (min_vec > vv);
#  ifdef USE_AVX2
      const Vec16thUint fail_bits = _pext_u32(veci32_movemask(fail_vec), 0x33333333);
#  else
      const __m128i vec_packed = _mm_packs_epi16(R_CAST(__m128i, fail_vec), R_CAST(__m128i, fail_vec));
      // unwanted high bits get truncated here.
      const Vec16thUint fail_bits = _mm_movemask_epi8(vec_packed);
#  endif
      genovec_alias[vidx] |= fail_bits;
    }
    const uint32_t remainder = sample_ct % kInt32PerVec;
    if (remainder) {
      const int32_t* trailing_start = R_CAST(const int32_t*, &(qual_main[fullvec_ct * kInt32PerVec]));
      const int32_t qual_min_i32 = S_CAST(int32_t, qual_min);
      Vec16thUint fail_bits = 0;
      for (uint32_t uii = 0; uii != remainder; ++uii) {
        if (qual_min_i32 > trailing_start[uii]) {
          fail_bits |= 3U << (2 * uii);
        }
      }
      genovec_alias[fullvec_ct] |= fail_bits;
    }
#else
    const int32_t* qual_iter = R_CAST(const int32_t*, qual_main);
    const int32_t qual_min_i32 = S_CAST(int32_t, qual_min);
    for (uint32_t widx = 0; ; ++widx) {
      if (widx >= word_ct_m1) {
        if (widx > word_ct_m1) {
          break;
        }
        loop_len = ModNz(sample_ct, kBitsPerWordD2);
      }
      uintptr_t geno_word = genovec[widx];
      for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits) {
        if (qual_min_i32 > qual_iter[sample_idx_lowbits]) {
          geno_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
        }
      }
      qual_iter = &(qual_iter[kBitsPerWordD2]);
      genovec[widx] = geno_word;
    }
#endif
  } else {
    return kBcfParseMalformedGeneric;
  }
  return kBcfParseOk;
}

BcfParseErr BcfParseGqDpUnaligned(const unsigned char* qual_type_start, uint32_t sample_ct, uint32_t qual_min, uint32_t qual_max, uintptr_t* genovec) {
  const unsigned char* qual_iter = qual_type_start;
  uint32_t qual_value_type;
  uint32_t qual_value_ct;
  if (unlikely(ScanBcfType(&qual_iter, &qual_value_type, &qual_value_ct) || (qual_value_ct > 1))) {
    return kBcfParseMalformedGeneric;
  }
  if (!qual_value_ct) {
    return kBcfParseOk;
  }
  return BcfParseGqDpMain(qual_iter, sample_ct, qual_min, qual_max, qual_value_type, genovec);
}

BcfParseErr BcfParseGqDpAligned(const unsigned char* qual_type_start, uint32_t sample_ct, uint32_t qual_min, uint32_t qual_max, uintptr_t* genovec) {
  const unsigned char* qual_iter = qual_type_start;
  uint32_t qual_value_type;
  uint32_t qual_value_ct;
  if (unlikely(ScanBcfTypeAligned(&qual_iter, &qual_value_type, &qual_value_ct) || (qual_value_ct > 1))) {
    return kBcfParseMalformedGeneric;
  }
  if (!qual_value_ct) {
    return kBcfParseOk;
  }
  return BcfParseGqDpMain(qual_iter, sample_ct, qual_min, qual_max, qual_value_type, genovec);
}

static_assert(kPglMaxAlleleCt == 255, "BcfScanGt() needs to be updated.");
BcfParseErr BcfScanGt(const BcfImportContext* bicp, const unsigned char* gt_start, const unsigned char** qual_starts, uint32_t* phase_or_dosage_found_ptr, uintptr_t* invfound_nypbuf) {
  // Just check for a phased het.
  const unsigned char* gt_main = gt_start;
  // Note that generic validation was already performed on these vectors.  We
  // just need to verify that e.g. GT values are int8/int16, and GQ/DP values
  // aren't characters.
  uint32_t gt_value_type;
  uint32_t gt_value_ct;
  if (unlikely(ScanBcfType(&gt_main, &gt_value_type, &gt_value_ct) || (gt_value_type > 2))) {
    return (gt_value_type == 3)? kBcfParseWideGt : kBcfParseMalformedGeneric;
  }
  if (gt_value_ct < 2) {
    // ploidy 0 or 1, phased genotypes are impossible
    return kBcfParseOk;
  }
  if (unlikely((gt_value_ct > 2) && bicp->bibc.error_on_polyploid)) {
    return kBcfParsePolyploidError;
  }
  const uint32_t sample_ct = bicp->bibc.sample_ct;
#ifdef __LP64__
  // Our search is over once we see a phased heterozygous call.  Usually, if
  // one exists at all, a filter-passing call will be present within the first
  // few variants; the important case to optimize is no-phased-genotypes,
  // all-diploid.
  // So, in the 64-bit diploid case, we perform a pre-scan for a set bottom bit
  // in all of the second-genotype-values in each pair.  If none exist, there
  // can't be any phased (or haploid/0-ploid) genotypes, and we can move on to
  // the next variant.
  // (Initially tried to make this an exhaustive scan, but then realized
  // halfcalls and SIMD don't get along that well.)
  if (gt_value_ct == 2) {
    if (gt_value_type == 1) {
      // int8
      if (sample_ct >= kInt16PerVec) {
        const uint32_t vec_ct_m1 = (sample_ct - 1) / kInt16PerVec;
        const VecU16* gt_alias = R_CAST(const VecU16*, gt_main);
        VecU16 found = vecu16_setzero();
        for (uint32_t vidx = 0; vidx != vec_ct_m1; ++vidx) {
          const VecU16 vv = vecu16_loadu(&(gt_alias[vidx]));
          found = found | vv;
        }
        const VecU16 last_vec = vecu16_loadu(&(gt_main[(sample_ct - kInt16PerVec) * 2]));
        found = found | last_vec;
        if (!vecu16_movemask(vecu16_srli(found, 1))) {
          // No phased (or haploid) genotypes at all.
          return kBcfParseOk;
        }
      }
    } else {
      // int16
      if (sample_ct >= kInt32PerVec) {
        const uint32_t vec_ct_m1 = (sample_ct - 1) / kInt32PerVec;
        const VecU32* gt_alias = R_CAST(const VecU32*, gt_main);
        VecU32 found = vecu32_setzero();
        for (uint32_t vidx = 0; vidx != vec_ct_m1; ++vidx) {
          const VecU32 vv = vecu32_loadu(&(gt_alias[vidx]));
          found = found | vv;
        }
        const VecU32 last_vec = vecu32_loadu(&(gt_main[(sample_ct - kInt32PerVec) * 4]));
        found = found | last_vec;
        if (!(vecu32_movemask(vecu32_srli(found, 1)) & 0x22222222U)) {
          return kBcfParseOk;
        }
      }
    }
  }
#endif
  const VcfHalfCall halfcall_mode = bicp->bibc.halfcall_mode;
  if ((!qual_starts[0]) && (!qual_starts[1])) {
    // No gq/dp check, so we can exit as soon as we find a match.
    if (gt_value_type == 1) {
      // int8
      if (gt_value_ct == 2) {
        // Usual case.
        // Low bit of each second byte is set iff the genotype is either phased
        // or non-diploid (0x81 "END_OF_VECTOR").
        const unsigned char* second_byte_stop = &(gt_main[sample_ct * 2 + 1]);
        for (const unsigned char* second_byte_iter = &(gt_main[1]); second_byte_iter != second_byte_stop; second_byte_iter = &(second_byte_iter[2])) {
          const uint32_t second_byte = *second_byte_iter;
          if ((second_byte & 0x81) == 1) {
            // phased
            const uint32_t first_allele_idx_p1 = second_byte_iter[-1] >> 1;
            const uint32_t second_allele_idx_p1 = second_byte >> 1;
            if (first_allele_idx_p1 != second_allele_idx_p1) {
              // phased het or halfcall
              // (or malformed, false positive ok there)
              if (first_allele_idx_p1 && second_allele_idx_p1) {
                *phase_or_dosage_found_ptr = 1;
                return kBcfParseOk;
              }
              if (halfcall_mode == kVcfHalfCallReference) {
                if (first_allele_idx_p1 + second_allele_idx_p1 != 1) {
                  *phase_or_dosage_found_ptr = 1;
                  return kBcfParseOk;
                }
              } else if (unlikely(halfcall_mode == kVcfHalfCallError)) {
                return kBcfParseHalfCallError;
              }
            }
          }
        }
      } else {
        // Ploidy > 2 treated as missing.
        // possible todo: Add SIMD code for the triploid case (advance 15 or 30
        // bytes on each iteration, etc.), since it does realistically come up
        // in large datasets, and it only takes one triploid sample to force
        // this code path to be taken for all samples.
        // (Tetraploidy is much rarer than triploidy.)
        const unsigned char* second_byte_stop = &(gt_main[sample_ct * gt_value_type + 1]);
        for (const unsigned char* second_byte_iter = &(gt_main[1]); second_byte_iter != second_byte_stop; second_byte_iter = &(second_byte_iter[gt_value_type])) {
          uint16_t second_and_third_bytes;
          memcpy(&second_and_third_bytes, second_byte_iter, 2);
          if ((second_and_third_bytes & 0x8181) == 0x8101) {
            // phased diploid
            const uint32_t first_allele_idx_p1 = second_byte_iter[-1] >> 1;
            const uint32_t second_allele_idx_p1 = (second_and_third_bytes & 0xff) >> 1;
            if (first_allele_idx_p1 != second_allele_idx_p1) {
              // phased het or halfcall
              if (first_allele_idx_p1 && second_allele_idx_p1) {
                *phase_or_dosage_found_ptr = 1;
                return kBcfParseOk;
              }
              if (halfcall_mode == kVcfHalfCallReference) {
                if (first_allele_idx_p1 + second_allele_idx_p1 != 1) {
                  *phase_or_dosage_found_ptr = 1;
                  return kBcfParseOk;
                }
              } else if (unlikely(halfcall_mode == kVcfHalfCallError)) {
                return kBcfParseHalfCallError;
              }
            }
          }
        }
      }
    } else {
      // int16
      if (gt_value_ct == 2) {
        const uint16_t* second_u16_stop = R_CAST(const uint16_t*, &(gt_main[sample_ct * 4 + 2]));
        for (const uint16_t* second_u16_iter = R_CAST(const uint16_t*, &(gt_main[2])); second_u16_iter != second_u16_stop; second_u16_iter = &(second_u16_iter[2])) {
          const uint32_t second_u16 = *second_u16_iter;
          if ((second_u16 & 0x8001) == 1) {
            // phased
            const uint32_t first_allele_idx_p1 = second_u16_iter[-1] >> 1;
            const uint32_t second_allele_idx_p1 = second_u16 >> 1;
            if (first_allele_idx_p1 != second_allele_idx_p1) {
              // phased het or halfcall
              if (first_allele_idx_p1 && second_allele_idx_p1) {
                *phase_or_dosage_found_ptr = 1;
                return kBcfParseOk;
              }
              if (halfcall_mode == kVcfHalfCallReference) {
                if (first_allele_idx_p1 + second_allele_idx_p1 != 1) {
                  *phase_or_dosage_found_ptr = 1;
                  return kBcfParseOk;
                }
              } else if (unlikely(halfcall_mode == kVcfHalfCallError)) {
                return kBcfParseHalfCallError;
              }
            }
          }
        }
      } else {
        // ploidy > 2
        const uint16_t* second_u16_stop = R_CAST(const uint16_t*, &(gt_main[sample_ct * 2 * gt_value_type + 2]));
        for (const uint16_t* second_u16_iter = R_CAST(const uint16_t*, &(gt_main[2])); second_u16_iter != second_u16_stop; second_u16_iter = &(second_u16_iter[gt_value_type])) {
          uint32_t second_and_third_u16s;
          memcpy(&second_and_third_u16s, second_u16_iter, 4);
          if ((second_and_third_u16s & 0x80018001U) == 0x80010001U) {
            // phased diploid
            const uint32_t first_allele_idx_p1 = second_u16_iter[-1] >> 1;
            const uint32_t second_allele_idx_p1 = (second_and_third_u16s & 0xffff) >> 1;
            if (first_allele_idx_p1 != second_allele_idx_p1) {
              // phased het or halfcall
              if (first_allele_idx_p1 && second_allele_idx_p1) {
                *phase_or_dosage_found_ptr = 1;
                return kBcfParseOk;
              }
              if (halfcall_mode == kVcfHalfCallReference) {
                if (first_allele_idx_p1 + second_allele_idx_p1 != 1) {
                  *phase_or_dosage_found_ptr = 1;
                  return kBcfParseOk;
                }
              } else if (unlikely(halfcall_mode == kVcfHalfCallError)) {
                return kBcfParseHalfCallError;
              }
            }
          }
        }
      }
    }
    return kBcfParseOk;
  }
  // We clear a bit when we find a phase/dosage we're keeping; then
  // BcfParseGqDpUnaligned() sets it iff the GQ or DP filter fails for that
  // sample; if any clear bits remain at the end, we're done with our scan.
  SetAllBits(sample_ct * 2, invfound_nypbuf);
  if (gt_value_type == 1) {
    if (gt_value_ct == 2) {
      const unsigned char* second_byte_start = &(gt_main[1]);
      for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        const uint32_t second_byte = second_byte_start[sample_idx * 2];
        if ((second_byte & 0x81) == 1) {
          // phased
          const uint32_t first_allele_idx_p1 = second_byte_start[sample_idx * 2 - 1] >> 1;
          const uint32_t second_allele_idx_p1 = second_byte >> 1;
          if (first_allele_idx_p1 != second_allele_idx_p1) {
            // phased het or halfcall
            if (first_allele_idx_p1 && second_allele_idx_p1) {
              ClearBit(sample_idx * 2, invfound_nypbuf);
            } else if (halfcall_mode == kVcfHalfCallReference) {
              if (first_allele_idx_p1 + second_allele_idx_p1 != 1) {
                ClearBit(sample_idx * 2, invfound_nypbuf);
              }
            } else if (unlikely(halfcall_mode == kVcfHalfCallError)) {
              return kBcfParseHalfCallError;
            }
          }
        }
      }
    } else {
      // ploidy > 2
      const unsigned char* second_byte_iter = &(gt_main[1]);
      for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx, second_byte_iter = &(second_byte_iter[gt_value_type])) {
        uint16_t second_and_third_bytes;
        memcpy(&second_and_third_bytes, second_byte_iter, 2);
        if ((second_and_third_bytes & 0x8181) == 0x8101) {
          // phased diploid
          const uint32_t first_allele_idx_p1 = second_byte_iter[-1] >> 1;
          const uint32_t second_allele_idx_p1 = (second_and_third_bytes & 0xff) >> 1;
          if (first_allele_idx_p1 != second_allele_idx_p1) {
            // phased het or halfcall
            if (first_allele_idx_p1 && second_allele_idx_p1) {
              ClearBit(sample_idx * 2, invfound_nypbuf);
            } else if (halfcall_mode == kVcfHalfCallReference) {
              if (first_allele_idx_p1 + second_allele_idx_p1 != 1) {
                ClearBit(sample_idx * 2, invfound_nypbuf);
              }
            } else if (unlikely(halfcall_mode == kVcfHalfCallError)) {
              return kBcfParseHalfCallError;
            }
          }
        }
      }
    }
  } else {
    // int16
    if (gt_value_ct == 2) {
      const uint16_t* second_u16_start = R_CAST(const uint16_t*, &(gt_main[2]));
      for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        const uint32_t second_u16 = second_u16_start[sample_idx * 2];
        if ((second_u16 & 0x8001) == 1) {
          // phased
          const uint32_t first_allele_idx_p1 = second_u16_start[sample_idx * 2 - 1] >> 1;
          const uint32_t second_allele_idx_p1 = second_u16 >> 1;
          if (first_allele_idx_p1 != second_allele_idx_p1) {
            // phased het or halfcall
            if (first_allele_idx_p1 && second_allele_idx_p1) {
              ClearBit(sample_idx * 2, invfound_nypbuf);
            } else if (halfcall_mode == kVcfHalfCallReference) {
              if (first_allele_idx_p1 + second_allele_idx_p1 != 1) {
                ClearBit(sample_idx * 2, invfound_nypbuf);
              }
            } else if (unlikely(halfcall_mode == kVcfHalfCallError)) {
              return kBcfParseHalfCallError;
            }
          }
        }
      }
    } else {
      // ploidy > 2
      const uint16_t* second_u16_iter = R_CAST(const uint16_t*, &(gt_main[2]));
      for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx, second_u16_iter = &(second_u16_iter[gt_value_type])) {
        uint32_t second_and_third_u16s;
        memcpy(&second_and_third_u16s, second_u16_iter, 4);
        if ((second_and_third_u16s & 0x80018001U) == 0x80010001U) {
          // phased diploid
          const uint32_t first_allele_idx_p1 = second_u16_iter[-1] >> 1;
          const uint32_t second_allele_idx_p1 = (second_and_third_u16s & 0xffff) >> 1;
          if (first_allele_idx_p1 != second_allele_idx_p1) {
            // phased het or halfcall
            if (first_allele_idx_p1 && second_allele_idx_p1) {
              ClearBit(sample_idx * 2, invfound_nypbuf);
            } else if (halfcall_mode == kVcfHalfCallReference) {
              if (first_allele_idx_p1 + second_allele_idx_p1 != 1) {
                ClearBit(sample_idx * 2, invfound_nypbuf);
              }
            } else if (unlikely(halfcall_mode == kVcfHalfCallError)) {
              return kBcfParseHalfCallError;
            }
          }
        }
      }
    }
  }
  if (AllBitsAreOne(invfound_nypbuf, sample_ct * 2)) {
    return kBcfParseOk;
  }
  for (uint32_t qual_idx = 0; qual_idx != 2; ++qual_idx) {
    if (!qual_starts[qual_idx]) {
      continue;
    }
    BcfParseErr bcf_parse_err = BcfParseGqDpUnaligned(qual_starts[qual_idx], sample_ct, bicp->bibc.qual_mins[qual_idx], bicp->bibc.qual_maxs[qual_idx], invfound_nypbuf);
    if (bcf_parse_err != kBcfParseOk) {
      return bcf_parse_err;
    }
    if (AllBitsAreOne(invfound_nypbuf, sample_ct * 2)) {
      return kBcfParseOk;
    }
  }
  *phase_or_dosage_found_ptr = 1;
  return kBcfParseOk;
}

// Dosage parsing is messy enough that we stick as closely as possible to the
// VCF logic, instead of trying to vectorize.
BoolErr ParseBcfBiallelicGp(const float* cur_gp_start, uint32_t is_haploid, double import_dosage_certainty, DosageParseResult* dpr_ptr, double* alt_dosage_ptr) {
  // See ParseVcfBiallelicGp().
  // P(0/0), P(0/1), P(1/1), etc.
  // assumes dpr initialized to kDosageParseOk
  // assumes *gp_iter is not missing
  // returns 1 if missing OR parsing error.  error: dpr still kDosageParseOk
  const float prob_0altf = cur_gp_start[0];
  if (unlikely((prob_0altf < S_CAST(float, 0.0)) || (prob_0altf > S_CAST(float, 1.0)))) {
    return 1;
  }
  const float prob_1altf = cur_gp_start[1];
  // second predicate written to be true on NaN
  if (unlikely((prob_1altf < S_CAST(float, 0.0)) || (!(prob_1altf <= S_CAST(float, 1.0))))) {
    return 1;
  }
  const double prob_0alt = S_CAST(double, prob_0altf);
  const double prob_1alt = S_CAST(double, prob_1altf);
  if (is_haploid) {
    const double denom = prob_0alt + prob_1alt;
    if (denom <= 2 * import_dosage_certainty) {
      if ((prob_0alt <= import_dosage_certainty) && (prob_1alt <= import_dosage_certainty)) {
        *dpr_ptr = kDosageParseForceMissing;
        return 1;
      }
    }
    *alt_dosage_ptr = 2 * prob_1alt / denom;
    return 0;
  }
  const float prob_2altf = cur_gp_start[2];
  if (unlikely((prob_2altf < S_CAST(float, 0.0)) || (!(prob_2altf <= S_CAST(float, 1.0))))) {
    return 1;
  }
  const double prob_2alt = S_CAST(double, prob_2altf);
  const double denom = prob_0alt + prob_1alt + prob_2alt;
  if (denom <= 3 * import_dosage_certainty) {
    if ((prob_0alt <= import_dosage_certainty) && (prob_1alt <= import_dosage_certainty) && (prob_2alt <= import_dosage_certainty)) {
      // force-missing
      // ok to use <= since we multiplied by (1 - epsilon)
      // during command-line parsing.  this lets us avoid
      // special-casing denom=0.
      *dpr_ptr = kDosageParseForceMissing;
      return 1;  // not really an error
    }
  }
  *alt_dosage_ptr = (prob_1alt + 2 * prob_2alt) / denom;
  return 0;
}

BoolErr ParseBcfBiallelicDosage(const float* cur_dosage_start, uint32_t is_haploid_or_0ploid, uint32_t dosage_is_gp, double import_dosage_certainty, DosageParseResult* dpr_ptr, uint32_t* dosage_int_ptr) {
  // See ParseVcfBiallelicDosage().
  // assumes dpr initialized to kDosageParseOk
  // returns 1 if missing OR parsing error.  error: dpr still kDosageParseOk.
  int32_t first_bits;
  memcpy(&first_bits, cur_dosage_start, 4);
  if (first_bits > 0x7f800000) {
    // ploidy-0, regular missing call, or NaN.
    *dpr_ptr = kDosageParseMissing;
    return 1;
  }
  double alt_dosage;
  if (dosage_is_gp) {
    if (ParseBcfBiallelicGp(cur_dosage_start, is_haploid_or_0ploid, import_dosage_certainty, dpr_ptr, &alt_dosage)) {
      return 1;
    }
  } else {
    alt_dosage = S_CAST(double, *cur_dosage_start);
    if (unlikely(alt_dosage < 0.0)) {
      return 1;
    }
    if (is_haploid_or_0ploid) {
      // possible todo: allow this to be suppressed (maybe upstream of this
      // function); 1000 Genomes phase 1 haploid dosages are still on 0..2
      // scale
      // right now the best approach for importing those files is commenting
      // out this line and recompiling...
      if (import_dosage_certainty != 0.0) {
        // quasi-bugfix (19 Feb 2019): dosage=DS import should respect
        // --import-dosage-certainty
        if (((1.0 - alt_dosage) <= import_dosage_certainty) && (alt_dosage <= import_dosage_certainty)) {
          *dpr_ptr = kDosageParseForceMissing;
          return 1;
        }
      }
      alt_dosage *= 2;
    } else {
      if (import_dosage_certainty != 0.0) {
        const double dist_from_1 = fabs(1.0 - alt_dosage);
        if ((1.0 - dist_from_1 <= import_dosage_certainty) && (dist_from_1 <= import_dosage_certainty)) {
          *dpr_ptr = kDosageParseForceMissing;
          return 1;
        }
      }
    }
    if (unlikely(alt_dosage > 2.0)) {
      return 1;
    }
  }
  *dosage_int_ptr = S_CAST(int32_t, alt_dosage * kDosageMid + 0.5);
  return 0;
}

BoolErr ParseBcfBiallelicHds(const float* dosage_main, const unsigned char* hds_main, uint32_t dosage_value_ct, uint32_t hds_value_ct, uint32_t sample_idx, uint32_t is_haploid_or_0ploid, uint32_t dosage_is_gp, double import_dosage_certainty, DosageParseResult* dpr_ptr, uint32_t* dosage_int_ptr, int32_t* cur_dphase_delta_ptr, uint32_t* hds_valid_ptr) {
  // See ParseVcfBiallelicHds().
  // assumes dpr initialized to kDosageParseOk
  // assumes cur_dphase_delta initialized to 0
  // assumes hds_valid initialized to 0
  // assumes dosage_main != nullptr and/or hds_main != nullptr
  // returns 1 if missing OR parsing error.  error: dpr still kDosageParseOk.
  if (hds_main) {
    const unsigned char* cur_hds_start = &(hds_main[sample_idx * hds_value_ct * sizeof(float)]);
    // search for HDS first, then DS
    int32_t first_bits;
    memcpy(&first_bits, cur_hds_start, 4);
    if (first_bits <= 0x7f800000) {
      float fxx;
      CopyFromUnalignedF(&fxx, cur_hds_start);
      double dosage1 = S_CAST(double, fxx);
      if (unlikely((dosage1 < 0.0) || (dosage1 > 1.0))) {
        return 1;
      }
      // if hds_valid and (cur_dphase_delta == 0), caller should override
      // hardcall-phase
      *hds_valid_ptr = 1;
      int32_t second_bits;
      CopyFromUnalignedOffsetI32(&second_bits, cur_hds_start, 1);
      if (second_bits > 0x7f800000) {
        // haploid ok, half-call not ok
        // 0x7f800002 == END_OF_VECTOR
        if (unlikely(second_bits != 0x7f800002)) {
          return 1;
        }
        if (import_dosage_certainty != 0.0) {
          if ((1.0 - dosage1 <= import_dosage_certainty) && (dosage1 <= import_dosage_certainty)) {
            *dpr_ptr = kDosageParseForceMissing;
            return 1;
          }
        }
        *dosage_int_ptr = S_CAST(int32_t, dosage1 * kDosageMax + 0.5);
        return 0;
      }
      CopyFromUnalignedOffsetF(&fxx, cur_hds_start, 1);
      double dosage2 = S_CAST(double, fxx);
      if (unlikely((dosage2 < 0.0) || (dosage2 > 1.0))) {
        return 1;
      }
      const double dosage_sum = dosage1 + dosage2;
      if (import_dosage_certainty != 0.0) {
        // Assume maximal het probability.
        const double dist_from_1 = fabs(1.0 - dosage_sum);
        if ((1.0 - dist_from_1 <= import_dosage_certainty) && (dist_from_1 <= import_dosage_certainty)) {
          *dpr_ptr = kDosageParseForceMissing;
          return 1;
        }
      }

      // force this to be nonnegative, since static_cast<int32_t> rounds
      // negative numbers toward zero
      const double dosage_diffp1 = 1.0 + dosage1 - dosage2;

      *dosage_int_ptr = S_CAST(int32_t, dosage_sum * kDosageMid + 0.5);
      *cur_dphase_delta_ptr = S_CAST(int32_t, dosage_diffp1 * kDosageMid + 0.5) - kDosageMid;
      return 0;
    }
    if (!dosage_main) {
      *dpr_ptr = kDosageParseMissing;
      return 1;
    }
  }
  const float* cur_dosage_start = &(dosage_main[sample_idx * dosage_value_ct]);
  return ParseBcfBiallelicDosage(cur_dosage_start, is_haploid_or_0ploid, dosage_is_gp, import_dosage_certainty, dpr_ptr, dosage_int_ptr);
}

uint32_t BcfGtIsPhasedHet(uint32_t gt_first, uint32_t gt_second, uint32_t gt_high_bit, VcfHalfCall halfcall_mode) {
  if (!(gt_second & 1)) {
    return 0;
  }
  const uint32_t first_allele_idx_p1 = (gt_first & (~gt_high_bit)) >> 1;
  const uint32_t second_allele_idx_p1 = (gt_second & (~gt_high_bit)) >> 1;
  if (first_allele_idx_p1 == second_allele_idx_p1) {
    return 0;
  }
  return (first_allele_idx_p1 && second_allele_idx_p1) || ((halfcall_mode == kVcfHalfCallReference) && (first_allele_idx_p1 + second_allele_idx_p1 != 1));
}

BcfParseErr BcfScanBiallelicHds(const BcfImportContext* bicp, const unsigned char* gt_start, const unsigned char** qual_starts, const unsigned char* dosage_start, const unsigned char* hds_start, uint32_t* phase_or_dosage_found_ptr, uintptr_t* __restrict invfound_nypbuf) {
  // See VcfScanBiallelicHdsLine().
  // DS+HDS
  // Only need to find phase *or* dosage.  We can expect this to happen
  // quickly (when we don't, it's essentially user error), so (unlike scanning
  // GT for a phased call) there's little point in spending much effort on
  // optimizing this.
  const unsigned char* hds_main = nullptr;
  uint32_t hds_value_ct = 0;
  if (hds_start) {
    const unsigned char* hds_main_raw = hds_start;
    uint32_t hds_value_type;
    if (unlikely(ScanBcfType(&hds_main_raw, &hds_value_type, &hds_value_ct))) {
      return kBcfParseMalformedGeneric;
    }
    if (unlikely(hds_value_type != 5)) {
      return kBcfParseNonfloatDosage;
    }
    if (hds_value_ct) {
      if (unlikely((hds_value_ct > 2) && bicp->bibc.error_on_polyploid)) {
        return kBcfParsePolyploidError;
      }
      hds_main = hds_main_raw;
    }
  }
  const float* dosage_main = nullptr;
  uint32_t dosage_value_ct = 0;
  if (dosage_start) {
    const unsigned char* dosage_main_raw = dosage_start;
    uint32_t dosage_value_type;
    if (unlikely(ScanBcfType(&dosage_main_raw, &dosage_value_type, &dosage_value_ct))) {
      return kBcfParseMalformedGeneric;
    }
    if (unlikely(dosage_value_type != 5)) {
      return kBcfParseNonfloatDosage;
    }
    if (dosage_value_ct) {
      dosage_main = R_CAST(const float*, dosage_main_raw);
    }
  }
  if ((!hds_value_ct) && (!dosage_value_ct)) {
    if (!gt_start) {
      // everything missing
      return kBcfParseOk;
    }
    return BcfScanGt(bicp, gt_start, qual_starts, phase_or_dosage_found_ptr, invfound_nypbuf);
  }
  const uint32_t sample_ct = bicp->bibc.sample_ct;
  ZeroWArr(NypCtToWordCt(sample_ct), invfound_nypbuf);
  for (uint32_t qual_idx = 0; qual_idx != 2; ++qual_idx) {
    if (!qual_starts[qual_idx]) {
      continue;
    }
    BcfParseErr bcf_parse_err = BcfParseGqDpUnaligned(qual_starts[qual_idx], sample_ct, bicp->bibc.qual_mins[qual_idx], bicp->bibc.qual_maxs[qual_idx], invfound_nypbuf);
    if (bcf_parse_err != kBcfParseOk) {
      return bcf_parse_err;
    }
    if (AllBitsAreOne(invfound_nypbuf, sample_ct * 2)) {
      return kBcfParseOk;
    }
  }
  const unsigned char* gt_main = gt_start;
  uint32_t gt_value_type = 0;
  uint32_t gt_value_ct = 0;
  if (gt_main) {
    if (unlikely(ScanBcfType(&gt_main, &gt_value_type, &gt_value_ct) || (gt_value_type > 3))) {
      return kBcfParseMalformedGeneric;
    }
    if (unlikely(gt_value_type > 1)) {
      return kBcfParseWideGt;
    }
    if (unlikely((gt_value_ct > 2) && bicp->bibc.error_on_polyploid)) {
      return kBcfParsePolyploidError;
    }
  }
  const uint32_t dosage_is_gp = bicp->dosage_is_gp;
  const uint32_t dosage_erase_halfdist = bicp->dosage_erase_halfdist;
  const double import_dosage_certainty = bicp->import_dosage_certainty;
  const VcfHalfCall halfcall_mode = bicp->bibc.halfcall_mode;
  uint32_t gt_first = 0;
  uint32_t gt_second = 0;
  uint32_t is_haploid_or_0ploid = (gt_value_ct == 1);
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    if (IsSet(invfound_nypbuf, sample_idx * 2)) {
      continue;
    }
    if (gt_value_ct > 1) {
      // no need to look at gt_first in haploid case.
      // only int8 supported, since variant is biallelic
      gt_first = gt_main[sample_idx * gt_value_ct];
      gt_second = gt_main[sample_idx * gt_value_ct + 1];
      is_haploid_or_0ploid = gt_second >> 7;
    }
    DosageParseResult dpr = kDosageParseOk;
    int32_t cur_dphase_delta = 0;
    uint32_t hds_valid = 0;
    uint32_t dosage_int;
    if (ParseBcfBiallelicHds(dosage_main, hds_main, dosage_value_ct, hds_value_ct, sample_idx, is_haploid_or_0ploid, dosage_is_gp, import_dosage_certainty, &dpr, &dosage_int, &cur_dphase_delta, &hds_valid)) {
      if (unlikely(!dpr)) {
        return kBcfParseInvalidDosage;
      }
      // if dpr != kDosageParseForceMissing, DS and HDS are missing.
      if ((dpr != kDosageParseForceMissing) && BcfGtIsPhasedHet(gt_first, gt_second, 0x80, halfcall_mode)) {
        *phase_or_dosage_found_ptr = 1;
        return kBcfParseOk;
      }
    } else if (cur_dphase_delta) {
      const uint32_t dphase_side1 = dosage_int + cur_dphase_delta;
      const uint32_t dphase_side2 = dosage_int - cur_dphase_delta;
      const uint32_t dphase_halfdist1 = DphaseHalfdist(dphase_side1);
      const uint32_t dphase_halfdist2 = DphaseHalfdist(dphase_side2);
      const uint32_t dphase_erase_halfdist = dosage_erase_halfdist + kDosage4th;
      if ((dphase_halfdist1 < dphase_erase_halfdist) ||
          (dphase_halfdist2 < dphase_erase_halfdist) ||
          (((dphase_side1 + kDosageMid) ^ (dphase_side2 + kDosageMid)) & kDosageMax)) {
        *phase_or_dosage_found_ptr = 1;
        return kBcfParseOk;
      }
    } else {
      const uint32_t cur_halfdist = BiallelicDosageHalfdist(dosage_int);
      if ((cur_halfdist < dosage_erase_halfdist) ||
          ((!hds_valid) && BcfGtIsPhasedHet(gt_first, gt_second, 0x80, halfcall_mode))) {
        *phase_or_dosage_found_ptr = 1;
        return kBcfParseOk;
      }
    }
  }
  return kBcfParseOk;
}

void BcfConvertBiallelicHaploidGt(const unsigned char* gt_main, uint32_t sample_ct, uintptr_t* __restrict genovec) {
  // Haploid, single byte per genotype.
  //   0    -> 3
  //   2    -> 0
  //   4    -> 2
  //   0x81 -> 3
  // Other values invalid, but we don't promise exhaustive error-checking, we
  // just need to (i) not blow up and (ii) try to have consistent behavior
  // between build types.
  // single-byte values, subtract 2 with wraparound, MIN(result, 3) works.
#ifdef __LP64__
  const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
  genovec[sample_ctl2 - 1] = 0;
  // overread ok since we perform no error detection
  const uint32_t vec_ct = DivUp(sample_ct, kBytesPerVec);
  const VecUc* gt_alias = R_CAST(const VecUc*, gt_main);
  Vec4thUint* genovec_alias = R_CAST(Vec4thUint*, genovec);
  const VecUc neg2 = vecuc_set1(0xfe);
  const VecUc three = vecuc_set1(3);
  const VecW zero = vecw_setzero();
  for (uint32_t vidx = 0; vidx != vec_ct; ++vidx) {
    const VecUc bcf_bytes = vecuc_loadu(&(gt_alias[vidx]));
    const VecUc bcf_bytes_m2 = bcf_bytes + neg2;
    VecUc converted_bytes = vecuc_min(bcf_bytes_m2, three);
    // Now gather bits 0,1,8,9,... via movemask.
    // todo: better ARM implementation
    converted_bytes = vecuc_permute0xd8_if_avx2(converted_bytes);
    VecW vec_lo = vecw_unpacklo8(R_CAST(VecW, converted_bytes), zero);
    VecW vec_hi = vecw_unpackhi8(R_CAST(VecW, converted_bytes), zero);
    vec_lo = vecw_slli(vec_lo, 7) | vecw_slli(vec_lo, 14);
    vec_hi = vecw_slli(vec_hi, 7) | vecw_slli(vec_hi, 14);
    const Vec4thUint lo_bits = vecw_movemask(vec_lo);
    const Vec4thUint hi_bits = vecw_movemask(vec_hi);
    genovec_alias[vidx] = lo_bits | (hi_bits << kBytesPerVec);
  }
  // Zero out the garbage at the end.
  const uint32_t inv_remainder = (-sample_ct) & (kBytesPerVec - 1);
  genovec_alias[vec_ct - 1] &= (~S_CAST(Vec4thUint, 0)) >> (2 * inv_remainder);
#else
  const unsigned char* gt_iter = gt_main;
  const uint32_t word_ct_m1 = (sample_ct - 1) / kBitsPerWordD2;
  uint32_t loop_len = kBitsPerWordD2;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= word_ct_m1) {
      if (widx > word_ct_m1) {
        break;
      }
      loop_len = ModNz(sample_ct, kBitsPerWordD2);
    }
    uintptr_t geno_word = 0;
    // possible todo: benchmark this vs. the forward loop, with and without
    // sparse-optimization
    for (uint32_t sample_idx_lowbits = loop_len; sample_idx_lowbits; ) {
      --sample_idx_lowbits;
      unsigned char raw_val_m2 = gt_iter[sample_idx_lowbits] - 2;
      if (raw_val_m2 > 3) {
        raw_val_m2 = 3;
      }
      geno_word |= raw_val_m2;
      geno_word = geno_word << 2;
    }
    genovec[widx] = geno_word;
    gt_iter = &(gt_iter[kBitsPerWordD2]);
  }
#endif
}

BcfParseErr BcfApplyGqDpFilters(const BcfImportBaseContext* bibcp, const GparseReadBcfMetadata* metap, const unsigned char* record_start, uintptr_t* __restrict genovec) {
  if ((metap->qual_vec_offsets[0] != UINT32_MAX) || (metap->qual_vec_offsets[1] != UINT32_MAX)) {
    const uint32_t sample_ct = bibcp->sample_ct;
    for (uint32_t qual_idx = 0; qual_idx != 2; ++qual_idx) {
      const uint32_t vec_offset = metap->qual_vec_offsets[qual_idx];
      if (vec_offset == UINT32_MAX) {
        continue;
      }
      const unsigned char* qual_start = &(record_start[vec_offset * S_CAST(uintptr_t, kBytesPerVec)]);
      BcfParseErr bcf_parse_err = BcfParseGqDpAligned(qual_start, sample_ct, bibcp->qual_mins[qual_idx], bibcp->qual_maxs[qual_idx], genovec);
      if (bcf_parse_err != kBcfParseOk) {
        return bcf_parse_err;
      }
    }
  }
  return kBcfParseOk;
}

// GT present, dosage/HDS not present, phase known to be irrelevant.
BcfParseErr BcfConvertUnphasedBiallelic(const BcfImportBaseContext* bibcp, const GparseReadBcfMetadata* metap, unsigned char* record_start, uintptr_t* __restrict genovec) {
  unsigned char* gt_type_start = &(record_start[metap->gt_vec_offset * S_CAST(uintptr_t, kBytesPerVec)]);
  const unsigned char* gt_main = gt_type_start;
  // this was partially validated in first pass
  uint32_t value_type;
  uint32_t value_ct;
  ScanBcfTypeAligned(&gt_main, &value_type, &value_ct);
  if (unlikely(value_type != 1)) {
    return (value_type <= 3)? kBcfParseWideGt : kBcfParseMalformedGeneric;
  }
  const uint32_t sample_ct = bibcp->sample_ct;
  // value_ct guaranteed to be positive
  if (value_ct == 1) {
    BcfConvertBiallelicHaploidGt(gt_main, sample_ct, genovec);
  } else {
    const unsigned char* biallelic_gt_lookup = bibcp->biallelic_gt_lookup;
#ifdef USE_SHUFFLE8
    const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
    genovec[sample_ctl2 - 1] = 0;
#else
    const uint32_t word_ct_m1 = (sample_ct - 1) / kBitsPerWordD2;
    uint32_t loop_len = kBitsPerWordD2;
#endif
    if (value_ct == 2) {
      // Regular diploid.
      // Possibly phased, we ignore that since the prescan told us no phased
      // calls are hets.
      // After right-shifting, the four cases for each byte are 0, 1, 2, and
      // 0x40.  This adds up to 16 cases for the byte pair, which works nicely
      // with _mm{256}_shuffle_epi8() parallel-lookup.  Conveniently, all
      // missing, haploid, and non-error half-call handling logic can be
      // encoded in the table.
      // As for erroring out on half-calls, that's handled by setting the
      // relevant table entries to 4 instead of 0..3.
#ifdef USE_SHUFFLE8
      const uint32_t vec_ct = DivUp(sample_ct, kInt16PerVec);
      const uint32_t inv_remainder = (-sample_ct) & (kInt16PerVec - 1);
      if (inv_remainder) {
        // Set all past-the-end bytes to 2.  This causes trailing genovec bits
        // to be zero, and doesn't create a half_call_error_bit2_vec false
        // positive.
        // (yes, it's a bit silly to use gt_type_start instead of gt_main
        // here...)
        memset(&(gt_type_start[kBytesPerVec + sample_ct * 2]), 2, 2 * inv_remainder);
      }
      const VecU16* gt_alias = R_CAST(const VecU16*, gt_main);
      Vec8thUint* genovec_alias = R_CAST(Vec8thUint*, genovec);

      const VecU16 hibit_mask = vecu16_set1(0x7f7f);
      const VecU16 three = vecu16_set1(0x303);
      const VecU16 mask_000f = vecu16_set1(0xf);
      const VecU16 lookup_vec = vecu16_loadu(biallelic_gt_lookup);
      VecU16 half_call_error_bit2_vec = vecu16_setzero();
      for (uint32_t vidx = 0; vidx != vec_ct; ++vidx) {
        const VecU16 bcf_bytes = vecu16_loadu(&(gt_alias[vidx]));
        const VecU16 shifted_bytes = vecu16_srli(bcf_bytes, 1) & hibit_mask;
        const VecU16 capped_bytes = vecu16_min8(shifted_bytes, three);
        const VecU16 ready_for_lookup = (capped_bytes | vecu16_srli(capped_bytes, 6)) & mask_000f;
        const VecU16 lookup_result = vecu16_shuffle8(lookup_vec, ready_for_lookup) & mask_000f;
        half_call_error_bit2_vec |= lookup_result;
        // todo: better ARM implementation
        const VecU16 ready_for_movemask = vecu16_slli(lookup_result, 7) | vecu16_slli(lookup_result, 14);
        const Vec8thUint geno_bits = vecu16_movemask(ready_for_movemask);
        genovec_alias[vidx] = geno_bits;
      }
      half_call_error_bit2_vec = vecu16_slli(half_call_error_bit2_vec, 5);
      if (unlikely(vecu16_movemask(half_call_error_bit2_vec))) {
        return kBcfParseHalfCallError;
      }
#else
      // todo: benchmark explicit vectorization of SSE2; I suspect that, with
      // no parallel-lookup operation, it would usually be worse than this
      // simple sparse-optimization.
      const uint16_t* gt_iter = R_CAST(const uint16_t*, gt_main);
      uintptr_t half_call_error_bit2 = 0;
      for (uint32_t widx = 0; ; ++widx) {
        if (widx >= word_ct_m1) {
          if (widx > word_ct_m1) {
            break;
          }
          loop_len = ModNz(sample_ct, kBitsPerWordD2);
        }
        uintptr_t geno_word = 0;
        for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits) {
          const uint32_t phaseless_val = gt_iter[sample_idx_lowbits] & 0xfefe;
          if (phaseless_val == 0x202) {
            continue;
          }
          uint32_t first_allele_idx_p1 = (phaseless_val & 0xff) >> 1;
          uint32_t second_allele_idx_p1 = phaseless_val >> 9;
          if (first_allele_idx_p1 > 3) {
            first_allele_idx_p1 = 3;
          }
          if (second_allele_idx_p1 > 3) {
            second_allele_idx_p1 = 3;
          }
          const uintptr_t result = biallelic_gt_lookup[first_allele_idx_p1 + second_allele_idx_p1 * 4];
          half_call_error_bit2 |= result;
          geno_word |= result << (2 * sample_idx_lowbits);
        }
        genovec[widx] = geno_word;
        gt_iter = &(gt_iter[kBitsPerWordD2]);
      }
      if (unlikely(half_call_error_bit2 & 4)) {
        return kBcfParseHalfCallError;
      }
#endif
    } else {
      // triploid, etc.
      if (unlikely(bibcp->error_on_polyploid)) {
        return kBcfParsePolyploidError;
      }
#ifdef USE_SHUFFLE8
      const uint32_t word_ct_m1 = (sample_ct - 1) / kBitsPerWordD2;
      uint32_t loop_len = kBitsPerWordD2;
#endif
      const unsigned char* gt_iter = gt_main;
      uintptr_t half_call_error_bit2 = 0;
      for (uint32_t widx = 0; ; ++widx) {
        if (widx >= word_ct_m1) {
          if (widx > word_ct_m1) {
            break;
          }
          loop_len = ModNz(sample_ct, kBitsPerWordD2);
        }
        uintptr_t geno_word = 0;
        for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits) {
          // may overread 1 byte
          uint32_t phaseless_val;
          memcpy(&phaseless_val, gt_iter, 4);
          gt_iter = &(gt_iter[value_ct]);
          phaseless_val &= 0xfffefe;
          if (phaseless_val == 0x810202) {
            continue;
          }
          uintptr_t result;
          if ((phaseless_val & 0x810000) != 0x810000) {
            // triploid+ (or malformed)
            result = 3;
          } else {
            unsigned char first_allele_idx_p1 = phaseless_val >> 1;
            unsigned char second_allele_idx_p1 = (phaseless_val >> 9) & 0x7f;
            if (first_allele_idx_p1 > 3) {
              first_allele_idx_p1 = 3;
            }
            if (second_allele_idx_p1 > 3) {
              second_allele_idx_p1 = 3;
            }
            result = biallelic_gt_lookup[first_allele_idx_p1 + second_allele_idx_p1 * 4];
            half_call_error_bit2 |= result;
          }
          geno_word |= result << (2 * sample_idx_lowbits);
        }
        genovec[widx] = geno_word;
      }
      if (unlikely(half_call_error_bit2 & 4)) {
        return kBcfParseHalfCallError;
      }
    }
  }
  return BcfApplyGqDpFilters(bibcp, metap, record_start, genovec);
}

static_assert(sizeof(AlleleCode) == 1, "BcfConvertMultiallelicHaploidInt8Gt() needs to be updated.");
// genovec assumed to be initialized by BcfApplyGqDpFilters()
void BcfConvertMultiallelicHaploidInt8Gt(const unsigned char* gt_main, uint32_t sample_ct, uint32_t allele_ct, uintptr_t* __restrict genovec, Halfword* __restrict patch_10_set_alias, AlleleCode** patch_10_iterp) {
  const unsigned char* gt_iter = gt_main;
  const uint32_t word_ct_m1 = (sample_ct - 1) / kBitsPerWordD2;
  DoubleAlleleCode* patch_10_vals_alias = R_CAST(DoubleAlleleCode*, *patch_10_iterp);
  DoubleAlleleCode* patch_10_iter = patch_10_vals_alias;
  uint32_t loop_len = kBitsPerWordD2;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= word_ct_m1) {
      if (widx > word_ct_m1) {
        break;
      }
      loop_len = ModNz(sample_ct, kBitsPerWordD2);
    }
    const uintptr_t gq_dp_fail_word = genovec[widx];
    uintptr_t geno_word = 0;
    uint32_t patch_10_set_hw = 0;
    for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits) {
      const unsigned char raw_val_m2 = gt_iter[sample_idx_lowbits] - 2;
      if ((!raw_val_m2) || ((gq_dp_fail_word >> (2 * sample_idx_lowbits)) & 1)) {
        continue;
      }
      const uint32_t allele_idx = raw_val_m2 >> 1;
      uintptr_t result;
      if (allele_idx >= allele_ct) {
        // missing and END_OF_VECTOR wind up here
        result = 3;
      } else {
        result = 2;
        // note that allele_idx == 0 is possible with malformed raw_val == 3,
        // we don't handle it "correctly" but we do need to avoid blowing up
        if (allele_idx >= 2) {
          patch_10_set_hw |= 1U << sample_idx_lowbits;
          *patch_10_iter++ = allele_idx * 0x101;
        }
      }
      geno_word |= result << (2 * sample_idx_lowbits);
    }
    genovec[widx] = geno_word | gq_dp_fail_word;
    patch_10_set_alias[widx] = patch_10_set_hw;
    gt_iter = &(gt_iter[kBitsPerWordD2]);
  }
  *patch_10_iterp = R_CAST(AlleleCode*, patch_10_iter);
}

static_assert(sizeof(AlleleCode) == 1, "BcfConvertMultiallelicHaploidInt8Gt() needs to be updated.");
// genovec assumed to be initialized by BcfApplyGqDpFilters()
void BcfConvertMultiallelicHaploidInt16Gt(const unsigned char* gt_main, uint32_t sample_ct, uint32_t allele_ct, uintptr_t* __restrict genovec, Halfword* __restrict patch_10_set_alias, AlleleCode** patch_10_iterp) {
  const uint16_t* gt_iter = R_CAST(const uint16_t*, gt_main);
  const uint32_t word_ct_m1 = (sample_ct - 1) / kBitsPerWordD2;
  DoubleAlleleCode* patch_10_vals_alias = R_CAST(DoubleAlleleCode*, *patch_10_iterp);
  DoubleAlleleCode* patch_10_iter = patch_10_vals_alias;
  uint32_t loop_len = kBitsPerWordD2;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= word_ct_m1) {
      if (widx > word_ct_m1) {
        break;
      }
      loop_len = ModNz(sample_ct, kBitsPerWordD2);
    }
    const uintptr_t gq_dp_fail_word = genovec[widx];
    uintptr_t geno_word = 0;
    uint32_t patch_10_set_hw = 0;
    for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits) {
      // intentional underflow
      const uint32_t raw_val_m2 = S_CAST(uint32_t, gt_iter[sample_idx_lowbits]) - 2;
      if ((!raw_val_m2) || ((gq_dp_fail_word >> (2 * sample_idx_lowbits)) & 1)) {
        continue;
      }
      const uint32_t allele_idx = raw_val_m2 >> 1;
      uintptr_t result;
      if (allele_idx >= allele_ct) {
        result = 3;
      } else {
        result = 2;
        if (allele_idx >= 2) {
          patch_10_set_hw |= 1U << sample_idx_lowbits;
          *patch_10_iter++ = allele_idx * 0x101;
        }
      }
      geno_word |= result << (2 * sample_idx_lowbits);
    }
    genovec[widx] = geno_word | gq_dp_fail_word;
    patch_10_set_alias[widx] = patch_10_set_hw;
    gt_iter = &(gt_iter[kBitsPerWordD2]);
  }
  *patch_10_iterp = R_CAST(AlleleCode*, patch_10_iter);
}

BcfParseErr BcfConvertUnphasedMultiallelic(const BcfImportBaseContext* bibcp, const GparseReadBcfMetadata* metap, const unsigned char* record_start, uint32_t allele_ct, uint32_t* __restrict patch_01_ctp, uint32_t* __restrict patch_10_ctp, uintptr_t* __restrict genovec, uintptr_t* __restrict patch_01_set, AlleleCode* __restrict patch_01_vals, uintptr_t* __restrict patch_10_set, AlleleCode* __restrict patch_10_vals) {
  const uint32_t sample_ct = bibcp->sample_ct;
  const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
  ZeroWArr(sample_ctl2, genovec);
  BcfParseErr bcf_parse_err = BcfApplyGqDpFilters(bibcp, metap, record_start, genovec);
  if (unlikely(bcf_parse_err != kBcfParseOk)) {
    return bcf_parse_err;
  }
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  const uint32_t word_ct_m1 = (sample_ct - 1) / kBitsPerWordD2;
  patch_01_set[sample_ctl - 1] = 0;
  patch_10_set[sample_ctl - 1] = 0;
  Halfword* patch_01_set_alias = R_CAST(Halfword*, patch_01_set);
  Halfword* patch_10_set_alias = R_CAST(Halfword*, patch_10_set);
  const VcfHalfCall halfcall_mode = bibcp->halfcall_mode;
  AlleleCode* patch_01_iter = patch_01_vals;
  AlleleCode* patch_10_iter = patch_10_vals;
  uint32_t loop_len = kBitsPerWordD2;
  const unsigned char* gt_main = &(record_start[metap->gt_vec_offset * S_CAST(uintptr_t, kBytesPerVec)]);
  uint32_t value_type;
  uint32_t value_ct;
  ScanBcfTypeAligned(&gt_main, &value_type, &value_ct);
  if (value_type == 1) {
    // int8
    if (allele_ct >= 64) {
      // don't misinterpret END_OF_VECTOR
      allele_ct = 63;
    }
    if (value_ct == 1) {
      ZeroWArr(sample_ctl, patch_01_set);
      BcfConvertMultiallelicHaploidInt8Gt(gt_main, sample_ct, allele_ct, genovec, patch_10_set_alias, &patch_10_iter);
    } else if (value_ct == 2) {
      const uint16_t* gt_iter = R_CAST(const uint16_t*, gt_main);
      for (uint32_t widx = 0; ; ++widx) {
        if (widx >= word_ct_m1) {
          if (widx > word_ct_m1) {
            break;
          }
          loop_len = ModNz(sample_ct, kBitsPerWordD2);
        }
        const uintptr_t gq_dp_fail_word = genovec[widx];
        uintptr_t geno_word = 0;
        uint32_t patch_01_set_hw = 0;
        uint32_t patch_10_set_hw = 0;
        for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits) {
          const uint32_t phaseless_val = gt_iter[sample_idx_lowbits] & 0xfefe;
          if ((phaseless_val == 0x202) || ((gq_dp_fail_word >> (2 * sample_idx_lowbits)) & 1)) {
            continue;
          }
          const uint32_t shifted_phaseless_val = phaseless_val >> 1;
          uint32_t first_allele_idx_p1 = shifted_phaseless_val & 0xff;
          uintptr_t result;
          if (first_allele_idx_p1 > allele_ct) {
            // assume end-of-vector
            result = 3;
          } else {
            uint32_t second_allele_idx_p1 = shifted_phaseless_val >> 8;
            const uint32_t cur_bit = 1U << sample_idx_lowbits;
            if (second_allele_idx_p1 > allele_ct) {
              // haploid
              if (!first_allele_idx_p1) {
                result = 3;
              } else {
                const uint32_t first_allele_idx = first_allele_idx_p1 - 1;
                if (first_allele_idx < 2) {
                  result = first_allele_idx * 2;
                } else {
                  result = 2;
                  patch_10_set_hw |= cur_bit;
                  *patch_10_iter++ = first_allele_idx;
                  *patch_10_iter++ = first_allele_idx;
                }
              }
            } else {
              // diploid
              if (second_allele_idx_p1 < first_allele_idx_p1) {
                const uint32_t uii = first_allele_idx_p1;
                first_allele_idx_p1 = second_allele_idx_p1;
                second_allele_idx_p1 = uii;
              }
              if (!first_allele_idx_p1) {
                // missing or half-call
                if ((!second_allele_idx_p1) || (halfcall_mode == kVcfHalfCallMissing)) {
                  result = 3;
                } else if (unlikely(halfcall_mode == kVcfHalfCallError)) {
                  return kBcfParseHalfCallError;
                } else {
                  const uint32_t second_allele_idx = second_allele_idx_p1 - 1;
                  if (second_allele_idx < 2) {
                    // kVcfHalfCallHaploid, kVcfHalfCallReference
                    result = second_allele_idx << halfcall_mode;
                  } else if (halfcall_mode == kVcfHalfCallReference) {
                    result = 1;
                    patch_01_set_hw |= cur_bit;
                    *patch_01_iter++ = second_allele_idx;
                  } else {
                    result = 2;
                    patch_10_set_hw |= cur_bit;
                    *patch_10_iter++ = second_allele_idx;
                    *patch_10_iter++ = second_allele_idx;
                  }
                }
              } else if (first_allele_idx_p1 == 1) {
                // note that we've handled the second_allele_idx_p1 == 0 and ==
                // 1 cases earlier
                result = 1;
                if (second_allele_idx_p1 > 2) {
                  patch_01_set_hw |= cur_bit;
                  *patch_01_iter++ = second_allele_idx_p1 - 1;
                }
              } else {
                result = 2;
                if (second_allele_idx_p1 > 2) {
                  patch_10_set_hw |= cur_bit;
                  *patch_10_iter++ = first_allele_idx_p1 - 1;
                  *patch_10_iter++ = second_allele_idx_p1 - 1;
                }
              }
            }
          }
          geno_word |= result << (2 * sample_idx_lowbits);
        }
        genovec[widx] = geno_word | gq_dp_fail_word;
        patch_01_set_alias[widx] = patch_01_set_hw;
        patch_10_set_alias[widx] = patch_10_set_hw;
        gt_iter = &(gt_iter[kBitsPerWordD2]);
      }
    } else {
      const unsigned char* gt_iter = gt_main;
      for (uint32_t widx = 0; ; ++widx) {
        if (widx >= word_ct_m1) {
          if (widx > word_ct_m1) {
            break;
          }
          loop_len = ModNz(sample_ct, kBitsPerWordD2);
        }
        const uintptr_t gq_dp_fail_word = genovec[widx];
        uintptr_t geno_word = 0;
        uint32_t patch_01_set_hw = 0;
        uint32_t patch_10_set_hw = 0;
        for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits) {
          // may overread 1 byte
          uint32_t phaseless_val;
          memcpy(&phaseless_val, gt_iter, 4);
          gt_iter = &(gt_iter[value_ct]);
          phaseless_val &= 0xfffefe;
          if ((phaseless_val == 0x810202) || ((gq_dp_fail_word >> (2 * sample_idx_lowbits)) & 1)) {
            continue;
          }
          uintptr_t result;
          if ((phaseless_val & 0x810000) != 0x810000) {
            // triploid+ (or malformed)
            result = 3;
          } else {
            phaseless_val &= 0xffff;
            const uint32_t shifted_phaseless_val = phaseless_val >> 1;
            // rest of this is identical to diploid case
            // todo: try making this a static inline function (didn't start
            // with that since the function would have 10+ parameters, many of
            // which may not be touched on a given function call)
            uint32_t first_allele_idx_p1 = shifted_phaseless_val;
            if (first_allele_idx_p1 > allele_ct) {
              // assume end-of-vector
              result = 3;
            } else {
              uint32_t second_allele_idx_p1 = shifted_phaseless_val >> 8;
              const uint32_t cur_bit = 1U << sample_idx_lowbits;
              if (second_allele_idx_p1 > allele_ct) {
                // haploid
                if (!first_allele_idx_p1) {
                  result = 3;
                } else {
                  const uint32_t first_allele_idx = first_allele_idx_p1 - 1;
                  if (first_allele_idx < 2) {
                    result = first_allele_idx * 2;
                  } else {
                    result = 2;
                    patch_10_set_hw |= cur_bit;
                    *patch_10_iter++ = first_allele_idx;
                    *patch_10_iter++ = first_allele_idx;
                  }
                }
              } else {
                // diploid
                if (second_allele_idx_p1 < first_allele_idx_p1) {
                  const uint32_t uii = first_allele_idx_p1;
                  first_allele_idx_p1 = second_allele_idx_p1;
                  second_allele_idx_p1 = uii;
                }
                if (!first_allele_idx_p1) {
                  // missing or half-call
                  if ((!second_allele_idx_p1) || (halfcall_mode == kVcfHalfCallMissing)) {
                    result = 3;
                  } else if (unlikely(halfcall_mode == kVcfHalfCallError)) {
                    return kBcfParseHalfCallError;
                  } else {
                    const uint32_t second_allele_idx = second_allele_idx_p1 - 1;
                    if (second_allele_idx < 2) {
                      // kVcfHalfCallHaploid, kVcfHalfCallReference
                      result = second_allele_idx << halfcall_mode;
                    } else if (halfcall_mode == kVcfHalfCallReference) {
                      result = 1;
                      patch_01_set_hw |= cur_bit;
                      *patch_01_iter++ = second_allele_idx;
                    } else {
                      result = 2;
                      patch_10_set_hw |= cur_bit;
                      *patch_10_iter++ = second_allele_idx;
                      *patch_10_iter++ = second_allele_idx;
                    }
                  }
                } else if (first_allele_idx_p1 == 1) {
                  result = 1;
                  if (second_allele_idx_p1 > 2) {
                    patch_01_set_hw |= cur_bit;
                    *patch_01_iter++ = second_allele_idx_p1 - 1;
                  }
                } else {
                  result = 2;
                  if (second_allele_idx_p1 > 2) {
                    patch_10_set_hw |= cur_bit;
                    *patch_10_iter++ = first_allele_idx_p1 - 1;
                    *patch_10_iter++ = second_allele_idx_p1 - 1;
                  }
                }
              }
            }
          }
          geno_word |= result << (2 * sample_idx_lowbits);
        }
        genovec[widx] = geno_word | gq_dp_fail_word;
        patch_01_set_alias[widx] = patch_01_set_hw;
        patch_10_set_alias[widx] = patch_10_set_hw;
      }
    }
  } else if (likely(value_type == 2)) {
    // int16
    if (value_ct == 1) {
      BcfConvertMultiallelicHaploidInt16Gt(gt_main, sample_ct, allele_ct, genovec, patch_10_set_alias, &patch_10_iter);
    } else if (value_ct == 2) {
      const uint32_t* gt_iter = R_CAST(const uint32_t*, gt_main);
      for (uint32_t widx = 0; ; ++widx) {
        if (widx >= word_ct_m1) {
          if (widx > word_ct_m1) {
            break;
          }
          loop_len = ModNz(sample_ct, kBitsPerWordD2);
        }
        const uintptr_t gq_dp_fail_word = genovec[widx];
        uintptr_t geno_word = 0;
        uint32_t patch_01_set_hw = 0;
        uint32_t patch_10_set_hw = 0;
        for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits) {
          const uint32_t phaseless_val = gt_iter[sample_idx_lowbits] & 0xfffefffeU;
          if ((phaseless_val == 0x20002) || ((gq_dp_fail_word >> (2 * sample_idx_lowbits)) & 1)) {
            continue;
          }
          const uint32_t shifted_phaseless_val = phaseless_val >> 1;
          uint32_t first_allele_idx_p1 = shifted_phaseless_val & 0xffff;
          uintptr_t result;
          if (first_allele_idx_p1 > allele_ct) {
            // assume end-of-vector
            result = 3;
          } else {
            uint32_t second_allele_idx_p1 = shifted_phaseless_val >> 16;
            const uint32_t cur_bit = 1U << sample_idx_lowbits;
            if (second_allele_idx_p1 > allele_ct) {
              // haploid
              if (!first_allele_idx_p1) {
                result = 3;
              } else {
                const uint32_t first_allele_idx = first_allele_idx_p1 - 1;
                if (first_allele_idx < 2) {
                  result = first_allele_idx * 2;
                } else {
                  result = 2;
                  patch_10_set_hw |= cur_bit;
                  *patch_10_iter++ = first_allele_idx;
                  *patch_10_iter++ = first_allele_idx;
                }
              }
            } else {
              // diploid
              if (second_allele_idx_p1 < first_allele_idx_p1) {
                const uint32_t uii = first_allele_idx_p1;
                first_allele_idx_p1 = second_allele_idx_p1;
                second_allele_idx_p1 = uii;
              }
              if (!first_allele_idx_p1) {
                // missing or half-call
                if ((!second_allele_idx_p1) || (halfcall_mode == kVcfHalfCallMissing)) {
                  result = 3;
                } else if (unlikely(halfcall_mode == kVcfHalfCallError)) {
                  return kBcfParseHalfCallError;
                } else {
                  const uint32_t second_allele_idx = second_allele_idx_p1 - 1;
                  if (second_allele_idx < 2) {
                    // kVcfHalfCallHaploid, kVcfHalfCallReference
                    result = second_allele_idx << halfcall_mode;
                  } else if (halfcall_mode == kVcfHalfCallReference) {
                    result = 1;
                    patch_01_set_hw |= cur_bit;
                    *patch_01_iter++ = second_allele_idx;
                  } else {
                    result = 2;
                    patch_10_set_hw |= cur_bit;
                    *patch_10_iter++ = second_allele_idx;
                    *patch_10_iter++ = second_allele_idx;
                  }
                }
              } else if (first_allele_idx_p1 == 1) {
                // note that we've handled the second_allele_idx_p1 == 0 and ==
                // 1 cases earlier
                result = 1;
                if (second_allele_idx_p1 > 2) {
                  patch_01_set_hw |= cur_bit;
                  *patch_01_iter++ = second_allele_idx_p1 - 1;
                }
              } else {
                result = 2;
                if (second_allele_idx_p1 > 2) {
                  patch_10_set_hw |= cur_bit;
                  *patch_10_iter++ = first_allele_idx_p1 - 1;
                  *patch_10_iter++ = second_allele_idx_p1 - 1;
                }
              }
            }
          }
          geno_word |= result << (2 * sample_idx_lowbits);
        }
        genovec[widx] = geno_word | gq_dp_fail_word;
        patch_01_set_alias[widx] = patch_01_set_hw;
        patch_10_set_alias[widx] = patch_10_set_hw;
        gt_iter = &(gt_iter[kBitsPerWordD2]);
      }
    } else {
      const uint16_t* gt_iter = R_CAST(const uint16_t*, gt_main);
      for (uint32_t widx = 0; ; ++widx) {
        if (widx >= word_ct_m1) {
          if (widx > word_ct_m1) {
            break;
          }
          loop_len = ModNz(sample_ct, kBitsPerWordD2);
        }
        const uintptr_t gq_dp_fail_word = genovec[widx];
        uintptr_t geno_word = 0;
        uint32_t patch_01_set_hw = 0;
        uint32_t patch_10_set_hw = 0;
        for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits) {
          // may overread 2 bytes
          uint64_t phaseless_val;
          memcpy(&phaseless_val, gt_iter, 8);
          gt_iter = &(gt_iter[value_ct]);
          phaseless_val &= 0xfffffffefffeLLU;
          if ((phaseless_val == 0x800100020002LLU) || ((gq_dp_fail_word >> (2 * sample_idx_lowbits)) & 1)) {
            continue;
          }
          uintptr_t result;
          if ((phaseless_val & 0x800100000000LLU) != 0x800100000000LLU) {
            // triploid+ (or malformed)
            result = 3;
          } else {
            const uint32_t shifted_phaseless_val = S_CAST(uint32_t, phaseless_val) >> 1;
            // rest of this is identical to diploid case
            uint32_t first_allele_idx_p1 = shifted_phaseless_val & 0xffff;
            if (first_allele_idx_p1 > allele_ct) {
              // assume end-of-vector
              result = 3;
            } else {
              uint32_t second_allele_idx_p1 = shifted_phaseless_val >> 16;
              const uint32_t cur_bit = 1U << sample_idx_lowbits;
              if (second_allele_idx_p1 > allele_ct) {
                // haploid
                if (!first_allele_idx_p1) {
                  result = 3;
                } else {
                  const uint32_t first_allele_idx = first_allele_idx_p1 - 1;
                  if (first_allele_idx < 2) {
                    result = first_allele_idx * 2;
                  } else {
                    result = 2;
                    patch_10_set_hw |= cur_bit;
                    *patch_10_iter++ = first_allele_idx;
                    *patch_10_iter++ = first_allele_idx;
                  }
                }
              } else {
                // diploid
                if (second_allele_idx_p1 < first_allele_idx_p1) {
                  const uint32_t uii = first_allele_idx_p1;
                  first_allele_idx_p1 = second_allele_idx_p1;
                  second_allele_idx_p1 = uii;
                }
                if (!first_allele_idx_p1) {
                  // missing or half-call
                  if ((!second_allele_idx_p1) || (halfcall_mode == kVcfHalfCallMissing)) {
                    result = 3;
                  } else if (unlikely(halfcall_mode == kVcfHalfCallError)) {
                    return kBcfParseHalfCallError;
                  } else {
                    const uint32_t second_allele_idx = second_allele_idx_p1 - 1;
                    if (second_allele_idx < 2) {
                      // kVcfHalfCallHaploid, kVcfHalfCallReference
                      result = second_allele_idx << halfcall_mode;
                    } else if (halfcall_mode == kVcfHalfCallReference) {
                      result = 1;
                      patch_01_set_hw |= cur_bit;
                      *patch_01_iter++ = second_allele_idx;
                    } else {
                      result = 2;
                      patch_10_set_hw |= cur_bit;
                      *patch_10_iter++ = second_allele_idx;
                      *patch_10_iter++ = second_allele_idx;
                    }
                  }
                } else if (first_allele_idx_p1 == 1) {
                  result = 1;
                  if (second_allele_idx_p1 > 2) {
                    patch_01_set_hw |= cur_bit;
                    *patch_01_iter++ = second_allele_idx_p1 - 1;
                  }
                } else {
                  result = 2;
                  if (second_allele_idx_p1 > 2) {
                    patch_10_set_hw |= cur_bit;
                    *patch_10_iter++ = first_allele_idx_p1 - 1;
                    *patch_10_iter++ = second_allele_idx_p1 - 1;
                  }
                }
              }
            }
          }
          geno_word |= result << (2 * sample_idx_lowbits);
        }
        genovec[widx] = geno_word | gq_dp_fail_word;
        patch_01_set_alias[widx] = patch_01_set_hw;
        patch_10_set_alias[widx] = patch_10_set_hw;
      }
    }
  } else {
    return (value_type == 3)? kBcfParseWideGt : kBcfParseMalformedGeneric;
  }
  *patch_01_ctp = patch_01_iter - patch_01_vals;
  *patch_10_ctp = S_CAST(uintptr_t, patch_10_iter - patch_10_vals) / 2;
  return kBcfParseOk;
}

// GT present, phase matters, dosage/HDS not present.
BcfParseErr BcfConvertPhasedBiallelic(const BcfImportBaseContext* bibcp, const GparseReadBcfMetadata* metap, unsigned char* record_start, uintptr_t* __restrict genovec, uintptr_t* __restrict phasepresent, uintptr_t* __restrict phaseinfo) {
  unsigned char* gt_type_start = &(record_start[metap->gt_vec_offset * S_CAST(uintptr_t, kBytesPerVec)]);
  const unsigned char* gt_main = gt_type_start;
  // this was partially validated in first pass
  uint32_t value_type;
  uint32_t value_ct;
  ScanBcfTypeAligned(&gt_main, &value_type, &value_ct);
  if (unlikely(value_type != 1)) {
    return (value_type <= 3)? kBcfParseWideGt : kBcfParseMalformedGeneric;
  }
  const uint32_t sample_ct = bibcp->sample_ct;
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  if (value_ct == 1) {
    ZeroWArr(sample_ctl, phasepresent);
    BcfConvertBiallelicHaploidGt(gt_main, sample_ct, genovec);
  } else {
    phasepresent[sample_ctl - 1] = 0;
    phaseinfo[sample_ctl - 1] = 0;
    const unsigned char* biallelic_gt_lookup = bibcp->biallelic_gt_lookup;
#ifdef USE_SHUFFLE8
    const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
    genovec[sample_ctl2 - 1] = 0;
#else
    const uint32_t word_ct_m1 = (sample_ct - 1) / kBitsPerWordD2;
    uint32_t loop_len = kBitsPerWordD2;
#endif
    if (value_ct == 2) {
      // Regular diploid.
      // Additional requirements beyond unphased case:
      // - Extract original bit 8 (low bit of second genotype value).  No need
      //   to do anything else (verify non-END_OF_VECTOR or het status) with
      //   it, here, that's taken care of after the GQ/DP filter.
      // - Extract original bit 9 (since the variant is biallelic, this bit
      //   should only be set when the second genotype value is 0, i.e. if the
      //   genotype is a phased het, it's 1|0).
#ifdef USE_SHUFFLE8
      const uint32_t vec_ct = DivUp(sample_ct, kInt16PerVec);
      const uint32_t inv_remainder = (-sample_ct) & (kInt16PerVec - 1);
      if (inv_remainder) {
        // Set all past-the-end bytes to 2.  This causes trailing genovec and
        // phasepresent bits to be zero, and doesn't create a
        // half_call_error_bit2_vec false positive.
        memset(&(gt_type_start[kBytesPerVec + sample_ct * 2]), 2, 2 * inv_remainder);
      }
      const VecU16* gt_alias = R_CAST(const VecU16*, gt_main);
      Vec8thUint* genovec_alias = R_CAST(Vec8thUint*, genovec);
      Vec16thUint* phasepresent_alias = R_CAST(Vec16thUint*, phasepresent);
      Vec16thUint* phaseinfo_alias = R_CAST(Vec16thUint*, phaseinfo);

      const VecU16 hibit_mask = vecu16_set1(0x7f7f);
      const VecU16 three = vecu16_set1(0x303);
      const VecU16 mask_000f = vecu16_set1(0xf);
#  ifndef USE_AVX2
      const VecU16 gather_even = vecu16_setr8(0, 2, 4, 6, 8, 10, 12, 14,
                                              -1, -1, -1, -1, -1, -1, -1, -1);
#  endif
      const VecU16 lookup_vec = vecu16_loadu(biallelic_gt_lookup);
      VecU16 half_call_error_bit2_vec = vecu16_setzero();
      for (uint32_t vidx = 0; vidx != vec_ct; ++vidx) {
        const VecU16 bcf_bytes = vecu16_loadu(&(gt_alias[vidx]));
        const VecU16 shifted_unmasked_bytes = vecu16_srli(bcf_bytes, 1);
        const VecU16 bit9s_at_bit7 = vecu16_srli(bcf_bytes, 2);
#  ifdef USE_AVX2
        const Vec16thUint cur_phasepresent = _pext_u32(vecu16_movemask(shifted_unmasked_bytes), 0x55555555);
        const Vec16thUint cur_phaseinfo = _pext_u32(vecu16_movemask(bit9s_at_bit7), 0x55555555);
#  else
        // todo: better ARM implementation
        const VecU16 phasepresent_vec = vecu16_shuffle8(shifted_unmasked_bytes, gather_even);
        const VecU16 phaseinfo_vec = vecu16_shuffle8(bit9s_at_bit7, gather_even);
        const Vec16thUint cur_phasepresent = vecu16_movemask(phasepresent_vec);
        const Vec16thUint cur_phaseinfo = vecu16_movemask(phaseinfo_vec);
#  endif
        phasepresent_alias[vidx] = cur_phasepresent;
        phaseinfo_alias[vidx] = cur_phaseinfo;

        const VecU16 shifted_bytes = shifted_unmasked_bytes & hibit_mask;
        const VecU16 capped_bytes = vecu16_min8(shifted_bytes, three);
        const VecU16 ready_for_lookup = (capped_bytes | vecu16_srli(capped_bytes, 6)) & mask_000f;
        const VecU16 lookup_result = vecu16_shuffle8(lookup_vec, ready_for_lookup) & mask_000f;
        half_call_error_bit2_vec |= lookup_result;
        const VecU16 ready_for_movemask = vecu16_slli(lookup_result, 7) | vecu16_slli(lookup_result, 14);
        const Vec8thUint geno_bits = vecu16_movemask(ready_for_movemask);
        genovec_alias[vidx] = geno_bits;
      }
      half_call_error_bit2_vec = vecu16_slli(half_call_error_bit2_vec, 5);
      if (unlikely(vecu16_movemask(half_call_error_bit2_vec))) {
        return kBcfParseHalfCallError;
      }
#else
      // todo: benchmark explicit vectorization of SSE2
      const uint16_t* gt_iter = R_CAST(const uint16_t*, gt_main);
      Halfword* phasepresent_alias = R_CAST(Halfword*, phasepresent);
      Halfword* phaseinfo_alias = R_CAST(Halfword*, phaseinfo);
      uintptr_t half_call_error_bit2 = 0;
      for (uint32_t widx = 0; ; ++widx) {
        if (widx >= word_ct_m1) {
          if (widx > word_ct_m1) {
            break;
          }
          loop_len = ModNz(sample_ct, kBitsPerWordD2);
        }
        uintptr_t geno_word = 0;
        uintptr_t phasepresent_hw_shifted = 0;
        Halfword phaseinfo_hw = 0;
        for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits) {
          const uintptr_t raw_val = gt_iter[sample_idx_lowbits];
          if ((raw_val & 0xfefe) == 0x202) {
            continue;
          }
          const uintptr_t bit_8 = raw_val & 0x100;
          phasepresent_hw_shifted |= bit_8 << sample_idx_lowbits;
          uint32_t second_allele_idx_p1 = raw_val >> 9;
          phaseinfo_hw |= (second_allele_idx_p1 & 1) << sample_idx_lowbits;
          uint32_t first_allele_idx_p1 = (raw_val & 0xff) >> 1;
          if (first_allele_idx_p1 > 3) {
            first_allele_idx_p1 = 3;
          }
          if (second_allele_idx_p1 > 3) {
            second_allele_idx_p1 = 3;
          }
          const uintptr_t result = biallelic_gt_lookup[first_allele_idx_p1 + second_allele_idx_p1 * 4];
          half_call_error_bit2 |= result;
          geno_word |= result << (2 * sample_idx_lowbits);
        }
        genovec[widx] = geno_word;
        phasepresent_alias[widx] = phasepresent_hw_shifted >> 8;
        phaseinfo_alias[widx] = phaseinfo_hw;
        gt_iter = &(gt_iter[kBitsPerWordD2]);
      }
      if (unlikely(half_call_error_bit2 & 4)) {
        return kBcfParseHalfCallError;
      }
#endif
    } else {
      // triploid, etc.
      if (unlikely(bibcp->error_on_polyploid)) {
        return kBcfParsePolyploidError;
      }
      Halfword* phasepresent_alias = R_CAST(Halfword*, phasepresent);
      Halfword* phaseinfo_alias = R_CAST(Halfword*, phaseinfo);
#ifdef USE_SHUFFLE8
      const uint32_t word_ct_m1 = (sample_ct - 1) / kBitsPerWordD2;
      uint32_t loop_len = kBitsPerWordD2;
#endif
      const unsigned char* gt_iter = gt_main;
      uintptr_t half_call_error_bit2 = 0;
      for (uint32_t widx = 0; ; ++widx) {
        if (widx >= word_ct_m1) {
          if (widx > word_ct_m1) {
            break;
          }
          loop_len = ModNz(sample_ct, kBitsPerWordD2);
        }
        uintptr_t geno_word = 0;
        uintptr_t phasepresent_hw_shifted = 0;
        Halfword phaseinfo_hw = 0;
        for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits) {
          // may overread 1 byte
          uint32_t raw_val;
          memcpy(&raw_val, gt_iter, 4);
          gt_iter = &(gt_iter[value_ct]);
          if ((raw_val & 0xfffefe) == 0x810202) {
            continue;
          }
          uintptr_t result;
          if ((raw_val & 0x810000) != 0x810000) {
            // triploid (or malformed)
            result = 3;
          } else {
            const uintptr_t bit_8 = raw_val & 0x100;
            phasepresent_hw_shifted |= bit_8 << sample_idx_lowbits;
            uint32_t second_allele_idx_p1 = (raw_val >> 9) & 0x7f;
            phaseinfo_hw |= (second_allele_idx_p1 & 1) << sample_idx_lowbits;
            uint32_t first_allele_idx_p1 = (raw_val & 0xff) >> 1;
            if (first_allele_idx_p1 > 3) {
              first_allele_idx_p1 = 3;
            }
            if (second_allele_idx_p1 > 3) {
              second_allele_idx_p1 = 3;
            }
            result = biallelic_gt_lookup[first_allele_idx_p1 + second_allele_idx_p1 * 4];
            half_call_error_bit2 |= result;
          }
          geno_word |= result << (2 * sample_idx_lowbits);
        }
        genovec[widx] = geno_word;
        phasepresent_alias[widx] = phasepresent_hw_shifted >> 8;
        phaseinfo_alias[widx] = phaseinfo_hw;
      }
      if (unlikely(half_call_error_bit2 & 4)) {
        return kBcfParseHalfCallError;
      }
    }
  }
  BcfParseErr bcf_parse_err = BcfApplyGqDpFilters(bibcp, metap, record_start, genovec);
  if (unlikely(bcf_parse_err != kBcfParseOk)) {
    return bcf_parse_err;
  }
  if (value_ct != 1) {
    const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
    // mask out all phasepresent bits which don't correspond to hets.
    MaskWordsToHalfwordsInvmatch(genovec, kMaskAAAA, sample_ctl2, phasepresent, phasepresent);
  }
  return kBcfParseOk;
}

BcfParseErr BcfConvertPhasedMultiallelic(const BcfImportBaseContext* bibcp, const GparseReadBcfMetadata* metap, unsigned char* record_start, uint32_t allele_ct, uint32_t* __restrict patch_01_ctp, uint32_t* __restrict patch_10_ctp, uintptr_t* __restrict genovec, uintptr_t* __restrict patch_01_set, AlleleCode* __restrict patch_01_vals, uintptr_t* __restrict patch_10_set, AlleleCode* __restrict patch_10_vals, uintptr_t* __restrict phasepresent, uintptr_t* __restrict phaseinfo) {
  // yes, there's lots of duplication with BcfConvertUnphasedMultiallelic...
  const uint32_t sample_ct = bibcp->sample_ct;
  const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
  ZeroWArr(sample_ctl2, genovec);
  BcfParseErr bcf_parse_err = BcfApplyGqDpFilters(bibcp, metap, record_start, genovec);
  if (unlikely(bcf_parse_err != kBcfParseOk)) {
    return bcf_parse_err;
  }
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  const uint32_t word_ct_m1 = (sample_ct - 1) / kBitsPerWordD2;
  patch_01_set[sample_ctl - 1] = 0;
  patch_10_set[sample_ctl - 1] = 0;
  Halfword* patch_01_set_alias = R_CAST(Halfword*, patch_01_set);
  Halfword* patch_10_set_alias = R_CAST(Halfword*, patch_10_set);
  phasepresent[sample_ctl - 1] = 0;
  Halfword* phasepresent_alias = R_CAST(Halfword*, phasepresent);
  Halfword* phaseinfo_alias = R_CAST(Halfword*, phaseinfo);
  const VcfHalfCall halfcall_mode = bibcp->halfcall_mode;
  AlleleCode* patch_01_iter = patch_01_vals;
  AlleleCode* patch_10_iter = patch_10_vals;
  uint32_t loop_len = kBitsPerWordD2;
  const unsigned char* gt_main = &(record_start[metap->gt_vec_offset * S_CAST(uintptr_t, kBytesPerVec)]);
  uint32_t value_type;
  uint32_t value_ct;
  ScanBcfTypeAligned(&gt_main, &value_type, &value_ct);
  if (unlikely((value_ct > 2) && bibcp->error_on_polyploid)) {
    return kBcfParsePolyploidError;
  }
  if (value_type == 1) {
    // int8
    if (allele_ct >= 64) {
      // don't misinterpret END_OF_VECTOR
      allele_ct = 63;
    }
    if (value_ct == 1) {
      ZeroWArr(sample_ctl, patch_01_set);
      ZeroWArr(sample_ctl, phasepresent);
      // phaseinfo doesn't matter
      BcfConvertMultiallelicHaploidInt8Gt(gt_main, sample_ct, allele_ct, genovec, patch_10_set_alias, &patch_10_iter);
    } else if (value_ct == 2) {
      const uint16_t* gt_iter = R_CAST(const uint16_t*, gt_main);
      for (uint32_t widx = 0; ; ++widx) {
        if (widx >= word_ct_m1) {
          if (widx > word_ct_m1) {
            break;
          }
          loop_len = ModNz(sample_ct, kBitsPerWordD2);
        }
        const uintptr_t gq_dp_fail_word = genovec[widx];
        uintptr_t geno_word = 0;
        uint32_t patch_01_set_hw = 0;
        uint32_t patch_10_set_hw = 0;
        uint32_t phasepresent_hw = 0;
        uint32_t phaseinfo_hw = 0;
        for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits) {
          const uint32_t raw_val = gt_iter[sample_idx_lowbits];
          const uint32_t phaseless_val = raw_val & 0xfefe;
          if ((phaseless_val == 0x202) || ((gq_dp_fail_word >> (2 * sample_idx_lowbits)) & 1)) {
            continue;
          }
          const uint32_t shifted_phaseless_val = phaseless_val >> 1;
          uint32_t first_allele_idx_p1 = shifted_phaseless_val & 0xff;
          uintptr_t result;
          if (first_allele_idx_p1 > allele_ct) {
            // assume end-of-vector
            result = 3;
          } else {
            uint32_t second_allele_idx_p1 = shifted_phaseless_val >> 8;
            const uint32_t cur_bit = 1U << sample_idx_lowbits;
            if (second_allele_idx_p1 > allele_ct) {
              // haploid
              if (!first_allele_idx_p1) {
                result = 3;
              } else {
                const uint32_t first_allele_idx = first_allele_idx_p1 - 1;
                if (first_allele_idx < 2) {
                  result = first_allele_idx * 2;
                } else {
                  result = 2;
                  patch_10_set_hw |= cur_bit;
                  *patch_10_iter++ = first_allele_idx;
                  *patch_10_iter++ = first_allele_idx;
                }
              }
            } else {
              // diploid
              if (second_allele_idx_p1 < first_allele_idx_p1) {
                const uint32_t uii = first_allele_idx_p1;
                first_allele_idx_p1 = second_allele_idx_p1;
                second_allele_idx_p1 = uii;
                phaseinfo_hw |= cur_bit;
              }
              if (!first_allele_idx_p1) {
                // missing or half-call
                if ((!second_allele_idx_p1) || (halfcall_mode == kVcfHalfCallMissing)) {
                  result = 3;
                } else if (unlikely(halfcall_mode == kVcfHalfCallError)) {
                  return kBcfParseHalfCallError;
                } else {
                  const uint32_t second_allele_idx = second_allele_idx_p1 - 1;
                  if (second_allele_idx < 2) {
                    // kVcfHalfCallHaploid, kVcfHalfCallReference
                    result = second_allele_idx << halfcall_mode;
                  } else if (halfcall_mode == kVcfHalfCallReference) {
                    result = 1;
                    patch_01_set_hw |= cur_bit;
                    *patch_01_iter++ = second_allele_idx;
                  } else {
                    result = 2;
                    patch_10_set_hw |= cur_bit;
                    *patch_10_iter++ = second_allele_idx;
                    *patch_10_iter++ = second_allele_idx;
                  }
                  if ((raw_val & 0x100) && (result == 1)) {
                    phasepresent_hw |= cur_bit;
                  }
                }
              } else {
                if ((raw_val & 0x100) && (first_allele_idx_p1 != second_allele_idx_p1)) {
                  phasepresent_hw |= cur_bit;
                }
                if (first_allele_idx_p1 == 1) {
                  // note that we've handled the second_allele_idx_p1 == 0 and
                  // == 1 cases earlier
                  result = 1;
                  if (second_allele_idx_p1 > 2) {
                    patch_01_set_hw |= cur_bit;
                    *patch_01_iter++ = second_allele_idx_p1 - 1;
                  }
                } else {
                  result = 2;
                  if (second_allele_idx_p1 > 2) {
                    patch_10_set_hw |= cur_bit;
                    *patch_10_iter++ = first_allele_idx_p1 - 1;
                    *patch_10_iter++ = second_allele_idx_p1 - 1;
                  }
                }
              }
            }
          }
          geno_word |= result << (2 * sample_idx_lowbits);
        }
        genovec[widx] = geno_word | gq_dp_fail_word;
        patch_01_set_alias[widx] = patch_01_set_hw;
        patch_10_set_alias[widx] = patch_10_set_hw;
        phasepresent_alias[widx] = phasepresent_hw;
        phaseinfo_alias[widx] = phaseinfo_hw;
        gt_iter = &(gt_iter[kBitsPerWordD2]);
      }
    } else {
      const unsigned char* gt_iter = gt_main;
      for (uint32_t widx = 0; ; ++widx) {
        if (widx >= word_ct_m1) {
          if (widx > word_ct_m1) {
            break;
          }
          loop_len = ModNz(sample_ct, kBitsPerWordD2);
        }
        const uintptr_t gq_dp_fail_word = genovec[widx];
        uintptr_t geno_word = 0;
        uint32_t patch_01_set_hw = 0;
        uint32_t patch_10_set_hw = 0;
        uint32_t phasepresent_hw = 0;
        uint32_t phaseinfo_hw = 0;
        for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits) {
          // may overread 1 byte
          uint32_t raw_val;
          memcpy(&raw_val, gt_iter, 4);
          gt_iter = &(gt_iter[value_ct]);
          const uint32_t phaseless_val = raw_val & 0xfffefe;
          if ((phaseless_val == 0x810202) || ((gq_dp_fail_word >> (2 * sample_idx_lowbits)) & 1)) {
            continue;
          }
          uintptr_t result;
          if ((phaseless_val & 0x810000) != 0x810000) {
            // triploid+ (or malformed)
            result = 3;
          } else {
            const uint32_t shifted_phaseless_val = (phaseless_val & 0xffff) >> 1;
            // rest of this is identical to diploid case
            uint32_t first_allele_idx_p1 = shifted_phaseless_val & 0xff;
            if (first_allele_idx_p1 > allele_ct) {
              // assume end-of-vector
              result = 3;
            } else {
              uint32_t second_allele_idx_p1 = shifted_phaseless_val >> 8;
              const uint32_t cur_bit = 1U << sample_idx_lowbits;
              if (second_allele_idx_p1 > allele_ct) {
                // haploid
                if (!first_allele_idx_p1) {
                  result = 3;
                } else {
                  const uint32_t first_allele_idx = first_allele_idx_p1 - 1;
                  if (first_allele_idx < 2) {
                    result = first_allele_idx * 2;
                  } else {
                    result = 2;
                    patch_10_set_hw |= cur_bit;
                    *patch_10_iter++ = first_allele_idx;
                    *patch_10_iter++ = first_allele_idx;
                  }
                }
              } else {
                // diploid
                if (second_allele_idx_p1 < first_allele_idx_p1) {
                  const uint32_t uii = first_allele_idx_p1;
                  first_allele_idx_p1 = second_allele_idx_p1;
                  second_allele_idx_p1 = uii;
                  phaseinfo_hw |= cur_bit;
                }
                if (!first_allele_idx_p1) {
                  // missing or half-call
                  if ((!second_allele_idx_p1) || (halfcall_mode == kVcfHalfCallMissing)) {
                    result = 3;
                  } else if (unlikely(halfcall_mode == kVcfHalfCallError)) {
                    return kBcfParseHalfCallError;
                  } else {
                    const uint32_t second_allele_idx = second_allele_idx_p1 - 1;
                    if (second_allele_idx < 2) {
                      // kVcfHalfCallHaploid, kVcfHalfCallReference
                      result = second_allele_idx << halfcall_mode;
                    } else if (halfcall_mode == kVcfHalfCallReference) {
                      result = 1;
                      patch_01_set_hw |= cur_bit;
                      *patch_01_iter++ = second_allele_idx;
                    } else {
                      result = 2;
                      patch_10_set_hw |= cur_bit;
                      *patch_10_iter++ = second_allele_idx;
                      *patch_10_iter++ = second_allele_idx;
                    }
                    if ((raw_val & 0x100) && (result == 1)) {
                      phasepresent_hw |= cur_bit;
                    }
                  }
                } else {
                  if ((raw_val & 0x100) && (first_allele_idx_p1 != second_allele_idx_p1)) {
                    phasepresent_hw |= cur_bit;
                  }
                  if (first_allele_idx_p1 == 1) {
                    result = 1;
                    if (second_allele_idx_p1 > 2) {
                      patch_01_set_hw |= cur_bit;
                      *patch_01_iter++ = second_allele_idx_p1 - 1;
                    }
                  } else {
                    result = 2;
                    if (second_allele_idx_p1 > 2) {
                      patch_10_set_hw |= cur_bit;
                      *patch_10_iter++ = first_allele_idx_p1 - 1;
                      *patch_10_iter++ = second_allele_idx_p1 - 1;
                    }
                  }
                }
              }
            }
          }
          geno_word |= result << (2 * sample_idx_lowbits);
        }
        genovec[widx] = geno_word | gq_dp_fail_word;
        patch_01_set_alias[widx] = patch_01_set_hw;
        patch_10_set_alias[widx] = patch_10_set_hw;
        phasepresent_alias[widx] = phasepresent_hw;
        phaseinfo_alias[widx] = phaseinfo_hw;
      }
    }
  } else if (likely(value_type == 2)) {
    // int16
    if (value_ct == 1) {
      ZeroWArr(sample_ctl, patch_01_set);
      ZeroWArr(sample_ctl, phasepresent);
      // phaseinfo doesn't matter
      BcfConvertMultiallelicHaploidInt16Gt(gt_main, sample_ct, allele_ct, genovec, patch_10_set_alias, &patch_10_iter);
    } else if (value_ct == 2) {
      const uint32_t* gt_iter = R_CAST(const uint32_t*, gt_main);
      for (uint32_t widx = 0; ; ++widx) {
        if (widx >= word_ct_m1) {
          if (widx > word_ct_m1) {
            break;
          }
          loop_len = ModNz(sample_ct, kBitsPerWordD2);
        }
        const uintptr_t gq_dp_fail_word = genovec[widx];
        uintptr_t geno_word = 0;
        uint32_t patch_01_set_hw = 0;
        uint32_t patch_10_set_hw = 0;
        uint32_t phasepresent_hw = 0;
        uint32_t phaseinfo_hw = 0;
        for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits) {
          const uint32_t raw_val = gt_iter[sample_idx_lowbits];
          const uint32_t phaseless_val = raw_val & 0xfffefffeU;
          if ((phaseless_val == 0x20002) || ((gq_dp_fail_word >> (2 * sample_idx_lowbits)) & 1)) {
            continue;
          }
          const uint32_t shifted_phaseless_val = phaseless_val >> 1;
          uint32_t first_allele_idx_p1 = shifted_phaseless_val & 0xffff;
          uintptr_t result;
          if (first_allele_idx_p1 > allele_ct) {
            // assume end-of-vector
            result = 3;
          } else {
            uint32_t second_allele_idx_p1 = shifted_phaseless_val >> 16;
            const uint32_t cur_bit = 1U << sample_idx_lowbits;
            if (second_allele_idx_p1 > allele_ct) {
              // haploid
              if (!first_allele_idx_p1) {
                result = 3;
              } else {
                const uint32_t first_allele_idx = first_allele_idx_p1 - 1;
                if (first_allele_idx < 2) {
                  result = first_allele_idx * 2;
                } else {
                  result = 2;
                  patch_10_set_hw |= cur_bit;
                  *patch_10_iter++ = first_allele_idx;
                  *patch_10_iter++ = first_allele_idx;
                }
              }
            } else {
              // diploid
              if (second_allele_idx_p1 < first_allele_idx_p1) {
                const uint32_t uii = first_allele_idx_p1;
                first_allele_idx_p1 = second_allele_idx_p1;
                second_allele_idx_p1 = uii;
                phaseinfo_hw |= cur_bit;
              }
              if (!first_allele_idx_p1) {
                // missing or half-call
                if ((!second_allele_idx_p1) || (halfcall_mode == kVcfHalfCallMissing)) {
                  result = 3;
                } else if (unlikely(halfcall_mode == kVcfHalfCallError)) {
                  return kBcfParseHalfCallError;
                } else {
                  const uint16_t second_allele_idx = second_allele_idx_p1 - 1;
                  if (second_allele_idx < 2) {
                    // kVcfHalfCallHaploid, kVcfHalfCallReference
                    result = second_allele_idx << halfcall_mode;
                  } else if (halfcall_mode == kVcfHalfCallReference) {
                    result = 1;
                    patch_01_set_hw |= cur_bit;
                    *patch_01_iter++ = second_allele_idx;
                  } else {
                    result = 2;
                    patch_10_set_hw |= cur_bit;
                    *patch_10_iter++ = second_allele_idx;
                    *patch_10_iter++ = second_allele_idx;
                  }
                  if ((raw_val & 0x10000) && (result == 1)) {
                    phasepresent_hw |= cur_bit;
                  }
                }
              } else {
                if ((raw_val & 0x10000) && (first_allele_idx_p1 != second_allele_idx_p1)) {
                  phasepresent_hw |= cur_bit;
                }
                if (first_allele_idx_p1 == 1) {
                  // note that we've handled the second_allele_idx_p1 == 0 and
                  // == 1 cases earlier
                  result = 1;
                  if (second_allele_idx_p1 > 2) {
                    patch_01_set_hw |= cur_bit;
                    *patch_01_iter++ = second_allele_idx_p1 - 1;
                  }
                } else {
                  result = 2;
                  if (second_allele_idx_p1 > 2) {
                    patch_10_set_hw |= cur_bit;
                    *patch_10_iter++ = first_allele_idx_p1 - 1;
                    *patch_10_iter++ = second_allele_idx_p1 - 1;
                  }
                }
              }
            }
          }
          geno_word |= result << (2 * sample_idx_lowbits);
        }
        genovec[widx] = geno_word | gq_dp_fail_word;
        patch_01_set_alias[widx] = patch_01_set_hw;
        patch_10_set_alias[widx] = patch_10_set_hw;
        phasepresent_alias[widx] = phasepresent_hw;
        phaseinfo_alias[widx] = phaseinfo_hw;
        gt_iter = &(gt_iter[kBitsPerWordD2]);
      }
    } else {
      const uint16_t* gt_iter = R_CAST(const uint16_t*, gt_main);
      for (uint32_t widx = 0; ; ++widx) {
        if (widx >= word_ct_m1) {
          if (widx > word_ct_m1) {
            break;
          }
          loop_len = ModNz(sample_ct, kBitsPerWordD2);
        }
        const uintptr_t gq_dp_fail_word = genovec[widx];
        uintptr_t geno_word = 0;
        uint32_t patch_01_set_hw = 0;
        uint32_t patch_10_set_hw = 0;
        uint32_t phasepresent_hw = 0;
        uint32_t phaseinfo_hw = 0;
        for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits) {
          // may overread 2 bytes
          uint64_t raw_val;
          memcpy(&raw_val, gt_iter, 8);
          gt_iter = &(gt_iter[value_ct]);
          const uint64_t phaseless_val = raw_val & 0xfffffffefffeLLU;
          if ((phaseless_val == 0x800100020002LLU) || ((gq_dp_fail_word >> (2 * sample_idx_lowbits)) & 1)) {
            continue;
          }
          uintptr_t result;
          if ((phaseless_val & 0x800100000000LLU) != 0x800100000000LLU) {
            // triploid+ (or malformed)
            result = 3;
          } else {
            const uint32_t shifted_phaseless_val = S_CAST(uint32_t, phaseless_val) >> 1;
            // rest of this is identical to diploid case
            uint32_t first_allele_idx_p1 = shifted_phaseless_val & 0xffff;
            if (first_allele_idx_p1 > allele_ct) {
              // assume end-of-vector
              result = 3;
            } else {
              uint32_t second_allele_idx_p1 = shifted_phaseless_val >> 16;
              const uint32_t cur_bit = 1U << sample_idx_lowbits;
              if (second_allele_idx_p1 > allele_ct) {
                // haploid
                if (!first_allele_idx_p1) {
                  result = 3;
                } else {
                  const uint32_t first_allele_idx = first_allele_idx_p1 - 1;
                  if (first_allele_idx < 2) {
                    result = first_allele_idx * 2;
                  } else {
                    result = 2;
                    patch_10_set_hw |= cur_bit;
                    *patch_10_iter++ = first_allele_idx;
                    *patch_10_iter++ = first_allele_idx;
                  }
                }
              } else {
                // diploid
                if (second_allele_idx_p1 < first_allele_idx_p1) {
                  const uint32_t uii = first_allele_idx_p1;
                  first_allele_idx_p1 = second_allele_idx_p1;
                  second_allele_idx_p1 = uii;
                  phaseinfo_hw |= cur_bit;
                }
                if (!first_allele_idx_p1) {
                  // missing or half-call
                  if ((!second_allele_idx_p1) || (halfcall_mode == kVcfHalfCallMissing)) {
                    result = 3;
                  } else if (unlikely(halfcall_mode == kVcfHalfCallError)) {
                    return kBcfParseHalfCallError;
                  } else {
                    const uint32_t second_allele_idx = second_allele_idx_p1 - 1;
                    if (second_allele_idx < 2) {
                      // kVcfHalfCallHaploid, kVcfHalfCallReference
                      result = second_allele_idx << halfcall_mode;
                    } else if (halfcall_mode == kVcfHalfCallReference) {
                      result = 1;
                      patch_01_set_hw |= cur_bit;
                      *patch_01_iter++ = second_allele_idx;
                    } else {
                      result = 2;
                      patch_10_set_hw |= cur_bit;
                      *patch_10_iter++ = second_allele_idx;
                      *patch_10_iter++ = second_allele_idx;
                    }
                    if ((raw_val & 0x10000) && (result == 1)) {
                      phasepresent_hw |= cur_bit;
                    }
                  }
                } else {
                  if ((raw_val & 0x10000) && (first_allele_idx_p1 != second_allele_idx_p1)) {
                    phasepresent_hw |= cur_bit;
                  }
                  if (first_allele_idx_p1 == 1) {
                    result = 1;
                    if (second_allele_idx_p1 > 2) {
                      patch_01_set_hw |= cur_bit;
                      *patch_01_iter++ = second_allele_idx_p1 - 1;
                    }
                  } else {
                    result = 2;
                    if (second_allele_idx_p1 > 2) {
                      patch_10_set_hw |= cur_bit;
                      *patch_10_iter++ = first_allele_idx_p1 - 1;
                      *patch_10_iter++ = second_allele_idx_p1 - 1;
                    }
                  }
                }
              }
            }
          }
        }
        genovec[widx] = geno_word | gq_dp_fail_word;
        patch_01_set_alias[widx] = patch_01_set_hw;
        patch_10_set_alias[widx] = patch_10_set_hw;
        phasepresent_alias[widx] = phasepresent_hw;
        phaseinfo_alias[widx] = phaseinfo_hw;
      }
    }
  } else {
    return (value_type == 3)? kBcfParseWideGt : kBcfParseMalformedGeneric;
  }
  *patch_01_ctp = patch_01_iter - patch_01_vals;
  *patch_10_ctp = S_CAST(uintptr_t, patch_10_iter - patch_10_vals) / 2;
  return kBcfParseOk;
}

BcfParseErr BcfConvertPhasedBiallelicDosage(const BcfImportContext* bicp, const GparseReadBcfMetadata* metap, unsigned char* record_start, uintptr_t* __restrict genovec, uintptr_t* __restrict phasepresent, uintptr_t* __restrict phaseinfo, uintptr_t* __restrict dosage_present, uintptr_t* __restrict dphase_present, Dosage** dosage_main_iter_ptr, SDosage** dphase_delta_iter_ptr) {
  // See VcfConvertPhasedBiallelicDosageLine().
  const BcfImportBaseContext* bibcp = &(bicp->bibc);
  const uint32_t sample_ct = bibcp->sample_ct;
  const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
  ZeroWArr(sample_ctl2, genovec);
  BcfParseErr bcf_parse_err = BcfApplyGqDpFilters(bibcp, metap, record_start, genovec);
  if (unlikely(bcf_parse_err != kBcfParseOk)) {
    return bcf_parse_err;
  }
  const unsigned char* hds_main = nullptr;
  uint32_t hds_value_ct = 0;
  if (metap->hds_vec_offset != UINT32_MAX) {
    hds_main = &(record_start[metap->hds_vec_offset * S_CAST(uintptr_t, kBytesPerVec)]);
    uint32_t hds_value_type;
    ScanBcfTypeAligned(&hds_main, &hds_value_type, &hds_value_ct);
    if (unlikely(hds_value_type != 5)) {
      return kBcfParseNonfloatDosage;
    }
    if (unlikely((hds_value_ct > 2) && bicp->bibc.error_on_polyploid)) {
      return kBcfParsePolyploidError;
    }
  }
  const float* dosage_main = nullptr;
  uint32_t dosage_value_ct = 0;
  if (metap->dosage_vec_offset != UINT32_MAX) {
    const unsigned char* dosage_main_raw = &(record_start[metap->dosage_vec_offset * S_CAST(uintptr_t, kBytesPerVec)]);
    uint32_t dosage_value_type;
    ScanBcfTypeAligned(&dosage_main_raw, &dosage_value_type, &dosage_value_ct);
    if (unlikely(dosage_value_type != 5)) {
      return kBcfParseNonfloatDosage;
    }
    dosage_main = R_CAST(const float*, dosage_main_raw);
  }
  const unsigned char* gt_main = nullptr;
  uint32_t gt_value_type = 0;
  uint32_t gt_value_ct = 0;
  if (metap->gt_vec_offset != UINT32_MAX) {
    gt_main = &(record_start[metap->gt_vec_offset * S_CAST(uintptr_t, kBytesPerVec)]);
    ScanBcfTypeAligned(&gt_main, &gt_value_type, &gt_value_ct);
    if (unlikely(gt_value_type > 1)) {
      return kBcfParseWideGt;
    }
    if (unlikely((gt_value_ct > 2) && bicp->bibc.error_on_polyploid)) {
      return kBcfParsePolyploidError;
    }
  }

  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  phasepresent[sample_ctl - 1] = 0;
  dosage_present[sample_ctl - 1] = 0;
  dphase_present[sample_ctl - 1] = 0;
  Halfword* phasepresent_alias = R_CAST(Halfword*, phasepresent);
  Halfword* phaseinfo_alias = R_CAST(Halfword*, phaseinfo);
  Halfword* dosage_present_alias = R_CAST(Halfword*, dosage_present);
  Halfword* dphase_present_alias = R_CAST(Halfword*, dphase_present);
  Dosage* dosage_main_iter = *dosage_main_iter_ptr;
  SDosage* dphase_delta_iter = *dphase_delta_iter_ptr;
  const unsigned char* biallelic_gt_lookup = bibcp->biallelic_gt_lookup;
  const uint32_t sample_ctl2_m1 = sample_ctl2 - 1;
  const uint32_t dosage_is_gp = bicp->dosage_is_gp;
  const uint32_t dosage_erase_halfdist = bicp->dosage_erase_halfdist;
  const uint32_t dphase_erase_halfdist = dosage_erase_halfdist + kDosage4th;
  const double import_dosage_certainty = bicp->import_dosage_certainty;
  uint32_t loop_len = kBitsPerWordD2;
  uint16_t gt_raw = 0;
  uint32_t is_haploid_or_0ploid = (gt_value_ct == 1);
  uintptr_t half_call_error_bit2 = 0;
  uint32_t sample_idx = 0;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= sample_ctl2_m1) {
      if (widx > sample_ctl2_m1) {
        break;
      }
      loop_len = ModNz(sample_ct, kBitsPerWordD2);
    }
    const uint32_t sample_idx_stop = sample_idx + loop_len;
    uintptr_t gq_dp_fail_word = genovec[widx];
    uintptr_t genovec_word = 0;
    uint32_t phasepresent_hw = 0;
    uint32_t phaseinfo_hw = 0;
    uint32_t dosage_present_hw = 0;
    uint32_t dphase_present_hw = 0;
    for (; sample_idx != sample_idx_stop; ++sample_idx, gq_dp_fail_word >>= 2) {
      if (gq_dp_fail_word & 1) {
        continue;
      }
      if (gt_value_ct) {
        memcpy(&gt_raw, &(gt_main[sample_idx * gt_value_ct]), 2);
        if (gt_value_ct > 1) {
          is_haploid_or_0ploid = gt_raw >> 15;
        }
      }
      const uint32_t sample_idx_lowbits = sample_idx % kBitsPerWordD2;
      const uint32_t shifted_bit = 1U << sample_idx_lowbits;
      DosageParseResult dpr = kDosageParseOk;
      uintptr_t cur_geno = 3;
      int32_t cur_dphase_delta = 0;
      uint32_t hds_valid = 0;
      uint32_t dosage_int;
      if (!ParseBcfBiallelicHds(dosage_main, hds_main, dosage_value_ct, hds_value_ct, sample_idx, is_haploid_or_0ploid, dosage_is_gp, import_dosage_certainty, &dpr, &dosage_int, &cur_dphase_delta, &hds_valid)) {
        if (hds_valid) {
          const uint32_t dphase_halfdist1 = DphaseHalfdist(dosage_int + cur_dphase_delta);
          const uint32_t dphase_halfdist2 = DphaseHalfdist(dosage_int - cur_dphase_delta);
          if ((dphase_halfdist1 < dphase_erase_halfdist) || (dphase_halfdist2 < dphase_erase_halfdist)) {
            // Usually no need to fill cur_geno here, since it'll get corrected
            // by --hard-call-threshold.
            if (cur_dphase_delta) {
              dphase_present_hw |= shifted_bit;
              *dphase_delta_iter++ = cur_dphase_delta;
            } else if (dosage_int == kDosageMid) {
              // quasi-bugfix (18 Jun 2023): if dosage_int == kDosageMid and
              // cur_dphase_delta == 0, we shouldn't save an explicit dosage at
              // all; this is just an ordinary unphased het.
              cur_geno = 1;
              goto BcfConvertPhasedBiallelicDosage_geno_done;
            }
            dosage_present_hw |= shifted_bit;
            *dosage_main_iter++ = dosage_int;
          } else {
            // Not saving dosage, since it's too close to an integer
            // (--dosage-erase-threshold).  Just directly synthesize the
            // hardcall we need.
            cur_geno = (dosage_int + kDosage4th) / kDosageMid;
            if (cur_geno == 1) {
              // Since dphase_erase_halfdist >= 8193, dphase_halfdist1 and
              // dphase_halfdist2 are both in [0, 8191] or [24577, 32768], so
              // it's always appropraite to save hardcall-phase.
              phasepresent_hw |= shifted_bit;
              if (cur_dphase_delta > 0) {
                phaseinfo_hw |= shifted_bit;
              }
            }
          }
          goto BcfConvertPhasedBiallelicDosage_geno_done;
        }
        // defer handling of unphased dosage
      } else if (unlikely(!dpr)) {
        return kBcfParseInvalidDosage;
      } else if (dpr == kDosageParseForceMissing) {
        goto BcfConvertPhasedBiallelicDosage_geno_done;
      }
      if (gt_main) {
        if (is_haploid_or_0ploid) {
          // gt_first == 0 -> missing (cur_geno == 3)
          // gt_first == 2 -> hom-ref (cur_geno == 0)
          // gt_first == 4 -> hom-alt (cur_geno == 2)
          cur_geno = (gt_raw & 0xff) - 2;
          // deliberate underflow
          if (cur_geno > 3) {
            cur_geno = 3;
          }
        } else if ((gt_value_ct == 2) || (gt_main[sample_idx * gt_value_ct + 2] == 0x81)) {
          // diploid
          const uint32_t bit_8 = gt_raw & 0x100;
          uint32_t first_allele_idx_p1 = (gt_raw & 0xff) >> 1;
          if (first_allele_idx_p1 > 3) {
            first_allele_idx_p1 = 3;
          }
          uint32_t second_allele_idx_p1 = gt_raw >> 9;
          if (second_allele_idx_p1 > 3) {
            second_allele_idx_p1 = 3;
          }
          cur_geno = biallelic_gt_lookup[first_allele_idx_p1 + second_allele_idx_p1 * 4];
          half_call_error_bit2 |= cur_geno;
          if ((cur_geno == 1) && bit_8) {
            phasepresent_hw |= shifted_bit;
            phaseinfo_hw |= (second_allele_idx_p1 & 1) << sample_idx_lowbits;
          }
        }
      }
      if (!dpr) {
        // now actually handle the unphased dosage
        const uint32_t cur_halfdist = BiallelicDosageHalfdist(dosage_int);
        if (cur_halfdist < dosage_erase_halfdist) {
          // ok for cur_geno to be 'wrong' for now, since it'll get corrected
          // by --hard-call-threshold
          dosage_present_hw |= shifted_bit;
          *dosage_main_iter++ = dosage_int;
        } else {
          // Not saving dosage, since it's too close to an integer, except
          // possibly in the implicit-phased-dosage edge case.
          // If that integer actually conflicts with the hardcall, we must
          // override the hardcall.
          cur_geno = (dosage_int + kDosage4th) / kDosageMid;
          if (phasepresent_hw & shifted_bit) {
            if (cur_geno != 1) {
              // Hardcall-phase no longer applies.
              phasepresent_hw ^= shifted_bit;
            } else if (cur_halfdist * 2 < dphase_erase_halfdist) {
              // Implicit phased-dosage, e.g. 0|0.99.  More stringent
              // dosage_erase_halfdist applies.
              dosage_present_hw |= shifted_bit;
              *dosage_main_iter++ = dosage_int;
            }
          }
        }
      }
    BcfConvertPhasedBiallelicDosage_geno_done:
      genovec_word |= cur_geno << (2 * sample_idx_lowbits);
    }
    genovec[widx] = genovec_word | gq_dp_fail_word;
    phasepresent_alias[widx] = phasepresent_hw;
    phaseinfo_alias[widx] = phaseinfo_hw;
    dosage_present_alias[widx] = dosage_present_hw;
    dphase_present_alias[widx] = dphase_present_hw;
  }
  if (unlikely(half_call_error_bit2 & 4)) {
    return kBcfParseHalfCallError;
  }
  *dosage_main_iter_ptr = dosage_main_iter;
  *dphase_delta_iter_ptr = dphase_delta_iter;
  return kBcfParseOk;
}

// 1: int8
// 2: int16
// 3: int32
// 5: float
// 7: char
static const unsigned char kBcfBytesPerElem[16] = {0, 1, 2, 4, 0, 4, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0};

typedef struct BcfGenoToPgenCtxStruct {
  BcfImportContext bic;
  uint32_t hard_call_halfdist;

  unsigned char** thread_wkspaces;

  uint32_t* thread_bidxs[2];
  GparseRecord* gparse[2];
  const uintptr_t* block_allele_idx_offsets[2];

  // PglErr set by main thread
  BcfParseErr* bcf_parse_errs;
  uintptr_t* err_vrec_idxs;
  uint32_t parse_failed;
} BcfGenoToPgenCtx;

THREAD_FUNC_DECL BcfGenoToPgenThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  BcfGenoToPgenCtx* ctx = S_CAST(BcfGenoToPgenCtx*, arg->sharedp->context);

  const BcfImportContext* bicp = &(ctx->bic);
  const BcfImportBaseContext* bibcp = &(bicp->bibc);
  const uint32_t sample_ct = bibcp->sample_ct;
  const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  const uint32_t hard_call_halfdist = ctx->hard_call_halfdist;
  unsigned char* thread_wkspace = ctx->thread_wkspaces[tidx];
  uintptr_t* patch_01_set = nullptr;
  AlleleCode* patch_01_vals = nullptr;
  uintptr_t* patch_10_set = nullptr;
  AlleleCode* patch_10_vals = nullptr;
  uintptr_t* phasepresent = nullptr;
  uintptr_t* phaseinfo = nullptr;
  uintptr_t* dosage_present = nullptr;
  Dosage* dosage_main = nullptr;
  uintptr_t* dphase_present = nullptr;
  SDosage* dphase_delta = nullptr;
  SDosage* tmp_dphase_delta = R_CAST(SDosage*, thread_wkspace);
  thread_wkspace = &(thread_wkspace[RoundUpPow2(sample_ct * sizeof(SDosage), kBytesPerVec)]);
  uintptr_t* write_patch_01_set = nullptr;
  AlleleCode* write_patch_01_vals = nullptr;
  uintptr_t* write_patch_10_set = nullptr;
  AlleleCode* write_patch_10_vals = nullptr;
  uintptr_t* write_phasepresent = nullptr;
  uintptr_t* write_phaseinfo = nullptr;
  uintptr_t* write_dosage_present = nullptr;
  Dosage* write_dosage_main = nullptr;
  uintptr_t* write_dphase_present = nullptr;
  SDosage* write_dphase_delta = nullptr;
  uint32_t cur_allele_ct = 2;
  uint32_t parity = 0;
  BcfParseErr bcf_parse_err = kBcfParseOk;
  uintptr_t vrec_idx = 0;
  do {
    const uintptr_t* block_allele_idx_offsets = ctx->block_allele_idx_offsets[parity];
    const uint32_t bidx_end = ctx->thread_bidxs[parity][tidx + 1];
    GparseRecord* cur_gparse = ctx->gparse[parity];

    for (uint32_t bidx = ctx->thread_bidxs[parity][tidx]; bidx != bidx_end; ++bidx) {
      GparseRecord* grp = &(cur_gparse[bidx]);
      uint32_t patch_01_ct = 0;
      uint32_t patch_10_ct = 0;
      uint32_t cur_phasepresent_exists = 0;
      uint32_t dosage_ct = 0;
      uint32_t dphase_ct = 0;
      unsigned char* record_start = grp->record_start;
      GparseFlags gparse_flags = grp->flags;
      if (gparse_flags == kfGparseNull) {
        SetAllBits(2 * sample_ct, R_CAST(uintptr_t*, record_start));
      } else {
        if (block_allele_idx_offsets) {
          cur_allele_ct = block_allele_idx_offsets[bidx + 1] - block_allele_idx_offsets[bidx];
        }
        uintptr_t* genovec = GparseGetPointers(thread_wkspace, sample_ct, cur_allele_ct, gparse_flags, &patch_01_set, &patch_01_vals, &patch_10_set, &patch_10_vals, &phasepresent, &phaseinfo, &dosage_present, &dosage_main, &dphase_present, &dphase_delta);
        uintptr_t* write_genovec = GparseGetPointers(record_start, sample_ct, cur_allele_ct, gparse_flags, &write_patch_01_set, &write_patch_01_vals, &write_patch_10_set, &write_patch_10_vals, &write_phasepresent, &write_phaseinfo, &write_dosage_present, &write_dosage_main, &write_dphase_present, &write_dphase_delta);
        const GparseReadBcfMetadata* metap = &(grp->metadata.read_bcf);
        if ((metap->hds_vec_offset == UINT32_MAX) && (metap->dosage_vec_offset == UINT32_MAX)) {
          if (!(gparse_flags & kfGparseHphase)) {
            if (cur_allele_ct == 2) {
              bcf_parse_err = BcfConvertUnphasedBiallelic(bibcp, metap, record_start, genovec);
            } else {
              bcf_parse_err = BcfConvertUnphasedMultiallelic(bibcp, metap, record_start, cur_allele_ct, &patch_01_ct, &patch_10_ct, genovec, patch_01_set, patch_01_vals, patch_10_set, patch_10_vals);
            }
          } else {
            if (cur_allele_ct == 2) {
              bcf_parse_err = BcfConvertPhasedBiallelic(bibcp, metap, record_start, genovec, phasepresent, phaseinfo);
            } else {
              bcf_parse_err = BcfConvertPhasedMultiallelic(bibcp, metap, record_start, cur_allele_ct, &patch_01_ct, &patch_10_ct, genovec, patch_01_set, patch_01_vals, patch_10_set, patch_10_vals, phasepresent, phaseinfo);
            }
            cur_phasepresent_exists = !AllWordsAreZero(phasepresent, sample_ctl);
          }
          if (unlikely(bcf_parse_err)) {
            vrec_idx = metap->rec_idx;
            goto BcfGenoToPgenThread_malformed;
          }
        } else {
          Dosage* dosage_main_iter = dosage_main;
          SDosage* dphase_delta_iter = dphase_delta;
          if (cur_allele_ct == 2) {
            bcf_parse_err = BcfConvertPhasedBiallelicDosage(bicp, metap, record_start, genovec, phasepresent, phaseinfo, dosage_present, dphase_present, &dosage_main_iter, &dphase_delta_iter);
          } else {
            // multiallelic dosage: shouldn't be possible to get here yet
            exit(S_CAST(int32_t, kPglRetInternalError));
          }
          if (unlikely(bcf_parse_err)) {
            vrec_idx = metap->rec_idx;
            goto BcfGenoToPgenThread_malformed;
          }
          dosage_ct = dosage_main_iter - dosage_main;
          if (dosage_ct) {
            dphase_ct = ApplyHardCallThreshPhased(dosage_present, dosage_main, dosage_ct, hard_call_halfdist, genovec, phasepresent, phaseinfo, dphase_present, dphase_delta, tmp_dphase_delta);
            memcpy(write_dosage_present, dosage_present, sample_ctl * sizeof(intptr_t));
            memcpy(write_dosage_main, dosage_main, dosage_ct * sizeof(Dosage));
            if (dphase_ct) {
              memcpy(write_dphase_present, dphase_present, sample_ctl * sizeof(intptr_t));
              memcpy(write_dphase_delta, dphase_delta, dphase_ct * sizeof(SDosage));
            }
          }
          cur_phasepresent_exists = !AllWordsAreZero(phasepresent, sample_ctl);
        }
        memcpy(write_genovec, genovec, sample_ctl2 * sizeof(intptr_t));
        if (patch_01_ct) {
          memcpy(write_patch_01_set, patch_01_set, sample_ctl * sizeof(intptr_t));
          memcpy(write_patch_01_vals, patch_01_vals, patch_01_ct * sizeof(AlleleCode));
        }
        if (patch_10_ct) {
          memcpy(write_patch_10_set, patch_10_set, sample_ctl * sizeof(intptr_t));
          memcpy(write_patch_10_vals, patch_10_vals, patch_10_ct * sizeof(AlleleCode) * 2);
        }
        if (cur_phasepresent_exists || dphase_ct) {
          memcpy(write_phasepresent, phasepresent, sample_ctl * sizeof(intptr_t));
          memcpy(write_phaseinfo, phaseinfo, sample_ctl * sizeof(intptr_t));
        }
      }

      grp->metadata.write.patch_01_ct = patch_01_ct;
      grp->metadata.write.patch_10_ct = patch_10_ct;
      grp->metadata.write.phasepresent_exists = cur_phasepresent_exists;
      grp->metadata.write.dosage_ct = dosage_ct;
      grp->metadata.write.multiallelic_dosage_ct = 0;
      grp->metadata.write.dphase_ct = dphase_ct;
      grp->metadata.write.multiallelic_dphase_ct = 0;
    }
    while (0) {
    BcfGenoToPgenThread_malformed:
      ctx->bcf_parse_errs[tidx] = bcf_parse_err;
      ctx->err_vrec_idxs[tidx] = vrec_idx;
      ctx->parse_failed = 1;
      break;
    }
    parity = 1 - parity;
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

// pgen_generated and psam_generated assumed to be initialized to 1.
PglErr BcfToPgen(const char* bcfname, const char* preexisting_psamname, const char* const_fid, const char* dosage_import_field, MiscFlags misc_flags, ImportFlags import_flags, LoadFilterLogFlags load_filter_log_import_flags, uint32_t no_samples_ok, uint32_t is_update_or_impute_sex, uint32_t is_splitpar, uint32_t is_sortvars, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, double import_dosage_certainty, char id_delim, char idspace_to, int32_t vcf_min_gq, int32_t vcf_min_dp, int32_t vcf_max_dp, VcfHalfCall halfcall_mode, FamCol fam_cols, uint32_t import_max_allele_ct, uint32_t max_thread_ct, char* outname, char* outname_end, ChrInfo* cip, uint32_t* pgen_generated_ptr, uint32_t* psam_generated_ptr) {
  // Yes, lots of this is copied-and-pasted from VcfToPgen(), but there are
  // enough differences that I don't think trying to handle them with the same
  // function is wise.

  // Possible todo: make this take proper advantage of an index file when a
  // chromosome filter has been specified.  (This requires an upgrade to
  // include/plink2_bgzf.)
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  FILE* bcffile = nullptr;
  const char* bgzf_errmsg = nullptr;
  char* pvar_cswritep = nullptr;
  uintptr_t vrec_idx = 0;
  const uint32_t half_call_explicit_error = (halfcall_mode == kVcfHalfCallError);
  PglErr reterr = kPglRetSuccess;
  BcfParseErr bcf_parse_err = kBcfParseOk;
  CompressStreamState pvar_css;
  PreinitCstream(&pvar_css);
  ThreadGroup tg;
  PreinitThreads(&tg);
  BcfGenoToPgenCtx ctx;
  BgzfRawMtDecompressStream bgzf;
  PreinitBgzfRawMtStream(&bgzf);
  STPgenWriter spgw;
  PreinitSpgw(&spgw);
  {
    // See TextFileOpenInternal().
    bcffile = fopen(bcfname, FOPEN_RB);
    if (unlikely(!bcffile)) {
      const uint32_t slen = strlen(bcfname);
      if (!StrEndsWith(bcfname, ".bcf", slen)) {
        logerrprintfww("Error: Failed to open %s : %s. (--bcf expects a complete filename; did you forget '.bcf' at the end?)\n", bcfname, strerror(errno));
      } else {
        logerrprintfww(kErrprintfFopen, bcfname, strerror(errno));
      }
      goto BcfToPgen_ret_OPEN_FAIL;
    }
    const uint32_t decompress_thread_ct = ClipU32(max_thread_ct - 1, 1, 4);
    uint32_t header_size;
    {
      char bgzf_header[16];
      uint32_t nbytes = fread_unlocked(bgzf_header, 1, 16, bcffile);
      if (unlikely(ferror_unlocked(bcffile))) {
        reterr = kPglRetReadFail;
        goto BcfToPgen_ret_BGZF_FAIL;
      }
      if (unlikely((nbytes != 16) || (!IsBgzfHeader(bgzf_header)))) {
        snprintf(g_logbuf, kLogbufSize, "Error: %s is not a BCF2 file.\n", bcfname);
        goto BcfToPgen_ret_MALFORMED_INPUT_WW;
      }
      reterr = BgzfRawMtStreamInit(bgzf_header, decompress_thread_ct, bcffile, nullptr, &bgzf, &bgzf_errmsg);
      if (unlikely(reterr)) {
        goto BcfToPgen_ret_BGZF_FAIL;
      }
      unsigned char prefix_buf[9];
      unsigned char* prefix_end = &(prefix_buf[9]);
      unsigned char* dummy = prefix_buf;
      reterr = BgzfRawMtStreamRead(prefix_end, &bgzf, &dummy, &bgzf_errmsg);
      if (unlikely(reterr)) {
        goto BcfToPgen_ret_BGZF_FAIL;
      }
      if (unlikely((dummy != prefix_end) || (!memequal_sk(prefix_buf, "BCF")))) {
        snprintf(g_logbuf, kLogbufSize, "Error: %s is not a BCF2 file.\n", bcfname);
        goto BcfToPgen_ret_MALFORMED_INPUT_WW;
      }
      if (unlikely(prefix_buf[3] != 2)) {
        if (prefix_buf[3] == 4) {
          snprintf(g_logbuf, kLogbufSize, "Error: %s appears to be a BCF1 file; --bcf only supports BCF2. Use 'bcftools view' to convert it to a PLINK-readable BCF.\n", bcfname);
        } else {
          snprintf(g_logbuf, kLogbufSize, "Error: %s is not a BCF2 file.\n", bcfname);
        }
        goto BcfToPgen_ret_MALFORMED_INPUT_WW;
      }
      if (unlikely(prefix_buf[4] > 2)) {
        snprintf(g_logbuf, kLogbufSize, "Error: %s appears to be formatted as BCFv2.%u; this " PROG_NAME_STR " build only supports v2.0-2.2. You may need to obtain an updated version of PLINK.\n", bcfname, prefix_buf[4]);
        goto BcfToPgen_ret_MALFORMED_INPUT_WW;
      }
      memcpy(&header_size, &(prefix_buf[5]), sizeof(int32_t));
      if (unlikely(header_size < 59)) {
        goto BcfToPgen_ret_MALFORMED_INPUT_GENERIC;
      }
    }
    // can't be const due to how --vcf-idspace-to is implemented
    // we could free this before the main conversion loop, but unlikely to be
    // more than a few MB so I won't bother
    // + 9 to simplify second pass
    char* vcf_header;
    {
      char* vcf_header_alloc;
      if (unlikely(bigstack_alloc_c(header_size + (9 * k1LU), &vcf_header_alloc))) {
        goto BcfToPgen_ret_NOMEM;
      }
      vcf_header = &(vcf_header_alloc[9]);
      unsigned char* header_load_iter = R_CAST(unsigned char*, vcf_header);
      unsigned char* header_load_end = &(header_load_iter[header_size]);
      reterr = BgzfRawMtStreamRead(header_load_end, &bgzf, &header_load_iter, &bgzf_errmsg);
      if (unlikely(reterr)) {
        goto BcfToPgen_ret_BGZF_FAIL;
      }
      if (unlikely(header_load_iter != header_load_end)) {
        goto BcfToPgen_ret_MALFORMED_INPUT_GENERIC;
      }
    }
    if (unlikely(!strequal_k_unsafe(&(vcf_header[header_size - 2]), "\n"))) {
      goto BcfToPgen_ret_MALFORMED_TEXT_HEADER;
    }
    char* vcf_header_end = &(vcf_header[header_size - 1]);
    const uint32_t prohibit_extra_chr = (misc_flags / kfMiscProhibitExtraChr) & 1;
    uint32_t dosage_import_field_slen = 0;

    uint32_t format_hds_search = 0;

    uint32_t unforced_gp = 0;
    if (dosage_import_field) {
      dosage_import_field_slen = strlen(dosage_import_field);
      ctx.bic.dosage_is_gp = 0;
      if (strequal_k(dosage_import_field, "HDS", dosage_import_field_slen)) {
        format_hds_search = 1;
        // special case: search for DS and HDS
        dosage_import_field = &(dosage_import_field[1]);
        dosage_import_field_slen = 2;
      } else if (strequal_k(dosage_import_field, "GP-force", dosage_import_field_slen)) {
        dosage_import_field = kGpText;
        dosage_import_field_slen = 2;
        ctx.bic.dosage_is_gp = 1;
      } else if (strequal_k(dosage_import_field, "GP", dosage_import_field_slen)) {
        unforced_gp = (import_dosage_certainty == 0.0);
        ctx.bic.dosage_is_gp = 1;
      }
    }

    ctx.bic.bibc.qual_mins[0] = vcf_min_gq;
    ctx.bic.bibc.qual_mins[1] = (vcf_min_dp == -1)? 0 : vcf_min_dp;
    ctx.bic.bibc.qual_maxs[0] = 0x7fffffff;
    ctx.bic.bibc.qual_maxs[1] = vcf_max_dp;
    // sample_ct set later
    if (halfcall_mode == kVcfHalfCallDefault) {
      halfcall_mode = kVcfHalfCallError;
    }
    ctx.bic.bibc.halfcall_mode = halfcall_mode;
    ctx.bic.bibc.error_on_polyploid = !(import_flags & kfImportPolyploidMissing);
    // Construct main biallelic GT-parsing lookup table.
    // Low 2 bits of index = first_allele_idx_p1 (or 3 for END_OF_VECTOR),
    //   high 2 bits = second_allele_idx_p1.
    // Table value is 2-bit genotype value, or 4 for half-call error.
    {
      unsigned char* biallelic_gt_lookup = ctx.bic.bibc.biallelic_gt_lookup;
      biallelic_gt_lookup[0] = 3;
      biallelic_gt_lookup[1] = 4;  // half-call
      biallelic_gt_lookup[2] = 4;
      biallelic_gt_lookup[3] = 3;

      biallelic_gt_lookup[5] = 0;
      biallelic_gt_lookup[6] = 1;
      // ignore second value if first value is end-of-vector
      biallelic_gt_lookup[7] = 3;

      biallelic_gt_lookup[10] = 2;
      biallelic_gt_lookup[11] = 3;

      biallelic_gt_lookup[12] = 3;
      biallelic_gt_lookup[13] = 0;
      biallelic_gt_lookup[14] = 2;
      biallelic_gt_lookup[15] = 3;

      if (halfcall_mode == kVcfHalfCallReference) {
        biallelic_gt_lookup[1] = 0;
        biallelic_gt_lookup[2] = 1;
      } else if (halfcall_mode == kVcfHalfCallHaploid) {
        biallelic_gt_lookup[1] = 0;
        biallelic_gt_lookup[2] = 2;
      } else if (halfcall_mode == kVcfHalfCallMissing) {
        biallelic_gt_lookup[1] = 3;
        biallelic_gt_lookup[2] = 3;
      }

      biallelic_gt_lookup[4] = biallelic_gt_lookup[1];
      biallelic_gt_lookup[8] = biallelic_gt_lookup[2];
      biallelic_gt_lookup[9] = biallelic_gt_lookup[6];
#ifdef USE_AVX2
      memcpy(&(biallelic_gt_lookup[16]), biallelic_gt_lookup, 16);
#endif
    }

    // always positive
    ctx.bic.dosage_erase_halfdist = kDosage4th - dosage_erase_thresh;

    ctx.bic.import_dosage_certainty = import_dosage_certainty;

    // Scan the text header block.
    // We perform the same validation as --vcf.  In addition, we need to
    // construct the contig and generic-string dictionaries; the latter may
    // require awareness of IDX keys.
    // (There is also a comment in the specification about explicitly
    // specifying a dictionary with a ##dictionary line, but AFAICT that is
    // left over from early brainstorming, is not consistent with the rest of
    // the design, and is not supported by bcftools so I assume it's a spec
    // error.)
    // hrec_add_idx() in htslib's vcf.c always puts IDX at the end of the
    // header line, so I'll be lazy and print a not-yet-supported error message
    // if it's positioned elsewhere for now.
    uint32_t chrset_exists = 0;
    uint32_t contig_string_idx_end = 0;
    uint32_t fif_string_idx_end = 1;  // fif = filter, info, format; initialized with FILTER/PASS
    uint32_t explicit_idx_keys = 2;  // 2 = unknown status
    uint32_t header_line_idx = 1;  // uint32 ok since uncompressed size < 2^32
    char* line_iter = vcf_header;
    for (; ; ++header_line_idx) {
      if (unlikely(line_iter == vcf_header_end)) {
        logerrputs("Error: No #CHROM header line in BCF text header block.\n");
        goto BcfToPgen_ret_MALFORMED_INPUT;
      }
      if (unlikely(*line_iter != '#')) {
        snprintf(g_logbuf, kLogbufSize, "Error: Line %u in BCF text header block does not start with '#'.\n", header_line_idx);
        goto BcfToPgen_ret_MALFORMED_INPUT_WW;
      }
      if (line_iter[1] != '#') {
        break;
      }
      char* line_main = &(line_iter[2]);
      char* line_end = AdvPastDelim(line_main, '\n');
      line_iter = line_end;
      // Recognized header lines:
      // ##fileformat: discard (regenerate; todo: conditionally error out)
      // ##fileDate: discard (regenerate)
      // ##source: discard (regenerate)
      // ##contig: ID added to contig dictionary, conditionally keep
      // ##FILTER: ID added to dictionary, otherwise passed through unchanged
      // ##INFO: note presence of INFO/PR, note presence of at least one non-PR
      //         field, keep data (though, if INFO/PR is the *only* field,
      //         omit it from the .pvar for consistency with --make-pgen
      //         default); ID added to dictionary
      //         update (8 Sep 2017): nonflag INFO/PR is noted, and not treated
      //         specially unless provisional-reference INFO/PR output would
      //         conflict with it
      // ##FORMAT: note presence of FORMAT/GT and FORMAT/GP, discard
      //           (regenerate); ID added to dictionary
      // ##chrSet: if recognized, perform consistency check and/or update
      //           chr_info
      //
      // Everything else (##reference, etc.) is passed through
      // unchanged.
      //
      // Because of how ##contig is handled (we only keep the lines which
      // correspond to chromosomes/contigs actually present in the BCF, and not
      // filtered out), we wait until the second pass to write the .pvar.
      if (StrStartsWithUnsafe(line_main, "chrSet=<")) {
        if (unlikely(chrset_exists)) {
          logerrputs("Error: Multiple ##chrSet header lines in BCF text header block.\n");
          goto BcfToPgen_ret_MALFORMED_INPUT;
        }
        chrset_exists = 1;
        // .pvar loader will print a warning if necessary
        reterr = ReadChrsetHeaderLine(&(line_main[8]), "--bcf file", misc_flags, header_line_idx, cip);
        if (unlikely(reterr)) {
          goto BcfToPgen_ret_1;
        }
        continue;
      }
      const uint32_t is_contig_line = StrStartsWithUnsafe(line_main, "contig=<");
      const uint32_t is_filter_line = StrStartsWithUnsafe(line_main, "FILTER=<");
      const uint32_t is_info_line = StrStartsWithUnsafe(line_main, "INFO=<");
      const uint32_t is_format_line = StrStartsWithUnsafe(line_main, "FORMAT=<");
      if (!(is_contig_line || is_filter_line || is_info_line || is_format_line)) {
        continue;
      }
      char* line_last_iter = &(line_end[-2]);
      // tolerate trailing spaces, not just \r
      while (ctou32(*line_last_iter) <= 32) {
        --line_last_iter;
      }
      if (unlikely(*line_last_iter != '>')) {
        snprintf(g_logbuf, kLogbufSize, "Error: Line %u in BCF text header block is malformed.\n", header_line_idx);
        goto BcfToPgen_ret_MALFORMED_INPUT_2;
      }
      --line_last_iter;
      // Assuming for now that the IDX key is at the end when present, it's
      // present iff:
      // - Last character inside <> is a digit.
      // - Last nondigit characters inside <> are ",IDX=".
      uint32_t cur_explicit_idx_key = 0;
      if (IsDigit(*line_last_iter)) {
        do {
          --line_last_iter;
        } while (IsDigit(*line_last_iter));
        cur_explicit_idx_key = StrStartsWithUnsafe(&(line_last_iter[-4]), ",IDX=");
      }
      // Points to leading '<' character.
      char* hkvline_iter = &(line_main[7 - 2 * is_info_line]);

      if (explicit_idx_keys == 2) {
        explicit_idx_keys = cur_explicit_idx_key;
        if (!cur_explicit_idx_key) {
          // The BCF specification does not require IDX= to be at the end of
          // the line, but it does require it to be present on all
          // contig/FILTER/INFO/FORMAT lines if it's present on any.
          // Thus, if a conforming writer puts IDX= in the middle of any lines,
          // we'll detect that either here, or in the next cur_explicit_idx_key
          // != explicit_idx_keys check.
          reterr = BcfHeaderLineIdxCheck(hkvline_iter, header_line_idx);
          if (unlikely(reterr)) {
            goto BcfToPgen_ret_1;
          }
        }
      } else {
        if (unlikely(cur_explicit_idx_key != explicit_idx_keys)) {
          reterr = BcfHeaderLineIdxCheck(hkvline_iter, header_line_idx);
          if (reterr != kPglRetSuccess) {
            goto BcfToPgen_ret_1;
          }
          logerrputs("Error: Some contig/FILTER/INFO/FORMAT line(s) in the BCF text header block have\nIDX keys, and others don't; this is not permitted by the specification.\n");
          goto BcfToPgen_ret_MALFORMED_INPUT;
        }
      }
      // Now points past leading '<' character.
      ++hkvline_iter;

      char* id_ptr;
      uint32_t id_slen;
      if (unlikely(HkvlineId(&hkvline_iter, &id_ptr, &id_slen))) {
        snprintf(g_logbuf, kLogbufSize, "Error: Line %u in BCF text header block is malformed.\n", header_line_idx);
        goto BcfToPgen_ret_MALFORMED_INPUT_2;
      }
      // Special case: FILTER/PASS IDX=0 line is implicit, we must be
      // indifferent to whether it's actually present.  (Though when it is
      // present, it's subject to the same IDX= always/never-present rule as
      // the other header lines, which is why we don't perform this check
      // earlier.)
      if (is_filter_line && strequal_k(id_ptr, "PASS", id_slen)) {
        continue;
      }
      if (!cur_explicit_idx_key) {
        if (is_contig_line) {
          ++contig_string_idx_end;
        } else {
          ++fif_string_idx_end;
        }
      } else {
        uint32_t val;
        if (unlikely(ScanUintDefcap(&(line_last_iter[1]), &val))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Invalid IDX= value on line %u of BCF text header block.\n", header_line_idx);
          goto BcfToPgen_ret_MALFORMED_INPUT_2;
        }
        if (is_contig_line) {
          if (val >= contig_string_idx_end) {
            contig_string_idx_end = val + 1;
          }
        } else {
          // Easier to check for duplicates on the second pass.
          if (val >= fif_string_idx_end) {
            fif_string_idx_end = val + 1;
          }
        }
      }
    }

    FinalizeChrset(load_filter_log_import_flags, cip);
    // don't call FinalizeChrInfo here, since this may be followed by --pmerge,
    // etc.

    if (unlikely(!StrStartsWithUnsafe(line_iter, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"))) {
      snprintf(g_logbuf, kLogbufSize, "Error: Header line %u of BCF text header block does not have expected field sequence after #CHROM.\n", header_line_idx);
      goto BcfToPgen_ret_MALFORMED_INPUT_WW;
    }
    uint32_t sample_ct = 0;
    {
      char* linebuf_iter = &(line_iter[strlen("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")]);
      if (StrStartsWithUnsafe(linebuf_iter, "\tFORMAT\t")) {
        reterr = VcfSampleLine(preexisting_psamname, const_fid, misc_flags, import_flags, fam_cols, id_delim, idspace_to, 'b', &(linebuf_iter[strlen("\tFORMAT\t")]), outname, outname_end, &sample_ct);
        if (unlikely(reterr)) {
          goto BcfToPgen_ret_1;
        }
      }
    }
    if (unlikely((!sample_ct) && (!no_samples_ok))) {
      logerrputs("Error: No samples in BCF text header block.  (This is only permitted when you\nhaven't specified another operation which requires genotype or sample\ninformation.)\n");
      goto BcfToPgen_ret_DEGENERATE_DATA;
    }
    if (unlikely(sample_ct >= (1 << 24))) {
      snprintf(g_logbuf, kLogbufSize, "Error: BCF text header block has %u sample IDs, which is larger than the BCF limit of 2^24 - 1.\n", sample_ct);
      goto BcfToPgen_ret_MALFORMED_INPUT_WW;
    }
    ctx.bic.bibc.sample_ct = sample_ct;
    const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);

    // Rescan header to save FILTER/INFO/FORMAT and contig ID dictionaries.
    // (We could also print a warning when an INFO key type declaration is
    // inconsistent with any variant entry, but that's not essential for
    // decoding; let's get this working first.)
    uintptr_t* bcf_contig_keep;
    const char** contig_names;
    uint32_t* contig_slens;
    const char** contig_out_names = nullptr; // spurious g++ 4.8 warning
    uint32_t* contig_out_slens = nullptr;
    char* contig_out_buf;
    const char** fif_strings;
    uint32_t* fif_slens;
    uintptr_t* bcf_contig_seen;
    uintptr_t* sample_nypbuf = nullptr;
    if (unlikely(bigstack_calloc_w(BitCtToWordCt(contig_string_idx_end), &bcf_contig_keep) ||
                 bigstack_calloc_kcp(contig_string_idx_end, &contig_names) ||
                 bigstack_calloc_u32(contig_string_idx_end, &contig_slens) ||
                 bigstack_alloc_kcp(contig_string_idx_end, &contig_out_names) ||
                 bigstack_alloc_u32(contig_string_idx_end, &contig_out_slens) ||
                 // no chrtoa() result is longer than (3 + kMaxChrTextnumSlen)
                 // characters long
                 bigstack_alloc_c((kMaxChrTextnum + 1 + kChrOffsetCt) * (3 + kMaxChrTextnumSlen), &contig_out_buf) ||
                 bigstack_calloc_kcp(fif_string_idx_end, &fif_strings) ||
                 bigstack_calloc_u32(fif_string_idx_end, &fif_slens) ||
                 bigstack_calloc_w(BitCtToWordCt(contig_string_idx_end), &bcf_contig_seen) ||
                 bigstack_alloc_w(sample_ctl2, &sample_nypbuf))) {
      goto BcfToPgen_ret_NOMEM;
    }
    const uint32_t header_line_ct = header_line_idx;
    // sidx = string index
    uint32_t gt_sidx = 0;  // replaces format_gt_exists
    uint32_t gq_sidx = 0;  // replaces format_gq_relevant
    uint32_t dp_sidx = 0;  // replaces format_dp_relevant
    uint32_t dosage_sidx = 0;  // replaces format_dosage_relevant
    uint32_t hds_sidx = 0;  // replaces later uses of format_hds_search

    // replaces info_pr_exists... except that UINT32_MAX is the not-present
    // value, just in case there's an INFO/PASS field for some reason.
    uint32_t pr_sidx = UINT32_MAX;

    uint32_t info_pr_nonflag_exists = 0;
    uint32_t info_nonpr_exists = 0;

    // we don't want to clutter the second pass with too many instances of
    // fwrite_ck(), so we compute FILTER-column, INFO-column and allele length
    // bounds, and then size the write buffer accordingly.
    uint32_t other_slen_ubound = kMaxFloatGSlen;  // QUAL, alleles, etc.
    uint64_t filter_info_slen_ubound = 1;
    {
      unsigned char* tmp_alloc_end = bigstack_end_mark;
      if (StoreStringAtEndK(g_bigstack_base, "PASS", strlen("PASS"), &tmp_alloc_end, &(fif_strings[0]))) {
        goto BcfToPgen_ret_NOMEM;
      }
      fif_slens[0] = 4;
      line_iter = vcf_header;
      uint32_t contig_idx = 0;
      uint32_t fif_idx = 1;
      uint32_t cur_header_idx = 0;
      for (header_line_idx = 1; header_line_idx != header_line_ct; ++header_line_idx) {
        char* line_main = &(line_iter[2]);
        char* line_end = AdvPastDelim(line_main, '\n');
        line_iter = line_end;
        const uint32_t is_contig_line = StrStartsWithUnsafe(line_main, "contig=<");
        const uint32_t is_filter_line = StrStartsWithUnsafe(line_main, "FILTER=<");
        const uint32_t is_info_line = StrStartsWithUnsafe(line_main, "INFO=<");
        const uint32_t is_format_line = StrStartsWithUnsafe(line_main, "FORMAT=<");
        if (!(is_contig_line || is_filter_line || is_info_line || is_format_line)) {
          continue;
        }
        const char* hkvline_iter = &(line_main[8 - 2 * is_info_line]);
        const char* id_ptr;
        uint32_t id_slen;
        // previously validated
        HkvlineIdK(&hkvline_iter, &id_ptr, &id_slen);

        if (is_filter_line && strequal_k(id_ptr, "PASS", id_slen)) {
          continue;
        }
        if (explicit_idx_keys) {
          char* line_last_iter = &(line_end[-2]);
          while (ctou32(*line_last_iter) <= 32) {
            --line_last_iter;
          }
          --line_last_iter;
          while (IsDigit(*line_last_iter)) {
            --line_last_iter;
          }
          // previously validated
          ScanUintDefcap(&(line_last_iter[1]), &cur_header_idx);
        }
        const char** target;
        uint32_t* slen_dst;
        if (is_contig_line) {
          if (!explicit_idx_keys) {
            cur_header_idx = contig_idx++;
          }
          target = &(contig_names[cur_header_idx]);
          slen_dst = &(contig_slens[cur_header_idx]);
        } else {
          if (!explicit_idx_keys) {
            cur_header_idx = fif_idx++;
          }
          target = &(fif_strings[cur_header_idx]);
          slen_dst = &(fif_slens[cur_header_idx]);
        }
        // May enforce id_slen <= kMaxIdSlen here later.  (Technically only
        // necessary if we use the key in e.g. error messages, which happens on
        // BCF export but not import.)
        if (*target) {
          // Duplicate IDX values are actually permitted, but only when the ID
          // string is identical.
          // (In the no-IDX case, if e.g. a FILTER and INFO key are identical,
          // their implicit IDX values are still different.)
          if (unlikely((id_slen != *slen_dst) || (!memequal(id_ptr, *target, id_slen)))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Multiple %s IDs in BCF text header block have IDX=%u.\n", is_contig_line? "contig" : "FILTER/INFO/FORMAT", cur_header_idx);
            goto BcfToPgen_ret_MALFORMED_INPUT_WW;
          }
        } else {
          // technically only need to null-terminate contig names
          if (StoreStringAtEndK(g_bigstack_base, id_ptr, id_slen, &tmp_alloc_end, target)) {
            goto BcfToPgen_ret_NOMEM;
          }
          *slen_dst = id_slen;
          if (id_slen > other_slen_ubound) {
            other_slen_ubound = id_slen;
          }
        }
        if (is_contig_line || is_filter_line) {
          continue;
        }
        const char* numstr;
        const char* typestr;
        uint32_t num_slen;
        uint32_t type_slen;
        if (unlikely(HkvlineNumTypeK(hkvline_iter, &numstr, &num_slen, &typestr, &type_slen))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Line %u in BCF text header block is malformed.\n", header_line_idx);
          goto BcfToPgen_ret_MALFORMED_INPUT_2;
        }
        if (is_info_line) {
          const uintptr_t is_flag = strequal_k(typestr, "Flag", type_slen);
          if (strequal_k(id_ptr, "PR", id_slen)) {
            if (unlikely((pr_sidx != UINT32_MAX) || info_pr_nonflag_exists)) {
              logerrputs("Error: Duplicate INFO/PR line in BCF text header block.\n");
              goto BcfToPgen_ret_MALFORMED_INPUT;
            }
            if (is_flag) {
              pr_sidx = cur_header_idx;
            } else {
              info_pr_nonflag_exists = 1;
              logerrprintfww("Warning: Line %u of BCF text header block has an unexpected definition of INFO/PR. This interferes with a few merge and liftover operations.\n", header_line_idx);
            }
          } else {
            info_nonpr_exists = 1;
          }
          continue;
        }
        // FORMAT
        if (strequal_k(id_ptr, "GT", id_slen)) {
          if (unlikely(gt_sidx)) {
            logerrputs("Error: Duplicate FORMAT/GT line in BCF text header block.\n");
            goto BcfToPgen_ret_MALFORMED_INPUT;
          }
          if (unlikely((!strequal_k(numstr, "1", num_slen)) || (!strequal_k(typestr, "String", type_slen)))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Line %u of BCF text header block does not have expected FORMAT/GT format.\n", header_line_idx);
            goto BcfToPgen_ret_MALFORMED_INPUT_WW;
          }
          gt_sidx = cur_header_idx;
        } else if ((vcf_min_gq != -1) && strequal_k(id_ptr, "GQ", id_slen)) {
          if (unlikely(gq_sidx)) {
            logerrputs("Error: Duplicate FORMAT/GQ header line in BCF text header block.\n");
            goto BcfToPgen_ret_MALFORMED_INPUT;
          }
          if (unlikely(!strequal_k(numstr, "1", num_slen))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Line %u of BCF text header block does not have expected FORMAT/GQ format.\n", header_line_idx);
            goto BcfToPgen_ret_MALFORMED_INPUT_WW;
          }
          gq_sidx = cur_header_idx;
        } else if (((vcf_min_dp != -1) || (vcf_max_dp != 0x7fffffff)) && strequal_k(id_ptr, "DP", id_slen)) {
          if (unlikely(dp_sidx)) {
            logerrputs("Error: Duplicate FORMAT/DP header line in BCF text header block.\n");
            goto BcfToPgen_ret_MALFORMED_INPUT;
          }
          if (unlikely(!strequal_k(numstr, "1", num_slen))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Line %u of BCF text header block does not have expected FORMAT/DP format.\n", header_line_idx);
            goto BcfToPgen_ret_MALFORMED_INPUT_WW;
          }
          dp_sidx = cur_header_idx;
        } else if (dosage_import_field) {
          if ((id_slen == dosage_import_field_slen) && memequal(id_ptr, dosage_import_field, dosage_import_field_slen)) {
            if (unlikely(dosage_sidx)) {
              snprintf(g_logbuf, kLogbufSize, "Error: Duplicate FORMAT/%s header line in BCF text header block.\n", dosage_import_field);
              goto BcfToPgen_ret_MALFORMED_INPUT_WW;
            }
            dosage_sidx = cur_header_idx;
          } else if (format_hds_search && strequal_k(id_ptr, "HDS", id_slen)) {
            if (unlikely(hds_sidx)) {
              logerrputs("Error: Duplicate FORMAT/HDS header line in BCF text header block.\n");
              goto BcfToPgen_ret_MALFORMED_INPUT;
            }
            hds_sidx = cur_header_idx;
          } else if (unforced_gp && strequal_k(id_ptr, "DS", id_slen)) {
            logerrputs("Error: --bcf dosage=GP specified, but --import-dosage-certainty was not and\nFORMAT/DS header line is present.\nSince " PROG_NAME_STR " collapses genotype probabilities down to dosages (even when\nperforming simple operations like \"" PROG_NAME_STR " --bcf ... --export bgen-1.2 ...\"),\n'dosage=GP' almost never makes sense in this situation.  Either change it to\n'dosage=DS' (if dosages are good enough for your analysis), or use another\nprogram to work with the genotype probabilities.\nThere is one notable exception to the preceding recommendation: you are writing\na script to process BCF files that are guaranteed to have FORMAT/GP, but may or\nmay not have FORMAT/DS.  You can use 'dosage=GP-force' to suppress this error\nin that situation.\n");
            goto BcfToPgen_ret_INCONSISTENT_INPUT;
          }
        }
      }
      BigstackEndSet(tmp_alloc_end);
    }
    const uint32_t require_gt = (load_filter_log_import_flags / kfLoadFilterLogVcfRequireGt) & 1;
    if (unlikely((!gt_sidx) && require_gt)) {
      logerrputs("Error: No FORMAT/GT key found in BCF text header block, when --vcf-require-gt\nwas specified.\n(If this header line is actually present, but with extra spaces or unusual\nfield ordering, standardize the header with e.g. bcftools.)\n");
      goto BcfToPgen_ret_INCONSISTENT_INPUT;
    }
    if ((!gq_sidx) && (vcf_min_gq != -1)) {
      logerrputs("Warning: No FORMAT/GQ key found in BCF text header block.  --vcf-min-gq\nignored.\n(If this header line is actually present, but with extra spaces or unusual\nfield ordering, standardize the header with e.g. bcftools.)\n");
    }
    if ((!dp_sidx) && ((vcf_max_dp != 0x7fffffff) || (vcf_min_dp != -1))) {
      logerrputs("Warning: No FORMAT/DP key found in BCF text header block.  --vcf-{max,min}-dp\nignored.\n(If this header line is actually present, but with extra spaces or unusual\nfield ordering, standardize the header with e.g. bcftools.)\n");
    }
    if (format_hds_search) {
      if (!hds_sidx) {
        if (!dosage_sidx) {
          logerrputs("Warning: No FORMAT/DS or :HDS key found in BCF text header block.  Dosages will\nnot be imported.\n(If these header line(s) are actually present, but with extra spaces or unusual\nfield ordering, standardize the header with e.g. bcftools.)\n");
        } else {
          logerrputs("Warning: No FORMAT/HDS key found in BCF text header block.  Dosages will be\nimported (from FORMAT/DS), but phase information will be limited or absent.\n");
        }
      }
    } else if ((!dosage_sidx) && dosage_import_field) {
      logerrprintfww("Warning: No FORMAT/%s key found in BCF text header block. Dosages will not be imported. (If this header line is actually present, but with extra spaces or unusual field ordering, standardize the header with e.g. bcftools.)\n", dosage_import_field);
    }
    const uint32_t ref_n_missing = (import_flags / kfImportVcfRefNMissing) & 1;
    if (unlikely(ref_n_missing && (pr_sidx == UINT32_MAX))) {
      logerrputs("Error: --vcf-ref-n-missing was specified, but the BCF does not have the\nINFO/PR header line that should be present in any .ped-derived BCF.\n");
      goto BcfToPgen_ret_INCONSISTENT_INPUT;
    }

    memcpy(contig_out_names, contig_names, contig_string_idx_end * sizeof(intptr_t));
    memcpy(contig_out_slens, contig_slens, contig_string_idx_end * sizeof(int32_t));

    unsigned char* bigstack_end_mark2 = g_bigstack_end;
    uintptr_t loadbuf_size = RoundDownPow2(bigstack_left() / 2, kEndAllocAlign);
#ifdef __LP64__
    // l_shared and l_indiv are of type uint32_t, so we never need to look at
    // more than 8 GiB at a time.
    if (loadbuf_size > (k1LU << 33)) {
      loadbuf_size = k1LU << 33;
    }
#endif
    // guaranteeing at least this many bytes allocated at bigstack_end
    // simplifies a bit of parsing logic
    if (loadbuf_size < (1 << 18)) {
      goto BcfToPgen_ret_NOMEM;
    }
    // Placed at end of arena so we can shrink it before the second pass
    // without fragmenting memory.
    unsigned char* loadbuf = S_CAST(unsigned char*, bigstack_end_alloc_raw(loadbuf_size));

    uint32_t variant_ct = 0;
    uintptr_t max_variant_ct = bigstack_left() / sizeof(intptr_t);
    if (pr_sidx != UINT32_MAX) {
      // nonref_flags
      max_variant_ct -= BitCtToAlignedWordCt(max_variant_ct) * kWordsPerVec;
    }
#ifdef __LP64__
    if (max_variant_ct > kPglMaxVariantCt) {
      max_variant_ct = kPglMaxVariantCt;
    }
#endif
    // max(32, l_shared + l_indiv - 24)
    // (there are cheaper ways to handle gt_exists, but the memory savings are
    // unlikely to be significant.)
    uintptr_t loadbuf_size_needed = 32;

    // Sum of GT, DS/GP, HDS, GQ, and DP vector lengths.
    uint32_t max_observed_rec_vecs = 0;

    const uint32_t max_variant_ctaw = BitCtToAlignedWordCt(max_variant_ct);
    // don't need dosage_flags or dphase_flags; dosage overrides GT so slow
    // parse needed
    uintptr_t* nonref_flags = nullptr;
    if (pr_sidx != UINT32_MAX) {
      nonref_flags = S_CAST(uintptr_t*, bigstack_alloc_raw_rd(max_variant_ctaw * sizeof(intptr_t)));
    }
    char* contig_out_buf_iter = contig_out_buf;
    uintptr_t* nonref_flags_iter = nonref_flags;
    uintptr_t* allele_idx_offsets = R_CAST(uintptr_t*, g_bigstack_base);
    uintptr_t nonref_word = 0;
    uintptr_t allele_idx_end = 0;
    uint32_t max_allele_ct = 0;
    uint32_t phase_or_dosage_found = 0;
    const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
    uint32_t par_warn_bcf_chrom = UINT32_MAX;
    while (1) {
      ++vrec_idx;  // 1-based since it's only used in error messages
      unsigned char* loadbuf_read_iter = loadbuf;
      reterr = BgzfRawMtStreamRead(&(loadbuf[32]), &bgzf, &loadbuf_read_iter, &bgzf_errmsg);
      if (unlikely(reterr)) {
        goto BcfToPgen_ret_BGZF_FAIL_N;
      }
      if (&(loadbuf[32]) != loadbuf_read_iter) {
        if (likely(loadbuf_read_iter == loadbuf)) {
          // EOF
          // possible todo: verify empty block is present at end
          break;
        }
        goto BcfToPgen_ret_VREC_GENERIC;
      }
      const uint32_t* vrec_header = R_CAST(uint32_t*, loadbuf);
      // IMPORTANT: Official specification is wrong about the ordering of these
      // fields as of Feb 2020!!  The correct ordering can be inferred from
      // bcf_read1_core() in htslib vcf.c:
      //   [0]: l_shared
      //   [1]: l_indiv
      //   [2]: chrom
      //   [3]: pos
      //   [4]: rlen
      //   [5]: qual, NOT n_allele_info
      //   [6]: low 16 bits n_info, high 16 bits n_allele
      //   [7]: low 24 bits n_sample, high 8 bits n_fmt
      const uint32_t l_shared = vrec_header[0];
      const uint32_t l_indiv = vrec_header[1];
      const uint32_t chrom = vrec_header[2];
      // Can ignore pos and qual on this pass.  (QUAL/FILTER enforcement is now
      // handled by the .pvar loader.)  Always ignore rlen for now, though we
      // may want to add a consistency check later.
      uint32_t n_allele = vrec_header[6] >> 16;
      const uint32_t n_info = vrec_header[6] & 0xffff;
      const uint32_t n_sample = vrec_header[7] & 0xffffff;
      const uint32_t n_fmt = vrec_header[7] >> 24;

      if (unlikely((l_shared < 24) || (chrom >= contig_string_idx_end) || (n_sample != sample_ct))) {
        goto BcfToPgen_ret_VREC_GENERIC;
      }
      const uint64_t second_load_size = l_shared + S_CAST(uint64_t, l_indiv) - 24;
      if (unlikely(second_load_size > loadbuf_size)) {
        goto BcfToPgen_ret_NOMEM;
      }
      if (second_load_size > loadbuf_size_needed) {
        loadbuf_size_needed = second_load_size;
      }
      loadbuf_read_iter = loadbuf;
      unsigned char* indiv_end = &(loadbuf[second_load_size]);
      reterr = BgzfRawMtStreamRead(indiv_end, &bgzf, &loadbuf_read_iter, &bgzf_errmsg);
      if (unlikely(reterr)) {
        goto BcfToPgen_ret_BGZF_FAIL_N;
      }
      if (unlikely(indiv_end != loadbuf_read_iter)) {
        goto BcfToPgen_ret_VREC_GENERIC;
      }

      // We've now decompressed to the end of the variant record, and can
      // safely skip the variant with "continue;".  There are currently three
      // cases where we might want to do this: chromosome filters,
      // --import-max-alleles, and --vcf-require-gt.
      if (n_allele > import_max_allele_ct) {
        continue;
      }
      if (IsSet(bcf_contig_seen, chrom)) {
        if (!IsSet(bcf_contig_keep, chrom)) {
          continue;
        }
      } else if (!require_gt) {
        const uint32_t contig_slen = contig_slens[chrom];
        if (unlikely(!contig_slen)) {
          putc_unlocked('\n', stdout);
          logerrputs("Error: Empty contig ID in --bcf file.");
          goto BcfToPgen_ret_MALFORMED_INPUT;
        }
        // Don't want to mutate cip in require_gt case until we know the contig
        // is being kept.
        uint32_t cur_chr_code;
        reterr = GetOrAddChrCode(contig_names[chrom], "--bcf file", 0, contig_slen, prohibit_extra_chr, cip, &cur_chr_code);
        if (unlikely(reterr)) {
          goto BcfToPgen_ret_1;
        }
        SetBit(chrom, bcf_contig_seen);
        if (!IsSet(cip->chr_mask, cur_chr_code)) {
          continue;
        }
        SetBit(chrom, bcf_contig_keep);
        if (cur_chr_code <= cip->max_code) {
          char* rendered_chr_name = contig_out_buf_iter;
          contig_out_buf_iter = chrtoa(cip, cur_chr_code, rendered_chr_name);
          contig_out_names[chrom] = rendered_chr_name;
          contig_out_slens[chrom] = contig_out_buf_iter - rendered_chr_name;
          if (cur_chr_code == x_code) {
            par_warn_bcf_chrom = chrom;
          }
        }
      }
      const unsigned char* shared_end = indiv_end - l_indiv;
      const unsigned char* parse_iter = shared_end;

      const unsigned char* gt_start = nullptr;
      const unsigned char* qual_starts[2];
      qual_starts[0] = nullptr;
      qual_starts[1] = nullptr;
      const unsigned char* dosage_start = nullptr;
      const unsigned char* hds_start = nullptr;
      uint32_t cur_observed_rec_vecs = 0;
      for (uint32_t fmt_idx = 0; fmt_idx != n_fmt; ++fmt_idx) {
        // 1. typed int indicating which FORMAT field
        // 2. shared type descriptor for each entry (usually a single byte)
        // 3. sample_ct entries
        uint32_t sidx;
        if (unlikely(ScanBcfTypedInt(&parse_iter, &sidx) || (parse_iter > indiv_end) || (sidx >= fif_string_idx_end))) {
          goto BcfToPgen_ret_VREC_GENERIC;
        }
        const uint32_t key_slen = fif_slens[sidx];
        const unsigned char* type_start = parse_iter;
        uint32_t value_type;
        uint32_t value_ct;
        if (unlikely((!key_slen) || ScanBcfType(&parse_iter, &value_type, &value_ct) || (parse_iter > indiv_end))) {
          goto BcfToPgen_ret_VREC_GENERIC;
        }
        const unsigned char* vec_start = parse_iter;
        if (value_ct) {
          const uint32_t bytes_per_elem = kBcfBytesPerElem[value_type];
          const uint64_t vec_byte_ct = bytes_per_elem * S_CAST(uint64_t, value_ct) * sample_ct;
          if (unlikely((!bytes_per_elem) || (S_CAST(uint64_t, indiv_end - parse_iter) < vec_byte_ct))) {
            goto BcfToPgen_ret_VREC_GENERIC;
          }
          parse_iter = &(parse_iter[vec_byte_ct]);
        }
        // defend against FORMAT/PASS edge case
        if (!sidx) {
          continue;
        }
        if (sidx == gt_sidx) {
          gt_start = type_start;
        } else if (sidx == gq_sidx) {
          qual_starts[0] = type_start;
        } else if (sidx == dp_sidx) {
          qual_starts[1] = type_start;
        } else if (sidx == dosage_sidx) {
          dosage_start = type_start;
        } else if (sidx == hds_sidx) {
          hds_start = type_start;
        } else {
          continue;
        }
#ifdef __LP64__
        // bcf-type guaranteed to fit in <= 16 bytes
        ++cur_observed_rec_vecs;
#else
        cur_observed_rec_vecs += DivUp(S_CAST(uintptr_t, vec_start - type_start), kBytesPerVec);
#endif
        cur_observed_rec_vecs += DivUp(S_CAST(uintptr_t, parse_iter - vec_start), kBytesPerVec);
      }
      if (require_gt) {
        if (!gt_start) {
          continue;
        }
        if (!IsSet(bcf_contig_seen, chrom)) {
          uint32_t cur_chr_code;
          reterr = GetOrAddChrCode(contig_names[chrom], "--bcf file", 0, strlen(contig_names[chrom]), prohibit_extra_chr, cip, &cur_chr_code);
          if (unlikely(reterr)) {
            goto BcfToPgen_ret_1;
          }
          SetBit(chrom, bcf_contig_seen);
          if (!IsSet(cip->chr_mask, cur_chr_code)) {
            continue;
          }
          SetBit(chrom, bcf_contig_keep);
          if (cur_chr_code <= cip->max_code) {
            char* rendered_chr_name = contig_out_buf_iter;
            contig_out_buf_iter = chrtoa(cip, cur_chr_code, rendered_chr_name);
            contig_out_names[chrom] = rendered_chr_name;
            contig_out_slens[chrom] = contig_out_buf_iter - rendered_chr_name;
            if (cur_chr_code == x_code) {
              par_warn_bcf_chrom = chrom;
            }
          }
        }
      }
      // We finally know for sure that we need to convert this variant.

      // Remainder of l_shared:
      //   ID (typed string)
      //   REF+ALT (n_allele typed strings)
      //   FILTER (typed vector of integers)
      //   INFO (n_info (typed integer, typed vector) pairs)
      // We scan these fields to compute an upper bound on the necessary .pvar
      // write-buffer size.  In the unlikely pr_sidx != UINT32_MAX case, we
      // also need to scan for presence of INFO/PR.
      parse_iter = loadbuf;
      // ID, REF, ALT
      const uint32_t str_ignore_ct = n_allele + 1;
      for (uint32_t uii = 0; uii != str_ignore_ct; ++uii) {
        const char* str_start;
        uint32_t slen;
        if (unlikely(ScanBcfTypedString(shared_end, &parse_iter, &str_start, &slen))) {
          goto BcfToPgen_ret_VREC_GENERIC;
        }
        if (slen > other_slen_ubound) {
          other_slen_ubound = slen;
        }
      }
      // bugfix (22 Feb 2020): need to treat n_allele as 2 after this point
      // when it's 1
      if (n_allele == 1) {
        n_allele = 2;
      }
      // FILTER
      uint32_t value_type;
      uint32_t value_ct;
      if (unlikely(ScanBcfType(&parse_iter, &value_type, &value_ct) || (parse_iter > shared_end))) {
        goto BcfToPgen_ret_VREC_GENERIC;
      }
      if (value_ct) {
        const uint32_t value_type_m1 = value_type - 1;
        const uintptr_t vec_byte_ct = S_CAST(uintptr_t, value_ct) << value_type_m1;
        if (unlikely((value_type_m1 > 2) || (S_CAST(uintptr_t, shared_end - parse_iter) < vec_byte_ct))) {
          goto BcfToPgen_ret_VREC_GENERIC;
        }
        uint64_t cur_filter_slen = value_ct - 1;  // delimiters
        if (!value_type_m1) {
          for (uint32_t filter_idx = 0; filter_idx != value_ct; ++filter_idx) {
            const uint32_t cur_val = parse_iter[filter_idx];
            if (cur_val == 0x80) {
              // missing
              ++cur_filter_slen;
            } else {
              const uint32_t key_slen = fif_slens[cur_val];
              if (unlikely((cur_val >= fif_string_idx_end) || (cur_val > 0x80) || (!key_slen))) {
                goto BcfToPgen_ret_VREC_GENERIC;
              }
              cur_filter_slen += key_slen;
            }
          }
        } else if (value_type_m1 == 1) {
          const uint16_t* filter_vec_alias = R_CAST(const uint16_t*, parse_iter);
          for (uint32_t filter_idx = 0; filter_idx != value_ct; ++filter_idx) {
            const uint32_t cur_val = filter_vec_alias[filter_idx];
            if (cur_val == 0x8000) {
              // missing
              ++cur_filter_slen;
            } else {
              // loadbuf_size >= 2^18, so fif_slens[65535] won't segfault
              const uint32_t key_slen = fif_slens[cur_val];
              if (unlikely((cur_val >= fif_string_idx_end) || (cur_val > 0x8000) || (!key_slen))) {
                goto BcfToPgen_ret_VREC_GENERIC;
              }
              cur_filter_slen += key_slen;
            }
          }
        } else {
          const uint32_t* filter_vec_alias = R_CAST(const uint32_t*, parse_iter);
          for (uint32_t filter_idx = 0; filter_idx != value_ct; ++filter_idx) {
            const uint32_t cur_val = filter_vec_alias[filter_idx];
            if (cur_val == 0x80000000U) {
              // missing
              ++cur_filter_slen;
            } else {
              if (unlikely(cur_val >= fif_string_idx_end)) {
                goto BcfToPgen_ret_VREC_GENERIC;
              }
              const uint32_t key_slen = fif_slens[cur_val];
              if (unlikely(!key_slen)) {
                goto BcfToPgen_ret_VREC_GENERIC;
              }
              cur_filter_slen += key_slen;
            }
          }
        }
        if (cur_filter_slen > filter_info_slen_ubound) {
          filter_info_slen_ubound = cur_filter_slen;
        }
        parse_iter = &(parse_iter[vec_byte_ct]);
      }
      // INFO
      uintptr_t info_pr_here = 0;
      if (n_info) {
        // bugfix (1 Jul 2021): avoid underflow when n_info == 0
        uint64_t cur_info_slen_ubound = n_info - 1;
        for (uint32_t uii = 0; uii != n_info; ++uii) {
          uint32_t sidx;
          if (unlikely(ScanBcfTypedInt(&parse_iter, &sidx) || (parse_iter > shared_end) || (sidx >= fif_string_idx_end))) {
            goto BcfToPgen_ret_VREC_GENERIC;
          }
          const uint32_t key_slen = fif_slens[sidx];
          if (unlikely((!key_slen) || ScanBcfType(&parse_iter, &value_type, &value_ct) || (parse_iter > shared_end))) {
            goto BcfToPgen_ret_VREC_GENERIC;
          }
          cur_info_slen_ubound += key_slen;
          if (value_ct) {
            const uint32_t bytes_per_elem = kBcfBytesPerElem[value_type];
            const uint64_t vec_byte_ct = bytes_per_elem * S_CAST(uint64_t, value_ct);
            if (unlikely((!bytes_per_elem) || (S_CAST(uint64_t, shared_end - parse_iter) < vec_byte_ct))) {
              goto BcfToPgen_ret_VREC_GENERIC;
            }
            if (value_type == 1) {
              // int8, max len 4 ("-127"), add 1 for comma
              cur_info_slen_ubound += 5 * S_CAST(uint64_t, value_ct);
            } else if (value_type == 2) {
              // int16, max len 6
              cur_info_slen_ubound += 7 * S_CAST(uint64_t, value_ct);
            } else if (value_type == 3) {
              // int32, max len 11
              cur_info_slen_ubound += 12 * S_CAST(uint64_t, value_ct);
            } else if (value_type == 5) {
              cur_info_slen_ubound += (kMaxFloatGSlen + 1) * S_CAST(uint64_t, value_ct);
            } else {
              // string
              cur_info_slen_ubound += value_ct + 1;
            }
            parse_iter = &(parse_iter[vec_byte_ct]);
          } else {
            // bugfix (8 Jan 2023): value_type == value_ct == 0
            // "untyped-missing" special case is always permitted, it's not
            // restricted to flags.
            //
            // value_ct == 0, value_type == 7 is also valid (string, missing).
            //
            // Otherwise, value_ct == 0 should not happen.
            if (value_type) {
              if (unlikely(value_type != 7)) {
                goto BcfToPgen_ret_VREC_GENERIC;
              }
              // Emit '=' with empty-string value.
              ++cur_info_slen_ubound;
            }
          }
          if (sidx == pr_sidx) {
            if (unlikely(info_pr_here)) {
              putc_unlocked('\n', stdout);
              snprintf(g_logbuf, kLogbufSize, "Error: Variant record #%" PRIuPTR " in --bcf file has multiple INFO/PR entries.\n", vrec_idx);
              goto BcfToPgen_ret_MALFORMED_INPUT_WW;
            }
            info_pr_here = 1;
          }
        }
        if (cur_info_slen_ubound > filter_info_slen_ubound) {
          filter_info_slen_ubound = cur_info_slen_ubound;
        }
      }
      if (unlikely(parse_iter != shared_end)) {
        goto BcfToPgen_ret_VREC_GENERIC;
      }

      if (max_allele_ct < n_allele) {
        if (n_allele > kPglMaxAlleleCt) {
          putc_unlocked('\n', stdout);
          logerrprintfww("Error: Variant record #%" PRIuPTR " in --bcf file has %u alleles; this build of " PROG_NAME_STR " is limited to " PGL_MAX_ALLELE_CT_STR ". (You can use \"--import-max-alleles " PGL_MAX_ALLELE_CT_STR "\" to filter out such variants.)\n", vrec_idx, n_allele);
          reterr = kPglRetNotYetSupported;
          goto BcfToPgen_ret_1;
        }
        max_allele_ct = n_allele;
      }
      if (max_observed_rec_vecs < cur_observed_rec_vecs) {
        max_observed_rec_vecs = cur_observed_rec_vecs;
      }
      allele_idx_offsets[variant_ct] = allele_idx_end;
      allele_idx_end += n_allele;
      const uint32_t variant_idx_lowbits = variant_ct % kBitsPerWord;
      if (pr_sidx != UINT32_MAX) {
        nonref_word |= info_pr_here << variant_idx_lowbits;
        if (variant_idx_lowbits == (kBitsPerWord - 1)) {
          *nonref_flags_iter++ = nonref_word;
          nonref_word = 0;
        }
      }
      if (sample_ct && (!phase_or_dosage_found)) {
        // Check if there's at least one phased het call, and/or at least one
        // relevant dosage.
        // Don't bother multithreading this since it's trivial.
        if (dosage_start || hds_start) {
          if (n_allele == 2) {
            bcf_parse_err = BcfScanBiallelicHds(&(ctx.bic), gt_start, qual_starts, dosage_start, hds_start, &phase_or_dosage_found, sample_nypbuf);
          } else {
            putc_unlocked('\n', stdout);
            logerrputs("Error: --bcf multiallelic dosage import is under development.\n");
            reterr = kPglRetNotYetSupported;
            goto BcfToPgen_ret_1;
          }
        } else if (gt_start) {
          // scan for phase
          bcf_parse_err = BcfScanGt(&(ctx.bic), gt_start, qual_starts, &phase_or_dosage_found, sample_nypbuf);
        }
        if (unlikely(bcf_parse_err)) {
          goto BcfToPgen_ret_PARSE;
        }
      }
      if (unlikely(variant_ct++ == max_variant_ct)) {
#ifdef __LP64__
        if (variant_ct == kPglMaxVariantCt) {
          putc_unlocked('\n', stdout);
          logerrputs("Error: " PROG_NAME_STR " does not support more than 2^31 - 3 variants.  We recommend using\nother software for very deep studies of small numbers of genomes.\n");
          goto BcfToPgen_ret_MALFORMED_INPUT;
        }
#endif
        goto BcfToPgen_ret_NOMEM;
      }
      if (!(variant_ct % 10000)) {
        printf("\r--bcf: %uk variants scanned.", variant_ct / 1000);
        fflush(stdout);
      }
    }
    const uintptr_t variant_skip_ct = vrec_idx - 1 - variant_ct;
    if (variant_ct % kBitsPerWord) {
      if (nonref_flags_iter) {
        *nonref_flags_iter = nonref_word;
      }
    } else if (unlikely(!variant_ct)) {
      if (!variant_skip_ct) {
        logerrputs("Error: No variants in --bcf file.\n");
        goto BcfToPgen_ret_DEGENERATE_DATA;
      }
      char* write_iter = strcpya_k(g_logbuf, "Error: All ");
      write_iter = wtoa(variant_skip_ct, write_iter);
      write_iter = strcpya_k(write_iter, " variant");
      if (variant_skip_ct != 1) {
        *write_iter++ = 's';
      }
      write_iter = strcpya_k(write_iter, " in --bcf file excluded by ");
      AppendLoadFilterFlagnames(load_filter_log_import_flags, &write_iter);
      strcpy_k(write_iter, ".\n");
      goto BcfToPgen_ret_INCONSISTENT_INPUT_WW;
    }

    putc_unlocked('\r', stdout);
    {
      char* write_iter = strcpya_k(g_logbuf, "--bcf: ");
      write_iter = wtoa(variant_ct + variant_skip_ct, write_iter);
      write_iter = strcpya_k(write_iter, " variant");
      if (variant_ct + variant_skip_ct != 1) {
        *write_iter++ = 's';
      }
      write_iter = strcpya_k(write_iter, " scanned");
      if (variant_skip_ct) {
        write_iter = strcpya_k(write_iter, "; ");
        write_iter = wtoa(variant_skip_ct, write_iter);
        write_iter = strcpya_k(write_iter, " excluded by ");
        AppendLoadFilterFlagnames(load_filter_log_import_flags, &write_iter);
        write_iter = strcpya_k(write_iter, ", ");
        write_iter = u32toa(variant_ct, write_iter);
        write_iter = strcpya_k(write_iter, " remaining");
      } else if (load_filter_log_import_flags) {
        write_iter = strcpya_k(write_iter, " (");
        AppendLoadFilterFlagnames(load_filter_log_import_flags, &write_iter);
        write_iter = strcpya_k(write_iter, " had no effect)");
      }
      strcpy_k(write_iter, ".\n");
      WordWrapB(0);
      logputsb();
    }

    if (allele_idx_end > 2 * variant_ct) {
      allele_idx_offsets[variant_ct] = allele_idx_end;
      BigstackFinalizeW(allele_idx_offsets, variant_ct + 1);
    } else {
      allele_idx_offsets = nullptr;
    }

    BigstackEndReset(bigstack_end_mark2);
    loadbuf_size = RoundUpPow2(loadbuf_size_needed, kEndAllocAlign);
    loadbuf = S_CAST(unsigned char*, bigstack_end_alloc_raw(loadbuf_size));

    reterr = BgzfRawMtStreamRewind(&bgzf, &bgzf_errmsg);
    if (unlikely(reterr)) {
      if (reterr == kPglRetDecompressFail) {
        DPrintf("DecompressFail during initial rewind");
        goto BcfToPgen_ret_REWIND_FAIL;
      }
      goto BcfToPgen_ret_BGZF_FAIL;
    }
    {
      unsigned char* bcf_header_buf = R_CAST(unsigned char*, &(vcf_header[-9]));
      unsigned char* bcf_iter = bcf_header_buf;
      unsigned char* bcf_header_end = &(bcf_header_buf[header_size + 9 * k1LU]);
      reterr = BgzfRawMtStreamRead(bcf_header_end, &bgzf, &bcf_iter, &bgzf_errmsg);
      if (unlikely(reterr)) {
        goto BcfToPgen_ret_BGZF_FAIL;
      }
      if (unlikely(bcf_iter != bcf_header_end)) {
        errno = 0;
        reterr = kPglRetReadFail;
        goto BcfToPgen_ret_BGZF_FAIL;
      }
    }

    uint32_t calc_thread_ct;
    // todo: tune this for BCF, these values were derived from VCF testing
    if (phase_or_dosage_found && (dosage_sidx || hds_sidx)) {
      calc_thread_ct = 1 + (sample_ct > 5) + (sample_ct > 12) + (sample_ct > 32) + (sample_ct > 512);
    } else {
      calc_thread_ct = 1 + (sample_ct > 40) + (sample_ct > 320);
    }
    if (calc_thread_ct + decompress_thread_ct > max_thread_ct) {
      calc_thread_ct = MAXV(1, max_thread_ct - decompress_thread_ct);
    }

    const uint32_t variant_ctl = BitCtToWordCt(variant_ct);
    PgenGlobalFlags phase_dosage_gflags = kfPgenGlobal0;
    GparseFlags gparse_flags = kfGparse0;  // yeah, this is a bit redundant
    if (phase_or_dosage_found || hds_sidx || dosage_sidx) {
      if ((!hds_sidx) && (!dosage_sidx)) {
        phase_dosage_gflags = kfPgenGlobalHardcallPhasePresent;
        gparse_flags = kfGparseHphase;
      } else {
        // thanks to automatic --hard-call-threshold, we may need to save
        // dosage-phase even when there's no HDS field
        phase_dosage_gflags = kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent;
        gparse_flags = kfGparseHphase | kfGparseDosage | kfGparseDphase;
      }
    }
    uint32_t nonref_flags_storage = 1;
    if (nonref_flags) {
      const uint32_t variant_ctl_m1 = variant_ctl - 1;
      const uintptr_t last_nonref_flags_word = nonref_flags[variant_ctl_m1];
      if (!last_nonref_flags_word) {
        for (uint32_t widx = 0; widx != variant_ctl_m1; ++widx) {
          if (nonref_flags[widx]) {
            nonref_flags_storage = 3;
            break;
          }
        }
      } else if (!((~last_nonref_flags_word) << ((-variant_ct) & (kBitsPerWord - 1)))) {
        nonref_flags_storage = 2;
        for (uint32_t widx = 0; widx != variant_ctl_m1; ++widx) {
          if (~nonref_flags[widx]) {
            nonref_flags_storage = 3;
            break;
          }
        }
      } else {
        nonref_flags_storage = 3;
      }
      if (nonref_flags_storage != 3) {
        // yeah, we may now have a temporary "memory leak" here (if
        // multiallelic variants are present, and thus allele_idx_offsets[]
        // must be kept around), but this array is typically only 1/64 the size
        // of allele_idx_offsets[].
        if (!allele_idx_offsets) {
          BigstackReset(nonref_flags);
        }
        nonref_flags = nullptr;
      }
    }

    const uint64_t write_chunk_ubound = 32 + MAXV(other_slen_ubound, filter_info_slen_ubound);
#ifndef __LP64__
    if (write_chunk_ubound > 0x7ff00000) {
      goto BcfToPgen_ret_NOMEM;
    }
#endif
    const uint32_t output_zst = ((import_flags & (kfImportKeepAutoconv | kfImportKeepAutoconvVzs)) != kfImportKeepAutoconv);
    snprintf(outname_end, kMaxOutfnameExtBlen, ".pvar");
    if (output_zst) {
      snprintf(&(outname_end[5]), kMaxOutfnameExtBlen - 5, ".zst");
    }
    const uintptr_t overflow_buf_size = kCompressStreamBlock + MAXV(write_chunk_ubound, kCompressStreamBlock);
    reterr = InitCstreamAlloc(outname, 0, output_zst, sample_ct? 1 : MAXV(1, max_thread_ct - decompress_thread_ct), overflow_buf_size, &pvar_css, &pvar_cswritep);
    if (unlikely(reterr)) {
      goto BcfToPgen_ret_1;
    }
    line_iter = vcf_header;
    uint32_t contig_idx = UINT32_MAX;  // deliberate overflow
    uint32_t idxeq_clipped = 0;
    for (header_line_idx = 1; header_line_idx != header_line_ct; ++header_line_idx) {
      char* line_start = line_iter;
      char* line_main = &(line_iter[2]);
      char* line_end = AdvPastDelim(line_main, '\n');
      line_iter = line_end;
      // chrSet skipped here since we call AppendChrsetLine after this loop
      if (StrStartsWithUnsafe(line_main, "fileformat=") || StrStartsWithUnsafe(line_main, "fileDate=") || StrStartsWithUnsafe(line_main, "source=") || StrStartsWithUnsafe(line_main, "FORMAT=") || StrStartsWithUnsafe(line_main, "chrSet=")) {
        continue;
      }
      // OS-agnostic newline clip.
      char* line_write_end = &(line_end[-1]);
      if (line_write_end[-1] == '\r') {
        --line_write_end;
      }
      const uint32_t is_contig_line = StrStartsWithUnsafe(line_main, "contig=<ID");
      if (explicit_idx_keys) {
        const uint32_t is_filter_line = StrStartsWithUnsafe(line_main, "FILTER=<ID");
        const uint32_t is_info_line = StrStartsWithUnsafe(line_main, "INFO=<ID");
        const uint32_t is_format_line = StrStartsWithUnsafe(line_main, "FORMAT=<ID");
        idxeq_clipped = is_contig_line || is_filter_line || is_info_line || is_format_line;
        if (idxeq_clipped) {
          // Clip ",IDX=123>".
          --line_write_end;
          do {
            --line_write_end;
          } while (IsDigit(*line_write_end));
          line_write_end = &(line_write_end[-4]);
          // line_write_end now points to the '=' in 'IDX='.
        }
      }
      if (is_contig_line) {
        if (explicit_idx_keys) {
          ScanUintDefcap(&(line_write_end[5]), &contig_idx);
        } else {
          ++contig_idx;
        }
        if (!IsSet(bcf_contig_keep, contig_idx)) {
          continue;
        }
      }
      if (unlikely(CsputsStd(line_start, line_write_end - line_start, &pvar_css, &pvar_cswritep))) {
        goto BcfToPgen_ret_WRITE_FAIL;
      }
      if (idxeq_clipped) {
        *pvar_cswritep++ = '>';
      }
      // NOT safe to use AppendBinaryEoln here.
#ifdef _WIN32
      pvar_cswritep = strcpya_k(pvar_cswritep, "\r\n");
#else
      *pvar_cswritep++ = '\n';
#endif
    }
    if (cip->chrset_source) {
      AppendChrsetLine(cip, &pvar_cswritep);
    }
    pvar_cswritep = strcpya_k(pvar_cswritep, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER");
    if (info_nonpr_exists) {
      pvar_cswritep = strcpya_k(pvar_cswritep, "\tINFO");
    }
    AppendBinaryEoln(&pvar_cswritep);
    if (unlikely(Cswrite(&pvar_css, &pvar_cswritep))) {
      goto BcfToPgen_ret_WRITE_FAIL;
    }

    // may as well have a functional progress meter in no-samples case
    uint32_t main_block_size = MINV(65536, variant_ct);
    uint32_t per_thread_block_limit = main_block_size;
    uint32_t cur_thread_block_vidx_limit = 1;
    uintptr_t per_thread_byte_limit = 0;
    unsigned char* geno_bufs[2];
    // defensive
    geno_bufs[0] = nullptr;
    geno_bufs[1] = nullptr;
    if (hard_call_thresh == UINT32_MAX) {
      hard_call_thresh = kDosageMid / 10;
    }
    uintptr_t* bcf_haploid_mask = nullptr;
    const uint32_t hard_call_halfdist = kDosage4th - hard_call_thresh;

    // See the analogous code in VcfToPgen().
    const uint32_t sex_info_avail = preexisting_psamname || is_update_or_impute_sex;
    uint32_t print_splitpar_warning = 0;
    if (par_warn_bcf_chrom != UINT32_MAX) {
      if ((import_flags & kfImportLaxChrX) || (!sample_ct)) {
        par_warn_bcf_chrom = UINT32_MAX;
      } else if (unlikely(!sex_info_avail)) {
        logerrputs("Error: chrX is present in the input file, but no sex information was provided;\nrerun this import with --psam, --update-sex, or --impute-sex.  --split-par may\nalso be appropriate.\n");
        goto BcfToPgen_ret_INCONSISTENT_INPUT;
      } else if ((!IsHumanChrset(cip)) || is_splitpar) {
        par_warn_bcf_chrom = UINT32_MAX;
      }
    }
    if (sample_ct) {
      if (dosage_import_field && strequal_k(dosage_import_field, "DS", dosage_import_field_slen) && (!require_gt)) {
        // Only need to initialize this when enforcing DS/haploid rule.
        const uint32_t word_ct = BitCtToWordCt(contig_string_idx_end);
        if (unlikely(bigstack_alloc_w(word_ct, &bcf_haploid_mask))) {
          goto BcfToPgen_ret_NOMEM;
        }
        if (cip->haploid_mask[0] & 1) {
          // all observed chromosomes are haploid
          memcpy(bcf_haploid_mask, bcf_contig_keep, word_ct * sizeof(intptr_t));
        } else {
          // chrX, chrY, and chrM are the affected chromosomes
          // bugfix (17 Jun 2022): we should be treating e.g. 23 and chrX
          // identically here, rather than only recognizing the latter
          ZeroWArr(word_ct, bcf_haploid_mask);
          const uint32_t y_code = cip->xymt_codes[kChrOffsetY];
          const uint32_t mt_code = cip->xymt_codes[kChrOffsetMT];
          uint32_t haploid_found = 0;
          for (uint32_t widx = 0; widx != word_ct; ++widx) {
            uintptr_t cur_word = bcf_contig_keep[widx];
            if (cur_word) {
              const uint32_t chrom_base = widx * kBitsPerWord;
              do {
                const uint32_t chrom = chrom_base + ctzw(cur_word);
                const uint32_t chr_code_raw = GetChrCodeRaw(contig_names[chrom]);
                if (((chr_code_raw < kChrRawX) && ((chr_code_raw == x_code) || (chr_code_raw == y_code) || (chr_code_raw == mt_code))) ||
                    (chr_code_raw == kChrRawX) || (chr_code_raw == kChrRawY) || (chr_code_raw == kChrRawMT)) {
                  haploid_found = 1;
                  SetBit(chrom, bcf_haploid_mask);
                }
                cur_word &= cur_word - 1;
              } while (cur_word);
            }
          }
          if (!haploid_found) {
            BigstackReset(bcf_haploid_mask);
            bcf_haploid_mask = nullptr;
          }
        }
      }
      snprintf(outname_end, kMaxOutfnameExtBlen, ".pgen");
      uintptr_t spgw_alloc_cacheline_ct;
      uint32_t max_vrec_len;
      reterr = SpgwInitPhase1(outname, allele_idx_offsets, nonref_flags, variant_ct, sample_ct, max_allele_ct, kPgenWriteBackwardSeek, phase_dosage_gflags, nonref_flags_storage, &spgw, &spgw_alloc_cacheline_ct, &max_vrec_len);
      if (unlikely(reterr)) {
        if (reterr == kPglRetOpenFail) {
          logerrprintfww(kErrprintfFopen, outname, strerror(errno));
        }
        goto BcfToPgen_ret_1;
      }
      unsigned char* spgw_alloc;
      if (unlikely(bigstack_alloc_uc(spgw_alloc_cacheline_ct * kCacheline, &spgw_alloc))) {
        goto BcfToPgen_ret_NOMEM;
      }
      SpgwInitPhase2(max_vrec_len, &spgw, spgw_alloc);

      if (unlikely(bigstack_alloc_ucp(calc_thread_ct, &ctx.thread_wkspaces) ||
                   bigstack_alloc_u32(calc_thread_ct + 1, &(ctx.thread_bidxs[0])) ||
                   bigstack_alloc_u32(calc_thread_ct + 1, &(ctx.thread_bidxs[1])) ||
                   bigstack_calloc_w(calc_thread_ct, &ctx.err_vrec_idxs))) {
        goto BcfToPgen_ret_NOMEM;
      }
      ctx.bcf_parse_errs = S_CAST(BcfParseErr*, bigstack_alloc_raw_rd(calc_thread_ct * sizeof(BcfParseErr)));
      if (unlikely(!ctx.bcf_parse_errs)) {
        goto BcfToPgen_ret_NOMEM;
      }
      ctx.hard_call_halfdist = hard_call_halfdist;
      ctx.parse_failed = 0;
      // defensive
      ctx.gparse[0] = nullptr;
      ctx.gparse[1] = nullptr;
      ctx.block_allele_idx_offsets[0] = nullptr;
      ctx.block_allele_idx_offsets[1] = nullptr;
      // Finished with all other memory allocations, so all remaining workspace
      // can be spent on multithreaded parsing.  Spend up to 1/6 on
      // g_thread_wkspaces (tune this fraction later).
      // Probable todo: factor out common parts with bgen-1.3 initialization
      // into separate function(s).
      uint64_t max_write_byte_ct = GparseWriteByteCt(sample_ct, max_allele_ct, gparse_flags);
      // always allocate tmp_dphase_delta for now
      uint64_t thread_wkspace_cl_ct = DivUp(max_write_byte_ct + sample_ct * sizeof(SDosage), kCacheline);
      uintptr_t cachelines_avail = bigstack_left() / (6 * kCacheline);
      if (calc_thread_ct * thread_wkspace_cl_ct > cachelines_avail) {
        if (unlikely(thread_wkspace_cl_ct > cachelines_avail)) {
          goto BcfToPgen_ret_NOMEM;
        }
        calc_thread_ct = cachelines_avail / thread_wkspace_cl_ct;
      }
      for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
        ctx.thread_wkspaces[tidx] = S_CAST(unsigned char*, bigstack_alloc_raw(thread_wkspace_cl_ct * kCacheline));
        ctx.bcf_parse_errs[tidx] = kBcfParseOk;
      }

      // be pessimistic re: rounding
      cachelines_avail = (bigstack_left() / kCacheline) - 4;
      const uint64_t max_bytes_req_per_variant = sizeof(GparseRecord) + MAXV(max_observed_rec_vecs * S_CAST(uintptr_t, kBytesPerVec), max_write_byte_ct) + calc_thread_ct;
      if (unlikely(cachelines_avail * kCacheline < 2 * max_bytes_req_per_variant)) {
        goto BcfToPgen_ret_NOMEM;
      }
      // use worst-case gparse_flags since lines will usually be similar
      uintptr_t min_bytes_req_per_variant = sizeof(GparseRecord) + GparseWriteByteCt(sample_ct, 2, gparse_flags);
      main_block_size = (cachelines_avail * kCacheline) / (min_bytes_req_per_variant * 2);
      // this is arbitrary, there's no connection to kPglVblockSize
      if (main_block_size > 65536) {
        main_block_size = 65536;
      }
      // divide by 2 for better parallelism in small-variant-count case
      // round up per_thread_block_limit so we only have two blocks
      if (main_block_size > DivUp(variant_ct, 2)) {
        main_block_size = DivUp(variant_ct, 2) + calc_thread_ct - 1;
      }
      // may as well guarantee divisibility
      per_thread_block_limit = main_block_size / calc_thread_ct;
      main_block_size = per_thread_block_limit * calc_thread_ct;
      if (unlikely(SetThreadCt(calc_thread_ct, &tg))) {
        goto BcfToPgen_ret_NOMEM;
      }
      ctx.gparse[0] = S_CAST(GparseRecord*, bigstack_alloc_raw_rd(main_block_size * sizeof(GparseRecord)));
      ctx.gparse[1] = S_CAST(GparseRecord*, bigstack_alloc_raw_rd(main_block_size * sizeof(GparseRecord)));
      SetThreadFuncAndData(BcfGenoToPgenThread, &ctx, &tg);
      cachelines_avail = bigstack_left() / (kCacheline * 2);
      geno_bufs[0] = S_CAST(unsigned char*, bigstack_alloc_raw(cachelines_avail * kCacheline));
      geno_bufs[1] = S_CAST(unsigned char*, bigstack_alloc_raw(cachelines_avail * kCacheline));
      // This is only used for comparison purposes, so it is unnecessary to
      // round it down to a multiple of kBytesPerVec even though every actual
      // record will be vector-aligned.
      per_thread_byte_limit = (cachelines_avail * kCacheline) / calc_thread_ct;
    }

    vrec_idx = 0;
    // Main workflow:
    // 1. Set n=0, load genotype data for first main_block_size variants
    //    while writing .pvar
    //
    // 2. Spawn threads processing batch n genotype data
    // 3. If n>0, write results for block (n-1)
    // 4. Increment n by 1
    // 5. Load/write-.pvar for batch (n+1) unless eof
    // 6. Join threads
    // 7. Goto step 2 unless eof
    //
    // 8. Write results for last block
    uint32_t prev_block_write_ct = 0;
    uintptr_t record_byte_ct = 0;
    uint32_t* thread_bidxs = nullptr;
    GparseRecord* cur_gparse = nullptr;
    unsigned char* geno_buf_iter = nullptr;
    unsigned char* cur_thread_byte_stop = nullptr;

    // Placed here since these values need to persist to the next block
    // iteration when we run out of geno_bufs[parity] space.
    uint32_t chrom = 0;
    uint32_t pos = 0;
    uint32_t qual_bits = 0;
    uint32_t n_allele = 0;
    uint32_t n_info = 0;
    uint32_t n_fmt = 0;
    const unsigned char* shared_end = nullptr;

    const unsigned char* gt_start = nullptr;
    const unsigned char* qual_starts[2];
    qual_starts[0] = nullptr;
    qual_starts[1] = nullptr;
    const unsigned char* dosage_start = nullptr;
    const unsigned char* hds_start = nullptr;
    uint32_t gt_type_blen = 0;
    uint32_t gt_main_blen = 0;
    uint32_t qual_type_blens[2];
    uint32_t qual_main_blens[2];
    qual_type_blens[0] = 0;
    qual_type_blens[1] = 0;
    qual_main_blens[0] = 0;
    qual_main_blens[1] = 0;
    uint32_t dosage_type_blen = 0;
    uint32_t dosage_main_blen = 0;
    uint32_t hds_type_blen = 0;
    uint32_t hds_main_blen = 0;
    uint32_t record_input_vec_ct = 0;
    uint32_t unflushed_line = 0;

    uint32_t parity = 0;
    for (uint32_t vidx_start = 0; ; ) {
      uint32_t cur_block_write_ct = 0;
      if (!IsLastBlock(&tg)) {
        const uint32_t block_vidx_limit = variant_ct - vidx_start;
        cur_thread_block_vidx_limit = MINV(block_vidx_limit, per_thread_block_limit);
        uint32_t cur_thread_fill_idx = 0;
        if (sample_ct) {
          thread_bidxs = ctx.thread_bidxs[parity];
          cur_gparse = ctx.gparse[parity];
          if (allele_idx_offsets) {
            ctx.block_allele_idx_offsets[parity] = &(allele_idx_offsets[vidx_start]);
          }
          geno_buf_iter = geno_bufs[parity];
          cur_thread_byte_stop = &(geno_buf_iter[per_thread_byte_limit]);
          thread_bidxs[0] = 0;
        }
        uint32_t block_vidx = 0;
        uint32_t sidx = 0;
        GparseRecord* grp;
        if (unflushed_line) {
          // We haven't actually written the current .pvar line or copied
          // genotype/quality data over.
          unflushed_line = 0;
          goto BcfToPgen_load_keep;
        }
        while (1) {
          {
            ++vrec_idx;
            unsigned char* loadbuf_read_iter = loadbuf;
            reterr = BgzfRawMtStreamRead(&(loadbuf[32]), &bgzf, &loadbuf_read_iter, &bgzf_errmsg);
            if (unlikely(reterr)) {
              goto BcfToPgen_ret_BGZF_FAIL_N;
            }
            if (&(loadbuf[32]) != loadbuf_read_iter) {
              DPrintf("read only %" PRIdPTR " out of 32 initial bytes, vrec_idx=%" PRIuPTR "\n", loadbuf_read_iter - loadbuf, vrec_idx);
              goto BcfToPgen_ret_REWIND_FAIL_N;
            }
            // [0]: l_shared
            // [1]: l_indiv
            // [2]: chrom
            // [3]: pos
            // [4]: rlen
            // [5]: qual, NOT n_allele_info
            // [6]: low 16 bits n_info, high 16 bits n_allele
            // [7]: low 24 bits n_sample, high 8 bits n_fmt
            const uint32_t* vrec_header = R_CAST(uint32_t*, loadbuf);
            const uint32_t l_shared = vrec_header[0];
            const uint32_t l_indiv = vrec_header[1];
            chrom = vrec_header[2];
            pos = vrec_header[3];
            qual_bits = vrec_header[5];
            n_allele = vrec_header[6] >> 16;
            n_info = vrec_header[6] & 0xffff;
            n_fmt = vrec_header[7] >> 24;
            // skip validation performed in first pass
            const uint64_t second_load_size = l_shared + S_CAST(uint64_t, l_indiv) - 24;
            loadbuf_read_iter = loadbuf;
            unsigned char* indiv_end = &(loadbuf[second_load_size]);
            reterr = BgzfRawMtStreamRead(indiv_end, &bgzf, &loadbuf_read_iter, &bgzf_errmsg);
            if (unlikely(reterr)) {
              goto BcfToPgen_ret_BGZF_FAIL_N;
            }
            if (unlikely(indiv_end != loadbuf_read_iter)) {
              DPrintf("read only %" PRIdPTR " out of %" PRIu64 " later bytes, vrec_idx=%" PRIuPTR "\n", loadbuf_read_iter - &(loadbuf[32]), second_load_size - 32, vrec_idx);
              goto BcfToPgen_ret_REWIND_FAIL_N;
            }

            // 1. check if we skip this variant.  chromosome filter,
            //    import_max_allele_ct, and require_gt can cause this.
            if ((n_allele > import_max_allele_ct) || (!IsSet(bcf_contig_keep, chrom))) {
              continue;
            }
            const uint32_t fail_on_ds_only = bcf_haploid_mask && IsSet(bcf_haploid_mask, chrom);

            // obvious todo: move duplicated code between first and second pass
            // into separate functions
            shared_end = indiv_end - l_indiv;
            const unsigned char* parse_iter = shared_end;

            gt_start = nullptr;
            qual_starts[0] = nullptr;
            qual_starts[1] = nullptr;
            dosage_start = nullptr;
            hds_start = nullptr;
            record_input_vec_ct = 0;
            uint32_t gt_exists = 0;
            for (uint32_t fmt_idx = 0; fmt_idx != n_fmt; ++fmt_idx) {
              // 1. typed int indicating which FORMAT field
              // 2. shared type descriptor for each entry (usually a single
              //    byte)
              // 3. sample_ct entries
              // previously validated
              ScanBcfTypedInt(&parse_iter, &sidx);
              const unsigned char* type_start = parse_iter;
              uint32_t value_type;
              uint32_t value_ct;
              ScanBcfType(&parse_iter, &value_type, &value_ct);
              const uint32_t bytes_per_elem = kBcfBytesPerElem[value_type];
              const uint32_t vec_byte_ct = bytes_per_elem * value_ct * sample_ct;
              const uint32_t type_blen = parse_iter - type_start;
              parse_iter = &(parse_iter[vec_byte_ct]);
              if ((!sidx) || (!value_ct)) {
                if (sidx == gt_sidx) {
                  gt_exists = 1;
                }
                continue;
              }
              if (sidx == gt_sidx) {
                // bugfix (25 Oct 2025)
                gt_exists = 1;
                gt_start = type_start;
                gt_type_blen = type_blen;
                gt_main_blen = vec_byte_ct;
              } else if (sidx == gq_sidx) {
                qual_starts[0] = type_start;
                qual_type_blens[0] = type_blen;
                qual_main_blens[0] = vec_byte_ct;
              } else if (sidx == dp_sidx) {
                qual_starts[1] = type_start;
                qual_type_blens[1] = type_blen;
                qual_main_blens[1] = vec_byte_ct;
              } else if (sidx == dosage_sidx) {
                dosage_start = type_start;
                dosage_type_blen = type_blen;
                dosage_main_blen = vec_byte_ct;
              } else if (sidx == hds_sidx) {
                hds_start = type_start;
                hds_type_blen = type_blen;
                hds_main_blen = vec_byte_ct;
              } else {
                continue;
              }
#ifdef __LP64__
              ++record_input_vec_ct;
#else
              record_input_vec_ct += DivUp(type_blen, kBytesPerVec);
#endif
              record_input_vec_ct += DivUp(vec_byte_ct, kBytesPerVec);
            }
            if (require_gt && (!gt_exists)) {
              continue;
            }
            // We already have enough information to determine
            // write_byte_ct_limit.
            if ((!gt_start) && (!dosage_sidx) && (!hds_sidx)) {
              gparse_flags = kfGparseNull;
            } else {
              if ((!phase_or_dosage_found) && (!dosage_sidx) && (!hds_sidx)) {
                gparse_flags = kfGparse0;
              } else {
                if (unlikely(fail_on_ds_only && (!gt_start) && dosage_start && (!hds_start))) {
                  snprintf(g_logbuf, kLogbufSize, "Error: Variant record #%" PRIuPTR " of --bcf file is for a chrX, chrM, or fully-haploid variant, and has a DS field without a companion GT field to clarify whether each DS value is on a 0..1 or 0..2 scale. This cannot be imported by " PROG_NAME_STR "; please e.g. regenerate the file with GT present.\n", vrec_idx);
                  goto BcfToPgen_ret_MALFORMED_INPUT_WWN;
                }
                gparse_flags = (dosage_start || hds_start)? (kfGparseHphase | kfGparseDosage | kfGparseDphase) : kfGparseHphase;
              }
            }
            const uintptr_t write_byte_ct_limit = GparseWriteByteCt(sample_ct, n_allele, gparse_flags);
            record_byte_ct = MAXV(record_input_vec_ct * S_CAST(uintptr_t, kBytesPerVec), write_byte_ct_limit);
            if ((block_vidx == cur_thread_block_vidx_limit) || (S_CAST(uintptr_t, cur_thread_byte_stop - geno_buf_iter) < record_byte_ct)) {
              thread_bidxs[++cur_thread_fill_idx] = block_vidx;
              if (cur_thread_fill_idx == calc_thread_ct) {
                unflushed_line = 1;
                break;
              }
              cur_thread_byte_stop = &(cur_thread_byte_stop[per_thread_byte_limit]);
              cur_thread_block_vidx_limit = MINV(cur_thread_block_vidx_limit + per_thread_block_limit, block_vidx_limit);
            }
          }
        BcfToPgen_load_keep:
          // CHROM, POS
          pvar_cswritep = memcpyax(pvar_cswritep, contig_out_names[chrom], contig_out_slens[chrom], '\t');
          const uint32_t pos1 = pos + 1;
          pvar_cswritep = u32toa_x(pos1, '\t', pvar_cswritep);
          if (unlikely(Cswrite(&pvar_css, &pvar_cswritep))) {
            goto BcfToPgen_ret_WRITE_FAIL;
          }
          if (par_warn_bcf_chrom == chrom) {
            if ((pos1 <= kPAR1IntersectionLast) || (pos1 >= kPAR2IntersectionFirst)) {
              if (unlikely(!is_sortvars)) {
                putc_unlocked('\n', stdout);
                logerrputs("Error: Human chrX pseudoautosomal variant(s) appear to be present in the input\nBCF, but --split-par was not specified.\n");
                goto BcfToPgen_ret_INCONSISTENT_INPUT;
              }
              print_splitpar_warning = 1;
              par_warn_bcf_chrom = UINT32_MAX;
            }
          }

          // ID
          const unsigned char* parse_iter = loadbuf;
          uint32_t slen;
          {
            const char* id_start;
            ScanBcfTypedString(shared_end, &parse_iter, &id_start, &slen);
            if (slen) {
              pvar_cswritep = memcpya(pvar_cswritep, id_start, slen);
            } else {
              *pvar_cswritep++ = '.';
            }
            *pvar_cswritep++ = '\t';
            if (unlikely(Cswrite(&pvar_css, &pvar_cswritep))) {
              goto BcfToPgen_ret_WRITE_FAIL;
            }

            // REF
            const char* ref_start;
            ScanBcfTypedString(shared_end, &parse_iter, &ref_start, &slen);
            pvar_cswritep = memcpyax(pvar_cswritep, ref_start, slen, '\t');
            if (ref_n_missing && (ref_start[0] == 'N') && (slen == 1)) {
              pvar_cswritep[-2] = '.';
            }

            // ALT
            if (n_allele == 1) {
              *pvar_cswritep++ = '.';
            } else {
              for (uint32_t allele_idx = 1; allele_idx != n_allele; ++allele_idx) {
                if (unlikely(Cswrite(&pvar_css, &pvar_cswritep))) {
                  goto BcfToPgen_ret_WRITE_FAIL;
                }
                const char* cur_alt_start;
                ScanBcfTypedString(shared_end, &parse_iter, &cur_alt_start, &slen);
                pvar_cswritep = memcpyax(pvar_cswritep, cur_alt_start, slen, ',');
              }
              --pvar_cswritep;
            }
            *pvar_cswritep++ = '\t';

            // QUAL
            if (S_CAST(int32_t, qual_bits) > 0x7f800000) {
              // NaN or missing
              // (could have fast path here for inf?)
              *pvar_cswritep++ = '.';
            } else {
              float qual_f;
              memcpy(&qual_f, &qual_bits, 4);
              pvar_cswritep = ftoa_g(qual_f, pvar_cswritep);
            }
            *pvar_cswritep++ = '\t';
            if (unlikely(Cswrite(&pvar_css, &pvar_cswritep))) {
              goto BcfToPgen_ret_WRITE_FAIL;
            }

            // FILTER
            uint32_t value_type;
            uint32_t value_ct;
            ScanBcfType(&parse_iter, &value_type, &value_ct);
            if (!value_ct) {
              *pvar_cswritep++ = '.';
            } else {
              const uint32_t value_type_m1 = value_type - 1;
              if (!value_type_m1) {
                for (uint32_t filter_idx = 0; filter_idx != value_ct; ++filter_idx) {
                  const uint32_t cur_sidx = parse_iter[filter_idx];
                  if (cur_sidx == 0x80) {
                    *pvar_cswritep++ = '.';
                  } else {
                    pvar_cswritep = memcpya(pvar_cswritep, fif_strings[cur_sidx], fif_slens[cur_sidx]);
                  }
                  *pvar_cswritep++ = ';';
                }
              } else if (value_type_m1 == 1) {
                for (uint32_t filter_idx = 0; filter_idx != value_ct; ++filter_idx) {
                  const uint32_t cur_sidx = CopyFromUnalignedOffsetU16ZX(parse_iter, filter_idx);
                  if (cur_sidx == 0x8000) {
                    *pvar_cswritep++ = '.';
                  } else {
                    pvar_cswritep = memcpya(pvar_cswritep, fif_strings[cur_sidx], fif_slens[cur_sidx]);
                  }
                  *pvar_cswritep++ = ';';
                }
              } else {
                for (uint32_t filter_idx = 0; filter_idx != value_ct; ++filter_idx) {
                  uint32_t cur_sidx;
                  CopyFromUnalignedOffsetU32(&cur_sidx, parse_iter, filter_idx);
                  if (cur_sidx == 0x80000000U) {
                    *pvar_cswritep++ = '.';
                  } else {
                    pvar_cswritep = memcpya(pvar_cswritep, fif_strings[cur_sidx], fif_slens[cur_sidx]);
                  }
                  *pvar_cswritep++ = ';';
                }
              }
              parse_iter = &(parse_iter[value_ct << value_type_m1]);
              --pvar_cswritep;
            }
            if (unlikely(Cswrite(&pvar_css, &pvar_cswritep))) {
              goto BcfToPgen_ret_WRITE_FAIL;
            }

            // INFO
            if (info_nonpr_exists) {
              *pvar_cswritep++ = '\t';
              if (!n_info) {
                *pvar_cswritep++ = '.';
              } else {
                for (uint32_t info_idx = 0; info_idx != n_info; ++info_idx) {
                  ScanBcfTypedInt(&parse_iter, &sidx);
                  pvar_cswritep = strcpya(pvar_cswritep, fif_strings[sidx]);

                  ScanBcfType(&parse_iter, &value_type, &value_ct);
                  const uint32_t bytes_per_elem = kBcfBytesPerElem[value_type];
                  const uint32_t vec_byte_ct = bytes_per_elem * value_ct;
                  const unsigned char* cur_vec_start = parse_iter;
                  parse_iter = &(parse_iter[vec_byte_ct]);
                  if (value_type == 7) {
                    // string

                    // Separated from other cases since we still write '=' when
                    // value_ct == 0.
                    *pvar_cswritep++ = '=';

                    // Unlike most other VCF/BCF fields, spaces are actually
                    // allowed by the VCF spec here, so we need to detect
                    // them.
                    // We error out on them for now (possible todo:
                    // autoconversion to "%20").
                    if (unlikely(memchr(cur_vec_start, ' ', value_ct))) {
                      snprintf(g_logbuf, kLogbufSize, "Error: INFO field in variant record #%" PRIuPTR " of --bcf file contains a space; this cannot be imported by " PROG_NAME_STR ". Remove or reformat the field before reattempting import.\n", vrec_idx);
                      goto BcfToPgen_ret_MALFORMED_INPUT_WWN;
                    }
                    pvar_cswritep = memcpya(pvar_cswritep, cur_vec_start, value_ct);
                  } else if (value_ct) {
                    *pvar_cswritep++ = '=';
                    // ugh.  not a separate function for now since no other
                    // code needs to do this
                    if (value_type == 1) {
                      // int8
                      for (uint32_t value_idx = 0; value_idx != value_ct; ++value_idx) {
                        const int8_t cur_val = S_CAST(int8_t, cur_vec_start[value_idx]);
                        if (cur_val != -128) {
                          pvar_cswritep = i32toa(cur_val, pvar_cswritep);
                        } else {
                          *pvar_cswritep++ = '.';
                        }
                        *pvar_cswritep++ = ',';
                      }
                    } else if (value_type == 2) {
                      // int16
                      for (uint32_t value_idx = 0; value_idx != value_ct; ++value_idx) {
                        int16_t cur_val;
                        CopyFromUnalignedOffsetI16(&cur_val, cur_vec_start, value_idx);
                        if (cur_val != -32768) {
                          pvar_cswritep = i32toa(cur_val, pvar_cswritep);
                        } else {
                          *pvar_cswritep++ = '.';
                        }
                        *pvar_cswritep++ = ',';
                      }
                    } else if (value_type == 3) {
                      // int32
                      for (uint32_t value_idx = 0; value_idx != value_ct; ++value_idx) {
                        int32_t cur_val;
                        CopyFromUnalignedOffsetI32(&cur_val, cur_vec_start, value_idx);
                        if (cur_val != (-2147483647 - 1)) {
                          pvar_cswritep = i32toa(cur_val, pvar_cswritep);
                        } else {
                          *pvar_cswritep++ = '.';
                        }
                        *pvar_cswritep++ = ',';
                      }
                    } else {
                      // float
                      for (uint32_t value_idx = 0; value_idx != value_ct; ++value_idx) {
                        uint32_t cur_bits;
                        CopyFromUnalignedOffsetU32(&cur_bits, cur_vec_start, value_idx);
                        if (cur_bits != 0x7f800001) {
                          float cur_float;
                          memcpy(&cur_float, &cur_bits, 4);
                          pvar_cswritep = ftoa_g(cur_float, pvar_cswritep);
                        } else {
                          *pvar_cswritep++ = '.';
                        }
                        *pvar_cswritep++ = ',';
                      }
                    }
                    --pvar_cswritep;
                  }
                  *pvar_cswritep++ = ';';
                }
                --pvar_cswritep;
              }
            }
            AppendBinaryEoln(&pvar_cswritep);
            if (unlikely(Cswrite(&pvar_css, &pvar_cswritep))) {
              goto BcfToPgen_ret_WRITE_FAIL;
            }
            if (!sample_ct) {
              if (++block_vidx == cur_thread_block_vidx_limit) {
                break;
              }
              continue;
            }
          }

          grp = &(cur_gparse[block_vidx]);
          grp->record_start = geno_buf_iter;
          grp->flags = gparse_flags;
          grp->metadata.read_bcf.qual_vec_offsets[0] = UINT32_MAX;
          grp->metadata.read_bcf.qual_vec_offsets[1] = UINT32_MAX;
          grp->metadata.read_bcf.gt_vec_offset = UINT32_MAX;
          grp->metadata.read_bcf.dosage_vec_offset = UINT32_MAX;
          grp->metadata.read_bcf.hds_vec_offset = UINT32_MAX;
          grp->metadata.read_bcf.rec_idx = vrec_idx;
          uintptr_t copy_vec_offset = 0;
          for (uint32_t qual_idx = 0; qual_idx != 2; ++qual_idx) {
            if (qual_starts[qual_idx]) {
              grp->metadata.read_bcf.qual_vec_offsets[qual_idx] = copy_vec_offset;
              const uint32_t type_blen = qual_type_blens[qual_idx];
              const unsigned char* src_iter = qual_starts[qual_idx];
              memcpy(&(geno_buf_iter[copy_vec_offset * kBytesPerVec]), src_iter, type_blen);
              src_iter = &(src_iter[type_blen]);
#ifdef __LP64__
              ++copy_vec_offset;
#else
              copy_vec_offset += DivUp(type_blen, kBytesPerVec);
#endif
              const uint32_t vec_blen = qual_main_blens[qual_idx];
              memcpy(&(geno_buf_iter[copy_vec_offset * kBytesPerVec]), src_iter, vec_blen);
              copy_vec_offset += DivUp(vec_blen, kBytesPerVec);
            }
          }
          if (gt_start) {
            grp->metadata.read_bcf.gt_vec_offset = copy_vec_offset;
            const unsigned char* src_iter = gt_start;
            memcpy(&(geno_buf_iter[copy_vec_offset * kBytesPerVec]), src_iter, gt_type_blen);
            src_iter = &(src_iter[gt_type_blen]);
#ifdef __LP64__
            ++copy_vec_offset;
#else
            copy_vec_offset += DivUp(gt_type_blen, kBytesPerVec);
#endif
            memcpy(&(geno_buf_iter[copy_vec_offset * kBytesPerVec]), src_iter, gt_main_blen);
            copy_vec_offset += DivUp(gt_main_blen, kBytesPerVec);
          }
          if (dosage_start) {
            grp->metadata.read_bcf.dosage_vec_offset = copy_vec_offset;
            const unsigned char* src_iter = dosage_start;
            memcpy(&(geno_buf_iter[copy_vec_offset * kBytesPerVec]), src_iter, dosage_type_blen);
            src_iter = &(src_iter[dosage_type_blen]);
#ifdef __LP64__
            ++copy_vec_offset;
#else
            copy_vec_offset += DivUp(dosage_type_blen, kBytesPerVec);
#endif
            memcpy(&(geno_buf_iter[copy_vec_offset * kBytesPerVec]), src_iter, dosage_main_blen);
            copy_vec_offset += DivUp(dosage_main_blen, kBytesPerVec);
          }
          if (hds_start) {
            grp->metadata.read_bcf.hds_vec_offset = copy_vec_offset;
            const unsigned char* src_iter = hds_start;
            memcpy(&(geno_buf_iter[copy_vec_offset * kBytesPerVec]), src_iter, hds_type_blen);
            src_iter = &(src_iter[hds_type_blen]);
#ifdef __LP64__
            ++copy_vec_offset;
#else
            copy_vec_offset += DivUp(hds_type_blen, kBytesPerVec);
#endif
            memcpy(&(geno_buf_iter[copy_vec_offset * kBytesPerVec]), src_iter, hds_main_blen);
            // copy_vec_offset += DivUp(hds_main_blen, kBytesPerVec);
          }
          geno_buf_iter = &(geno_buf_iter[record_byte_ct]);
          ++block_vidx;

          // true iff this is the last variant we're keeping in the entire
          // file
          if (block_vidx == block_vidx_limit) {
            for (; cur_thread_fill_idx != calc_thread_ct; ) {
              // save endpoint for current thread, and tell any leftover
              // threads to do nothing
              thread_bidxs[++cur_thread_fill_idx] = block_vidx;
            }
            break;
          }
        }
        cur_block_write_ct = block_vidx;
      }
      if (sample_ct) {
        if (vidx_start) {
          JoinThreads(&tg);
          if (unlikely(ctx.parse_failed)) {
            goto BcfToPgen_ret_THREAD_PARSE;
          }
        }
        if (!IsLastBlock(&tg)) {
          if (vidx_start + cur_block_write_ct == variant_ct) {
            DeclareLastThreadBlock(&tg);
          }
          if (unlikely(SpawnThreads(&tg))) {
            goto BcfToPgen_ret_THREAD_CREATE_FAIL;
          }
        }
        parity = 1 - parity;
        if (vidx_start) {
          // write *previous* block results
          reterr = GparseFlush(ctx.gparse[parity], allele_idx_offsets, prev_block_write_ct, &spgw);
          if (unlikely(reterr)) {
            goto BcfToPgen_ret_1;
          }
        }
      } else if (vidx_start + cur_block_write_ct == variant_ct) {
        break;
      }
      if (vidx_start == variant_ct) {
        break;
      }
      if (vidx_start) {
        printf("\r--bcf: %uk variants converted.", vidx_start / 1000);
        if (vidx_start <= main_block_size) {
          fputs("    \b\b\b\b", stdout);
        }
        fflush(stdout);
      }
      vidx_start += cur_block_write_ct;
      prev_block_write_ct = cur_block_write_ct;
    }
    if (unlikely(CswriteCloseNull(&pvar_css, pvar_cswritep))) {
      goto BcfToPgen_ret_WRITE_FAIL;
    }
    if (sample_ct) {
      reterr = SpgwFinish(&spgw);
      if (unlikely(reterr)) {
        goto BcfToPgen_ret_1;
      }
    }
    putc_unlocked('\r', stdout);
    char* write_iter = strcpya_k(g_logbuf, "--bcf: ");
    const uint32_t outname_base_slen = outname_end - outname;
    if (sample_ct) {
      write_iter = memcpya(write_iter, outname, outname_base_slen + 5);
      write_iter = strcpya_k(write_iter, " + ");
    } else {
      *pgen_generated_ptr = 0;
    }
    write_iter = memcpya(write_iter, outname, outname_base_slen);
    write_iter = strcpya_k(write_iter, ".pvar");
    if (output_zst) {
      write_iter = strcpya_k(write_iter, ".zst");
    }
    if (sample_ct && (!preexisting_psamname)) {
      write_iter = strcpya_k(write_iter, " + ");
      write_iter = memcpya(write_iter, outname, outname_base_slen);
      write_iter = strcpya_k(write_iter, ".psam");
    } else {
      *psam_generated_ptr = 0;
    }
    write_iter = strcpya_k(write_iter, " written");
    if (!sample_ct) {
      write_iter = strcpya_k(write_iter, " (no samples present)");
    }
    strcpy_k(write_iter, ".\n");
    WordWrapB(0);
    logputsb();
    if (print_splitpar_warning) {
      logerrputs("Warning: Human chrX pseudoautosomal variant(s) appear to be present in the\ninput BCF.  You probably want to include --split-par in your next command.\n");
    }
  }
  while (0) {
  BcfToPgen_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  BcfToPgen_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  BcfToPgen_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  BcfToPgen_ret_BGZF_FAIL_N:
    putc_unlocked('\n', stdout);
  BcfToPgen_ret_BGZF_FAIL:
    // ReadFail, DecompressFail, ThreadCreateFail, and RewindFail possible.
    if (reterr == kPglRetReadFail) {
      logerrprintfww(kErrprintfFread, bcfname, rstrerror(errno));
    } else if (reterr == kPglRetDecompressFail) {
      logerrprintfww(kErrprintfDecompress, bcfname, bgzf_errmsg);
    } else if (reterr == kPglRetRewindFail) {
      // forgot this case (13 Mar 2024)
      logerrprintfww(kErrprintfRewind, "--bcf file");
    }
    break;
  BcfToPgen_ret_REWIND_FAIL_N:
    putc_unlocked('\n', stdout);
  BcfToPgen_ret_REWIND_FAIL:
    logerrprintfww(kErrprintfRewind, "--bcf file");
    reterr = kPglRetRewindFail;
    break;
  BcfToPgen_ret_THREAD_PARSE:
    {
      // Doesn't really cost us anything to report the first error if multiple
      // converter threads error out in the same block, though we can't
      // generally guarantee that we'll report the first error in the file.
      for (uint32_t tidx = 0; ; ++tidx) {
        bcf_parse_err = ctx.bcf_parse_errs[tidx];
        if (bcf_parse_err) {
          vrec_idx = ctx.err_vrec_idxs[tidx];
          break;
        }
      }
    }
  BcfToPgen_ret_PARSE:
    if (bcf_parse_err == kBcfParseHalfCallError) {
      putc_unlocked('\n', stdout);
      logerrprintf("Error: Variant record #%" PRIuPTR " of --bcf file has a GT half-call.\n", vrec_idx);
      if (!half_call_explicit_error) {
        logerrputs("Use --vcf-half-call to specify how these should be processed.\n");
      }
      reterr = kPglRetMalformedInput;
      break;
    } else if (bcf_parse_err == kBcfParseInvalidDosage) {
      // probable todo: distinguish HDS errors (right now, it just prints
      // "invalid DS field" on all HDS errors).
      putc_unlocked('\n', stdout);
      logerrprintfww("Error: Variant record #%" PRIuPTR " of --bcf file has an invalid %s field.\n", vrec_idx, dosage_import_field);
      reterr = kPglRetInconsistentInput;
      break;
    } else if (bcf_parse_err == kBcfParsePolyploidError) {
      putc_unlocked('\n', stdout);
      logerrprintfww("Error: Variant record #%" PRIuPTR " of --bcf file has a polyploid genotype.%s\n", vrec_idx, (import_flags & kfImportPolyploidExplicitError)? "" : " (Use '--polyploid-mode missing' to treat these as missing values.)");
      reterr = kPglRetInconsistentInput;
    } else if (bcf_parse_err == kBcfParseWideGt) {
      putc_unlocked('\n', stdout);
      logerrprintfww("Error: Variant record #%" PRIuPTR " of --bcf file uses unexpectedly-wide integers for the GT field; only int8s are currently supported for n_allele < 64, and int16s for n_allele in [64, 16383].\n", vrec_idx);
      reterr = kPglRetNotYetSupported;
      break;
    } else if (bcf_parse_err == kBcfParseFloatDp) {
      putc_unlocked('\n', stdout);
      logerrprintfww("Error: Variant record #%" PRIuPTR " of --bcf file has a floating-point DP field; only integers are supported for now.\n", vrec_idx);
      reterr = kPglRetNotYetSupported;
      break;
    } else if (bcf_parse_err == kBcfParseNonfloatDosage) {
      putc_unlocked('\n', stdout);
      logerrprintfww("Error: Variant record #%" PRIuPTR " of --bcf file has a dosage field that isn't of Float type; this isn't currently supported.\n", vrec_idx);
      reterr = kPglRetNotYetSupported;
      break;
    }
    // kBcfParseMalformedGeneric
  BcfToPgen_ret_VREC_GENERIC:
    putc_unlocked('\n', stdout);
    logerrprintf("Error: Variant record #%" PRIuPTR " of --bcf file is malformed.\n", vrec_idx);
    reterr = kPglRetMalformedInput;
    break;
  BcfToPgen_ret_MALFORMED_INPUT_2:
    logerrputsb();
  BcfToPgen_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  BcfToPgen_ret_MALFORMED_INPUT_WWN:
    putc_unlocked('\n', stdout);
  BcfToPgen_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetMalformedInput;
    break;
  BcfToPgen_ret_MALFORMED_INPUT_GENERIC:
    logerrputs("Error: Malformed BCF file.\n");
    reterr = kPglRetMalformedInput;
    break;
  BcfToPgen_ret_MALFORMED_TEXT_HEADER:
    logerrputs("Error: Malformed BCF text header block.\n");
    reterr = kPglRetMalformedInput;
    break;
  BcfToPgen_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  BcfToPgen_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  BcfToPgen_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  BcfToPgen_ret_DEGENERATE_DATA:
    reterr = kPglRetDegenerateData;
    break;
  }
 BcfToPgen_ret_1:
  CleanupSpgw(&spgw, &reterr);
  CleanupThreads(&tg);
  CswriteCloseCond(&pvar_css, pvar_cswritep);
  CleanupBgzfRawMtStream(&bgzf);
  fclose_cond(bcffile);
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}


ENUM_U31_DEF_START()
  kOxSampleColSkip, // processed during first pass
  kOxSampleColParent,
  kOxSampleColCatUnknown,
  kOxSampleColSex,
  kOxSampleColBinary,
  kOxSampleColQuantitative,
  kOxSampleColCatNumeric,
  kOxSampleColCatNames,
ENUM_U31_DEF_END(OxSampleCol);

static_assert(5 * kMaxIdBlen <= kDecompressChunkSize, "OxSampleToPsam needs FID+IID+SID+PAT+MAT to fit in writebuf.");
PglErr OxSampleToPsam(const char* samplename, const char* const_fid, const char* ox_missing_code, const char* missing_catname, MiscFlags misc_flags, ImportFlags import_flags, uint32_t psam_01, char id_delim, char* outname, char* outname_end, uint32_t* sample_ct_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  FILE* psamfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  uintptr_t line_idx = 0;
  TextStream sample_txs;
  PreinitTextStream(&sample_txs);
  {
    uint32_t missing_catname_slen = strlen(missing_catname);

    uint32_t max_line_blen;
    if (unlikely(StandardizeMaxLineBlen(bigstack_left() / 4, &max_line_blen))) {
      goto OxSampleToPsam_ret_NOMEM;
    }
    reterr = InitTextStream(samplename, max_line_blen, 1, &sample_txs);
    if (unlikely(reterr)) {
      if (reterr == kPglRetOpenFail) {
        const uint32_t slen = strlen(samplename);
        if ((!StrEndsWith(samplename, ".sample", slen)) &&
            (!StrEndsWith(samplename, ".sample.gz", slen))) {
          logerrprintfww("Error: Failed to open %s : %s. (--sample expects a complete filename; did you forget '.sample' at the end?)\n", samplename, strerror(errno));
          goto OxSampleToPsam_ret_1;
        }
      }
      goto OxSampleToPsam_ret_TSTREAM_FAIL;
    }
    uint32_t mc_ct = 0;
    uintptr_t max_mc_blen = 1;
    char* sorted_mc = nullptr;
    if (!ox_missing_code) {
      if (unlikely(bigstack_alloc_c(3, &sorted_mc))) {
        goto OxSampleToPsam_ret_NOMEM;
      }
      strcpy_k(sorted_mc, "NA");
      mc_ct = 1;
      max_mc_blen = 3;
    } else {
      // er, this should use something like
      // CountAndMeasureMultistrReverseAlloc()...
      const char* missing_code_iter = ox_missing_code;
      while (*missing_code_iter) {
        while (*missing_code_iter == ',') {
          ++missing_code_iter;
        }
        if (!(*missing_code_iter)) {
          break;
        }
        ++mc_ct;
        const char* token_end = Strchrnul(missing_code_iter, ',');
        uintptr_t token_slen = token_end - missing_code_iter;
        if (token_slen >= max_mc_blen) {
          max_mc_blen = token_slen + 1;
        }
        missing_code_iter = token_end;
      }
      if (mc_ct) {
        if (unlikely(bigstack_alloc_c(mc_ct * max_mc_blen, &sorted_mc))) {
          goto OxSampleToPsam_ret_NOMEM;
        }
        missing_code_iter = ox_missing_code;
        for (uintptr_t mc_idx = 0; mc_idx != mc_ct; ++mc_idx) {
          while (*missing_code_iter == ',') {
            ++missing_code_iter;
          }
          const char* token_end = Strchrnul(missing_code_iter, ',');
          uintptr_t token_slen = token_end - missing_code_iter;
          memcpyx(&(sorted_mc[mc_idx * max_mc_blen]), missing_code_iter, token_slen, '\0');
          missing_code_iter = token_end;
        }
        // this was temporarily broken in June 2018 due to first
        // strcmp_overread() implementation returning 0 instead of -1 on
        // less-than
        qsort(sorted_mc, mc_ct, max_mc_blen, strcmp_overread_casted);
      }
    }

    // First pass: Validate header rows, load ID and father/mother columns if
    // present, check for all-missing phenotype columns.
    line_idx = 1;
    char* line_start = TextGet(&sample_txs);
    if (unlikely(!line_start)) {
      if (TextStreamErrcode2(&sample_txs, &reterr)) {
        goto OxSampleToPsam_ret_TSTREAM_FAIL;
      }
      logerrputs("Error: Empty .sample file.\n");
      goto OxSampleToPsam_ret_MALFORMED_INPUT;
    }
    char* token_end = CurTokenEnd(line_start);
    // Update (31 Aug 2020): First column is now treated as ID, with any column
    // name allowed, to support QCTOOLv2's .sample dialect.
    // Tolerate tab delimiter for now, though .sample spec technically
    // prohibits that.
    char* linebuf_iter = FirstNonTspace(token_end);
    uint32_t token_slen = strlen_se(linebuf_iter);
    const uint32_t id2_exists = strequal_k(linebuf_iter, "ID_2", token_slen);
    uint32_t col_ct = 1;
    if (id2_exists) {
      linebuf_iter = FirstNonTspace(&(linebuf_iter[token_slen]));
      id_delim = ' '; // this guarantees no conflict
      ++col_ct;
    }
    // 0 = not present, otherwise zero-based index; this is fine since first
    //     column has to be (part of) sample ID
    uint32_t missing_col = 0;
    uint32_t father_col = 0;
    uint32_t mother_col = 0;
    uint32_t sex_col = 0;
    while (!IsEolnKns(*linebuf_iter)) {
      token_end = CurTokenEnd(linebuf_iter);
      token_slen = token_end - linebuf_iter;
      if (token_slen == 3) {
        if (MatchUpperK(linebuf_iter, "SEX")) {
          if (unlikely(sex_col)) {
            logerrputs("Error: Multiple sex columns in .sample file.\n");
            goto OxSampleToPsam_ret_MALFORMED_INPUT;
          }
          sex_col = col_ct;
        }
      } else if (token_slen == 4) {
        if (unlikely(memequal(linebuf_iter, "ID_2", 4))) {
          if (id2_exists) {
            logerrputs("Error: Multiple 'ID_2' columns in .sample file.\n");
          } else {
            logerrputs("Error: Improperly positioned 'ID_2' column in .sample file (must be second).\n");
          }
          goto OxSampleToPsam_ret_MALFORMED_INPUT;
        }
      } else if (token_slen == 6) {
        if (MatchUpperK(linebuf_iter, "FATHER")) {
          if (unlikely(father_col)) {
            logerrputs("Error: Multiple 'father' columns in .sample file.\n");
            goto OxSampleToPsam_ret_MALFORMED_INPUT;
          }
          father_col = col_ct;
        } else if (MatchUpperK(linebuf_iter, "MOTHER")) {
          if (unlikely(mother_col)) {
            logerrputs("Error: Multiple 'mother' columns in .sample file.\n");
            goto OxSampleToPsam_ret_MALFORMED_INPUT;
          }
          mother_col = col_ct;
        }
      } else if (MatchUpperKLen(linebuf_iter, "MISSING", token_slen)) {
        if (unlikely(missing_col)) {
          logerrputs("Error: Multiple 'missing' columns in .sample file.\n");
          goto OxSampleToPsam_ret_MALFORMED_INPUT;
        }
        missing_col = col_ct;
      }
      ++col_ct;
      linebuf_iter = FirstNonTspace(token_end);
    }
    const uint32_t col_ctl = BitCtToWordCt(col_ct);
    uintptr_t* col_first_pass_remaining;
    uintptr_t* col_nm;
    unsigned char* col_types;
    // bugfix (21 Sep 2023): forgot to 0-initialize col_first_pass_remaining
    if (unlikely(bigstack_calloc_w(col_ctl, &col_first_pass_remaining) ||
                 bigstack_calloc_w(col_ctl, &col_nm) ||
                 bigstack_alloc_uc(col_ct, &col_types))) {
      goto OxSampleToPsam_ret_NOMEM;
    }
    if (col_ct > 1 + id2_exists) {
      FillBitsNz(1 + id2_exists, col_ct, col_first_pass_remaining);
    }
    ++line_idx;
    linebuf_iter = TextGet(&sample_txs);
    if (unlikely(!linebuf_iter)) {
      if (!TextStreamErrcode2(&sample_txs, &reterr)) {
        logerrputs("Error: Only one line in .sample file.\n");
        goto OxSampleToPsam_ret_MALFORMED_INPUT;
      }
      goto OxSampleToPsam_ret_TSTREAM_FAIL;
    }
    uint32_t at_least_one_binary_pheno = 0;
    uint32_t initial_skip_col_ct = 0;
    for (uint32_t col_idx = 0; col_idx != col_ct; ++col_idx) {
      linebuf_iter = FirstNonTspace(linebuf_iter);
      const unsigned char col_type_char = *linebuf_iter;
      if (unlikely(IsEolnKns(col_type_char))) {
        logerrputs("Error: Second .sample header line has fewer tokens than the first.\n");
        goto OxSampleToPsam_ret_MALFORMED_INPUT;
      }
      if (unlikely(!IsEolnKns(linebuf_iter[1]))) {
        *CurTokenEnd(linebuf_iter) = '\0';
        snprintf(g_logbuf, kLogbufSize, "Error: Unrecognized .sample variable type '%s'.\n", linebuf_iter);
        goto OxSampleToPsam_ret_MALFORMED_INPUT_WW;
      }
      // can make this table-driven
      OxSampleCol cur_col_type;
      if (col_type_char == 'D') {
        cur_col_type = kOxSampleColCatUnknown;
      } else if (col_type_char == 'B') {
        at_least_one_binary_pheno = 1;
        cur_col_type = kOxSampleColBinary;
      } else if ((col_type_char == 'C') || (col_type_char == 'P')) {
        cur_col_type = kOxSampleColQuantitative;
      } else if (likely(col_type_char == '0')) {
        ++initial_skip_col_ct;
        cur_col_type = kOxSampleColSkip;
      } else {
        snprintf(g_logbuf, kLogbufSize, "Error: Unrecognized .sample variable type '%c'.\n", col_type_char);
        goto OxSampleToPsam_ret_MALFORMED_INPUT_2;
      }
      col_types[col_idx] = cur_col_type;
      ++linebuf_iter;
    }
    linebuf_iter = FirstNonTspace(linebuf_iter);
    if (unlikely(!IsEolnKns(*linebuf_iter))) {
      logerrputs("Error: Second .sample header line has more tokens than the first.\n");
      goto OxSampleToPsam_ret_MALFORMED_INPUT;
    }
    if (unlikely(col_types[0] != kOxSampleColSkip)) {
      logerrputs("Error: First column of second .sample header line must have type '0'.\n");
      goto OxSampleToPsam_ret_MALFORMED_INPUT;
    }
    if (id2_exists) {
      if (unlikely(col_types[1] != kOxSampleColSkip)) {
        logerrputs("Error: ID_2 column in .sample file doesn't have type '0'.\n");
        goto OxSampleToPsam_ret_MALFORMED_INPUT;
      }
    }
    if (missing_col) {
      if (unlikely(col_types[missing_col] != kOxSampleColSkip)) {
        logerrputs("Error: 'missing' column in .sample file doesn't have type '0'.\n");
        goto OxSampleToPsam_ret_MALFORMED_INPUT;
      }
      ClearBit(missing_col, col_first_pass_remaining);
    }
    if (initial_skip_col_ct != 1 + id2_exists + (missing_col != 0)) {
      logerrputs("Error: Only ID and 'missing' columns are permitted to have type '0' in .sample\nfiles.\n");
      goto OxSampleToPsam_ret_MALFORMED_INPUT;
    }
    if (father_col) {
      if (unlikely(col_types[father_col] != kOxSampleColCatUnknown)) {
        logerrputs("Error: .sample 'father' column is not of type 'D'.\n");
        goto OxSampleToPsam_ret_MALFORMED_INPUT;
      }
      col_types[father_col] = kOxSampleColParent;
    }
    if (mother_col) {
      if (unlikely(col_types[mother_col] != kOxSampleColCatUnknown)) {
        logerrputs("Error: .sample 'mother' column is not of type 'D'.\n");
        goto OxSampleToPsam_ret_MALFORMED_INPUT;
      }
      col_types[mother_col] = kOxSampleColParent;
    }
    if (sex_col) {
      if (unlikely(col_types[sex_col] != kOxSampleColCatUnknown)) {
        logerrputs("Error: .sample sex column is not of type 'D'.\n");
        goto OxSampleToPsam_ret_MALFORMED_INPUT;
      }
      col_types[sex_col] = kOxSampleColSex;
      SetBit(sex_col, col_nm);
      ClearBit(sex_col, col_first_pass_remaining);
    }
    if (at_least_one_binary_pheno) {
      if (unlikely((bsearch_strbox("0", sorted_mc, 1, max_mc_blen, mc_ct) != -1) || (bsearch_strbox("1", sorted_mc, 1, max_mc_blen, mc_ct) != -1))) {
        logerrputs("Error: '0' and '1' are unacceptable missing case/control phenotype codes.\n");
        goto OxSampleToPsam_ret_INCONSISTENT_INPUT;
      }
    }
    ImportSampleIdContext isic;
    InitImportSampleIdContext(const_fid, import_flags, id_delim, &isic);
    unsigned char* tmp_alloc_base = g_bigstack_base;
    unsigned char* tmp_alloc_end = BigstackEndRoundedDown();
    char* all_ids_start = R_CAST(char*, tmp_alloc_base);
    const uint32_t parental_col_exists = father_col || mother_col;
    // doesn't include ID
    const uint32_t first_pass_main_col_ct = PopcountWords(col_first_pass_remaining, col_ctl);
    uint32_t first_pass_finished_col_ct = 0;
    char* father_start = nullptr;
    char* mother_start = nullptr;
    uint32_t father_slen = 0;
    uint32_t mother_slen = 0;
    while (1) {
      ++line_idx;
      line_start = TextGet(&sample_txs);
      if (!line_start) {
        if (likely(!TextStreamErrcode2(&sample_txs, &reterr))) {
          break;
        }
        goto OxSampleToPsam_ret_TSTREAM_FAIL;
      }
      token_end = CurTokenEnd(line_start);
      const uintptr_t col1_slen = token_end - line_start;
      if (id2_exists) {
        char* iid_start = FirstNonTspace(token_end);
        if (unlikely(IsEolnKns(*iid_start))) {
          goto OxSampleToPsam_ret_MISSING_TOKENS;
        }
        char* iid_end = CurTokenEnd(iid_start);
        const uintptr_t col2_slen = iid_end - iid_start;
        if (unlikely(col1_slen + 1 + col2_slen >= S_CAST(uintptr_t, tmp_alloc_end - tmp_alloc_base))) {
          goto OxSampleToPsam_ret_NOMEM;
        }
        tmp_alloc_base = memcpyua(tmp_alloc_base, line_start, col1_slen);
        *tmp_alloc_base++ = ' ';
        tmp_alloc_base = memcpyua(tmp_alloc_base, iid_start, col2_slen);
        token_end = iid_end;
      } else {
        if (unlikely(col1_slen >= S_CAST(uintptr_t, tmp_alloc_end - tmp_alloc_base))) {
          goto OxSampleToPsam_ret_NOMEM;
        }
        tmp_alloc_base = memcpyua(tmp_alloc_base, line_start, col1_slen);
      }
      *tmp_alloc_base++ = '\0';
      uintptr_t col_uidx_base = 0;
      uintptr_t cur_bits = col_first_pass_remaining[0];
      uint32_t prev_col_uidx = id2_exists;
      for (uint32_t uii = first_pass_finished_col_ct; uii != first_pass_main_col_ct; ++uii) {
        const uint32_t col_uidx = BitIter1(col_first_pass_remaining, &col_uidx_base, &cur_bits);
        linebuf_iter = NextTokenMult(token_end, col_uidx - prev_col_uidx);
        if (unlikely(!linebuf_iter)) {
          goto OxSampleToPsam_ret_MISSING_TOKENS;
        }
        prev_col_uidx = col_uidx;
        token_end = CurTokenEnd(linebuf_iter);
        token_slen = token_end - linebuf_iter;
        const OxSampleCol cur_col_type = S_CAST(OxSampleCol, col_types[col_uidx]);
        if (bsearch_strbox(linebuf_iter, sorted_mc, token_slen, max_mc_blen, mc_ct) != -1) {
          // Missing value.
          if (cur_col_type == kOxSampleColParent) {
            if (col_uidx == father_col) {
              father_start = nullptr;
            } else {
              mother_start = nullptr;
            }
          }
          continue;
        }
        if (cur_col_type == kOxSampleColCatUnknown) {
          if (IsCategoricalPhenostrNocsv(linebuf_iter)) {
            col_types[col_uidx] = kOxSampleColCatNames;
          } else {
            // Assume this is a column of positive integers.  If this
            // assumption is false, we'll detect the problem during the second
            // pass.
            col_types[col_uidx] = kOxSampleColCatNumeric;
          }
          SetBit(col_uidx, col_nm);
          ClearBit(col_uidx, col_first_pass_remaining);
          ++first_pass_finished_col_ct;
        } else if ((cur_col_type == kOxSampleColBinary) || (cur_col_type == kOxSampleColQuantitative)) {
          SetBit(col_uidx, col_nm);
          ClearBit(col_uidx, col_first_pass_remaining);
          ++first_pass_finished_col_ct;
        } else if (cur_col_type == kOxSampleColParent) {
          if (col_uidx == father_col) {
            father_start = linebuf_iter;
            father_slen = token_end - linebuf_iter;
          } else {
            mother_start = linebuf_iter;
            mother_slen = token_end - linebuf_iter;
          }
        }
      }
      if (parental_col_exists) {
        if (father_start) {
          tmp_alloc_base = memcpyua(tmp_alloc_base, father_start, father_slen);
          *tmp_alloc_base++ = '\0';
        } else {
          tmp_alloc_base = memcpyua(tmp_alloc_base, "0", 2);
        }
        if (mother_start) {
          tmp_alloc_base = memcpyua(tmp_alloc_base, mother_start, mother_slen);
          *tmp_alloc_base++ = '\0';
        } else {
          tmp_alloc_base = memcpyua(tmp_alloc_base, "0", 2);
        }
      }
    }
    BigstackBaseSet(tmp_alloc_base);
    if (unlikely(line_idx > 0x80000001U)) {
      logerrputs("Error: " PROG_NAME_STR " does not support more than 2^31 - 2 samples.\n");
      goto OxSampleToPsam_ret_MALFORMED_INPUT;
    }
    const uint32_t sample_ct = line_idx - 3;
    if (unlikely(!sample_ct)) {
      logerrputs("Error: No samples in .sample file.\n");
      goto OxSampleToPsam_ret_DEGENERATE_DATA;
    }
    const char* all_ids_iter = all_ids_start;
    const uint32_t iid_sid = (misc_flags / kfMiscIidSid) & 1;
    uint32_t write_fid = 0;
    uint32_t write_sid = 0;
    uint32_t nz_parent_exists = 0;
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      if (id_delim) {
        if (id2_exists) {
          // fast path for original case
          if (!write_fid) {
            write_fid = !memequal(all_ids_iter, "0 ", 2);
          }
          all_ids_iter = strnul(all_ids_iter);
        } else {
          const char* first_part_end = Strchrnul(all_ids_iter, id_delim);
          if (unlikely(!(*first_part_end))) {
            logerrprintfww("Error: No instances of --id-delim argument '%c' in sample ID '%s'.\n", id_delim, all_ids_iter);
            goto OxSampleToPsam_ret_INCONSISTENT_INPUT;
          }
          const uint32_t first_slen = first_part_end - all_ids_iter;
          const char* second_part_start = &(first_part_end[1]);
          const char* second_part_end = Strchrnul(second_part_start, id_delim);
          if (*second_part_end) {
            const char* third_part_start = &(second_part_end[1]);
            const char* third_part_end = Strchrnul(third_part_start, id_delim);
            if (unlikely(*third_part_end)) {
              logerrprintfww("Error: Too many instances of --id-delim argument '%c' in sample ID '%s'.\n", id_delim, all_ids_iter);
              goto OxSampleToPsam_ret_INCONSISTENT_INPUT;
            }
            write_fid |= (first_slen != 1) || (*all_ids_iter != '0');
            write_sid |= (S_CAST(uintptr_t, third_part_end - third_part_start) != 1) || (*third_part_start != '0');
            all_ids_iter = third_part_end;
          } else {
            if (!iid_sid) {
              write_fid |= (first_slen != 1) || (*all_ids_iter != '0');
            } else {
              write_sid |= (S_CAST(uintptr_t, second_part_end - second_part_start) != 1) || (*second_part_start != '0');
            }
            all_ids_iter = second_part_end;
          }
        }
      } else {
        all_ids_iter = strnul(all_ids_iter);
      }
      ++all_ids_iter;
      if (parental_col_exists) {
        const char* paternal_id_end = strnul(all_ids_iter);
        const char* maternal_id_start = &(paternal_id_end[1]);
        const char* maternal_id_end = strnul(maternal_id_start);
        if (!nz_parent_exists) {
          // octal
          nz_parent_exists = !memequal(all_ids_iter, "\60\0\60", 4);
        }
        all_ids_iter = &(maternal_id_end[1]);
      }
    }
    if (id_delim) {
      if (!iid_sid) {
        isic.fid_delim_mode = write_fid? kImportFidDelimModeAlwaysCopy : kImportFidDelimModeAlwaysOmit;
        isic.sid_delim_mode = write_sid? kImportSidDelimModeCopyOr0 : kImportSidDelimModeNonexistOrOmit;
      } else {
        isic.fid_delim_mode = write_fid? kImportFidDelimModeCopyOr0 : kImportFidDelimModeOmitWhenPresent;
        isic.sid_delim_mode = write_sid? kImportSidDelimModeAlwaysCopy : kImportSidDelimModeNonexistOrOmit;
      }
    } else if (isic.const_fid || isic.double_id) {
      write_fid = 1;
    }

    snprintf(outname_end, kMaxOutfnameExtBlen, ".psam");
    if (unlikely(fopen_checked(outname, FOPEN_WB, &psamfile))) {
      goto OxSampleToPsam_ret_OPEN_FAIL;
    }
    // We add col_ct, since positive-integer categorical phenotypes are
    // lengthened by 1 character ('C' added in front).  (This bound can easily
    // be tightened, but it should practically never matter.)
    uintptr_t linebuf_size = max_line_blen + col_ct + kDecompressChunkSize;
    linebuf_size = MINV(linebuf_size, kMaxLongLine);
    char* writebuf;
    if (unlikely(bigstack_alloc_c(linebuf_size, &writebuf))) {
      goto OxSampleToPsam_ret_NOMEM;
    }
    char* write_iter = writebuf;
    *write_iter++ = '#';
    if (write_fid) {
      write_iter = strcpya_k(write_iter, "FID\t");
    }
    write_iter = strcpya_k(write_iter, "IID");
    if (write_sid) {
      write_iter = strcpya_k(write_iter, "\tSID");
    }
    if (nz_parent_exists) {
      write_iter = strcpya_k(write_iter, "\tPAT\tMAT");
    }
    write_iter = strcpya_k(write_iter, "\tSEX");

    reterr = TextRewind(&sample_txs);
    if (unlikely(reterr)) {
      goto OxSampleToPsam_ret_TSTREAM_FAIL;
    }
    line_idx = 1;
    line_start = TextGet(&sample_txs);
    if (unlikely(!line_start)) {
      goto OxSampleToPsam_ret_TSTREAM_REWIND_FAIL;
    }
    token_end = line_start;
    for (uint32_t col_idx = 0; col_idx != col_ct; ++col_idx) {
      linebuf_iter = FirstNonTspace(token_end);
      if (unlikely(IsEolnKns(*linebuf_iter))) {
        goto OxSampleToPsam_ret_REWIND_FAIL;
      }
      token_end = CurTokenEnd(linebuf_iter);
      const OxSampleCol cur_col_type = S_CAST(OxSampleCol, col_types[col_idx]);
      if (cur_col_type > kOxSampleColSex) {
        if (IsSet(col_nm, col_idx)) {
          *write_iter++ = '\t';
          write_iter = memcpya(write_iter, linebuf_iter, token_end - linebuf_iter);
        }
      }
    }
    AppendBinaryEoln(&write_iter);
    if (unlikely(fwrite_checked(writebuf, write_iter - writebuf, psamfile))) {
      goto OxSampleToPsam_ret_WRITE_FAIL;
    }
    ++line_idx;
    line_start = TextGet(&sample_txs);
    if (unlikely(!line_start)) {
      goto OxSampleToPsam_ret_TSTREAM_REWIND_FAIL;
    }
    all_ids_iter = all_ids_start;
    const char ctrl_outchar = '1' - psam_01;
    const char case_outchar = ctrl_outchar + 1;
    line_start = TextGet(&sample_txs);
    for (; line_start; line_start = TextGet(&sample_txs)) {
      ++line_idx;
      write_iter = writebuf;
      const char* sample_id_end = strnul(all_ids_iter);
      reterr = ImportSampleId(all_ids_iter, sample_id_end, &isic, &write_iter);
      if (unlikely(reterr)) {
        goto OxSampleToPsam_ret_1;
      }
      all_ids_iter = &(sample_id_end[1]);
      if (parental_col_exists) {
        const char* paternal_id_end = strnul(all_ids_iter);
        if (nz_parent_exists) {
          *write_iter++ = '\t';
          write_iter = memcpya(write_iter, all_ids_iter, paternal_id_end - all_ids_iter);
        }
        all_ids_iter = &(paternal_id_end[1]);
        const char* maternal_id_end = strnul(all_ids_iter);
        if (nz_parent_exists) {
          *write_iter++ = '\t';
          write_iter = memcpya(write_iter, all_ids_iter, maternal_id_end - all_ids_iter);
        }
        all_ids_iter = &(maternal_id_end[1]);
      }
      // flush now since backfilled sex is variable-length ("NA" vs. "1"/"2")
      if (unlikely(fwrite_checked(writebuf, write_iter - writebuf, psamfile))) {
        goto OxSampleToPsam_ret_WRITE_FAIL;
      }

      char* cur_writebuf_start = writebuf;
      write_iter = strcpya_k(writebuf, "\tNA");
      token_end = line_start;
      for (uint32_t col_idx = 0; col_idx != col_ct; ++col_idx) {
        linebuf_iter = FirstNonTspace(token_end);
        if (unlikely(IsEolnKns(*linebuf_iter))) {
          goto OxSampleToPsam_ret_MISSING_TOKENS;
        }
        token_end = CurTokenEnd(linebuf_iter);
        const OxSampleCol cur_col_type = S_CAST(OxSampleCol, col_types[col_idx]);
        if ((cur_col_type < kOxSampleColSex) || (!IsSet(col_nm, col_idx))) {
          continue;
        }
        token_slen = token_end - linebuf_iter;
        const uint32_t is_missing = (bsearch_strbox(linebuf_iter, sorted_mc, token_slen, max_mc_blen, mc_ct) != -1);
        if (cur_col_type >= kOxSampleColCatNumeric) {
          *write_iter++ = '\t';
          if (!is_missing) {
            double dxx = 0.0;
            char* num_end = ScanadvDouble(linebuf_iter, &dxx);
            if (cur_col_type == kOxSampleColCatNumeric) {
              *write_iter++ = 'C';
              // .sample files are relatively small, so let's go ahead and
              // (i) validate we have a positive integer < 2^31
              // (ii) convert e.g. 9000000, 9000000., 9.0e6 all to 9000000
              int32_t ii = S_CAST(int32_t, dxx);
              if (unlikely((num_end != token_end) || (ii <= 0) || (S_CAST(double, ii) != dxx))) {
                *token_end = '\0';
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid categorical phenotype value '%s' on line %" PRIuPTR ", column %u of .sample file (positive integer < 2^31 or --missing-code value expected).\n", linebuf_iter, line_idx, col_idx + 1);
                goto OxSampleToPsam_ret_INCONSISTENT_INPUT_WW;
              }
              write_iter = u32toa(ii, write_iter);
            } else {
              if (unlikely(num_end)) {
                logerrputs("Error: Mixed numeric/non-numeric categorical phenotype in .sample file.\n");
                goto OxSampleToPsam_ret_MALFORMED_INPUT;
              }
              write_iter = memcpya(write_iter, linebuf_iter, token_slen);
            }
          } else {
            write_iter = memcpya(write_iter, missing_catname, missing_catname_slen);
          }
        } else if (cur_col_type != kOxSampleColSex) {
          *write_iter++ = '\t';
          if (!is_missing) {
            if (cur_col_type == kOxSampleColBinary) {
              // tolerate "control"/"case" as well as 0/1
              if (strequal_k(linebuf_iter, "0", token_slen) || strequal_k(linebuf_iter, "control", token_slen)) {
                *write_iter++ = ctrl_outchar;
              } else if (likely(strequal_k(linebuf_iter, "1", token_slen) || strequal_k(linebuf_iter, "case", token_slen))) {
                *write_iter++ = case_outchar;
              } else {
                *token_end = '\0';
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid binary phenotype value '%s' on line %" PRIuPTR ", column %u of .sample file ('0', '1', or --missing-code value expected).\n", linebuf_iter, line_idx, col_idx + 1);
                goto OxSampleToPsam_ret_INCONSISTENT_INPUT_WW;
              }
            } else {
              assert(cur_col_type == kOxSampleColQuantitative);
              double dxx = 0.0;
              if (unlikely(!ScantokDouble(linebuf_iter, &dxx))) {
                *token_end = '\0';
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid quantitative phenotype value '%s' on line %" PRIuPTR ", column %u of .sample file (non-infinite number or --missing-code value expected).\n", linebuf_iter, line_idx, col_idx + 1);
                goto OxSampleToPsam_ret_INCONSISTENT_INPUT_WW;
              }
              // used over memcpy to make --data and --data --make-pgen the
              // same (could make that conditional on keep_autoconv?)
              write_iter = dtoa_g(dxx, write_iter);
            }
          } else {
            // bugfix (7 Jul 2023): forgot this branch
            write_iter = strcpya_k(write_iter, "NA");
          }
        } else {
          // sex
          if (!is_missing) {
            const unsigned char sex_ucc = *linebuf_iter;
            if ((token_slen == 1) && ((S_CAST(uint32_t, sex_ucc) - 49) < 2)) {
              ++cur_writebuf_start;
              cur_writebuf_start[0] = '\t';
              cur_writebuf_start[1] = sex_ucc;
            } else if (unlikely((token_slen != 1) || (sex_ucc != '0'))) {
              // tolerate '0' as a sex-only missing code even when not
              // explicitly specified
              *token_end = '\0';
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid sex code '%s' on line %" PRIuPTR ", column %u of .sample file ('0', '1', '2', or --missing-code value expected).\n", linebuf_iter, line_idx, col_idx + 1);
              goto OxSampleToPsam_ret_INCONSISTENT_INPUT_WW;
            }
          }
        }
      }
      AppendBinaryEoln(&write_iter);
      if (unlikely(fwrite_checked(cur_writebuf_start, write_iter - cur_writebuf_start, psamfile))) {
        goto OxSampleToPsam_ret_WRITE_FAIL;
      }
    }
    if (unlikely(TextStreamErrcode2(&sample_txs, &reterr) ||
                 (line_idx != sample_ct + 2))) {
      goto OxSampleToPsam_ret_TSTREAM_REWIND_FAIL;
    }

    // no final writebuf flush since we didn't use usual manual-streaming
    // strategy
    // (probably change this?)
    if (unlikely(fclose_null(&psamfile))) {
      goto OxSampleToPsam_ret_WRITE_FAIL;
    }
    logprintfww("%u sample%s imported from .sample file to %s .\n", sample_ct, (sample_ct == 1)? "" : "s", outname);
    *sample_ct_ptr = sample_ct;
  }
  while (0) {
  OxSampleToPsam_ret_TSTREAM_FAIL:
    TextStreamErrPrint(".sample file", &sample_txs);
    break;
  OxSampleToPsam_ret_TSTREAM_REWIND_FAIL:
    TextStreamErrPrintRewind(".sample file", &sample_txs, &reterr);
    break;
  OxSampleToPsam_ret_REWIND_FAIL:
    logerrprintfww(kErrprintfRewind, ".sample file");
    reterr = kPglRetRewindFail;
    break;
  OxSampleToPsam_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  OxSampleToPsam_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  OxSampleToPsam_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  OxSampleToPsam_ret_MISSING_TOKENS:
    logerrprintf("Error: Line %" PRIuPTR " of .sample file has fewer tokens than expected.\n", line_idx);
    reterr = kPglRetMalformedInput;
    break;
  OxSampleToPsam_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
  OxSampleToPsam_ret_MALFORMED_INPUT_2:
    logerrputsb();
  OxSampleToPsam_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  OxSampleToPsam_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  OxSampleToPsam_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  OxSampleToPsam_ret_DEGENERATE_DATA:
    reterr = kPglRetDegenerateData;
    break;
  }
 OxSampleToPsam_ret_1:
  CleanupTextStream2(".sample file", &sample_txs, &reterr);
  fclose_cond(psamfile);
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

uint32_t Bgen11DosageImportCheck(uint32_t dosage_int_sum_thresh, uint32_t import_dosage_certainty_int, uint32_t dosage_erase_halfdist, uint32_t dosage_int0, uint32_t dosage_int1, uint32_t dosage_int2) {
  const uint32_t dosage_int_sum = dosage_int0 + dosage_int1 + dosage_int2;
  if ((dosage_int_sum <= dosage_int_sum_thresh) && (dosage_int0 < import_dosage_certainty_int) && (dosage_int1 < import_dosage_certainty_int) && (dosage_int2 < import_dosage_certainty_int)) {
    return 0;
  }
  // ties realistically happen, use banker's rounding
  // 1/65536 -> 0/32768
  // 3/65536 -> 2/32768
  // 5/65536 -> 2/32768
  const DosageProd write_dosage_int_numer = S_CAST(DosageProd, kDosageMid) * dosage_int1 + S_CAST(DosageProd, kDosageMax) * dosage_int2;
  uint32_t write_dosage_int;
  if (dosage_int_sum == kDosageMax) {
    // optimize common case
    write_dosage_int = ((write_dosage_int_numer + kDosageMid) / kDosageMax) - ((write_dosage_int_numer % (2 * S_CAST(DosageProd, kDosageMax))) == kDosageMid);
  } else {
    write_dosage_int = (write_dosage_int_numer + (dosage_int_sum / 2)) / dosage_int_sum;
    write_dosage_int -= (2 * (write_dosage_int_numer - write_dosage_int * dosage_int_sum) == dosage_int_sum) * (write_dosage_int % 2);
  }
  const uint32_t halfdist = BiallelicDosageHalfdist(write_dosage_int);
  return (halfdist < dosage_erase_halfdist);
}

void Bgen11DosageImportUpdate(uint32_t dosage_int_sum_thresh, uint32_t import_dosage_certainty_int, uint32_t hard_call_halfdist, uint32_t dosage_erase_halfdist, uint32_t sample_idx_lowbits, uint32_t dosage_int0, uint32_t dosage_int1, uint32_t dosage_int2, uintptr_t* genovec_word_ptr, uint32_t* dosage_present_hw_ptr, Dosage** dosage_main_iterp) {
  const uint32_t dosage_int_sum = dosage_int0 + dosage_int1 + dosage_int2;
  if (dosage_int_sum <= dosage_int_sum_thresh) {
    if ((dosage_int0 < import_dosage_certainty_int) && (dosage_int1 < import_dosage_certainty_int) && (dosage_int2 < import_dosage_certainty_int)) {
      *genovec_word_ptr |= (3 * k1LU) << (2 * sample_idx_lowbits);
      return;
    }
  }
  const DosageProd write_dosage_int_numer = S_CAST(DosageProd, kDosageMid) * dosage_int1 + S_CAST(DosageProd, kDosageMax) * dosage_int2;
  uint32_t write_dosage_int;
  if (dosage_int_sum == kDosageMax) {
    write_dosage_int = ((write_dosage_int_numer + kDosageMid) / kDosageMax) - ((write_dosage_int_numer % (2 * S_CAST(DosageProd, kDosageMax))) == kDosageMid);
  } else {
    write_dosage_int = (write_dosage_int_numer + (dosage_int_sum / 2)) / dosage_int_sum;
    write_dosage_int -= (2 * (write_dosage_int_numer - write_dosage_int * dosage_int_sum) == dosage_int_sum) * (write_dosage_int % 2);
  }
  const uint32_t halfdist = BiallelicDosageHalfdist(write_dosage_int);
  if (halfdist < hard_call_halfdist) {
    *genovec_word_ptr |= (3 * k1LU) << (2 * sample_idx_lowbits);
  } else {
    *genovec_word_ptr |= ((write_dosage_int + (kDosage4th * k1LU)) / kDosageMid) << (2 * sample_idx_lowbits);
    if (halfdist >= dosage_erase_halfdist) {
      return;
    }
  }
  *dosage_present_hw_ptr |= 1U << sample_idx_lowbits;
  **dosage_main_iterp = write_dosage_int;
  *dosage_main_iterp += 1;
}

// This may make a (tiny) allocation from the bottom of bigstack.
// single_chr_str_ptr and single_chr_slen_ptr must be nullptr or non-null
// together.
// Ok for file_descrip to be nullptr if cur_chr_code_ptr is nullptr.
PglErr InitOxfordSingleChr(const char* ox_single_chr_str, const char* file_descrip, const char** single_chr_str_ptr, uint32_t* single_chr_slen_ptr, uint32_t* cur_chr_code_ptr, ChrInfo* cip) {
  const uint32_t chr_code_raw = GetChrCodeRaw(ox_single_chr_str);
  if (chr_code_raw == UINT32_MAX) {
    // command-line parser guarantees that prohibit_extra_chr is false here
    const uint32_t chr_slen = strlen(ox_single_chr_str);
    if (single_chr_str_ptr) {
      *single_chr_str_ptr = ox_single_chr_str;
      *single_chr_slen_ptr = chr_slen;
    }
    PglErr reterr = kPglRetSuccess;
    if (cur_chr_code_ptr) {
      if (TryToAddChrName(ox_single_chr_str, file_descrip, 0, chr_slen, 0, cur_chr_code_ptr, cip) != kPglRetSuccess) {
        reterr = kPglRetInvalidCmdline;
      }
    }
    return reterr;
  }
  uint32_t chr_code = chr_code_raw;
  if (chr_code > cip->max_code) {
    if (chr_code < kMaxContigs) {
      logerrputs("Error: --oxford-single-chr chromosome code is not in the chromosome set.\n");
      return kPglRetInvalidCmdline;
    }
    chr_code = cip->xymt_codes[chr_code - kMaxContigs];
    if (IsI32Neg(chr_code)) {
      logerrputs("Error: --oxford-single-chr chromosome code is not in the chromosome set.\n");
      return kPglRetInvalidCmdline;
    }
  }
  if (!IsSet(cip->chr_mask, chr_code)) {
    logerrputs("Error: --oxford-single-chr chromosome code is excluded by chromosome filter.\n");
    return kPglRetInvalidCmdline;
  }
  if (cur_chr_code_ptr) {
    *cur_chr_code_ptr = chr_code;
  }
  if (single_chr_str_ptr) {
    char* chr_buf;
    if (unlikely(bigstack_alloc_c(kCacheline, &chr_buf))) {
      return kPglRetNomem;
    }
    char* chr_name_end = chrtoa(cip, chr_code, chr_buf);
    *single_chr_str_ptr = chr_buf;
    *single_chr_slen_ptr = chr_name_end - chr_buf;
  }
  return kPglRetSuccess;
}

static_assert(sizeof(Dosage) == 2, "OxGenToPgen() needs to be updated.");
PglErr OxGenToPgen(const char* genname, const char* samplename, const char* const_fid, const char* ox_single_chr_str, const char* ox_missing_code, const char* missing_catname, MiscFlags misc_flags, ImportFlags import_flags, LoadFilterLogFlags load_filter_log_import_flags, OxfordImportFlags oxford_import_flags, uint32_t psam_01, uint32_t is_splitpar, uint32_t is_sortvars, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, double import_dosage_certainty, char id_delim, uint32_t max_thread_ct, char* outname, char* outname_end, ChrInfo* cip) {
  // Experimented with making this single-pass; that actually benchmarked a bit
  // slower than the current 2-pass algorithm.  Yes, that might be because the
  // Mac I tested on has a crappy filesystem, but it's still enough reason to
  // leave this alone.
  unsigned char* bigstack_mark = g_bigstack_base;
  char* pvar_cswritep = nullptr;
  PglErr reterr = kPglRetSuccess;
  uintptr_t line_idx = 0;
  CompressStreamState pvar_css;
  PreinitCstream(&pvar_css);
  TextStream gen_txs;
  STPgenWriter spgw;
  PreinitTextStream(&gen_txs);
  PreinitSpgw(&spgw);
  {
    uint32_t sample_ct;
    reterr = OxSampleToPsam(samplename, const_fid, ox_missing_code, missing_catname, misc_flags, import_flags, psam_01, id_delim, outname, outname_end, &sample_ct);
    if (unlikely(reterr)) {
      goto OxGenToPgen_ret_1;
    }
    if (unlikely(sample_ct > (kMaxLongLine / 6))) {
      // impossible for a valid .gen line to fit in maximum-length load buffer
      logerrputs("Error: Too many samples for .gen file converter.\n");
      reterr = kPglRetNotYetSupported;
      goto OxGenToPgen_ret_1;
    }
    // Two passes:
    // 1. Count # of (non-chromosome-filtered) variants, write .pvar file, and
    //    check if at least one non-hardcall needs to be saved.
    // 2. Write the .pgen.
    uint32_t max_line_blen;
    if (unlikely(StandardizeMaxLineBlen(bigstack_left() / 4, &max_line_blen))) {
      goto OxGenToPgen_ret_NOMEM;
    }
    const uint32_t decompress_thread_ct = 1 + (max_thread_ct > 2);
    reterr = ForceNonFifo(genname);
    if (unlikely(reterr)) {
      if (reterr == kPglRetOpenFail) {
        const uint32_t slen = strlen(genname);
        if ((!StrEndsWith(genname, ".gen", slen)) &&
            (!StrEndsWith(genname, ".gen.gz", slen))) {
          logerrprintfww("Error: Failed to open %s : %s. (--gen expects a complete filename; did you forget '.gen' at the end?)\n", genname, strerror(errno));
        } else {
          logerrprintfww(kErrprintfFopen, genname, strerror(errno));
        }
      } else {
        logerrprintfww(kErrprintfRewind, genname);
      }
      goto OxGenToPgen_ret_1;
    }
    reterr = InitTextStream(genname, max_line_blen, decompress_thread_ct, &gen_txs);
    if (unlikely(reterr)) {
      goto OxGenToPgen_ret_TSTREAM_FAIL;
    }
    const uint32_t prohibit_extra_chr = (misc_flags / kfMiscProhibitExtraChr) & 1;
    FinalizeChrset(load_filter_log_import_flags, cip);

    snprintf(outname_end, kMaxOutfnameExtBlen, ".pvar");
    const uint32_t output_zst = (import_flags / kfImportKeepAutoconvVzs) & 1;
    if (output_zst) {
      snprintf(&(outname_end[5]), kMaxOutfnameExtBlen - 5, ".zst");
    }
    const uintptr_t overflow_buf_size = kCompressStreamBlock + max_line_blen;
    unsigned char* bigstack_mark2 = g_bigstack_base;
    reterr = InitCstreamAlloc(outname, 0, output_zst, 1, overflow_buf_size, &pvar_css, &pvar_cswritep);
    if (unlikely(reterr)) {
      goto OxGenToPgen_ret_1;
    }

    const char* single_chr_str = nullptr;
    uint32_t single_chr_slen = 0;
    uint32_t cur_chr_code = 0;
    if (ox_single_chr_str) {
      reterr = InitOxfordSingleChr(ox_single_chr_str, ".gen file", &single_chr_str, &single_chr_slen, &cur_chr_code, cip);
      if (unlikely(reterr)) {
        goto OxGenToPgen_ret_1;
      }
    }

    if (cip->chrset_source) {
      AppendChrsetLine(cip, &pvar_cswritep);
    }
    pvar_cswritep = strcpya_k(pvar_cswritep, "#CHROM\tPOS\tID\tREF\tALT" EOLN_STR);

    // Explicit 32768 instead of kDosageMax since this is driven by the BGEN
    // 1.1 format, not plink2's dosage representation.
    // Note that command-line parser multiplies import_dosage_certainty by
    // (1 - kSmallEpsilon), and we want import_dosage_certainty_int to be 1
    // when import_dosage_certainty is zero.
    uint32_t import_dosage_certainty_int = 1 + S_CAST(int32_t, import_dosage_certainty * 32768);
    const uint32_t dosage_int_sum_thresh = 3 * (import_dosage_certainty_int - 1);

    const uint32_t dosage_erase_halfdist = kDosage4th - dosage_erase_thresh;
    const uint32_t prov_ref_allele_second = !(oxford_import_flags & kfOxfordImportRefFirst);
    uint32_t dosage_exists = 0;
    uint32_t variant_ct = 0;
    uintptr_t variant_skip_ct = 0;
    char* line_iter = TextLineEnd(&gen_txs);
    ++line_idx;
    if (unlikely(!TextGetUnsafe2(&gen_txs, &line_iter))) {
      if (TextStreamErrcode2(&gen_txs, &reterr)) {
        goto OxGenToPgen_ret_TSTREAM_FAIL;
      }
      logerrputs("Error: Empty .gen file.\n");
      goto OxGenToPgen_ret_DEGENERATE_DATA;
    }
    uint32_t is_v2 = 0;
    {
      const uint32_t token_ct = CountTokens(line_iter);
      const uint32_t expected_v2_token_ct = 3 * sample_ct + 6;
      if (token_ct == expected_v2_token_ct) {
        is_v2 = 1;
      } else if (unlikely(token_ct != expected_v2_token_ct - 1)) {
        logerrprintfww("Error: Unexpected number of columns in .gen file (%u or %u expected).\n", expected_v2_token_ct - 1, expected_v2_token_ct);
        goto OxGenToPgen_ret_INCONSISTENT_INPUT;
      }
    }
    // sex_info_avail guaranteed to be true since .sample is required
    uint32_t par_warn_code = UINT32_MAX;
    uint32_t print_splitpar_warning = 0;
    if (!(import_flags & kfImportLaxChrX)) {
      const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
      if ((!IsI32Neg(x_code)) && IsHumanChrset(cip) && (!is_splitpar)) {
        par_warn_code = x_code;
      }
    }
    goto OxGenToPgen_loop_start;
    for (; TextGetUnsafe2(&gen_txs, &line_iter); ++line_idx) {
    OxGenToPgen_loop_start:
      ;  // make this work with plain C, as opposed to just C++, compilers
      char* chr_code_str = line_iter;
      char* chr_code_end = CurTokenEnd(chr_code_str);
      const char* variant_id_str = FirstNonTspace(chr_code_end);
      if (unlikely(IsEolnKns(*variant_id_str))) {
        goto OxGenToPgen_ret_MISSING_TOKENS;
      }
      if (is_v2) {
        variant_id_str = FirstNonTspace(CurTokenEnd(variant_id_str));
        if (unlikely(IsEolnKns(*variant_id_str))) {
          goto OxGenToPgen_ret_MISSING_TOKENS;
        }
      }

      if (!single_chr_str) {
        reterr = GetOrAddChrCodeDestructive(".gen file", line_idx, prohibit_extra_chr, chr_code_str, chr_code_end, cip, &cur_chr_code);
        if (unlikely(reterr)) {
          if (strequal_k(chr_code_str, "---", chr_code_end - chr_code_str)) {
            logerrputs("(Did you forget --oxford-single-chr?)\n");
          }
          goto OxGenToPgen_ret_1;
        }
        if (!IsSet(cip->chr_mask, cur_chr_code)) {
          ++variant_skip_ct;
          continue;
        }
        pvar_cswritep = chrtoa(cip, cur_chr_code, pvar_cswritep);
      } else {
        pvar_cswritep = memcpya(pvar_cswritep, single_chr_str, single_chr_slen);
      }
      *pvar_cswritep++ = '\t';
      ++variant_ct;

      const char* variant_id_end = CurTokenEnd(variant_id_str);
      const char* pos_str = FirstNonTspace(variant_id_end);
      if (unlikely(IsEolnKns(*pos_str))) {
        goto OxGenToPgen_ret_MISSING_TOKENS;
      }
      const char* pos_end = CurTokenEnd(pos_str);
      uint32_t cur_bp;
      if (unlikely(ScanUintDefcap(pos_str, &cur_bp))) {
        putc_unlocked('\n', stdout);
        snprintf(g_logbuf, kLogbufSize, "Error: Invalid bp coordinate on line %" PRIuPTR " of %s.\n", line_idx, genname);
        goto OxGenToPgen_ret_MALFORMED_INPUT_WW;
      }

      pvar_cswritep = u32toa_x(cur_bp, '\t', pvar_cswritep);
      if (par_warn_code == cur_chr_code) {
        if ((cur_bp <= kPAR1IntersectionLast) || (cur_bp >= kPAR2IntersectionFirst)) {
          if (unlikely(!is_sortvars)) {
            putc_unlocked('\n', stdout);
            logerrputs("Error: Human chrX pseudoautosomal variant(s) appear to be present in the input\n.gen, but --split-par was not specified.\n");
            goto OxGenToPgen_ret_INCONSISTENT_INPUT;
          }
          print_splitpar_warning = 1;
          par_warn_code = UINT32_MAX;
        }
      }
      const uint32_t variant_id_slen = variant_id_end - variant_id_str;
      if (unlikely(variant_id_slen > kMaxIdSlen)) {
        putc_unlocked('\n', stdout);
        logerrputs("Error: Variant names are limited to " MAX_ID_SLEN_STR " characters.\n");
        goto OxGenToPgen_ret_MALFORMED_INPUT;
      }
      pvar_cswritep = memcpyax(pvar_cswritep, variant_id_str, variant_id_slen, '\t');

      // .gen specification does not define which column should be expected to
      // be the reference allele, and which the alternate.  plink 1.9 assumed
      // alt was usually first, but the reverse seems to be more common now.
      // So:
      //   If 'ref-first' or 'ref-last' was specified, we know what to do.
      //   If not, we treat the second allele as the provisional reference, for
      //     backward compatibility.
      const char* first_allele_str = FirstNonTspace(pos_end);
      if (unlikely(IsEolnKns(*first_allele_str))) {
        goto OxGenToPgen_ret_MISSING_TOKENS;
      }
      const char* first_allele_end = CurTokenEnd(first_allele_str);
      const char* second_allele_str = FirstNonTspace(first_allele_end);
      if (unlikely(IsEolnKns(*second_allele_str))) {
        goto OxGenToPgen_ret_MISSING_TOKENS;
      }
      const char* linebuf_iter = CurTokenEnd(second_allele_str);
      if (!prov_ref_allele_second) {
        pvar_cswritep = memcpyax(pvar_cswritep, first_allele_str, first_allele_end - first_allele_str, '\t');
        pvar_cswritep = memcpya(pvar_cswritep, second_allele_str, linebuf_iter - second_allele_str);
      } else {
        pvar_cswritep = memcpyax(pvar_cswritep, second_allele_str, linebuf_iter - second_allele_str, '\t');
        pvar_cswritep = memcpya(pvar_cswritep, first_allele_str, first_allele_end - first_allele_str);
      }
      AppendBinaryEoln(&pvar_cswritep);
      if (unlikely(Cswrite(&pvar_css, &pvar_cswritep))) {
        goto OxGenToPgen_ret_WRITE_FAIL;
      }

      if (!dosage_exists) {
        for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
          linebuf_iter = FirstNonTspace(linebuf_iter);
          const char cc = *linebuf_iter;
          if (unlikely(IsEolnKns(cc))) {
            goto OxGenToPgen_ret_MISSING_TOKENS;
          }
          char cc2 = linebuf_iter[1];
          if ((cc2 == ' ') || (cc2 == '\t')) {
            // fast handling of "1 0 0 ", "0 1 0 ", "0 0 1 ", "0 0 0 " cases
            cc2 = linebuf_iter[3];
            if (((cc2 == ' ') || (cc2 == '\t')) && (ctou32(linebuf_iter[5]) <= 32)) {
              const uint32_t uii = ctou32(cc) - 48;
              const uint32_t ujj = ctou32(linebuf_iter[2]) - 48;
              const uint32_t ukk = ctou32(linebuf_iter[4]) - 48;
              if (((uii | ujj | ukk) < 2) && (uii + ujj + ukk < 2)) {
                linebuf_iter = &(linebuf_iter[5]);
                continue;
              }
            }
          }
          double prob_0alt;
          const char* first_dosage_str_end = ScanadvDouble(linebuf_iter, &prob_0alt);
          if (!first_dosage_str_end) {
            // triple-NA, etc. ok; treat as missing value
            linebuf_iter = NextTokenMult(linebuf_iter, 2);
            if (unlikely(!linebuf_iter)) {
              goto OxGenToPgen_ret_MISSING_TOKENS;
            }
            linebuf_iter = CurTokenEnd(linebuf_iter);
            continue;
          }
          if (unlikely(ctou32(*first_dosage_str_end) > ' ')) {
            goto OxGenToPgen_ret_INVALID_DOSAGE;
          }
          linebuf_iter = FirstNonTspace(first_dosage_str_end);
          if (unlikely(IsEolnKns(*linebuf_iter))) {
            goto OxGenToPgen_ret_MISSING_TOKENS;
          }
          double prob_1alt;
          linebuf_iter = ScantokDouble(linebuf_iter, &prob_1alt);
          if (unlikely(!linebuf_iter)) {
            goto OxGenToPgen_ret_INVALID_DOSAGE;
          }
          linebuf_iter = FirstNonTspace(linebuf_iter);
          if (unlikely(IsEolnKns(*linebuf_iter))) {
            goto OxGenToPgen_ret_MISSING_TOKENS;
          }
          double prob_2alt;
          linebuf_iter = ScantokDouble(linebuf_iter, &prob_2alt);
          if (unlikely(!linebuf_iter)) {
            goto OxGenToPgen_ret_INVALID_DOSAGE;
          }
          // bugfix: forgot the "multiply by 32768" part of "multiply by 32768
          // and round" .gen -> .bgen conversion.
          prob_0alt *= 32768;
          prob_1alt *= 32768;
          prob_2alt *= 32768;

          // now treat this identically to bgen-1.1
          // Compare with 65535.4999999999 instead of 65535.5 since 0.5 +
          // <first floating point number below 65535.5> may evaluate to 65536.
          if (unlikely((prob_0alt < 0.0) || (prob_0alt >= 65535.4999999999) || (prob_1alt < 0.0) || (prob_1alt >= 65535.4999999999) || (prob_2alt < 0.0) || (prob_2alt >= 65535.4999999999))) {
            goto OxGenToPgen_ret_INVALID_DOSAGE;
          }
          const uint32_t dosage_int0 = S_CAST(int32_t, prob_0alt + 0.5);
          const uint32_t dosage_int1 = S_CAST(int32_t, prob_1alt + 0.5);
          const uint32_t dosage_int2 = S_CAST(int32_t, prob_2alt + 0.5);
          dosage_exists = Bgen11DosageImportCheck(dosage_int_sum_thresh, import_dosage_certainty_int, dosage_erase_halfdist, dosage_int0, dosage_int1, dosage_int2);
          if (dosage_exists) {
            break;
          }
        }
      }
      if (!(variant_ct % 1000)) {
        printf("\r--data/--gen: %uk variants scanned.", variant_ct / 1000);
        fflush(stdout);
      }
      line_iter = AdvPastDelim(K_CAST(char*, linebuf_iter), '\n');
    }
    if (unlikely(TextStreamErrcode2(&gen_txs, &reterr))) {
      goto OxGenToPgen_ret_TSTREAM_FAIL;
    }
    putc_unlocked('\r', stdout);
    if (unlikely(CswriteCloseNull(&pvar_css, pvar_cswritep))) {
      goto OxGenToPgen_ret_WRITE_FAIL;
    }
    if (unlikely(!variant_ct)) {
      if (!variant_skip_ct) {
        logerrputs("Error: No variants in .gen file.\n");
        goto OxGenToPgen_ret_INCONSISTENT_INPUT;
      }
      char* write_iter = strcpya_k(g_logbuf, "Error: All ");
      write_iter = wtoa(variant_skip_ct, write_iter);
      write_iter = strcpya_k(write_iter, " variant");
      if (variant_skip_ct != 1) {
        *write_iter++ = 's';
      }
      write_iter = strcpya_k(write_iter, " in .gen file excluded by ");
      AppendLoadFilterFlagnames(load_filter_log_import_flags, &write_iter);
      strcpy_k(write_iter, ".\n");
      goto OxGenToPgen_ret_INCONSISTENT_INPUT_WW;
    }
    {
      char* write_iter = strcpya_k(g_logbuf, "--data/--gen: ");
      write_iter = wtoa(variant_ct + variant_skip_ct, write_iter);
      write_iter = strcpya_k(write_iter, " variant");
      if (variant_ct + variant_skip_ct != 1) {
        *write_iter++ = 's';
      }
      write_iter = strcpya_k(write_iter, " scanned");
      if (!dosage_exists) {
        write_iter = strcpya_k(write_iter, " (all hardcalls)");
      }
      if (variant_skip_ct) {
        write_iter = strcpya_k(write_iter, "; ");
        write_iter = wtoa(variant_skip_ct, write_iter);
        write_iter = strcpya_k(write_iter, " excluded by ");
        AppendLoadFilterFlagnames(load_filter_log_import_flags, &write_iter);
        write_iter = strcpya_k(write_iter, ", ");
        write_iter = u32toa(variant_ct, write_iter);
        write_iter = strcpya_k(write_iter, " remaining");
      } else if (load_filter_log_import_flags) {
        write_iter = strcpya_k(write_iter, " (");
        AppendLoadFilterFlagnames(load_filter_log_import_flags, &write_iter);
        write_iter = strcpya_k(write_iter, " had no effect)");
      }
      strcpy_k(write_iter, ".\n");
      WordWrapB(0);
      logputsb();
    }

    // second pass
    BigstackReset(bigstack_mark2);
    if (TextIsMt(&gen_txs) && (decompress_thread_ct > 1)) {
      // Close and reopen, so that we can reduce decompress_thread_ct to 1.
      char* dst = R_CAST(char*, TextStreamMemStart(&gen_txs));
      if (unlikely(CleanupTextStream(&gen_txs, &reterr))) {
        logerrprintfww(kErrprintfFread, ".gen file", rstrerror(errno));
        goto OxGenToPgen_ret_1;
      }
      reterr = TextStreamOpenEx(genname, kMaxLongLine, max_line_blen, 1, nullptr, dst, &gen_txs);
      if (unlikely(reterr)) {
        goto OxGenToPgen_ret_TSTREAM_FAIL;
      }
    } else {
      reterr = TextRewind(&gen_txs);
      if (unlikely(reterr)) {
        goto OxGenToPgen_ret_TSTREAM_FAIL;
      }
    }

    snprintf(outname_end, kMaxOutfnameExtBlen, ".pgen");
    uintptr_t spgw_alloc_cacheline_ct;
    uint32_t max_vrec_len;
    reterr = SpgwInitPhase1(outname, nullptr, nullptr, variant_ct, sample_ct, 0, kPgenWriteBackwardSeek, dosage_exists? kfPgenGlobalDosagePresent : kfPgenGlobal0, (oxford_import_flags & kfOxfordImportRefUnknown)? 2 : 1, &spgw, &spgw_alloc_cacheline_ct, &max_vrec_len);
    if (unlikely(reterr)) {
      if (reterr == kPglRetOpenFail) {
        logerrprintfww(kErrprintfFopen, outname, strerror(errno));
      }
      goto OxGenToPgen_ret_1;
    }
    unsigned char* spgw_alloc;
    if (unlikely(bigstack_alloc_uc(spgw_alloc_cacheline_ct * kCacheline, &spgw_alloc))) {
      goto OxGenToPgen_ret_NOMEM;
    }
    SpgwInitPhase2(max_vrec_len, &spgw, spgw_alloc);

    const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
    const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
    uintptr_t* genovec;
    Halfword* dosage_present_hwarr;
    // if we weren't using bigstack_alloc, this would need to be sample_ctaw2
    if (unlikely(bigstack_alloc_w(sample_ctl2, &genovec) ||
                 bigstack_alloc_hw(2 * sample_ctl, &dosage_present_hwarr))) {
      goto OxGenToPgen_ret_NOMEM;
    }
    Dosage* dosage_main = nullptr;
    if (dosage_exists) {
      if (unlikely(bigstack_alloc_dosage(sample_ct, &dosage_main))) {
        goto OxGenToPgen_ret_NOMEM;
      }
    }
    if (hard_call_thresh == UINT32_MAX) {
      hard_call_thresh = kDosageMid / 10;
    }
    const uint32_t hard_call_halfdist = kDosage4th - hard_call_thresh;
    const uint32_t sample_ctl2_m1 = sample_ctl2 - 1;
    const uintptr_t line_ct = line_idx - 1;
    uint32_t vidx = 0;
    line_iter = TextLineEnd(&gen_txs);
    for (line_idx = 1; line_idx <= line_ct; ++line_idx) {
      reterr = TextGetUnsafe(&gen_txs, &line_iter);
      if (unlikely(reterr)) {
        goto OxGenToPgen_ret_TSTREAM_FAIL;
      }
      char* chr_code_str = line_iter;
      char* chr_code_end = CurTokenEnd(chr_code_str);
      if (variant_skip_ct) {
        *chr_code_end = '\0';
        const uint32_t chr_code = GetChrCode(chr_code_str, cip, chr_code_end - chr_code_str);
        if (!IsSet(cip->chr_mask, chr_code)) {
          line_iter = AdvPastDelim(chr_code_end, '\n');
          continue;
        }
      }
      const char* linebuf_iter = NextTokenMult(FirstNonTspace(&(chr_code_end[1])), 4 + is_v2);
      uint32_t inner_loop_last = kBitsPerWordD2 - 1;
      Dosage* dosage_main_iter = dosage_main;
      for (uint32_t widx = 0; ; ++widx) {
        if (widx >= sample_ctl2_m1) {
          if (widx > sample_ctl2_m1) {
            break;
          }
          inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
        }
        uintptr_t genovec_word = 0;
        uint32_t dosage_present_hw = 0;
        for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
          linebuf_iter = FirstNonTspace(linebuf_iter);
          const char cc = *linebuf_iter;
          if (unlikely(IsEolnKns(cc))) {
            goto OxGenToPgen_ret_MISSING_TOKENS;
          }
          char cc2 = linebuf_iter[1];
          if ((cc2 == ' ') || (cc2 == '\t')) {
            cc2 = linebuf_iter[3];
            if (((cc2 == ' ') || (cc2 == '\t')) && (ctou32(linebuf_iter[5]) <= 32)) {
              const uint32_t uii = ctou32(cc) - 48;
              const uint32_t ujj = ctou32(linebuf_iter[2]) - 48;
              const uint32_t ukk = ctou32(linebuf_iter[4]) - 48;
              const uint32_t all_or = uii | ujj | ukk;
              if (all_or < 2) {
                if (!all_or) {
                  genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
                  linebuf_iter = &(linebuf_iter[5]);
                  continue;
                }
                if (uii + ujj + ukk == 1) {
                  uintptr_t cur_geno = ukk * 2 + ujj;
                  genovec_word |= cur_geno << (2 * sample_idx_lowbits);
                  linebuf_iter = &(linebuf_iter[5]);
                  continue;
                }
              }
            }
          }
          double prob_0alt;
          const char* first_dosage_str_end = ScanadvDouble(linebuf_iter, &prob_0alt);
          if (!first_dosage_str_end) {
            // ignore next two tokens if first token in triplet is not numeric,
            // since we treat this as missing regardless
            linebuf_iter = NextTokenMult(linebuf_iter, 2);
            if (unlikely(!linebuf_iter)) {
              goto OxGenToPgen_ret_MISSING_TOKENS;
            }
            genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
            linebuf_iter = CurTokenEnd(linebuf_iter);
            continue;
          }
          if (unlikely(ctou32(*first_dosage_str_end) > ' ')) {
            goto OxGenToPgen_ret_INVALID_DOSAGE;
          }
          linebuf_iter = FirstNonTspace(first_dosage_str_end);
          if (unlikely(IsEolnKns(*linebuf_iter))) {
            goto OxGenToPgen_ret_MISSING_TOKENS;
          }
          double prob_1alt;
          linebuf_iter = ScantokDouble(linebuf_iter, &prob_1alt);
          if (unlikely(!linebuf_iter)) {
            goto OxGenToPgen_ret_INVALID_DOSAGE;
          }
          linebuf_iter = FirstNonTspace(linebuf_iter);
          if (unlikely(IsEolnKns(*linebuf_iter))) {
            goto OxGenToPgen_ret_MISSING_TOKENS;
          }
          double prob_2alt;
          linebuf_iter = ScantokDouble(linebuf_iter, &prob_2alt);
          if (unlikely(!linebuf_iter)) {
            goto OxGenToPgen_ret_INVALID_DOSAGE;
          }
          // bugfix
          prob_0alt *= 32768;
          prob_1alt *= 32768;
          prob_2alt *= 32768;

          if (unlikely((prob_0alt < 0.0) || (prob_0alt >= 65535.4999999999) || (prob_1alt < 0.0) || (prob_1alt >= 65535.4999999999) || (prob_2alt < 0.0) || (prob_2alt >= 65535.4999999999))) {
            goto OxGenToPgen_ret_INVALID_DOSAGE;
          }
          const uint32_t dosage_int0 = S_CAST(int32_t, prob_0alt + 0.5);
          const uint32_t dosage_int1 = S_CAST(int32_t, prob_1alt + 0.5);
          const uint32_t dosage_int2 = S_CAST(int32_t, prob_2alt + 0.5);
          Bgen11DosageImportUpdate(dosage_int_sum_thresh, import_dosage_certainty_int, hard_call_halfdist, dosage_erase_halfdist, sample_idx_lowbits, dosage_int0, dosage_int1, dosage_int2, &genovec_word, &dosage_present_hw, &dosage_main_iter);
        }
        genovec[widx] = genovec_word;
        dosage_present_hwarr[widx] = dosage_present_hw;
      }
      if (prov_ref_allele_second) {
        GenovecInvertUnsafe(sample_ct, genovec);
        ZeroTrailingNyps(sample_ct, genovec);
      }
      if (dosage_main_iter != dosage_main) {
        const uint32_t dosage_ct = dosage_main_iter - dosage_main;
        if (prov_ref_allele_second) {
          BiallelicDosage16Invert(dosage_ct, dosage_main);
        }
        uintptr_t* __attribute__((may_alias)) dosage_present_warr = R_CAST(uintptr_t*, dosage_present_hwarr);
        reterr = SpgwAppendBiallelicGenovecDosage16(genovec, dosage_present_warr, dosage_main, dosage_ct, &spgw);
        if (unlikely(reterr)) {
          goto OxGenToPgen_ret_1;
        }
      } else {
        if (unlikely(SpgwAppendBiallelicGenovec(genovec, &spgw))) {
          goto OxGenToPgen_ret_WRITE_FAIL;
        }
      }
      ++vidx;
      if (!(vidx % 1000)) {
        printf("\r--data/--gen: %uk variants converted.", vidx / 1000);
        fflush(stdout);
      }
      line_iter = AdvPastDelim(K_CAST(char*, linebuf_iter), '\n');
    }
    reterr = SpgwFinish(&spgw);
    if (unlikely(reterr)) {
      goto OxGenToPgen_ret_1;
    }
    putc_unlocked('\r', stdout);
    char* write_iter = strcpya_k(g_logbuf, "--data/--gen: ");
    const uint32_t outname_base_slen = outname_end - outname;
    write_iter = memcpya(write_iter, outname, outname_base_slen + 5);
    write_iter = strcpya_k(write_iter, " + ");
    write_iter = memcpya(write_iter, outname, outname_base_slen);
    write_iter = strcpya_k(write_iter, ".pvar");
    if (output_zst) {
      write_iter = strcpya_k(write_iter, ".zst");
    }
    snprintf(write_iter, kLogbufSize - 2 * kPglFnamesize - 64, " written.\n");
    WordWrapB(0);
    logputsb();
    if (print_splitpar_warning) {
      logerrputs("Warning: Human chrX pseudoautosomal variant(s) appear to be present in the\ninput .gen.  You probably want to include --split-par in your next command.\n");
    }
  }
  while (0) {
  OxGenToPgen_ret_TSTREAM_FAIL:
    TextStreamErrPrint(".gen file", &gen_txs);
    break;
  OxGenToPgen_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  OxGenToPgen_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  OxGenToPgen_ret_INVALID_DOSAGE:
    putc_unlocked('\n', stdout);
    logerrprintf("Error: Line %" PRIuPTR " of .gen file has an invalid dosage value.\n", line_idx);
    reterr = kPglRetMalformedInput;
    break;
  OxGenToPgen_ret_MISSING_TOKENS:
    putc_unlocked('\n', stdout);
    logerrprintf("Error: Line %" PRIuPTR " of .gen file has fewer tokens than expected.\n", line_idx);
    reterr = kPglRetMalformedInput;
    break;
  OxGenToPgen_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  OxGenToPgen_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  OxGenToPgen_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  OxGenToPgen_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  OxGenToPgen_ret_DEGENERATE_DATA:
    reterr = kPglRetDegenerateData;
    break;
  }
 OxGenToPgen_ret_1:
  CleanupSpgw(&spgw, &reterr);
  CleanupTextStream2(".gen file", &gen_txs, &reterr);
  CswriteCloseCond(&pvar_css, pvar_cswritep);
  BigstackReset(bigstack_mark);
  return reterr;
}

typedef struct BgenImportCommonStruct {
  uint32_t sample_ct;
  uint32_t dosage_erase_halfdist;
  uint32_t compression_mode;

  struct libdeflate_decompressor** libdeflate_decompressors;

  unsigned char** compressed_geno_starts[2];

  uint32_t cur_block_size;
} BgenImportCommon;

typedef struct Bgen11DosageScanCtxStruct {
  BgenImportCommon* common;

  uint16_t** bgen_geno_bufs;
  uint32_t import_dosage_certainty_int;

  uint32_t dosage_exists;
  PglErr reterr;  // only kPglRetMalformedInput possible for now
} Bgen11DosageScanCtx;

THREAD_FUNC_DECL Bgen11DosageScanThread(void* raw_arg) {
  // this bails as soon as a single non-hardcall is detected.  still
  // multithreaded due to relatively low speed of zlib_decompress() call,
  // practical value of handling the all-hardcall case efficiently, and reduced
  // code complexity (locally more complex, but globally cleaner due to overlap
  // with Bgen11GenoToPgenThread()).
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  Bgen11DosageScanCtx* ctx = S_CAST(Bgen11DosageScanCtx*, arg->sharedp->context);
  BgenImportCommon* bicp = ctx->common;

  const uint32_t sample_ct = bicp->sample_ct;
  uint16_t* bgen_geno_buf = ctx->bgen_geno_bufs[tidx];
  struct libdeflate_decompressor* decompressor = bicp->libdeflate_decompressors[tidx];
  const uint32_t calc_thread_ct = GetThreadCt(arg->sharedp);
  // hard_call_halfdist irrelevant here
  const uint32_t dosage_erase_halfdist = bicp->dosage_erase_halfdist;
  const uint32_t import_dosage_certainty_int = ctx->import_dosage_certainty_int;
  const uint32_t dosage_int_sum_thresh = 3 * (import_dosage_certainty_int - 1);
  const uint32_t compression_mode = bicp->compression_mode;
  // uint32_t vidx_base = 0;
  uint32_t parity = 0;
  do {
    const uintptr_t cur_block_size = bicp->cur_block_size;
    const uint32_t vidx_end = ((tidx + 1) * cur_block_size) / calc_thread_ct;
    uint32_t vidx = (tidx * cur_block_size) / calc_thread_ct;
    unsigned char* compressed_geno_iter = bicp->compressed_geno_starts[parity][vidx];
    uint16_t* bgen_probs = bgen_geno_buf;
    for (; vidx != vidx_end; ++vidx) {
      if (compression_mode) {
        uint32_t compressed_block_byte_ct;
        memcpy(&compressed_block_byte_ct, compressed_geno_iter, 4);
        compressed_geno_iter = &(compressed_geno_iter[4]);
        if (unlikely(libdeflate_zlib_decompress(decompressor, compressed_geno_iter, compressed_block_byte_ct, bgen_probs, 6 * sample_ct, nullptr) != LIBDEFLATE_SUCCESS)) {
          break;
        }
        compressed_geno_iter = &(compressed_geno_iter[compressed_block_byte_ct]);
      } else {
        bgen_probs = R_CAST(uint16_t*, compressed_geno_iter);
        compressed_geno_iter = &(compressed_geno_iter[6 * sample_ct]);
      }
      const uint16_t* bgen_probs_iter = bgen_probs;
      for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        const uint32_t dosage_int0 = *bgen_probs_iter++;
        const uint32_t dosage_int1 = *bgen_probs_iter++;
        const uint32_t dosage_int2 = *bgen_probs_iter++;
        const uint32_t dosage_int_sum = dosage_int0 + dosage_int1 + dosage_int2;
        if ((dosage_int_sum > dosage_int_sum_thresh) || (dosage_int0 >= import_dosage_certainty_int) || (dosage_int1 >= import_dosage_certainty_int) || (dosage_int2 >= import_dosage_certainty_int)) {
          // ties realistically happen, use banker's rounding
          // 1/65536 -> 0/32768
          // 3/65536 -> 2/32768
          // 5/65536 -> 2/32768
          const DosageProd write_dosage_int_numer = S_CAST(DosageProd, kDosageMid) * dosage_int1 + S_CAST(DosageProd, kDosageMax) * dosage_int2;
          uint32_t write_dosage_int;
          if (dosage_int_sum == kDosageMax) {
            // optimize common case
            write_dosage_int = ((write_dosage_int_numer + kDosageMid) / kDosageMax) - ((write_dosage_int_numer % (2 * S_CAST(DosageProd, kDosageMax))) == kDosageMid);
          } else {
            write_dosage_int = (write_dosage_int_numer + (dosage_int_sum / 2)) / dosage_int_sum;
            write_dosage_int -= (2 * (write_dosage_int_numer - write_dosage_int * dosage_int_sum) == dosage_int_sum) * (write_dosage_int % 2);
          }
          const uint32_t halfdist = BiallelicDosageHalfdist(write_dosage_int);
          if (halfdist < dosage_erase_halfdist) {
            goto Bgen11DosageScanThread_dosage_found;
          }
        }
      }
    }
    if (unlikely(vidx != vidx_end)) {
      // g_error_vidxs[tidx] = vidx + vidx_base;
      ctx->reterr = kPglRetMalformedInput;
    }
    while (0) {
    Bgen11DosageScanThread_dosage_found:
      ctx->dosage_exists = 1;
    }
    // vidx_base += cur_block_size;
    parity = 1 - parity;
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

typedef struct Bgen11GenoToPgenCtxStruct {
  BgenImportCommon* common;

  uint16_t** bgen_geno_bufs;
  uint32_t hard_call_halfdist;
  uint32_t import_dosage_certainty_int;
  uint32_t prov_ref_allele_second;

  uintptr_t* write_genovecs[2];
  uint32_t* write_dosage_cts[2];
  uintptr_t* write_dosage_presents[2];
  Dosage* write_dosage_mains[2];

  PglErr reterr;  // only kPglRetMalformedInput possible for now
} Bgen11GenoToPgenCtx;

static_assert(sizeof(Dosage) == 2, "Bgen11GenoToPgenThread() needs to be updated.");
THREAD_FUNC_DECL Bgen11GenoToPgenThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  Bgen11GenoToPgenCtx* ctx = S_CAST(Bgen11GenoToPgenCtx*, arg->sharedp->context);
  BgenImportCommon* bicp = ctx->common;

  const uintptr_t sample_ct = bicp->sample_ct;
  uint16_t* bgen_geno_buf = ctx->bgen_geno_bufs[tidx];
  struct libdeflate_decompressor* decompressor = bicp->libdeflate_decompressors[tidx];
  const uint32_t calc_thread_ct = GetThreadCt(arg->sharedp);
  const uint32_t hard_call_halfdist = ctx->hard_call_halfdist;
  const uint32_t dosage_erase_halfdist = bicp->dosage_erase_halfdist;
  const uint32_t import_dosage_certainty_int = ctx->import_dosage_certainty_int;
  const uint32_t dosage_int_sum_thresh = 3 * (import_dosage_certainty_int - 1);
  const uint32_t compression_mode = bicp->compression_mode;
  const uint32_t prov_ref_allele_second = ctx->prov_ref_allele_second;
  const uintptr_t sample_ctaw2 = NypCtToAlignedWordCt(sample_ct);
  const uint32_t sample_ctl2_m1 = (sample_ct - 1) / kBitsPerWordD2;
  const uintptr_t sample_ctaw = BitCtToAlignedWordCt(sample_ct);
  // uint32_t vidx_base = 0;
  uint32_t parity = 0;
  do {
    const uintptr_t cur_block_write_ct = bicp->cur_block_size;
    uint32_t vidx = (tidx * cur_block_write_ct) / calc_thread_ct;
    const uint32_t vidx_end = ((tidx + 1) * cur_block_write_ct) / calc_thread_ct;
    unsigned char* compressed_geno_iter = bicp->compressed_geno_starts[parity][vidx];
    uintptr_t* write_genovec_iter = &(ctx->write_genovecs[parity][vidx * sample_ctaw2]);
    uint32_t* write_dosage_ct_iter = &(ctx->write_dosage_cts[parity][vidx]);
    Halfword* write_dosage_present_iter = DowncastWToHW(&(ctx->write_dosage_presents[parity][vidx * sample_ctaw]));
    Dosage* write_dosage_main_iter = &(ctx->write_dosage_mains[parity][vidx * sample_ct]);
    uint16_t* bgen_probs = bgen_geno_buf;
    for (; vidx != vidx_end; ++vidx) {
      if (compression_mode) {
        uint32_t compressed_block_byte_ct;
        memcpy(&compressed_block_byte_ct, compressed_geno_iter, 4);
        compressed_geno_iter = &(compressed_geno_iter[4]);
        if (unlikely(libdeflate_zlib_decompress(decompressor, compressed_geno_iter, compressed_block_byte_ct, bgen_probs, 6 * sample_ct, nullptr) != LIBDEFLATE_SUCCESS)) {
          break;
        }
        compressed_geno_iter = &(compressed_geno_iter[compressed_block_byte_ct]);
      } else {
        bgen_probs = R_CAST(uint16_t*, compressed_geno_iter);
        compressed_geno_iter = &(compressed_geno_iter[6 * sample_ct]);
      }
      const uint16_t* bgen_probs_iter = bgen_probs;
      Dosage* cur_dosage_main_iter = write_dosage_main_iter;
      uint32_t inner_loop_last = kBitsPerWordD2 - 1;
      for (uint32_t widx = 0; ; ++widx) {
        if (widx >= sample_ctl2_m1) {
          if (widx > sample_ctl2_m1) {
            break;
          }
          inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
        }
        uintptr_t genovec_word = 0;
        uint32_t dosage_present_hw = 0;
        for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
          const uint32_t dosage_int0 = *bgen_probs_iter++;
          const uint32_t dosage_int1 = *bgen_probs_iter++;
          const uint32_t dosage_int2 = *bgen_probs_iter++;
          const uint32_t dosage_int_sum = dosage_int0 + dosage_int1 + dosage_int2;
          if (dosage_int_sum <= dosage_int_sum_thresh) {
            if ((dosage_int0 < import_dosage_certainty_int) && (dosage_int1 < import_dosage_certainty_int) && (dosage_int2 < import_dosage_certainty_int)) {
              genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
              continue;
            }
          }
          const DosageProd write_dosage_int_numer = S_CAST(DosageProd, kDosageMid) * dosage_int1 + S_CAST(DosageProd, kDosageMax) * dosage_int2;
          uint32_t write_dosage_int;
          if (dosage_int_sum == kDosageMax) {
            write_dosage_int = ((write_dosage_int_numer + kDosageMid) / kDosageMax) - ((write_dosage_int_numer % (2 * S_CAST(DosageProd, kDosageMax))) == kDosageMid);
          } else {
            write_dosage_int = (write_dosage_int_numer + (dosage_int_sum / 2)) / dosage_int_sum;
            write_dosage_int -= (2 * (write_dosage_int_numer - write_dosage_int * dosage_int_sum) == dosage_int_sum) * (write_dosage_int % 2);
          }
          const uint32_t halfdist = BiallelicDosageHalfdist(write_dosage_int);
          if (halfdist < hard_call_halfdist) {
            genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
          } else {
            genovec_word |= ((write_dosage_int + (kDosage4th * k1LU)) / kDosageMid) << (2 * sample_idx_lowbits);
            if (halfdist >= dosage_erase_halfdist) {
              continue;
            }
          }
          dosage_present_hw |= 1U << sample_idx_lowbits;
          *cur_dosage_main_iter++ = write_dosage_int;
        }
        write_genovec_iter[widx] = genovec_word;
        write_dosage_present_iter[widx] = dosage_present_hw;
      }
      const uint32_t dosage_ct = cur_dosage_main_iter - write_dosage_main_iter;
      if (prov_ref_allele_second) {
        GenovecInvertUnsafe(sample_ct, write_genovec_iter);
        ZeroTrailingNyps(sample_ct, write_genovec_iter);
        if (dosage_ct) {
          BiallelicDosage16Invert(dosage_ct, write_dosage_main_iter);
        }
      }
      *write_dosage_ct_iter++ = dosage_ct;
      write_genovec_iter = &(write_genovec_iter[sample_ctaw2]);
      write_dosage_present_iter = &(write_dosage_present_iter[2 * sample_ctaw]);
      write_dosage_main_iter = &(write_dosage_main_iter[sample_ct]);
    }
    if (unlikely(vidx != vidx_end)) {
      // g_error_vidxs[tidx] = vidx + vidx_base;
      ctx->reterr = kPglRetMalformedInput;
    }
    // vidx_base += cur_block_write_ct;
    parity = 1 - parity;
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

HEADER_INLINE uintptr_t Bgen13GetOneVal(const unsigned char* prob_start, uint64_t prob_offset, uint32_t bit_precision, uintptr_t numer_mask) {
  const uint64_t bit_offset = prob_offset * bit_precision;
  uint64_t relevant_bits;
  // This can read slightly past the end of the buffer.
  memcpy(&relevant_bits, &(prob_start[bit_offset / CHAR_BIT]), sizeof(int64_t));
  return (relevant_bits >> (bit_offset % CHAR_BIT)) & numer_mask;
}

HEADER_INLINE void Bgen13GetTwoVals(const unsigned char* prob_start, uint64_t prob_offset, uint32_t bit_precision, uintptr_t numer_mask, uintptr_t* first_val_ptr, uintptr_t* second_val_ptr) {
  const uint64_t bit_offset = prob_offset * bit_precision;
  uint64_t relevant_bits;
  // This can read slightly past the end of the buffer.
  // Note that with bit_precision=29 and variable ploidy,
  // (bit_offset % CHAR_BIT) == 7 is possible, so we may only get 57 bits when
  // we need 58; thus we don't support 29-31 bits for now.
  memcpy(&relevant_bits, &(prob_start[bit_offset / CHAR_BIT]), sizeof(int64_t));
  relevant_bits = relevant_bits >> (bit_offset % CHAR_BIT);
  *first_val_ptr = relevant_bits & numer_mask;
  *second_val_ptr = (relevant_bits >> bit_precision) & numer_mask;
}

// returns nonzero if phase or dosage found.
uint32_t Bgen13ScanBiallelicPhased(uintptr_t numer_a1, uintptr_t numer_a2, uintptr_t numer_mask, uint32_t numer_certainty_min, uint64_t magic_preadd, uint64_t magic_mult, uint32_t magic_postshift, uint32_t dosage_erase_halfdist) {
  // we're almost certainly exiting quickly, so hardcall bypass is unimportant
  // here.

  // --import-dosage-certainty: Original plan was to treat the haplotype
  // dosages as independent, but unfortunately that produces a gross
  // export-reimport discontinuity.  To avoid that discontinuity, we should
  // assume maximal het probability.
  const uintptr_t numer_sum = numer_a1 + numer_a2;
  if (numer_certainty_min) {
    const uint32_t dist_from_1 = abs_i32(numer_sum - numer_mask);
    if ((dist_from_1 < numer_certainty_min) && (numer_mask - dist_from_1 < numer_certainty_min)) {
      return 0;
    }
  }
  const uint32_t write_dosage_int = (magic_preadd + magic_mult * numer_sum) >> magic_postshift;
  if (numer_a1 != numer_a2) {
    // 'Found' if we save phase OR dosage info.
    // * Dosage is erased if phasedist1 and phasedist2 are both >=
    //   dphase_erase_halfdist.  But we still save phase if the two sides don't
    //   round to the same integer.
    // * Otherwise, we usually save dosage.  However, there's one exception:
    //   write_dphase_delta_val == 0 and the dosage_erase_halfdist inequality
    //   holds.  Then we just save an unphased hardcall.

    // +numer_mask to force this to be nonnegative
    const uint32_t write_dphase_delta_p1 = (magic_preadd + magic_mult * (numer_a1 + numer_mask - numer_a2)) >> magic_postshift;
    const int32_t write_dphase_delta_val = write_dphase_delta_p1 - kDosageMid;

    // now on 0..32768 scale
    const uint32_t dphase_side1 = write_dosage_int + write_dphase_delta_val;
    const uint32_t dphase_side2 = write_dosage_int - write_dphase_delta_val;

    const uint32_t dphase_halfdist1 = DphaseHalfdist(dphase_side1);
    const uint32_t dphase_halfdist2 = DphaseHalfdist(dphase_side2);
    const uint32_t dphase_erase_halfdist = dosage_erase_halfdist + kDosage4th;
    if ((dphase_halfdist1 >= dphase_erase_halfdist) && (dphase_halfdist2 >= dphase_erase_halfdist)) {
      return (((dphase_side1 + kDosageMid) ^ (dphase_side2 + kDosageMid)) & kDosageMax);
    }
    if (write_dphase_delta_val != 0) {
      return 1;
    }
  }
  const uint32_t halfdist = BiallelicDosageHalfdist(write_dosage_int);
  return (halfdist < dosage_erase_halfdist);
}

/*
// Reliably fast division by constants of the form 2^n - 1; see
//   http://ridiculousfish.com/blog/posts/labor-of-division-episode-iii.html
// The general case also requires a preshift parameter, but it's always zero
// for the odd .bgen denominators.
typedef struct BgenMagicNumStruct {
  uint64_t totq_magic;
  uint32_t totq_postshift;
  uint32_t totq_incr;
} BgenMagicNum;

static const BgenMagicNum kBgenMagicNums[25] = {
  {0, 0, 0},
  {1, 0, 0},
  {2863311531U, 33, 0},
  {1227133513U, 33, 1},
  {2290649225U, 35, 0},
  {1108378657U, 35, 1},
  {1090785345U, 36, 1},
  {270549121, 35, 1},
  {2155905153U, 39, 0},
  {134480385, 36, 1},
  {1074791425U, 40, 1},
  {4196353, 33, 1},
  {16781313, 36, 1},
  {67117057, 39, 1},
  {268451841, 42, 1},
  {1073774593U, 45, 1},
  {2147516417U, 47, 0},
  {131073, 34, 1},
  {262145, 36, 1},
  {524289, 38, 1},
  {1048577, 40, 1},
  {2097153, 42, 1},
  {4194305, 44, 1},
  {8388609, 46, 1},
  {16777217, 48, 1}
  // need to switch to a different algorithm past this point thanks to overflow
  // issues
};
*/

// bgen-1.3 diploid import requires the following operation (haploid is
// identical except b is always zero):
//   round((32768a + 16384b)/(2^{bit precision} - 1))
//   floor((32768a + 16384b)/(2^{bit_precision} - 1) + 0.5)
// = floor((32768a + 16384b + 2^{bit_precision - 1})
//     / (2^{bit_precision} - 1))
// = (totq_magic * (32768a + 16384b + 2^{bits-1} + totq_incr))
//     >> totq_postshift
//
// This works fine for bit_precision <= 16, anyway.  There are two issues which
// come up with higher precision:
// 1. The ridiculous_fish magic numbers assume a 32-bit dividend.  Our dividend
//    is guaranteed to be divisible by 2^14, but it can be as large as
//      (2^{bits} - 1) * 2^15 + 2^{bits-1}.
// 2. Relatedly, the current sequence of operations multiplies totq_magic by
//    (dividend + totq_incr) (where totq_incr is zero or one); this
//    intermediate result must not overflow a uint64_t.
// It turns out that neither of these problems are relevant for bit_precision
// in [17, 24], but there's a failure at bits=25.
//
// However, the following rearrangement works for bits in [2, 31]:
//   (preadd + (totq_magic * (2a + b))) >> (totq_postshift - 14)
//   where preadd := (totq_magic * (2^{bits-1} + totq_incr)) >> 14
// and setting {preadd = 0, mult = 2^14, postshift = 0} works well enough for
// bits=1.
// It fails for bits=32, 2a+b=2147614719.
//
// The following code was used to generate and validate the body of the table
// below:
// for (uint32_t uii = 2; uii <= 31; ++uii) {
//   const uint32_t divisor = (1LLU << uii) - 1U;
//   uint64_t mult;
//   uint32_t preshift;
//   uint32_t postshift;
//   uint32_t incr;
//   DivisionMagicNums(divisor, &mult, &preshift, &postshift, &incr);
//   const uint64_t half = 1U << (uii - 1);
//   const uint64_t preadd = (mult * (half + incr)) >> 14;
//   const uint32_t postshift_m14 = postshift - 14;
//   printf("{%llu%s, %llu%s, %u},\n", preadd, (preadd > 0x7fffffff)? "LLU" : "", mult, (mult > 0x7fffffff)? "U" : "", postshift_m14);
//   for (uint64_t ullii = 0; ullii <= 2 * divisor; ++ullii) {
//     uint64_t numer = ullii * 16384;
//     uint64_t true_result = (numer + half) / divisor;
//     uint64_t my_result = (preadd + mult * ullii) >> postshift_m14;
//     if (true_result != my_result) {
//       printf("failure: bits=%u, 2a+b=%llu\n", uii, ullii);
//       exit(1);
//     }
//   }
// }

static_assert(kDosageMid == 16384, "bgen-1.3 import magic numbers must be changed.");
typedef struct BgenMagicNumStruct {
  uint64_t preadd;
  uint32_t mult;  // must copy this into a uint64_t
  uint32_t postshift;
} BgenMagicNum;

// We throw a not-yet-supported error on bits>28 for now.
CONSTI32(kMaxBgenImportBits, 28);

static const BgenMagicNum kBgenMagicNums[kMaxBgenImportBits + 1] = {
  {0, 0, 0},
  {0, 16384, 0},
  {349525, 2863311531U, 19},
  {374491, 1227133513, 19},
  {1118481, 2290649225U, 21},
  {1150051, 1108378657, 21},
  {2197016, 1090785345, 22},
  {1073345, 270549121, 21},
  {16843009, 2155905153U, 25},
  {2109464, 134480385, 22},
  {33652832, 1074791425, 26},
  {262528, 4196353, 19},
  {2098688, 16781313, 22},
  {16783360, 67117057, 25},
  {134242305, 268451841, 28},
  {1073840131, 1073774593, 31},
  {4295032834LLU, 2147516417U, 33},
  {524300, 131073, 20},
  {2097176, 262145, 22},
  {8388656, 524289, 24},
  {33554528, 1048577, 26},
  {134217920, 2097153, 28},
  {536871296, 4194305, 30},
  {2147484416LLU, 8388609, 32},
  {8589936128LLU, 16777217, 34},
  {34359741440LLU, 33554433, 36},
  {137438959616LLU, 67108865, 38},
  {549755826176LLU, 134217729, 40},
  {2199023280128LLU, 268435457, 42}
  // ,{8796093071360LLU, 536870913, 44}
  // ,{35184372187136LLU, 1073741825, 46}
  // ,{140737488551936LLU, 2147483649U, 48}
};

typedef struct Bgen13DosageOrPhaseScanCtxStruct {
    BgenImportCommon* common;
  // For each bit precision level, how large must
  //   max(numerators, 2^{bit_precision} - 1 - [sum of numerators])
  // be to avoid throwing out the genotype?
  // Idempotence (consistently get same data back after export and reimport) is
  // not possible for --import-dosage-certainty, so we may as well operate on
  // the pre-conversion numerators.
  uint32_t* bgen_import_dosage_certainty_thresholds;

  unsigned char** thread_wkspaces;

  uint32_t* thread_bidxs[2];
  uint16_t* bgen_allele_cts[2];
  uint32_t* uncompressed_genodata_byte_cts[2];

  const char** err_extra;

  // high 32 bits = bidx, earlier one takes precedence
  // next 16 bits = tidx if relevant (to look up err_extra)
  // next 8 bits = uint8_t(BgenImportErrSubtype)
  // low 8 bits = uint8_t(PglErr)
  uint64_t err_info;

  uint32_t error_on_polyploid;
  uint32_t vidx_start;

  uint32_t dosage_exists;
} Bgen13DosageOrPhaseScanCtx;

ENUM_U31_DEF_START()
  kBgenImportErrSubtype0,
  kBgenImportErrSubtypeLibdeflateZdecompress,
  kBgenImportErrSubtypeZstdDecompress,
  kBgenImportErrSubtypeSampleCtMismatch,
  kBgenImportErrSubtypePloidyOutOfRange,
  kBgenImportErrSubtypePolyploid,
  kBgenImportErrSubtypeIsPhasedOutOfRange,
  kBgenImportErrSubtypeBitPrecisionOutOfRange,
  kBgenImportErrSubtypeBitPrecisionUnsupported,
  kBgenImportErrSubtypeAlleleCtUnsupported,
  kBgenImportErrSubtypeUncompressedByteCtMismatch,
  kBgenImportErrSubtypeNumeratorOverflow,
  kBgenImportErrSubtypeProbOffsetMismatch,
  kBgenImportErrSubtypeInvalidMissingPloidy,
ENUM_U31_DEF_END(BgenImportErrSubtype);

void PrintBgenImportErr(const char** err_extra, uint64_t err_info, uint32_t vidx_base) {
  const BgenImportErrSubtype err_subtype = S_CAST(BgenImportErrSubtype, (err_info >> 8) & 255);
  const uint32_t vidx = (err_info >> 32) + vidx_base;
  logputs("\n");
  switch (err_subtype) {
  case kBgenImportErrSubtypeLibdeflateZdecompress:
    snprintf(g_logbuf, kLogbufSize, "Error: .bgen libdeflate_zlib_decompress() failure, vidx=%u.\n", vidx);
    break;
  case kBgenImportErrSubtypeZstdDecompress:
    {
      const uint32_t tidx = (err_info >> 16) & 65535;
      snprintf(g_logbuf, kLogbufSize, "Error: .bgen ZSTD_decompress() failure, vidx=%u: %s\n", vidx, err_extra[tidx]);
    }
    break;
  case kBgenImportErrSubtypeSampleCtMismatch:
    snprintf(g_logbuf, kLogbufSize, "Error: .bgen sample_ct mismatch, vidx=%u.\n", vidx);
    break;
  case kBgenImportErrSubtypePloidyOutOfRange:
    snprintf(g_logbuf, kLogbufSize, "Error: .bgen ploidy out of range, vidx=%u.\n", vidx);
    break;
  case kBgenImportErrSubtypePolyploid:
    snprintf(g_logbuf, kLogbufSize, "Error: Polyploid genotype in .bgen file, vidx=%u.\n", vidx);
    break;
  case kBgenImportErrSubtypeIsPhasedOutOfRange:
    snprintf(g_logbuf, kLogbufSize, "Error: .bgen is_phased out of range, vidx=%u.\n", vidx);
    break;
  case kBgenImportErrSubtypeBitPrecisionOutOfRange:
    snprintf(g_logbuf, kLogbufSize, "Error: .bgen bit_precision out of range, vidx=%u.\n", vidx);
    break;
  case kBgenImportErrSubtypeBitPrecisionUnsupported:
    snprintf(g_logbuf, kLogbufSize, "Error: .bgen bit_precision unsupported, vidx=%u.\n", vidx);
    break;
  case kBgenImportErrSubtypeAlleleCtUnsupported:
    snprintf(g_logbuf, kLogbufSize, "Error: .bgen allele_ct unsupported, vidx=%u.\n", vidx);
    break;
  case kBgenImportErrSubtypeUncompressedByteCtMismatch:
    snprintf(g_logbuf, kLogbufSize, "Error: .bgen uncompressed_byte_ct mismatch, vidx=%u.\n", vidx);
    break;
  case kBgenImportErrSubtypeNumeratorOverflow:
    snprintf(g_logbuf, kLogbufSize, "Error: .bgen numerator overflow, vidx=%u.\n", vidx);
    break;
  case kBgenImportErrSubtypeProbOffsetMismatch:
    snprintf(g_logbuf, kLogbufSize, "Error: .bgen prob_offset mismatch, vidx=%u.\n", vidx);
    break;
  case kBgenImportErrSubtypeInvalidMissingPloidy:
    snprintf(g_logbuf, kLogbufSize, "Error: Invalid .bgen missing_and_ploidy byte, vidx=%u.\n", vidx);
    break;
  default:
    ;
    const PglErr reterr = S_CAST(PglErr, err_info & 255);
    if (reterr == kPglRetNomem) {
      return;
    }
    assert(0);
    snprintf(g_logbuf, kLogbufSize, "Error: .bgen import internal error, vidx=%u.\n", vidx);
  }
  WordWrapB(0);
  logerrputsb();
}

static_assert(sizeof(Dosage) == 2, "Bgen13DosageOrPhaseScanThread() needs to be updated.");
THREAD_FUNC_DECL Bgen13DosageOrPhaseScanThread(void* raw_arg) {
  // This bails as soon as a single phased or dosage call is detected; that's
  // enough to determine the header format.
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  Bgen13DosageOrPhaseScanCtx* ctx = S_CAST(Bgen13DosageOrPhaseScanCtx*, arg->sharedp->context);
  BgenImportCommon* bicp = ctx->common;

  const uint32_t sample_ct = bicp->sample_ct;
  const uint32_t dosage_erase_halfdist = bicp->dosage_erase_halfdist;
  const uint32_t* bgen_import_dosage_certainty_thresholds = ctx->bgen_import_dosage_certainty_thresholds;
  const uint32_t compression_mode = bicp->compression_mode;
  const unsigned char* cur_uncompressed_geno = nullptr;
  struct libdeflate_decompressor* decompressor = bicp->libdeflate_decompressors[tidx];
  if (compression_mode) {
    cur_uncompressed_geno = ctx->thread_wkspaces[tidx];
  }
  uint32_t parity = 0;
  uint64_t new_err_info = 0;
  do {
    {
      const uint32_t bidx_end = ctx->thread_bidxs[parity][tidx + 1];
      uint32_t bidx = ctx->thread_bidxs[parity][tidx];
      unsigned char** compressed_geno_starts = bicp->compressed_geno_starts[parity];
      const uint16_t* bgen_allele_cts = ctx->bgen_allele_cts[parity];
      const uint32_t* uncompressed_genodata_byte_cts = ctx->uncompressed_genodata_byte_cts[parity];
      for (; bidx != bidx_end; ++bidx) {
        const unsigned char* compressed_geno_start = compressed_geno_starts[bidx];
        const unsigned char* compressed_geno_end = compressed_geno_starts[bidx + 1];
        uint32_t compressed_byte_ct = compressed_geno_end - compressed_geno_start;
        uint32_t uncompressed_byte_ct;
        if (compression_mode) {
          uncompressed_byte_ct = uncompressed_genodata_byte_cts[bidx];
          if (compression_mode == 1) {
            if (unlikely(libdeflate_zlib_decompress(decompressor, compressed_geno_start, compressed_byte_ct, K_CAST(unsigned char*, cur_uncompressed_geno), uncompressed_byte_ct, nullptr) != LIBDEFLATE_SUCCESS)) {
              new_err_info = (S_CAST(uint64_t, bidx) << 32) | (S_CAST(uint32_t, kBgenImportErrSubtypeLibdeflateZdecompress) << 8) | S_CAST(uint32_t, kPglRetMalformedInput);
              goto Bgen13DosageOrPhaseScanThread_err;
            }
          } else {
            const uintptr_t extracted_byte_ct = ZSTD_decompress(K_CAST(unsigned char*, cur_uncompressed_geno), uncompressed_byte_ct, compressed_geno_start, compressed_byte_ct);
            if (unlikely(extracted_byte_ct != uncompressed_byte_ct)) {
              assert(ZSTD_isError(extracted_byte_ct));
              ctx->err_extra[tidx] = ZSTD_getErrorName(extracted_byte_ct);
              new_err_info = (S_CAST(uint64_t, bidx) << 32) | (tidx << 16) | (S_CAST(uint32_t, kBgenImportErrSubtypeZstdDecompress) << 8) | S_CAST(uint32_t, kPglRetMalformedInput);
              goto Bgen13DosageOrPhaseScanThread_err;
            }
          }
        } else {
          cur_uncompressed_geno = compressed_geno_start;
          uncompressed_byte_ct = compressed_byte_ct;
        }
        // 4 bytes: sample_ct
        // 2 bytes: # of alleles, must match bgen_allele_cts[bidx]
        // 1 byte: min ploidy
        // 1 byte: max ploidy
        // sample_ct bytes: low 6 bits = ploidy, top bit = missingness
        // 1 byte: 1 if phased, 0 if not
        // 1 byte: # of bits of probability precision
        uint32_t stored_sample_ct;
        memcpy(&stored_sample_ct, cur_uncompressed_geno, sizeof(int32_t));
        if (unlikely((uncompressed_byte_ct < 10 + sample_ct) || (sample_ct != stored_sample_ct))) {
          new_err_info = (S_CAST(uint64_t, bidx) << 32) | (S_CAST(uint32_t, kBgenImportErrSubtypeSampleCtMismatch) << 8) | S_CAST(uint32_t, kPglRetMalformedInput);
          goto Bgen13DosageOrPhaseScanThread_err;
        }
        const uint32_t cur_allele_ct = bgen_allele_cts[bidx];
        /*
        uint16_t stored_allele_ct;
        memcpy(&stored_allele_ct, &(cur_uncompressed_geno[4]), sizeof(int16_t));
        if (unlikely(stored_allele_ct != cur_allele_ct)) {
          goto Bgen13DosageOrPhaseScanThread_malformed;
        }
        */
        const uint32_t min_ploidy = cur_uncompressed_geno[6];
        const uint32_t max_ploidy = cur_uncompressed_geno[7];
        if (unlikely((min_ploidy > max_ploidy) || (max_ploidy > 63))) {
          new_err_info = (S_CAST(uint64_t, bidx) << 32) | (S_CAST(uint32_t, kBgenImportErrSubtypePloidyOutOfRange) << 8) | S_CAST(uint32_t, kPglRetMalformedInput);
          goto Bgen13DosageOrPhaseScanThread_err;
        }
        if (unlikely(max_ploidy > 2)) {
          if (ctx->error_on_polyploid) {
            // you can't fire me, I quit!
            new_err_info = (S_CAST(uint64_t, bidx) << 32) | (S_CAST(uint32_t, kBgenImportErrSubtypePolyploid) << 8) | S_CAST(uint32_t, kPglRetInconsistentInput);
            goto Bgen13DosageOrPhaseScanThread_err;
          }
          new_err_info = (S_CAST(uint64_t, bidx) << 32) | (S_CAST(uint32_t, kBgenImportErrSubtypePolyploid) << 8) | S_CAST(uint32_t, kPglRetNotYetSupported);
          goto Bgen13DosageOrPhaseScanThread_err;
        }
        const unsigned char* missing_and_ploidy_info = &(cur_uncompressed_geno[8]);
        const unsigned char* probs_start = &(cur_uncompressed_geno[10 + sample_ct]);
        const uint32_t is_phased = probs_start[-2];
        if (unlikely(is_phased > 1)) {
          new_err_info = (S_CAST(uint64_t, bidx) << 32) | (S_CAST(uint32_t, kBgenImportErrSubtypeIsPhasedOutOfRange) << 8) | S_CAST(uint32_t, kPglRetMalformedInput);
          goto Bgen13DosageOrPhaseScanThread_err;
        }
        const uint32_t bit_precision = probs_start[-1];
        if (unlikely((!bit_precision) || (bit_precision > 32))) {
          new_err_info = (S_CAST(uint64_t, bidx) << 32) | (S_CAST(uint32_t, kBgenImportErrSubtypeBitPrecisionOutOfRange) << 8) | S_CAST(uint32_t, kPglRetMalformedInput);
          goto Bgen13DosageOrPhaseScanThread_err;
        }
        if (unlikely(bit_precision > kMaxBgenImportBits)) {
          new_err_info = (S_CAST(uint64_t, bidx) << 32) | (S_CAST(uint32_t, kBgenImportErrSubtypeBitPrecisionUnsupported) << 8) | S_CAST(uint32_t, kPglRetNotYetSupported);
          goto Bgen13DosageOrPhaseScanThread_err;
        }

        if (cur_allele_ct != 2) {
          // shouldn't be possible to get here for now
          assert(0);
          new_err_info = (S_CAST(uint64_t, bidx) << 32) | (S_CAST(uint32_t, kBgenImportErrSubtypeAlleleCtUnsupported) << 8) | S_CAST(uint32_t, kPglRetNotYetSupported);
          goto Bgen13DosageOrPhaseScanThread_err;
        }
        const uint64_t magic_preadd = kBgenMagicNums[bit_precision].preadd;
        const uint64_t magic_mult = kBgenMagicNums[bit_precision].mult;
        const uint32_t magic_postshift = kBgenMagicNums[bit_precision].postshift;

        // also equal to denominator
        const uintptr_t numer_mask = (1U << bit_precision) - 1;

        uint32_t numer_certainty_min = 0;
        if (bgen_import_dosage_certainty_thresholds) {
          numer_certainty_min = bgen_import_dosage_certainty_thresholds[bit_precision];
        }

        if (min_ploidy == max_ploidy) {
          // faster handling of common cases (no need to keep checking if
          // we've read past the end)
          if (unlikely(uncompressed_byte_ct != 10 + sample_ct + DivUp(S_CAST(uint64_t, bit_precision) * max_ploidy * sample_ct, CHAR_BIT))) {
            new_err_info = (S_CAST(uint64_t, bidx) << 32) | (S_CAST(uint32_t, kBgenImportErrSubtypeUncompressedByteCtMismatch) << 8) | S_CAST(uint32_t, kPglRetMalformedInput);
            goto Bgen13DosageOrPhaseScanThread_err;
          }
          if (max_ploidy < 2) {
            if (!max_ploidy) {
              // don't need to do anything in all-ploidy-0 case
              continue;
            }
            // biallelic, haploid
            // make sample_idx uint64_t so that when we multiply it by
            // bit_precision to get bit_offset, we don't have to worry about
            // overflow
            for (uintptr_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
              const uint32_t missing_and_ploidy = missing_and_ploidy_info[sample_idx];
              // treat anything else as missing
              if (missing_and_ploidy != 1) {
                continue;
              }
              const uintptr_t numer_a = Bgen13GetOneVal(probs_start, sample_idx, bit_precision, numer_mask);
              if ((numer_a < numer_certainty_min) && (numer_mask - numer_certainty_min < numer_a)) {
                continue;
              }
              const uint32_t write_dosage_int = (magic_preadd + magic_mult * 2 * numer_a) >> magic_postshift;
              const uint32_t halfdist = BiallelicDosageHalfdist(write_dosage_int);
              if (halfdist < dosage_erase_halfdist) {
                goto Bgen13DosageOrPhaseScanThread_found;
              }
            }
            continue;
          }
          // could combine the diploid phased and unphased cases into one loop,
          // but much of the similarity is superficial: ((2 + 1) choose 2) - 1
          // just happens to be equal to (2 - 1) * 2
          if (!is_phased) {
            // It's likely only one thread can execute this due to memory
            // limitations (since we don't know whether multiallelic variants
            // are present in advance, we try to use 4 GiB buffers), so this is
            // frequently more of a bottleneck than the full-decode loop when
            // no phase or dosage info is in the file.
            uintptr_t sample_idx = 0;
#ifdef __LP64__
            // Fast paths for common cases.
            if (bit_precision == 8) {
              // 1x2 bytes per entry
              const uint32_t full_vec_ct = sample_ct / kInt16PerVec;
              // If all bytes are equal to 0 or numer_mask, all calls in this
              // block must be hardcall/missing.
              const VecUc vec0 = vecuc_setzero();
              const VecUc vecmax = vecuc_set1(numer_mask);
              uint32_t vec_idx = 0;
              for (; vec_idx != full_vec_ct; ++vec_idx) {
                const VecUc cur_vec = vecuc_loadu(&(probs_start[vec_idx * kBytesPerVec]));
                const VecUc safe_bytes = (cur_vec == vec0) | (cur_vec == vecmax);
                if (!vec0255_is_all_set(safe_bytes)) {
                  break;
                }
              }
              sample_idx = vec_idx * kInt16PerVec;
            } else if (bit_precision == 16) {
              // 2x2 bytes per entry
              const uint32_t full_vec_ct = sample_ct / (kBytesPerVec / 4);
              // If all uint16s are equal to 0 or numer_mask, all calls in this
              // block must be hardcall/missing.
              const VecU16 vec0 = vecu16_setzero();
              const VecU16 vecmax = vecu16_set1(numer_mask);
              uint32_t vec_idx = 0;
              for (; vec_idx != full_vec_ct; ++vec_idx) {
                const VecU16 cur_vec = vecu16_loadu(&(probs_start[vec_idx * kBytesPerVec]));
                const VecU16 safe_u16s = (cur_vec == vec0) | (cur_vec == vecmax);
                if (!vec0255u16_is_all_set(safe_u16s)) {
                  break;
                }
              }
              sample_idx = vec_idx * (kBytesPerVec / 4);
            }
#endif
            for (; sample_idx != sample_ct; ++sample_idx) {
              const uint32_t missing_and_ploidy = missing_and_ploidy_info[sample_idx];
              // treat anything else as missing
              if (missing_and_ploidy != 2) {
                continue;
              }
              uintptr_t numer_aa;
              uintptr_t numer_ab;
              Bgen13GetTwoVals(probs_start, sample_idx * 2, bit_precision, numer_mask, &numer_aa, &numer_ab);
              // common trivial cases
              if (!numer_aa) {
                if ((!numer_ab) || (numer_ab == numer_mask)) {
                  continue;
                }
              } else if ((numer_aa == numer_mask) && (!numer_ab)) {
                continue;
              }
              if ((numer_aa < numer_certainty_min) && (numer_ab < numer_certainty_min) && (numer_mask - numer_certainty_min < numer_aa + numer_ab)) {
                // treat as missing
                continue;
              }
              const uint32_t write_dosage_int = (magic_preadd + magic_mult * (2 * numer_aa + numer_ab)) >> magic_postshift;
              const uint32_t halfdist = BiallelicDosageHalfdist(write_dosage_int);
              if (halfdist < dosage_erase_halfdist) {
                goto Bgen13DosageOrPhaseScanThread_found;
              }
            }
            continue;
          }
          // biallelic, phased, diploid
          for (uintptr_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
            const uint32_t missing_and_ploidy = missing_and_ploidy_info[sample_idx];
            if (missing_and_ploidy != 2) {
              continue;
            }
            uintptr_t numer_a1;
            uintptr_t numer_a2;
            Bgen13GetTwoVals(probs_start, sample_idx * 2, bit_precision, numer_mask, &numer_a1, &numer_a2);
            if (Bgen13ScanBiallelicPhased(numer_a1, numer_a2, numer_mask, numer_certainty_min, magic_preadd, magic_mult, magic_postshift, dosage_erase_halfdist)) {
              goto Bgen13DosageOrPhaseScanThread_found;
            }
          }
          continue;
        }
        // biallelic, variable ploidy
        const uint64_t remaining_bit_ct = 8LLU * (uncompressed_byte_ct - S_CAST(uintptr_t, probs_start - cur_uncompressed_geno));
        const uint64_t prob_offset_end = remaining_bit_ct / bit_precision;
        uintptr_t prob_offset = 0;
        if (!is_phased) {
          for (uintptr_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
            uint32_t missing_and_ploidy = missing_and_ploidy_info[sample_idx];
            if (missing_and_ploidy == 2) {
              if (unlikely(prob_offset + 2 > prob_offset_end)) {
                new_err_info = (S_CAST(uint64_t, bidx) << 32) | (S_CAST(uint32_t, kBgenImportErrSubtypeProbOffsetMismatch) << 8) | S_CAST(uint32_t, kPglRetMalformedInput);
                goto Bgen13DosageOrPhaseScanThread_err;
              }
              uintptr_t numer_aa;
              uintptr_t numer_ab;
              Bgen13GetTwoVals(probs_start, prob_offset, bit_precision, numer_mask, &numer_aa, &numer_ab);
              prob_offset += 2;
              if ((numer_aa < numer_certainty_min) && (numer_ab < numer_certainty_min) && (numer_mask - numer_certainty_min < numer_aa + numer_ab)) {
                // treat as missing
                continue;
              }
              const uint32_t write_dosage_int = (magic_preadd + magic_mult * (2 * numer_aa + numer_ab)) >> magic_postshift;
              const uint32_t halfdist = BiallelicDosageHalfdist(write_dosage_int);
              if (halfdist < dosage_erase_halfdist) {
                goto Bgen13DosageOrPhaseScanThread_found;
              }
            } else {
              if (missing_and_ploidy == 1) {
                // this is slightly liberal (we don't care to validate in
                // missing case); implementation would be different in a bgen
                // validator
                // note that this can't be moved before the if(), since
                // missing-ploidy could be zero
                if (unlikely(prob_offset >= prob_offset_end)) {
                  new_err_info = (S_CAST(uint64_t, bidx) << 32) | (S_CAST(uint32_t, kBgenImportErrSubtypeProbOffsetMismatch) << 8) | S_CAST(uint32_t, kPglRetMalformedInput);
                  goto Bgen13DosageOrPhaseScanThread_err;
                }
                const uintptr_t numer_a = Bgen13GetOneVal(probs_start, prob_offset, bit_precision, numer_mask);
                ++prob_offset;
                if ((numer_a < numer_certainty_min) && (numer_mask - numer_certainty_min < numer_a)) {
                  continue;
                }
                const uint32_t write_dosage_int = (magic_preadd + magic_mult * numer_a * 2) >> magic_postshift;
                const uint32_t halfdist = BiallelicDosageHalfdist(write_dosage_int);
                if (halfdist < dosage_erase_halfdist) {
                  goto Bgen13DosageOrPhaseScanThread_found;
                }
              } else {
                // treat as missing
                missing_and_ploidy &= 127;
                if (unlikely(missing_and_ploidy > 2)) {
                  new_err_info = (S_CAST(uint64_t, bidx) << 32) | (S_CAST(uint32_t, kBgenImportErrSubtypeInvalidMissingPloidy) << 8) | S_CAST(uint32_t, kPglRetMalformedInput);
                  goto Bgen13DosageOrPhaseScanThread_err;
                }
                prob_offset += missing_and_ploidy;
              }
            }
          }
          continue;
        }
        // biallelic, variable ploidy, phased
        for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
          uint32_t missing_and_ploidy = missing_and_ploidy_info[sample_idx];
          if (missing_and_ploidy == 2) {
            if (unlikely(prob_offset + 2 > prob_offset_end)) {
              new_err_info = (S_CAST(uint64_t, bidx) << 32) | (S_CAST(uint32_t, kBgenImportErrSubtypeProbOffsetMismatch) << 8) | S_CAST(uint32_t, kPglRetMalformedInput);
              goto Bgen13DosageOrPhaseScanThread_err;
            }
            uintptr_t numer_a1;
            uintptr_t numer_a2;
            Bgen13GetTwoVals(probs_start, prob_offset, bit_precision, numer_mask, &numer_a1, &numer_a2);
            prob_offset += 2;
            if (Bgen13ScanBiallelicPhased(numer_a1, numer_a2, numer_mask, numer_certainty_min, magic_preadd, magic_mult, magic_postshift, dosage_erase_halfdist)) {
              goto Bgen13DosageOrPhaseScanThread_found;
            }
          } else {
            if (missing_and_ploidy == 1) {
              if (unlikely(prob_offset >= prob_offset_end)) {
                new_err_info = (S_CAST(uint64_t, bidx) << 32) | (S_CAST(uint32_t, kBgenImportErrSubtypeProbOffsetMismatch) << 8) | S_CAST(uint32_t, kPglRetMalformedInput);
                goto Bgen13DosageOrPhaseScanThread_err;
              }
              const uintptr_t numer_a = Bgen13GetOneVal(probs_start, prob_offset, bit_precision, numer_mask);
              ++prob_offset;
              if ((numer_a < numer_certainty_min) && (numer_mask - numer_certainty_min < numer_a)) {
                continue;
              }
              const uint32_t write_dosage_int = (magic_preadd + magic_mult * numer_a * 2) >> magic_postshift;
              const uint32_t halfdist = BiallelicDosageHalfdist(write_dosage_int);
              if (halfdist < dosage_erase_halfdist) {
                goto Bgen13DosageOrPhaseScanThread_found;
              }
            } else {
              // treat as missing
              missing_and_ploidy &= 127;
              if (unlikely(missing_and_ploidy > 2)) {
                new_err_info = (S_CAST(uint64_t, bidx) << 32) | (S_CAST(uint32_t, kBgenImportErrSubtypeInvalidMissingPloidy) << 8) | S_CAST(uint32_t, kPglRetMalformedInput);
                goto Bgen13DosageOrPhaseScanThread_err;
              }
              prob_offset += missing_and_ploidy;
            }
          }
        }
      }
    }
    while (0) {
    Bgen13DosageOrPhaseScanThread_found:
      ctx->dosage_exists = 1;
      break;
    }
    parity = 1 - parity;
  } while (!THREAD_BLOCK_FINISH(arg));
  while (0) {
  Bgen13DosageOrPhaseScanThread_err:
    UpdateU64IfSmaller(new_err_info, &ctx->err_info);
    THREAD_BLOCK_FINISH(arg);
    break;
  }
  THREAD_RETURN;
}


// two trivial cases manually checked before this function is called, compiler
// is currently too dumm to handle them appropriately quickly if they're the
// first lines of this function, even when it's inline
// returns 1 if no dosage
// todo: more tests re: whether this should be inlined
uint32_t Bgen13ConvertBiallelicPhased(uint32_t sample_idx_lowbits, uintptr_t numer_a1, uintptr_t numer_a2, uintptr_t numer_mask, uint32_t numer_certainty_min, uint64_t magic_preadd, uint64_t magic_mult, uint32_t magic_postshift, uint32_t hard_call_halfdist, uint32_t dosage_erase_halfdist, uintptr_t* genovec_word_ptr, uint32_t* __restrict phasepresent_hw_ptr, uint32_t* __restrict phaseinfo_hw_ptr, uint32_t* __restrict dphase_present_hw_ptr, SDosage** dphase_delta_iterp, uint32_t* write_dosage_int_ptr) {
  const uint32_t sample_idx_lowbits_x2 = sample_idx_lowbits * 2;
  if (numer_certainty_min) {
    // bugfix (29 Apr 2019): This check was not implemented correctly.
    //
    // Correct VCF HDS-import code:
    //   dist_from_1 = fabs(1.0 - dosage_sum);
    //   if ((1.0 - dist_from_1 <= import_dosage_certainty) && (dist_from_1 <= import_dosage_certainty)) {
    //     (force missing)
    //   }
    const uintptr_t numer_sum = numer_a1 + numer_a2;
    const uint32_t dist_from_1 = abs_i32(numer_sum - numer_mask);
    if ((dist_from_1 < numer_certainty_min) && (numer_mask - dist_from_1 < numer_certainty_min)) {
      *genovec_word_ptr |= (3 * k1LU) << sample_idx_lowbits_x2;
      return 1;
    }
  }
  const uint32_t write_dosage_int = (magic_preadd + magic_mult * (numer_a1 + numer_a2)) >> magic_postshift;
  if (numer_a1 == numer_a2) {
    const uint32_t halfdist = BiallelicDosageHalfdist(write_dosage_int);
    if (halfdist < hard_call_halfdist) {
      *genovec_word_ptr |= (3 * k1LU) << sample_idx_lowbits_x2;
    } else {
      *genovec_word_ptr |= ((write_dosage_int + (kDosage4th * k1LU)) / kDosageMid) << sample_idx_lowbits_x2;
      if (halfdist >= dosage_erase_halfdist) {
        return 1;
      }
    }
    *write_dosage_int_ptr = write_dosage_int;
    return 0;
  }
  const uint32_t write_dphase_delta_p1 = (magic_preadd + magic_mult * (numer_a1 + numer_mask - numer_a2)) >> magic_postshift;
  int32_t write_dphase_delta_val = write_dphase_delta_p1 - kDosageMid;

  // now on 0..32768 scale
  const uint32_t dphase_side1 = write_dosage_int + write_dphase_delta_val;
  const uint32_t dphase_side2 = write_dosage_int - write_dphase_delta_val;

  const uint32_t dphase_halfdist1 = DphaseHalfdist(dphase_side1);
  const uint32_t dphase_halfdist2 = DphaseHalfdist(dphase_side2);
  const uint32_t dphase_erase_halfdist = dosage_erase_halfdist + kDosage4th;
  if ((dphase_halfdist1 >= dphase_erase_halfdist) && (dphase_halfdist2 >= dphase_erase_halfdist)) {
    // --dosage-erase-threshold only applies to phased dosage when both sides
    // are within [threshold/2] of an integer; otherwise we have a
    // discontinuity between
    //   dosage=0.103=0.0515|0.0515 (dphase omitted) and
    //   dosage=0.103=0.051|0.052.
    // Note that, since we don't allow the --dosage-erase-threshold value
    // to exceed the --hard-call-threshold value, the dosage component of
    // 0|1:0.85 cannot be erased by any --dosage-erase-threshold setting when
    // importing with --hard-call-threshold 0.2; it is necessary to use e.g.
    // "--make-pgen erase-dosage" instead.
    const uintptr_t geno1 = (dphase_side1 + kDosageMid) / kDosageMax;
    const uintptr_t geno2 = (dphase_side2 + kDosageMid) / kDosageMax;
    uintptr_t cur_geno = geno1 + geno2;
    // This duplicate code isn't strictly necessary, but I expect it to pay off
    // since all phased-het hardcalls go here.
    *genovec_word_ptr |= cur_geno << sample_idx_lowbits_x2;
    if (cur_geno == 1) {
      *phasepresent_hw_ptr |= 1U << sample_idx_lowbits;
      *phaseinfo_hw_ptr |= geno1 << sample_idx_lowbits;
    }
    return 1;
  }
  const uint32_t halfdist = BiallelicDosageHalfdist(write_dosage_int);
  // reference allele assumed to be second here, so haplotype dosages are for
  // alt1.
  if (halfdist < hard_call_halfdist) {
    // --hard-call-threshold must not care about phased dosage; otherwise we
    // have a discontinuity between
    //   dosage=0.998=0.499|0.499 (dphase omitted) and
    //   0.998=0.498|0.5.
    *genovec_word_ptr |= (3 * k1LU) << sample_idx_lowbits_x2;
  } else {
    const uintptr_t cur_geno = (write_dosage_int + (kDosage4th * k1LU)) / kDosageMid;
    *genovec_word_ptr |= cur_geno << sample_idx_lowbits_x2;
    if (cur_geno == 1) {
      // generate hardcall-phase from dosage-phase iff delta > 0.5.
      // bugfix (22 Apr 2018): (write_dphase_delta_val >> 31) is actually
      // undefined behavior for an int32_t; must cast to uint32_t first
      const uint32_t neg_sign_bit = -(S_CAST(uint32_t, write_dphase_delta_val) >> 31);
      const uint32_t abs_write_dphase_delta_val = (S_CAST(uint32_t, write_dphase_delta_val) ^ neg_sign_bit) - neg_sign_bit;
      if (abs_write_dphase_delta_val > kDosage4th) {
        *phasepresent_hw_ptr |= 1U << sample_idx_lowbits;
        // phaseinfo_hw bit set for 1|0
        // todo: check if branch is faster
        *phaseinfo_hw_ptr |= (neg_sign_bit + 1) << sample_idx_lowbits;
        if ((abs_write_dphase_delta_val == write_dosage_int) || (abs_write_dphase_delta_val + write_dosage_int == kDosageMax)) {
          // can omit explicit dphase_delta in this case
          write_dphase_delta_val = 0;
        }
      } else if ((halfdist >= dosage_erase_halfdist) && (!write_dphase_delta_val)) {
        // unphased het hardcall special case
        return 1;
      }
    }
  }
  // we should never get here in bit_precision == 1 case
  if (write_dphase_delta_val != 0) {
    *dphase_present_hw_ptr |= 1U << sample_idx_lowbits;
    **dphase_delta_iterp = write_dphase_delta_val;
    *dphase_delta_iterp += 1;
  }
  *write_dosage_int_ptr = write_dosage_int;
  return 0;
}

typedef struct Bgen13GenoToPgenCtxStruct {
  BgenImportCommon* common;
  uint32_t error_on_polyploid;
  uint32_t hard_call_halfdist;
  uint32_t* bgen_import_dosage_certainty_thresholds;
  uint32_t prov_ref_allele_second;

  unsigned char** thread_wkspaces;
  uint32_t* thread_bidxs[2];
  GparseRecord* gparse[2];
  const uintptr_t* block_allele_idx_offsets[2];

  const char** err_extra;

  uint64_t err_info;
} Bgen13GenoToPgenCtx;

static_assert(sizeof(Dosage) == 2, "Bgen13GenoToPgenThread() needs to be updated.");
THREAD_FUNC_DECL Bgen13GenoToPgenThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  Bgen13GenoToPgenCtx* ctx = S_CAST(Bgen13GenoToPgenCtx*, arg->sharedp->context);
  BgenImportCommon* bicp = ctx->common;

  const uintptr_t sample_ct = bicp->sample_ct;
  const uint32_t hard_call_halfdist = ctx->hard_call_halfdist;
  const uint32_t dosage_erase_halfdist = bicp->dosage_erase_halfdist;
  const uint32_t* bgen_import_dosage_certainty_thresholds = ctx->bgen_import_dosage_certainty_thresholds;
  const uint32_t compression_mode = bicp->compression_mode;
  const uint32_t prov_ref_allele_second = ctx->prov_ref_allele_second;
  const uint32_t sample_ctl2_m1 = (sample_ct - 1) / kBitsPerWordD2;
  const uintptr_t sample_ctl = BitCtToWordCt(sample_ct);
  const unsigned char* cur_uncompressed_geno = ctx->thread_wkspaces[tidx];
  struct libdeflate_decompressor* decompressor = bicp->libdeflate_decompressors[tidx];
  uintptr_t* patch_01_set = nullptr;
  AlleleCode* patch_01_vals = nullptr;
  uintptr_t* patch_10_set = nullptr;
  AlleleCode* patch_10_vals = nullptr;
  uintptr_t* __attribute__((may_alias)) phasepresent = nullptr;
  uintptr_t* __attribute__((may_alias)) phaseinfo = nullptr;
  uintptr_t* __attribute__((may_alias)) dosage_present = nullptr;
  Dosage* dosage_main = nullptr;
  uintptr_t* __attribute__((may_alias)) dphase_present = nullptr;
  SDosage* dphase_delta = nullptr;
  uint32_t cur_allele_ct = 2;
  uint32_t parity = 0;
  uint64_t new_err_info = 0;
  do {
    {
      const uintptr_t* block_allele_idx_offsets = ctx->block_allele_idx_offsets[parity];
      const uint32_t bidx_end = ctx->thread_bidxs[parity][tidx + 1];
      uint32_t bidx = ctx->thread_bidxs[parity][tidx];
      GparseRecord* cur_gparse = ctx->gparse[parity];

      for (; bidx != bidx_end; ++bidx) {
        GparseRecord* grp = &(cur_gparse[bidx]);
        const uint32_t compressed_byte_ct = grp->metadata.read_bgen.input_byte_ct;
        uint32_t uncompressed_byte_ct;
        if (compression_mode) {
          uncompressed_byte_ct = grp->metadata.read_bgen.uncompressed_byte_ct;
          if (compression_mode == 1) {
            if (unlikely(libdeflate_zlib_decompress(decompressor, grp->record_start, compressed_byte_ct, K_CAST(unsigned char*, cur_uncompressed_geno), uncompressed_byte_ct, nullptr) != LIBDEFLATE_SUCCESS)) {
              new_err_info = (S_CAST(uint64_t, bidx) << 32) | (S_CAST(uint32_t, kBgenImportErrSubtypeLibdeflateZdecompress) << 8) | S_CAST(uint32_t, kPglRetMalformedInput);
              goto Bgen13GenoToPgenThread_err;
            }
          } else {
            const uintptr_t extracted_byte_ct = ZSTD_decompress(K_CAST(unsigned char*, cur_uncompressed_geno), uncompressed_byte_ct, grp->record_start, compressed_byte_ct);
            if (unlikely(extracted_byte_ct != uncompressed_byte_ct)) {
              assert(ZSTD_isError(extracted_byte_ct));
              ctx->err_extra[tidx] = ZSTD_getErrorName(extracted_byte_ct);
              new_err_info = (S_CAST(uint64_t, bidx) << 32) | (tidx << 16) | (S_CAST(uint32_t, kBgenImportErrSubtypeZstdDecompress) << 8) | S_CAST(uint32_t, kPglRetMalformedInput);
              goto Bgen13GenoToPgenThread_err;
            }
          }
        } else {
          uncompressed_byte_ct = compressed_byte_ct;
          memcpy(K_CAST(unsigned char*, cur_uncompressed_geno), grp->record_start, compressed_byte_ct);
        }
        // 4 bytes: sample_ct
        // 2 bytes: # of alleles, must match bgen_allele_cts[bidx]
        // 1 byte: min ploidy
        // 1 byte: max ploidy
        // sample_ct bytes: low 6 bits = ploidy, top bit = missingness
        // 1 byte: 1 if phased, 0 if not
        // 1 byte: # of bits of probability precision (we just support 8 and 16
        //         for now, add others later)
        uint32_t stored_sample_ct;
        memcpy(&stored_sample_ct, cur_uncompressed_geno, sizeof(int32_t));
        if (unlikely((uncompressed_byte_ct < 10 + sample_ct) || (sample_ct != stored_sample_ct))) {
          new_err_info = (S_CAST(uint64_t, bidx) << 32) | (S_CAST(uint32_t, kBgenImportErrSubtypeSampleCtMismatch) << 8) | S_CAST(uint32_t, kPglRetMalformedInput);
          goto Bgen13GenoToPgenThread_err;
        }
        if (block_allele_idx_offsets) {
          cur_allele_ct = block_allele_idx_offsets[bidx + 1] - block_allele_idx_offsets[bidx];
        }
        /*
        uint16_t stored_allele_ct;
        memcpy_k(&stored_allele_ct, &(cur_uncompressed_geno[4]), sizeof(int16_t));
        if (unlikely(stored_allele_ct != cur_allele_ct)) {
          goto Bgen13GenoToPgenThread_malformed;
        }
        */
        const uint32_t min_ploidy = cur_uncompressed_geno[6];
        const uint32_t max_ploidy = cur_uncompressed_geno[7];
        if (unlikely((min_ploidy > max_ploidy) || (max_ploidy > 63))) {
          new_err_info = (S_CAST(uint64_t, bidx) << 32) | (S_CAST(uint32_t, kBgenImportErrSubtypePloidyOutOfRange) << 8) | S_CAST(uint32_t, kPglRetMalformedInput);
          goto Bgen13GenoToPgenThread_err;
        }
        if (unlikely(max_ploidy > 2)) {
          if (ctx->error_on_polyploid) {
            new_err_info = (S_CAST(uint64_t, bidx) << 32) | (S_CAST(uint32_t, kBgenImportErrSubtypePolyploid) << 8) | S_CAST(uint32_t, kPglRetInconsistentInput);
            goto Bgen13GenoToPgenThread_err;
          }
          new_err_info = (S_CAST(uint64_t, bidx) << 32) | (S_CAST(uint32_t, kBgenImportErrSubtypePolyploid) << 8) | S_CAST(uint32_t, kPglRetNotYetSupported);
          goto Bgen13GenoToPgenThread_err;
        }
        const unsigned char* missing_and_ploidy_iter = &(cur_uncompressed_geno[8]);
        const unsigned char* probs_start = &(cur_uncompressed_geno[10 + sample_ct]);
        const uint32_t is_phased = probs_start[-2];
        if (unlikely(is_phased > 1)) {
          new_err_info = (S_CAST(uint64_t, bidx) << 32) | (S_CAST(uint32_t, kBgenImportErrSubtypeIsPhasedOutOfRange) << 8) | S_CAST(uint32_t, kPglRetMalformedInput);
          goto Bgen13GenoToPgenThread_err;
        }
        const uint32_t bit_precision = probs_start[-1];
        if (unlikely((!bit_precision) || (bit_precision > 32))) {
          new_err_info = (S_CAST(uint64_t, bidx) << 32) | (S_CAST(uint32_t, kBgenImportErrSubtypeBitPrecisionOutOfRange) << 8) | S_CAST(uint32_t, kPglRetMalformedInput);
          goto Bgen13GenoToPgenThread_err;
        }
        if (unlikely(bit_precision > kMaxBgenImportBits)) {
          new_err_info = (S_CAST(uint64_t, bidx) << 32) | (S_CAST(uint32_t, kBgenImportErrSubtypeBitPrecisionUnsupported) << 8) | S_CAST(uint32_t, kPglRetNotYetSupported);
          goto Bgen13GenoToPgenThread_err;
        }
        if (cur_allele_ct != 2) {
          // shouldn't be possible to get here for now
          assert(0);
          new_err_info = (S_CAST(uint64_t, bidx) << 32) | (S_CAST(uint32_t, kBgenImportErrSubtypeAlleleCtUnsupported) << 8) | S_CAST(uint32_t, kPglRetNotYetSupported);
          goto Bgen13GenoToPgenThread_err;
        }
        uintptr_t* genovec = GparseGetPointers(grp->record_start, sample_ct, cur_allele_ct, grp->flags, &patch_01_set, &patch_01_vals, &patch_10_set, &patch_10_vals, &phasepresent, &phaseinfo, &dosage_present, &dosage_main, &dphase_present, &dphase_delta);
        Dosage* dosage_main_iter = dosage_main;
        SDosage* dphase_delta_iter = dphase_delta;
        uint32_t cur_phasepresent_exists = 0;
        // turns out UKB haplotype files still use bit_precision == 16, so
        // don't bother implementing movemask-based bit_precision == 1
        // optimization for now
        const uint64_t magic_preadd = kBgenMagicNums[bit_precision].preadd;
        const uint64_t magic_mult = kBgenMagicNums[bit_precision].mult;
        const uint32_t magic_postshift = kBgenMagicNums[bit_precision].postshift;

        // also equal to denominator
        const uintptr_t numer_mask = (1U << bit_precision) - 1;

        uint32_t numer_certainty_min = 0;
        if (bgen_import_dosage_certainty_thresholds) {
          numer_certainty_min = bgen_import_dosage_certainty_thresholds[bit_precision];
        }

        uint32_t inner_loop_last = kBitsPerWordD2 - 1;
        if (min_ploidy == max_ploidy) {
          // faster handling of common cases (no need to keep checking if
          // we've read past the end)
          if (unlikely(uncompressed_byte_ct != 10 + sample_ct + DivUp(S_CAST(uint64_t, bit_precision) * max_ploidy * sample_ct, CHAR_BIT))) {
            new_err_info = (S_CAST(uint64_t, bidx) << 32) | (S_CAST(uint32_t, kBgenImportErrSubtypeUncompressedByteCtMismatch) << 8) | S_CAST(uint32_t, kPglRetMalformedInput);
            goto Bgen13GenoToPgenThread_err;
          }
          if (max_ploidy == 2) {
            if (!is_phased) {
              // Could start with vectorized scan which proceeds until first
              // inexact dosage?  Logic needs to be slightly different from
              // scan thread since we need to look at missing_and_ploidy.
              // On the other hand, decompression already seems to be ~3/4 of
              // compute time here, so I won't worry about that optimization
              // for now.
              for (uint32_t widx = 0; ; ++widx) {
                if (widx >= sample_ctl2_m1) {
                  if (widx > sample_ctl2_m1) {
                    break;
                  }
                  inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
                }
                const uintptr_t prob_offset_base = widx * kBitsPerWord;
                uintptr_t genovec_word = 0;
                uint32_t dosage_present_hw = 0;
                for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
                  const uint32_t missing_and_ploidy = *missing_and_ploidy_iter++;
                  if (missing_and_ploidy != 2) {
                    // (could also validate that missing_and_ploidy == 130)
                  Bgen13GenoToPgenThread_diploid_unphased_missing:
                    genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
                    continue;
                  }
                  uintptr_t numer_aa;
                  uintptr_t numer_ab;
                  Bgen13GetTwoVals(probs_start, prob_offset_base + 2 * sample_idx_lowbits, bit_precision, numer_mask, &numer_aa, &numer_ab);
                  // common trivial cases
                  if (!numer_aa) {
                    if (!numer_ab) {
                      continue;
                    }
                    if (numer_ab == numer_mask) {
                      genovec_word |= k1LU << (2 * sample_idx_lowbits);
                      continue;
                    }
                  } else if ((numer_aa == numer_mask) && (!numer_ab)) {
                    genovec_word |= (2 * k1LU) << (2 * sample_idx_lowbits);
                    continue;
                  }
                  if (unlikely(numer_aa + numer_ab > numer_mask)) {
                    new_err_info = (S_CAST(uint64_t, bidx) << 32) | (S_CAST(uint32_t, kBgenImportErrSubtypeNumeratorOverflow) << 8) | S_CAST(uint32_t, kPglRetMalformedInput);
                    goto Bgen13GenoToPgenThread_err;
                  }
                  if ((numer_aa < numer_certainty_min) && (numer_ab < numer_certainty_min) && (numer_mask - numer_certainty_min < numer_aa + numer_ab)) {
                    // missing due to --import-dosage-certainty
                    goto Bgen13GenoToPgenThread_diploid_unphased_missing;
                  }
                  const uint32_t write_dosage_int = (magic_preadd + magic_mult * (2 * numer_aa + numer_ab)) >> magic_postshift;
                  const uint32_t halfdist = BiallelicDosageHalfdist(write_dosage_int);
                  if (halfdist < hard_call_halfdist) {
                    genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
                  } else {
                    genovec_word |= ((write_dosage_int + (kDosage4th * k1LU)) / kDosageMid) << (2 * sample_idx_lowbits);
                    if (halfdist >= dosage_erase_halfdist) {
                      continue;
                    }
                  }
                  dosage_present_hw |= 1U << sample_idx_lowbits;
                  *dosage_main_iter++ = write_dosage_int;
                }
                genovec[widx] = genovec_word;
                R_CAST(Halfword*, dosage_present)[widx] = dosage_present_hw;
              }
            } else {
              // biallelic, phased, diploid
              // may want a fast path for bit_precision == 1
              phasepresent[sample_ctl - 1] = 0;
              for (uint32_t widx = 0; ; ++widx) {
                if (widx >= sample_ctl2_m1) {
                  if (widx > sample_ctl2_m1) {
                    break;
                  }
                  inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
                }
                const uintptr_t prob_offset_base = widx * kBitsPerWord;
                uintptr_t genovec_word = 0;
                uint32_t phasepresent_hw = 0;
                uint32_t phaseinfo_hw = 0;
                uint32_t dosage_present_hw = 0;
                uint32_t dphase_present_hw = 0;
                for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
                  const uint32_t missing_and_ploidy = *missing_and_ploidy_iter++;
                  if (missing_and_ploidy != 2) {
                    genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
                    continue;
                  }
                  uintptr_t numer_a1;
                  uintptr_t numer_a2;
                  Bgen13GetTwoVals(probs_start, prob_offset_base + 2 * sample_idx_lowbits, bit_precision, numer_mask, &numer_a1, &numer_a2);
                  if ((!numer_a1) && (!numer_a2)) {
                    continue;
                  }
                  if ((numer_a1 == numer_mask) && (numer_a2 == numer_mask)) {
                    genovec_word |= (2 * k1LU) << (2 * sample_idx_lowbits);
                    continue;
                  }
                  uint32_t write_dosage_int;
                  if (Bgen13ConvertBiallelicPhased(sample_idx_lowbits, numer_a1, numer_a2, numer_mask, numer_certainty_min, magic_preadd, magic_mult, magic_postshift, hard_call_halfdist, dosage_erase_halfdist, &genovec_word, &phasepresent_hw, &phaseinfo_hw, &dphase_present_hw, &dphase_delta_iter, &write_dosage_int)) {
                    continue;
                  }
                  dosage_present_hw |= 1U << sample_idx_lowbits;
                  *dosage_main_iter++ = write_dosage_int;
                }
                genovec[widx] = genovec_word;
                R_CAST(Halfword*, phasepresent)[widx] = phasepresent_hw;
                R_CAST(Halfword*, phaseinfo)[widx] = phaseinfo_hw;
                R_CAST(Halfword*, dosage_present)[widx] = dosage_present_hw;
                R_CAST(Halfword*, dphase_present)[widx] = dphase_present_hw;
              }
              cur_phasepresent_exists = !AllWordsAreZero(phasepresent, sample_ctl);
            }
          } else if (max_ploidy) {
            // biallelic, haploid
            for (uint32_t widx = 0; ; ++widx) {
              if (widx >= sample_ctl2_m1) {
                if (widx > sample_ctl2_m1) {
                  break;
                }
                inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
              }
              const uintptr_t prob_offset_base = widx * kBitsPerWordD2;
              uintptr_t genovec_word = 0;
              uint32_t dosage_present_hw = 0;
              for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
                const uint32_t missing_and_ploidy = *missing_and_ploidy_iter++;
                if (missing_and_ploidy != 1) {
                Bgen13GenoToPgenThread_haploid_missing:
                  genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
                  continue;
                }
                const uintptr_t numer_a = Bgen13GetOneVal(probs_start, prob_offset_base + sample_idx_lowbits, bit_precision, numer_mask);
                if ((numer_a < numer_certainty_min) && (numer_mask - numer_certainty_min < numer_a)) {
                  goto Bgen13GenoToPgenThread_haploid_missing;
                }
                const uint32_t write_dosage_int = (magic_preadd + magic_mult * numer_a * 2) >> magic_postshift;
                const uint32_t halfdist = BiallelicDosageHalfdist(write_dosage_int);
                if (halfdist < hard_call_halfdist) {
                  genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
                } else {
                  genovec_word |= ((write_dosage_int + (kDosage4th * k1LU)) / kDosageMid) << (2 * sample_idx_lowbits);
                  if (halfdist >= dosage_erase_halfdist) {
                    continue;
                  }
                }
                dosage_present_hw |= 1U << sample_idx_lowbits;
                *dosage_main_iter++ = write_dosage_int;
              }
              genovec[widx] = genovec_word;
              R_CAST(Halfword*, dosage_present)[widx] = dosage_present_hw;
            }
          } else {
            // all-ploidy-0, set everything to missing
            SetAllBits(sample_ct * 2, genovec);
          }
        } else {
          // biallelic, variable ploidy
          const uint64_t remaining_bit_ct = 8LLU * (uncompressed_byte_ct - S_CAST(uintptr_t, probs_start - cur_uncompressed_geno));
          const uint64_t prob_offset_end = remaining_bit_ct / bit_precision;
          uintptr_t prob_offset = 0;
          if (!is_phased) {
            for (uint32_t widx = 0; ; ++widx) {
              if (widx >= sample_ctl2_m1) {
                if (widx > sample_ctl2_m1) {
                  break;
                }
                inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
              }
              uintptr_t genovec_word = 0;
              uint32_t dosage_present_hw = 0;
              for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
                uint32_t missing_and_ploidy = *missing_and_ploidy_iter++;
                uint32_t write_dosage_int;
                if (missing_and_ploidy == 2) {
                  if (unlikely(prob_offset + 2 > prob_offset_end)) {
                    new_err_info = (S_CAST(uint64_t, bidx) << 32) | (S_CAST(uint32_t, kBgenImportErrSubtypeProbOffsetMismatch) << 8) | S_CAST(uint32_t, kPglRetMalformedInput);
                    goto Bgen13GenoToPgenThread_err;
                  }
                  uintptr_t numer_aa;
                  uintptr_t numer_ab;
                  Bgen13GetTwoVals(probs_start, prob_offset, bit_precision, numer_mask, &numer_aa, &numer_ab);
                  prob_offset += 2;
                  if (unlikely(numer_aa + numer_ab > numer_mask)) {
                    new_err_info = (S_CAST(uint64_t, bidx) << 32) | (S_CAST(uint32_t, kBgenImportErrSubtypeNumeratorOverflow) << 8) | S_CAST(uint32_t, kPglRetMalformedInput);
                    goto Bgen13GenoToPgenThread_err;
                  }
                  if ((numer_aa < numer_certainty_min) && (numer_ab < numer_certainty_min) && (numer_mask - numer_certainty_min < numer_aa + numer_ab)) {
                    // missing due to --import-dosage-certainty
                    goto Bgen13GenoToPgenThread_generic_unphased_missing;
                  }
                  write_dosage_int = (magic_preadd + magic_mult * (2 * numer_aa + numer_ab)) >> magic_postshift;
                } else {
                  if (missing_and_ploidy != 1) {
                    missing_and_ploidy &= 127;
                    if (unlikely(missing_and_ploidy > 2)) {
                      new_err_info = (S_CAST(uint64_t, bidx) << 32) | (S_CAST(uint32_t, kBgenImportErrSubtypeInvalidMissingPloidy) << 8) | S_CAST(uint32_t, kPglRetMalformedInput);
                      goto Bgen13GenoToPgenThread_err;
                    }
                    prob_offset += missing_and_ploidy;
                  Bgen13GenoToPgenThread_generic_unphased_missing:
                    genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
                    continue;
                  }
                  if (unlikely(prob_offset >= prob_offset_end)) {
                    new_err_info = (S_CAST(uint64_t, bidx) << 32) | (S_CAST(uint32_t, kBgenImportErrSubtypeProbOffsetMismatch) << 8) | S_CAST(uint32_t, kPglRetMalformedInput);
                    goto Bgen13GenoToPgenThread_err;
                  }
                  const uintptr_t numer_a = Bgen13GetOneVal(probs_start, prob_offset, bit_precision, numer_mask);
                  ++prob_offset;
                  if ((numer_a < numer_certainty_min) && (numer_mask - numer_certainty_min < numer_a)) {
                    goto Bgen13GenoToPgenThread_generic_unphased_missing;
                  }
                  write_dosage_int = (magic_preadd + magic_mult * numer_a * 2) >> magic_postshift;
                }
                const uint32_t halfdist = BiallelicDosageHalfdist(write_dosage_int);
                if (halfdist < hard_call_halfdist) {
                  genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
                } else {
                  genovec_word |= ((write_dosage_int + (kDosage4th * k1LU)) / kDosageMid) << (2 * sample_idx_lowbits);
                  if (halfdist >= dosage_erase_halfdist) {
                    continue;
                  }
                }
                dosage_present_hw |= 1U << sample_idx_lowbits;
                *dosage_main_iter++ = write_dosage_int;
              }
              genovec[widx] = genovec_word;
              R_CAST(Halfword*, dosage_present)[widx] = dosage_present_hw;
            }
          } else {
            // biallelic, phased, variable ploidy
            phasepresent[sample_ctl - 1] = 0;
            for (uint32_t widx = 0; ; ++widx) {
              if (widx >= sample_ctl2_m1) {
                if (widx > sample_ctl2_m1) {
                  break;
                }
                inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
              }
              uintptr_t genovec_word = 0;
              uint32_t phasepresent_hw = 0;
              uint32_t phaseinfo_hw = 0;
              uint32_t dosage_present_hw = 0;
              uint32_t dphase_present_hw = 0;
              uint32_t write_dosage_int;
              for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
                uint32_t missing_and_ploidy = *missing_and_ploidy_iter++;
                if (missing_and_ploidy == 2) {
                  if (unlikely(prob_offset + 2 > prob_offset_end)) {
                    new_err_info = (S_CAST(uint64_t, bidx) << 32) | (S_CAST(uint32_t, kBgenImportErrSubtypeProbOffsetMismatch) << 8) | S_CAST(uint32_t, kPglRetMalformedInput);
                    goto Bgen13GenoToPgenThread_err;
                  }
                  uintptr_t numer_a1;
                  uintptr_t numer_a2;
                  Bgen13GetTwoVals(probs_start, prob_offset, bit_precision, numer_mask, &numer_a1, &numer_a2);
                  prob_offset += 2;
                  if ((!numer_a1) && (!numer_a2)) {
                    continue;
                  }
                  if ((numer_a1 == numer_mask) && (numer_a2 == numer_mask)) {
                    genovec_word |= (2 * k1LU) << (2 * sample_idx_lowbits);
                    continue;
                  }
                  if (Bgen13ConvertBiallelicPhased(sample_idx_lowbits, numer_a1, numer_a2, numer_mask, numer_certainty_min, magic_preadd, magic_mult, magic_postshift, hard_call_halfdist, dosage_erase_halfdist, &genovec_word, &phasepresent_hw, &phaseinfo_hw, &dphase_present_hw, &dphase_delta_iter, &write_dosage_int)) {
                    continue;
                  }
                } else {
                  if (missing_and_ploidy != 1) {
                    missing_and_ploidy &= 127;
                    if (unlikely(missing_and_ploidy > 2)) {
                      new_err_info = (S_CAST(uint64_t, bidx) << 32) | (S_CAST(uint32_t, kBgenImportErrSubtypeInvalidMissingPloidy) << 8) | S_CAST(uint32_t, kPglRetMalformedInput);
                      goto Bgen13GenoToPgenThread_err;
                    }
                    prob_offset += missing_and_ploidy;
                  Bgen13GenoToPgenThread_generic_phased_missing:
                    genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
                    continue;
                  }
                  if (unlikely(prob_offset >= prob_offset_end)) {
                    new_err_info = (S_CAST(uint64_t, bidx) << 32) | (S_CAST(uint32_t, kBgenImportErrSubtypeProbOffsetMismatch) << 8) | S_CAST(uint32_t, kPglRetMalformedInput);
                    goto Bgen13GenoToPgenThread_err;
                  }
                  const uintptr_t numer_a = Bgen13GetOneVal(probs_start, prob_offset, bit_precision, numer_mask);
                  ++prob_offset;
                  if ((numer_a < numer_certainty_min) && (numer_mask - numer_certainty_min < numer_a)) {
                    goto Bgen13GenoToPgenThread_generic_phased_missing;
                  }
                  write_dosage_int = (magic_preadd + magic_mult * numer_a * 2) >> magic_postshift;
                  const uint32_t halfdist = BiallelicDosageHalfdist(write_dosage_int);
                  if (halfdist < hard_call_halfdist) {
                    genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
                  } else {
                    genovec_word |= ((write_dosage_int + (kDosage4th * k1LU)) / kDosageMid) << (2 * sample_idx_lowbits);
                    // bugfix (22 Apr 2018): forgot this check
                    if (halfdist >= dosage_erase_halfdist) {
                      continue;
                    }
                  }
                }
                dosage_present_hw |= 1U << sample_idx_lowbits;
                *dosage_main_iter++ = write_dosage_int;
              }
              genovec[widx] = genovec_word;
              // bugfix (22 Apr 2018): didn't save phasepresent, etc.
              R_CAST(Halfword*, phasepresent)[widx] = phasepresent_hw;
              R_CAST(Halfword*, phaseinfo)[widx] = phaseinfo_hw;
              R_CAST(Halfword*, dosage_present)[widx] = dosage_present_hw;
              R_CAST(Halfword*, dphase_present)[widx] = dphase_present_hw;
            }
            cur_phasepresent_exists = !AllWordsAreZero(phasepresent, sample_ctl);
          }
        }
        if (!(sample_ctl2_m1 % 2)) {
          // needed for dosage_present.  don't think it's needed for
          // dphase_present, but play that safe for now
          R_CAST(Halfword*, dosage_present)[sample_ctl2_m1 + 1] = 0;
          R_CAST(Halfword*, dphase_present)[sample_ctl2_m1 + 1] = 0;
        }
        const uint32_t dosage_ct = dosage_main_iter - dosage_main;
        uint32_t dphase_ct = dphase_delta_iter - dphase_delta;
        // note that this is inverted from bgen-1.1
        if (!prov_ref_allele_second) {
          GenovecInvertUnsafe(sample_ct, genovec);
          ZeroTrailingNyps(sample_ct, genovec);
          if (dosage_ct) {
            BiallelicDosage16Invert(dosage_ct, dosage_main);
            // currently no code path here where dosage_ct < dphase_ct
            if (dphase_ct) {
              BiallelicDphase16Invert(dphase_ct, dphase_delta);
            }
          }
        }
        grp->metadata.write.patch_01_ct = 0;
        grp->metadata.write.patch_10_ct = 0;
        grp->metadata.write.phasepresent_exists = cur_phasepresent_exists;
        grp->metadata.write.dosage_ct = dosage_ct;
        grp->metadata.write.multiallelic_dosage_ct = 0;
        grp->metadata.write.dphase_ct = dphase_ct;
        grp->metadata.write.multiallelic_dphase_ct = 0;
      }
    }
    parity = 1 - parity;
  } while (!THREAD_BLOCK_FINISH(arg));
  while (0) {
  Bgen13GenoToPgenThread_err:
    UpdateU64IfSmaller(new_err_info, &ctx->err_info);
    THREAD_BLOCK_FINISH(arg);
    break;
  }
  THREAD_RETURN;
}

static_assert(sizeof(Dosage) == 2, "OxBgenToPgen() needs to be updated.");
PglErr OxBgenToPgen(const char* bgenname, const char* samplename, const char* const_fid, const char* ox_single_chr_str, const char* ox_missing_code, const char* missing_catname, MiscFlags misc_flags, ImportFlags import_flags, LoadFilterLogFlags load_filter_log_import_flags, OxfordImportFlags oxford_import_flags, uint32_t psam_01, uint32_t is_update_or_impute_sex, uint32_t is_splitpar, uint32_t is_sortvars, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, double import_dosage_certainty, char id_delim, char idspace_to, uint32_t import_max_allele_ct, uint32_t max_thread_ct, char* outname, char* outname_end, ChrInfo* cip) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  FILE* bgenfile = nullptr;

  // only if no sample file specified, and .bgen has sample IDs.  (possible
  // todo: consistency check when both sources of sample IDs are present?)
  FILE* psamfile = nullptr;

  char* pvar_cswritep = nullptr;
  const uint32_t sex_info_avail = samplename[0] || is_update_or_impute_sex;
  CompressStreamState pvar_css;
  PreinitCstream(&pvar_css);
  ThreadGroup tg;
  PreinitThreads(&tg);
  STPgenWriter spgw;
  PglErr reterr = kPglRetSuccess;
  PreinitSpgw(&spgw);
  BgenImportCommon common;
  common.libdeflate_decompressors = nullptr;
  {
    // Pass 1: Determine whether there's at least one non-hardcall needs to be
    //         saved, and if a chromosome filter was specified, count the
    //         number of variants which pass the filter.
    //         For bgen-1.2/1.3, the .pvar is also written in this pass.
    //         For bgen-1.1, we can usually early-bail when no chromosome
    //         filter is involved, so .pvar writing is postponed till the
    //         second pass.
    // Pass 2: Write .pgen file.
    //
    // Possible todo: reduce this to 1-pass in some circumstances.  Punt for
    // now, though, since the current workflow makes good use of the scanning
    // pass to bound conversion-pass per-thread buffer sizes.
    if (unlikely(fopen_checked(bgenname, FOPEN_RB, &bgenfile))) {
      goto OxBgenToPgen_ret_OPEN_FAIL;
    }
    uint32_t initial_uints[5];
    if (unlikely(!fread_unlocked(&initial_uints[0], 20, 1, bgenfile))) {
      // this could be malformed input as well; could distinguish later?
      goto OxBgenToPgen_ret_READ_FAIL;
    }
    if (unlikely(initial_uints[1] > initial_uints[0])) {
      logerrputs("Error: Invalid .bgen header.\n");
      goto OxBgenToPgen_ret_MALFORMED_INPUT;
    }
    const uint32_t sample_ct = initial_uints[3];
    if (unlikely(initial_uints[4] && (initial_uints[4] != 0x6e656762))) {
      logerrputs("Error: Invalid .bgen magic number.\n");
      goto OxBgenToPgen_ret_MALFORMED_INPUT;
    }
    const uint32_t header_variant_ct = initial_uints[2];
    if (unlikely(!header_variant_ct)) {
      logerrputs("Error: Empty .bgen file.\n");
      goto OxBgenToPgen_ret_DEGENERATE_DATA;
    }

    if (unlikely(fseeko(bgenfile, initial_uints[1], SEEK_SET))) {
      goto OxBgenToPgen_ret_READ_FAIL;
    }
    uint32_t header_flags;
    if (unlikely(!fread_unlocked(&header_flags, 4, 1, bgenfile))) {
      goto OxBgenToPgen_ret_READ_FAIL;
    }
    const uint32_t compression_mode = header_flags & 3;
    const uint32_t layout = (header_flags >> 2) & 15;
    if (unlikely(!layout)) {
      logerrputs("Error: BGEN v1.0 files are not supported by " PROG_NAME_STR ".\n");
      goto OxBgenToPgen_ret_MALFORMED_INPUT;
    }
    if (unlikely((compression_mode == 3) || (layout > 2))) {
      logerrputs("Error: Unrecognized BGEN version.  Use gen-convert or a similar tool to\ndowncode to BGEN v1.3 if you want to process this data with " PROG_NAME_STR ".\n");
      goto OxBgenToPgen_ret_MALFORMED_INPUT;
    }
    if (unlikely((compression_mode == 2) && (layout == 1))) {
      logerrputs("Error: Invalid .bgen header.\n");
      goto OxBgenToPgen_ret_MALFORMED_INPUT;
    }
    logprintf("--bgen: %u variant%s declared in header, format v1.%c.\n", header_variant_ct, (header_variant_ct == 1)? "" : "s", (layout == 1)? '1' : ((compression_mode == 2)? '3' : '2'));
    uint32_t raw_variant_ct = header_variant_ct;
    if (samplename[0]) {
      uint32_t sfile_sample_ct;
      reterr = OxSampleToPsam(samplename, const_fid, ox_missing_code, missing_catname, misc_flags, import_flags, psam_01, id_delim, outname, outname_end, &sfile_sample_ct);
      if (unlikely(reterr)) {
        goto OxBgenToPgen_ret_1;
      }
      if (unlikely(sfile_sample_ct != sample_ct)) {
        logerrprintf("Error: .sample file has %u sample%s, while .bgen file has %u.\n", sfile_sample_ct, (sfile_sample_ct == 1)? "" : "s", sample_ct);
        goto OxBgenToPgen_ret_INCONSISTENT_INPUT;
      }
      if (header_flags >> 31) {
        uint32_t sample_id_block_byte_ct;
        uint32_t sample_id_block_entry_ct;
        if (unlikely((!fread_unlocked(&sample_id_block_byte_ct, 4, 1, bgenfile)) ||
                     (!fread_unlocked(&sample_id_block_entry_ct, 4, 1, bgenfile)))) {
          goto OxBgenToPgen_ret_READ_FAIL;
        }
        if (unlikely((S_CAST(uint64_t, sample_id_block_byte_ct) + initial_uints[1] > initial_uints[0]) ||
                     (sample_id_block_entry_ct != sample_ct))) {
          logerrputs("Error: Invalid .bgen header.\n");
          goto OxBgenToPgen_ret_MALFORMED_INPUT;
        }
      }
    } else {
      if (unlikely(!(header_flags >> 31))) {
        logerrputs("Error: .bgen file does not contain sample IDs, and no .sample file was\nspecified.\n");
        goto OxBgenToPgen_ret_INCONSISTENT_INPUT;
      }
      // possible todo: optionally error out if sample IDs aren't consistent
      // between .bgen and .sample, using ImportIidFromSampleId()

      // see VcfSampleLine()
      uint32_t sample_id_block_byte_ct;
      uint32_t sample_id_block_entry_ct;
      if (unlikely((!fread_unlocked(&sample_id_block_byte_ct, 4, 1, bgenfile)) ||
                   (!fread_unlocked(&sample_id_block_entry_ct, 4, 1, bgenfile)))) {
        goto OxBgenToPgen_ret_READ_FAIL;
      }
      if (unlikely((sample_id_block_byte_ct < 8) ||
                   (S_CAST(uint64_t, sample_id_block_byte_ct) + initial_uints[1] > initial_uints[0]) ||
                   (sample_id_block_entry_ct != sample_ct))) {
        logerrputs("Error: Invalid .bgen header.\n");
        goto OxBgenToPgen_ret_MALFORMED_INPUT;
      }
      sample_id_block_byte_ct -= 8;
      char* sample_id_block_main;
      if (unlikely(bigstack_alloc_c(sample_id_block_byte_ct, &sample_id_block_main))) {
        goto OxBgenToPgen_ret_NOMEM;
      }
      char* sample_id_block_end = &(sample_id_block_main[sample_id_block_byte_ct]);
      if (unlikely(fread_checked(sample_id_block_main, sample_id_block_byte_ct, bgenfile))) {
        goto OxBgenToPgen_ret_READ_FAIL;
      }

      // Always check if any tab/eoln characters are present, and error out if
      // so.
      // if id_delim != ' ', also check if spaces are present; if so, replace
      // with --idspace-to character or error out
      char* sample_id_block_iter = sample_id_block_main;
      for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        uint16_t input_id_slen;
        memcpy_k(&input_id_slen, sample_id_block_iter, sizeof(int16_t));

        // need to check this to avoid read-past-the-end indeterminate
        // behavior
        if (unlikely(S_CAST(uintptr_t, sample_id_block_end - sample_id_block_iter) < input_id_slen + (2 * k1LU))) {
          logerrputs("Error: Invalid .bgen header.\n");
          goto OxBgenToPgen_ret_MALFORMED_INPUT;
        }
        char* sample_id_iter = &(sample_id_block_iter[2]);
        char* sample_id_end = &(sample_id_iter[input_id_slen]);
        uint32_t char_code_min = 32 + (id_delim != ' ');
        for (; sample_id_iter != sample_id_end; ++sample_id_iter) {
          const uint32_t char_code = ctou32(*sample_id_iter);
          if (char_code < char_code_min) {
            if (unlikely(char_code < 32)) {
              logerrputs("Error: .bgen sample ID contains tabs, newlines, and/or nonprinting characters.\n");
              goto OxBgenToPgen_ret_MALFORMED_INPUT;
            }
            if (unlikely(!idspace_to)) {
              logerrputs("Error: .bgen sample ID contains space(s).  Use --idspace-to to convert them to\nanother character, or \"--id-delim ' '\" to interpret the spaces as FID/IID or\nIID/SID delimiters.\n");
              goto OxBgenToPgen_ret_INCONSISTENT_INPUT;
            }
            *sample_id_iter = idspace_to;
          }
        }
        sample_id_block_iter = sample_id_end;
      }
      snprintf(outname_end, kMaxOutfnameExtBlen, ".psam");
      if (unlikely(fopen_checked(outname, FOPEN_WB, &psamfile))) {
        goto OxBgenToPgen_ret_OPEN_FAIL;
      }
      ImportSampleIdContext isic;
      InitImportSampleIdContext(const_fid, import_flags, id_delim, &isic);
      uint32_t write_fid = 0;
      uint32_t write_sid = 0;
      if (id_delim) {
        sample_id_block_iter = sample_id_block_main;
        const uint32_t iid_sid = (misc_flags / kfMiscIidSid) & 1;
        for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
          uint16_t input_id_slen;
          memcpy_k(&input_id_slen, sample_id_block_iter, sizeof(int16_t));
          char* sample_id_iter = &(sample_id_block_iter[2]);
          // previously verified that this is in-bounds
          char* sample_id_end = &(sample_id_iter[input_id_slen]);
          char* first_delim = S_CAST(char*, memchr(sample_id_iter, ctou32(id_delim), sample_id_end - sample_id_iter));
          if (unlikely(!first_delim)) {
            snprintf(g_logbuf, kLogbufSize, "Error: No '%c' in sample ID.\n", id_delim);
            goto OxBgenToPgen_ret_INCONSISTENT_INPUT_2;
          }
          const uint32_t first_slen = first_delim - sample_id_iter;
          char* second_part_start = &(first_delim[1]);
          char* maybe_second_part_end = S_CAST(char*, memchr(sample_id_iter, ctou32(id_delim), sample_id_end - second_part_start));
          if (maybe_second_part_end == nullptr) {
            if (!iid_sid) {
              write_fid |= (first_slen != 1) || (sample_id_iter[0] != '0');
            } else {
              const uint32_t sid_slen = sample_id_end - second_part_start;
              write_sid |= (sid_slen != 1) || (second_part_start[0] != '0');
            }
          } else {
            write_fid |= (first_slen != 1) || (sample_id_iter[0] != '0');
            char* sid_start = &(maybe_second_part_end[1]);
            const uint32_t sid_slen = sample_id_end - sid_start;
            if (unlikely(memchr(sid_start, ctou32(id_delim), sid_slen))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Too many instances of --id-delim argument '%c' in sample ID.\n", id_delim);
              goto OxBgenToPgen_ret_INCONSISTENT_INPUT_2;
            }
            write_sid |= (sid_slen != 1) || (sid_start[0] != '0');
          }
          sample_id_block_iter = sample_id_end;
        }
        if (!iid_sid) {
          isic.fid_delim_mode = write_fid? kImportFidDelimModeAlwaysCopy : kImportFidDelimModeAlwaysOmit;
          isic.sid_delim_mode = write_sid? kImportSidDelimModeCopyOr0 : kImportSidDelimModeNonexistOrOmit;
        } else {
          isic.fid_delim_mode = write_fid? kImportFidDelimModeCopyOr0 : kImportFidDelimModeOmitWhenPresent;
          isic.sid_delim_mode = write_sid? kImportSidDelimModeAlwaysCopy : kImportSidDelimModeNonexistOrOmit;
        }
      } else if (isic.const_fid || isic.double_id) {
        write_fid = 1;
      }
      char* textbuf = g_textbuf;
      char* write_iter = textbuf;
      *write_iter++ = '#';
      if (write_fid) {
        write_iter = strcpya_k(write_iter, "FID\t");
      }
      write_iter = strcpya_k(write_iter, "IID");
      if (write_sid) {
        write_iter = strcpya_k(write_iter, "\tSID");
      }
      write_iter = strcpya_k(write_iter, "\tSEX");
      AppendBinaryEoln(&write_iter);
      char* textbuf_flush = &(textbuf[kMaxMediumLine]);
      sample_id_block_iter = sample_id_block_main;
      for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        uint16_t input_id_slen;
        memcpy_k(&input_id_slen, sample_id_block_iter, sizeof(int16_t));
        char* sample_id_start = &(sample_id_block_iter[2]);
        char* sample_id_end = &(sample_id_start[input_id_slen]);
        reterr = ImportSampleId(sample_id_start, sample_id_end, &isic, &write_iter);
        if (unlikely(reterr)) {
          goto OxBgenToPgen_ret_1;
        }
        // SEX
        write_iter = strcpya_k(write_iter, "\tNA");
        AppendBinaryEoln(&write_iter);
        if (write_iter >= textbuf_flush) {
          if (unlikely(fwrite_checked(textbuf, write_iter - textbuf, psamfile))) {
            goto OxBgenToPgen_ret_WRITE_FAIL;
          }
          write_iter = textbuf;
        }
        sample_id_block_iter = sample_id_end;
      }
      if (unlikely(sample_id_block_iter != &(sample_id_block_main[sample_id_block_byte_ct]))) {
        logerrputs("Error: Invalid .bgen header.\n");
        goto OxBgenToPgen_ret_MALFORMED_INPUT;
      }
      if (write_iter != textbuf) {
        if (unlikely(fwrite_checked(textbuf, write_iter - textbuf, psamfile))) {
          goto OxBgenToPgen_ret_WRITE_FAIL;
        }
      }
      BigstackReset(sample_id_block_main);
      if (unlikely(fclose_null(&psamfile))) {
        goto OxBgenToPgen_ret_WRITE_FAIL;
      }
      logprintfww("--bgen: %u sample ID%s written to %s .\n", sample_ct, (sample_ct == 1)? "" : "s", outname);
    }
    if (unlikely(fseeko(bgenfile, initial_uints[0] + 4, SEEK_SET))) {
      goto OxBgenToPgen_ret_READ_FAIL;
    }
    const uint32_t prohibit_extra_chr = (misc_flags / kfMiscProhibitExtraChr) & 1;
    FinalizeChrset(load_filter_log_import_flags, cip);
    const uint32_t autosome_ct_p1 = cip->autosome_ct + 1;
    uint32_t chr_filter_exists = (PopcountBitRange(cip->chr_mask, 0, autosome_ct_p1) != autosome_ct_p1) || ((!prohibit_extra_chr) && (cip->is_include_stack || cip->incl_excl_name_stack));
    if (!chr_filter_exists) {
      for (uint32_t xymt_idx = 0; xymt_idx != kChrOffsetCt; ++xymt_idx) {
        if (!IsI32Neg(cip->xymt_codes[xymt_idx])) {
          if (!IsSet(cip->chr_mask, autosome_ct_p1 + xymt_idx)) {
            chr_filter_exists = 1;
            break;
          }
        }
      }
    }

    if (unlikely(BIGSTACK_ALLOC_X(struct libdeflate_decompressor*, max_thread_ct, &common.libdeflate_decompressors))) {
      goto OxBgenToPgen_ret_NOMEM;
    }
    ZeroPtrArr(max_thread_ct, common.libdeflate_decompressors);
    // make libdeflate_alloc_decompressor() calls later, when we know how many
    // decompressor threads we're using

    uint32_t cur_chr_code = 0;
    if (ox_single_chr_str) {
      reterr = InitOxfordSingleChr(ox_single_chr_str, "--bgen file", nullptr, nullptr, &cur_chr_code, cip);
      if (unlikely(reterr)) {
        goto OxBgenToPgen_ret_1;
      }
    }
    snprintf(outname_end, kMaxOutfnameExtBlen, ".pvar");
    const uint32_t output_zst = (import_flags / kfImportKeepAutoconvVzs) & 1;
    if (output_zst) {
      snprintf(&(outname_end[5]), kMaxOutfnameExtBlen - 5, ".zst");
    }
    const uintptr_t overflow_buf_size = 2 * kCompressStreamBlock + kCacheline;
    reterr = InitCstreamAlloc(outname, 0, output_zst, 1, overflow_buf_size, &pvar_css, &pvar_cswritep);
    if (unlikely(reterr)) {
      goto OxBgenToPgen_ret_1;
    }
    if (cip->chrset_source) {
      AppendChrsetLine(cip, &pvar_cswritep);
    }
    pvar_cswritep = strcpya_k(pvar_cswritep, "#CHROM\tPOS\tID\tREF\tALT" EOLN_STR);

    const uint32_t snpid_chr = (oxford_import_flags & kfOxfordImportBgenSnpIdChr);

    // true for both provisional-reference and real-reference second
    const uint32_t prov_ref_allele_second = !(oxford_import_flags & kfOxfordImportRefFirst);

    const uint32_t allow_overstated_variant_ct = (import_flags / kfImportLaxBgen) & 1;

    if (hard_call_thresh == UINT32_MAX) {
      hard_call_thresh = kDosageMid / 10;
    }
    const uint32_t sample_ctaw2 = NypCtToAlignedWordCt(sample_ct);
    const uint32_t sample_ctaw = BitCtToAlignedWordCt(sample_ct);
    uint32_t dosage_exists = 0;
    common.sample_ct = sample_ct;
    common.dosage_erase_halfdist = kDosage4th - dosage_erase_thresh;
    common.compression_mode = compression_mode;
    uint32_t par_warn_code = UINT32_MAX;
    uint32_t print_splitpar_warning = 0;
    if (!(import_flags & kfImportLaxChrX)) {
      const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
      if ((!IsI32Neg(x_code)) && ((!sex_info_avail) || (IsHumanChrset(cip) && (!is_splitpar)))) {
        par_warn_code = x_code;
      }
    }
    uint32_t variant_ct = 0;
    // temporary kludge
    uint32_t multiallelic_tmp_skip_ct = 0;
    if (layout == 1) {
      // v1.1
      // this block belongs in its own function...
      Bgen11DosageScanCtx scan_ctx;
      scan_ctx.common = &common;
      scan_ctx.dosage_exists = 0;
      scan_ctx.reterr = kPglRetSuccess;
      uintptr_t loadbuf_size = RoundDownPow2(bigstack_left() / 4, kCacheline);
#ifdef __LP64__
      if (loadbuf_size > kMaxLongLine) {
        loadbuf_size = kMaxLongLine;
      }
#endif
      // must have enough space for chromosome and variant IDs
      if (unlikely(loadbuf_size < 2 * 65536)) {
        goto OxBgenToPgen_ret_NOMEM;
      }
      unsigned char* loadbuf = S_CAST(unsigned char*, bigstack_alloc_raw(loadbuf_size));
      scan_ctx.import_dosage_certainty_int = 1 + S_CAST(int32_t, import_dosage_certainty * 32768);
      uintptr_t bgen_geno_max_byte_ct = 6LU * sample_ct;
      if (compression_mode) {
        bgen_geno_max_byte_ct = libdeflate_deflate_compress_bound(nullptr, bgen_geno_max_byte_ct);
      }
      if (unlikely(bgen_geno_max_byte_ct > UINT32_MAX)) {
        logerrputs("Error: Too many samples for .bgen format.\n");
        goto OxBgenToPgen_ret_MALFORMED_INPUT;
      }
      bgen_geno_max_byte_ct += compression_mode * 4;
      // thread-count-independent:
      //   (everything after "2 *" rounded up to cacheline)
      //   compressed_geno_bufs: 2 * bgen_geno_max_byte_ct * main_block_size
      //   common.compressed_geno_starts: 2 * sizeof(intptr_t) *
      //                                  main_block_size
      //   ctx.write_genovecs: 2 * sample_ctaw2 * sizeof(intptr_t) *
      //                       main_block_size
      //   ctx.write_dosage_cts: 2 * sizeof(int32_t) * main_block_size
      //   ctx.write_dosage_presents: 2 * sample_ctaw * sizeof(intptr_t) *
      //                              main_block_size
      //   ctx.write_dosage_mains (main bottleneck): 2 * sample_ct *
      //                                           sizeof(Dosage)
      // additional requirement per thread:
      //   ctx.bgen_geno_bufs: sample_ct * 3 * sizeof(int16_t)

      uint32_t calc_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
      if ((!compression_mode) && (calc_thread_ct > 2)) {
        // computation doesn't seem to saturate when decompression is involved
        // (well, it didn't before libdeflate, anyway.  recheck this.)
        calc_thread_ct = 2;
      }
      if (calc_thread_ct > header_variant_ct) {
        calc_thread_ct = header_variant_ct;
      }
      if (unlikely(bigstack_alloc_u16p(calc_thread_ct, &scan_ctx.bgen_geno_bufs))) {
        goto OxBgenToPgen_ret_NOMEM;
      }
      const uint32_t sample_ct_x3 = sample_ct * 3;
      for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
        if (unlikely(bigstack_alloc_u16(sample_ct_x3, &(scan_ctx.bgen_geno_bufs[tidx])))) {
          goto OxBgenToPgen_ret_NOMEM;
        }
      }
      uintptr_t cachelines_avail_m12 = bigstack_left() / kCacheline;
      // reserve 1/8 of remaining memory for writer
      cachelines_avail_m12 -= cachelines_avail_m12 / 8;
      if (unlikely(cachelines_avail_m12 < 12)) {
        goto OxBgenToPgen_ret_NOMEM;
      }
      // we're making 12 allocations; be pessimistic re: rounding
      cachelines_avail_m12 -= 12;
      const uintptr_t bytes_req_per_in_block_variant = 2 * (bgen_geno_max_byte_ct + sizeof(intptr_t) + sample_ctaw2 * sizeof(intptr_t) + sizeof(int32_t) + sample_ctaw * sizeof(intptr_t) + sample_ct * sizeof(Dosage));
      uintptr_t main_block_size = (cachelines_avail_m12 * kCacheline) / bytes_req_per_in_block_variant;
      if (main_block_size > 65536) {
        main_block_size = 65536;
      } else if (unlikely(main_block_size < 8)) {
        // this threshold is arbitrary
        goto OxBgenToPgen_ret_NOMEM;
      }
      if (calc_thread_ct > main_block_size / 8) {
        calc_thread_ct = main_block_size / 8;
      }
      if (unlikely(SetThreadCt(calc_thread_ct, &tg))) {
        goto OxBgenToPgen_ret_NOMEM;
      }
      for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
        common.libdeflate_decompressors[tidx] = libdeflate_alloc_decompressor();
        if (unlikely(!common.libdeflate_decompressors[tidx])) {
          goto OxBgenToPgen_ret_NOMEM;
        }
      }
      Bgen11GenoToPgenCtx ctx;
      ctx.common = &common;
      ctx.bgen_geno_bufs = scan_ctx.bgen_geno_bufs;
      ctx.hard_call_halfdist = kDosage4th - hard_call_thresh;
      ctx.import_dosage_certainty_int = scan_ctx.import_dosage_certainty_int;
      ctx.prov_ref_allele_second = prov_ref_allele_second;
      ctx.reterr = kPglRetSuccess;
      unsigned char* compressed_geno_bufs[2];
      if (unlikely(bigstack_alloc_uc(bgen_geno_max_byte_ct * main_block_size, &(compressed_geno_bufs[0])) ||
                   bigstack_alloc_uc(bgen_geno_max_byte_ct * main_block_size, &(compressed_geno_bufs[1])) ||
                   bigstack_alloc_ucp(main_block_size, &(common.compressed_geno_starts[0])) ||
                   bigstack_alloc_ucp(main_block_size, &(common.compressed_geno_starts[1])) ||
                   bigstack_alloc_w(sample_ctaw2 * main_block_size, &(ctx.write_genovecs[0])) ||
                   bigstack_alloc_w(sample_ctaw2 * main_block_size, &(ctx.write_genovecs[1])) ||
                   bigstack_alloc_u32(main_block_size, &(ctx.write_dosage_cts[0])) ||
                   bigstack_alloc_u32(main_block_size, &(ctx.write_dosage_cts[1])) ||
                   bigstack_alloc_w(sample_ctaw * main_block_size, &(ctx.write_dosage_presents[0])) ||
                   bigstack_alloc_w(sample_ctaw * main_block_size, &(ctx.write_dosage_presents[1])) ||
                   bigstack_alloc_dosage(sample_ct * main_block_size, &(ctx.write_dosage_mains[0])) ||
                   bigstack_alloc_dosage(sample_ct * main_block_size, &(ctx.write_dosage_mains[1])))) {
        // this should be impossible
        assert(0);
        goto OxBgenToPgen_ret_NOMEM;
      }
      SetThreadFuncAndData(Bgen11DosageScanThread, &scan_ctx, &tg);

      // likely cases are (i) non-hardcall near top of the file, and (ii) no
      // non-hardcalls at all.  to handle the first case efficiently, we want
      // the first blocks to be small so we bail quickly; to handle the second
      // case efficiently, we want large blocks on average.  so we start with
      // a minimal block size and then repeatedly double.
      uint32_t block_vidx = 0;
      uint32_t cur_block_size = calc_thread_ct;
      uint32_t parity = 0;
      uintptr_t compressed_block_byte_ct = 6LU * sample_ct;
      unsigned char** compressed_geno_starts = common.compressed_geno_starts[0];
      unsigned char* bgen_geno_iter = compressed_geno_bufs[0];
      uint32_t skip = 0;
      for (uint32_t variant_uidx = 0; variant_uidx != header_variant_ct; ) {
        uint32_t uii;
        {
          const uintptr_t bytes_read = fread_unlocked(&uii, 1, 4, bgenfile);
          if (bytes_read != 4) {
            // don't know of a bgen-1.1 use case, but may as well make
            // --lax-bgen-import have a consistent effect.
            if (likely(allow_overstated_variant_ct && (!bytes_read))) {
              raw_variant_ct = variant_uidx;
              putc_unlocked('\n', stdout);
              logprintf("--lax-bgen-import: .bgen file actually contains %u variant%s.\n", raw_variant_ct, (raw_variant_ct == 1)? "" : "s");
              break;
            }
            goto OxBgenToPgen_ret_READ_FAIL;
          }
        }
        if (unlikely(uii != sample_ct)) {
          putc_unlocked('\n', stdout);
          logerrputs("Error: Unexpected number of samples specified in SNP block header.\n");
          goto OxBgenToPgen_ret_MALFORMED_INPUT;
        }
        uint16_t snpid_slen;
        if (unlikely(!fread_unlocked(&snpid_slen, 2, 1, bgenfile))) {
          goto OxBgenToPgen_ret_READ_FAIL;
        }
        if (!snpid_chr) {
          if (unlikely(fseeko(bgenfile, snpid_slen, SEEK_CUR))) {
            goto OxBgenToPgen_ret_READ_FAIL;
          }
        } else {
          if (unlikely(!snpid_slen)) {
            putc_unlocked('\n', stdout);
            logerrputs("Error: Length-0 SNP ID in .bgen file.\n");
            goto OxBgenToPgen_ret_INCONSISTENT_INPUT;
          }
          if (unlikely(!fread_unlocked(loadbuf, snpid_slen, 1, bgenfile))) {
            goto OxBgenToPgen_ret_READ_FAIL;
          }
          // loadbuf[snpid_slen] = '\0';
        }
        uint16_t rsid_slen;
        if (unlikely(!fread_unlocked(&rsid_slen, 2, 1, bgenfile))) {
          goto OxBgenToPgen_ret_READ_FAIL;
        }
        if (unlikely(fseeko(bgenfile, rsid_slen, SEEK_CUR))) {
          goto OxBgenToPgen_ret_READ_FAIL;
        }
        uint16_t chr_name_slen;
        if (unlikely(!fread_unlocked(&chr_name_slen, 2, 1, bgenfile))) {
          goto OxBgenToPgen_ret_READ_FAIL;
        }
        if (ox_single_chr_str) {
          if (unlikely(fseeko(bgenfile, chr_name_slen, SEEK_CUR))) {
            goto OxBgenToPgen_ret_READ_FAIL;
          }
        } else {
          if (!snpid_chr) {
            if (unlikely(!chr_name_slen)) {
              putc_unlocked('\n', stdout);
              logerrputs("Error: Length-0 chromosome ID in .bgen file.\n");
              goto OxBgenToPgen_ret_INCONSISTENT_INPUT;
            }
            if (unlikely(!fread_unlocked(loadbuf, chr_name_slen, 1, bgenfile))) {
              goto OxBgenToPgen_ret_READ_FAIL;
            }
            if (strequal_k(R_CAST(char*, loadbuf), "NA", chr_name_slen)) {
              loadbuf[0] = '0';
              chr_name_slen = 1;
            }
          } else {
            chr_name_slen = snpid_slen;
          }
          loadbuf[chr_name_slen] = '\0';
          reterr = GetOrAddChrCode(R_CAST(char*, loadbuf), "--bgen file", 0, chr_name_slen, prohibit_extra_chr, cip, &cur_chr_code);
          if (unlikely(reterr)) {
            goto OxBgenToPgen_ret_1;
          }
          skip = !IsSet(cip->chr_mask, cur_chr_code);
        }

        uint32_t cur_bp;  // ignore in this pass
        if (unlikely(!fread_unlocked(&cur_bp, 4, 1, bgenfile))) {
          goto OxBgenToPgen_ret_READ_FAIL;
        }

        // allele count always 2 and not stored when layout=1
        for (uint32_t allele_idx = 0; allele_idx != 2; ++allele_idx) {
          uint32_t allele_slen;
          if (unlikely(!fread_unlocked(&allele_slen, 4, 1, bgenfile))) {
            goto OxBgenToPgen_ret_READ_FAIL;
          }
          if (unlikely(fseeko(bgenfile, allele_slen, SEEK_CUR))) {
            goto OxBgenToPgen_ret_READ_FAIL;
          }
        }

        if (compression_mode) {
#ifdef __LP64__
          compressed_block_byte_ct = 0;
#endif
          if (unlikely(!fread_unlocked(&compressed_block_byte_ct, 4, 1, bgenfile))) {
            goto OxBgenToPgen_ret_READ_FAIL;
          }
        }
        ++variant_uidx;
        if (!(variant_uidx % 1000)) {
          printf("\r--bgen: %uk variants scanned.", variant_uidx / 1000);
          fflush(stdout);
        }
        if (dosage_exists || skip) {
          if (unlikely(fseeko(bgenfile, compressed_block_byte_ct, SEEK_CUR))) {
            goto OxBgenToPgen_ret_READ_FAIL;
          }
          // bugfix (25 Jun 2017): block_vidx should be left unchanged here
          variant_ct += 1 - skip;
          continue;
        }
        compressed_geno_starts[block_vidx] = bgen_geno_iter;
        if (compression_mode) {
          memcpy(bgen_geno_iter, &compressed_block_byte_ct, 4);
          bgen_geno_iter = &(bgen_geno_iter[4]);
        }
        if (unlikely(fread_checked(bgen_geno_iter, compressed_block_byte_ct, bgenfile))) {
          goto OxBgenToPgen_ret_READ_FAIL;
        }
        bgen_geno_iter = &(bgen_geno_iter[compressed_block_byte_ct]);
        ++block_vidx;
        if (block_vidx == cur_block_size) {
          parity = 1 - parity;
          if (ThreadsAreActive(&tg)) {
            // process *previous* block results
            JoinThreads(&tg);
            reterr = scan_ctx.reterr;
            if (unlikely(reterr)) {
              goto OxBgenToPgen_ret_bgen11_thread_fail;
            }
            dosage_exists = scan_ctx.dosage_exists;
            if (dosage_exists) {
              // don't need to scan for any more dosages
              StopThreads(&tg);
              if ((!chr_filter_exists) && (!allow_overstated_variant_ct)) {
                break;
              }
              continue;
            }
          }
          common.cur_block_size = cur_block_size;
          if (unlikely(SpawnThreads(&tg))) {
            goto OxBgenToPgen_ret_THREAD_CREATE_FAIL;
          }
          compressed_geno_starts = common.compressed_geno_starts[parity];
          bgen_geno_iter = compressed_geno_bufs[parity];
          block_vidx = 0;
          variant_ct += cur_block_size;
          if (cur_block_size < main_block_size) {
            cur_block_size *= 2;
            if (cur_block_size > main_block_size) {
              cur_block_size = main_block_size;
            }
          }
        }
      }

      if (!chr_filter_exists) {
        variant_ct = raw_variant_ct;
      } else {
        variant_ct += block_vidx;
        if (variant_ct < calc_thread_ct) {
          if (unlikely(!variant_ct)) {
            putc_unlocked('\n', stdout);
            char* write_iter = strcpya_k(g_logbuf, "Error: All ");
            write_iter = u32toa(raw_variant_ct, write_iter);
            write_iter = strcpya_k(write_iter, " variant");
            if (raw_variant_ct != 1) {
              *write_iter++ = 's';
            }
            write_iter = strcpya_k(write_iter, " in .bgen file excluded by ");
            AppendLoadFilterFlagnames(load_filter_log_import_flags, &write_iter);
            strcpy_k(write_iter, ".\n");
            goto OxBgenToPgen_ret_INCONSISTENT_INPUT_WW;
          }
          // bugfix (7 Oct 2017): with fewer variants than threads, need to
          // force initial launch here
          common.cur_block_size = variant_ct;
          if (unlikely(SpawnThreads(&tg))) {
            goto OxBgenToPgen_ret_THREAD_CREATE_FAIL;
          }
          block_vidx = 0;
        }
      }
      if (ThreadsAreActive(&tg)) {
        JoinThreads(&tg);
        reterr = scan_ctx.reterr;
        if (unlikely(reterr)) {
          goto OxBgenToPgen_ret_bgen11_thread_fail;
        }
        if (block_vidx && (!scan_ctx.dosage_exists)) {
          common.cur_block_size = block_vidx;
        } else {
          common.cur_block_size = 0;
        }
        DeclareLastThreadBlock(&tg);
        if (unlikely(SpawnThreads(&tg))) {
          goto OxBgenToPgen_ret_THREAD_CREATE_FAIL;
        }
        JoinThreads(&tg);
        reterr = scan_ctx.reterr;
        if (unlikely(reterr)) {
          goto OxBgenToPgen_ret_bgen11_thread_fail;
        }
        dosage_exists = scan_ctx.dosage_exists;
      }

      if (unlikely(fseeko(bgenfile, initial_uints[0] + 4, SEEK_SET))) {
        goto OxBgenToPgen_ret_READ_FAIL;
      }
      snprintf(outname_end, kMaxOutfnameExtBlen, ".pgen");
      uintptr_t spgw_alloc_cacheline_ct;
      uint32_t max_vrec_len;
      reterr = SpgwInitPhase1(outname, nullptr, nullptr, variant_ct, sample_ct, 0, kPgenWriteBackwardSeek, dosage_exists? kfPgenGlobalDosagePresent : kfPgenGlobal0, (oxford_import_flags & kfOxfordImportRefUnknown)? 2 : 1, &spgw, &spgw_alloc_cacheline_ct, &max_vrec_len);
      if (unlikely(reterr)) {
        if (reterr == kPglRetOpenFail) {
          putc_unlocked('\n', stdout);
          logerrprintfww(kErrprintfFopen, outname, strerror(errno));
        }
        goto OxBgenToPgen_ret_1;
      }
      unsigned char* spgw_alloc;
      if (unlikely(bigstack_alloc_uc(spgw_alloc_cacheline_ct * kCacheline, &spgw_alloc))) {
        goto OxBgenToPgen_ret_NOMEM;
      }
      SpgwInitPhase2(max_vrec_len, &spgw, spgw_alloc);

      // Main workflow:
      // 1. Set n=0, load genotype data for first main_block_size variants
      //    while writing .pvar
      //
      // 2. Spawn threads processing batch n genotype data
      // 3. If n>0, write results for block (n-1)
      // 4. Increment n by 1
      // 5. Load/write-.pvar for batch (n+1) unless eof
      // 6. Join threads
      // 7. Goto step 2 unless eof
      //
      // 8. Write results for last block
      //
      // (May be better to change this to use one output buffer instead of 2.)
      uint32_t prev_block_write_ct = 0;
      parity = 0;
      SetThreadFuncAndData(Bgen11GenoToPgenThread, &ctx, &tg);
      for (uint32_t vidx_start = 0; ; ) {
        uint32_t cur_block_write_ct = 0;
        if (!IsLastBlock(&tg)) {
          cur_block_write_ct = MINV(variant_ct - vidx_start, main_block_size);
          compressed_geno_starts = common.compressed_geno_starts[parity];
          bgen_geno_iter = compressed_geno_bufs[parity];
          for (block_vidx = 0; block_vidx != cur_block_write_ct; ) {
            uint32_t uii;
            if (unlikely(!fread_unlocked(&uii, 4, 1, bgenfile))) {
              goto OxBgenToPgen_ret_READ_FAIL;
            }
            if (unlikely(uii != sample_ct)) {
              putc_unlocked('\n', stdout);
              logerrputs("Error: Unexpected number of samples specified in SNP block header.\n");
              goto OxBgenToPgen_ret_MALFORMED_INPUT;
            }
            uint16_t snpid_slen;
            if (unlikely(!fread_unlocked(&snpid_slen, 2, 1, bgenfile))) {
              goto OxBgenToPgen_ret_READ_FAIL;
            }
            char* rsid_start = R_CAST(char*, loadbuf);
            if (!snpid_chr) {
              if (unlikely(fseeko(bgenfile, snpid_slen, SEEK_CUR))) {
                goto OxBgenToPgen_ret_READ_FAIL;
              }
            } else {
              if (unlikely(!snpid_slen)) {
                putc_unlocked('\n', stdout);
                logerrputs("Error: Length-0 SNP ID in .bgen file.\n");
                goto OxBgenToPgen_ret_INCONSISTENT_INPUT;
              }
              if (unlikely(!fread_unlocked(loadbuf, snpid_slen, 1, bgenfile))) {
                goto OxBgenToPgen_ret_READ_FAIL;
              }
              // loadbuf[snpid_slen] = '\0';
              rsid_start = R_CAST(char*, &(loadbuf[snpid_slen + 1]));
            }
            uint16_t rsid_slen;
            if (unlikely(!fread_unlocked(&rsid_slen, 2, 1, bgenfile))) {
              goto OxBgenToPgen_ret_READ_FAIL;
            }
            if (unlikely(!rsid_slen)) {
              putc_unlocked('\n', stdout);
              logerrputs("Error: Length-0 rsID in .bgen file.\n");
              goto OxBgenToPgen_ret_MALFORMED_INPUT;
            }
            if (unlikely(!fread_unlocked(rsid_start, rsid_slen, 1, bgenfile))) {
              goto OxBgenToPgen_ret_READ_FAIL;
            }
            char* loadbuf_iter = &(rsid_start[rsid_slen]);
            char* chr_name_start = loadbuf_iter;
            uint16_t chr_name_slen;
            if (unlikely(!fread_unlocked(&chr_name_slen, 2, 1, bgenfile))) {
              goto OxBgenToPgen_ret_READ_FAIL;
            }
            if (ox_single_chr_str) {
              if (unlikely(fseeko(bgenfile, chr_name_slen, SEEK_CUR))) {
                goto OxBgenToPgen_ret_READ_FAIL;
              }
            } else {
              if (!snpid_chr) {
                if (unlikely(!chr_name_slen)) {
                  putc_unlocked('\n', stdout);
                  logerrputs("Error: Length-0 chromosome ID in .bgen file.\n");
                  goto OxBgenToPgen_ret_INCONSISTENT_INPUT;
                }
                if (unlikely(!fread_unlocked(chr_name_start, chr_name_slen, 1, bgenfile))) {
                  goto OxBgenToPgen_ret_READ_FAIL;
                }
                if (strequal_k(chr_name_start, "NA", chr_name_slen)) {
                  chr_name_start[0] = '0';
                  chr_name_slen = 1;
                }
              } else {
                if (unlikely(fseeko(bgenfile, chr_name_slen, SEEK_CUR))) {
                  goto OxBgenToPgen_ret_READ_FAIL;
                }
                chr_name_start = R_CAST(char*, loadbuf);
                chr_name_slen = snpid_slen;
              }
              chr_name_start[chr_name_slen] = '\0';
              cur_chr_code = GetChrCode(chr_name_start, cip, chr_name_slen);
              skip = !IsSet(cip->chr_mask, cur_chr_code);
            }

            uint32_t cur_bp;
            uint32_t a1_slen;
            if (unlikely((!fread_unlocked(&cur_bp, 4, 1, bgenfile)) ||
                         (!fread_unlocked(&a1_slen, 4, 1, bgenfile)))) {
              goto OxBgenToPgen_ret_READ_FAIL;
            }
            if (skip) {
              uint32_t a2_slen;
              if (unlikely(fseeko(bgenfile, a1_slen, SEEK_CUR) ||
                           (!fread_unlocked(&a2_slen, 4, 1, bgenfile)) ||
                           fseeko(bgenfile, a2_slen, SEEK_CUR))) {
                goto OxBgenToPgen_ret_READ_FAIL;
              }
              if (compression_mode) {
#ifdef __LP64__
                compressed_block_byte_ct = 0;
#endif
                if (unlikely(!fread_unlocked(&compressed_block_byte_ct, 4, 1, bgenfile))) {
                  goto OxBgenToPgen_ret_READ_FAIL;
                }
              }
              if (unlikely(fseeko(bgenfile, compressed_block_byte_ct, SEEK_CUR))) {
                goto OxBgenToPgen_ret_READ_FAIL;
              }
              continue;
            }
            char* a1_ptr = loadbuf_iter;
            if (unlikely(!a1_slen)) {
              putc_unlocked('\n', stdout);
              logerrputs("Error: Empty allele code in .bgen file.\n");
              goto OxBgenToPgen_ret_MALFORMED_INPUT;
            }
            // TODO: enforce a consistent, configurable limit across the
            // program (2^26 default?)
            if (unlikely(a1_slen > 1000000000)) {
              putc_unlocked('\n', stdout);
              logerrputs("Error: Allele code in .bgen file has more than 1 billion characters.\n");
              goto OxBgenToPgen_ret_MALFORMED_INPUT;
            }
            if (unlikely(a1_slen + S_CAST(uintptr_t, a1_ptr - R_CAST(char*, loadbuf)) > loadbuf_size)) {
              goto OxBgenToPgen_ret_NOMEM;
            }
            if (unlikely(!fread_unlocked(a1_ptr, a1_slen, 1, bgenfile))) {
              goto OxBgenToPgen_ret_READ_FAIL;
            }
            char* a2_ptr = &(a1_ptr[a1_slen]);
            uint32_t a2_slen;
            if (unlikely(!fread_unlocked(&a2_slen, 4, 1, bgenfile))) {
              goto OxBgenToPgen_ret_READ_FAIL;
            }
            if (unlikely(!a2_slen)) {
              putc_unlocked('\n', stdout);
              logerrputs("Error: Empty allele code in .bgen file.\n");
              goto OxBgenToPgen_ret_MALFORMED_INPUT;
            }
            if (unlikely(a2_slen > 1000000000)) {
              putc_unlocked('\n', stdout);
              logerrputs("Error: Allele code in .bgen file has more than 1 billion characters.\n");
              goto OxBgenToPgen_ret_MALFORMED_INPUT;
            }
            if (unlikely(a2_slen + S_CAST(uintptr_t, a2_ptr - R_CAST(char*, loadbuf)) > loadbuf_size)) {
              goto OxBgenToPgen_ret_NOMEM;
            }
            if (unlikely(!fread_unlocked(a2_ptr, a2_slen, 1, bgenfile))) {
              goto OxBgenToPgen_ret_READ_FAIL;
            }
            if (compression_mode) {
#ifdef __LP64__
              compressed_block_byte_ct = 0;
#endif
              if (unlikely(!fread_unlocked(&compressed_block_byte_ct, 4, 1, bgenfile))) {
                goto OxBgenToPgen_ret_READ_FAIL;
              }
            }
            pvar_cswritep = chrtoa(cip, cur_chr_code, pvar_cswritep);
            *pvar_cswritep++ = '\t';
            if (unlikely(cur_bp > kMaxBp)) {
              putc_unlocked('\n', stdout);
              logerrputs("Error: Invalid bp coordinate (> 2^31 - 2) in .bgen file\n");
              goto OxBgenToPgen_ret_MALFORMED_INPUT;
            }
            if (par_warn_code == cur_chr_code) {
              if (unlikely(!sex_info_avail)) {
                goto OxBgenToPgen_ret_SLOPPY_CHRX_1;
              }
              if ((cur_bp <= kPAR1IntersectionLast) || (cur_bp >= kPAR2IntersectionFirst)) {
                if (unlikely(!is_sortvars)) {
                  goto OxBgenToPgen_ret_SLOPPY_CHRX_2;
                }
                print_splitpar_warning = 1;
                par_warn_code = UINT32_MAX;
              }
            }
            pvar_cswritep = u32toa_x(cur_bp, '\t', pvar_cswritep);
            pvar_cswritep = memcpyax(pvar_cswritep, rsid_start, rsid_slen, '\t');
            if (prov_ref_allele_second) {
              uint32_t swap_slen = a1_slen;
              a1_slen = a2_slen;
              a2_slen = swap_slen;
              char* swap_ptr = a1_ptr;
              a1_ptr = a2_ptr;
              a2_ptr = swap_ptr;
            }
            if (unlikely(Cswrite(&pvar_css, &pvar_cswritep))) {
              goto OxBgenToPgen_ret_WRITE_FAIL;
            }
            if (a1_slen < kMaxMediumLine) {
              pvar_cswritep = memcpya(pvar_cswritep, a1_ptr, a1_slen);
            } else {
              if (unlikely(CsputsStd(a1_ptr, a1_slen, &pvar_css, &pvar_cswritep))) {
                goto OxBgenToPgen_ret_WRITE_FAIL;
              }
            }
            *pvar_cswritep++ = '\t';
            if (unlikely(Cswrite(&pvar_css, &pvar_cswritep))) {
              goto OxBgenToPgen_ret_WRITE_FAIL;
            }
            if (a2_slen < kMaxMediumLine) {
              pvar_cswritep = memcpya(pvar_cswritep, a2_ptr, a2_slen);
            } else {
              if (unlikely(CsputsStd(a2_ptr, a2_slen, &pvar_css, &pvar_cswritep))) {
                goto OxBgenToPgen_ret_WRITE_FAIL;
              }
            }
            AppendBinaryEoln(&pvar_cswritep);

            compressed_geno_starts[block_vidx] = bgen_geno_iter;
            if (compression_mode) {
              memcpy(bgen_geno_iter, &compressed_block_byte_ct, 4);
              bgen_geno_iter = &(bgen_geno_iter[4]);
            }
            if (unlikely(fread_checked(bgen_geno_iter, compressed_block_byte_ct, bgenfile))) {
              goto OxBgenToPgen_ret_READ_FAIL;
            }
            bgen_geno_iter = &(bgen_geno_iter[compressed_block_byte_ct]);
            ++block_vidx;
          }
        }
        if (vidx_start) {
          JoinThreads(&tg);
          reterr = ctx.reterr;
          if (unlikely(reterr)) {
            goto OxBgenToPgen_ret_bgen11_thread_fail;
          }
        }
        if (!IsLastBlock(&tg)) {
          common.cur_block_size = cur_block_write_ct;
          if (vidx_start + cur_block_write_ct == variant_ct) {
            DeclareLastThreadBlock(&tg);
          }
          if (unlikely(SpawnThreads(&tg))) {
            goto OxBgenToPgen_ret_THREAD_CREATE_FAIL;
          }
        }
        parity = 1 - parity;
        if (vidx_start) {
          // write *previous* block results
          uintptr_t* write_genovec_iter = ctx.write_genovecs[parity];
          uint32_t* write_dosage_ct_iter = ctx.write_dosage_cts[parity];
          uintptr_t* write_dosage_present_iter = ctx.write_dosage_presents[parity];
          Dosage* write_dosage_main_iter = ctx.write_dosage_mains[parity];
          for (uint32_t vidx = vidx_start - prev_block_write_ct; vidx != vidx_start; ++vidx) {
            const uint32_t cur_dosage_ct = *write_dosage_ct_iter++;
            if (!cur_dosage_ct) {
              if (unlikely(SpgwAppendBiallelicGenovec(write_genovec_iter, &spgw))) {
                goto OxBgenToPgen_ret_WRITE_FAIL;
              }
            } else {
              reterr = SpgwAppendBiallelicGenovecDosage16(write_genovec_iter, write_dosage_present_iter, write_dosage_main_iter, cur_dosage_ct, &spgw);
              if (unlikely(reterr)) {
                goto OxBgenToPgen_ret_1;
              }
            }
            write_genovec_iter = &(write_genovec_iter[sample_ctaw2]);
            write_dosage_present_iter = &(write_dosage_present_iter[sample_ctaw]);
            write_dosage_main_iter = &(write_dosage_main_iter[sample_ct]);
          }
        }
        if (vidx_start == variant_ct) {
          break;
        }
        if (vidx_start) {
          printf("\r--bgen: %uk variants converted.", vidx_start / 1000);
          if (vidx_start <= main_block_size) {
            fputs("    \b\b\b\b", stdout);
          }
          fflush(stdout);
        }
        vidx_start += cur_block_write_ct;
        prev_block_write_ct = cur_block_write_ct;
      }
    } else {
      // v1.2-1.3
      uintptr_t* allele_idx_offsets;
      if (unlikely(bigstack_end_alloc_w(raw_variant_ct + 1, &allele_idx_offsets))) {
        goto OxBgenToPgen_ret_NOMEM;
      }

      Bgen13DosageOrPhaseScanCtx scan_ctx;
      scan_ctx.common = &common;
      scan_ctx.bgen_import_dosage_certainty_thresholds = nullptr;
      if (import_dosage_certainty > (1.0 - kSmallEpsilon) / 3.0) {
        scan_ctx.bgen_import_dosage_certainty_thresholds = S_CAST(uint32_t*, bigstack_alloc_raw_rd((kMaxBgenImportBits + 1) * sizeof(int32_t)));
        for (uint32_t bit_precision = 1; bit_precision != (kMaxBgenImportBits + 1); ++bit_precision) {
          const uint32_t denom = (1U << bit_precision) - 1;
          const double denom_d = u31tod(denom);
          scan_ctx.bgen_import_dosage_certainty_thresholds[bit_precision] = 1 + S_CAST(int32_t, import_dosage_certainty * denom_d);
        }
      }
      // bugfix (2 Jul 2017): if max_thread_ct == 1 but there's >12GiB memory,
      //   limit to 1 thread rather than (max_thread_ct - 1)...
      uint32_t calc_thread_ct_limit = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
      if (calc_thread_ct_limit > raw_variant_ct) {
        calc_thread_ct_limit = raw_variant_ct;
      }

      scan_ctx.thread_wkspaces = S_CAST(unsigned char**, bigstack_alloc_raw_rd(calc_thread_ct_limit * sizeof(intptr_t)));
      scan_ctx.thread_bidxs[0] = S_CAST(uint32_t*, bigstack_alloc_raw_rd((calc_thread_ct_limit + 1) * sizeof(int32_t)));
      scan_ctx.thread_bidxs[1] = S_CAST(uint32_t*, bigstack_alloc_raw_rd((calc_thread_ct_limit + 1) * sizeof(int32_t)));
      uintptr_t main_block_size = 65536;
      if (unlikely(bigstack_alloc_kcp(calc_thread_ct_limit, &(scan_ctx.err_extra)) ||
                   // ***** all bigstack allocations from this point on are
                   //       reset before pass 2 *****
                   // probably want to change this to use Gparse...
                   bigstack_alloc_u16(main_block_size, &(scan_ctx.bgen_allele_cts[0])) ||
                   bigstack_alloc_u16(main_block_size, &(scan_ctx.bgen_allele_cts[1])) ||
                   bigstack_alloc_ucp(main_block_size + 1, &(common.compressed_geno_starts[0])) ||
                   bigstack_alloc_ucp(main_block_size + 1, &(common.compressed_geno_starts[1])))) {
        goto OxBgenToPgen_ret_NOMEM;
      }
      if (compression_mode) {
        if (unlikely(bigstack_alloc_u32(main_block_size, &(scan_ctx.uncompressed_genodata_byte_cts[0])) ||
                     bigstack_alloc_u32(main_block_size, &(scan_ctx.uncompressed_genodata_byte_cts[1])))) {
          goto OxBgenToPgen_ret_NOMEM;
        }
      } else {
        // defensive
        scan_ctx.uncompressed_genodata_byte_cts[0] = nullptr;
        scan_ctx.uncompressed_genodata_byte_cts[1] = nullptr;
      }

      // ploidy >2 is not supported by PLINK 2.  (A future build may have code
      // to treat those calls as missing instead of erroring out, as is done
      // with VCF ploidy >2.  But I'll wait until this case actually comes up
      // in the wild...)
      // But even without that, the diploid worst case of 65535 alleles,
      // unphased 32-bit probabilities blows past the 4GiB uncompressed record
      // size limit with just 1 sample!  Consequences:
      // * A simple way to avoid unnecessary NOMEM errors is to give each
      //   thread 4GiB of decompression workspace on the first pass.  This may
      //   greatly reduce the number of decompression worker threads we can
      //   deploy, but for the first pass that's acceptable: the worker threads
      //   will usually all exit almost immediately (since we just need to
      //   determine whether *any* phase/dosage info needs to be saved).
      // * Even 1 thread x 4GiB won't always be available, especially since we
      //   have a double-buffering workflow which requires additional
      //   allocations summing to more than twice the decompression workspace.
      //   So we need to be able to fall back to a smaller decompression
      //   workspace size, and throw NOMEM when it's insufficient.
      // * Of course, records will almost always be far smaller than 4GiB.
      //   During the first pass, we'll see every uncompressed record size
      //   (even if the decompression worker threads terminate early), so we
      //   can usually increase the number of worker threads before the second
      //   pass.
      // Overall memory allocation for first pass:
      //   loadbuf_size (~1/7, up to 2GiB) : Chromosome code/variant ID/allele
      //                                     code load buffer.
      //   mainbuf_size (~2/7) : Compressed genotype data buffer 0, up to 4GiB
      //                         per decompression thread
      //   mainbuf_size        : Compressed genotype data buffer 1
      //   mainbuf_size        : Decompression thread workspace(s)
      // Second pass:
      //   mainbuf_size (~1/6) : Decompression thread workspaces.
      //   16K                 : .bgen chromosome code load buffer.
      //   remainder (~5/6)    : Single-threaded .pgen writer, and Gparse
      //                         double-buffer.
      uintptr_t loadbuf_size = RoundDownPow2(bigstack_left() / 7, kCacheline);
      if (loadbuf_size > kMaxLongLine) {
        loadbuf_size = kMaxLongLine;
      } else if (unlikely(loadbuf_size < 2 * 65536)) {
        // don't want to worry about chromosome/variant ID buffer space checks
        // in inner loop
        goto OxBgenToPgen_ret_NOMEM;
      }
      unsigned char* loadbuf = S_CAST(unsigned char*, bigstack_alloc_raw(loadbuf_size));

      uintptr_t mainbuf_size = RoundDownPow2(bigstack_left() / 3, kCacheline);
      uint32_t calc_thread_ct = 1;
      uintptr_t thread_wkspace_size;
#ifdef __LP64__
      // hard compressed and uncompressed record length limits of 2^31 - 1
      // bytes, since these are represented as uint32s in the file.
      if (mainbuf_size > 0x100000000LLU) {
        thread_wkspace_size = 0x100000000LLU;
        mainbuf_size &= 0xffffffff00000000LLU;
        calc_thread_ct = mainbuf_size >> 32;
        if (calc_thread_ct > calc_thread_ct_limit) {
          calc_thread_ct = calc_thread_ct_limit;
          mainbuf_size = S_CAST(uintptr_t, calc_thread_ct_limit) << 32;
        }
      } else {
        thread_wkspace_size = mainbuf_size;
      }
#else
      thread_wkspace_size = mainbuf_size;
#endif
      // note that thread_wkspace_size is the size limit for a compressed
      // variant record *and* the uncompressed form

      if (main_block_size > raw_variant_ct + calc_thread_ct - 1) {
        main_block_size = raw_variant_ct + calc_thread_ct - 1;
      }
      uint32_t per_thread_block_limit = main_block_size / calc_thread_ct;
      // may as well guarantee divisibility
      main_block_size = per_thread_block_limit * calc_thread_ct;
      if (unlikely(SetThreadCt(calc_thread_ct, &tg))) {
        goto OxBgenToPgen_ret_NOMEM;
      }
      if (compression_mode == 1) {
        for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
          common.libdeflate_decompressors[tidx] = libdeflate_alloc_decompressor();
          if (unlikely(!common.libdeflate_decompressors[tidx])) {
            goto OxBgenToPgen_ret_NOMEM;
          }
        }
      }
      unsigned char* compressed_geno_bufs[2];
      compressed_geno_bufs[0] = S_CAST(unsigned char*, bigstack_alloc_raw(mainbuf_size));
      compressed_geno_bufs[1] = S_CAST(unsigned char*, bigstack_alloc_raw(mainbuf_size));
      for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
        scan_ctx.thread_wkspaces[tidx] = S_CAST(unsigned char*, bigstack_alloc_raw(thread_wkspace_size));
      }
      scan_ctx.err_info = (~0LLU) << 32;
      scan_ctx.error_on_polyploid = !(import_flags & kfImportPolyploidMissing);
      SetThreadFuncAndData(Bgen13DosageOrPhaseScanThread, &scan_ctx, &tg);

      uint32_t block_vidx = 0;

      // bgen-1.2 and -1.3 records can vary wildly in size, so we're a bit more
      // careful with load balancing here.  (possible todo: bgzf-style atomic
      // operations)
      uint32_t cur_per_thread_block_limit = 1;
      uint32_t cur_thread_block_vidx_limit = 1;
      uint32_t cur_thread_fill_idx = 0;

      uint32_t parity = 0;
      uint32_t* thread_bidxs = scan_ctx.thread_bidxs[0];
      uint16_t* bgen_allele_cts = scan_ctx.bgen_allele_cts[0];
      unsigned char** compressed_geno_starts = common.compressed_geno_starts[0];
      uint32_t* uncompressed_genodata_byte_cts = scan_ctx.uncompressed_genodata_byte_cts[0];
      unsigned char* bgen_geno_iter = compressed_geno_bufs[0];
      unsigned char* cur_geno_buf_end = &(bgen_geno_iter[thread_wkspace_size]);
      thread_bidxs[0] = 0;
      compressed_geno_starts[0] = bgen_geno_iter;
      uintptr_t* allele_idx_offsets_iter = allele_idx_offsets;
      uintptr_t tot_allele_ct = 0;
      uint32_t max_allele_ct = 2;
      uint32_t max_compressed_geno_blen = 0;
      uint32_t max_uncompressed_geno_blen = 0;
      uint32_t uncompressed_genodata_byte_ct = 0;
      uint32_t skip_chr = 0;
      // uint32_t import_max_alleles_skip_ct = 0;

      for (uint32_t variant_uidx = 0; variant_uidx != header_variant_ct; ) {
        // format is mostly identical to bgen 1.1; but there's no sample count,
        // and there is an allele count
        // logic is more similar to the second bgen 1.1 pass since we write the
        // .pvar here.
        uint16_t snpid_slen;
        {
          const uintptr_t bytes_read = fread_unlocked(&snpid_slen, 1, 2, bgenfile);
          if (bytes_read != 2) {
            if (likely(!bytes_read)) {
              putc_unlocked('\n', stdout);
              if (likely(allow_overstated_variant_ct)) {
                raw_variant_ct = variant_uidx;
                logprintf("--lax-bgen-import: .bgen file actually contains %u variant%s.\n", raw_variant_ct, (raw_variant_ct == 1)? "" : "s");
                break;
              }
              logerrprintfww("Error: .bgen file actually contains %u variant%s; header is incorrect. Add --lax-bgen-import to force import.\n", variant_uidx, (variant_uidx == 1)? "" : "s");
              goto OxBgenToPgen_ret_INCONSISTENT_INPUT;
            }
            goto OxBgenToPgen_ret_READ_FAIL;
          }
        }
        char* rsid_start = R_CAST(char*, loadbuf);
        if (!snpid_chr) {
          if (unlikely(fseeko(bgenfile, snpid_slen, SEEK_CUR))) {
            goto OxBgenToPgen_ret_READ_FAIL;
          }
        } else {
          if (unlikely(!snpid_slen)) {
            putc_unlocked('\n', stdout);
            logerrputs("Error: Length-0 SNP ID in .bgen file.\n");
            goto OxBgenToPgen_ret_INCONSISTENT_INPUT;
          }
          if (unlikely(!fread_unlocked(loadbuf, snpid_slen, 1, bgenfile))) {
            goto OxBgenToPgen_ret_READ_FAIL;
          }
          // loadbuf[snpid_slen] = '\0';
          rsid_start = R_CAST(char*, &(loadbuf[snpid_slen + 1]));
        }
        uint16_t rsid_slen;
        if (unlikely(!fread_unlocked(&rsid_slen, 2, 1, bgenfile))) {
          goto OxBgenToPgen_ret_READ_FAIL;
        }
        if (unlikely(!rsid_slen)) {
          putc_unlocked('\n', stdout);
          logerrputs("Error: Length-0 rsID in .bgen file.\n");
          goto OxBgenToPgen_ret_MALFORMED_INPUT;
        }
        if (unlikely(!fread_unlocked(rsid_start, rsid_slen, 1, bgenfile))) {
          goto OxBgenToPgen_ret_READ_FAIL;
        }
        char* loadbuf_iter = &(rsid_start[rsid_slen]);
        char* chr_name_start = loadbuf_iter;
        uint16_t chr_name_slen;
        if (unlikely(!fread_unlocked(&chr_name_slen, 2, 1, bgenfile))) {
          goto OxBgenToPgen_ret_READ_FAIL;
        }
        if (ox_single_chr_str) {
          if (unlikely(fseeko(bgenfile, chr_name_slen, SEEK_CUR))) {
            goto OxBgenToPgen_ret_READ_FAIL;
          }
        } else {
          if (!snpid_chr) {
            if (unlikely(!chr_name_slen)) {
              putc_unlocked('\n', stdout);
              logerrputs("Error: Length-0 chromosome ID in .bgen file.\n");
              goto OxBgenToPgen_ret_INCONSISTENT_INPUT;
            }
            if (unlikely(!fread_unlocked(chr_name_start, chr_name_slen, 1, bgenfile))) {
              goto OxBgenToPgen_ret_READ_FAIL;
            }
            if (strequal_k(chr_name_start, "NA", chr_name_slen)) {
              chr_name_start[0] = '0';
              chr_name_slen = 1;
            }
          } else {
            if (unlikely(fseeko(bgenfile, chr_name_slen, SEEK_CUR))) {
              goto OxBgenToPgen_ret_READ_FAIL;
            }
            chr_name_start = R_CAST(char*, loadbuf);
            chr_name_slen = snpid_slen;
          }
          // chromosome ID length restriction enforced here, so we don't check
          // earlier
          chr_name_start[chr_name_slen] = '\0';
          reterr = GetOrAddChrCode(chr_name_start, "--bgen file", 0, chr_name_slen, prohibit_extra_chr, cip, &cur_chr_code);
          if (unlikely(reterr)) {
            goto OxBgenToPgen_ret_1;
          }
          skip_chr = !IsSet(cip->chr_mask, cur_chr_code);
        }

        uint32_t cur_bp;
        uint32_t cur_allele_ct = 0;
        if (unlikely((!fread_unlocked(&cur_bp, 4, 1, bgenfile)) ||
                     (!fread_unlocked(&cur_allele_ct, 2, 1, bgenfile)))) {
          goto OxBgenToPgen_ret_READ_FAIL;
        }
        if (unlikely(cur_allele_ct < 2)) {
          // this is undefined in the 1.3 standard; prohibit for now
          putc_unlocked('\n', stdout);
          logerrputs("Error: .bgen variant has fewer than two alleles.\n");
          goto OxBgenToPgen_ret_MALFORMED_INPUT;
        }
        ++variant_uidx;
        if (!(variant_uidx % 1000)) {
          printf("\r--bgen: %uk variants scanned.", variant_uidx / 1000);
          fflush(stdout);
        }

        // the "cur_allele_ct > 2" part is a temporary kludge
        if (skip_chr || (cur_allele_ct > 2)) {
          if (!skip_chr) {
            if (cur_allele_ct > import_max_allele_ct) {
              // ++import_max_alleles_skip_ct;
            } else {
              ++multiallelic_tmp_skip_ct;
            }
          } else {
            chr_filter_exists = 1;
          }
          for (uint32_t allele_idx = 0; allele_idx != cur_allele_ct; ++allele_idx) {
            uint32_t cur_allele_slen;
            if (unlikely((!fread_unlocked(&cur_allele_slen, 4, 1, bgenfile)) ||
                         fseeko(bgenfile, cur_allele_slen, SEEK_CUR))) {
              goto OxBgenToPgen_ret_READ_FAIL;
            }
          }
          uint32_t genodata_byte_ct;
          if (unlikely((!fread_unlocked(&genodata_byte_ct, 4, 1, bgenfile)) ||
                       fseeko(bgenfile, genodata_byte_ct, SEEK_CUR))) {
            goto OxBgenToPgen_ret_READ_FAIL;
          }
          continue;
        }
        if (unlikely(rsid_slen > kMaxIdSlen)) {
          // enforce this iff we aren't skipping
          putc_unlocked('\n', stdout);
          logerrputs("Error: Variant names are limited to " MAX_ID_SLEN_STR " characters.\n");
          goto OxBgenToPgen_ret_MALFORMED_INPUT;
        }
        // special handling of first two alleles since either may be
        // reference, so we may need to swap order
        char* a1_ptr = loadbuf_iter;
        uint32_t a1_slen;
        if (unlikely(!fread_unlocked(&a1_slen, 4, 1, bgenfile))) {
          goto OxBgenToPgen_ret_READ_FAIL;
        }
        if (unlikely(!a1_slen)) {
          putc_unlocked('\n', stdout);
          logerrputs("Error: Empty allele code in .bgen file.\n");
          goto OxBgenToPgen_ret_MALFORMED_INPUT;
        }
        if (unlikely(a1_slen > 1000000000)) {
          putc_unlocked('\n', stdout);
          logerrputs("Error: Allele code in .bgen file has more than 1 billion characters.\n");
          goto OxBgenToPgen_ret_MALFORMED_INPUT;
        }
        if (unlikely(a1_slen + S_CAST(uintptr_t, a1_ptr - R_CAST(char*, loadbuf)) > loadbuf_size)) {
          goto OxBgenToPgen_ret_NOMEM;
        }
        if (unlikely(!fread_unlocked(a1_ptr, a1_slen, 1, bgenfile))) {
          goto OxBgenToPgen_ret_READ_FAIL;
        }
        char* a2_ptr = &(a1_ptr[a1_slen]);
        uint32_t a2_slen;
        if (unlikely(!fread_unlocked(&a2_slen, 4, 1, bgenfile))) {
          goto OxBgenToPgen_ret_READ_FAIL;
        }
        if (unlikely(!a2_slen)) {
          putc_unlocked('\n', stdout);
          logerrputs("Error: Empty allele code in .bgen file.\n");
          goto OxBgenToPgen_ret_MALFORMED_INPUT;
        }
        if (unlikely(a2_slen > 1000000000)) {
          putc_unlocked('\n', stdout);
          logerrputs("Error: Allele code in .bgen file has more than 1 billion characters.\n");
          goto OxBgenToPgen_ret_MALFORMED_INPUT;
        }
        if (unlikely(a2_slen + S_CAST(uintptr_t, a2_ptr - R_CAST(char*, loadbuf)) > loadbuf_size)) {
          goto OxBgenToPgen_ret_NOMEM;
        }
        if (unlikely(!fread_unlocked(a2_ptr, a2_slen, 1, bgenfile))) {
          goto OxBgenToPgen_ret_READ_FAIL;
        }
        pvar_cswritep = chrtoa(cip, cur_chr_code, pvar_cswritep);
        *pvar_cswritep++ = '\t';
        if (unlikely(cur_bp > kMaxBp)) {
          putc_unlocked('\n', stdout);
          logerrputs("Error: Invalid bp coordinate (> 2^31 - 2) in .bgen file\n");
          goto OxBgenToPgen_ret_MALFORMED_INPUT;
        }
        if (par_warn_code == cur_chr_code) {
          if (unlikely(!sex_info_avail)) {
            goto OxBgenToPgen_ret_SLOPPY_CHRX_1;
          }
          if ((cur_bp <= kPAR1IntersectionLast) || (cur_bp >= kPAR2IntersectionFirst)) {
            if (unlikely(!is_sortvars)) {
              goto OxBgenToPgen_ret_SLOPPY_CHRX_2;
            }
            print_splitpar_warning = 1;
            par_warn_code = UINT32_MAX;
          }
        }
        pvar_cswritep = u32toa_x(cur_bp, '\t', pvar_cswritep);
        pvar_cswritep = memcpyax(pvar_cswritep, rsid_start, rsid_slen, '\t');
        if (prov_ref_allele_second) {
          const uint32_t swap_slen = a1_slen;
          a1_slen = a2_slen;
          a2_slen = swap_slen;
          char* swap_ptr = a1_ptr;
          a1_ptr = a2_ptr;
          a2_ptr = swap_ptr;
        }
        // allele codes may be too large for write buffer
        if (unlikely(Cswrite(&pvar_css, &pvar_cswritep))) {
          goto OxBgenToPgen_ret_WRITE_FAIL;
        }
        if (a1_slen < kMaxMediumLine) {
          pvar_cswritep = memcpya(pvar_cswritep, a1_ptr, a1_slen);
        } else {
          if (unlikely(CsputsStd(a1_ptr, a1_slen, &pvar_css, &pvar_cswritep))) {
            goto OxBgenToPgen_ret_WRITE_FAIL;
          }
        }
        *pvar_cswritep++ = '\t';
        if (unlikely(Cswrite(&pvar_css, &pvar_cswritep))) {
          goto OxBgenToPgen_ret_WRITE_FAIL;
        }
        if (a2_slen < kMaxMediumLine) {
          pvar_cswritep = memcpya(pvar_cswritep, a2_ptr, a2_slen);
        } else {
          if (unlikely(CsputsStd(a2_ptr, a2_slen, &pvar_css, &pvar_cswritep))) {
            goto OxBgenToPgen_ret_WRITE_FAIL;
          }
        }
        for (uint32_t allele_idx = 2; allele_idx != cur_allele_ct; ++allele_idx) {
          // (can't actually reach here yet since we're skipping multiallelics
          // for now)
          // safe to use entire loadbuf for this
          assert(0);
          uint32_t cur_allele_slen;
          if (unlikely(!fread_unlocked(&cur_allele_slen, 4, 1, bgenfile))) {
            goto OxBgenToPgen_ret_READ_FAIL;
          }
          if (unlikely(!cur_allele_slen)) {
            putc_unlocked('\n', stdout);
            logerrputs("Error: Empty allele code in .bgen file.\n");
            goto OxBgenToPgen_ret_MALFORMED_INPUT;
          }
          if (unlikely(cur_allele_slen > 1000000000)) {
            putc_unlocked('\n', stdout);
            logerrputs("Error: Allele code in .bgen file has more than 1 billion characters.\n");
            goto OxBgenToPgen_ret_MALFORMED_INPUT;
          }
          if (unlikely(cur_allele_slen > loadbuf_size)) {
            goto OxBgenToPgen_ret_NOMEM;
          }
          if (unlikely(!fread_unlocked(loadbuf, cur_allele_slen, 1, bgenfile))) {
            goto OxBgenToPgen_ret_READ_FAIL;
          }
          *pvar_cswritep++ = ',';
          if (unlikely(Cswrite(&pvar_css, &pvar_cswritep))) {
            goto OxBgenToPgen_ret_WRITE_FAIL;
          }
          if (cur_allele_slen < kMaxMediumLine) {
            pvar_cswritep = memcpya(pvar_cswritep, loadbuf, cur_allele_slen);
          } else {
            if (unlikely(CsputsStd(R_CAST(char*, loadbuf), cur_allele_slen, &pvar_css, &pvar_cswritep))) {
              goto OxBgenToPgen_ret_WRITE_FAIL;
            }
          }
        }

        AppendBinaryEoln(&pvar_cswritep);
        *allele_idx_offsets_iter++ = tot_allele_ct;
        tot_allele_ct += cur_allele_ct;
        if (cur_allele_ct > max_allele_ct) {
          max_allele_ct = cur_allele_ct;
        }
        uint32_t genodata_byte_ct;
        if (unlikely(!fread_unlocked(&genodata_byte_ct, 4, 1, bgenfile))) {
          goto OxBgenToPgen_ret_READ_FAIL;
        }
        if (genodata_byte_ct > max_compressed_geno_blen) {
          max_compressed_geno_blen = genodata_byte_ct;
        }
        if (uncompressed_genodata_byte_cts) {
          if (unlikely(genodata_byte_ct < 4)) {
            logerrputs("Error: Invalid compressed block length in .bgen file.\n");
            goto OxBgenToPgen_ret_MALFORMED_INPUT;
          }
          if (unlikely(!fread_unlocked(&uncompressed_genodata_byte_ct, 4, 1, bgenfile))) {
            goto OxBgenToPgen_ret_READ_FAIL;
          }
          if (uncompressed_genodata_byte_ct > max_uncompressed_geno_blen) {
            max_uncompressed_geno_blen = uncompressed_genodata_byte_ct;
          }
          genodata_byte_ct -= 4;
        }
        if (dosage_exists) {
          if (unlikely(fseeko(bgenfile, genodata_byte_ct, SEEK_CUR))) {
            goto OxBgenToPgen_ret_READ_FAIL;
          }
          ++variant_ct;
          continue;
        }

        if ((block_vidx == cur_thread_block_vidx_limit) || (S_CAST(uintptr_t, cur_geno_buf_end - bgen_geno_iter) < genodata_byte_ct)) {
          if (unlikely(!block_vidx)) {
            goto OxBgenToPgen_ret_NOMEM;
          }
          thread_bidxs[++cur_thread_fill_idx] = block_vidx;
          if (cur_thread_fill_idx == calc_thread_ct) {
            parity = 1 - parity;
            if (ThreadsAreActive(&tg)) {
              // process *previous* block results
              JoinThreads(&tg);
              reterr = S_CAST(PglErr, scan_ctx.err_info & 255);
              if (unlikely(reterr)) {
                PrintBgenImportErr(scan_ctx.err_extra, scan_ctx.err_info, scan_ctx.vidx_start);
                goto OxBgenToPgen_ret_1;
              }
              dosage_exists = scan_ctx.dosage_exists;
              if (dosage_exists) {
                // don't need to scan for any more dosages
                StopThreads(&tg);

                // however, unlike bgen-1.1 case, we can never do full
                // early-exit since we have to scan for multiallelic variants:
                // writer must be initialized with (i) an accurate variant
                // count, which is temporarily affected by skipped multiallelic
                // variants, and (ii) when we no longer skip them, the
                // PgenWriter constructor still needs a maximum allele count so
                // it can allocate properly-sized buffers.
                if (unlikely(fseeko(bgenfile, genodata_byte_ct, SEEK_CUR))) {
                  goto OxBgenToPgen_ret_READ_FAIL;
                }
                ++variant_ct;
                continue;
              }
            }
            scan_ctx.vidx_start = variant_ct;
            if (unlikely(SpawnThreads(&tg))) {
              goto OxBgenToPgen_ret_THREAD_CREATE_FAIL;
            }
            compressed_geno_starts = common.compressed_geno_starts[parity];
            uncompressed_genodata_byte_cts = scan_ctx.uncompressed_genodata_byte_cts[parity];
            thread_bidxs = scan_ctx.thread_bidxs[parity];
            bgen_allele_cts = scan_ctx.bgen_allele_cts[parity];
            bgen_geno_iter = compressed_geno_bufs[parity];
            thread_bidxs[0] = 0;
            compressed_geno_starts[0] = bgen_geno_iter;
            variant_ct += block_vidx;
            block_vidx = 0;
            if (cur_per_thread_block_limit < per_thread_block_limit) {
              cur_per_thread_block_limit *= 2;
              if (cur_per_thread_block_limit > per_thread_block_limit) {
                cur_per_thread_block_limit = per_thread_block_limit;
              }
            }
            cur_thread_block_vidx_limit = 0;
            cur_thread_fill_idx = 0;
          }
          cur_geno_buf_end = &(bgen_geno_iter[thread_wkspace_size]);
          cur_thread_block_vidx_limit += cur_per_thread_block_limit;
        }
        bgen_allele_cts[block_vidx] = cur_allele_ct;
        if (uncompressed_genodata_byte_cts) {
          uncompressed_genodata_byte_cts[block_vidx] = uncompressed_genodata_byte_ct;
        }
        if (unlikely(fread_checked(bgen_geno_iter, genodata_byte_ct, bgenfile))) {
          goto OxBgenToPgen_ret_READ_FAIL;
        }
        bgen_geno_iter = &(bgen_geno_iter[genodata_byte_ct]);
        compressed_geno_starts[++block_vidx] = bgen_geno_iter;
      }
      variant_ct += block_vidx;
      if (multiallelic_tmp_skip_ct) {
        putc_unlocked('\n', stdout);
        logerrprintfww("Warning: %u multiallelic variant%s skipped%s (not yet supported).\n", multiallelic_tmp_skip_ct, (multiallelic_tmp_skip_ct == 1)? "" : "s", (import_max_allele_ct < 0x7ffffffe)? ", on top of --import-max-alleles filter" : "");
      }
      if (unlikely(!variant_ct)) {
        putc_unlocked('\n', stdout);
        char* write_iter = strcpya_k(g_logbuf, "Error: All ");
        if (multiallelic_tmp_skip_ct) {
          write_iter = strcpya_k(write_iter, "remaining ");
        }
        write_iter = u32toa(raw_variant_ct - multiallelic_tmp_skip_ct, write_iter);
        write_iter = strcpya_k(write_iter, " variant");
        if (raw_variant_ct - multiallelic_tmp_skip_ct != 1) {
          *write_iter++ = 's';
        }
        write_iter = strcpya_k(write_iter, " in .bgen file excluded by ");
        AppendLoadFilterFlagnames(load_filter_log_import_flags, &write_iter);
        strcpy_k(write_iter, ".\n");
        goto OxBgenToPgen_ret_INCONSISTENT_INPUT_WW;
      }
      if (variant_ct == block_vidx) {
        // with multiple threads, there's no guarantee that even the first
        // decompression job has launched (e.g. there's only 1 variant on the
        // relevant chromosome in the entire .bgen, and calc_thread_ct == 2).
        thread_bidxs[cur_thread_fill_idx + 1] = block_vidx;
        scan_ctx.vidx_start = 0;
        if (unlikely(SpawnThreads(&tg))) {
          goto OxBgenToPgen_ret_THREAD_CREATE_FAIL;
        }
        block_vidx = 0;
      }
      if (ThreadsAreActive(&tg)) {
        JoinThreads(&tg);
        reterr = S_CAST(PglErr, scan_ctx.err_info & 255);
        if (unlikely(reterr)) {
          PrintBgenImportErr(scan_ctx.err_extra, scan_ctx.err_info, scan_ctx.vidx_start);
          goto OxBgenToPgen_ret_1;
        }
        if ((!block_vidx) || scan_ctx.dosage_exists) {
          // ignore thread_bidxs[] in this case
          StopThreads(&tg);
        } else {
          for (; cur_thread_fill_idx != calc_thread_ct; ) {
            // save endpoint for current thread, and tell any leftover threads
            // to do nothing
            thread_bidxs[++cur_thread_fill_idx] = block_vidx;
          }
          DeclareLastThreadBlock(&tg);
          scan_ctx.vidx_start = variant_ct - block_vidx;
          SpawnThreads(&tg);
          JoinThreads(&tg);
          reterr = S_CAST(PglErr, scan_ctx.err_info & 255);
          if (unlikely(reterr)) {
            PrintBgenImportErr(scan_ctx.err_extra, scan_ctx.err_info, scan_ctx.vidx_start);
            goto OxBgenToPgen_ret_bgen13_thread_fail;
          }
        }
        dosage_exists = scan_ctx.dosage_exists;
      }

      if (max_allele_ct == 2) {
        allele_idx_offsets = nullptr;
        BigstackEndReset(bigstack_end_mark);
      } else {
        // not yet possible
        reterr = kPglRetNotYetSupported;
        goto OxBgenToPgen_ret_1;
        *allele_idx_offsets_iter = tot_allele_ct;
      }
      if (unlikely(fseeko(bgenfile, initial_uints[0] + 4, SEEK_SET))) {
        goto OxBgenToPgen_ret_READ_FAIL;
      }
      snprintf(outname_end, kMaxOutfnameExtBlen, ".pgen");
      uintptr_t spgw_alloc_cacheline_ct;
      uint32_t max_vrec_len;
      reterr = SpgwInitPhase1(outname, allele_idx_offsets, nullptr, variant_ct, sample_ct, max_allele_ct, kPgenWriteBackwardSeek, dosage_exists? (kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent) : kfPgenGlobal0, (oxford_import_flags & kfOxfordImportRefUnknown)? 2 : 1, &spgw, &spgw_alloc_cacheline_ct, &max_vrec_len);
      if (unlikely(reterr)) {
        if (reterr == kPglRetOpenFail) {
          logerrprintfww(kErrprintfFopen, outname, strerror(errno));
        }
        goto OxBgenToPgen_ret_1;
      }

      BigstackReset(scan_ctx.bgen_allele_cts[0]);

      // only needs to fit chromosome codes in second pass
      loadbuf = S_CAST(unsigned char*, bigstack_alloc_raw_rd(kMaxIdBlen));
      unsigned char* spgw_alloc;
      if (unlikely(bigstack_alloc_uc(spgw_alloc_cacheline_ct * kCacheline, &spgw_alloc))) {
        goto OxBgenToPgen_ret_NOMEM;
      }
      if (!uncompressed_genodata_byte_cts) {
        max_uncompressed_geno_blen = max_compressed_geno_blen;
      }
      // Now that we know max_uncompressed_geno_blen, try to increase
      // calc_thread_ct, and resize ctx.thread_wkspaces[tidx] (and also resize
      // compressed_geno_bufs[] in next step).
      // Additional *6 in denominator since we want to limit these allocations
      // to 1/6 of remaining workspace.
      thread_wkspace_size = RoundUpPow2(max_uncompressed_geno_blen, kCacheline);
      // bugfix (16 Jul 2017): was computing cachelines_avail, not bytes_avail
      uintptr_t bytes_avail = RoundDownPow2(bigstack_left() / 6, kCacheline);
      uint32_t old_calc_thread_ct = calc_thread_ct;
      if (calc_thread_ct_limit * thread_wkspace_size <= bytes_avail) {
        calc_thread_ct = calc_thread_ct_limit;
      } else {
        calc_thread_ct = bytes_avail / thread_wkspace_size;
        if (unlikely(!calc_thread_ct)) {
          goto OxBgenToPgen_ret_NOMEM;
        }
      }
      if (unlikely(SetThreadCt(calc_thread_ct, &tg))) {
        goto OxBgenToPgen_ret_NOMEM;
      }
      Bgen13GenoToPgenCtx ctx;
      ctx.common = &common;
      ctx.err_extra = scan_ctx.err_extra;
      ctx.error_on_polyploid = !(import_flags & kfImportPolyploidMissing);
      ctx.hard_call_halfdist = kDosage4th - hard_call_thresh;
      ctx.bgen_import_dosage_certainty_thresholds = scan_ctx.bgen_import_dosage_certainty_thresholds;
      ctx.prov_ref_allele_second = prov_ref_allele_second;
      ctx.thread_wkspaces = scan_ctx.thread_wkspaces;
      for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
        ctx.thread_wkspaces[tidx] = S_CAST(unsigned char*, bigstack_alloc_raw(thread_wkspace_size));
      }
      ctx.thread_bidxs[0] = scan_ctx.thread_bidxs[0];
      ctx.thread_bidxs[1] = scan_ctx.thread_bidxs[1];
      if (compression_mode == 1) {
        for (uint32_t tidx = old_calc_thread_ct; tidx < calc_thread_ct; ++tidx) {
          common.libdeflate_decompressors[tidx] = libdeflate_alloc_decompressor();
          if (unlikely(!common.libdeflate_decompressors[tidx])) {
            goto OxBgenToPgen_ret_NOMEM;
          }
        }
      }
      bytes_avail -= thread_wkspace_size * calc_thread_ct;

      // const GparseFlags gparse_flags = dosage_exists? (kfGparseHphase | kfGparseDosage | kfGparseDphase) : kfGparse0;
      // unconditionally reserve space for everything but multiallelics for now
      const GparseFlags gparse_flags = kfGparseHphase | kfGparseDosage | kfGparseDphase;
      const uint64_t max_write_byte_ct = GparseWriteByteCt(sample_ct, max_allele_ct, gparse_flags);
      uintptr_t cachelines_avail = bigstack_left() / kCacheline;
      if (unlikely(cachelines_avail < 4)) {
        goto OxBgenToPgen_ret_NOMEM;
      }
      // we're making 4 allocations; be pessimistic re: rounding
      cachelines_avail -= 4;
      // may as well include calc_thread_ct thanks to potential
      // per_thread_byte_limit adverse rounding
      const uint64_t max_bytes_req_per_variant = sizeof(GparseRecord) + MAXV(max_compressed_geno_blen, max_write_byte_ct) + calc_thread_ct;
      if (unlikely(cachelines_avail * kCacheline < 2 * max_bytes_req_per_variant)) {
        goto OxBgenToPgen_ret_NOMEM;
      }
      uintptr_t min_bytes_req_per_variant = sizeof(GparseRecord) + GparseWriteByteCt(sample_ct, 2, gparse_flags);
      main_block_size = (cachelines_avail * kCacheline) / (min_bytes_req_per_variant * 2);
      // this is arbitrary, there's no connection to kPglVblockSize
      if (main_block_size > 65536) {
        main_block_size = 65536;
      }
      // divide by 2 for better parallelism in small-variant-count case
      // round up per_thread_block_limit so we only have two blocks
      if (main_block_size > DivUp(raw_variant_ct, 2)) {
        main_block_size = DivUp(raw_variant_ct, 2) + calc_thread_ct - 1;
      }
      // may as well guarantee divisibility
      per_thread_block_limit = main_block_size / calc_thread_ct;
      main_block_size = per_thread_block_limit * calc_thread_ct;
      ctx.gparse[0] = S_CAST(GparseRecord*, bigstack_alloc_raw_rd(main_block_size * sizeof(GparseRecord)));
      ctx.gparse[1] = S_CAST(GparseRecord*, bigstack_alloc_raw_rd(main_block_size * sizeof(GparseRecord)));
      ctx.block_allele_idx_offsets[0] = nullptr;  // defensive
      ctx.block_allele_idx_offsets[1] = nullptr;
      ctx.err_info = (~0LLU) << 32;
      SetThreadFuncAndData(Bgen13GenoToPgenThread, &ctx, &tg);
      cachelines_avail = bigstack_left() / (kCacheline * 2);
      if (unlikely(bigstack_alloc_uc(cachelines_avail * kCacheline, &(compressed_geno_bufs[0])) ||
                   bigstack_alloc_uc(cachelines_avail * kCacheline, &(compressed_geno_bufs[1])))) {
        assert(0);
        goto OxBgenToPgen_ret_NOMEM;
      }
      const uintptr_t per_thread_byte_limit = (cachelines_avail * kCacheline) / calc_thread_ct;

      SpgwInitPhase2(max_vrec_len, &spgw, spgw_alloc);

      // Main workflow:
      // 1. Set n=0, load genotype data for first main_block_size variants
      //    while writing .pvar
      //
      // 2. Spawn threads processing batch n genotype data
      // 3. If n>0, write results for block (n-1)
      // 4. Increment n by 1
      // 5. Load/write-.pvar for batch (n+1) unless eof
      // 6. Join threads
      // 7. Goto step 2 unless eof
      //
      // 8. Write results for last block
      uint32_t prev_block_write_ct = 0;
      uint32_t prev_genodata_byte_ct = 0;
      uintptr_t prev_record_byte_ct = 0;
      uint32_t prev_allele_ct = 0;
      parity = 0;
      for (uint32_t vidx_start = 0; ; ) {
        uint32_t cur_block_write_ct = 0;
        if (!IsLastBlock(&tg)) {
          const uint32_t block_vidx_limit = variant_ct - vidx_start;
          cur_thread_block_vidx_limit = MINV(block_vidx_limit, per_thread_block_limit);
          cur_thread_fill_idx = 0;
          thread_bidxs = ctx.thread_bidxs[parity];
          GparseRecord* cur_gparse = ctx.gparse[parity];
          if (allele_idx_offsets) {
            ctx.block_allele_idx_offsets[parity] = &(allele_idx_offsets[vidx_start]);
          }
          bgen_geno_iter = compressed_geno_bufs[parity];
          unsigned char* cur_thread_byte_stop = &(bgen_geno_iter[per_thread_byte_limit]);
          thread_bidxs[0] = 0;
          block_vidx = 0;
          // strictly speaking, prev_genodata_byte_ct and genodata_byte_ct
          // can be collapsed into one variable, as well as
          // {block_vidx, cur_block_write_ct}, but not a big deal if the
          // compiler fails to see this
          uint32_t genodata_byte_ct = prev_genodata_byte_ct;
          uint32_t cur_allele_ct = prev_allele_ct;
          uintptr_t record_byte_ct = prev_record_byte_ct;
          GparseRecord* grp;
          if (!genodata_byte_ct) {
            goto OxBgenToPgen_load13_start;
          }
          // we may stop before main_block_size due to insufficient space in
          // compressed_geno_buf.  if so, the file pointer is right before the
          // genotype data, rather than at the beginning of a variant record.
          skip_chr = 0;  // defensive
          while (1) {
            grp = &(cur_gparse[block_vidx]);
            grp->record_start = bgen_geno_iter;
            grp->flags = gparse_flags;
            grp->metadata.read_bgen.input_byte_ct = genodata_byte_ct;
            if (compression_mode) {
              grp->metadata.read_bgen.uncompressed_byte_ct = uncompressed_genodata_byte_ct;
            }
            if (unlikely(fread_checked(bgen_geno_iter, genodata_byte_ct, bgenfile))) {
              goto OxBgenToPgen_ret_READ_FAIL;
            }
            bgen_geno_iter = &(bgen_geno_iter[record_byte_ct]);
            ++block_vidx;

            uint16_t snpid_slen;
            // true iff this is the last variant we're keeping in the entire
            // file
            if (block_vidx == block_vidx_limit) {
              for (; cur_thread_fill_idx != calc_thread_ct; ) {
                // save endpoint for current thread, and tell any leftover
                // threads to do nothing
                thread_bidxs[++cur_thread_fill_idx] = block_vidx;
              }
              break;
            }
          OxBgenToPgen_load13_start:
            if (unlikely(!fread_unlocked(&snpid_slen, 2, 1, bgenfile))) {
              goto OxBgenToPgen_ret_READ_FAIL;
            }

            if (!snpid_chr) {
              if (unlikely(fseeko(bgenfile, snpid_slen, SEEK_CUR))) {
                goto OxBgenToPgen_ret_READ_FAIL;
              }
            } else {
              if (unlikely(!fread_unlocked(loadbuf, snpid_slen, 1, bgenfile))) {
                goto OxBgenToPgen_ret_READ_FAIL;
              }
              loadbuf[snpid_slen] = '\0';
            }
            uint16_t rsid_slen;
            if (unlikely(!fread_unlocked(&rsid_slen, 2, 1, bgenfile))) {
              goto OxBgenToPgen_ret_READ_FAIL;
            }
            if (unlikely(fseeko(bgenfile, rsid_slen, SEEK_CUR))) {
              goto OxBgenToPgen_ret_READ_FAIL;
            }
            uint16_t chr_name_slen;
            if (unlikely(!fread_unlocked(&chr_name_slen, 2, 1, bgenfile))) {
              goto OxBgenToPgen_ret_READ_FAIL;
            }
            if (ox_single_chr_str) {
              if (unlikely(fseeko(bgenfile, chr_name_slen, SEEK_CUR))) {
                goto OxBgenToPgen_ret_READ_FAIL;
              }
            } else {
              if (!snpid_chr) {
                if (unlikely(!fread_unlocked(loadbuf, chr_name_slen, 1, bgenfile))) {
                  goto OxBgenToPgen_ret_READ_FAIL;
                }
                if (strequal_k(R_CAST(char*, loadbuf), "NA", chr_name_slen)) {
                  memcpy_k(loadbuf, "0", 2);
                  chr_name_slen = 1;
                } else {
                  loadbuf[chr_name_slen] = '\0';
                }
              } else {
                if (unlikely(fseeko(bgenfile, chr_name_slen, SEEK_CUR))) {
                  goto OxBgenToPgen_ret_READ_FAIL;
                }
                chr_name_slen = snpid_slen;
              }
              if (chr_filter_exists) {
                const uint32_t cur_chr_code2 = GetChrCode(R_CAST(char*, loadbuf), cip, chr_name_slen);

                // we scanned all the variants
                assert(!IsI32Neg(cur_chr_code2));

                skip_chr = !IsSet(cip->chr_mask, cur_chr_code2);
              }
            }

            uint32_t cur_bp;  // ignore in this pass
            if (unlikely(!fread_unlocked(&cur_bp, 4, 1, bgenfile))) {
              goto OxBgenToPgen_ret_READ_FAIL;
            }

            cur_allele_ct = 0;
            if (unlikely(!fread_unlocked(&cur_allele_ct, 2, 1, bgenfile))) {
              goto OxBgenToPgen_ret_READ_FAIL;
            }
            for (uint32_t allele_idx = 0; allele_idx != cur_allele_ct; ++allele_idx) {
              uint32_t allele_slen;
              if (unlikely(!fread_unlocked(&allele_slen, 4, 1, bgenfile))) {
                goto OxBgenToPgen_ret_READ_FAIL;
              }
              if (unlikely(fseeko(bgenfile, allele_slen, SEEK_CUR))) {
                goto OxBgenToPgen_ret_READ_FAIL;
              }
            }
            if (unlikely(!fread_unlocked(&genodata_byte_ct, 4, 1, bgenfile))) {
              goto OxBgenToPgen_ret_READ_FAIL;
            }

            // "cur_allele_ct > 2" is temporary kludge
            if (skip_chr || (cur_allele_ct > 2)) {
              if (unlikely(fseeko(bgenfile, genodata_byte_ct, SEEK_CUR))) {
                goto OxBgenToPgen_ret_READ_FAIL;
              }
              goto OxBgenToPgen_load13_start;
            }
            if (compression_mode) {
              if (unlikely(!fread_unlocked(&uncompressed_genodata_byte_ct, 4, 1, bgenfile))) {
                goto OxBgenToPgen_ret_READ_FAIL;
              }
              genodata_byte_ct -= 4;
            }
            const uintptr_t write_byte_ct_limit = GparseWriteByteCt(sample_ct, cur_allele_ct, gparse_flags);
            record_byte_ct = MAXV(RoundUpPow2(genodata_byte_ct, kBytesPerVec), write_byte_ct_limit);

            if ((block_vidx == cur_thread_block_vidx_limit) || (S_CAST(uintptr_t, cur_thread_byte_stop - bgen_geno_iter) < record_byte_ct)) {
              thread_bidxs[++cur_thread_fill_idx] = block_vidx;
              if (cur_thread_fill_idx == calc_thread_ct) {
                prev_allele_ct = cur_allele_ct;
                prev_genodata_byte_ct = genodata_byte_ct;
                prev_record_byte_ct = record_byte_ct;
                break;
              }
              cur_thread_byte_stop = &(cur_thread_byte_stop[per_thread_byte_limit]);
              cur_thread_block_vidx_limit = MINV(cur_thread_block_vidx_limit + per_thread_block_limit, block_vidx_limit);
            }
          }
          cur_block_write_ct = block_vidx;
        }
        if (vidx_start) {
          JoinThreads(&tg);
          reterr = S_CAST(PglErr, ctx.err_info);
          if (unlikely(reterr)) {
            PrintBgenImportErr(ctx.err_extra, ctx.err_info, vidx_start - prev_block_write_ct);
            goto OxBgenToPgen_ret_1;
          }
        }
        if (!IsLastBlock(&tg)) {
          if (vidx_start + cur_block_write_ct == variant_ct) {
            DeclareLastThreadBlock(&tg);
          }
          if (unlikely(SpawnThreads(&tg))) {
            goto OxBgenToPgen_ret_THREAD_CREATE_FAIL;
          }
        }
        parity = 1 - parity;
        if (vidx_start) {
          // write *previous* block results
          reterr = GparseFlush(ctx.gparse[parity], allele_idx_offsets, prev_block_write_ct, &spgw);
          if (unlikely(reterr)) {
            goto OxBgenToPgen_ret_1;
          }
        }
        if (vidx_start == variant_ct) {
          break;
        }
        if (vidx_start) {
          printf("\r--bgen: %uk variants converted.", vidx_start / 1000);
          if (vidx_start <= main_block_size) {
            fputs("    \b\b\b\b", stdout);
          }
          fflush(stdout);
        }
        vidx_start += cur_block_write_ct;
        prev_block_write_ct = cur_block_write_ct;
      }
    }
    if (unlikely(CswriteCloseNull(&pvar_css, pvar_cswritep))) {
      goto OxBgenToPgen_ret_WRITE_FAIL;
    }

    reterr = SpgwFinish(&spgw);
    if (unlikely(reterr)) {
      goto OxBgenToPgen_ret_1;
    }
    putc_unlocked('\r', stdout);
    char* write_iter = strcpya_k(g_logbuf, "--bgen: ");
    write_iter = u32toa(variant_ct, write_iter);
    write_iter = strcpya_k(write_iter, " variant");
    if (variant_ct != 1) {
      *write_iter++ = 's';
    }
    write_iter = strcpya_k(write_iter, " remaining");
    if (load_filter_log_import_flags) {
      write_iter = strcpya_k(write_iter, " (");
      if (raw_variant_ct - multiallelic_tmp_skip_ct != variant_ct) {
        write_iter = u32toa(raw_variant_ct - multiallelic_tmp_skip_ct - variant_ct, write_iter);
        write_iter = strcpya_k(write_iter, " excluded by ");
        AppendLoadFilterFlagnames(load_filter_log_import_flags, &write_iter);
        *write_iter++ = ')';
      } else {
        AppendLoadFilterFlagnames(load_filter_log_import_flags, &write_iter);
        write_iter = strcpya_k(write_iter, " had no effect)");
      }
    }
    write_iter = strcpya_k(write_iter, "; ");
    const uint32_t outname_base_slen = outname_end - outname;
    write_iter = memcpya(write_iter, outname, outname_base_slen + 5);
    write_iter = strcpya_k(write_iter, " + ");
    write_iter = memcpya(write_iter, outname, outname_base_slen);
    write_iter = strcpya_k(write_iter, ".pvar");
    if (output_zst) {
      write_iter = strcpya_k(write_iter, ".zst");
    }
    write_iter = strcpya_k(write_iter, " written");
    if (!dosage_exists) {
      write_iter = strcpya_k(write_iter, " (only unphased hardcalls)");
    }
    snprintf(write_iter, kLogbufSize - 2 * kPglFnamesize - 512, ".\n");
    WordWrapB(0);
    logputsb();
    if (print_splitpar_warning) {
      logerrputs("Warning: Human chrX pseudoautosomal variant(s) appear to be present in the\ninput .bgen.  You probably want to include --split-par in your next command.\n");
    }
  }
  while (0) {
  OxBgenToPgen_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  OxBgenToPgen_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  OxBgenToPgen_ret_READ_FAIL:
    if (feof_unlocked(bgenfile)) {
      errno = 0;
    }
    putc_unlocked('\n', stdout);
    logerrprintfww(kErrprintfFread, bgenname, rstrerror(errno));
    reterr = kPglRetReadFail;
    break;
  OxBgenToPgen_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  OxBgenToPgen_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  OxBgenToPgen_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
  OxBgenToPgen_ret_INCONSISTENT_INPUT_2:
    logerrputsb();
  OxBgenToPgen_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  OxBgenToPgen_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  OxBgenToPgen_ret_DEGENERATE_DATA:
    reterr = kPglRetDegenerateData;
    break;
  OxBgenToPgen_ret_SLOPPY_CHRX_1:
    putc_unlocked('\n', stdout);
    logerrputs("Error: chrX is present in the input file, but no sex information was provided;\nrerun this import with --sample, --update-sex, or --impute-sex.  --split-par\nmay also be appropriate.\n");
    reterr = kPglRetInconsistentInput;
    break;
  OxBgenToPgen_ret_SLOPPY_CHRX_2:
    putc_unlocked('\n', stdout);
    logerrputs("Error: Human chrX pseudoautosomal variant(s) appear to be present in the\ninput .bgen, but --split-par was not specified.\n");
    reterr = kPglRetInconsistentInput;
    break;
  OxBgenToPgen_ret_bgen13_thread_fail:
    if (reterr == kPglRetInconsistentInput) {
      // --polyploid-mode note doesn't help here
      putc_unlocked('\n', stdout);
      logerrputs("Error: Polyploid genotype in .bgen file.\n");
    } else if (reterr == kPglRetMalformedInput) {
    OxBgenToPgen_ret_bgen11_thread_fail:
      putc_unlocked('\n', stdout);
      logerrputs("Error: Invalid compressed SNP block in .bgen file.\n");
    } else if (reterr == kPglRetNotYetSupported) {
      putc_unlocked('\n', stdout);
      logerrputs("Error: BGEN import doesn't currently support multiallelic variants, 29-32 bit\nprobability precision, or ploidy > 2.\n");
    }
    // note that nomem is also possible here
  }
 OxBgenToPgen_ret_1:
  if (common.libdeflate_decompressors) {
    for (uint32_t tidx = 0; tidx != max_thread_ct; ++tidx) {
      if (!common.libdeflate_decompressors[tidx]) {
        break;
      }
      libdeflate_free_decompressor(common.libdeflate_decompressors[tidx]);
    }
    // common.libdeflate_decompressors = nullptr;
  }
  CleanupSpgw(&spgw, &reterr);
  CleanupThreads(&tg);
  fclose_cond(bgenfile);
  fclose_cond(psamfile);
  CswriteCloseCond(&pvar_css, pvar_cswritep);
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

PglErr OxHapslegendToPgen(const char* hapsname, const char* legendname, const char* samplename, const char* const_fid, const char* ox_single_chr_str, const char* ox_missing_code, const char* missing_catname, MiscFlags misc_flags, ImportFlags import_flags, LoadFilterLogFlags load_filter_log_import_flags, OxfordImportFlags oxford_import_flags, uint32_t psam_01, uint32_t is_update_or_impute_sex, uint32_t is_splitpar, uint32_t is_sortvars, char id_delim, uint32_t max_thread_ct, char* outname, char* outname_end, ChrInfo* cip, uint32_t* pgi_generated_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* psamfile = nullptr;
  uintptr_t line_idx_haps = 0;
  uintptr_t line_idx_legend = 0;
  PglErr reterr = kPglRetSuccess;
  char* pvar_cswritep = nullptr;
  CompressStreamState pvar_css;
  PreinitCstream(&pvar_css);
  TextStream haps_txs;
  TextStream legend_txs;
  STPgenWriter spgw;
  PreinitTextStream(&haps_txs);
  PreinitTextStream(&legend_txs);
  PreinitSpgw(&spgw);
  {
    uint32_t sfile_sample_ct = 0;
    if (samplename[0]) {
      reterr = OxSampleToPsam(samplename, const_fid, ox_missing_code, missing_catname, misc_flags, import_flags, psam_01, id_delim, outname, outname_end, &sfile_sample_ct);
      if (unlikely(reterr)) {
        goto OxHapslegendToPgen_ret_1;
      }
      if (unlikely(sfile_sample_ct > (kMaxLongLine / 4))) {
        logerrputs("Error: Too many samples for .haps file converter.\n");
        reterr = kPglRetNotYetSupported;
        goto OxHapslegendToPgen_ret_1;
      }
    }

    const uint32_t output_zst = (import_flags / kfImportKeepAutoconvVzs) & 1;
    uint32_t max_line_blen;
    if (StandardizeMaxLineBlen(bigstack_left() / (4 + output_zst), &max_line_blen)) {
      goto OxHapslegendToPgen_ret_NOMEM;
    }
    reterr = InitTextStream(hapsname, max_line_blen, ClipU32(max_thread_ct - 1, 1, 4), &haps_txs);
    if (unlikely(reterr)) {
      if (reterr == kPglRetOpenFail) {
        const uint32_t slen = strlen(hapsname);
        if ((!StrEndsWith(hapsname, ".haps", slen)) &&
            (!StrEndsWith(hapsname, ".haps.gz", slen))) {
          logerrprintfww("Error: Failed to open %s : %s. (--haps expects a complete filename; did you forget '.haps' at the end?)\n", hapsname, strerror(errno));
          goto OxHapslegendToPgen_ret_1;
        }
      }
      goto OxHapslegendToPgen_ret_TSTREAM_FAIL_HAPS;
    }
    snprintf(outname_end, kMaxOutfnameExtBlen, ".pvar");
    if (output_zst) {
      snprintf(&(outname_end[5]), kMaxOutfnameExtBlen - 5, ".zst");
    }
    const uintptr_t overflow_buf_size = kCompressStreamBlock + max_line_blen;
    reterr = InitCstreamAlloc(outname, 0, output_zst, 1, overflow_buf_size, &pvar_css, &pvar_cswritep);
    if (unlikely(reterr)) {
      goto OxHapslegendToPgen_ret_1;
    }
    pvar_cswritep = strcpya_k(pvar_cswritep, "#CHROM\tPOS\tID\tREF\tALT" EOLN_STR);
    ++line_idx_haps;
    char* haps_line_iter = TextGet(&haps_txs);
    if (unlikely(!haps_line_iter)) {
      if (!TextStreamErrcode2(&haps_txs, &reterr)) {
        snprintf(g_logbuf, kLogbufSize, "Error: %s is empty.\n", hapsname);
        goto OxHapslegendToPgen_ret_INCONSISTENT_INPUT_WW;
      }
      goto OxHapslegendToPgen_ret_TSTREAM_FAIL_HAPS;
    }
    const uint32_t token_ct = CountTokens(haps_line_iter);
    FinalizeChrset(load_filter_log_import_flags, cip);
    const char* single_chr_str = nullptr;
    char* legend_line_start = nullptr;
    uint32_t single_chr_slen = 0;
    uint32_t cur_chr_code = 0;
    uint32_t is_haploid = 0;
    // Count number of samples in the .haps file.
    // Note that the user is not required to provide a .sample file.  So, if a
    // .sample file is provided, verify the count is as expected; if not,
    // generate a .psam with the correct sample count once we know what it is.
    uint32_t sample_ct;
    if (legendname[0]) {
      if (unlikely(token_ct % 2)) {
        snprintf(g_logbuf, kLogbufSize, "Error: %s has an odd number of tokens in the first line. (With --haps + --legend, the .haps file is expected to have no header columns.)\n", hapsname);
        goto OxHapslegendToPgen_ret_MALFORMED_INPUT_WW;
      }
      sample_ct = token_ct / 2;
      // Other .legend-specific initialization.
      reterr = InitOxfordSingleChr(ox_single_chr_str, "--legend argument", &single_chr_str, &single_chr_slen, &cur_chr_code, cip);
      if (unlikely(reterr)) {
        goto OxHapslegendToPgen_ret_1;
      }
      is_haploid = IsSet(cip->haploid_mask, cur_chr_code);
      reterr = InitTextStream(legendname, max_line_blen, 1, &legend_txs);
      if (unlikely(reterr)) {
        if (reterr == kPglRetOpenFail) {
          const uint32_t slen = strlen(legendname);
          if (!StrEndsWith(legendname, ".legend", slen)) {
            logerrprintfww("Error: Failed to open %s : %s. (--legend expects a complete filename; did you forget '.legend' at the end?)\n", legendname, strerror(errno));
            goto OxHapslegendToPgen_ret_1;
          }
        }
        goto OxHapslegendToPgen_ret_TSTREAM_FAIL_LEGEND;
      }
      ++line_idx_legend;
      legend_line_start = TextGet(&legend_txs);
      if (unlikely(!legend_line_start)) {
        if (!TextStreamErrcode2(&legend_txs, &reterr)) {
          snprintf(g_logbuf, kLogbufSize, "Error: %s is empty.\n", legendname);
          goto OxHapslegendToPgen_ret_MALFORMED_INPUT_WW;
        }
        goto OxHapslegendToPgen_ret_TSTREAM_FAIL_LEGEND;
      }
      // require at least 4 columns, in ID/pos/A1/A2 order; header text is
      // permitted to vary.  tolerate and ignore extra columns.
      if (unlikely(!NextTokenMult(legend_line_start, 3))) {
        goto OxHapslegendToPgen_ret_MISSING_TOKENS_LEGEND;
      }
    } else {
      // !legendname[0]
      if (unlikely((token_ct < 7) || (!(token_ct % 2)))) {
        snprintf(g_logbuf, kLogbufSize, "Error: Unexpected token count in line %" PRIuPTR " of %s (should be odd, >5).\n", line_idx_haps, hapsname);
        goto OxHapslegendToPgen_ret_MALFORMED_INPUT_WW;
      }
      sample_ct = (token_ct - 5) / 2;
    }
    if (sfile_sample_ct) {
      if (unlikely(sfile_sample_ct != sample_ct)) {
        snprintf(g_logbuf, kLogbufSize, "Error: .sample file has %u sample%s, while %s has %u.\n", sfile_sample_ct, (sfile_sample_ct == 1)? "" : "s", hapsname, sample_ct);
        goto OxHapslegendToPgen_ret_INCONSISTENT_INPUT_WW;
      }
    } else {
      // create a dummy .psam file with "per0", "per1", etc. IDs, matching
      // --dummy
      unsigned char* bigstack_mark2 = g_bigstack_base;
      char* writebuf = S_CAST(char*, bigstack_alloc_raw(kMaxMediumLine + max_line_blen));
      char* writebuf_flush = &(writebuf[kMaxMediumLine]);
      snprintf(outname_end, kMaxOutfnameExtBlen, ".psam");
      if (unlikely(fopen_checked(outname, FOPEN_WB, &psamfile))) {
        goto OxHapslegendToPgen_ret_OPEN_FAIL;
      }
      char* write_iter = strcpya_k(writebuf, "#IID\tSEX" EOLN_STR);
      for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        write_iter = strcpya_k(write_iter, "per");
        write_iter = u32toa(sample_idx, write_iter);
        write_iter = strcpya_k(write_iter, "\tNA" EOLN_STR);
        if (unlikely(fwrite_ck(writebuf_flush, psamfile, &write_iter))) {
          goto OxHapslegendToPgen_ret_WRITE_FAIL;
        }
      }
      if (unlikely(fclose_flush_null(writebuf_flush, write_iter, &psamfile))) {
        goto OxHapslegendToPgen_ret_WRITE_FAIL;
      }
      BigstackReset(bigstack_mark2);
    }

    const uint32_t variant_ct_limit = MINV(kPglMaxVariantCt, bigstack_left() / 8);
    const uint32_t keep_pgi = !(import_flags & kfImportKeepAutoconv);
    snprintf(outname_end, kMaxOutfnameExtBlen, ".pgen");
    if (keep_pgi) {
      *pgi_generated_ptr = 1;
    }
    uintptr_t spgw_alloc_cacheline_ct;
    uint32_t max_vrec_len;
    reterr = SpgwInitPhase1(outname, nullptr, nullptr, variant_ct_limit, sample_ct, 0, keep_pgi? kPgenWriteSeparateIndex : kPgenWriteAndCopy, kfPgenGlobalHardcallPhasePresent, (oxford_import_flags & (kfOxfordImportRefFirst | kfOxfordImportRefLast))? 1 : 2, &spgw, &spgw_alloc_cacheline_ct, &max_vrec_len);
    if (unlikely(reterr)) {
      if (reterr == kPglRetOpenFail) {
        logerrprintfww(kErrprintfFopen, outname, strerror(errno));
      }
      goto OxHapslegendToPgen_ret_1;
    }
    unsigned char* spgw_alloc;
    if (unlikely(bigstack_alloc_uc(spgw_alloc_cacheline_ct * kCacheline, &spgw_alloc))) {
      goto OxHapslegendToPgen_ret_NOMEM;
    }
    SpgwInitPhase2(max_vrec_len, &spgw, spgw_alloc);
    const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
    const uint32_t sample_ctl2_m1 = sample_ctl2 - 1;
    const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
    const uint32_t prohibit_extra_chr = (misc_flags / kfMiscProhibitExtraChr) & 1;
    const uint32_t prov_ref_allele_second = !(oxford_import_flags & kfOxfordImportRefFirst);
    const uint32_t phaseinfo_match_4char = prov_ref_allele_second? 0x20312030 : 0x20302031;
    const uint32_t phaseinfo_match = 1 + prov_ref_allele_second;
    uint32_t variant_ct = 0;
    uintptr_t variant_skip_ct = 0;
    const uint32_t sex_info_avail = sfile_sample_ct || is_update_or_impute_sex;
    uint32_t par_warn_code = UINT32_MAX;
    uint32_t print_splitpar_warning = 0;
    if (!(import_flags & kfImportLaxChrX)) {
      const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
      if ((!IsI32Neg(x_code)) && ((!sex_info_avail) || (IsHumanChrset(cip) && (!is_splitpar)))) {
        par_warn_code = x_code;
      }
    }
    uintptr_t* genovec;
    Halfword* phaseinfo_hwarr;
    if (unlikely(bigstack_alloc_w(sample_ctl2, &genovec) ||
                 bigstack_alloc_hw(2 * sample_ctl, &phaseinfo_hwarr))) {
      goto OxHapslegendToPgen_ret_NOMEM;
    }
    goto OxHapslegendToPgen_first_line;
    while (1) {
      ++line_idx_haps;
      reterr = TextGetUnsafe(&haps_txs, &haps_line_iter);
      if (reterr) {
        if (unlikely(reterr != kPglRetEof)) {
          goto OxHapslegendToPgen_ret_TSTREAM_FAIL_HAPS;
        }
        if (legend_line_start) {
          if (unlikely(TextGet(&legend_txs) != nullptr)) {
            snprintf(g_logbuf, kLogbufSize, "Error: %s has fewer nonheader lines than %s.\n", hapsname, legendname);
            goto OxHapslegendToPgen_ret_INCONSISTENT_INPUT_WW;
          }
        }
        reterr = kPglRetSuccess;
        break;
      }
    OxHapslegendToPgen_first_line:
      if (legend_line_start) {
        ++line_idx_legend;
        legend_line_start = TextGet(&legend_txs);
        if (unlikely(legend_line_start == nullptr)) {
          snprintf(g_logbuf, kLogbufSize, "Error: %s has fewer nonheader lines than %s.\n", legendname, hapsname);
          goto OxHapslegendToPgen_ret_INCONSISTENT_INPUT_WW;
        }
      }
      char* linebuf_iter = haps_line_iter;
      char* variant_id_start;
      char* variant_id_end;
      char* bp_start;
      char* bp_end;
      char* allele0_start;
      char* allele0_end;
      char* allele1_start;
      char* allele1_end;
      if (!legend_line_start) {
        char* chr_code_end = CurTokenEnd(haps_line_iter);
        variant_id_start = FirstNonTspace(chr_code_end);
        variant_id_end = FirstSpaceOrEoln(variant_id_start);
        bp_start = FirstNonTspace(variant_id_end);
        bp_end = FirstSpaceOrEoln(bp_start);
        allele0_start = FirstNonTspace(bp_end);
        allele0_end = FirstSpaceOrEoln(allele0_start);
        allele1_start = FirstNonTspace(allele0_end);
        allele1_end = FirstSpaceOrEoln(allele1_start);
        linebuf_iter = FirstNonTspace(allele1_end);
        if (unlikely(IsEolnKns(*linebuf_iter))) {
          goto OxHapslegendToPgen_ret_MISSING_TOKENS_HAPS;
        }
        reterr = GetOrAddChrCodeDestructive("--haps file", line_idx_haps, prohibit_extra_chr, haps_line_iter, chr_code_end, cip, &cur_chr_code);
        if (unlikely(reterr)) {
          goto OxHapslegendToPgen_ret_1;
        }
        if (!IsSet(cip->chr_mask, cur_chr_code)) {
          ++variant_skip_ct;
          haps_line_iter = AdvPastDelim(linebuf_iter, '\n');
          continue;
        }
        pvar_cswritep = chrtoa(cip, cur_chr_code, pvar_cswritep);
        is_haploid = IsSet(cip->haploid_mask, cur_chr_code);
      } else {
        variant_id_start = legend_line_start;
        variant_id_end = FirstSpaceOrEoln(variant_id_start);
        bp_start = FirstNonTspace(variant_id_end);
        bp_end = FirstSpaceOrEoln(bp_start);
        allele0_start = FirstNonTspace(bp_end);
        allele0_end = FirstSpaceOrEoln(allele0_start);
        allele1_start = FirstNonTspace(allele0_end);
        if (unlikely(IsEolnKns(*allele1_start))) {
          goto OxHapslegendToPgen_ret_MISSING_TOKENS_LEGEND;
        }
        allele1_end = FirstSpaceOrEoln(allele1_start);

        pvar_cswritep = memcpya(pvar_cswritep, single_chr_str, single_chr_slen);
      }
      *pvar_cswritep++ = '\t';
      const uint32_t id_slen = variant_id_end - variant_id_start;
      if (unlikely(id_slen > kMaxIdSlen)) {
        putc_unlocked('\n', stdout);
        logerrputs("Error: Variant names are limited to " MAX_ID_SLEN_STR " characters.\n");
        goto OxHapslegendToPgen_ret_MALFORMED_INPUT;
      }
      uint32_t cur_bp;
      if (unlikely(ScanUintDefcap(bp_start, &cur_bp))) {
        putc_unlocked('\n', stdout);
        if (legend_line_start) {
          logprintfww("Error: Invalid bp coordinate on line %" PRIuPTR " of %s.\n", line_idx_legend, legendname);
        } else {
          logprintfww("Error: Invalid bp coordinate on line %" PRIuPTR " of %s.\n", line_idx_haps, hapsname);
        }
        goto OxHapslegendToPgen_ret_MALFORMED_INPUT;
      }
      if (par_warn_code == cur_chr_code) {
        if (unlikely(!sex_info_avail)) {
          putc_unlocked('\n', stdout);
          logerrputs("Error: chrX is present in the input, but no sex information was provided; rerun\nthis import with --sample, --update-sex, or --impute-sex.  --split-par may also\nbe appropriate.\n");
          goto OxHapslegendToPgen_ret_INCONSISTENT_INPUT;
        }
        if ((cur_bp <= kPAR1IntersectionLast) || (cur_bp >= kPAR2IntersectionFirst)) {
          if (unlikely(!is_sortvars)) {
            putc_unlocked('\n', stdout);
            logerrputs("Error: Human chrX pseudoautosomal variant(s) appear to be present in the input,\nbut --split-par was not specified.\n");
            goto OxHapslegendToPgen_ret_INCONSISTENT_INPUT;
          }
          print_splitpar_warning = 1;
          par_warn_code = UINT32_MAX;
        }
      }
      pvar_cswritep = u32toa_x(cur_bp, '\t', pvar_cswritep);
      pvar_cswritep = memcpyax(pvar_cswritep, variant_id_start, id_slen, '\t');
      if (!prov_ref_allele_second) {
        pvar_cswritep = memcpyax(pvar_cswritep, allele0_start, allele0_end - allele0_start, '\t');
        pvar_cswritep = memcpya(pvar_cswritep, allele1_start, allele1_end - allele1_start);
      } else {
        pvar_cswritep = memcpyax(pvar_cswritep, allele1_start, allele1_end - allele1_start, '\t');
        pvar_cswritep = memcpya(pvar_cswritep, allele0_start, allele0_end - allele0_start);
      }
      AppendBinaryEoln(&pvar_cswritep);
      if (unlikely(Cswrite(&pvar_css, &pvar_cswritep))) {
        goto OxHapslegendToPgen_ret_WRITE_FAIL;
      }
      uintptr_t genovec_word_or = 0;
      uint32_t inner_loop_last = kBitsPerWordD2 - 1;
      // optimize common case: autosomal diploid, always exactly one space
      // this loop is time-critical; all my attemps to merge in the haploid
      // case have caused >10% slowdowns
      if ((!is_haploid) && IsEoln(linebuf_iter[sample_ct * 4 - 1])) {
        haps_line_iter = AdvPastDelim(&(linebuf_iter[sample_ct * 4 - 1]), '\n');
        linebuf_iter[sample_ct * 4 - 1] = ' ';
#ifdef USE_SSE2
        const VecU16 all0 = vecu16_set1(0x2030);
        const VecU16 all1 = vecu16_set1(0x2031);
        const uint32_t fullword_ct = sample_ct / kBitsPerWordD2;
        for (uint32_t widx = 0; widx != fullword_ct; ++widx) {
          uintptr_t geno_first = 0;
          for (uint32_t uii = 0; uii != 2; ++uii) {
            VecU16 cur_chars = vecu16_loadu(linebuf_iter);
            linebuf_iter += kBytesPerVec;
            // todo: better ARM implementation
            uintptr_t zero_mm = vecu16_movemask(cur_chars == all0);
            uintptr_t one_mm = vecu16_movemask(cur_chars == all1);
            cur_chars = vecu16_loadu(linebuf_iter);
            linebuf_iter += kBytesPerVec;
            zero_mm |= S_CAST(uintptr_t, vecu16_movemask(cur_chars == all0)) << kBytesPerVec;
            one_mm |= S_CAST(uintptr_t, vecu16_movemask(cur_chars == all1)) << kBytesPerVec;
#  ifndef USE_AVX2
            cur_chars = vecu16_loadu(linebuf_iter);
            linebuf_iter += kBytesPerVec;
            zero_mm |= S_CAST(uintptr_t, vecu16_movemask(cur_chars == all0)) << 32;
            one_mm |= S_CAST(uintptr_t, vecu16_movemask(cur_chars == all1)) << 32;
            cur_chars = vecu16_loadu(linebuf_iter);
            linebuf_iter += kBytesPerVec;
            zero_mm |= S_CAST(uintptr_t, vecu16_movemask(cur_chars == all0)) << 48;
            one_mm |= S_CAST(uintptr_t, vecu16_movemask(cur_chars == all1)) << 48;
#  endif
            if (unlikely(~(zero_mm | one_mm))) {
              // todo: other error messages
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid token on line %" PRIuPTR " of %s.\n", line_idx_haps, hapsname);
              goto OxHapslegendToPgen_ret_MALFORMED_INPUT_WW;
            }
            geno_first |= S_CAST(uintptr_t, PackWordToHalfwordMask5555(one_mm)) << (uii * 32);
          }
          const uintptr_t geno_second = (geno_first >> 1) & kMask5555;
          geno_first &= kMask5555;
          const uintptr_t geno_sum = geno_first + geno_second;
          genovec[widx] = geno_sum;
          genovec_word_or |= geno_sum;
          uintptr_t phaseinfo_hw = geno_second & (~geno_first);
          if (!prov_ref_allele_second) {
            phaseinfo_hw = geno_first & (~geno_second);
          }
          phaseinfo_hw = PackWordToHalfword(phaseinfo_hw);
          phaseinfo_hwarr[widx] = phaseinfo_hw;
        }
        const uint32_t remainder = sample_ct % kBitsPerWordD2;
        if (remainder) {
          const unsigned char* linebuf_iter_uc = R_CAST(const unsigned char*, linebuf_iter);
          uintptr_t genovec_word = 0;
          Halfword phaseinfo_hw = 0;
          for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != remainder; ++sample_idx_lowbits) {
            uint32_t cur_hap_4char;
            CopyFromUnalignedIncrU32(&cur_hap_4char, &linebuf_iter_uc);
            if (unlikely((cur_hap_4char & 0xfffefffeU) != 0x20302030)) {
              // todo: other error messages
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid token on line %" PRIuPTR " of %s.\n", line_idx_haps, hapsname);
              goto OxHapslegendToPgen_ret_MALFORMED_INPUT_WW;
            }
            // bugfix (28 Apr 2019): this needs to be uintptr_t for the next
            // left-shift to work
            const uintptr_t new_geno = (cur_hap_4char + (cur_hap_4char >> 16)) & 3;
            genovec_word |= new_geno << (2 * sample_idx_lowbits);
            if (cur_hap_4char == phaseinfo_match_4char) {
              phaseinfo_hw |= 1U << sample_idx_lowbits;
            }
          }
          genovec[fullword_ct] = genovec_word;
          genovec_word_or |= genovec_word;
          phaseinfo_hwarr[fullword_ct] = phaseinfo_hw;
        }
#else  // !USE_SSE2
        const unsigned char* linebuf_iter_uc = R_CAST(const unsigned char*, linebuf_iter);
        for (uint32_t widx = 0; ; ++widx) {
          if (widx >= sample_ctl2_m1) {
            if (widx > sample_ctl2_m1) {
              break;
            }
            inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
          }
          uintptr_t genovec_word = 0;
          uint32_t phaseinfo_hw = 0;
          for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
            // assumes little-endian
            uint32_t cur_hap_4char;
            CopyFromUnalignedIncrU32(&cur_hap_4char, &linebuf_iter_uc);
            if (unlikely((cur_hap_4char & 0xfffefffeU) != 0x20302030)) {
              if ((cur_hap_4char & 0xfffffffeU) == 0x202d2030) {
                // "0 - ", "1 - "
                goto OxHapslegendToPgen_ret_HAPLOID_TOKEN;
              }
              // any character < 32?
              if ((((cur_hap_4char & 0xe0e0e0e0U) * 7) & 0x80808080U) != 0x80808080U) {
                goto OxHapslegendToPgen_ret_MISSING_TOKENS_HAPS;
              }
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid token on line %" PRIuPTR " of %s.\n", line_idx_haps, hapsname);
              goto OxHapslegendToPgen_ret_MALFORMED_INPUT_WW;
            }
            const uintptr_t new_geno = (cur_hap_4char + (cur_hap_4char >> 16)) & 3;
            genovec_word |= new_geno << (2 * sample_idx_lowbits);
            if (cur_hap_4char == phaseinfo_match_4char) {
              phaseinfo_hw |= 1U << sample_idx_lowbits;
            }
          }
          genovec[widx] = genovec_word;
          genovec_word_or |= genovec_word;
          phaseinfo_hwarr[widx] = phaseinfo_hw;
        }
#endif  // !USE_SSE2
      } else {
        for (uint32_t widx = 0; ; ++widx) {
          if (widx >= sample_ctl2_m1) {
            if (widx > sample_ctl2_m1) {
              break;
            }
            inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
          }
          uintptr_t genovec_word = 0;
          uint32_t phaseinfo_hw = 0;
          for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
            const uint32_t first_hap_char_code = ctou32(*linebuf_iter);
            const uint32_t first_hap_int = first_hap_char_code - 48;
            char* post_first_hap = &(linebuf_iter[1]);
            if (unlikely((first_hap_int >= 2) || (*post_first_hap != ' '))) {
              if ((first_hap_char_code <= 32) || (ctou32(*post_first_hap) < 32)) {
                goto OxHapslegendToPgen_ret_MISSING_TOKENS_HAPS;
              }
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid token on line %" PRIuPTR " of %s.\n", line_idx_haps, hapsname);
              goto OxHapslegendToPgen_ret_MALFORMED_INPUT_WW;
            }
            char* second_hap = post_first_hap;
            uint32_t second_hap_char_code;
            do {
              second_hap_char_code = ctou32(*(++second_hap));
            } while (second_hap_char_code == 32);
            char* post_second_hap = &(second_hap[1]);
            const uint32_t post_second_hap_char_code = ctou32(*post_second_hap);
            uint32_t second_hap_int = second_hap_char_code - 48;
            if ((second_hap_int >= 2) || (post_second_hap_char_code > 32)) {
              if (likely(is_haploid && (second_hap_char_code == 45))) {
                // could require --sample, and require this sample to be male
                // in this case?
                second_hap_int = first_hap_int;
              } else {
                if (second_hap_char_code <= 32) {
                  goto OxHapslegendToPgen_ret_MISSING_TOKENS_HAPS;
                }
                if ((second_hap_char_code == 45) && (post_second_hap_char_code <= 32)) {
                  goto OxHapslegendToPgen_ret_HAPLOID_TOKEN;
                }
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid token on line %" PRIuPTR " of %s.\n", line_idx_haps, hapsname);
                goto OxHapslegendToPgen_ret_MALFORMED_INPUT_WW;
              }
            }
            genovec_word |= S_CAST(uintptr_t, first_hap_int + second_hap_int) << (2 * sample_idx_lowbits);
            if (first_hap_int + 2 * second_hap_int == phaseinfo_match) {
              phaseinfo_hw |= 1U << sample_idx_lowbits;
            }
            linebuf_iter = FirstNonChar(post_second_hap, ' ');
          }
          genovec[widx] = genovec_word;
          genovec_word_or |= genovec_word;
          phaseinfo_hwarr[widx] = phaseinfo_hw;
        }
        haps_line_iter = AdvPastDelim(linebuf_iter, '\n');
      }
      if (prov_ref_allele_second) {
        GenovecInvertUnsafe(sample_ct, genovec);
        ZeroTrailingNyps(sample_ct, genovec);
      }
      if (genovec_word_or & kMask5555) {
        uintptr_t* __attribute__((may_alias)) phaseinfo_warr = R_CAST(uintptr_t*, phaseinfo_hwarr);
        if (unlikely(SpgwAppendBiallelicGenovecHphase(genovec, nullptr, phaseinfo_warr, &spgw))) {
          goto OxHapslegendToPgen_ret_WRITE_FAIL;
        }
      } else {
        if (unlikely(SpgwAppendBiallelicGenovec(genovec, &spgw))) {
          goto OxHapslegendToPgen_ret_WRITE_FAIL;
        }
      }
      if (!(++variant_ct % 1000)) {
        printf("\r--haps%s: %uk variants converted.", legendname[0]? " + --legend" : "", variant_ct / 1000);
        fflush(stdout);
      }
    }
    if (unlikely(CswriteCloseNull(&pvar_css, pvar_cswritep))) {
      goto OxHapslegendToPgen_ret_WRITE_FAIL;
    }
    // bugfix (25 Oct 2025): this check needs to happen before SpgwFinish()
    if (unlikely(!variant_ct)) {
      char* write_iter = strcpya_k(g_logbuf, "Error: All ");
      write_iter = wtoa(variant_skip_ct, write_iter);
      write_iter = strcpya_k(write_iter, " variant");
      if (variant_skip_ct != 1) {
        *write_iter++ = 's';
      }
      write_iter = strcpya_k(write_iter, " in ");
      write_iter = strcpya(write_iter, hapsname);
      write_iter = strcpya_k(write_iter, " excluded by ");
      AppendLoadFilterFlagnames(load_filter_log_import_flags, &write_iter);
      strcpy_k(write_iter, ".\n");
      goto OxHapslegendToPgen_ret_INCONSISTENT_INPUT_WW;
    }
    reterr = SpgwFinish(&spgw);
    if (unlikely(reterr)) {
      PgenWriteFinishErrPrint(reterr, outname, outname_end);
      goto OxHapslegendToPgen_ret_1;
    }
    putc_unlocked('\r', stdout);
    char* write_iter = strcpya_k(g_logbuf, "--haps");
    if (legend_line_start) {
      write_iter = strcpya_k(write_iter, " + --legend");
    }
    write_iter = strcpya_k(write_iter, ": ");
    write_iter = u32toa(variant_ct, write_iter);
    write_iter = strcpya_k(write_iter, " variant");
    if (variant_ct != 1) {
      *write_iter++ = 's';
    }
    if (variant_skip_ct) {
      write_iter = strcpya_k(write_iter, " remaining (");
      write_iter = wtoa(variant_ct + variant_skip_ct, write_iter);
      write_iter = strcpya_k(write_iter, " processed, ");
      write_iter = wtoa(variant_skip_ct, write_iter);
      write_iter = strcpya_k(write_iter, " excluded by ");
      AppendLoadFilterFlagnames(load_filter_log_import_flags, &write_iter);
      *write_iter++ = ')';
    } else {
      write_iter = strcpya_k(write_iter, " processed");
      if (load_filter_log_import_flags) {
        write_iter = strcpya_k(write_iter, " (");
        AppendLoadFilterFlagnames(load_filter_log_import_flags, &write_iter);
        write_iter = strcpya_k(write_iter, " had no effect)");
      }
    }
    write_iter = strcpya_k(write_iter, "; ");
    const uint32_t outname_base_slen = outname_end - outname;
    write_iter = memcpya(write_iter, outname, outname_base_slen);
    write_iter = strcpya_k(write_iter, ".pgen + ");
    if (!sfile_sample_ct) {
      write_iter = memcpya(write_iter, outname, outname_base_slen);
      write_iter = strcpya_k(write_iter, ".psam + ");
    }
    write_iter = memcpya(write_iter, outname, outname_base_slen);
    write_iter = strcpya_k(write_iter, ".pvar");
    if (output_zst) {
      write_iter = strcpya_k(write_iter, ".zst");
    }
    snprintf(write_iter, kLogbufSize - 3 * kPglFnamesize - 64, " written.\n");
    WordWrapB(0);
    logputsb();
    if (print_splitpar_warning) {
      logerrputs("Warning: Human chrX pseudoautosomal variant(s) appear to be present in the\ninput.  You probably want to include --split-par in your next command.\n");
    }
  }
  while (0) {
  OxHapslegendToPgen_ret_TSTREAM_FAIL_HAPS:
    putc_unlocked('\n', stdout);
    TextStreamErrPrint(hapsname, &haps_txs);
    break;
  OxHapslegendToPgen_ret_TSTREAM_FAIL_LEGEND:
    TextStreamErrPrint(legendname, &legend_txs);
    break;
  OxHapslegendToPgen_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  OxHapslegendToPgen_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  OxHapslegendToPgen_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  OxHapslegendToPgen_ret_MISSING_TOKENS_HAPS:
    snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx_haps, hapsname);
  OxHapslegendToPgen_ret_MALFORMED_INPUT_WW:
    putc_unlocked('\n', stdout);
    WordWrapB(0);
    logerrputsb();
  OxHapslegendToPgen_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  OxHapslegendToPgen_ret_MISSING_TOKENS_LEGEND:
    putc_unlocked('\n', stdout);
    snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx_legend, legendname);
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetMalformedInput;
    break;
  OxHapslegendToPgen_ret_HAPLOID_TOKEN:
    putc_unlocked('\n', stdout);
    snprintf(g_logbuf, kLogbufSize, "Error: Haploid-only token on line %" PRIuPTR " of %s.\n", line_idx_haps, hapsname);
    WordWrapB(0);
    logerrputsb();
    reterr = legendname[0]? kPglRetInconsistentInput : kPglRetMalformedInput;
    break;
  OxHapslegendToPgen_ret_INCONSISTENT_INPUT_WW:
    putc_unlocked('\n', stdout);
    WordWrapB(0);
    logerrputsb();
  OxHapslegendToPgen_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 OxHapslegendToPgen_ret_1:
  CleanupSpgw(&spgw, &reterr);
  CleanupTextStream2(legendname, &legend_txs, &reterr);
  CleanupTextStream2(hapsname, &haps_txs, &reterr);
  CswriteCloseCond(&pvar_css, pvar_cswritep);
  fclose_cond(psamfile);
  BigstackReset(bigstack_mark);
  return reterr;
}


// probable todo: use the VCF/BGEN-import parallelization strategy to speed
// this up.
static_assert(sizeof(Dosage) == 2, "Plink1DosageToPgen() needs to be updated.");
PglErr Plink1DosageToPgen(const char* dosagename, const char* famname, const char* mapname, const char* import_single_chr_str, const Plink1DosageInfo* pdip, const char* missing_catname, MiscFlags misc_flags, ImportFlags import_flags, LoadFilterLogFlags load_filter_log_import_flags, uint32_t psam_01, FamCol fam_cols, int32_t missing_pheno, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, double import_dosage_certainty, uint32_t max_thread_ct, char* outname, char* outname_end, ChrInfo* cip) {
  // Tried making this single-pass, doesn't appear to be an improvement.
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;

  // these are not allocated on bigstack, and must be explicitly freed
  PhenoCol* pheno_cols = nullptr;
  char* pheno_names = nullptr;
  uint32_t pheno_ct = 0;

  FILE* psamfile = nullptr;
  uintptr_t line_idx = 0;
  Plink1DosageFlags flags = pdip->flags;
  PglErr reterr = kPglRetSuccess;
  char* pvar_cswritep = nullptr;
  CompressStreamState pvar_css;
  PreinitCstream(&pvar_css);
  textFILE dosage_txf;
  TextStream dosage_txs;
  STPgenWriter spgw;
  PreinitTextFile(&dosage_txf);
  PreinitTextStream(&dosage_txs);
  PreinitSpgw(&spgw);
  {
    // 1. Read .fam file.  (May as well support most .psam files too, since
    //    it's the same driver function.  However, unless 'noheader' modifier
    //    is present, SID field cannot be used for disambiguation.)
    uint32_t raw_sample_ct = 0;
    uintptr_t* sample_include = nullptr;
    PedigreeIdInfo pii;
    InitPedigreeIdInfo(misc_flags, &pii);
    uintptr_t* sex_nm = nullptr;
    uintptr_t* sex_male = nullptr;
    uintptr_t* founder_info = nullptr;
    uintptr_t max_pheno_name_blen = 0;
    reterr = LoadPsam(famname, nullptr, missing_catname, fam_cols, 0x7fffffff, missing_pheno, (misc_flags / kfMiscAffection01) & 1, (misc_flags / kfMiscNoCategorical) & 1, (misc_flags / kfMiscNeg9PhenoReallyMissing) & 1, max_thread_ct, &pii, &sample_include, &founder_info, &sex_nm, &sex_male, &pheno_cols, &pheno_names, &raw_sample_ct, &pheno_ct, &max_pheno_name_blen);
    if (unlikely(reterr)) {
      goto Plink1DosageToPgen_ret_1;
    }

    // 2. Read dosage-file header line if it exists, then write new .psam.
    const uint32_t first_data_col_idx = pdip->skips[0] + pdip->skips[1] + pdip->skips[2] + 3;
    uint32_t sample_ct = 0;
    uint32_t* dosage_sample_idx_to_fam_uidx;
    if (unlikely(bigstack_end_alloc_u32(raw_sample_ct, &dosage_sample_idx_to_fam_uidx))) {
      goto Plink1DosageToPgen_ret_NOMEM;
    }
    if (flags & kfPlink1DosageNoheader) {
      sample_ct = raw_sample_ct;
      for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        dosage_sample_idx_to_fam_uidx[sample_idx] = sample_idx;
      }
    } else {
      ZeroWArr(BitCtToWordCt(raw_sample_ct), sample_include);
      const uint32_t tmp_htable_size = GetHtableFastSize(raw_sample_ct);
      uint32_t* htable_tmp;
      char* idbuf;
      if (unlikely(bigstack_end_alloc_u32(tmp_htable_size, &htable_tmp) ||
                   bigstack_end_alloc_c(pii.sii.max_sample_id_blen, &idbuf))) {
        goto Plink1DosageToPgen_ret_NOMEM;
      }
      const uint32_t duplicate_idx = PopulateStrboxHtable(pii.sii.sample_ids, raw_sample_ct, pii.sii.max_sample_id_blen, tmp_htable_size, htable_tmp);
      if (unlikely(duplicate_idx)) {
        char* duplicate_sample_id = &(pii.sii.sample_ids[duplicate_idx * pii.sii.max_sample_id_blen]);
        char* duplicate_fid_end = AdvToDelim(duplicate_sample_id, '\t');
        *duplicate_fid_end = ' ';
        snprintf(g_logbuf, kLogbufSize, "Error: Duplicate sample ID \"%s\" in .fam file.\n", duplicate_sample_id);
        goto Plink1DosageToPgen_ret_MALFORMED_INPUT_WW;
      }

      uintptr_t dst_capacity = bigstack_left();
#ifdef __LP64__
      if (dst_capacity > S_CAST(uintptr_t, kMaxLongLine) + S_CAST(uintptr_t, kDecompressChunkSize)) {
        dst_capacity = S_CAST(uintptr_t, kMaxLongLine) + S_CAST(uintptr_t, kDecompressChunkSize);
      } else if (unlikely(dst_capacity < kDecompressMinCapacity)) {
        goto Plink1DosageToPgen_ret_NOMEM;
      }
#else
      if (unlikely(dst_capacity < kDecompressMinCapacity)) {
        goto Plink1DosageToPgen_ret_NOMEM;
      }
#endif
      // not formally allocated
      char* dst = R_CAST(char*, g_bigstack_base);
      reterr = TextFileOpenEx(dosagename, kMaxLongLine, dst_capacity, dst, &dosage_txf);
      if (unlikely(reterr)) {
        goto Plink1DosageToPgen_ret_TFILE_FAIL2;
      }
      line_idx = 1;
      char* line_start = TextFileGet(&dosage_txf);
      if (unlikely(!line_start)) {
        reterr = TextFileRawErrcode(&dosage_txf);
        goto Plink1DosageToPgen_ret_TFILE_FAIL2;
      }
      char* loadbuf_iter = NextTokenMult(line_start, first_data_col_idx);
      if (unlikely(!loadbuf_iter)) {
        goto Plink1DosageToPgen_ret_MISSING_TOKENS;
      }
      const char id_delim = pdip->id_delim;
      do {
        char* fid_end = CurTokenEnd(loadbuf_iter);
        char* iid_start;
        char* iid_end;
        uint32_t iid_slen;
        if (id_delim) {
          iid_end = fid_end;
          fid_end = S_CAST(char*, memchr(loadbuf_iter, ctou32(id_delim), iid_end - loadbuf_iter));
          if (unlikely(!fid_end)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Sample ID in --import-dosage file does not contain '%c' delimiter.\n", id_delim);
            goto Plink1DosageToPgen_ret_INCONSISTENT_INPUT_2;
          }
          iid_start = &(fid_end[1]);
          iid_slen = iid_end - iid_start;
          if (unlikely(memchr(iid_start, ctou32(id_delim), iid_slen))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Sample ID in --import-dosage file has multiple instances of '%c'.\n", id_delim);
            goto Plink1DosageToPgen_ret_INCONSISTENT_INPUT_2;
          }
        } else {
          iid_start = FirstNonTspace(fid_end);
          if (unlikely(IsEolnKns(*iid_start))) {
            goto Plink1DosageToPgen_ret_MISSING_TOKENS;
          }
          iid_end = CurTokenEnd(iid_start);
          iid_slen = iid_end - iid_start;
        }
        const uint32_t fid_slen = fid_end - loadbuf_iter;
        const uint32_t cur_id_slen = fid_slen + iid_slen + 1;
        if (unlikely(cur_id_slen >= pii.sii.max_sample_id_blen)) {
          logerrputs("Error: .fam file does not contain all sample IDs in dosage file.\n");
          goto Plink1DosageToPgen_ret_INCONSISTENT_INPUT;
        }
        char* idbuf_iid = memcpyax(idbuf, loadbuf_iter, fid_slen, '\t');
        memcpyx(idbuf_iid, iid_start, iid_slen, '\0');
        uint32_t sample_uidx = StrboxHtableFind(idbuf, pii.sii.sample_ids, htable_tmp, pii.sii.max_sample_id_blen, cur_id_slen, tmp_htable_size);
        if (unlikely(sample_uidx == UINT32_MAX)) {
          logerrputs("Error: .fam file does not contain all sample IDs in dosage file.\n");
          goto Plink1DosageToPgen_ret_INCONSISTENT_INPUT;
        }
        if (unlikely(IsSet(sample_include, sample_uidx))) {
          idbuf_iid[-1] = ' ';
          snprintf(g_logbuf, kLogbufSize, "Error: Duplicate sample ID \"%s\" in dosage file.\n", idbuf);
          goto Plink1DosageToPgen_ret_MALFORMED_INPUT_WW;
        }
        SetBit(sample_uidx, sample_include);
        dosage_sample_idx_to_fam_uidx[sample_ct++] = sample_uidx;
        loadbuf_iter = FirstNonTspace(iid_end);
      } while (!IsEolnKns(*loadbuf_iter));
      // Not worth the trouble of move-constructing dosage_txs from dosage_txf,
      // since we may need to load the .map in between, and we'll need to
      // rewind again anyway.
      if (unlikely(CleanupTextFile2(dosagename, &dosage_txf, &reterr))) {
        goto Plink1DosageToPgen_ret_1;
      }
      line_idx = 0;
    }

    snprintf(outname_end, kMaxOutfnameExtBlen, ".psam");
    if (unlikely(fopen_checked(outname, FOPEN_WB, &psamfile))) {
      goto Plink1DosageToPgen_ret_OPEN_FAIL;
    }
    char* writebuf = g_textbuf;
    char* writebuf_flush = &(writebuf[kMaxMediumLine]);
    char* write_iter = writebuf;
    *write_iter++ = '#';
    const uint32_t write_fid = DataFidColIsRequired(sample_include, &pii.sii, sample_ct, 1);
    if (write_fid) {
      write_iter = strcpya_k(write_iter, "FID\t");
    }
    write_iter = strcpya_k(write_iter, "IID");
    const uint32_t write_sid = DataSidColIsRequired(sample_include, pii.sii.sids, sample_ct, pii.sii.max_sid_blen, 1);
    if (write_sid) {
      write_iter = strcpya_k(write_iter, "\tSID");
    }
    const uint32_t write_parents = DataParentalColsAreRequired(sample_include, &pii.sii, &pii.parental_id_info, sample_ct, 1);
    if (write_parents) {
      write_iter = strcpya_k(write_iter, "\tPAT\tMAT");
    }
    write_iter = strcpya_k(write_iter, "\tSEX");
    for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
      *write_iter++ = '\t';
      write_iter = strcpya(write_iter, &(pheno_names[pheno_idx * max_pheno_name_blen]));
      if (unlikely(fwrite_ck(writebuf_flush, psamfile, &write_iter))) {
        goto Plink1DosageToPgen_ret_WRITE_FAIL;
      }
    }
    AppendBinaryEoln(&write_iter);
    const char ctrl_char = '1' - psam_01;
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      const uint32_t sample_uidx = dosage_sample_idx_to_fam_uidx[sample_idx];
      const char* cur_sample_id = &(pii.sii.sample_ids[sample_uidx * pii.sii.max_sample_id_blen]);
      if (!write_fid) {
        cur_sample_id = AdvPastDelim(cur_sample_id, '\t');
      }
      write_iter = strcpya(write_iter, cur_sample_id);
      if (write_sid) {
        *write_iter++ = '\t';
        write_iter = strcpya(write_iter, &(pii.sii.sids[sample_uidx * pii.sii.max_sid_blen]));
      }
      if (write_parents) {
        *write_iter++ = '\t';
        write_iter = strcpyax(write_iter, &(pii.parental_id_info.paternal_ids[pii.parental_id_info.max_paternal_id_blen * sample_uidx]), '\t');
        write_iter = strcpya(write_iter, &(pii.parental_id_info.maternal_ids[pii.parental_id_info.max_maternal_id_blen * sample_uidx]));
      }
      *write_iter++ = '\t';
      if (IsSet(sex_nm, sample_uidx)) {
        *write_iter++ = '2' - IsSet(sex_male, sample_uidx);
      } else {
        write_iter = strcpya_k(write_iter, "NA");
      }
      for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
        if (unlikely(fwrite_ck(writebuf_flush, psamfile, &write_iter))) {
          goto Plink1DosageToPgen_ret_WRITE_FAIL;
        }
        *write_iter++ = '\t';
        write_iter = AppendPhenoStrEx(&(pheno_cols[pheno_idx]), "NA", 2, sample_uidx, ctrl_char, write_iter);
      }
      AppendBinaryEoln(&write_iter);
      if (unlikely(fwrite_ck(writebuf_flush, psamfile, &write_iter))) {
        goto Plink1DosageToPgen_ret_WRITE_FAIL;
      }
    }
    if (unlikely(fclose_flush_null(writebuf_flush, write_iter, &psamfile))) {
      goto Plink1DosageToPgen_ret_WRITE_FAIL;
    }
    // Don't need sample info any more.
    BigstackEndReset(bigstack_end_mark);

    // 3. Read .map file if it exists.
    uint32_t max_variant_id_slen = 1;
    ChrIdx* variant_chr_codes = nullptr;
    uint32_t* variant_bps = nullptr;
    char** variant_ids = nullptr;
    double* variant_cms = nullptr;
    uint32_t* variant_id_htable = nullptr;
    uintptr_t* variant_already_seen = nullptr;
    uint32_t variant_id_htable_size = 0;
    uint32_t map_variant_ct = 0;
    FinalizeChrset(load_filter_log_import_flags, cip);
    if (mapname) {
      reterr = LoadMap(mapname, misc_flags, load_filter_log_import_flags, cip, &max_variant_id_slen, &variant_chr_codes, &variant_bps, &variant_ids, &variant_cms, &map_variant_ct);
      if (unlikely(reterr)) {
        goto Plink1DosageToPgen_ret_1;
      }
      const uint32_t map_variant_ctl = BitCtToWordCt(map_variant_ct);
      if (unlikely(bigstack_alloc_w(map_variant_ctl, &variant_already_seen))) {
        goto Plink1DosageToPgen_ret_NOMEM;
      }
      SetAllBits(map_variant_ct, variant_already_seen);
      unsigned char* bigstack_end_mark2 = g_bigstack_end;

      // allow hash table to only use half of available memory
      g_bigstack_end = &(g_bigstack_base[RoundDownPow2(bigstack_left() / 2, kEndAllocAlign)]);

      reterr = AllocAndPopulateIdHtableMt(variant_already_seen, TO_CONSTCPCONSTP(variant_ids), map_variant_ct, bigstack_left() / 2, max_thread_ct, &variant_id_htable, nullptr, &variant_id_htable_size, nullptr);
      if (unlikely(reterr)) {
        goto Plink1DosageToPgen_ret_1;
      }
      g_bigstack_end = bigstack_end_mark2;
      ZeroWArr(map_variant_ctl, variant_already_seen);
    }

    // 4. Dosage file pass 1: count variants, check whether any decimal dosages
    //    need to be saved, write .pvar.
    //
    // Lots of overlap with OxGenToPgen().

    const uint32_t output_zst = (import_flags / kfImportKeepAutoconvVzs) & 1;
    // This is overly conservative in the .pvar.zst output case, but it
    // shouldn't matter.
    uintptr_t ulii = bigstack_left() >> (1 + output_zst);
    // Main write buffer needs to be a bit larger than the corresponding read
    // buffer; also need a small chromosome buffer.
    if (unlikely(ulii < RoundUpPow2((kMaxMediumLine + kCacheline) / 2, kCacheline))) {
      goto Plink1DosageToPgen_ret_NOMEM;
    }
    ulii -= RoundUpPow2((kMaxMediumLine + kCacheline) / 2, kCacheline);
    uint32_t max_line_blen;
    if (unlikely(StandardizeMaxLineBlen(ulii, &max_line_blen))) {
      goto Plink1DosageToPgen_ret_NOMEM;
    }
    reterr = InitTextStream(dosagename, max_line_blen, ClipU32(max_thread_ct - 1, 1, 3), &dosage_txs);
    if (unlikely(reterr)) {
      goto Plink1DosageToPgen_ret_TSTREAM_REWIND1_FAIL;
    }
    snprintf(outname_end, kMaxOutfnameExtBlen, ".pvar");
    if (output_zst) {
      snprintf(&(outname_end[5]), kMaxOutfnameExtBlen - 5, ".zst");
    }
    unsigned char* bigstack_mark2 = g_bigstack_base;
    const uintptr_t overflow_buf_size = kCompressStreamBlock + kDecompressChunkSize + max_line_blen;
    reterr = InitCstreamAlloc(outname, 0, output_zst, 1, overflow_buf_size, &pvar_css, &pvar_cswritep);
    if (unlikely(reterr)) {
      goto Plink1DosageToPgen_ret_1;
    }
    const uint32_t prohibit_extra_chr = (misc_flags / kfMiscProhibitExtraChr) & 1;
    const char* single_chr_str = nullptr;
    uint32_t single_chr_slen = 0;
    const uint32_t chr_col_idx = pdip->chr_col_idx;
    const uint32_t check_chr_col = (chr_col_idx != UINT32_MAX);
    if (!check_chr_col) {
      if (import_single_chr_str) {
        uint32_t chr_code_raw = GetChrCodeRaw(import_single_chr_str);
        if (chr_code_raw == UINT32_MAX) {
          // command-line parser guarantees that prohibit_extra_chr is false
          // here
          single_chr_str = import_single_chr_str;
          single_chr_slen = strlen(import_single_chr_str);
        } else {
          uint32_t chr_code = chr_code_raw;
          if (chr_code > cip->max_code) {
            if (unlikely(chr_code < kMaxContigs)) {
              logerrputs("Error: --import-dosage single-chr= code is not in the chromosome set.\n");
              goto Plink1DosageToPgen_ret_INVALID_CMDLINE;
            }
            chr_code = cip->xymt_codes[chr_code - kMaxContigs];
            if (unlikely(IsI32Neg(chr_code))) {
              logerrputs("Error: --import-dosage single-chr= code is not in the chromosome set.\n");
              goto Plink1DosageToPgen_ret_INVALID_CMDLINE;
            }
          }
          if (unlikely(!IsSet(cip->chr_mask, chr_code))) {
            logerrputs("Error: --import-dosage single-chr= code is excluded by chromosome filter.\n");
            goto Plink1DosageToPgen_ret_INVALID_CMDLINE;
          }
          // check this allocation?...
          char* chr_buf = S_CAST(char*, bigstack_alloc_raw(kCacheline));
          char* chr_name_end = chrtoa(cip, chr_code, chr_buf);
          single_chr_str = chr_buf;
          single_chr_slen = chr_name_end - chr_buf;
        }
      } else {
        // default to "chr0"
        if (unlikely(!IsSet(cip->chr_mask, 0))) {
          logerrputs("Error: No --import-dosage chromosome information specified, and chr0 excluded.\n");
          goto Plink1DosageToPgen_ret_INVALID_CMDLINE;
        }
        char* chr_buf = S_CAST(char*, bigstack_alloc_raw(kCacheline));
        char* chr_name_end = chrtoa(cip, 0, chr_buf);
        single_chr_str = chr_buf;
        single_chr_slen = chr_name_end - chr_buf;
      }
    }
    pvar_cswritep = strcpya_k(pvar_cswritep, "#CHROM\tPOS\tID\tREF\tALT");
    if (variant_cms) {
      pvar_cswritep = strcpya_k(pvar_cswritep, "\tCM");
    }
    AppendBinaryEoln(&pvar_cswritep);
    // types:
    // 0 = #CHROM
    // 1 = POS
    // 2 = ID
    // 3 = REF
    // 4 = ALT
    // 5 = first data column
    // (command-line parser verifies that CHROM/POS don't collide with
    // anything else)
    uint64_t parse_table[6];
    // high bits = col index, low bits = col type
    const uint32_t id_col_idx = pdip->skips[0];
    const uint32_t prov_ref_allele_second = !(flags & kfPlink1DosageRefFirst);
    uint32_t ref_col_idx = id_col_idx + pdip->skips[1] + 1;
    const uint32_t alt_col_idx = ref_col_idx + (!prov_ref_allele_second);
    ref_col_idx += prov_ref_allele_second;
    parse_table[0] = (S_CAST(uint64_t, id_col_idx) << 32) + 2;
    parse_table[1] = (S_CAST(uint64_t, ref_col_idx) << 32) + 3;
    parse_table[2] = (S_CAST(uint64_t, alt_col_idx) << 32) + 4;
    uint32_t relevant_initial_col_ct = 3;
    if (check_chr_col) {
      parse_table[relevant_initial_col_ct++] = S_CAST(uint64_t, chr_col_idx) << 32;
    }
    const uint32_t check_pos_col = (pdip->pos_col_idx != UINT32_MAX);
    if (check_pos_col) {
      parse_table[relevant_initial_col_ct++] = (S_CAST(uint64_t, pdip->pos_col_idx) << 32) + 1;
    }
#if (__GNUC__ >= 12) && (__GNUC__ <= 15)
    // https://github.com/cms-sw/cmssw/issues/44582
    // confirmed this fires for gcc 12, 14, and 15
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Warray-bounds"
#endif
    STD_SORT(relevant_initial_col_ct, u64cmp, parse_table);
#if (__GNUC__ >= 12) && (__GNUC__ <= 15)
#  pragma GCC diagnostic pop
#endif
    uint32_t col_skips[6];
    uint32_t col_types[6];
    for (uint32_t uii = 0; uii != relevant_initial_col_ct; ++uii) {
      const uint64_t parse_table_entry = parse_table[uii];
      col_skips[uii] = parse_table_entry >> 32;
      col_types[uii] = S_CAST(uint32_t, parse_table_entry);
    }
    col_skips[relevant_initial_col_ct] = first_data_col_idx;
    col_types[relevant_initial_col_ct++] = 5;
    for (uint32_t uii = relevant_initial_col_ct - 1; uii; --uii) {
      col_skips[uii] -= col_skips[uii - 1];
    }

    double dosage_multiplier = kDosageMid;
    double dosage_ceil = 32767.5 / 16384.0;
    if (flags & kfPlink1DosageFormatSingle01) {
      dosage_multiplier = kDosageMax;
      dosage_ceil = 32767.5 / 32768.0;
    }
    uint32_t format_triple = (flags / kfPlink1DosageFormatTriple) & 1;
    uint32_t format_infer = !(flags & (kfPlink1DosageFormatSingle | kfPlink1DosageFormatDouble | kfPlink1DosageFormatTriple));
    // dosage_erase_halfdist corresponds to very high confidence (keep only
    // hardcall), while force_missing_halfdist corresponds to very low
    // confidence.
    const uint32_t dosage_erase_halfdist = kDosage4th - dosage_erase_thresh;
    const uint32_t force_missing_halfdist_p1 = (import_dosage_certainty == 0.0)? 0 : (1 + S_CAST(uint32_t, (import_dosage_certainty - 0.5) * kDosageMid));
    uint32_t dosage_exists = 0;
    uint32_t variant_ct = 0;
    uintptr_t variant_skip_ct = 0;
    uint32_t variant_uidx = 0;
    char* line_iter = TextLineEnd(&dosage_txs);
    if (!(flags & kfPlink1DosageNoheader)) {
      // Skip header line.
      line_idx = 1;
      reterr = TextGetUnsafe(&dosage_txs, &line_iter);
      if (unlikely(reterr)) {
        goto Plink1DosageToPgen_ret_TSTREAM_REWIND1_FAIL;
      }
      line_iter = AdvPastDelim(line_iter, '\n');
    }
    ++line_idx;
    for (; TextGetUnsafe2(&dosage_txs, &line_iter); line_iter = AdvPastDelim(line_iter, '\n'), ++line_idx) {
      char* token_ptrs[6];
      uint32_t token_slens[6];
      line_iter = TokenLex0(line_iter, col_types, col_skips, relevant_initial_col_ct, token_ptrs, token_slens);
      if (unlikely(!line_iter)) {
        goto Plink1DosageToPgen_ret_MISSING_TOKENS;
      }
      // ID
      const char* variant_id = token_ptrs[2];
      const uint32_t variant_id_slen = token_slens[2];
      if (map_variant_ct) {
        variant_uidx = VariantIdDupflagHtableFind(variant_id, TO_CONSTCPCONSTP(variant_ids), variant_id_htable, variant_id_slen, variant_id_htable_size, max_variant_id_slen);
        if (variant_uidx >> 31) {
          if (likely(variant_uidx == UINT32_MAX)) {
            ++variant_skip_ct;
            continue;
          }
          snprintf(g_logbuf, kLogbufSize, "Error: Variant ID '%s' appears multiple times in .map file.\n", variant_ids[variant_uidx & 0x7fffffff]);
          goto Plink1DosageToPgen_ret_MALFORMED_INPUT_WW;
        }
        if (unlikely(IsSet(variant_already_seen, variant_uidx))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Variant ID '%s' appears multiple times in --import-dosage file.\n", variant_ids[variant_uidx]);
          goto Plink1DosageToPgen_ret_MALFORMED_INPUT_WW;
        }
        // already performed chromosome filtering
        pvar_cswritep = chrtoa(cip, variant_chr_codes[variant_uidx], pvar_cswritep);
        *pvar_cswritep++ = '\t';
        pvar_cswritep = u32toa_x(variant_bps[variant_uidx], '\t', pvar_cswritep);
        pvar_cswritep = memcpya(pvar_cswritep, variant_id, variant_id_slen);
      } else {
        if (unlikely(variant_id_slen > kMaxIdSlen)) {
          putc_unlocked('\n', stdout);
          logerrputs("Error: Variant names are limited to " MAX_ID_SLEN_STR " characters.\n");
          goto Plink1DosageToPgen_ret_MALFORMED_INPUT;
        }
        // #CHROM
        if (check_chr_col) {
          char* chr_code_str = token_ptrs[0];
          char* chr_code_end = &(chr_code_str[token_slens[0]]);
          uint32_t cur_chr_code;
          reterr = GetOrAddChrCodeDestructive("--import-dosage file", line_idx, prohibit_extra_chr, chr_code_str, chr_code_end, cip, &cur_chr_code);
          if (unlikely(reterr)) {
            goto Plink1DosageToPgen_ret_1;
          }
          if (!IsSet(cip->chr_mask, cur_chr_code)) {
            ++variant_skip_ct;
            continue;
          }
          pvar_cswritep = chrtoa(cip, cur_chr_code, pvar_cswritep);
        } else {
          pvar_cswritep = memcpya(pvar_cswritep, single_chr_str, single_chr_slen);
        }
        *pvar_cswritep++ = '\t';
        // POS
        if (check_pos_col) {
          const char* pos_str = token_ptrs[1];
          // no need to support negative values here
          uint32_t cur_bp;
          if (unlikely(ScanUintDefcap(pos_str, &cur_bp))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid bp coordinate on line %" PRIuPTR " of %s.\n", line_idx, dosagename);
            goto Plink1DosageToPgen_ret_MALFORMED_INPUT_WW;
          }
          pvar_cswritep = u32toa(cur_bp, pvar_cswritep);
        } else {
          *pvar_cswritep++ = '0';
        }
        *pvar_cswritep++ = '\t';
        pvar_cswritep = memcpya(pvar_cswritep, variant_id, variant_id_slen);
      }
      ++variant_ct;
      *pvar_cswritep++ = '\t';
      // REF, ALT
      pvar_cswritep = memcpyax(pvar_cswritep, token_ptrs[3], token_slens[3], '\t');
      pvar_cswritep = memcpya(pvar_cswritep, token_ptrs[4], token_slens[4]);
      if (variant_cms) {
        *pvar_cswritep++ = '\t';
        pvar_cswritep = dtoa_g_p8(variant_cms[variant_uidx], pvar_cswritep);
      }
      AppendBinaryEoln(&pvar_cswritep);
      if (unlikely(Cswrite(&pvar_css, &pvar_cswritep))) {
        goto Plink1DosageToPgen_ret_WRITE_FAIL;
      }
      if (!dosage_exists) {
        char* linebuf_iter = token_ptrs[5];
        if (format_infer) {
          const uint32_t remaining_col_ct = CountTokens(linebuf_iter);
          if (remaining_col_ct == sample_ct) {
            flags |= kfPlink1DosageFormatSingle;
          } else if (remaining_col_ct == sample_ct * 3) {
            format_triple = 1;
          } else if (unlikely(remaining_col_ct != sample_ct * 2)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Unexpected format=infer column count in --import-dosage file (%u; should be %u, %u, or %u).\n", remaining_col_ct, sample_ct, sample_ct * 2, sample_ct * 3);
            goto Plink1DosageToPgen_ret_INCONSISTENT_INPUT_WW;
          }
          format_infer = 0;
        }
        if (flags & kfPlink1DosageFormatSingle) {
          // todo: modify these loops to cleanly update line_iter
          for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
            if (unlikely(!linebuf_iter)) {
              goto Plink1DosageToPgen_ret_MISSING_TOKENS;
            }
            double a1_dosage;
            char* str_end = ScanadvDouble(linebuf_iter, &a1_dosage);
            if ((!str_end) || (a1_dosage < (0.5 / 32768.0)) || (a1_dosage >= dosage_ceil)) {
              linebuf_iter = NextToken(linebuf_iter);
              continue;
            }
            if (unlikely(!IsSpaceOrEoln(*str_end))) {
              str_end = CurTokenEnd(str_end);
              *str_end = '\0';
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid numeric token '%s' on line %" PRIuPTR " of --import-dosage file.\n", linebuf_iter, line_idx);
              goto Plink1DosageToPgen_ret_MALFORMED_INPUT_WW;
            }
            a1_dosage *= dosage_multiplier;
            const uint32_t dosage_int = S_CAST(int32_t, a1_dosage + 0.5);
            const uint32_t halfdist = BiallelicDosageHalfdist(dosage_int);
            if ((halfdist < dosage_erase_halfdist) && (halfdist >= force_missing_halfdist_p1)) {
              dosage_exists = 1;
              break;
            }
            linebuf_iter = FirstNonTspace(str_end);
          }
        } else {
          // for compatibility with plink 1.x, do not actually parse third
          // value of each triplet if format=3
          for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
            if (unlikely(!linebuf_iter)) {
              goto Plink1DosageToPgen_ret_MISSING_TOKENS;
            }
            double prob_2a1;
            char* str_end = ScanadvDouble(linebuf_iter, &prob_2a1);
            if (!str_end) {
              linebuf_iter = NextTokenMult(linebuf_iter, 2 + format_triple);
              continue;
            }
            if (unlikely(!IsSpaceOrEoln(*str_end))) {
              str_end = CurTokenEnd(str_end);
              *str_end = '\0';
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid numeric token '%s' on line %" PRIuPTR " of --import-dosage file.\n", linebuf_iter, line_idx);
              goto Plink1DosageToPgen_ret_MALFORMED_INPUT_WW;
            }
            linebuf_iter = FirstNonTspace(str_end);
            if (unlikely(IsEolnKns(*linebuf_iter))) {
              goto Plink1DosageToPgen_ret_MISSING_TOKENS;
            }
            double prob_1a1;
            str_end = ScanadvDouble(linebuf_iter, &prob_1a1);
            if (!str_end) {
              linebuf_iter = NextTokenMult(linebuf_iter, 1 + format_triple);
              continue;
            }
            if (unlikely(!IsSpaceOrEoln(*str_end))) {
              str_end = CurTokenEnd(str_end);
              *str_end = '\0';
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid numeric token '%s' on line %" PRIuPTR " of --import-dosage file.\n", linebuf_iter, line_idx);
              goto Plink1DosageToPgen_ret_MALFORMED_INPUT_WW;
            }
            linebuf_iter = NextTokenMult(str_end, 1 + format_triple);
            double prob_one_or_two_a1 = prob_2a1 + prob_1a1;
            if ((prob_2a1 < 0.0) || (prob_1a1 < 0.0) || (prob_one_or_two_a1 > 1.01 * (1 + kSmallEpsilon))) {
              continue;
            }
            if (prob_one_or_two_a1 > 1.0) {
              const double rescale = 1.0 / prob_one_or_two_a1;
              prob_2a1 *= rescale;
              prob_1a1 *= rescale;
              prob_one_or_two_a1 = 1.0;
            }
            const uint32_t dosage_int = S_CAST(int32_t, prob_2a1 * 32768 + prob_1a1 * 16384 + 0.5);
            const uint32_t halfdist = BiallelicDosageHalfdist(dosage_int);
            if ((halfdist < dosage_erase_halfdist) && ((prob_2a1 >= import_dosage_certainty) || (prob_1a1 >= import_dosage_certainty) || (prob_one_or_two_a1 <= 1.0 - import_dosage_certainty))) {
              dosage_exists = 1;
              break;
            }
          }
        }
      }
      if (!(variant_ct % 1000)) {
        printf("\r--import-dosage: %uk variants scanned.", variant_ct / 1000);
        fflush(stdout);
      }
    }
    if (TextStreamErrcode2(&dosage_txs, &reterr)) {
      goto Plink1DosageToPgen_ret_TSTREAM_REWIND1_FAIL;
    }
    putc_unlocked('\r', stdout);
    if (unlikely(CswriteCloseNull(&pvar_css, pvar_cswritep))) {
      goto Plink1DosageToPgen_ret_WRITE_FAIL;
    }
    if (unlikely(!variant_ct)) {
      if (!variant_skip_ct) {
        logerrputs("Error: Empty --import-dosage file.\n");
        goto Plink1DosageToPgen_ret_DEGENERATE_DATA;
      }
      char* log_write_iter = strcpya_k(g_logbuf, "Error: All ");
      log_write_iter = wtoa(variant_skip_ct, log_write_iter);
      log_write_iter = strcpya_k(log_write_iter, " variant");
      if (variant_skip_ct != 1) {
        *log_write_iter++ = 's';
      }
      log_write_iter = strcpya_k(log_write_iter, " in --import-dosage file excluded by ");
      AppendLoadFilterFlagnames(load_filter_log_import_flags, &log_write_iter);
      strcpy_k(log_write_iter, ".\n");
      goto Plink1DosageToPgen_ret_INCONSISTENT_INPUT_WW;
    }
    {
      char* log_write_iter = strcpya_k(g_logbuf, "--import-dosage: ");
      log_write_iter = wtoa(variant_ct + variant_skip_ct, log_write_iter);
      log_write_iter = strcpya_k(log_write_iter, " variant");
      if (variant_ct + variant_skip_ct != 1) {
        *log_write_iter++ = 's';
      }
      log_write_iter = strcpya_k(log_write_iter, " scanned");
      if (!dosage_exists) {
        log_write_iter = strcpya_k(log_write_iter, " (all hardcalls)");
      }
      if (variant_skip_ct) {
        log_write_iter = strcpya_k(log_write_iter, "; ");
        log_write_iter = wtoa(variant_skip_ct, log_write_iter);
        log_write_iter = strcpya_k(log_write_iter, " excluded by ");
        AppendLoadFilterFlagnames(load_filter_log_import_flags, &log_write_iter);
        log_write_iter = strcpya_k(log_write_iter, ", ");
        log_write_iter = u32toa(variant_ct, log_write_iter);
        log_write_iter = strcpya_k(log_write_iter, " remaining");
      } else if (load_filter_log_import_flags) {
        log_write_iter = strcpya_k(log_write_iter, " (");
        AppendLoadFilterFlagnames(load_filter_log_import_flags, &log_write_iter);
        log_write_iter = strcpya_k(log_write_iter, " had no effect)");
      }
      strcpy_k(log_write_iter, ".\n");
      WordWrapB(0);
      logputsb();
    }

    // 5. Dosage file pass 2: write .pgen.
    BigstackReset(bigstack_mark2);
    reterr = TextRewind(&dosage_txs);
    if (unlikely(reterr)) {
      goto Plink1DosageToPgen_ret_TSTREAM_REWIND1_FAIL;
    }
    line_iter = TextLineEnd(&dosage_txs);
    const uintptr_t line_ct = line_idx - 1;
    line_idx = 0;
    if (!(flags & kfPlink1DosageNoheader)) {
      // skip header line again
      ++line_idx;
      reterr = TextGetUnsafe(&dosage_txs, &line_iter);
      if (unlikely(reterr)) {
        goto Plink1DosageToPgen_ret_TSTREAM_REWIND2_FAIL;
      }
      line_iter = AdvPastDelim(line_iter, '\n');
    }
    snprintf(outname_end, kMaxOutfnameExtBlen, ".pgen");
    uintptr_t spgw_alloc_cacheline_ct;
    uint32_t max_vrec_len;
    reterr = SpgwInitPhase1(outname, nullptr, nullptr, variant_ct, sample_ct, 0, kPgenWriteBackwardSeek, dosage_exists? kfPgenGlobalDosagePresent: kfPgenGlobal0, (flags & (kfPlink1DosageRefFirst | kfPlink1DosageRefLast))? 1 : 2, &spgw, &spgw_alloc_cacheline_ct, &max_vrec_len);
    if (unlikely(reterr)) {
      if (reterr == kPglRetOpenFail) {
        logerrprintfww(kErrprintfFopen, outname, strerror(errno));
      }
      goto Plink1DosageToPgen_ret_1;
    }
    unsigned char* spgw_alloc;
    if (unlikely(bigstack_alloc_uc(spgw_alloc_cacheline_ct * kCacheline, &spgw_alloc))) {
      goto Plink1DosageToPgen_ret_NOMEM;
    }
    SpgwInitPhase2(max_vrec_len, &spgw, spgw_alloc);

    const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
    const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
    uintptr_t* genovec;
    Halfword* dosage_present_hwarr;
    if (unlikely(bigstack_alloc_w(sample_ctl2, &genovec) ||
                 bigstack_alloc_hw(2 * sample_ctl, &dosage_present_hwarr))) {
      goto Plink1DosageToPgen_ret_NOMEM;
    }
    Dosage* dosage_main = nullptr;
    if (dosage_exists) {
      if (unlikely(bigstack_alloc_dosage(sample_ct, &dosage_main))) {
        goto Plink1DosageToPgen_ret_NOMEM;
      }
    }
    if (hard_call_thresh == UINT32_MAX) {
      hard_call_thresh = kDosageMid / 10;
    }
    dosage_ceil = 2.02 * (1 + kSmallEpsilon);
    if (flags & kfPlink1DosageFormatSingle01) {
      dosage_ceil = 1.01 * (1 + kSmallEpsilon);
    }
    const uint32_t hard_call_halfdist = kDosage4th - hard_call_thresh;
    const uint32_t sample_ctl2_m1 = sample_ctl2 - 1;
    uint32_t vidx = 0;
    for (; line_idx < line_ct; line_iter = AdvPastDelim(line_iter, '\n')) {
      ++line_idx;
      reterr = TextGetUnsafe(&dosage_txs, &line_iter);
      if (unlikely(reterr)) {
        goto Plink1DosageToPgen_ret_TSTREAM_REWIND2_FAIL;
      }
      if (variant_skip_ct) {
        if (map_variant_ct) {
          char* variant_id = NextTokenMult0(line_iter, id_col_idx);
          const uint32_t variant_id_slen = strlen_se(variant_id);
          if (VariantIdDupflagHtableFind(variant_id, TO_CONSTCPCONSTP(variant_ids), variant_id_htable, variant_id_slen, variant_id_htable_size, max_variant_id_slen) == UINT32_MAX) {
            continue;
          }
          line_iter = NextTokenMult(variant_id, first_data_col_idx - id_col_idx);
        } else {
          char* chr_code_str = NextTokenMult0(line_iter, chr_col_idx);
          char* chr_code_end = CurTokenEnd(chr_code_str);
          line_iter = NextTokenMult(chr_code_end, first_data_col_idx - chr_col_idx);
          *chr_code_end = '\0';
          const uint32_t chr_code = GetChrCode(chr_code_str, cip, chr_code_end - chr_code_str);
          if (!IsSet(cip->chr_mask, chr_code)) {
            continue;
          }
        }
      } else {
        line_iter = NextTokenMult(line_iter, first_data_col_idx);
      }
      if (!line_iter) {
        goto Plink1DosageToPgen_ret_MISSING_TOKENS;
      }
      uint32_t inner_loop_last = kBitsPerWordD2 - 1;
      Dosage* dosage_main_iter = dosage_main;
      char* linebuf_iter = line_iter;
      for (uint32_t widx = 0; ; ++widx) {
        if (widx >= sample_ctl2_m1) {
          if (widx > sample_ctl2_m1) {
            break;
          }
          inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
        }
        uintptr_t genovec_word = 0;
        uint32_t dosage_present_hw = 0;
        if (flags & kfPlink1DosageFormatSingle) {
          for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
            if (unlikely(IsEolnKns(*linebuf_iter))) {
              goto Plink1DosageToPgen_ret_MISSING_TOKENS;
            }
            double a1_dosage;
            char* str_end = ScanadvDouble(linebuf_iter, &a1_dosage);
            if ((!str_end) || (a1_dosage < 0.0) || (a1_dosage > dosage_ceil)) {
              genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
              linebuf_iter = FirstNonTspace(CurTokenEnd(linebuf_iter));
              continue;
            }
            if (unlikely(!IsSpaceOrEoln(*str_end))) {
              str_end = CurTokenEnd(str_end);
              *str_end = '\0';
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid numeric token '%s' on line %" PRIuPTR " of --import-dosage file.\n", linebuf_iter, line_idx);
              goto Plink1DosageToPgen_ret_MALFORMED_INPUT_WW;
            }
            linebuf_iter = FirstNonTspace(str_end);
            uint32_t dosage_int = S_CAST(int32_t, a1_dosage * dosage_multiplier + 0.5);
            if (dosage_int > kDosageMax) {
              dosage_int = kDosageMax;
            }
            const uint32_t cur_halfdist = BiallelicDosageHalfdist(dosage_int);
            if (cur_halfdist < hard_call_halfdist) {
              genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
              if (cur_halfdist < force_missing_halfdist_p1) {
                continue;
              }
            } else {
              genovec_word |= ((dosage_int + (kDosage4th * k1LU)) / kDosageMid) << (2 * sample_idx_lowbits);
              if (cur_halfdist >= dosage_erase_halfdist) {
                continue;
              }
            }
            dosage_present_hw |= 1U << sample_idx_lowbits;
            *dosage_main_iter++ = dosage_int;
          }
        } else {
          for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
            if (unlikely(!linebuf_iter)) {
              goto Plink1DosageToPgen_ret_MISSING_TOKENS;
            }
            double prob_2a1;
            char* str_end = ScanadvDouble(linebuf_iter, &prob_2a1);
            if (!str_end) {
              genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
              linebuf_iter = NextTokenMult(linebuf_iter, 2 + format_triple);
              continue;
            }
            if (unlikely(!IsSpaceOrEoln(*str_end))) {
              str_end = CurTokenEnd(str_end);
              *str_end = '\0';
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid numeric token '%s' on line %" PRIuPTR " of --import-dosage file.\n", linebuf_iter, line_idx);
              goto Plink1DosageToPgen_ret_MALFORMED_INPUT_WW;
            }
            linebuf_iter = FirstNonTspace(str_end);
            if (unlikely(IsEolnKns(*linebuf_iter))) {
              goto Plink1DosageToPgen_ret_MISSING_TOKENS;
            }
            double prob_1a1;
            str_end = ScanadvDouble(linebuf_iter, &prob_1a1);
            if (!str_end) {
              genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
              linebuf_iter = NextTokenMult(linebuf_iter, 1 + format_triple);
              continue;
            }
            if (unlikely(!IsSpaceOrEoln(*str_end))) {
              str_end = CurTokenEnd(str_end);
              *str_end = '\0';
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid numeric token '%s' on line %" PRIuPTR " of --import-dosage file.\n", linebuf_iter, line_idx);
              goto Plink1DosageToPgen_ret_MALFORMED_INPUT_WW;
            }
            linebuf_iter = NextTokenMult(str_end, 1 + format_triple);
            double prob_one_or_two_a1 = prob_2a1 + prob_1a1;
            if ((prob_2a1 < 0.0) || (prob_1a1 < 0.0) || (prob_one_or_two_a1 > 1.01 * (1 + kSmallEpsilon))) {
              genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
              continue;
            }
            if (prob_one_or_two_a1 > 1.0) {
              const double rescale = 1.0 / prob_one_or_two_a1;
              prob_2a1 *= rescale;
              prob_1a1 *= rescale;
              prob_one_or_two_a1 = 1.0;
            }
            if ((prob_2a1 < import_dosage_certainty) && (prob_1a1 < import_dosage_certainty) && (prob_one_or_two_a1 > 1.0 - import_dosage_certainty)) {
              genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
            }
            const uint32_t dosage_int = S_CAST(int32_t, prob_2a1 * 32768 + prob_1a1 * 16384 + 0.5);
            const uint32_t cur_halfdist = BiallelicDosageHalfdist(dosage_int);
            if (cur_halfdist < hard_call_halfdist) {
              genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
            } else {
              genovec_word |= ((dosage_int + (kDosage4th * k1LU)) / kDosageMid) << (2 * sample_idx_lowbits);
              if (cur_halfdist >= dosage_erase_halfdist) {
                continue;
              }
            }
            dosage_present_hw |= 1U << sample_idx_lowbits;
            *dosage_main_iter++ = dosage_int;
          }
        }
        genovec[widx] = genovec_word;
        dosage_present_hwarr[widx] = dosage_present_hw;
      }
      if (!prov_ref_allele_second) {
        GenovecInvertUnsafe(sample_ct, genovec);
        ZeroTrailingNyps(sample_ct, genovec);
      }
      if (dosage_main_iter != dosage_main) {
        const uint32_t dosage_ct = dosage_main_iter - dosage_main;
        if (!prov_ref_allele_second) {
          BiallelicDosage16Invert(dosage_ct, dosage_main);
        }
        uintptr_t* __attribute__((may_alias)) dosage_present_warr = R_CAST(uintptr_t*, dosage_present_hwarr);
        reterr = SpgwAppendBiallelicGenovecDosage16(genovec, dosage_present_warr, dosage_main, dosage_ct, &spgw);
        if (unlikely(reterr)) {
          goto Plink1DosageToPgen_ret_1;
        }
      } else {
        if (unlikely(SpgwAppendBiallelicGenovec(genovec, &spgw))) {
          goto Plink1DosageToPgen_ret_WRITE_FAIL;
        }
      }
      ++vidx;
      if (!(vidx % 1000)) {
        printf("\r--import-dosage: %uk variants converted.", vidx / 1000);
        fflush(stdout);
      }
    }
    reterr = SpgwFinish(&spgw);
    if (unlikely(reterr)) {
      goto Plink1DosageToPgen_ret_1;
    }
    putc_unlocked('\r', stdout);
    write_iter = strcpya_k(g_logbuf, "--import-dosage: ");
    const uint32_t outname_base_slen = outname_end - outname;
    write_iter = memcpya(write_iter, outname, outname_base_slen + 5);
    write_iter = strcpya_k(write_iter, " + ");
    write_iter = memcpya(write_iter, outname, outname_base_slen);
    write_iter = strcpya_k(write_iter, ".pvar");
    if (output_zst) {
      write_iter = strcpya_k(write_iter, ".zst");
    }
    write_iter = strcpya_k(write_iter, " + ");
    write_iter = memcpya(write_iter, outname, outname_base_slen);
    strcpy_k(write_iter, ".psam written.\n");
    WordWrapB(0);
    logputsb();
  }
  while (0) {
  Plink1DosageToPgen_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  Plink1DosageToPgen_ret_TFILE_FAIL2:
    if (TextFileEof(&dosage_txf)) {
      logerrprintfww("Error: %s is empty.\n", dosagename);
      reterr = kPglRetInconsistentInput;
    } else {
      TextFileErrPrint(dosagename, &dosage_txf);
    }
    break;
  Plink1DosageToPgen_ret_TSTREAM_REWIND1_FAIL:
    if (!(flags & kfPlink1DosageNoheader)) {
    Plink1DosageToPgen_ret_TSTREAM_REWIND2_FAIL:
      if ((reterr == kPglRetOpenFail) || (reterr == kPglRetEof)) {
        reterr = kPglRetRewindFail;
      }
    }
    TextStreamErrPrintRewind(dosagename, &dosage_txs, &reterr);
    break;
  Plink1DosageToPgen_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  Plink1DosageToPgen_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  Plink1DosageToPgen_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  Plink1DosageToPgen_ret_MISSING_TOKENS:
    snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, dosagename);
  Plink1DosageToPgen_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    putc_unlocked('\n', stdout);
    logerrputsb();
  Plink1DosageToPgen_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  Plink1DosageToPgen_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
  Plink1DosageToPgen_ret_INCONSISTENT_INPUT_2:
    logerrputsb();
  Plink1DosageToPgen_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  Plink1DosageToPgen_ret_DEGENERATE_DATA:
    reterr = kPglRetDegenerateData;
    break;
  }
 Plink1DosageToPgen_ret_1:
  CleanupSpgw(&spgw, &reterr);
  fclose_cond(psamfile);
  CswriteCloseCond(&pvar_css, pvar_cswritep);
  CleanupTextFile2(dosagename, &dosage_txf, &reterr);
  CleanupTextStream2(dosagename, &dosage_txs, &reterr);
  free_cond(pheno_names);
  CleanupPhenoCols(pheno_ct, pheno_cols);
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}


typedef struct GenerateDummyCtxStruct {
  uint32_t sample_ct;
  // binary search over cdf is faster than (int)(log(drand)/log(q)) for
  // truncated geometric distribution
  STD_ARRAY_DECL(uint64_t, kBitsPerWordD2, phase_geomdist);
  STD_ARRAY_DECL(uint64_t, kBitsPerWordD2, dosage_geomdist);
  uint32_t phase_invert;
  uint32_t dosage_geomdist_max;
  uint32_t hard_call_halfdist;
  uint32_t dosage_erase_halfdist;

  uint32_t geno3_thresh_ct;
  // Usually (<nonmissing geno freq> * 2^32).
  // If <nonmissing geno freq> is in (0, 1e-7), it's rounded up to 1e-7 so that
  // discretization doesn't really interfere with Hardy-Weinberg equilibrium.
  uint64_t* geno3_thresh_arr;

  sfmt_t** sfmtp_arr;

  uint32_t cur_block_write_ct;

  uintptr_t* write_genovecs[2];
  uintptr_t* write_phasepresents[2];
  uintptr_t* write_phaseinfos[2];
  uint32_t* write_dosage_cts[2];
  uintptr_t* write_dosage_presents[2];
  Dosage* write_dosage_mains[2];
  uint32_t* write_dphase_cts[2];
  uintptr_t* write_dphase_presents[2];
  SDosage* write_dphase_deltas[2];
} GenerateDummyCtx;

static_assert(sizeof(Dosage) == 2, "GenerateDummyThread() needs to be updated.");
THREAD_FUNC_DECL GenerateDummyThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  GenerateDummyCtx* ctx = S_CAST(GenerateDummyCtx*, arg->sharedp->context);

  const uint32_t sample_ct = ctx->sample_ct;
  const uint32_t calc_thread_ct = GetThreadCt(arg->sharedp);
  STD_ARRAY_KREF(uint64_t, kBitsPerWordD2) phase_geomdist = ctx->phase_geomdist;
  STD_ARRAY_KREF(uint64_t, kBitsPerWordD2) dosage_geomdist = ctx->dosage_geomdist;
  const uint32_t geno3_thresh_ct = ctx->geno3_thresh_ct;
  const uint64_t* geno3_thresh_arr = ctx->geno3_thresh_arr;
  const uint32_t phase_invert = ctx->phase_invert;
  const uint32_t phase_is_variable = (phase_geomdist[kBitsPerWordD2 - 1] != 0);
  const uint32_t phase_exists = phase_invert || phase_is_variable;
  const uint32_t dosage_geomdist_max = ctx->dosage_geomdist_max;
  const uint32_t dosage_exists = (dosage_geomdist_max != kBitsPerWord);
  const uint32_t hard_call_halfdist = ctx->hard_call_halfdist;
  const uint32_t dosage_erase_halfdist = ctx->dosage_erase_halfdist;
  const uintptr_t sample_ctaw2 = NypCtToAlignedWordCt(sample_ct);
  const uint32_t sample_ctl2_m1 = (sample_ct - 1) / kBitsPerWordD2;
  const uintptr_t sample_ctaw = BitCtToAlignedWordCt(sample_ct);
  sfmt_t* sfmtp = ctx->sfmtp_arr[tidx];
  uint64_t u64rand = sfmt_genrand_uint64(sfmtp);
  uint32_t rand16_left = 4;
  uint32_t ld_denom = 0;
  uint32_t ld_invert = 0;
  double afreq = sfmt_to_res53(sfmt_genrand_uint64(sfmtp));
  uint64_t geno3_thresh = 1LLU << 32;
  double geno3_thresh_d = 4294967296.0;
  if (geno3_thresh_ct) {
    uint32_t uii = 0;
    if (geno3_thresh_ct > 1) {
      // assumes rand16_left != 0
      uii = ((u64rand & 65535) * geno3_thresh_ct) >> 16;
      --rand16_left;
    }
    geno3_thresh = geno3_thresh_arr[uii];
    geno3_thresh_d = u63tod(geno3_thresh);
  }
  uint32_t parity = 0;
  do {
    const uintptr_t cur_block_write_ct = ctx->cur_block_write_ct;
    uint32_t vidx = (tidx * cur_block_write_ct) / calc_thread_ct;
    const uint32_t vidx_end = ((tidx + 1) * cur_block_write_ct) / calc_thread_ct;
    uintptr_t* write_genovec_iter = &(ctx->write_genovecs[parity][vidx * sample_ctaw2]);
    Halfword* write_phasepresent_iter = DowncastWToHW(&(ctx->write_phasepresents[parity][vidx * sample_ctaw]));
    Halfword* write_phaseinfo_iter = DowncastWToHW(&(ctx->write_phaseinfos[parity][vidx * sample_ctaw]));
    uint32_t* write_dosage_ct_iter = &(ctx->write_dosage_cts[parity][vidx]);
    Halfword* write_dosage_present_iter = DowncastWToHW(&(ctx->write_dosage_presents[parity][vidx * sample_ctaw]));
    Dosage* write_dosage_main_iter = &(ctx->write_dosage_mains[parity][vidx * sample_ct]);
    uint32_t* write_dphase_ct_iter = &(ctx->write_dphase_cts[parity][vidx]);
    Halfword* write_dphase_present_iter = DowncastWToHW(&(ctx->write_dphase_presents[parity][vidx * sample_ctaw]));
    SDosage* write_dphase_delta_iter = &(ctx->write_dphase_deltas[parity][vidx * sample_ct]);
    uintptr_t* prev_genovec = nullptr;
    for (; vidx != vidx_end; ++vidx) {
      Dosage* cur_dosage_main_iter = write_dosage_main_iter;
      SDosage* cur_dphase_delta_iter = write_dphase_delta_iter;
      uint64_t geno1_thresh;
      uint64_t geno2_thresh;
      {
        const uint64_t cur_rand64 = sfmt_genrand_uint64(sfmtp);
        // ~50% chance of high LD
        if (cur_rand64 & 1) {
          prev_genovec = nullptr;
          // Only update allele and missing frequencies if not simulating LD.

          // Could make this follow e.g. a beta(0.2, 0.2) distribution, but
          // that's unimportant for software testing, so I won't bother
          // unless/until that inverse-CDF is already implemented in this
          // codebase for some reason.
          afreq = sfmt_to_res53(cur_rand64);
          if (geno3_thresh_ct > 1) {
            if (!rand16_left) {
              u64rand = sfmt_genrand_uint64(sfmtp);
              rand16_left = 4;
            }
            const uint32_t uii = ((u64rand & 65535) * geno3_thresh_ct) >> 16;
            u64rand >>= 16;
            --rand16_left;
            geno3_thresh = geno3_thresh_arr[uii];
            geno3_thresh_d = u63tod(geno3_thresh);
          }
        } else {
          // When generating a genovec_word in high LD with the previous one,
          // we use the following process to select positions to generate a
          // new genotype:
          // 1. Draw a uint16.
          // 2. Multiply by ld_denom, and right-shift 16 (Lemire's fast
          //    alternative to the modulo reduction).
          // 3. If this is >= loop_len, stop selecting.  Otherwise, generate a
          //    new genotype for this value of sample_idx_lowbits.
          // Our chance of continuing the loop is r=(loop_len/ld_denom) on each
          // iteration, so the expected number of iterations is
          //   r/(1-r) = loop_len / (ld_denom - loop_len).
          // Given a target genotype-reselection rate of R, we want this last
          // expectation to be near R*loop_len.  (Technically a bit larger,
          // since we're sampling with replacement, but we only need this
          // calculation to hit the right ballpark.)  Solving this equation:
          //   R * loop_len = loop_len / (ld_denom - loop_len)
          //   R = 1 / (ld_denom - loop_len)
          //   ld_denom - loop_len = 1 / R
          //   ld_denom = loop_len + (1 / R)
          // Since the chance of regenerating the same genotype as before is
          // sometimes smaller than 1/2 but never by that much, choosing R =
          // (2 / kPglMaxDifflistLenDivisor) should test a mix of
          // LD-compression and not-quite-LD-compression cases, while smaller
          // values of R test the former more heavily.
          const uint32_t rand_shift = (cur_rand64 >> 1) & 3;
          ld_denom = kBitsPerWordD2 + ((kPglMaxDifflistLenDivisor / 2) << rand_shift);
          // 1/16 chance of LD-invert.
          ld_invert = ((cur_rand64 & 120) == 120);
        }
        const double inv_afreq = 1.0 - afreq;
        // respect Hardy-Weinberg equilibrium
        geno1_thresh = S_CAST(int64_t, geno3_thresh_d * afreq * afreq);
        geno2_thresh = S_CAST(int64_t, geno3_thresh_d * (1.0 - inv_afreq * inv_afreq));
      }
      uint32_t loop_len = kBitsPerWordD2;
      uint32_t loop_len_d2 = kBitsPerWordD4;
      for (uint32_t widx = 0; ; ++widx) {
        if (widx >= sample_ctl2_m1) {
          if (widx > sample_ctl2_m1) {
            break;
          }
          loop_len = ModNz(sample_ct, kBitsPerWordD2);
          loop_len_d2 = loop_len / 2;
          if (prev_genovec) {
            ld_denom = ld_denom + loop_len - kBitsPerWordD2;
          }
        }
        uintptr_t genovec_word;
        if (!prev_genovec) {
          genovec_word = 0;
          for (uint32_t rand_idx = 0; rand_idx != loop_len_d2; ++rand_idx) {

            // sfmt_genrand_uint64 calls can't be mixed with
            // sfmt_genrand_uint32 calls, so use it here even in 32-bit build
            const uint64_t cur_rand = sfmt_genrand_uint64(sfmtp);
            const uintptr_t rand_lowbits = cur_rand & UINT32_MAX;
            uintptr_t cur_geno;
            if (rand_lowbits >= geno2_thresh) {
              cur_geno = 2 + (rand_lowbits >= geno3_thresh);
            } else {
              cur_geno = (rand_lowbits >= geno1_thresh);
            }
            genovec_word = genovec_word << 2;
            genovec_word |= cur_geno;
            const uintptr_t rand_highbits = cur_rand >> 32;
            if (rand_highbits >= geno2_thresh) {
              cur_geno = 2 + (rand_highbits >= geno3_thresh);
            } else {
              cur_geno = (rand_highbits >= geno1_thresh);
            }
            genovec_word = genovec_word << 2;
            genovec_word |= cur_geno;
          }
          if (loop_len % 2) {
            const uintptr_t rand_lowbits = sfmt_genrand_uint64(sfmtp) & UINT32_MAX;
            uintptr_t cur_geno;
            if (rand_lowbits >= geno2_thresh) {
              cur_geno = 2 + (rand_lowbits >= geno3_thresh);
            } else {
              cur_geno = (rand_lowbits >= geno1_thresh);
            }
            genovec_word = genovec_word << 2;
            genovec_word |= cur_geno;
          }
        } else {
          genovec_word = prev_genovec[widx];
          if (ld_invert) {
            genovec_word = InvertGenoWordUnsafe(genovec_word);
            genovec_word = bzhi_max(genovec_word, loop_len * 2);
          }
          while (1) {
            if (!rand16_left) {
              u64rand = sfmt_genrand_uint64(sfmtp);
              rand16_left = 4;
            }
            const uint32_t sample_idx_lowbits = ((u64rand & 65535) * ld_denom) >> 16;
            u64rand >>= 16;
            --rand16_left;
            if (sample_idx_lowbits >= loop_len) {
              break;
            }
            if (rand16_left < 2) {
              u64rand = sfmt_genrand_uint64(sfmtp);
              rand16_left = 4;
            }
            const uintptr_t rand_lowbits = u64rand & UINT32_MAX;
            u64rand >>= 32;
            rand16_left -= 2;
            uintptr_t cur_geno;
            if (rand_lowbits >= geno2_thresh) {
              cur_geno = 2 + (rand_lowbits >= geno3_thresh);
            } else {
              cur_geno = (rand_lowbits >= geno1_thresh);
            }
            const uint32_t bit_shift_ct = sample_idx_lowbits * 2;
            genovec_word &= ~((3 * k1LU) << bit_shift_ct);
            genovec_word |= cur_geno << bit_shift_ct;
          }
        }
        // set phasepresent bits may correspond to missing or homozygous calls;
        // that's resolved later.
        uint32_t phasepresent_possible_hw = 0;
        uint32_t phaseinfo_hw = 0;
        if (phase_exists) {
          if (phase_is_variable) {
            uint32_t sample_idx_lowbits = 0;
            while (1) {
              sample_idx_lowbits += UpperBoundNonemptyU64(&phase_geomdist[0], kBitsPerWordD2, sfmt_genrand_uint64(sfmtp));
              if (sample_idx_lowbits >= loop_len) {
                break;
              }
              phasepresent_possible_hw |= 1U << sample_idx_lowbits;
              ++sample_idx_lowbits;
            }
          }
          if (phase_invert) {
            phasepresent_possible_hw = S_CAST(Halfword, ~phasepresent_possible_hw);
          }
          const uint32_t kInt16PerHalfword = (kBitsPerWord / 2) / 16;
          if (rand16_left < kInt16PerHalfword) {
            u64rand = sfmt_genrand_uint64(sfmtp);
            rand16_left = 4;
          }
          phaseinfo_hw = S_CAST(Halfword, u64rand);
          u64rand >>= kBitsPerWord / 2;
          rand16_left -= kInt16PerHalfword;
        }
        uint32_t dosage_present_hw = 0;
        uint32_t dphase_present_hw = 0;
        if (dosage_exists) {
          // deliberate overflow
          uint32_t sample_idx_lowbits = UINT32_MAX;
          while (1) {
            ++sample_idx_lowbits;
            if (dosage_geomdist_max) {
              sample_idx_lowbits += UpperBoundNonemptyU64(&dosage_geomdist[0], dosage_geomdist_max, sfmt_genrand_uint64(sfmtp));
            }
            if (sample_idx_lowbits >= loop_len) {
              break;
            }
            if (((genovec_word >> (2 * sample_idx_lowbits)) & 3) == 3) {
              continue;
            }
            if (!rand16_left) {
              u64rand = sfmt_genrand_uint64(sfmtp);
              rand16_left = 4;
            }
            // possible todo: make this respect simulated allele frequency
            //
            // repeated squaring provides a quick and dirty way to do this,
            // e.g. square a uniform 0..1-scaled value 4 times to get an
            // expectation of 1/17.
            const uint32_t dosage_int = ((u64rand & 65535) + 1) / 2;
            u64rand >>= 16;
            --rand16_left;
            const uint32_t halfdist = BiallelicDosageHalfdist(dosage_int);
            if (halfdist < dosage_erase_halfdist) {
              *cur_dosage_main_iter++ = dosage_int;
              const uint32_t shifted_bit = 1U << sample_idx_lowbits;
              dosage_present_hw |= shifted_bit;
              if (phasepresent_possible_hw & shifted_bit) {
                dphase_present_hw |= shifted_bit;
                int32_t cur_dphase_delta = DosageHomdist(dosage_int);
                // To test both the maximum-difference and the
                // nonmaximal-difference cases, we subtract 1 from this number
                // if it's even.
                //
                // Note that, as of this writing, pgenlib_write
                // AppendDphase16() does *not* omit a dphase_delta entry when
                // the implicit value inferred from hardcall-phase and dosage
                // would be correct.  It is the caller's responsibility to
                // exploit this case.  (MakePgenThread doesn't do this, maybe
                // it should.  Exporting to VCF/BCF and re-importing does do
                // the trick.)
                cur_dphase_delta -= 1 - (cur_dphase_delta & 1);
                if (!(phaseinfo_hw & shifted_bit)) {
                  cur_dphase_delta = -cur_dphase_delta;
                }
                *cur_dphase_delta_iter++ = cur_dphase_delta;
              }
              if (halfdist < hard_call_halfdist) {
                genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
                continue;
              }
            }
            genovec_word &= ~((3 * k1LU) << (2 * sample_idx_lowbits));
            genovec_word |= ((dosage_int + (kDosage4th * k1LU)) / kDosageMid) << (2 * sample_idx_lowbits);
          }
        }
        write_genovec_iter[widx] = genovec_word;
        const uint32_t cur_hets = Pack01ToHalfword(genovec_word);
        const uint32_t phasepresent_hw = phasepresent_possible_hw & cur_hets;
        write_phasepresent_iter[widx] = phasepresent_hw;
        write_phaseinfo_iter[widx] = phaseinfo_hw & phasepresent_hw;
        write_dosage_present_iter[widx] = dosage_present_hw;
        write_dphase_present_iter[widx] = dphase_present_hw;
      }
      ZeroTrailingNyps(sample_ct, write_genovec_iter);
      uintptr_t* __attribute__((may_alias)) write_phasepresent_warr = R_CAST(uintptr_t*, write_phasepresent_iter);
      uintptr_t* __attribute__((may_alias)) write_phaseinfo_warr = R_CAST(uintptr_t*, write_phaseinfo_iter);
      uintptr_t* __attribute__((may_alias)) write_dosage_present_warr = R_CAST(uintptr_t*, write_dosage_present_iter);
      ZeroTrailingBits(sample_ct, write_phasepresent_warr);
      ZeroTrailingBits(sample_ct, write_phaseinfo_warr);
      ZeroTrailingBits(sample_ct, write_dosage_present_warr);
      const uint32_t dosage_ct = cur_dosage_main_iter - write_dosage_main_iter;
      *write_dosage_ct_iter++ = dosage_ct;
      const uint32_t dphase_ct = cur_dphase_delta_iter - write_dphase_delta_iter;
      *write_dphase_ct_iter++ = dphase_ct;
      prev_genovec = write_genovec_iter;
      write_genovec_iter = &(write_genovec_iter[sample_ctaw2]);
      write_phasepresent_iter = &(write_phasepresent_iter[2 * sample_ctaw]);
      write_phaseinfo_iter = &(write_phaseinfo_iter[2 * sample_ctaw]);
      write_dosage_present_iter = &(write_dosage_present_iter[2 * sample_ctaw]);
      write_dosage_main_iter = &(write_dosage_main_iter[sample_ct]);
      write_dphase_present_iter = &(write_dphase_present_iter[2 * sample_ctaw]);
      write_dphase_delta_iter = &(write_dphase_delta_iter[sample_ct]);
    }
    parity = 1 - parity;
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

PglErr EigIndToPsam(const char* indname, const char* const_fid, const char* missing_catname, MiscFlags misc_flags, ImportFlags import_flags, uint32_t psam_01, char id_delim, char* outname, char* outname_end, uint32_t* sample_ct_ptr, uint32_t* hash_ptr) {
  // Some overlap with OxSampleToPsam().
  unsigned char* bigstack_mark = g_bigstack_base;
  TextStream ind_txs;
  PreinitTextStream(&ind_txs);
  uintptr_t line_idx = 0;
  FILE* psamfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  {
    // Scanning pass: determine write_fid and write_sid, distinguish phenotype
    // cases.
    // Phenotype cases:
    // 0. all Ignore: no output phenotype column
    // 1. all Control/Case/Ignore: Control -> '1' - psam_01,
    //    Case -> '2' - psam_01, Ignore -> 'NA'
    // 2. otherwise, if all numeric/Ignore, Ignore -> 'NA'
    // 3. otherwise, if all catname/Ignore, Ignore -> missing category name
    // otherwise, discard phenotype
    // (yes, this doesn't correspond to how CONVERTF treats Ignore)
    reterr = InitTextStreamEx(indname, 0, kTextStreamBlenFast, kTextStreamBlenFast, 1, &ind_txs);
    if (unlikely(reterr)) {
      goto EigIndToPsam_ret_TSTREAM_FAIL;
    }
    ImportSampleIdContext isic;
    InitImportSampleIdContext(const_fid, import_flags, id_delim, &isic);
    const uint32_t iid_sid = (misc_flags / kfMiscIidSid) & 1;
    uint32_t write_fid = 0;
    uint32_t write_sid = 0;
    uint32_t pheno_casectrl_seen = 0;
    uint32_t pheno_numeric_seen = 0;
    uint32_t pheno_catname_seen = 0;
    while (1) {
      const char* sample_id_start = TextGet(&ind_txs);
      if (!sample_id_start) {
        break;
      }
      ++line_idx;
      // 1. sample ID
      if (unlikely(IsEolnKns(*sample_id_start))) {
        goto EigIndToPsam_ret_MISSING_TOKENS;
      }
      const char* sample_id_end = CurTokenEnd(sample_id_start);
      if (id_delim) {
        const char* first_delim = S_CAST(const char*, memchr(sample_id_start, ctou32(id_delim), sample_id_end - sample_id_start));
        if (unlikely(!first_delim)) {
          snprintf(g_logbuf, kLogbufSize, "Error: No '%c' in sample ID.\n", id_delim);
          goto EigIndToPsam_ret_INCONSISTENT_INPUT_2;
        }
        const uint32_t first_slen = first_delim - sample_id_start;
        const char* second_part_start = &(first_delim[1]);
        const char* maybe_second_part_end = S_CAST(const char*, memchr(second_part_start, ctou32(id_delim), sample_id_end - second_part_start));
        if (maybe_second_part_end == nullptr) {
          if (!iid_sid) {
            write_fid |= (first_slen != 1) || (sample_id_start[0] != '0');
          } else {
            const uint32_t sid_slen = sample_id_end - second_part_start;
            write_sid |= (sid_slen != 1) || (second_part_start[0] != '0');
          }
        } else {
          write_fid |= (first_slen != 1) || (sample_id_start[0] != '0');
          const char* sid_start = &(maybe_second_part_end[1]);
          const uint32_t sid_slen = sample_id_end - sid_start;
          if (unlikely(memchr(sid_start, ctou32(id_delim), sid_slen))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Too many instances of --id-delim argument '%c' in sample ID.\n", id_delim);
            goto EigIndToPsam_ret_INCONSISTENT_INPUT_2;
          }
          write_sid |= (sid_slen != 1) || (sid_start[0] != '0');
        }
      }
      // 2. sex (ignore for now)
      // 3. phenotype
      const char* pheno_start = NextTokenMult(sample_id_end, 2);
      if (unlikely(!pheno_start)) {
        goto EigIndToPsam_ret_MISSING_TOKENS;
      }
      const uint32_t pheno_slen = strlen_se(pheno_start);
      if (strequal_k(pheno_start, "Ignore", pheno_slen)) {
        continue;
      }
      if (strequal_k(pheno_start, "Case", pheno_slen) ||
          strequal_k(pheno_start, "Control", pheno_slen)) {
        pheno_casectrl_seen = 1;
        continue;
      }
      if (IsCategoricalPhenostrNocsv(pheno_start)) {
        if (unlikely(pheno_slen > kMaxIdSlen)) {
          logerrputs("Error: Categorical phenotypes are limited to " MAX_ID_SLEN_STR " characters.\n");
          goto EigIndToPsam_ret_MALFORMED_INPUT;
        }
        pheno_catname_seen = 1;
      } else {
        pheno_numeric_seen = 1;
      }
    }
    if (unlikely(TextStreamErrcode2(&ind_txs, &reterr))) {
      goto EigIndToPsam_ret_TSTREAM_FAIL;
    }
    if (unlikely(line_idx > 0x7ffffffe)) {
      logerrputs("Error: " PROG_NAME_STR " does not support more than 2^31 - 2 samples.\n");
      goto EigIndToPsam_ret_MALFORMED_INPUT;
    }
    const uint32_t sample_ct = line_idx;
    if (unlikely(!sample_ct)) {
      logerrputs("Error: No samples in .ind file.\n");
      goto EigIndToPsam_ret_DEGENERATE_DATA;
    }
    if (id_delim) {
      if (!iid_sid) {
        isic.fid_delim_mode = write_fid? kImportFidDelimModeAlwaysCopy : kImportFidDelimModeAlwaysOmit;
        isic.sid_delim_mode = write_sid? kImportSidDelimModeCopyOr0 : kImportSidDelimModeNonexistOrOmit;
      } else {
        isic.fid_delim_mode = write_fid? kImportFidDelimModeCopyOr0 : kImportFidDelimModeOmitWhenPresent;
        isic.sid_delim_mode = write_sid? kImportSidDelimModeAlwaysCopy : kImportSidDelimModeNonexistOrOmit;
      }
    } else if (isic.const_fid || isic.double_id) {
      write_fid = 1;
    }
    PhenoDtype pheno_dtype = kPhenoDtypeOther;
    {
      const uint32_t seen_sum = pheno_casectrl_seen + pheno_numeric_seen + pheno_catname_seen;
      if (seen_sum > 1) {
        logerrputs("Warning: Mixed label types in .ind file; no phenotype column will be generated.\n");
      } else if (pheno_casectrl_seen) {
        pheno_dtype = kPhenoDtypeCc;
      } else if (pheno_numeric_seen) {
        pheno_dtype = kPhenoDtypeQt;
      } else if (pheno_catname_seen) {
        pheno_dtype = kPhenoDtypeCat;
      }
    }

    snprintf(outname_end, kMaxOutfnameExtBlen, ".psam");
    if (unlikely(fopen_checked(outname, FOPEN_WB, &psamfile))) {
      goto EigIndToPsam_ret_OPEN_FAIL;
    }
    char* write_iter = g_textbuf;
    char* textbuf_flush = &(write_iter[kMaxMediumLine]);
    *write_iter++ = '#';
    if (write_fid) {
      write_iter = strcpya_k(write_iter, "FID\t");
    }
    write_iter = strcpya_k(write_iter, "IID");
    if (write_sid) {
      write_iter = strcpya_k(write_iter, "\tSID");
    }
    write_iter = strcpya_k(write_iter, "\tSEX");
    if (pheno_dtype != kPhenoDtypeOther) {
      write_iter = strcpya_k(write_iter, "\tPHENO1");
    }
    AppendBinaryEoln(&write_iter);

    reterr = TextRewind(&ind_txs);
    if (unlikely(reterr)) {
      goto EigIndToPsam_ret_TSTREAM_FAIL;
    }
    const char ctrl_char = '1' - psam_01;
    const uint32_t missing_catname_slen = strlen(missing_catname);
    uint32_t eighash = 0;
    uint32_t max_sample_id_slen = 0;
    for (line_idx = 1; line_idx <= sample_ct; ++line_idx) {
      const char* sample_id_start = TextGet(&ind_txs);
      if (unlikely(!sample_id_start)) {
        goto EigIndToPsam_ret_TSTREAM_REWIND_FAIL;
      }
      const char* sample_id_end = FirstSpaceOrEoln(sample_id_start);
      const uint32_t sample_id_slen = sample_id_end - sample_id_start;
      if (sample_id_slen > max_sample_id_slen) {
        max_sample_id_slen = sample_id_slen;
      }
      UpdateEighash(sample_id_start, sample_id_slen, &eighash);
      reterr = ImportSampleId(sample_id_start, sample_id_end, &isic, &write_iter);
      if (unlikely(reterr)) {
        goto EigIndToPsam_ret_1;
      }
      const char* sex_start = FirstNonTspace(sample_id_end);
      const char* sex_end = FirstSpaceOrEoln(sex_start);
      const char* pheno_start = FirstNonTspace(sex_end);
      if (unlikely(IsEolnKns(*pheno_start))) {
        goto EigIndToPsam_ret_REWIND_FAIL;
      }
      const uint32_t sex_slen = sex_end - sex_start;
      if (unlikely(sex_slen != 1)) {
        goto EigIndToPsam_ret_INVALID_SEX;
      }
      *write_iter++ = '\t';
      const char sexchar = *sex_start;
      if (sexchar == 'M') {
        *write_iter++ = '1';
      } else if (sexchar == 'F') {
        *write_iter++ = '2';
      } else if (likely(sexchar == 'U')) {
        write_iter = strcpya_k(write_iter, "NA");
      } else {
        goto EigIndToPsam_ret_INVALID_SEX;
      }
      if (pheno_dtype != kPhenoDtypeOther) {
        *write_iter++ = '\t';
        const uint32_t pheno_slen = CurTokenEnd(pheno_start) - pheno_start;
        if (pheno_dtype == kPhenoDtypeCc) {
          // has to be Control / Case / Ignore
          if (pheno_slen == 7) {
            *write_iter++ = ctrl_char;
          } else if (pheno_slen == 4) {
            *write_iter++ = ctrl_char + 1;
          } else {
            write_iter = strcpya_k(write_iter, "NA");
          }
        } else {
          if (strequal_k(pheno_start, "Ignore", pheno_slen)) {
            if (pheno_dtype == kPhenoDtypeQt) {
              write_iter = strcpya_k(write_iter, "NA");
            } else {
              write_iter = memcpya(write_iter, missing_catname, missing_catname_slen);
            }
          } else {
            write_iter = memcpya(write_iter, pheno_start, pheno_slen);
          }
        }
      }
      AppendBinaryEoln(&write_iter);
      if (unlikely(fwrite_ck(textbuf_flush, psamfile, &write_iter))) {
        goto EigIndToPsam_ret_WRITE_FAIL;
      }
    }
    // We should be at end of file.
    if (unlikely(TextGet(&ind_txs))) {
      goto EigIndToPsam_ret_TSTREAM_REWIND_FAIL;
    }
    if (unlikely(fclose_flush_null(textbuf_flush, write_iter, &psamfile))) {
      goto EigIndToPsam_ret_WRITE_FAIL;
    }
    logprintfww("--eigind: %u sample%s imported to %s .\n", sample_ct, (sample_ct == 1)? "" : "s", outname);
    if (max_sample_id_slen > 39) {
      logerrprintfww("Warning: Longest sample ID had length %u; this violates the 39 character limit in the specification.\n", max_sample_id_slen);
    }
    *sample_ct_ptr = sample_ct;
    *hash_ptr = eighash;
  }
  while (0) {
  EigIndToPsam_ret_TSTREAM_FAIL:
    TextStreamErrPrint(".ind file", &ind_txs);
    break;
  EigIndToPsam_ret_TSTREAM_REWIND_FAIL:
    TextStreamErrPrintRewind(".ind file", &ind_txs, &reterr);
    break;
  EigIndToPsam_ret_REWIND_FAIL:
    logerrprintfww(kErrprintfRewind, ".ind file");
    reterr = kPglRetRewindFail;
    break;
  EigIndToPsam_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  EigIndToPsam_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  EigIndToPsam_ret_MISSING_TOKENS:
    logerrprintfww("Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, indname);
  EigIndToPsam_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  EigIndToPsam_ret_INVALID_SEX:
    logerrprintfww("Error: Line %" PRIuPTR " of %s has an invalid sex entry ('M', 'F', or 'U' expected).\n", line_idx, indname);
    reterr = kPglRetMalformedInput;
    break;
  EigIndToPsam_ret_INCONSISTENT_INPUT_2:
    logerrputsb();
    reterr = kPglRetInconsistentInput;
    break;
  EigIndToPsam_ret_DEGENERATE_DATA:
    reterr = kPglRetDegenerateData;
    break;
  }
 EigIndToPsam_ret_1:
  CleanupTextStream2(".ind file", &ind_txs, &reterr);
  fclose_cond(psamfile);
  BigstackReset(bigstack_mark);
  return reterr;
}

static inline uint32_t GetEigChrCode(const ChrInfo* cip, const char* chr_code_start) {
  const uint32_t autosome_ct = cip->autosome_ct;
  const uint32_t chr_code_raw = GetChrCodeRaw(chr_code_start);
  if (chr_code_raw <= autosome_ct) {
    return chr_code_raw;
  }
  if (chr_code_raw <= autosome_ct + 2) {
    return cip->xymt_codes[chr_code_raw - autosome_ct - 1];
  }
  if ((chr_code_raw == 90) || (chr_code_raw == 91)) {
    // 90 -> MT, 91 -> XY
    return cip->xymt_codes[kChrOffsetMT + 90 - chr_code_raw];
  }
  if (likely((chr_code_raw >= kMaxContigs) && (chr_code_raw < UINT32_MAXM1))) {
    return cip->xymt_codes[chr_code_raw - kMaxContigs];
  }
  // EIGENSOFT also maps 2A and 2B to 2, but I don't see any need to support
  // that.
  return UINT32_MAX;
}

PglErr EigSnpToPvar(const char* snpname, const ChrInfo* cip, ImportFlags import_flags, uint32_t max_thread_ct, char* outname, char* outname_end, uintptr_t** variant_include_ptr, uint32_t* raw_variant_ct_ptr, uint32_t* hash_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  TextStream snp_txs;
  PreinitTextStream(&snp_txs);
  uintptr_t line_idx = 0;

  char* cswritep = nullptr;
  CompressStreamState pvar_css;
  PreinitCstream(&pvar_css);
  PglErr reterr = kPglRetSuccess;
  {
    if (unlikely(cip->autosome_ct > 87)) {
      logerrputs("Error: EIGENSOFT .snp format does not support autosome_ct > 87.\n");
      goto EigSnpToPvar_ret_INCONSISTENT_INPUT;
    }
    // ok to temporarily spend 1/2 of memory on variant_include, don't have to
    // worry about long allele codes
    const uint32_t raw_variant_ct_limit = MINV(0x7ffffffd, bigstack_left() * 4);
    if (unlikely(bigstack_alloc_w(BitCtToWordCt(raw_variant_ct_limit), variant_include_ptr))) {
      goto EigSnpToPvar_ret_NOMEM;
    }
    uintptr_t* variant_include = *variant_include_ptr;
    snprintf(outname_end, kMaxOutfnameExtBlen, ".pvar");
    const uint32_t output_zst = (import_flags / kfImportKeepAutoconvVzs) & 1;
    if (output_zst) {
      snprintf(&(outname_end[5]), kMaxOutfnameExtBlen - 5, ".zst");
    }
    const uintptr_t overflow_buf_size = kCompressStreamBlock + kMaxIdSlen + 64;
    const uint32_t io_thread_ct = MAXV(1, max_thread_ct - 1);
    // io_thread_ct ignored by InitCstreamAlloc() when output_zst false
    reterr = InitCstreamAlloc(outname, 0, output_zst, io_thread_ct, overflow_buf_size, &pvar_css, &cswritep);
    if (unlikely(reterr)) {
      goto EigSnpToPvar_ret_1;
    }

    // First pass: just check for nonzero centimorgan value (after applying
    // chromosome filter).
    reterr = SizeAndInitTextStream(snpname, bigstack_left(), output_zst? 1 : io_thread_ct, &snp_txs);
    if (unlikely(reterr)) {
      goto EigSnpToPvar_ret_TSTREAM_FAIL;
    }
    char* line_iter = TextLineEnd(&snp_txs);
    uint32_t at_least_one_nzero_cm = 0;
    for (; ; line_iter = AdvPastDelim(line_iter, '\n')) {
      reterr = TextGetUnsafe(&snp_txs, &line_iter);
      if (reterr) {
        if (likely(reterr == kPglRetEof)) {
          reterr = kPglRetSuccess;
          break;
        }
        goto EigSnpToPvar_ret_TSTREAM_FAIL;
      }
      ++line_idx;
      char* chr_code_start = FirstNonTspace(FirstSpaceOrEoln(line_iter));
      char* chr_code_end = FirstSpaceOrEoln(chr_code_start);
      char* cm_start = FirstNonTspace(chr_code_end);
      if (unlikely(IsEolnKns(*cm_start))) {
        goto EigSnpToPvar_ret_MISSING_TOKENS;
      }
      char* cm_end = FirstSpaceOrEoln(cm_start);
      line_iter = cm_end;
      // Can't use GetChrCode[Counted]() since MT/XY order is reversed.
      const uint32_t chr_code = GetEigChrCode(cip, chr_code_start);
      if (unlikely(IsI32Neg(chr_code))) {
        *chr_code_end = '\0';
        ChrError(chr_code_start, ".snp file", cip, line_idx, UINT32_MAXM1);
        goto EigSnpToPvar_ret_MALFORMED_INPUT;
      }
      if (!IsSet(cip->chr_mask, chr_code)) {
        continue;
      }
      double cur_cm;
      if (unlikely(!ScantokDouble(cm_start, &cur_cm))) {
        snprintf(g_logbuf, kLogbufSize, "Error: Invalid centimorgan position on line %" PRIuPTR " of .snp file.\n", line_idx);
        goto EigSnpToPvar_ret_MALFORMED_INPUT;
      }
      if (cur_cm != 0.0) {
        at_least_one_nzero_cm = 1;
        break;
      }
    }
    if (unlikely(line_idx == 0)) {
      logerrputs("Error: No variants in .snp file.\n");
      goto EigSnpToPvar_ret_DEGENERATE_DATA;
    }

    cswritep = strcpya_k(cswritep, "#CHROM\tPOS\tID\tREF\tALT");
    if (at_least_one_nzero_cm) {
      cswritep = strcpya_k(cswritep, "\tCM");
    }
    AppendBinaryEoln(&cswritep);

    reterr = TextRewind(&snp_txs);
    if (unlikely(reterr)) {
      goto EigSnpToPvar_ret_TSTREAM_FAIL;
    }
    line_idx = 0;
    line_iter = TextLineEnd(&snp_txs);
    uintptr_t variant_include_word = 0;
    uint32_t eighash = 0;
    for (; TextGetUnsafe2(&snp_txs, &line_iter); line_iter = AdvPastDelim(line_iter, '\n')) {
      const uint32_t variant_uidx = line_idx;
      if (unlikely(variant_uidx == raw_variant_ct_limit)) {
        if (raw_variant_ct_limit == 0x7ffffffd) {
          logerrputs("Error: " PROG_NAME_STR " does not support more than 2^31 - 3 variants.  We recommend using\nother software for very deep studies of small numbers of genomes.\n");
          goto EigSnpToPvar_ret_MALFORMED_INPUT;
        }
        goto EigSnpToPvar_ret_NOMEM;
      }
      ++line_idx;
      char* variant_id_start = line_iter;
      char* variant_id_end = FirstSpaceOrEoln(variant_id_start);
      char* chr_code_start = FirstNonTspace(variant_id_end);
      char* chr_code_end = FirstSpaceOrEoln(chr_code_start);
      char* cm_start = FirstNonTspace(chr_code_end);
      char* cm_end = FirstSpaceOrEoln(cm_start);
      char* bp_start = FirstNonTspace(cm_end);
      char* bp_end = FirstSpaceOrEoln(bp_start);
      char* ref_start = FirstNonTspace(bp_end);
      char* ref_end = FirstSpaceOrEoln(ref_start);
      char* alt_start = FirstNonTspace(ref_end);
      if (unlikely(IsEolnKns(*alt_start))) {
        goto EigSnpToPvar_ret_MISSING_TOKENS;
      }
      char* alt_end = FirstSpaceOrEoln(alt_start);
      line_iter = alt_end;

      const uint32_t variant_id_slen = variant_id_end - variant_id_start;
      UpdateEighash(variant_id_start, variant_id_slen, &eighash);

      const uint32_t chr_code = GetEigChrCode(cip, chr_code_start);
      if (unlikely(IsI32Neg(chr_code))) {
        *chr_code_end = '\0';
        ChrError(chr_code_start, ".snp file", cip, line_idx, UINT32_MAXM1);
        goto EigSnpToPvar_ret_MALFORMED_INPUT;
      }
      const uintptr_t keep_variant = IsSet(cip->chr_mask, chr_code);
      const uint32_t variant_uidx_lowbits = variant_uidx % kBitsPerWord;
      variant_include_word |= keep_variant << (variant_uidx % kBitsPerWord);
      if (variant_uidx_lowbits == kBitsPerWord - 1) {
        variant_include[variant_uidx / kBitsPerWord] = variant_include_word;
        variant_include_word = 0;
      }
      if (!keep_variant) {
        continue;
      }
      cswritep = chrtoa(cip, chr_code, cswritep);
      *cswritep++ = '\t';

      uint32_t cur_bp;
      if (unlikely(ScanUintDefcap(bp_start, &cur_bp))) {
        logerrprintf("Error: Invalid bp coordinate on line %" PRIuPTR " of .snp file.\n", line_idx);
        goto EigSnpToPvar_ret_MALFORMED_INPUT;
      }
      cswritep = u32toa_x(cur_bp, '\t', cswritep);

      if (unlikely(variant_id_slen > kMaxIdSlen)) {
        logerrputs("Error: Variant names are limited to " MAX_ID_SLEN_STR " characters.\n");
        goto EigSnpToPvar_ret_MALFORMED_INPUT;
      }
      cswritep = memcpyax(cswritep, variant_id_start, variant_id_slen, '\t');

      const uint32_t ref_slen = ref_end - ref_start;
      const uint32_t alt_slen = alt_end - alt_start;
      if (unlikely((ref_slen != 1) || (alt_slen != 1))) {
        logerrprintf("Error: Multi-character allele code on line %" PRIuPTR " of .snp file.\n", line_idx);
        goto EigSnpToPvar_ret_MALFORMED_INPUT;
      }
      *cswritep++ = *ref_start;
      *cswritep++ = '\t';
      char alt_char = *alt_start;
      if (alt_char == 'X') {
        alt_char = '.';
      }
      *cswritep++ = alt_char;

      if (at_least_one_nzero_cm) {
        double cur_cm;
        char* cm_parse_end = ScanadvDouble(cm_start, &cur_cm);
        if (unlikely(cm_parse_end != cm_end)) {
          logerrprintf("Error: Invalid centimorgan position on line %" PRIuPTR " of .snp file.\n", line_idx);
          goto EigSnpToPvar_ret_MALFORMED_INPUT;
        }
        *cswritep++ = '\t';
        cswritep = dtoa_g_p8(cur_cm, cswritep);
      }

      AppendBinaryEoln(&cswritep);
      if (unlikely(Cswrite(&pvar_css, &cswritep))) {
        goto EigSnpToPvar_ret_WRITE_FAIL;
      }
    }
    if (unlikely(TextStreamErrcode2(&snp_txs, &reterr))) {
      goto EigSnpToPvar_ret_TSTREAM_REWIND_FAIL;
    }
    const uint32_t raw_variant_ct = line_idx;
    if (unlikely(raw_variant_ct == 0)) {
      // this isn't what we saw the first time around
      goto EigSnpToPvar_ret_TSTREAM_REWIND_FAIL;
    }
    uintptr_t* last_variant_include_word_ptr = &(variant_include[(raw_variant_ct - 1) / kBitsPerWord]);
    if (raw_variant_ct % kBitsPerWord) {
      *last_variant_include_word_ptr = variant_include_word;
    }
    if (unlikely(CswriteCloseNull(&pvar_css, cswritep))) {
      goto EigSnpToPvar_ret_WRITE_FAIL;
    }
    bigstack_mark = R_CAST(unsigned char*, RoundUpPow2(R_CAST(uintptr_t, &(last_variant_include_word_ptr[1])), kCacheline));
    *raw_variant_ct_ptr = raw_variant_ct;
    *hash_ptr = eighash;
  }
  while (0) {
  EigSnpToPvar_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  EigSnpToPvar_ret_TSTREAM_FAIL:
    TextStreamErrPrint(".snp file", &snp_txs);
    break;
  EigSnpToPvar_ret_TSTREAM_REWIND_FAIL:
    TextStreamErrPrintRewind(".snp file", &snp_txs, &reterr);
    break;
  EigSnpToPvar_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  EigSnpToPvar_ret_MISSING_TOKENS:
    logerrprintfww("Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, snpname);
  EigSnpToPvar_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  EigSnpToPvar_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  EigSnpToPvar_ret_DEGENERATE_DATA:
    reterr = kPglRetDegenerateData;
    break;
  }
 EigSnpToPvar_ret_1:
  CswriteCloseCond(&pvar_css, cswritep);
  CleanupTextStream2(snpname, &snp_txs, &reterr);
  BigstackReset(bigstack_mark);
  return reterr;
}

typedef struct EigGenoToPgenCtxStruct {
  uint32_t sample_ct;

  // loadbufs[i] transformed to writebufs[i] by thread 0, then
  // compressed/written by thread 1
  // Must be safe to overread loadbufs by a vector-length.
  unsigned char* loadbufs[2];
  uint32_t block_transform_cts[2];
  uintptr_t* writebufs[2];
  uint32_t block_write_cts[2];

  STPgenWriter* spgwp;

  PglErr write_reterr;
  int32_t write_errno;
} EigGenoToPgenCtx;

THREAD_FUNC_DECL EigGenoToPgenThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uint32_t tidx = arg->tidx;
  EigGenoToPgenCtx* ctx = S_CAST(EigGenoToPgenCtx*, arg->sharedp->context);

  const uint32_t sample_ct = ctx->sample_ct;
  const uint32_t sample_ctaw2 = NypCtToAlignedWordCt(sample_ct);
  uint32_t parity = 0;
  if (tidx == 0) {
    // tidx 0 = transform thread
    const uint32_t rec_copy_blen = DivUp(sample_ct, 4);
    const uint32_t rec_blen = MAXV(48, rec_copy_blen);
    // 1. Reverse nyp order.
    // 2. Perform the following nyp mapping:
    //      0 -> 2
    //      1 -> 1
    //      2 -> 0
    //      3 -> 3
#ifdef USE_SSE2
    const uint32_t vec_ct = sample_ctaw2 / kWordsPerVec;
    const VecW mask_0f0f = vecw_set1(kMask0F0F);
#  ifdef USE_SHUFFLE8
    // Use lookup table to finish reversing nyp order, perform the 0 <-> 2
    // nyp mapping, and (for the originally-low nybble) perform the final
    // left-shift-4 at the same time.
    //    0 ->  0 -> 10
    //    1 ->  4 ->  6
    //    2 ->  8 ->  2
    //    3 -> 12 -> 14
    //    4 ->  1 ->  9
    //    5 ->  5 ->  5
    //    6 ->  9 ->  1
    //    7 -> 13 -> 13
    //    8 ->  2 ->  8
    //    9 ->  6 ->  4
    //   10 -> 10 ->  0
    //   11 -> 14 -> 12
    //   12 ->  3 -> 11
    //   13 ->  7 ->  7
    //   14 -> 11 ->  3
    //   15 -> 15 -> 15
    const VecW eig_lookup_high = vecw_setr8(10, 6, 2, 14, 9, 5, 1, 13,
                                            8, 4, 0, 12, 11, 7, 3, 15);
    const VecW eig_lookup_low = vecw_slli(eig_lookup_high, 4);
#  else
    const VecW mask_3333 = vecw_set1(kMask3333);
    const VecW mask_5555 = vecw_set1(kMask5555);
#  endif
#else  // !USE_SSE2
    const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
#endif
    do {
      const uint32_t cur_transform_ct = ctx->block_transform_cts[parity];
      unsigned char* loadbuf_iter = ctx->loadbufs[parity];
      uintptr_t* writebuf_iter = ctx->writebufs[parity];
      for (uint32_t bidx = 0; bidx != cur_transform_ct; ++bidx) {
#ifdef USE_SSE2
        unsigned char* loadbuf_viter = loadbuf_iter;
        VecW* write_viter = R_CAST(VecW*, writebuf_iter);
        for (uint32_t vidx = 0; vidx != vec_ct; ++vidx) {
          const VecW input_vec = vecw_loadu(loadbuf_viter);
          loadbuf_viter = &(loadbuf_viter[kBytesPerVec]);
          const VecW high_nybbles_shifted = vecw_srli(input_vec, 4) & mask_0f0f;
          const VecW low_nybbles = input_vec & mask_0f0f;
#  ifdef USE_SHUFFLE8
          const VecW low_transformed = vecw_shuffle8(eig_lookup_low, low_nybbles);
          const VecW high_transformed = vecw_shuffle8(eig_lookup_high, high_nybbles_shifted);
          const VecW result = low_transformed | high_transformed;
#  else
          const VecW nybble_swapped = high_nybbles_shifted | vecw_slli(low_nybbles, 4);
          const VecW high_nyps = vecw_and_notfirst(mask_3333, nybble_swapped);
          const VecW low_nyps = mask_3333 & nybble_swapped;
          const VecW nyp_reversed = vecw_srli(high_nyps, 2) | vecw_slli(low_nyps, 2);
          const VecW is_even = vecw_and_notfirst(nyp_reversed, mask_5555);
          const VecW result = nyp_reversed ^ vecw_slli(is_even, 1);
#  endif
          *write_viter++ = result;
        }
#else
        for (uint32_t widx = 0; widx != sample_ctl2; ++widx) {
          uintptr_t input_word;
          CopyFromUnalignedOffsetW(&input_word, loadbuf_iter, widx);
          const uintptr_t high_nybbles = input_word & (~kMask0F0F);
          const uintptr_t low_nybbles = input_word & kMask0F0F;
          const uintptr_t nybble_swapped = (high_nybbles >> 4) | (low_nybbles << 4);
          const uintptr_t high_nyps = nybble_swapped & (~kMask3333);
          const uintptr_t low_nyps = nybble_swapped & kMask3333;
          const uintptr_t nyp_reversed = (high_nyps >> 2) | (low_nyps << 2);
          const uintptr_t is_even = (~nyp_reversed) & kMask5555;
          writebuf_iter[widx] = nyp_reversed ^ (is_even << 1);
        }
#endif
        ZeroTrailingNyps(sample_ct, writebuf_iter);
        loadbuf_iter = &(loadbuf_iter[rec_blen]);
        writebuf_iter = &(writebuf_iter[sample_ctaw2]);
      }
      ctx->block_write_cts[parity] = cur_transform_ct;
      parity = 1 - parity;
    } while (!THREAD_BLOCK_FINISH(arg));
    THREAD_RETURN;
  }
  // tidx 1 = compress/write thread (expected to be the bottleneck)
  // note that this does nothing during the first timestep
  STPgenWriter* spgwp = ctx->spgwp;
  while (!THREAD_BLOCK_FINISH(arg)) {
    const uint32_t cur_write_ct = ctx->block_write_cts[parity];
    uintptr_t* writebuf_iter = ctx->writebufs[parity];
    for (uint32_t bidx = 0; bidx != cur_write_ct; ++bidx) {
      PglErr reterr = SpgwAppendBiallelicGenovec(writebuf_iter, spgwp);
      if (unlikely(reterr)) {
        ctx->write_reterr = reterr;
        ctx->write_errno = errno;
        break;
      }
      writebuf_iter = &(writebuf_iter[sample_ctaw2]);
    }
    parity = 1 - parity;
  }
  THREAD_RETURN;
}

PglErr EigGenoToPgen(const char* outname, const char* genoname, const uintptr_t* variant_include,  uint32_t variant_ct, uint32_t sample_ct, FILE* genofile) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  // This is a decent fit for MTPgenWriter -- .pgen compression represents a
  // relatively high share of the computational workload -- but STPgenWriter
  // should be fast enough, and doesn't have MTPgenWriter's memory usage
  // problem for million+ sample datasets.
  STPgenWriter spgw;
  PreinitSpgw(&spgw);
  ThreadGroup tg;
  PreinitThreads(&tg);
  EigGenoToPgenCtx ctx;
  {
    uintptr_t spgw_alloc_cacheline_ct;
    uint32_t max_vrec_len;
    reterr = SpgwInitPhase1(outname, nullptr, nullptr, variant_ct, sample_ct, 0, kPgenWriteBackwardSeek, kfPgenGlobal0, 1, &spgw, &spgw_alloc_cacheline_ct, &max_vrec_len);
    if (unlikely(reterr)) {
      if (reterr == kPglRetOpenFail) {
        logerrprintfww(kErrprintfFopen, outname, strerror(errno));
      }
      goto EigGenoToPgen_ret_1;
    }
    unsigned char* spgw_alloc;
    if (unlikely(bigstack_alloc_uc(spgw_alloc_cacheline_ct * kCacheline, &spgw_alloc))) {
      goto EigGenoToPgen_ret_NOMEM;
    }
    SpgwInitPhase2(max_vrec_len, &spgw, spgw_alloc);

    const uint32_t sample_ctaw2 = NypCtToAlignedWordCt(sample_ct);
    const uint32_t input_rec_blen = MAXV(48, DivUp(sample_ct, 4));
    uintptr_t cachelines_avail_m4 = bigstack_left() / kCacheline;
    if (unlikely(cachelines_avail_m4 < 4)) {
      goto EigGenoToPgen_ret_NOMEM;
    }
    cachelines_avail_m4 -= 4;
    const uintptr_t bytes_req_per_in_block_variant = 2 * (input_rec_blen + sample_ctaw2 * sizeof(intptr_t));
    uintptr_t main_block_size = (cachelines_avail_m4 * kCacheline) / bytes_req_per_in_block_variant;
    if (main_block_size > 65536) {
      main_block_size = 65536;
    } else if (unlikely(main_block_size < 8)) {
      // this threshold is arbitrary
      goto EigGenoToPgen_ret_NOMEM;
    }

    if (unlikely(SetThreadCt(2, &tg))) {
      goto EigGenoToPgen_ret_NOMEM;
    }
    ctx.sample_ct = sample_ct;
    unsigned char* loadbufs[2];
    if (unlikely(bigstack_alloc_uc(input_rec_blen * main_block_size, &(loadbufs[0])) ||
                 bigstack_alloc_uc(input_rec_blen * main_block_size, &(loadbufs[1])) ||
                 bigstack_alloc_w(sample_ctaw2 * main_block_size, &(ctx.writebufs[0])) ||
                 bigstack_alloc_w(sample_ctaw2 * main_block_size, &(ctx.writebufs[1])))) {
      // this should be impossible
      assert(0);
      goto EigGenoToPgen_ret_NOMEM;
    }
    ctx.loadbufs[0] = loadbufs[0];
    ctx.loadbufs[1] = loadbufs[1];
    ctx.spgwp = &spgw;
    ctx.write_reterr = kPglRetSuccess;
    ctx.write_errno = 0;
    SetThreadFuncAndData(EigGenoToPgenThread, &ctx, &tg);

    fputs("--eiggeno: 0%", stdout);
    uint32_t variant_uidx = 0;
    uint32_t parity = 0;
    uint32_t pct = 0;
    uint32_t next_print_idx = (variant_ct + 99) / 100;
    for (uint32_t variant_idx = 0; ; ) {
      const uint32_t cur_block_size = MINV(main_block_size, variant_ct - variant_idx);
      if (cur_block_size) {
        unsigned char* loadbuf_iter = loadbufs[parity];
        uint32_t bidx = 0;
        do {
          const uint32_t variant_uidx_start = AdvTo1Bit(variant_include, variant_uidx);
          const uint32_t variant_uidx_stop = AdvBoundedTo0Bit(variant_include, variant_uidx_start, variant_uidx_start + cur_block_size - bidx);
          const uintptr_t load_vct = variant_uidx_stop - variant_uidx_start;
          const uintptr_t load_bct = load_vct * input_rec_blen;
          if (unlikely(fseeko(genofile, (variant_uidx_start + 1LLU) * input_rec_blen, SEEK_SET) ||
                       fread_checked(loadbuf_iter, load_bct, genofile))) {
            goto EigGenoToPgen_ret_READ_FAIL;
          }
          variant_uidx = variant_uidx_stop;
          bidx += load_vct;
          loadbuf_iter = &(loadbuf_iter[load_bct]);
        } while (bidx != cur_block_size);
      }
      ctx.block_transform_cts[parity] = cur_block_size;
      if (variant_idx) {
        JoinThreads(&tg);
        reterr = ctx.write_reterr;
        if (unlikely(reterr)) {
          if (reterr == kPglRetWriteFail) {
            errno = ctx.write_errno;
          }
          goto EigGenoToPgen_ret_1;
        }
      }
      if (!cur_block_size) {
        DeclareLastThreadBlock(&tg);
      }
      if (unlikely(SpawnThreads(&tg))) {
        goto EigGenoToPgen_ret_THREAD_CREATE_FAIL;
      }
      if (!cur_block_size) {
        break;
      }
      if (variant_idx >= next_print_idx) {
        if (pct > 10) {
          putc_unlocked('\b', stdout);
        }
        pct = (variant_idx * 100LLU) / variant_ct;
        printf("\b\b%u%%", pct++);
        fflush(stdout);
        next_print_idx = (pct * S_CAST(uint64_t, variant_ct) + 99) / 100;
      }
      variant_idx += cur_block_size;
      parity = 1 - parity;
    }
    JoinThreads(&tg);
    reterr = ctx.write_reterr;
    if (unlikely(reterr)) {
      if (reterr == kPglRetWriteFail) {
        errno = ctx.write_errno;
      }
      goto EigGenoToPgen_ret_1;
    }
    reterr = SpgwFinish(&spgw);
    if (unlikely(reterr)) {
      goto EigGenoToPgen_ret_1;
    }
  }
  while (0) {
  EigGenoToPgen_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  EigGenoToPgen_ret_READ_FAIL:
    if (feof_unlocked(genofile)) {
      errno = 0;
    }
    putc_unlocked('\n', stdout);
    logerrprintfww(kErrprintfFread, genoname, rstrerror(errno));
    reterr = kPglRetReadFail;
    break;
  EigGenoToPgen_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 EigGenoToPgen_ret_1:
  CleanupSpgw(&spgw, &reterr);
  CleanupThreads(&tg);
  BigstackReset(bigstack_mark);
  return reterr;
}

typedef struct EigTgenoRecodeCtxStruct {
  uint32_t sample_ct;
  uint32_t raw_load_batch_size;
  uint32_t raw_load_batch_ct;
  uint32_t cur_vidx_ctaw2;
  uintptr_t* tgeno_loadbuf;
} EigTgenoRecodeCtx;

void EigTgenoRecodeMain(uintptr_t tidx, uintptr_t thread_ct, EigTgenoRecodeCtx* ctx) {
  // division of labor: (2 * max_thread_ct - 1) pieces.  Extra threads launched
  // when all but last piece has been loaded.  Each extra thread recodes two
  // pieces.  Original thread then loads last piece and recodes it.
  const uint32_t sample_ct = ctx->sample_ct;
  const uint32_t raw_load_batch_ct = ctx->raw_load_batch_ct;
  const uint32_t piece_ct = thread_ct * 2 - 1;
  const uint32_t batch_idx_start = (raw_load_batch_ct * S_CAST(uint64_t, tidx) * 2) / piece_ct;
  const uint32_t raw_load_batch_size = ctx->raw_load_batch_size;
  const uint32_t sample_idx_start = batch_idx_start * raw_load_batch_size;
  const uint32_t sample_idx_end = (tidx == thread_ct - 1)? sample_ct : (((raw_load_batch_ct * (S_CAST(uint64_t, tidx) + 1) * 2) / piece_ct) * raw_load_batch_size);
  const uintptr_t cur_vidx_ctaw2 = ctx->cur_vidx_ctaw2;
  uintptr_t* tgeno_loadbuf_iter = &(ctx->tgeno_loadbuf[sample_idx_start * cur_vidx_ctaw2]);
#ifdef USE_SSE2
  const uintptr_t vec_ct = (sample_idx_end - sample_idx_start) * (cur_vidx_ctaw2 / kWordsPerVec);
  const VecW mask_0f0f = vecw_set1(kMask0F0F);
#  ifdef USE_SHUFFLE8
  const VecW eig_lookup_high = vecw_setr8(10, 6, 2, 14, 9, 5, 1, 13,
                                          8, 4, 0, 12, 11, 7, 3, 15);
  const VecW eig_lookup_low = vecw_slli(eig_lookup_high, 4);
#  else
  const VecW mask_3333 = vecw_set1(kMask3333);
  const VecW mask_5555 = vecw_set1(kMask5555);
#  endif
  VecW* tgeno_loadbuf_viter = R_CAST(VecW*, tgeno_loadbuf_iter);
  for (uintptr_t vidx = 0; vidx != vec_ct; ++vidx) {
    const VecW input_vec = *tgeno_loadbuf_viter;
    const VecW high_nybbles_shifted = vecw_srli(input_vec, 4) & mask_0f0f;
    const VecW low_nybbles = input_vec & mask_0f0f;
#  ifdef USE_SHUFFLE8
    const VecW low_transformed = vecw_shuffle8(eig_lookup_low, low_nybbles);
    const VecW high_transformed = vecw_shuffle8(eig_lookup_high, high_nybbles_shifted);
    const VecW result = low_transformed | high_transformed;
#  else
    const VecW nybble_swapped = high_nybbles_shifted | vecw_slli(low_nybbles, 4);
    const VecW high_nyps = vecw_and_notfirst(mask_3333, nybble_swapped);
    const VecW low_nyps = mask_3333 & nybble_swapped;
    const VecW nyp_reversed = vecw_srli(high_nyps, 2) | vecw_slli(low_nyps, 2);
    const VecW is_even = vecw_and_notfirst(nyp_reversed, mask_5555);
    const VecW result = nyp_reversed ^ vecw_slli(is_even, 1);
#  endif
    *tgeno_loadbuf_viter++ = result;
  }
#else
  const uintptr_t word_ct = (sample_idx_end - sample_idx_start) * cur_vidx_ctaw2;
  for (uintptr_t widx = 0; widx != word_ct; ++widx) {
    const uintptr_t input_word = *tgeno_loadbuf_iter;
    const uintptr_t high_nybbles = input_word & (~kMask0F0F);
    const uintptr_t low_nybbles = input_word & kMask0F0F;
    const uintptr_t nybble_swapped = (high_nybbles >> 4) | (low_nybbles << 4);
    const uintptr_t high_nyps = nybble_swapped & (~kMask3333);
    const uintptr_t low_nyps = nybble_swapped & kMask3333;
    const uintptr_t nyp_reversed = (high_nyps >> 2) | (low_nyps << 2);
    const uintptr_t is_even = (~nyp_reversed) & kMask5555;
    *tgeno_loadbuf_iter++ = nyp_reversed ^ (is_even << 1);
  }
#endif
}

THREAD_FUNC_DECL EigTgenoRecodeThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uint32_t tidx = arg->tidx;
  const uint32_t thread_ct = GetThreadCt(arg->sharedp) + 1;
  EigTgenoRecodeCtx* ctx = S_CAST(EigTgenoRecodeCtx*, arg->sharedp->context);
  do {
    EigTgenoRecodeMain(tidx, thread_ct, ctx);
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

typedef struct EigTgenoTransposeCtxStruct {
  uint32_t sample_ct;
  uint32_t loadbuf_ul_stride;

  uintptr_t* tgeno_loadbuf_iter;

  VecW** thread_vecaligned_bufs;
  uintptr_t** thread_write_genovecs;
  PgenWriterCommon** pwcs;

  uint32_t cur_block_write_ct;
} EigTgenoTransposeCtx;

void EigTgenoTransposeMain(uint32_t tidx, EigTgenoTransposeCtx* ctx) {
  // nearly identical to Plink1SmajTransposeMain().
  const uint32_t sample_ct = ctx->sample_ct;
  const uint32_t sample_ctaw2 = NypCtToAlignedWordCt(sample_ct);
  const uintptr_t transpose_block_ct_m1 = (sample_ct - 1) / kPglNypTransposeBatch;
  PgenWriterCommon* pwcp = ctx->pwcs[tidx];
  VecW* vecaligned_buf = ctx->thread_vecaligned_bufs[tidx];
  uintptr_t* write_genovec = ctx->thread_write_genovecs[tidx];
  const uintptr_t cur_block_write_ct = ctx->cur_block_write_ct;
  const uintptr_t loadbuf_ul_stride = ctx->loadbuf_ul_stride;
  uint32_t write_idx = tidx * kPglVblockSize;
  uintptr_t* read_iter = &(ctx->tgeno_loadbuf_iter[write_idx / kBitsPerWordD2]);
  const uint32_t write_idx_end = MINV(write_idx + kPglVblockSize, cur_block_write_ct);
  while (write_idx < write_idx_end) {
    const uintptr_t* read_iter2 = read_iter;
    // uintptr_t* write_iter = write_genovec;
    const uint32_t vblock_size = MINV(kPglNypTransposeBatch, write_idx_end - write_idx);
    uint32_t read_batch_size = kPglNypTransposeBatch;
    for (uintptr_t transpose_block_idx = 0; ; ++transpose_block_idx) {
      if (transpose_block_idx >= transpose_block_ct_m1) {
        if (transpose_block_idx > transpose_block_ct_m1) {
          break;
        }
        read_batch_size = ModNz(sample_ct, kPglNypTransposeBatch);
      }
      TransposeNypblock(read_iter2, loadbuf_ul_stride, sample_ctaw2, read_batch_size, vblock_size, &(write_genovec[transpose_block_idx * kPglNypTransposeWords]), vecaligned_buf);
      read_iter2 = &(read_iter2[kPglNypTransposeBatch * loadbuf_ul_stride]);
    }
    uintptr_t* cur_write_genovec = write_genovec;
    for (uint32_t uii = 0; uii != vblock_size; ++uii, cur_write_genovec = &(cur_write_genovec[sample_ctaw2])) {
      ZeroTrailingNyps(sample_ct, cur_write_genovec);
      PwcAppendBiallelicGenovec(cur_write_genovec, pwcp);
    }
    write_idx += vblock_size;
    read_iter = &(read_iter[kPglNypTransposeWords]);
  }
}

THREAD_FUNC_DECL EigTgenoTransposeThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uint32_t tidx = arg->tidx;
  EigTgenoTransposeCtx* ctx = S_CAST(EigTgenoTransposeCtx*, arg->sharedp->context);
  do {
    EigTgenoTransposeMain(tidx, ctx);
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

static_assert((kPglNypTransposeBatch % kNypsPerVec == 0) && (kPglNypTransposeBatch % kBitsPerWord == 0) && ((kPglNypTransposeBatch & (kPglNypTransposeBatch - 1)) == 0), "EigTgenoToPgen() needs to be updated.");
PglErr EigTgenoToPgen(const char* outname, const char* genoname, const uintptr_t* variant_include, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t sample_ct, uint32_t max_thread_ct, FILE* genofile) {
  // See Plink1SampleMajorToPgen().
  // Wouldn't be ridiculous to merge into a single function.
  unsigned char* bigstack_mark = g_bigstack_base;
  MTPgenWriter* mpgwp = nullptr;
  PglErr reterr = kPglRetSuccess;
  ThreadGroup recode_tg;
  PreinitThreads(&recode_tg);
  ThreadGroup write_tg;
  PreinitThreads(&write_tg);
  EigTgenoRecodeCtx recode_ctx;
  EigTgenoTransposeCtx write_ctx;
  STPgenWriter spgw;
  PreinitSpgw(&spgw);
  {
    const uint32_t input_rec_blen = MAXV(48, NypCtToByteCt(raw_variant_ct));
    unsigned char* raw_loadbuf = nullptr;
    uint32_t raw_load_batch_size = 1;
#ifndef __APPLE__
    if (input_rec_blen < 5120) {
      // assuming 4K block size, fseek won't let us avoid reading many
      // unnecessary disk blocks
      raw_load_batch_size += 131071 / input_rec_blen;
      if (unlikely(bigstack_alloc_uc(raw_load_batch_size * input_rec_blen, &raw_loadbuf))) {
        goto EigTgenoToPgen_ret_NOMEM;
      }
      if (unlikely(fseeko(genofile, 48, SEEK_SET))) {
        goto EigTgenoToPgen_ret_READ_FAIL;
      }
    }
#else
    // macOS seems to suck at seek/read interleaving
    if (input_rec_blen < 1048576) {
      raw_load_batch_size += 2097151 / input_rec_blen;
      if (unlikely(bigstack_alloc_uc(raw_load_batch_size * input_rec_blen, &raw_loadbuf))) {
        goto EigTgenoToPgen_ret_NOMEM;
      }
      if (unlikely(fseeko(genofile, 48, SEEK_SET))) {
        goto EigTgenoToPgen_ret_READ_FAIL;
      }
    }
#endif
    const uint32_t raw_load_batch_ct_m1 = (sample_ct - 1) / raw_load_batch_size;
    if (!raw_load_batch_ct_m1) {
      raw_load_batch_size = sample_ct;
    }

    const uintptr_t sample_ctaw2 = NypCtToAlignedWordCt(sample_ct);
    uint32_t cur_vidx_ct = raw_variant_ct;
    uintptr_t cur_vidx_ctaw2 = NypCtToAlignedWordCt(cur_vidx_ct);
    uintptr_t* tgeno_loadbuf = nullptr;
    uint32_t write_thread_ct = 0;
    unsigned char* bigstack_mark2 = g_bigstack_base;
    while ((max_thread_ct > 1) && (variant_ct == raw_variant_ct) && (variant_ct > kPglVblockSize * 2)) {
      if (bigstack_alloc_w(sample_ct * cur_vidx_ctaw2, &tgeno_loadbuf)) {
        break;
      }
      uintptr_t alloc_base_cacheline_ct;
      uint64_t mpgw_per_thread_cacheline_ct;
      uint32_t vrec_len_byte_ct;
      uint64_t vblock_cacheline_ct;
      MpgwInitPhase1(nullptr, variant_ct, sample_ct, kfPgenGlobal0, &alloc_base_cacheline_ct, &mpgw_per_thread_cacheline_ct, &vrec_len_byte_ct, &vblock_cacheline_ct);
#ifndef __LP64__
      if ((mpgw_per_thread_cacheline_ct > (0x7fffffff / kCacheline)) || (vblock_cacheline_ct > (0x7fffffff / kCacheline))) {
        break;
      }
#endif

      write_thread_ct = DivUp(variant_ct, kPglVblockSize);
      if (write_thread_ct >= max_thread_ct) {
        write_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
      }
      mpgwp = S_CAST(MTPgenWriter*, bigstack_alloc((write_thread_ct + DivUp(sizeof(MTPgenWriter), kBytesPerWord)) * sizeof(intptr_t)));
      if (!mpgwp) {
        break;
      }
      PreinitMpgw(mpgwp);
      if (bigstack_alloc_vp(write_thread_ct, &write_ctx.thread_vecaligned_bufs) ||
          bigstack_alloc_wp(write_thread_ct, &write_ctx.thread_write_genovecs)) {
        mpgwp = nullptr;
        break;
      }
      write_ctx.pwcs = &(mpgwp->pwcs[0]);
      uintptr_t cachelines_avail = bigstack_left() / kCacheline;
      // inner loop transposes kPglNypTransposeBatch variants at a time
      const uintptr_t transpose_thread_cacheline_ct = kPglNypTransposeBufbytes / kCacheline + NypCtToVecCt(sample_ct) * (kPglNypTransposeBatch / kVecsPerCacheline);
      if (cachelines_avail < write_thread_ct * S_CAST(uint64_t, transpose_thread_cacheline_ct)) {
        mpgwp = nullptr;
        break;
      }
      for (uint32_t tidx = 0; tidx != write_thread_ct; ++tidx) {
        write_ctx.thread_vecaligned_bufs[tidx] = S_CAST(VecW*, bigstack_alloc_raw(kPglNypTransposeBufbytes));
        write_ctx.thread_write_genovecs[tidx] = S_CAST(uintptr_t*, bigstack_alloc_raw(NypCtToVecCt(sample_ct) * kBytesPerVec * kPglNypTransposeBatch));
      }
      cachelines_avail = bigstack_left() / kCacheline;
      // We already allocated tgeno_loadbuf, so compression write buffers are
      // free to use all of remaining workspace.
      if (cachelines_avail < alloc_base_cacheline_ct + mpgw_per_thread_cacheline_ct * write_thread_ct) {
        if (cachelines_avail >= alloc_base_cacheline_ct + 2 * mpgw_per_thread_cacheline_ct) {
          write_thread_ct = (cachelines_avail - alloc_base_cacheline_ct) / mpgw_per_thread_cacheline_ct;
        } else {
          mpgwp = nullptr;
          break;
        }
      }
      // Ok, we have enough memory to support at least 2 compressor threads.
      unsigned char* mpgw_alloc = S_CAST(unsigned char*, bigstack_alloc_raw((alloc_base_cacheline_ct + mpgw_per_thread_cacheline_ct * write_thread_ct) * kCacheline));
      reterr = MpgwInitPhase2(outname, nullptr, variant_ct, sample_ct, kPgenWriteBackwardSeek, kfPgenGlobal0, 1, vrec_len_byte_ct, vblock_cacheline_ct, write_thread_ct, mpgw_alloc, mpgwp);
      if (unlikely(reterr)) {
        if (reterr == kPglRetOpenFail) {
          logerrprintfww(kErrprintfFopen, outname, strerror(errno));
        }
        goto EigTgenoToPgen_ret_1;
      }
      write_ctx.sample_ct = sample_ct;
      write_ctx.loadbuf_ul_stride = cur_vidx_ctaw2;
      if (unlikely(SetThreadCt(write_thread_ct - 1, &write_tg))) {
        goto EigTgenoToPgen_ret_NOMEM;
      }
      SetThreadFuncAndData(EigTgenoTransposeThread, &write_ctx, &write_tg);
      break;
    }
    VecW* vecaligned_buf = nullptr;
    uintptr_t* write_genovec = nullptr;
    if (!mpgwp) {
      BigstackReset(bigstack_mark2);
      uintptr_t spgw_alloc_cacheline_ct;
      uint32_t max_vrec_len;
      reterr = SpgwInitPhase1(outname, nullptr, nullptr, variant_ct, sample_ct, 0, kPgenWriteBackwardSeek, kfPgenGlobal0, 1, &spgw, &spgw_alloc_cacheline_ct, &max_vrec_len);
      if (unlikely(reterr)) {
        if (reterr == kPglRetOpenFail) {
          logerrprintfww(kErrprintfFopen, outname, strerror(errno));
        }
        goto EigTgenoToPgen_ret_1;
      }
      unsigned char* spgw_alloc;
      if (unlikely(bigstack_alloc_uc(spgw_alloc_cacheline_ct * kCacheline, &spgw_alloc) ||
                   bigstack_alloc_v(kPglNypTransposeBufbytes / kBytesPerVec, &vecaligned_buf) ||
                   bigstack_alloc_w(sample_ctaw2 * kPglNypTransposeBatch, &write_genovec))) {
        goto EigTgenoToPgen_ret_NOMEM;
      }
      SpgwInitPhase2(max_vrec_len, &spgw, spgw_alloc);

      uintptr_t cachelines_avail = bigstack_left() / kCacheline;
      const uint64_t full_load_vecs_req = sample_ct * S_CAST(uint64_t, NypCtToVecCt(raw_variant_ct));
      if (full_load_vecs_req > cachelines_avail * kVecsPerCacheline) {
        // Load the largest multiple of kPglNypTransposeBatch variants at a
        // time that fits into the remaining workspace.
        const uint64_t min_load_cl_req = DivUpU64(sample_ct * NypCtToVecCt(kPglNypTransposeBatch), kVecsPerCacheline);
        const uint32_t vbatches_per_load = cachelines_avail / min_load_cl_req;
        if (unlikely(!vbatches_per_load)) {
          g_failed_alloc_attempt_size = min_load_cl_req * kCacheline;
          goto EigTgenoToPgen_ret_NOMEM;
        }
        cur_vidx_ct = vbatches_per_load * kPglNypTransposeBatch;
      }
      if (cur_vidx_ct > variant_ct) {
        // Possible with --chr.  Avoid loading way too much in this case.
        // (Further optimizations possible here, but this should be a rare case
        // so I'm just trying to get down to the right order of magnitude.)
        cur_vidx_ct = MINV(RoundUpPow2(variant_ct, kPglNypTransposeBatch), raw_variant_ct);
      }
      cur_vidx_ctaw2 = NypCtToAlignedWordCt(cur_vidx_ct);
      if (unlikely(bigstack_alloc_w(sample_ct * cur_vidx_ctaw2, &tgeno_loadbuf))) {
        // this shouldn't be possible
        assert(0);
        goto EigTgenoToPgen_ret_NOMEM;
      }
    }

    uint32_t cur_vidx_ct4 = NypCtToByteCt(cur_vidx_ct);
    const uint32_t raw_load_batch_ct = raw_load_batch_ct_m1 + 1;
    recode_ctx.sample_ct = sample_ct;
    recode_ctx.raw_load_batch_size = raw_load_batch_size;
    recode_ctx.raw_load_batch_ct = raw_load_batch_ct;
    recode_ctx.cur_vidx_ctaw2 = cur_vidx_ctaw2;
    recode_ctx.tgeno_loadbuf = tgeno_loadbuf;
    const uint32_t thread_ct_m1 = MINV(max_thread_ct, raw_load_batch_ct) - 1;
    const uint32_t recode_launch_batch_idx = (raw_load_batch_ct * 2LLU * thread_ct_m1) / (2 * thread_ct_m1 + 1);
    if (unlikely(SetThreadCt0(thread_ct_m1, &recode_tg))) {
      goto EigTgenoToPgen_ret_NOMEM;
    }
    SetThreadFuncAndData(EigTgenoRecodeThread, &recode_ctx, &recode_tg);
    const uintptr_t transpose_block_ct_m1 = (sample_ct - 1) / kPglNypTransposeBatch;
    const uint32_t raw_pass_ct = 1 + (raw_variant_ct - 1) / cur_vidx_ct;
    uint32_t pass_ct = raw_pass_ct;
    if ((variant_ct < raw_variant_ct) && (raw_pass_ct > 1)) {
      // may be able to skip some raw passes.
      const uint32_t raw_pass_ct_m1 = raw_pass_ct - 1;
      assert(cur_vidx_ct % kBitsPerWord == 0);
      const uintptr_t* variant_include_iter = variant_include;
      const uintptr_t* variant_include_end = &(variant_include[BitCtToWordCt(raw_variant_ct)]);
      uint32_t word_ct = cur_vidx_ct / kBitsPerWord;
      pass_ct = 0;
      for (uint32_t pass_uidx = 0; ; ++pass_uidx) {
        if (pass_uidx >= raw_pass_ct_m1) {
          if (pass_uidx > raw_pass_ct_m1) {
            break;
          }
          word_ct = variant_include_end - variant_include_iter;
        }
        if (!AllWordsAreZero(variant_include_iter, word_ct)) {
          ++pass_ct;
        }
        variant_include_iter = &(variant_include_iter[word_ct]);
      }
    }
    uint32_t cur_vidx_base = 0;
    uint32_t pct = 0;
    for (uint32_t pass_idx1 = 1; ; cur_vidx_base += cur_vidx_ct) {
      if (raw_variant_ct - cur_vidx_base <= cur_vidx_ct) {
        cur_vidx_ct = raw_variant_ct - cur_vidx_base;
        cur_vidx_ct4 = NypCtToByteCt(cur_vidx_ct);
        cur_vidx_ctaw2 = NypCtToAlignedWordCt(cur_vidx_ct);
        recode_ctx.cur_vidx_ctaw2 = cur_vidx_ctaw2;
      }
      if (AllWordsAreZero(&(variant_include[cur_vidx_base / kBitsPerWord]), BitCtToWordCt(cur_vidx_ct))) {
        continue;
      }
      uint32_t cur_raw_load_batch_size = raw_load_batch_size;
      uintptr_t* tgeno_loadbuf_iter = tgeno_loadbuf;
      putc_unlocked('\r', stdout);
      printf("--eiggeno pass %u/%u: loading and recoding... 0%%", pass_idx1, pass_ct);
      fflush(stdout);
      pct = 0;
      uint32_t next_print_idx = (raw_load_batch_ct + 99) / 100;
      const uint64_t seek_addl_offset = 48 + cur_vidx_base / 4;
      for (uint32_t raw_load_batch_idx = 0; ; ) {
        if ((raw_load_batch_idx == recode_launch_batch_idx) && thread_ct_m1) {
          if (pass_idx1 == pass_ct) {
            DeclareLastThreadBlock(&recode_tg);
          }
          if (unlikely(SpawnThreads(&recode_tg))) {
            goto EigTgenoToPgen_ret_THREAD_CREATE_FAIL;
          }
        }
        // possible todo: check if multithreaded reading is faster if we know
        // we're reading from a SSD
        if (raw_load_batch_size == 1) {
          if (unlikely(fseeko(genofile, seek_addl_offset + raw_load_batch_idx * S_CAST(uint64_t, input_rec_blen), SEEK_SET) ||
                       (!fread_unlocked(tgeno_loadbuf_iter, cur_vidx_ct4, 1, genofile)))) {
            goto EigTgenoToPgen_ret_READ_FAIL;
          }
          tgeno_loadbuf_iter = &(tgeno_loadbuf_iter[cur_vidx_ctaw2]);
        } else {
          if (unlikely(!fread_unlocked(raw_loadbuf, cur_raw_load_batch_size * input_rec_blen, 1, genofile))) {
            goto EigTgenoToPgen_ret_READ_FAIL;
          }
          unsigned char* raw_loadbuf_iter = &(raw_loadbuf[cur_vidx_base / 4]);
          for (uint32_t uii = 0; uii != cur_raw_load_batch_size; ++uii) {
            memcpy(tgeno_loadbuf_iter, raw_loadbuf_iter, cur_vidx_ct4);
            raw_loadbuf_iter = &(raw_loadbuf_iter[input_rec_blen]);
            tgeno_loadbuf_iter = &(tgeno_loadbuf_iter[cur_vidx_ctaw2]);
          }
        }
        ++raw_load_batch_idx;
        if (raw_load_batch_idx >= raw_load_batch_ct_m1) {
          if (raw_load_batch_idx > raw_load_batch_ct_m1) {
            break;
          }
          cur_raw_load_batch_size = sample_ct - raw_load_batch_idx * raw_load_batch_size;
        }
        if (raw_load_batch_idx >= next_print_idx) {
          if (pct > 10) {
            putc_unlocked('\b', stdout);
          }
          pct = (raw_load_batch_idx * 100LLU) / raw_load_batch_ct;
          printf("\b\b%u%%", pct++);
          fflush(stdout);
          next_print_idx = (pct * S_CAST(uint64_t, raw_load_batch_ct) + 99) / 100;
        }
      }
      EigTgenoRecodeMain(thread_ct_m1, thread_ct_m1 + 1, &recode_ctx);
      JoinThreads0(&recode_tg);
      putc_unlocked('\r', stdout);
      printf("--eiggeno pass %u/%u: transposing and compressing... 0%%", pass_idx1, pass_ct);
      fflush(stdout);
      pct = 0;
      if (!mpgwp) {
        uintptr_t* read_iter = tgeno_loadbuf;
        next_print_idx = (cur_vidx_ct + 99) / 100;
        for (uint32_t write_uidx = 0; write_uidx < cur_vidx_ct; ) {
          const uintptr_t* read_iter2 = read_iter;
          const uint32_t vblock_size = MINV(kPglNypTransposeBatch, cur_vidx_ct - write_uidx);
          const uint32_t vblock_start = cur_vidx_base + write_uidx;
          const uint32_t vblock_end = vblock_start + vblock_size;
          // With --chr, we may skip many batches.
          for (uint32_t variant_uidx = vblock_start; ; ) {
            variant_uidx = AdvBoundedTo1Bit(variant_include, variant_uidx, vblock_end);
            if (variant_uidx == vblock_end) {
              break;
            }
            uint32_t read_batch_size = kPglNypTransposeBatch;
            for (uintptr_t transpose_block_idx = 0; ; ++transpose_block_idx) {
              if (transpose_block_idx >= transpose_block_ct_m1) {
                if (transpose_block_idx > transpose_block_ct_m1) {
                  break;
                }
                read_batch_size = ModNz(sample_ct, kPglNypTransposeBatch);
              }
              TransposeNypblock(read_iter2, cur_vidx_ctaw2, sample_ctaw2, read_batch_size, vblock_size, &(write_genovec[transpose_block_idx * kPglNypTransposeWords]), vecaligned_buf);
              read_iter2 = &(read_iter2[kPglNypTransposeBatch * cur_vidx_ctaw2]);
            }
            uintptr_t* cur_write_genovec = &(write_genovec[(variant_uidx - vblock_start) * sample_ctaw2]);
            const uint32_t variant_uidx_stop = AdvBoundedTo0Bit(variant_include, variant_uidx, vblock_end);
            for (; variant_uidx != variant_uidx_stop; ++variant_uidx, cur_write_genovec = &(cur_write_genovec[sample_ctaw2])) {
              ZeroTrailingNyps(sample_ct, cur_write_genovec);
              reterr = SpgwAppendBiallelicGenovec(cur_write_genovec, &spgw);
              if (unlikely(reterr)) {
                goto EigTgenoToPgen_ret_WRITE_FAIL;
              }
            }
          }
          if (write_uidx >= next_print_idx) {
            if (pct > 10) {
              putc_unlocked('\b', stdout);
            }
            pct = (write_uidx * 100LLU) / cur_vidx_ct;
            printf("\b\b%u%%", pct++);
            fflush(stdout);
            next_print_idx = (pct * S_CAST(uint64_t, cur_vidx_ct) + 99) / 100;
          }
          write_uidx += vblock_size;
          read_iter = &(read_iter[kPglNypTransposeWords]);
        }
      } else {
        const uint32_t vblock_group_ct = 1 + ((variant_ct - 1) / (write_thread_ct * kPglVblockSize));
        next_print_idx = (vblock_group_ct + 99) / 100;
        write_ctx.cur_block_write_ct = write_thread_ct * kPglVblockSize;
        for (uint32_t load_idx = 0; load_idx != vblock_group_ct; ++load_idx) {
          if (load_idx >= next_print_idx) {
            if (pct > 10) {
              putc_unlocked('\b', stdout);
            }
            pct = (load_idx * 100LLU) / vblock_group_ct;
            printf("\b\b%u%%", pct++);
            fflush(stdout);
            next_print_idx = (pct * S_CAST(uint64_t, vblock_group_ct) + 99) / 100;
          }
          write_ctx.tgeno_loadbuf_iter = &(tgeno_loadbuf[load_idx * write_thread_ct * (kPglVblockSize / kBitsPerWordD2)]);
          if (load_idx == vblock_group_ct - 1) {
            DeclareLastThreadBlock(&write_tg);
            write_ctx.cur_block_write_ct = variant_ct - load_idx * write_thread_ct * kPglVblockSize;
          }
          if (unlikely(SpawnThreads(&write_tg))) {
            goto EigTgenoToPgen_ret_THREAD_CREATE_FAIL;
          }
          EigTgenoTransposeMain(write_thread_ct - 1, &write_ctx);
          JoinThreads(&write_tg);
          reterr = MpgwFlush(mpgwp);
          if (unlikely(reterr)) {
            goto EigTgenoToPgen_ret_WRITE_FAIL;
          }
        }
      }
      if (pass_idx1 == pass_ct) {
        break;
      }
      ++pass_idx1;
      fputs("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b                     ", stdout);
      if (unlikely(fseeko(genofile, 48, SEEK_SET))) {
        goto EigTgenoToPgen_ret_READ_FAIL;
      }
    }
    if (!mpgwp) {
      reterr = SpgwFinish(&spgw);
      if (unlikely(reterr)) {
        goto EigTgenoToPgen_ret_1;
      }
    } else {
      mpgwp = nullptr;
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\bdone.\n", stdout);
  }
  while (0) {
  EigTgenoToPgen_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  EigTgenoToPgen_ret_READ_FAIL:
    if (feof_unlocked(genofile)) {
      errno = 0;
    }
    logputs("\n");
    logerrprintfww(kErrprintfFread, genoname, rstrerror(errno));
    reterr = kPglRetReadFail;
    break;
  EigTgenoToPgen_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  EigTgenoToPgen_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 EigTgenoToPgen_ret_1:
  CleanupThreads(&write_tg);
  CleanupThreads(&recode_tg);
  CleanupMpgw(mpgwp, &reterr);
  CleanupSpgw(&spgw, &reterr);
  BigstackReset(bigstack_mark);
  return reterr;
}

// Move this to a more central location if we ever need it elsewhere.
const char* ScanadvHexU32(const char* hex_start, uint32_t* valp) {
  // Mirror scanf("%x"), except prohibit leading -, and error out on overflow.
  if (*hex_start == '+') {
    ++hex_start;
  }
  if ((hex_start[0] == '0') && ((hex_start[1] & 0xdf) == 'X')) {
    hex_start = &(hex_start[2]);
  }
  uint32_t val = 0;
  for (const char* hex_iter = hex_start; ; ++hex_iter) {
    uint32_t cur_digit = ctou32(*hex_iter) - 48;
    if (cur_digit >= 10) {
      // 'A'/'a' = 0, 'B'/'b' = 1, etc.
      const uint32_t letter_idx = (cur_digit & 0xffffffdfU) - 17;
      if (letter_idx > 5) {
        *valp = val;
        return (hex_iter == hex_start)? nullptr : hex_iter;
      }
      cur_digit = letter_idx + 10;
    }
    if (unlikely(val >= (1U << 28))) {
      // overflow
      return nullptr;
    }
    val = (val * 16) + cur_digit;
  }
}

PglErr EigfileToPgen(const char* genoname, const char* indname, const char* snpname, const char* const_fid, const char* missing_catname, MiscFlags misc_flags, ImportFlags import_flags, LoadFilterLogFlags load_filter_log_import_flags, uint32_t psam_01, char id_delim, uint32_t max_thread_ct, char* outname, char* outname_end, ChrInfo* cip) {
  FILE* genofile = nullptr;
  PglErr reterr = kPglRetSuccess;
  {
    uint32_t sample_ct;
    uint32_t s_hash;
    reterr = EigIndToPsam(indname, const_fid, missing_catname, misc_flags, import_flags, psam_01, id_delim, outname, outname_end, &sample_ct, &s_hash);
    if (unlikely(reterr)) {
      goto EigfileToPgen_ret_1;
    }
    FinalizeChrset(load_filter_log_import_flags, cip);

    uintptr_t* variant_include;
    uint32_t raw_variant_ct;
    uint32_t v_hash;
    reterr = EigSnpToPvar(snpname, cip, import_flags, max_thread_ct, outname, outname_end, &variant_include, &raw_variant_ct, &v_hash);
    if (unlikely(reterr)) {
      goto EigfileToPgen_ret_1;
    }
    const uint32_t variant_ct = PopcountWords(variant_include, BitCtToWordCt(raw_variant_ct));
    if (unlikely(!variant_ct)) {
      char* write_iter = strcpya_k(g_logbuf, "Error: All ");
      write_iter = u32toa(raw_variant_ct, write_iter);
      write_iter = strcpya_k(write_iter, " variant");
      if (raw_variant_ct != 1) {
        *write_iter++ = 's';
      }
      write_iter = strcpya_k(write_iter, " in .snp file excluded by ");
      AppendLoadFilterFlagnames(load_filter_log_import_flags, &write_iter);
      strcpy_k(write_iter, ".\n");
      WordWrapB(0);
      logputsb();
      goto EigfileToPgen_ret_INCONSISTENT_INPUT;
    }
    {
      char* write_iter = strcpya_k(g_logbuf, "--eigsnp: ");
      write_iter = u32toa(raw_variant_ct, write_iter);
      write_iter = strcpya_k(write_iter, " variant");
      if (raw_variant_ct != 1) {
        *write_iter++ = 's';
      }
      if (variant_ct != raw_variant_ct) {
        write_iter = strcpya_k(write_iter, " scanned; ");
        write_iter = u32toa(raw_variant_ct - variant_ct, write_iter);
        write_iter = strcpya_k(write_iter, " excluded by ");
        AppendLoadFilterFlagnames(load_filter_log_import_flags, &write_iter);
        write_iter = strcpya_k(write_iter, ", ");
        write_iter = u32toa(variant_ct, write_iter);
      }
      write_iter = strcpya_k(write_iter, " imported to ");
      write_iter = strcpyax(write_iter, outname, ' ');
      if (load_filter_log_import_flags && (variant_ct == raw_variant_ct)) {
        *write_iter++ = '(';
        AppendLoadFilterFlagnames(load_filter_log_import_flags, &write_iter);
        write_iter = strcpya_k(write_iter, " had no effect)");
      }
      strcpy_k(write_iter, ".\n");
      WordWrapB(0);
      logputsb();
    }

    // open genofile, check header
    if (unlikely(fopen_checked(genoname, FOPEN_RB, &genofile))) {
      goto EigfileToPgen_ret_OPEN_FAIL;
    }
    if (unlikely(fseeko(genofile, 0, SEEK_END))) {
      goto EigfileToPgen_ret_READ_FAIL;
    }
    const uint64_t fsize = ftello(genofile);
    if (unlikely(fsize < 49)) {
      snprintf(g_logbuf, kLogbufSize, "Error: %s is too small to be a valid EIGENSOFT PACKEDANCESTRYMAP or TGENO file.\n", genoname);
      goto EigfileToPgen_ret_INCONSISTENT_INPUT_WW;
    }
    rewind(genofile);
    char header[48];
    if (unlikely(!fread_unlocked(header, 48, 1, genofile))) {
      goto EigfileToPgen_ret_READ_FAIL;
    }
    uint32_t is_tgeno = 0;
    if ((!memequal(header, "GENO", 4)) || (ctou32(header[4]) > 32)) {
      if (unlikely((!memequal(header, "TGENO", 5)) || (ctou32(header[5]) > 32))) {
        snprintf(g_logbuf, kLogbufSize, "Error: %s does not start with 'GENO ' or 'TGENO '.\n", genoname);
        goto EigfileToPgen_ret_INCONSISTENT_INPUT_WW;
      }
      is_tgeno = 1;
    }
    char* sample_ct_start = &(header[5]);
    if (unlikely(!memchr(sample_ct_start, 0, 43))) {
      snprintf(g_logbuf, kLogbufSize, "Error: %s header is not null-terminated.\n", genoname);
      goto EigfileToPgen_ret_INCONSISTENT_INPUT_WW;
    }
    sample_ct_start = FirstNonTspace(sample_ct_start);
    char* sample_ct_end = FirstSpaceOrEoln(sample_ct_start);
    char* variant_ct_start = FirstNonTspace(sample_ct_end);
    char* variant_ct_end = FirstSpaceOrEoln(variant_ct_start);
    char* sample_hash_start = FirstNonTspace(variant_ct_end);
    char* sample_hash_end = FirstSpaceOrEoln(sample_hash_start);
    char* variant_hash_start = FirstNonTspace(sample_hash_end);
    if (unlikely(IsEolnKns(*variant_hash_start))) {
      snprintf(g_logbuf, kLogbufSize, "Error: %s header has fewer tokens than expected.\n", genoname);
      goto EigfileToPgen_ret_INCONSISTENT_INPUT_WW;
    }
    char* variant_hash_end = FirstSpaceOrEoln(variant_hash_start);
    if (unlikely(!IsEolnKns(*FirstNonTspace(variant_hash_end)))) {
      snprintf(g_logbuf, kLogbufSize, "Error: %s header has more tokens than expected.\n", genoname);
      goto EigfileToPgen_ret_INCONSISTENT_INPUT_WW;
    }

    uint32_t header_val;
    if (unlikely(ScanUintDefcap(sample_ct_start, &header_val) || (header_val != sample_ct))) {
      *sample_ct_end = '\0';
      snprintf(g_logbuf, kLogbufSize, "Error: Sample count mismatch between %s header (%s) and .ind file.\n", genoname, sample_ct_start);
      goto EigfileToPgen_ret_INCONSISTENT_INPUT_WW;
    }
    if (unlikely(ScanUintDefcap(variant_ct_start, &header_val) || (header_val != raw_variant_ct))) {
      *sample_ct_end = '\0';
      snprintf(g_logbuf, kLogbufSize, "Error: Variant count mismatch between %s header (%s) and .snp file.\n", genoname, variant_ct_start);
      goto EigfileToPgen_ret_INCONSISTENT_INPUT_WW;
    }
    if (!(import_flags & kfImportEigNohash)) {
      uint32_t observed_hash;
      if (unlikely((!ScanadvHexU32(sample_hash_start, &observed_hash)) || (observed_hash != s_hash))) {
        *sample_hash_end = '\0';
        snprintf(g_logbuf, kLogbufSize, "Error: Sample-ID hash mismatch between %s header (%s) and .ind file (%x). (If you intentionally modified the IDs, use the 'nohash' modifier to suppress this error.)\n", genoname, sample_hash_start, s_hash);
        goto EigfileToPgen_ret_INCONSISTENT_INPUT_WW;
      }
      if (unlikely((!ScanadvHexU32(variant_hash_start, &observed_hash)) || (observed_hash != v_hash))) {
        *variant_hash_end = '\0';
        snprintf(g_logbuf, kLogbufSize, "Error: Variant-ID hash mismatch between %s header (%s) and .snp file (%x). (If you intentionally modified the IDs, use the 'nohash' modifier to suppress this error.)\n", genoname, variant_hash_start, v_hash);
        goto EigfileToPgen_ret_INCONSISTENT_INPUT_WW;
      }
    }
    snprintf(outname_end, kMaxOutfnameExtBlen, ".pgen");
    if (!is_tgeno) {
      const uintptr_t input_rec_blen = MAXV(48, DivUp(sample_ct, 4));
      const uint64_t fsize_expected = (raw_variant_ct + 1LLU) * input_rec_blen;
      if (fsize != fsize_expected) {
        snprintf(g_logbuf, kLogbufSize, "Error: Unexpected %s file size (%u sample%s from .ind file, %u variant%s from .snp file, thus %" PRIu64 " bytes expected; observed %" PRIu64 ").\n", genoname, sample_ct, (sample_ct == 1)? "" : "s", raw_variant_ct, (raw_variant_ct == 1)? "" : "s", fsize_expected, fsize);
        goto EigfileToPgen_ret_INCONSISTENT_INPUT_WW;
      }
      reterr = EigGenoToPgen(outname, genoname, variant_include, variant_ct, sample_ct, genofile);
    } else {
      const uint64_t input_rec_blen = MAXV(48, DivUp(raw_variant_ct, 4));
      const uint64_t fsize_expected = 48 + sample_ct * input_rec_blen;
      if (fsize != fsize_expected) {
        snprintf(g_logbuf, kLogbufSize, "Error: Unexpected %s file size (%u sample%s from .ind file, %u variant%s from .snp file, thus %" PRIu64 " bytes expected; observed %" PRIu64 ").\n", genoname, sample_ct, (sample_ct == 1)? "" : "s", raw_variant_ct, (raw_variant_ct == 1)? "" : "s", fsize_expected, fsize);
        goto EigfileToPgen_ret_INCONSISTENT_INPUT_WW;
      }
      reterr = EigTgenoToPgen(outname, genoname, variant_include, raw_variant_ct, variant_ct, sample_ct, max_thread_ct, genofile);
    }
    if (unlikely(reterr)) {
      goto EigfileToPgen_ret_1;
    }

    putc_unlocked('\r', stdout);
    logprintfww("--eiggeno: %s written.\n", outname);
  }
  while (0) {
  EigfileToPgen_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  EigfileToPgen_ret_READ_FAIL:
    if (feof_unlocked(genofile)) {
      errno = 0;
    }
    putc_unlocked('\n', stdout);
    logerrprintfww(kErrprintfFread, genoname, rstrerror(errno));
    reterr = kPglRetReadFail;
    break;
  EigfileToPgen_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  EigfileToPgen_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 EigfileToPgen_ret_1:
  fclose_cond(genofile);
  return reterr;
}

static_assert(sizeof(Dosage) == 2, "GenerateDummy() needs to be updated.");
PglErr GenerateDummy(const GenDummyInfo* gendummy_info_ptr, ImportFlags import_flags, LoadFilterLogFlags load_filter_log_import_flags, uint32_t psam_01, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, uint32_t max_thread_ct, sfmt_t* sfmtp, char* outname, char* outname_end, ChrInfo* cip) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* psamfile = nullptr;
  char* pvar_cswritep = nullptr;
  CompressStreamState pvar_css;
  PreinitCstream(&pvar_css);
  ThreadGroup tg;
  PreinitThreads(&tg);
  STPgenWriter spgw;
  PglErr reterr = kPglRetSuccess;
  PreinitSpgw(&spgw);
  {
    FinalizeChrset(load_filter_log_import_flags, cip);
    if (unlikely(!IsSet(cip->chr_mask, 1))) {
      logerrputs("Error: --dummy cannot be used when chromosome 1 is excluded.\n");
      goto GenerateDummy_ret_INVALID_CMDLINE;
    }
    if (unlikely(IsSet(cip->haploid_mask, 1))) {
      logerrputs("Error: --dummy cannot be used to generate haploid data.\n");
      goto GenerateDummy_ret_INVALID_CMDLINE;
    }
    char chr1_name_buf[5];
    char* chr1_name_end = chrtoa(cip, 1, chr1_name_buf);
    *chr1_name_end = '\t';
    const uint32_t chr1_name_blen = 1 + S_CAST(uintptr_t, chr1_name_end - chr1_name_buf);
    const uint32_t sample_ct = gendummy_info_ptr->sample_ct;
    const uint32_t variant_ct = gendummy_info_ptr->variant_ct;
    // missing pheno string is always "NA"
    const GenDummyFlags flags = gendummy_info_ptr->flags;
    uint16_t alleles[13];
    uint32_t four_alleles = 0;
    if (flags & kfGenDummyAcgt) {
      memcpy_k(alleles, "\tA\tC\tA\tG\tA\tT\tC\tG\tC\tT\tG\tT\tA", 26);
      four_alleles = 1;
    } else if (flags & kfGenDummy1234) {
      memcpy_k(alleles, "\t1\t2\t1\t3\t1\t4\t2\t3\t2\t4\t3\t4\t1", 26);
      four_alleles = 1;
    } else if (flags & kfGenDummy12) {
      memcpy_k(alleles, "\t1\t2\t1", 6);
    } else {
      memcpy_k(alleles, "\tA\tB\tA", 6);
    }

    snprintf(outname_end, kMaxOutfnameExtBlen, ".pvar");
    const uint32_t output_zst = (import_flags / kfImportKeepAutoconvVzs) & 1;
    if (output_zst) {
      snprintf(&(outname_end[5]), kMaxOutfnameExtBlen - 5, ".zst");
    }
    const uintptr_t overflow_buf_size = 2 * kCompressStreamBlock;
    reterr = InitCstreamAlloc(outname, 0, output_zst, max_thread_ct, overflow_buf_size, &pvar_css, &pvar_cswritep);
    if (unlikely(reterr)) {
      goto GenerateDummy_ret_1;
    }
    pvar_cswritep = strcpya_k(pvar_cswritep, "#CHROM\tPOS\tID\tREF\tALT");
    AppendBinaryEoln(&pvar_cswritep);
    if (four_alleles) {
      uint32_t urand = 0;
      for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
        if (!(variant_idx % 8)) {
          if (unlikely(Cswrite(&pvar_css, &pvar_cswritep))) {
            goto GenerateDummy_ret_WRITE_FAIL;
          }
          do {
            urand = sfmt_genrand_uint32(sfmtp);
          } while (urand < 425132032U);  // 2^32 - 12^8
        }
        const uint32_t quotient = urand / 12;
        const uint32_t remainder = urand - (quotient * 12U);
        urand = quotient;
        pvar_cswritep = memcpya(pvar_cswritep, chr1_name_buf, chr1_name_blen);
        pvar_cswritep = u32toa(variant_idx, pvar_cswritep);
        pvar_cswritep = strcpya_k(pvar_cswritep, "\tsnp");
        pvar_cswritep = u32toa(variant_idx, pvar_cswritep);
        pvar_cswritep = memcpya(pvar_cswritep, &(alleles[remainder]), 4);
        AppendBinaryEoln(&pvar_cswritep);
      }
    } else {
      uint32_t urand = 0;
      for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
        if (!(variant_idx % 32)) {
          if (unlikely(Cswrite(&pvar_css, &pvar_cswritep))) {
            goto GenerateDummy_ret_WRITE_FAIL;
          }
          urand = sfmt_genrand_uint32(sfmtp);
        }
        const uint32_t remainder = urand & 1;
        urand >>= 1;
        pvar_cswritep = memcpya(pvar_cswritep, chr1_name_buf, chr1_name_blen);
        pvar_cswritep = u32toa(variant_idx, pvar_cswritep);
        pvar_cswritep = strcpya_k(pvar_cswritep, "\tsnp");
        pvar_cswritep = u32toa(variant_idx, pvar_cswritep);
        pvar_cswritep = memcpya(pvar_cswritep, &(alleles[remainder]), 4);
        AppendBinaryEoln(&pvar_cswritep);
      }
    }
    if (unlikely(CswriteCloseNull(&pvar_css, pvar_cswritep))) {
      goto GenerateDummy_ret_WRITE_FAIL;
    }
    BigstackReset(bigstack_mark);

    snprintf(outname_end, kMaxOutfnameExtBlen, ".psam");
    if (unlikely(fopen_checked(outname, FOPEN_WB, &psamfile))) {
      goto GenerateDummy_ret_OPEN_FAIL;
    }
    const uint32_t pheno_ct = gendummy_info_ptr->pheno_ct;
    char* writebuf;
    if (unlikely(bigstack_alloc_c(kMaxMediumLine + 48 + pheno_ct * MAXV(kMaxMissingPhenostrBlen, 16), &writebuf))) {
      goto GenerateDummy_ret_NOMEM;
    }
    char* writebuf_flush = &(writebuf[kMaxMediumLine]);
    // Alpha 2 change: no more FID column
    char* write_iter = strcpya_k(writebuf, "#IID\tSEX");
    for (uint32_t pheno_idx_p1 = 1; pheno_idx_p1 <= pheno_ct; ++pheno_idx_p1) {
      write_iter = strcpya_k(write_iter, "\tPHENO");
      write_iter = u32toa(pheno_idx_p1, write_iter);
    }
    AppendBinaryEoln(&write_iter);
    const uint32_t pheno_m_check = (gendummy_info_ptr->pheno_mfreq >= k2m32 * 0.5);
    const uint32_t pheno_m32 = S_CAST(uint32_t, gendummy_info_ptr->pheno_mfreq * 4294967296.0 - 0.5);
    if ((flags & kfGenDummyScalarPheno) && pheno_ct) {
      uint32_t saved_rnormal = 0;
      double saved_rnormal_val = 0.0;
      for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        if (unlikely(fwrite_ck(writebuf_flush, psamfile, &write_iter))) {
          goto GenerateDummy_ret_WRITE_FAIL;
        }
        write_iter = strcpya_k(write_iter, "per");
        write_iter = u32toa(sample_idx, write_iter);
        // could add option to add some males/unknown gender
        write_iter = strcpya_k(write_iter, "\t2");
        for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
          *write_iter++ = '\t';
          if (pheno_m_check && (sfmt_genrand_uint32(sfmtp) <= pheno_m32)) {
            write_iter = strcpya_k(write_iter, "NA");
          } else {
            double dxx;
            if (saved_rnormal) {
              dxx = saved_rnormal_val;
            } else {
              dxx = RandNormal(sfmtp, &saved_rnormal_val);
            }
            saved_rnormal_val = 1 - saved_rnormal_val;
            write_iter = dtoa_g(dxx, write_iter);
          }
        }
        AppendBinaryEoln(&write_iter);
      }
    } else {
      uint32_t urand = sfmt_genrand_uint32(sfmtp);
      uint32_t urand_bits_left = 32;
      const char ctrl_char = '1' - psam_01;
      for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        if (unlikely(fwrite_ck(writebuf_flush, psamfile, &write_iter))) {
          goto GenerateDummy_ret_WRITE_FAIL;
        }
        // bugfix (9 Mar 2018): forgot to remove FID column here
        write_iter = strcpya_k(write_iter, "per");
        write_iter = u32toa(sample_idx, write_iter);
        write_iter = strcpya_k(write_iter, "\t2");
        for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
          *write_iter++ = '\t';
          if (pheno_m_check && (sfmt_genrand_uint32(sfmtp) <= pheno_m32)) {
            write_iter = strcpya_k(write_iter, "NA");
          } else {
            if (!urand_bits_left) {
              urand = sfmt_genrand_uint32(sfmtp);
              urand_bits_left = 32;
            }
            *write_iter++ = (urand & 1) + ctrl_char;
            urand >>= 1;
            --urand_bits_left;
          }
        }
        AppendBinaryEoln(&write_iter);
      }
    }
    if (unlikely(fclose_flush_null(writebuf_flush, write_iter, &psamfile))) {
      goto GenerateDummy_ret_WRITE_FAIL;
    }

    BigstackReset(writebuf);
    snprintf(outname_end, kMaxOutfnameExtBlen, ".pgen");
    GenerateDummyCtx ctx;
    {
      const uint32_t geno_mfreq_ct = gendummy_info_ptr->geno_mfreq_ct;
      ctx.geno3_thresh_ct = geno_mfreq_ct;
      ctx.geno3_thresh_arr = nullptr;
      if (geno_mfreq_ct) {
        if (unlikely(geno_mfreq_ct > 65536)) {
          logerrputs("Error: --dummy is limited to 65536 missing-dosage frequencies.\n");
          goto GenerateDummy_ret_INVALID_CMDLINE;
        }
        uint32_t warning_9999999 = 0;
        if (unlikely(bigstack_alloc_u64(geno_mfreq_ct, &ctx.geno3_thresh_arr))) {
          goto GenerateDummy_ret_NOMEM;
        }
        const double* geno_mfreqs = gendummy_info_ptr->geno_mfreqs;
        for (uint32_t uii = 0; uii != geno_mfreq_ct; ++uii) {
          double nm_freq = 1.0 - geno_mfreqs[uii];
          if ((nm_freq < 1e-7) && (nm_freq > 0.0)) {
            warning_9999999 = 1;
            nm_freq = 1e-7;
          }
          ctx.geno3_thresh_arr[uii] = S_CAST(int64_t, nm_freq * 4294967296.0);
        }
        if (warning_9999999) {
          logerrputs("Warning: reducing missing dosage freq(s) in (0.9999999, 1) to 0.9999999.\n");
        }
      }
    }
    const double phase_freq = gendummy_info_ptr->phase_freq;
    ctx.phase_invert = 0;
    const uint32_t is_phased = (phase_freq >= k2m53);
    if (!is_phased) {
      ctx.phase_geomdist[kBitsPerWordD2 - 1] = 0;
    } else {
      double remaining_prob = 1.0;
      ctx.phase_invert = (phase_freq > 0.5);
      if (ctx.phase_invert) {
        for (uint32_t uii = 0; uii != kBitsPerWordD2; ++uii) {
          remaining_prob *= phase_freq;
          ctx.phase_geomdist[uii] = -S_CAST(uint64_t, remaining_prob * k2p64);
        }
      } else {
        const double phase_nfreq = 1.0 - phase_freq;
        for (uint32_t uii = 0; uii != kBitsPerWordD2; ++uii) {
          remaining_prob *= phase_nfreq;
          ctx.phase_geomdist[uii] = -S_CAST(uint64_t, remaining_prob * k2p64);
        }
      }
    }
    const double dosage_nfreq = 1.0 - gendummy_info_ptr->dosage_freq;
    if (dosage_nfreq >= 1.0) {
      ctx.dosage_geomdist_max = kBitsPerWord;  // used as a flag
    } else {
      double remaining_prob = 1.0;
      for (uint32_t uii = 0; uii != kBitsPerWordD2; ++uii) {
        remaining_prob *= dosage_nfreq;
        ctx.dosage_geomdist[uii] = -S_CAST(uint64_t, remaining_prob * k2p64);
      }
      uint32_t dosage_geomdist_max = kBitsPerWordD2;
      for (; dosage_geomdist_max; --dosage_geomdist_max) {
        if (ctx.dosage_geomdist[dosage_geomdist_max - 1] != 0) {
          break;
        }
      }
      ctx.dosage_geomdist_max = dosage_geomdist_max;
    }
    PgenGlobalFlags gflags = kfPgenGlobal0;
    if (is_phased) {
      gflags |= kfPgenGlobalHardcallPhasePresent;
    }
    if (dosage_nfreq < 1.0) {
      gflags |= kfPgenGlobalDosagePresent;
    }
    uintptr_t spgw_alloc_cacheline_ct;
    uint32_t max_vrec_len;
    reterr = SpgwInitPhase1(outname, nullptr, nullptr, variant_ct, sample_ct, 0, kPgenWriteBackwardSeek, gflags, 1, &spgw, &spgw_alloc_cacheline_ct, &max_vrec_len);
    if (unlikely(reterr)) {
      if (reterr == kPglRetOpenFail) {
        logerrprintfww(kErrprintfFopen, outname, strerror(errno));
      }
      goto GenerateDummy_ret_1;
    }
    unsigned char* spgw_alloc;
    if (unlikely(bigstack_alloc_uc(spgw_alloc_cacheline_ct * kCacheline, &spgw_alloc))) {
      goto GenerateDummy_ret_NOMEM;
    }
    SpgwInitPhase2(max_vrec_len, &spgw, spgw_alloc);

    // thread-count-independent:
    //   (everything after "2 *" rounded up to cacheline)
    //   ctx.write_genovecs: 2 * sample_ctaw2 * sizeof(intptr_t) *
    //                       main_block_size
    //   ctx.write_dosage_cts: 2 * sizeof(int32_t) * main_block_size
    //   ctx.write_dosage_presents: 2 * sample_ctaw * sizeof(intptr_t) *
    //                              main_block_size
    //   ctx.write_dosage_mains: 2 * sample_ct * sizeof(Dosage)
    uint32_t calc_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
    // saturates around 4 compute threads, both with and without dosage
    // (todo: test this on something other than a MacBook Pro, could just be a
    // hyperthreading artifact)
    if (calc_thread_ct > 4) {
      calc_thread_ct = 4;
    }
    if (unlikely(InitAllocSfmtpArr(calc_thread_ct, 0, sfmtp, &ctx.sfmtp_arr))) {
      goto GenerateDummy_ret_NOMEM;
    }
    const uint32_t sample_ctaw2 = NypCtToAlignedWordCt(sample_ct);
    const uint32_t sample_ctaw = BitCtToAlignedWordCt(sample_ct);
    // we're making 18 allocations; be pessimistic re: rounding
    uintptr_t cachelines_avail_m18 = bigstack_left() / kCacheline;
    if (unlikely(cachelines_avail_m18 < 18)) {
      goto GenerateDummy_ret_NOMEM;
    }
    cachelines_avail_m18 -= 18;
    const uintptr_t bytes_req_per_in_block_variant = 2 * (sample_ctaw2 * sizeof(intptr_t) + 2 * sizeof(int32_t) + sample_ctaw * 4 * sizeof(intptr_t) + sample_ct * 2 * sizeof(Dosage));
    uintptr_t main_block_size = (cachelines_avail_m18 * kCacheline) / bytes_req_per_in_block_variant;
    if (main_block_size > 65536) {
      main_block_size = 65536;
    } else if (unlikely(main_block_size < 8)) {
      // this threshold is arbitrary
      goto GenerateDummy_ret_NOMEM;
    }
    if (calc_thread_ct > main_block_size / 8) {
      calc_thread_ct = main_block_size / 8;
    }
    if (unlikely(SetThreadCt(calc_thread_ct, &tg))) {
      goto GenerateDummy_ret_NOMEM;
    }
    ctx.sample_ct = sample_ct;
    if (unlikely(bigstack_alloc_w(sample_ctaw2 * main_block_size, &(ctx.write_genovecs[0])) ||
                 bigstack_alloc_w(sample_ctaw2 * main_block_size, &(ctx.write_genovecs[1])) ||
                 bigstack_alloc_w(sample_ctaw * main_block_size, &(ctx.write_phasepresents[0])) ||
                 bigstack_alloc_w(sample_ctaw * main_block_size, &(ctx.write_phasepresents[1])) ||
                 bigstack_alloc_w(sample_ctaw * main_block_size, &(ctx.write_phaseinfos[0])) ||
                 bigstack_alloc_w(sample_ctaw * main_block_size, &(ctx.write_phaseinfos[1])) ||
                 bigstack_alloc_u32(main_block_size, &(ctx.write_dosage_cts[0])) ||
                 bigstack_alloc_u32(main_block_size, &(ctx.write_dosage_cts[1])) ||
                 bigstack_alloc_w(sample_ctaw * main_block_size, &(ctx.write_dosage_presents[0])) ||
                 bigstack_alloc_w(sample_ctaw * main_block_size, &(ctx.write_dosage_presents[1])) ||
                 bigstack_alloc_dosage(sample_ct * main_block_size, &(ctx.write_dosage_mains[0])) ||
                 bigstack_alloc_dosage(sample_ct * main_block_size, &(ctx.write_dosage_mains[1])) ||
                 bigstack_alloc_u32(main_block_size, &(ctx.write_dphase_cts[0])) ||
                 bigstack_alloc_u32(main_block_size, &(ctx.write_dphase_cts[1])) ||
                 bigstack_alloc_w(sample_ctaw * main_block_size, &(ctx.write_dphase_presents[0])) ||
                 bigstack_alloc_w(sample_ctaw * main_block_size, &(ctx.write_dphase_presents[1])) ||
                 bigstack_alloc_dphase(sample_ct * main_block_size, &(ctx.write_dphase_deltas[0])) ||
                 bigstack_alloc_dphase(sample_ct * main_block_size, &(ctx.write_dphase_deltas[1])))) {
      // this should be impossible
      assert(0);
      goto GenerateDummy_ret_NOMEM;
    }
    // bugfix (3 Nov 2017): forgot to handle hard_call_thresh default value
    if (hard_call_thresh == UINT32_MAX) {
      hard_call_thresh = kDosageMid / 10;
    }
    ctx.hard_call_halfdist = kDosage4th - hard_call_thresh;
    ctx.dosage_erase_halfdist = kDosage4th - dosage_erase_thresh;
    SetThreadFuncAndData(GenerateDummyThread, &ctx, &tg);

    // Main workflow:
    // 1. Set n=0
    //
    // 2. Spawn threads generating batch n genotype data
    // 3. If n>0, write results for block (n-1)
    // 4. Increment n by 1
    // 5. Join threads
    // 6. Goto step 2 unless eof
    //
    // 7. Write results for last block
    uint32_t prev_block_write_ct = 0;
    uint32_t parity = 0;
    for (uint32_t vidx_start = 0; ; ) {
      uint32_t cur_block_write_ct = 0;
      if (!IsLastBlock(&tg)) {
        cur_block_write_ct = MINV(variant_ct - vidx_start, main_block_size);
      }
      if (vidx_start) {
        JoinThreads(&tg);
        // GenerateDummyThread() never errors out
      }
      if (!IsLastBlock(&tg)) {
        ctx.cur_block_write_ct = cur_block_write_ct;
        if (vidx_start + cur_block_write_ct == variant_ct) {
          DeclareLastThreadBlock(&tg);
        }
        if (unlikely(SpawnThreads(&tg))) {
          goto GenerateDummy_ret_THREAD_CREATE_FAIL;
        }
      }
      parity = 1 - parity;
      if (vidx_start) {
        // write *previous* block results
        uintptr_t* write_genovec_iter = ctx.write_genovecs[parity];
        uintptr_t* write_phasepresent_iter = ctx.write_phasepresents[parity];
        uintptr_t* write_phaseinfo_iter = ctx.write_phaseinfos[parity];
        uint32_t* write_dosage_ct_iter = ctx.write_dosage_cts[parity];
        uintptr_t* write_dosage_present_iter = ctx.write_dosage_presents[parity];
        Dosage* write_dosage_main_iter = ctx.write_dosage_mains[parity];
        uint32_t* write_dphase_ct_iter = ctx.write_dphase_cts[parity];
        uintptr_t* write_dphase_present_iter = ctx.write_dphase_presents[parity];
        SDosage* write_dphase_delta_iter = ctx.write_dphase_deltas[parity];
        for (uint32_t vidx = vidx_start - prev_block_write_ct; vidx != vidx_start; ++vidx) {
          const uint32_t cur_dosage_ct = *write_dosage_ct_iter++;
          const uint32_t cur_dphase_ct = *write_dphase_ct_iter++;
          if (!cur_dosage_ct) {
            if (!is_phased) {
              if (unlikely(SpgwAppendBiallelicGenovec(write_genovec_iter, &spgw))) {
                goto GenerateDummy_ret_WRITE_FAIL;
              }
            } else {
              if (unlikely(SpgwAppendBiallelicGenovecHphase(write_genovec_iter, write_phasepresent_iter, write_phaseinfo_iter, &spgw))) {
                goto GenerateDummy_ret_WRITE_FAIL;
              }
            }
          } else {
            if (!is_phased) {
              reterr = SpgwAppendBiallelicGenovecDosage16(write_genovec_iter, write_dosage_present_iter, write_dosage_main_iter, cur_dosage_ct, &spgw);
              if (unlikely(reterr)) {
                goto GenerateDummy_ret_1;
              }
            } else {
              reterr = SpgwAppendBiallelicGenovecDphase16(write_genovec_iter, write_phasepresent_iter, write_phaseinfo_iter, write_dosage_present_iter, write_dphase_present_iter, write_dosage_main_iter, write_dphase_delta_iter, cur_dosage_ct, cur_dphase_ct, &spgw);
              if (unlikely(reterr)) {
                goto GenerateDummy_ret_1;
              }
            }
          }
          write_genovec_iter = &(write_genovec_iter[sample_ctaw2]);
          write_phasepresent_iter = &(write_phasepresent_iter[sample_ctaw]);
          write_phaseinfo_iter = &(write_phaseinfo_iter[sample_ctaw]);
          write_dosage_present_iter = &(write_dosage_present_iter[sample_ctaw]);
          write_dosage_main_iter = &(write_dosage_main_iter[sample_ct]);
          write_dphase_present_iter = &(write_dphase_present_iter[sample_ctaw]);
          write_dphase_delta_iter = &(write_dphase_delta_iter[sample_ct]);
        }
      }
      if (vidx_start == variant_ct) {
        break;
      }
      if (vidx_start) {
        printf("\r--dummy: %uk variants written.", vidx_start / 1000);
        fflush(stdout);
      }
      vidx_start += cur_block_write_ct;
      prev_block_write_ct = cur_block_write_ct;
    }
    reterr = SpgwFinish(&spgw);
    if (unlikely(reterr)) {
      goto GenerateDummy_ret_1;
    }

    putc_unlocked('\r', stdout);
    *outname_end = '\0';
    logprintfww("Dummy data (%u sample%s, %u SNP%s) written to %s.pgen + %s.pvar%s + %s.psam .\n", sample_ct, (sample_ct == 1)? "" : "s", variant_ct, (variant_ct == 1)? "" : "s", outname, outname, output_zst? ".zst" : "", outname);
  }
  while (0) {
  GenerateDummy_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  GenerateDummy_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  GenerateDummy_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  GenerateDummy_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  GenerateDummy_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 GenerateDummy_ret_1:
  CleanupSpgw(&spgw, &reterr);
  CleanupThreads(&tg);
  fclose_cond(psamfile);
  CswriteCloseCond(&pvar_css, pvar_cswritep);
  BigstackReset(bigstack_mark);
  return reterr;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
