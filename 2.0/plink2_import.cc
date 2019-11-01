// This file is part of PLINK 2.00, copyright (C) 2005-2019 Shaun Purcell,
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


#include "pgenlib_write.h"
#include "plink2_import.h"
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
  gendummy_info_ptr->geno_mfreq = 0.0;
  gendummy_info_ptr->pheno_mfreq = 0.0;
  gendummy_info_ptr->dosage_freq = 0.0;
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
  uint16_t gt_present;
  uint16_t qual_present;
  uint32_t dosage_field_idx;
  uint32_t hds_field_idx;
  uintptr_t line_idx;  // for error reporting
} GparseReadVcfMetadata;

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

PglErr GparseFlush(const GparseRecord* grp, uint32_t write_block_size, STPgenWriter* spgwp) {
  PglErr reterr = kPglRetSuccess;
  {
    const uintptr_t* allele_idx_offsets = SpgwGetAlleleIdxOffsets(spgwp);
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
          reterr = SpgwAppendMultiallelicSparse(genovec, patch_01_set, patch_01_vals, patch_10_set, patch_10_vals, cur_gwmp->patch_01_ct, cur_gwmp->patch_10_ct, spgwp);
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
                reterr = SpgwAppendMultiallelicSparse(genovec, patch_01_set, patch_01_vals, patch_10_set, patch_10_vals, cur_gwmp->patch_01_ct, cur_gwmp->patch_10_ct, spgwp);
              } else {
                reterr = SpgwAppendMultiallelicGenovecHphase(genovec, patch_01_set, patch_01_vals, patch_10_set, patch_10_vals, phasepresent, phaseinfo, cur_gwmp->patch_01_ct, cur_gwmp->patch_10_ct, spgwp);
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


typedef struct ImportSampleIdContextStruct {
  const char* const_fid;
  uint32_t const_fid_slen;
  uint32_t double_id;
  uint32_t id_delim_sid;
  uint32_t two_part_null_fid;
  uint32_t two_part_null_sid;
  uint32_t omit_fid_0;
  char id_delim;
} ImportSampleIdContext;

void InitImportSampleIdContext(const char* const_fid, MiscFlags misc_flags, ImportFlags import_flags, char id_delim, ImportSampleIdContext* isicp) {
  isicp->const_fid = nullptr;
  isicp->const_fid_slen = 0;
  if (const_fid && (!memequal_k(const_fid, "0", 2))) {
    isicp->const_fid_slen = strlen(const_fid);
    isicp->const_fid = const_fid;
  }
  isicp->double_id = (import_flags / kfImportDoubleId) & 1;
  isicp->id_delim_sid = (misc_flags / kfMiscIidSid) & 1;
  isicp->two_part_null_fid = 0;
  isicp->two_part_null_sid = 0;
  isicp->omit_fid_0 = 0;
  isicp->id_delim = id_delim;
}

PglErr ImportSampleId(const char* input_id_iter, const char* input_id_end, const ImportSampleIdContext* isicp, char** write_iterp) {
  PglErr reterr = kPglRetSuccess;
  {
    const char id_delim = isicp->id_delim;
    char* write_iter = *write_iterp;
    if (id_delim) {
      const char* first_delim = S_CAST(const char*, memchr(input_id_iter, ctou32(id_delim), input_id_end - input_id_iter));
      if (unlikely(!first_delim)) {
        snprintf(g_logbuf, kLogbufSize, "Error: No '%c' in sample ID.\n", id_delim);
        goto ImportSampleId_ret_INCONSISTENT_INPUT_2;
      }
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
      const char* second_part_end = S_CAST(const char*, memchr(second_part_start, ctou32(id_delim), input_id_end - second_part_start));
      const char* third_part_start = nullptr;
      uint32_t third_part_slen = 0;
      if (second_part_end) {
        if (unlikely(second_part_start == second_part_end)) {
          snprintf(g_logbuf, kLogbufSize, "Error: Consecutive instances of '%c' in sample ID.\n", id_delim);
          goto ImportSampleId_ret_INCONSISTENT_INPUT_2;
        }
        third_part_start = &(second_part_end[1]);
        third_part_slen = input_id_end - third_part_start;
        if (unlikely(memchr(third_part_start, ctou32(id_delim), third_part_slen))) {
          snprintf(g_logbuf, kLogbufSize, "Error: More than two instances of '%c' in sample ID.\n", id_delim);
          goto ImportSampleId_ret_INCONSISTENT_INPUT_2;
        }
        if (unlikely(third_part_slen > kMaxIdSlen)) {
          goto ImportSampleId_ret_MALFORMED_INPUT_LONG_ID;
        }
      } else {
        second_part_end = input_id_end;
        if (isicp->id_delim_sid) {
          if (unlikely((first_part_slen == 1) && (*input_id_iter == '0'))) {
            logerrputs("Error: Sample ID induces an invalid IID of '0'.\n");
            goto ImportSampleId_ret_INCONSISTENT_INPUT;
          }
          if (isicp->two_part_null_fid) {
            write_iter = strcpya_k(write_iter, "0\t");
          }
        }
      }
      if (!isicp->omit_fid_0) {
        write_iter = memcpyax(write_iter, input_id_iter, first_part_slen, '\t');
      }
      const uint32_t second_part_slen = second_part_end - second_part_start;
      if (unlikely((third_part_start || (!isicp->id_delim_sid)) && (second_part_slen == 1) && (*second_part_end == '0'))) {
        logerrputs("Error: Sample ID induces an invalid IID of '0'.\n");
        goto ImportSampleId_ret_INCONSISTENT_INPUT;
      }
      if (unlikely(second_part_slen > kMaxIdSlen)) {
        goto ImportSampleId_ret_MALFORMED_INPUT_LONG_ID;
      }
      write_iter = memcpya(write_iter, second_part_start, second_part_slen);
      if (third_part_start) {
        *write_iter++ = '\t';
        write_iter = memcpya(write_iter, third_part_start, third_part_slen);
      } else if (isicp->two_part_null_sid) {
        write_iter = strcpya_k(write_iter, "\t0");
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
  return kPglRetSuccess;
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
      if (isicp->id_delim_sid && (!isicp->two_part_null_fid)) {
        // only possible if there's never a third field
        iid_slen = first_delim - input_id_iter;
      } else {
        const char* second_part_start = &(first_delim[1]);
        const char* second_part_end = S_CAST(const char*, memchr(second_part_start, ctou32(id_delim), input_id_end - second_part_start));
        if ((!isicp->id_delim_sid) || second_part_end) {
          iid_start = second_part_start;
          if (!second_part_end) {
            second_part_end = input_id_end;
          }
          iid_slen = second_part_end - second_part_start;
        } else {
          iid_slen = first_delim - iid_start;
        }
      }
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
  return kPglRetSuccess;
}

PglErr VcfSampleLine(const char* preexisting_psamname, const char* const_fid, MiscFlags misc_flags, ImportFlags import_flags, FamCol fam_cols, char id_delim, char idspace_to, char flag_char, char* sample_line_first_id, char* outname, char* outname_end, uint32_t* sample_ct_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  uintptr_t line_idx = 0;
  PglErr reterr = kPglRetSuccess;
  TextStream psam_txs;
  PreinitTextStream(&psam_txs);
  {
    ImportSampleIdContext isic;
    InitImportSampleIdContext(const_fid, misc_flags, import_flags, id_delim, &isic);
    uint32_t write_fid = 1;
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
      uint32_t nonzero_first_field_observed = 0;
      while (ctou32(sample_line_iter[0]) >= ' ') {
        const char* token_end = strchrnul_n(sample_line_iter, '\t');
        const char* first_delim = S_CAST(const char*, memchr(sample_line_iter, ctou32(id_delim), token_end - sample_line_iter));
        if (unlikely(!first_delim)) {
          snprintf(g_logbuf, kLogbufSize, "Error: No '%c' in sample ID.\n", id_delim);
          goto VcfSampleLine_ret_INCONSISTENT_INPUT_2;
        }
        if (!nonzero_first_field_observed) {
          nonzero_first_field_observed = (first_delim != &(sample_line_iter[1])) || (sample_line_iter[0] != '0');
        }
        sample_line_iter = &(first_delim[1]);
        if (memchr(sample_line_iter, ctou32(id_delim), token_end - sample_line_iter)) {
          write_fid = 1;
          write_sid = 1;
          if (!nonzero_first_field_observed) {
            // bugfix (22 Feb 2018): this was != instead of ==
            while (*token_end == '\t') {
              sample_line_iter = &(token_end[1]);
              token_end = strchrnul_n(sample_line_iter, '\t');
              if ((sample_line_iter[0] != '0') || (sample_line_iter[1] != id_delim)) {
                nonzero_first_field_observed = 1;
                break;
              }
            }
          }
          break;
        }
        if (*token_end != '\t') {
          if (isic.id_delim_sid) {
            write_fid = 0;
            write_sid = 1;
          }
          break;
        }
        sample_line_iter = &(token_end[1]);
      }
      if (!nonzero_first_field_observed) {
        write_fid = 0;
        isic.omit_fid_0 = 1;
      }
    } else if ((!isic.const_fid) && (!isic.double_id)) {
      write_fid = 0;
    }
    const char* sample_line_iter = sample_line_first_id;
    uint32_t sample_ct = 0;
    if (!preexisting_psamname) {
      snprintf(outname_end, kMaxOutfnameExtBlen, ".psam");
      if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
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
        if (write_fid || isic.omit_fid_0) {
          isic.two_part_null_fid = isic.id_delim_sid;
          isic.two_part_null_sid = 1 - isic.id_delim_sid;
        }
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
        if (unlikely(fwrite_ck(textbuf_flush, outfile, &write_iter))) {
          goto VcfSampleLine_ret_WRITE_FAIL;
        }
        if (*token_end != '\t') {
          break;
        }
        sample_line_iter = &(token_end[1]);
      }
      if (unlikely(fclose_flush_null(textbuf_flush, write_iter, &outfile))) {
        goto VcfSampleLine_ret_WRITE_FAIL;
      }
    } else {
      // check consistency of IIDs between VCF and .psam file.
      reterr = SizeAndInitTextStream(preexisting_psamname, bigstack_left(), 1, &psam_txs);
      if (unlikely(reterr)) {
        goto VcfSampleLine_ret_TSTREAM_FAIL;
      }
      uint32_t sample_line_eoln = (ctou32(sample_line_iter[0]) < 32);
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
      uint32_t fid_present;
      if (psam_line_start[0] == '#') {
        // only check for matching IIDs for now.
        fid_present = (psam_line_start[1] == 'F');
        // bugfix (12 Apr 2018): forgot to skip this line
        ++line_idx;
        psam_line_start = TextGet(&psam_txs);
      } else {
        fid_present = fam_cols & kfFamCol1;
      }
      for (; psam_line_start; ++line_idx, psam_line_start = TextGet(&psam_txs)) {
        if (unlikely(psam_line_start[0] == '#')) {
          snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of %s starts with a '#'. (This is only permitted before the first nonheader line, and a #FID/IID header line is present it must denote the end of the header block.)\n", line_idx, preexisting_psamname);
          goto VcfSampleLine_ret_MALFORMED_INPUT_WW;
        }
        const char* psam_iid_start = psam_line_start;
        if (fid_present) {
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
  fclose_cond(outfile);
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
  uint32_t qual_present = 0;
  qual_field_idxs[0] = UINT32_MAX;
  qual_field_idxs[1] = UINT32_MAX;
  if (vcf_min_gq >= 0) {
    qual_field_idxs[0] = GetVcfFormatPosition("GQ", format_start, format_end, 2);
    qual_present = (qual_field_idxs[0] != UINT32_MAX);
  }
  if ((vcf_min_dp >= 0) || (vcf_max_dp != 0x7fffffff)) {
    qual_field_idxs[1] = GetVcfFormatPosition("DP", format_start, format_end, 2);
    if (qual_field_idxs[1] != UINT32_MAX) {
      qual_present = 1;
    }
  }
  return qual_present;
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

  // line-specific
  uint32_t gt_present;
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

ENUM_U31_DEF_START()
  kDosageParseOk,
  kDosageParseMissing,
  kDosageParseForceMissing
ENUM_U31_DEF_END(DosageParseResult);

BoolErr ParseVcfGp(const char* gp_iter, uint32_t is_haploid, double import_dosage_certainty, DosageParseResult* dpr_ptr, double* alt_dosage_ptr) {
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
    if (ParseVcfGp(gtext_iter, is_haploid, import_dosage_certainty, dpr_ptr, &alt_dosage)) {
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
  // returns 1 if missing OR parsing error.  error: no_dosage_here still 0.

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
      if (unlikely((!ScanadvDouble(hds_gtext_iter, &dosage2)) || (dosage2 < 0.0) || (dosage2 > 1.0))) {
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
  kVcfParseInvalidDosage
ENUM_U31_DEF_END(VcfParseErr);

VcfParseErr VcfScanBiallelicHdsLine(const VcfImportContext* vicp, const char* format_end, uint32_t* phase_or_dosage_found_ptr, char** line_iter_ptr) {
  // DS+HDS
  // only need to find phase *or* dosage
  const uint32_t sample_ct = vicp->vibc.sample_ct;
  const VcfHalfCall halfcall_mode = vicp->vibc.halfcall_mode;
  const uint32_t gt_present = vicp->vibc.gt_present;
  STD_ARRAY_KREF(uint32_t, 2) qual_field_skips = vicp->vibc.qual_field_skips;
  STD_ARRAY_KREF(int32_t, 2) qual_line_mins = vicp->vibc.qual_line_mins;
  STD_ARRAY_KREF(int32_t, 2) qual_line_maxs = vicp->vibc.qual_line_maxs;
  const uint32_t qual_field_ct = vicp->vibc.qual_field_ct;

  const uint32_t dosage_erase_halfdist = vicp->dosage_erase_halfdist;
  const double import_dosage_certainty = vicp->import_dosage_certainty;
  const uint32_t dosage_field_idx = vicp->dosage_field_idx;
  const uint32_t hds_field_idx = vicp->hds_field_idx;

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
      // minor bugfix: this needs to be continue, not break
      continue;
    }
    if (gt_present) {
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
      if ((dpr != kDosageParseForceMissing) && cur_gt_phased && VcfIsHetShort(cur_gtext_start, halfcall_mode)) {
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
  const uint32_t gt_present = vicp->vibc.gt_present;
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
        if (gt_present) {
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
          // --import-dosage-certainty, gt_present must be ignored
          goto VcfConvertPhasedBiallelicDosageLine_geno_done;
        }
        if (gt_present) {
          cur_geno = ctow(*linebuf_iter) - 48;
          // this nasty code should be spelled out in as few places as possible
          if (cur_geno <= 1) {
            if (is_haploid) {
              cur_geno *= 2;
            } else {
              const char cc = linebuf_iter[3];
              if (((cc != '/') && (cc != '|')) || (linebuf_iter[4] == '.')) {
                const uintptr_t second_allele_idx = ctow(linebuf_iter[2]) - 48;
                if (second_allele_idx <= 1) {
                  cur_geno += second_allele_idx;
                  if (cur_gt_phased && (cur_geno == 1)) {
                    phasepresent_hw |= shifted_bit;
#ifdef USE_AVX2
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
              } else {
                // code triploids, etc. as missing
                // (probable todo: make the scanning pass take this into
                // account)
                // might want to subject handling of 0/0/. to --vcf-half-call
                // control
                cur_geno = 3;
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
            const char cc = linebuf_iter[3];
            if (((cc != '/') && (cc != '|')) || (linebuf_iter[4] == '.')) {
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
            } else {
              // code triploids, etc. as missing
              // might want to subject handling of 0/0/. to --vcf-half-call
              // control
              cur_geno = 3;
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
            // probable todo: check for triploid
            // kVcfHalfCallHaploid, kVcfHalfCallReference
            cur_geno <<= halfcall_mode;
          } else {
            cur_geno = 3;
          }
        } else {
          cur_geno = 3;
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
          cur_geno = 3;
        } else {
          const uint32_t is_haploid = (linebuf_iter[1] != '/') && (linebuf_iter[1] != '|');
          cur_geno = ctow(*linebuf_iter) - 48;
          if (cur_geno < allele_ct) {
            if (cur_geno <= 1) {
              if (is_haploid) {
                cur_geno *= 2;
              } else {
                const char cc = linebuf_iter[3];
                if (((cc != '/') && (cc != '|')) || (linebuf_iter[4] == '.')) {
                  const uintptr_t second_allele_idx = ctow(linebuf_iter[2]) - 48;
                  if (second_allele_idx <= 1) {
                    cur_geno += second_allele_idx;
                  } else if (second_allele_idx < allele_ct) {
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
                } else {
                  // code triploids, etc. as missing
                  cur_geno = 3;
                }
              }
            } else {
              const uint32_t shifted_bit = 1U << sample_idx_lowbits;
              if (is_haploid) {
                patch_10_hw |= shifted_bit;
                *patch_10_vals_iter++ = cur_geno;
                *patch_10_vals_iter++ = cur_geno;
                cur_geno = 2;
              } else {
                const char cc = linebuf_iter[3];
                if (((cc != '/') && (cc != '|')) || (linebuf_iter[4] == '.')) {
                  uintptr_t second_allele_idx = ctow(linebuf_iter[2]) - 48;
                  if (second_allele_idx < allele_ct) {
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
                  } else if (unlikely(second_allele_idx != (~k0LU) * 2)) {
                    // not '.'
                    return kVcfParseInvalidGt;
                  } else if (halfcall_mode == kVcfHalfCallMissing) {
                    cur_geno = 3;
                  } else if (unlikely(halfcall_mode == kVcfHalfCallError)) {
                    return kVcfParseHalfCallError;
                  } else if (halfcall_mode == kVcfHalfCallReference) {
                    patch_01_hw |= shifted_bit;
                    *patch_01_vals_iter++ = cur_geno;
                    cur_geno = 1;
                  } else {
                    // kVcfHalfCallHaploid
                    patch_10_hw |= shifted_bit;
                    *patch_10_vals_iter++ = cur_geno;
                    *patch_10_vals_iter++ = cur_geno;
                    cur_geno = 2;
                  }
                } else {
                  // code triploids, etc. as missing
                  cur_geno = 3;
                }
              }
            }
          } else if (unlikely(cur_geno != (~k0LU) * 2)) {
            // not '.'
            return kVcfParseInvalidGt;
          } else if (halfcall_mode != kVcfHalfCallMissing) {
            const char second_allele_char = linebuf_iter[2];
            if ((second_allele_char != '.') && (!is_haploid)) {
              cur_geno = ctow(second_allele_char) - 48;
              if (unlikely(cur_geno >= allele_ct)) {
                return kVcfParseInvalidGt;
              }
              if (unlikely(halfcall_mode == kVcfHalfCallError)) {
                return kVcfParseHalfCallError;
              }
              // probable todo: check for triploid
              if (cur_geno <= 1) {
                // works for both kVcfHalfCallReference and kVcfHalfCallHaploid
                cur_geno <<= halfcall_mode;
              } else {
                const uint32_t shifted_bit = 1U << sample_idx_lowbits;
                if (halfcall_mode == kVcfHalfCallReference) {
                  patch_01_hw |= shifted_bit;
                  *patch_01_vals_iter++ = cur_geno;
                  cur_geno = 1;
                } else {
                  // kVcfHalfCallHaploid
                  patch_10_hw |= shifted_bit;
                  *patch_10_vals_iter++ = cur_geno;
                  *patch_10_vals_iter++ = cur_geno;
                  cur_geno = 2;
                }
              }
            } else {
              cur_geno = 3;
            }
          } else {
            cur_geno = 3;
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
            // '/' = ascii 47, '|' = ascii 124
            const uint32_t is_haploid = (next_digit != (~k0LU)) && (next_digit != 76);
            if (cur_geno <= 1) {
              if (is_haploid) {
                cur_geno *= 2;
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
                  if (((next_digit != ~k0LU) && (next_digit != 76)) || second_allele_iter[1] == '.') {
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
                    // code triploids, etc. as missing
                    cur_geno = 3;
                  }
                } else if (unlikely(second_allele_idx != (~k0LU) * 2)) {
                  // not '.'
                  return kVcfParseInvalidGt;
                } else if (halfcall_mode == kVcfHalfCallMissing) {
                  cur_geno = 3;
                } else if (unlikely(halfcall_mode == kVcfHalfCallError)) {
                  return kVcfParseHalfCallError;
                } else {
                  // kVcfHalfCallReference, kVcfHalfCallHaploid
                  cur_geno <<= halfcall_mode;
                }
              }
            } else {
              const uint32_t shifted_bit = 1U << sample_idx_lowbits;
              if (is_haploid) {
                patch_10_hw |= shifted_bit;
                *patch_10_vals_iter++ = cur_geno;
                *patch_10_vals_iter++ = cur_geno;
                cur_geno = 2;
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
                  if (((next_digit != ~k0LU) && (next_digit != 76)) || second_allele_iter[1] == '.') {
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
                  } else {
                    // code triploids, etc. as missing
                    cur_geno = 3;
                  }
                } else if (unlikely(second_allele_idx != (~k0LU) * 2)) {
                  // not '.'
                  return kVcfParseInvalidGt;
                } else if (halfcall_mode == kVcfHalfCallMissing) {
                  cur_geno = 3;
                } else if (unlikely(halfcall_mode == kVcfHalfCallError)) {
                  return kVcfParseHalfCallError;
                } else if (halfcall_mode == kVcfHalfCallReference) {
                  patch_01_hw |= shifted_bit;
                  *patch_01_vals_iter++ = cur_geno;
                  cur_geno = 1;
                } else {
                  // kVcfHalfCallHaploid
                  patch_10_hw |= shifted_bit;
                  *patch_10_vals_iter++ = cur_geno;
                  *patch_10_vals_iter++ = cur_geno;
                  cur_geno = 2;
                }
              }
            }
          } else if (unlikely(cur_geno != (~k0LU) * 2)) {
            // not '.'
            return kVcfParseInvalidGt;
          } else if (halfcall_mode != kVcfHalfCallMissing) {
            const char second_allele_first_char = linebuf_iter[2];
            if ((second_allele_first_char != '.') && ((linebuf_iter[1] == '/') || (linebuf_iter[1] == '|'))) {
              cur_geno = ctow(second_allele_first_char) - 48;
              if (unlikely(cur_geno >= 10)) {
                return kVcfParseInvalidGt;
              }
              for (const char* second_allele_iter = &(linebuf_iter[3]); ; ++second_allele_iter) {
                uintptr_t next_digit = ctow(*second_allele_iter) - 48;
                if (next_digit >= 10) {
                  break;
                }
                cur_geno = cur_geno * 10 + next_digit;
                if (cur_geno >= allele_ct) {
                  return kVcfParseInvalidGt;
                }
              }
              if (unlikely(halfcall_mode == kVcfHalfCallError)) {
                return kVcfParseHalfCallError;
              }
              if (cur_geno <= 1) {
                // works for both kVcfHalfCallReference and kVcfHalfCallHaploid
                cur_geno <<= halfcall_mode;
              } else {
                const uint32_t shifted_bit = 1U << sample_idx_lowbits;
                if (halfcall_mode == kVcfHalfCallReference) {
                  patch_01_hw |= shifted_bit;
                  *patch_01_vals_iter++ = cur_geno;
                  cur_geno = 1;
                } else {
                  // kVcfHalfCallHaploid
                  patch_10_hw |= shifted_bit;
                  *patch_10_vals_iter++ = cur_geno;
                  *patch_10_vals_iter++ = cur_geno;
                  cur_geno = 2;
                }
              }
            } else {
              cur_geno = 3;
            }
          } else {
            cur_geno = 3;
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
        cur_geno = 3;
      } else {
        const uint32_t is_phased = (linebuf_iter[1] == '|');
        const uint32_t is_haploid = (!is_phased) && (linebuf_iter[1] != '/');
        cur_geno = ctow(*linebuf_iter) - 48;
        if (cur_geno <= 1) {
          if (is_haploid) {
            cur_geno *= 2;
          } else {
            const char cc = linebuf_iter[3];
            if (((cc != '/') && (cc != '|')) || (linebuf_iter[4] == '.')) {
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
            } else {
              // code triploids, etc. as missing
              // might want to subject handling of 0/0/. to --vcf-half-call
              // control
              cur_geno = 3;
            }
          }
        } else if (unlikely(cur_geno != (~k0LU) * 2)) {
          // not '.'
          return kVcfParseInvalidGt;
        } else if (halfcall_mode != kVcfHalfCallMissing) {
          const char second_allele_char = linebuf_iter[2];
          if ((second_allele_char != '.') && (!is_haploid)) {
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
          } else {
            cur_geno = 3;
          }
        } else {
          cur_geno = 3;
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
          cur_geno = 3;
        } else {
          const uint32_t is_phased = (linebuf_iter[1] == '|');
          const uint32_t is_haploid = (!is_phased) && (linebuf_iter[1] != '/');
          cur_geno = ctow(*linebuf_iter) - 48;
          if (cur_geno < allele_ct) {
            if (cur_geno <= 1) {
              if (is_haploid) {
                cur_geno *= 2;
              } else {
                const char cc = linebuf_iter[3];
                if (((cc != '/') && (cc != '|')) || (linebuf_iter[4] == '.')) {
                  const uintptr_t second_allele_idx = ctow(linebuf_iter[2]) - 48;
                  if (second_allele_idx < allele_ct) {
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
                  } else if (unlikely(second_allele_idx != (~k0LU) * 2)) {
                    // not '.'
                    return kVcfParseInvalidGt;
                  } else if (halfcall_mode == kVcfHalfCallMissing) {
                    cur_geno = 3;
                  } else if (unlikely(halfcall_mode == kVcfHalfCallError)) {
                    return kVcfParseHalfCallError;
                  } else {
                    // kVcfHalfCallHaploid, kVcfHalfCallReference
                    if (is_phased && (halfcall_mode == kVcfHalfCallReference) && cur_geno) {
                      phasepresent_hw |= shifted_bit;
                      phaseinfo_hw |= shifted_bit;
                    } else {
                      cur_geno <<= halfcall_mode;
                    }
                  }
                } else {
                  // code triploids, etc. as missing
                  cur_geno = 3;
                }
              }
            } else {
              // first allele index > 1
              if (is_haploid) {
                patch_10_hw |= shifted_bit;
                *patch_10_vals_iter++ = cur_geno;
                *patch_10_vals_iter++ = cur_geno;
                cur_geno = 2;
              } else {
                const char cc = linebuf_iter[3];
                if (((cc != '/') && (cc != '|')) || (linebuf_iter[4] == '.')) {
                  uintptr_t second_allele_idx = ctow(linebuf_iter[2]) - 48;
                  if (second_allele_idx < allele_ct) {
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
                  } else if (unlikely(second_allele_idx != (~k0LU) * 2)) {
                    // not '.'
                    return kVcfParseInvalidGt;
                  } else if (halfcall_mode == kVcfHalfCallMissing) {
                    cur_geno = 3;
                  } else if (unlikely(halfcall_mode == kVcfHalfCallError)) {
                    return kVcfParseHalfCallError;
                  } else if (halfcall_mode == kVcfHalfCallReference) {
                    if (is_phased) {
                      phasepresent_hw |= shifted_bit;
                      phaseinfo_hw |= shifted_bit;
                    }
                    patch_01_hw |= shifted_bit;
                    *patch_01_vals_iter++ = cur_geno;
                    cur_geno = 1;
                  } else {
                    // kVcfHalfCallHaploid
                    patch_10_hw |= shifted_bit;
                    *patch_10_vals_iter++ = cur_geno;
                    *patch_10_vals_iter++ = cur_geno;
                    cur_geno = 2;
                  }
                } else {
                  // code triploids, etc. as missing
                  cur_geno = 3;
                }
              }
            }
          } else if (unlikely(cur_geno != (~k0LU) * 2)) {
            // not '.'
            return kVcfParseInvalidGt;
          } else if (halfcall_mode != kVcfHalfCallMissing) {
            const char second_allele_char = linebuf_iter[2];
            if ((second_allele_char != '.') && (!is_haploid)) {
              cur_geno = ctow(second_allele_char) - 48;
              if (unlikely(cur_geno >= allele_ct)) {
                return kVcfParseInvalidGt;
              }
              if (unlikely(halfcall_mode == kVcfHalfCallError)) {
                return kVcfParseHalfCallError;
              }
              // probable todo: check for triploid
              if (is_phased && (halfcall_mode == kVcfHalfCallReference) && cur_geno) {
                phasepresent_hw |= shifted_bit;
              }
              if (cur_geno <= 1) {
                // works for both kVcfHalfCallReference and kVcfHalfCallHaploid
                cur_geno <<= halfcall_mode;
              } else {
                if (halfcall_mode == kVcfHalfCallReference) {
                  patch_01_hw |= shifted_bit;
                  *patch_01_vals_iter++ = cur_geno;
                  cur_geno = 1;
                } else {
                  // kVcfHalfCallHaploid
                  patch_10_hw |= shifted_bit;
                  *patch_10_vals_iter++ = cur_geno;
                  *patch_10_vals_iter++ = cur_geno;
                  cur_geno = 2;
                }
              }
            } else {
              cur_geno = 3;
            }
          } else {
            cur_geno = 3;
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
            // '/' = ascii 47, '|' = ascii 124
            const uint32_t is_phased = (next_digit == 76);
            const uint32_t is_haploid = (!is_phased) && (next_digit != (~k0LU));
            if (cur_geno <= 1) {
              if (is_haploid) {
                cur_geno *= 2;
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
                  if (((next_digit != ~k0LU) && (next_digit != 76)) || second_allele_iter[1] == '.') {
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
                    // code triploids, etc. as missing
                    cur_geno = 3;
                  }
                } else if (unlikely(second_allele_idx != (~k0LU) * 2)) {
                  // not '.'
                  return kVcfParseInvalidGt;
                } else if (halfcall_mode == kVcfHalfCallMissing) {
                  cur_geno = 3;
                } else if (unlikely(halfcall_mode == kVcfHalfCallError)) {
                  return kVcfParseHalfCallError;
                } else {
                  // kVcfHalfCallReference, kVcfHalfCallHaploid
                  if (is_phased && (halfcall_mode == kVcfHalfCallReference) && cur_geno) {
                    phasepresent_hw |= shifted_bit;
                    phaseinfo_hw |= shifted_bit;
                  } else {
                    cur_geno <<= halfcall_mode;
                  }
                }
              }
            } else {
              // first allele index > 1
              if (is_haploid) {
                patch_10_hw |= shifted_bit;
                *patch_10_vals_iter++ = cur_geno;
                *patch_10_vals_iter++ = cur_geno;
                cur_geno = 2;
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
                  if (((next_digit != ~k0LU) && (next_digit != 76)) || second_allele_iter[1] == '.') {
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
                  } else {
                    // code triploids, etc. as missing
                    cur_geno = 3;
                  }
                } else if (unlikely(second_allele_idx != (~k0LU) * 2)) {
                  // not '.'
                  return kVcfParseInvalidGt;
                } else if (halfcall_mode == kVcfHalfCallMissing) {
                  cur_geno = 3;
                } else if (unlikely(halfcall_mode == kVcfHalfCallError)) {
                  return kVcfParseHalfCallError;
                } else if (halfcall_mode == kVcfHalfCallReference) {
                  if (is_phased) {
                    phasepresent_hw |= shifted_bit;
                    phaseinfo_hw |= shifted_bit;
                  }
                  patch_01_hw |= shifted_bit;
                  *patch_01_vals_iter++ = cur_geno;
                  cur_geno = 1;
                } else {
                  // kVcfHalfCallHaploid
                  patch_10_hw |= shifted_bit;
                  *patch_10_vals_iter++ = cur_geno;
                  *patch_10_vals_iter++ = cur_geno;
                  cur_geno = 2;
                }
              }
            }
          } else if (unlikely(cur_geno != (~k0LU) * 2)) {
            // not '.'
            return kVcfParseInvalidGt;
          } else if (halfcall_mode != kVcfHalfCallMissing) {
            const char second_allele_first_char = linebuf_iter[2];
            if ((second_allele_first_char != '.') && ((linebuf_iter[1] == '/') || (linebuf_iter[1] == '|'))) {
              cur_geno = ctow(second_allele_first_char) - 48;
              if (unlikely(cur_geno >= 10)) {
                return kVcfParseInvalidGt;
              }
              for (const char* second_allele_iter = &(linebuf_iter[3]); ; ++second_allele_iter) {
                uintptr_t next_digit = ctow(*second_allele_iter) - 48;
                if (next_digit >= 10) {
                  break;
                }
                cur_geno = cur_geno * 10 + next_digit;
                if (cur_geno >= allele_ct) {
                  return kVcfParseInvalidGt;
                }
              }
              if (unlikely(halfcall_mode == kVcfHalfCallError)) {
                return kVcfParseHalfCallError;
              }
              // probable todo: check for triploid
              if ((linebuf_iter[1] == '|') && (halfcall_mode == kVcfHalfCallReference) && cur_geno) {
                phasepresent_hw |= shifted_bit;
              }
              if (cur_geno <= 1) {
                // works for both kVcfHalfCallReference and kVcfHalfCallHaploid
                cur_geno <<= halfcall_mode;
              } else {
                if (halfcall_mode == kVcfHalfCallReference) {
                  patch_01_hw |= shifted_bit;
                  *patch_01_vals_iter++ = cur_geno;
                  cur_geno = 1;
                } else {
                  // kVcfHalfCallHaploid
                  patch_10_hw |= shifted_bit;
                  *patch_10_vals_iter++ = cur_geno;
                  *patch_10_vals_iter++ = cur_geno;
                  cur_geno = 2;
                }
              }
            } else {
              cur_geno = 3;
            }
          } else {
            cur_geno = 3;
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
        uint32_t qual_field_ct = grp->metadata.read_vcf.qual_present;
        if (qual_field_ct) {
          qual_field_ct = VcfQualScanInit2(grp->metadata.read_vcf.qual_field_idxs, qual_mins, qual_maxs, vic.vibc.qual_field_skips, vic.vibc.qual_line_mins, vic.vibc.qual_line_maxs);
        }
        vic.vibc.gt_present = grp->metadata.read_vcf.gt_present;
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

static_assert(!kVcfHalfCallReference, "VcfToPgen() assumes kVcfHalfCallReference == 0.");
static_assert(kVcfHalfCallHaploid == 1, "VcfToPgen() assumes kVcfHalfCallHaploid == 1.");
PglErr VcfToPgen(const char* vcfname, const char* preexisting_psamname, const char* const_fid, const char* dosage_import_field, MiscFlags misc_flags, ImportFlags import_flags, uint32_t no_samples_ok, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, double import_dosage_certainty, char id_delim, char idspace_to, int32_t vcf_min_gq, int32_t vcf_min_dp, int32_t vcf_max_dp, VcfHalfCall halfcall_mode, FamCol fam_cols, uint32_t max_thread_ct, char* outname, char* outname_end, ChrInfo* cip, uint32_t* pgen_generated_ptr, uint32_t* psam_generated_ptr) {
  // Now performs a 2-pass load.  Yes, this can be slower than plink 1.9, but
  // it's necessary to use the Pgen_writer classes for now (since we need to
  // know upfront how many variants there are, and whether phase/dosage is
  // present).
  // preexisting_psamname should be nullptr if no such file was specified.
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  FILE* pvarfile = nullptr;
  uintptr_t line_idx = 1;
  const uint32_t half_call_explicit_error = (halfcall_mode == kVcfHalfCallError);
  PglErr reterr = kPglRetSuccess;
  VcfParseErr vcf_parse_err = kVcfParseOk;
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

    reterr = InitTextStreamEx(vcfname, 1, kMaxLongLine, max_line_blen, ClipU32(max_thread_ct - 1, 1, 4), &vcf_txs);
    if (unlikely(reterr)) {
      if (reterr == kPglRetOpenFail) {
        const uint32_t slen = strlen(vcfname);
        if ((!StrEndsWith(vcfname, ".vcf", slen)) &&
            (!StrEndsWith(vcfname, ".vcf.gz", slen))) {
          logerrprintfww("Error: Failed to open %s : %s. (--vcf expects a complete filename; did you forget '.vcf' at the end?)\n", vcfname, strerror(errno));
          goto VcfToPgen_ret_1;
        }
      }
      goto VcfToPgen_ret_TSTREAM_FAIL;
    }
    const uint32_t allow_extra_chrs = (misc_flags / kfMiscAllowExtraChrs) & 1;
    uint32_t dosage_import_field_slen = 0;

    // == 2 when searched for and found in header
    // then becomes a boolean
    uint32_t format_hds_search = 0;

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
    vic.vibc.halfcall_mode = halfcall_mode;

    vic.dosage_is_gp = ctx.dosage_is_gp;

    // always positive
    vic.dosage_erase_halfdist = kDosage4th - dosage_erase_thresh;

    vic.import_dosage_certainty = import_dosage_certainty;

    vic.dosage_field_idx = UINT32_MAX;
    vic.hds_field_idx = UINT32_MAX;

    max_line_blen = 1;  // Now means "max_observed_line_blen".
    uint32_t format_gt_present = 0;
    uint32_t format_gq_relevant = 0;
    uint32_t format_dp_relevant = 0;
    uint32_t format_dosage_relevant = 0;
    uint32_t info_pr_present = 0;
    uint32_t info_pr_nonflag_present = 0;
    uint32_t info_nonpr_present = 0;
    uint32_t chrset_present = 0;
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
        if ((line_idx == 1) && (memequal_k(line_iter, "BCF", 3))) {
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
      // ##INFO: note presence of INFO:PR, note presence of at least one non-PR
      //         field, keep data (though, if INFO:PR is the *only* field,
      //         omit it from the .pvar for consistency with --make-pgen
      //         default)
      //         update (8 Sep 2017): nonflag INFO:PR is noted, and not treated
      //         specially unless provisional-reference INFO:PR output would
      //         conflict with it
      // ##FORMAT: note presence of FORMAT:GT and FORMAT:GP, discard
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
      // filtered out), we wait until second pass to write the .pvar.
      if (StrStartsWithUnsafe(&(line_iter[2]), "chrSet=<")) {
        if (unlikely(chrset_present)) {
          logerrputs("Error: Multiple ##chrSet header lines in --vcf file.\n");
          goto VcfToPgen_ret_MALFORMED_INPUT;
        }
        chrset_present = 1;
        // .pvar loader will print a warning if necessary
        reterr = ReadChrsetHeaderLine(&(line_iter[10]), "--vcf file", misc_flags, line_idx, cip);
        if (unlikely(reterr)) {
          goto VcfToPgen_ret_1;
        }
      } else if (StrStartsWithUnsafe(&(line_iter[2]), "FORMAT=<ID=GT,Number=")) {
        if (unlikely(format_gt_present)) {
          logerrputs("Error: Duplicate FORMAT:GT header line in --vcf file.\n");
          goto VcfToPgen_ret_MALFORMED_INPUT;
        }
        if (unlikely(!StrStartsWithUnsafe(&(line_iter[2 + strlen("FORMAT=<ID=GT,Number=")]), "1,Type=String,Description="))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Header line %" PRIuPTR " of --vcf file does not have expected FORMAT:GT format.\n", line_idx);
          goto VcfToPgen_ret_MALFORMED_INPUT_WW;
        }
        format_gt_present = 1;
      } else if ((vcf_min_gq != -1) && StrStartsWithUnsafe(&(line_iter[2]), "FORMAT=<ID=GQ,Number=1,Type=")) {
        if (unlikely(format_gq_relevant)) {
          logerrputs("Error: Duplicate FORMAT:GQ header line in --vcf file.\n");
          goto VcfToPgen_ret_MALFORMED_INPUT;
        }
        format_gq_relevant = 1;
      } else if (((vcf_min_dp != -1) || (vcf_max_dp != 0x7fffffff)) && StrStartsWithUnsafe(&(line_iter[2]), "FORMAT=<ID=DP,Number=1,Type=")) {
        if (unlikely(format_dp_relevant)) {
          logerrputs("Error: Duplicate FORMAT:DP header line in --vcf file.\n");
          goto VcfToPgen_ret_MALFORMED_INPUT;
        }
        format_dp_relevant = 1;
      } else if (dosage_import_field && StrStartsWithUnsafe(&(line_iter[2]), "FORMAT=<ID=")) {
        // strlen("##FORMAT=<ID=") == 13
        if (memequal(&(line_iter[13]), dosage_import_field, dosage_import_field_slen) && (line_iter[13 + dosage_import_field_slen] == ',')) {
          if (unlikely(format_dosage_relevant)) {
            logerrprintfww("Error: Duplicate FORMAT:%s header line in --vcf file.\n", dosage_import_field);
            goto VcfToPgen_ret_MALFORMED_INPUT_WW;
          }
          format_dosage_relevant = 1;
        } else if (format_hds_search && memequal_k(&(line_iter[13]), "HDS,", 4)) {
          if (unlikely(format_hds_search == 2)) {
            logerrputs("Error: Duplicate FORMAT:HDS header line in --vcf file.\n");
            goto VcfToPgen_ret_MALFORMED_INPUT;
          }
          format_hds_search = 2;
        } else if (unforced_gp && memequal_k(&(line_iter[13]), "DS,", 3)) {
          logerrputs("Error: --vcf dosage=GP specified, but --import-dosage-certainty was not and\nFORMAT:DS header line is present.\nSince " PROG_NAME_STR " collapses genotype probabilities down to dosages (even when\nperforming simple operations like \"" PROG_NAME_STR " --vcf ... --export bgen-1.2 ...\"),\n'dosage=GP' almost never makes sense in this situation.  Either change it to\n'dosage=DS' (if dosages are good enough for your analysis), or use another\nprogram to work with the genotype probabilities.\nThere is one notable exception to the preceding recommendation: you are writing\na script to process VCF files that are guaranteed to have FORMAT:GP, but may or\nmay not have FORMAT:DS.  You can use 'dosage=GP-force' to suppress this error\nin that situation.\n");
          goto VcfToPgen_ret_INCONSISTENT_INPUT;
        }
      } else if (StrStartsWithUnsafe(&(line_iter[2]), "INFO=<ID=")) {
        if (StrStartsWithUnsafe(&(line_iter[2 + strlen("INFO=<ID=")]), "PR,Number=")) {
          if (unlikely(info_pr_present || info_pr_nonflag_present)) {
            logerrputs("Error: Duplicate INFO:PR header line in --vcf file.\n");
            goto VcfToPgen_ret_MALFORMED_INPUT;
          }
          info_pr_nonflag_present = !StrStartsWithUnsafe(&(line_iter[2 + strlen("INFO=<ID=PR,Number=")]), "0,Type=Flag,Description=");
          info_pr_present = 1 - info_pr_nonflag_present;
          if (info_pr_nonflag_present) {
            logerrprintfww("Warning: Header line %" PRIuPTR " of --vcf file has an unexpected definition of INFO:PR. This interferes with a few merge and liftover operations.\n", line_idx);
          }
        } else {
          info_nonpr_present = 1;
        }
      }
      line_iter = AdvPastDelim(line_iter, '\n');
      const uint32_t prev_line_blen = line_iter - prev_line_start;
      if (prev_line_blen > max_line_blen) {
        max_line_blen = prev_line_blen;
      }
    }
    const uint32_t require_gt = (import_flags / kfImportVcfRequireGt) & 1;
    if (unlikely((!format_gt_present) && require_gt)) {
      logerrputs("Error: No FORMAT:GT field in --vcf file header, when --vcf-require-gt was\nspecified.\n");
      goto VcfToPgen_ret_INCONSISTENT_INPUT;
    }
    if ((!format_gq_relevant) && (vcf_min_gq != -1)) {
      logerrputs("Warning: No FORMAT:GQ field in --vcf file header.  --vcf-min-gq ignored.\n");
      vcf_min_gq = -1;
    }
    if ((!format_dp_relevant) && ((vcf_max_dp != 0x7fffffff) || (vcf_min_dp != -1))) {
      logerrputs("Warning: No FORMAT:DP field in --vcf file header.  --vcf-{max,min}-dp ignored.\n");
      vcf_max_dp = 0x7fffffff;
      vcf_min_dp = -1;
    }
    const uint32_t format_gq_or_dp_relevant = format_gq_relevant || format_dp_relevant;
    if (format_hds_search) {
      --format_hds_search;
      if (!format_hds_search) {
        if (!format_dosage_relevant) {
          logerrputs("Warning: No FORMAT:DS or :HDS field in --vcf file header. Dosages will not be\nimported.\n");
        } else {
          logerrputs("Warning: No FORMAT:HDS field in --vcf file header.  Dosages will be imported\n(from FORMAT:DS), but phase information will be limited or absent.\n");
        }
      }
    } else if ((!format_dosage_relevant) && dosage_import_field) {
      logerrprintfww("Warning: No FORMAT:%s field in --vcf file header. Dosages will not be imported.\n", dosage_import_field);
    }
    FinalizeChrset(misc_flags, cip);
    // don't call FinalizeChrInfo here, since this may be followed by
    // --pmerge, etc.

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
      goto VcfToPgen_ret_INCONSISTENT_INPUT;
    }
    vic.vibc.sample_ct = sample_ct;
    // bugfix (5 Jun 2018): must initialize qual_field_ct to zero
    vic.vibc.qual_field_ct = 0;

    uint32_t variant_ct = 0;
    uintptr_t max_variant_ct = bigstack_left() / sizeof(intptr_t);
    if (info_pr_present) {
      // nonref_flags
      max_variant_ct -= BitCtToAlignedWordCt(max_variant_ct) * kWordsPerVec;
    }
#ifdef __LP64__
    if (max_variant_ct > 0x7ffffffd) {
      max_variant_ct = 0x7ffffffd;
    }
#endif
    uintptr_t base_chr_present[kChrExcludeWords];
    ZeroWArr(kChrExcludeWords, base_chr_present);

    const uintptr_t header_line_ct = line_idx;
    const uint32_t max_variant_ctaw = BitCtToAlignedWordCt(max_variant_ct);
    // don't need dosage_flags or dphase_flags; dosage overrides GT so slow
    // parse needed
    uintptr_t* nonref_flags = nullptr;
    if (info_pr_present) {
      nonref_flags = S_CAST(uintptr_t*, bigstack_alloc_raw_rd(max_variant_ctaw * sizeof(intptr_t)));
    }
    uintptr_t* nonref_flags_iter = nonref_flags;
    uintptr_t* allele_idx_offsets = R_CAST(uintptr_t*, g_bigstack_base);
    if (halfcall_mode == kVcfHalfCallDefault) {
      halfcall_mode = kVcfHalfCallError;
    }
    uintptr_t max_postformat_blen = 1;  // starting from tab at end of FORMAT
    uintptr_t variant_skip_ct = 0;
    uintptr_t nonref_word = 0;
    uintptr_t allele_idx_end = 0;
    uint32_t max_alt_ct = 1;
    uint32_t max_allele_slen = 1;
    uint32_t max_qualfilterinfo_slen = 6;
    uint32_t phase_or_dosage_found = 0;

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
        goto VcfToPgen_ret_MALFORMED_INPUT_WW;
      }

      // note REF length
      char* ref_allele_start = &(id_end[1]);
      linebuf_iter = FirstPrespace(ref_allele_start);
      if (unlikely(*linebuf_iter != '\t')) {
        goto VcfToPgen_ret_MISSING_TOKENS;
      }
      uint32_t cur_max_allele_slen = linebuf_iter - ref_allele_start;

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
        vic.vibc.gt_present = memequal_k(linebuf_iter, "GT", 2) && ((linebuf_iter[2] == ':') || (linebuf_iter[2] == '\t'));
        if (require_gt && (!vic.vibc.gt_present)) {
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
      reterr = GetOrAddChrCodeDestructive("--vcf file", line_idx, allow_extra_chrs, line_iter, chr_code_end, cip, &cur_chr_code);
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
      if (info_pr_present) {
        if (PrInInfoToken(info_end - info_start, info_start)) {
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
        if (phase_or_dosage_found) {
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
        // Don't bother multithreading this since it's trivial.
        if ((vic.hds_field_idx != UINT32_MAX) || (vic.dosage_field_idx != UINT32_MAX)) {
          if (alt_ct == 1) {
            vcf_parse_err = VcfScanBiallelicHdsLine(&vic, format_end, &phase_or_dosage_found, &line_iter);
          } else {
            logerrputs("Error: --vcf multiallelic dosage import is under development.\n");
            reterr = kPglRetNotYetSupported;
            goto VcfToPgen_ret_1;
          }
          if (unlikely(vcf_parse_err)) {
            goto VcfToPgen_ret_PARSE;
          }
        } else {
          if (!vic.vibc.gt_present) {
            goto VcfToPgen_linescan_done;
          }
          if (alt_ct < 10) {
            phase_or_dosage_found = VcfScanShortallelicLine(&(vic.vibc), format_end, &line_iter);
          } else {
            phase_or_dosage_found = VcfScanLongallelicLine(&(vic.vibc), format_end, &line_iter);
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
        if (variant_ct == 0x7ffffffd) {
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
      logerrputs("Error: No variants in --vcf file.\n");
      goto VcfToPgen_ret_INCONSISTENT_INPUT;
    }

    putc_unlocked('\r', stdout);
    if (!variant_skip_ct) {
      logprintf("--vcf: %u variant%s scanned.\n", variant_ct, (variant_ct == 1)? "" : "s");
    } else {
      logprintf("--vcf: %u variant%s scanned (%" PRIuPTR " skipped).\n", variant_ct, (variant_ct == 1)? "" : "s", variant_skip_ct);
    }

    if (allele_idx_end > 2 * variant_ct) {
      allele_idx_offsets[variant_ct] = allele_idx_end;
      BigstackFinalizeW(allele_idx_offsets, variant_ct + 1);
    } else {
      allele_idx_offsets = nullptr;
    }

    uint32_t calc_thread_ct;
    {
      // Close file, then reopen with a smaller line-load buffer and (if bgzf)
      // reduce decompression thread count.  2 is good in the simplest cases
      // (no GQ/DP filter, no dosage), otherwise limit to 1.
      uint32_t decompress_thread_ct = 1;
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
      if (CleanupTextStream2(vcfname, &vcf_txs, &reterr)) {
        goto VcfToPgen_ret_1;
      }
      BigstackEndReset(bigstack_end_mark);
      reterr = InitTextStreamEx(vcfname, 1, kMaxLongLine, MAXV(max_line_blen, kTextStreamBlenFast), decompress_thread_ct, &vcf_txs);
      if (unlikely(reterr)) {
        goto VcfToPgen_ret_TSTREAM_REWIND_FAIL;
      }
      if (calc_thread_ct + decompress_thread_ct > max_thread_ct) {
        calc_thread_ct = MAXV(1, max_thread_ct - decompress_thread_ct);
      }
    }

    snprintf(outname_end, kMaxOutfnameExtBlen, ".pvar");
    if (unlikely(fopen_checked(outname, FOPEN_WB, &pvarfile))) {
      goto VcfToPgen_ret_OPEN_FAIL;
    }
    for (line_idx = 1, line_iter = TextLineEnd(&vcf_txs); ; ++line_idx, line_iter = AdvPastDelim(line_iter, '\n')) {
      reterr = TextNextLineUnsafe(&vcf_txs, &line_iter);
      if (unlikely(reterr)) {
        goto VcfToPgen_ret_TSTREAM_REWIND_FAIL;
      }
      if (line_idx == header_line_ct) {
        break;
      }
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
      char* line_end = AdvToDelim(line_iter, '\n');
#ifdef _WIN32
      if (line_end[-1] == '\r') {
        --line_end;
      }
      // NOT safe to use AppendBinaryEoln here.
      if (unlikely(fwrite_checked(line_iter, line_end - line_iter, pvarfile))) {
        goto VcfToPgen_ret_WRITE_FAIL;
      }
      if (unlikely(fputs_checked("\r\n", pvarfile))) {
        goto VcfToPgen_ret_WRITE_FAIL;
      }
#else
      char* line_write_end;
      if (line_end[-1] == '\r') {
        line_write_end = line_end;
        line_end[-1] = '\n';
      } else {
        line_write_end = &(line_end[1]);
      }
      if (unlikely(fwrite_checked(line_iter, line_write_end - line_iter, pvarfile))) {
        goto VcfToPgen_ret_WRITE_FAIL;
      }
#endif
      line_iter = line_end;
    }
    char* write_iter = g_textbuf;
    if (cip->chrset_source) {
      AppendChrsetLine(cip, &write_iter);
    }
    write_iter = strcpya_k(write_iter, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER");
    if (info_nonpr_present) {
      write_iter = strcpya_k(write_iter, "\tINFO");
    }
    AppendBinaryEoln(&write_iter);
    if (unlikely(fwrite_checked(g_textbuf, write_iter - g_textbuf, pvarfile))) {
      goto VcfToPgen_ret_WRITE_FAIL;
    }

    const uint32_t variant_ctl = BitCtToWordCt(variant_ct);
    PgenGlobalFlags phase_dosage_gflags = kfPgenGlobal0;
    GparseFlags gparse_flags = kfGparse0;  // yeah, this is a bit redundant
    // bugfix (22 Jun 2018): if dosage= was specified, we need to reserve
    // phasepresent/dosage_present buffers for each variant even when we know
    // they'll be irrelevant, since they're needed during conversion.
    if (phase_or_dosage_found || format_hds_search || format_dosage_relevant) {
      if ((!format_hds_search) && (!format_dosage_relevant)) {
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
    char* writebuf;
    if (unlikely(bigstack_alloc_c(2 * max_allele_slen + max_qualfilterinfo_slen + kMaxMediumLine + kMaxIdSlen + 32, &writebuf))) {
      goto VcfToPgen_ret_NOMEM;
    }
    write_iter = writebuf;
    char* writebuf_flush = &(writebuf[kMaxMediumLine]);
    if (max_alt_ct > kPglMaxAltAlleleCt) {
      logerrprintfww("Error: VCF file has a variant with %u ALT alleles; this build of " PROG_NAME_STR " is limited to %u.\n", max_alt_ct, kPglMaxAltAlleleCt);
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
      reterr = SpgwInitPhase1(outname, allele_idx_offsets, nonref_flags, variant_ct, sample_ct, phase_dosage_gflags, nonref_flags_storage, &spgw, &spgw_alloc_cacheline_ct, &max_vrec_len);
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

      if (unlikely(
              bigstack_alloc_ucp(calc_thread_ct, &ctx.thread_wkspaces) ||
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
            ctx.block_allele_idx_offsets[parity] = &(SpgwGetAlleleIdxOffsets(&spgw)[vidx_start]);
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
          grp->metadata.read_vcf.gt_present = vic.vibc.gt_present;
          grp->metadata.read_vcf.qual_present = vic.vibc.qual_field_ct;
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
            goto VcfToPgen_ret_TSTREAM_REWIND_FAIL;
          }

          // 1. check if we skip this variant.  chromosome filter and
          //    require_gt can cause this.
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
          } else {
            if (chr_code_base >= kMaxContigs) {
              chr_code_base = cip->xymt_codes[chr_code_base - kMaxContigs];
            }
            if (IsI32Neg(chr_code_base) || (!IsSet(base_chr_present, chr_code_base))) {
              assert(variant_skip_ct);
              line_iter = chr_code_end;
              goto VcfToPgen_load_start;
            }
          }
          // chr_code_base is now a proper numeric chromosome index for
          // non-contigs, and UINT32_MAX if it's a contig name
          char* pos_str = &(chr_code_end[1]);
          char* pos_str_end = AdvToDelim(pos_str, '\t');
          // copy ID, REF verbatim
          linebuf_iter = AdvToNthDelim(&(pos_str_end[1]), 2, '\t');

          // ALT, QUAL, FILTER, INFO
          char* filter_end = AdvToNthDelim(&(linebuf_iter[1]), 3, '\t');
          char* format_start = nullptr;
          char* info_end;
          if (sample_ct) {
            info_end = AdvToDelim(&(filter_end[1]), '\t');
            format_start = &(info_end[1]);
            vic.vibc.gt_present = memequal_k(format_start, "GT", 2) && ((format_start[2] == ':') || (format_start[2] == '\t'));
            if (require_gt && (!vic.vibc.gt_present)) {
              line_iter = format_start;
              goto VcfToPgen_load_start;
            }
          } else {
            info_end = NextPrespace(filter_end);
          }

          // make sure POS starts with an integer, apply --output-chr setting
          uint32_t cur_bp;
          if (unlikely(ScanUintDefcap(pos_str, &cur_bp))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid POS on line %" PRIuPTR " of --vcf file.\n", line_idx);
            goto VcfToPgen_ret_MALFORMED_INPUT_2N;
          }

          if (chr_code_base == UINT32_MAX) {
            write_iter = memcpya(write_iter, line_iter, chr_code_end - line_iter);
          } else {
            write_iter = chrtoa(cip, chr_code_base, write_iter);
          }
          *write_iter++ = '\t';
          write_iter = u32toa(cur_bp, write_iter);

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
            write_iter = memcpya(write_iter, copy_start, linebuf_iter - copy_start);
            if (fwrite_ck(writebuf_flush, pvarfile, &write_iter)) {
              goto VcfToPgen_ret_WRITE_FAIL;
            }
            if (ucc != ',') {
              break;
            }
            copy_start = linebuf_iter;
          }

          if (info_nonpr_present) {
            // VCF specification permits whitespace in INFO field, while PVAR
            // does not.  Check for whitespace and error out if necessary.
            if (unlikely(memchr(filter_end, ' ', info_end - filter_end))) {
              snprintf(g_logbuf, kLogbufSize, "Error: INFO field on line %" PRIuPTR " of --vcf file contains a space; this cannot be imported by " PROG_NAME_STR ".  Remove or reformat the field before reattempting import.\n", line_idx);
              goto VcfToPgen_ret_MALFORMED_INPUT_2N;
            }
            write_iter = memcpya(write_iter, linebuf_iter, info_end - linebuf_iter);
          } else {
            write_iter = memcpya(write_iter, linebuf_iter, filter_end - linebuf_iter);
          }
          AppendBinaryEoln(&write_iter);
          if (!sample_ct) {
            if (++block_vidx == cur_thread_block_vidx_limit) {
              break;
            }
            line_iter = info_end;
            goto VcfToPgen_load_start;
          }
          if ((!vic.vibc.gt_present) && (!format_dosage_relevant) && (!format_hds_search)) {
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
          reterr = GparseFlush(ctx.gparse[parity], prev_block_write_ct, &spgw);
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
    if (unlikely(fclose_flush_null(writebuf_flush, write_iter, &pvarfile))) {
      goto VcfToPgen_ret_WRITE_FAIL;
    }
    if (sample_ct) {
      SpgwFinish(&spgw);
    }
    putc_unlocked('\r', stdout);
    write_iter = strcpya_k(g_logbuf, "--vcf: ");
    const uint32_t outname_base_slen = outname_end - outname;
    if (sample_ct) {
      write_iter = memcpya(write_iter, outname, outname_base_slen + 5);
      write_iter = strcpya_k(write_iter, " + ");
    } else {
      *pgen_generated_ptr = 0;
    }
    write_iter = memcpya(write_iter, outname, outname_base_slen);
    write_iter = strcpya_k(write_iter, ".pvar");
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
  }
  while (0) {
  VcfToPgen_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  VcfToPgen_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  VcfToPgen_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  VcfToPgen_ret_TSTREAM_FAIL:
    putc_unlocked('\n', stdout);
    TextStreamErrPrint("--vcf file", &vcf_txs);
    break;
  VcfToPgen_ret_TSTREAM_REWIND_FAIL:
    putc_unlocked('\n', stdout);
    TextStreamErrPrintRewind("--vcf file", &vcf_txs, &reterr);
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
      putc_unlocked('\n', stdout);
      logerrprintfww("Error: Line %" PRIuPTR " of --vcf file has an invalid %s field.\n", line_idx, dosage_import_field);
      reterr = kPglRetInconsistentInput;
      break;
    }
  VcfToPgen_ret_MISSING_TOKENS:
    putc_unlocked('\n', stdout);
    logerrprintf("Error: Line %" PRIuPTR " of --vcf file has fewer tokens than expected.\n", line_idx);
    reterr = kPglRetMalformedInput;
    break;
  VcfToPgen_ret_MALFORMED_INPUT_2N:
    logputs("\n");
    logerrputsb();
  VcfToPgen_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  VcfToPgen_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetMalformedInput;
    break;
  VcfToPgen_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  VcfToPgen_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 VcfToPgen_ret_1:
  CleanupSpgw(&spgw, &reterr);
  CleanupThreads(&tg);
  CleanupTextStream2("--vcf file", &vcf_txs, &reterr);
  fclose_cond(pvarfile);
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}


PglErr OxSampleToPsam(const char* samplename, const char* ox_missing_code, ImportFlags import_flags, char* outname, char* outname_end, uint32_t* sample_ct_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* psamfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  uintptr_t line_idx = 0;
  TextStream sample_txs;
  PreinitTextStream(&sample_txs);
  {
    uint32_t omp_slen = 2;
    char output_missing_pheno[kMaxMissingPhenostrBlen];
    if (import_flags & kfImportKeepAutoconv) {
      // must use --output-missing-phenotype parameter, which we've validated
      // to be consistent with --input-missing-phenotype
      omp_slen = strlen(g_output_missing_pheno);
      memcpy(output_missing_pheno, g_output_missing_pheno, omp_slen);
    } else {
      // use "NA" since that's always safe
      memcpy_k(output_missing_pheno, "NA", 2);
    }
    const char* missing_catname = g_missing_catname;
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
        const char* token_end = strchrnul(missing_code_iter, ',');
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
          const char* token_end = strchrnul(missing_code_iter, ',');
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

    // New first pass: check whether, from the third line on, all first tokens
    // are '0'.  If yes, we omit FID from the output.
    uint32_t write_fid = 0;
    uint32_t header_lines_left = 2;
    while (1) {
      ++line_idx;
      char* line_start = TextGet(&sample_txs);
      if (!line_start) {
        // bugfix (16 Feb 2018): don't check this if we break out of the loop
        // on non-0 FID
        if (unlikely(TextStreamErrcode2(&sample_txs, &reterr))) {
          goto OxSampleToPsam_ret_TSTREAM_FAIL;
        }
        if (unlikely(header_lines_left)) {
          logerrputs("Error: Empty .sample file.\n");
          goto OxSampleToPsam_ret_MALFORMED_INPUT;
        }
        break;
      }
      if (header_lines_left) {
        --header_lines_left;
        continue;
      }
      if ((line_start[0] != '0') || (!IsSpaceOrEoln(line_start[1]))) {
        write_fid = 1;
        break;
      }
    }
    reterr = TextRewind(&sample_txs);
    if (unlikely(reterr)) {
      goto OxSampleToPsam_ret_TSTREAM_FAIL;
    }
    line_idx = 1;
    char* line_start = TextGet(&sample_txs);
    if (unlikely(!line_start)) {
      reterr = TextStreamRawErrcode(&sample_txs);
      goto OxSampleToPsam_ret_TSTREAM_REWIND_FAIL;
    }
    char* token_end = CurTokenEnd(line_start);
    // todo: allow arbitrary first column name, and ID_2/missing to be absent,
    // to match SNPTEST 2.5.4's de facto change to the .sample file format.
    // (Probably wait till 2.5.4 is a stable version, though.)
    if (unlikely(!strequal_k(line_start, "ID_1", token_end - line_start))) {
      goto OxSampleToPsam_ret_INVALID_SAMPLE_HEADER_1;
    }
    // currently accepts tab as delimiter, though .sample spec technically
    // prohibits that
    char* linebuf_iter = FirstNonTspace(token_end);
    uint32_t token_slen = strlen_se(linebuf_iter);
    if (unlikely(!strequal_k(linebuf_iter, "ID_2", token_slen))) {
      goto OxSampleToPsam_ret_INVALID_SAMPLE_HEADER_1;
    }
    linebuf_iter = FirstNonTspace(&(linebuf_iter[token_slen]));
    token_slen = strlen_se(linebuf_iter);
    if (unlikely(!MatchUpperKLen(linebuf_iter, "MISSING", token_slen))) {
      goto OxSampleToPsam_ret_INVALID_SAMPLE_HEADER_1;
    }
    linebuf_iter = FirstNonTspace(&(linebuf_iter[token_slen]));

    snprintf(outname_end, kMaxOutfnameExtBlen, ".psam");
    if (unlikely(fopen_checked(outname, FOPEN_WB, &psamfile))) {
      goto OxSampleToPsam_ret_OPEN_FAIL;
    }
    // categorical phenotypes are lengthened by 1 character ('C' added in
    // front), so this needs to be 50% larger than maximum possible line length
    // to handle worst case
    uintptr_t linebuf_size = max_line_blen + kDecompressChunkSize;
    linebuf_size += linebuf_size / 2;
    // bugfix (19 Mar 2018): this needs to be min, not max...
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
    write_iter = strcpya_k(write_iter, "IID\tSEX");

    // 0 = not present, otherwise zero-based index (this is fine since first
    //     column has to be FID)
    uint32_t sex_col = 0;

    uint32_t col_ct = 3;

    while (!IsEolnKns(*linebuf_iter)) {
      token_end = CurTokenEnd(linebuf_iter);
      token_slen = token_end - linebuf_iter;
      if (MatchUpperKLen(linebuf_iter, "SEX", token_slen)) {
        if (unlikely(sex_col)) {
          logerrputs("Error: Multiple sex columns in .sample file.\n");
          goto OxSampleToPsam_ret_MALFORMED_INPUT;
        }
        sex_col = col_ct;
      }
      ++col_ct;
      linebuf_iter = FirstNonTspace(token_end);
    }

    ++line_idx;
    line_start = TextGet(&sample_txs);
    if (unlikely(!line_start)) {
      if (!TextStreamErrcode2(&sample_txs, &reterr)) {
        logerrputs("Error: Only one line in .sample file.\n");
        goto OxSampleToPsam_ret_MALFORMED_INPUT;
      }
      goto OxSampleToPsam_ret_TSTREAM_FAIL;
    }
    token_end = CurTokenEnd(line_start);
    if (unlikely((S_CAST(uintptr_t, token_end - line_start) != 1) || (*line_start != '0'))) {
      goto OxSampleToPsam_ret_INVALID_SAMPLE_HEADER_2;
    }
    linebuf_iter = FirstNonTspace(token_end);
    token_slen = strlen_se(linebuf_iter);
    if (unlikely((token_slen != 1) || (*linebuf_iter != '0'))) {
      goto OxSampleToPsam_ret_INVALID_SAMPLE_HEADER_2;
    }
    linebuf_iter = FirstNonTspace(&(linebuf_iter[1]));
    token_slen = strlen_se(linebuf_iter);
    if (unlikely((token_slen != 1) || (*linebuf_iter != '0'))) {
      goto OxSampleToPsam_ret_INVALID_SAMPLE_HEADER_2;
    }
    linebuf_iter = &(linebuf_iter[1]);

    const uint32_t col_ctl = BitCtToWordCt(col_ct);
    uintptr_t* col_is_categorical;
    uintptr_t* col_is_qt;
    if (unlikely(
            bigstack_calloc_w(col_ctl, &col_is_categorical) ||
            bigstack_calloc_w(col_ctl, &col_is_qt))) {
      goto OxSampleToPsam_ret_NOMEM;
    }
    uint32_t at_least_one_binary_pheno = 0;
    for (uint32_t col_idx = 3; col_idx != col_ct; ++col_idx) {
      linebuf_iter = FirstNonTspace(linebuf_iter);
      unsigned char col_type_char = *linebuf_iter;
      if (unlikely(IsEolnKns(col_type_char))) {
        logerrputs("Error: Second .sample header line has fewer tokens than the first.\n");
        goto OxSampleToPsam_ret_MALFORMED_INPUT;
      }
      if (unlikely(linebuf_iter[1] > ' ')) {
        goto OxSampleToPsam_ret_INVALID_SAMPLE_HEADER_2;
      }
      if (col_idx == sex_col) {
        if (unlikely(col_type_char != 'D')) {
          logerrputs("Error: .sample sex column is not of type 'D'.\n");
          goto OxSampleToPsam_ret_MALFORMED_INPUT;
        }
      } else {
        if ((col_type_char == 'C') || (col_type_char == 'P')) {
          SetBit(col_idx, col_is_qt);
        } else if (col_type_char == 'D') {
          SetBit(col_idx, col_is_categorical);
        } else {
          at_least_one_binary_pheno = 1;
          if (unlikely(col_type_char != 'B')) {
            snprintf(g_logbuf, kLogbufSize, "Error: Unrecognized .sample variable type '%c'.\n", col_type_char);
            goto OxSampleToPsam_ret_MALFORMED_INPUT_2;
          }
        }
      }
      ++linebuf_iter;
    }
    if (at_least_one_binary_pheno) {
      // check for pathological case
      if (unlikely((bsearch_str("0", sorted_mc, 1, max_mc_blen, mc_ct) != -1) || (bsearch_str("1", sorted_mc, 1, max_mc_blen, mc_ct) != -1))) {
        logerrputs("Error: '0' and '1' are unacceptable missing case/control phenotype codes.\n");
        goto OxSampleToPsam_ret_INCONSISTENT_INPUT;
      }
    }
    // to make --data and --data --make-pgen consistent, we do a two-pass load,
    // checking for empty phenotypes in the first pass.
    uintptr_t* col_keep;
    if (unlikely(bigstack_alloc_w(col_ctl, &col_keep))) {
      goto OxSampleToPsam_ret_NOMEM;
    }
    col_keep[0] = 7;
    ZeroWArr(col_ctl - 1, &(col_keep[1]));
    uint32_t uncertain_col_ct = col_ct - 3;
    if (sex_col) {
      // we don't care if sex column is all-NA
      SetBit(sex_col, col_keep);
      --uncertain_col_ct;
    }
    while (uncertain_col_ct) {
      ++line_idx;
      char* line_iter = TextGet(&sample_txs);
      if (!line_iter) {
        if (likely(!TextStreamErrcode2(&sample_txs, &reterr))) {
          break;
        }
        goto OxSampleToPsam_ret_TSTREAM_FAIL;
      }

      const uint32_t old_uncertain_col_ct = uncertain_col_ct;
      uintptr_t col_uidx_base = 0;
      uintptr_t col_keep_inv_bits = ~col_keep[0];
      uint32_t old_col_uidx = 0;
      for (uint32_t uncertain_col_idx = 0; uncertain_col_idx != old_uncertain_col_ct; ++uncertain_col_idx) {
        const uint32_t col_uidx = BitIter0(col_keep, &col_uidx_base, &col_keep_inv_bits);
        line_iter = NextTokenMult(line_iter, col_uidx - old_col_uidx);
        if (unlikely(!line_iter)) {
          goto OxSampleToPsam_ret_MISSING_TOKENS;
        }
        token_end = CurTokenEnd(line_iter);
        token_slen = token_end - line_iter;
        if (bsearch_str(line_iter, sorted_mc, token_slen, max_mc_blen, mc_ct) == -1) {
          SetBit(col_uidx, col_keep);
          --uncertain_col_ct;
        }
        line_iter = token_end;
        old_col_uidx = col_uidx;
      }
    }

    reterr = TextRewind(&sample_txs);
    if (unlikely(reterr)) {
      goto OxSampleToPsam_ret_TSTREAM_FAIL;
    }
    line_idx = 0;

    uint32_t sample_ct_p2 = 0;
    while (1) {
      ++line_idx;
      line_start = TextGet(&sample_txs);
      if (!line_start) {
        if (likely(!TextStreamErrcode2(&sample_txs, &reterr))) {
          break;
        }
        goto OxSampleToPsam_ret_TSTREAM_FAIL;
      }
      ++sample_ct_p2;
      if (sample_ct_p2 < 3) {
        // header lines
        if (sample_ct_p2 == 1) {
          // bugfix (21 Mar 2018): can't set line_iter to nullptr here when
          // col_ct == 3
          char* line_iter = NextTokenMult(line_start, 2);
          line_iter = FirstNonTspace(CurTokenEnd(line_iter));
          for (uint32_t col_idx = 3; col_idx != col_ct; ++col_idx) {
            token_end = CurTokenEnd(line_iter);
            if (IsSet(col_keep, col_idx) && (col_idx != sex_col)) {
              *write_iter++ = '\t';
              write_iter = memcpya(write_iter, line_iter, token_end - line_iter);
            }
            line_iter = FirstNonTspace(token_end);
          }
          AppendBinaryEoln(&write_iter);
          if (unlikely(fwrite_checked(writebuf, write_iter - writebuf, psamfile))) {
            goto OxSampleToPsam_ret_WRITE_FAIL;
          }
          write_iter = writebuf;
        }
        continue;
      }
      if (unlikely(sample_ct_p2 == 0x80000001U)) {
        logerrputs("Error: " PROG_NAME_STR " does not support more than 2^31 - 2 samples.\n");
        goto OxSampleToPsam_ret_MALFORMED_INPUT;
      }

      // FID
      token_end = CurTokenEnd(line_start);
      write_iter = writebuf;
      if (write_fid) {
        write_iter = memcpyax(write_iter, line_start, token_end - line_start, '\t');
      }

      // IID
      linebuf_iter = FirstNonTspace(token_end);
      if (unlikely(IsEolnKns(*linebuf_iter))) {
        goto OxSampleToPsam_ret_MISSING_TOKENS;
      }
      token_end = CurTokenEnd(linebuf_iter);
      write_iter = memcpya(write_iter, linebuf_iter, token_end - linebuf_iter);

      // MISSING
      linebuf_iter = FirstNonTspace(token_end);
      if (unlikely(IsEolnKns(*linebuf_iter))) {
        goto OxSampleToPsam_ret_MISSING_TOKENS;
      }
      token_end = CurTokenEnd(linebuf_iter);

      // flush now since backfilled sex is variable-length ("NA" vs. "1"/"2")
      if (unlikely(fwrite_checked(writebuf, write_iter - writebuf, psamfile))) {
        goto OxSampleToPsam_ret_WRITE_FAIL;
      }
      char* cur_writebuf_start = writebuf;
      write_iter = strcpya_k(writebuf, "\tNA");
      for (uint32_t col_idx = 3; col_idx != col_ct; ++col_idx) {
        linebuf_iter = FirstNonTspace(token_end);
        if (unlikely(IsEolnKns(*linebuf_iter))) {
          goto OxSampleToPsam_ret_MISSING_TOKENS;
        }
        token_end = CurTokenEnd(linebuf_iter);
        if (!IsSet(col_keep, col_idx)) {
          continue;
        }
        token_slen = token_end - linebuf_iter;
        const uint32_t is_missing = (bsearch_str(linebuf_iter, sorted_mc, token_slen, max_mc_blen, mc_ct) != -1);
        if (col_idx == sex_col) {
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
        } else {
          *write_iter++ = '\t';
          if (IsSet(col_is_categorical, col_idx)) {
            if (!is_missing) {
              *write_iter++ = 'C';
              // .sample files are relatively small, so let's go ahead and
              // (i) validate we have a positive integer < 2^31
              // (ii) convert e.g. 9000000, 9000000., 9.0e6 all to 9000000
              double dxx = 0.0;
              char* num_end = ScanadvDouble(linebuf_iter, &dxx);
              int32_t ii = S_CAST(int32_t, dxx);
              if (unlikely((num_end != token_end) || (ii <= 0) || (S_CAST(double, ii) != dxx))) {
                *token_end = '\0';
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid categorical phenotype '%s' on line %" PRIuPTR ", column %u of .sample file (positive integer < 2^31 or --missing-code value expected).\n", linebuf_iter, line_idx, col_idx + 1);
                goto OxSampleToPsam_ret_INCONSISTENT_INPUT_WW;
              }
              write_iter = u32toa(ii, write_iter);
            } else {
              write_iter = memcpya(write_iter, missing_catname, missing_catname_slen);
            }
          } else if (!is_missing) {
            if (IsSet(col_is_qt, col_idx)) {
              double dxx = 0.0;
              if (unlikely(!ScantokDouble(linebuf_iter, &dxx))) {
                *token_end = '\0';
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid quantitative phenotype '%s' on line %" PRIuPTR ", column %u of .sample file (non-infinite number or --missing-code value expected).\n", linebuf_iter, line_idx, col_idx + 1);
                goto OxSampleToPsam_ret_INCONSISTENT_INPUT_WW;
              }
              // used over memcpy to make --data and --data --make-pgen the
              // same (could make that conditional on keep_autoconv?)
              write_iter = dtoa_g(dxx, write_iter);
            } else {
              const uint32_t cc_char_m48 = ctou32(*linebuf_iter) - 48;
              if (likely((token_slen == 1) && (cc_char_m48 < 2))) {
                *write_iter++ = cc_char_m48 + '1';
              } else {
                *token_end = '\0';
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid binary phenotype value '%s' on line %" PRIuPTR ", column %u of .sample file ('0', '1', or --missing-code value expected).\n", linebuf_iter, line_idx, col_idx + 1);
                goto OxSampleToPsam_ret_INCONSISTENT_INPUT_WW;
              }
            }
          } else {
            write_iter = memcpya(write_iter, output_missing_pheno, omp_slen);
          }
        }
      }
      AppendBinaryEoln(&write_iter);
      if (unlikely(fwrite_checked(cur_writebuf_start, write_iter - cur_writebuf_start, psamfile))) {
        goto OxSampleToPsam_ret_WRITE_FAIL;
      }
    }

    // no final writebuf flush since we didn't use usual manual-streaming
    // strategy
    if (unlikely(fclose_null(&psamfile))) {
      goto OxSampleToPsam_ret_WRITE_FAIL;
    }
    const uint32_t sample_ct = sample_ct_p2 - 2;
    if (unlikely(!sample_ct)) {
      logerrputs("Error: No samples in .sample file.\n");
      goto OxSampleToPsam_ret_INCONSISTENT_INPUT;
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
  OxSampleToPsam_ret_MALFORMED_INPUT_2:
    logerrputsb();
  OxSampleToPsam_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  OxSampleToPsam_ret_INVALID_SAMPLE_HEADER_1:
    logerrputs("Error: Invalid first header line in .sample file.\n");
    reterr = kPglRetMalformedInput;
    break;
  OxSampleToPsam_ret_INVALID_SAMPLE_HEADER_2:
    logerrputs("Error: Invalid second header line in .sample file.\n");
    reterr = kPglRetMalformedInput;
    break;
  OxSampleToPsam_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  OxSampleToPsam_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 OxSampleToPsam_ret_1:
  CleanupTextStream2(".sample file", &sample_txs, &reterr);
  fclose_cond(psamfile);
  BigstackReset(bigstack_mark);
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

BoolErr InitOxfordSingleChr(const char* ox_single_chr_str, const char** single_chr_str_ptr, uint32_t* single_chr_slen_ptr, uint32_t* cur_chr_code_ptr, ChrInfo* cip) {
  const uint32_t chr_code_raw = GetChrCodeRaw(ox_single_chr_str);
  if (chr_code_raw == UINT32_MAX) {
    // command-line parser guarantees that allow_extra_chrs is true here
    const uint32_t chr_slen = strlen(ox_single_chr_str);
    if (single_chr_str_ptr) {
      *single_chr_str_ptr = ox_single_chr_str;
      *single_chr_slen_ptr = chr_slen;
      return 0;
    }
    return (TryToAddChrName(ox_single_chr_str, "--bgen file", 0, chr_slen, 1, cur_chr_code_ptr, cip) != kPglRetSuccess);
  }
  uint32_t chr_code = chr_code_raw;
  if (chr_code > cip->max_code) {
    if (chr_code < kMaxContigs) {
      logerrputs("Error: --oxford-single-chr chromosome code is not in the chromosome set.\n");
      return 1;
    }
    chr_code = cip->xymt_codes[chr_code - kMaxContigs];
    if (IsI32Neg(chr_code)) {
      logerrputs("Error: --oxford-single-chr chromosome code is not in the chromosome set.\n");
      return 1;
    }
  }
  if (!IsSet(cip->chr_mask, chr_code)) {
    logerrputs("Error: --oxford-single-chr chromosome code is excluded by chromosome filter.\n");
    return 1;
  }
  // bugfix (19 Mar 2018): need a not here...
  if (!single_chr_str_ptr) {
    *cur_chr_code_ptr = chr_code;
    return 0;
  }
  // can't fail due to OxGenToPgen()'s writebuf allocation logic
  char* chr_buf = S_CAST(char*, bigstack_alloc_raw(kCacheline));
  char* chr_name_end = chrtoa(cip, chr_code, chr_buf);
  *single_chr_str_ptr = chr_buf;
  *single_chr_slen_ptr = chr_name_end - chr_buf;
  return 0;
}

static_assert(sizeof(Dosage) == 2, "OxGenToPgen() needs to be updated.");
PglErr OxGenToPgen(const char* genname, const char* samplename, const char* ox_single_chr_str, const char* ox_missing_code, MiscFlags misc_flags, ImportFlags import_flags, OxfordImportFlags oxford_import_flags, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, double import_dosage_certainty, uint32_t max_thread_ct, char* outname, char* outname_end, ChrInfo* cip) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* pvarfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  uintptr_t line_idx = 0;
  TextStream gen_txs;
  STPgenWriter spgw;
  PreinitTextStream(&gen_txs);
  PreinitSpgw(&spgw);
  {
    uint32_t sample_ct;
    reterr = OxSampleToPsam(samplename, ox_missing_code, import_flags, outname, outname_end, &sample_ct);
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
    reterr = InitTextStream(genname, max_line_blen, decompress_thread_ct, &gen_txs);
    if (unlikely(reterr)) {
      if (reterr == kPglRetOpenFail) {
        const uint32_t slen = strlen(genname);
        if ((!StrEndsWith(genname, ".gen", slen)) &&
            (!StrEndsWith(genname, ".gen.gz", slen))) {
          logerrprintfww("Error: Failed to open %s : %s. (--gen expects a complete filename; did you forget '.gen' at the end?)\n", genname, strerror(errno));
          goto OxGenToPgen_ret_1;
        }
      }
      goto OxGenToPgen_ret_TSTREAM_FAIL;
    }
    const uint32_t allow_extra_chrs = (misc_flags / kfMiscAllowExtraChrs) & 1;
    FinalizeChrset(misc_flags, cip);

    char* writebuf;
    if (unlikely(bigstack_alloc_c(kMaxMediumLine + max_line_blen, &writebuf))) {
      // shouldn't actually be possible, but may as well defend against changes
      // to how RLstream allocation works
      goto OxGenToPgen_ret_NOMEM;
    }
    char* writebuf_flush = &(writebuf[kMaxMediumLine]);

    const char* single_chr_str = nullptr;
    uint32_t single_chr_slen = 0;
    if (ox_single_chr_str) {
      if (unlikely(InitOxfordSingleChr(ox_single_chr_str, &single_chr_str, &single_chr_slen, nullptr, cip))) {
        goto OxGenToPgen_ret_INVALID_CMDLINE;
      }
    }

    snprintf(outname_end, kMaxOutfnameExtBlen, ".pvar");
    if (unlikely(fopen_checked(outname, FOPEN_WB, &pvarfile))) {
      goto OxGenToPgen_ret_OPEN_FAIL;
    }
    char* write_iter = writebuf;
    if (cip->chrset_source) {
      AppendChrsetLine(cip, &write_iter);
    }
    write_iter = strcpya_k(write_iter, "#CHROM\tPOS\tID\tREF\tALT" EOLN_STR);

    // Explicit 32768 instead of kDosageMax since this is driven by the BGEN
    // 1.1 format, not plink2's dosage representation.
    // Note that command-line parser multiplies import_dosage_certainty by
    // (1 - kSmallEpsilon), and we want import_dosage_certainty_int to be 1
    // when import_dosage_certainty is zero.
    uint32_t import_dosage_certainty_int = 1 + S_CAST(int32_t, import_dosage_certainty * 32768);
    const uint32_t dosage_int_sum_thresh = 3 * (import_dosage_certainty_int - 1);

    const uint32_t dosage_erase_halfdist = kDosage4th - dosage_erase_thresh;
    const uint32_t prov_ref_allele_second = !(oxford_import_flags & kfOxfordImportRefFirst);
    uint32_t dosage_is_present = 0;
    uint32_t variant_ct = 0;
    uintptr_t variant_skip_ct = 0;
    char* line_iter = TextLineEnd(&gen_txs);
    ++line_idx;
    for (; TextGetUnsafe2(&gen_txs, &line_iter); ++line_idx) {
      char* chr_code_str = line_iter;
      char* chr_code_end = CurTokenEnd(chr_code_str);
      const char* variant_id_str = FirstNonTspace(chr_code_end);
      if (unlikely(IsEolnKns(*variant_id_str))) {
        goto OxGenToPgen_ret_MISSING_TOKENS;
      }

      if (!single_chr_str) {
        uint32_t cur_chr_code;
        reterr = GetOrAddChrCodeDestructive(".gen file", line_idx, allow_extra_chrs, chr_code_str, chr_code_end, cip, &cur_chr_code);
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
        write_iter = chrtoa(cip, cur_chr_code, write_iter);
      } else {
        write_iter = memcpya(write_iter, single_chr_str, single_chr_slen);
      }
      *write_iter++ = '\t';
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

      write_iter = u32toa_x(cur_bp, '\t', write_iter);
      const uint32_t variant_id_slen = variant_id_end - variant_id_str;
      if (unlikely(variant_id_slen > kMaxIdSlen)) {
        putc_unlocked('\n', stdout);
        logerrputs("Error: Variant names are limited to " MAX_ID_SLEN_STR " characters.\n");
        goto OxGenToPgen_ret_MALFORMED_INPUT;
      }
      write_iter = memcpyax(write_iter, variant_id_str, variant_id_slen, '\t');

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
        write_iter = memcpyax(write_iter, first_allele_str, first_allele_end - first_allele_str, '\t');
        write_iter = memcpya(write_iter, second_allele_str, linebuf_iter - second_allele_str);
      } else {
        write_iter = memcpyax(write_iter, second_allele_str, linebuf_iter - second_allele_str, '\t');
        write_iter = memcpya(write_iter, first_allele_str, first_allele_end - first_allele_str);
      }
      AppendBinaryEoln(&write_iter);
      if (unlikely(fwrite_ck(writebuf_flush, pvarfile, &write_iter))) {
        goto OxGenToPgen_ret_WRITE_FAIL;
      }

      if (!dosage_is_present) {
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
          dosage_is_present = Bgen11DosageImportCheck(dosage_int_sum_thresh, import_dosage_certainty_int, dosage_erase_halfdist, dosage_int0, dosage_int1, dosage_int2);
          if (dosage_is_present) {
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
    if (unlikely(fclose_flush_null(writebuf_flush, write_iter, &pvarfile))) {
      goto OxGenToPgen_ret_WRITE_FAIL;
    }
    if (unlikely(!variant_ct)) {
      if (!variant_skip_ct) {
        logerrputs("Error: Empty .gen file.\n");
        goto OxGenToPgen_ret_INCONSISTENT_INPUT;
      }
      logerrprintfww("Error: All %" PRIuPTR " variant%s in .gen file skipped due to chromosome filter.\n", variant_skip_ct, (variant_skip_ct == 1)? "" : "s");
      goto OxGenToPgen_ret_INCONSISTENT_INPUT;
    }
    logprintf("--data/--gen: %u variant%s scanned%s.\n", variant_ct, (variant_ct == 1)? "" : "s", dosage_is_present? "" : " (all hardcalls)");

    // second pass
    BigstackReset(writebuf);
    if (TextIsMt(&gen_txs) && (decompress_thread_ct > 1)) {
      // Close and reopen, so that we can reduce decompress_thread_ct to 1.
      char* dst = R_CAST(char*, TextStreamMemStart(&gen_txs));
      if (unlikely(CleanupTextStream(&gen_txs, &reterr))) {
        logerrprintfww(kErrprintfFread, ".gen file", strerror(errno));
        goto OxGenToPgen_ret_1;
      }
      reterr = TextStreamOpenEx(genname, kMaxLongLine, max_line_blen, 1, nullptr, dst, &gen_txs);
      if (unlikely(reterr)) {
        goto OxGenToPgen_ret_TSTREAM_REWIND_FAIL;
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
    reterr = SpgwInitPhase1(outname, nullptr, nullptr, variant_ct, sample_ct, dosage_is_present? kfPgenGlobalDosagePresent : kfPgenGlobal0, (oxford_import_flags & (kfOxfordImportRefFirst | kfOxfordImportRefLast))? 1 : 2, &spgw, &spgw_alloc_cacheline_ct, &max_vrec_len);
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
    uintptr_t* dosage_present;
    // if we weren't using bigstack_alloc, this would need to be sample_ctaw2
    if (unlikely(
            bigstack_alloc_w(sample_ctl2, &genovec) ||
            bigstack_alloc_w(sample_ctl, &dosage_present))) {
      goto OxGenToPgen_ret_NOMEM;
    }
    Dosage* dosage_main = nullptr;
    if (dosage_is_present) {
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
        goto OxGenToPgen_ret_TSTREAM_REWIND_FAIL;
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
      const char* linebuf_iter = NextTokenMult(FirstNonTspace(&(chr_code_end[1])), 4);
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
        R_CAST(Halfword*, dosage_present)[widx] = dosage_present_hw;
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
        reterr = SpgwAppendBiallelicGenovecDosage16(genovec, dosage_present, dosage_main, dosage_ct, &spgw);
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
    SpgwFinish(&spgw);
    putc_unlocked('\r', stdout);
    write_iter = strcpya_k(g_logbuf, "--data/--gen: ");
    const uint32_t outname_base_slen = outname_end - outname;
    write_iter = memcpya(write_iter, outname, outname_base_slen + 5);
    write_iter = strcpya_k(write_iter, " + ");
    write_iter = memcpya(write_iter, outname, outname_base_slen);
    write_iter = strcpya_k(write_iter, ".pvar");
    snprintf(write_iter, kLogbufSize - 2 * kPglFnamesize - 64, " written.\n");
    WordWrapB(0);
    logputsb();
  }
  while (0) {
  OxGenToPgen_ret_TSTREAM_FAIL:
    TextStreamErrPrint(".gen file", &gen_txs);
    break;
  OxGenToPgen_ret_TSTREAM_REWIND_FAIL:
    TextStreamErrPrintRewind(".gen file", &gen_txs, &reterr);
    break;
  OxGenToPgen_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  OxGenToPgen_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  OxGenToPgen_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  OxGenToPgen_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
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
  OxGenToPgen_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 OxGenToPgen_ret_1:
  CleanupSpgw(&spgw, &reterr);
  CleanupTextStream2(".gen file", &gen_txs, &reterr);
  fclose_cond(pvarfile);
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

  uint32_t dosage_is_present;
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
      ctx->dosage_is_present = 1;
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
    uintptr_t* write_dosage_present_iter = &(ctx->write_dosage_presents[parity][vidx * sample_ctaw]);
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
        R_CAST(Halfword*, write_dosage_present_iter)[widx] = dosage_present_hw;
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
      write_dosage_present_iter = &(write_dosage_present_iter[sample_ctaw]);
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

#ifdef __arm__
#  error "Unaligned accesses in Bgen13GetOneVal()."
#endif
HEADER_INLINE uintptr_t Bgen13GetOneVal(const unsigned char* prob_start, uint64_t prob_offset, uint32_t bit_precision, uintptr_t numer_mask) {
  const uint64_t bit_offset = prob_offset * bit_precision;
  uint64_t relevant_bits;
  // This can read slightly past the end of the buffer.
  memcpy(&relevant_bits, &(prob_start[bit_offset / CHAR_BIT]), sizeof(int64_t));
  return (relevant_bits >> (bit_offset % CHAR_BIT)) & numer_mask;
}

#ifdef __arm__
#  error "Unaligned accesses in Bgen13GetTwoVals()."
#endif
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

  // high 32 bits = tidx for now (smaller takes precedence)
  // low 32 bits = uint32_t(PglErr)
  uint64_t err_info;

  uint32_t dosage_is_present;
} Bgen13DosageOrPhaseScanCtx;

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
              // possible todo: report variant index
              goto Bgen13DosageOrPhaseScanThread_malformed;
            }
          } else {
            const uintptr_t extracted_byte_ct = ZSTD_decompress(K_CAST(unsigned char*, cur_uncompressed_geno), uncompressed_byte_ct, compressed_geno_start, compressed_byte_ct);
            if (unlikely(extracted_byte_ct != uncompressed_byte_ct)) {
              // possible todo: inspect error code
              goto Bgen13DosageOrPhaseScanThread_malformed;
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
          goto Bgen13DosageOrPhaseScanThread_malformed;
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
          goto Bgen13DosageOrPhaseScanThread_malformed;
        }
        if (unlikely(max_ploidy > 2)) {
          goto Bgen13DosageOrPhaseScanThread_not_yet_supported;
        }
        const unsigned char* missing_and_ploidy_info = &(cur_uncompressed_geno[8]);
        const unsigned char* probs_start = &(cur_uncompressed_geno[10 + sample_ct]);
        const uint32_t is_phased = probs_start[-2];
        if (unlikely(is_phased > 1)) {
          goto Bgen13DosageOrPhaseScanThread_malformed;
        }
        const uint32_t bit_precision = probs_start[-1];
        if (unlikely((!bit_precision) || (bit_precision > 32))) {
          goto Bgen13DosageOrPhaseScanThread_malformed;
        }
        if (unlikely(bit_precision > kMaxBgenImportBits)) {
          goto Bgen13DosageOrPhaseScanThread_not_yet_supported;
        }

        if (cur_allele_ct != 2) {
          // shouldn't be possible to get here for now
          assert(0);
          goto Bgen13DosageOrPhaseScanThread_not_yet_supported;
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
            goto Bgen13DosageOrPhaseScanThread_malformed;
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
              const uint32_t full_vec_ct = sample_ct / (kBytesPerVec / 2);
              // If all bytes are equal to 0 or numer_mask, all calls in this
              // block must be hardcall/missing.
              const VecUc vec0 = vecuc_setzero();
              const VecUc vecmax = vecuc_set1(numer_mask);
              uint32_t vec_idx = 0;
              for (; vec_idx != full_vec_ct; ++vec_idx) {
                const VecUc cur_vec = vecuc_loadu(&(probs_start[vec_idx * kBytesPerVec]));
                const VecUc safe_bytes = (cur_vec == vec0) | (cur_vec == vecmax);
                if (vecuc_movemask(safe_bytes) != kVec8thUintMax) {
                  break;
                }
              }
              sample_idx = vec_idx * (kBytesPerVec / 2);
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
                if (vecu16_movemask(safe_u16s) != kVec8thUintMax) {
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
                goto Bgen13DosageOrPhaseScanThread_malformed;
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
                  goto Bgen13DosageOrPhaseScanThread_malformed;
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
                  goto Bgen13DosageOrPhaseScanThread_malformed;
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
              goto Bgen13DosageOrPhaseScanThread_malformed;
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
                goto Bgen13DosageOrPhaseScanThread_malformed;
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
                goto Bgen13DosageOrPhaseScanThread_malformed;
              }
              prob_offset += missing_and_ploidy;
            }
          }
        }
      }
    }
    while (0) {
    Bgen13DosageOrPhaseScanThread_malformed:
      {
        const uint64_t new_err_info = (S_CAST(uint64_t, tidx) << 32) | S_CAST(uint32_t, kPglRetMalformedInput);
        UpdateU64IfSmaller(new_err_info, &ctx->err_info);
      }
      break;
    Bgen13DosageOrPhaseScanThread_not_yet_supported:
      {
        const uint64_t new_err_info = (S_CAST(uint64_t, tidx) << 32) | S_CAST(uint32_t, kPglRetNotYetSupported);
        UpdateU64IfSmaller(new_err_info, &ctx->err_info);
      }
      break;
    Bgen13DosageOrPhaseScanThread_found:
      ctx->dosage_is_present = 1;
      break;
    }
    parity = 1 - parity;
  } while (!THREAD_BLOCK_FINISH(arg));
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
  uint32_t hard_call_halfdist;
  uint32_t* bgen_import_dosage_certainty_thresholds;
  uint32_t prov_ref_allele_second;

  unsigned char** thread_wkspaces;
  uint32_t* thread_bidxs[2];
  GparseRecord* gparse[2];
  const uintptr_t* block_allele_idx_offsets[2];

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
  uintptr_t* phasepresent = nullptr;
  uintptr_t* phaseinfo = nullptr;
  uintptr_t* dosage_present = nullptr;
  Dosage* dosage_main = nullptr;
  uintptr_t* dphase_present = nullptr;
  SDosage* dphase_delta = nullptr;
  uint32_t cur_allele_ct = 2;
  uint32_t parity = 0;
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
              // possible todo: report variant index
              goto Bgen13GenoToPgenThread_malformed;
            }
          } else {
            const uintptr_t extracted_byte_ct = ZSTD_decompress(K_CAST(unsigned char*, cur_uncompressed_geno), uncompressed_byte_ct, grp->record_start, compressed_byte_ct);
            if (unlikely(extracted_byte_ct != uncompressed_byte_ct)) {
              // possible todo: inspect error code
              goto Bgen13GenoToPgenThread_malformed;
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
          goto Bgen13GenoToPgenThread_malformed;
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
          goto Bgen13GenoToPgenThread_malformed;
        }
        if (unlikely(max_ploidy > 2)) {
          goto Bgen13GenoToPgenThread_not_yet_supported;
        }
        const unsigned char* missing_and_ploidy_iter = &(cur_uncompressed_geno[8]);
        const unsigned char* probs_start = &(cur_uncompressed_geno[10 + sample_ct]);
        const uint32_t is_phased = probs_start[-2];
        if (unlikely(is_phased > 1)) {
          goto Bgen13GenoToPgenThread_malformed;
        }
        const uint32_t bit_precision = probs_start[-1];
        if (unlikely((!bit_precision) || (bit_precision > 32))) {
          goto Bgen13GenoToPgenThread_malformed;
        }
        if (unlikely(bit_precision > kMaxBgenImportBits)) {
          goto Bgen13GenoToPgenThread_not_yet_supported;
        }
        if (cur_allele_ct != 2) {
          // shouldn't be possible to get here for now
          assert(0);
          goto Bgen13GenoToPgenThread_not_yet_supported;
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
            goto Bgen13GenoToPgenThread_malformed;
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
                    goto Bgen13GenoToPgenThread_malformed;
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
                    goto Bgen13GenoToPgenThread_malformed;
                  }
                  uintptr_t numer_aa;
                  uintptr_t numer_ab;
                  Bgen13GetTwoVals(probs_start, prob_offset, bit_precision, numer_mask, &numer_aa, &numer_ab);
                  prob_offset += 2;
                  if (unlikely(numer_aa + numer_ab > numer_mask)) {
                    goto Bgen13GenoToPgenThread_malformed;
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
                      goto Bgen13GenoToPgenThread_malformed;
                    }
                    prob_offset += missing_and_ploidy;
                  Bgen13GenoToPgenThread_generic_unphased_missing:
                    genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
                    continue;
                  }
                  if (unlikely(prob_offset >= prob_offset_end)) {
                    goto Bgen13GenoToPgenThread_malformed;
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
                    goto Bgen13GenoToPgenThread_malformed;
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
                      goto Bgen13GenoToPgenThread_malformed;
                    }
                    prob_offset += missing_and_ploidy;
                  Bgen13GenoToPgenThread_generic_phased_missing:
                    genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
                    continue;
                  }
                  if (unlikely(prob_offset >= prob_offset_end)) {
                    goto Bgen13GenoToPgenThread_malformed;
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
          // do we actually need this?  well, play it safe for now
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
    while (0) {
    Bgen13GenoToPgenThread_malformed:
      {
        const uint64_t new_err_info = (S_CAST(uint64_t, tidx) << 32) | S_CAST(uint32_t, kPglRetMalformedInput);
        UpdateU64IfSmaller(new_err_info, &ctx->err_info);
      }
      break;
    Bgen13GenoToPgenThread_not_yet_supported:
      {
        const uint64_t new_err_info = (S_CAST(uint64_t, tidx) << 32) | S_CAST(uint32_t, kPglRetNotYetSupported);
        UpdateU64IfSmaller(new_err_info, &ctx->err_info);
      }
      break;
    }
    parity = 1 - parity;
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

static_assert(sizeof(Dosage) == 2, "OxBgenToPgen() needs to be updated.");
PglErr OxBgenToPgen(const char* bgenname, const char* samplename, const char* const_fid, const char* ox_single_chr_str, const char* ox_missing_code, MiscFlags misc_flags, ImportFlags import_flags, OxfordImportFlags oxford_import_flags, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, double import_dosage_certainty, char id_delim, char idspace_to, uint32_t max_thread_ct, char* outname, char* outname_end, ChrInfo* cip) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  FILE* bgenfile = nullptr;

  // only if no sample file specified, and .bgen has sample IDs.  (possible
  // todo: consistency check when both sources of sample IDs are present?)
  FILE* psamfile = nullptr;

  FILE* pvarfile = nullptr;
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
    const uint32_t raw_variant_ct = initial_uints[2];
    if (unlikely(!raw_variant_ct)) {
      logerrputs("Error: Empty .bgen file.\n");
      goto OxBgenToPgen_ret_INCONSISTENT_INPUT;
    }
    const uint32_t sample_ct = initial_uints[3];
    if (unlikely(initial_uints[4] && (initial_uints[4] != 0x6e656762))) {
      logerrputs("Error: Invalid .bgen magic number.\n");
      goto OxBgenToPgen_ret_MALFORMED_INPUT;
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
    logprintf("--bgen: %u variant%s detected, format v1.%c.\n", raw_variant_ct, (raw_variant_ct == 1)? "" : "s", (layout == 1)? '1' : ((compression_mode == 2)? '3' : '2'));
    if (samplename[0]) {
      uint32_t sfile_sample_ct;
      reterr = OxSampleToPsam(samplename, ox_missing_code, import_flags, outname, outname_end, &sfile_sample_ct);
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
        if (unlikely(
                (!fread_unlocked(&sample_id_block_byte_ct, 4, 1, bgenfile)) ||
                (!fread_unlocked(&sample_id_block_entry_ct, 4, 1, bgenfile)))) {
          goto OxBgenToPgen_ret_READ_FAIL;
        }
        if (unlikely(
                (S_CAST(uint64_t, sample_id_block_byte_ct) + initial_uints[1] > initial_uints[0]) ||
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
      if (unlikely(
              (!fread_unlocked(&sample_id_block_byte_ct, 4, 1, bgenfile)) ||
              (!fread_unlocked(&sample_id_block_entry_ct, 4, 1, bgenfile)))) {
        goto OxBgenToPgen_ret_READ_FAIL;
      }
      if (unlikely(
              (sample_id_block_byte_ct < 8) ||
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
      InitImportSampleIdContext(const_fid, misc_flags, import_flags, id_delim, &isic);
      uint32_t write_fid = 1;
      uint32_t write_sid = 0;
      if (id_delim) {
        sample_id_block_iter = sample_id_block_main;
        uint32_t nonzero_first_field_observed = 0;
        for (uint32_t sample_idx = 0; ; ) {
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
          if (!nonzero_first_field_observed) {
            nonzero_first_field_observed = (first_delim != &(sample_id_iter[1])) || (sample_id_iter[0] != '0');
          }
          sample_id_iter = &(first_delim[1]);
          if (memchr(sample_id_iter, ctou32(id_delim), sample_id_end - sample_id_iter)) {
            isic.two_part_null_fid = isic.id_delim_sid;
            isic.two_part_null_sid = 1 - isic.id_delim_sid;
            write_fid = 1;
            write_sid = 1;
            if (!nonzero_first_field_observed) {
              sample_id_block_iter = sample_id_end;
              while (++sample_idx != sample_ct) {
                uint16_t new_input_id_slen;
                memcpy_k(&new_input_id_slen, sample_id_block_iter, sizeof(int16_t));
                sample_id_iter = &(sample_id_block_iter[2]);
                // We actually need to error out on new_input_id_slen < 2, but
                // that'll happen in the next loop.
                if ((new_input_id_slen < 2) || (sample_id_iter[0] != '0') || (sample_id_iter[1] != id_delim)) {
                  nonzero_first_field_observed = 1;
                  break;
                }
                sample_id_block_iter = &(sample_id_iter[new_input_id_slen]);
              }
            }
            break;
          }
          sample_id_block_iter = sample_id_end;
          if (++sample_idx == sample_ct) {
            if (isic.id_delim_sid) {
              write_fid = 0;
              write_sid = 1;
            }
            break;
          }
        }
        if (!nonzero_first_field_observed) {
          write_fid = 0;
          isic.omit_fid_0 = 1;
        }
      } else if ((!isic.const_fid) && (!isic.double_id)) {
        write_fid = 0;
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
    const uint32_t allow_extra_chrs = (misc_flags / kfMiscAllowExtraChrs) & 1;
    FinalizeChrset(misc_flags, cip);
    const uint32_t autosome_ct_p1 = cip->autosome_ct + 1;
    uint32_t chr_filter_present = (PopcountBitRange(cip->chr_mask, 0, autosome_ct_p1) != autosome_ct_p1) || (allow_extra_chrs && (cip->is_include_stack || cip->incl_excl_name_stack));
    if (!chr_filter_present) {
      for (uint32_t xymt_idx = 0; xymt_idx != kChrOffsetCt; ++xymt_idx) {
        if (!IsI32Neg(cip->xymt_codes[xymt_idx])) {
          if (!IsSet(cip->chr_mask, autosome_ct_p1 + xymt_idx)) {
            chr_filter_present = 1;
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

    char* writebuf;
    if (unlikely(bigstack_alloc_c(2 * kMaxMediumLine + kCacheline, &writebuf))) {
      goto OxBgenToPgen_ret_NOMEM;
    }
    char* writebuf_flush = &(writebuf[kMaxMediumLine]);

    uint32_t cur_chr_code = 0;
    if (ox_single_chr_str) {
      if (unlikely(InitOxfordSingleChr(ox_single_chr_str, nullptr, nullptr, &cur_chr_code, cip))) {
        goto OxBgenToPgen_ret_INVALID_CMDLINE;
      }
    }
    snprintf(outname_end, kMaxOutfnameExtBlen, ".pvar");
    if (unlikely(fopen_checked(outname, FOPEN_WB, &pvarfile))) {
      goto OxBgenToPgen_ret_OPEN_FAIL;
    }
    char* write_iter = writebuf;
    if (cip->chrset_source) {
      AppendChrsetLine(cip, &write_iter);
    }
    write_iter = strcpya_k(write_iter, "#CHROM\tPOS\tID\tREF\tALT" EOLN_STR);

    const uint32_t snpid_chr = (oxford_import_flags & kfOxfordImportBgenSnpIdChr);

    // true for both provisional-reference and real-reference second
    const uint32_t prov_ref_allele_second = !(oxford_import_flags & kfOxfordImportRefFirst);

    if (hard_call_thresh == UINT32_MAX) {
      hard_call_thresh = kDosageMid / 10;
    }
    const uint32_t sample_ctaw2 = NypCtToAlignedWordCt(sample_ct);
    const uint32_t sample_ctaw = BitCtToAlignedWordCt(sample_ct);
    uint32_t dosage_is_present = 0;
    common.sample_ct = sample_ct;
    common.dosage_erase_halfdist = kDosage4th - dosage_erase_thresh;
    common.compression_mode = compression_mode;
    if (layout == 1) {
      // v1.1
      // this block belongs in its own function...
      Bgen11DosageScanCtx scan_ctx;
      scan_ctx.common = &common;
      scan_ctx.dosage_is_present = 0;
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
      if (calc_thread_ct > raw_variant_ct) {
        calc_thread_ct = raw_variant_ct;
      }
      if (unlikely(
              bigstack_alloc_u16p(calc_thread_ct, &scan_ctx.bgen_geno_bufs))) {
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
      if (unlikely(
              bigstack_alloc_uc(bgen_geno_max_byte_ct * main_block_size, &(compressed_geno_bufs[0])) ||
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
      uint32_t variant_ct = 0;
      uint32_t block_vidx = 0;
      uint32_t cur_block_size = calc_thread_ct;
      uint32_t parity = 0;
      uintptr_t compressed_block_byte_ct = 6LU * sample_ct;
      unsigned char** compressed_geno_starts = common.compressed_geno_starts[0];
      unsigned char* bgen_geno_iter = compressed_geno_bufs[0];
      uint32_t skip = 0;
      for (uint32_t variant_uidx = 0; variant_uidx != raw_variant_ct; ) {
        uint32_t uii;
        if (unlikely(!fread_unlocked(&uii, 4, 1, bgenfile))) {
          goto OxBgenToPgen_ret_READ_FAIL;
        }
        if (unlikely(uii != sample_ct)) {
          logputs("\n");
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
            logputs("\n");
            logerrputs("Error: Length-0 SNP ID in .bgen file.\n");
            goto OxBgenToPgen_ret_INCONSISTENT_INPUT;
          }
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
            if (unlikely(!chr_name_slen)) {
              logputs("\n");
              logerrputs("Error: Length-0 chromosome ID in .bgen file.\n");
              goto OxBgenToPgen_ret_INCONSISTENT_INPUT;
            }
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
            chr_name_slen = snpid_slen;
          }
          reterr = GetOrAddChrCodeDestructive("--bgen file", 0, allow_extra_chrs, R_CAST(char*, loadbuf), R_CAST(char*, &(loadbuf[chr_name_slen])), cip, &cur_chr_code);
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
        if (dosage_is_present || skip) {
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
            dosage_is_present = scan_ctx.dosage_is_present;
            if (dosage_is_present) {
              // don't need to scan for any more dosages
              StopThreads(&tg);
              if (!chr_filter_present) {
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

      if (!chr_filter_present) {
        variant_ct = raw_variant_ct;
      } else {
        variant_ct += block_vidx;
        if (variant_ct < calc_thread_ct) {
          if (unlikely(!variant_ct)) {
            logputs("\n");
            logerrprintfww("Error: All %u variant%s in .bgen file skipped due to chromosome filter.\n", raw_variant_ct, (raw_variant_ct == 1)? "" : "s");
            goto OxBgenToPgen_ret_INCONSISTENT_INPUT;
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
        if (block_vidx && (!scan_ctx.dosage_is_present)) {
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
        dosage_is_present = scan_ctx.dosage_is_present;
      }

      if (unlikely(fseeko(bgenfile, initial_uints[0] + 4, SEEK_SET))) {
        goto OxBgenToPgen_ret_READ_FAIL;
      }
      snprintf(outname_end, kMaxOutfnameExtBlen, ".pgen");
      uintptr_t spgw_alloc_cacheline_ct;
      uint32_t max_vrec_len;
      reterr = SpgwInitPhase1(outname, nullptr, nullptr, variant_ct, sample_ct, dosage_is_present? kfPgenGlobalDosagePresent : kfPgenGlobal0, (oxford_import_flags & (kfOxfordImportRefFirst | kfOxfordImportRefLast))? 1 : 2, &spgw, &spgw_alloc_cacheline_ct, &max_vrec_len);
      if (unlikely(reterr)) {
        if (reterr == kPglRetOpenFail) {
          logputs("\n");
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
              logputs("\n");
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
                logputs("\n");
                logerrputs("Error: Length-0 SNP ID in .bgen file.\n");
                goto OxBgenToPgen_ret_INCONSISTENT_INPUT;
              }
              if (unlikely(!fread_unlocked(loadbuf, snpid_slen, 1, bgenfile))) {
                goto OxBgenToPgen_ret_READ_FAIL;
              }
              loadbuf[snpid_slen] = '\0';
              rsid_start = R_CAST(char*, &(loadbuf[snpid_slen + 1]));
            }
            uint16_t rsid_slen;
            if (unlikely(!fread_unlocked(&rsid_slen, 2, 1, bgenfile))) {
              goto OxBgenToPgen_ret_READ_FAIL;
            }
            if (unlikely(!rsid_slen)) {
              logputs("\n");
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
                  logputs("\n");
                  logerrputs("Error: Length-0 chromosome ID in .bgen file.\n");
                  goto OxBgenToPgen_ret_INCONSISTENT_INPUT;
                }
                if (unlikely(!fread_unlocked(chr_name_start, chr_name_slen, 1, bgenfile))) {
                  goto OxBgenToPgen_ret_READ_FAIL;
                }
                if (strequal_k(chr_name_start, "NA", chr_name_slen)) {
                  strcpy_k(chr_name_start, "0");
                  chr_name_slen = 1;
                } else {
                  chr_name_start[chr_name_slen] = '\0';
                }
              } else {
                if (unlikely(fseeko(bgenfile, chr_name_slen, SEEK_CUR))) {
                  goto OxBgenToPgen_ret_READ_FAIL;
                }
                chr_name_start = R_CAST(char*, loadbuf);
                chr_name_slen = snpid_slen;
              }
              reterr = GetOrAddChrCodeDestructive("--bgen file", 0, allow_extra_chrs, chr_name_start, &(chr_name_start[chr_name_slen]), cip, &cur_chr_code);
              if (unlikely(reterr)) {
                goto OxBgenToPgen_ret_1;
              }
              skip = !IsSet(cip->chr_mask, cur_chr_code);
            }

            uint32_t cur_bp;
            uint32_t a1_slen;
            if (unlikely(
                    (!fread_unlocked(&cur_bp, 4, 1, bgenfile)) ||
                    (!fread_unlocked(&a1_slen, 4, 1, bgenfile)))) {
              goto OxBgenToPgen_ret_READ_FAIL;
            }
            if (skip) {
              uint32_t a2_slen;
              if (unlikely(
                      fseeko(bgenfile, a1_slen, SEEK_CUR) ||
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
              logputs("\n");
              logerrputs("Error: Empty allele code in .bgen file.\n");
              goto OxBgenToPgen_ret_MALFORMED_INPUT;
            }
            if (unlikely(a1_slen > 1000000000)) {
              logputs("\n");
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
              logputs("\n");
              logerrputs("Error: Empty allele code in .bgen file.\n");
              goto OxBgenToPgen_ret_MALFORMED_INPUT;
            }
            if (unlikely(a2_slen > 1000000000)) {
              logputs("\n");
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
            write_iter = chrtoa(cip, cur_chr_code, write_iter);
            *write_iter++ = '\t';
            if (unlikely(cur_bp > 0x7ffffffe)) {
              logputs("\n");
              logerrputs("Error: Invalid bp coordinate (> 2^31 - 2) in .bgen file\n");
              goto OxBgenToPgen_ret_MALFORMED_INPUT;
            }
            write_iter = u32toa_x(cur_bp, '\t', write_iter);
            write_iter = memcpyax(write_iter, rsid_start, rsid_slen, '\t');
            if (prov_ref_allele_second) {
              uint32_t swap_slen = a1_slen;
              a1_slen = a2_slen;
              a2_slen = swap_slen;
              char* swap_ptr = a1_ptr;
              a1_ptr = a2_ptr;
              a2_ptr = swap_ptr;
            }
            if ((write_iter >= writebuf_flush) || (a1_slen >= kMaxMediumLine)) {
              if (unlikely(fwrite_checked(writebuf, write_iter - writebuf, pvarfile))) {
                goto OxBgenToPgen_ret_WRITE_FAIL;
              }
              write_iter = writebuf;
            }
            if (a1_slen < kMaxMediumLine) {
              write_iter = memcpya(write_iter, a1_ptr, a1_slen);
            } else {
              if (unlikely(fwrite_checked(a1_ptr, a1_slen, pvarfile))) {
                goto OxBgenToPgen_ret_WRITE_FAIL;
              }
            }
            *write_iter++ = '\t';
            if ((write_iter >= writebuf_flush) || (a2_slen >= kMaxMediumLine)) {
              if (unlikely(fwrite_checked(writebuf, write_iter - writebuf, pvarfile))) {
                goto OxBgenToPgen_ret_WRITE_FAIL;
              }
              write_iter = writebuf;
            }
            if (a2_slen < kMaxMediumLine) {
              write_iter = memcpya(write_iter, a2_ptr, a2_slen);
            } else {
              if (unlikely(fwrite_checked(a2_ptr, a2_slen, pvarfile))) {
                goto OxBgenToPgen_ret_WRITE_FAIL;
              }
            }
            AppendBinaryEoln(&write_iter);

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
      // ***** all bigstack allocations from this point on are reset before
      //       pass 2 *****
      // probably want to change this to use Gparse...
      uintptr_t main_block_size = 65536;
      if (unlikely(
              bigstack_alloc_u16(main_block_size, &(scan_ctx.bgen_allele_cts[0])) ||
              bigstack_alloc_u16(main_block_size, &(scan_ctx.bgen_allele_cts[1])) ||
              bigstack_alloc_ucp(main_block_size + 1, &(common.compressed_geno_starts[0])) ||
              bigstack_alloc_ucp(main_block_size + 1, &(common.compressed_geno_starts[1])))) {
        goto OxBgenToPgen_ret_NOMEM;
      }
      if (compression_mode) {
        if (unlikely(
                bigstack_alloc_u32(main_block_size, &(scan_ctx.uncompressed_genodata_byte_cts[0])) ||
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
      SetThreadFuncAndData(Bgen13DosageOrPhaseScanThread, &scan_ctx, &tg);

      uint32_t variant_ct = 0;

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
      uint32_t skip = 0;

      // temporary kludge
      uint32_t multiallelic_skip_ct = 0;

      for (uint32_t variant_uidx = 0; variant_uidx != raw_variant_ct; ) {
        // format is mostly identical to bgen 1.1; but there's no sample count,
        // and there is an allele count
        // logic is more similar to the second bgen 1.1 pass since we write the
        // .pvar here.
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
            logputs("\n");
            logerrputs("Error: Length-0 SNP ID in .bgen file.\n");
            goto OxBgenToPgen_ret_INCONSISTENT_INPUT;
          }
          if (unlikely(!fread_unlocked(loadbuf, snpid_slen, 1, bgenfile))) {
            goto OxBgenToPgen_ret_READ_FAIL;
          }
          loadbuf[snpid_slen] = '\0';
          rsid_start = R_CAST(char*, &(loadbuf[snpid_slen + 1]));
        }
        uint16_t rsid_slen;
        if (unlikely(!fread_unlocked(&rsid_slen, 2, 1, bgenfile))) {
          goto OxBgenToPgen_ret_READ_FAIL;
        }
        if (unlikely(!rsid_slen)) {
          logputs("\n");
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
              logputs("\n");
              logerrputs("Error: Length-0 chromosome ID in .bgen file.\n");
              goto OxBgenToPgen_ret_INCONSISTENT_INPUT;
            }
            if (unlikely(!fread_unlocked(chr_name_start, chr_name_slen, 1, bgenfile))) {
              goto OxBgenToPgen_ret_READ_FAIL;
            }
            if (strequal_k(chr_name_start, "NA", chr_name_slen)) {
              strcpy_k(chr_name_start, "0");
              chr_name_slen = 1;
            } else {
              chr_name_start[chr_name_slen] = '\0';
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
          reterr = GetOrAddChrCodeDestructive("--bgen file", 0, allow_extra_chrs, chr_name_start, &(chr_name_start[chr_name_slen]), cip, &cur_chr_code);
          if (unlikely(reterr)) {
            goto OxBgenToPgen_ret_1;
          }
          skip = !IsSet(cip->chr_mask, cur_chr_code);
        }

        uint32_t cur_bp;
        uint32_t cur_allele_ct = 0;
        if (unlikely(
                (!fread_unlocked(&cur_bp, 4, 1, bgenfile)) ||
                (!fread_unlocked(&cur_allele_ct, 2, 1, bgenfile)))) {
          goto OxBgenToPgen_ret_READ_FAIL;
        }
        if (unlikely(cur_allele_ct < 2)) {
          // this is undefined in the 1.3 standard; prohibit for now
          logputs("\n");
          logerrputs("Error: .bgen variant has fewer than two alleles.\n");
          goto OxBgenToPgen_ret_MALFORMED_INPUT;
        }
        ++variant_uidx;
        if (!(variant_uidx % 1000)) {
          printf("\r--bgen: %uk variants scanned.", variant_uidx / 1000);
          fflush(stdout);
        }

        // the "cur_allele_ct > 2" part is a temporary kludge
        if (skip || (cur_allele_ct > 2)) {
          if (!skip) {
            ++multiallelic_skip_ct;
          }
          for (uint32_t allele_idx = 0; allele_idx != cur_allele_ct; ++allele_idx) {
            uint32_t cur_allele_slen;
            if (unlikely(
                    (!fread_unlocked(&cur_allele_slen, 4, 1, bgenfile)) ||
                    fseeko(bgenfile, cur_allele_slen, SEEK_CUR))) {
              goto OxBgenToPgen_ret_READ_FAIL;
            }
          }
          uint32_t genodata_byte_ct;
          if (unlikely(
                  (!fread_unlocked(&genodata_byte_ct, 4, 1, bgenfile)) ||
                  fseeko(bgenfile, genodata_byte_ct, SEEK_CUR))) {
            goto OxBgenToPgen_ret_READ_FAIL;
          }
          continue;
        }
        if (unlikely(rsid_slen > kMaxIdSlen)) {
          // enforce this iff we aren't skipping
          logputs("\n");
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
          logputs("\n");
          logerrputs("Error: Empty allele code in .bgen file.\n");
          goto OxBgenToPgen_ret_MALFORMED_INPUT;
        }
        if (unlikely(a1_slen > 1000000000)) {
          logputs("\n");
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
          logputs("\n");
          logerrputs("Error: Empty allele code in .bgen file.\n");
          goto OxBgenToPgen_ret_MALFORMED_INPUT;
        }
        if (unlikely(a2_slen > 1000000000)) {
          logputs("\n");
          logerrputs("Error: Allele code in .bgen file has more than 1 billion characters.\n");
          goto OxBgenToPgen_ret_MALFORMED_INPUT;
        }
        if (unlikely(a2_slen + S_CAST(uintptr_t, a2_ptr - R_CAST(char*, loadbuf)) > loadbuf_size)) {
          goto OxBgenToPgen_ret_NOMEM;
        }
        if (unlikely(!fread_unlocked(a2_ptr, a2_slen, 1, bgenfile))) {
          goto OxBgenToPgen_ret_READ_FAIL;
        }
        write_iter = chrtoa(cip, cur_chr_code, write_iter);
        *write_iter++ = '\t';
        if (unlikely(cur_bp > 0x7ffffffe)) {
          logputs("\n");
          logerrputs("Error: Invalid bp coordinate (> 2^31 - 2) in .bgen file\n");
          goto OxBgenToPgen_ret_MALFORMED_INPUT;
        }
        write_iter = u32toa_x(cur_bp, '\t', write_iter);
        write_iter = memcpyax(write_iter, rsid_start, rsid_slen, '\t');
        if (prov_ref_allele_second) {
          const uint32_t swap_slen = a1_slen;
          a1_slen = a2_slen;
          a2_slen = swap_slen;
          char* swap_ptr = a1_ptr;
          a1_ptr = a2_ptr;
          a2_ptr = swap_ptr;
        }
        // allele codes may be too large for write buffer, so we special-case
        // this instead of using fwrite_ck()
        if ((write_iter >= writebuf_flush) || (a1_slen >= kMaxMediumLine)) {
          if (unlikely(fwrite_checked(writebuf, write_iter - writebuf, pvarfile))) {
            goto OxBgenToPgen_ret_WRITE_FAIL;
          }
          write_iter = writebuf;
        }
        if (a1_slen < kMaxMediumLine) {
          write_iter = memcpya(write_iter, a1_ptr, a1_slen);
        } else {
          if (unlikely(fwrite_checked(a1_ptr, a1_slen, pvarfile))) {
            goto OxBgenToPgen_ret_WRITE_FAIL;
          }
        }
        *write_iter++ = '\t';
        if ((write_iter >= writebuf_flush) || (a2_slen >= kMaxMediumLine)) {
          if (unlikely(fwrite_checked(writebuf, write_iter - writebuf, pvarfile))) {
            goto OxBgenToPgen_ret_WRITE_FAIL;
          }
          write_iter = writebuf;
        }
        if (a2_slen < kMaxMediumLine) {
          write_iter = memcpya(write_iter, a2_ptr, a2_slen);
        } else {
          if (unlikely(fwrite_checked(a2_ptr, a2_slen, pvarfile))) {
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
            logputs("\n");
            logerrputs("Error: Empty allele code in .bgen file.\n");
            goto OxBgenToPgen_ret_MALFORMED_INPUT;
          }
          if (unlikely(cur_allele_slen > 1000000000)) {
            logputs("\n");
            logerrputs("Error: Allele code in .bgen file has more than 1 billion characters.\n");
            goto OxBgenToPgen_ret_MALFORMED_INPUT;
          }
          if (unlikely(cur_allele_slen > loadbuf_size)) {
            goto OxBgenToPgen_ret_NOMEM;
          }
          if (unlikely(!fread_unlocked(loadbuf, cur_allele_slen, 1, bgenfile))) {
            goto OxBgenToPgen_ret_READ_FAIL;
          }
          *write_iter++ = ',';
          if ((write_iter >= writebuf_flush) || (cur_allele_slen >= kMaxMediumLine)) {
            if (unlikely(fwrite_checked(writebuf, write_iter - writebuf, pvarfile))) {
              goto OxBgenToPgen_ret_WRITE_FAIL;
            }
            write_iter = writebuf;
          }
          if (cur_allele_slen < kMaxMediumLine) {
            write_iter = memcpya(write_iter, loadbuf, cur_allele_slen);
          } else {
            if (unlikely(fwrite_checked(loadbuf, cur_allele_slen, pvarfile))) {
              goto OxBgenToPgen_ret_WRITE_FAIL;
            }
          }
        }

        AppendBinaryEoln(&write_iter);
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
        if (dosage_is_present) {
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
              reterr = S_CAST(PglErr, scan_ctx.err_info);
              if (unlikely(reterr)) {
                goto OxBgenToPgen_ret_bgen13_thread_fail;
              }
              dosage_is_present = scan_ctx.dosage_is_present;
              if (dosage_is_present) {
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
      if (multiallelic_skip_ct) {
        logputs("\n");
        logerrprintfww("Warning: %u multiallelic variant%s skipped (not yet supported).\n", multiallelic_skip_ct, (multiallelic_skip_ct == 1)? "" : "s");
      }
      if (unlikely(!variant_ct)) {
        logputs("\n");
        logerrprintf("Error: All %u variant%s in .bgen file skipped.\n", raw_variant_ct, (raw_variant_ct == 1)? "" : "s");
        goto OxBgenToPgen_ret_INCONSISTENT_INPUT;
      }
      if (variant_ct == block_vidx) {
        // with multiple threads, there's no guarantee that even the first
        // decompression job has launched (e.g. there's only 1 variant on the
        // relevant chromosome in the entire .bgen, and calc_thread_ct == 2).
        thread_bidxs[cur_thread_fill_idx + 1] = block_vidx;
        if (unlikely(SpawnThreads(&tg))) {
          goto OxBgenToPgen_ret_THREAD_CREATE_FAIL;
        }
        block_vidx = 0;
      }
      chr_filter_present = (variant_ct + multiallelic_skip_ct != raw_variant_ct);
      if (ThreadsAreActive(&tg)) {
        JoinThreads(&tg);
        reterr = S_CAST(PglErr, scan_ctx.err_info);
        if (unlikely(reterr)) {
          goto OxBgenToPgen_ret_bgen13_thread_fail;
        }
        if ((!block_vidx) || scan_ctx.dosage_is_present) {
          // ignore thread_bidxs[] in this case
          StopThreads(&tg);
        } else {
          for (; cur_thread_fill_idx != calc_thread_ct; ) {
            // save endpoint for current thread, and tell any leftover threads
            // to do nothing
            thread_bidxs[++cur_thread_fill_idx] = block_vidx;
          }
          DeclareLastThreadBlock(&tg);
          SpawnThreads(&tg);
          JoinThreads(&tg);
          reterr = S_CAST(PglErr, scan_ctx.err_info);
          if (unlikely(reterr)) {
            goto OxBgenToPgen_ret_bgen13_thread_fail;
          }
        }
        dosage_is_present = scan_ctx.dosage_is_present;
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
      reterr = SpgwInitPhase1(outname, allele_idx_offsets, nullptr, variant_ct, sample_ct, dosage_is_present? (kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent) : kfPgenGlobal0, (oxford_import_flags & (kfOxfordImportRefFirst | kfOxfordImportRefLast))? 1 : 2, &spgw, &spgw_alloc_cacheline_ct, &max_vrec_len);
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

      // const GparseFlags gparse_flags = dosage_is_present? (kfGparseHphase | kfGparseDosage | kfGparseDphase) : kfGparse0;
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
      if (unlikely(
              bigstack_alloc_uc(cachelines_avail * kCacheline, &(compressed_geno_bufs[0])) ||
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
            ctx.block_allele_idx_offsets[parity] = &(SpgwGetAlleleIdxOffsets(&spgw)[vidx_start]);
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
              if (chr_filter_present) {
                const uint32_t cur_chr_code2 = GetChrCode(R_CAST(char*, loadbuf), cip, chr_name_slen);

                // we scanned all the variants
                assert(!IsI32Neg(cur_chr_code2));

                skip = !IsSet(cip->chr_mask, cur_chr_code2);
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
            if (skip || (cur_allele_ct > 2)) {
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
            goto OxBgenToPgen_ret_bgen13_thread_fail;
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
          reterr = GparseFlush(ctx.gparse[parity], prev_block_write_ct, &spgw);
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
    if (unlikely(fclose_flush_null(writebuf_flush, write_iter, &pvarfile))) {
      goto OxBgenToPgen_ret_WRITE_FAIL;
    }

    SpgwFinish(&spgw);
    putc_unlocked('\r', stdout);
    write_iter = strcpya_k(g_logbuf, "--bgen: ");
    const uint32_t outname_base_slen = outname_end - outname;
    write_iter = memcpya(write_iter, outname, outname_base_slen + 5);
    write_iter = strcpya_k(write_iter, " + ");
    write_iter = memcpya(write_iter, outname, outname_base_slen);
    write_iter = strcpya_k(write_iter, ".pvar");
    write_iter = strcpya_k(write_iter, " written");
    if (!dosage_is_present) {
      write_iter = strcpya_k(write_iter, " (only unphased hardcalls)");
    }
    snprintf(write_iter, kLogbufSize - 2 * kPglFnamesize - 64, ".\n");
    WordWrapB(0);
    logputsb();
  }
  while (0) {
  OxBgenToPgen_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  OxBgenToPgen_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  OxBgenToPgen_ret_READ_FAIL:
    logerrprintfww(kErrprintfFread, bgenname, strerror(errno));
    reterr = kPglRetReadFail;
    break;
  OxBgenToPgen_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  OxBgenToPgen_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  OxBgenToPgen_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  OxBgenToPgen_ret_INCONSISTENT_INPUT_2:
    logerrputsb();
  OxBgenToPgen_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  OxBgenToPgen_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  OxBgenToPgen_ret_bgen13_thread_fail:
    if (reterr == kPglRetMalformedInput) {
    OxBgenToPgen_ret_bgen11_thread_fail:
      logputs("\n");
      logerrputs("Error: Invalid compressed SNP block in .bgen file.\n");
    } else if (reterr == kPglRetNotYetSupported) {
      logputs("\n");
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
  fclose_cond(pvarfile);
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

BoolErr ImportLegendCols(const char* fname, uintptr_t line_idx, uint32_t prov_ref_allele_second, const char** loadbuf_iter_ptr, char** write_iter_ptr, uint32_t* variant_ct_ptr) {
  {
    if (*variant_ct_ptr == 0x7ffffffd) {
      logerrputs("Error: " PROG_NAME_STR " does not support more than 2^31 - 3 variants.  We recommend using\nother software for very deep studies of small numbers of genomes.\n");
      return 1;
    }
    *variant_ct_ptr += 1;
    char* write_iter = *write_iter_ptr;
    *write_iter++ = '\t';
    const char* id_start = *loadbuf_iter_ptr;
    const char* id_end = CurTokenEnd(id_start);
    const uint32_t id_slen = id_end - id_start;
    if (id_slen > kMaxIdSlen) {
      logerrputs("Error: Variant names are limited to " MAX_ID_SLEN_STR " characters.\n");
      return 1;
    }
    const char* pos_str = FirstNonTspace(id_end);
    if (unlikely(!pos_str)) {
      goto ImportLegendCols_ret_MISSING_TOKENS;
    }
    const char* pos_end = CurTokenEnd(pos_str);
    uint32_t cur_bp;
    if (ScanUintDefcap(pos_str, &cur_bp)) {
      logpreprintfww("Error: Invalid bp coordinate on line %" PRIuPTR " of %s.\n", line_idx, fname);
      return 1;
    }
    write_iter = u32toa_x(cur_bp, '\t', write_iter);
    write_iter = memcpyax(write_iter, id_start, id_slen, '\t');
    const char* first_allele_str = FirstNonTspace(pos_end);
    if (unlikely(IsEolnKns(*first_allele_str))) {
      goto ImportLegendCols_ret_MISSING_TOKENS;
    }
    const char* first_allele_end = CurTokenEnd(first_allele_str);
    const char* second_allele_str = FirstNonTspace(first_allele_end);
    if (unlikely(IsEolnKns(*second_allele_str))) {
      goto ImportLegendCols_ret_MISSING_TOKENS;
    }
    const char* second_allele_end = CurTokenEnd(second_allele_str);
    if (!prov_ref_allele_second) {
      write_iter = memcpyax(write_iter, first_allele_str, first_allele_end - first_allele_str, '\t');
      write_iter = memcpya(write_iter, second_allele_str, second_allele_end - second_allele_str);
    } else {
      write_iter = memcpyax(write_iter, second_allele_str, second_allele_end - second_allele_str, '\t');
      write_iter = memcpya(write_iter, first_allele_str, first_allele_end - first_allele_str);
    }
    *write_iter_ptr = write_iter;
    AppendBinaryEoln(write_iter_ptr);
    *loadbuf_iter_ptr = second_allele_end;
    return 0;
  }
  {
  ImportLegendCols_ret_MISSING_TOKENS:
    logpreprintfww("Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, fname);
    return 1;
  }
}

PglErr ScanHapsForHet(const char* loadbuf_iter, const char* hapsname, uint32_t sample_ct, uint32_t is_haploid, uintptr_t line_idx_haps, uint32_t* at_least_one_het_ptr) {
  PglErr reterr = kPglRetSuccess;
  {
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      const uint32_t first_hap_char_code = ctou32(*loadbuf_iter);
      const uint32_t first_hap_int = first_hap_char_code - 48;
      // will .haps files ever support triallelic variants?  don't worry about
      // that for now
      const char* post_first_hap = &(loadbuf_iter[1]);
      if (unlikely((first_hap_int >= 2) || (ctou32(*post_first_hap) > 32))) {
        if (first_hap_char_code <= 32) {
          goto ScanHapsForHet_ret_MISSING_TOKENS;
        }
        goto ScanHapsForHet_ret_INVALID_TOKEN;
      }
      const char* second_hap = FirstNonTspace(post_first_hap);
      const char* post_second_hap = &(second_hap[1]);
      const uint32_t second_hap_char_code = ctou32(*second_hap);
      const uint32_t second_hap_int = second_hap_char_code - 48;
      const uint32_t post_second_hap_char_code = ctou32(*post_second_hap);
      if ((second_hap_int >= 2) || (post_second_hap_char_code > 32)) {
        // if haploid, permit '-' in second column
        if (unlikely((!is_haploid) || (second_hap_char_code != 45))) {
          if (second_hap_char_code <= 32) {
            goto ScanHapsForHet_ret_MISSING_TOKENS;
          }
          if ((second_hap_char_code == 45) && (post_second_hap_char_code <= 32)) {
            goto ScanHapsForHet_ret_HAPLOID_TOKEN;
          }
          goto ScanHapsForHet_ret_INVALID_TOKEN;
        }
      } else if (first_hap_int != second_hap_int) {
        *at_least_one_het_ptr = 1;
        break;
      }
      loadbuf_iter = FirstNonTspace(post_second_hap);
    }
  }
  while (0) {
  ScanHapsForHet_ret_HAPLOID_TOKEN:
    snprintf(g_logbuf, kLogbufSize, "Error: Haploid-only token on line %" PRIuPTR " of %s.\n", line_idx_haps, hapsname);
    reterr = kPglRetInconsistentInput;
    break;
  ScanHapsForHet_ret_INVALID_TOKEN:
    snprintf(g_logbuf, kLogbufSize, "Error: Invalid token on line %" PRIuPTR " of %s.\n", line_idx_haps, hapsname);
    reterr = kPglRetMalformedInput;
    break;
  ScanHapsForHet_ret_MISSING_TOKENS:
    snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx_haps, hapsname);
    reterr = kPglRetMalformedInput;
    break;
  }
  return reterr;
}

#ifdef __arm__
#  error "Unaligned accesses in OxHapslegendToPgen()."
#endif
PglErr OxHapslegendToPgen(const char* hapsname, const char* legendname, const char* samplename, const char* ox_single_chr_str, const char* ox_missing_code, MiscFlags misc_flags, ImportFlags import_flags, OxfordImportFlags oxford_import_flags, uint32_t max_thread_ct, char* outname, char* outname_end, ChrInfo* cip) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  uintptr_t line_idx_haps = 0;
  uintptr_t line_idx_legend = 0;
  PglErr reterr = kPglRetSuccess;
  TextStream haps_txs;
  TextStream legend_txs;
  STPgenWriter spgw;
  PreinitTextStream(&haps_txs);
  PreinitTextStream(&legend_txs);
  PreinitSpgw(&spgw);
  {
    uint32_t sfile_sample_ct = 0;
    if (samplename[0]) {
      reterr = OxSampleToPsam(samplename, ox_missing_code, import_flags, outname, outname_end, &sfile_sample_ct);
      if (unlikely(reterr)) {
        goto OxHapslegendToPgen_ret_1;
      }
      if (unlikely(sfile_sample_ct > (kMaxLongLine / 4))) {
        logerrputs("Error: Too many samples for .haps file converter.\n");
        reterr = kPglRetNotYetSupported;
        goto OxHapslegendToPgen_ret_1;
      }
    }

    uint32_t max_line_blen;
    if (StandardizeMaxLineBlen(bigstack_left() / 4, &max_line_blen)) {
      goto OxHapslegendToPgen_ret_NOMEM;
    }
    reterr = InitTextStream(hapsname, max_line_blen, ClipU32(max_thread_ct - 1, 1, 4), &haps_txs);
    if (unlikely(reterr)) {
      goto OxHapslegendToPgen_ret_TSTREAM_FAIL_HAPS;
    }
    uintptr_t writebuf_size = max_line_blen + kMaxMediumLine;
    char* writebuf = S_CAST(char*, bigstack_alloc_raw(writebuf_size));
    char* writebuf_flush = &(writebuf[kMaxMediumLine]);
    snprintf(outname_end, kMaxOutfnameExtBlen, ".pvar");
    if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
      goto OxHapslegendToPgen_ret_OPEN_FAIL;
    }
    char* write_iter = strcpya_k(writebuf, "#CHROM\tPOS\tID\tREF\tALT" EOLN_STR);
    ++line_idx_haps;
    char* haps_line_start = TextGet(&haps_txs);
    if (unlikely(!haps_line_start)) {
      if (!TextStreamErrcode2(&haps_txs, &reterr)) {
        snprintf(g_logbuf, kLogbufSize, "Error: %s is empty.\n", hapsname);
        goto OxHapslegendToPgen_ret_INCONSISTENT_INPUT_WW;
      }
      goto OxHapslegendToPgen_ret_TSTREAM_FAIL_HAPS;
    }
    const uint32_t token_ct = CountTokens(haps_line_start);
    // pass 1: count variants, write .pvar file, may as well also verify
    // there's at least one heterozygous call
    FinalizeChrset(misc_flags, cip);
    const uint32_t allow_extra_chrs = (misc_flags / kfMiscAllowExtraChrs) & 1;
    const uint32_t prov_ref_allele_second = !(oxford_import_flags & kfOxfordImportRefFirst);
    uint32_t at_least_one_het = 0;
    uintptr_t variant_skip_ct = 0;
    uint32_t variant_ct = 0;
    uint32_t is_haploid = 0;
    uint32_t sample_ct;
    // support both .haps + .legend (.haps expected to contain no header
    // columns), and pure .haps
    if (legendname[0]) {
      assert(ox_single_chr_str);
      if (unlikely(token_ct % 2)) {
        snprintf(g_logbuf, kLogbufSize, "Error: %s has an odd number of tokens in the first line. (With --haps + --legend, the .haps file is expected to have no header columns.)\n", hapsname);
        goto OxHapslegendToPgen_ret_MALFORMED_INPUT_WW;
      }
      sample_ct = token_ct / 2;
      if (unlikely(sfile_sample_ct && (sfile_sample_ct != sample_ct))) {
        snprintf(g_logbuf, kLogbufSize, "Error: .sample file has %u sample%s, while %s has %u.\n", sfile_sample_ct, (sfile_sample_ct == 1)? "" : "s", hapsname, sample_ct);
        goto OxHapslegendToPgen_ret_INCONSISTENT_INPUT_WW;
      }
      reterr = TextRewind(&haps_txs);
      if (unlikely(reterr)) {
        goto OxHapslegendToPgen_ret_TSTREAM_FAIL_HAPS;
      }
      line_idx_haps = 0;
      const uint32_t chr_code_raw = GetChrCodeRaw(ox_single_chr_str);
      const char* single_chr_str = nullptr;
      uint32_t single_chr_slen;
      char chr_buf[8];  // nothing longer than e.g. "chrMT" for now
      if (chr_code_raw == UINT32_MAX) {
        // command-line parser guarantees that allow_extra_chrs is true here
        single_chr_str = ox_single_chr_str;
        single_chr_slen = strlen(ox_single_chr_str);
      } else {
        uint32_t chr_code = chr_code_raw;
        if (chr_code > cip->max_code) {
          if (unlikely(chr_code < kMaxContigs)) {
            logerrputs("Error: --legend chromosome code is not in the chromosome set.\n");
            goto OxHapslegendToPgen_ret_INVALID_CMDLINE;
          }
          chr_code = cip->xymt_codes[chr_code - kMaxContigs];
          if (unlikely(IsI32Neg(chr_code))) {
            logerrputs("Error: --legend chromosome code is not in the chromosome set.\n");
            goto OxHapslegendToPgen_ret_INVALID_CMDLINE;
          }
        }
        if (unlikely(!IsSet(cip->chr_mask, chr_code))) {
          logerrputs("Error: --legend chromosome code is excluded by chromosome filter.\n");
          goto OxHapslegendToPgen_ret_INVALID_CMDLINE;
        }
        is_haploid = IsSet(cip->haploid_mask, chr_code);
        char* chr_name_end = chrtoa(cip, chr_code, chr_buf);
        single_chr_str = chr_buf;
        single_chr_slen = chr_name_end - chr_buf;
      }
      reterr = InitTextStream(legendname, max_line_blen, 1, &legend_txs);
      if (unlikely(reterr)) {
        goto OxHapslegendToPgen_ret_1;
      }
      ++line_idx_legend;
      char* legend_line_start = TextGet(&legend_txs);
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
      while (1) {
        ++line_idx_legend;
        legend_line_start = TextGet(&legend_txs);
        if (!legend_line_start) {
          if (likely(!TextStreamErrcode2(&legend_txs, &reterr))) {
            break;
          }
          goto OxHapslegendToPgen_ret_TSTREAM_FAIL_LEGEND;
        }
        const char* linebuf_iter = legend_line_start;
        write_iter = memcpya(write_iter, single_chr_str, single_chr_slen);
        if (unlikely(ImportLegendCols(legendname, line_idx_legend, prov_ref_allele_second, &linebuf_iter, &write_iter, &variant_ct))) {
          goto OxHapslegendToPgen_ret_MALFORMED_INPUT;
        }
        if (unlikely(fwrite_ck(writebuf_flush, outfile, &write_iter))) {
          goto OxHapslegendToPgen_ret_WRITE_FAIL;
        }
        if (!at_least_one_het) {
          ++line_idx_haps;
          haps_line_start = TextGet(&haps_txs);
          if (unlikely(!haps_line_start)) {
            if (!TextStreamErrcode2(&haps_txs, &reterr)) {
              if (line_idx_haps == 1) {
                goto OxHapslegendToPgen_ret_TSTREAM_REWIND_FAIL_HAPS;
              }
              snprintf(g_logbuf, kLogbufSize, "Error: %s has fewer nonheader lines than %s.\n", hapsname, legendname);
              goto OxHapslegendToPgen_ret_INCONSISTENT_INPUT_WW;
            }
            goto OxHapslegendToPgen_ret_TSTREAM_FAIL_HAPS;
          }
          reterr = ScanHapsForHet(haps_line_start, hapsname, sample_ct, is_haploid, line_idx_haps, &at_least_one_het);
          if (unlikely(reterr)) {
            WordWrapB(0);
            logerrputsb();
            goto OxHapslegendToPgen_ret_1;
          }
        }
      }
      BigstackReset(TextStreamMemStart(&legend_txs));
      if (CleanupTextStream2(legendname, &legend_txs, &reterr)) {
        goto OxHapslegendToPgen_ret_1;
      }
    } else {
      if (unlikely((token_ct < 7) || (!(token_ct % 2)))) {
        snprintf(g_logbuf, kLogbufSize, "Error: Unexpected token count in line %" PRIuPTR " of %s (should be odd, >5).\n", line_idx_haps, hapsname);
        goto OxHapslegendToPgen_ret_MALFORMED_INPUT_WW;
      }
      sample_ct = (token_ct - 5) / 2;
      if (unlikely(sfile_sample_ct && (sfile_sample_ct != sample_ct))) {
        snprintf(g_logbuf, kLogbufSize, "Error: .sample file has %u sample%s, while %s has %u.\n", sfile_sample_ct, (sfile_sample_ct == 1)? "" : "s", hapsname, sample_ct);
        goto OxHapslegendToPgen_ret_INCONSISTENT_INPUT_WW;
      }
      char* haps_line_iter = haps_line_start;
      while (1) {
        char* chr_code_end = CurTokenEnd(haps_line_iter);
        char* linebuf_iter = FirstNonTspace(chr_code_end);
        if (unlikely(IsEolnKns(*linebuf_iter))) {
          goto OxHapslegendToPgen_ret_MISSING_TOKENS_HAPS;
        }
        uint32_t cur_chr_code;
        reterr = GetOrAddChrCodeDestructive("--haps file", line_idx_haps, allow_extra_chrs, haps_line_iter, chr_code_end, cip, &cur_chr_code);
        if (unlikely(reterr)) {
          goto OxHapslegendToPgen_ret_1;
        }
        if (!IsSet(cip->chr_mask, cur_chr_code)) {
          ++variant_skip_ct;
        } else {
          is_haploid = IsSet(cip->haploid_mask, cur_chr_code);
          write_iter = chrtoa(cip, cur_chr_code, write_iter);
          if (unlikely(ImportLegendCols(hapsname, line_idx_haps, prov_ref_allele_second, K_CAST(const char**, &linebuf_iter), &write_iter, &variant_ct))) {
            putc_unlocked('\n', stdout);
            goto OxHapslegendToPgen_ret_MALFORMED_INPUT;
          }
          if (unlikely(fwrite_ck(writebuf_flush, outfile, &write_iter))) {
            goto OxHapslegendToPgen_ret_WRITE_FAIL;
          }
          if (!at_least_one_het) {
            linebuf_iter = FirstNonTspace(linebuf_iter);
            if (unlikely(ScanHapsForHet(linebuf_iter, hapsname, sample_ct, is_haploid, line_idx_haps, &at_least_one_het))) {
              // override InconsistentInput return code since chromosome info
              // was also gathered from .haps file
              goto OxHapslegendToPgen_ret_MALFORMED_INPUT_WW;
            }
          }
        }
        haps_line_iter = AdvPastDelim(linebuf_iter, '\n');
        ++line_idx_haps;
        if (!TextGetUnsafe2(&haps_txs, &haps_line_iter)) {
          if (likely(!TextStreamErrcode2(&haps_txs, &reterr))) {
            break;
          }
          goto OxHapslegendToPgen_ret_TSTREAM_FAIL_HAPS;
        }
      }
      if (unlikely(!variant_ct)) {
        snprintf(g_logbuf, kLogbufSize, "Error: All %" PRIuPTR " variant%s in %s skipped due to chromosome filter.\n", variant_skip_ct, (variant_skip_ct == 1)? "" : "s", hapsname);
        goto OxHapslegendToPgen_ret_INCONSISTENT_INPUT_WW;
      }
    }
    if (unlikely(fclose_flush_null(writebuf_flush, write_iter, &outfile))) {
      goto OxHapslegendToPgen_ret_WRITE_FAIL;
    }
    if (!sfile_sample_ct) {
      // create a dummy .psam file with "per0", "per1", etc. IDs, matching
      // --dummy
      snprintf(outname_end, kMaxOutfnameExtBlen, ".psam");
      if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
        goto OxHapslegendToPgen_ret_OPEN_FAIL;
      }
      write_iter = strcpya_k(writebuf, "#IID\tSEX" EOLN_STR);
      for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        write_iter = strcpya_k(write_iter, "per");
        write_iter = u32toa(sample_idx, write_iter);
        write_iter = strcpya_k(write_iter, "\tNA" EOLN_STR);
        if (unlikely(fwrite_ck(writebuf_flush, outfile, &write_iter))) {
          goto OxHapslegendToPgen_ret_WRITE_FAIL;
        }
      }
      if (unlikely(fclose_flush_null(writebuf_flush, write_iter, &outfile))) {
        goto OxHapslegendToPgen_ret_WRITE_FAIL;
      }
    }
    reterr = TextRewind(&haps_txs);
    if (unlikely(reterr)) {
      goto OxHapslegendToPgen_ret_TSTREAM_FAIL_HAPS;
    }
    line_idx_haps = 0;
    BigstackReset(writebuf);
    logprintf("--haps%s: %u variant%s scanned.\n", legendname[0]? " + --legend" : "", variant_ct, (variant_ct == 1)? "" : "s");
    snprintf(outname_end, kMaxOutfnameExtBlen, ".pgen");
    uintptr_t spgw_alloc_cacheline_ct;
    uint32_t max_vrec_len;
    reterr = SpgwInitPhase1(outname, nullptr, nullptr, variant_ct, sample_ct, at_least_one_het? kfPgenGlobalHardcallPhasePresent : kfPgenGlobal0, (oxford_import_flags & (kfOxfordImportRefFirst | kfOxfordImportRefLast))? 1 : 2, &spgw, &spgw_alloc_cacheline_ct, &max_vrec_len);
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
    const uint32_t phaseinfo_match_4char = prov_ref_allele_second? 0x20312030 : 0x20302031;
    const uint32_t phaseinfo_match = 1 + prov_ref_allele_second;
    uintptr_t* genovec;
    uintptr_t* phaseinfo;
    if (unlikely(
            bigstack_alloc_w(sample_ctl2, &genovec) ||
            bigstack_alloc_w(sample_ctl, &phaseinfo))) {
      goto OxHapslegendToPgen_ret_NOMEM;
    }

    char* haps_line_iter = TextLineEnd(&haps_txs);
    for (uint32_t vidx = 0; vidx != variant_ct; ) {
      ++line_idx_haps;
      reterr = TextGetUnsafe(&haps_txs, &haps_line_iter);
      if (unlikely(reterr)) {
        if ((reterr == kPglRetEof) && (line_idx_haps > 1) && (legendname[0])) {
          snprintf(g_logbuf, kLogbufSize, "Error: %s has fewer nonheader lines than %s.\n", hapsname, legendname);
          goto OxHapslegendToPgen_ret_INCONSISTENT_INPUT_WW;
        }
        goto OxHapslegendToPgen_ret_TSTREAM_REWIND_FAIL_HAPS;
      }
      char* linebuf_iter = haps_line_iter;
      if (!legendname[0]) {
        char* chr_code_end = CurTokenEnd(haps_line_iter);
        const uint32_t cur_chr_code = GetChrCodeCounted(cip, chr_code_end - haps_line_iter, haps_line_iter);
        if (!IsSet(cip->chr_mask, cur_chr_code)) {
          haps_line_iter = AdvPastDelim(haps_line_iter, '\n');
          continue;
        }
        is_haploid = IsSet(cip->haploid_mask, cur_chr_code);
        linebuf_iter = NextTokenMult(FirstNonTspace(chr_code_end), 4);
        if (unlikely(!linebuf_iter)) {
          goto OxHapslegendToPgen_ret_MISSING_TOKENS_HAPS;
        }
      }
      uintptr_t genovec_word_or = 0;
      uint32_t inner_loop_last = kBitsPerWordD2 - 1;
      // optimize common case: autosomal diploid, always exactly one space
      // this loop is time-critical; all my attemps to merge in the haploid
      // case have caused >10% slowdowns
      if ((!is_haploid) && (ctou32(linebuf_iter[sample_ct * 4 - 1]) < 32)) {
        haps_line_iter = AdvPastDelim(&(linebuf_iter[sample_ct * 4 - 1]), '\n');
        linebuf_iter[sample_ct * 4 - 1] = ' ';
#ifdef __LP64__
        const VecU16* linebuf_viter = R_CAST(const VecU16*, linebuf_iter);
        const VecU16 all0 = vecu16_set1(0x2030);
        const VecU16 all1 = vecu16_set1(0x2031);
        const uint32_t fullword_ct = sample_ct / kBitsPerWordD2;
        Halfword* phaseinfo_alias = R_CAST(Halfword*, phaseinfo);
        for (uint32_t widx = 0; widx != fullword_ct; ++widx) {
          uintptr_t geno_first = 0;
          for (uint32_t uii = 0; uii != 2; ++uii) {
            VecU16 cur_chars = vecu16_loadu(linebuf_viter);
            ++linebuf_viter;
            uintptr_t zero_mm = vecu16_movemask(cur_chars == all0);
            uintptr_t one_mm = vecu16_movemask(cur_chars == all1);
            cur_chars = vecu16_loadu(linebuf_viter);
            ++linebuf_viter;
            zero_mm |= S_CAST(uintptr_t, vecu16_movemask(cur_chars == all0)) << kBytesPerVec;
            one_mm |= S_CAST(uintptr_t, vecu16_movemask(cur_chars == all1)) << kBytesPerVec;
#  ifndef USE_AVX2
            cur_chars = vecu16_loadu(linebuf_viter);
            ++linebuf_viter;
            zero_mm |= S_CAST(uintptr_t, vecu16_movemask(cur_chars == all0)) << 32;
            one_mm |= S_CAST(uintptr_t, vecu16_movemask(cur_chars == all1)) << 32;
            cur_chars = vecu16_loadu(linebuf_viter);
            ++linebuf_viter;
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
          phaseinfo_alias[widx] = phaseinfo_hw;
        }
        const uint32_t remainder = sample_ct % kBitsPerWordD2;
        if (remainder) {
          const uint32_t* linebuf_alias32_iter = R_CAST(const uint32_t*, linebuf_viter);
          uintptr_t genovec_word = 0;
          Halfword phaseinfo_hw = 0;
          for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != remainder; ++sample_idx_lowbits) {
            uint32_t cur_hap_4char = *linebuf_alias32_iter++;
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
          phaseinfo_alias[fullword_ct] = phaseinfo_hw;
        }
#else  // !__LP64__
        const uint32_t* linebuf_alias32_iter = R_CAST(const uint32_t*, linebuf_iter);
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
            uint32_t cur_hap_4char = *linebuf_alias32_iter++;
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
          R_CAST(Halfword*, phaseinfo)[widx] = phaseinfo_hw;
        }
#endif  // !__LP64__
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
          R_CAST(Halfword*, phaseinfo)[widx] = phaseinfo_hw;
        }
        haps_line_iter = AdvPastDelim(linebuf_iter, '\n');
      }
      if (prov_ref_allele_second) {
        GenovecInvertUnsafe(sample_ct, genovec);
        ZeroTrailingNyps(sample_ct, genovec);
      }
      if (genovec_word_or & kMask5555) {
        if (unlikely(SpgwAppendBiallelicGenovecHphase(genovec, nullptr, phaseinfo, &spgw))) {
          goto OxHapslegendToPgen_ret_WRITE_FAIL;
        }
      } else {
        if (unlikely(SpgwAppendBiallelicGenovec(genovec, &spgw))) {
          goto OxHapslegendToPgen_ret_WRITE_FAIL;
        }
      }
      if (!(++vidx % 1000)) {
        printf("\r--haps%s: %uk variants converted.", legendname[0]? " + --legend" : "", vidx / 1000);
        fflush(stdout);
      }
    }
    SpgwFinish(&spgw);
    putc_unlocked('\r', stdout);
    write_iter = strcpya_k(g_logbuf, "--haps");
    if (legendname[0]) {
      write_iter = strcpya_k(write_iter, " + --legend");
    }
    write_iter = strcpya_k(write_iter, ": ");
    const uint32_t outname_base_slen = outname_end - outname;
    write_iter = memcpya(write_iter, outname, outname_base_slen + 5);
    write_iter = strcpya_k(write_iter, " + ");
    if (!sfile_sample_ct) {
      write_iter = memcpya(write_iter, outname, outname_base_slen);
      write_iter = strcpya_k(write_iter, ".psam + ");
    }
    write_iter = memcpya(write_iter, outname, outname_base_slen);
    snprintf(write_iter, kLogbufSize - 3 * kPglFnamesize - 64, ".pvar written.\n");
    WordWrapB(0);
    logputsb();
  }
  while (0) {
  OxHapslegendToPgen_ret_TSTREAM_REWIND_FAIL_HAPS:
    putc_unlocked('\n', stdout);
    TextStreamErrPrintRewind(hapsname, &haps_txs, &reterr);
    break;
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
    reterr = kPglRetOpenFail;
    break;
  OxHapslegendToPgen_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
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
    reterr = kPglRetInconsistentInput;
    break;
  }
 OxHapslegendToPgen_ret_1:
  CleanupSpgw(&spgw, &reterr);
  CleanupTextStream2(legendname, &legend_txs, &reterr);
  CleanupTextStream2(hapsname, &haps_txs, &reterr);
  fclose_cond(outfile);
  BigstackReset(bigstack_mark);
  return reterr;
}


// could add an option to LoadPvar() to not require allele columns, but .map
// is easy enough to write a separate loader for...
CONSTI32(kLoadMapBlockSize, 65536);

// assumes FinalizeChrset() has already been called.
static_assert(kMaxContigs <= 65536, "LoadMap() needs to be updated.");
PglErr LoadMap(const char* mapname, MiscFlags misc_flags, ChrInfo* cip, uint32_t* max_variant_id_slen_ptr, uint16_t** variant_chr_codes_ptr, uint32_t** variant_bps_ptr, char*** variant_ids_ptr, double** variant_cms_ptr, uint32_t* variant_ct_ptr) {
  // caller should call ForgetExtraChrNames(1, cip) after finishing import.
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  uintptr_t line_idx = 0;
  PglErr reterr = kPglRetSuccess;
  TextStream map_txs;
  PreinitTextStream(&map_txs);
  {
    // Workspace used as follows:
    // |--linebuf--|--temp-->----|----<- variant IDs --|
    //            1/4                                 end
    // linebuf is overwritten with the main return arrays at the end.
    reterr = SizeAndInitTextStream(mapname, bigstack_left() / 4, 1, &map_txs);
    if (unlikely(reterr)) {
      goto LoadMap_ret_TSTREAM_FAIL;
    }
    char* line_start;
    do {
      ++line_idx;
      line_start = TextGet(&map_txs);
      if (unlikely(!line_start)) {
        if (!TextStreamErrcode2(&map_txs, &reterr)) {
          logerrputs("Error: Empty .map file.\n");
          goto LoadMap_ret_INCONSISTENT_INPUT;
        }
        goto LoadMap_ret_TSTREAM_FAIL;
      }
    } while (*line_start == '#');
    uint32_t map_cols = 3;
    {
      const char* linebuf_iter = NextTokenMult(line_start, 2);
      if (unlikely(!linebuf_iter)) {
        goto LoadMap_ret_MISSING_TOKENS;
      }
      linebuf_iter = NextToken(linebuf_iter);
      if (linebuf_iter) {
        linebuf_iter = NextToken(linebuf_iter);
        if (!linebuf_iter) {
          map_cols = 4;
        } else {
          linebuf_iter = NextToken(linebuf_iter);
          if (linebuf_iter) {
            if (unlikely(NextToken(linebuf_iter))) {
              // do NOT permit >6 columns, .bim is ok but .pvar is not
              // (pointless to support .pvar for legacy formats)
              snprintf(g_logbuf, kLogbufSize, "Error: %s is not a .map/.bim file (too many columns).\n", mapname);
              goto LoadMap_ret_MALFORMED_INPUT_WW;
            }
            map_cols = 4;
          }
        }
      }
    }

    const uint32_t allow_extra_chrs = (misc_flags / kfMiscAllowExtraChrs) & 1;
    uint32_t max_variant_id_slen = *max_variant_id_slen_ptr;
    unsigned char* tmp_alloc_base = g_bigstack_base;
    unsigned char* tmp_alloc_end = bigstack_end_mark;
    uint16_t* cur_chr_codes = nullptr;
    uint32_t* cur_bps = nullptr;
    char** cur_ids = nullptr;
    double* cur_cms = nullptr;
    double cur_cm = 0.0;
    uint32_t at_least_one_nzero_cm = 0;
    uint32_t variant_ct = 0;
    char* line_iter = line_start;
    while (1) {
      {
        // chrom, id, (cm?), pos
        char* chr_code_end = CurTokenEnd(line_iter);
        char* linebuf_iter = FirstNonTspace(chr_code_end);
        if (unlikely(IsEolnKns(*linebuf_iter))) {
          goto LoadMap_ret_MISSING_TOKENS;
        }
        uint32_t cur_chr_code;
        reterr = GetOrAddChrCodeDestructive(".map file", line_idx, allow_extra_chrs, line_iter, chr_code_end, cip, &cur_chr_code);
        if (unlikely(reterr)) {
          goto LoadMap_ret_1;
        }
        if (!IsSet(cip->chr_mask, cur_chr_code)) {
          line_iter = linebuf_iter;
          goto LoadMap_skip_variant;
        }
        char* token_end = CurTokenEnd(linebuf_iter);
        uint32_t id_slen = token_end - linebuf_iter;
        if (id_slen > max_variant_id_slen) {
          max_variant_id_slen = id_slen;
        }
        tmp_alloc_end -= id_slen + 1;
        if (unlikely(tmp_alloc_end < tmp_alloc_base)) {
          goto LoadMap_ret_NOMEM;
        }
        memcpyx(tmp_alloc_end, linebuf_iter, id_slen, '\0');
        linebuf_iter = FirstNonTspace(token_end);
        if (unlikely(IsEolnKns(*linebuf_iter))) {
          goto LoadMap_ret_MISSING_TOKENS;
        }

        if (map_cols == 4) {
          char* cm_end = ScantokDouble(linebuf_iter, &cur_cm);
          if (unlikely(!cm_end)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid centimorgan position on line %" PRIuPTR " of %s.\n", line_idx, mapname);
            goto LoadMap_ret_MALFORMED_INPUT_WW;
          }
          at_least_one_nzero_cm = (cur_cm != 0.0);
          linebuf_iter = NextToken(cm_end);
          if (unlikely(!linebuf_iter)) {
            goto LoadMap_ret_MISSING_TOKENS;
          }
        }
        int32_t cur_bp;
        if (unlikely(ScanIntAbsDefcap(linebuf_iter, &cur_bp))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Invalid bp coordinate on line %" PRIuPTR " of %s.\n", line_idx, mapname);
          goto LoadMap_ret_MALFORMED_INPUT_WW;
        }
        line_iter = linebuf_iter;
        if (cur_bp < 0) {
          goto LoadMap_skip_variant;
        }

        const uint32_t variant_idx_lowbits = variant_ct % kLoadMapBlockSize;
        if (!variant_idx_lowbits) {
          if (unlikely(S_CAST(uintptr_t, tmp_alloc_end - tmp_alloc_base) <= kLoadMapBlockSize * (sizeof(int16_t) + sizeof(int32_t) + sizeof(intptr_t) + sizeof(double)))) {
            goto LoadMap_ret_NOMEM;
          }
          cur_chr_codes = R_CAST(uint16_t*, tmp_alloc_base);
          tmp_alloc_base = R_CAST(unsigned char*, &(cur_chr_codes[kLoadMapBlockSize]));
          cur_bps = R_CAST(uint32_t*, tmp_alloc_base);
          tmp_alloc_base = R_CAST(unsigned char*, &(cur_bps[kLoadMapBlockSize]));
          cur_ids = R_CAST(char**, tmp_alloc_base);
          tmp_alloc_base = R_CAST(unsigned char*, &(cur_ids[kLoadMapBlockSize]));
          cur_cms = R_CAST(double*, tmp_alloc_base);
          tmp_alloc_base = R_CAST(unsigned char*, &(cur_cms[kLoadMapBlockSize]));
        }
        cur_chr_codes[variant_idx_lowbits] = cur_chr_code;
        cur_ids[variant_idx_lowbits] = R_CAST(char*, tmp_alloc_end);
        cur_cms[variant_idx_lowbits] = cur_cm;
        cur_bps[variant_idx_lowbits] = cur_bp;
        ++variant_ct;
      }
    LoadMap_skip_variant:
      ++line_idx;
      line_iter = AdvPastDelim(line_iter, '\n');
      if (!TextGetUnsafe2(&map_txs, &line_iter)) {
        if (!TextStreamErrcode2(&map_txs, &reterr)) {
          break;
        }
        goto LoadMap_ret_TSTREAM_FAIL;
      }
      if (unlikely(line_iter[0] == '#')) {
        snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of %s starts with a '#'. (This is only permitted before the first nonheader line.)\n", line_idx, mapname);
        goto LoadMap_ret_MALFORMED_INPUT_WW;
      }
    }
    if (unlikely(max_variant_id_slen > kMaxIdSlen)) {
      logerrputs("Error: Variant names are limited to " MAX_ID_SLEN_STR " characters.\n");
      goto LoadMap_ret_MALFORMED_INPUT;
    }

    if (unlikely(!variant_ct)) {
      logerrputs("Error: All variants in .map/.bim file skipped due to chromosome filter.\n");
      goto LoadMap_ret_INCONSISTENT_INPUT;
    }
    // true requirement is weaker, but whatever
    g_bigstack_end = g_bigstack_base;
    g_bigstack_base = TextStreamMemStart(&map_txs);
    if (CleanupTextStream2(mapname, &map_txs, &reterr)) {
      goto LoadMap_ret_1;
    }

    if (unlikely(
            bigstack_alloc_u16(variant_ct, variant_chr_codes_ptr) ||
            bigstack_alloc_u32(variant_ct, variant_bps_ptr) ||
            bigstack_alloc_cp(variant_ct, variant_ids_ptr))) {
      goto LoadMap_ret_NOMEM;
    }
    uint16_t* variant_chr_codes = *variant_chr_codes_ptr;
    uint32_t* variant_bps = *variant_bps_ptr;
    char** variant_ids = *variant_ids_ptr;
    double* variant_cms = nullptr;
    if (at_least_one_nzero_cm) {
      if (unlikely(bigstack_alloc_d(variant_ct, variant_cms_ptr))) {
        goto LoadMap_ret_NOMEM;
      }
      variant_cms = *variant_cms_ptr;
    } else {
      *variant_cms_ptr = nullptr;
    }
    *max_variant_id_slen_ptr = max_variant_id_slen;
    *variant_ct_ptr = variant_ct;
    const uint32_t full_block_ct = variant_ct / kLoadMapBlockSize;
    bigstack_mark = g_bigstack_base;
    unsigned char* read_iter = g_bigstack_end;  // bugfix (30 May 2018)
    BigstackEndSet(tmp_alloc_end);
    bigstack_end_mark = g_bigstack_end;

    for (uint32_t block_idx = 0; block_idx != full_block_ct; ++block_idx) {
      memcpy(&(variant_chr_codes[block_idx * kLoadMapBlockSize]), read_iter, kLoadMapBlockSize * sizeof(int16_t));
      read_iter = &(read_iter[kLoadMapBlockSize * sizeof(int16_t)]);
      memcpy(&(variant_bps[block_idx * kLoadMapBlockSize]), read_iter, kLoadMapBlockSize * sizeof(int32_t));
      read_iter = &(read_iter[kLoadMapBlockSize * sizeof(int32_t)]);
      memcpy(&(variant_ids[block_idx * kLoadMapBlockSize]), read_iter, kLoadMapBlockSize * sizeof(intptr_t));
      read_iter = &(read_iter[kLoadMapBlockSize * sizeof(intptr_t)]);
      if (at_least_one_nzero_cm) {
        memcpy(&(variant_cms[block_idx * kLoadMapBlockSize]), read_iter, kLoadMapBlockSize * sizeof(double));
      }
      read_iter = &(read_iter[kLoadMapBlockSize * sizeof(double)]);
    }
    const uint32_t variant_ct_lowbits = variant_ct % kLoadMapBlockSize;
    memcpy(&(variant_chr_codes[full_block_ct * kLoadMapBlockSize]), read_iter, variant_ct_lowbits * sizeof(int16_t));
    read_iter = &(read_iter[kLoadMapBlockSize * sizeof(int16_t)]);
    memcpy(&(variant_bps[full_block_ct * kLoadMapBlockSize]), read_iter, variant_ct_lowbits * sizeof(int32_t));
    read_iter = &(read_iter[kLoadMapBlockSize * sizeof(int32_t)]);
    memcpy(&(variant_ids[full_block_ct * kLoadMapBlockSize]), read_iter, variant_ct_lowbits * sizeof(intptr_t));
    if (at_least_one_nzero_cm) {
      read_iter = &(read_iter[kLoadMapBlockSize * sizeof(intptr_t)]);
      memcpy(&(variant_cms[full_block_ct * kLoadMapBlockSize]), read_iter, variant_ct_lowbits * sizeof(double));
    }
  }
  while (0) {
  LoadMap_ret_TSTREAM_FAIL:
    TextStreamErrPrint(mapname, &map_txs);
    break;
  LoadMap_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  LoadMap_ret_MISSING_TOKENS:
    snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, mapname);
  LoadMap_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  LoadMap_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  LoadMap_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 LoadMap_ret_1:
  // ForgetExtraChrNames(1, cip);
  CleanupTextStream2(mapname, &map_txs, &reterr);
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

// probable todo: use the VCF/BGEN-import parallelization strategy to speed
// this up.
static_assert(sizeof(Dosage) == 2, "Plink1DosageToPgen() needs to be updated.");
PglErr Plink1DosageToPgen(const char* dosagename, const char* famname, const char* mapname, const char* import_single_chr_str, const Plink1DosageInfo* pdip, MiscFlags misc_flags, ImportFlags import_flags, FamCol fam_cols, int32_t missing_pheno, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, double import_dosage_certainty, uint32_t max_thread_ct, char* outname, char* outname_end, ChrInfo* cip) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;

  // these are not allocated on bigstack, and must be explicitly freed
  PhenoCol* pheno_cols = nullptr;
  char* pheno_names = nullptr;
  uint32_t pheno_ct = 0;

  FILE* outfile = nullptr;
  uintptr_t line_idx = 0;
  Plink1DosageFlags flags = pdip->flags;
  PglErr reterr = kPglRetSuccess;
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
    reterr = LoadPsam(famname, nullptr, fam_cols, 0x7fffffff, missing_pheno, (misc_flags / kfMiscAffection01) & 1, max_thread_ct, &pii, &sample_include, &founder_info, &sex_nm, &sex_male, &pheno_cols, &pheno_names, &raw_sample_ct, &pheno_ct, &max_pheno_name_blen);
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
      if (unlikely(
              bigstack_end_alloc_u32(tmp_htable_size, &htable_tmp) ||
              bigstack_end_alloc_c(pii.sii.max_sample_id_blen, &idbuf))) {
        goto Plink1DosageToPgen_ret_NOMEM;
      }
      const uint32_t duplicate_idx = PopulateStrboxHtable(pii.sii.sample_ids, raw_sample_ct, pii.sii.max_sample_id_blen, tmp_htable_size, htable_tmp);
      if (unlikely(duplicate_idx)) {
        char* duplicate_sample_id = &(pii.sii.sample_ids[duplicate_idx * pii.sii.max_sample_id_blen]);
        char* duplicate_fid_end = AdvToDelim(duplicate_sample_id, '\t');
        *duplicate_fid_end = ' ';
        snprintf(g_logbuf, kLogbufSize, "Error: Duplicate sample ID '%s' in .fam file.\n", duplicate_sample_id);
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
          snprintf(g_logbuf, kLogbufSize, "Error: Duplicate sample ID '%s' in dosage file.\n", idbuf);
          goto Plink1DosageToPgen_ret_MALFORMED_INPUT_WW;
        }
        SetBit(sample_uidx, sample_include);
        dosage_sample_idx_to_fam_uidx[sample_ct++] = sample_uidx;
        loadbuf_iter = FirstNonTspace(iid_end);
      } while (!IsEolnKns(*loadbuf_iter));
      // Not worth the trouble of move-constructing dosage_txs from dosage_txf,
      // since we may need to load the .map in between, and we'll need to
      // rewind again anyway.
      if (CleanupTextFile2(dosagename, &dosage_txf, &reterr)) {
        goto Plink1DosageToPgen_ret_1;
      }
      line_idx = 0;
    }

    snprintf(outname_end, kMaxOutfnameExtBlen, ".psam");
    if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
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
    const uint32_t write_parents = DataParentalColsAreRequired(sample_include, &pii, sample_ct, 1);
    if (write_parents) {
      write_iter = strcpya_k(write_iter, "\tPAT\tMAT");
    }
    write_iter = strcpya_k(write_iter, "\tSEX");
    for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
      *write_iter++ = '\t';
      write_iter = strcpya(write_iter, &(pheno_names[pheno_idx * max_pheno_name_blen]));
      if (unlikely(fwrite_ck(writebuf_flush, outfile, &write_iter))) {
        goto Plink1DosageToPgen_ret_WRITE_FAIL;
      }
    }
    AppendBinaryEoln(&write_iter);
    uint32_t omp_slen = 2;
    char output_missing_pheno[kMaxMissingPhenostrBlen];
    if (import_flags & kfImportKeepAutoconv) {
      omp_slen = strlen(g_output_missing_pheno);
      memcpy(output_missing_pheno, g_output_missing_pheno, omp_slen);
    } else {
      memcpy_k(output_missing_pheno, "NA", 2);
    }
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
        if (unlikely(fwrite_ck(writebuf_flush, outfile, &write_iter))) {
          goto Plink1DosageToPgen_ret_WRITE_FAIL;
        }
        *write_iter++ = '\t';
        write_iter = AppendPhenoStr(&(pheno_cols[pheno_idx]), output_missing_pheno, omp_slen, sample_uidx, write_iter);
      }
      AppendBinaryEoln(&write_iter);
      if (unlikely(fwrite_ck(writebuf_flush, outfile, &write_iter))) {
        goto Plink1DosageToPgen_ret_WRITE_FAIL;
      }
    }
    if (unlikely(fclose_flush_null(writebuf_flush, write_iter, &outfile))) {
      goto Plink1DosageToPgen_ret_WRITE_FAIL;
    }
    // Don't need sample info any more.
    BigstackEndReset(bigstack_end_mark);

    // 3. Read .map file if it exists.
    uint32_t max_variant_id_slen = 1;
    uint16_t* variant_chr_codes = nullptr;
    uint32_t* variant_bps = nullptr;
    char** variant_ids = nullptr;
    double* variant_cms = nullptr;
    uint32_t* variant_id_htable = nullptr;
    uintptr_t* variant_already_seen = nullptr;
    uint32_t variant_id_htable_size = 0;
    uint32_t map_variant_ct = 0;
    FinalizeChrset(misc_flags, cip);
    if (mapname) {
      reterr = LoadMap(mapname, misc_flags, cip, &max_variant_id_slen, &variant_chr_codes, &variant_bps, &variant_ids, &variant_cms, &map_variant_ct);
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

    uintptr_t ulii = bigstack_left() / 2;
    // writebuf needs kMaxMediumLine more bytes, and we also need to allocate a
    // chromosome buffer.
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
    if (unlikely(bigstack_alloc_c(kMaxMediumLine + kDecompressChunkSize + max_line_blen, &writebuf))) {
      // shouldn't be possible
      goto Plink1DosageToPgen_ret_NOMEM;
    }
    writebuf_flush = &(writebuf[kMaxMediumLine]);
    const uint32_t allow_extra_chrs = (misc_flags / kfMiscAllowExtraChrs) & 1;
    const char* single_chr_str = nullptr;
    uint32_t single_chr_slen = 0;
    const uint32_t chr_col_idx = pdip->chr_col_idx;
    const uint32_t check_chr_col = (chr_col_idx != UINT32_MAX);
    if (!check_chr_col) {
      if (import_single_chr_str) {
        uint32_t chr_code_raw = GetChrCodeRaw(import_single_chr_str);
        if (chr_code_raw == UINT32_MAX) {
          // command-line parser guarantees that allow_extra_chrs is true here
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
    snprintf(outname_end, kMaxOutfnameExtBlen, ".pvar");
    if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
      goto Plink1DosageToPgen_ret_OPEN_FAIL;
    }
    write_iter = strcpya_k(writebuf, "#CHROM\tPOS\tID\tREF\tALT");
    if (variant_cms) {
      write_iter = strcpya_k(write_iter, "\tCM");
    }
    AppendBinaryEoln(&write_iter);
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

    STD_SORT(relevant_initial_col_ct, u64cmp, parse_table);
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
    uint32_t dosage_is_present = 0;
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
        write_iter = chrtoa(cip, variant_chr_codes[variant_uidx], write_iter);
        *write_iter++ = '\t';
        write_iter = u32toa_x(variant_bps[variant_uidx], '\t', write_iter);
        write_iter = memcpya(write_iter, variant_id, variant_id_slen);
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
          reterr = GetOrAddChrCodeDestructive("--import-dosage file", line_idx, allow_extra_chrs, chr_code_str, chr_code_end, cip, &cur_chr_code);
          if (unlikely(reterr)) {
            goto Plink1DosageToPgen_ret_1;
          }
          if (!IsSet(cip->chr_mask, cur_chr_code)) {
            ++variant_skip_ct;
            continue;
          }
          write_iter = chrtoa(cip, cur_chr_code, write_iter);
        } else {
          write_iter = memcpya(write_iter, single_chr_str, single_chr_slen);
        }
        *write_iter++ = '\t';
        // POS
        if (check_pos_col) {
          const char* pos_str = token_ptrs[1];
          // no need to support negative values here
          uint32_t cur_bp;
          if (unlikely(ScanUintDefcap(pos_str, &cur_bp))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid bp coordinate on line %" PRIuPTR " of %s.\n", line_idx, dosagename);
            goto Plink1DosageToPgen_ret_MALFORMED_INPUT_WW;
          }
          write_iter = u32toa(cur_bp, write_iter);
        } else {
          *write_iter++ = '0';
        }
        *write_iter++ = '\t';
        write_iter = memcpya(write_iter, variant_id, variant_id_slen);
      }
      ++variant_ct;
      *write_iter++ = '\t';
      // REF, ALT
      write_iter = memcpyax(write_iter, token_ptrs[3], token_slens[3], '\t');
      write_iter = memcpya(write_iter, token_ptrs[4], token_slens[4]);
      if (variant_cms) {
        *write_iter++ = '\t';
        write_iter = dtoa_g(variant_cms[variant_uidx], write_iter);
      }
      AppendBinaryEoln(&write_iter);
      if (unlikely(fwrite_ck(writebuf_flush, outfile, &write_iter))) {
        goto Plink1DosageToPgen_ret_WRITE_FAIL;
      }
      if (!dosage_is_present) {
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
            if ((!linebuf_iter) || (a1_dosage < (0.5 / 32768.0)) || (a1_dosage >= dosage_ceil)) {
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
              dosage_is_present = 1;
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
              dosage_is_present = 1;
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
    if (unlikely(fclose_flush_null(writebuf_flush, write_iter, &outfile))) {
      goto Plink1DosageToPgen_ret_WRITE_FAIL;
    }
    if (unlikely(!variant_ct)) {
      if (!variant_skip_ct) {
        logerrputs("Error: Empty --import-dosage file.\n");
        goto Plink1DosageToPgen_ret_INCONSISTENT_INPUT;
      }
      logerrprintfww("Error: All %" PRIuPTR " variant%s in --import-dosage file skipped.\n", variant_skip_ct, (variant_skip_ct == 1)? "" : "s");
      goto Plink1DosageToPgen_ret_INCONSISTENT_INPUT;
    }
    logprintf("--import-dosage: %u variant%s scanned%s.\n", variant_ct, (variant_ct == 1)? "" : "s", dosage_is_present? "" : " (all hardcalls)");

    // 5. Dosage file pass 2: write .pgen.
    BigstackReset(writebuf);
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
    reterr = SpgwInitPhase1(outname, nullptr, nullptr, variant_ct, sample_ct, dosage_is_present? kfPgenGlobalDosagePresent: kfPgenGlobal0, (flags & (kfPlink1DosageRefFirst | kfPlink1DosageRefLast))? 1 : 2, &spgw, &spgw_alloc_cacheline_ct, &max_vrec_len);
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
    uintptr_t* dosage_present;
    if (unlikely(
            bigstack_alloc_w(sample_ctl2, &genovec) ||
            bigstack_alloc_w(sample_ctl, &dosage_present))) {
      goto Plink1DosageToPgen_ret_NOMEM;
    }
    Dosage* dosage_main = nullptr;
    if (dosage_is_present) {
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
        R_CAST(Halfword*, dosage_present)[widx] = dosage_present_hw;
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
        reterr = SpgwAppendBiallelicGenovecDosage16(genovec, dosage_present, dosage_main, dosage_ct, &spgw);
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
    SpgwFinish(&spgw);
    putc_unlocked('\r', stdout);
    write_iter = strcpya_k(g_logbuf, "--import-dosage: ");
    const uint32_t outname_base_slen = outname_end - outname;
    write_iter = memcpya(write_iter, outname, outname_base_slen + 5);
    write_iter = strcpya_k(write_iter, " + ");
    write_iter = memcpya(write_iter, outname, outname_base_slen);
    write_iter = strcpya_k(write_iter, ".pvar + ");
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
  }
 Plink1DosageToPgen_ret_1:
  CleanupSpgw(&spgw, &reterr);
  ForgetExtraChrNames(1, cip);
  fclose_cond(outfile);
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
  STD_ARRAY_DECL(uint64_t, kBitsPerWordD2, geno_missing_geomdist);
  STD_ARRAY_DECL(uint64_t, kBitsPerWordD2, dosage_geomdist);
  uint32_t geno_missing_invert;
  uint32_t dosage_geomdist_max;
  uint32_t hard_call_halfdist;
  uint32_t dosage_erase_halfdist;

  sfmt_t** sfmtp_arr;

  uint32_t cur_block_write_ct;

  uintptr_t* write_genovecs[2];
  uint32_t* write_dosage_cts[2];
  uintptr_t* write_dosage_presents[2];
  Dosage* write_dosage_mains[2];
} GenerateDummyCtx;

static_assert(sizeof(Dosage) == 2, "GenerateDummyThread() needs to be updated.");
THREAD_FUNC_DECL GenerateDummyThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  GenerateDummyCtx* ctx = S_CAST(GenerateDummyCtx*, arg->sharedp->context);

  const uint32_t sample_ct = ctx->sample_ct;
  const uint32_t calc_thread_ct = GetThreadCt(arg->sharedp);
  STD_ARRAY_KREF(uint64_t, kBitsPerWordD2) geno_missing_geomdist = ctx->geno_missing_geomdist;
  STD_ARRAY_KREF(uint64_t, kBitsPerWordD2) dosage_geomdist = ctx->dosage_geomdist;
  const uint32_t geno_missing_invert = ctx->geno_missing_invert;
  const uint32_t geno_missing_check = geno_missing_invert || (geno_missing_geomdist[kBitsPerWordD2 - 1] != 0);
  const uint32_t dosage_geomdist_max = ctx->dosage_geomdist_max;
  const uint32_t dosage_is_present = (dosage_geomdist_max != kBitsPerWord);
  const uint32_t hard_call_halfdist = ctx->hard_call_halfdist;
  const uint32_t dosage_erase_halfdist = ctx->dosage_erase_halfdist;
  const uintptr_t sample_ctaw2 = NypCtToAlignedWordCt(sample_ct);
  const uint32_t sample_ctl2_m1 = (sample_ct - 1) / kBitsPerWordD2;
  const uintptr_t sample_ctaw = BitCtToAlignedWordCt(sample_ct);
  sfmt_t* sfmtp = ctx->sfmtp_arr[tidx];
  uint64_t u64rand = sfmt_genrand_uint64(sfmtp);
  uint32_t rand16_left = 4;
  uint32_t parity = 0;
  do {
    const uintptr_t cur_block_write_ct = ctx->cur_block_write_ct;
    uint32_t vidx = (tidx * cur_block_write_ct) / calc_thread_ct;
    const uint32_t vidx_end = ((tidx + 1) * cur_block_write_ct) / calc_thread_ct;
    uintptr_t* write_genovec_iter = &(ctx->write_genovecs[parity][vidx * sample_ctaw2]);
    uint32_t* write_dosage_ct_iter = &(ctx->write_dosage_cts[parity][vidx]);
    uintptr_t* write_dosage_present_iter = &(ctx->write_dosage_presents[parity][vidx * sample_ctaw]);

    // bugfix (23 Jul 2017): multiply by sample_ct, not sample_ctaw
    Dosage* write_dosage_main_iter = &(ctx->write_dosage_mains[parity][vidx * sample_ct]);
    for (; vidx != vidx_end; ++vidx) {
      Dosage* cur_dosage_main_iter = write_dosage_main_iter;
      uint32_t loop_len = kBitsPerWordD2;
      for (uint32_t widx = 0; ; ++widx) {
        if (widx >= sample_ctl2_m1) {
          if (widx > sample_ctl2_m1) {
            break;
          }
          loop_len = ModNz(sample_ct, kBitsPerWordD2);
        }
        // sfmt_genrand_uint64 calls can't be mixed with sfmt_genrand_uint32
        // calls, so use it here even in 32-bit build
        uintptr_t genovec_word = sfmt_genrand_uint64(sfmtp);
        genovec_word = genovec_word - ((genovec_word >> 1) & kMask5555);
        if (geno_missing_check) {
          uintptr_t missing_mask = 0;
          uint32_t sample_idx_lowbits = 0;
          while (1) {
            sample_idx_lowbits += CountSortedLeqU64(&geno_missing_geomdist[0], kBitsPerWordD2, sfmt_genrand_uint64(sfmtp));
            if (sample_idx_lowbits >= loop_len) {
              break;
            }
            missing_mask |= (3 * k1LU) << (2 * sample_idx_lowbits);
            ++sample_idx_lowbits;
          }
          if (geno_missing_invert) {
            missing_mask = ~missing_mask;
          }
          genovec_word |= missing_mask;
        }
        uint32_t dosage_present_hw = 0;
        if (dosage_is_present) {
          // deliberate overflow
          uint32_t sample_idx_lowbits = UINT32_MAX;
          while (1) {
            ++sample_idx_lowbits;
            if (dosage_geomdist_max) {
              sample_idx_lowbits += CountSortedLeqU64(&dosage_geomdist[0], dosage_geomdist_max, sfmt_genrand_uint64(sfmtp));
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
            const uint32_t dosage_int = ((u64rand & 65535) + 1) / 2;
            u64rand >>= 16;
            --rand16_left;
            const uint32_t halfdist = BiallelicDosageHalfdist(dosage_int);
            if (halfdist < dosage_erase_halfdist) {
              *cur_dosage_main_iter++ = dosage_int;
              dosage_present_hw |= 1U << sample_idx_lowbits;
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
        R_CAST(Halfword*, write_dosage_present_iter)[widx] = dosage_present_hw;
      }
      ZeroTrailingNyps(sample_ct, write_genovec_iter);
      const uint32_t dosage_ct = cur_dosage_main_iter - write_dosage_main_iter;
      *write_dosage_ct_iter++ = dosage_ct;
      write_genovec_iter = &(write_genovec_iter[sample_ctaw2]);
      write_dosage_present_iter = &(write_dosage_present_iter[sample_ctaw]);
      write_dosage_main_iter = &(write_dosage_main_iter[sample_ct]);
    }
    parity = 1 - parity;
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

static_assert(sizeof(Dosage) == 2, "GenerateDummy() needs to be updated.");
PglErr GenerateDummy(const GenDummyInfo* gendummy_info_ptr, MiscFlags misc_flags, ImportFlags import_flags, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, uint32_t max_thread_ct, sfmt_t* sfmtp, char* outname, char* outname_end, ChrInfo* cip) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  ThreadGroup tg;
  PreinitThreads(&tg);
  STPgenWriter spgw;
  PglErr reterr = kPglRetSuccess;
  PreinitSpgw(&spgw);
  {
    FinalizeChrset(misc_flags, cip);
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
    char* textbuf = g_textbuf;
    char* textbuf_flush = &(textbuf[kMaxMediumLine]);
    snprintf(outname_end, kMaxOutfnameExtBlen, ".pvar");
    if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
      goto GenerateDummy_ret_OPEN_FAIL;
    }
    char* write_iter = strcpya_k(textbuf, "#CHROM\tPOS\tID\tREF\tALT");
    AppendBinaryEoln(&write_iter);
    if (four_alleles) {
      uint32_t urand = 0;
      for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
        if (!(variant_idx % 8)) {
          if (unlikely(fwrite_ck(textbuf_flush, outfile, &write_iter))) {
            goto GenerateDummy_ret_WRITE_FAIL;
          }
          do {
            urand = sfmt_genrand_uint32(sfmtp);
          } while (urand < 425132032U);  // 2^32 - 12^8
        }
        const uint32_t quotient = urand / 12;
        const uint32_t remainder = urand - (quotient * 12U);
        urand = quotient;
        write_iter = memcpya(write_iter, chr1_name_buf, chr1_name_blen);
        write_iter = u32toa(variant_idx, write_iter);
        write_iter = strcpya_k(write_iter, "\tsnp");
        write_iter = u32toa(variant_idx, write_iter);
        write_iter = memcpya(write_iter, &(alleles[remainder]), 4);
        AppendBinaryEoln(&write_iter);
      }
    } else {
      uint32_t urand = 0;
      for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
        if (!(variant_idx % 32)) {
          if (unlikely(fwrite_ck(textbuf_flush, outfile, &write_iter))) {
            goto GenerateDummy_ret_WRITE_FAIL;
          }
          urand = sfmt_genrand_uint32(sfmtp);
        }
        const uint32_t remainder = urand & 1;
        urand >>= 1;
        write_iter = memcpya(write_iter, chr1_name_buf, chr1_name_blen);
        write_iter = u32toa(variant_idx, write_iter);
        write_iter = strcpya_k(write_iter, "\tsnp");
        write_iter = u32toa(variant_idx, write_iter);
        write_iter = memcpya(write_iter, &(alleles[remainder]), 4);
        AppendBinaryEoln(&write_iter);
      }
    }
    if (unlikely(fclose_flush_null(textbuf_flush, write_iter, &outfile))) {
      goto GenerateDummy_ret_WRITE_FAIL;
    }

    snprintf(outname_end, kMaxOutfnameExtBlen, ".psam");
    if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
      goto GenerateDummy_ret_OPEN_FAIL;
    }
    const uint32_t pheno_ct = gendummy_info_ptr->pheno_ct;
    char* writebuf;
    if (unlikely(bigstack_alloc_c(kMaxMediumLine + 48 + pheno_ct * MAXV(kMaxMissingPhenostrBlen, 16), &writebuf))) {
      goto GenerateDummy_ret_NOMEM;
    }
    char* writebuf_flush = &(writebuf[kMaxMediumLine]);
    uint32_t omp_slen = 2;
    char output_missing_pheno[kMaxMissingPhenostrBlen];
    if (import_flags & kfImportKeepAutoconv) {
      // must use --output-missing-phenotype parameter, which we've validated
      // to be consistent with --input-missing-phenotype
      omp_slen = strlen(g_output_missing_pheno);
      memcpy(output_missing_pheno, g_output_missing_pheno, omp_slen);
    } else {
      // use "NA" since that's always safe
      memcpy_k(output_missing_pheno, "NA", 2);
    }
    // Alpha 2 change: no more FID column
    write_iter = strcpya_k(writebuf, "#IID\tSEX");
    for (uint32_t pheno_idx_p1 = 1; pheno_idx_p1 <= pheno_ct; ++pheno_idx_p1) {
      write_iter = strcpya_k(write_iter, "\tPHENO");
      write_iter = u32toa(pheno_idx_p1, write_iter);
    }
    AppendBinaryEoln(&write_iter);
    const uint32_t pheno_m_check = (gendummy_info_ptr->pheno_mfreq >= kRecip2m32 * 0.5);
    const uint32_t pheno_m32 = S_CAST(uint32_t, gendummy_info_ptr->pheno_mfreq * 4294967296.0 - 0.5);
    if ((flags & kfGenDummyScalarPheno) && pheno_ct) {
      uint32_t saved_rnormal = 0;
      double saved_rnormal_val = 0.0;
      for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        if (unlikely(fwrite_ck(writebuf_flush, outfile, &write_iter))) {
          goto GenerateDummy_ret_WRITE_FAIL;
        }
        write_iter = strcpya_k(write_iter, "per");
        write_iter = u32toa(sample_idx, write_iter);
        // could add option to add some males/unknown gender
        write_iter = strcpya_k(write_iter, "\t2");
        for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
          *write_iter++ = '\t';
          if (pheno_m_check && (sfmt_genrand_uint32(sfmtp) <= pheno_m32)) {
            write_iter = memcpya(write_iter, output_missing_pheno, omp_slen);
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
      for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        if (unlikely(fwrite_ck(writebuf_flush, outfile, &write_iter))) {
          goto GenerateDummy_ret_WRITE_FAIL;
        }
        // bugfix (9 Mar 2018): forgot to remove FID column here
        write_iter = strcpya_k(write_iter, "per");
        write_iter = u32toa(sample_idx, write_iter);
        write_iter = strcpya_k(write_iter, "\t2");
        for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
          *write_iter++ = '\t';
          if (pheno_m_check && (sfmt_genrand_uint32(sfmtp) <= pheno_m32)) {
            write_iter = memcpya(write_iter, output_missing_pheno, omp_slen);
          } else {
            if (!urand_bits_left) {
              urand = sfmt_genrand_uint32(sfmtp);
              urand_bits_left = 32;
            }
            *write_iter++ = (urand & 1) + '1';
            urand >>= 1;
            --urand_bits_left;
          }
        }
        AppendBinaryEoln(&write_iter);
      }
    }
    if (unlikely(fclose_flush_null(writebuf_flush, write_iter, &outfile))) {
      goto GenerateDummy_ret_WRITE_FAIL;
    }

    BigstackReset(writebuf);
    snprintf(outname_end, kMaxOutfnameExtBlen, ".pgen");
    const double geno_mfreq = gendummy_info_ptr->geno_mfreq;
    GenerateDummyCtx ctx;
    ctx.geno_missing_invert = 0;
    if (geno_mfreq < kRecip2m53) {
      // beyond this point, 1-x may just be 1
      ctx.geno_missing_geomdist[kBitsPerWordD2 - 1] = 0;
    } else {
      double remaining_prob = 1.0;
      ctx.geno_missing_invert = (geno_mfreq > 0.5);
      if (ctx.geno_missing_invert) {
        for (uint32_t uii = 0; uii != kBitsPerWordD2; ++uii) {
          remaining_prob *= geno_mfreq;
          ctx.geno_missing_geomdist[uii] = -S_CAST(uint64_t, remaining_prob * k2m64);
        }
      } else {
        const double geno_nmfreq = 1.0 - geno_mfreq;
        for (uint32_t uii = 0; uii != kBitsPerWordD2; ++uii) {
          remaining_prob *= geno_nmfreq;
          ctx.geno_missing_geomdist[uii] = -S_CAST(uint64_t, remaining_prob * k2m64);
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
        ctx.dosage_geomdist[uii] = -S_CAST(uint64_t, remaining_prob * k2m64);
      }
      uint32_t dosage_geomdist_max = kBitsPerWordD2;
      for (; dosage_geomdist_max; --dosage_geomdist_max) {
        if (ctx.dosage_geomdist[dosage_geomdist_max - 1] != 0) {
          break;
        }
      }
      ctx.dosage_geomdist_max = dosage_geomdist_max;
    }
    uintptr_t spgw_alloc_cacheline_ct;
    uint32_t max_vrec_len;
    reterr = SpgwInitPhase1(outname, nullptr, nullptr, variant_ct, sample_ct, (dosage_nfreq >= 1.0)? kfPgenGlobal0 : kfPgenGlobalDosagePresent, 1, &spgw, &spgw_alloc_cacheline_ct, &max_vrec_len);
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
    uintptr_t cachelines_avail_m8 = bigstack_left() / kCacheline;
    if (unlikely(cachelines_avail_m8 < 8)) {
      goto GenerateDummy_ret_NOMEM;
    }
    // we're making 8 allocations; be pessimistic re: rounding
    cachelines_avail_m8 -= 8;
    const uintptr_t bytes_req_per_in_block_variant = 2 * (sample_ctaw2 * sizeof(intptr_t) + sizeof(int32_t) + sample_ctaw * sizeof(intptr_t) + sample_ct * sizeof(Dosage));
    uintptr_t main_block_size = (cachelines_avail_m8 * kCacheline) / bytes_req_per_in_block_variant;
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
    if (unlikely(
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
        uint32_t* write_dosage_ct_iter = ctx.write_dosage_cts[parity];
        uintptr_t* write_dosage_present_iter = ctx.write_dosage_presents[parity];
        Dosage* write_dosage_main_iter = ctx.write_dosage_mains[parity];
        for (uint32_t vidx = vidx_start - prev_block_write_ct; vidx != vidx_start; ++vidx) {
          const uint32_t cur_dosage_ct = *write_dosage_ct_iter++;
          if (!cur_dosage_ct) {
            if (unlikely(SpgwAppendBiallelicGenovec(write_genovec_iter, &spgw))) {
              goto GenerateDummy_ret_WRITE_FAIL;
            }
          } else {
            reterr = SpgwAppendBiallelicGenovecDosage16(write_genovec_iter, write_dosage_present_iter, write_dosage_main_iter, cur_dosage_ct, &spgw);
            if (unlikely(reterr)) {
              goto GenerateDummy_ret_1;
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
        printf("\r--dummy: %uk variants written.", vidx_start / 1000);
        fflush(stdout);
      }
      vidx_start += cur_block_write_ct;
      prev_block_write_ct = cur_block_write_ct;
    }
    SpgwFinish(&spgw);

    putc_unlocked('\r', stdout);
    *outname_end = '\0';
    logprintfww("Dummy data (%u sample%s, %u SNP%s) written to %s.pgen + %s.pvar + %s.psam .\n", sample_ct, (sample_ct == 1)? "" : "s", variant_ct, (variant_ct == 1)? "" : "s", outname, outname, outname);
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
  fclose_cond(outfile);
  BigstackReset(bigstack_mark);
  return reterr;
}


typedef struct Plink1SmajTransposeCtxStruct {
  uint32_t sample_ct;
  uint32_t loadbuf_ul_stride;

  uintptr_t* plink1_smaj_loadbuf_iter;

  VecW** thread_vecaligned_bufs;
  uintptr_t** thread_write_genovecs;
  PgenWriterCommon** pwcs;

  uint32_t cur_block_write_ct;
} Plink1SmajTransposeCtx;

void Plink1SmajTransposeMain(uint32_t tidx, Plink1SmajTransposeCtx* ctx) {
  const uint32_t sample_ct = ctx->sample_ct;
  const uint32_t sample_ctaw2 = NypCtToAlignedWordCt(sample_ct);
  const uint32_t write_batch_ct_m1 = (sample_ct - 1) / kPglNypTransposeBatch;
  PgenWriterCommon* pwcp = ctx->pwcs[tidx];
  VecW* vecaligned_buf = ctx->thread_vecaligned_bufs[tidx];
  uintptr_t* write_genovec = ctx->thread_write_genovecs[tidx];
  const uintptr_t cur_block_write_ct = ctx->cur_block_write_ct;
  const uint32_t loadbuf_ul_stride = ctx->loadbuf_ul_stride;
  uint32_t write_idx = tidx * kPglVblockSize;
  uintptr_t* read_iter = &(ctx->plink1_smaj_loadbuf_iter[write_idx / kBitsPerWordD2]);
  const uint32_t write_idx_end = MINV(write_idx + kPglVblockSize, cur_block_write_ct);
  while (write_idx < write_idx_end) {
    const uintptr_t* read_iter2 = read_iter;
    // uintptr_t* write_iter = write_genovec;
    const uint32_t vblock_size = MINV(kPglNypTransposeBatch, write_idx_end - write_idx);
    uint32_t read_batch_size = kPglNypTransposeBatch;
    for (uint32_t write_batch_idx = 0; ; ++write_batch_idx) {
      if (write_batch_idx >= write_batch_ct_m1) {
        if (write_batch_idx > write_batch_ct_m1) {
          break;
        }
        read_batch_size = ModNz(sample_ct, kPglNypTransposeBatch);
      }
      TransposeNypblock(read_iter2, loadbuf_ul_stride, sample_ctaw2, read_batch_size, vblock_size, &(write_genovec[write_batch_idx * kPglNypTransposeWords]), vecaligned_buf);
      read_iter2 = &(read_iter2[kPglNypTransposeBatch * loadbuf_ul_stride]);
    }
    for (uint32_t uii = 0; uii != vblock_size; ++uii) {
      uintptr_t* cur_write_genovec = &(write_genovec[uii * sample_ctaw2]);
      PgrPlink1ToPlink2InplaceUnsafe(sample_ct, cur_write_genovec);
      ZeroTrailingNyps(sample_ct, cur_write_genovec);
      PwcAppendBiallelicGenovec(cur_write_genovec, pwcp);
    }
    write_idx += vblock_size;
    read_iter = &(read_iter[kPglNypTransposeWords]);
  }
}

THREAD_FUNC_DECL Plink1SmajTransposeThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uint32_t tidx = arg->tidx;
  Plink1SmajTransposeCtx* ctx = S_CAST(Plink1SmajTransposeCtx*, arg->sharedp->context);
  do {
    Plink1SmajTransposeMain(tidx, ctx);
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

PglErr Plink1SampleMajorToPgen(const char* pgenname, uintptr_t variant_ct, uintptr_t sample_ct, uint32_t real_ref_alleles, uint32_t max_thread_ct, FILE* infile) {
  unsigned char* bigstack_mark = g_bigstack_base;
  MTPgenWriter* mpgwp = nullptr;
  ThreadGroup tg;
  PreinitThreads(&tg);
  Plink1SmajTransposeCtx ctx;
  PglErr reterr = kPglRetSuccess;
  {
    // file size already validated by PgfiInitPhase1()
    logprintfww("Sample-major .bed file detected.  Transposing to %s .\n", pgenname);
    fputs("0%", stdout);
    fflush(stdout);
    if (unlikely((!variant_ct) || (!sample_ct))) {
      logputs("\n");
      logerrputs("Error: Zero-variant/zero-sample .pgen writing is not supported.\n");
      goto Plink1SampleMajorToPgen_ret_INCONSISTENT_INPUT;
    }
    const uint32_t variant_ct4 = NypCtToByteCt(variant_ct);
    unsigned char* raw_loadbuf = nullptr;
    uint32_t raw_load_batch_size = 1;
    if (variant_ct4 < 5120) {
      // assuming 4K block size, fseek won't let us avoid reading many
      // unnecessary disk blocks
      raw_load_batch_size += 131071 / variant_ct4;
      if (unlikely(bigstack_alloc_uc(raw_load_batch_size * variant_ct4, &raw_loadbuf))) {
        goto Plink1SampleMajorToPgen_ret_NOMEM;
      }
    }
    const uint32_t raw_load_batch_ct_m1 = (sample_ct - 1) / raw_load_batch_size;
    if (!raw_load_batch_ct_m1) {
      raw_load_batch_size = sample_ct;
    }
    const uint32_t raw_load_batch_ct = raw_load_batch_ct_m1 + 1;
    uintptr_t alloc_base_cacheline_ct;
    uint64_t mpgw_per_thread_cacheline_ct;
    uint32_t vrec_len_byte_ct;
    uint64_t vblock_cacheline_ct;
    MpgwInitPhase1(nullptr, variant_ct, sample_ct, kfPgenGlobal0, &alloc_base_cacheline_ct, &mpgw_per_thread_cacheline_ct, &vrec_len_byte_ct, &vblock_cacheline_ct);
#ifndef __LP64__
    if (unlikely((mpgw_per_thread_cacheline_ct > (0x7fffffff / kCacheline)) || (vblock_cacheline_ct > (0x7fffffff / kCacheline)))) {
      goto Plink1SampleMajorToPgen_ret_NOMEM;
    }
#endif

    uint32_t calc_thread_ct = DivUp(variant_ct, kPglVblockSize);
    if (calc_thread_ct >= max_thread_ct) {
      calc_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
    }
    // note that BIGSTACK_ALLOC_X() doesn't work here due to variable-size
    // array at end
    mpgwp = S_CAST(MTPgenWriter*, bigstack_alloc((calc_thread_ct + DivUp(sizeof(MTPgenWriter), kBytesPerWord)) * sizeof(intptr_t)));
    if (unlikely(!mpgwp)) {
      goto Plink1SampleMajorToPgen_ret_NOMEM;
    }
    mpgwp->pgen_outfile = nullptr;
    if (unlikely(
            bigstack_alloc_vp(calc_thread_ct, &ctx.thread_vecaligned_bufs) ||
            bigstack_alloc_wp(calc_thread_ct, &ctx.thread_write_genovecs))) {
      goto Plink1SampleMajorToPgen_ret_NOMEM;
    }
    ctx.pwcs = &(mpgwp->pwcs[0]);
    uintptr_t cachelines_avail = bigstack_left() / kCacheline;
    // inner loop transposes kPglNypTransposeBatch variants at a time
    const uintptr_t transpose_thread_cacheline_ct = kPglNypTransposeBufbytes / kCacheline + NypCtToVecCt(sample_ct) * (kPglNypTransposeBatch / kVecsPerCacheline);
    if (unlikely(cachelines_avail < calc_thread_ct * S_CAST(uint64_t, transpose_thread_cacheline_ct))) {
      goto Plink1SampleMajorToPgen_ret_NOMEM;
    }
    for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
      ctx.thread_vecaligned_bufs[tidx] = S_CAST(VecW*, bigstack_alloc_raw(kPglNypTransposeBufbytes));
      ctx.thread_write_genovecs[tidx] = S_CAST(uintptr_t*, bigstack_alloc_raw(NypCtToVecCt(sample_ct) * kBytesPerVec * kPglNypTransposeBatch));
    }
    cachelines_avail = bigstack_left() / kCacheline;
    // Main workflow:
    // 1. Load next calc_thread_ct * load_multiplier * kPglVblockSize
    //    variants.
    //    calc_thread_ct is reduced as necessary to ensure the compression
    //    write buffers use <= 1/8 of total workspace.
    //    with calc_thread_ct determined, load_multiplier is then chosen to use
    //    as much of the remaining workspace as possible.
    // 2. Repeat load_multiplier times:
    //    a. Spawn threads processing calc_thread_ct vblocks
    //    b. Join threads
    //    c. Flush results
    // 3. Goto step 1 unless eof.  (load_multiplier may be smaller on last
    //    iteration.)
    // No double-buffering here since main bottleneck is how many variants we
    // can load at once.
    if ((cachelines_avail / 8) < alloc_base_cacheline_ct + mpgw_per_thread_cacheline_ct * calc_thread_ct) {
      if ((cachelines_avail / 8) >= alloc_base_cacheline_ct + mpgw_per_thread_cacheline_ct) {
        calc_thread_ct = ((cachelines_avail / 8) - alloc_base_cacheline_ct) / mpgw_per_thread_cacheline_ct;
      } else if (likely((cachelines_avail / 3) >= alloc_base_cacheline_ct + mpgw_per_thread_cacheline_ct)) {
        calc_thread_ct = 1;
      } else {
        // possible todo: simple single-threaded fallback
        // report this value since it can plausibly come up
        g_failed_alloc_attempt_size = (alloc_base_cacheline_ct + mpgw_per_thread_cacheline_ct) * kCacheline * 3;
        goto Plink1SampleMajorToPgen_ret_NOMEM;
      }
    }
    // todo: determine appropriate calc_thread_ct limit.  (should not be less
    // than 7-8.)
    unsigned char* mpgw_alloc = S_CAST(unsigned char*, bigstack_alloc_raw((alloc_base_cacheline_ct + mpgw_per_thread_cacheline_ct * calc_thread_ct) * kCacheline));
    reterr = MpgwInitPhase2(pgenname, nullptr, nullptr, variant_ct, sample_ct, kfPgenGlobal0, 2 - real_ref_alleles, vrec_len_byte_ct, vblock_cacheline_ct, calc_thread_ct, mpgw_alloc, mpgwp);
    if (unlikely(reterr)) {
      if (reterr == kPglRetOpenFail) {
        logputs("\n");
        logerrprintfww(kErrprintfFopen, pgenname, strerror(errno));
      }
      goto Plink1SampleMajorToPgen_ret_1;
    }

    cachelines_avail = bigstack_left() / kCacheline;
    const uint64_t full_load_vecs_req = sample_ct * S_CAST(uint64_t, NypCtToAlignedWordCt(variant_ct));
    uintptr_t* plink1_smaj_loadbuf;
    uint32_t load_multiplier;
    uint32_t cur_vidx_ct;
    if (full_load_vecs_req > cachelines_avail * kVecsPerCacheline) {
      // each iteration requires ((kPglVblockSize / 4) * calc_thread_ct *
      //   sample_ct) bytes to be loaded
      load_multiplier = cachelines_avail / ((kPglVblockSize / (4 * kCacheline)) * calc_thread_ct * S_CAST(uintptr_t, sample_ct));
      assert(load_multiplier);
      cur_vidx_ct = load_multiplier * calc_thread_ct * kPglVblockSize;
      plink1_smaj_loadbuf = S_CAST(uintptr_t*, bigstack_alloc_raw_rd((cur_vidx_ct / 4) * S_CAST(uintptr_t, sample_ct)));
      // bugfix (18 Nov 2017): this may be larger than variant_ct
      if (cur_vidx_ct > variant_ct) {
        cur_vidx_ct = variant_ct;
        load_multiplier = 1 + (cur_vidx_ct - 1) / (kPglVblockSize * calc_thread_ct);
      }
    } else {
      load_multiplier = 1 + ((variant_ct - 1) / (calc_thread_ct * kPglVblockSize));
      cur_vidx_ct = variant_ct;
      plink1_smaj_loadbuf = S_CAST(uintptr_t*, bigstack_alloc_raw_rd(full_load_vecs_req * kBytesPerVec));
    }
    uint32_t cur_vidx_ct4 = NypCtToByteCt(cur_vidx_ct);
    uint32_t cur_vidx_ctaw2 = NypCtToAlignedWordCt(cur_vidx_ct);
    const uint32_t pass_ct = 1 + (variant_ct - 1) / cur_vidx_ct;
    ctx.sample_ct = sample_ct;
    ctx.loadbuf_ul_stride = NypCtToVecCt(cur_vidx_ct) * kWordsPerVec;
    if (unlikely(SetThreadCt0(calc_thread_ct - 1, &tg))) {
      goto Plink1SampleMajorToPgen_ret_NOMEM;
    }
    SetThreadFuncAndData(Plink1SmajTransposeThread, &ctx, &tg);
    uint32_t pass_idx1 = 0;
    for (uint32_t cur_vidx_base = 0; ; ) {
      uint32_t cur_raw_load_batch_size = raw_load_batch_size;
      uintptr_t* smaj_loadbuf_iter = plink1_smaj_loadbuf;
      putc_unlocked('\r', stdout);
      ++pass_idx1;
      printf("Pass %u/%u: loading... 0%%", pass_idx1, pass_ct);
      fflush(stdout);
      uint32_t pct = 0;
      uint32_t next_print_idx = raw_load_batch_ct / 100;
      const uint64_t seek_addl_offset = 3 + cur_vidx_base / 4;
      for (uint32_t raw_load_batch_idx = 0; ; ) {
        if (raw_load_batch_size == 1) {
          if (unlikely(fseeko(infile, seek_addl_offset + raw_load_batch_idx * S_CAST(uint64_t, variant_ct4), SEEK_SET))) {
            goto Plink1SampleMajorToPgen_ret_READ_FAIL;
          }
          if (unlikely(!fread_unlocked(smaj_loadbuf_iter, cur_vidx_ct4, 1, infile))) {
            goto Plink1SampleMajorToPgen_ret_READ_FAIL;
          }
          smaj_loadbuf_iter = &(smaj_loadbuf_iter[cur_vidx_ctaw2]);
        } else {
          if (unlikely(!fread_unlocked(raw_loadbuf, cur_raw_load_batch_size * variant_ct4, 1, infile))) {
            goto Plink1SampleMajorToPgen_ret_READ_FAIL;
          }
          unsigned char* raw_loadbuf_iter = &(raw_loadbuf[cur_vidx_base / 4]);
          for (uint32_t uii = 0; uii != cur_raw_load_batch_size; ++uii) {
            memcpy(smaj_loadbuf_iter, raw_loadbuf_iter, cur_vidx_ct4);
            raw_loadbuf_iter = &(raw_loadbuf_iter[variant_ct4]);
            smaj_loadbuf_iter = &(smaj_loadbuf_iter[cur_vidx_ctaw2]);
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
          next_print_idx = (pct * S_CAST(uint64_t, raw_load_batch_ct)) / 100;
        }
      }
      const uintptr_t last_tidx = calc_thread_ct - 1;
      uint32_t load_idx = 0;
      ctx.cur_block_write_ct = calc_thread_ct * kPglVblockSize;
      putc_unlocked('\r', stdout);
      printf("Pass %u/%u: transposing and compressing... 0%%", pass_idx1, pass_ct);
      ReinitThreads(&tg);
      pct = 0;
      next_print_idx = load_idx / 100;
      do {
        if (load_idx >= next_print_idx) {
          if (pct > 10) {
            putc_unlocked('\b', stdout);
          }
          pct = (load_idx * 100LLU) / load_multiplier;
          printf("\b\b%u%%", pct++);
          fflush(stdout);
          next_print_idx = (pct * S_CAST(uint64_t, load_multiplier)) / 100;
        }
        ctx.plink1_smaj_loadbuf_iter = &(plink1_smaj_loadbuf[load_idx * calc_thread_ct * (kPglVblockSize / kBitsPerWordD2)]);
        if (++load_idx == load_multiplier) {
          DeclareLastThreadBlock(&tg);
          ctx.cur_block_write_ct = cur_vidx_ct - (load_idx - 1) * calc_thread_ct * kPglVblockSize;
        }
        if (last_tidx) {
          if (unlikely(SpawnThreads(&tg))) {
            goto Plink1SampleMajorToPgen_ret_THREAD_CREATE_FAIL;
          }
        }
        Plink1SmajTransposeMain(last_tidx, &ctx);
        if (last_tidx) {
          JoinThreads(&tg);
          // Plink1SmajTransposeThread() never errors out
        }
        reterr = MpgwFlush(mpgwp);
        if (unlikely(reterr)) {
          goto Plink1SampleMajorToPgen_ret_WRITE_FAIL;
        }
      } while (!IsLastBlock(&tg));
      cur_vidx_base += cur_vidx_ct;
      if (cur_vidx_base == variant_ct) {
        if (pct > 10) {
          putc_unlocked('\b', stdout);
        }
        break;
      }
      fputs("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b                     ", stdout);
      // assumes PgfiInitPhase1() leaves file pointer at byte 3; otherwise,
      // necessary to put this at top of main loop
      if (unlikely(fseeko(infile, 3, SEEK_SET))) {
        goto Plink1SampleMajorToPgen_ret_READ_FAIL;
      }
      if (variant_ct - cur_vidx_base <= cur_vidx_ct) {
        cur_vidx_ct = variant_ct - cur_vidx_base;
        cur_vidx_ct4 = NypCtToByteCt(cur_vidx_ct);
        cur_vidx_ctaw2 = NypCtToAlignedWordCt(cur_vidx_ct);
        ctx.loadbuf_ul_stride = NypCtToVecCt(cur_vidx_ct) * kWordsPerVec;
        load_multiplier = 1 + (cur_vidx_ct - 1) / (kPglVblockSize * calc_thread_ct);
      }
    }
    mpgwp = nullptr;
    fputs("\b\bdone.\n", stdout);
    logprintf("Transpose complete.\n");
  }
  while (0) {
  Plink1SampleMajorToPgen_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  Plink1SampleMajorToPgen_ret_READ_FAIL:
    logerrprintfww(kErrprintfFread, ".bed file", strerror(errno));
    reterr = kPglRetReadFail;
    break;
  Plink1SampleMajorToPgen_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  Plink1SampleMajorToPgen_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  Plink1SampleMajorToPgen_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 Plink1SampleMajorToPgen_ret_1:
  CleanupThreads(&tg);
  CleanupMpgw(mpgwp, &reterr);
  BigstackReset(bigstack_mark);
  return reterr;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
