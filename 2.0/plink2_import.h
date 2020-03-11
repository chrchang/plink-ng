#ifndef __PLINK2_IMPORT_H__
#define __PLINK2_IMPORT_H__

// This file is part of PLINK 2.00, copyright (C) 2005-2020 Shaun Purcell,
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


#include "plink2_data.h"
#include "include/SFMT.h"

#ifdef __cplusplus
namespace plink2 {
#endif

FLAGSET_DEF_START()
  kfImport0,
  kfImportKeepAutoconv = (1 << 0),
  kfImportKeepAutoconvVzs = (1 << 1),
  kfImportDoubleId = (1 << 2),
  kfImportVcfRequireGt = (1 << 3)
FLAGSET_DEF_END(ImportFlags);

ENUM_U31_DEF_START()
  kVcfHalfCallReference,
  kVcfHalfCallHaploid,
  kVcfHalfCallMissing,
  kVcfHalfCallError,
  // gets converted to kVcfHalfCallError, but with a different error message
  kVcfHalfCallDefault
ENUM_U31_DEF_END(VcfHalfCall);

FLAGSET_DEF_START()
  kfOxfordImport0,
  kfOxfordImportRefFirst = (1 << 0),
  kfOxfordImportRefLast = (1 << 1),
  kfOxfordImportRefUnknown = (1 << 2),
  kfOxfordImportRefAll = ((kfOxfordImportRefUnknown * 2) - kfOxfordImportRefFirst),
  kfOxfordImportBgenSnpIdChr = (1 << 3)
FLAGSET_DEF_END(OxfordImportFlags);

FLAGSET_DEF_START()
  kfPlink1Dosage0,
  kfPlink1DosageNoheader = (1 << 0),
  kfPlink1DosageFormatSingle = (1 << 1),
  kfPlink1DosageFormatSingle01 = (1 << 2),
  kfPlink1DosageFormatDouble = (1 << 3),
  kfPlink1DosageFormatTriple = (1 << 4),
  kfPlink1DosageRefFirst = (1 << 5),
  kfPlink1DosageRefLast = (1 << 6)
FLAGSET_DEF_END(Plink1DosageFlags);

FLAGSET_DEF_START()
  kfGenDummy0,
  kfGenDummyAcgt = (1 << 0),
  kfGenDummy1234 = (1 << 1),
  kfGenDummy12 = (1 << 2),
  kfGenDummyScalarPheno = (1 << 3)
FLAGSET_DEF_END(GenDummyFlags);

typedef struct Plink1DosageInfoStruct {
  Plink1DosageFlags flags;
  STD_ARRAY_DECL(uint32_t, 3, skips);
  uint32_t chr_col_idx;  // 0-based
  uint32_t pos_col_idx;
  char id_delim;
} Plink1DosageInfo;

typedef struct GenDummyInfoStruct {
  GenDummyFlags flags;
  uint32_t sample_ct;
  uint32_t variant_ct;
  uint32_t pheno_ct;
  double geno_mfreq;
  double pheno_mfreq;
  double dosage_freq;
} GenDummyInfo;

void InitPlink1Dosage(Plink1DosageInfo* plink1_dosage_info_ptr);

void InitGenDummy(GenDummyInfo* gendummy_info_ptr);

PglErr VcfToPgen(const char* vcfname, const char* preexisting_psamname, const char* const_fid, const char* dosage_import_field, MiscFlags misc_flags, ImportFlags import_flags, uint32_t no_samples_ok, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, double import_dosage_certainty, char id_delim, char idspace_to, int32_t vcf_min_gq, int32_t vcf_min_dp, int32_t vcf_max_dp, VcfHalfCall halfcall_mode, FamCol fam_cols, uint32_t max_thread_ct, char* outname, char* outname_end, ChrInfo* cip, uint32_t* pgen_generated_ptr, uint32_t* psam_generated_ptr);

PglErr BcfToPgen(const char* bcfname, const char* preexisting_psamname, const char* const_fid, const char* dosage_import_field, MiscFlags misc_flags, ImportFlags import_flags, uint32_t no_samples_ok, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, double import_dosage_certainty, char id_delim, char idspace_to, int32_t vcf_min_gq, int32_t vcf_min_dp, int32_t vcf_max_dp, VcfHalfCall halfcall_mode, FamCol fam_cols, uint32_t max_thread_ct, char* outname, char* outname_end, ChrInfo* cip, uint32_t* pgen_generated_ptr, uint32_t* psam_generated_ptr);

PglErr OxGenToPgen(const char* genname, const char* samplename, const char* ox_single_chr_str, const char* ox_missing_code, MiscFlags misc_flags, ImportFlags import_flags, OxfordImportFlags oxford_import_flags, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, double import_dosage_certainty, uint32_t max_thread_ct, char* outname, char* outname_end, ChrInfo* cip);

PglErr OxBgenToPgen(const char* bgenname, const char* samplename, const char* const_fid, const char* ox_single_chr_str, const char* ox_missing_code, MiscFlags misc_flags, ImportFlags import_flags, OxfordImportFlags oxford_import_flags, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, double import_dosage_certainty, char id_delim, char idspace_to, uint32_t max_thread_ct, char* outname, char* outname_end, ChrInfo* cip);

PglErr OxHapslegendToPgen(const char* hapsname, const char* legendname, const char* samplename, const char* ox_single_chr_str, const char* ox_missing_code, MiscFlags misc_flags, ImportFlags import_flags, OxfordImportFlags oxford_import_flags, uint32_t max_thread_ct, char* outname, char* outname_end, ChrInfo* cip);

PglErr Plink1DosageToPgen(const char* dosagename, const char* famname, const char* mapname, const char* import_single_chr_str, const Plink1DosageInfo* pdip, MiscFlags misc_flags, ImportFlags import_flags, FamCol fam_cols, int32_t missing_pheno, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, double import_dosage_certainty, uint32_t max_thread_ct, char* outname, char* outname_end, ChrInfo* cip);

PglErr GenerateDummy(const GenDummyInfo* gendummy_info_ptr, MiscFlags misc_flags, ImportFlags import_flags, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, uint32_t max_thread_ct, sfmt_t* sfmtp, char* outname, char* outname_end, ChrInfo* cip);

PglErr Plink1SampleMajorToPgen(const char* pgenname, uintptr_t variant_ct, uintptr_t sample_ct, uint32_t real_ref_alleles, uint32_t max_thread_ct, FILE* infile);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_IMPORT_H__
