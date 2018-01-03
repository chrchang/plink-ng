#ifndef __PLINK2_IMPORT_H__
#define __PLINK2_IMPORT_H__

// This file is part of PLINK 2.00, copyright (C) 2005-2018 Shaun Purcell,
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

#ifdef __cplusplus
namespace plink2 {
#endif

ENUM_U31_DEF_START()
  kVcfHalfCallReference,
  kVcfHalfCallHaploid,
  kVcfHalfCallMissing,
  kVcfHalfCallError,
  // gets converted to kVcfHalfCallError, but with a different error message
  kVcfHalfCallDefault
ENUM_U31_DEF_END(vcf_half_call_t);

FLAGSET_DEF_START()
  kfOxfordImport0,
  kfOxfordImportRefFirst = (1 << 0),
  kfOxfordImportRefSecond = (1 << 1),
  kfOxfordImportBgenSnpIdChr = (1 << 2)
FLAGSET_DEF_END(oxford_import_t);

FLAGSET_DEF_START()
  kfPlink1Dosage0,
  kfPlink1DosageNoheader = (1 << 0),
  kfPlink1DosageFormatSingle = (1 << 1),
  kfPlink1DosageFormatSingle01 = (1 << 2),
  kfPlink1DosageFormatDouble = (1 << 3),
  kfPlink1DosageFormatTriple = (1 << 4),
  kfPlink1DosageRefFirst = (1 << 5),
  kfPlink1DosageRefSecond = (1 << 6)
FLAGSET_DEF_END(plink1_dosage_flags_t);

FLAGSET_DEF_START()
  kfGenDummy0,
  kfGenDummyAcgt = (1 << 0),
  kfGenDummy1234 = (1 << 1),
  kfGenDummy12 = (1 << 2),
  kfGenDummyScalarPheno = (1 << 3)
FLAGSET_DEF_END(gendummy_flags_t);

typedef struct plink1_dosage_info_struct {
  plink1_dosage_flags_t flags;
  uint32_t skips[3];
  uint32_t chr_col_idx; // 0-based
  uint32_t pos_col_idx;
  char id_delim;
} plink1_dosage_info_t;

typedef struct gendummy_info_struct {
  gendummy_flags_t flags;
  uint32_t sample_ct;
  uint32_t variant_ct;
  uint32_t pheno_ct;
  double geno_mfreq;
  double pheno_mfreq;
  double dosage_freq;
} gendummy_info_t;

void init_plink1_dosage(plink1_dosage_info_t* plink1_dosage_info_ptr);

void init_gendummy(gendummy_info_t* gendummy_info_ptr);

pglerr_t vcf_to_pgen(const char* vcfname, const char* preexisting_psamname, const char* const_fid, const char* dosage_import_field, misc_flags_t misc_flags, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, double import_dosage_certainty, char id_delim, char idspace_to, int32_t vcf_min_gq, int32_t vcf_min_dp, vcf_half_call_t vcf_half_call, fam_col_t fam_cols, char* outname, char* outname_end, chr_info_t* cip);

pglerr_t ox_gen_to_pgen(const char* genname, const char* samplename, const char* ox_single_chr_str, const char* ox_missing_code, misc_flags_t misc_flags, oxford_import_t oxford_import_flags, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, double import_dosage_certainty, char* outname, char* outname_end, chr_info_t* cip);

pglerr_t ox_bgen_to_pgen(const char* bgenname, const char* samplename, const char* const_fid, const char* ox_missing_code, misc_flags_t misc_flags, oxford_import_t oxford_import_flags, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, double import_dosage_certainty, char id_delim, char idspace_to, uint32_t max_thread_ct, char* outname, char* outname_end, chr_info_t* cip);

pglerr_t ox_hapslegend_to_pgen(const char* hapsname, const char* legendname, const char* samplename, const char* ox_single_chr_str, const char* ox_missing_code, misc_flags_t misc_flags, oxford_import_t oxford_import_flags, char* outname, char* outname_end, chr_info_t* cip);

pglerr_t plink1_dosage_to_pgen(const char* dosagename, const char* famname, const char* mapname, const char* import_single_chr_str, const plink1_dosage_info_t* pdip, misc_flags_t misc_flags, fam_col_t fam_cols, int32_t missing_pheno, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, double import_dosage_certainty, uint32_t max_thread_ct, char* outname, char* outname_end, chr_info_t* cip);

pglerr_t generate_dummy(const gendummy_info_t* gendummy_info_ptr, misc_flags_t misc_flags, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, uint32_t max_thread_ct, char* outname, char* outname_end, chr_info_t* cip);

pglerr_t plink1_sample_major_to_pgen(const char* pgenname, uintptr_t variant_ct, uintptr_t sample_ct, uint32_t real_ref_alleles, uint32_t max_thread_ct, FILE* infile);

#ifdef __cplusplus
} // namespace plink2
#endif

#endif // __PLINK2_IMPORT_H__
