#ifndef __PLINK2_PVAR_H__
#define __PLINK2_PVAR_H__

// This library is part of PLINK 2.00, copyright (C) 2005-2020 Shaun Purcell,
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

#include "plink2_common.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// Key .pvar properties: (i) .bim files are valid .pvar files; (ii) .vcf files,
// with or without post-INFO columns removed, are valid .pvar files; (iii)
// these files remain valid when {b}gzipped.
//
// File starts with an arbitrary (possibly zero) number of header lines marked
// by a leading '#'.  All lines which don't start with '#CHROM' are currently
// ignored.  The #CHROM line specifies the columns in the .pvar file; the
// following column headers are recognized:
//   POS (bp coordinate)
//   ID (variant ID)
//   REF (reference allele)
//   ALT (alternate alleles, comma-separated)
//   QUAL (phred-scaled quality score for whether the locus is variable at all)
//   FILTER (PASS, ., or semicolon-separated list of failing filter codes)
//   INFO (semicolon-separated list of flags and key-value pairs, with types
//     declared in header)
//   FORMAT (terminates header line parsing, so any VCF sample data is ignored)
//   CM (centimorgan position)
// The loader will error out if multiple #CHROM lines are in the header for
// some bizarre reason.  If no #CHROM line is present, it is assumed to be
// "#CHROM ID CM POS ALT REF" for .bim compatibility, or
// "#CHROM ID POS ALT REF" in the case where the first nonheader line has
// exactly 5 columns.


// allele_idx_offsets[] is a length-(variant_ct + 1) array of reference allele
// indices in allele_storage[].  The alt1 allele index is the reference allele
// index plus 1, etc.  The number of alt alleles for (0-based) variant n is
// (allele_idx_offsets[n+1] - allele_idx_offsets[n] - 1).

// Memory for all the return arrays is allocated off the bottom of g_bigstack.
// (I.e. it's assumed that you'll want to keep these arrays until your program
// is about to terminate.)


PglErr ReadChrsetHeaderLine(const char* chrset_iter, const char* file_descrip, MiscFlags misc_flags, uintptr_t line_idx, ChrInfo* cip);

typedef struct VaridTemplateStruct {
  NONCOPYABLE(VaridTemplateStruct);
  const char* missing_id_match;
  char* chr_output_name_buf;
  STD_ARRAY_DECL(const char*, 5, segs);
  STD_ARRAY_DECL(uint32_t, 5, seg_lens);
  STD_ARRAY_DECL(uint32_t, 4, insert_types);
  uint32_t insert_ct;
  uint32_t base_len;
  uint32_t chr_slen;
  uint32_t alleles_needed;

  // next two values (and the first two) are common between templates
  uint32_t new_id_max_allele_slen;
  uint32_t overflow_substitute_blen;  // only nonzero if overflow -> missing
} VaridTemplate;

void VaridTemplateInit(const char* varid_template_str, const char* missing_id_match, char* chr_output_name_buf, uint32_t new_id_max_allele_slen, uint32_t overflow_substitute_blen, VaridTemplate* vtp);

char* VaridTemplateWrite(const VaridTemplate* vtp, const char* ref_start, const char* alt1_start, uint32_t cur_bp, uint32_t ref_token_slen, uint32_t extra_alt_ct, uint32_t alt_token_slen, char* dst);

// assumes info_token[-1] is safe to read
// may set info_token[info_slen] to \0, since it needs to use strstr()
// (todo: try memmem()... except it isn't available on 32-bit mingw?)
char* PrInInfoToken(uint32_t info_slen, char* info_token);

// cip, max_variant_id_slen, and info_reload are in/out parameters.
// Chromosome filtering is performed if cip requests it.
PglErr LoadPvar(const char* pvarname, const char* var_filter_exceptions_flattened, const char* varid_template_str, const char* varid_multi_template_str, const char* varid_multi_nonsnp_template_str, const char* missing_varid_match, const char* require_info_flattened, const char* require_no_info_flattened, const CmpExpr* extract_if_info_exprp, const CmpExpr* exclude_if_info_exprp, MiscFlags misc_flags, PvarPsamFlags pvar_psam_flags, uint32_t xheader_needed, uint32_t qualfilter_needed, float var_min_qual, uint32_t splitpar_bound1, uint32_t splitpar_bound2, uint32_t new_variant_id_max_allele_slen, uint32_t snps_only, uint32_t split_chr_ok, uint32_t filter_min_allele_ct, uint32_t filter_max_allele_ct, uint32_t max_thread_ct, ChrInfo* cip, uint32_t* max_variant_id_slen_ptr, uint32_t* info_reload_slen_ptr, UnsortedVar* vpos_sortstatus_ptr, char** xheader_ptr, uintptr_t** variant_include_ptr, uint32_t** variant_bps_ptr, char*** variant_ids_ptr, uintptr_t** allele_idx_offsets_ptr, const char*** allele_storage_ptr, uintptr_t** qual_present_ptr, float** quals_ptr, uintptr_t** filter_present_ptr, uintptr_t** filter_npass_ptr, char*** filter_storage_ptr, uintptr_t** nonref_flags_ptr, double** variant_cms_ptr, ChrIdx** chr_idxs_ptr, uint32_t* raw_variant_ct_ptr, uint32_t* variant_ct_ptr, uint32_t* max_allele_ct_ptr, uint32_t* max_allele_slen_ptr, uintptr_t* xheader_blen_ptr, InfoFlags* info_flags_ptr, uint32_t* max_filter_slen_ptr);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_PVAR_H__
