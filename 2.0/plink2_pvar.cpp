// This file is part of PLINK 2.00, copyright (C) 2005-2017 Shaun Purcell,
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
#include "plink2_decompress.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// this used to employ a backward-growing linked list, with variable-size
// elements, but demultiplexing was relatively expensive.  now we allocate
// size-64k pos[], allele_idxs[], ids[], cms[], etc. blocks, and just memcpy
// those chunks at the end.  (cms[] is lazy-initialized.)
CONSTU31(kLoadPvarBlockSize, 65536);
static_assert(!(kLoadPvarBlockSize & (kLoadPvarBlockSize - 1)), "kLoadPvarBlockSize must be a power of 2.");
static_assert(kLoadPvarBlockSize >= (kMaxMediumLine / 8), "kLoadPvarBlockSize cannot be smaller than kMaxMediumLine / 8.");


static_assert(!kChrOffsetX, "read_chrset_header_line() assumes kChrOffsetX == 0.");
static_assert(kChrOffsetY == 1, "read_chrset_header_line() assumes kChrOffsetY == 1.");
static_assert(kChrOffsetPAR1 == 4, "read_chrset_header_line() assumes kChrOffsetPAR1 == 4.");
pglerr_t read_chrset_header_line(char* chrset_iter, const char* file_descrip, misc_flags_t misc_flags, uintptr_t line_idx, chr_info_t* cip) {
  // chrset_iter is expected to point to first character after
  // "##chrSet=<".
  pglerr_t reterr = kPglRetSuccess;
  {
    uint32_t cmdline_autosome_ct = 0;
    uint32_t cmdline_haploid = 0;
    int32_t cmdline_xymt_codes[kChrOffsetCt];
    if (cip->chrset_source == kChrsetSourceCmdline) {
      if (misc_flags & kfMiscChrOverrideCmdline) {
	goto read_chrset_header_line_ret_1;
      }
      if (!(misc_flags & kfMiscChrOverrideFile)) {
	// save off info we need for consistency check
	cmdline_autosome_ct = cip->autosome_ct;
	cmdline_haploid = cip->haploid_mask[0] & 1;
	memcpy(cmdline_xymt_codes, cip->xymt_codes, kChrOffsetCt * sizeof(int32_t));
      }
      fill_ulong_zero(kChrMaskWords, cip->haploid_mask);
    }
    for (uint32_t uii = 0; uii < kChrOffsetCt; ++uii) {
      cip->xymt_codes[uii] = -2;
    }
    if (!memcmp(chrset_iter, "haploidAutosomeCt=", 18)) {
      uint32_t explicit_haploid_ct;
      if (scan_posint_capped(&(chrset_iter[18]), kMaxChrTextnum, &explicit_haploid_ct)) {
	sprintf(g_logbuf, "Error: Header line %" PRIuPTR " of %s has an invalid ##chrSet haploid count (max %u).\n", line_idx, file_descrip, kMaxChrTextnum);
	goto read_chrset_header_line_ret_MALFORMED_INPUT_WW;
      }
      // could verify that X, Y, etc. are not present?
      if (cmdline_autosome_ct) {
	if (!cmdline_haploid) {
	  sprintf(g_logbuf, "Error: Header line %" PRIuPTR " of %s specifies a haploid genome, while a diploid genome was specified on the command line.\n", line_idx, file_descrip);
	  goto read_chrset_header_line_ret_INCONSISTENT_INPUT_WW;
	}
	if (explicit_haploid_ct != cmdline_autosome_ct) {
	  sprintf(g_logbuf, "Error: Header line %" PRIuPTR " of %s specifies %u autosome%s, while the command line specified %u.\n", line_idx, file_descrip, explicit_haploid_ct, (explicit_haploid_ct == 1)? "" : "s", cmdline_autosome_ct);
	  goto read_chrset_header_line_ret_INCONSISTENT_INPUT_WW;
	}
      }
      cip->autosome_ct = explicit_haploid_ct;
      fill_all_bits(explicit_haploid_ct + 1, cip->haploid_mask);
    } else {
      if (memcmp(chrset_iter, "autosomePairCt=", 15)) {
	sprintf(g_logbuf, "Error: Header line %" PRIuPTR " of %s does not have expected ##chrSet format.\n", line_idx, file_descrip);
	goto read_chrset_header_line_ret_MALFORMED_INPUT_WW;
      }
      chrset_iter = &(chrset_iter[15]);
      uint32_t explicit_autosome_ct;
      if (scanadv_posint_capped(kMaxChrTextnum, &chrset_iter, &explicit_autosome_ct)) {
	sprintf(g_logbuf, "Error: Header line %" PRIuPTR " of %s has an invalid ##chrSet autosome count (max %u).\n", line_idx, file_descrip, kMaxChrTextnum);
	goto read_chrset_header_line_ret_MALFORMED_INPUT_WW;
      }
      cip->autosome_ct = explicit_autosome_ct;
      if (*chrset_iter != '>') {
	if (*chrset_iter != ',') {
	  sprintf(g_logbuf, "Error: Header line %" PRIuPTR " of %s does not have expected ##chrSet format.\n", line_idx, file_descrip);
	  goto read_chrset_header_line_ret_MALFORMED_INPUT_WW;
	}
	// this can theoretically be confused by e.g. a Description="..." field
	// containing commas not followed by spaces
	while (1) {
	  ++chrset_iter;
	  // uppercase
	  uint32_t first_char_ui = ((unsigned char)(*chrset_iter)) & 0xdf;

	  uint32_t second_char_ui = (unsigned char)chrset_iter[1];
	  // 44 is ',', 62 is '>'
	  if ((second_char_ui == 44) || (second_char_ui == 62)) {
	    if (first_char_ui == 77) {
	      // M
	      cip->xymt_codes[kChrOffsetMT] = explicit_autosome_ct + 1 + kChrOffsetMT;
	    } else {
	      first_char_ui -= 88; // X = 0, Y = 1, everything else larger
	      if (first_char_ui < 2) {
		cip->xymt_codes[first_char_ui] = explicit_autosome_ct + 1 + first_char_ui;
	      }
	    }
	  } else {
	    second_char_ui &= 0xdf;
	    const uint32_t third_char_ui = (unsigned char)chrset_iter[2];
	    if ((third_char_ui == 44) || (third_char_ui == 62)) {
	      if ((first_char_ui == 88) && (second_char_ui == 89)) {
		// XY
		cip->xymt_codes[kChrOffsetXY] = explicit_autosome_ct + 1 + kChrOffsetXY;
	      } else if ((first_char_ui == 77) && (second_char_ui == 84)) {
		// MT
		cip->xymt_codes[kChrOffsetMT] = explicit_autosome_ct + 1 + kChrOffsetMT;
	      }
	    } else if ((first_char_ui == 80) && (second_char_ui == 65) && ((third_char_ui & 0xdf) == 82)) {	    
	      // PAR1, PAR2
	      const uint32_t par_idx_m1 = ((unsigned char)chrset_iter[3]) - '1';
	      if ((par_idx_m1 < 2) && ((chrset_iter[4] == ',') || (chrset_iter[4] == '>'))) {
		cip->xymt_codes[kChrOffsetPAR1] = explicit_autosome_ct + 1 + kChrOffsetPAR1 + par_idx_m1;
	      }
	    }
	  }
	  chrset_iter = strchr(chrset_iter, ',');
	  if (!chrset_iter) {
	    break;
	  }
	}
      }
      if (cmdline_autosome_ct) {
	if (cmdline_haploid) {
	  sprintf(g_logbuf, "Error: Header line %" PRIuPTR " of %s specifies a diploid genome, while a haploid genome was specified on the command line.\n", line_idx, file_descrip);
	  goto read_chrset_header_line_ret_INCONSISTENT_INPUT_WW;
	}
	if (explicit_autosome_ct != cmdline_autosome_ct) {
	  sprintf(g_logbuf, "Error: Header line %" PRIuPTR " of %s specifies %u autosome%s, while the command line specified %u.\n", line_idx, file_descrip, explicit_autosome_ct, (explicit_autosome_ct == 1)? "" : "s", cmdline_autosome_ct);
	  goto read_chrset_header_line_ret_INCONSISTENT_INPUT_WW;
	}
	for (uint32_t xymt_idx = 0; xymt_idx < kChrOffsetPAR1; ++xymt_idx) {
	  // it's okay if the command line doesn't explicitly exclude e.g. chrX
	  // while for whatever reason it is excluded from ##chrSet; but the
	  // reverse can create problems
	  if ((cmdline_xymt_codes[xymt_idx] < 0) && (cip->xymt_codes[xymt_idx] >= 0)) {
	    sprintf(g_logbuf, "Error: Header line %" PRIuPTR " of %s specifies a chromosome set including %s, while the command line excludes it.\n", line_idx, file_descrip, g_xymt_log_names[xymt_idx]);
	    goto read_chrset_header_line_ret_INCONSISTENT_INPUT_WW;
	  }
	}
      }
    }
    cip->chrset_source = kChrsetSourceFile;
  }
  while (0) {
  read_chrset_header_line_ret_MALFORMED_INPUT_WW:
    wordwrapb(0);
    logerrprintb();
    reterr = kPglRetMalformedInput;
    break;
  read_chrset_header_line_ret_INCONSISTENT_INPUT_WW:
    wordwrapb(0);
    logerrprintb();
    reterr = kPglRetInconsistentInput;
    break;
  }
 read_chrset_header_line_ret_1:
  return reterr;
}

void varid_template_init(const char* varid_template, uint32_t* template_insert_ct_ptr, uint32_t* template_base_len_ptr, uint32_t* alleles_needed_ptr, const char* varid_template_segs[5], uint32_t* varid_template_seg_lens, uint32_t* varid_template_types) {
  // template string was previously validated
  // varid_template is only input, everything else is output values
  const char* varid_template_iter = varid_template;
  uint32_t template_insert_ct = 0;
  uint32_t template_base_len = 0;
  unsigned char ucc = (unsigned char)(*varid_template_iter);
  uint32_t alleles_needed = 0; // bit 0 = ref, bit 1 = alt, bit 2 = ascii sort
  varid_template_segs[0] = varid_template_iter;
  do {
    if (ucc <= '@') {
      uint32_t seg_len;
      uint32_t insert_type;
      if (ucc == '@') {
	seg_len = (uintptr_t)(varid_template_iter - varid_template_segs[template_insert_ct]);
	insert_type = 0;
	goto varid_template_init_match;
      }
      if (ucc == '#') {
	seg_len = (uintptr_t)(varid_template_iter - varid_template_segs[template_insert_ct]);
	insert_type = 1;
	goto varid_template_init_match;
      }
      if (ucc == '$') {
	seg_len = (uintptr_t)(varid_template_iter - varid_template_segs[template_insert_ct]);
	{
	  const uint32_t uii = (unsigned char)(*(++varid_template_iter));
	  if (uii <= '2') {
	    alleles_needed += 2; // this happens twice
	    insert_type = uii - 48; // '1' -> type 2, '2' -> type 3
	  } else {
	    // 'r' -> type 2, 'a' -> type 3
	    insert_type = 1 + ((uii & 0xdf) == 'A');
	  }
	  alleles_needed += insert_type;
	  ++insert_type;
	}
      varid_template_init_match:
	varid_template_seg_lens[template_insert_ct] = seg_len;
	template_base_len += seg_len;
	varid_template_types[template_insert_ct++] = insert_type;
	varid_template_segs[template_insert_ct] = &(varid_template_iter[1]);
      }
    }
    ucc = (unsigned char)(*(++varid_template_iter));
  } while (ucc);
  const uint32_t seg_len = (uintptr_t)(varid_template_iter - varid_template_segs[template_insert_ct]);
  varid_template_seg_lens[template_insert_ct] = seg_len;
  *template_insert_ct_ptr = template_insert_ct;
  *template_base_len_ptr = template_base_len + seg_len;
  *alleles_needed_ptr = alleles_needed;
}

void backfill_chr_idxs(const chr_info_t* cip, uint32_t chrs_encountered_m1, uint32_t offset, uint32_t end_vidx, chr_idx_t* chr_idxs) {
  uint32_t chr_fo_idx = chrs_encountered_m1;
  while (1) {
    uint32_t start_vidx = cip->chr_fo_vidx_start[chr_fo_idx];
    if (start_vidx < offset) {
      start_vidx = offset;
    }
    chr_idx_t* chr_idxs_write_base = &(chr_idxs[start_vidx - offset]);
    const uint32_t vidx_ct = end_vidx - start_vidx;
    const chr_idx_t cur_chr_idx = (uint32_t)(cip->chr_file_order[chr_fo_idx]);
    for (uint32_t uii = 0; uii < vidx_ct; ++uii) {
      chr_idxs_write_base[uii] = cur_chr_idx;
    }
    if (start_vidx == offset) {
      return;
    }
    end_vidx = start_vidx;
    --chr_fo_idx;
  }
}

char* pr_in_info_token(uint32_t info_slen, char* info_token) {
  if ((!memcmp(info_token, "PR", 2)) && ((info_slen == 2) || (info_token[2] == ';'))) {
    return info_token;
  }
  if (!memcmp(&(info_token[((int32_t)info_slen) - 3]), ";PR", 3)) {
    return &(info_token[info_slen - 2]);
  }
  info_token[info_slen] = '\0';
  char* first_info_end = strchr(info_token, ';');
  if (!first_info_end) {
    return nullptr;
  }
  char* pr_prestart = strstr(first_info_end, ";PR;");
  return pr_prestart? nullptr : (&(pr_prestart[1]));
}

pglerr_t splitpar(const uint32_t* variant_bps, unsorted_var_t vpos_sortstatus, uint32_t splitpar_bound1, uint32_t splitpar_bound2, uintptr_t* variant_include, uintptr_t* loaded_chr_mask, chr_info_t* cip, uint32_t* chrs_encountered_m1_ptr, uint32_t* exclude_ct_ptr) {
  const int32_t x_code = cip->xymt_codes[kChrOffsetX];
  if ((x_code < 0) || (!is_set(loaded_chr_mask, x_code))) {
    logerrprint("Warning: --split-par had no effect (no X chromosome in dataset).\n");
    return kPglRetSuccess;
  }
  const int32_t par1_code = cip->xymt_codes[kChrOffsetPAR1];
  const int32_t par2_code = cip->xymt_codes[kChrOffsetPAR2];
  if (par2_code < 0) {
    // may want to remove this restriction later
    logerrprint("Error: --split-par cannot currently be used with a custom chromosome set.\n");
    return kPglRetInvalidCmdline;
  }
  if (is_set(loaded_chr_mask, par1_code) || is_set(loaded_chr_mask, par2_code)) {
    logerrprint("Error: --split-par cannot be used on a dataset which already contains a PAR1 or\nPAR2 region.\n");
    return kPglRetInvalidCmdline;
  }
  if (vpos_sortstatus & kfUnsortedVarBp) {
    logerrprint("Error: --split-par cannot be used with an unsorted .bim/.pvar file.\n");
    return kPglRetInvalidCmdline;
  }
  const uint32_t orig_xchr_fo_idx = cip->chr_idx_to_foidx[(uint32_t)x_code];
  const uint32_t orig_x_start = cip->chr_fo_vidx_start[orig_xchr_fo_idx];
  const uint32_t orig_x_end = cip->chr_fo_vidx_start[orig_xchr_fo_idx + 1];
  const uint32_t par1_end = orig_x_start + uint32arr_greater_than(&(variant_bps[orig_x_start]), orig_x_end - orig_x_start, splitpar_bound1 + 1);
  const uint32_t par2_start = par1_end + uint32arr_greater_than(&(variant_bps[par1_end]), orig_x_end - par1_end, splitpar_bound2);
  uint32_t tot_codes_changed = (par1_end - orig_x_start) + (orig_x_end - par2_start);
  if (!tot_codes_changed) {
    logerrprint("Warning: --split-par had no effect (no X variants were in the PARs).\n");
    return kPglRetSuccess;
  }
  // one of the PARs, and/or the main chrX body, may be empty; that's not a big
  // deal
  *chrs_encountered_m1_ptr += 2;
  const uint32_t chrs_encountered_m1 = *chrs_encountered_m1_ptr;
  cip->chr_fo_vidx_start[chrs_encountered_m1 + 1] = cip->chr_fo_vidx_start[chrs_encountered_m1 - 1];
  for (uint32_t chr_fo_idx = chrs_encountered_m1 - 2; chr_fo_idx > orig_xchr_fo_idx; --chr_fo_idx) {
    cip->chr_fo_vidx_start[chr_fo_idx + 2] = cip->chr_fo_vidx_start[chr_fo_idx];
    const int32_t cur_chr_idx = cip->chr_file_order[chr_fo_idx];
    cip->chr_file_order[chr_fo_idx + 2] = cur_chr_idx;
    cip->chr_idx_to_foidx[cur_chr_idx] = chr_fo_idx + 2;
  }
  cip->chr_fo_vidx_start[orig_xchr_fo_idx + 1] = par1_end;
  cip->chr_fo_vidx_start[orig_xchr_fo_idx + 2] = par2_start;
  cip->chr_file_order[orig_xchr_fo_idx] = par1_code;
  cip->chr_file_order[orig_xchr_fo_idx + 1] = x_code;
  cip->chr_file_order[orig_xchr_fo_idx + 2] = par2_code;
  cip->chr_idx_to_foidx[(uint32_t)par1_code] = orig_xchr_fo_idx;
  cip->chr_idx_to_foidx[(uint32_t)x_code] = orig_xchr_fo_idx + 1;
  cip->chr_idx_to_foidx[(uint32_t)par2_code] = orig_xchr_fo_idx + 2;
  uintptr_t* chr_mask = cip->chr_mask;
  if (par1_end > orig_x_start) {
    if (!is_set(chr_mask, par1_code)) {
      *exclude_ct_ptr += popcount_bit_idx(variant_include, orig_x_start, par1_end);
      clear_bits_nz(orig_x_start, par1_end, variant_include);
    } else {
      set_bit(par1_code, loaded_chr_mask);
    }
  }
  if (par1_end == par2_start) {
    clear_bit(x_code, chr_mask);
  } else if (!is_set(chr_mask, x_code)) {
    clear_bit(x_code, chr_mask);
    *exclude_ct_ptr += popcount_bit_idx(variant_include, par1_end, par2_start);
    clear_bits_nz(par1_end, par2_start, variant_include);
  }
  if (par2_start < orig_x_end) {
    if (!is_set(chr_mask, par2_code)) {
      *exclude_ct_ptr += popcount_bit_idx(variant_include, par2_start, orig_x_end);
      clear_bits_nz(par2_start, orig_x_end, variant_include);
    } else {
      set_bit(par2_code, loaded_chr_mask);
    }
  }
  LOGPRINTF("--split-par: %u chromosome code%s changed.\n", tot_codes_changed, (tot_codes_changed == 1)? "" : "s");
  return kPglRetSuccess;
}

// --input-missing-genotype code set to 1 by load_pvar()
static uint8_t acgtm_bool_table[256] = {
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

static inline uint32_t is_acgtm(unsigned char ucc) {
  return (uint32_t)(acgtm_bool_table[ucc]);
}

static_assert((!(kMaxIdSlen % kCacheline)), "load_pvar() must be updated.");
pglerr_t load_pvar(const char* pvarname, char* var_filter_exceptions_flattened, const char* varid_template, const char* missing_varid_match, misc_flags_t misc_flags, pvar_psam_t pvar_psam_modifier, exportf_flags_t exportf_modifier, float var_min_qual, uint32_t splitpar_bound1, uint32_t splitpar_bound2, uint32_t new_variant_id_max_allele_slen, uint32_t snps_only, uint32_t split_chr_ok, chr_info_t* cip, uint32_t* max_variant_id_slen_ptr, uint32_t* info_reload_slen_ptr, unsorted_var_t* vpos_sortstatus_ptr, char** xheader_ptr, uintptr_t** variant_include_ptr, uint32_t** variant_bps_ptr, char*** variant_ids_ptr, uintptr_t** variant_allele_idxs_ptr, char*** allele_storage_ptr, uintptr_t** qual_present_ptr, float** quals_ptr, uintptr_t** filter_present_ptr, uintptr_t** filter_npass_ptr, char*** filter_storage_ptr, uintptr_t** nonref_flags_ptr, double** variant_cms_ptr, chr_idx_t** chr_idxs_ptr, uint32_t* raw_variant_ct_ptr, uint32_t* variant_ct_ptr, uint32_t* max_allele_slen_ptr, uintptr_t* xheader_blen_ptr, uint32_t* xheader_info_pr_ptr, uint32_t* max_filter_slen_ptr) {
  // chr_info, max_variant_id_slen, and info_reload_slen are in/out; just
  // outparameters after them.  (Due to its large size in some VCFs, INFO is
  // not kept in memory for now.  This has a speed penalty, of course; maybe
  // it's worthwhile to conditionally load it later.)

  // variant_allele_idxs currently assumed to be initialized to nullptr

  // should handle raw_variant_ct == 0 properly

  // todo: upgrade this to handle split chromosomes, unsorted chromosomes/bp
  //   coordinates, maybe skipping of allele code loading
  // probable todo: load INFO:END.  (does this allow the CNV module to be
  //   unified with the rest of the program?)  but this will probably wait
  //   until I need to analyze some sort of CNV data, and that day keeps
  //   getting postponed...
  // possible todo: require FILTER to only contain values declared in header,
  //   and modify its storage accordingly?  (pointless for now, but worthwhile
  //   to keep an eye on what typical VCF files look like.)
  
  // Workspace is used as follows:
  // |--header, allele_storage->----|--other return arrays---|--loadbuf--|-
  //                                                        1/4
  //
  // -temp-->----|----<- filter failures, variant IDs, long alleles--|
  //                                                                end
  // I.e. on successful return, both bigstack_base and bigstack_end will move.
  // This is designed to be called near the start of a program, at a time when
  // no large temporary buffer is needed.
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;

  gzFile gz_infile = nullptr;
  uintptr_t line_idx = 0;
  uint32_t max_allele_slen = 1;
  pglerr_t reterr = kPglRetSuccess;
  {
    reterr = gzopen_read_checked(pvarname, &gz_infile);
    if (reterr) {
      goto load_pvar_ret_1;
    }
    const uintptr_t initial_bigstack_size = bigstack_left();
    uintptr_t loadbuf_size = round_down_pow2(initial_bigstack_size / 4, kCacheline);
    char* loadbuf = (char*)(&(bigstack_mark[loadbuf_size]));
    if (loadbuf_size > kMaxLongLine) {
      loadbuf_size = kMaxLongLine;
    } else if (loadbuf_size <= kLoadPvarBlockSize * 2 * sizeof(intptr_t)) {
      goto load_pvar_ret_NOMEM;
    }
    loadbuf[loadbuf_size - 1] = ' ';

    char* xheader_end = ((pvar_psam_modifier & kfPvarColXheader) || (exportf_modifier & kfExportfVcf))? ((char*)bigstack_mark) : nullptr;
    uint32_t chrset_present = 0;
    uint32_t info_pr_present = 0;
    uint32_t info_nonpr_present = 0;
    char* loadbuf_first_token;
    while (1) {
      ++line_idx;
      // strangely, gzgets tends to be more than twice as fast as fgets on my
      // dev machine.  may as well support gzipped input files everywhere...
      // (update: now using zstd's zlibWrapper gzgets.  todo: verify that this
      // wrapper has negligible performance cost.)
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
	if (!gzeof(gz_infile)) {
	  goto load_pvar_ret_READ_FAIL;
	}
	loadbuf_first_token = loadbuf;
	loadbuf_first_token[0] = '\0';
	break;
      }
      if (!loadbuf[loadbuf_size - 1]) {
	if (loadbuf_size == kMaxLongLine) {
	  goto load_pvar_ret_LONG_LINE;
	}
	goto load_pvar_ret_NOMEM;
      }
      loadbuf_first_token = skip_initial_spaces(loadbuf);
      if ((loadbuf_first_token[0] != '#') || (!strcmp_se(loadbuf_first_token, "#CHROM", 6))) {
	break;
      }
      if (!memcmp(loadbuf_first_token, "##INFO=<ID=PR,Number=", 21)) {
	if (info_pr_present) {
	  sprintf(g_logbuf, "Error: Duplicate INFO:PR header line in %s.\n", pvarname);
	  goto load_pvar_ret_MALFORMED_INPUT_WW;
	}
	if (memcmp(&(loadbuf_first_token[21]), "0,Type=Flag,Description=", 24)) {
	  sprintf(g_logbuf, "Error: Header line %" PRIuPTR " of %s does not have expected INFO:PR format.\n", line_idx, pvarname);
	  goto load_pvar_ret_MALFORMED_INPUT_WW;
	}
	info_pr_present = 1;
      } else if ((!info_nonpr_present) && (!memcmp(loadbuf_first_token, "##INFO=<ID=", 11))) {
	info_nonpr_present = 1;
      }
      if (!memcmp(loadbuf_first_token, "##chrSet=<", 10)) {
	if (chrset_present) {
	  sprintf(g_logbuf, "Error: Multiple ##chrSet header lines in %s.\n", pvarname);
	  goto load_pvar_ret_MALFORMED_INPUT_WW;
	}
	chrset_present = 1;
	const uint32_t cmdline_chrset = (cip->chrset_source == kChrsetSourceCmdline) && (!(misc_flags & kfMiscChrOverrideFile));
	reterr = read_chrset_header_line(&(loadbuf_first_token[10]), pvarname, misc_flags, line_idx, cip);
	if (reterr) {
	  goto load_pvar_ret_1;
	}
	if (!cmdline_chrset) {
	  const uint32_t autosome_ct = cip->autosome_ct;
	  if (cip->haploid_mask[0] & 1) {
	    LOGPRINTF("chrSet header line: %u autosome%s (haploid).\n", autosome_ct, (autosome_ct == 1)? "" : "s");
	  } else {
	    LOGPRINTF("chrSet header line: %u autosome pair%s.\n", autosome_ct, (autosome_ct == 1)? "" : "s");
	  }
	}
      } else if (xheader_end) {
	// if the "pvar file" was actually a VCF, suppress the same lines we'd
	// suppress when importing with --vcf.
	if (memcmp(loadbuf_first_token, "##fileformat=", 13) && memcmp(loadbuf_first_token, "##fileDate=", 11) && memcmp(loadbuf_first_token, "##source=", 9) && memcmp(loadbuf_first_token, "##FORMAT=", 9)) {
	  uint32_t line_slen = strlen(loadbuf_first_token);
	  if (loadbuf_first_token[line_slen - 1] == '\n') {
	    --line_slen;
	    if (loadbuf_first_token[line_slen - 1] == '\r') {
	      --line_slen;
	    }
	  }
	  if ((uintptr_t)(loadbuf - xheader_end) < line_slen + 2) {
	    goto load_pvar_ret_NOMEM;
	  }
	  xheader_end = memcpya(xheader_end, loadbuf_first_token, line_slen);
	  append_binary_eoln(&xheader_end);
	}
      }
    }
    *xheader_info_pr_ptr = info_pr_present;
    if (xheader_end) {
      *xheader_ptr = (char*)bigstack_mark;
      *xheader_blen_ptr = (uintptr_t)(xheader_end - (*xheader_ptr));
      g_bigstack_base = (unsigned char*)round_up_pow2((uintptr_t)xheader_end, kCacheline);
    }
    finalize_chrset(misc_flags, cip);
    char** allele_storage = (char**)g_bigstack_base;
    char** allele_storage_iter = allele_storage;

    uint32_t col_skips[8];
    uint32_t col_types[8];
    uint32_t relevant_postchr_col_ct = 5;
    uint32_t alt_col_idx = 4;
    uint32_t load_qual_col = 0;
    uint32_t load_filter_col = 0;
    uint32_t info_col_present = 0;
    uint32_t cm_col_present = 0;
    if (loadbuf_first_token[0] == '#') {
      // parse header
      // [-1] = #CHROM (must be first column)
      // [0] = POS
      // [1] = ID
      // [2] = REF
      // [3] = ALT
      // [4] = QUAL
      // [5] = FILTER
      // [6] = INFO
      // [7] = CM (usually absent)
      
      // code is similar to plink 1.9 annotate() and gene_report(), but they
      // don't have a forced first column
      // might want to write plink2_common library functions for this...
      uint32_t col_idx = 0;
      char* token_end = &(loadbuf_first_token[6]);
      uint32_t found_header_bitset = 0;
      relevant_postchr_col_ct = 0;
      char* loadbuf_iter;
      while (1) {
        loadbuf_iter = skip_initial_spaces(token_end);
	if (is_eoln_kns(*loadbuf_iter)) {
	  break;
	}
	++col_idx;
	token_end = token_endnn(loadbuf_iter);
        const uint32_t token_slen = (uintptr_t)(token_end - loadbuf_iter);
	uint32_t cur_col_type;
	if (token_slen <= 3) {
	  if (token_slen == 3) {
	    if (!memcmp(loadbuf_iter, "POS", 3)) {
	      cur_col_type = 0;
	    } else if (!memcmp(loadbuf_iter, "REF", 3)) {
	      cur_col_type = 2;
	    } else if (!memcmp(loadbuf_iter, "ALT", 3)) {
	      cur_col_type = 3;
	      alt_col_idx = col_idx;
	    } else {
	      continue;
	    }
	  } else if (token_slen == 2) {
	    if (!memcmp(loadbuf_iter, "ID", 2)) {
	      cur_col_type = 1;
	    } else if (!memcmp(loadbuf_iter, "CM", 2)) {
	      cur_col_type = 7;
	      cm_col_present = 1;
	    } else {
	      continue;
	    }
	  } else {
	    continue;
	  }
	} else if ((token_slen == 4) && (!memcmp(loadbuf_iter, "QUAL", 4))) {
	  load_qual_col = 2 * ((pvar_psam_modifier & (kfPvarColMaybequal | kfPvarColQual)) || (exportf_modifier & kfExportfVcf)) + (var_min_qual != -1);
	  if (!load_qual_col) {
	    continue;
	  }
	  cur_col_type = 4;
	} else if ((token_slen == 4) && (!memcmp(loadbuf_iter, "INFO", 4))) {
	  cur_col_type = 6;
	  info_col_present = 1;
	} else if (token_slen == 6) {
	  if (!memcmp(loadbuf_iter, "FILTER", 6)) {
	    load_filter_col = 2 * ((pvar_psam_modifier & (kfPvarColMaybefilter | kfPvarColFilter)) || (exportf_modifier & kfExportfVcf)) + ((misc_flags / kfMiscExcludePvarFilterFail) & 1);
	    if (!load_filter_col) {
	      continue;
	    }
	    cur_col_type = 5;
	  } else if (!memcmp(loadbuf_iter, "FORMAT", 6)) {
	    break;
	  } else {
	    continue;
	  }
	} else {
	  continue;
	}
	const uint32_t cur_col_type_shifted = 1 << cur_col_type;
	if (found_header_bitset & cur_col_type_shifted) {
	  *token_end = '\0';
	  sprintf(g_logbuf, "Error: Duplicate column header '%s' on line %" PRIuPTR " of %s.\n", loadbuf_iter, line_idx, pvarname);
	  goto load_pvar_ret_MALFORMED_INPUT_WW;
	}
	found_header_bitset |= cur_col_type_shifted;
	col_skips[relevant_postchr_col_ct] = col_idx;
	col_types[relevant_postchr_col_ct++] = cur_col_type;
      }
      if ((found_header_bitset & 0x0f) != 0x0f) {
	sprintf(g_logbuf, "Error: Missing column header(s) on line %" PRIuPTR " of %s. (POS, ID, REF, and ALT are required.)\n", line_idx, pvarname);
	goto load_pvar_ret_MALFORMED_INPUT_WW;
      }
      for (uint32_t rpc_col_idx = relevant_postchr_col_ct - 1; rpc_col_idx; --rpc_col_idx) {
	col_skips[rpc_col_idx] -= col_skips[rpc_col_idx - 1];
      }
      loadbuf_first_token[0] = '\0'; // forces line to be skipped by main loop
    } else if (loadbuf_first_token[0]) {
      col_skips[0] = 1;
      col_skips[1] = 1;
      col_skips[2] = 1;
      col_skips[3] = 1;
      col_types[0] = 1;
      // CM column is formally optional in headerless .pvar files (and it was
      // "secretly" optional for the standard plink 1.9 standard .bim loader).
      // If the line has exactly 5 columns, assume CM is omitted.
      char* loadbuf_iter = next_token_mult(loadbuf_first_token, 4);
      if (!loadbuf_iter) {
	goto load_pvar_ret_MISSING_TOKENS;
      }
      loadbuf_iter = next_token(loadbuf_iter);
      if (!loadbuf_iter) {
        // #CHROM ID POS ALT REF
	relevant_postchr_col_ct = 4;
	col_types[1] = 0;
	col_types[2] = 3;
	col_types[3] = 2;
	alt_col_idx = 3;
      } else {
        // #CHROM ID CM POS ALT REF
        col_skips[4] = 1;
	col_types[1] = 7;
	col_types[2] = 0;
	col_types[3] = 3;
	col_types[4] = 2;
	// alt_col_idx = 4;
	cm_col_present = 1;
      }
    }
    uint32_t info_reload_slen = *info_reload_slen_ptr;
    if (!info_col_present) {
      info_pr_present = 0;
      info_reload_slen = 0;
    } else if ((!info_pr_present) && (!info_reload_slen)) {
      info_col_present = 0;
    }
    // done with header, loadbuf_first_token now points to beginning of first
    // real line.
    uint32_t max_variant_id_slen = *max_variant_id_slen_ptr;
    uint32_t chrs_encountered_m1 = 0xffffffffU; // intentional overflow
    uint32_t prev_chr_code = 0xffffffffU; // force initial mismatch
    uint32_t raw_variant_ct = 0;
    uintptr_t* chr_mask = cip->chr_mask;
    const char* missing_allele_str = &(g_one_char_strs[92]);
    double last_cm = -DBL_MAX;
    int32_t last_bp = 0;

    // this way, we only need to check allele_storage_iter against this (i)
    // when processing a multiallelic variant or (ii) at the end of a block
    char** allele_storage_limit = (char**)(&(loadbuf[kLoadPvarBlockSize * (-2) * sizeof(intptr_t)]));
    
    unsigned char* tmp_alloc_base = (unsigned char*)(&(loadbuf[loadbuf_size]));
    uintptr_t* loaded_chr_mask = (uintptr_t*)tmp_alloc_base;
    // bugfix (2 Jun 2017): forgot to zero-initialize loaded_chr_mask
    fill_ulong_zero(kChrMaskWords, loaded_chr_mask);
    tmp_alloc_base = &(tmp_alloc_base[round_up_pow2(kChrMaskWords * sizeof(intptr_t), kCacheline)]);
    unsigned char* tmp_alloc_end = bigstack_end_mark;
    uint32_t fexcept_ct = 0;
    uintptr_t max_fexcept_blen = 2;
    char* sorted_fexcepts = nullptr;
    if (var_filter_exceptions_flattened) {
      char** strptr_arr = (char**)tmp_alloc_end;
      if (count_and_measure_multistr_reverse_alloc(var_filter_exceptions_flattened, ((uintptr_t)(tmp_alloc_end - tmp_alloc_base)) / sizeof(intptr_t), &fexcept_ct, &max_fexcept_blen, &strptr_arr)) {
	goto load_pvar_ret_NOMEM;
      }
      if ((uintptr_t)(((unsigned char*)strptr_arr) - tmp_alloc_base) < fexcept_ct * max_fexcept_blen) {
	goto load_pvar_ret_NOMEM;
      }
      strptr_arr_sort(fexcept_ct, strptr_arr);
      sorted_fexcepts = (char*)tmp_alloc_base;
      fexcept_ct = copy_and_dedup_sorted_strptrs_to_strbox(strptr_arr, fexcept_ct, max_fexcept_blen, sorted_fexcepts);
      tmp_alloc_base = &(tmp_alloc_base[round_up_pow2(fexcept_ct * max_fexcept_blen, kCacheline)]);
    }
    char* chr_output_name_buf = nullptr;
    const char* varid_template_segs[5];
    uint32_t insert_slens[4];
    uint32_t varid_template_seg_lens[5];
    uint32_t varid_template_insert_types[4];
    uint32_t varid_template_insert_ct = 0;
    uint32_t varid_template_base_len = 0;
    uint32_t varid_alleles_needed = 0;
    uint32_t missing_varid_blen = 0;
    uint32_t missing_varid_match_slen = 0;
    fill_uint_zero(4, insert_slens);
    if (varid_template) {
      if ((uintptr_t)(tmp_alloc_end - tmp_alloc_base) < kMaxIdSlen) {
	goto load_pvar_ret_NOMEM;
      }
      chr_output_name_buf = (char*)tmp_alloc_base;
      tmp_alloc_base = &(tmp_alloc_base[kMaxIdSlen]);
      if (!missing_varid_match) {
	missing_varid_match = &(g_one_char_strs[92]); // '.'
      }
      missing_varid_blen = strlen(missing_varid_match);
      if (misc_flags & kfMiscSetMissingVarIds) {
	missing_varid_match_slen = missing_varid_blen;
      }
      ++missing_varid_blen;
      varid_template_init(varid_template, &varid_template_insert_ct, &varid_template_base_len, &varid_alleles_needed, varid_template_segs, varid_template_seg_lens, varid_template_insert_types);
    }

    // prevent later return-array allocations from overlapping with temporary
    // storage
    g_bigstack_end = tmp_alloc_base;

    // prevent variant_id_htable_find from breaking
    if (((const char*)tmp_alloc_end) > (&(g_one_char_strs[512 - kMaxIdSlen]))) {
      // const_cast
      tmp_alloc_end = (unsigned char*)((uintptr_t)(&(g_one_char_strs[512 - kMaxIdSlen])));
    }
    const uint32_t allow_extra_chrs = (misc_flags / kfMiscAllowExtraChrs) & 1;
    const uint32_t merge_par = (misc_flags / kfMiscMergePar) & 1;
    const int32_t x_code = cip->xymt_codes[kChrOffsetX];
    const int32_t par2_code = cip->xymt_codes[kChrOffsetPAR2];
    int32_t parx_code = cip->xymt_codes[kChrOffsetPAR1];
    uint32_t merge_par_ct = 0;

    // Corner case: with --split-par + --not-chr x, we should keep the
    // pseudoautosomal regions.  To facilitate this, we temporarily don't mask
    // out chrX; splitpar() handles this properly later.
    const uint32_t splitpar_and_exclude_x = splitpar_bound2 && (x_code >= 0) && (!is_set(cip->chr_mask, x_code));
    if (splitpar_and_exclude_x) {
      set_bit(x_code, cip->chr_mask);
    }

    if (snps_only > 1) {
      acgtm_bool_table[(unsigned char)(*g_input_missing_geno_ptr)] = 1;
    }
    
    uint32_t* cur_bps = nullptr;
    uintptr_t* cur_allele_idxs = nullptr;
    char** cur_ids = nullptr;
    uintptr_t* cur_include = nullptr;
    uintptr_t* cur_qual_present = nullptr;
    float* cur_quals = nullptr;
    uintptr_t* cur_filter_present = nullptr;
    uintptr_t* cur_filter_npass = nullptr;
    char** cur_filter_storage = nullptr;
    uintptr_t* cur_nonref_flags = nullptr;
    uint32_t max_filter_slen = 0;
    uint32_t exclude_ct = 0;

    // only allocated when necessary
    // if we want to scale this approach to more fields, we'll need to add a
    // few pointers to the start of each block.  right now, we force cur_cms[]
    // to be allocated before cur_chr_idxs[] when both are present, but this
    // is error-prone.
    uint32_t at_least_one_npass_filter = 0;
    uint32_t at_least_one_nzero_cm = 0;
    const uint32_t new_variant_id_overflow_missing = (misc_flags / kfMiscNewVarIdOverflowMissing) & 1;
    uintptr_t new_variant_id_allele_len_overflow = 0;
    double* cur_cms = nullptr;
    uint32_t cms_start_block = 0xffffffffU;

    chr_idx_t* cur_chr_idxs = nullptr;
    uint32_t chr_idxs_start_block = 0xffffffffU;
    uint32_t is_split_chr = 0;
    unsorted_var_t vpos_sortstatus = kfUnsortedVar0;
    
    while (1) {
      if (!is_eoln_kns(*loadbuf_first_token)) {
#ifdef __LP64__
	// maximum prime < 2^32 is 4294967291; quadratic hashing guarantee
	// breaks down past that divided by 2.
	if (raw_variant_ct == 0x7ffffffd) {
	  logerrprint("Error: " PROG_NAME_STR " does not support more than 2^31 - 3 variants.  We recommend other\nsoftware, such as PLINK/SEQ, for very deep studies of small numbers of genomes.\n");
	  goto load_pvar_ret_MALFORMED_INPUT;
	}
#endif
	const uint32_t variant_idx_lowbits = raw_variant_ct % kLoadPvarBlockSize;
	if (!variant_idx_lowbits) {
	  if (((uintptr_t)(tmp_alloc_end - tmp_alloc_base) <= kLoadPvarBlockSize * (sizeof(int32_t) + 2 * sizeof(intptr_t) + at_least_one_nzero_cm * sizeof(double)) + is_split_chr * sizeof(chr_idx_t) + (1 + info_pr_present) * (kLoadPvarBlockSize / CHAR_BIT) + (load_qual_col? ((kLoadPvarBlockSize / CHAR_BIT) + kLoadPvarBlockSize * sizeof(float)) : 0) + (load_filter_col? (2 * (kLoadPvarBlockSize / CHAR_BIT) + kLoadPvarBlockSize * sizeof(intptr_t)) : 0)) || (allele_storage_iter >= allele_storage_limit)) {
	    goto load_pvar_ret_NOMEM;
	  }
	  cur_bps = (uint32_t*)tmp_alloc_base;
	  cur_allele_idxs = (uintptr_t*)(&(tmp_alloc_base[kLoadPvarBlockSize * sizeof(int32_t)]));
	  cur_ids = (char**)(&(tmp_alloc_base[kLoadPvarBlockSize * (sizeof(int32_t) + sizeof(intptr_t))]));
	  cur_include = (uintptr_t*)(&(tmp_alloc_base[kLoadPvarBlockSize * (sizeof(int32_t) + 2 * sizeof(intptr_t))]));
	  fill_ulong_one(kLoadPvarBlockSize / kBitsPerWord, cur_include);
	  tmp_alloc_base = &(tmp_alloc_base[kLoadPvarBlockSize * (sizeof(int32_t) + 2 * sizeof(intptr_t)) + (kLoadPvarBlockSize / CHAR_BIT)]);
	  if (load_qual_col > 1) {
	    cur_qual_present = (uintptr_t*)tmp_alloc_base;
	    fill_ulong_zero(kLoadPvarBlockSize / kBitsPerWord, cur_qual_present);
	    cur_quals = (float*)(&(tmp_alloc_base[kLoadPvarBlockSize / CHAR_BIT]));
	    tmp_alloc_base = &(tmp_alloc_base[kLoadPvarBlockSize * sizeof(float) + (kLoadPvarBlockSize / CHAR_BIT)]);
	  }
	  if (load_filter_col > 1) {
	    cur_filter_present = (uintptr_t*)tmp_alloc_base;
	    cur_filter_npass = (uintptr_t*)(&(tmp_alloc_base[kLoadPvarBlockSize / CHAR_BIT]));
	    cur_filter_storage = (char**)(&(tmp_alloc_base[2 * (kLoadPvarBlockSize / CHAR_BIT)]));
	    fill_ulong_zero(kLoadPvarBlockSize / kBitsPerWord, cur_filter_present);
	    fill_ulong_zero(kLoadPvarBlockSize / kBitsPerWord, cur_filter_npass);
	    tmp_alloc_base = &(tmp_alloc_base[2 * (kLoadPvarBlockSize / CHAR_BIT) + kLoadPvarBlockSize * sizeof(intptr_t)]);
	  }
	  if (info_pr_present) {
	    cur_nonref_flags = (uintptr_t*)tmp_alloc_base;
	    fill_ulong_zero(kLoadPvarBlockSize / kBitsPerWord, cur_nonref_flags);
	    tmp_alloc_base = &(tmp_alloc_base[kLoadPvarBlockSize / CHAR_BIT]);
	  }
	  if (at_least_one_nzero_cm) {
	    cur_cms = (double*)tmp_alloc_base;
	    fill_double_zero(kLoadPvarBlockSize, cur_cms);
	    tmp_alloc_base = (unsigned char*)(&(cur_cms[kLoadPvarBlockSize]));
	  }
	  if (is_split_chr) {
	    cur_chr_idxs = (chr_idx_t*)tmp_alloc_base;
	    tmp_alloc_base = (unsigned char*)(&(cur_chr_idxs[kLoadPvarBlockSize]));
	  }
	}
	char* loadbuf_iter = token_endnn(loadbuf_first_token);	
	// #CHROM
	if (!(*loadbuf_iter)) {
	  goto load_pvar_ret_MISSING_TOKENS;
	}
	int32_t cur_chr_code;
	reterr = get_or_add_chr_code_destructive(".pvar file", line_idx, allow_extra_chrs, loadbuf_first_token, loadbuf_iter, cip, &cur_chr_code);
	if (reterr) {
	  goto load_pvar_ret_1;
	}
	if (merge_par) {
	  if (cur_chr_code == par2_code) {
	    // don't permit PAR1 variants after PAR2
	    parx_code = par2_code;
	  }
	  if (cur_chr_code == parx_code) {
	    ++merge_par_ct;
	    cur_chr_code = x_code;
	  }
	}
	if (((uint32_t)cur_chr_code) != prev_chr_code) {
	  prev_chr_code = cur_chr_code;
	  if (!is_split_chr) {
	    if (is_set(loaded_chr_mask, cur_chr_code)) {
	      if (!split_chr_ok) {
		sprintf(g_logbuf, "Error: %s has a split chromosome. Use --make-pgen by itself to remedy this.\n", pvarname);
		goto load_pvar_ret_MALFORMED_INPUT_WW;
	      }
	      if ((uintptr_t)(tmp_alloc_end - tmp_alloc_base) < kLoadPvarBlockSize * sizeof(chr_idx_t)) {
		goto load_pvar_ret_NOMEM;
	      }
	      cur_chr_idxs = (chr_idx_t*)tmp_alloc_base;
	      tmp_alloc_base = (unsigned char*)(&(cur_chr_idxs[kLoadPvarBlockSize]));
	      // may want to track the first problem variant index
	      // cip->chr_fo_vidx_start[chrs_encountered_m1] = raw_variant_ct;
	      backfill_chr_idxs(cip, chrs_encountered_m1, round_down_pow2(raw_variant_ct, kLoadPvarBlockSize), raw_variant_ct, cur_chr_idxs);
	      chr_idxs_start_block = raw_variant_ct / kLoadPvarBlockSize;
	      is_split_chr = 1;
	      vpos_sortstatus |= kfUnsortedVarBp | kfUnsortedVarCm | kfUnsortedVarSplitChr;
	    } else {
	      // how much of this do we need in split-chrom case?
	      cip->chr_file_order[++chrs_encountered_m1] = cur_chr_code;
	      cip->chr_fo_vidx_start[chrs_encountered_m1] = raw_variant_ct;
	      cip->chr_idx_to_foidx[(uint32_t)cur_chr_code] = chrs_encountered_m1;
	      last_bp = 0;
	      last_cm = -DBL_MAX;
	    }
	  }
	  set_bit(cur_chr_code, loaded_chr_mask);
	  if (chr_output_name_buf) {
	    varid_template_base_len -= insert_slens[0];
	    char* chr_name_end = chr_name_write(cip, (uint32_t)cur_chr_code, chr_output_name_buf);
	    insert_slens[0] = (uintptr_t)(chr_name_end - chr_output_name_buf);
	    varid_template_base_len += insert_slens[0];
	  }
	}
	*loadbuf_iter = '\t';

	// could make this store (and cur_allele_idxs[] allocation) conditional
	// on a multiallelic variant being sighted, but unlike the CM column
	// this should become common
	cur_allele_idxs[variant_idx_lowbits] = (uintptr_t)(allele_storage_iter - allele_storage);

        char* token_ptrs[8];
	uint32_t token_slens[8];	
	if (is_set(chr_mask, cur_chr_code) || info_pr_present) {
	  for (uint32_t rpc_col_idx = 0; rpc_col_idx < relevant_postchr_col_ct; ++rpc_col_idx) {
	    const uint32_t cur_col_type = col_types[rpc_col_idx];
	    loadbuf_iter = next_token_mult(loadbuf_iter, col_skips[rpc_col_idx]);
	    if (!loadbuf_iter) {
	      goto load_pvar_ret_MISSING_TOKENS;
	    }
	    token_ptrs[cur_col_type] = loadbuf_iter;
	    char* token_end = token_endnn(loadbuf_iter);
	    token_slens[cur_col_type] = (uintptr_t)(token_end - loadbuf_iter);
	    loadbuf_iter = token_end;
	  }
	  if (info_col_present) {
	    const uint32_t info_slen = token_slens[6];
	    if (info_slen > info_reload_slen) {
	      info_reload_slen = info_slen;
	    }
	    if (info_pr_present) {
	      // always load all nonref_flags entries so they can be compared
	      // against .pgen for now.

	      // (todo: general INFO filtering code)
	      char* info_token = token_ptrs[6];
	      if (((!memcmp(info_token, "PR", 2)) && ((info_slen == 2) || (info_token[2] == ';'))) || (!memcmp(&(info_token[((int32_t)info_slen) - 3]), ";PR", 3))) {
		SET_BIT(variant_idx_lowbits, cur_nonref_flags);
	      } else {
		info_token[info_slen] = '\0';
		char* first_info_end = strchr(info_token, ';');
		if (first_info_end && strstr(first_info_end, ";PR;")) {
		  SET_BIT(variant_idx_lowbits, cur_nonref_flags);
		}
	      }
	      if (!is_set(chr_mask, cur_chr_code)) {
		goto load_pvar_skip_variant;
	      }
	    }
	  }
	  
	  // POS
	  int32_t cur_bp;
	  if (scan_int_abs_defcap(token_ptrs[0], &cur_bp)) {
	    sprintf(g_logbuf, "Error: Invalid bp coordinate on line %" PRIuPTR " of %s.\n", line_idx, pvarname);
	    goto load_pvar_ret_MALFORMED_INPUT_WW;
	  }

	  if (cur_bp < 0) {
	    goto load_pvar_skip_variant;
	  }
	  
	  // QUAL
	  if (load_qual_col) {
	    char* qual_token = token_ptrs[4];
	    if ((qual_token[0] != '.') || (qual_token[1] > ' ')) {
	      float cur_qual;
	      if (scan_float(qual_token, &cur_qual)) {
		sprintf(g_logbuf, "Error: Invalid QUAL value on line %" PRIuPTR " of %s.\n", line_idx, pvarname);
		goto load_pvar_ret_MALFORMED_INPUT_WW;
	      }
	      if ((load_qual_col & 1) && (cur_qual < var_min_qual)) {
		goto load_pvar_skip_variant;
	      }
	      if (load_qual_col > 1) {
	        SET_BIT(variant_idx_lowbits, cur_qual_present);
		// possible todo: optimize all-quals-same case
		// possible todo: conditionally allocate, like cur_cms
		cur_quals[variant_idx_lowbits] = cur_qual;
	      }
	    } else if (load_qual_col & 1) {
	      goto load_pvar_skip_variant;
	    }
	  }

	  // avoid repeating the ALT string split in --set-...-var-ids case
	  loadbuf_iter = token_ptrs[3];
	  uint32_t remaining_alt_char_ct = token_slens[3];
	  // handle --snps-only here instead of later, since it reduces the
	  // amount of data we need to load
	  if (snps_only) {
	    if ((token_slens[2] != 1) || (!(remaining_alt_char_ct % 2))) {
	      goto load_pvar_skip_variant;
	    }
	    const uint32_t extra_alt_ct = remaining_alt_char_ct / 2;
	    for (uint32_t uii = 0; uii < extra_alt_ct; ++uii) {
	      // no need to check for empty allele code here, that'll be
	      // caught later
	      if (loadbuf_iter[2 * uii + 1] != ',') {
		goto load_pvar_skip_variant;
	      }
	    }
	    if (snps_only > 1) {
	      // just-acgt
	      if (!is_acgtm(token_ptrs[2][0])) {
		goto load_pvar_skip_variant;
	      }
	      for (uint32_t uii = 0; uii <= extra_alt_ct; ++uii) {
		if (!is_acgtm(loadbuf_iter[2 * uii])) {
		  goto load_pvar_skip_variant;
		}
	      }
	    }
	  }
	  
	  // FILTER
	  if (load_filter_col) {
	    char* filter_token = token_ptrs[5];
	    const uint32_t filter_slen = token_slens[5];
	    if ((filter_slen > 1) || (filter_token[0] != '.')) {
	      if ((filter_slen != 4) || memcmp(filter_token, "PASS", 4)) {
		if (load_filter_col & 1) {
		  if (!fexcept_ct) {
		    goto load_pvar_skip_variant;
		  }
		  char* filter_token_iter = filter_token;
		  uint32_t remaining_byte_ct = filter_slen;
		  while (1) {
		    char* cur_filter_name_end = (char*)memchr(filter_token_iter, ';', remaining_byte_ct);
		    uint32_t cur_slen = remaining_byte_ct;
		    if (cur_filter_name_end) {
		      cur_slen = (uintptr_t)(cur_filter_name_end - filter_token_iter);
		    }
		    // possible todo: error out on "PASS", since that
		    // shouldn't coexist with other filters
		    // possible todo: maintain a dictionary of FILTER
		    // strings, analogous to what BCF2 does on disk
		    if (bsearch_str(filter_token_iter, sorted_fexcepts, cur_slen, max_fexcept_blen, fexcept_ct) == -1) {
		      goto load_pvar_skip_variant;
		    }
		    const uint32_t cur_blen = cur_slen + 1;
		    if (cur_blen >= remaining_byte_ct) {
		      break;
		    }
		    filter_token_iter = &(filter_token_iter[cur_blen]);
		    remaining_byte_ct -= cur_blen;
		  }
		}
		if (load_filter_col > 1) {
		  SET_BIT(variant_idx_lowbits, cur_filter_npass);
		  at_least_one_npass_filter = 1;
		  // possible todo: detect repeated filter values, store more
		  // compactly
		  if (filter_slen > max_filter_slen) {
		    max_filter_slen = filter_slen;
		  }
		  tmp_alloc_end -= filter_slen + 1;
		  if (tmp_alloc_end < tmp_alloc_base) {
		    goto load_pvar_ret_NOMEM;
		  }
		  cur_filter_storage[variant_idx_lowbits] = (char*)tmp_alloc_end;
		  memcpyx(tmp_alloc_end, filter_token, filter_slen, '\0');
		}
	      }
	      if (load_filter_col > 1) {
		SET_BIT(variant_idx_lowbits, cur_filter_present);
	      }
	    }
	  }

	  if (cur_chr_idxs) {
	    cur_chr_idxs[variant_idx_lowbits] = (uint32_t)cur_chr_code;
	  }
	  if (cur_bp < last_bp) {
	    vpos_sortstatus |= kfUnsortedVarBp;
	  }
	  cur_bps[variant_idx_lowbits] = cur_bp;
	  last_bp = cur_bp;
	  char* alt_allele_iter = (char*)memchr(loadbuf_iter, ',', remaining_alt_char_ct);
	  uint32_t id_slen;
	  if ((!varid_template) || (missing_varid_match_slen && ((token_slens[1] != missing_varid_match_slen) || memcmp(token_ptrs[1], missing_varid_match, missing_varid_match_slen)))) {
	    id_slen = token_slens[1];
	    tmp_alloc_end -= id_slen + 1;
	    if (tmp_alloc_end < tmp_alloc_base) {
	      goto load_pvar_ret_NOMEM;
	    }
	    memcpyx(tmp_alloc_end, token_ptrs[1], id_slen, '\0');
	  } else {
	    insert_slens[1] = int_slen(cur_bp);
	    uint32_t ref_slen = 0;
	    uint32_t cur_overflow = 0;
	    char* tmp_allele_ptrs[2];
	    if (varid_alleles_needed & 1) {
	      ref_slen = token_slens[2];
	      if (ref_slen > new_variant_id_max_allele_slen) {
		ref_slen = new_variant_id_max_allele_slen;
		cur_overflow = 1;
	      }
	      insert_slens[2] = ref_slen;
	      tmp_allele_ptrs[0] = token_ptrs[2];
	    }
	    if (varid_alleles_needed > 1) {
	      uint32_t alt1_slen;
	      if (!alt_allele_iter) {
		alt1_slen = remaining_alt_char_ct;
	      } else {
		alt1_slen = (uintptr_t)(alt_allele_iter - loadbuf_iter);
	      }
	      if (alt1_slen > new_variant_id_max_allele_slen) {
		alt1_slen = new_variant_id_max_allele_slen;
		++cur_overflow;
	      }
	      if (varid_alleles_needed <= 3) {
	      load_pvar_keep_allele_ascii_order:
		insert_slens[3] = alt1_slen;
		tmp_allele_ptrs[1] = loadbuf_iter;
	      } else {
		uint32_t smaller_slen = alt1_slen;
		const int32_t ref_slen_geq = (ref_slen >= alt1_slen);
		if (!ref_slen_geq) {
		  smaller_slen = ref_slen;
		}
		int32_t memcmp_result = memcmp(token_ptrs[2], loadbuf_iter, smaller_slen);
		if (!memcmp_result) {
		  memcmp_result = ref_slen_geq;
		}
		if (memcmp_result <= 0) {
		  goto load_pvar_keep_allele_ascii_order;
		}
		insert_slens[3] = ref_slen;
		tmp_allele_ptrs[1] = tmp_allele_ptrs[0];
		insert_slens[2] = alt1_slen;
		tmp_allele_ptrs[0] = loadbuf_iter;
	      }
	    }
	    id_slen = varid_template_base_len + insert_slens[1] + insert_slens[2] + insert_slens[3];
	    if (new_variant_id_overflow_missing && cur_overflow) {
	      tmp_alloc_end -= missing_varid_blen;
	      if (tmp_alloc_end < tmp_alloc_base) {
		goto load_pvar_ret_NOMEM;
	      }
	      memcpy(tmp_alloc_end, missing_varid_match, missing_varid_blen);
	      id_slen = 0;
	      cur_overflow = 1;
	    } else {
	      tmp_alloc_end -= id_slen + 1;
	      if (tmp_alloc_end < tmp_alloc_base) {
		goto load_pvar_ret_NOMEM;
	      }
	      char* id_iter = (char*)tmp_alloc_end;
	      char* insert_ptrs[4];
	      for (uint32_t insert_idx = 0; insert_idx < varid_template_insert_ct; ++insert_idx) {
		id_iter = memcpya(id_iter, varid_template_segs[insert_idx], varid_template_seg_lens[insert_idx]);
		const uint32_t cur_insert_type = varid_template_insert_types[insert_idx];
		insert_ptrs[cur_insert_type] = id_iter;
		id_iter = &(id_iter[insert_slens[cur_insert_type]]);
	      }
	      memcpyx(id_iter, varid_template_segs[varid_template_insert_ct], varid_template_seg_lens[varid_template_insert_ct], '\0');

	      memcpy(insert_ptrs[0], chr_output_name_buf, insert_slens[0]);
	      uint32toa(cur_bp, insert_ptrs[1]);
	      for (uint32_t insert_type_idx = 2; insert_type_idx < varid_template_insert_ct; ++insert_type_idx) {
		memcpy(insert_ptrs[insert_type_idx], tmp_allele_ptrs[insert_type_idx - 2], insert_slens[insert_type_idx]);
	      }
	    }
            new_variant_id_allele_len_overflow += cur_overflow;
	  }
	  if (id_slen > max_variant_id_slen) {
	    max_variant_id_slen = id_slen;
	  }
	  cur_ids[variant_idx_lowbits] = (char*)tmp_alloc_end;

	  // REF
	  char* ref_allele = token_ptrs[2];
	  const uint32_t ref_slen = token_slens[2];
	  if (ref_slen == 1) {
	    // const_cast
	    *allele_storage_iter = (char*)((uintptr_t)(&(g_one_char_strs[2 * ref_allele[0]])));
	  } else {
	    tmp_alloc_end -= ref_slen + 1;
	    if (tmp_alloc_end < tmp_alloc_base) {
	      goto load_pvar_ret_NOMEM;
	    }
	    memcpyx(tmp_alloc_end, ref_allele, ref_slen, '\0');
	    *allele_storage_iter = (char*)tmp_alloc_end;
	    if (ref_slen > max_allele_slen) {
	      max_allele_slen = ref_slen;
	    }
	  }
	  ++allele_storage_iter;

	  // ALT
	  if (alt_allele_iter) {
	    do {
	      if (allele_storage_iter >= allele_storage_limit) {
		goto load_pvar_ret_NOMEM;
	      }
	      const uint32_t cur_allele_slen = (uintptr_t)(alt_allele_iter - loadbuf_iter);
	      // possible todo: convert '0' to '.' here?
	      if (cur_allele_slen == 1) {
		// const_cast
		*allele_storage_iter = (char*)((uintptr_t)(&(g_one_char_strs[2 * loadbuf_iter[0]])));
	      } else {
		if (!cur_allele_slen) {
		  goto load_pvar_ret_EMPTY_ALLELE_CODE;
		}
		tmp_alloc_end -= cur_allele_slen + 1;
		if (tmp_alloc_end < tmp_alloc_base) {
		  goto load_pvar_ret_NOMEM;
		}
		memcpyx(tmp_alloc_end, loadbuf_iter, cur_allele_slen, '\0');
		*allele_storage_iter = (char*)tmp_alloc_end;
		if (cur_allele_slen > max_allele_slen) {
		  max_allele_slen = cur_allele_slen;
		}
	      }
	      ++allele_storage_iter;
	      remaining_alt_char_ct -= cur_allele_slen + 1;
	      loadbuf_iter = alt_allele_iter;
	      alt_allele_iter = (char*)memchr(loadbuf_iter, ',', remaining_alt_char_ct);
	    } while (alt_allele_iter);
	    if (!remaining_alt_char_ct) {
	      goto load_pvar_ret_EMPTY_ALLELE_CODE;
	    }
	  }
	  if (remaining_alt_char_ct == 1) {
	    // const_cast
	    *allele_storage_iter = (char*)((uintptr_t)(&(g_one_char_strs[2 * loadbuf_iter[0]])));
	  } else {
	    tmp_alloc_end -= remaining_alt_char_ct + 1;
	    if (tmp_alloc_end < tmp_alloc_base) {
	      goto load_pvar_ret_NOMEM;
	    }
	    memcpyx(tmp_alloc_end, loadbuf_iter, remaining_alt_char_ct, '\0');
	    *allele_storage_iter = (char*)tmp_alloc_end;
	    if (remaining_alt_char_ct > max_allele_slen) {
	      max_allele_slen = remaining_alt_char_ct;
	    }
	  }
	  ++allele_storage_iter;

	  // CM
	  if (cm_col_present) {
	    char* cm_token = token_ptrs[7];
	    if ((cm_token[0] != '0') || (cm_token[1] > ' ')) {
	      double cur_cm;
	      if (!scanadv_double(cm_token, &cur_cm)) {
		sprintf(g_logbuf, "Error: Invalid centimorgan position on line %" PRIuPTR " of %s.\n", line_idx, pvarname);
		goto load_pvar_ret_MALFORMED_INPUT_WW;
	      }
	      if (cur_cm < last_cm) {
		vpos_sortstatus |= kfUnsortedVarCm;
	      } else {
		last_cm = cur_cm;
	      }
	      if (cur_cm != 0.0) {
		if (!at_least_one_nzero_cm) {
		  if ((uintptr_t)(tmp_alloc_end - tmp_alloc_base) < kLoadPvarBlockSize * sizeof(double)) {
		    goto load_pvar_ret_NOMEM;
		  }
		  if (cur_chr_idxs) {
		    // reposition cur_chr_idxs[] after cur_cms[]
		    cur_cms = (double*)cur_chr_idxs;
		    cur_chr_idxs = (chr_idx_t*)(&(cur_cms[kLoadPvarBlockSize]));
		    memcpy(cur_chr_idxs, cur_cms, kLoadPvarBlockSize * sizeof(chr_idx_t));
		    tmp_alloc_base = (unsigned char*)(&(cur_chr_idxs[kLoadPvarBlockSize]));
		  } else {
		    cur_cms = (double*)tmp_alloc_base;
		    tmp_alloc_base = (unsigned char*)(&(cur_cms[kLoadPvarBlockSize]));
		  }
		  fill_double_zero(kLoadPvarBlockSize, cur_cms);
		  cms_start_block = raw_variant_ct / kLoadPvarBlockSize;
		  at_least_one_nzero_cm = 1;
		}
	        cur_cms[variant_idx_lowbits] = cur_cm;
	      }
	    }
	  }
	} else {
	  token_ptrs[3] = next_token_mult(loadbuf_iter, alt_col_idx);
	  if (!token_ptrs[3]) {
	    goto load_pvar_ret_MISSING_TOKENS;
	  }
	  token_slens[3] = strlen_se(token_ptrs[3]);
	load_pvar_skip_variant:
	  ++exclude_ct;
	  clear_bit(variant_idx_lowbits, cur_include);
	  cur_bps[variant_idx_lowbits] = last_bp;
	  // necessary to check alt allele count
	  // const_cast
	  *allele_storage_iter++ = (char*)((uintptr_t)missing_allele_str);
	  *allele_storage_iter++ = (char*)((uintptr_t)missing_allele_str);
	  loadbuf_iter = token_ptrs[3];
	  char* token_end = &(loadbuf_iter[token_slens[3]]);
	  while (1) {
	    loadbuf_iter = (char*)memchr(loadbuf_iter, ',', (uintptr_t)(token_end - loadbuf_iter));
	    if (!loadbuf_iter) {
	      break;
	    }
	    if (allele_storage_iter >= allele_storage_limit) {
	      goto load_pvar_ret_NOMEM;
	    }
	    ++loadbuf_iter;
	    // const_cast
	    *allele_storage_iter++ = (char*)((uintptr_t)missing_allele_str);
	  }
	}
	++raw_variant_ct;
      }
      ++line_idx;
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
	if (!gzeof(gz_infile)) {
	  goto load_pvar_ret_READ_FAIL;
	}
	break;
      }
      if (!loadbuf[loadbuf_size - 1]) {
	if (loadbuf_size == kMaxLongLine) {
	  goto load_pvar_ret_LONG_LINE;
	}
	goto load_pvar_ret_NOMEM;
      }
      loadbuf_first_token = skip_initial_spaces(loadbuf);
      if (loadbuf_first_token[0] == '#') {
	sprintf(g_logbuf, "Error: Line %" PRIuPTR " of %s starts with a '#'. (This is only permitted before the first nonheader line, and if a #CHROM header line is present it must denote the end of the header block.)\n", line_idx, pvarname);
	goto load_pvar_ret_MALFORMED_INPUT_WW;
      }
    }
    if (max_variant_id_slen > kMaxIdSlen) {
      logerrprint("Error: Variant names are limited to " MAX_ID_SLEN_STR " characters.\n");
      goto load_pvar_ret_MALFORMED_INPUT;
    }
    if (new_variant_id_allele_len_overflow) {
      if (new_variant_id_overflow_missing) {
	LOGERRPRINTFWW("Warning: %" PRIuPTR " variant ID%s %s due to allele code length.\n", new_variant_id_allele_len_overflow, (new_variant_id_allele_len_overflow == 1)? "" : "s", missing_varid_match_slen? "unchanged by --set-missing-var-ids" : "erased by --set-all-var-ids");
	if (max_variant_id_slen < missing_varid_blen - 1) {
	  max_variant_id_slen = missing_varid_blen - 1;
	}
      } else if (misc_flags & kfMiscNewVarIdOverflowTruncate) {
	LOGERRPRINTF("Warning: %" PRIuPTR " allele code%s truncated by --set-%s-var-ids.\n", new_variant_id_allele_len_overflow, (new_variant_id_allele_len_overflow == 1)? "" : "s", missing_varid_match_slen? "missing" : "all");
      } else {
	LOGERRPRINTFWW("Error: %" PRIuPTR " allele code%s too long for --set-%s-var-ids. You should either switch to a different allele/variant naming scheme for long indels, or use --new-id-max-allele-len to raise the length limit.\n", new_variant_id_allele_len_overflow, (new_variant_id_allele_len_overflow == 1)? "" : "s", missing_varid_match_slen? "missing" : "all");
	goto load_pvar_ret_INCONSISTENT_INPUT;
      }
    }
    if (gzclose_null(&gz_infile)) {
      goto load_pvar_ret_READ_FAIL;
    }
    *max_variant_id_slen_ptr = max_variant_id_slen;
    *max_allele_slen_ptr = max_allele_slen;
    *max_filter_slen_ptr = max_filter_slen;
    *raw_variant_ct_ptr = raw_variant_ct;
    uintptr_t allele_idx_end = (uintptr_t)(allele_storage_iter - allele_storage);
    bigstack_finalize_ul((uintptr_t*)allele_storage, allele_idx_end);
    uintptr_t* variant_allele_idxs = nullptr;
    const uint32_t full_block_ct = raw_variant_ct / kLoadPvarBlockSize;
    const uintptr_t raw_variant_ct_lowbits = raw_variant_ct % kLoadPvarBlockSize;
    // todo: determine whether we want variant_include to be guaranteed to be
    // terminated by a zero bit
    const uint32_t raw_variant_ctl = BITCT_TO_WORDCT(raw_variant_ct);
    if (bigstack_alloc_ul(raw_variant_ctl, variant_include_ptr) ||
	bigstack_alloc_ui(raw_variant_ct, variant_bps_ptr) ||
	bigstack_alloc_cp(raw_variant_ct, variant_ids_ptr)) {
      goto load_pvar_ret_NOMEM;
    }
    uintptr_t* qual_present = nullptr;
    float* quals = nullptr;
    if (load_qual_col > 1) {
      if (bigstack_alloc_ul(raw_variant_ctl, qual_present_ptr) ||
	  bigstack_alloc_f(raw_variant_ct, quals_ptr)) {
	goto load_pvar_ret_NOMEM;
      }
      qual_present = *qual_present_ptr;
      quals = *quals_ptr;
    }
    uintptr_t* filter_present = nullptr;
    uintptr_t* filter_npass = nullptr;
    char** filter_storage = nullptr;
    if (load_filter_col > 1) {
      if (bigstack_alloc_ul(raw_variant_ctl, filter_present_ptr) ||
	  bigstack_alloc_ul(raw_variant_ctl, filter_npass_ptr)) {
	goto load_pvar_ret_NOMEM;
      }
      filter_present = *filter_present_ptr;
      filter_npass = *filter_npass_ptr;
      if (at_least_one_npass_filter) {
	// possible todo: store this in a sparse manner
	if (bigstack_alloc_cp(raw_variant_ct, filter_storage_ptr)) {
	  goto load_pvar_ret_NOMEM;
	}
	filter_storage = *filter_storage_ptr;
      }
    }
    uintptr_t* nonref_flags = nullptr;
    if (info_pr_present) {
      if (bigstack_alloc_ul(raw_variant_ctl, nonref_flags_ptr)) {
	goto load_pvar_ret_NOMEM;
      }
      nonref_flags = *nonref_flags_ptr;
    }
    // load_qual_col > 1:
    //   kLoadPvarBlockSize / CHAR_BIT for qual_present
    //   kLoadPvarBlockSize * sizeof(float) for quals
    // load_filter_col > 1:
    //   2 * (kLoadPvarBlockSize / CHAR_BIT) for filter_present, filter_npass
    //   kLoadPvarBlockSize * sizeof(intptr_t) for filter_storage
    // at_least_one_nzero_cm:
    //   kLoadPvarBlockSize * sizeof(double)
    // is_split_chr:
    //   kLoadPvarBlockSize * sizeof(chr_idx_t)
    unsigned char* read_iter = g_bigstack_end;
    uint32_t* variant_bps = *variant_bps_ptr;
    char** variant_ids = *variant_ids_ptr;
    uintptr_t* variant_include = *variant_include_ptr;
    for (uint32_t block_idx = 0; block_idx < full_block_ct; ++block_idx) {
      memcpy(&(variant_bps[block_idx * kLoadPvarBlockSize]), read_iter, kLoadPvarBlockSize * sizeof(int32_t));
      // skip over variant_allele_idxs
      read_iter = &(read_iter[kLoadPvarBlockSize * (sizeof(int32_t) + sizeof(intptr_t))]);
      memcpy(&(variant_ids[block_idx * kLoadPvarBlockSize]), read_iter, kLoadPvarBlockSize * sizeof(intptr_t));
      read_iter = &(read_iter[kLoadPvarBlockSize * sizeof(intptr_t)]);
      memcpy(&(variant_include[block_idx * (kLoadPvarBlockSize / kBitsPerWord)]), read_iter, kLoadPvarBlockSize / CHAR_BIT);
      read_iter = &(read_iter[kLoadPvarBlockSize / CHAR_BIT]);
      if (qual_present) {
	memcpy(&(qual_present[block_idx * (kLoadPvarBlockSize / kBitsPerWord)]), read_iter, kLoadPvarBlockSize / CHAR_BIT);
	read_iter = &(read_iter[kLoadPvarBlockSize / CHAR_BIT]);
	memcpy(&(quals[block_idx * kLoadPvarBlockSize]), read_iter, kLoadPvarBlockSize * sizeof(float));
	read_iter = &(read_iter[kLoadPvarBlockSize * sizeof(float)]);
      }
      if (filter_present) {
	memcpy(&(filter_present[block_idx * (kLoadPvarBlockSize / kBitsPerWord)]), read_iter, kLoadPvarBlockSize / CHAR_BIT);
	read_iter = &(read_iter[kLoadPvarBlockSize / CHAR_BIT]);
	memcpy(&(filter_npass[block_idx * (kLoadPvarBlockSize / kBitsPerWord)]), read_iter, kLoadPvarBlockSize / CHAR_BIT);
	read_iter = &(read_iter[kLoadPvarBlockSize / CHAR_BIT]);
	if (filter_storage) {
	  memcpy(&(filter_storage[block_idx * kLoadPvarBlockSize]), read_iter, kLoadPvarBlockSize * sizeof(intptr_t));
	}
	read_iter = &(read_iter[kLoadPvarBlockSize * sizeof(intptr_t)]);
      }
      if (info_pr_present) {
	memcpy(&(nonref_flags[block_idx * (kLoadPvarBlockSize / kBitsPerWord)]), read_iter, kLoadPvarBlockSize / CHAR_BIT);
	read_iter = &(read_iter[kLoadPvarBlockSize / CHAR_BIT]);
      }
      // skip over cms
      if (block_idx >= cms_start_block) {
        read_iter = &(read_iter[kLoadPvarBlockSize * sizeof(double)]);
      }
      // skip over chr_idxs
      if (block_idx >= chr_idxs_start_block) {
	read_iter = &(read_iter[kLoadPvarBlockSize * sizeof(chr_idx_t)]);
      }
    }
    memcpy(&(variant_bps[full_block_ct * kLoadPvarBlockSize]), read_iter, raw_variant_ct_lowbits * sizeof(int32_t));
    read_iter = &(read_iter[kLoadPvarBlockSize * (sizeof(int32_t) + sizeof(intptr_t))]);
    memcpy(&(variant_ids[full_block_ct * kLoadPvarBlockSize]), read_iter, raw_variant_ct_lowbits * sizeof(intptr_t));
    read_iter = &(read_iter[kLoadPvarBlockSize * sizeof(intptr_t)]);
    const uint32_t last_bitblock_size = DIV_UP(raw_variant_ct_lowbits, CHAR_BIT);
    memcpy(&(variant_include[full_block_ct * (kLoadPvarBlockSize / kBitsPerWord)]), read_iter, last_bitblock_size);
    zero_trailing_bits(raw_variant_ct, variant_include);
    read_iter = &(read_iter[kLoadPvarBlockSize / CHAR_BIT]);
    if (qual_present) {
      memcpy(&(qual_present[full_block_ct * (kLoadPvarBlockSize / kBitsPerWord)]), read_iter, last_bitblock_size);
      zero_trailing_bits(raw_variant_ct, qual_present);
      read_iter = &(read_iter[kLoadPvarBlockSize / CHAR_BIT]);
      memcpy(&(quals[full_block_ct * kLoadPvarBlockSize]), read_iter, raw_variant_ct_lowbits * sizeof(float));
      read_iter = &(read_iter[kLoadPvarBlockSize * sizeof(float)]);
    }
    if (filter_present) {
      memcpy(&(filter_present[full_block_ct * (kLoadPvarBlockSize / kBitsPerWord)]), read_iter, last_bitblock_size);
      zero_trailing_bits(raw_variant_ct, filter_present);
      read_iter = &(read_iter[kLoadPvarBlockSize / CHAR_BIT]);
      memcpy(&(filter_npass[full_block_ct * (kLoadPvarBlockSize / kBitsPerWord)]), read_iter, last_bitblock_size);
      zero_trailing_bits(raw_variant_ct, filter_npass);
      read_iter = &(read_iter[kLoadPvarBlockSize / CHAR_BIT]);
      if (filter_storage) {
	memcpy(&(filter_storage[full_block_ct * kLoadPvarBlockSize]), read_iter, raw_variant_ct_lowbits * sizeof(intptr_t));
	read_iter = &(read_iter[kLoadPvarBlockSize * sizeof(intptr_t)]);
      }
    }
    if (info_pr_present) {
      memcpy(&(nonref_flags[full_block_ct * (kLoadPvarBlockSize / kBitsPerWord)]), read_iter, last_bitblock_size);
      zero_trailing_bits(raw_variant_ct, nonref_flags);
      // read_iter = &(read_iter[kLoadPvarBlockSize / CHAR_BIT]);
    }
    const uintptr_t read_iter_stride_base = kLoadPvarBlockSize * (sizeof(int32_t) + 2 * sizeof(intptr_t)) + (kLoadPvarBlockSize / CHAR_BIT) + (load_qual_col > 1) * ((kLoadPvarBlockSize / CHAR_BIT) + kLoadPvarBlockSize * sizeof(float)) + (load_filter_col > 1) * (2 * (kLoadPvarBlockSize / CHAR_BIT) + kLoadPvarBlockSize * sizeof(intptr_t)) + info_pr_present * (kLoadPvarBlockSize / CHAR_BIT);
    if (allele_idx_end > 2 * ((uintptr_t)raw_variant_ct)) {
      if (bigstack_alloc_ul(raw_variant_ct + 1, variant_allele_idxs_ptr)) {
	goto load_pvar_ret_NOMEM;
      }
      variant_allele_idxs = *variant_allele_idxs_ptr;
      uintptr_t* allele_idx_read_iter = (uintptr_t*)(&(g_bigstack_end[kLoadPvarBlockSize * sizeof(int32_t)]));
      for (uint32_t block_idx = 0; block_idx < full_block_ct; ++block_idx) {
	memcpy(&(variant_allele_idxs[block_idx * kLoadPvarBlockSize]), allele_idx_read_iter, kLoadPvarBlockSize * sizeof(intptr_t));
	allele_idx_read_iter = (uintptr_t*)(((uintptr_t)allele_idx_read_iter) + read_iter_stride_base + (block_idx >= cms_start_block) * kLoadPvarBlockSize * sizeof(double) + (block_idx >= chr_idxs_start_block) * kLoadPvarBlockSize * sizeof(chr_idx_t));
      }
      memcpy(&(variant_allele_idxs[full_block_ct * kLoadPvarBlockSize]), allele_idx_read_iter, raw_variant_ct_lowbits * sizeof(intptr_t));
      variant_allele_idxs[raw_variant_ct] = allele_idx_end;
    }
    if (at_least_one_nzero_cm) {
      if (bigstack_alloc_d(raw_variant_ct, variant_cms_ptr)) {
	goto load_pvar_ret_NOMEM;
      }
      double* variant_cms = *variant_cms_ptr;
      fill_double_zero(cms_start_block * kLoadPvarBlockSize, variant_cms);
      double* cms_read_iter = (double*)(&(g_bigstack_end[read_iter_stride_base * (cms_start_block + 1)]));
      if (cms_start_block > chr_idxs_start_block) {
	cms_read_iter = (double*)(((uintptr_t)cms_read_iter) + kLoadPvarBlockSize * sizeof(chr_idx_t) * (cms_start_block - chr_idxs_start_block));
      }
      for (uint32_t block_idx = cms_start_block; block_idx < full_block_ct; ++block_idx) {
	memcpy(&(variant_cms[block_idx * kLoadPvarBlockSize]), cms_read_iter, kLoadPvarBlockSize * sizeof(double));
	cms_read_iter = (double*)(((uintptr_t)cms_read_iter) + read_iter_stride_base + kLoadPvarBlockSize * sizeof(double) + (block_idx >= chr_idxs_start_block) * kLoadPvarBlockSize * sizeof(chr_idx_t));
      }
      memcpy(&(variant_cms[full_block_ct * kLoadPvarBlockSize]), cms_read_iter, raw_variant_ct_lowbits * sizeof(double));
    } else {
      *variant_cms_ptr = nullptr;
    }
    if (!is_split_chr) {
      cip->chr_fo_vidx_start[chrs_encountered_m1 + 1] = raw_variant_ct;
      if (splitpar_bound2) {
	if (splitpar_and_exclude_x) {
	  clear_bit(x_code, chr_mask);
	}
	reterr = splitpar(variant_bps, *vpos_sortstatus_ptr, splitpar_bound1, splitpar_bound2, variant_include, loaded_chr_mask, cip, &chrs_encountered_m1, &exclude_ct);
	if (reterr) {
	  goto load_pvar_ret_1;
	}
      } else if (merge_par) {
	if (merge_par_ct) {
	  LOGPRINTF("--merge-par: %u chromosome code%s changed.\n", merge_par_ct, (merge_par_ct == 1)? "" : "s");
	} else {
	  logerrprint("Warning: --merge-par had no effect (no PAR1/PAR2 chromosome codes present).\n");
	}
      }
      cip->chr_ct = chrs_encountered_m1 + 1;
    } else {
      chr_idx_t* chr_idxs = (chr_idx_t*)bigstack_alloc(raw_variant_ct * sizeof(chr_idx_t));
      if (!chr_idxs) {
	goto load_pvar_ret_NOMEM;
      }
      *chr_idxs_ptr = chr_idxs;
      if (chr_idxs_start_block) {
	const uint32_t end_vidx = chr_idxs_start_block * kLoadPvarBlockSize;
	uint32_t chr_fo_idx = chrs_encountered_m1;
	while (cip->chr_fo_vidx_start[chr_fo_idx] >= end_vidx) {
	  --chr_fo_idx;
	}
	backfill_chr_idxs(cip, chr_fo_idx, 0, end_vidx, chr_idxs);
      }
      chr_idx_t* chr_idxs_read_iter = (chr_idx_t*)(&(g_bigstack_end[read_iter_stride_base * (chr_idxs_start_block + 1)]));
      if (chr_idxs_start_block >= cms_start_block) {
	chr_idxs_read_iter = (chr_idx_t*)(((uintptr_t)chr_idxs_read_iter) + kLoadPvarBlockSize * sizeof(double) * (chr_idxs_start_block + 1 - cms_start_block));
      }
      for (uint32_t block_idx = chr_idxs_start_block; block_idx < full_block_ct;) {
	memcpy(&(chr_idxs[block_idx * kLoadPvarBlockSize]), chr_idxs_read_iter, kLoadPvarBlockSize * sizeof(chr_idx_t));
	++block_idx;
	chr_idxs_read_iter = (chr_idx_t*)(((uintptr_t)chr_idxs_read_iter) + read_iter_stride_base + kLoadPvarBlockSize * sizeof(chr_idx_t) + (block_idx >= cms_start_block) * kLoadPvarBlockSize * sizeof(double));
      }
      memcpy(&(chr_idxs[full_block_ct * kLoadPvarBlockSize]), chr_idxs_read_iter, raw_variant_ct_lowbits * sizeof(chr_idx_t));
    }
    const uint32_t last_chr_code = cip->max_code + cip->name_ct;
    const uint32_t chr_word_ct = BITCT_TO_WORDCT(last_chr_code + 1);
    bitvec_and(loaded_chr_mask, chr_word_ct, chr_mask);
    bigstack_end_set(tmp_alloc_end);
    *variant_ct_ptr = raw_variant_ct - exclude_ct;
    *vpos_sortstatus_ptr = vpos_sortstatus;
    *allele_storage_ptr = allele_storage;
    // if only INFO:PR present, no need to reload
    *info_reload_slen_ptr = info_nonpr_present? info_reload_slen : 0;
  }

  while (0) {
  load_pvar_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  load_pvar_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  load_pvar_ret_EMPTY_ALLELE_CODE:
    LOGERRPRINTFWW("Error: Empty allele code on line %" PRIuPTR " of %s.\n", line_idx, pvarname);
    reterr = kPglRetMalformedInput;
    break;
  load_pvar_ret_LONG_LINE:
    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of %s is pathologically long.\n", line_idx, pvarname);
  load_pvar_ret_MALFORMED_INPUT_WW:
    wordwrapb(0);
    logerrprintb();
  load_pvar_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  load_pvar_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  load_pvar_ret_MISSING_TOKENS:
    LOGERRPRINTFWW("Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, pvarname);
    reterr = kPglRetMalformedInput;
    break;
  }
 load_pvar_ret_1:
  if (reterr) {
    bigstack_double_reset(bigstack_mark, bigstack_end_mark);
  }
  gzclose_cond(gz_infile);
  return reterr;
}

#ifdef __cplusplus
} // namespace plink2
#endif
