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


#include "plink2_import.h"
#include "plink2_psam.h"
#include "plink2_pvar.h"
#include "plink2_random.h"

#include "zstd/lib/zstd.h"

#ifdef __cplusplus
namespace plink2 {
#endif

void init_plink1_dosage(plink1_dosage_info_t* plink1_dosage_info_ptr) {
  plink1_dosage_info_ptr->flags = kfPlink1Dosage0;
  fill_uint_zero(3, plink1_dosage_info_ptr->skips);
  plink1_dosage_info_ptr->chr_col_idx = UINT32_MAX;
  plink1_dosage_info_ptr->pos_col_idx = UINT32_MAX;
  plink1_dosage_info_ptr->id_delim = '\0';
}

void init_gendummy(gendummy_info_t* gendummy_info_ptr) {
  gendummy_info_ptr->flags = kfGenDummy0;
  gendummy_info_ptr->pheno_ct = 1;
  gendummy_info_ptr->geno_mfreq = 0.0;
  gendummy_info_ptr->pheno_mfreq = 0.0;
  gendummy_info_ptr->dosage_freq = 0.0;
}

pglerr_t vcf_sample_line(const char* preexisting_psamname, const char* const_fid, uint32_t double_id, fam_col_t fam_cols, char id_delim, char idspace_to, char flag_char, char* sample_line_first_id, char* outname, char* outname_end, uintptr_t* sample_ct_ptr) {
  gzFile gz_infile = nullptr;
  FILE* outfile = nullptr;
  uintptr_t line_idx = 0;
  pglerr_t reterr = kPglRetSuccess;
  {
    uintptr_t const_fid_len = 0;
    if (const_fid) {
      const_fid_len = strlen(const_fid);
    } else if ((!double_id) && (!id_delim)) {
      // default: --double-id + --id-delim
      double_id = 1;
      id_delim = '_';
    }
    const uint32_t double_or_const_fid = double_id || const_fid;
    if (id_delim != ' ') {
      char* sample_line_iter = strchr(sample_line_first_id, ' ');
      if (sample_line_iter) {
        if (!idspace_to) {
          logerrprint("Error: VCF/BCF2 sample ID contains space(s).  Use --idspace-to to convert them\nto another character, or \"--id-delim ' '\" to interpret the spaces as FID/IID\ndelimiters.\n");
          goto vcf_sample_line_ret_INCONSISTENT_INPUT;
        }
        do {
          *sample_line_iter = idspace_to;
          sample_line_iter = strchr(&(sample_line_iter[1]), ' ');
        } while (sample_line_iter);
      }
    }
    char* sample_line_iter = sample_line_first_id;
    uintptr_t sample_ct = 0;
    if (!preexisting_psamname) {
      strcpy(outname_end, ".psam");
      if (fopen_checked(outname, FOPEN_WB, &outfile)) {
        goto vcf_sample_line_ret_OPEN_FAIL;
      }
      char* write_iter = g_textbuf;
      char* textbuf_flush = &(write_iter[kMaxMediumLine]);
      write_iter = strcpya(write_iter, "#FID\tIID");
      uint32_t sid_present = 0;
      if (id_delim) {
        while (((unsigned char)sample_line_iter[0]) >= ' ') {
          char* token_end = strchr(sample_line_iter, '\t');
          if (!token_end) {
            token_end = next_prespace(sample_line_iter);
          }
          char* first_delim = (char*)memchr(sample_line_iter, (unsigned char)id_delim, (uintptr_t)(token_end - sample_line_iter));
          if (first_delim) {
            sample_line_iter = &(first_delim[1]);
            if (memchr(sample_line_iter, (unsigned char)id_delim, (uintptr_t)(token_end - sample_line_iter)) != nullptr) {
              sid_present = 1;
              write_iter = strcpya(write_iter, "\tSID");
              break;
            }
          }
          if (*token_end != '\t') {
            break;
          }
          sample_line_iter = &(token_end[1]);
        }
        sample_line_iter = sample_line_first_id;
      }
      write_iter = strcpya(write_iter, "\tSEX");
      append_binary_eoln(&write_iter);
      while (((unsigned char)sample_line_iter[0]) >= ' ') {
        ++sample_ct;
        char* token_end = strchr(sample_line_iter, '\t');
        if (!token_end) {
          token_end = next_prespace(sample_line_iter);
        }
        const uint32_t token_slen = (uintptr_t)(token_end - sample_line_iter);
        if ((*sample_line_iter == '0') && (token_slen == 1)) {
          logerrprint("Error: Sample ID cannot be '0'.\n");
          goto vcf_sample_line_ret_MALFORMED_INPUT;
        }
        if (id_delim) {
          if (*sample_line_iter == id_delim) {
            sprintf(g_logbuf, "Error: '%c' at beginning of sample ID.\n", id_delim);
            goto vcf_sample_line_ret_INCONSISTENT_INPUT_2;
          }
          if (sample_line_iter[token_slen - 1] == id_delim) {
            sprintf(g_logbuf, "Error: '%c' at end of sample ID.\n", id_delim);
            goto vcf_sample_line_ret_INCONSISTENT_INPUT_2;
          }
          char* first_delim = (char*)memchr(sample_line_iter, (unsigned char)id_delim, token_slen);
          if (!first_delim) {
            if (double_or_const_fid) {
              goto vcf_sample_line_nopsam_one_id;
            }
            sprintf(g_logbuf, "Error: No '%c' in sample ID.\n", id_delim);
            goto vcf_sample_line_ret_INCONSISTENT_INPUT_2;
          }
          char* iid_start = &(first_delim[1]);
          char* iid_end = (char*)memchr(iid_start, (unsigned char)id_delim, (uintptr_t)(token_end - iid_start));
          const char* sid_start = &(g_one_char_strs[96]);
          uint32_t sid_slen = 1;
          if (iid_end) {
            if (iid_start == iid_end) {
              sprintf(g_logbuf, "Error: Consecutive instances of '%c' in sample ID.\n", id_delim);
              goto vcf_sample_line_ret_INCONSISTENT_INPUT_DELIM;
            }
            sid_start = &(iid_end[1]);
            sid_slen = (uintptr_t)(token_end - sid_start);
            if (memchr(sid_start, (unsigned char)id_delim, sid_slen)) {
              sprintf(g_logbuf, "Error: More than two instances of '%c' in sample ID.\n", id_delim);
              goto vcf_sample_line_ret_INCONSISTENT_INPUT_DELIM;
            }
            if (sid_slen > kMaxIdSlen) {
              logerrprint("Error: SIDs are limited to " MAX_ID_SLEN_STR " characters.\n");
              goto vcf_sample_line_ret_MALFORMED_INPUT;
            }
          } else {
            iid_end = token_end;
          }
          const uint32_t fid_slen = (uintptr_t)(first_delim - sample_line_iter);
          if (fid_slen > kMaxIdSlen) {
            // strictly speaking, you could have e.g. a 20k char ID which
            // splits into a valid FID/IID pair with the right delimiter, but
            // you're not supposed to have sample IDs anywhere near that length
            // so I'll classify this as MalformedInput.
            goto vcf_sample_line_ret_MALFORMED_INPUT_LONG_ID;
          }
          write_iter = memcpyax(write_iter, sample_line_iter, fid_slen, '\t');
          const uint32_t iid_slen = (uintptr_t)(iid_end - iid_start);
          if ((*iid_start == '0') && (iid_slen == 1)) {
            logerrprint("Error: Sample ID induces an invalid IID of '0'.\n");
            goto vcf_sample_line_ret_INCONSISTENT_INPUT;
          }
          if (iid_slen > kMaxIdSlen) {
            goto vcf_sample_line_ret_MALFORMED_INPUT_LONG_ID;
          }
          write_iter = memcpya(write_iter, iid_start, iid_slen);
          if (sid_present) {
            *write_iter++ = '\t';
            write_iter = memcpya(write_iter, sid_start, sid_slen);
          }
        } else {
        vcf_sample_line_nopsam_one_id:
          if (token_slen > kMaxIdSlen) {
            goto vcf_sample_line_ret_MALFORMED_INPUT_LONG_ID;
          }
          if (double_id) {
            write_iter = memcpya(write_iter, sample_line_iter, token_slen);
          } else {
            write_iter = memcpya(write_iter, const_fid, const_fid_len);
          }
          *write_iter++ = '\t';
          write_iter = memcpya(write_iter, sample_line_iter, token_slen);
          if (sid_present) {
            write_iter = strcpya(write_iter, "\t0");
          }
        }
        // PAT/MAT/PHENO1 not required in .psam file
        // SEX now included, so that --vcf + --out has the same effect as --vcf
        // + --make-pgen + --out
        write_iter = memcpyl3a(write_iter, "\tNA");
        append_binary_eoln(&write_iter);
        if (fwrite_ck(textbuf_flush, outfile, &write_iter)) {
          goto vcf_sample_line_ret_WRITE_FAIL;
        }
        if (*token_end != '\t') {
          break;
        }
        sample_line_iter = &(token_end[1]);
      }
      if (fclose_flush_null(textbuf_flush, write_iter, &outfile)) {
        goto vcf_sample_line_ret_WRITE_FAIL;
      }
    } else {
      // check consistency of IIDs between VCF and .psam file.
      reterr = gzopen_read_checked(preexisting_psamname, &gz_infile);
      if (reterr) {
        goto vcf_sample_line_ret_1;
      }
      uintptr_t loadbuf_size = bigstack_left();
      if (loadbuf_size > kMaxLongLine) {
        loadbuf_size = kMaxLongLine;
      } else if (loadbuf_size <= kMaxMediumLine) {
        goto vcf_sample_line_ret_NOMEM;
      }
      char* loadbuf = (char*)g_bigstack_base;
      // not formally allocated for now
      loadbuf[loadbuf_size - 1] = ' ';
      char* loadbuf_first_token;
      do {
        ++line_idx;
        if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
          if (!gzeof(gz_infile)) {
            goto vcf_sample_line_ret_READ_FAIL;
          }
          loadbuf_first_token = loadbuf;
          loadbuf_first_token[0] = '\0';
          break;
        }
        if (!loadbuf[loadbuf_size - 1]) {
          if (loadbuf_size == kMaxLongLine) {
            goto vcf_sample_line_ret_LONG_LINE;
          }
          goto vcf_sample_line_ret_NOMEM;
        }
        loadbuf_first_token = skip_initial_spaces(loadbuf);
      } while (is_eoln_kns(*loadbuf_first_token) || ((loadbuf_first_token[0] == '#') && (!tokequal_k(&(loadbuf_first_token[1]), "FID")) && (!tokequal_k(&(loadbuf_first_token[1]), "IID"))));
      uint32_t fid_present;
      if (loadbuf_first_token[0] == '#') {
        // only care about position of IID column
        fid_present = (loadbuf_first_token[1] == 'F');
      } else {
        fid_present = fam_cols & kfFamCol1;
      }
      while (1) {
        if (!is_eoln_kns(*loadbuf_first_token)) {
          char* psam_iid_start = loadbuf_first_token;
          if (fid_present) {
            psam_iid_start = skip_initial_spaces(token_endnn(psam_iid_start));
            if (is_eoln_kns(*psam_iid_start)) {
              goto vcf_sample_line_ret_MISSING_TOKENS;
            }
          }
          if (((unsigned char)sample_line_iter[0]) < ' ') {
            sprintf(g_logbuf, "Error: --%ccf file contains fewer sample IDs than %s.\n", flag_char, preexisting_psamname);
            goto vcf_sample_line_ret_INCONSISTENT_INPUT_WW;
          }
          ++sample_ct;
          char* sample_line_token_end = strchr(sample_line_iter, '\t');
          if (!sample_line_token_end) {
            sample_line_token_end = next_prespace(sample_line_iter);
          }
          uint32_t sample_line_token_slen = (uintptr_t)(sample_line_token_end - sample_line_iter);
          if ((*sample_line_iter == '0') && (sample_line_token_slen == 1)) {
            logerrprint("Error: Sample ID cannot be '0'.\n");
            goto vcf_sample_line_ret_MALFORMED_INPUT;
          }
          char* sample_line_iid_start = sample_line_iter;
          if (id_delim) {
            char* first_delim = (char*)memchr(sample_line_iter, (unsigned char)id_delim, sample_line_token_slen);
            if (!first_delim) {
              if (!double_or_const_fid) {
                sprintf(g_logbuf, "Error: No '%c' in sample ID.\n", id_delim);
                goto vcf_sample_line_ret_INCONSISTENT_INPUT_2;
              }
            } else {
              sample_line_iid_start = &(first_delim[1]);
              sample_line_token_slen = (uintptr_t)(sample_line_token_end - sample_line_iid_start);
              char* sample_line_iid_end = (char*)memchr(sample_line_iid_start, (unsigned char)id_delim, sample_line_token_slen);
              if (sample_line_iid_end) {
                // don't bother erroring out on >2 instances of delimiter for
                // now
                sample_line_token_slen = (uintptr_t)(sample_line_iid_end - sample_line_iid_start);
              }
              if ((*sample_line_iid_start == '0') && (sample_line_token_slen == 1)) {
                logerrprint("Error: Sample ID induces an invalid IID of '0'.\n");
                goto vcf_sample_line_ret_INCONSISTENT_INPUT;
              }
            }
          }
          if (sample_line_token_slen > kMaxIdSlen) {
            goto vcf_sample_line_ret_MALFORMED_INPUT_LONG_ID;
          }
          if (memcmp(sample_line_iid_start, psam_iid_start, sample_line_token_slen) || (((unsigned char)psam_iid_start[sample_line_token_slen]) > 32)) {
            sprintf(g_logbuf, "Error: Mismatched IDs between --%ccf file and %s.\n", flag_char, preexisting_psamname);
            goto vcf_sample_line_ret_INCONSISTENT_INPUT_WW;
          }
          sample_line_iter = &(sample_line_token_end[1]);
        }
        ++line_idx;
        if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
          if (!gzeof(gz_infile)) {
            goto vcf_sample_line_ret_READ_FAIL;
          }
          break;
        }
        if (!loadbuf[loadbuf_size - 1]) {
          if (loadbuf_size == kMaxLongLine) {
            goto vcf_sample_line_ret_LONG_LINE;
          }
          goto vcf_sample_line_ret_NOMEM;
        }
        loadbuf_first_token = skip_initial_spaces(loadbuf);
        if (loadbuf_first_token[0] == '#') {
          sprintf(g_logbuf, "Error: Line %" PRIuPTR " of %s starts with a '#'. (This is only permitted before the first nonheader line, and a #FID/IID header line is present it must denote the end of the header block.)\n", line_idx, preexisting_psamname);
          goto vcf_sample_line_ret_MALFORMED_INPUT_WW;
        }
      }
      if (gzclose_null(&gz_infile)) {
        goto vcf_sample_line_ret_READ_FAIL;
      }
      if (((unsigned char)sample_line_iter[0]) >= ' ') {
        sprintf(g_logbuf, "Error: --%ccf file contains more sample IDs than %s.\n", flag_char, preexisting_psamname);
        goto vcf_sample_line_ret_INCONSISTENT_INPUT_WW;
      }
    }
    *sample_ct_ptr = sample_ct;
  }
  while (0) {
  vcf_sample_line_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  vcf_sample_line_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  vcf_sample_line_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  vcf_sample_line_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  vcf_sample_line_ret_LONG_LINE:
    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of %s is pathologically long.\n", line_idx, preexisting_psamname);
  vcf_sample_line_ret_MALFORMED_INPUT_WW:
    wordwrapb(0);
    logerrprintb();
  vcf_sample_line_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  vcf_sample_line_ret_MISSING_TOKENS:
    LOGERRPRINTFWW("Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, preexisting_psamname);
    reterr = kPglRetMalformedInput;
    break;
  vcf_sample_line_ret_INCONSISTENT_INPUT_DELIM:
    logerrprintb();
    if (id_delim == '_') {
      logerrprint("If you do not want '_' to be treated as a FID/IID delimiter, use --double-id or\n--const-fid to choose a different method of converting VCF sample IDs to PLINK\nIDs, or --id-delim to change the FID/IID delimiter.\n");
    }
    reterr = kPglRetInconsistentInput;
    break;
  vcf_sample_line_ret_MALFORMED_INPUT_LONG_ID:
    logerrprint("Error: FIDs and IIDs are limited to " MAX_ID_SLEN_STR " characters.\n");
    reterr = kPglRetMalformedInput;
    break;
  vcf_sample_line_ret_INCONSISTENT_INPUT_WW:
    wordwrapb(0);
  vcf_sample_line_ret_INCONSISTENT_INPUT_2:
    logerrprintb();
  vcf_sample_line_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 vcf_sample_line_ret_1:
  gzclose_cond(gz_infile);
  fclose_cond(outfile);
  return reterr;
}

uint32_t vcf_is_het_short(const char* first_gchar_ptr, vcf_half_call_t vcf_half_call) {
  // '.' == ascii 46, '0' == ascii 48
  // if kVcfHalfCallReference, ./0 is not phased, but ./1 is
  const uint32_t first_gchar = (unsigned char)first_gchar_ptr[0];
  const uint32_t second_gchar = (unsigned char)first_gchar_ptr[2];
  return (first_gchar != second_gchar) && (((first_gchar != 46) && (second_gchar != 46)) || ((vcf_half_call == kVcfHalfCallReference) && ((first_gchar > 48) || (second_gchar > 48))));
}

uint32_t get_vcf_format_position(const char* __restrict needle, uint32_t needle_slen, char* format_start, char* format_end) {
  *format_end = '\0';
  // assumes first field is GT
  char* token_start = &(format_start[3]);
  uint32_t field_idx = 0;
  while (1) {
    char* token_end = strchr(token_start, ':');
    ++field_idx;
    if (!token_end) {
      if (((uintptr_t)(format_end - token_start) != needle_slen) || memcmp(token_start, needle, needle_slen)) {
        field_idx = 0;
      }
      break;
    }
    if (((uintptr_t)(token_end - token_start) == needle_slen) && (!memcmp(token_start, needle, needle_slen))) {
      break;
    }
    token_start = &(token_end[1]);
  }
  *format_end = '\t';
  return field_idx;
}

uint32_t vcf_qual_scan_init(int32_t vcf_min_gq, int32_t vcf_min_dp, char* format_start, char* format_end, uint32_t* qual_field_skips, int32_t* qual_thresholds) {
  uint32_t gq_field_idx = 0;
  if (vcf_min_gq >= 0) {
    gq_field_idx = get_vcf_format_position("GQ", 2, format_start, format_end);
  }
  uint32_t qual_field_ct = 0;
  uint32_t dp_field_idx = 0;
  if (vcf_min_dp >= 0) {
    dp_field_idx = get_vcf_format_position("DP", 2, format_start, format_end);
    if (dp_field_idx && ((!gq_field_idx) || (dp_field_idx < gq_field_idx))) {
      qual_field_skips[0] = dp_field_idx;
      qual_thresholds[0] = vcf_min_dp;
      qual_field_ct = 1;
      dp_field_idx = 0;
    }
  }
  if (gq_field_idx) {
    qual_field_skips[qual_field_ct] = gq_field_idx;
    qual_thresholds[qual_field_ct++] = vcf_min_gq;
    if (dp_field_idx) {
      qual_field_skips[qual_field_ct] = dp_field_idx;
      qual_thresholds[qual_field_ct++] = vcf_min_dp;
    }
  }
  if (qual_field_ct == 2) {
    qual_field_skips[1] -= qual_field_skips[0];
  }
  return qual_field_ct;
}

char* vcf_field_advance(char* gtext_iter, char* gtext_end, uint32_t advance_ct) {
  // assumes advance_ct is nonzero
  do {
    char* field_end = (char*)memchr(gtext_iter, ':', gtext_end - gtext_iter);
    if (!field_end) {
      return nullptr;
    }
    gtext_iter = &(field_end[1]);
  } while (--advance_ct);
  return gtext_iter;
}

// returns 1 if a quality check failed
// assumes either 1 or 2 qual fields, otherwise change this to a loop
uint32_t vcf_check_quals(const uint32_t* qual_field_skips, const int32_t* qual_thresholds, char* gtext_iter, char* gtext_end, uint32_t qual_field_ct) {
  gtext_iter = vcf_field_advance(gtext_iter, gtext_end, qual_field_skips[0]);
  if (!gtext_iter) {
    return 0;
  }
  int32_t ii;
  if ((!scan_int32(gtext_iter, &ii)) && (ii < qual_thresholds[0])) {
    return 1;
  }
  if (qual_field_ct == 1) {
    return 0;
  }
  gtext_iter = vcf_field_advance(gtext_iter, gtext_end, qual_field_skips[1]);
  if (!gtext_iter) {
    return 0;
  }
  return (!scan_int32(gtext_iter, &ii)) && (ii < qual_thresholds[1]);
}

boolerr_t parse_vcf_gp(char* gp_iter, uint32_t is_haploid, double import_dosage_certainty, uint32_t* is_missing_ptr, double* alt_dosage_ptr) {
  // P(0/0), P(0/1), P(1/1), etc.
  // assumes is_missing initialized to 0
  double prob_0alt;
  gp_iter = scanadv_double(gp_iter, &prob_0alt);
  if ((!gp_iter) || (prob_0alt < 0.0) || (prob_0alt > 1.0) || (*gp_iter != ',')) {
    return 1;
  }
  double prob_1alt;
  gp_iter = scanadv_double(&(gp_iter[1]), &prob_1alt);
  if ((!gp_iter) || (prob_1alt < 0.0) || (prob_1alt > 1.0)) {
    return 1;
  }
  if (is_haploid) {
    const double denom = prob_0alt + prob_1alt;
    if (denom <= 2 * import_dosage_certainty) {
      if ((prob_0alt <= import_dosage_certainty) && (prob_1alt <= import_dosage_certainty)) {
        *is_missing_ptr = 1;
        return 1;
      }
    }
    *alt_dosage_ptr = 2 * prob_1alt / denom;
    return 0;
  }
  double prob_2alt;
  if ((*gp_iter != ',') || (!scanadv_double(&(gp_iter[1]), &prob_2alt)) || (prob_2alt < 0.0) || (prob_2alt > 1.0)) {
    return 1;
  }
  const double denom = prob_0alt + prob_1alt + prob_2alt;
  if (denom <= 3 * import_dosage_certainty) {
    if ((prob_0alt <= import_dosage_certainty) && (prob_1alt <= import_dosage_certainty) && (prob_2alt <= import_dosage_certainty)) {
      // treat as missing
      // ok to use <= since we multiplied by (1 - epsilon)
      // during command-line parsing.  this lets us avoid
      // special-casing denom=0.
      *is_missing_ptr = 1;
      return 1; // not really an error
    }
  }
  *alt_dosage_ptr = (prob_1alt + 2 * prob_2alt) / denom;
  return 0;
}

boolerr_t parse_vcf_dosage(char* gtext_iter, char* gtext_end, uint32_t dosage_field_idx, uint32_t is_haploid, uint32_t dosage_is_gp, double import_dosage_certainty, uint32_t* is_missing_ptr, uint32_t* dosage_int_ptr) {
  // assumes is_missing initialized to 0
  // assumes dosage_field_idx != 0
  // returns 1 if missing OR parsing error.
  uint32_t field_idx = 0;
  do {
    char* field_end = (char*)memchr(gtext_iter, ':', gtext_end - gtext_iter);
    if (!field_end) {
      *is_missing_ptr = 1;
      return 1;
    }
    gtext_iter = &(field_end[1]);
  } while (++field_idx < dosage_field_idx);
  if (((gtext_iter[0] == '.') || (gtext_iter[0] == '?')) && (((uint32_t)((unsigned char)gtext_iter[1])) - 48 >= 10)) {
    // missing field (dot/'?' followed by non-digit)
    // could enforce gtext_iter[1] == colon, comma, etc.?
    *is_missing_ptr = 1;
    return 1;
  }
  double alt_dosage;
  if (dosage_is_gp) {
    if (parse_vcf_gp(gtext_iter, is_haploid, import_dosage_certainty, is_missing_ptr, &alt_dosage)) {
      return 1;
    }
  } else {
    if ((!scanadv_double(gtext_iter, &alt_dosage)) || (alt_dosage < 0.0)) {
      return 1;
    }
    if (is_haploid) {
      // possible todo: allow this to be suppressed (maybe upstream of this
      // function); 1000 Genomes phase 1 haploid dosages are still on 0..2
      // scale
      alt_dosage *= 2;
    }
    if (alt_dosage > 2.0) {
      return 1;
    }
  }
  *dosage_int_ptr = (int32_t)(alt_dosage * kDosageMid + 0.5);
  return 0;
}

static_assert(!kVcfHalfCallReference, "vcf_to_pgen() assumes kVcfHalfCallReference == 0.");
static_assert(kVcfHalfCallHaploid == 1, "vcf_to_pgen() assumes kVcfHalfCallHaploid == 1.");
pglerr_t vcf_to_pgen(const char* vcfname, const char* preexisting_psamname, const char* const_fid, const char* dosage_import_field, misc_flags_t misc_flags, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, double import_dosage_certainty, char id_delim, char idspace_to, int32_t vcf_min_gq, int32_t vcf_min_dp, vcf_half_call_t vcf_half_call, fam_col_t fam_cols, char* outname, char* outname_end, chr_info_t* cip) {
  // Now performs a 2-pass load.  Yes, this can be slower than plink 1.9, but
  // it's necessary to use the Pgen_writer classes for now (since we need to
  // know upfront how many variants there are, and whether phase/dosage is
  // present).
  // preexisting_psamname should be nullptr if no such file was specified.
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  FILE* pvarfile = nullptr;
  uintptr_t line_idx = 0;
  const uint32_t vcf_half_call_explicit_error = (vcf_half_call == kVcfHalfCallError);
  st_pgen_writer_t spgw;
  pglerr_t reterr = kPglRetSuccess;
  gzFile gz_infile;
  spgw_preinit(&spgw);
  {
    // don't use gzopen_read_checked() since we want to customize the error
    // message
    gz_infile = gzopen(vcfname, FOPEN_RB);
    if (!gz_infile) {
      const uint32_t slen = strlen(vcfname);
      if (str_endswith(vcfname, ".vcf", slen) ||
          str_endswith(vcfname, ".vcf.gz", slen)) {
        LOGERRPRINTFWW(g_errstr_fopen, vcfname);
      } else {
        LOGERRPRINTFWW("Error: Failed to open %s. (--vcf expects a complete filename; did you forget '.vcf' at the end?)\n", vcfname);
      }
      goto vcf_to_pgen_ret_OPEN_FAIL;
    }
    if (gzbuffer(gz_infile, 131072)) {
      goto vcf_to_pgen_ret_NOMEM;
    }
    uintptr_t loadbuf_size = bigstack_left() / 4;
    if (loadbuf_size > kMaxLongLine) {
      loadbuf_size = kMaxLongLine;
    } else if (loadbuf_size <= kMaxMediumLine) {
      goto vcf_to_pgen_ret_NOMEM;
    } else {
      loadbuf_size = round_up_pow2(loadbuf_size, kEndAllocAlign);
    }
    char* loadbuf = (char*)bigstack_end_alloc_raw(loadbuf_size);
    loadbuf[loadbuf_size - 1] = ' ';
    char* loadbuf_iter;
    const uint32_t allow_extra_chrs = (misc_flags / kfMiscAllowExtraChrs) & 1;
    uint32_t dosage_import_field_slen = 0;
    if (dosage_import_field) {
      dosage_import_field_slen = strlen(dosage_import_field);
    }
    const uint32_t dosage_is_gp = strequal_k(dosage_import_field, "GP", dosage_import_field_slen);
    uint32_t format_gt_present = 0;
    uint32_t format_gq_relevant = 0;
    uint32_t format_dp_relevant = 0;
    uint32_t format_dosage_relevant = 0;
    uint32_t info_pr_present = 0;
    uint32_t info_pr_nonflag_present = 0;
    uint32_t info_nonpr_present = 0;
    uint32_t chrset_present = 0;
    while (1) {
      ++line_idx;
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
        if (!gzeof(gz_infile)) {
          goto vcf_to_pgen_ret_READ_FAIL;
        }
        logerrprint("Error: No #CHROM header line or variant records in --vcf file.\n");
        goto vcf_to_pgen_ret_MALFORMED_INPUT;
      }
      if ((line_idx == 1) && (!memcmp(loadbuf, "BCF", 3))) {
        // this is more informative than "missing header line"...
        if (loadbuf[3] == 2) {
          sprintf(g_logbuf, "Error: %s appears to be a BCF2 file. Try --bcf instead of --vcf.\n", vcfname);
          goto vcf_to_pgen_ret_MALFORMED_INPUT_WW;
        }
        if (loadbuf[3] == 4) {
          sprintf(g_logbuf, "Error: %s appears to be a BCF1 file. Use 'bcftools view' to convert it to a PLINK-readable VCF.\n", vcfname);
          goto vcf_to_pgen_ret_MALFORMED_INPUT_WW;
        }
      }
      if (!loadbuf[loadbuf_size - 1]) {
        if (loadbuf_size == kMaxLongLine) {
          goto vcf_to_pgen_ret_LONG_LINE;
        }
        goto vcf_to_pgen_ret_NOMEM;
      }
      // don't tolerate leading spaces
      loadbuf_iter = loadbuf;
      if (*loadbuf_iter != '#') {
        logerrprint("Error: No #CHROM header line in --vcf file.\n");
        goto vcf_to_pgen_ret_MALFORMED_INPUT;
      }
      if (loadbuf_iter[1] != '#') {
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
      if (str_startswith2(&(loadbuf_iter[2]), "chrSet=<")) {
        if (chrset_present) {
          logerrprint("Error: Multiple ##chrSet header lines in --vcf file.\n");
          goto vcf_to_pgen_ret_MALFORMED_INPUT;
        }
        chrset_present = 1;
        // .pvar loader will print a warning if necessary
        reterr = read_chrset_header_line(&(loadbuf_iter[10]), "--vcf file", misc_flags, line_idx, cip);
        if (reterr) {
          goto vcf_to_pgen_ret_1;
        }
      } else if (str_startswith2(&(loadbuf_iter[2]), "FORMAT=<ID=GT,Number=")) {
        if (format_gt_present) {
          logerrprint("Error: Duplicate FORMAT:GT header line in --vcf file.\n");
          goto vcf_to_pgen_ret_MALFORMED_INPUT;
        }
        if (!str_startswith2(&(loadbuf_iter[2 + strlen("FORMAT=<ID=GT,Number=")]), "1,Type=String,Description=")) {
          sprintf(g_logbuf, "Error: Header line %" PRIuPTR " of --vcf file does not have expected FORMAT:GT format.\n", line_idx);
          goto vcf_to_pgen_ret_MALFORMED_INPUT_WW;
        }
        format_gt_present = 1;
      } else if ((vcf_min_gq != -1) && str_startswith2(&(loadbuf_iter[2]), "FORMAT=<ID=GQ,Number=1,Type=")) {
        if (format_gq_relevant) {
          logerrprint("Error: Duplicate FORMAT:GQ header line in --vcf file.\n");
          goto vcf_to_pgen_ret_MALFORMED_INPUT;
        }
        format_gq_relevant = 1;
      } else if ((vcf_min_dp != -1) && str_startswith2(&(loadbuf_iter[2]), "FORMAT=<ID=DP,Number=1,Type=")) {
        if (format_dp_relevant) {
          logerrprint("Error: Duplicate FORMAT:DP header line in --vcf file.\n");
          goto vcf_to_pgen_ret_MALFORMED_INPUT;
        }
        format_dp_relevant = 1;
      } else if (dosage_import_field && str_startswith2(&(loadbuf_iter[2]), "FORMAT=<ID=") && (!memcmp(&(loadbuf_iter[2 + strlen("FORMAT=<ID=")]), dosage_import_field, dosage_import_field_slen)) && (loadbuf_iter[2 + strlen("FORMAT=<ID=") + dosage_import_field_slen] == ',')) {
        if (format_dosage_relevant) {
          LOGERRPRINTFWW("Error: Duplicate FORMAT:%s header line in --vcf file.\n", dosage_import_field);
          goto vcf_to_pgen_ret_MALFORMED_INPUT_WW;
        }
        format_dosage_relevant = 1;
      } else if (str_startswith2(&(loadbuf_iter[2]), "INFO=<ID=")) {
        if (str_startswith2(&(loadbuf_iter[2 + strlen("INFO=<ID=")]), "PR,Number=")) {
          if (info_pr_present || info_pr_nonflag_present) {
            logerrprint("Error: Duplicate INFO:PR header line in --vcf file.\n");
            goto vcf_to_pgen_ret_MALFORMED_INPUT;
          }
          info_pr_nonflag_present = !str_startswith2(&(loadbuf_iter[2 + strlen("INFO=<ID=PR,Number=")]), "0,Type=Flag,Description=");
          info_pr_present = 1 - info_pr_nonflag_present;
          if (info_pr_nonflag_present) {
            LOGERRPRINTFWW("Warning: Header line %" PRIuPTR " of --vcf file has an unexpected definition of INFO:PR. This interferes with a few merge and liftover operations.\n", line_idx);
          }
        } else {
          info_nonpr_present = 1;
        }
      }
    }
    const uint32_t require_gt = (misc_flags / kfMiscVcfRequireGt) & 1;
    if ((!format_gt_present) && require_gt) {
      // todo: allow_no_variants exception
      logerrprint("Error: No GT field in --vcf file header, when --vcf-require-gt was specified.\n");
      goto vcf_to_pgen_ret_INCONSISTENT_INPUT;
    }
    if ((!format_gq_relevant) && (vcf_min_gq != -1)) {
      logerrprint("Warning: No GQ field in --vcf file header.  --vcf-min-gq ignored.\n");
      vcf_min_gq = -1;
    }
    if ((!format_dp_relevant) && (vcf_min_dp != -1)) {
      logerrprint("Warning: No DP field in --vcf file header.  --vcf-min-dp ignored.\n");
      vcf_min_dp = -1;
    }
    const uint32_t format_gq_or_dp_relevant = format_gq_relevant || format_dp_relevant;
    if ((!format_dosage_relevant) && dosage_import_field) {
      LOGERRPRINTFWW("Warning: No %s field in --vcf file header. Dosages will not be imported.\n", dosage_import_field);
    }
    finalize_chrset(misc_flags, cip);
    // don't call finalize_chr_info here, since this may be followed by
    // --pmerge, etc.

    if (!str_startswith2(loadbuf_iter, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")) {
      sprintf(g_logbuf, "Error: Header line %" PRIuPTR " of --vcf file does not have expected field sequence after #CHROM.\n", line_idx);
      goto vcf_to_pgen_ret_MALFORMED_INPUT_WW;
    }
    loadbuf_iter = &(loadbuf_iter[strlen("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")]);
    const uint32_t double_id = (misc_flags / kfMiscDoubleId) & 1;
    uintptr_t sample_ct = 0;
    if (str_startswith2(loadbuf_iter, "\tFORMAT\t")) {
      reterr = vcf_sample_line(preexisting_psamname, const_fid, double_id, fam_cols, id_delim, idspace_to, 'v', &(loadbuf_iter[strlen("\tFORMAT\t")]), outname, outname_end, &sample_ct);
      if (reterr) {
        goto vcf_to_pgen_ret_1;
      }
      /*
    } else if (allow_no_samples) {
      strcpy(outname_end, ".psam");
      if (fopen_checked(outname, "w", &outfile)) {
        goto vcf_to_pgen_ret_OPEN_FAIL;
      }
      fputs("#FID\tIID\n", outfile);
      if (fclose_null(&outfile)) {
        goto vcf_to_pgen_ret_WRITE_FAIL;
      }
      */
    }
    // todo: allow_no_samples exception
    if (!sample_ct) {
      logerrprint("Error: No samples in --vcf file.\n");
      goto vcf_to_pgen_ret_INCONSISTENT_INPUT;
    }

    uint32_t variant_ct = 0;
    uint32_t max_alt_ct = 1;
    uintptr_t* variant_allele_idxs = (uintptr_t*)g_bigstack_base;
    uintptr_t max_variant_ct = (uintptr_t)(((uintptr_t*)g_bigstack_end) - variant_allele_idxs);
    max_variant_ct -= BITCT_TO_ALIGNED_WORDCT(max_variant_ct) * kWordsPerVec;
    if (format_dosage_relevant) {
      max_variant_ct -= BITCT_TO_ALIGNED_WORDCT(max_variant_ct) * kWordsPerVec;
    }
    if (info_pr_present) {
      max_variant_ct -= BITCT_TO_ALIGNED_WORDCT(max_variant_ct) * kWordsPerVec;
    }
#ifdef __LP64__
    if (max_variant_ct > 0x7ffffffd) {
      max_variant_ct = 0x7ffffffd;
    }
#endif
    uintptr_t base_chr_present[kChrExcludeWords];
    fill_ulong_zero(kChrExcludeWords, base_chr_present);

    const uintptr_t header_line_ct = line_idx;
    const uint32_t max_variant_ctaw = BITCT_TO_ALIGNED_WORDCT(max_variant_ct);
    uintptr_t* phasing_flags = (uintptr_t*)bigstack_end_alloc_raw_rd(max_variant_ctaw * sizeof(intptr_t));
    uintptr_t* phasing_flags_iter = phasing_flags;
    uintptr_t* dosage_flags = nullptr;
    if (format_dosage_relevant) {
      dosage_flags = (uintptr_t*)bigstack_end_alloc_raw_rd(max_variant_ctaw * sizeof(intptr_t));
    }
    uintptr_t* dosage_flags_iter = dosage_flags;
    uintptr_t* nonref_flags = nullptr;
    if (info_pr_present) {
      nonref_flags = (uintptr_t*)bigstack_end_alloc_raw_rd(max_variant_ctaw * sizeof(intptr_t));
    }
    uintptr_t* nonref_flags_iter = nonref_flags;
    if (vcf_half_call == kVcfHalfCallDefault) {
      vcf_half_call = kVcfHalfCallError;
    }
    uintptr_t variant_skip_ct = 0;
    uintptr_t phasing_word = 0;
    uintptr_t dosage_word = 0;
    uintptr_t nonref_word = 0;
    uintptr_t allele_idx_end = 0;
    uint32_t max_allele_slen = 1;
    uint32_t max_qualfilterinfo_slen = 6;
    uint32_t qual_field_ct = 0;

    const uint32_t dosage_erase_halfdist = kDosage4th - dosage_erase_thresh;

    // temporary kludge
    uintptr_t multiallelic_skip_ct = 0;

    while (1) {
      ++line_idx;
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
        if (!gzeof(gz_infile)) {
          goto vcf_to_pgen_ret_READ_FAIL;
        }
        break;
      }
      if (!loadbuf[loadbuf_size - 1]) {
        if (loadbuf_size == kMaxLongLine) {
          goto vcf_to_pgen_ret_LONG_LINE_N;
        }
        goto vcf_to_pgen_ret_NOMEM;
      }
      // do tolerate trailing newlines
      if ((unsigned char)(*loadbuf) <= 32) {
        if (*loadbuf == ' ') {
          sprintf(g_logbuf, "Error: Leading space on line %" PRIuPTR " of --vcf file.\n", line_idx);
          goto vcf_to_pgen_ret_MALFORMED_INPUT_2N;
        }
        continue;
      }
      loadbuf_iter = loadbuf;
      char* chr_code_end = strchr(loadbuf, '\t');
      if (!chr_code_end) {
        goto vcf_to_pgen_ret_MISSING_TOKENS;
      }
      // QUAL/FILTER enforcement is now postponed till .pvar loading.  only
      // other things we do during the scanning pass are (i) count alt alleles,
      // and (ii) check whether any phased genotype calls are present.

      char* pos_end = strchr(&(chr_code_end[1]), '\t');
      if (!pos_end) {
        goto vcf_to_pgen_ret_MISSING_TOKENS;
      }

      // may as well check ID length here
      // postpone POS validation till second pass so we only have to parse it
      // once
      char* id_end = strchr(&(pos_end[1]), '\t');
      if (!id_end) {
        goto vcf_to_pgen_ret_MISSING_TOKENS;
      }
      if ((uintptr_t)(id_end - pos_end) > kMaxIdBlen) {
        sprintf(g_logbuf, "Error: Invalid ID on line %" PRIuPTR " of --vcf file (max " MAX_ID_SLEN_STR " chars).\n", line_idx);
        goto vcf_to_pgen_ret_MALFORMED_INPUT_WW;
      }

      // note REF length
      char* ref_allele_start = &(id_end[1]);
      loadbuf_iter = strchr(ref_allele_start, '\t');
      if (!loadbuf_iter) {
        goto vcf_to_pgen_ret_MISSING_TOKENS;
      }
      uint32_t cur_max_allele_slen = (uintptr_t)(loadbuf_iter - ref_allele_start);

      uint32_t alt_ct = 1;
      unsigned char ucc;
      // treat ALT=. as if it were an actual allele for now
      while (1) {
        char* cur_allele_start = ++loadbuf_iter;
        ucc = (unsigned char)(*loadbuf_iter);
        if ((ucc <= ',') && (ucc != '*')) {
          sprintf(g_logbuf, "Error: Invalid alternate allele on line %" PRIuPTR " of --vcf file.\n", line_idx);
          goto vcf_to_pgen_ret_MALFORMED_INPUT_2N;
        }
        do {
          ucc = (unsigned char)(*(++loadbuf_iter));
          // allow GATK 3.4 <*:DEL> symbolic allele
        } while ((ucc > ',') || (ucc == '*'));
        const uint32_t cur_allele_slen = (uintptr_t)(loadbuf_iter - cur_allele_start);
        if (cur_allele_slen > cur_max_allele_slen) {
          cur_max_allele_slen = cur_allele_slen;
        }
        if (ucc != ',') {
          break;
        }
        ++alt_ct;
      }

      // temporary kludge
      if (alt_ct > 1) {
        ++multiallelic_skip_ct;
        continue;
      }

      if (ucc != '\t') {
        sprintf(g_logbuf, "Error: Malformed ALT field on line %" PRIuPTR " of --vcf file.\n", line_idx);
        goto vcf_to_pgen_ret_MALFORMED_INPUT_2N;
      }
      if (alt_ct > max_alt_ct) {
        max_alt_ct = alt_ct;
      }

      // skip QUAL, FILTER
      char* qual_start_m1 = loadbuf_iter;
      for (uint32_t uii = 0; uii < 2; ++uii) {
        loadbuf_iter = strchr(&(loadbuf_iter[1]), '\t');
        if (!loadbuf_iter) {
          goto vcf_to_pgen_ret_MISSING_TOKENS;
        }
      }

      // possibly check for FORMAT:GT before proceeding
      char* info_start = &(loadbuf_iter[1]);
      char* info_end = strchr(info_start, '\t');
      if (!info_end) {
        goto vcf_to_pgen_ret_MISSING_TOKENS;
      }
      loadbuf_iter = &(info_end[1]);
      const uint32_t gt_missing = memcmp(loadbuf_iter, "GT", 2) || ((loadbuf_iter[2] != ':') && (loadbuf_iter[2] != '\t'));
      if (require_gt && gt_missing) {
        ++variant_skip_ct;
        continue;
      }
      const uint32_t cur_qualfilterinfo_slen = (uintptr_t)(info_end - qual_start_m1);

      // all converters *do* respect chromosome filters
      // wait till this point to apply it, since we don't want to
      // add a contig name to the hash table unless at least one variant on
      // that contig wasn't filtered out for other reasons.
      int32_t cur_chr_code;
      reterr = get_or_add_chr_code_destructive("--vcf file", line_idx, allow_extra_chrs, loadbuf, chr_code_end, cip, &cur_chr_code);
      if (reterr) {
        goto vcf_to_pgen_ret_1;
      }
      if (!is_set(cip->chr_mask, cur_chr_code)) {
        ++variant_skip_ct;
        continue;
      }
      if (cur_max_allele_slen > max_allele_slen) {
        max_allele_slen = cur_max_allele_slen;
      }
      if (cur_qualfilterinfo_slen > max_qualfilterinfo_slen) {
        max_qualfilterinfo_slen = cur_qualfilterinfo_slen;
      }
      if ((uint32_t)cur_chr_code <= cip->max_code) {
        set_bit(cur_chr_code, base_chr_present);
      }

      variant_allele_idxs[variant_ct] = allele_idx_end;
      allele_idx_end += alt_ct + 1;
      const uint32_t variant_idx_lowbits = variant_ct % kBitsPerWord;
      if (info_pr_present) {
        if (pr_in_info_token((uintptr_t)(info_end - info_start), info_start)) {
          nonref_word |= k1LU << variant_idx_lowbits;
        }
        if (variant_idx_lowbits == (kBitsPerWord - 1)) {
          *nonref_flags_iter++ = nonref_word;
          nonref_word = 0;
        }
      }
      if (!gt_missing) {
        // todo: sample_ct == 0 case
        // possible todo: import dosages when GT missing

        // loadbuf_iter currently points to beginning of FORMAT field
        char* format_end = strchr(loadbuf_iter, '\t');
        if (!format_end) {
          goto vcf_to_pgen_ret_MISSING_TOKENS;
        }

        uint32_t qual_field_skips[2];
        int32_t qual_thresholds[2];
        if (format_gq_or_dp_relevant) {
          qual_field_ct = vcf_qual_scan_init(vcf_min_gq, vcf_min_dp, loadbuf_iter, format_end, qual_field_skips, qual_thresholds);
        }
        // if nonzero, 0-based index of dosage field
        uint32_t dosage_field_idx = 0;
        if (format_dosage_relevant) {
          dosage_field_idx = get_vcf_format_position(dosage_import_field, dosage_import_field_slen, loadbuf_iter, format_end);
        }

        // check if there's at least one phased het call, and/or at least one
        // relevant dosage
        if (alt_ct < 10) {
          // always check for a phased het
          char* phasescan_iter = format_end;
          while (1) {
            // this should quickly fail if there are no phased calls at all.
            phasescan_iter = strchr(&(phasescan_iter[1]), '|');
            if (!phasescan_iter) {
              break;
            }
            if (phasescan_iter[-2] != '\t') {
              // at least one other gdata field uses the '|' character.
              // switch to iterating over tabs.
              while (1) {
                phasescan_iter = strchr(&(phasescan_iter[1]), '\t');
                if (!phasescan_iter) {
                  break;
                }
                if (phasescan_iter[2] == '|') {
                  if (vcf_is_het_short(&(phasescan_iter[1]), vcf_half_call)) {
                    if (qual_field_ct) {
                      char* cur_gtext_end = strchr(phasescan_iter, '\t');
                      if (!cur_gtext_end) {
                        cur_gtext_end = next_prespace(phasescan_iter);
                      }
                      if (vcf_check_quals(qual_field_skips, qual_thresholds, phasescan_iter, cur_gtext_end, qual_field_ct)) {
                        break;
                      }
                    }
                    phasing_word |= k1LU << variant_idx_lowbits;
                    break;
                  }
                }
              }
              break;
            }
            if (vcf_is_het_short(&(phasescan_iter[-1]), vcf_half_call)) {
              if (qual_field_ct) {
                char* cur_gtext_end = strchr(phasescan_iter, '\t');
                if (!cur_gtext_end) {
                  cur_gtext_end = next_prespace(phasescan_iter);
                }
                if (vcf_check_quals(qual_field_skips, qual_thresholds, phasescan_iter, cur_gtext_end, qual_field_ct)) {
                  break;
                }
              }
              phasing_word |= k1LU << variant_idx_lowbits;
              break;
            }
            phasescan_iter = strchr(&(phasescan_iter[1]), '\t');
            if (!phasescan_iter) {
              break;
            }
          }
          if (dosage_field_idx) {
            char* dosagescan_iter = format_end;
            for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
              char* cur_gtext_start = ++dosagescan_iter;
              char* cur_gtext_end = strchr(dosagescan_iter, '\t');
              if (!cur_gtext_end) {
                if (sample_idx + 1 != sample_ct) {
                  goto vcf_to_pgen_ret_MISSING_TOKENS;
                }
                cur_gtext_end = next_prespace(dosagescan_iter);
              }
              dosagescan_iter = cur_gtext_end;
              if (qual_field_ct) {
                if (vcf_check_quals(qual_field_skips, qual_thresholds, cur_gtext_start, cur_gtext_end, qual_field_ct)) {
                  break;
                }
              }
              const uint32_t is_haploid = (cur_gtext_start[1] != '/') && (cur_gtext_start[1] != '|');
              uint32_t is_missing = 0;
              uint32_t dosage_int;
              if (parse_vcf_dosage(cur_gtext_start, cur_gtext_end, dosage_field_idx, is_haploid, dosage_is_gp, import_dosage_certainty, &is_missing, &dosage_int)) {
                if (is_missing) {
                  continue;
                }
                goto vcf_to_pgen_ret_INVALID_DOSAGE;
              }
              const uint32_t cur_halfdist = biallelic_dosage_halfdist(dosage_int);
              if (cur_halfdist < dosage_erase_halfdist) {
                goto vcf_to_pgen_dosagescan_hit;
              }
            }
          }
          if (0) {
          vcf_to_pgen_dosagescan_hit:
            dosage_word |= k1LU << variant_idx_lowbits;
          }
        } else {
          // alt_ct >= 10
          // todo
        }
      }
      if (variant_idx_lowbits == (kBitsPerWord - 1)) {
        *phasing_flags_iter++ = phasing_word;
        phasing_word = 0;
        if (dosage_flags_iter) {
          *dosage_flags_iter++ = dosage_word;
          dosage_word = 0;
        }
      }
      if (variant_ct++ == max_variant_ct) {
#ifdef __LP64__
        if (variant_ct == 0x7ffffffd) {
          logerrprint("Error: " PROG_NAME_STR " does not support more than 2^31 - 3 variants.  We recommend other\nsoftware, such as PLINK/SEQ, for very deep studies of small numbers of genomes.\n");
          goto vcf_to_pgen_ret_MALFORMED_INPUT;
        }
#endif
        goto vcf_to_pgen_ret_NOMEM;
      }
      if (!(variant_ct % 1000)) {
        printf("\r--vcf: %uk variants scanned.", variant_ct / 1000);
        fflush(stdout);
      }
    }
    if (variant_ct % kBitsPerWord) {
      *phasing_flags_iter = phasing_word;
      if (dosage_flags_iter) {
        *dosage_flags_iter = dosage_word;
      }
      if (nonref_flags_iter) {
        *nonref_flags_iter = nonref_word;
      }
    } else if (!variant_ct) {
      // todo: allow_no_variants exception
      logerrprint("Error: No variants in --vcf file.\n");
      goto vcf_to_pgen_ret_INCONSISTENT_INPUT;
    }

    putc_unlocked('\r', stdout);
    if (!variant_skip_ct) {
      LOGPRINTF("--vcf: %u variant%s scanned.\n", variant_ct, (variant_ct == 1)? "" : "s");
    } else {
      LOGPRINTF("--vcf: %u variant%s scanned (%" PRIuPTR " skipped).\n", variant_ct, (variant_ct == 1)? "" : "s", variant_skip_ct);
    }

    // temporary kludge
    if (multiallelic_skip_ct) {
      LOGERRPRINTFWW("Warning: %" PRIuPTR " multiallelic variant%s %sskipped (not yet supported).\n", multiallelic_skip_ct, (multiallelic_skip_ct == 1)? "" : "s", variant_skip_ct? "also " : "");
    }

    if (gzrewind(gz_infile)) {
      goto vcf_to_pgen_ret_READ_FAIL;
    }
    const uintptr_t line_ct = line_idx - 1;

    if (allele_idx_end > 2 * variant_ct) {
      variant_allele_idxs[variant_ct] = allele_idx_end;
      bigstack_finalize_ul(variant_allele_idxs, variant_ct + 1);
    } else {
      variant_allele_idxs = nullptr;
    }

    strcpy(outname_end, ".pvar");
    if (fopen_checked(outname, FOPEN_WB, &pvarfile)) {
      goto vcf_to_pgen_ret_OPEN_FAIL;
    }
    line_idx = 0;
    while (1) {
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
        goto vcf_to_pgen_ret_READ_FAIL;
      }
      if (++line_idx == header_line_ct) {
        break;
      }
      // don't use textbuf here, since header line length could theoretically
      // exceed kMaxMediumLine bytes
      if (str_startswith2(loadbuf, "##fileformat=") || str_startswith2(loadbuf, "##fileDate=") || str_startswith2(loadbuf, "##source=") || str_startswith2(loadbuf, "##FORMAT=") || str_startswith2(loadbuf, "##chrSet=")) {
        continue;
      }
      const uint32_t line_slen = strlen(loadbuf);
      if (str_startswith2(loadbuf, "##contig=<ID=")) {
        char* contig_name_start = &(loadbuf[strlen("##contig=<ID=")]);
        char* contig_name_end = strchr(contig_name_start, ',');
        if (!contig_name_end) {
          // could search backwards from end of line in this case, but if
          // contig names are long enough for that to matter we have other
          // problems...
          contig_name_end = strchr(contig_name_start, '>');
          if (!contig_name_end) {
            sprintf(g_logbuf, "Error: Header line %" PRIuPTR " of --vcf file does not have expected ##contig format.\n", line_idx);
            goto vcf_to_pgen_ret_MALFORMED_INPUT_WW;
          }
        }
        const int32_t cur_chr_code = get_chr_code_counted(cip, (uintptr_t)(contig_name_end - contig_name_start), contig_name_start);
        if (cur_chr_code < 0) {
          continue;
        }
        if ((uint32_t)cur_chr_code <= cip->max_code) {
          if (!is_set(base_chr_present, cur_chr_code)) {
            continue;
          }
        } else {
          if (!is_set(cip->chr_mask, cur_chr_code)) {
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
      // don't need to check for \n since this is pass #2, we already validated
      // file contents, and we can't be on the last line of a file lacking a
      // final eoln since #CHROM is still to come
      char* line_end = &(loadbuf[line_slen - 1]);
      if (line_end[-1] == '\r') {
        --line_end;
      }
      append_binary_eoln(&line_end);
      if (fwrite_checked(loadbuf, (uintptr_t)(line_end - loadbuf), pvarfile)) {
        goto vcf_to_pgen_ret_WRITE_FAIL;
      }
    }
    char* write_iter = g_textbuf;
    if (cip->chrset_source) {
      append_chrset_line(cip, &write_iter);
    }
    write_iter = strcpya(write_iter, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER");
    if (info_nonpr_present) {
      write_iter = strcpya(write_iter, "\tINFO");
    }
    append_binary_eoln(&write_iter);
    if (fwrite_checked(g_textbuf, write_iter - g_textbuf, pvarfile)) {
      goto vcf_to_pgen_ret_WRITE_FAIL;
    }

    const uint32_t variant_ctl = BITCT_TO_WORDCT(variant_ct);
    pgen_global_flags_t phase_dosage_gflags = are_all_words_zero(phasing_flags, variant_ctl)? kfPgenGlobal0 : kfPgenGlobalHardcallPhasePresent;
    if (format_dosage_relevant) {
      if (are_all_words_zero(dosage_flags, variant_ctl)) {
        format_dosage_relevant = 0;
      } else {
        phase_dosage_gflags |= kfPgenGlobalDosagePresent;
      }
    }
    uint32_t nonref_flags_storage = 1;
    if (nonref_flags) {
      const uint32_t variant_ctl_m1 = variant_ctl - 1;
      const uintptr_t last_nonref_flags_word = nonref_flags[variant_ctl_m1];
      if (!last_nonref_flags_word) {
        for (uint32_t widx = 0; widx < variant_ctl_m1; ++widx) {
          if (nonref_flags[widx]) {
            nonref_flags_storage = 3;
            break;
          }
        }
        // todo: replace kBitsPerWord - MOD_NZ...
      } else if (!((~last_nonref_flags_word) << (kBitsPerWord - MOD_NZ(variant_ct, kBitsPerWord)))) {
        nonref_flags_storage = 2;
        for (uint32_t widx = 0; widx < variant_ctl_m1; ++widx) {
          if (~nonref_flags[widx]) {
            nonref_flags_storage = 3;
            break;
          }
        }
      } else {
        nonref_flags_storage = 3;
      }
      if (nonref_flags_storage != 3) {
        nonref_flags = nullptr;
        bigstack_end_reset(phasing_flags);
      }
    }
    strcpy(outname_end, ".pgen");
    uintptr_t spgw_alloc_cacheline_ct;
    uint32_t max_vrec_len;
    reterr = spgw_init_phase1(outname, variant_allele_idxs, nonref_flags, variant_ct, sample_ct, phase_dosage_gflags, nonref_flags_storage, &spgw, &spgw_alloc_cacheline_ct, &max_vrec_len);
    if (reterr) {
      goto vcf_to_pgen_ret_1;
    }
    unsigned char* spgw_alloc;
    if (bigstack_alloc_uc(spgw_alloc_cacheline_ct * kCacheline, &spgw_alloc)) {
      goto vcf_to_pgen_ret_NOMEM;
    }
    spgw_init_phase2(max_vrec_len, &spgw, spgw_alloc);

    const uint32_t sample_ctl2 = QUATERCT_TO_WORDCT(sample_ct);
    uintptr_t* genovec;
    // if we weren't using bigstack_alloc, this would need to be sample_ctaw2
    if (bigstack_alloc_ul(sample_ctl2, &genovec)) {
      goto vcf_to_pgen_ret_NOMEM;
    }
    // nothing should go wrong if trailing word is garbage, but keep an eye on
    // this
    // fill_ulong_zero(sample_ctaw2 - sample_ctl2, &(genovec[sample_ctl2]));
    const uint32_t sample_ctl = BITCT_TO_WORDCT(sample_ct);
    uintptr_t* phasepresent = nullptr;
    uintptr_t* phaseinfo = nullptr;
    if (phase_dosage_gflags & kfPgenGlobalHardcallPhasePresent) {
      if (bigstack_alloc_ul(sample_ctl, &phasepresent) ||
          bigstack_alloc_ul(sample_ctl, &phaseinfo)) {
        goto vcf_to_pgen_ret_NOMEM;
      }
      // bugfix (7 Sep 2017): phasepresent can't have nonzero trailing bits
      phasepresent[sample_ctl - 1] = 0;
    }
    uintptr_t* dosage_present = nullptr;
    dosage_t* dosage_vals = nullptr;
    if (phase_dosage_gflags & kfPgenGlobalDosagePresent) {
      if (bigstack_alloc_ul(sample_ctl, &dosage_present) ||
          bigstack_alloc_dosage(sample_ct, &dosage_vals)) {
        goto vcf_to_pgen_ret_NOMEM;
      }
      dosage_present[sample_ctl - 1] = 0;
    }

    char* writebuf;
    if (bigstack_alloc_c(2 * max_allele_slen + max_qualfilterinfo_slen + kMaxMediumLine + kMaxIdSlen + 32, &writebuf)) {
      goto vcf_to_pgen_ret_NOMEM;
    }
    write_iter = writebuf;
    char* writebuf_flush = &(writebuf[kMaxMediumLine]);

    if (hard_call_thresh == UINT32_MAX) {
      hard_call_thresh = kDosageMid / 10;
    }
    const uint32_t hard_call_halfdist = kDosage4th - hard_call_thresh;

    uint32_t vidx = 0;
    for (line_idx = header_line_ct + 1; line_idx <= line_ct; ++line_idx) {
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
        goto vcf_to_pgen_ret_READ_FAIL;
      }
      if ((unsigned char)(*loadbuf) < 32) {
        continue;
      }
      // 1. check if we skip this variant.  chromosome filter, require_gt, and
      //    (temporarily) multiple alt alleles can cause this.
      char* chr_code_end = (char*)rawmemchr(loadbuf, '\t');
      int32_t chr_code_base = get_chr_code_raw(loadbuf);
      if (chr_code_base == -1) {
        // skip hash table lookup if we know we aren't skipping the variant
        if (variant_skip_ct) {
          *chr_code_end = '\0';
          const uint32_t chr_code = id_htable_find(loadbuf, cip->nonstd_names, cip->nonstd_id_htable, (uintptr_t)(chr_code_end - loadbuf), kChrHtableSize);
          if ((chr_code == UINT32_MAX) || (!IS_SET(cip->chr_mask, chr_code))) {
            continue;
          }
          *chr_code_end = '\t';
        }
      } else {
        if (chr_code_base >= ((int32_t)kMaxContigs)) {
          chr_code_base = cip->xymt_codes[chr_code_base - kMaxContigs];
        }
        if ((chr_code_base < 0) || (!is_set(base_chr_present, chr_code_base))) {
          assert(variant_skip_ct);
          continue;
        }
      }
      // chr_code_base is now a proper numeric chromosome index for
      // non-contigs, and -1 if it's a contig name
      char* pos_str = &(chr_code_end[1]);
      char* pos_str_end = (char*)rawmemchr(pos_str, '\t');
      loadbuf_iter = pos_str_end;
      // copy ID, REF verbatim
      for (uint32_t uii = 0; uii < 2; ++uii) {
        loadbuf_iter = (char*)rawmemchr(&(loadbuf_iter[1]), '\t');
      }

      // ALT, QUAL, FILTER, INFO
      char* filter_end = loadbuf_iter;
      for (uint32_t uii = 0; uii < 3; ++uii) {
        filter_end = (char*)rawmemchr(&(filter_end[1]), '\t');
      }
      char* info_end = (char*)rawmemchr(&(filter_end[1]), '\t');
      char* format_start = &(info_end[1]);
      const uint32_t gt_missing = memcmp(format_start, "GT", 2) || ((format_start[2] != ':') && (format_start[2] != '\t'));
      if (require_gt && gt_missing) {
        continue;
      }

      // make sure POS starts with an integer, apply --output-chr setting
      uint32_t cur_bp;
      if (scan_uint_defcap(pos_str, &cur_bp)) {
        sprintf(g_logbuf, "Error: Invalid POS on line %" PRIuPTR " of --vcf file.\n", line_idx);
        goto vcf_to_pgen_ret_MALFORMED_INPUT_2N;
      }

      // temporary kludge
      char* write_line_start = write_iter;

      if (chr_code_base == -1) {
        write_iter = memcpya(write_iter, loadbuf, (uintptr_t)(chr_code_end - loadbuf));
      } else {
        write_iter = chr_name_write(cip, chr_code_base, write_iter);
      }
      *write_iter++ = '\t';
      write_iter = uint32toa(cur_bp, write_iter);

      uint32_t alt_ct = 1;
      char* copy_start = pos_str_end;
      while (1) {
        ++loadbuf_iter;
        unsigned char ucc;
        do {
          ucc = (unsigned char)(*(++loadbuf_iter));
          // allow GATK 3.4 <*:DEL> symbolic allele
        } while ((ucc > ',') || (ucc == '*'));

        // temporary kludge
        if (ucc == ',') {
          alt_ct = 2;
          break;
        }

        write_iter = memcpya(write_iter, copy_start, (uintptr_t)(loadbuf_iter - copy_start));
        // unsafe to flush for now due to multiallelic kludge
        /*
        if (fwrite_ck(writebuf_flush, pvarfile, &write_iter)) {
          goto vcf_to_pgen_ret_WRITE_FAIL;
        }
        */
        if (ucc != ',') {
          break;
        }
        copy_start = loadbuf_iter;
        ++alt_ct;
      }

      // temporary kludge
      if (alt_ct > 1) {
        write_iter = write_line_start;
        continue;
      }
      if (fwrite_ck(writebuf_flush, pvarfile, &write_iter)) {
        goto vcf_to_pgen_ret_WRITE_FAIL;
      }

      write_iter = memcpya(write_iter, loadbuf_iter, (uintptr_t)((info_nonpr_present ? info_end : filter_end) - loadbuf_iter));
      append_binary_eoln(&write_iter);

      if (gt_missing) {
        fill_all_bits(2 * sample_ct, genovec);
        if (spgw_append_biallelic_genovec(genovec, &spgw)) {
          goto vcf_to_pgen_ret_WRITE_FAIL;
        }
      } else {
        loadbuf_iter = (char*)rawmemchr(format_start, '\t');
        uint32_t qual_field_skips[2];
        int32_t qual_thresholds[2];
        if (format_gq_or_dp_relevant) {
          qual_field_ct = vcf_qual_scan_init(vcf_min_gq, vcf_min_dp, format_start, loadbuf_iter, qual_field_skips, qual_thresholds);
        }

        uint32_t dosage_field_idx = 0;
        dosage_t* dosage_vals_iter = dosage_vals;
        if (dosage_flags && IS_SET(dosage_flags, vidx)) {
          dosage_field_idx = get_vcf_format_position(dosage_import_field, dosage_import_field_slen, format_start, loadbuf_iter);
        }

        // todo: multiallelic variants
        ++loadbuf_iter;
        const uint32_t sample_ctl2_m1 = sample_ctl2 - 1;
        uint32_t inner_loop_last = kBitsPerWordD2 - 1;
        uint32_t widx = 0;
        if (!IS_SET(phasing_flags, vidx)) {
          while (1) {
            if (widx >= sample_ctl2_m1) {
              if (widx > sample_ctl2_m1) {
                break;
              }
              inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
            }
            uintptr_t genovec_word = 0;
            uint32_t dosage_present_hw = 0;
            for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
              char* cur_gtext_end = strchr(loadbuf_iter, '\t');
              if (!cur_gtext_end) {
                if ((sample_idx_lowbits != inner_loop_last) || (widx != sample_ctl2_m1)) {
                  goto vcf_to_pgen_ret_MISSING_TOKENS;
                }
                cur_gtext_end = next_prespace(loadbuf_iter);
              }
              uintptr_t cur_geno;
              if (qual_field_ct) {
                if (vcf_check_quals(qual_field_skips, qual_thresholds, loadbuf_iter, cur_gtext_end, qual_field_ct)) {
                  goto vcf_to_pgen_force_missing_1;
                }
              }
              {
                const uint32_t is_haploid = (loadbuf_iter[1] != '/') && (loadbuf_iter[1] != '|');
                cur_geno = (unsigned char)(*loadbuf_iter) - 48;
                if (cur_geno <= 1) {
                  if (is_haploid) {
                    cur_geno *= 2;
                  } else {
                    const char cc = loadbuf_iter[3];
                    if (((cc != '/') && (cc != '|')) || (loadbuf_iter[4] == '.')) {
                      // code triploids, etc. as missing
                      // might want to subject handling of 0/0/. to
                      // --vcf-half-call control
                      const uintptr_t second_allele_idx = ((unsigned char)loadbuf_iter[2]) - 48;
                      if (second_allele_idx <= 1) {
                        cur_geno += second_allele_idx;
                      } else if (second_allele_idx != (uintptr_t)((intptr_t)(-2))) {
                        // not '.'
                        goto vcf_to_pgen_ret_INVALID_GT;
                      } else if (vcf_half_call == kVcfHalfCallMissing) {
                        cur_geno = 3;
                      } else if (vcf_half_call == kVcfHalfCallError) {
                        goto vcf_to_pgen_ret_HALF_CALL_ERROR;
                      } else {
                        // kVcfHalfCallHaploid, kVcfHalfCallReference
                        cur_geno <<= vcf_half_call;
                      }
                    }
                  }
                } else if (cur_geno != (uintptr_t)((intptr_t)(-2))) {
                  // not '.'
                  goto vcf_to_pgen_ret_INVALID_GT;
                } else if (vcf_half_call != kVcfHalfCallMissing) {
                  const char second_allele_char = loadbuf_iter[2];
                  if ((second_allele_char != '.') && ((loadbuf_iter[1] == '/') || (loadbuf_iter[1] == '|'))) {
                    cur_geno = ((unsigned char)second_allele_char) - 48;
                    if (cur_geno > 1) {
                      goto vcf_to_pgen_ret_INVALID_GT;
                    }
                    if (vcf_half_call == kVcfHalfCallError) {
                      goto vcf_to_pgen_ret_HALF_CALL_ERROR;
                    }
                    // kVcfHalfCallHaploid, kVcfHalfCallReference
                    cur_geno <<= vcf_half_call;
                  } else {
                    cur_geno = 3;
                  }
                } else {
                  cur_geno = 3;
                }
                if (dosage_field_idx) {
                  uint32_t is_missing = 0;
                  uint32_t dosage_int;
                  if (!parse_vcf_dosage(loadbuf_iter, cur_gtext_end, dosage_field_idx, is_haploid, dosage_is_gp, import_dosage_certainty, &is_missing, &dosage_int)) {
                    const uint32_t cur_halfdist = biallelic_dosage_halfdist(dosage_int);
                    if ((cur_geno == 3) && (cur_halfdist >= hard_call_halfdist)) {
                      // only overwrite GT if (i) it was missing, and (ii) the
                      // dosage's distance from the nearest hardcall doesn't
                      // exceed the --hard-call-threshold value.
                      // (possible todo: warn or error out if dosage and GT are
                      // inconsistent)
                      cur_geno = (dosage_int + (kDosage4th * k1LU)) / kDosageMid;
                    }
                    if (cur_halfdist < dosage_erase_halfdist) {
                      dosage_present_hw |= 1U << sample_idx_lowbits;
                      *dosage_vals_iter++ = dosage_int;
                    }
                  } else if (!is_missing) {
                    goto vcf_to_pgen_ret_INVALID_DOSAGE;
                  }
                }
              }
              while (0) {
              vcf_to_pgen_force_missing_1:
                cur_geno = 3;
              }
              genovec_word |= cur_geno << (2 * sample_idx_lowbits);
              loadbuf_iter = &(cur_gtext_end[1]);
            }
            genovec[widx] = genovec_word;
            if (dosage_field_idx) {
              ((halfword_t*)dosage_present)[widx] = dosage_present_hw;
            }
            ++widx;
          }
          if (!dosage_field_idx) {
            if (spgw_append_biallelic_genovec(genovec, &spgw)) {
              goto vcf_to_pgen_ret_WRITE_FAIL;
            }
          } else {
            assert(dosage_vals_iter != dosage_vals);
            if (spgw_append_biallelic_genovec_dosage16(genovec, dosage_present, dosage_vals, dosage_vals_iter - dosage_vals, &spgw)) {
              goto vcf_to_pgen_ret_WRITE_FAIL;
            }
          }
        } else {
          while (1) {
            if (widx >= sample_ctl2_m1) {
              if (widx > sample_ctl2_m1) {
                break;
              }
              inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
            }
            uintptr_t genovec_word = 0;
            uint32_t phasepresent_halfword = 0;
            uint32_t phaseinfo_halfword = 0;
            uint32_t dosage_present_hw = 0;
            for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
              char* cur_gtext_end = strchr(loadbuf_iter, '\t');
              if (!cur_gtext_end) {
                if ((sample_idx_lowbits != inner_loop_last) || (widx != sample_ctl2_m1)) {
                  goto vcf_to_pgen_ret_MISSING_TOKENS;
                }
                cur_gtext_end = next_prespace(loadbuf_iter);
              }
              uintptr_t cur_geno;
              if (qual_field_ct) {
                if (vcf_check_quals(qual_field_skips, qual_thresholds, loadbuf_iter, cur_gtext_end, qual_field_ct)) {
                  goto vcf_to_pgen_force_missing_2;
                }
              }
              {
                const uint32_t is_phased = (loadbuf_iter[1] == '|');
                const uint32_t is_haploid = (!is_phased) && (loadbuf_iter[1] != '/');
                cur_geno = (unsigned char)(*loadbuf_iter) - 48;
                if (cur_geno <= 1) {
                  if (is_haploid) {
                    cur_geno *= 2;
                  } else {
                    const char cc = loadbuf_iter[3];
                    if (((cc != '/') && (cc != '|')) || (loadbuf_iter[4] == '.')) {
                      // code triploids, etc. as missing
                      // might want to subject handling of 0/0/. to
                      // --vcf-half-call control
                      const uintptr_t second_allele_idx = ((unsigned char)loadbuf_iter[2]) - 48;
                      if (second_allele_idx <= 1) {
                        cur_geno += second_allele_idx;
                        // todo: check if this should be less branchy
                        if (is_phased && (cur_geno == 1)) {
                          const uint32_t shifted_bit = 1U << sample_idx_lowbits;
                          phasepresent_halfword |= shifted_bit;
                          if (!second_allele_idx) {
                            // 1|0
                            phaseinfo_halfword |= shifted_bit;
                          }
                        }
                      } else if (second_allele_idx != (uintptr_t)((intptr_t)(-2))) {
                        // not '.'
                        goto vcf_to_pgen_ret_INVALID_GT;
                      } else if (vcf_half_call == kVcfHalfCallMissing) {
                        cur_geno = 3;
                      } else if (vcf_half_call == kVcfHalfCallError) {
                        goto vcf_to_pgen_ret_HALF_CALL_ERROR;
                      } else {
                        // kVcfHalfCallHaploid, kVcfHalfCallReference
                        cur_geno <<= vcf_half_call;
                      }
                    }
                  }
                } else if (cur_geno != (uintptr_t)((intptr_t)(-2))) {
                  // not '.'
                  goto vcf_to_pgen_ret_INVALID_GT;
                } else if (vcf_half_call != kVcfHalfCallMissing) {
                  const char second_allele_char = loadbuf_iter[2];
                  if ((second_allele_char != '.') && ((loadbuf_iter[1] == '/') || (loadbuf_iter[1] == '|'))) {
                    cur_geno = ((unsigned char)second_allele_char) - 48;
                    if (cur_geno > 1) {
                      goto vcf_to_pgen_ret_INVALID_GT;
                    }
                    if (vcf_half_call == kVcfHalfCallError) {
                      goto vcf_to_pgen_ret_HALF_CALL_ERROR;
                    }
                    // kVcfHalfCallHaploid, kVcfHalfCallReference
                    cur_geno <<= vcf_half_call;
                  } else {
                    cur_geno = 3;
                  }
                } else {
                  cur_geno = 3;
                }
                if (dosage_field_idx) {
                  uint32_t is_missing = 0;
                  uint32_t dosage_int;
                  if (!parse_vcf_dosage(loadbuf_iter, cur_gtext_end, dosage_field_idx, is_haploid, dosage_is_gp, import_dosage_certainty, &is_missing, &dosage_int)) {
                    const uint32_t cur_halfdist = biallelic_dosage_halfdist(dosage_int);
                    if ((cur_geno == 3) && (cur_halfdist >= hard_call_halfdist)) {
                      cur_geno = (dosage_int + (kDosage4th * k1LU)) / kDosageMid;
                    }
                    if (cur_halfdist < dosage_erase_halfdist) {
                      dosage_present_hw |= 1U << sample_idx_lowbits;
                      *dosage_vals_iter++ = dosage_int;
                    }
                  } else if (!is_missing) {
                    goto vcf_to_pgen_ret_INVALID_DOSAGE;
                  }
                }
              }
              while (0) {
              vcf_to_pgen_force_missing_2:
                cur_geno = 3;
              }
              genovec_word |= cur_geno << (2 * sample_idx_lowbits);
              loadbuf_iter = &(cur_gtext_end[1]);
            }
            genovec[widx] = genovec_word;
            ((halfword_t*)phasepresent)[widx] = (halfword_t)phasepresent_halfword;
            ((halfword_t*)phaseinfo)[widx] = (halfword_t)phaseinfo_halfword;
            if (dosage_field_idx) {
              ((halfword_t*)dosage_present)[widx] = dosage_present_hw;
            }
            ++widx;
          }
          if (!dosage_field_idx) {
            if (spgw_append_biallelic_genovec_hphase(genovec, phasepresent, phaseinfo, &spgw)) {
              goto vcf_to_pgen_ret_WRITE_FAIL;
            }
          } else {
            if (spgw_append_biallelic_genovec_hphase_dosage16(genovec, phasepresent, phaseinfo, dosage_present, dosage_vals, dosage_vals_iter - dosage_vals, &spgw)) {
              goto vcf_to_pgen_ret_WRITE_FAIL;
            }
          }
        }
      }
      if (!(++vidx % 1000)) {
        printf("\r--vcf: %uk variants converted.", vidx / 1000);
        fflush(stdout);
      }
    }
    if (fclose_flush_null(writebuf_flush, write_iter, &pvarfile)) {
      goto vcf_to_pgen_ret_WRITE_FAIL;
    }
    spgw_finish(&spgw);
    putc_unlocked('\r', stdout);
    write_iter = strcpya(g_logbuf, "--vcf: ");
    const uint32_t outname_base_slen = (uintptr_t)(outname_end - outname);
    write_iter = memcpya(write_iter, outname, outname_base_slen + 5);
    write_iter = memcpyl3a(write_iter, " + ");
    write_iter = memcpya(write_iter, outname, outname_base_slen);
    write_iter = strcpya(write_iter, ".pvar");
    if (!preexisting_psamname) {
      write_iter = memcpyl3a(write_iter, " + ");
      write_iter = memcpya(write_iter, outname, outname_base_slen);
      write_iter = strcpya(write_iter, ".psam");
    }
    strcpy(write_iter, " written.\n");
    wordwrapb(0);
    logprintb();
  }
  while (0) {
  vcf_to_pgen_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  vcf_to_pgen_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  vcf_to_pgen_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  vcf_to_pgen_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  vcf_to_pgen_ret_HALF_CALL_ERROR:
    putc_unlocked('\n', stdout);
    LOGERRPRINTF("Error: Line %" PRIuPTR " of --vcf file has a GT half-call.\n", line_idx);
    if (!vcf_half_call_explicit_error) {
      logerrprint("Use --vcf-half-call to specify how these should be processed.\n");
    }
    reterr = kPglRetMalformedInput;
    break;
  vcf_to_pgen_ret_INVALID_GT:
    putc_unlocked('\n', stdout);
    LOGERRPRINTF("Error: Line %" PRIuPTR " of --vcf file has an invalid GT field.\n", line_idx);
    reterr = kPglRetMalformedInput;
    break;
  vcf_to_pgen_ret_MISSING_TOKENS:
    putc_unlocked('\n', stdout);
    LOGERRPRINTF("Error: Line %" PRIuPTR " of --vcf file has fewer tokens than expected.\n", line_idx);
    reterr = kPglRetMalformedInput;
    break;
  vcf_to_pgen_ret_LONG_LINE_N:
    putc_unlocked('\n', stdout);
  vcf_to_pgen_ret_LONG_LINE:
    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of --vcf file is pathologically long.\n", line_idx);
  vcf_to_pgen_ret_MALFORMED_INPUT_2N:
    logprint("\n");
    logerrprintb();
  vcf_to_pgen_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  vcf_to_pgen_ret_MALFORMED_INPUT_WW:
    wordwrapb(0);
    logerrprintb();
    reterr = kPglRetMalformedInput;
    break;
  vcf_to_pgen_ret_INVALID_DOSAGE:
    putc_unlocked('\n', stdout);
    LOGERRPRINTFWW("Error: Line %" PRIuPTR " of --vcf file has an invalid %s field.\n", line_idx, dosage_import_field);
  vcf_to_pgen_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 vcf_to_pgen_ret_1:
  if (spgw_cleanup(&spgw) && (!reterr)) {
    reterr = kPglRetWriteFail;
  }
  gzclose_cond(gz_infile);
  fclose_cond(pvarfile);
  bigstack_double_reset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

pglerr_t ox_sample_to_psam(const char* samplename, const char* ox_missing_code, misc_flags_t misc_flags, char* outname, char* outname_end, uint32_t* sample_ct_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  gzFile gz_infile = nullptr;
  FILE* psamfile = nullptr;
  pglerr_t reterr = kPglRetSuccess;
  uintptr_t loadbuf_size = 0;
  uintptr_t line_idx = 0;
  {
    uint32_t omp_slen = 2;
    char output_missing_pheno[kMaxMissingPhenostrBlen];
    if (misc_flags & kfMiscKeepAutoconv) {
      // must use --output-missing-phenotype parameter, which we've validated
      // to be consistent with --input-missing-phenotype
      omp_slen = strlen(g_output_missing_pheno);
      memcpy(output_missing_pheno, g_output_missing_pheno, omp_slen);
    } else {
      // use "NA" since that's always safe
      memcpy(output_missing_pheno, "NA", 2);
    }
    const char* missing_catname = g_missing_catname;
    uint32_t missing_catname_slen = strlen(missing_catname);

    gz_infile = gzopen(samplename, FOPEN_RB);
    if (!gz_infile) {
      const uint32_t slen = strlen(samplename);
      if (str_endswith(samplename, ".sample", slen) ||
          str_endswith(samplename, ".sample.gz", slen)) {
        LOGERRPRINTFWW(g_errstr_fopen, samplename);
      } else {
        LOGERRPRINTFWW("Error: Failed to open %s. (--sample expects a complete filename; did you forget '.sample' at the end?)\n", samplename);
      }
      goto ox_sample_to_psam_ret_OPEN_FAIL;
    }
    if (gzbuffer(gz_infile, 131072)) {
      goto ox_sample_to_psam_ret_NOMEM;
    }
    uint32_t mc_ct = 0;
    uintptr_t max_mc_blen = 1;
    char* sorted_mc = nullptr;
    if (!ox_missing_code) {
      if (bigstack_alloc_c(3, &sorted_mc)) {
        goto ox_sample_to_psam_ret_NOMEM;
      }
      memcpy(sorted_mc, "NA", 3);
      mc_ct = 1;
      max_mc_blen = 3;
    } else {
      // er, this should use something like
      // count_and_measure_multistr_reverse_alloc()...
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
        uintptr_t token_slen = (uintptr_t)(token_end - missing_code_iter);
        if (token_slen >= max_mc_blen) {
          max_mc_blen = token_slen + 1;
        }
        missing_code_iter = token_end;
      }
      if (mc_ct) {
        if (bigstack_alloc_c(mc_ct * max_mc_blen, &sorted_mc)) {
          goto ox_sample_to_psam_ret_NOMEM;
        }
        missing_code_iter = ox_missing_code;
        for (uintptr_t mc_idx = 0; mc_idx < mc_ct; ++mc_idx) {
          while (*missing_code_iter == ',') {
            ++missing_code_iter;
          }
          const char* token_end = strchrnul(missing_code_iter, ',');
          uintptr_t token_slen = (uintptr_t)(token_end - missing_code_iter);
          memcpyx(&(sorted_mc[mc_idx * max_mc_blen]), missing_code_iter, token_slen, '\0');
          missing_code_iter = token_end;
        }
        qsort(sorted_mc, mc_ct, max_mc_blen, strcmp_casted);
      }
    }
    loadbuf_size = bigstack_left() / 4;
    if (loadbuf_size > kMaxLongLine) {
      loadbuf_size = kMaxLongLine;
    } else if (loadbuf_size <= kMaxMediumLine) {
      goto ox_sample_to_psam_ret_NOMEM;
    } else {
      loadbuf_size = round_up_pow2(loadbuf_size, kCacheline);
    }
    char* loadbuf = (char*)bigstack_alloc_raw(loadbuf_size);
    loadbuf[loadbuf_size - 1] = ' ';
    char* loadbuf_first_token;
    do {
      ++line_idx;
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
        if (!gzeof(gz_infile)) {
          goto ox_sample_to_psam_ret_READ_FAIL;
        }
        logerrprint("Error: Empty .sample file.\n");
        goto ox_sample_to_psam_ret_MALFORMED_INPUT;
      }
      if (!loadbuf[loadbuf_size - 1]) {
        goto ox_sample_to_psam_ret_LONG_LINE;
      }
      loadbuf_first_token = skip_initial_spaces(loadbuf);
    } while (is_eoln_kns(*loadbuf_first_token));
    char* token_end = token_endnn(loadbuf_first_token);
    if (!strequal_k(loadbuf_first_token, "ID_1", (uintptr_t)(token_end - loadbuf_first_token))) {
      goto ox_sample_to_psam_ret_INVALID_SAMPLE_HEADER_1;
    }
    // currently accepts tab as delimiter, though .sample spec technically
    // prohibits that
    char* loadbuf_iter = skip_initial_spaces(token_end);
    uint32_t token_slen = strlen_se(loadbuf_iter);
    if (!strequal_k(loadbuf_iter, "ID_2", token_slen)) {
      goto ox_sample_to_psam_ret_INVALID_SAMPLE_HEADER_1;
    }
    loadbuf_iter = skip_initial_spaces(&(loadbuf_iter[token_slen]));
    token_slen = strlen_se(loadbuf_iter);
    if (!match_upper_k(loadbuf_iter, "MISSING", token_slen)) {
      goto ox_sample_to_psam_ret_INVALID_SAMPLE_HEADER_1;
    }
    loadbuf_iter = skip_initial_spaces(&(loadbuf_iter[token_slen]));

    strcpy(outname_end, ".psam");
    if (fopen_checked(outname, FOPEN_WB, &psamfile)) {
      goto ox_sample_to_psam_ret_OPEN_FAIL;
    }
    // categorical phenotypes are lengthened by 1 character ('C' added in
    // front), so this needs to be 50% larger than loadbuf to handle worst case
    char* writebuf = (char*)bigstack_alloc_raw_rd(loadbuf_size + (loadbuf_size / 2));
    char* write_iter = strcpya(writebuf, "#FID\tIID\tSEX");

    // 0 = not present, otherwise zero-based index (this is fine since first
    //     column has to be FID)
    uint32_t sex_col = 0;

    uint32_t col_ct = 3;

    while (!is_eoln_kns(*loadbuf_iter)) {
      token_end = token_endnn(loadbuf_iter);
      token_slen = (uintptr_t)(token_end - loadbuf_iter);
      if (match_upper_k(loadbuf_iter, "SEX", token_slen)) {
        if (sex_col) {
          logerrprint("Error: Multiple sex columns in .sample file.\n");
          goto ox_sample_to_psam_ret_MALFORMED_INPUT;
        }
        sex_col = col_ct;
      }
      ++col_ct;
      loadbuf_iter = skip_initial_spaces(token_end);
    }

    do {
      ++line_idx;
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
        if (!gzeof(gz_infile)) {
          goto ox_sample_to_psam_ret_READ_FAIL;
        }
        logerrprint("Error: Only one nonempty line in .sample file.\n");
        goto ox_sample_to_psam_ret_MALFORMED_INPUT;
      }
      if (!loadbuf[loadbuf_size - 1]) {
        goto ox_sample_to_psam_ret_LONG_LINE;
      }
      loadbuf_first_token = skip_initial_spaces(loadbuf);
    } while (is_eoln_kns(*loadbuf_first_token));

    token_end = token_endnn(loadbuf_first_token);
    if ((((uintptr_t)(token_end - loadbuf_first_token)) != 1) || (*loadbuf_first_token != '0')) {
      goto ox_sample_to_psam_ret_INVALID_SAMPLE_HEADER_2;
    }
    loadbuf_iter = skip_initial_spaces(token_end);
    token_slen = strlen_se(loadbuf_iter);
    if ((token_slen != 1) || (*loadbuf_iter != '0')) {
      goto ox_sample_to_psam_ret_INVALID_SAMPLE_HEADER_2;
    }
    loadbuf_iter = skip_initial_spaces(&(loadbuf_iter[1]));
    token_slen = strlen_se(loadbuf_iter);
    if ((token_slen != 1) || (*loadbuf_iter != '0')) {
      goto ox_sample_to_psam_ret_INVALID_SAMPLE_HEADER_2;
    }
    loadbuf_iter = &(loadbuf_iter[1]);

    const uint32_t col_ctl = BITCT_TO_WORDCT(col_ct);
    uintptr_t* col_is_categorical = (uintptr_t*)bigstack_alloc_raw_rd(col_ctl * sizeof(intptr_t));
    uintptr_t* col_is_qt = (uintptr_t*)bigstack_alloc_raw_rd(col_ctl * sizeof(intptr_t));
    fill_ulong_zero(col_ctl, col_is_categorical);
    fill_ulong_zero(col_ctl, col_is_qt);
    uint32_t at_least_one_binary_pheno = 0;
    for (uint32_t col_idx = 3; col_idx < col_ct; ++col_idx) {
      loadbuf_iter = skip_initial_spaces(loadbuf_iter);
      unsigned char col_type_char = *loadbuf_iter;
      if (is_eoln_kns(col_type_char)) {
        logerrprint("Error: Second .sample header line has fewer tokens than the first.\n");
        goto ox_sample_to_psam_ret_MALFORMED_INPUT;
      }
      if (loadbuf_iter[1] > ' ') {
        goto ox_sample_to_psam_ret_INVALID_SAMPLE_HEADER_2;
      }
      if (col_idx == sex_col) {
        if (col_type_char != 'D') {
          logerrprint("Error: .sample sex column is not of type 'D'.\n");
          goto ox_sample_to_psam_ret_MALFORMED_INPUT;
        }
      } else {
        if ((col_type_char == 'C') || (col_type_char == 'P')) {
          set_bit(col_idx, col_is_qt);
        } else if (col_type_char == 'D') {
          set_bit(col_idx, col_is_categorical);
        } else {
          at_least_one_binary_pheno = 1;
          if (col_type_char != 'B') {
            sprintf(g_logbuf, "Error: Unrecognized .sample variable type '%c'.\n", col_type_char);
            goto ox_sample_to_psam_ret_MALFORMED_INPUT_2;
          }
        }
      }
      ++loadbuf_iter;
    }
    if (at_least_one_binary_pheno) {
      // check for pathological case
      if ((bsearch_str("0", sorted_mc, 1, max_mc_blen, mc_ct) != -1) || (bsearch_str("1", sorted_mc, 1, max_mc_blen, mc_ct) != -1)) {
        logerrprint("Error: '0' and '1' are unacceptable missing case/control phenotype codes.\n");
        goto ox_sample_to_psam_ret_INCONSISTENT_INPUT;
      }
    }
    // to make --data and --data --make-pgen consistent, we do a two-pass load,
    // checking for empty phenotypes in the first pass.
    uintptr_t* col_keep = (uintptr_t*)bigstack_alloc_raw(round_up_pow2(col_ct, kCacheline * CHAR_BIT) / CHAR_BIT);
    col_keep[0] = 7;
    fill_ulong_zero(col_ctl - 1, &(col_keep[1]));
    uint32_t uncertain_col_ct = col_ct - 3;
    if (sex_col) {
      // we don't care if sex column is all-NA
      set_bit(sex_col, col_keep);
      --uncertain_col_ct;
    }
    while (uncertain_col_ct) {
      ++line_idx;
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
        break;
      }
      if (!loadbuf[loadbuf_size - 1]) {
        goto ox_sample_to_psam_ret_LONG_LINE;
      }
      loadbuf_iter = skip_initial_spaces(loadbuf);
      if (is_eoln_kns(*loadbuf_iter)) {
        continue;
      }

      const uint32_t old_uncertain_col_ct = uncertain_col_ct;
      uint32_t old_col_uidx = 0;
      uint32_t col_uidx = 0;
      for (uint32_t uncertain_col_idx = 0; uncertain_col_idx < old_uncertain_col_ct; ++uncertain_col_idx, ++col_uidx) {
        next_unset_unsafe_ck(col_keep, &col_uidx);
        loadbuf_iter = next_token_mult(loadbuf_iter, col_uidx - old_col_uidx);
        if (!loadbuf_iter) {
          goto ox_sample_to_psam_ret_MISSING_TOKENS;
        }
        token_end = token_endnn(loadbuf_iter);
        token_slen = (uintptr_t)(token_end - loadbuf_iter);
        if (bsearch_str(loadbuf_iter, sorted_mc, token_slen, max_mc_blen, mc_ct) == -1) {
          set_bit(col_uidx, col_keep);
          --uncertain_col_ct;
        }
        loadbuf_iter = token_end;
        old_col_uidx = col_uidx;
      }
    }

    if (gzrewind(gz_infile)) {
      goto ox_sample_to_psam_ret_READ_FAIL;
    }
    line_idx = 0;

    uint32_t sample_ct_p2 = 0;
    while (1) {
      ++line_idx;
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
        break;
      }
      if (!loadbuf[loadbuf_size - 1]) {
        goto ox_sample_to_psam_ret_LONG_LINE;
      }
      loadbuf_first_token = skip_initial_spaces(loadbuf);
      if (is_eoln_kns(*loadbuf_first_token)) {
        continue;
      }
      ++sample_ct_p2;
      if (sample_ct_p2 < 3) {
        // header lines
        if (sample_ct_p2 == 1) {
          loadbuf_iter = next_token_mult(loadbuf_first_token, 3);
          for (uint32_t col_idx = 3; col_idx < col_ct; ++col_idx) {
            token_end = token_endnn(loadbuf_iter);
            if (is_set(col_keep, col_idx) && (col_idx != sex_col)) {
              *write_iter++ = '\t';
              write_iter = memcpya(write_iter, loadbuf_iter, (uintptr_t)(token_end - loadbuf_iter));
            }
            loadbuf_iter = skip_initial_spaces(token_end);
          }
          append_binary_eoln(&write_iter);
          if (fwrite_checked(writebuf, (uintptr_t)(write_iter - writebuf), psamfile)) {
            goto ox_sample_to_psam_ret_WRITE_FAIL;
          }
          write_iter = writebuf;
        }
        continue;
      }
      if (sample_ct_p2 == 0x80000001U) {
        logerrprint("Error: " PROG_NAME_STR " does not support more than 2^31 - 2 samples.\n");
        goto ox_sample_to_psam_ret_MALFORMED_INPUT;
      }

      // FID
      token_end = token_endnn(loadbuf_first_token);
      write_iter = memcpyax(writebuf, loadbuf_first_token, (uintptr_t)(token_end - loadbuf_first_token), '\t');

      // IID
      loadbuf_iter = skip_initial_spaces(token_end);
      if (is_eoln_kns(*loadbuf_iter)) {
        goto ox_sample_to_psam_ret_MISSING_TOKENS;
      }
      token_end = token_endnn(loadbuf_iter);
      write_iter = memcpya(write_iter, loadbuf_iter, (uintptr_t)(token_end - loadbuf_iter));

      // MISSING
      loadbuf_iter = skip_initial_spaces(token_end);
      if (is_eoln_kns(*loadbuf_iter)) {
        goto ox_sample_to_psam_ret_MISSING_TOKENS;
      }
      token_end = token_endnn(loadbuf_iter);

      // flush now since backfilled sex is variable-length ("NA" vs. "1"/"2")
      if (fwrite_checked(writebuf, (uintptr_t)(write_iter - writebuf), psamfile)) {
        goto ox_sample_to_psam_ret_WRITE_FAIL;
      }
      char* cur_writebuf_start = writebuf;
      write_iter = memcpyl3a(writebuf, "\tNA");
      for (uint32_t col_idx = 3; col_idx < col_ct; ++col_idx) {
        loadbuf_iter = skip_initial_spaces(token_end);
        if (is_eoln_kns(*loadbuf_iter)) {
          goto ox_sample_to_psam_ret_MISSING_TOKENS;
        }
        token_end = token_endnn(loadbuf_iter);
        if (!is_set(col_keep, col_idx)) {
          continue;
        }
        token_slen = (uintptr_t)(token_end - loadbuf_iter);
        const uint32_t is_missing = (bsearch_str(loadbuf_iter, sorted_mc, token_slen, max_mc_blen, mc_ct) != -1);
        if (col_idx == sex_col) {
          if (!is_missing) {
            const unsigned char sex_ucc = *loadbuf_iter;
            if ((token_slen == 1) && ((((uint32_t)sex_ucc) - 49) < 2)) {
              ++cur_writebuf_start;
              cur_writebuf_start[0] = '\t';
              cur_writebuf_start[1] = sex_ucc;
            } else if ((token_slen != 1) || (sex_ucc != '0')) {
              // tolerate '0' as a sex-only missing code even when not
              // explicitly specified
              *token_end = '\0';
              sprintf(g_logbuf, "Error: Invalid sex code '%s' on line %" PRIuPTR ", column %u of .sample file ('0', '1', '2', or --missing-code value expected).\n", loadbuf_iter, line_idx, col_idx + 1);
              goto ox_sample_to_psam_ret_INCONSISTENT_INPUT_WW;
            }
          }
        } else {
          *write_iter++ = '\t';
          if (is_set(col_is_categorical, col_idx)) {
            if (!is_missing) {
              *write_iter++ = 'C';
              // .sample files are relatively small, so let's go ahead and
              // (i) validate we have a positive integer < 2^31
              // (ii) convert e.g. 9000000, 9000000., 9.0e6 all to 9000000
              double dxx = 0.0;
              char* num_end = scanadv_double(loadbuf_iter, &dxx);
              int32_t ii = (int32_t)dxx;
              if ((num_end != token_end) || (ii <= 0) || (((double)ii) != dxx)) {
                *token_end = '\0';
                sprintf(g_logbuf, "Error: Invalid categorical phenotype '%s' on line %" PRIuPTR ", column %u of .sample file (positive integer < 2^31 or --missing-code value expected).\n", loadbuf_iter, line_idx, col_idx + 1);
                goto ox_sample_to_psam_ret_INCONSISTENT_INPUT_WW;
              }
              write_iter = uint32toa(ii, write_iter);
            } else {
              write_iter = memcpya(write_iter, missing_catname, missing_catname_slen);
            }
          } else if (!is_missing) {
            if (is_set(col_is_qt, col_idx)) {
              double dxx = 0.0;
              char* num_end = scanadv_double(loadbuf_iter, &dxx);
              if (num_end != token_end) {
                *token_end = '\0';
                sprintf(g_logbuf, "Error: Invalid quantitative phenotype '%s' on line %" PRIuPTR ", column %u of .sample file (non-infinite number or --missing-code value expected).\n", loadbuf_iter, line_idx, col_idx + 1);
                goto ox_sample_to_psam_ret_INCONSISTENT_INPUT_WW;
              }
              // used over memcpy to make --data and --data --make-pgen the
              // same (could make that conditional on keep_autoconv?)
              write_iter = dtoa_g(dxx, write_iter);
            } else {
              const uint32_t cc_char_m48 = (uint32_t)((unsigned char)(*loadbuf_iter)) - 48;
              if ((token_slen == 1) && (cc_char_m48 < 2)) {
                *write_iter++ = cc_char_m48 + '1';
              } else {
                *token_end = '\0';
                sprintf(g_logbuf, "Error: Invalid binary phenotype value '%s' on line %" PRIuPTR ", column %u of .sample file ('0', '1', or --missing-code value expected).\n", loadbuf_iter, line_idx, col_idx + 1);
                goto ox_sample_to_psam_ret_INCONSISTENT_INPUT_WW;
              }
            }
          } else {
            write_iter = memcpya(write_iter, output_missing_pheno, omp_slen);
          }
        }
      }
      append_binary_eoln(&write_iter);
      if (fwrite_checked(cur_writebuf_start, (uintptr_t)(write_iter - cur_writebuf_start), psamfile)) {
        goto ox_sample_to_psam_ret_WRITE_FAIL;
      }
    }
    if (!gzeof(gz_infile)) {
      goto ox_sample_to_psam_ret_READ_FAIL;
    }

    // no final writebuf flush since we didn't use usual manual-streaming
    // strategy
    if (fclose_null(&psamfile)) {
      goto ox_sample_to_psam_ret_WRITE_FAIL;
    }
    const uint32_t sample_ct = sample_ct_p2 - 2;
    if ((!sample_ct) && (!(misc_flags & kfMiscAllowNoSamples))) {
      logerrprint("Error: No samples in .sample file.\n");
      goto ox_sample_to_psam_ret_INCONSISTENT_INPUT;
    }
    LOGPRINTFWW("%u sample%s imported from .sample file to %s .\n", sample_ct, (sample_ct == 1)? "" : "s", outname);
    *sample_ct_ptr = sample_ct;
  }
  while (0) {
  ox_sample_to_psam_ret_LONG_LINE:
    if (loadbuf_size == kMaxLongLine) {
      LOGERRPRINTF("Error: Line %" PRIuPTR " of .sample file is pathologically long.\n", line_idx);
      reterr = kPglRetMalformedInput;
      break;
    }
  ox_sample_to_psam_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ox_sample_to_psam_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  ox_sample_to_psam_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  ox_sample_to_psam_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  ox_sample_to_psam_ret_MISSING_TOKENS:
    LOGERRPRINTF("Error: Line %" PRIuPTR " of .sample file has fewer tokens than expected.\n", line_idx);
    reterr = kPglRetMalformedInput;
    break;
  ox_sample_to_psam_ret_MALFORMED_INPUT_2:
    logerrprintb();
  ox_sample_to_psam_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  ox_sample_to_psam_ret_INVALID_SAMPLE_HEADER_1:
    logerrprint("Error: Invalid first header line in .sample file.\n");
    reterr = kPglRetMalformedInput;
    break;
  ox_sample_to_psam_ret_INVALID_SAMPLE_HEADER_2:
    logerrprint("Error: Invalid second header line in .sample file.\n");
    reterr = kPglRetMalformedInput;
    break;
  ox_sample_to_psam_ret_INCONSISTENT_INPUT_WW:
    wordwrapb(0);
    logerrprintb();
  ox_sample_to_psam_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
  gzclose_cond(gz_infile);
  fclose_cond(psamfile);
  bigstack_reset(bigstack_mark);
  return reterr;
}

uint32_t bgen11_dosage_import_check(uint32_t dosage_int_sum_thresh, uint32_t import_dosage_certainty_int, uint32_t dosage_erase_halfdist, uint32_t dosage_int0, uint32_t dosage_int1, uint32_t dosage_int2) {
  const uint32_t dosage_int_sum = dosage_int0 + dosage_int1 + dosage_int2;
  if ((dosage_int_sum <= dosage_int_sum_thresh) && (dosage_int0 < import_dosage_certainty_int) && (dosage_int1 < import_dosage_certainty_int) && (dosage_int2 < import_dosage_certainty_int)) {
    return 0;
  }
  // ties realistically happen, use banker's rounding
  // 1/65536 -> 0/32768
  // 3/65536 -> 2/32768
  // 5/65536 -> 2/32768
  const dosage_prod_t write_dosage_int_numer = ((dosage_prod_t)kDosageMid) * dosage_int1 + ((dosage_prod_t)kDosageMax) * dosage_int2;
  uint32_t write_dosage_int;
  if (dosage_int_sum == kDosageMax) {
    // optimize common case
    write_dosage_int = ((write_dosage_int_numer + kDosageMid) / kDosageMax) - ((write_dosage_int_numer % (2 * ((dosage_prod_t)kDosageMax))) == kDosageMid);
  } else {
    write_dosage_int = (write_dosage_int_numer + (dosage_int_sum / 2)) / dosage_int_sum;
    write_dosage_int -= (2 * (write_dosage_int_numer - write_dosage_int * dosage_int_sum) == dosage_int_sum) * (write_dosage_int % 2);
  }
  const uint32_t halfdist = biallelic_dosage_halfdist(write_dosage_int);
  return (halfdist < dosage_erase_halfdist);
}

void bgen11_dosage_import_update(uint32_t dosage_int_sum_thresh, uint32_t import_dosage_certainty_int, uint32_t hard_call_halfdist, uint32_t dosage_erase_halfdist, uint32_t sample_idx_lowbits, uint32_t dosage_int0, uint32_t dosage_int1, uint32_t dosage_int2, uintptr_t* genovec_word_ptr, uint32_t* dosage_present_hw_ptr, dosage_t** dosage_vals_iterp) {
  const uint32_t dosage_int_sum = dosage_int0 + dosage_int1 + dosage_int2;
  if (dosage_int_sum <= dosage_int_sum_thresh) {
    if ((dosage_int0 < import_dosage_certainty_int) && (dosage_int1 < import_dosage_certainty_int) && (dosage_int2 < import_dosage_certainty_int)) {
      *genovec_word_ptr |= (3 * k1LU) << (2 * sample_idx_lowbits);
      return;
    }
  }
  const dosage_prod_t write_dosage_int_numer = ((dosage_prod_t)kDosageMid) * dosage_int1 + ((dosage_prod_t)kDosageMax) * dosage_int2;
  uint32_t write_dosage_int;
  if (dosage_int_sum == kDosageMax) {
    write_dosage_int = ((write_dosage_int_numer + kDosageMid) / kDosageMax) - ((write_dosage_int_numer % (2 * ((dosage_prod_t)kDosageMax))) == kDosageMid);
  } else {
    write_dosage_int = (write_dosage_int_numer + (dosage_int_sum / 2)) / dosage_int_sum;
    write_dosage_int -= (2 * (write_dosage_int_numer - write_dosage_int * dosage_int_sum) == dosage_int_sum) * (write_dosage_int % 2);
  }
  const uint32_t halfdist = biallelic_dosage_halfdist(write_dosage_int);
  if (halfdist < hard_call_halfdist) {
    *genovec_word_ptr |= (3 * k1LU) << (2 * sample_idx_lowbits);
  } else {
    *genovec_word_ptr |= ((write_dosage_int + (kDosage4th * k1LU)) / kDosageMid) << (2 * sample_idx_lowbits);
    if (halfdist >= dosage_erase_halfdist) {
      return;
    }
  }
  *dosage_present_hw_ptr |= 1U << sample_idx_lowbits;
  **dosage_vals_iterp = write_dosage_int;
  *dosage_vals_iterp += 1;
}

static_assert(sizeof(dosage_t) == 2, "ox_gen_to_pgen() needs to be updated.");
pglerr_t ox_gen_to_pgen(const char* genname, const char* samplename, const char* ox_single_chr_str, const char* ox_missing_code, misc_flags_t misc_flags, oxford_import_t oxford_import_flags, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, double import_dosage_certainty, char* outname, char* outname_end, chr_info_t* cip) {
  unsigned char* bigstack_mark = g_bigstack_base;
  gzFile gz_infile = nullptr;
  FILE* pvarfile = nullptr;
  st_pgen_writer_t spgw;
  pglerr_t reterr = kPglRetSuccess;
  uintptr_t loadbuf_size = 0;
  uintptr_t line_idx = 0;
  spgw_preinit(&spgw);
  {
    uint32_t sample_ct;
    reterr = ox_sample_to_psam(samplename, ox_missing_code, misc_flags, outname, outname_end, &sample_ct);
    if (reterr) {
      goto ox_gen_to_pgen_ret_1;
    }
    if (sample_ct > (kMaxLongLine / 6)) {
      // impossible for a valid .gen line to fit in maximum-length load buffer
      logerrprint("Error: Too many samples for .gen file converter.\n");
      reterr = kPglRetNotYetSupported;
      goto ox_gen_to_pgen_ret_1;
    }
    // Two passes:
    // 1. Count # of (non-chromosome-filtered) variants, write .pvar file, and
    //    check if at least one non-hardcall needs to be saved.
    // 2. Write the .pgen.
    gz_infile = gzopen(genname, FOPEN_RB);
    if (!gz_infile) {
      const uint32_t slen = strlen(genname);
      if (str_endswith(genname, ".gen", slen) ||
          str_endswith(genname, ".gen.gz", slen)) {
        LOGERRPRINTFWW(g_errstr_fopen, genname);
      } else {
        LOGERRPRINTFWW("Error: Failed to open %s. (--gen expects a complete filename; did you forget '.gen' at the end?)\n", genname);
      }
      goto ox_gen_to_pgen_ret_OPEN_FAIL;
    }
    if (gzbuffer(gz_infile, 131072)) {
      goto ox_gen_to_pgen_ret_NOMEM;
    }
    loadbuf_size = bigstack_left() / 4;
    if (loadbuf_size > kMaxLongLine) {
      loadbuf_size = kMaxLongLine;
    } else if (loadbuf_size <= kMaxMediumLine) {
      goto ox_gen_to_pgen_ret_NOMEM;
    } else {
      loadbuf_size = round_up_pow2(loadbuf_size, kCacheline);
    }
    char* loadbuf = (char*)bigstack_alloc_raw(loadbuf_size);
    loadbuf[loadbuf_size - 1] = ' ';
    const uint32_t allow_extra_chrs = (misc_flags / kfMiscAllowExtraChrs) & 1;
    finalize_chrset(misc_flags, cip);

    char* writebuf = (char*)bigstack_alloc_raw(kMaxMediumLine + loadbuf_size);
    char* writebuf_flush = &(writebuf[kMaxMediumLine]);

    const char* single_chr_str = nullptr;
    uint32_t single_chr_slen = 0;
    if (ox_single_chr_str) {
      int32_t chr_code_raw = get_chr_code_raw(ox_single_chr_str);
      if (chr_code_raw == -1) {
        // command-line parser guarantees that allow_extra_chrs is true here
        single_chr_str = ox_single_chr_str;
        single_chr_slen = strlen(ox_single_chr_str);
      } else {
        uint32_t chr_code = chr_code_raw;
        if (chr_code > cip->max_code) {
          if (chr_code < kMaxContigs) {
            logerrprint("Error: --oxford-single-chr chromosome code is not in the chromosome set.\n");
            goto ox_gen_to_pgen_ret_INVALID_CMDLINE;
          }
          chr_code = cip->xymt_codes[chr_code - kMaxContigs];
          if (((int32_t)chr_code) < 0) {
            logerrprint("Error: --oxford-single-chr chromosome code is not in the chromosome set.\n");
            goto ox_gen_to_pgen_ret_INVALID_CMDLINE;
          }
        }
        if (!is_set(cip->chr_mask, chr_code)) {
          // could permit this in --allow-no-vars case, but it's silly
          logerrprint("Error: --oxford-single-chr chromosome code is excluded by chromosome filter.\n");
          goto ox_gen_to_pgen_ret_INVALID_CMDLINE;
        }
        char* chr_buf = (char*)bigstack_alloc_raw(kCacheline);
        char* chr_name_end = chr_name_write(cip, chr_code, chr_buf);
        single_chr_str = chr_buf;
        single_chr_slen = (uintptr_t)(chr_name_end - chr_buf);
      }
    }

    strcpy(outname_end, ".pvar");
    if (fopen_checked(outname, FOPEN_WB, &pvarfile)) {
      goto ox_gen_to_pgen_ret_OPEN_FAIL;
    }
    char* write_iter = writebuf;
    if (cip->chrset_source) {
      append_chrset_line(cip, &write_iter);
    }
    write_iter = strcpya(write_iter, "#CHROM\tPOS\tID\tREF\tALT" EOLN_STR);

    // Explicit 32768 instead of kDosageMax since this is driven by the BGEN
    // 1.1 format, not plink2's dosage representation.
    // Note that command-line parser multiplies import_dosage_certainty by
    // (1 - kSmallEpsilon), and we want import_dosage_certainty_int to be 1
    // when import_dosage_certainty is zero.
    uint32_t import_dosage_certainty_int = 1 + (int32_t)(import_dosage_certainty * 32768);
    const uint32_t dosage_int_sum_thresh = 3 * (import_dosage_certainty_int - 1);

    const uint32_t dosage_erase_halfdist = kDosage4th - dosage_erase_thresh;
    const uint32_t prov_ref_allele_second = !(oxford_import_flags & kfOxfordImportRefFirst);
    uint32_t dosage_is_present = 0;
    uint32_t variant_ct = 0;
    uintptr_t variant_skip_ct = 0;
    while (1) {
      ++line_idx;
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
        if (!gzeof(gz_infile)) {
          goto ox_gen_to_pgen_ret_READ_FAIL;
        }
        break;
      }
      if (!loadbuf[loadbuf_size - 1]) {
        goto ox_gen_to_pgen_ret_LONG_LINE_N;
      }
      char* chr_code_str = skip_initial_spaces(loadbuf);
      if (is_eoln_kns(*chr_code_str)) {
        continue;
      }
      char* chr_code_end = token_endnn(chr_code_str);
      char* variant_id_str = skip_initial_spaces(chr_code_end);
      if (is_eoln_kns(*variant_id_str)) {
        goto ox_gen_to_pgen_ret_MISSING_TOKENS;
      }

      if (!single_chr_str) {
        int32_t cur_chr_code;
        reterr = get_or_add_chr_code_destructive(".gen file", line_idx, allow_extra_chrs, chr_code_str, chr_code_end, cip, &cur_chr_code);
        if (reterr) {
          if (strequal_k(chr_code_str, "---", (uintptr_t)(chr_code_end - chr_code_str))) {
            logerrprint("(Did you forget --oxford-single-chr?)\n");
          }
          goto ox_gen_to_pgen_ret_1;
        }
        if (!is_set(cip->chr_mask, cur_chr_code)) {
          ++variant_skip_ct;
          continue;
        }
        write_iter = chr_name_write(cip, cur_chr_code, write_iter);
      } else {
        write_iter = memcpya(write_iter, single_chr_str, single_chr_slen);
      }
      *write_iter++ = '\t';
      ++variant_ct;

      char* variant_id_end = token_endnn(variant_id_str);
      char* pos_str = skip_initial_spaces(variant_id_end);
      if (is_eoln_kns(*pos_str)) {
        goto ox_gen_to_pgen_ret_MISSING_TOKENS;
      }
      char* pos_end = token_endnn(pos_str);
      uint32_t cur_bp;
      if (scan_uint_defcap(pos_str, &cur_bp)) {
        putc_unlocked('\n', stdout);
        sprintf(g_logbuf, "Error: Invalid bp coordinate on line %" PRIuPTR " of %s.\n", line_idx, genname);
        goto ox_gen_to_pgen_ret_MALFORMED_INPUT_WW;
      }

      write_iter = uint32toa_x(cur_bp, '\t', write_iter);
      const uint32_t variant_id_slen = (uintptr_t)(variant_id_end - variant_id_str);
      if (variant_id_slen > kMaxIdSlen) {
        putc_unlocked('\n', stdout);
        logerrprint("Error: Variant names are limited to " MAX_ID_SLEN_STR " characters.\n");
        goto ox_gen_to_pgen_ret_MALFORMED_INPUT;
      }
      write_iter = memcpyax(write_iter, variant_id_str, variant_id_slen, '\t');

      // .gen specification does not define which column should be expected to
      // be the reference allele, and which the alternate.  plink 1.9 assumed
      // alt was usually first, but the reverse seems to be more common now.
      // So:
      //   If 'ref-first' or 'ref-second' was specified, we know what to do.
      //   If not, we treat the second allele as the provisional reference, for
      //     backward compatibility.
      char* first_allele_str = skip_initial_spaces(pos_end);
      if (is_eoln_kns(*first_allele_str)) {
        goto ox_gen_to_pgen_ret_MISSING_TOKENS;
      }
      char* first_allele_end = token_endnn(first_allele_str);
      char* second_allele_str = skip_initial_spaces(first_allele_end);
      if (is_eoln_kns(*second_allele_str)) {
        goto ox_gen_to_pgen_ret_MISSING_TOKENS;
      }
      char* loadbuf_iter = token_endnn(second_allele_str);
      if (!prov_ref_allele_second) {
        write_iter = memcpyax(write_iter, first_allele_str, (uintptr_t)(first_allele_end - first_allele_str), '\t');
        write_iter = memcpya(write_iter, second_allele_str, (uintptr_t)(loadbuf_iter - second_allele_str));
      } else {
        write_iter = memcpyax(write_iter, second_allele_str, (uintptr_t)(loadbuf_iter - second_allele_str), '\t');
        write_iter = memcpya(write_iter, first_allele_str, (uintptr_t)(first_allele_end - first_allele_str));
      }
      append_binary_eoln(&write_iter);
      if (fwrite_ck(writebuf_flush, pvarfile, &write_iter)) {
        goto ox_gen_to_pgen_ret_WRITE_FAIL;
      }

      if (!dosage_is_present) {
        for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
          loadbuf_iter = skip_initial_spaces(loadbuf_iter);
          const char cc = *loadbuf_iter;
          if (is_eoln_kns(cc)) {
            goto ox_gen_to_pgen_ret_MISSING_TOKENS;
          }
          char cc2 = loadbuf_iter[1];
          if ((cc2 == ' ') || (cc2 == '\t')) {
            // fast handling of "1 0 0 ", "0 1 0 ", "0 0 1 ", "0 0 0 " cases
            cc2 = loadbuf_iter[3];
            if (((cc2 == ' ') || (cc2 == '\t')) && ((unsigned char)(loadbuf_iter[5]) <= 32)) {
              const uint32_t uii = ((uint32_t)((unsigned char)cc)) - 48;
              const uint32_t ujj = ((uint32_t)((unsigned char)loadbuf_iter[2])) - 48;
              const uint32_t ukk = ((uint32_t)((unsigned char)loadbuf_iter[4])) - 48;
              if (((uii | ujj | ukk) < 2) && (uii + ujj + ukk < 2)) {
                loadbuf_iter = &(loadbuf_iter[5]);
                continue;
              }
            }
          }
          double prob_0alt;
          char* first_dosage_str_end = scanadv_double(loadbuf_iter, &prob_0alt);
          if (!first_dosage_str_end) {
            // triple-NA, etc. ok; treat as missing value
            loadbuf_iter = next_token_mult(loadbuf_iter, 2);
            if (!loadbuf_iter) {
              goto ox_gen_to_pgen_ret_MISSING_TOKENS;
            }
            loadbuf_iter = token_endnn(loadbuf_iter);
            continue;
          }
          if ((unsigned char)(*first_dosage_str_end) > ' ') {
            goto ox_gen_to_pgen_ret_INVALID_DOSAGE;
          }
          loadbuf_iter = skip_initial_spaces(first_dosage_str_end);
          if (is_eoln_kns(*loadbuf_iter)) {
            goto ox_gen_to_pgen_ret_MISSING_TOKENS;
          }
          double prob_1alt;
          loadbuf_iter = scanadv_double(loadbuf_iter, &prob_1alt);
          if ((!loadbuf_iter) || ((unsigned char)(*loadbuf_iter) > ' ')) {
            goto ox_gen_to_pgen_ret_INVALID_DOSAGE;
          }
          loadbuf_iter = skip_initial_spaces(loadbuf_iter);
          if (is_eoln_kns(*loadbuf_iter)) {
            goto ox_gen_to_pgen_ret_MISSING_TOKENS;
          }
          double prob_2alt;
          loadbuf_iter = scanadv_double(loadbuf_iter, &prob_2alt);
          if ((!loadbuf_iter) || ((unsigned char)(*loadbuf_iter) > ' ')) {
            goto ox_gen_to_pgen_ret_INVALID_DOSAGE;
          }
          // bugfix: forgot the "multiply by 32768" part of "multiply by 32768
          // and round" .gen -> .bgen conversion.
          prob_0alt *= 32768;
          prob_1alt *= 32768;
          prob_2alt *= 32768;

          // now treat this identically to bgen-1.1
          // Compare with 65535.4999999999 instead of 65535.5 since 0.5 +
          // [first floating point number below 65535.5] may evaluate to 65536.
          if ((prob_0alt < 0.0) || (prob_0alt >= 65535.4999999999) || (prob_1alt < 0.0) || (prob_1alt >= 65535.4999999999) || (prob_2alt < 0.0) || (prob_2alt >= 65535.4999999999)) {
            goto ox_gen_to_pgen_ret_INVALID_DOSAGE;
          }
          const uint32_t dosage_int0 = (int32_t)(prob_0alt + 0.5);
          const uint32_t dosage_int1 = (int32_t)(prob_1alt + 0.5);
          const uint32_t dosage_int2 = (int32_t)(prob_2alt + 0.5);
          dosage_is_present = bgen11_dosage_import_check(dosage_int_sum_thresh, import_dosage_certainty_int, dosage_erase_halfdist, dosage_int0, dosage_int1, dosage_int2);
          if (dosage_is_present) {
            break;
          }
        }
      }
      if (!(variant_ct % 1000)) {
        printf("\r--data/--gen: %uk variants scanned.", variant_ct / 1000);
        fflush(stdout);
      }
    }
    putc_unlocked('\r', stdout);
    if (fclose_flush_null(writebuf_flush, write_iter, &pvarfile)) {
      goto ox_gen_to_pgen_ret_WRITE_FAIL;
    }
    if (!variant_ct) {
      if (!variant_skip_ct) {
        // permit this in --allow-no-vars case?
        logerrprint("Error: Empty .gen file.\n");
        goto ox_gen_to_pgen_ret_INCONSISTENT_INPUT;
      }
      LOGERRPRINTFWW("Error: All %" PRIuPTR " variant%s in .gen file skipped due to chromosome filter.\n", variant_skip_ct, (variant_skip_ct == 1)? "" : "s");
      goto ox_gen_to_pgen_ret_INCONSISTENT_INPUT;
    }
    LOGPRINTF("--data/--gen: %u variant%s scanned%s.\n", variant_ct, (variant_ct == 1)? "" : "s", dosage_is_present? "" : " (all hardcalls)");

    // second pass
    bigstack_reset(writebuf);
    if (gzrewind(gz_infile)) {
      goto ox_gen_to_pgen_ret_READ_FAIL;
    }
    strcpy(outname_end, ".pgen");
    uintptr_t spgw_alloc_cacheline_ct;
    uint32_t max_vrec_len;
    reterr = spgw_init_phase1(outname, nullptr, nullptr, variant_ct, sample_ct, dosage_is_present? kfPgenGlobalDosagePresent : kfPgenGlobal0, (oxford_import_flags & (kfOxfordImportRefFirst | kfOxfordImportRefSecond))? 1 : 2, &spgw, &spgw_alloc_cacheline_ct, &max_vrec_len);
    if (reterr) {
      goto ox_gen_to_pgen_ret_1;
    }
    unsigned char* spgw_alloc;
    if (bigstack_alloc_uc(spgw_alloc_cacheline_ct * kCacheline, &spgw_alloc)) {
      goto ox_gen_to_pgen_ret_NOMEM;
    }
    spgw_init_phase2(max_vrec_len, &spgw, spgw_alloc);

    const uint32_t sample_ctl2 = QUATERCT_TO_WORDCT(sample_ct);
    const uint32_t sample_ctl = BITCT_TO_WORDCT(sample_ct);
    uintptr_t* genovec;
    uintptr_t* dosage_present;
    // if we weren't using bigstack_alloc, this would need to be sample_ctaw2
    if (bigstack_alloc_ul(sample_ctl2, &genovec) ||
        bigstack_alloc_ul(sample_ctl, &dosage_present)) {
      goto ox_gen_to_pgen_ret_NOMEM;
    }
    dosage_t* dosage_vals = nullptr;
    if (dosage_is_present) {
      if (bigstack_alloc_dosage(sample_ct, &dosage_vals)) {
        goto ox_gen_to_pgen_ret_NOMEM;
      }
    }
    if (hard_call_thresh == UINT32_MAX) {
      hard_call_thresh = kDosageMid / 10;
    }
    const uint32_t hard_call_halfdist = kDosage4th - hard_call_thresh;
    const uint32_t sample_ctl2_m1 = sample_ctl2 - 1;
    const uintptr_t line_ct = line_idx - 1;
    uint32_t vidx = 0;
    for (line_idx = 1; line_idx <= line_ct; ++line_idx) {
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
        goto ox_gen_to_pgen_ret_READ_FAIL;
      }
      char* chr_code_str = skip_initial_spaces(loadbuf);
      if (is_eoln_kns(*chr_code_str)) {
        continue;
      }
      char* chr_code_end = token_endnn(chr_code_str);
      char* loadbuf_iter = skip_initial_spaces(chr_code_end);
      if (variant_skip_ct) {
        *chr_code_end = '\0';
        const uint32_t chr_code = get_chr_code(chr_code_str, cip, (uintptr_t)(chr_code_end - chr_code_str));
        if (!is_set(cip->chr_mask, chr_code)) {
          continue;
        }
      }
      loadbuf_iter = next_token_mult(loadbuf_iter, 4);
      uint32_t inner_loop_last = kBitsPerWordD2 - 1;
      uint32_t widx = 0;
      dosage_t* dosage_vals_iter = dosage_vals;
      while (1) {
        if (widx >= sample_ctl2_m1) {
          if (widx > sample_ctl2_m1) {
            break;
          }
          inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
        }
        uintptr_t genovec_word = 0;
        uint32_t dosage_present_hw = 0;
        for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
          loadbuf_iter = skip_initial_spaces(loadbuf_iter);
          const char cc = *loadbuf_iter;
          if (is_eoln_kns(cc)) {
            goto ox_gen_to_pgen_ret_MISSING_TOKENS;
          }
          char cc2 = loadbuf_iter[1];
          if ((cc2 == ' ') || (cc2 == '\t')) {
            cc2 = loadbuf_iter[3];
            if (((cc2 == ' ') || (cc2 == '\t')) && ((unsigned char)(loadbuf_iter[5]) <= 32)) {
              const uint32_t uii = ((uint32_t)((unsigned char)cc)) - 48;
              const uint32_t ujj = ((uint32_t)((unsigned char)loadbuf_iter[2])) - 48;
              const uint32_t ukk = ((uint32_t)((unsigned char)loadbuf_iter[4])) - 48;
              const uint32_t all_or = uii | ujj | ukk;
              if (all_or < 2) {
                if (!all_or) {
                  genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
                  loadbuf_iter = &(loadbuf_iter[5]);
                  continue;
                } else if (uii + ujj + ukk == 1) {
                  uintptr_t cur_geno = ukk * 2 + ujj;
                  genovec_word |= cur_geno << (2 * sample_idx_lowbits);
                  loadbuf_iter = &(loadbuf_iter[5]);
                  continue;
                }
              }
            }
          }
          double prob_0alt;
          char* first_dosage_str_end = scanadv_double(loadbuf_iter, &prob_0alt);
          if (!first_dosage_str_end) {
            // ignore next two tokens if first token in triplet is not numeric,
            // since we treat this as missing regardless
            loadbuf_iter = next_token_mult(loadbuf_iter, 2);
            if (!loadbuf_iter) {
              goto ox_gen_to_pgen_ret_MISSING_TOKENS;
            }
            genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
            loadbuf_iter = token_endnn(loadbuf_iter);
            continue;
          }
          if ((unsigned char)(*first_dosage_str_end) > ' ') {
            goto ox_gen_to_pgen_ret_INVALID_DOSAGE;
          }
          loadbuf_iter = skip_initial_spaces(first_dosage_str_end);
          if (is_eoln_kns(*loadbuf_iter)) {
            goto ox_gen_to_pgen_ret_MISSING_TOKENS;
          }
          double prob_1alt;
          loadbuf_iter = scanadv_double(loadbuf_iter, &prob_1alt);
          if ((!loadbuf_iter) || ((unsigned char)(*loadbuf_iter) > ' ')) {
            goto ox_gen_to_pgen_ret_INVALID_DOSAGE;
          }
          loadbuf_iter = skip_initial_spaces(loadbuf_iter);
          if (is_eoln_kns(*loadbuf_iter)) {
            goto ox_gen_to_pgen_ret_MISSING_TOKENS;
          }
          double prob_2alt;
          loadbuf_iter = scanadv_double(loadbuf_iter, &prob_2alt);
          if ((!loadbuf_iter) || ((unsigned char)(*loadbuf_iter) > ' ')) {
            goto ox_gen_to_pgen_ret_INVALID_DOSAGE;
          }
          // bugfix
          prob_0alt *= 32768;
          prob_1alt *= 32768;
          prob_2alt *= 32768;

          if ((prob_0alt < 0.0) || (prob_0alt >= 65535.4999999999) || (prob_1alt < 0.0) || (prob_1alt >= 65535.4999999999) || (prob_2alt < 0.0) || (prob_2alt >= 65535.4999999999)) {
            goto ox_gen_to_pgen_ret_INVALID_DOSAGE;
          }
          const uint32_t dosage_int0 = (int32_t)(prob_0alt + 0.5);
          const uint32_t dosage_int1 = (int32_t)(prob_1alt + 0.5);
          const uint32_t dosage_int2 = (int32_t)(prob_2alt + 0.5);
          bgen11_dosage_import_update(dosage_int_sum_thresh, import_dosage_certainty_int, hard_call_halfdist, dosage_erase_halfdist, sample_idx_lowbits, dosage_int0, dosage_int1, dosage_int2, &genovec_word, &dosage_present_hw, &dosage_vals_iter);
        }
        genovec[widx] = genovec_word;
        ((halfword_t*)dosage_present)[widx] = (halfword_t)dosage_present_hw;
        ++widx;
      }
      if (prov_ref_allele_second) {
        genovec_invert_unsafe(sample_ct, genovec);
        zero_trailing_quaters(sample_ct, genovec);
      }
      if (dosage_vals_iter != dosage_vals) {
        const uint32_t dosage_ct = (uintptr_t)(dosage_vals_iter - dosage_vals);
        if (prov_ref_allele_second) {
          biallelic_dosage16_invert(dosage_ct, dosage_vals);
        }
        if (spgw_append_biallelic_genovec_dosage16(genovec, dosage_present, dosage_vals, dosage_ct, &spgw)) {
          goto ox_gen_to_pgen_ret_WRITE_FAIL;
        }
      } else {
        if (spgw_append_biallelic_genovec(genovec, &spgw)) {
          goto ox_gen_to_pgen_ret_WRITE_FAIL;
        }
      }
      ++vidx;
      if (!(vidx % 1000)) {
        printf("\r--data/--gen: %uk variants converted.", vidx / 1000);
        fflush(stdout);
      }
    }
    spgw_finish(&spgw);
    putc_unlocked('\r', stdout);
    write_iter = strcpya(g_logbuf, "--data/--gen: ");
    const uint32_t outname_base_slen = (uintptr_t)(outname_end - outname);
    write_iter = memcpya(write_iter, outname, outname_base_slen + 5);
    write_iter = memcpyl3a(write_iter, " + ");
    write_iter = memcpya(write_iter, outname, outname_base_slen);
    write_iter = strcpya(write_iter, ".pvar");
    strcpy(write_iter, " written.\n");
    wordwrapb(0);
    logprintb();
  }
  while (0) {
  ox_gen_to_pgen_ret_LONG_LINE_N:
    putc_unlocked('\n', stdout);
    if (loadbuf_size == kMaxLongLine) {
      LOGERRPRINTF("Error: Line %" PRIuPTR " of .gen file is pathologically long.\n", line_idx);
      reterr = kPglRetMalformedInput;
      break;
    }
  ox_gen_to_pgen_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ox_gen_to_pgen_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  ox_gen_to_pgen_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  ox_gen_to_pgen_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  ox_gen_to_pgen_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  ox_gen_to_pgen_ret_INVALID_DOSAGE:
    putc_unlocked('\n', stdout);
    LOGERRPRINTF("Error: Line %" PRIuPTR " of .gen file has an invalid dosage value.\n", line_idx);
    reterr = kPglRetMalformedInput;
    break;
  ox_gen_to_pgen_ret_MISSING_TOKENS:
    putc_unlocked('\n', stdout);
    LOGERRPRINTF("Error: Line %" PRIuPTR " of .gen file has fewer tokens than expected.\n", line_idx);
    reterr = kPglRetMalformedInput;
    break;
  ox_gen_to_pgen_ret_MALFORMED_INPUT_WW:
    wordwrapb(0);
    logerrprintb();
  ox_gen_to_pgen_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  ox_gen_to_pgen_ret_INCONSISTENT_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  }
 ox_gen_to_pgen_ret_1:
  if (spgw_cleanup(&spgw) && (!reterr)) {
    reterr = kPglRetWriteFail;
  }
  gzclose_cond(gz_infile);
  fclose_cond(pvarfile);
  bigstack_reset(bigstack_mark);
  return reterr;
}

// a few multithread globals
static uint16_t** g_bgen_geno_bufs = nullptr; // per-thread


// per-variant (could make compressed_geno_starts per-thread)
static unsigned char** g_compressed_geno_starts[2] = {nullptr, nullptr};
static uintptr_t* g_write_genovecs[2] = {nullptr, nullptr};
static uint32_t* g_write_dosage_cts[2] = {nullptr, nullptr};
static uintptr_t* g_write_dosage_presents[2] = {nullptr, nullptr};
static dosage_t* g_write_dosage_val_bufs[2] = {nullptr, nullptr};

static uint32_t g_sample_ct = 0;
static uint32_t g_calc_thread_ct = 0;
static uint32_t g_cur_block_write_ct = 0;
static uint32_t g_hard_call_halfdist = 0;
static uint32_t g_dosage_erase_halfdist = 0;
static uint32_t g_import_dosage_certainty_int = 0;
static uint32_t g_compression_mode = 0;
static uint32_t g_dosage_is_present = 0;
static uint32_t g_prov_ref_allele_second = 0;
static pglerr_t g_error_ret = kPglRetSuccess;

// static uint32_t* g_error_vidxs = nullptr; // per-thread

THREAD_FUNC_DECL bgen11_dosage_scan_thread(void* arg) {
  // this bails as soon as a single non-hardcall is detected.  still
  // multithreaded due to low speed of uncompress() call, practical value of
  // handling the all-hardcall case efficiently, and reduced code complexity
  // (locally more complex, but globally cleaner due to overlap with
  // bgen11_geno_to_pgen_thread()).
  const uintptr_t tidx = (uintptr_t)arg;
  const uint32_t sample_ct = g_sample_ct;
  uint16_t* bgen_geno_buf = g_bgen_geno_bufs[tidx];
  const uint32_t calc_thread_ct = g_calc_thread_ct;
  // hard_call_halfdist irrelevant here
  const uint32_t dosage_erase_halfdist = g_dosage_erase_halfdist;
  const uint32_t import_dosage_certainty_int = g_import_dosage_certainty_int;
  const uint32_t dosage_int_sum_thresh = 3 * (import_dosage_certainty_int - 1);
  const uint32_t compression_mode = g_compression_mode;
  // uint32_t vidx_base = 0;
  uint32_t parity = 0;
  while (1) {
    const uint32_t is_last_block = g_is_last_thread_block;
    const uintptr_t cur_block_write_ct = g_cur_block_write_ct;
    const uint32_t vidx_end = ((tidx + 1) * cur_block_write_ct) / calc_thread_ct;
    uint32_t vidx = (tidx * cur_block_write_ct) / calc_thread_ct;
    unsigned char* compressed_geno_iter = g_compressed_geno_starts[parity][vidx];
    uint16_t* bgen_probs = bgen_geno_buf;
    for (; vidx < vidx_end; ++vidx) {
      if (compression_mode) {
        uint32_t compressed_block_byte_ct;
        memcpy(&compressed_block_byte_ct, compressed_geno_iter, 4);
        compressed_geno_iter = &(compressed_geno_iter[4]);
        uLongf zlib_ulongf = 6 * sample_ct;
        if (uncompress((Bytef*)bgen_probs, &zlib_ulongf, (Bytef*)compressed_geno_iter, compressed_block_byte_ct) != Z_OK) {
          break;
        }
        compressed_geno_iter = &(compressed_geno_iter[compressed_block_byte_ct]);
      } else {
        bgen_probs = (uint16_t*)compressed_geno_iter;
        compressed_geno_iter = &(compressed_geno_iter[6 * sample_ct]);
      }
      const uint16_t* bgen_probs_iter = bgen_probs;
      for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
        const uint32_t dosage_int0 = (uint32_t)(*bgen_probs_iter++);
        const uint32_t dosage_int1 = (uint32_t)(*bgen_probs_iter++);
        const uint32_t dosage_int2 = (uint32_t)(*bgen_probs_iter++);
        const uint32_t dosage_int_sum = dosage_int0 + dosage_int1 + dosage_int2;
        if ((dosage_int_sum > dosage_int_sum_thresh) || (dosage_int0 >= import_dosage_certainty_int) || (dosage_int1 >= import_dosage_certainty_int) || (dosage_int2 >= import_dosage_certainty_int)) {
          // ties realistically happen, use banker's rounding
          // 1/65536 -> 0/32768
          // 3/65536 -> 2/32768
          // 5/65536 -> 2/32768
          const dosage_prod_t write_dosage_int_numer = ((dosage_prod_t)kDosageMid) * dosage_int1 + ((dosage_prod_t)kDosageMax) * dosage_int2;
          uint32_t write_dosage_int;
          if (dosage_int_sum == kDosageMax) {
            // optimize common case
            write_dosage_int = ((write_dosage_int_numer + kDosageMid) / kDosageMax) - ((write_dosage_int_numer % (2 * ((dosage_prod_t)kDosageMax))) == kDosageMid);
          } else {
            write_dosage_int = (write_dosage_int_numer + (dosage_int_sum / 2)) / dosage_int_sum;
            write_dosage_int -= (2 * (write_dosage_int_numer - write_dosage_int * dosage_int_sum) == dosage_int_sum) * (write_dosage_int % 2);
          }
          const uint32_t halfdist = biallelic_dosage_halfdist(write_dosage_int);
          if (halfdist < dosage_erase_halfdist) {
            goto bgen11_dosage_scan_thread_dosage_found;
          }
        }
      }
    }
    if (vidx != vidx_end) {
      // g_error_vidxs[tidx] = vidx + vidx_base;
      g_error_ret = kPglRetMalformedInput;
    }
    while (0) {
    bgen11_dosage_scan_thread_dosage_found:
      g_dosage_is_present = 1;
    }
    // vidx_base += cur_block_write_ct;
    if (is_last_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
    parity = 1 - parity;
  }
}

static_assert(sizeof(dosage_t) == 2, "bgen11_geno_to_pgen_thread() needs to be updated.");
THREAD_FUNC_DECL bgen11_geno_to_pgen_thread(void* arg) {
  const uintptr_t tidx = (uintptr_t)arg;
  const uintptr_t sample_ct = g_sample_ct;
  uint16_t* bgen_geno_buf = g_bgen_geno_bufs[tidx];
  const uint32_t calc_thread_ct = g_calc_thread_ct;
  const uint32_t hard_call_halfdist = g_hard_call_halfdist;
  const uint32_t dosage_erase_halfdist = g_dosage_erase_halfdist;
  const uint32_t import_dosage_certainty_int = g_import_dosage_certainty_int;
  const uint32_t dosage_int_sum_thresh = 3 * (import_dosage_certainty_int - 1);
  const uint32_t compression_mode = g_compression_mode;
  const uint32_t prov_ref_allele_second = g_prov_ref_allele_second;
  const uintptr_t sample_ctaw2 = QUATERCT_TO_ALIGNED_WORDCT(sample_ct);
  const uint32_t sample_ctl2_m1 = (sample_ct - 1) / kBitsPerWordD2;
  const uintptr_t sample_ctaw = BITCT_TO_ALIGNED_WORDCT(sample_ct);
  // uint32_t vidx_base = 0;
  uint32_t parity = 0;
  while (1) {
    const uint32_t is_last_block = g_is_last_thread_block;
    const uintptr_t cur_block_write_ct = g_cur_block_write_ct;
    uint32_t vidx = (tidx * cur_block_write_ct) / calc_thread_ct;
    const uint32_t vidx_end = ((tidx + 1) * cur_block_write_ct) / calc_thread_ct;
    unsigned char* compressed_geno_iter = g_compressed_geno_starts[parity][vidx];
    uintptr_t* write_genovec_iter = &(g_write_genovecs[parity][vidx * sample_ctaw2]);
    uint32_t* write_dosage_ct_iter = &(g_write_dosage_cts[parity][vidx]);
    uintptr_t* write_dosage_present_iter = &(g_write_dosage_presents[parity][vidx * sample_ctaw]);
    dosage_t* write_dosage_vals_iter = &(g_write_dosage_val_bufs[parity][vidx * sample_ct]);
    uint16_t* bgen_probs = bgen_geno_buf;
    for (; vidx < vidx_end; ++vidx) {
      if (compression_mode) {
        uint32_t compressed_block_byte_ct;
        memcpy(&compressed_block_byte_ct, compressed_geno_iter, 4);
        compressed_geno_iter = &(compressed_geno_iter[4]);
        uLongf zlib_ulongf = 6 * sample_ct;
        if (uncompress((Bytef*)bgen_probs, &zlib_ulongf, (Bytef*)compressed_geno_iter, compressed_block_byte_ct) != Z_OK) {
          break;
        }
        compressed_geno_iter = &(compressed_geno_iter[compressed_block_byte_ct]);
      } else {
        bgen_probs = (uint16_t*)compressed_geno_iter;
        compressed_geno_iter = &(compressed_geno_iter[6 * sample_ct]);
      }
      const uint16_t* bgen_probs_iter = bgen_probs;
      dosage_t* cur_dosage_vals_iter = write_dosage_vals_iter;
      uint32_t inner_loop_last = kBitsPerWordD2 - 1;
      uint32_t widx = 0;
      while (1) {
        if (widx >= sample_ctl2_m1) {
          if (widx > sample_ctl2_m1) {
            break;
          }
          inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
        }
        uintptr_t genovec_word = 0;
        uint32_t dosage_present_hw = 0;
        for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
          const uint32_t dosage_int0 = (uint32_t)(*bgen_probs_iter++);
          const uint32_t dosage_int1 = (uint32_t)(*bgen_probs_iter++);
          const uint32_t dosage_int2 = (uint32_t)(*bgen_probs_iter++);
          const uint32_t dosage_int_sum = dosage_int0 + dosage_int1 + dosage_int2;
          if (dosage_int_sum <= dosage_int_sum_thresh) {
            if ((dosage_int0 < import_dosage_certainty_int) && (dosage_int1 < import_dosage_certainty_int) && (dosage_int2 < import_dosage_certainty_int)) {
              genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
              continue;
            }
          }
          const dosage_prod_t write_dosage_int_numer = ((dosage_prod_t)kDosageMid) * dosage_int1 + ((dosage_prod_t)kDosageMax) * dosage_int2;
          uint32_t write_dosage_int;
          if (dosage_int_sum == kDosageMax) {
            write_dosage_int = ((write_dosage_int_numer + kDosageMid) / kDosageMax) - ((write_dosage_int_numer % (2 * ((dosage_prod_t)kDosageMax))) == kDosageMid);
          } else {
            write_dosage_int = (write_dosage_int_numer + (dosage_int_sum / 2)) / dosage_int_sum;
            write_dosage_int -= (2 * (write_dosage_int_numer - write_dosage_int * dosage_int_sum) == dosage_int_sum) * (write_dosage_int % 2);
          }
          const uint32_t halfdist = biallelic_dosage_halfdist(write_dosage_int);
          if (halfdist < hard_call_halfdist) {
            genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
          } else {
            genovec_word |= ((write_dosage_int + (kDosage4th * k1LU)) / kDosageMid) << (2 * sample_idx_lowbits);
            if (halfdist >= dosage_erase_halfdist) {
              continue;
            }
          }
          dosage_present_hw |= 1U << sample_idx_lowbits;
          *cur_dosage_vals_iter++ = write_dosage_int;
        }
        write_genovec_iter[widx] = genovec_word;
        ((halfword_t*)write_dosage_present_iter)[widx] = (halfword_t)dosage_present_hw;
        ++widx;
      }
      const uint32_t dosage_ct = (uintptr_t)(cur_dosage_vals_iter - write_dosage_vals_iter);
      if (prov_ref_allele_second) {
        genovec_invert_unsafe(sample_ct, write_genovec_iter);
        zero_trailing_quaters(sample_ct, write_genovec_iter);
        if (dosage_ct) {
          biallelic_dosage16_invert(dosage_ct, write_dosage_vals_iter);
        }
      }
      *write_dosage_ct_iter++ = dosage_ct;
      write_genovec_iter = &(write_genovec_iter[sample_ctaw2]);
      write_dosage_present_iter = &(write_dosage_present_iter[sample_ctaw]);
      write_dosage_vals_iter = &(write_dosage_vals_iter[sample_ct]);
    }
    if (vidx != vidx_end) {
      // g_error_vidxs[tidx] = vidx + vidx_base;
      g_error_ret = kPglRetMalformedInput;
    }
    // vidx_base += cur_block_write_ct;
    if (is_last_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
    parity = 1 - parity;
  }
}

static unsigned char** g_thread_wkspaces = nullptr;
static uint32_t* g_thread_bidxs[2] = {nullptr, nullptr};
static uint16_t* g_bgen_allele_cts[2] = {nullptr, nullptr};
static uint32_t* g_uncompressed_genodata_byte_cts[2] = {nullptr, nullptr};

// for each bit precision level, how large must
//   max(numerators, 2^{bit_precision} - 1 - [sum of numerators])
// be to avoid throwing out the genotype?
static uint32_t* g_bgen_import_dosage_certainty_thresholds = nullptr;

// Reliably fast division by constants of the form 2^n - 1; see
//   http://ridiculousfish.com/blog/posts/labor-of-division-episode-iii.html
// The general case also requires a preshift parameter, but it's always zero
// for the odd .bgen denominators.
typedef struct bgen_magic_num_struct {
  uint64_t totq_magic;
  uint32_t totq_postshift;
  uint32_t totq_incr;
} bgen_magic_num_t;

static const bgen_magic_num_t kBgenMagicNums[25] = {
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
  {2147516417U, 47, 0}
  // todo: check whether something similar works for 17-32 bit cases
  /*
  ,{131073, 34, 1},
  {262145, 36, 1},
  {524289, 38, 1},
  {1048577, 40, 1},
  {2097153, 42, 1},
  {4194305, 44, 1},
  {8388609, 46, 1},
  {16777217, 48, 1},
  {33554433, 50, 1},
  {67108865, 52, 1},
  {134217729, 54, 1},
  {268435457, 56, 1},
  {536870913, 58, 1},
  {1073741825U, 60, 1},
  {2147483649U, 62, 1},
  {2147483649U, 63, 0}
  */
};

static_assert(sizeof(dosage_t) == 2, "bgen13_dosage_or_phase_scan_thread() needs to be updated.");
THREAD_FUNC_DECL bgen13_dosage_or_phase_scan_thread(void* arg) {
  // This bails as soon as a single phased or dosage call is detected.  We
  // provisionally assume e.g. phased calls are also present when dosages are,
  // and clean up the relevant header bytes when the assumption is untrue.
  // (Well, that's how it'll work after phased dosages are implemented,
  // anyway.)
  const uintptr_t tidx = (uintptr_t)arg;
  const uint32_t sample_ct = g_sample_ct;
  const uint32_t dosage_erase_halfdist = g_dosage_erase_halfdist;
  const uint32_t* bgen_import_dosage_certainty_thresholds = g_bgen_import_dosage_certainty_thresholds;
  const uint32_t compression_mode = g_compression_mode;
  const unsigned char* cur_uncompressed_geno = nullptr;
  if (compression_mode) {
    cur_uncompressed_geno = g_thread_wkspaces[tidx];
  }
  uint32_t parity = 0;
  while (1) {
    const uint32_t is_last_block = g_is_last_thread_block;

    // this is just used as an no-error flag
    const uintptr_t cur_block_write_ct = g_cur_block_write_ct;
    if (cur_block_write_ct) {
      const uint32_t bidx_end = g_thread_bidxs[parity][tidx + 1];
      uint32_t bidx = g_thread_bidxs[parity][tidx];
      unsigned char** compressed_geno_starts = g_compressed_geno_starts[parity];
      const uint16_t* bgen_allele_cts = g_bgen_allele_cts[parity];
      const uint32_t* uncompressed_genodata_byte_cts = g_uncompressed_genodata_byte_cts[parity];
      for (; bidx < bidx_end; ++bidx) {
        const unsigned char* compressed_geno_start = compressed_geno_starts[bidx];
        const unsigned char* compressed_geno_end = compressed_geno_starts[bidx + 1];
        uint32_t compressed_byte_ct = (uintptr_t)(compressed_geno_end - compressed_geno_start);
        uint32_t uncompressed_byte_ct;
        if (compression_mode) {
          uncompressed_byte_ct = uncompressed_genodata_byte_cts[bidx];
          if (compression_mode == 1) {
            uLongf zlib_ulongf = uncompressed_byte_ct;
            // const_cast
            if (uncompress((Bytef*)((uintptr_t)cur_uncompressed_geno), &zlib_ulongf, (const Bytef*)compressed_geno_start, compressed_byte_ct) != Z_OK) {
              // possible todo: report variant index
              goto bgen13_dosage_or_phase_scan_thread_malformed;
            }
          } else {
            // const_cast
            const uintptr_t extracted_byte_ct = ZSTD_decompress((void*)((uintptr_t)cur_uncompressed_geno), uncompressed_byte_ct, compressed_geno_start, compressed_byte_ct);
            if (extracted_byte_ct != uncompressed_byte_ct) {
              // possible todo: inspect error code
              goto bgen13_dosage_or_phase_scan_thread_malformed;
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
        // 1 byte: # of bits of probability precision (we just support 8 and 16
        //         for now, add others later)
#ifdef __arm__
  #error "Unaligned accesses in bgen13_dosage_or_phase_scan_thread()."
#endif
        if ((uncompressed_byte_ct < 10 + sample_ct) || (sample_ct != *((const uint32_t*)cur_uncompressed_geno))) {
          goto bgen13_dosage_or_phase_scan_thread_malformed;
        }
        const uint32_t cur_allele_ct = bgen_allele_cts[bidx];
        if (*((const uint16_t*)(&(cur_uncompressed_geno[4]))) != cur_allele_ct) {
          goto bgen13_dosage_or_phase_scan_thread_malformed;
        }
        const uint32_t min_ploidy = cur_uncompressed_geno[6];
        const uint32_t max_ploidy = cur_uncompressed_geno[7];
        if ((min_ploidy > max_ploidy) || (max_ploidy > 63)) {
          goto bgen13_dosage_or_phase_scan_thread_malformed;
        }
        if (max_ploidy > 2) {
          goto bgen13_dosage_or_phase_scan_thread_not_yet_supported;
        }
        const unsigned char* missing_and_ploidy_info = &(cur_uncompressed_geno[8]);
        const unsigned char* uncompressed_geno_iter = &(cur_uncompressed_geno[8 + sample_ct]);
        const uint32_t is_phased = *uncompressed_geno_iter++;
        if (is_phased > 1) {
          goto bgen13_dosage_or_phase_scan_thread_malformed;
        }
        const uint32_t bit_precision = *uncompressed_geno_iter++;
        if ((!bit_precision) || (bit_precision > 32)) {
          goto bgen13_dosage_or_phase_scan_thread_malformed;
        }
        if (bit_precision > 16) {
          goto bgen13_dosage_or_phase_scan_thread_not_yet_supported;
        }
        const uint64_t totq_magic = kBgenMagicNums[bit_precision].totq_magic;
        const uint32_t totq_postshift = kBgenMagicNums[bit_precision].totq_postshift;
        uint32_t totq_incr = kBgenMagicNums[bit_precision].totq_incr;
        const uint32_t bytes_per_prob = DIV_UP(bit_precision, CHAR_BIT);

        // also equal to denominator
        const uintptr_t numer_mask = (1U << bit_precision) - 1;

        // diploid (haploid is identical except b is always zero):
        //   round((32768a + 16384b)/(2^{bit precision} - 1))
        //   floor((32768a + 16384b)/(2^{bit_precision} - 1) + 0.5)
        // = floor((32768a + 16384b + 2^{bit_precision - 1})
        //     / (2^{bit_precision} - 1))
        // = (totq_magic * (32768a + 16384b + 2^{bits-1} + totq_incr))
        //     >> totq_postshift
        //
        // This works fine for bit_precision <= 16, anyway.  There are two
        // issues which come up with higher precision:
        // 1. The ridiculous_fish magic numbers assume a 32-bit dividend.  Our
        //    dividend is guaranteed to be divisible by 2^14, but it can be as
        //    large as
        //      (2^{bits} - 1) * 2^15 + 2^{bits-1}.
        //    I would not be surprised if a similar approach still works with
        //    bits > 16, but I'm pretty sure the magic-number-generating
        //    function would need to be different.
        // 2. Relatedly, the current sequence of operations multiplies
        //    totq_magic by (dividend + totq_incr) (where totq_incr is zero or
        //    one); this intermediate result must not overflow a uint64_t.
        //
        // Meanwhile, idempotence is not possible for --import-dosage-certainty
        // anyway, so we apply that check to the pre-conversion numerators.
        totq_incr += 1U << (bit_precision - 1);
        uint32_t numer_certainty_min = 0;
        if (bgen_import_dosage_certainty_thresholds) {
          numer_certainty_min = bgen_import_dosage_certainty_thresholds[bit_precision];
        }

        if (is_phased) {
          // todo
          goto bgen13_dosage_or_phase_scan_thread_not_yet_supported;
        } else {
          if (cur_allele_ct == 2) {
            if (min_ploidy == max_ploidy) {
              // faster handling of common cases (no need to keep checking if
              // we've read past the end)
              if (uncompressed_byte_ct != (1 + bytes_per_prob * (max_ploidy * k1LU)) * sample_ct + 10) {
                goto bgen13_dosage_or_phase_scan_thread_malformed;
              }
              if (max_ploidy == 2) {
                for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
                  const uint32_t missing_and_ploidy = missing_and_ploidy_info[sample_idx];
                  // treat anything else as missing
                  if (missing_and_ploidy == 2) {
                    const unsigned char* sample_probs_start = &(uncompressed_geno_iter[sample_idx * 2 * bytes_per_prob]);
                    // this can read slightly past the end of the buffer
                    const uintptr_t numer_aa = (*((const uint32_t*)sample_probs_start)) & numer_mask;
                    const uintptr_t numer_ab = (*((const uint32_t*)(&(sample_probs_start[bytes_per_prob])))) & numer_mask;
                    if ((numer_aa < numer_certainty_min) && (numer_ab < numer_certainty_min) && (numer_mask - numer_certainty_min < numer_aa + numer_ab)) {
                      // treat as missing
                      continue;
                    }
                    const uint32_t write_dosage_int = (totq_magic * (kDosageMax * ((uint64_t)numer_aa) + kDosageMid * ((uint64_t)numer_ab) + totq_incr)) >> totq_postshift;
                    const uint32_t halfdist = biallelic_dosage_halfdist(write_dosage_int);
                    if (halfdist < dosage_erase_halfdist) {
                      goto bgen13_dosage_or_phase_scan_thread_found;
                    }
                  }
                }
              } else if (max_ploidy == 1) {
                for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
                  const uint32_t missing_and_ploidy = missing_and_ploidy_info[sample_idx];
                  if (missing_and_ploidy == 1) {
                    const unsigned char* sample_probs_start = &(uncompressed_geno_iter[sample_idx * bytes_per_prob]);
                    const uintptr_t numer_a = (*((const uint32_t*)sample_probs_start)) & numer_mask;
                    if ((numer_a < numer_certainty_min) && (numer_mask - numer_certainty_min < numer_a)) {
                      continue;
                    }
                    const uint32_t write_dosage_int = (totq_magic * (kDosageMax * ((uint64_t)numer_a) + totq_incr)) >> totq_postshift;
                    const uint32_t halfdist = biallelic_dosage_halfdist(write_dosage_int);
                    if (halfdist < dosage_erase_halfdist) {
                      goto bgen13_dosage_or_phase_scan_thread_found;
                    }
                  }
                }
              }
              // don't need to do anything in all-ploidy-0 case
            } else {
              const unsigned char* uncompressed_geno_end = &(cur_uncompressed_geno[uncompressed_byte_ct]);
              for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
                if (uncompressed_geno_iter > uncompressed_geno_end) {
                  goto bgen13_dosage_or_phase_scan_thread_malformed;
                }
                uint32_t missing_and_ploidy = missing_and_ploidy_info[sample_idx];
                if (missing_and_ploidy == 2) {
                  const uintptr_t numer_aa = (*((const uint32_t*)uncompressed_geno_iter)) & numer_mask;
                  const uintptr_t numer_ab = (*((const uint32_t*)(&(uncompressed_geno_iter[bytes_per_prob])))) & numer_mask;
                  uncompressed_geno_iter = &(uncompressed_geno_iter[2 * bytes_per_prob]);
                  if ((numer_aa < numer_certainty_min) && (numer_ab < numer_certainty_min) && (numer_mask - numer_certainty_min < numer_aa + numer_ab)) {
                    // treat as missing
                    continue;
                  }
                  const uint32_t write_dosage_int = (totq_magic * (kDosageMax * ((uint64_t)numer_aa) + kDosageMid * ((uint64_t)numer_ab) + totq_incr)) >> totq_postshift;
                  const uint32_t halfdist = biallelic_dosage_halfdist(write_dosage_int);
                  if (halfdist < dosage_erase_halfdist) {
                    goto bgen13_dosage_or_phase_scan_thread_found;
                  }
                } else if (missing_and_ploidy == 1) {
                  const uintptr_t numer_a = (*((const uint32_t*)uncompressed_geno_iter)) & numer_mask;
                  uncompressed_geno_iter = &(uncompressed_geno_iter[bytes_per_prob]);
                  if ((numer_a < numer_certainty_min) && (numer_mask - numer_certainty_min < numer_a)) {
                    continue;
                  }
                  const uint32_t write_dosage_int = (totq_magic * (kDosageMax * ((uint64_t)numer_a) + totq_incr)) >> totq_postshift;
                  const uint32_t halfdist = biallelic_dosage_halfdist(write_dosage_int);
                  if (halfdist < dosage_erase_halfdist) {
                    goto bgen13_dosage_or_phase_scan_thread_found;
                  }
                } else {
                  // treat as missing
                  missing_and_ploidy &= 127;
                  if (missing_and_ploidy > 2) {
                    goto bgen13_dosage_or_phase_scan_thread_malformed;
                  }
                  uncompressed_geno_iter = &(uncompressed_geno_iter[missing_and_ploidy * bytes_per_prob]);
                }
              }
            }
          } else {
            // todo: unphased multiallelic variants
            // (shouldn't currently be possible to reach here, I/O thread skips
            // multiallelics for now)
            assert(0);
            goto bgen13_dosage_or_phase_scan_thread_not_yet_supported;
          }
        }
      }
    }
    while (0) {
    bgen13_dosage_or_phase_scan_thread_malformed:
      g_error_ret = kPglRetMalformedInput;
      break;
    bgen13_dosage_or_phase_scan_thread_not_yet_supported:
      g_error_ret = kPglRetNotYetSupported;
      break;
    bgen13_dosage_or_phase_scan_thread_found:
      g_dosage_is_present = 1;
      break;
    }
    if (is_last_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
    parity = 1 - parity;
  }
}

static uintptr_t* g_write_phasepresents[2] = {nullptr, nullptr};
static uintptr_t* g_write_phaseinfos[2] = {nullptr, nullptr};
static uintptr_t* g_write_dphase_presents[2] = {nullptr, nullptr};
static uint32_t* g_write_dphase_cts[2] = {nullptr, nullptr};

static_assert(sizeof(dosage_t) == 2, "bgen13_geno_to_pgen_thread() needs to be updated.");
THREAD_FUNC_DECL bgen13_geno_to_pgen_thread(void* arg) {
  const uintptr_t tidx = (uintptr_t)arg;
  const uintptr_t sample_ct = g_sample_ct;
  const uint32_t hard_call_halfdist = g_hard_call_halfdist;
  const uint32_t dosage_erase_halfdist = g_dosage_erase_halfdist;
  const uint32_t* bgen_import_dosage_certainty_thresholds = g_bgen_import_dosage_certainty_thresholds;
  const uint32_t compression_mode = g_compression_mode;
  const uint32_t prov_ref_allele_second = g_prov_ref_allele_second;
  const uintptr_t sample_ctaw2 = QUATERCT_TO_ALIGNED_WORDCT(sample_ct);
  const uint32_t sample_ctl2_m1 = (sample_ct - 1) / kBitsPerWordD2;
  const uintptr_t sample_ctaw = BITCT_TO_ALIGNED_WORDCT(sample_ct);
  const unsigned char* cur_uncompressed_geno = nullptr;
  if (compression_mode) {
    cur_uncompressed_geno = g_thread_wkspaces[tidx];
  }
  uint32_t parity = 0;
  while (1) {
    const uint32_t is_last_block = g_is_last_thread_block;

    // this is just used as an no-error flag
    const uintptr_t cur_block_write_ct = g_cur_block_write_ct;
    if (cur_block_write_ct) {
      const uint32_t bidx_end = g_thread_bidxs[parity][tidx + 1];
      uint32_t bidx = g_thread_bidxs[parity][tidx];
      unsigned char** compressed_geno_starts = g_compressed_geno_starts[parity];
      uintptr_t* write_genovec_iter = &(g_write_genovecs[parity][bidx * sample_ctaw2]);
      uint32_t* write_dosage_ct_iter = &(g_write_dosage_cts[parity][bidx]);
      uint32_t* write_dphase_ct_iter = &(g_write_dphase_cts[parity][bidx]);
      uintptr_t* write_dosage_present_iter = &(g_write_dosage_presents[parity][bidx * sample_ctaw]);
      uintptr_t* write_dphase_present_iter = &(g_write_dphase_presents[parity][bidx * sample_ctaw]);
      dosage_t* write_dosage_vals_iter = &(g_write_dosage_val_bufs[parity][bidx * sample_ct * 2]);
      const uint16_t* bgen_allele_ct_iter = &(g_bgen_allele_cts[parity][bidx]);
      const uint32_t* uncompressed_genodata_byte_ct_iter = &(g_uncompressed_genodata_byte_cts[parity][bidx]);
      for (; bidx < bidx_end; ++bidx) {
        const unsigned char* compressed_geno_start = compressed_geno_starts[bidx];
        const unsigned char* compressed_geno_end = compressed_geno_starts[bidx + 1];
        uint32_t compressed_byte_ct = (uintptr_t)(compressed_geno_end - compressed_geno_start);
        uint32_t uncompressed_byte_ct;
        if (compression_mode) {
          uncompressed_byte_ct = *uncompressed_genodata_byte_ct_iter++;
          if (compression_mode == 1) {
            uLongf zlib_ulongf = uncompressed_byte_ct;
            // const_cast
            if (uncompress((Bytef*)((uintptr_t)cur_uncompressed_geno), &zlib_ulongf, (const Bytef*)compressed_geno_start, compressed_byte_ct) != Z_OK) {
              // possible todo: report variant index
              goto bgen13_geno_to_pgen_thread_malformed;
            }
          } else {
            // const_cast
            const uintptr_t extracted_byte_ct = ZSTD_decompress((void*)((uintptr_t)cur_uncompressed_geno), uncompressed_byte_ct, compressed_geno_start, compressed_byte_ct);
            if (extracted_byte_ct != uncompressed_byte_ct) {
              // possible todo: inspect error code
              goto bgen13_geno_to_pgen_thread_malformed;
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
        // 1 byte: # of bits of probability precision (we just support 8 and 16
        //         for now, add others later)
        if ((uncompressed_byte_ct < 10 + sample_ct) || (sample_ct != *((const uint32_t*)cur_uncompressed_geno))) {
          goto bgen13_geno_to_pgen_thread_malformed;
        }
        const uint32_t cur_allele_ct = *bgen_allele_ct_iter++;
        if (*((const uint16_t*)(&(cur_uncompressed_geno[4]))) != cur_allele_ct) {
          goto bgen13_geno_to_pgen_thread_malformed;
        }
        const uint32_t min_ploidy = cur_uncompressed_geno[6];
        const uint32_t max_ploidy = cur_uncompressed_geno[7];
        if ((min_ploidy > max_ploidy) || (max_ploidy > 63)) {
          goto bgen13_geno_to_pgen_thread_malformed;
        }
        if (max_ploidy > 2) {
          goto bgen13_geno_to_pgen_thread_not_yet_supported;
        }
        const unsigned char* missing_and_ploidy_iter = &(cur_uncompressed_geno[8]);
        const unsigned char* uncompressed_geno_iter = &(cur_uncompressed_geno[8 + sample_ct]);
        const uint32_t is_phased = *uncompressed_geno_iter++;
        if (is_phased > 1) {
          goto bgen13_geno_to_pgen_thread_malformed;
        }
        const uint32_t bit_precision = *uncompressed_geno_iter++;
        if ((!bit_precision) || (bit_precision > 32)) {
          goto bgen13_geno_to_pgen_thread_malformed;
        }
        if (bit_precision > 16) {
          goto bgen13_geno_to_pgen_thread_not_yet_supported;
        }
        const uint64_t totq_magic = kBgenMagicNums[bit_precision].totq_magic;
        const uint32_t totq_postshift = kBgenMagicNums[bit_precision].totq_postshift;
        uint32_t totq_incr = kBgenMagicNums[bit_precision].totq_incr;
        const uint32_t bytes_per_prob = DIV_UP(bit_precision, CHAR_BIT);

        // also equal to denominator
        const uintptr_t numer_mask = (1U << bit_precision) - 1;

        totq_incr += 1U << (bit_precision - 1);
        uint32_t numer_certainty_min = 0;
        if (bgen_import_dosage_certainty_thresholds) {
          numer_certainty_min = bgen_import_dosage_certainty_thresholds[bit_precision];
        }

        dosage_t* cur_dosage_vals_iter = write_dosage_vals_iter;
        uint32_t inner_loop_last = kBitsPerWordD2 - 1;
        uint32_t widx = 0;
        if (is_phased) {
          // todo
          goto bgen13_geno_to_pgen_thread_not_yet_supported;
        } else {
          // fill_ulong_zero(sample_ctaw, write_dphase_present_iter);
          if (cur_allele_ct == 2) {
            if (min_ploidy == max_ploidy) {
              // faster handling of common cases (no need to keep checking if
              // we've read past the end)
              if (uncompressed_byte_ct != (bytes_per_prob * (max_ploidy * k1LU) + 1) * sample_ct + 10) {
                goto bgen13_geno_to_pgen_thread_malformed;
              }
              if (max_ploidy == 2) {
                while (1) {
                  if (widx >= sample_ctl2_m1) {
                    if (widx > sample_ctl2_m1) {
                      break;
                    }
                    inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
                  }
                  uintptr_t genovec_word = 0;
                  uint32_t dosage_present_hw = 0;
                  for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits, uncompressed_geno_iter = &(uncompressed_geno_iter[2 * bytes_per_prob])) {
                    const uint32_t missing_and_ploidy = *missing_and_ploidy_iter++;
                    if (missing_and_ploidy == 2) {
#ifdef __arm__
  #error "Unaligned accesses in bgen13_geno_to_pgen_thread()."
#endif
                      const uintptr_t numer_aa = (*((const uint32_t*)uncompressed_geno_iter)) & numer_mask;
                      const uintptr_t numer_ab = (*((const uint32_t*)(&(uncompressed_geno_iter[bytes_per_prob])))) & numer_mask;
                      if (numer_aa + numer_ab > numer_mask) {
                        goto bgen13_geno_to_pgen_thread_malformed;
                      }
                      if ((numer_aa < numer_certainty_min) && (numer_ab < numer_certainty_min) && (numer_mask - numer_certainty_min < numer_aa + numer_ab)) {
                        // missing due to --import-dosage-certainty
                        goto bgen13_geno_to_pgen_thread_diploid_missing;
                      }
                      const uint32_t write_dosage_int = (totq_magic * (kDosageMax * ((uint64_t)numer_aa) + kDosageMid * ((uint64_t)numer_ab) + totq_incr)) >> totq_postshift;
                      const uint32_t halfdist = biallelic_dosage_halfdist(write_dosage_int);
                      if (halfdist < hard_call_halfdist) {
                        genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
                      } else {
                        genovec_word |= ((write_dosage_int + (kDosage4th * k1LU)) / kDosageMid) << (2 * sample_idx_lowbits);
                        if (halfdist >= dosage_erase_halfdist) {
                          continue;
                        }
                      }
                      dosage_present_hw |= 1U << sample_idx_lowbits;
                      *cur_dosage_vals_iter++ = write_dosage_int;
                    } else {
                      // (could also validate that missing_and_ploidy == 130)
                    bgen13_geno_to_pgen_thread_diploid_missing:
                      genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
                    }
                  }
                  write_genovec_iter[widx] = genovec_word;
                  ((halfword_t*)write_dosage_present_iter)[widx] = (halfword_t)dosage_present_hw;
                  ++widx;
                }
              } else if (max_ploidy == 1) {
                while (1) {
                  if (widx >= sample_ctl2_m1) {
                    if (widx > sample_ctl2_m1) {
                      break;
                    }
                    inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
                  }
                  uintptr_t genovec_word = 0;
                  uint32_t dosage_present_hw = 0;
                  for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits, uncompressed_geno_iter = &(uncompressed_geno_iter[bytes_per_prob])) {
                    const uint32_t missing_and_ploidy = *missing_and_ploidy_iter++;
                    if (missing_and_ploidy == 1) {
                      const uintptr_t numer_a = (*((const uint32_t*)uncompressed_geno_iter)) & numer_mask;
                      if ((numer_a < numer_certainty_min) && (numer_mask - numer_certainty_min < numer_a)) {
                        goto bgen13_geno_to_pgen_thread_haploid_missing;
                      }
                      const uint32_t write_dosage_int = (totq_magic * (kDosageMax * ((uint64_t)numer_a) + totq_incr)) >> totq_postshift;
                      const uint32_t halfdist = biallelic_dosage_halfdist(write_dosage_int);
                      if (halfdist < hard_call_halfdist) {
                        genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
                      } else {
                        genovec_word |= ((write_dosage_int + (kDosage4th * k1LU)) / kDosageMid) << (2 * sample_idx_lowbits);
                        if (halfdist >= dosage_erase_halfdist) {
                          continue;
                        }
                      }
                      dosage_present_hw |= 1U << sample_idx_lowbits;
                      *cur_dosage_vals_iter++ = write_dosage_int;
                    } else {
                    bgen13_geno_to_pgen_thread_haploid_missing:
                      genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
                    }
                  }
                  write_genovec_iter[widx] = genovec_word;
                  ((halfword_t*)write_dosage_present_iter)[widx] = (halfword_t)dosage_present_hw;
                  ++widx;
                }
              }
              // don't need to do anything in all-ploidy-0 case
            } else {
              const unsigned char* uncompressed_geno_end = &(cur_uncompressed_geno[uncompressed_byte_ct]);
              while (1) {
                if (widx >= sample_ctl2_m1) {
                  if (widx > sample_ctl2_m1) {
                    break;
                  }
                  inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
                }
                uintptr_t genovec_word = 0;
                uint32_t dosage_present_hw = 0;
                for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
                  if (uncompressed_geno_iter > uncompressed_geno_end) {
                    goto bgen13_geno_to_pgen_thread_malformed;
                  }
                  uint32_t missing_and_ploidy = *missing_and_ploidy_iter++;
                  uint32_t write_dosage_int;
                  if (missing_and_ploidy == 2) {
                    const uintptr_t numer_aa = (*((const uint32_t*)uncompressed_geno_iter)) & numer_mask;
                    const uintptr_t numer_ab = (*((const uint32_t*)(&(uncompressed_geno_iter[bytes_per_prob])))) & numer_mask;
                    uncompressed_geno_iter = &(uncompressed_geno_iter[2 * bytes_per_prob]);
                    if (numer_aa + numer_ab > numer_mask) {
                      goto bgen13_geno_to_pgen_thread_malformed;
                    }
                    if ((numer_aa < numer_certainty_min) && (numer_ab < numer_certainty_min) && (numer_mask - numer_certainty_min < numer_aa + numer_ab)) {
                      // missing due to --import-dosage-certainty
                      goto bgen13_geno_to_pgen_thread_generic_missing;
                    }
                    write_dosage_int = (totq_magic * (kDosageMax * ((uint64_t)numer_aa) + kDosageMid * ((uint64_t)numer_ab) + totq_incr)) >> totq_postshift;
                  } else if (missing_and_ploidy == 1) {
                    const uintptr_t numer_a = (*((const uint32_t*)uncompressed_geno_iter)) & numer_mask;
                    uncompressed_geno_iter = &(uncompressed_geno_iter[bytes_per_prob]);
                    if ((numer_a < numer_certainty_min) && (numer_mask - numer_certainty_min < numer_a)) {
                      goto bgen13_geno_to_pgen_thread_generic_missing;
                    }
                    write_dosage_int = (totq_magic * (kDosageMax * ((uint64_t)numer_a) + totq_incr)) >> totq_postshift;
                  } else {
                    missing_and_ploidy &= 127;
                    if (missing_and_ploidy > 2) {
                      goto bgen13_geno_to_pgen_thread_malformed;
                    }
                    uncompressed_geno_iter = &(uncompressed_geno_iter[missing_and_ploidy * bytes_per_prob]);
                  bgen13_geno_to_pgen_thread_generic_missing:
                    genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
                    continue;
                  }
                  const uint32_t halfdist = biallelic_dosage_halfdist(write_dosage_int);
                  if (halfdist < hard_call_halfdist) {
                    genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
                  } else {
                    genovec_word |= ((write_dosage_int + (kDosage4th * k1LU)) / kDosageMid) << (2 * sample_idx_lowbits);
                    if (halfdist >= dosage_erase_halfdist) {
                      continue;
                    }
                  }
                  dosage_present_hw |= 1U << sample_idx_lowbits;
                  *cur_dosage_vals_iter++ = write_dosage_int;
                }
                write_genovec_iter[widx] = genovec_word;
                ((halfword_t*)write_dosage_present_iter)[widx] = (halfword_t)dosage_present_hw;
                ++widx;
              }
            }
            const uint32_t dosage_ct = (uintptr_t)(cur_dosage_vals_iter - write_dosage_vals_iter);
            // note that this is inverted from bgen-1.1
            if (!prov_ref_allele_second) {
              genovec_invert_unsafe(sample_ct, write_genovec_iter);
              zero_trailing_quaters(sample_ct, write_genovec_iter);
              if (dosage_ct) {
                biallelic_dosage16_invert(dosage_ct, write_dosage_vals_iter);
              }
            }
            *write_dosage_ct_iter++ = dosage_ct;
            *write_dphase_ct_iter++ = 0;
            write_genovec_iter = &(write_genovec_iter[sample_ctaw2]);
            write_dosage_present_iter = &(write_dosage_present_iter[sample_ctaw]);
            write_dphase_present_iter = &(write_dphase_present_iter[sample_ctaw]);
            write_dosage_vals_iter = &(write_dosage_vals_iter[sample_ct * 2]);
          } else {
            // todo: unphased multiallelic variants
            assert(0);
            goto bgen13_geno_to_pgen_thread_not_yet_supported;
          }
        }
      }
    }
    while (0) {
    bgen13_geno_to_pgen_thread_malformed:
      g_error_ret = kPglRetMalformedInput;
      break;
    bgen13_geno_to_pgen_thread_not_yet_supported:
      g_error_ret = kPglRetNotYetSupported;
      break;
    }
    if (is_last_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
    parity = 1 - parity;
  }
}

static_assert(sizeof(dosage_t) == 2, "ox_bgen_to_pgen() needs to be updated.");
pglerr_t ox_bgen_to_pgen(const char* bgenname, const char* samplename, const char* const_fid, const char* ox_missing_code, misc_flags_t misc_flags, oxford_import_t oxford_import_flags, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, double import_dosage_certainty, char id_delim, char idspace_to, uint32_t max_thread_ct, char* outname, char* outname_end, chr_info_t* cip) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  FILE* bgenfile = nullptr;

  // only if no sample file specified, and .bgen has sample IDs.  (possible
  // todo: consistency check when both sources of sample IDs are present?)
  FILE* psamfile = nullptr;

  FILE* pvarfile = nullptr;
  threads_state_t ts;
  init_threads3z(&ts);
  st_pgen_writer_t spgw;
  pglerr_t reterr = kPglRetSuccess;
  spgw_preinit(&spgw);
  {
    // Pass 1: Determine whether there's at least one non-hardcall needs to be
    //         saved, and if a chromosome filter was specified, count the
    //         number of variants which pass the filter.
    //         For bgen-1.2/1.3, the .pvar is also written in this pass.
    //         For bgen-1.1, we can usually early-bail when no chromosome
    //         filter is involved, so .pvar writing is postponed till the
    //         second pass.
    // Pass 2: Write .pgen file.

    if (fopen_checked(bgenname, FOPEN_RB, &bgenfile)) {
      goto ox_bgen_to_pgen_ret_OPEN_FAIL;
    }
    uint32_t initial_uints[5];
    if (!fread_unlocked(initial_uints, 20, 1, bgenfile)) {
      // this could be malformed input as well; could distinguish later?
      goto ox_bgen_to_pgen_ret_READ_FAIL;
    }
    if (initial_uints[1] > initial_uints[0]) {
      logerrprint("Error: Invalid .bgen header.\n");
      goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
    }
    const uint32_t raw_variant_ct = initial_uints[2];
    if (!raw_variant_ct) {
      // permit this in --allow-no-vars case?
      logerrprint("Error: Empty .bgen file.\n");
      goto ox_bgen_to_pgen_ret_INCONSISTENT_INPUT;
    }
    const uint32_t sample_ct = initial_uints[3];
    if (initial_uints[4] && (initial_uints[4] != 0x6e656762)) {
      logerrprint("Error: Invalid .bgen magic number.\n");
      goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
    }

    if (fseeko(bgenfile, initial_uints[1], SEEK_SET)) {
      goto ox_bgen_to_pgen_ret_READ_FAIL;
    }
    uint32_t header_flags;
    if (!fread_unlocked(&header_flags, 4, 1, bgenfile)) {
      goto ox_bgen_to_pgen_ret_READ_FAIL;
    }
    const uint32_t compression_mode = header_flags & 3;
    const uint32_t layout = (header_flags >> 2) & 15;
    if (!layout) {
      logerrprint("Error: BGEN v1.0 files are not supported by " PROG_NAME_STR ".\n");
      goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
    }
    if ((compression_mode == 3) || (layout > 2)) {
      logerrprint("Error: Unrecognized BGEN version.  Use gen-convert or a similar tool to\ndowncode to BGEN v1.3 if you want to process this data with " PROG_NAME_STR ".\n");
      goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
    }
    if ((compression_mode == 2) && (layout == 1)) {
      logerrprint("Error: Invalid .bgen header.\n");
      goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
    }
    LOGPRINTF("--bgen: %u variant%s detected, format v1.%c.\n", raw_variant_ct, (raw_variant_ct == 1)? "" : "s", (layout == 1)? '1' : ((compression_mode == 2)? '3' : '2'));
    if (samplename[0]) {
      uint32_t sfile_sample_ct;
      reterr = ox_sample_to_psam(samplename, ox_missing_code, misc_flags, outname, outname_end, &sfile_sample_ct);
      if (reterr) {
        goto ox_bgen_to_pgen_ret_1;
      }
      if (sfile_sample_ct != sample_ct) {
        LOGERRPRINTF("Error: .sample file has %u sample%s, while .bgen file has %u.\n", sfile_sample_ct, (sfile_sample_ct == 1)? "" : "s", sample_ct);
        goto ox_bgen_to_pgen_ret_INCONSISTENT_INPUT;
      }
      if (header_flags >> 31) {
        uint32_t sample_id_block_byte_ct;
        uint32_t sample_id_block_entry_ct;
        if ((!fread_unlocked(&sample_id_block_byte_ct, 4, 1, bgenfile)) ||
            (!fread_unlocked(&sample_id_block_entry_ct, 4, 1, bgenfile))) {
          goto ox_bgen_to_pgen_ret_READ_FAIL;
        }
        if ((((uint64_t)sample_id_block_byte_ct) + initial_uints[1] > initial_uints[0]) ||
            (sample_id_block_entry_ct != sample_ct)) {
          logerrprint("Error: Invalid .bgen header.\n");
          goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
        }
      }
    } else {
      if (!(header_flags >> 31)) {
        logerrprint("Error: .bgen file does not contain sample IDs, and no .sample file was\nspecified.\n");
        goto ox_bgen_to_pgen_ret_INCONSISTENT_INPUT;
      }
      // possible todo: optionally error out if sample IDs aren't consistent
      // between .bgen and .sample

      // see vcf_sample_line()
      // probable todo: wrap much of this in its own function
      uint32_t double_id = (misc_flags / kfMiscDoubleId) & 1;
      uintptr_t const_fid_len = 0;
      if (const_fid) {
        const_fid_len = strlen(const_fid);
      } else if ((!double_id) && (!id_delim)) {
        // default: --double-id + --id-delim
        double_id = 1;
        id_delim = '_';
      }
      const uint32_t double_or_const_fid = double_id || const_fid;
      uint32_t sample_id_block_byte_ct;
      uint32_t sample_id_block_entry_ct;
      if ((!fread_unlocked(&sample_id_block_byte_ct, 4, 1, bgenfile)) ||
          (!fread_unlocked(&sample_id_block_entry_ct, 4, 1, bgenfile))) {
        goto ox_bgen_to_pgen_ret_READ_FAIL;
      }
      if ((sample_id_block_byte_ct < 8) ||
          (((uint64_t)sample_id_block_byte_ct) + initial_uints[1] > initial_uints[0]) ||
          (sample_id_block_entry_ct != sample_ct)) {
        logerrprint("Error: Invalid .bgen header.\n");
        goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
      }
      sample_id_block_byte_ct -= 8;
      unsigned char* sample_id_block_main = bigstack_alloc(sample_id_block_byte_ct);
      if (!sample_id_block_main) {
        goto ox_bgen_to_pgen_ret_NOMEM;
      }
      unsigned char* sample_id_block_end = &(sample_id_block_main[sample_id_block_byte_ct]);
      if (fread_checked(sample_id_block_main, sample_id_block_byte_ct, bgenfile)) {
        goto ox_bgen_to_pgen_ret_READ_FAIL;
      }

      // always check if any tab/eoln characters are present, and error out if
      // so
      // if id_delim != ' ', also check if spaces are present; if so, replace
      // with --idspace-to character or error out
      unsigned char* sample_id_block_iter = sample_id_block_main;
      for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
        // don't use uint16_t here since we add 2
        const uint32_t input_id_slen = *((uint16_t*)sample_id_block_iter);

        // need to check this to avoid read-past-the-end indeterminate
        // behavior
        if ((uintptr_t)(sample_id_block_end - sample_id_block_iter) < input_id_slen + 2) {
          logerrprint("Error: Invalid .bgen header.\n");
          goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
        }
        unsigned char* sample_id_iter = &(sample_id_block_iter[2]);
        unsigned char* sample_id_end = &(sample_id_iter[input_id_slen]);
        uint32_t char_code_min = 32 + (id_delim != ' ');
        for (; sample_id_iter != sample_id_end; ++sample_id_iter) {
          const uint32_t char_code = *sample_id_iter;
          if (char_code < char_code_min) {
            if (char_code < 32) {
              logerrprint("Error: .bgen sample ID contains tabs, newlines, and/or nonprinting characters.\n");
              goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
            }
            if (!idspace_to) {
              logerrprint("Error: .bgen sample ID contains space(s).  Use --idspace-to to convert them to\nanother character, or \"--id-delim ' '\" to interpret the spaces as FID/IID\ndelimiters.\n");
              goto ox_bgen_to_pgen_ret_INCONSISTENT_INPUT;
            }
            *sample_id_iter = idspace_to;
          }
        }
        sample_id_block_iter = sample_id_end;
      }
      strcpy(outname_end, ".psam");
      if (fopen_checked(outname, FOPEN_WB, &psamfile)) {
        goto ox_bgen_to_pgen_ret_OPEN_FAIL;
      }
      char* textbuf = g_textbuf;
      char* write_iter = strcpya(textbuf, "#FID\tIID");
      uint32_t sid_present = 0;
      if (id_delim) {
        // check if three-part IDs are present
        sample_id_block_iter = sample_id_block_main;
        for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
          const uint32_t input_id_slen = *((uint16_t*)sample_id_block_iter);
          if ((uintptr_t)(sample_id_block_end - sample_id_block_iter) < input_id_slen + 2) {
            logerrprint("Error: Invalid .bgen header.\n");
            goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
          }
          unsigned char* sample_id_start = &(sample_id_block_iter[2]);
          unsigned char* sample_id_end = &(sample_id_start[input_id_slen]);
          unsigned char* first_delim = (unsigned char*)memchr(sample_id_start, (unsigned char)id_delim, (uintptr_t)(sample_id_end - sample_id_start));
          if (first_delim) {
            unsigned char* iid_start = &(first_delim[1]);
            if (memchr(iid_start, (unsigned char)id_delim, (uintptr_t)(sample_id_end - iid_start)) != nullptr) {
              sid_present = 1;
              write_iter = strcpya(write_iter, "\tSID");
              break;
            }
          }
          sample_id_block_iter = sample_id_end;
        }
      }
      write_iter = strcpya(write_iter, "\tSEX");
      append_binary_eoln(&write_iter);
      char* textbuf_flush = &(textbuf[kMaxMediumLine]);
      sample_id_block_iter = sample_id_block_main;
      for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
        const uint32_t input_id_slen = *((uint16_t*)sample_id_block_iter);
        if ((uintptr_t)(sample_id_block_end - sample_id_block_iter) < input_id_slen + 2) {
          logerrprint("Error: Invalid .bgen header.\n");
          goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
        }
        unsigned char* sample_id_start = &(sample_id_block_iter[2]);
        if (input_id_slen <= 1) {
          if (!input_id_slen) {
            logerrprint("Error: Empty sample ID in .bgen file.\n");
            goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
          }
          if (*sample_id_start == '0') {
            logerrprint("Error: Sample ID cannot be '0'.\n");
            goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
          }
        }
        unsigned char* sample_id_end = &(sample_id_start[input_id_slen]);
        if (id_delim) {
          if (*sample_id_start == id_delim) {
            sprintf(g_logbuf, "Error: '%c' at beginning of sample ID.\n", id_delim);
            goto ox_bgen_to_pgen_ret_INCONSISTENT_INPUT_2;
          }
          unsigned char* first_delim = (unsigned char*)memchr(sample_id_start, (unsigned char)id_delim, input_id_slen);
          if (!first_delim) {
            if (double_or_const_fid) {
              goto ox_bgen_to_pgen_one_sample_id;
            }
            sprintf(g_logbuf, "Error: No '%c' in sample ID.\n", id_delim);
            goto ox_bgen_to_pgen_ret_INCONSISTENT_INPUT_2;
          }
          unsigned char* iid_start = &(first_delim[1]);
          unsigned char* iid_end = (unsigned char*)memchr(iid_start, (unsigned char)id_delim, (uintptr_t)(sample_id_end - iid_start));
          const unsigned char* sid_start = (const unsigned char*)(&(g_one_char_strs[96]));
          uint32_t sid_slen = 1;
          if (iid_end) {
            if (iid_start == iid_end) {
              sprintf(g_logbuf, "Error: Consecutive instances of '%c' in sample ID.\n", id_delim);
              goto ox_bgen_to_pgen_ret_INCONSISTENT_INPUT_DELIM;
            }
            sid_start = &(iid_end[1]);
            sid_slen = (uintptr_t)(sample_id_end - sid_start);
            if (memchr(sid_start, (unsigned char)id_delim, sid_slen)) {
              sprintf(g_logbuf, "Error: More than two instances of '%c' in sample ID.\n", id_delim);
              goto ox_bgen_to_pgen_ret_INCONSISTENT_INPUT_DELIM;
            }
          } else {
            iid_end = sample_id_end;
          }
          const uint32_t fid_slen = (uintptr_t)(first_delim - sample_id_start);
          if (fid_slen > kMaxIdSlen) {
            goto ox_bgen_to_pgen_ret_MALFORMED_INPUT_LONG_ID;
          }
          write_iter = memcpyax(write_iter, sample_id_start, fid_slen, '\t');
          const uint32_t iid_slen = (uintptr_t)(iid_end - iid_start);
          if ((*iid_start == '0') && (iid_slen == 1)) {
            logerrprint("Error: Sample ID induces an invalid IID of '0'.\n");
            goto ox_bgen_to_pgen_ret_INCONSISTENT_INPUT;
          }
          if (iid_slen > kMaxIdSlen) {
            goto ox_bgen_to_pgen_ret_MALFORMED_INPUT_LONG_ID;
          }
          write_iter = memcpya(write_iter, iid_start, iid_slen);
          if (sid_present) {
            *write_iter++ = '\t';
            write_iter = memcpya(write_iter, sid_start, sid_slen);
          }
        } else {
        ox_bgen_to_pgen_one_sample_id:
          if (input_id_slen > kMaxIdSlen) {
            goto ox_bgen_to_pgen_ret_MALFORMED_INPUT_LONG_ID;
          }
          if (double_id) {
            write_iter = memcpya(write_iter, sample_id_start, input_id_slen);
          } else {
            write_iter = memcpya(write_iter, const_fid, const_fid_len);
          }
          *write_iter++ = '\t';
          write_iter = memcpya(write_iter, sample_id_start, input_id_slen);
          if (sid_present) {
            write_iter = strcpya(write_iter, "\t0");
          }
        }
        // SEX
        write_iter = memcpyl3a(write_iter, "\tNA");
        append_binary_eoln(&write_iter);
        if (write_iter >= textbuf_flush) {
          if (fwrite_checked(textbuf, (uintptr_t)(write_iter - textbuf), psamfile)) {
            goto ox_bgen_to_pgen_ret_WRITE_FAIL;
          }
          write_iter = textbuf;
        }
        sample_id_block_iter = sample_id_end;
      }
      if (sample_id_block_iter != &(sample_id_block_main[sample_id_block_byte_ct])) {
        logerrprint("Error: Invalid .bgen header.\n");
        goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
      }
      if (write_iter != textbuf) {
        if (fwrite_checked(textbuf, (uintptr_t)(write_iter - textbuf), psamfile)) {
          goto ox_bgen_to_pgen_ret_WRITE_FAIL;
        }
      }
      bigstack_reset(sample_id_block_main);
      if (fclose_null(&psamfile)) {
        goto ox_bgen_to_pgen_ret_WRITE_FAIL;
      }
      LOGPRINTFWW("--bgen: %u sample ID%s written to %s .\n", sample_ct, (sample_ct == 1)? "" : "s", outname);
    }
    if (fseeko(bgenfile, initial_uints[0] + 4, SEEK_SET)) {
      goto ox_bgen_to_pgen_ret_READ_FAIL;
    }
    const uint32_t allow_extra_chrs = (misc_flags / kfMiscAllowExtraChrs) & 1;
    finalize_chrset(misc_flags, cip);
    const uint32_t autosome_ct_p1 = cip->autosome_ct + 1;
    uint32_t chr_filter_present = (popcount_bit_idx(cip->chr_mask, 0, autosome_ct_p1) != autosome_ct_p1) || (allow_extra_chrs && (cip->is_include_stack || cip->incl_excl_name_stack));
    if (!chr_filter_present) {
      for (uint32_t xymt_idx = 0; xymt_idx < kChrOffsetCt; ++xymt_idx) {
        if (cip->xymt_codes[xymt_idx] >= 0) {
          if (!is_set(cip->chr_mask, autosome_ct_p1 + xymt_idx)) {
            chr_filter_present = 1;
            break;
          }
        }
      }
    }

    char* writebuf = (char*)bigstack_alloc_raw(2 * kMaxMediumLine + kCacheline);
    char* writebuf_flush = &(writebuf[kMaxMediumLine]);
    strcpy(outname_end, ".pvar");
    if (fopen_checked(outname, FOPEN_WB, &pvarfile)) {
      goto ox_bgen_to_pgen_ret_OPEN_FAIL;
    }
    char* write_iter = writebuf;
    if (cip->chrset_source) {
      append_chrset_line(cip, &write_iter);
    }
    write_iter = strcpya(write_iter, "#CHROM\tPOS\tID\tREF\tALT" EOLN_STR);

    const uint32_t snpid_chr = (oxford_import_flags & kfOxfordImportBgenSnpIdChr);

    // true for both provisional-reference and real-reference second
    const uint32_t prov_ref_allele_second = !(oxford_import_flags & kfOxfordImportRefFirst);

    if (hard_call_thresh == UINT32_MAX) {
      hard_call_thresh = kDosageMid / 10;
    }
    const uint32_t sample_ctaw2 = QUATERCT_TO_ALIGNED_WORDCT(sample_ct);
    const uint32_t sample_ctaw = BITCT_TO_ALIGNED_WORDCT(sample_ct);
    uint32_t dosage_is_present = 0;
    g_sample_ct = sample_ct;
    g_hard_call_halfdist = kDosage4th - hard_call_thresh;
    g_dosage_erase_halfdist = kDosage4th - dosage_erase_thresh;
    g_compression_mode = compression_mode;
    g_prov_ref_allele_second = prov_ref_allele_second;
    g_error_ret = kPglRetSuccess;
    g_dosage_is_present = 0;
    if (layout == 1) {
      // v1.1
      uintptr_t loadbuf_size = round_down_pow2(bigstack_left() / 4, kCacheline);
#ifdef __LP64__
      if (loadbuf_size > kMaxLongLine) {
        loadbuf_size = kMaxLongLine;
      }
#endif
      // must have enough space for chromosome and variant IDs
      if (loadbuf_size < 2 * 65536) {
        goto ox_bgen_to_pgen_ret_NOMEM;
      }
      unsigned char* loadbuf = (unsigned char*)bigstack_alloc_raw(loadbuf_size);
      g_import_dosage_certainty_int = 1 + (int32_t)(import_dosage_certainty * 32768);
      uintptr_t bgen_geno_max_byte_ct = 6LU * sample_ct;
      if (compression_mode) {
        bgen_geno_max_byte_ct = compressBound(bgen_geno_max_byte_ct);
      }
      if (bgen_geno_max_byte_ct > UINT32_MAX) {
        logerrprint("Error: Too many samples for .bgen format.\n");
        goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
      }
      bgen_geno_max_byte_ct += compression_mode * 4;
      // thread-count-independent:
      //   (everything after "2 *" rounded up to cacheline)
      //   compressed_geno_bufs: 2 * bgen_geno_max_byte_ct * main_block_size
      //   g_compressed_geno_starts: 2 * sizeof(intptr_t) * main_block_size
      //   g_write_genovecs: 2 * sample_ctaw2 * sizeof(intptr_t) *
      //                     main_block_size
      //   g_write_dosage_cts: 2 * sizeof(int32_t) * main_block_size
      //   g_write_dosage_presents: 2 * sample_ctaw * sizeof(intptr_t) *
      //                            main_block_size
      //   g_write_dosage_val_bufs (main bottleneck): 2 * sample_ct *
      //                                              sizeof(dosage_t)
      // additional requirement per thread:
      //   g_bgen_geno_bufs: sample_ct * 3 * sizeof(int16_t)

      uint32_t calc_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
      if ((!compression_mode) && (calc_thread_ct > 2)) {
        // computation doesn't seem to saturate when decompression is involved
        calc_thread_ct = 2;
      }
      if (calc_thread_ct > raw_variant_ct) {
        calc_thread_ct = raw_variant_ct;
      }
      if (bigstack_alloc_thread(calc_thread_ct, &ts.threads) ||
          bigstack_alloc_usip(calc_thread_ct, &g_bgen_geno_bufs)) {
        goto ox_bgen_to_pgen_ret_NOMEM;
      }
      const uint32_t sample_ct_x3 = sample_ct * 3;
      for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
        if (bigstack_alloc_usi(sample_ct_x3, &(g_bgen_geno_bufs[tidx]))) {
          goto ox_bgen_to_pgen_ret_NOMEM;
        }
      }
      uintptr_t cachelines_avail_m12 = bigstack_left() / kCacheline;
      // reserve 1/8 of remaining memory for writer
      cachelines_avail_m12 -= cachelines_avail_m12 / 8;
      if (cachelines_avail_m12 < 12) {
        goto ox_bgen_to_pgen_ret_NOMEM;
      }
      // we're making 12 allocations; be pessimistic re: rounding
      cachelines_avail_m12 -= 12;
      const uintptr_t bytes_req_per_in_block_variant = 2 * (bgen_geno_max_byte_ct + sizeof(intptr_t) + sample_ctaw2 * sizeof(intptr_t) + sizeof(int32_t) + sample_ctaw * sizeof(intptr_t) + sample_ct * sizeof(dosage_t));
      uintptr_t main_block_size = (cachelines_avail_m12 * kCacheline) / bytes_req_per_in_block_variant;
      if (main_block_size > 65536) {
        main_block_size = 65536;
      } else if (main_block_size < 8) {
        // this threshold is arbitrary
        goto ox_bgen_to_pgen_ret_NOMEM;
      }
      if (calc_thread_ct > main_block_size / 8) {
        calc_thread_ct = main_block_size / 8;
      }
      ts.calc_thread_ct = calc_thread_ct;
      g_calc_thread_ct = calc_thread_ct;
      unsigned char* compressed_geno_bufs[2];
      if (bigstack_alloc_uc(bgen_geno_max_byte_ct * main_block_size, &(compressed_geno_bufs[0])) ||
          bigstack_alloc_uc(bgen_geno_max_byte_ct * main_block_size, &(compressed_geno_bufs[1])) ||
          bigstack_alloc_ucp(main_block_size, &(g_compressed_geno_starts[0])) ||
          bigstack_alloc_ucp(main_block_size, &(g_compressed_geno_starts[1])) ||
          bigstack_alloc_ul(sample_ctaw2 * main_block_size, &(g_write_genovecs[0])) ||
          bigstack_alloc_ul(sample_ctaw2 * main_block_size, &(g_write_genovecs[1])) ||
          bigstack_alloc_ui(main_block_size, &(g_write_dosage_cts[0])) ||
          bigstack_alloc_ui(main_block_size, &(g_write_dosage_cts[1])) ||
          bigstack_alloc_ul(sample_ctaw * main_block_size, &(g_write_dosage_presents[0])) ||
          bigstack_alloc_ul(sample_ctaw * main_block_size, &(g_write_dosage_presents[1])) ||
          bigstack_alloc_dosage(sample_ct * main_block_size, &(g_write_dosage_val_bufs[0])) ||
          bigstack_alloc_dosage(sample_ct * main_block_size, &(g_write_dosage_val_bufs[1]))) {
        // this should be impossible
        assert(0);
        goto ox_bgen_to_pgen_ret_NOMEM;
      }

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
      unsigned char** compressed_geno_starts = g_compressed_geno_starts[0];
      unsigned char* bgen_geno_iter = compressed_geno_bufs[0];
      for (uint32_t variant_uidx = 0; variant_uidx < raw_variant_ct; ) {
        uint32_t uii;
        if (!fread_unlocked(&uii, 4, 1, bgenfile)) {
          goto ox_bgen_to_pgen_ret_READ_FAIL;
        }
        if (uii != sample_ct) {
          logprint("\n");
          logerrprint("Error: Unexpected number of samples specified in SNP block header.\n");
          goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
        }
        uint16_t snpid_slen;
        if (!fread_unlocked(&snpid_slen, 2, 1, bgenfile)) {
          goto ox_bgen_to_pgen_ret_READ_FAIL;
        }
        if (!snpid_chr) {
          if (fseeko(bgenfile, snpid_slen, SEEK_CUR)) {
            goto ox_bgen_to_pgen_ret_READ_FAIL;
          }
        } else {
          if (!snpid_slen) {
            logprint("\n");
            logerrprint("Error: Length-0 SNP ID in .bgen file.\n");
            goto ox_bgen_to_pgen_ret_INCONSISTENT_INPUT;
          }
          if (!fread_unlocked(loadbuf, snpid_slen, 1, bgenfile)) {
            goto ox_bgen_to_pgen_ret_READ_FAIL;
          }
          loadbuf[snpid_slen] = '\0';
        }
        uint16_t rsid_slen;
        if (!fread_unlocked(&rsid_slen, 2, 1, bgenfile)) {
          goto ox_bgen_to_pgen_ret_READ_FAIL;
        }
        if (fseeko(bgenfile, rsid_slen, SEEK_CUR)) {
          goto ox_bgen_to_pgen_ret_READ_FAIL;
        }
        uint16_t chr_name_slen;
        if (!fread_unlocked(&chr_name_slen, 2, 1, bgenfile)) {
          goto ox_bgen_to_pgen_ret_READ_FAIL;
        }
        if (!snpid_chr) {
          if (!chr_name_slen) {
            logprint("\n");
            logerrprint("Error: Length-0 chromosome ID in .bgen file.\n");
            goto ox_bgen_to_pgen_ret_INCONSISTENT_INPUT;
          }
          if (!fread_unlocked(loadbuf, chr_name_slen, 1, bgenfile)) {
            goto ox_bgen_to_pgen_ret_READ_FAIL;
          }
          if (strequal_k((char*)loadbuf, "NA", chr_name_slen)) {
            strcpy((char*)loadbuf, "0");
            chr_name_slen = 1;
          } else {
            loadbuf[chr_name_slen] = '\0';
          }
        } else {
          if (fseeko(bgenfile, chr_name_slen, SEEK_CUR)) {
            goto ox_bgen_to_pgen_ret_READ_FAIL;
          }
          chr_name_slen = snpid_slen;
        }
        int32_t cur_chr_code;
        reterr = get_or_add_chr_code_destructive("--bgen file", 0, allow_extra_chrs, (char*)loadbuf, (char*)(&(loadbuf[chr_name_slen])), cip, &cur_chr_code);
        if (reterr) {
          goto ox_bgen_to_pgen_ret_1;
        }
        const uint32_t skip = !is_set(cip->chr_mask, cur_chr_code);

        uint32_t cur_bp; // ignore in this pass
        if (!fread_unlocked(&cur_bp, 4, 1, bgenfile)) {
          goto ox_bgen_to_pgen_ret_READ_FAIL;
        }

        // allele count always 2 and not stored when layout=1
        for (uint32_t allele_idx = 0; allele_idx < 2; ++allele_idx) {
          uint32_t allele_slen;
          if (!fread_unlocked(&allele_slen, 4, 1, bgenfile)) {
            goto ox_bgen_to_pgen_ret_READ_FAIL;
          }
          if (fseeko(bgenfile, allele_slen, SEEK_CUR)) {
            goto ox_bgen_to_pgen_ret_READ_FAIL;
          }
        }

        if (compression_mode) {
#ifdef __LP64__
          compressed_block_byte_ct = 0;
#endif
          if (!fread_unlocked(&compressed_block_byte_ct, 4, 1, bgenfile)) {
            goto ox_bgen_to_pgen_ret_READ_FAIL;
          }
        }
        ++variant_uidx;
        if (!(variant_uidx % 1000)) {
          printf("\r--bgen: %uk variants scanned.", variant_uidx / 1000);
          fflush(stdout);
        }
        if (dosage_is_present || skip) {
          if (fseeko(bgenfile, compressed_block_byte_ct, SEEK_CUR)) {
            goto ox_bgen_to_pgen_ret_READ_FAIL;
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
        if (fread_checked(bgen_geno_iter, compressed_block_byte_ct, bgenfile)) {
          goto ox_bgen_to_pgen_ret_READ_FAIL;
        }
        bgen_geno_iter = &(bgen_geno_iter[compressed_block_byte_ct]);
        ++block_vidx;
        if (block_vidx == cur_block_size) {
          parity = 1 - parity;
          if (ts.thread_func_ptr) {
            // process *previous* block results
            join_threads3z(&ts);
            reterr = g_error_ret;
            if (reterr) {
              logprint("\n");
              logerrprint("Error: Invalid compressed SNP block in .bgen file.\n");
              goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
            }
            dosage_is_present = g_dosage_is_present;
            if (dosage_is_present) {
              // don't need to scan for any more dosages
              stop_threads3z(&ts, &g_cur_block_write_ct);
              if (!chr_filter_present) {
                break;
              }
              continue;
            }
          }
          g_cur_block_write_ct = cur_block_size;
          ts.thread_func_ptr = bgen11_dosage_scan_thread;
          if (spawn_threads3z(variant_ct, &ts)) {
            goto ox_bgen_to_pgen_ret_THREAD_CREATE_FAIL;
          }
          compressed_geno_starts = g_compressed_geno_starts[parity];
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
          if (!variant_ct) {
            logprint("\n");
            LOGERRPRINTFWW("Error: All %u variant%s in .bgen file skipped due to chromosome filter.\n", raw_variant_ct, (raw_variant_ct == 1)? "" : "s");
            goto ox_bgen_to_pgen_ret_INCONSISTENT_INPUT;
          }
          // bugfix (7 Oct 2017): with fewer variants than threads, need to
          // force initial launch here
          g_cur_block_write_ct = variant_ct;
          ts.thread_func_ptr = bgen11_dosage_scan_thread;
          if (spawn_threads3z(0, &ts)) {
            goto ox_bgen_to_pgen_ret_THREAD_CREATE_FAIL;
          }
          block_vidx = 0;
        }
      }
      if (ts.thread_func_ptr) {
        join_threads3z(&ts);
        reterr = g_error_ret;
        if (reterr) {
          logprint("\n");
          logerrprint("Error: Invalid compressed SNP block in .bgen file.\n");
          goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
        }
        if (block_vidx && (!g_dosage_is_present)) {
          g_cur_block_write_ct = block_vidx;
        } else {
          g_cur_block_write_ct = 0;
        }
        ts.is_last_block = 1;
        if (spawn_threads3z(1, &ts)) {
          goto ox_bgen_to_pgen_ret_THREAD_CREATE_FAIL;
        }
        join_threads3z(&ts);
        dosage_is_present = g_dosage_is_present;
      }

      if (fseeko(bgenfile, initial_uints[0] + 4, SEEK_SET)) {
        goto ox_bgen_to_pgen_ret_READ_FAIL;
      }
      strcpy(outname_end, ".pgen");
      uintptr_t spgw_alloc_cacheline_ct;
      uint32_t max_vrec_len;
      reterr = spgw_init_phase1(outname, nullptr, nullptr, variant_ct, sample_ct, dosage_is_present? kfPgenGlobalDosagePresent : kfPgenGlobal0, (oxford_import_flags & (kfOxfordImportRefFirst | kfOxfordImportRefSecond))? 1 : 2, &spgw, &spgw_alloc_cacheline_ct, &max_vrec_len);
      if (reterr) {
        goto ox_bgen_to_pgen_ret_1;
      }
      unsigned char* spgw_alloc;
      if (bigstack_alloc_uc(spgw_alloc_cacheline_ct * kCacheline, &spgw_alloc)) {
        goto ox_bgen_to_pgen_ret_NOMEM;
      }
      spgw_init_phase2(max_vrec_len, &spgw, spgw_alloc);

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
      uint32_t vidx_start = 0;
      uint32_t prev_block_write_ct = 0;
      parity = 0;
      reinit_threads3z(&ts);
      while (1) {
        uint32_t cur_block_write_ct = 0;
        if (!ts.is_last_block) {
          cur_block_write_ct = MINV(variant_ct - vidx_start, main_block_size);
          compressed_geno_starts = g_compressed_geno_starts[parity];
          bgen_geno_iter = compressed_geno_bufs[parity];
          for (block_vidx = 0; block_vidx < cur_block_write_ct;) {
            uint32_t uii;
            if (!fread_unlocked(&uii, 4, 1, bgenfile)) {
              goto ox_bgen_to_pgen_ret_READ_FAIL;
            }
            if (uii != sample_ct) {
              logprint("\n");
              logerrprint("Error: Unexpected number of samples specified in SNP block header.\n");
              goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
            }
            uint16_t snpid_slen;
            if (!fread_unlocked(&snpid_slen, 2, 1, bgenfile)) {
              goto ox_bgen_to_pgen_ret_READ_FAIL;
            }
            char* rsid_start = (char*)loadbuf;
            if (!snpid_chr) {
              if (fseeko(bgenfile, snpid_slen, SEEK_CUR)) {
                goto ox_bgen_to_pgen_ret_READ_FAIL;
              }
            } else {
              if (!snpid_slen) {
                logprint("\n");
                logerrprint("Error: Length-0 SNP ID in .bgen file.\n");
                goto ox_bgen_to_pgen_ret_INCONSISTENT_INPUT;
              }
              if (!fread_unlocked(loadbuf, snpid_slen, 1, bgenfile)) {
                goto ox_bgen_to_pgen_ret_READ_FAIL;
              }
              loadbuf[snpid_slen] = '\0';
              rsid_start = (char*)(&(loadbuf[snpid_slen + 1]));
            }
            uint16_t rsid_slen;
            if (!fread_unlocked(&rsid_slen, 2, 1, bgenfile)) {
              goto ox_bgen_to_pgen_ret_READ_FAIL;
            }
            if (!rsid_slen) {
              logprint("\n");
              logerrprint("Error: Length-0 rsID in .bgen file.\n");
              goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
            }
            if (!fread_unlocked(rsid_start, rsid_slen, 1, bgenfile)) {
              goto ox_bgen_to_pgen_ret_READ_FAIL;
            }
            char* loadbuf_iter = &(rsid_start[rsid_slen]);
            char* chr_name_start = loadbuf_iter;
            uint16_t chr_name_slen;
            if (!fread_unlocked(&chr_name_slen, 2, 1, bgenfile)) {
              goto ox_bgen_to_pgen_ret_READ_FAIL;
            }
            if (!snpid_chr) {
              if (!chr_name_slen) {
                logprint("\n");
                logerrprint("Error: Length-0 chromosome ID in .bgen file.\n");
                goto ox_bgen_to_pgen_ret_INCONSISTENT_INPUT;
              }
              if (!fread_unlocked(chr_name_start, chr_name_slen, 1, bgenfile)) {
                goto ox_bgen_to_pgen_ret_READ_FAIL;
              }
              if (strequal_k(chr_name_start, "NA", chr_name_slen)) {
                strcpy(chr_name_start, "0");
                chr_name_slen = 1;
              } else {
                chr_name_start[chr_name_slen] = '\0';
              }
            } else {
              if (fseeko(bgenfile, chr_name_slen, SEEK_CUR)) {
                goto ox_bgen_to_pgen_ret_READ_FAIL;
              }
              chr_name_start = (char*)loadbuf;
              chr_name_slen = snpid_slen;
            }
            int32_t cur_chr_code;
            reterr = get_or_add_chr_code_destructive("--bgen file", 0, allow_extra_chrs, (char*)chr_name_start, &(chr_name_start[chr_name_slen]), cip, &cur_chr_code);
            if (reterr) {
              goto ox_bgen_to_pgen_ret_1;
            }
            const uint32_t skip = !is_set(cip->chr_mask, cur_chr_code);

            uint32_t cur_bp;
            uint32_t a1_slen;
            if ((!fread_unlocked(&cur_bp, 4, 1, bgenfile)) ||
                (!fread_unlocked(&a1_slen, 4, 1, bgenfile))) {
              goto ox_bgen_to_pgen_ret_READ_FAIL;
            }
            if (skip) {
              uint32_t a2_slen;
              if (fseeko(bgenfile, a1_slen, SEEK_CUR) ||
                  (!fread_unlocked(&a2_slen, 4, 1, bgenfile)) ||
                  fseeko(bgenfile, a2_slen, SEEK_CUR)) {
                goto ox_bgen_to_pgen_ret_READ_FAIL;
              }
              if (compression_mode) {
#ifdef __LP64__
                compressed_block_byte_ct = 0;
#endif
                if (!fread_unlocked(&compressed_block_byte_ct, 4, 1, bgenfile)) {
                  goto ox_bgen_to_pgen_ret_READ_FAIL;
                }
              }
              if (fseeko(bgenfile, compressed_block_byte_ct, SEEK_CUR)) {
                goto ox_bgen_to_pgen_ret_READ_FAIL;
              }
              continue;
            }
            char* a1_ptr = loadbuf_iter;
            if (!a1_slen) {
              logprint("\n");
              logerrprint("Error: Empty allele code in .bgen file.\n");
              goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
            }
            if (a1_slen > 1000000000) {
              logprint("\n");
              logerrprint("Error: Allele code in .bgen file has more than 1 billion characters.\n");
              goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
            }
            if (a1_slen + (uintptr_t)(a1_ptr - ((char*)loadbuf)) > loadbuf_size) {
              goto ox_bgen_to_pgen_ret_NOMEM;
            }
            if (!fread_unlocked(a1_ptr, a1_slen, 1, bgenfile)) {
              goto ox_bgen_to_pgen_ret_READ_FAIL;
            }
            char* a2_ptr = &(a1_ptr[a1_slen]);
            uint32_t a2_slen;
            if (!fread_unlocked(&a2_slen, 4, 1, bgenfile)) {
              goto ox_bgen_to_pgen_ret_READ_FAIL;
            }
            if (!a2_slen) {
              logprint("\n");
              logerrprint("Error: Empty allele code in .bgen file.\n");
              goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
            }
            if (a2_slen > 1000000000) {
              logprint("\n");
              logerrprint("Error: Allele code in .bgen file has more than 1 billion characters.\n");
              goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
            }
            if (a2_slen + (uintptr_t)(a2_ptr - ((char*)loadbuf)) > loadbuf_size) {
              goto ox_bgen_to_pgen_ret_NOMEM;
            }
            if (!fread_unlocked(a2_ptr, a2_slen, 1, bgenfile)) {
              goto ox_bgen_to_pgen_ret_READ_FAIL;
            }
            if (compression_mode) {
#ifdef __LP64__
              compressed_block_byte_ct = 0;
#endif
              if (!fread_unlocked(&compressed_block_byte_ct, 4, 1, bgenfile)) {
                goto ox_bgen_to_pgen_ret_READ_FAIL;
              }
            }
            write_iter = chr_name_write(cip, cur_chr_code, write_iter);
            *write_iter++ = '\t';
            if (cur_bp > 0x7ffffffe) {
              logprint("\n");
              logerrprint("Error: Invalid bp coordinate (> 2^31 - 2) in .bgen file\n");
              goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
            }
            write_iter = uint32toa_x(cur_bp, '\t', write_iter);
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
              if (fwrite_checked(writebuf, (uintptr_t)(write_iter - writebuf), pvarfile)) {
                goto ox_bgen_to_pgen_ret_WRITE_FAIL;
              }
              write_iter = writebuf;
            }
            if (a1_slen < kMaxMediumLine) {
              write_iter = memcpya(write_iter, a1_ptr, a1_slen);
            } else {
              if (fwrite_checked(a1_ptr, a1_slen, pvarfile)) {
                goto ox_bgen_to_pgen_ret_WRITE_FAIL;
              }
            }
            *write_iter++ = '\t';
            if ((write_iter >= writebuf_flush) || (a2_slen >= kMaxMediumLine)) {
              if (fwrite_checked(writebuf, (uintptr_t)(write_iter - writebuf), pvarfile)) {
                goto ox_bgen_to_pgen_ret_WRITE_FAIL;
              }
              write_iter = writebuf;
            }
            if (a2_slen < kMaxMediumLine) {
              write_iter = memcpya(write_iter, a2_ptr, a2_slen);
            } else {
              if (fwrite_checked(a2_ptr, a2_slen, pvarfile)) {
                goto ox_bgen_to_pgen_ret_WRITE_FAIL;
              }
            }
            append_binary_eoln(&write_iter);

            compressed_geno_starts[block_vidx] = bgen_geno_iter;
            if (compression_mode) {
              memcpy(bgen_geno_iter, &compressed_block_byte_ct, 4);
              bgen_geno_iter = &(bgen_geno_iter[4]);
            }
            if (fread_checked(bgen_geno_iter, compressed_block_byte_ct, bgenfile)) {
              goto ox_bgen_to_pgen_ret_READ_FAIL;
            }
            bgen_geno_iter = &(bgen_geno_iter[compressed_block_byte_ct]);
            ++block_vidx;
          }
        }
        if (vidx_start) {
          join_threads3z(&ts);
          reterr = g_error_ret;
          if (reterr) {
            logprint("\n");
            logerrprint("Error: Invalid compressed SNP block in .bgen file.\n");
            goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
          }
        }
        if (!ts.is_last_block) {
          g_cur_block_write_ct = cur_block_write_ct;
          ts.is_last_block = (vidx_start + cur_block_write_ct == variant_ct);
          ts.thread_func_ptr = bgen11_geno_to_pgen_thread;
          if (spawn_threads3z(vidx_start, &ts)) {
            goto ox_bgen_to_pgen_ret_THREAD_CREATE_FAIL;
          }
        }
        parity = 1 - parity;
        if (vidx_start) {
          // write *previous* block results
          uintptr_t* write_genovec_iter = g_write_genovecs[parity];
          uint32_t* write_dosage_ct_iter = g_write_dosage_cts[parity];
          uintptr_t* write_dosage_present_iter = g_write_dosage_presents[parity];
          dosage_t* write_dosage_vals_iter = g_write_dosage_val_bufs[parity];
          for (uint32_t vidx = vidx_start - prev_block_write_ct; vidx < vidx_start; ++vidx) {
            const uint32_t cur_dosage_ct = *write_dosage_ct_iter++;
            if (!cur_dosage_ct) {
              if (spgw_append_biallelic_genovec(write_genovec_iter, &spgw)) {
                goto ox_bgen_to_pgen_ret_WRITE_FAIL;
              }
            } else {
              if (spgw_append_biallelic_genovec_dosage16(write_genovec_iter, write_dosage_present_iter, write_dosage_vals_iter, cur_dosage_ct, &spgw)) {
                goto ox_bgen_to_pgen_ret_WRITE_FAIL;
              }
            }
            write_genovec_iter = &(write_genovec_iter[sample_ctaw2]);
            write_dosage_present_iter = &(write_dosage_present_iter[sample_ctaw]);
            write_dosage_vals_iter = &(write_dosage_vals_iter[sample_ct]);
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
      if (bigstack_end_alloc_ul(raw_variant_ct + 1, &allele_idx_offsets)) {
        goto ox_bgen_to_pgen_ret_NOMEM;
      }

      g_bgen_import_dosage_certainty_thresholds = nullptr;
      if (import_dosage_certainty > (1.0 - kSmallEpsilon) / 3.0) {
        g_bgen_import_dosage_certainty_thresholds = (uint32_t*)bigstack_alloc_raw_rd(25 * sizeof(int32_t));
        for (uint32_t bit_precision = 1; bit_precision <= 16; ++bit_precision) {
          const uint32_t denom = (1U << bit_precision) - 1;
          g_bgen_import_dosage_certainty_thresholds[bit_precision] = 1 + (int32_t)(import_dosage_certainty * ((int32_t)denom));
        }
      }
      // bugfix (2 Jul 2017): if max_thread_ct == 1 but there's >12GB memory,
      //   limit to 1 thread rather than (max_thread_ct - 1)...
      uint32_t calc_thread_ct_limit = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
      if (calc_thread_ct_limit > raw_variant_ct) {
        calc_thread_ct_limit = raw_variant_ct;
      }
      if (bigstack_alloc_thread(calc_thread_ct_limit, &ts.threads)) {
        goto ox_bgen_to_pgen_ret_NOMEM;
      }

      g_thread_wkspaces = (unsigned char**)bigstack_alloc_raw_rd(calc_thread_ct_limit * sizeof(intptr_t));
      g_thread_bidxs[0] = (uint32_t*)bigstack_alloc_raw_rd((calc_thread_ct_limit + 1) * sizeof(int32_t));
      g_thread_bidxs[1] = (uint32_t*)bigstack_alloc_raw_rd((calc_thread_ct_limit + 1) * sizeof(int32_t));
      // ***** all allocations from this point on are reset before pass 2 *****
      uintptr_t main_block_size = 65536;
      if (bigstack_alloc_usi(main_block_size, &(g_bgen_allele_cts[0])) ||
          bigstack_alloc_usi(main_block_size, &(g_bgen_allele_cts[1])) ||
          bigstack_alloc_ucp(main_block_size + 1, &(g_compressed_geno_starts[0])) ||
          bigstack_alloc_ucp(main_block_size + 1, &(g_compressed_geno_starts[1]))) {
        goto ox_bgen_to_pgen_ret_NOMEM;
      }
      if (compression_mode) {
        if (bigstack_alloc_ui(main_block_size, &(g_uncompressed_genodata_byte_cts[0])) ||
            bigstack_alloc_ui(main_block_size, &(g_uncompressed_genodata_byte_cts[1]))) {
          goto ox_bgen_to_pgen_ret_NOMEM;
        }
      } else {
        // defensive
        g_uncompressed_genodata_byte_cts[0] = nullptr;
        g_uncompressed_genodata_byte_cts[1] = nullptr;
      }

      // ploidy >2 is not supported by PLINK 2.  (A future build may have code
      // to treat those calls as missing instead of erroring out, as is done
      // with VCF ploidy >2.  But I'll wait until this case actually comes up
      // in the wild...)
      // But even without that, the diploid worst case of 65535 alleles,
      // unphased 32-bit probabilities blows past the 4GB uncompressed record
      // size limit with just 1 sample!  Consequences:
      // * A simple way to avoid unnecessary NOMEM errors is to give each
      //   thread 4GB of decompression workspace on the first pass.  This may
      //   greatly reduce the number of decompression worker threads we can
      //   deploy, but for the first pass that's acceptable: the worker threads
      //   will usually all exit almost immediately (since we just need to
      //   determine whether *any* phase/dosage info needs to be saved).
      // * Even 1 thread x 4GB won't always be available, especially since we
      //   have a double-buffering workflow which requires additional
      //   allocations summing to more than twice the decompression workspace.
      //   So we need to be able to fall back to a smaller decompression
      //   workspace size, and throw NOMEM when it's insufficient.
      // * Of course, records will almost always be far smaller than 4GB.
      //   During the first pass, we'll see every uncompressed record size
      //   (even if the decompression worker threads terminate early), so we
      //   can usually increase the number of worker threads before the second
      //   pass.
      // Overall memory allocation for first pass:
      //   loadbuf_size (~1/7, up to 2GB) : Chromosome code/variant ID/allele
      //                                    code load buffer.
      //   mainbuf_size (~2/7) : Compressed genotype data buffer 0, up to 4GB
      //                         per decompression thread
      //   mainbuf_size        : Compressed genotype data buffer 1
      //   mainbuf_size        : Decompression thread workspace(s)
      // Second pass:
      //   mainbuf_size (~1/6) : Decompression thread workspaces.
      //   16K                 : .bgen chromosome code load buffer.
      //   remainder (~5/6)    : Compressed genotype data buffers, writer, and
      //                         write buffers.
      uintptr_t loadbuf_size = round_down_pow2(bigstack_left() / 7, kCacheline);
      if (loadbuf_size > kMaxLongLine) {
        loadbuf_size = kMaxLongLine;
      } else if (loadbuf_size < 2 * 65536) {
        // don't want to worry about chromosome/variant ID buffer space checks
        // in inner loop
        goto ox_bgen_to_pgen_ret_NOMEM;
      }
      unsigned char* loadbuf = bigstack_alloc_raw(loadbuf_size);

      uintptr_t mainbuf_size = round_down_pow2(bigstack_left() / 3, kCacheline);
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
          mainbuf_size = ((uintptr_t)calc_thread_ct_limit) << 32;
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
      ts.calc_thread_ct = calc_thread_ct;
      g_calc_thread_ct = calc_thread_ct;
      unsigned char* compressed_geno_bufs[2];
      compressed_geno_bufs[0] = bigstack_alloc_raw(mainbuf_size);
      compressed_geno_bufs[1] = bigstack_alloc_raw(mainbuf_size);
      for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
        g_thread_wkspaces[tidx] = bigstack_alloc_raw(thread_wkspace_size);
      }

      uint32_t variant_ct = 0;

      uint32_t block_vidx = 0;

      // bgen-1.2 and -1.3 records can vary wildly in size, so we're a bit more
      // careful with load balancing here.
      uint32_t cur_per_thread_block_limit = 1;
      uint32_t cur_thread_block_vidx_limit = 1;
      uint32_t cur_thread_fill_idx = 0;

      uint32_t parity = 0;
      uint32_t* thread_bidxs = g_thread_bidxs[0];
      uint16_t* bgen_allele_cts = g_bgen_allele_cts[0];
      unsigned char** compressed_geno_starts = g_compressed_geno_starts[0];
      uint32_t* uncompressed_genodata_byte_cts = g_uncompressed_genodata_byte_cts[0];
      unsigned char* bgen_geno_iter = compressed_geno_bufs[0];
      unsigned char* cur_geno_buf_end = &(bgen_geno_iter[thread_wkspace_size]);
      thread_bidxs[0] = 0;
      compressed_geno_starts[0] = bgen_geno_iter;
      uintptr_t* allele_idx_offsets_iter = allele_idx_offsets;
      uintptr_t tot_allele_ct = 0;
      uint32_t max_geno_blen = 0;
      uint32_t uncompressed_genodata_byte_ct = 0;

      // temporary kludge
      uint32_t multiallelic_skip_ct = 0;

      g_cur_block_write_ct = 1; // just used as a flag

      for (uint32_t variant_uidx = 0; variant_uidx < raw_variant_ct; ) {
        // format is mostly identical to bgen 1.1; but there's no sample count,
        // and there is an allele count
        // logic is more similar to the second bgen 1.1 pass since we write the
        // .pvar here.
        uint16_t snpid_slen;
        if (!fread_unlocked(&snpid_slen, 2, 1, bgenfile)) {
          goto ox_bgen_to_pgen_ret_READ_FAIL;
        }
        char* rsid_start = (char*)loadbuf;
        if (!snpid_chr) {
          if (fseeko(bgenfile, snpid_slen, SEEK_CUR)) {
            goto ox_bgen_to_pgen_ret_READ_FAIL;
          }
        } else {
          if (!snpid_slen) {
            logprint("\n");
            logerrprint("Error: Length-0 SNP ID in .bgen file.\n");
            goto ox_bgen_to_pgen_ret_INCONSISTENT_INPUT;
          }
          if (!fread_unlocked(loadbuf, snpid_slen, 1, bgenfile)) {
            goto ox_bgen_to_pgen_ret_READ_FAIL;
          }
          loadbuf[snpid_slen] = '\0';
          rsid_start = (char*)(&(loadbuf[snpid_slen + 1]));
        }
        uint16_t rsid_slen;
        if (!fread_unlocked(&rsid_slen, 2, 1, bgenfile)) {
          goto ox_bgen_to_pgen_ret_READ_FAIL;
        }
        if (!rsid_slen) {
          logprint("\n");
          logerrprint("Error: Length-0 rsID in .bgen file.\n");
          goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
        }
        if (!fread_unlocked(rsid_start, rsid_slen, 1, bgenfile)) {
          goto ox_bgen_to_pgen_ret_READ_FAIL;
        }
        char* loadbuf_iter = &(rsid_start[rsid_slen]);
        char* chr_name_start = loadbuf_iter;
        uint16_t chr_name_slen;
        if (!fread_unlocked(&chr_name_slen, 2, 1, bgenfile)) {
          goto ox_bgen_to_pgen_ret_READ_FAIL;
        }
        if (!snpid_chr) {
          if (!chr_name_slen) {
            logprint("\n");
            logerrprint("Error: Length-0 chromosome ID in .bgen file.\n");
            goto ox_bgen_to_pgen_ret_INCONSISTENT_INPUT;
          }
          if (!fread_unlocked(chr_name_start, chr_name_slen, 1, bgenfile)) {
            goto ox_bgen_to_pgen_ret_READ_FAIL;
          }
          if (strequal_k(chr_name_start, "NA", chr_name_slen)) {
            strcpy(chr_name_start, "0");
            chr_name_slen = 1;
          } else {
            chr_name_start[chr_name_slen] = '\0';
          }
        } else {
          if (fseeko(bgenfile, chr_name_slen, SEEK_CUR)) {
            goto ox_bgen_to_pgen_ret_READ_FAIL;
          }
          chr_name_start = (char*)loadbuf;
          chr_name_slen = snpid_slen;
        }
        // chromosome ID length restriction enforced here, so we don't check
        // earlier
        int32_t cur_chr_code;
        reterr = get_or_add_chr_code_destructive("--bgen file", 0, allow_extra_chrs, (char*)chr_name_start, &(chr_name_start[chr_name_slen]), cip, &cur_chr_code);
        if (reterr) {
          goto ox_bgen_to_pgen_ret_1;
        }
        uint32_t skip = !is_set(cip->chr_mask, cur_chr_code);

        uint32_t cur_bp;
        uint32_t cur_allele_ct = 0;
        if ((!fread_unlocked(&cur_bp, 4, 1, bgenfile)) ||
            (!fread_unlocked(&cur_allele_ct, 2, 1, bgenfile))) {
          goto ox_bgen_to_pgen_ret_READ_FAIL;
        }
        if (cur_allele_ct < 2) {
          // this is undefined in the 1.3 standard; prohibit for now
          logprint("\n");
          logerrprint("Error: .bgen variant has fewer than two alleles.\n");
          goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
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
          for (uint32_t allele_idx = 0; allele_idx < cur_allele_ct; ++allele_idx) {
            uint32_t cur_allele_slen;
            if ((!fread_unlocked(&cur_allele_slen, 4, 1, bgenfile)) ||
                fseeko(bgenfile, cur_allele_slen, SEEK_CUR)) {
              goto ox_bgen_to_pgen_ret_READ_FAIL;
            }
          }
          uint32_t genodata_byte_ct;
          if ((!fread_unlocked(&genodata_byte_ct, 4, 1, bgenfile)) ||
              fseeko(bgenfile, genodata_byte_ct, SEEK_CUR)) {
            goto ox_bgen_to_pgen_ret_READ_FAIL;
          }
          continue;
        }
        if (rsid_slen > kMaxIdSlen) {
          // enforce this iff we aren't skipping
          logprint("\n");
          logerrprint("Error: Variant names are limited to " MAX_ID_SLEN_STR " characters.\n");
          goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
        }
        // special handling of first two alleles since either may be
        // reference, so we may need to swap order
        char* a1_ptr = loadbuf_iter;
        uint32_t a1_slen;
        if (!fread_unlocked(&a1_slen, 4, 1, bgenfile)) {
          goto ox_bgen_to_pgen_ret_READ_FAIL;
        }
        if (!a1_slen) {
          logprint("\n");
          logerrprint("Error: Empty allele code in .bgen file.\n");
          goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
        }
        if (a1_slen > 1000000000) {
          logprint("\n");
          logerrprint("Error: Allele code in .bgen file has more than 1 billion characters.\n");
          goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
        }
        if (a1_slen + (uintptr_t)(a1_ptr - ((char*)loadbuf)) > loadbuf_size) {
          goto ox_bgen_to_pgen_ret_NOMEM;
        }
        if (!fread_unlocked(a1_ptr, a1_slen, 1, bgenfile)) {
          goto ox_bgen_to_pgen_ret_READ_FAIL;
        }
        char* a2_ptr = &(a1_ptr[a1_slen]);
        uint32_t a2_slen;
        if (!fread_unlocked(&a2_slen, 4, 1, bgenfile)) {
          goto ox_bgen_to_pgen_ret_READ_FAIL;
        }
        if (!a2_slen) {
          logprint("\n");
          logerrprint("Error: Empty allele code in .bgen file.\n");
          goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
        }
        if (a2_slen > 1000000000) {
          logprint("\n");
          logerrprint("Error: Allele code in .bgen file has more than 1 billion characters.\n");
          goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
        }
        if (a2_slen + (uintptr_t)(a2_ptr - ((char*)loadbuf)) > loadbuf_size) {
          goto ox_bgen_to_pgen_ret_NOMEM;
        }
        if (!fread_unlocked(a2_ptr, a2_slen, 1, bgenfile)) {
          goto ox_bgen_to_pgen_ret_READ_FAIL;
        }
        write_iter = chr_name_write(cip, cur_chr_code, write_iter);
        *write_iter++ = '\t';
        if (cur_bp > 0x7ffffffe) {
          logprint("\n");
          logerrprint("Error: Invalid bp coordinate (> 2^31 - 2) in .bgen file\n");
          goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
        }
        write_iter = uint32toa_x(cur_bp, '\t', write_iter);
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
          if (fwrite_checked(writebuf, (uintptr_t)(write_iter - writebuf), pvarfile)) {
            goto ox_bgen_to_pgen_ret_WRITE_FAIL;
          }
          write_iter = writebuf;
        }
        if (a1_slen < kMaxMediumLine) {
          write_iter = memcpya(write_iter, a1_ptr, a1_slen);
        } else {
          if (fwrite_checked(a1_ptr, a1_slen, pvarfile)) {
            goto ox_bgen_to_pgen_ret_WRITE_FAIL;
          }
        }
        *write_iter++ = '\t';
        if ((write_iter >= writebuf_flush) || (a2_slen >= kMaxMediumLine)) {
          if (fwrite_checked(writebuf, (uintptr_t)(write_iter - writebuf), pvarfile)) {
            goto ox_bgen_to_pgen_ret_WRITE_FAIL;
          }
          write_iter = writebuf;
        }
        if (a2_slen < kMaxMediumLine) {
          write_iter = memcpya(write_iter, a2_ptr, a2_slen);
        } else {
          if (fwrite_checked(a2_ptr, a2_slen, pvarfile)) {
            goto ox_bgen_to_pgen_ret_WRITE_FAIL;
          }
        }
        for (uint32_t allele_idx = 2; allele_idx < cur_allele_ct; ++allele_idx) {
          // (can't actually reach here yet since we're skipping multiallelics
          // for now)
          // safe to use entire loadbuf for this
          assert(0);
          uint32_t cur_allele_slen;
          if (!fread_unlocked(&cur_allele_slen, 4, 1, bgenfile)) {
            goto ox_bgen_to_pgen_ret_READ_FAIL;
          }
          if (!cur_allele_slen) {
            logprint("\n");
            logerrprint("Error: Empty allele code in .bgen file.\n");
            goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
          }
          if (cur_allele_slen > 1000000000) {
            logprint("\n");
            logerrprint("Error: Allele code in .bgen file has more than 1 billion characters.\n");
            goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
          }
          if (cur_allele_slen > loadbuf_size) {
            goto ox_bgen_to_pgen_ret_NOMEM;
          }
          if (!fread_unlocked(loadbuf, cur_allele_slen, 1, bgenfile)) {
            goto ox_bgen_to_pgen_ret_READ_FAIL;
          }
          *write_iter++ = ',';
          if ((write_iter >= writebuf_flush) || (cur_allele_slen >= kMaxMediumLine)) {
            if (fwrite_checked(writebuf, (uintptr_t)(write_iter - writebuf), pvarfile)) {
              goto ox_bgen_to_pgen_ret_WRITE_FAIL;
            }
            write_iter = writebuf;
          }
          if (cur_allele_slen < kMaxMediumLine) {
            write_iter = memcpya(write_iter, loadbuf, cur_allele_slen);
          } else {
            if (fwrite_checked(loadbuf, cur_allele_slen, pvarfile)) {
              goto ox_bgen_to_pgen_ret_WRITE_FAIL;
            }
          }
        }

        append_binary_eoln(&write_iter);
        *allele_idx_offsets_iter++ = tot_allele_ct;
        tot_allele_ct += cur_allele_ct;
        uint32_t genodata_byte_ct;
        if (!fread_unlocked(&genodata_byte_ct, 4, 1, bgenfile)) {
          goto ox_bgen_to_pgen_ret_READ_FAIL;
        }
        if (genodata_byte_ct > max_geno_blen) {
          max_geno_blen = genodata_byte_ct;
        }
        if (uncompressed_genodata_byte_cts) {
          if (genodata_byte_ct < 4) {
            logerrprint("Error: Invalid compressed block length in .bgen file.\n");
            goto ox_bgen_to_pgen_ret_MALFORMED_INPUT;
          }
          if (!fread_unlocked(&uncompressed_genodata_byte_ct, 4, 1, bgenfile)) {
            goto ox_bgen_to_pgen_ret_READ_FAIL;
          }
          if (uncompressed_genodata_byte_ct > max_geno_blen) {
            max_geno_blen = uncompressed_genodata_byte_ct;
          }
          genodata_byte_ct -= 4;
        }
        if (dosage_is_present) {
          if (fseeko(bgenfile, genodata_byte_ct, SEEK_CUR)) {
            goto ox_bgen_to_pgen_ret_READ_FAIL;
          }
          ++variant_ct;
          continue;
        }

        if ((block_vidx == cur_thread_block_vidx_limit) || ((uintptr_t)(cur_geno_buf_end - bgen_geno_iter) < genodata_byte_ct)) {
          if (!block_vidx) {
            goto ox_bgen_to_pgen_ret_NOMEM;
          }
          thread_bidxs[++cur_thread_fill_idx] = block_vidx;
          if (cur_thread_fill_idx == calc_thread_ct) {
            parity = 1 - parity;
            if (ts.thread_func_ptr) {
              // process *previous* block results
              join_threads3z(&ts);
              reterr = g_error_ret;
              if (reterr) {
                goto ox_bgen_to_pgen_ret_bgen13_thread_fail;
              }
              dosage_is_present = g_dosage_is_present;
              if (dosage_is_present) {
                // don't need to scan for any more dosages
                stop_threads3z(&ts, &g_cur_block_write_ct);

                // however, unlike bgen-1.1 case, we can never do full
                // early-exit since we have to scan for multiallelic variants:
                // writer must be initialized with (i) an accurate variant
                // count, which is affected by skipped multiallelic variants,
                // and (ii) when we no longer skip them, the PgenWriter
                // constructor still needs a maximum allele count so it can
                // allocate properly-sized buffers.
                if (fseeko(bgenfile, genodata_byte_ct, SEEK_CUR)) {
                  goto ox_bgen_to_pgen_ret_READ_FAIL;
                }
                ++variant_ct;
                continue;
              }
            }
            ts.thread_func_ptr = bgen13_dosage_or_phase_scan_thread;
            if (spawn_threads3z(variant_ct, &ts)) {
              goto ox_bgen_to_pgen_ret_THREAD_CREATE_FAIL;
            }
            compressed_geno_starts = g_compressed_geno_starts[parity];
            uncompressed_genodata_byte_cts = g_uncompressed_genodata_byte_cts[parity];
            thread_bidxs = g_thread_bidxs[parity];
            bgen_allele_cts = g_bgen_allele_cts[parity];
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
        if (fread_checked(bgen_geno_iter, genodata_byte_ct, bgenfile)) {
          goto ox_bgen_to_pgen_ret_READ_FAIL;
        }
        bgen_geno_iter = &(bgen_geno_iter[genodata_byte_ct]);
        compressed_geno_starts[++block_vidx] = bgen_geno_iter;
      }
      variant_ct += block_vidx;
      if (multiallelic_skip_ct) {
        logprint("\n");
        LOGERRPRINTFWW("Warning: %u multiallelic variant%s skipped (not yet supported).\n", multiallelic_skip_ct, (multiallelic_skip_ct == 1)? "" : "s");
      }
      if (!variant_ct) {
        logprint("\n");
        LOGERRPRINTF("Error: All %u variant%s in .bgen file skipped.\n", raw_variant_ct, (raw_variant_ct == 1)? "" : "s");
        goto ox_bgen_to_pgen_ret_INCONSISTENT_INPUT;
      }
      if (variant_ct == block_vidx) {
        // with multiple threads, there's no guarantee that even the first
        // decompression job has launched (e.g. there's only 1 variant on the
        // relevant chromosome in the entire .bgen, and calc_thread_ct == 2).
        thread_bidxs[cur_thread_fill_idx + 1] = block_vidx;
        ts.thread_func_ptr = bgen13_dosage_or_phase_scan_thread;
        if (spawn_threads3z(variant_ct, &ts)) {
          goto ox_bgen_to_pgen_ret_THREAD_CREATE_FAIL;
        }
        block_vidx = 0;
      }
      chr_filter_present = (variant_ct + multiallelic_skip_ct != raw_variant_ct);
      if (ts.thread_func_ptr) {
        join_threads3z(&ts);
        reterr = g_error_ret;
        if (reterr) {
          goto ox_bgen_to_pgen_ret_bgen13_thread_fail;
        }
        if ((!block_vidx) || g_dosage_is_present) {
          // ignore thread_bidxs[] in this case
          g_cur_block_write_ct = 0;
        } else {
          for (; cur_thread_fill_idx < calc_thread_ct; ) {
            // save endpoint for current thread, and tell any leftover threads
            // to do nothing
            thread_bidxs[++cur_thread_fill_idx] = block_vidx;
          }
        }
        ts.is_last_block = 1;
        if (spawn_threads3z(1, &ts)) {
          goto ox_bgen_to_pgen_ret_THREAD_CREATE_FAIL;
        }
        join_threads3z(&ts);
        dosage_is_present = g_dosage_is_present;
      }

      if (tot_allele_ct == variant_ct * 2) {
        allele_idx_offsets = nullptr;
        bigstack_end_reset(bigstack_end_mark);
      } else {
        // not yet possible
        assert(0);
        *allele_idx_offsets_iter = tot_allele_ct;
      }
      if (fseeko(bgenfile, initial_uints[0] + 4, SEEK_SET)) {
        goto ox_bgen_to_pgen_ret_READ_FAIL;
      }
      strcpy(outname_end, ".pgen");
      uintptr_t spgw_alloc_cacheline_ct;
      uint32_t max_vrec_len;
      reterr = spgw_init_phase1(outname, allele_idx_offsets, nullptr, variant_ct, sample_ct, dosage_is_present? (kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent) : kfPgenGlobal0, (oxford_import_flags & (kfOxfordImportRefFirst | kfOxfordImportRefSecond))? 1 : 2, &spgw, &spgw_alloc_cacheline_ct, &max_vrec_len);
      if (reterr) {
        goto ox_bgen_to_pgen_ret_1;
      }

      bigstack_reset(g_bgen_allele_cts[0]);

      // only needs to fit chromosome codes in second pass
      loadbuf = bigstack_alloc_raw_rd(kMaxIdBlen);
      unsigned char* spgw_alloc;
      if (bigstack_alloc_uc(spgw_alloc_cacheline_ct * kCacheline, &spgw_alloc)) {
        goto ox_bgen_to_pgen_ret_NOMEM;
      }
      // Now that we know max_geno_blen, try to increase calc_thread_ct, and
      // resize g_thread_wkspaces[tidx] (and also resize compressed_geno_bufs[]
      // in next step).
      // Additional *6 in denominator since we want to limit these allocations
      // to 1/6 of remaining workspace.
      thread_wkspace_size = round_up_pow2(max_geno_blen, kCacheline);
      // bugfix (16 Jul 2017): was computing cachelines_avail, not bytes_avail
      uintptr_t bytes_avail = round_down_pow2(bigstack_left() / 6, kCacheline);
      if (calc_thread_ct_limit * thread_wkspace_size <= bytes_avail) {
        calc_thread_ct = calc_thread_ct_limit;
      } else {
        calc_thread_ct = bytes_avail / thread_wkspace_size;
        if (!calc_thread_ct) {
          goto ox_bgen_to_pgen_ret_NOMEM;
        }
      }
      ts.calc_thread_ct = calc_thread_ct;
      for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
        g_thread_wkspaces[tidx] = bigstack_alloc_raw(thread_wkspace_size);
      }
      bytes_avail -= thread_wkspace_size * calc_thread_ct;
      // Per-write-buffer-variant allocations:
      //   g_bgen_allele_cts: 2 * sizeof(int16_t)
      //   g_compressed_geno_starts: 2 * sizeof(intptr_t)
      //   g_uncompressed_genodata_byte_cts: 2 * sizeof(int32_t)
      //     (unless compression_mode == 0)
      //   g_write_genovecs: 2 * sample_ctaw2 * sizeof(intptr_t)
      //   g_write_phasepresents: 2 * sample_ctaw * sizeof(intptr_t)
      //   g_write_phaseinfos: 2 * sample_ctaw * sizeof(intptr_t)
      //   g_write_dosage_presents: 2 * sample_ctaw * sizeof(intptr_t)
      //   g_write_dphase_presents: 2 * sample_ctaw * sizeof(intptr_t)
      //   g_write_dosage_cts: 2 * sizeof(int32_t)
      //   g_write_dphase_cts: 2 * sizeof(int32_t)
      //   g_write_dosage_val_bufs (the big one): 2 * sample_ct * 2 *
      //                                          sizeof(dosage_t)
      //     additional factor of 2 here is due to phased-dosage support.  (not
      //     actually implemented yet, but will be soon.)
      //   g_compressed_geno_bufs (the other big one): up to 2 * max_geno_blen
      //     The "up to" here is due to the possibility that a few variants
      //     require much more space than the rest; unlikely now, but will be
      //     important when multiallelic support is added.  To defend against
      //     that possibility, we limit g_compressed_geno_bufs[] to 50% of the
      //     total allocation here.
      uintptr_t cachelines_avail_m24 = bigstack_left() / kCacheline;
      if (cachelines_avail_m24 < 24) {
        goto ox_bgen_to_pgen_ret_NOMEM;
      }
      // we're making up to 24 allocations; be pessimistic re: rounding
      // (g_compressed_geno_starts has +1, but we have enough room for error)
      cachelines_avail_m24 -= 24;
      uintptr_t bytes_req_per_in_block_variant = 2 * (sizeof(int16_t) + sizeof(intptr_t) + sample_ctaw2 * sizeof(intptr_t) + sample_ctaw * 4 * sizeof(intptr_t) + 2 * sizeof(int32_t) + sample_ct * 2 * sizeof(dosage_t));
      if (compression_mode) {
        bytes_req_per_in_block_variant += 2 * sizeof(int32_t);
      }
      // 50% cap
      mainbuf_size = MINV(max_geno_blen, bytes_req_per_in_block_variant / 2);
      // bugfix (16 Jul 2017): forgot to include this term
      // (17 Jul 2017): forgot to multiply by 2
      bytes_req_per_in_block_variant += 2 * mainbuf_size;
      main_block_size = (cachelines_avail_m24 * kCacheline) / bytes_req_per_in_block_variant;
      if (main_block_size > 65536) {
        main_block_size = 65536;
      }
      if (main_block_size > raw_variant_ct + calc_thread_ct - 1) {
        main_block_size = raw_variant_ct + calc_thread_ct - 1;
      }
      per_thread_block_limit = main_block_size / calc_thread_ct;
      // may as well guarantee divisibility
      main_block_size = per_thread_block_limit * calc_thread_ct;
      mainbuf_size *= main_block_size;
      if (mainbuf_size < max_geno_blen) {
        // bugfix (2 Jul 2017): don't error out here if the entire .bgen has
        // e.g. only one variant
        goto ox_bgen_to_pgen_ret_NOMEM;
      }
      if (bigstack_alloc_usi(main_block_size, &(g_bgen_allele_cts[0])) ||
          bigstack_alloc_usi(main_block_size, &(g_bgen_allele_cts[1])) ||
          bigstack_alloc_ucp(main_block_size + 1, &(g_compressed_geno_starts[0])) ||
          bigstack_alloc_ucp(main_block_size + 1, &(g_compressed_geno_starts[1])) ||
          bigstack_alloc_uc(mainbuf_size, &(compressed_geno_bufs[0])) ||
          bigstack_alloc_uc(mainbuf_size, &(compressed_geno_bufs[1])) ||
          bigstack_alloc_ul(sample_ctaw2 * main_block_size, &(g_write_genovecs[0])) ||
          bigstack_alloc_ul(sample_ctaw2 * main_block_size, &(g_write_genovecs[1])) ||
          bigstack_alloc_ul(sample_ctaw * main_block_size, &(g_write_phasepresents[0])) ||
          bigstack_alloc_ul(sample_ctaw * main_block_size, &(g_write_phasepresents[1])) ||
          bigstack_alloc_ul(sample_ctaw * main_block_size, &(g_write_phaseinfos[0])) ||
          bigstack_alloc_ul(sample_ctaw * main_block_size, &(g_write_phaseinfos[1])) ||
          bigstack_alloc_ul(sample_ctaw * main_block_size, &(g_write_dosage_presents[0])) ||
          bigstack_alloc_ul(sample_ctaw * main_block_size, &(g_write_dosage_presents[1])) ||
          bigstack_alloc_ul(sample_ctaw * main_block_size, &(g_write_dphase_presents[0])) ||
          bigstack_alloc_ul(sample_ctaw * main_block_size, &(g_write_dphase_presents[1])) ||
          bigstack_alloc_ui(main_block_size, &(g_write_dosage_cts[0])) ||
          bigstack_alloc_ui(main_block_size, &(g_write_dosage_cts[1])) ||
          bigstack_alloc_ui(main_block_size, &(g_write_dphase_cts[0])) ||
          bigstack_alloc_ui(main_block_size, &(g_write_dphase_cts[1])) ||
          bigstack_alloc_dosage(sample_ct * 2 * main_block_size, &(g_write_dosage_val_bufs[0])) ||
          bigstack_alloc_dosage(sample_ct * 2 * main_block_size, &(g_write_dosage_val_bufs[1]))) {
        // this should be impossible
        assert(0);
        goto ox_bgen_to_pgen_ret_NOMEM;
      }
      if (compression_mode) {
        if (bigstack_alloc_ui(main_block_size, &(g_uncompressed_genodata_byte_cts[0])) ||
            bigstack_alloc_ui(main_block_size, &(g_uncompressed_genodata_byte_cts[1]))) {
          assert(0);
          goto ox_bgen_to_pgen_ret_NOMEM;
        }
      }
      spgw_init_phase2(max_vrec_len, &spgw, spgw_alloc);

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
      // (May be better to change this to use one output buffer instead of 2,
      // due to high memory requirement.)
      uint32_t vidx_start = 0;
      uint32_t prev_block_write_ct = 0;
      uint32_t prev_genodata_byte_ct = 0;
      uint32_t prev_allele_ct = 0;
      uint32_t skip = 0;
      parity = 0;
      reinit_threads3z(&ts);
      g_cur_block_write_ct = 1;
      while (1) {
        uint32_t cur_block_write_ct = 0;
        if (!ts.is_last_block) {
          const uint32_t block_vidx_limit = variant_ct - vidx_start;
          cur_thread_block_vidx_limit = MINV(block_vidx_limit, per_thread_block_limit);
          cur_thread_fill_idx = 0;
          thread_bidxs = g_thread_bidxs[parity];
          bgen_allele_cts = g_bgen_allele_cts[parity];
          compressed_geno_starts = g_compressed_geno_starts[parity];
          uncompressed_genodata_byte_cts = g_uncompressed_genodata_byte_cts[parity];
          bgen_geno_iter = compressed_geno_bufs[parity];
          cur_geno_buf_end = &(bgen_geno_iter[thread_wkspace_size]);
          thread_bidxs[0] = 0;
          compressed_geno_starts[0] = bgen_geno_iter;
          block_vidx = 0;
          // strictly speaking, prev_genodata_byte_ct and genodata_byte_ct
          // can be collapsed into one variable, as well as
          // {block_vidx, cur_block_write_ct}, but not a big deal if the
          // compiler fails to see this
          uint32_t genodata_byte_ct = prev_genodata_byte_ct;
          uint32_t cur_allele_ct = prev_allele_ct;
          if (!genodata_byte_ct) {
            goto ox_bgen_to_pgen_load13_start;
          }
          // we may stop before main_block_size due to insufficient space in
          // compressed_geno_buf.  if so, the file pointer is right before the
          // genotype data, rather than at the beginning of a variant record.
          while (1) {
            bgen_allele_cts[block_vidx] = cur_allele_ct;
            if (uncompressed_genodata_byte_cts) {
              uncompressed_genodata_byte_cts[block_vidx] = uncompressed_genodata_byte_ct;
            }
            if (fread_checked(bgen_geno_iter, genodata_byte_ct, bgenfile)) {
              goto ox_bgen_to_pgen_ret_READ_FAIL;
            }
            bgen_geno_iter = &(bgen_geno_iter[genodata_byte_ct]);
            compressed_geno_starts[++block_vidx] = bgen_geno_iter;

            uint16_t snpid_slen;
            // true iff this is the last variant we're keeping in the entire
            // file
            if (block_vidx == block_vidx_limit) {
              for (; cur_thread_fill_idx < calc_thread_ct; ) {
                // save endpoint for current thread, and tell any leftover
                // threads to do nothing
                thread_bidxs[++cur_thread_fill_idx] = block_vidx;
              }
              break;
            }
          ox_bgen_to_pgen_load13_start:
            if (!fread_unlocked(&snpid_slen, 2, 1, bgenfile)) {
              goto ox_bgen_to_pgen_ret_READ_FAIL;
            }

            if (!snpid_chr) {
              if (fseeko(bgenfile, snpid_slen, SEEK_CUR)) {
                goto ox_bgen_to_pgen_ret_READ_FAIL;
              }
            } else {
              if (!fread_unlocked(loadbuf, snpid_slen, 1, bgenfile)) {
                goto ox_bgen_to_pgen_ret_READ_FAIL;
              }
              loadbuf[snpid_slen] = '\0';
            }
            uint16_t rsid_slen;
            if (!fread_unlocked(&rsid_slen, 2, 1, bgenfile)) {
              goto ox_bgen_to_pgen_ret_READ_FAIL;
            }
            if (fseeko(bgenfile, rsid_slen, SEEK_CUR)) {
              goto ox_bgen_to_pgen_ret_READ_FAIL;
            }
            uint16_t chr_name_slen;
            if (!fread_unlocked(&chr_name_slen, 2, 1, bgenfile)) {
              goto ox_bgen_to_pgen_ret_READ_FAIL;
            }
            if (!snpid_chr) {
              if (!fread_unlocked(loadbuf, chr_name_slen, 1, bgenfile)) {
                goto ox_bgen_to_pgen_ret_READ_FAIL;
              }
              if (strequal_k((char*)loadbuf, "NA", chr_name_slen)) {
                strcpy((char*)loadbuf, "0");
                chr_name_slen = 1;
              } else {
                loadbuf[chr_name_slen] = '\0';
              }
            } else {
              if (fseeko(bgenfile, chr_name_slen, SEEK_CUR)) {
                goto ox_bgen_to_pgen_ret_READ_FAIL;
              }
              chr_name_slen = snpid_slen;
            }
            if (chr_filter_present) {
              const int32_t cur_chr_code = get_chr_code((char*)loadbuf, cip, chr_name_slen);
              assert(cur_chr_code >= 0); // we scanned all the variants
              skip = !is_set(cip->chr_mask, cur_chr_code);
            }

            uint32_t cur_bp; // ignore in this pass
            if (!fread_unlocked(&cur_bp, 4, 1, bgenfile)) {
              goto ox_bgen_to_pgen_ret_READ_FAIL;
            }

            cur_allele_ct = 0;
            if (!fread_unlocked(&cur_allele_ct, 2, 1, bgenfile)) {
              goto ox_bgen_to_pgen_ret_READ_FAIL;
            }
            for (uint32_t allele_idx = 0; allele_idx < cur_allele_ct; ++allele_idx) {
              uint32_t allele_slen;
              if (!fread_unlocked(&allele_slen, 4, 1, bgenfile)) {
                goto ox_bgen_to_pgen_ret_READ_FAIL;
              }
              if (fseeko(bgenfile, allele_slen, SEEK_CUR)) {
                goto ox_bgen_to_pgen_ret_READ_FAIL;
              }
            }
            if (!fread_unlocked(&genodata_byte_ct, 4, 1, bgenfile)) {
              goto ox_bgen_to_pgen_ret_READ_FAIL;
            }

            // "cur_allele_ct > 2" is temporary kludge
            if (skip || (cur_allele_ct > 2)) {
              if (fseeko(bgenfile, genodata_byte_ct, SEEK_CUR)) {
                goto ox_bgen_to_pgen_ret_READ_FAIL;
              }
              goto ox_bgen_to_pgen_load13_start;
            }
            if (uncompressed_genodata_byte_cts) {
              if (!fread_unlocked(&uncompressed_genodata_byte_ct, 4, 1, bgenfile)) {
                goto ox_bgen_to_pgen_ret_READ_FAIL;
              }
              genodata_byte_ct -= 4;
            }

            if ((block_vidx == cur_thread_block_vidx_limit) || ((uintptr_t)(cur_geno_buf_end - bgen_geno_iter) < genodata_byte_ct)) {
              thread_bidxs[++cur_thread_fill_idx] = block_vidx;
              if (cur_thread_fill_idx == calc_thread_ct) {
                prev_allele_ct = cur_allele_ct;
                prev_genodata_byte_ct = genodata_byte_ct;
                break;
              }
              cur_geno_buf_end = &(bgen_geno_iter[thread_wkspace_size]);
              cur_thread_block_vidx_limit = MINV(cur_thread_block_vidx_limit + per_thread_block_limit, block_vidx_limit);
            }
          }
          cur_block_write_ct = block_vidx;
        }
        if (vidx_start) {
          join_threads3z(&ts);
          reterr = g_error_ret;
          if (reterr) {
            goto ox_bgen_to_pgen_ret_bgen13_thread_fail;
          }
        }
        if (!ts.is_last_block) {
          ts.is_last_block = (vidx_start + cur_block_write_ct == variant_ct);
          ts.thread_func_ptr = bgen13_geno_to_pgen_thread;
          if (spawn_threads3z(vidx_start, &ts)) {
            goto ox_bgen_to_pgen_ret_THREAD_CREATE_FAIL;
          }
        }
        parity = 1 - parity;
        if (vidx_start) {
          // write *previous* block results
          const uintptr_t* write_genovec_iter = g_write_genovecs[parity];
          // const uintptr_t* write_phasepresents = g_write_phasepresents[parity];
          // const uintptr_t* write_phaseinfos = g_write_phaseinfos[parity];
          const uintptr_t* write_dosage_presents = g_write_dosage_presents[parity];
          // const uintptr_t* write_dphase_presents = g_write_dphase_presents[parity];
          const uint32_t* write_dosage_cts = g_write_dosage_cts[parity];
          const uint32_t* write_dphase_cts = g_write_dphase_cts[parity];
          const dosage_t* write_dosage_val_bufs = g_write_dosage_val_bufs[parity];
          for (uintptr_t block_vidx = 0; block_vidx < prev_block_write_ct; ++block_vidx) {
            const uint32_t cur_dosage_ct = write_dosage_cts[block_vidx];
            const uint32_t cur_dphase_ct = write_dphase_cts[block_vidx];
            if (!cur_dphase_ct) {
              if (!cur_dosage_ct) {
                if (spgw_append_biallelic_genovec(write_genovec_iter, &spgw)) {
                  goto ox_bgen_to_pgen_ret_WRITE_FAIL;
                }
              } else {
                if (spgw_append_biallelic_genovec_dosage16(write_genovec_iter, &(write_dosage_presents[block_vidx * sample_ctaw]), &(write_dosage_val_bufs[block_vidx * 2 * sample_ct]), cur_dosage_ct, &spgw)) {
                  goto ox_bgen_to_pgen_ret_WRITE_FAIL;
                }
              }
            } else {
              // todo
            }
            write_genovec_iter = &(write_genovec_iter[sample_ctaw2]);
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
    if (fclose_flush_null(writebuf_flush, write_iter, &pvarfile)) {
      goto ox_bgen_to_pgen_ret_WRITE_FAIL;
    }

    spgw_finish(&spgw);
    putc_unlocked('\r', stdout);
    write_iter = strcpya(g_logbuf, "--bgen: ");
    const uint32_t outname_base_slen = (uintptr_t)(outname_end - outname);
    write_iter = memcpya(write_iter, outname, outname_base_slen + 5);
    write_iter = memcpyl3a(write_iter, " + ");
    write_iter = memcpya(write_iter, outname, outname_base_slen);
    write_iter = strcpya(write_iter, ".pvar");
    write_iter = strcpya(write_iter, " written");
    if (!dosage_is_present) {
      write_iter = strcpya(write_iter, " (only hardcalls)");
    }
    strcpy(write_iter, ".\n");
    wordwrapb(0);
    logprintb();
  }
  while (0) {
  ox_bgen_to_pgen_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ox_bgen_to_pgen_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  ox_bgen_to_pgen_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  ox_bgen_to_pgen_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  ox_bgen_to_pgen_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  ox_bgen_to_pgen_ret_INCONSISTENT_INPUT_DELIM:
    logerrprintb();
    if (id_delim == '_') {
      logerrprint("If you do not want '_' to be treated as a FID/IID delimiter, use --double-id or\n--const-fid to choose a different method of converting .bgen sample IDs to\nPLINK IDs, or --id-delim to change the FID/IID delimiter.\n");
    }
    reterr = kPglRetInconsistentInput;
    break;
  ox_bgen_to_pgen_ret_MALFORMED_INPUT_LONG_ID:
    logerrprint("Error: FIDs and IIDs are limited to " MAX_ID_SLEN_STR " characters.\n");
    reterr = kPglRetMalformedInput;
    break;
  ox_bgen_to_pgen_ret_INCONSISTENT_INPUT_2:
    logerrprintb();
  ox_bgen_to_pgen_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  ox_bgen_to_pgen_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  ox_bgen_to_pgen_ret_bgen13_thread_fail:
    if (reterr == kPglRetMalformedInput) {
      logprint("\n");
      logerrprint("Error: Invalid compressed SNP block in .bgen file.\n");
    } else if (reterr == kPglRetNotYetSupported) {
      logprint("\n");
      logerrprint("Error: BGEN import doesn't currently support phased variants, >16-bit\nprobability precision, or ploidy > 2.\n");
    }
  }
 ox_bgen_to_pgen_ret_1:
  if (spgw_cleanup(&spgw) && (!reterr)) {
    reterr = kPglRetWriteFail;
  }
  threads3z_cleanup(&ts, &g_cur_block_write_ct);
  fclose_cond(bgenfile);
  fclose_cond(psamfile);
  fclose_cond(pvarfile);
  bigstack_double_reset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

boolerr_t import_legend_cols(const char* fname, uintptr_t line_idx, uint32_t prov_ref_allele_second, char** loadbuf_iter_ptr, char** write_iter_ptr, uint32_t* variant_ct_ptr) {
  {
    if (*variant_ct_ptr == 0x7ffffffd) {
      logerrprint("Error: " PROG_NAME_STR " does not support more than 2^31 - 3 variants.  We recommend other\nsoftware, such as PLINK/SEQ, for very deep studies of small numbers of genomes.\n");
      return 1;
    }
    *variant_ct_ptr += 1;
    char* write_iter = *write_iter_ptr;
    *write_iter++ = '\t';
    char* id_start = *loadbuf_iter_ptr;
    char* id_end = token_endnn(id_start);
    const uint32_t id_slen = (uintptr_t)(id_end - id_start);
    if (id_slen > kMaxIdSlen) {
      logerrprint("Error: Variant names are limited to " MAX_ID_SLEN_STR " characters.\n");
      return 1;
    }
    char* pos_str = skip_initial_spaces(id_end);
    if (!pos_str) {
      goto import_legend_cols_ret_MISSING_TOKENS;
    }
    char* pos_end = token_endnn(pos_str);
    uint32_t cur_bp;
    if (scan_uint_defcap(pos_str, &cur_bp)) {
      LOGPREPRINTFWW("Error: Invalid bp coordinate on line %" PRIuPTR " of %s.\n", line_idx, fname);
      return 1;
    }
    write_iter = uint32toa_x(cur_bp, '\t', write_iter);
    write_iter = memcpyax(write_iter, id_start, id_slen, '\t');
    char* first_allele_str = skip_initial_spaces(pos_end);
    if (is_eoln_kns(*first_allele_str)) {
      goto import_legend_cols_ret_MISSING_TOKENS;
    }
    char* first_allele_end = token_endnn(first_allele_str);
    char* second_allele_str = skip_initial_spaces(first_allele_end);
    if (is_eoln_kns(*second_allele_str)) {
      goto import_legend_cols_ret_MISSING_TOKENS;
    }
    char* second_allele_end = token_endnn(second_allele_str);
    if (!prov_ref_allele_second) {
      write_iter = memcpyax(write_iter, first_allele_str, (uintptr_t)(first_allele_end - first_allele_str), '\t');
      write_iter = memcpya(write_iter, second_allele_str, (uintptr_t)(second_allele_end - second_allele_str));
    } else {
      write_iter = memcpyax(write_iter, second_allele_str, (uintptr_t)(second_allele_end - second_allele_str), '\t');
      write_iter = memcpya(write_iter, first_allele_str, (uintptr_t)(first_allele_end - first_allele_str));
    }
    *write_iter_ptr = write_iter;
    append_binary_eoln(write_iter_ptr);
    *loadbuf_iter_ptr = second_allele_end;
    return 0;
  }
  {
  import_legend_cols_ret_MISSING_TOKENS:
    LOGPREPRINTFWW("Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, fname);
    return 1;
  }
}

pglerr_t scan_haps_for_het(char* loadbuf_iter, const char* hapsname, uint32_t sample_ct, uint32_t is_haploid_or_mt, uintptr_t line_idx_haps, uint32_t* at_least_one_het_ptr) {
  pglerr_t reterr = kPglRetSuccess;
  {
    for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
      const uint32_t first_hap_char_code = (uint32_t)((unsigned char)(*loadbuf_iter));
      const uint32_t first_hap_int = first_hap_char_code - 48;
      // will .haps files ever support triallelic variants?  don't worry about
      // that for now
      char* post_first_hap = &(loadbuf_iter[1]);
      if ((first_hap_int >= 2) || ((unsigned char)(*post_first_hap) > 32)) {
        if (first_hap_char_code <= 32) {
          goto scan_haps_for_het_ret_MISSING_TOKENS;
        }
        goto scan_haps_for_het_ret_INVALID_TOKEN;
      }
      char* second_hap = skip_initial_spaces(post_first_hap);
      char* post_second_hap = &(second_hap[1]);
      const uint32_t second_hap_char_code = (uint32_t)((unsigned char)(*second_hap));
      const uint32_t second_hap_int = second_hap_char_code - 48;
      const uint32_t post_second_hap_char_code = (uint32_t)((unsigned char)(*post_second_hap));
      if ((second_hap_int >= 2) || (post_second_hap_char_code > 32)) {
        // if haploid or MT, permit '-' in second column
        if ((!is_haploid_or_mt) || (second_hap_char_code != 45)) {
          if (second_hap_char_code <= 32) {
            goto scan_haps_for_het_ret_MISSING_TOKENS;
          }
          if ((second_hap_char_code == 45) && (post_second_hap_char_code <= 32)) {
            goto scan_haps_for_het_ret_HAPLOID_TOKEN;
          }
          goto scan_haps_for_het_ret_INVALID_TOKEN;
        }
      } else if (first_hap_int != second_hap_int) {
        *at_least_one_het_ptr = 1;
        break;
      }
      loadbuf_iter = skip_initial_spaces(post_second_hap);
    }
  }
  while (0) {
  scan_haps_for_het_ret_HAPLOID_TOKEN:
    sprintf(g_logbuf, "Error: Haploid/MT-only token on line %" PRIuPTR " of %s.\n", line_idx_haps, hapsname);
    reterr = kPglRetInconsistentInput;
    break;
  scan_haps_for_het_ret_INVALID_TOKEN:
    sprintf(g_logbuf, "Error: Invalid token on line %" PRIuPTR " of %s.\n", line_idx_haps, hapsname);
    reterr = kPglRetMalformedInput;
    break;
  scan_haps_for_het_ret_MISSING_TOKENS:
    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx_haps, hapsname);
    reterr = kPglRetMalformedInput;
    break;
  }
  return reterr;
}

#ifdef __arm__
  #error "Unaligned accesses in ox_hapslegend_to_pgen()."
#endif
pglerr_t ox_hapslegend_to_pgen(const char* hapsname, const char* legendname, const char* samplename, const char* ox_single_chr_str, const char* ox_missing_code, misc_flags_t misc_flags, oxford_import_t oxford_import_flags, char* outname, char* outname_end, chr_info_t* cip) {
  unsigned char* bigstack_mark = g_bigstack_base;
  gzFile gz_hapsfile = nullptr;
  gzFile gz_legendfile = nullptr;
  FILE* outfile = nullptr;
  uintptr_t line_idx_haps = 0;
  uintptr_t line_idx_legend = 0;
  pglerr_t reterr = kPglRetSuccess;
  st_pgen_writer_t spgw;
  uintptr_t loadbuf_size;
  spgw_preinit(&spgw);
  {
    uint32_t sfile_sample_ct = 0;
    if (samplename[0]) {
      reterr = ox_sample_to_psam(samplename, ox_missing_code, misc_flags, outname, outname_end, &sfile_sample_ct);
      if (reterr) {
        goto ox_hapslegend_to_pgen_ret_1;
      }
      if (sfile_sample_ct > (kMaxLongLine / 4)) {
        logerrprint("Error: Too many samples for .haps file converter.\n");
        reterr = kPglRetNotYetSupported;
        goto ox_hapslegend_to_pgen_ret_1;
      }
    }

    reterr = gzopen_read_checked(hapsname, &gz_hapsfile);
    if (reterr) {
      goto ox_hapslegend_to_pgen_ret_1;
    }
    uintptr_t writebuf_size = bigstack_left() / 2;
    if (writebuf_size < 2 * kMaxMediumLine + kCacheline) {
      return kPglRetNomem;
#ifdef __LP64__
      // in 32-bit case, kMaxLongLine + kMaxMediumLine overflows
    } else if (writebuf_size > kMaxLongLine + kMaxMediumLine) {
      writebuf_size = kMaxLongLine + kMaxMediumLine;
#endif
    } else {
      writebuf_size &= ~(kCacheline - 1);
    }
    loadbuf_size = writebuf_size - kMaxMediumLine;
    char* loadbuf = (char*)bigstack_alloc_raw(loadbuf_size);
    loadbuf[loadbuf_size - 1] = ' ';
    char* writebuf = (char*)bigstack_alloc_raw(writebuf_size);
    char* writebuf_flush = &(writebuf[kMaxMediumLine]);
    strcpy(outname_end, ".pvar");
    if (fopen_checked(outname, FOPEN_WB, &outfile)) {
      goto ox_hapslegend_to_pgen_ret_OPEN_FAIL;
    }
    char* write_iter = strcpya(writebuf, "#CHROM\tPOS\tID\tREF\tALT" EOLN_STR);
    char* loadbuf_first_token;
    do {
      ++line_idx_haps;
      if (!gzgets(gz_hapsfile, loadbuf, loadbuf_size)) {
        if (!gzeof(gz_hapsfile)) {
          goto ox_hapslegend_to_pgen_ret_READ_FAIL;
        }
        sprintf(g_logbuf, "Error: %s is empty.\n", hapsname);
        goto ox_hapslegend_to_pgen_ret_INCONSISTENT_INPUT_WW;
      }
      if (!loadbuf[loadbuf_size - 1]) {
        goto ox_hapslegend_to_pgen_ret_LONG_LINE_HAP;
      }
      loadbuf_first_token = skip_initial_spaces(loadbuf);
    } while (is_eoln_kns(*loadbuf_first_token));
    const uint32_t token_ct = count_tokens(loadbuf_first_token);
    // pass 1: count variants, write .pvar file, may as well also verify
    // there's at least one heterozygous call
    finalize_chrset(misc_flags, cip);
    const uint32_t allow_extra_chrs = (misc_flags / kfMiscAllowExtraChrs) & 1;
    const int32_t mt_code = cip->xymt_codes[kChrOffsetMT];
    const uint32_t prov_ref_allele_second = !(oxford_import_flags & kfOxfordImportRefFirst);
    uint32_t at_least_one_het = 0;
    uintptr_t variant_skip_ct = 0;
    uint32_t variant_ct = 0;
    uint32_t is_haploid_or_mt = 0;
    uint32_t sample_ct;
    // support both .haps + .legend (.haps expected to contain no header
    // columns), and pure .haps
    if (legendname[0]) {
      assert(ox_single_chr_str);
      if (token_ct % 2) {
        sprintf(g_logbuf, "Error: %s has an odd number of tokens in the first line. (With --haps + --legend, the .haps file is expected to have no header columns.)\n", hapsname);
        goto ox_hapslegend_to_pgen_ret_MALFORMED_INPUT_WW;
      }
      sample_ct = token_ct / 2;
      if (sfile_sample_ct && (sfile_sample_ct != sample_ct)) {
        sprintf(g_logbuf, "Error: .sample file has %u sample%s, while %s has %u.\n", sfile_sample_ct, (sfile_sample_ct == 1)? "" : "s", hapsname, sample_ct);
        goto ox_hapslegend_to_pgen_ret_INCONSISTENT_INPUT_WW;
      }
      if (gzrewind(gz_hapsfile)) {
        goto ox_hapslegend_to_pgen_ret_READ_FAIL;
      }
      line_idx_haps = 0;
      int32_t chr_code_raw = get_chr_code_raw(ox_single_chr_str);
      const char* single_chr_str = nullptr;
      uint32_t single_chr_slen;
      char chr_buf[8]; // nothing longer than e.g. "chrMT" for now
      if (chr_code_raw == -1) {
        // command-line parser guarantees that allow_extra_chrs is true here
        single_chr_str = ox_single_chr_str;
        single_chr_slen = strlen(ox_single_chr_str);
      } else {
        uint32_t chr_code = chr_code_raw;
        if (chr_code > cip->max_code) {
          if (chr_code < kMaxContigs) {
            logerrprint("Error: --legend chromosome code is not in the chromosome set.\n");
            goto ox_hapslegend_to_pgen_ret_INVALID_CMDLINE;
          }
          chr_code = cip->xymt_codes[chr_code - kMaxContigs];
          if (((int32_t)chr_code) < 0) {
            logerrprint("Error: --legend chromosome code is not in the chromosome set.\n");
            goto ox_hapslegend_to_pgen_ret_INVALID_CMDLINE;
          }
        }
        if (!is_set(cip->chr_mask, chr_code)) {
          logerrprint("Error: --legend chromosome code is excluded by chromosome filter.\n");
          goto ox_hapslegend_to_pgen_ret_INVALID_CMDLINE;
        }
        is_haploid_or_mt = is_set(cip->haploid_mask, chr_code) || (((int32_t)chr_code) == mt_code);
        char* chr_name_end = chr_name_write(cip, chr_code, chr_buf);
        single_chr_str = chr_buf;
        single_chr_slen = (uintptr_t)(chr_name_end - chr_buf);
      }
      reterr = gzopen_read_checked(legendname, &gz_legendfile);
      if (reterr) {
        goto ox_hapslegend_to_pgen_ret_1;
      }
      do {
        ++line_idx_legend;
        if (!gzgets(gz_legendfile, loadbuf, loadbuf_size)) {
          if (!gzeof(gz_legendfile)) {
            goto ox_hapslegend_to_pgen_ret_READ_FAIL;
          }
          sprintf(g_logbuf, "Error: %s is empty.\n", legendname);
          goto ox_hapslegend_to_pgen_ret_MALFORMED_INPUT_WW;
        }
        if (!loadbuf[loadbuf_size - 1]) {
          goto ox_hapslegend_to_pgen_ret_LONG_LINE_LEGEND;
        }
        loadbuf_first_token = skip_initial_spaces(loadbuf);
      } while (is_eoln_kns(*loadbuf_first_token));
      // require at least 4 columns, in ID/pos/A1/A2 order; header text is
      // permitted to vary.  tolerate and ignore extra columns.
      if (!next_token_mult(loadbuf_first_token, 3)) {
        goto ox_hapslegend_to_pgen_ret_MISSING_TOKENS_LEGEND;
      }
      while (1) {
        ++line_idx_legend;
        if (!gzgets(gz_legendfile, loadbuf, loadbuf_size)) {
          if (!gzeof(gz_legendfile)) {
            goto ox_hapslegend_to_pgen_ret_READ_FAIL;
          }
          break;
        }
        if (!loadbuf[loadbuf_size - 1]) {
          goto ox_hapslegend_to_pgen_ret_LONG_LINE_LEGEND;
        }
        loadbuf_first_token = skip_initial_spaces(loadbuf);
        if (is_eoln_kns(*loadbuf_first_token)) {
          continue;
        }
        char* loadbuf_iter = loadbuf_first_token;
        write_iter = memcpya(write_iter, single_chr_str, single_chr_slen);
        if (import_legend_cols(legendname, line_idx_legend, prov_ref_allele_second, &loadbuf_iter, &write_iter, &variant_ct)) {
          putc_unlocked('\n', stdout);
          goto ox_hapslegend_to_pgen_ret_MALFORMED_INPUT;
        }
        if (fwrite_ck(writebuf_flush, outfile, &write_iter)) {
          goto ox_hapslegend_to_pgen_ret_WRITE_FAIL;
        }
        if (!at_least_one_het) {
          do {
            ++line_idx_haps;
            if (!gzgets(gz_hapsfile, loadbuf, loadbuf_size)) {
              if (!gzeof(gz_hapsfile)) {
                goto ox_hapslegend_to_pgen_ret_READ_FAIL;
              }
              sprintf(g_logbuf, "Error: %s has fewer nonheader lines than %s.\n", hapsname, legendname);
              goto ox_hapslegend_to_pgen_ret_INCONSISTENT_INPUT_WW;
            }
            if (!loadbuf[loadbuf_size - 1]) {
              goto ox_hapslegend_to_pgen_ret_LONG_LINE_HAP;
            }
            loadbuf_first_token = skip_initial_spaces(loadbuf);
          } while (is_eoln_kns(*loadbuf_first_token));
          reterr = scan_haps_for_het(loadbuf_first_token, hapsname, sample_ct, is_haploid_or_mt, line_idx_haps, &at_least_one_het);
          if (reterr) {
            putc_unlocked('\n', stdout);
            wordwrapb(0);
            logerrprintb();
            goto ox_hapslegend_to_pgen_ret_1;
          }
        }
      }
      if (gzrewind(gz_legendfile)) {
        goto ox_hapslegend_to_pgen_ret_READ_FAIL;
      }
    } else {
      assert(!legendname[0]);
      if ((token_ct < 7) || (!(token_ct % 2))) {
        sprintf(g_logbuf, "Error: Unexpected token count in line %" PRIuPTR " of %s (should be odd, >5).\n", line_idx_haps, hapsname);
        goto ox_hapslegend_to_pgen_ret_MALFORMED_INPUT_WW;
      }
      sample_ct = (token_ct - 5) / 2;
      if (sfile_sample_ct && (sfile_sample_ct != sample_ct)) {
        sprintf(g_logbuf, "Error: .sample file has %u sample%s, while %s has %u.\n", sfile_sample_ct, (sfile_sample_ct == 1)? "" : "s", hapsname, sample_ct);
        goto ox_hapslegend_to_pgen_ret_INCONSISTENT_INPUT_WW;
      }
      while (1) {
        if (!is_eoln_kns(*loadbuf_first_token)) {
          char* chr_code_end = token_endnn(loadbuf_first_token);
          char* loadbuf_iter = skip_initial_spaces(chr_code_end);
          if (is_eoln_kns(*loadbuf_iter)) {
            goto ox_hapslegend_to_pgen_ret_MISSING_TOKENS_HAPS;
          }
          int32_t cur_chr_code;
          reterr = get_or_add_chr_code_destructive("--haps file", line_idx_haps, allow_extra_chrs, loadbuf_first_token, chr_code_end, cip, &cur_chr_code);
          if (reterr) {
            goto ox_hapslegend_to_pgen_ret_1;
          }
          if (!is_set(cip->chr_mask, cur_chr_code)) {
            ++variant_skip_ct;
          } else {
            is_haploid_or_mt = is_set(cip->haploid_mask, cur_chr_code) || (cur_chr_code == mt_code);
            write_iter = chr_name_write(cip, cur_chr_code, write_iter);
            if (import_legend_cols(hapsname, line_idx_haps, prov_ref_allele_second, &loadbuf_iter, &write_iter, &variant_ct)) {
              putc_unlocked('\n', stdout);
              goto ox_hapslegend_to_pgen_ret_MALFORMED_INPUT;
            }
            if (fwrite_ck(writebuf_flush, outfile, &write_iter)) {
              goto ox_hapslegend_to_pgen_ret_WRITE_FAIL;
            }
            if (!at_least_one_het) {
              loadbuf_iter = skip_initial_spaces(loadbuf_iter);
              if (scan_haps_for_het(loadbuf_iter, hapsname, sample_ct, is_haploid_or_mt, line_idx_haps, &at_least_one_het)) {
                // override InconsistentInput return code since chromosome info
                // was also gathered from .haps file
                goto ox_hapslegend_to_pgen_ret_MALFORMED_INPUT_WW;
              }
            }
          }
        }
        ++line_idx_haps;
        if (!gzgets(gz_hapsfile, loadbuf, loadbuf_size)) {
          if (!gzeof(gz_hapsfile)) {
            goto ox_hapslegend_to_pgen_ret_READ_FAIL;
          }
          break;
        }
        if (!loadbuf[loadbuf_size - 1]) {
          goto ox_hapslegend_to_pgen_ret_LONG_LINE_HAP;
        }
        loadbuf_first_token = skip_initial_spaces(loadbuf);
      }
      if (!variant_ct) {
        sprintf(g_logbuf, "Error: All %" PRIuPTR " variant%s in %s skipped due to chromosome filter.\n", variant_skip_ct, (variant_skip_ct == 1)? "" : "s", hapsname);
        goto ox_hapslegend_to_pgen_ret_INCONSISTENT_INPUT_WW;
      }
    }
    if (fclose_flush_null(writebuf_flush, write_iter, &outfile)) {
      goto ox_hapslegend_to_pgen_ret_WRITE_FAIL;
    }
    if (!sfile_sample_ct) {
      // create a dummy .psam file with "per0/per0", "per1/per1", etc. IDs,
      // matching --dummy
      strcpy(outname_end, ".psam");
      if (fopen_checked(outname, FOPEN_WB, &outfile)) {
        goto ox_hapslegend_to_pgen_ret_OPEN_FAIL;
      }
      write_iter = strcpya(writebuf, "#FID\tIID\tSEX" EOLN_STR);
      for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
        write_iter = memcpyl3a(write_iter, "per");
        write_iter = uint32toa(sample_idx, write_iter);
        write_iter = strcpya(write_iter, "\tper");
        write_iter = uint32toa(sample_idx, write_iter);
        write_iter = strcpya(write_iter, "\tNA" EOLN_STR);
        if (fwrite_ck(writebuf_flush, outfile, &write_iter)) {
          goto ox_hapslegend_to_pgen_ret_WRITE_FAIL;
        }
      }
      if (fclose_flush_null(writebuf_flush, write_iter, &outfile)) {
        goto ox_hapslegend_to_pgen_ret_WRITE_FAIL;
      }
    }
    if (gzrewind(gz_hapsfile)) {
      goto ox_hapslegend_to_pgen_ret_READ_FAIL;
    }
    line_idx_haps = 0;
    bigstack_reset(writebuf);
    putc_unlocked('\r', stdout);
    LOGPRINTF("--haps%s: %u variant%s scanned.\n", legendname[0]? " + --legend" : "", variant_ct, (variant_ct == 1)? "" : "s");
    strcpy(outname_end, ".pgen");
    uintptr_t spgw_alloc_cacheline_ct;
    uint32_t max_vrec_len;
    reterr = spgw_init_phase1(outname, nullptr, nullptr, variant_ct, sample_ct, at_least_one_het? kfPgenGlobalHardcallPhasePresent : kfPgenGlobal0, (oxford_import_flags & (kfOxfordImportRefFirst | kfOxfordImportRefSecond))? 1 : 2, &spgw, &spgw_alloc_cacheline_ct, &max_vrec_len);
    if (reterr) {
      goto ox_hapslegend_to_pgen_ret_1;
    }
    unsigned char* spgw_alloc;
    if (bigstack_alloc_uc(spgw_alloc_cacheline_ct * kCacheline, &spgw_alloc)) {
      goto ox_hapslegend_to_pgen_ret_NOMEM;
    }
    spgw_init_phase2(max_vrec_len, &spgw, spgw_alloc);
    const uint32_t sample_ctl2 = QUATERCT_TO_WORDCT(sample_ct);
    const uint32_t sample_ctl2_m1 = sample_ctl2 - 1;
    const uint32_t sample_ctl = BITCT_TO_WORDCT(sample_ct);
    const uint32_t phaseinfo_match_4char = prov_ref_allele_second? 0x20312030 : 0x20302031;
    const uint32_t phaseinfo_match = 1 + prov_ref_allele_second;
    uintptr_t* genovec;
    uintptr_t* phaseinfo;
    if (bigstack_alloc_ul(sample_ctl2, &genovec) ||
        bigstack_alloc_ul(sample_ctl, &phaseinfo)) {
      goto ox_hapslegend_to_pgen_ret_NOMEM;
    }

    for (uint32_t vidx = 0; vidx < variant_ct;) {
      ++line_idx_haps;
      if (!gzgets(gz_hapsfile, loadbuf, loadbuf_size)) {
        if ((!gzeof(gz_hapsfile)) || (!legendname[0])) {
          goto ox_hapslegend_to_pgen_ret_READ_FAIL;
        }
        sprintf(g_logbuf, "Error: %s has fewer nonheader lines than %s.\n", hapsname, legendname);
        goto ox_hapslegend_to_pgen_ret_INCONSISTENT_INPUT_WW;
      }
      if (!loadbuf[loadbuf_size - 1]) {
        goto ox_hapslegend_to_pgen_ret_LONG_LINE_HAP;
      }
      loadbuf_first_token = skip_initial_spaces(loadbuf);
      if (is_eoln_kns(*loadbuf_first_token)) {
        continue;
      }
      char* loadbuf_iter = loadbuf_first_token;
      if (!legendname[0]) {
        char* chr_code_end = token_endnn(loadbuf_first_token);
        const int32_t cur_chr_code = get_chr_code_counted(cip, (uintptr_t)(chr_code_end - loadbuf_first_token), loadbuf_first_token);
        if (!is_set(cip->chr_mask, cur_chr_code)) {
          continue;
        }
        is_haploid_or_mt = is_set(cip->haploid_mask, cur_chr_code) || (cur_chr_code == mt_code);
        loadbuf_iter = next_token_mult(skip_initial_spaces(chr_code_end), 4);
        if (!loadbuf_iter) {
          goto ox_hapslegend_to_pgen_ret_MISSING_TOKENS_HAPS;
        }
      }
      uintptr_t genovec_word_or = 0;
      uint32_t inner_loop_last = kBitsPerWordD2 - 1;
      uint32_t widx = 0;
      // optimize common case: autosomal diploid, always exactly one space
      // this loop is time-critical; all my attemps to merge in the haploid/MT
      // case have caused >10% slowdowns
      if ((!is_haploid_or_mt) && ((unsigned char)loadbuf_iter[sample_ct * 4 - 1] < 32)) {
        loadbuf_iter[sample_ct * 4 - 1] = ' ';
        const uint32_t* loadbuf_alias32_iter = (const uint32_t*)loadbuf_iter;
        while (1) {
          if (widx >= sample_ctl2_m1) {
            if (widx > sample_ctl2_m1) {
              break;
            }
            inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
          }
          uintptr_t genovec_word = 0;
          uint32_t phaseinfo_halfword = 0;
          for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
            // assumes little-endian
            uint32_t cur_hap_4char = *loadbuf_alias32_iter++;
            if ((cur_hap_4char & 0xfffefffeU) != 0x20302030) {
              if ((cur_hap_4char & 0xfffffffeU) == 0x202d2030) {
                // "0 - ", "1 - "
                goto ox_hapslegend_to_pgen_ret_HAPLOID_TOKEN;
              }
              // any character < 32?
              if ((((cur_hap_4char & 0xe0e0e0e0U) * 7) & 0x80808080U) != 0x80808080U) {
                goto ox_hapslegend_to_pgen_ret_MISSING_TOKENS_HAPS;
              }
              sprintf(g_logbuf, "Error: Invalid token on line %" PRIuPTR " of %s.\n", line_idx_haps, hapsname);
              goto ox_hapslegend_to_pgen_ret_MALFORMED_INPUT_WW;
            }
            const uintptr_t new_geno = (cur_hap_4char + (cur_hap_4char >> 16)) & 3;
            genovec_word |= new_geno << (2 * sample_idx_lowbits);
            if (cur_hap_4char == phaseinfo_match_4char) {
              phaseinfo_halfword |= 1U << sample_idx_lowbits;
            }
          }
          genovec[widx] = genovec_word;
          genovec_word_or |= genovec_word;
          ((halfword_t*)phaseinfo)[widx] = phaseinfo_halfword;
          ++widx;
        }
      } else {
        while (1) {
          if (widx >= sample_ctl2_m1) {
            if (widx > sample_ctl2_m1) {
              break;
            }
            inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
          }
          uintptr_t genovec_word = 0;
          uint32_t phaseinfo_halfword = 0;
          for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
            const uint32_t first_hap_char_code = (uint32_t)((unsigned char)(*loadbuf_iter));
            const uint32_t first_hap_int = first_hap_char_code - 48;
            char* post_first_hap = &(loadbuf_iter[1]);
            if ((first_hap_int >= 2) || ((unsigned char)(*post_first_hap) > 32)) {
              if (first_hap_char_code <= 32) {
                goto ox_hapslegend_to_pgen_ret_MISSING_TOKENS_HAPS;
              }
              sprintf(g_logbuf, "Error: Invalid token on line %" PRIuPTR " of %s.\n", line_idx_haps, hapsname);
              goto ox_hapslegend_to_pgen_ret_MALFORMED_INPUT_WW;
            }
            char* second_hap = skip_initial_spaces(post_first_hap);
            char* post_second_hap = &(second_hap[1]);
            const uint32_t second_hap_char_code = (uint32_t)((unsigned char)(*second_hap));
            const uint32_t post_second_hap_char_code = (uint32_t)((unsigned char)(*post_second_hap));
            uint32_t second_hap_int = second_hap_char_code - 48;
            if ((second_hap_int >= 2) || (post_second_hap_char_code > 32)) {
              if (is_haploid_or_mt && (second_hap_char_code == 45)) {
                // could require --sample, and require this sample to be male
                // in this case?
                second_hap_int = first_hap_int;
              } else {
                if (second_hap_char_code <= 32) {
                  goto ox_hapslegend_to_pgen_ret_MISSING_TOKENS_HAPS;
                }
                if ((second_hap_char_code == 45) && (post_second_hap_char_code <= 32)) {
                  goto ox_hapslegend_to_pgen_ret_HAPLOID_TOKEN;
                }
                sprintf(g_logbuf, "Error: Invalid token on line %" PRIuPTR " of %s.\n", line_idx_haps, hapsname);
                goto ox_hapslegend_to_pgen_ret_MALFORMED_INPUT_WW;
              }
            }
            genovec_word |= ((uintptr_t)(first_hap_int + second_hap_int)) << (2 * sample_idx_lowbits);
            if (first_hap_int + 2 * second_hap_int == phaseinfo_match) {
              phaseinfo_halfword |= 1U << sample_idx_lowbits;
            }
            loadbuf_iter = skip_initial_spaces(post_second_hap);
          }
          genovec[widx] = genovec_word;
          genovec_word_or |= genovec_word;
          ((halfword_t*)phaseinfo)[widx] = phaseinfo_halfword;
          ++widx;
        }
      }
      if (prov_ref_allele_second) {
        genovec_invert_unsafe(sample_ct, genovec);
        zero_trailing_quaters(sample_ct, genovec);
      }
      if (genovec_word_or & kMask5555) {
        if (spgw_append_biallelic_genovec_hphase(genovec, nullptr, phaseinfo, &spgw)) {
          goto ox_hapslegend_to_pgen_ret_WRITE_FAIL;
        }
      } else {
        if (spgw_append_biallelic_genovec(genovec, &spgw)) {
          goto ox_hapslegend_to_pgen_ret_WRITE_FAIL;
        }
      }
      if (!(++vidx % 1000)) {
        printf("\r--haps%s: %uk variants converted.", legendname[0]? " + --legend" : "", vidx / 1000);
        fflush(stdout);
      }
    }
    spgw_finish(&spgw);
    putc_unlocked('\r', stdout);
    write_iter = strcpya(g_logbuf, "--haps");
    if (legendname[0]) {
      write_iter = strcpya(write_iter, " + --legend");
    }
    write_iter = strcpya(write_iter, ": ");
    const uint32_t outname_base_slen = (uintptr_t)(outname_end - outname);
    write_iter = memcpya(write_iter, outname, outname_base_slen + 5);
    write_iter = memcpyl3a(write_iter, " + ");
    if (!sfile_sample_ct) {
      write_iter = memcpya(write_iter, outname, outname_base_slen);
      write_iter = strcpya(write_iter, ".psam + ");
    }
    write_iter = memcpya(write_iter, outname, outname_base_slen);
    strcpy(write_iter, ".pvar written.\n");
    wordwrapb(0);
    logprintb();
  }
  while (0) {
  ox_hapslegend_to_pgen_ret_LONG_LINE_HAP:
    if (loadbuf_size == kMaxLongLine) {
      sprintf(g_logbuf, "Line %" PRIuPTR " of %s is pathologically long.\n", line_idx_haps, hapsname);
      goto ox_hapslegend_to_pgen_ret_MALFORMED_INPUT_WW;
    }
  ox_hapslegend_to_pgen_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ox_hapslegend_to_pgen_ret_LONG_LINE_LEGEND:
    if (loadbuf_size == kMaxLongLine) {
      sprintf(g_logbuf, "Line %" PRIuPTR " of %s is pathologically long.\n", line_idx_legend, legendname);
      goto ox_hapslegend_to_pgen_ret_MALFORMED_INPUT_WW;
    }
    reterr = kPglRetNomem;
    break;
  ox_hapslegend_to_pgen_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  ox_hapslegend_to_pgen_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  ox_hapslegend_to_pgen_ret_WRITE_FAIL:
    reterr = kPglRetOpenFail;
    break;
  ox_hapslegend_to_pgen_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  ox_hapslegend_to_pgen_ret_MISSING_TOKENS_HAPS:
    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx_haps, hapsname);
  ox_hapslegend_to_pgen_ret_MALFORMED_INPUT_WW:
    putc_unlocked('\n', stdout);
    wordwrapb(0);
    logerrprintb();
  ox_hapslegend_to_pgen_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  ox_hapslegend_to_pgen_ret_MISSING_TOKENS_LEGEND:
    putc_unlocked('\n', stdout);
    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx_legend, legendname);
    wordwrapb(0);
    logerrprintb();
    reterr = kPglRetMalformedInput;
    break;
  ox_hapslegend_to_pgen_ret_HAPLOID_TOKEN:
    putc_unlocked('\n', stdout);
    sprintf(g_logbuf, "Error: Haploid/MT-only token on line %" PRIuPTR " of %s.\n", line_idx_haps, hapsname);
    wordwrapb(0);
    logerrprintb();
    reterr = legendname[0]? kPglRetInconsistentInput : kPglRetMalformedInput;
    break;
  ox_hapslegend_to_pgen_ret_INCONSISTENT_INPUT_WW:
    putc_unlocked('\n', stdout);
    wordwrapb(0);
    logerrprintb();
    reterr = kPglRetInconsistentInput;
    break;
  }
 ox_hapslegend_to_pgen_ret_1:
  spgw_cleanup(&spgw);
  fclose_cond(outfile);
  gzclose_cond(gz_legendfile);
  gzclose_cond(gz_hapsfile);
  bigstack_reset(bigstack_mark);
  return reterr;
}


// could add an option to load_pvar() to not require allele columns, but .map
// is easy enough to write a separate loader for...
CONSTU31(kLoadMapBlockSize, 65536);

// assumes finalize_chrset() has already been called.
static_assert(kMaxContigs <= 65536, "load_map() needs to be updated.");
pglerr_t load_map(const char* mapname, misc_flags_t misc_flags, chr_info_t* cip, uint32_t* max_variant_id_slen_ptr, uint16_t** variant_chr_codes_ptr, uint32_t** variant_bps_ptr, char*** variant_ids_ptr, double** variant_cms_ptr, uint32_t* variant_ct_ptr) {
  // caller should call forget_extra_chr_names(1, cip) after finishing import.
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  gzFile gz_infile = nullptr;
  uintptr_t loadbuf_size = 0;
  uintptr_t line_idx = 0;
  pglerr_t reterr = kPglRetSuccess;
  {
    reterr = gzopen_read_checked(mapname, &gz_infile);
    if (reterr) {
      goto load_map_ret_1;
    }
    // Workspace used as follows:
    // |--loadbuf--|--temp-->----|----<- variant IDs --|
    //            1/4                                 end
    // loadbuf is overwritten with the main return arrays at the end.
    loadbuf_size = round_down_pow2(bigstack_left() / 4, kCacheline);
    char* loadbuf = (char*)g_bigstack_base;
    if (loadbuf_size > kMaxLongLine) {
      loadbuf_size = kMaxLongLine;
    } else if (loadbuf_size <= kMaxMediumLine) {
      goto load_map_ret_NOMEM;
    }
    loadbuf[loadbuf_size - 1] = ' ';
    char* loadbuf_first_token;
    do {
      ++line_idx;
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
        if (!gzeof(gz_infile)) {
          goto load_map_ret_READ_FAIL;
        }
        logerrprint("Error: Empty .map file.\n");
        goto load_map_ret_INCONSISTENT_INPUT;
      }
      if (!loadbuf[loadbuf_size - 1]) {
        goto load_map_ret_LONG_LINE;
      }
      loadbuf_first_token = skip_initial_spaces(loadbuf);
    } while (is_eoln_kns(*loadbuf_first_token) || (*loadbuf_first_token == '#'));
    char* loadbuf_iter = next_token_mult(loadbuf_first_token, 2);
    if (!loadbuf_iter) {
      goto load_map_ret_MISSING_TOKENS;
    }
    uint32_t map_cols = 3;
    loadbuf_iter = next_token(loadbuf_iter);
    if (loadbuf_iter) {
      loadbuf_iter = next_token(loadbuf_iter);
      if (!loadbuf_iter) {
        map_cols = 4;
      } else {
        loadbuf_iter = next_token(loadbuf_iter);
        if (loadbuf_iter) {
          if (next_token(loadbuf_iter)) {
            // do NOT permit >6 columns, .bim is ok but .pvar is not
            // (pointless to support .pvar for legacy formats)
            sprintf(g_logbuf, "Error: %s is not a .map/.bim file (too many columns).\n", mapname);
            goto load_map_ret_MALFORMED_INPUT_WW;
          }
          map_cols = 4;
        }
      }
    }

    const uint32_t allow_extra_chrs = (misc_flags / kfMiscAllowExtraChrs) & 1;
    uint32_t max_variant_id_slen = *max_variant_id_slen_ptr;
    unsigned char* tmp_alloc_base = (unsigned char*)(&(loadbuf[loadbuf_size]));
    unsigned char* tmp_alloc_end = bigstack_end_mark;
    uint16_t* cur_chr_codes = nullptr;
    uint32_t* cur_bps = nullptr;
    char** cur_ids = nullptr;
    double* cur_cms = nullptr;
    double cur_cm = 0.0;
    uint32_t at_least_one_nzero_cm = 0;
    uint32_t variant_ct = 0;
    while (1) {
      if (!is_eoln_kns(*loadbuf_first_token)) {
        // chrom, id, (cm?), pos
        char* loadbuf_iter = token_endnn(loadbuf_first_token);
        if (!(*loadbuf_iter)) {
          goto load_map_ret_MISSING_TOKENS;
        }
        int32_t cur_chr_code;
        reterr = get_or_add_chr_code_destructive(".map file", line_idx, allow_extra_chrs, loadbuf_first_token, loadbuf_iter, cip, &cur_chr_code);
        if (reterr) {
          goto load_map_ret_1;
        }
        if (!is_set(cip->chr_mask, cur_chr_code)) {
          goto load_map_skip_variant;
        }
        loadbuf_iter = skip_initial_spaces(&(loadbuf_iter[1]));
        if (is_eoln_kns(*loadbuf_iter)) {
          goto load_map_ret_MISSING_TOKENS;
        }
        char* token_end = token_endnn(loadbuf_iter);
        uint32_t id_slen = (uintptr_t)(token_end - loadbuf_iter);
        if (id_slen > max_variant_id_slen) {
          max_variant_id_slen = id_slen;
        }
        tmp_alloc_end -= id_slen + 1;
        if (tmp_alloc_end < tmp_alloc_base) {
          goto load_map_ret_NOMEM;
        }
        memcpyx(tmp_alloc_end, loadbuf_iter, id_slen, '\0');
        loadbuf_iter = skip_initial_spaces(token_end);
        if (is_eoln_kns(*loadbuf_iter)) {
          goto load_map_ret_MISSING_TOKENS;
        }

        if (map_cols == 4) {
          char* cm_end = scanadv_double(loadbuf_iter, &cur_cm);
          if (!cm_end) {
            sprintf(g_logbuf, "Error: Invalid centimorgan position on line %" PRIuPTR " of %s.\n", line_idx, mapname);
            goto load_map_ret_MALFORMED_INPUT_WW;
          }
          at_least_one_nzero_cm = (cur_cm != 0.0);
          loadbuf_iter = next_token(cm_end);
          if (!loadbuf_iter) {
            goto load_map_ret_MISSING_TOKENS;
          }
        }
        int32_t cur_bp;
        if (scan_int_abs_defcap(loadbuf_iter, &cur_bp)) {
          sprintf(g_logbuf, "Error: Invalid bp coordinate on line %" PRIuPTR " of %s.\n", line_idx, mapname);
          goto load_map_ret_MALFORMED_INPUT_WW;
        }
        if (cur_bp < 0) {
          goto load_map_skip_variant;
        }

        const uint32_t variant_idx_lowbits = variant_ct % kLoadMapBlockSize;
        if (!variant_idx_lowbits) {
          if ((uintptr_t)(tmp_alloc_end - tmp_alloc_base) <= kLoadMapBlockSize * (sizeof(int16_t) + sizeof(int32_t) + sizeof(intptr_t) + sizeof(double))) {
            goto load_map_ret_NOMEM;
          }
          cur_chr_codes = (uint16_t*)tmp_alloc_base;
          tmp_alloc_base = (unsigned char*)(&(cur_chr_codes[kLoadMapBlockSize]));
          cur_bps = (uint32_t*)tmp_alloc_base;
          tmp_alloc_base = (unsigned char*)(&(cur_bps[kLoadMapBlockSize]));
          cur_ids = (char**)tmp_alloc_base;
          tmp_alloc_base = (unsigned char*)(&(cur_ids[kLoadMapBlockSize]));
          cur_cms = (double*)tmp_alloc_base;
          tmp_alloc_base = (unsigned char*)(&(cur_cms[kLoadMapBlockSize]));
        }
        cur_chr_codes[variant_idx_lowbits] = (uint32_t)cur_chr_code;
        cur_ids[variant_idx_lowbits] = (char*)tmp_alloc_end;
        cur_cms[variant_idx_lowbits] = cur_cm;
        cur_bps[variant_idx_lowbits] = (uint32_t)cur_bp;
        ++variant_ct;
      }
    load_map_skip_variant:
      ++line_idx;
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
        if (!gzeof(gz_infile)) {
          goto load_map_ret_READ_FAIL;
        }
        break;
      }
      if (!loadbuf[loadbuf_size - 1]) {
        goto load_map_ret_LONG_LINE;
      }
      loadbuf_first_token = skip_initial_spaces(loadbuf);
      if (loadbuf_first_token[0] == '#') {
        sprintf(g_logbuf, "Error: Line %" PRIuPTR " of %s starts with a '#'. (This is only permitted before the first nonheader line.)\n", line_idx, mapname);
        goto load_map_ret_MALFORMED_INPUT_WW;
      }
    }
    if (max_variant_id_slen > kMaxIdSlen) {
      logerrprint("Error: Variant names are limited to " MAX_ID_SLEN_STR " characters.\n");
      goto load_map_ret_MALFORMED_INPUT;
    }

    if (!variant_ct) {
      logerrprint("Error: All variants in .map/.bim file skipped due to chromosome filter.\n");
      goto load_map_ret_INCONSISTENT_INPUT;
    }
    tmp_alloc_base = (unsigned char*)(&(loadbuf[loadbuf_size]));
    // true requirement is weaker, but whatever
    g_bigstack_end = tmp_alloc_base;

    if (bigstack_alloc_usi(variant_ct, variant_chr_codes_ptr) ||
        bigstack_alloc_ui(variant_ct, variant_bps_ptr) ||
        bigstack_alloc_cp(variant_ct, variant_ids_ptr)) {
      goto load_map_ret_NOMEM;
    }
    uint16_t* variant_chr_codes = *variant_chr_codes_ptr;
    uint32_t* variant_bps = *variant_bps_ptr;
    char** variant_ids = *variant_ids_ptr;
    double* variant_cms = nullptr;
    if (at_least_one_nzero_cm) {
      if (bigstack_alloc_d(variant_ct, variant_cms_ptr)) {
        goto load_map_ret_NOMEM;
      }
      variant_cms = *variant_cms_ptr;
    } else {
      *variant_cms_ptr = nullptr;
    }
    *max_variant_id_slen_ptr = max_variant_id_slen;
    *variant_ct_ptr = variant_ct;
    const uint32_t full_block_ct = variant_ct / kLoadMapBlockSize;
    bigstack_mark = g_bigstack_base;
    bigstack_end_set(tmp_alloc_end);
    bigstack_end_mark = g_bigstack_end;

    unsigned char* read_iter = tmp_alloc_base;
    for (uint32_t block_idx = 0; block_idx < full_block_ct; ++block_idx) {
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
  load_map_ret_LONG_LINE:
    if (loadbuf_size == kMaxLongLine) {
      LOGERRPRINTFWW("Error: Line %" PRIuPTR " of %s is pathologically long.\n", line_idx, mapname);
      reterr = kPglRetMalformedInput;
      break;
    }
  load_map_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  load_map_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  load_map_ret_MISSING_TOKENS:
    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, mapname);
  load_map_ret_MALFORMED_INPUT_WW:
    wordwrapb(0);
    logerrprintb();
  load_map_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  load_map_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 load_map_ret_1:
  // forget_extra_chr_names(1, cip);
  gzclose_cond(gz_infile);
  bigstack_double_reset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

static_assert(sizeof(dosage_t) == 2, "plink1_dosage_to_pgen() needs to be updated.");
pglerr_t plink1_dosage_to_pgen(const char* dosagename, const char* famname, const char* mapname, const char* import_single_chr_str, const plink1_dosage_info_t* pdip, misc_flags_t misc_flags, fam_col_t fam_cols, int32_t missing_pheno, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, double import_dosage_certainty, uint32_t max_thread_ct, char* outname, char* outname_end, chr_info_t* cip) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;

  // these are not allocated on bigstack, and must be explicitly freed
  pheno_col_t* pheno_cols = nullptr;
  char* pheno_names = nullptr;
  uint32_t pheno_ct = 0;

  gzFile gz_infile = nullptr;
  FILE* outfile = nullptr;
  uintptr_t loadbuf_size = 0;
  uintptr_t line_idx = 0;
  st_pgen_writer_t spgw;
  pglerr_t reterr = kPglRetSuccess;
  spgw_preinit(&spgw);
  {
    // 1. Read .fam file.  (May as well support most .psam files too, since
    //    it's the same driver function.  However, unless 'noheader' modifier
    //    is present, SID field cannot be used for disambiguation.)
    uintptr_t max_sample_id_blen = 4;
    uintptr_t max_sid_blen = 0;
    uintptr_t max_paternal_id_blen = 2;
    uintptr_t max_maternal_id_blen = 2;
    uint32_t raw_sample_ct = 0;
    uintptr_t* sample_include = nullptr;
    char* sample_ids = nullptr;
    char* sids = nullptr;
    char* paternal_ids = nullptr;
    char* maternal_ids = nullptr;
    uintptr_t* sex_nm = nullptr;
    uintptr_t* sex_male = nullptr;
    uintptr_t* founder_info = nullptr;
    uintptr_t max_pheno_name_blen = 0;
    reterr = load_psam(famname, nullptr, fam_cols, 0x7fffffff, missing_pheno, (misc_flags / kfMiscAffection01) & 1, &max_sample_id_blen, &max_sid_blen, &max_paternal_id_blen, &max_maternal_id_blen, &sample_include, &sample_ids, &sids, &paternal_ids, &maternal_ids, &founder_info, &sex_nm, &sex_male, &pheno_cols, &pheno_names, &raw_sample_ct, &pheno_ct, &max_pheno_name_blen);
    if (reterr) {
      goto plink1_dosage_to_pgen_ret_1;
    }

    // 2. Read dosage-file header line if it exists, then write new .psam.
    reterr = gzopen_read_checked(dosagename, &gz_infile);
    if (reterr) {
      goto plink1_dosage_to_pgen_ret_1;
    }

    const uint32_t first_data_col_idx = pdip->skips[0] + pdip->skips[1] + pdip->skips[2] + 3;
    uint32_t sample_ct = 0;
    uint32_t* dosage_sample_idx_to_fam_uidx;
    if (bigstack_end_alloc_ui(raw_sample_ct, &dosage_sample_idx_to_fam_uidx)) {
      goto plink1_dosage_to_pgen_ret_NOMEM;
    }
    plink1_dosage_flags_t flags = pdip->flags;
    if (flags & kfPlink1DosageNoheader) {
      sample_ct = raw_sample_ct;
      for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
        dosage_sample_idx_to_fam_uidx[sample_idx] = sample_idx;
      }
    } else {
      fill_ulong_zero(BITCT_TO_WORDCT(raw_sample_ct), sample_include);
      const uint32_t tmp_htable_size = get_htable_fast_size(raw_sample_ct);
      uint32_t* htable_tmp;
      char* idbuf;
      if (bigstack_end_alloc_ui(tmp_htable_size, &htable_tmp) ||
          bigstack_end_alloc_c(max_sample_id_blen, &idbuf)) {
        goto plink1_dosage_to_pgen_ret_NOMEM;
      }
      const uint32_t duplicate_idx = populate_strbox_htable(sample_ids, raw_sample_ct, max_sample_id_blen, tmp_htable_size, htable_tmp);
      if (duplicate_idx) {
        char* duplicate_sample_id = &(sample_ids[duplicate_idx * max_sample_id_blen]);
        char* duplicate_fid_end = (char*)rawmemchr(duplicate_sample_id, '\t');
        *duplicate_fid_end = ' ';
        sprintf(g_logbuf, "Error: Duplicate sample ID '%s' in .fam file.\n", duplicate_sample_id);
        goto plink1_dosage_to_pgen_ret_MALFORMED_INPUT_WW;
      }

      loadbuf_size = bigstack_left();
      if (loadbuf_size > kMaxLongLine) {
        loadbuf_size = kMaxLongLine;
      } else if (loadbuf_size <= kMaxMediumLine) {
        goto plink1_dosage_to_pgen_ret_NOMEM;
      }
      // not formally allocated
      char* loadbuf = (char*)g_bigstack_base;
      loadbuf[loadbuf_size - 1] = ' ';
      char* loadbuf_first_token;
      do {
        ++line_idx;
        if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
          if (!gzeof(gz_infile)) {
            goto plink1_dosage_to_pgen_ret_READ_FAIL;
          }
          sprintf(g_logbuf, "Error: %s is empty.\n", dosagename);
          goto plink1_dosage_to_pgen_ret_INCONSISTENT_INPUT_WW;
        }
        if (!loadbuf[loadbuf_size - 1]) {
          goto plink1_dosage_to_pgen_ret_LONG_LINE;
        }
        loadbuf_first_token = skip_initial_spaces(loadbuf);
      } while (is_eoln_kns(*loadbuf_first_token));
      char* loadbuf_iter = next_token_mult(loadbuf_first_token, first_data_col_idx);
      if (!loadbuf_iter) {
        goto plink1_dosage_to_pgen_ret_MISSING_TOKENS;
      }
      const char id_delim = pdip->id_delim;
      do {
        char* fid_end = token_endnn(loadbuf_iter);
        char* iid_start;
        char* iid_end;
        uint32_t iid_slen;
        if (id_delim) {
          iid_end = fid_end;
          fid_end = (char*)memchr(loadbuf_iter, (unsigned char)id_delim, (uintptr_t)(iid_end - loadbuf_iter));
          if (!fid_end) {
            sprintf(g_logbuf, "Error: Sample ID in --import-dosage file does not contain '%c' delimiter.\n", id_delim);
            goto plink1_dosage_to_pgen_ret_INCONSISTENT_INPUT_2;
          }
          iid_start = &(fid_end[1]);
          iid_slen = (uintptr_t)(iid_end - iid_start);
          if (memchr(iid_start, (unsigned char)id_delim, (uintptr_t)(iid_end - iid_start))) {
            sprintf(g_logbuf, "Error: Sample ID in --import-dosage file has multiple instances of '%c'.\n", id_delim);
            goto plink1_dosage_to_pgen_ret_INCONSISTENT_INPUT_2;
          }
        } else {
          iid_start = skip_initial_spaces(fid_end);
          if (is_eoln_kns(*iid_start)) {
            goto plink1_dosage_to_pgen_ret_MISSING_TOKENS;
          }
          iid_end = token_endnn(iid_start);
          iid_slen = (uintptr_t)(iid_end - iid_start);
        }
        const uint32_t fid_slen = (uintptr_t)(fid_end - loadbuf_iter);
        const uint32_t cur_id_slen = fid_slen + iid_slen + 1;
        if (cur_id_slen >= max_sample_id_blen) {
          logerrprint("Error: .fam file does not contain all sample IDs in dosage file.\n");
          goto plink1_dosage_to_pgen_ret_INCONSISTENT_INPUT;
        }
        char* idbuf_iid = memcpyax(idbuf, loadbuf_iter, fid_slen, '\t');
        memcpyx(idbuf_iid, iid_start, iid_slen, '\0');
        uint32_t sample_uidx = strbox_htable_find(idbuf, sample_ids, htable_tmp, max_sample_id_blen, cur_id_slen, tmp_htable_size);
        if (sample_uidx == UINT32_MAX) {
          logerrprint("Error: .fam file does not contain all sample IDs in dosage file.\n");
          goto plink1_dosage_to_pgen_ret_INCONSISTENT_INPUT;
        }
        if (is_set(sample_include, sample_uidx)) {
          idbuf_iid[-1] = ' ';
          sprintf(g_logbuf, "Error: Duplicate sample ID '%s' in dosage file.\n", idbuf);
          goto plink1_dosage_to_pgen_ret_MALFORMED_INPUT_WW;
        }
        set_bit(sample_uidx, sample_include);
        dosage_sample_idx_to_fam_uidx[sample_ct++] = sample_uidx;
        loadbuf_iter = skip_initial_spaces(iid_end);
      } while (!is_eoln_kns(*loadbuf_iter));
    }

    strcpy(outname_end, ".psam");
    if (fopen_checked(outname, FOPEN_WB, &outfile)) {
      goto plink1_dosage_to_pgen_ret_OPEN_FAIL;
    }
    char* writebuf = g_textbuf;
    char* writebuf_flush = &(writebuf[kMaxMediumLine]);
    char* write_iter = strcpya(writebuf, "#FID\tIID");
    const uint32_t write_sid = sid_col_required(sample_include, sids, sample_ct, max_sid_blen, 1);
    if (write_sid) {
      write_iter = strcpya(write_iter, "\tSID");
    }
    const uint32_t write_parents = is_parental_info_present(sample_include, paternal_ids, maternal_ids, sample_ct, max_paternal_id_blen, max_maternal_id_blen);
    if (write_parents) {
      write_iter = strcpya(write_iter, "\tPAT\tMAT");
    }
    write_iter = strcpya(write_iter, "\tSEX");
    for (uint32_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
      *write_iter++ = '\t';
      write_iter = strcpya(write_iter, &(pheno_names[pheno_idx * max_pheno_name_blen]));
      if (fwrite_ck(writebuf_flush, outfile, &write_iter)) {
        goto plink1_dosage_to_pgen_ret_WRITE_FAIL;
      }
    }
    append_binary_eoln(&write_iter);
    uint32_t omp_slen = 2;
    char output_missing_pheno[kMaxMissingPhenostrBlen];
    if (misc_flags & kfMiscKeepAutoconv) {
      omp_slen = strlen(g_output_missing_pheno);
      memcpy(output_missing_pheno, g_output_missing_pheno, omp_slen);
    } else {
      memcpy(output_missing_pheno, "NA", 2);
    }
    for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
      const uint32_t sample_uidx = dosage_sample_idx_to_fam_uidx[sample_idx];
      write_iter = strcpya(write_iter, &(sample_ids[sample_uidx * max_sample_id_blen]));
      if (write_sid) {
        *write_iter++ = '\t';
        write_iter = strcpya(write_iter, &(sids[sample_uidx * max_sid_blen]));
      }
      if (write_parents) {
        *write_iter++ = '\t';
        write_iter = strcpyax(write_iter, &(paternal_ids[max_paternal_id_blen * sample_uidx]), '\t');
        write_iter = strcpya(write_iter, &(maternal_ids[max_maternal_id_blen * sample_uidx]));
      }
      *write_iter++ = '\t';
      if (is_set(sex_nm, sample_uidx)) {
        *write_iter++ = '2' - is_set(sex_male, sample_uidx);
      } else {
        write_iter = strcpya(write_iter, "NA");
      }
      for (uint32_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
        if (fwrite_ck(writebuf_flush, outfile, &write_iter)) {
          goto plink1_dosage_to_pgen_ret_WRITE_FAIL;
        }
        *write_iter++ = '\t';
        write_iter = append_pheno_str(&(pheno_cols[pheno_idx]), output_missing_pheno, omp_slen, sample_uidx, write_iter);
      }
      append_binary_eoln(&write_iter);
      if (fwrite_ck(writebuf_flush, outfile, &write_iter)) {
        goto plink1_dosage_to_pgen_ret_WRITE_FAIL;
      }
    }
    if (fclose_flush_null(writebuf_flush, write_iter, &outfile)) {
      goto plink1_dosage_to_pgen_ret_WRITE_FAIL;
    }
    // Don't need sample info any more.
    bigstack_end_reset(bigstack_end_mark);

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
    finalize_chrset(misc_flags, cip);
    if (mapname) {
      reterr = load_map(mapname, misc_flags, cip, &max_variant_id_slen, &variant_chr_codes, &variant_bps, &variant_ids, &variant_cms, &map_variant_ct);
      if (reterr) {
        goto plink1_dosage_to_pgen_ret_1;
      }
      const uint32_t map_variant_ctl = BITCT_TO_WORDCT(map_variant_ct);
      if (bigstack_alloc_ul(map_variant_ctl, &variant_already_seen)) {
        goto plink1_dosage_to_pgen_ret_NOMEM;
      }
      fill_all_bits(map_variant_ct, variant_already_seen);
      unsigned char* bigstack_end_mark2 = g_bigstack_end;
      g_bigstack_end = &(g_bigstack_base[round_down_pow2(bigstack_left() / 2, kEndAllocAlign)]); // allow hash table to only use half of available memory
      reterr = alloc_and_populate_id_htable_mt(variant_already_seen, variant_ids, map_variant_ct, max_thread_ct, &variant_id_htable, nullptr, &variant_id_htable_size);
      if (reterr) {
        goto plink1_dosage_to_pgen_ret_1;
      }
      g_bigstack_end = bigstack_end_mark2;
      fill_ulong_zero(map_variant_ctl, variant_already_seen);
    }

    // 4. Dosage file pass 1: count variants, check whether any decimal dosages
    //    need to be saved, write .pvar.
    //
    // Lots of overlap with ox_gen_to_pgen().
    loadbuf_size = bigstack_left() / 2;
    if (loadbuf_size <= kMaxMediumLine) {
      goto plink1_dosage_to_pgen_ret_NOMEM;
    }
    loadbuf_size -= kMaxMediumLine;
    if (loadbuf_size > kMaxLongLine) {
      loadbuf_size = kMaxLongLine;
    } else {
      loadbuf_size = round_up_pow2(loadbuf_size, kCacheline);
    }
    char* loadbuf = (char*)bigstack_alloc_raw(loadbuf_size);
    loadbuf[loadbuf_size - 1] = ' ';
    writebuf = (char*)bigstack_alloc_raw(kMaxMediumLine + loadbuf_size);
    writebuf_flush = &(writebuf[kMaxMediumLine]);
    const uint32_t allow_extra_chrs = (misc_flags / kfMiscAllowExtraChrs) & 1;
    const char* single_chr_str = nullptr;
    uint32_t single_chr_slen = 0;
    const uint32_t chr_col_idx = pdip->chr_col_idx;
    const uint32_t check_chr_col = (chr_col_idx != UINT32_MAX);
    if (!check_chr_col) {
      if (import_single_chr_str) {
        int32_t chr_code_raw = get_chr_code_raw(import_single_chr_str);
        if (chr_code_raw == -1) {
          // command-line parser guarantees that allow_extra_chrs is true here
          single_chr_str = import_single_chr_str;
          single_chr_slen = strlen(import_single_chr_str);
        } else {
          uint32_t chr_code = chr_code_raw;
          if (chr_code > cip->max_code) {
            if (chr_code < kMaxContigs) {
              logerrprint("Error: --import-dosage single-chr= code is not in the chromosome set.\n");
              goto plink1_dosage_to_pgen_ret_INVALID_CMDLINE;
            }
            chr_code = cip->xymt_codes[chr_code - kMaxContigs];
            if (((int32_t)chr_code) < 0) {
              logerrprint("Error: --import-dosage single-chr= code is not in the chromosome set.\n");
              goto plink1_dosage_to_pgen_ret_INVALID_CMDLINE;
            }
          }
          if (!is_set(cip->chr_mask, chr_code)) {
            // could permit this in --allow-no-vars case, but it's silly
            logerrprint("Error: --import-dosage single-chr= code is excluded by chromosome filter.\n");
            goto plink1_dosage_to_pgen_ret_INVALID_CMDLINE;
          }
          char* chr_buf = (char*)bigstack_alloc_raw(kCacheline);
          char* chr_name_end = chr_name_write(cip, chr_code, chr_buf);
          single_chr_str = chr_buf;
          single_chr_slen = (uintptr_t)(chr_name_end - chr_buf);
        }
      } else {
        // default to "chr0"
        if (!is_set(cip->chr_mask, 0)) {
          logerrprint("Error: No --import-dosage chromosome information specified, and chr0 excluded.\n");
          goto plink1_dosage_to_pgen_ret_INVALID_CMDLINE;
        }
        char* chr_buf = (char*)bigstack_alloc_raw(kCacheline);
        char* chr_name_end = chr_name_write(cip, 0, chr_buf);
        single_chr_str = chr_buf;
        single_chr_slen = (uintptr_t)(chr_name_end - chr_buf);
      }
    }
    strcpy(outname_end, ".pvar");
    if (fopen_checked(outname, FOPEN_WB, &outfile)) {
      goto plink1_dosage_to_pgen_ret_OPEN_FAIL;
    }
    write_iter = strcpya(writebuf, "#CHROM\tPOS\tID\tREF\tALT");
    if (variant_cms) {
      write_iter = memcpyl3a(write_iter, "\tCM");
    }
    append_binary_eoln(&write_iter);
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
    parse_table[0] = (((uint64_t)id_col_idx) << 32) + 2;
    parse_table[1] = (((uint64_t)ref_col_idx) << 32) + 3;
    parse_table[2] = (((uint64_t)alt_col_idx) << 32) + 4;
    uint32_t relevant_initial_col_ct = 3;
    if (check_chr_col) {
      parse_table[relevant_initial_col_ct++] = ((uint64_t)chr_col_idx) << 32;
    }
    const uint32_t check_pos_col = (pdip->pos_col_idx != UINT32_MAX);
    if (check_pos_col) {
      parse_table[relevant_initial_col_ct++] = (((uint64_t)(pdip->pos_col_idx)) << 32) + 1;
    }
    qsort(parse_table, relevant_initial_col_ct, sizeof(int64_t), uint64cmp);
    uint32_t col_skips[6];
    uint32_t col_types[6];
    for (uint32_t uii = 0; uii < relevant_initial_col_ct; ++uii) {
      const uint64_t parse_table_entry = parse_table[uii];
      col_skips[uii] = parse_table_entry >> 32;
      col_types[uii] = (uint32_t)parse_table_entry;
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
    const uint32_t dosage_erase_halfdist = kDosage4th - dosage_erase_thresh;
    uint32_t dosage_is_present = 0;
    uint32_t variant_ct = 0;
    uintptr_t variant_skip_ct = 0;
    uint32_t variant_uidx = 0;
    while (1) {
      ++line_idx;
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
        if (!gzeof(gz_infile)) {
          goto plink1_dosage_to_pgen_ret_READ_FAIL;
        }
        break;
      }
      if (!loadbuf[loadbuf_size - 1]) {
        goto plink1_dosage_to_pgen_ret_LONG_LINE_N;
      }
      char* loadbuf_iter = skip_initial_spaces(loadbuf);
      if (is_eoln_kns(*loadbuf_iter)) {
        continue;
      }
      char* token_ptrs[6];
      uint32_t token_slens[6];
      for (uint32_t ric_col_idx = 0; ric_col_idx < relevant_initial_col_ct; ++ric_col_idx) {
        const uint32_t cur_col_type = col_types[ric_col_idx];
        loadbuf_iter = next_token_multz(loadbuf_iter, col_skips[ric_col_idx]);
        if (!loadbuf_iter) {
          goto plink1_dosage_to_pgen_ret_MISSING_TOKENS;
        }
        token_ptrs[cur_col_type] = loadbuf_iter;
        char* token_end = token_endnn(loadbuf_iter);
        token_slens[cur_col_type] = (uintptr_t)(token_end - loadbuf_iter);
        loadbuf_iter = token_end;
      }
      // ID
      const char* variant_id = token_ptrs[2];
      const uint32_t variant_id_slen = token_slens[2];
      if (map_variant_ct) {
        variant_uidx = variant_id_dupflag_htable_find(variant_id, variant_ids, variant_id_htable, variant_id_slen, variant_id_htable_size, max_variant_id_slen);
        if (variant_uidx >> 31) {
          if (variant_uidx == UINT32_MAX) {
            ++variant_skip_ct;
            continue;
          }
          sprintf(g_logbuf, "Error: Variant ID '%s' appears multiple times in .map file.\n", variant_ids[variant_uidx & 0x7fffffff]);
          goto plink1_dosage_to_pgen_ret_MALFORMED_INPUT_WW;
        }
        if (is_set(variant_already_seen, variant_uidx)) {
          sprintf(g_logbuf, "Error: Variant ID '%s' appears multiple times in --import-dosage file.\n", variant_ids[variant_uidx]);
          goto plink1_dosage_to_pgen_ret_MALFORMED_INPUT_WW;
        }
        // already performed chromosome filtering
        write_iter = chr_name_write(cip, (uint32_t)variant_chr_codes[variant_uidx], write_iter);
        *write_iter++ = '\t';
        write_iter = uint32toa_x(variant_bps[variant_uidx], '\t', write_iter);
        write_iter = memcpya(write_iter, variant_id, variant_id_slen);
      } else {
        if (variant_id_slen > kMaxIdSlen) {
          putc_unlocked('\n', stdout);
          logerrprint("Error: Variant names are limited to " MAX_ID_SLEN_STR " characters.\n");
          goto plink1_dosage_to_pgen_ret_MALFORMED_INPUT;
        }
        // #CHROM
        if (check_chr_col) {
          char* chr_code_str = token_ptrs[0];
          char* chr_code_end = &(chr_code_str[token_slens[0]]);
          int32_t cur_chr_code;
          reterr = get_or_add_chr_code_destructive("--import-dosage file", line_idx, allow_extra_chrs, chr_code_str, chr_code_end, cip, &cur_chr_code);
          if (reterr) {
            goto plink1_dosage_to_pgen_ret_1;
          }
          if (!is_set(cip->chr_mask, cur_chr_code)) {
            ++variant_skip_ct;
            continue;
          }
          write_iter = chr_name_write(cip, cur_chr_code, write_iter);
        } else {
          write_iter = memcpya(write_iter, single_chr_str, single_chr_slen);
        }
        *write_iter++ = '\t';
        // POS
        if (check_pos_col) {
          char* pos_str = token_ptrs[1];
          // no need to support negative values here
          uint32_t cur_bp;
          if (scan_uint_defcap(pos_str, &cur_bp)) {
            sprintf(g_logbuf, "Error: Invalid bp coordinate on line %" PRIuPTR " of %s.\n", line_idx, dosagename);
            goto plink1_dosage_to_pgen_ret_MALFORMED_INPUT_WW;
          }
          write_iter = uint32toa(cur_bp, write_iter);
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
      append_binary_eoln(&write_iter);
      if (fwrite_ck(writebuf_flush, outfile, &write_iter)) {
        goto plink1_dosage_to_pgen_ret_WRITE_FAIL;
      }
      if (!dosage_is_present) {
        loadbuf_iter = token_ptrs[5];
        if (format_infer) {
          const uint32_t remaining_col_ct = count_tokens(loadbuf_iter);
          if (remaining_col_ct == sample_ct) {
            flags |= kfPlink1DosageFormatSingle;
          } else if (remaining_col_ct == sample_ct * 3) {
            format_triple = 1;
          } else if (remaining_col_ct != sample_ct * 2) {
            sprintf(g_logbuf, "Error: Unexpected format=infer column count in --import-dosage file (%u; should be %u, %u, or %u).\n", remaining_col_ct, sample_ct, sample_ct * 2, sample_ct * 3);
            goto plink1_dosage_to_pgen_ret_INCONSISTENT_INPUT_WW;
          }
          format_infer = 0;
        }
        if (flags & kfPlink1DosageFormatSingle) {
          for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
            if (!loadbuf_iter) {
              goto plink1_dosage_to_pgen_ret_MISSING_TOKENS;
            }
            double a1_dosage;
            char* str_end = scanadv_double(loadbuf_iter, &a1_dosage);
            if ((!loadbuf_iter) || (a1_dosage < (0.5 / 32768.0)) || (a1_dosage >= dosage_ceil)) {
              loadbuf_iter = next_token(loadbuf_iter);
              continue;
            }
            a1_dosage *= dosage_multiplier;
            const uint32_t dosage_int = (uint32_t)(a1_dosage + 0.5);
            const uint32_t halfdist = biallelic_dosage_halfdist(dosage_int);
            if (halfdist < dosage_erase_halfdist) {
              dosage_is_present = 1;
              break;
            }
            loadbuf_iter = next_token(str_end);
          }
        } else {
          // for compatibility with plink 1.x, do not actually parse third
          // value of each triplet if format=3
          for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
            if (!loadbuf_iter) {
              goto plink1_dosage_to_pgen_ret_MISSING_TOKENS;
            }
            double prob_2a1;
            char* str_end = scanadv_double(loadbuf_iter, &prob_2a1);
            if (!str_end) {
              loadbuf_iter = next_token_mult(loadbuf_iter, 2 + format_triple);
              continue;
            }
            loadbuf_iter = next_token(str_end);
            if (!loadbuf_iter) {
              goto plink1_dosage_to_pgen_ret_MISSING_TOKENS;
            }
            double prob_1a1;
            str_end = scanadv_double(loadbuf_iter, &prob_1a1);
            if (!str_end) {
              loadbuf_iter = next_token_mult(loadbuf_iter, 1 + format_triple);
              continue;
            }
            loadbuf_iter = next_token_mult(str_end, 1 + format_triple);
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
            const uint32_t dosage_int = (uint32_t)(prob_2a1 * 32768 + prob_1a1 * 16384 + 0.5);
            const uint32_t halfdist = biallelic_dosage_halfdist(dosage_int);
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
    putc_unlocked('\r', stdout);
    if (fclose_flush_null(writebuf_flush, write_iter, &outfile)) {
      goto plink1_dosage_to_pgen_ret_WRITE_FAIL;
    }
    if (!variant_ct) {
      if (!variant_skip_ct) {
        logerrprint("Error: Empty --import-dosage file.\n");
        goto plink1_dosage_to_pgen_ret_INCONSISTENT_INPUT;
      }
      LOGERRPRINTFWW("Error: All %" PRIuPTR " variant%s in --import-dosage file skipped.\n", variant_skip_ct, (variant_skip_ct == 1)? "" : "s");
      goto plink1_dosage_to_pgen_ret_INCONSISTENT_INPUT;
    }
    LOGPRINTF("--import-dosage: %u variant%s scanned%s.\n", variant_ct, (variant_ct == 1)? "" : "s", dosage_is_present? "" : " (all hardcalls)");

    // 5. Dosage file pass 2: write .pgen.
    bigstack_reset(writebuf);
    if (gzrewind(gz_infile)) {
      goto plink1_dosage_to_pgen_ret_READ_FAIL;
    }
    const uintptr_t line_ct = line_idx - 1;
    line_idx = 0;
    if (!(flags & kfPlink1DosageNoheader)) {
      // skip header line again
      char* loadbuf_first_token;
      do {
        ++line_idx;
        if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
          goto plink1_dosage_to_pgen_ret_READ_FAIL;
        }
        loadbuf_first_token = skip_initial_spaces(loadbuf);
      } while (is_eoln_kns(*loadbuf_first_token));
    }
    strcpy(outname_end, ".pgen");
    uintptr_t spgw_alloc_cacheline_ct;
    uint32_t max_vrec_len;
    reterr = spgw_init_phase1(outname, nullptr, nullptr, variant_ct, sample_ct, dosage_is_present? kfPgenGlobalDosagePresent: kfPgenGlobal0, (flags & (kfPlink1DosageRefFirst | kfPlink1DosageRefSecond))? 1 : 2, &spgw, &spgw_alloc_cacheline_ct, &max_vrec_len);
    if (reterr) {
      goto plink1_dosage_to_pgen_ret_1;
    }
    unsigned char* spgw_alloc;
    if (bigstack_alloc_uc(spgw_alloc_cacheline_ct * kCacheline, &spgw_alloc)) {
      goto plink1_dosage_to_pgen_ret_NOMEM;
    }
    spgw_init_phase2(max_vrec_len, &spgw, spgw_alloc);

    const uint32_t sample_ctl2 = QUATERCT_TO_WORDCT(sample_ct);
    const uint32_t sample_ctl = BITCT_TO_WORDCT(sample_ct);
    uintptr_t* genovec;
    uintptr_t* dosage_present;
    if (bigstack_alloc_ul(sample_ctl2, &genovec) ||
        bigstack_alloc_ul(sample_ctl, &dosage_present)) {
      goto plink1_dosage_to_pgen_ret_NOMEM;
    }
    dosage_t* dosage_vals = nullptr;
    if (dosage_is_present) {
      if (bigstack_alloc_dosage(sample_ct, &dosage_vals)) {
        goto plink1_dosage_to_pgen_ret_NOMEM;
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
    do {
      ++line_idx;
      if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
        goto plink1_dosage_to_pgen_ret_READ_FAIL;
      }
      char* loadbuf_iter = skip_initial_spaces(loadbuf);
      if (is_eoln_kns(*loadbuf_iter)) {
        continue;
      }
      if (variant_skip_ct) {
        if (map_variant_ct) {
          char* variant_id = next_token_multz(loadbuf_iter, id_col_idx);
          const uint32_t variant_id_slen = strlen_se(variant_id);
          if (variant_id_dupflag_htable_find(variant_id, variant_ids, variant_id_htable, variant_id_slen, variant_id_htable_size, max_variant_id_slen) == UINT32_MAX) {
            continue;
          }
          loadbuf_iter = next_token_mult(variant_id, first_data_col_idx - id_col_idx);
        } else {
          char* chr_code_str = next_token_multz(loadbuf_iter, chr_col_idx);
          char* chr_code_end = token_endnn(chr_code_str);
          loadbuf_iter = next_token_mult(chr_code_end, first_data_col_idx - chr_col_idx);
          *chr_code_end = '\0';
          const uint32_t chr_code = get_chr_code(chr_code_str, cip, (uintptr_t)(chr_code_end - chr_code_str));
          if (!is_set(cip->chr_mask, chr_code)) {
            continue;
          }
        }
      } else {
        loadbuf_iter = next_token_mult(loadbuf_iter, first_data_col_idx);
      }
      uint32_t inner_loop_last = kBitsPerWordD2 - 1;
      uint32_t widx = 0;
      dosage_t* dosage_vals_iter = dosage_vals;
      while (1) {
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
            if (!loadbuf_iter) {
              goto plink1_dosage_to_pgen_ret_MISSING_TOKENS;
            }
            double a1_dosage;
            char* str_end = scanadv_double(loadbuf_iter, &a1_dosage);
            if ((!str_end) || (a1_dosage < 0.0) || (a1_dosage > dosage_ceil)) {
              genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
              loadbuf_iter = next_token(loadbuf_iter);
              continue;
            }
            loadbuf_iter = next_token(str_end);
            uint32_t dosage_int = (uint32_t)(a1_dosage * dosage_multiplier + 0.5);
            if (dosage_int > kDosageMax) {
              dosage_int = kDosageMax;
            }
            const uint32_t cur_halfdist = biallelic_dosage_halfdist(dosage_int);
            if (cur_halfdist < hard_call_halfdist) {
              genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
            } else {
              genovec_word |= ((dosage_int + (kDosage4th * k1LU)) / kDosageMid) << (2 * sample_idx_lowbits);
              if (cur_halfdist >= dosage_erase_halfdist) {
                continue;
              }
            }
            dosage_present_hw |= 1U << sample_idx_lowbits;
            *dosage_vals_iter++ = dosage_int;
          }
        } else {
          for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
            if (!loadbuf_iter) {
              goto plink1_dosage_to_pgen_ret_MISSING_TOKENS;
            }
            double prob_2a1;
            char* str_end = scanadv_double(loadbuf_iter, &prob_2a1);
            if (!str_end) {
              genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
              loadbuf_iter = next_token_mult(loadbuf_iter, 2 + format_triple);
              continue;
            }
            loadbuf_iter = next_token(str_end);
            if (!loadbuf_iter) {
              goto plink1_dosage_to_pgen_ret_MISSING_TOKENS;
            }
            double prob_1a1;
            str_end = scanadv_double(loadbuf_iter, &prob_1a1);
            if (!str_end) {
              genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
              loadbuf_iter = next_token_mult(loadbuf_iter, 1 + format_triple);
              continue;
            }
            loadbuf_iter = next_token_mult(str_end, 1 + format_triple);
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
            const uint32_t dosage_int = (uint32_t)(prob_2a1 * 32768 + prob_1a1 * 16384 + 0.5);
            const uint32_t cur_halfdist = biallelic_dosage_halfdist(dosage_int);
            if (cur_halfdist < hard_call_halfdist) {
              genovec_word |= (3 * k1LU) << (2 * sample_idx_lowbits);
            } else {
              genovec_word |= ((dosage_int + (kDosage4th * k1LU)) / kDosageMid) << (2 * sample_idx_lowbits);
              if (cur_halfdist >= dosage_erase_halfdist) {
                continue;
              }
            }
            dosage_present_hw |= 1U << sample_idx_lowbits;
            *dosage_vals_iter++ = dosage_int;
          }
        }
        genovec[widx] = genovec_word;
        ((halfword_t*)dosage_present)[widx] = (halfword_t)dosage_present_hw;
        ++widx;
      }
      if (!prov_ref_allele_second) {
        genovec_invert_unsafe(sample_ct, genovec);
        zero_trailing_quaters(sample_ct, genovec);
      }
      if (dosage_vals_iter != dosage_vals) {
        const uint32_t dosage_ct = (uintptr_t)(dosage_vals_iter - dosage_vals);
        if (!prov_ref_allele_second) {
          biallelic_dosage16_invert(dosage_ct, dosage_vals);
        }
        if (spgw_append_biallelic_genovec_dosage16(genovec, dosage_present, dosage_vals, dosage_ct, &spgw)) {
          goto plink1_dosage_to_pgen_ret_WRITE_FAIL;
        }
      } else {
        if (spgw_append_biallelic_genovec(genovec, &spgw)) {
          goto plink1_dosage_to_pgen_ret_WRITE_FAIL;
        }
      }
      ++vidx;
      if (!(vidx % 1000)) {
        printf("\r--import-dosage: %uk variants converted.", vidx / 1000);
        fflush(stdout);
      }
    } while (line_idx < line_ct);
    spgw_finish(&spgw);
    putc_unlocked('\r', stdout);
    write_iter = strcpya(g_logbuf, "--import-dosage: ");
    const uint32_t outname_base_slen = (uintptr_t)(outname_end - outname);
    write_iter = memcpya(write_iter, outname, outname_base_slen + 5);
    write_iter = memcpyl3a(write_iter, " + ");
    write_iter = memcpya(write_iter, outname, outname_base_slen);
    write_iter = strcpya(write_iter, ".pvar + ");
    write_iter = memcpya(write_iter, outname, outname_base_slen);
    write_iter = strcpya(write_iter, ".psam written.\n");
    wordwrapb(0);
    logprintb();
  }
  while (0) {
  plink1_dosage_to_pgen_ret_LONG_LINE_N:
    putc_unlocked('\n', stdout);
  plink1_dosage_to_pgen_ret_LONG_LINE:
    if (loadbuf_size == kMaxLongLine) {
      LOGERRPRINTFWW("Error: Line %" PRIuPTR " of %s is pathologically long.\n", line_idx, dosagename);
      reterr = kPglRetMalformedInput;
      break;
    }
  plink1_dosage_to_pgen_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  plink1_dosage_to_pgen_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  plink1_dosage_to_pgen_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  plink1_dosage_to_pgen_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  plink1_dosage_to_pgen_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  plink1_dosage_to_pgen_ret_MISSING_TOKENS:
    sprintf(g_logbuf, "Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, dosagename);
  plink1_dosage_to_pgen_ret_MALFORMED_INPUT_WW:
    wordwrapb(0);
    putc_unlocked('\n', stdout);
    logerrprintb();
  plink1_dosage_to_pgen_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  plink1_dosage_to_pgen_ret_INCONSISTENT_INPUT_WW:
    wordwrapb(0);
  plink1_dosage_to_pgen_ret_INCONSISTENT_INPUT_2:
    logerrprintb();
  plink1_dosage_to_pgen_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 plink1_dosage_to_pgen_ret_1:
  if (spgw_cleanup(&spgw) && (!reterr)) {
    reterr = kPglRetWriteFail;
  }
  forget_extra_chr_names(1, cip);
  fclose_cond(outfile);
  gzclose_cond(gz_infile);
  free_cond(pheno_names);
  cleanup_pheno_cols(pheno_ct, pheno_cols);
  bigstack_double_reset(bigstack_mark, bigstack_end_mark);
  return reterr;
}


// binary search over cdf is faster than (int)(log(drand)/log(q)) for truncated
// geometric distribution
static uint64_t g_geno_missing_geomdist[kBitsPerWordD2];
static uint64_t g_dosage_geomdist[kBitsPerWordD2];
static uint32_t g_geno_missing_invert = 0;
static uint32_t g_dosage_geomdist_max = 0;

static_assert(sizeof(dosage_t) == 2, "generate_dummy_thread() needs to be updated.");
THREAD_FUNC_DECL generate_dummy_thread(void* arg) {
  const uintptr_t tidx = (uintptr_t)arg;
  const uint32_t sample_ct = g_sample_ct;
  const uint32_t calc_thread_ct = g_calc_thread_ct;
  const uint64_t* geno_missing_geomdist = g_geno_missing_geomdist;
  const uint64_t* dosage_geomdist = g_dosage_geomdist;
  const uint32_t geno_missing_invert = g_geno_missing_invert;
  const uint32_t geno_missing_check = geno_missing_invert || (geno_missing_geomdist[kBitsPerWordD2 - 1] != 0);
  const uint32_t dosage_geomdist_max = g_dosage_geomdist_max;
  const uint32_t dosage_is_present = (dosage_geomdist_max != kBitsPerWord);
  const uint32_t hard_call_halfdist = g_hard_call_halfdist;
  const uint32_t dosage_erase_halfdist = g_dosage_erase_halfdist;
  const uintptr_t sample_ctaw2 = QUATERCT_TO_ALIGNED_WORDCT(sample_ct);
  const uint32_t sample_ctl2_m1 = (sample_ct - 1) / kBitsPerWordD2;
  const uintptr_t sample_ctaw = BITCT_TO_ALIGNED_WORDCT(sample_ct);
  sfmt_t* sfmtp = g_sfmtp_arr[tidx];
  uint64_t ullrand = sfmt_genrand_uint64(sfmtp);
  uint32_t rand16_left = 4;
  uint32_t parity = 0;
  while (1) {
    const uint32_t is_last_block = g_is_last_thread_block;
    const uintptr_t cur_block_write_ct = g_cur_block_write_ct;
    uint32_t vidx = (tidx * cur_block_write_ct) / calc_thread_ct;
    const uint32_t vidx_end = ((tidx + 1) * cur_block_write_ct) / calc_thread_ct;
    uintptr_t* write_genovec_iter = &(g_write_genovecs[parity][vidx * sample_ctaw2]);
    uint32_t* write_dosage_ct_iter = &(g_write_dosage_cts[parity][vidx]);
    uintptr_t* write_dosage_present_iter = &(g_write_dosage_presents[parity][vidx * sample_ctaw]);

    // bugfix (23 Jul 2017): multiply by sample_ct, not sample_ctaw
    dosage_t* write_dosage_vals_iter = &(g_write_dosage_val_bufs[parity][vidx * sample_ct]);
    for (; vidx < vidx_end; ++vidx) {
      dosage_t* cur_dosage_vals_iter = write_dosage_vals_iter;
      uint32_t loop_len = kBitsPerWordD2;
      uint32_t widx = 0;
      while (1) {
        if (widx >= sample_ctl2_m1) {
          if (widx > sample_ctl2_m1) {
            break;
          }
          loop_len = MOD_NZ(sample_ct, kBitsPerWordD2);
        }
        // sfmt_genrand_uint64 calls can't be mixed with sfmt_genrand_uint32
        // calls, so use it here even in 32-bit build
        uintptr_t genovec_word = sfmt_genrand_uint64(sfmtp);
        genovec_word = genovec_word - ((genovec_word >> 1) & kMask5555);
        if (geno_missing_check) {
          uintptr_t missing_mask = 0;
          uint32_t sample_idx_lowbits = 0;
          while (1) {
            sample_idx_lowbits += uint64arr_geq(geno_missing_geomdist, kBitsPerWordD2, sfmt_genrand_uint64(sfmtp));
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
              sample_idx_lowbits += uint64arr_geq(dosage_geomdist, dosage_geomdist_max, sfmt_genrand_uint64(sfmtp));
            }
            if (sample_idx_lowbits >= loop_len) {
              break;
            }
            if (((genovec_word >> (2 * sample_idx_lowbits)) & 3) == 3) {
              continue;
            }
            if (!rand16_left) {
              ullrand = sfmt_genrand_uint64(sfmtp);
              rand16_left = 4;
            }
            const uint32_t dosage_int = ((ullrand & 65535) + 1) / 2;
            ullrand >>= 16;
            --rand16_left;
            const uint32_t halfdist = biallelic_dosage_halfdist(dosage_int);
            if (halfdist < dosage_erase_halfdist) {
              *cur_dosage_vals_iter++ = dosage_int;
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
        ((halfword_t*)write_dosage_present_iter)[widx] = (halfword_t)dosage_present_hw;
        ++widx;
      }
      zero_trailing_quaters(sample_ct, write_genovec_iter);
      const uint32_t dosage_ct = (uintptr_t)(cur_dosage_vals_iter - write_dosage_vals_iter);
      *write_dosage_ct_iter++ = dosage_ct;
      write_genovec_iter = &(write_genovec_iter[sample_ctaw2]);
      write_dosage_present_iter = &(write_dosage_present_iter[sample_ctaw]);
      write_dosage_vals_iter = &(write_dosage_vals_iter[sample_ct]);
    }
    if (is_last_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
    parity = 1 - parity;
  }
}

static_assert(sizeof(dosage_t) == 2, "generate_dummy() needs to be updated.");
pglerr_t generate_dummy(const gendummy_info_t* gendummy_info_ptr, misc_flags_t misc_flags, uint32_t hard_call_thresh, uint32_t dosage_erase_thresh, uint32_t max_thread_ct, char* outname, char* outname_end, chr_info_t* cip) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  threads_state_t ts;
  init_threads3z(&ts);
  st_pgen_writer_t spgw;
  pglerr_t reterr = kPglRetSuccess;
  spgw_preinit(&spgw);
  {
    finalize_chrset(misc_flags, cip);
    if (!is_set(cip->chr_mask, 1)) {
      logerrprint("Error: --dummy cannot be used when chromosome 1 is excluded.\n");
      goto generate_dummy_ret_INVALID_CMDLINE;
    }
    if (is_set(cip->haploid_mask, 1)) {
      logerrprint("Error: --dummy cannot be used to generate haploid data.\n");
      goto generate_dummy_ret_INVALID_CMDLINE;
    }
    char chr1_name_buf[5];
    char* chr1_name_end = chr_name_write(cip, 1, chr1_name_buf);
    *chr1_name_end = '\t';
    const uint32_t chr1_name_blen = 1 + (uintptr_t)(chr1_name_end - chr1_name_buf);
    const uint32_t sample_ct = gendummy_info_ptr->sample_ct;
    const uint32_t variant_ct = gendummy_info_ptr->variant_ct;
    // missing pheno string is always "NA"
    const gendummy_flags_t flags = gendummy_info_ptr->flags;
    uint16_t alleles[13];
    uint32_t four_alleles = 0;
    if (flags & kfGenDummyAcgt) {
      memcpy(alleles, "\tA\tC\tA\tG\tA\tT\tC\tG\tC\tT\tG\tT\tA", 26);
      four_alleles = 1;
    } else if (flags & kfGenDummy1234) {
      memcpy(alleles, "\t1\t2\t1\t3\t1\t4\t2\t3\t2\t4\t3\t4\t1", 26);
      four_alleles = 1;
    } else if (flags & kfGenDummy12) {
      memcpy(alleles, "\t1\t2\t1", 6);
    } else {
      memcpy(alleles, "\tA\tB\tA", 6);
    }
    char* textbuf = g_textbuf;
    char* textbuf_flush = &(textbuf[kMaxMediumLine]);
    strcpy(outname_end, ".pvar");
    if (fopen_checked(outname, FOPEN_WB, &outfile)) {
      goto generate_dummy_ret_OPEN_FAIL;
    }
    char* write_iter = strcpya(textbuf, "#CHROM\tPOS\tID\tREF\tALT");
    append_binary_eoln(&write_iter);
    if (four_alleles) {
      uint32_t urand = 0;
      for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx) {
        if (!(variant_idx % 8)) {
          if (fwrite_ck(textbuf_flush, outfile, &write_iter)) {
            goto generate_dummy_ret_WRITE_FAIL;
          }
          do {
            urand = sfmt_genrand_uint32(&g_sfmt);
          } while (urand < 425132032U); // 2^32 - 12^8
        }
        const uint32_t quotient = urand / 12;
        const uint32_t remainder = urand - (quotient * 12U);
        urand = quotient;
        write_iter = memcpya(write_iter, chr1_name_buf, chr1_name_blen);
        write_iter = uint32toa(variant_idx, write_iter);
        write_iter = strcpya(write_iter, "\tsnp");
        write_iter = uint32toa(variant_idx, write_iter);
        write_iter = memcpya(write_iter, &(alleles[remainder]), 4);
        append_binary_eoln(&write_iter);
      }
    } else {
      uint32_t urand = 0;
      for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx) {
        if (!(variant_idx % 32)) {
          if (fwrite_ck(textbuf_flush, outfile, &write_iter)) {
            goto generate_dummy_ret_WRITE_FAIL;
          }
          urand = sfmt_genrand_uint32(&g_sfmt);
        }
        const uint32_t remainder = urand & 1;
        urand >>= 1;
        write_iter = memcpya(write_iter, chr1_name_buf, chr1_name_blen);
        write_iter = uint32toa(variant_idx, write_iter);
        write_iter = strcpya(write_iter, "\tsnp");
        write_iter = uint32toa(variant_idx, write_iter);
        write_iter = memcpya(write_iter, &(alleles[remainder]), 4);
        append_binary_eoln(&write_iter);
      }
    }
    if (fclose_flush_null(textbuf_flush, write_iter, &outfile)) {
      goto generate_dummy_ret_WRITE_FAIL;
    }

    strcpy(outname_end, ".psam");
    if (fopen_checked(outname, FOPEN_WB, &outfile)) {
      goto generate_dummy_ret_OPEN_FAIL;
    }
    const uint32_t pheno_ct = gendummy_info_ptr->pheno_ct;
    char* writebuf;
    if (bigstack_alloc_c(kMaxMediumLine + 48 + pheno_ct * MAXV(kMaxMissingPhenostrBlen, 16), &writebuf)) {
      goto generate_dummy_ret_NOMEM;
    }
    char* writebuf_flush = &(writebuf[kMaxMediumLine]);
    uint32_t omp_slen = 2;
    char output_missing_pheno[kMaxMissingPhenostrBlen];
    if (misc_flags & kfMiscKeepAutoconv) {
      // must use --output-missing-phenotype parameter, which we've validated
      // to be consistent with --input-missing-phenotype
      omp_slen = strlen(g_output_missing_pheno);
      memcpy(output_missing_pheno, g_output_missing_pheno, omp_slen);
    } else {
      // use "NA" since that's always safe
      memcpy(output_missing_pheno, "NA", 2);
    }
    write_iter = strcpya(writebuf, "#FID\tIID\tSEX");
    for (uint32_t pheno_idx_p1 = 1; pheno_idx_p1 <= pheno_ct; ++pheno_idx_p1) {
      write_iter = strcpya(write_iter, "\tPHENO");
      write_iter = uint32toa(pheno_idx_p1, write_iter);
    }
    append_binary_eoln(&write_iter);
    const uint32_t pheno_m_check = (gendummy_info_ptr->pheno_mfreq >= kRecip2m32 * 0.5);
    const uint32_t pheno_m32 = (uint32_t)(gendummy_info_ptr->pheno_mfreq * 4294967296.0 - 0.5);
    if ((flags & kfGenDummyScalarPheno) && pheno_ct) {
      uint32_t saved_rnormal = 0;
      double saved_rnormal_val = 0.0;
      for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
        if (fwrite_ck(writebuf_flush, outfile, &write_iter)) {
          goto generate_dummy_ret_WRITE_FAIL;
        }
        write_iter = memcpyl3a(write_iter, "per");
        write_iter = uint32toa(sample_idx, write_iter);
        write_iter = strcpya(write_iter, "\tper");
        write_iter = uint32toa(sample_idx, write_iter);
        // could add option to add some males/unknown gender
        write_iter = strcpya(write_iter, "\t2");
        for (uint32_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
          *write_iter++ = '\t';
          if (pheno_m_check && (sfmt_genrand_uint32(&g_sfmt) <= pheno_m32)) {
            write_iter = memcpya(write_iter, output_missing_pheno, omp_slen);
          } else {
            double dxx;
            if (saved_rnormal) {
              dxx = saved_rnormal_val;
            } else {
              dxx = rand_normal(&g_sfmt, &saved_rnormal_val);
            }
            saved_rnormal_val = 1 - saved_rnormal_val;
            write_iter = dtoa_g(dxx, write_iter);
          }
        }
        append_binary_eoln(&write_iter);
      }
    } else {
      uint32_t urand = sfmt_genrand_uint32(&g_sfmt);
      uint32_t urand_bits_left = 32;
      for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
        if (fwrite_ck(writebuf_flush, outfile, &write_iter)) {
          goto generate_dummy_ret_WRITE_FAIL;
        }
        write_iter = memcpyl3a(write_iter, "per");
        write_iter = uint32toa(sample_idx, write_iter);
        write_iter = strcpya(write_iter, "\tper");
        write_iter = uint32toa(sample_idx, write_iter);
        write_iter = strcpya(write_iter, "\t2");
        for (uint32_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
          *write_iter++ = '\t';
          if (pheno_m_check && (sfmt_genrand_uint32(&g_sfmt) <= pheno_m32)) {
            write_iter = memcpya(write_iter, output_missing_pheno, omp_slen);
          } else {
            if (!urand_bits_left) {
              urand = sfmt_genrand_uint32(&g_sfmt);
              urand_bits_left = 32;
            }
            *write_iter++ = (char)((urand & 1) + '1');
            urand >>= 1;
            --urand_bits_left;
          }
        }
        append_binary_eoln(&write_iter);
      }
    }
    if (fclose_flush_null(writebuf_flush, write_iter, &outfile)) {
      goto generate_dummy_ret_WRITE_FAIL;
    }

    bigstack_reset(writebuf);
    strcpy(outname_end, ".pgen");
    const double geno_mfreq = gendummy_info_ptr->geno_mfreq;
    if (geno_mfreq < kRecip2m53) {
      // beyond this point, 1-x may just be 1
      g_geno_missing_geomdist[kBitsPerWordD2 - 1] = 0;
    } else {
      double remaining_prob = 1.0;
      g_geno_missing_invert = (geno_mfreq > 0.5);
      if (g_geno_missing_invert) {
        for (uint32_t uii = 0; uii < kBitsPerWordD2; ++uii) {
          remaining_prob *= geno_mfreq;
          g_geno_missing_geomdist[uii] = -((uint64_t)(remaining_prob * k2m64));
        }
      } else {
        const double geno_nmfreq = 1.0 - geno_mfreq;
        for (uint32_t uii = 0; uii < kBitsPerWordD2; ++uii) {
          remaining_prob *= geno_nmfreq;
          g_geno_missing_geomdist[uii] = -((uint64_t)(remaining_prob * k2m64));
        }
      }
    }
    const double dosage_nfreq = 1.0 - gendummy_info_ptr->dosage_freq;
    if (dosage_nfreq >= 1.0) {
      g_dosage_geomdist_max = kBitsPerWord; // used as a flag
    } else {
      double remaining_prob = 1.0;
      for (uint32_t uii = 0; uii < kBitsPerWordD2; ++uii) {
        remaining_prob *= dosage_nfreq;
        g_dosage_geomdist[uii] = -((uint64_t)(remaining_prob * k2m64));
      }
      uint32_t dosage_geomdist_max = kBitsPerWordD2;
      for (; dosage_geomdist_max; --dosage_geomdist_max) {
        if (g_dosage_geomdist[dosage_geomdist_max - 1] != 0) {
          break;
        }
      }
      g_dosage_geomdist_max = dosage_geomdist_max;
    }
    uintptr_t spgw_alloc_cacheline_ct;
    uint32_t max_vrec_len;
    reterr = spgw_init_phase1(outname, nullptr, nullptr, variant_ct, sample_ct, (dosage_nfreq >= 1.0)? kfPgenGlobal0 : kfPgenGlobalDosagePresent, 1, &spgw, &spgw_alloc_cacheline_ct, &max_vrec_len);
    if (reterr) {
      goto generate_dummy_ret_1;
    }
    unsigned char* spgw_alloc;
    if (bigstack_alloc_uc(spgw_alloc_cacheline_ct * kCacheline, &spgw_alloc)) {
      goto generate_dummy_ret_NOMEM;
    }
    spgw_init_phase2(max_vrec_len, &spgw, spgw_alloc);

    // thread-count-independent:
    //   (everything after "2 *" rounded up to cacheline)
    //   g_write_genovecs: 2 * sample_ctaw2 * sizeof(intptr_t) *
    //                     main_block_size
    //   g_write_dosage_cts: 2 * sizeof(int32_t) * main_block_size
    //   g_write_dosage_presents: 2 * sample_ctaw * sizeof(intptr_t) *
    //                            main_block_size
    //   g_write_dosage_val_bufs: 2 * sample_ct * sizeof(dosage_t)
    uint32_t calc_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
    // saturates around 4 compute threads, both with and without dosage
    // (todo: test this on something other than a MacBook Pro, could just be a
    // hyperthreading artifact)
    if (calc_thread_ct > 4) {
      calc_thread_ct = 4;
    }
    if (bigstack_init_sfmtp(calc_thread_ct, 0)) {
      goto generate_dummy_ret_NOMEM;
    }
    if (bigstack_alloc_thread(calc_thread_ct, &ts.threads)) {
      goto generate_dummy_ret_NOMEM;
    }
    const uint32_t sample_ctaw2 = QUATERCT_TO_ALIGNED_WORDCT(sample_ct);
    const uint32_t sample_ctaw = BITCT_TO_ALIGNED_WORDCT(sample_ct);
    uintptr_t cachelines_avail_m8 = bigstack_left() / kCacheline;
    if (cachelines_avail_m8 < 8) {
      goto generate_dummy_ret_NOMEM;
    }
    // we're making 8 allocations; be pessimistic re: rounding
    cachelines_avail_m8 -= 8;
    const uintptr_t bytes_req_per_in_block_variant = 2 * (sample_ctaw2 * sizeof(intptr_t) + sizeof(int32_t) + sample_ctaw * sizeof(intptr_t) + sample_ct * sizeof(dosage_t));
    uintptr_t main_block_size = (cachelines_avail_m8 * kCacheline) / bytes_req_per_in_block_variant;
    if (main_block_size > 65536) {
      main_block_size = 65536;
    } else if (main_block_size < 8) {
      // this threshold is arbitrary
      goto generate_dummy_ret_NOMEM;
    }
    if (calc_thread_ct > main_block_size / 8) {
      calc_thread_ct = main_block_size / 8;
    }
    ts.calc_thread_ct = calc_thread_ct;
    g_calc_thread_ct = calc_thread_ct;
    g_sample_ct = sample_ct;
    if (bigstack_alloc_ul(sample_ctaw2 * main_block_size, &(g_write_genovecs[0])) ||
        bigstack_alloc_ul(sample_ctaw2 * main_block_size, &(g_write_genovecs[1])) ||
        bigstack_alloc_ui(main_block_size, &(g_write_dosage_cts[0])) ||
        bigstack_alloc_ui(main_block_size, &(g_write_dosage_cts[1])) ||
        bigstack_alloc_ul(sample_ctaw * main_block_size, &(g_write_dosage_presents[0])) ||
        bigstack_alloc_ul(sample_ctaw * main_block_size, &(g_write_dosage_presents[1])) ||
        bigstack_alloc_dosage(sample_ct * main_block_size, &(g_write_dosage_val_bufs[0])) ||
        bigstack_alloc_dosage(sample_ct * main_block_size, &(g_write_dosage_val_bufs[1]))) {
      // this should be impossible
      assert(0);
      goto generate_dummy_ret_NOMEM;
    }
    // bugfix (3 Nov 2017): forgot to handle hard_call_thresh default value
    if (hard_call_thresh == UINT32_MAX) {
      hard_call_thresh = kDosageMid / 10;
    }
    g_hard_call_halfdist = kDosage4th - hard_call_thresh;
    g_dosage_erase_halfdist = kDosage4th - dosage_erase_thresh;

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
    uint32_t vidx_start = 0;
    uint32_t prev_block_write_ct = 0;
    uint32_t parity = 0;
    while (1) {
      uint32_t cur_block_write_ct = 0;
      if (!ts.is_last_block) {
        cur_block_write_ct = MINV(variant_ct - vidx_start, main_block_size);
      }
      if (vidx_start) {
        join_threads3z(&ts);
      }
      if (!ts.is_last_block) {
        g_cur_block_write_ct = cur_block_write_ct;
        ts.is_last_block = (vidx_start + cur_block_write_ct == variant_ct);
        ts.thread_func_ptr = generate_dummy_thread;
        if (spawn_threads3z(vidx_start, &ts)) {
          goto generate_dummy_ret_THREAD_CREATE_FAIL;
        }
      }
      parity = 1 - parity;
      if (vidx_start) {
        // write *previous* block results
        uintptr_t* write_genovec_iter = g_write_genovecs[parity];
        uint32_t* write_dosage_ct_iter = g_write_dosage_cts[parity];
        uintptr_t* write_dosage_present_iter = g_write_dosage_presents[parity];
        dosage_t* write_dosage_vals_iter = g_write_dosage_val_bufs[parity];
        for (uint32_t vidx = vidx_start - prev_block_write_ct; vidx < vidx_start; ++vidx) {
          const uint32_t cur_dosage_ct = *write_dosage_ct_iter++;
          if (!cur_dosage_ct) {
            if (spgw_append_biallelic_genovec(write_genovec_iter, &spgw)) {
              goto generate_dummy_ret_WRITE_FAIL;
            }
          } else {
            if (spgw_append_biallelic_genovec_dosage16(write_genovec_iter, write_dosage_present_iter, write_dosage_vals_iter, cur_dosage_ct, &spgw)) {
              goto generate_dummy_ret_WRITE_FAIL;
            }
          }
          write_genovec_iter = &(write_genovec_iter[sample_ctaw2]);
          write_dosage_present_iter = &(write_dosage_present_iter[sample_ctaw]);
          write_dosage_vals_iter = &(write_dosage_vals_iter[sample_ct]);
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
    spgw_finish(&spgw);

    putc_unlocked('\r', stdout);
    *outname_end = '\0';
    LOGPRINTFWW("Dummy data (%u sample%s, %u SNP%s) written to %s.pgen + %s.pvar + %s.psam .\n", sample_ct, (sample_ct == 1)? "" : "s", variant_ct, (variant_ct == 1)? "" : "s", outname, outname, outname);
  }
  while (0) {
  generate_dummy_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  generate_dummy_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  generate_dummy_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  generate_dummy_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  generate_dummy_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 generate_dummy_ret_1:
  if (spgw_cleanup(&spgw) && (!reterr)) {
    reterr = kPglRetWriteFail;
  }
  threads3z_cleanup(&ts, &g_cur_block_write_ct);
  fclose_cond(outfile);
  bigstack_reset(bigstack_mark);
  return reterr;
}


// more multithread globals
// g_sample_ct, g_calc_thread_ct, g_cur_block_write_ct, g_hard_call_halfdist,
// g_dosage_erase_halfdist, g_error_ret declared earlier

static vul_t** g_thread_vecaligned_bufs = nullptr;
static uintptr_t** g_thread_write_genovecs = nullptr;

static uintptr_t* g_plink1_smaj_loadbuf_iter = nullptr;
static pgen_writer_common_t** g_pwcs = nullptr;
static uint32_t g_stride = 0;

THREAD_FUNC_DECL plink1_smaj_transpose_thread(void* arg) {
  const uintptr_t tidx = (uintptr_t)arg;
  const uint32_t sample_ct = g_sample_ct;
  const uint32_t sample_ctaw2 = QUATERCT_TO_ALIGNED_WORDCT(sample_ct);
  const uint32_t write_batch_ct_m1 = (sample_ct - 1) / kPglQuaterTransposeBatch;
  const uint32_t calc_thread_ct = g_calc_thread_ct;
  pgen_writer_common_t* pwcp = g_pwcs[tidx];
  vul_t* vecaligned_buf = g_thread_vecaligned_bufs[tidx];
  uintptr_t* write_genovec = g_thread_write_genovecs[tidx];
  while (1) {
    const uint32_t is_last_block = g_is_last_thread_block;
    const uintptr_t cur_block_write_ct = g_cur_block_write_ct;
    const uint32_t loadbuf_ul_stride = g_stride;
    uint32_t write_idx = tidx * kPglVblockSize;
    uintptr_t* read_iter = &(g_plink1_smaj_loadbuf_iter[write_idx / kBitsPerWordD2]);
    const uint32_t write_idx_end = MINV(write_idx + kPglVblockSize, cur_block_write_ct);
    while (write_idx < write_idx_end) {
      const uintptr_t* read_iter2 = read_iter;
      // uintptr_t* write_iter = write_genovec;
      const uint32_t vblock_size = MINV(kPglQuaterTransposeBatch, write_idx_end - write_idx);
      uint32_t write_batch_idx = 0;
      uint32_t read_batch_size = kPglQuaterTransposeBatch;
      while (1) {
        if (write_batch_idx >= write_batch_ct_m1) {
          if (write_batch_idx > write_batch_ct_m1) {
            break;
          }
          read_batch_size = MOD_NZ(sample_ct, kPglQuaterTransposeBatch);
        }
        transpose_quaterblock(read_iter2, loadbuf_ul_stride, sample_ctaw2, read_batch_size, vblock_size, &(write_genovec[write_batch_idx * kPglQuaterTransposeWords]), vecaligned_buf);
        read_iter2 = &(read_iter2[kPglQuaterTransposeBatch * loadbuf_ul_stride]);
        ++write_batch_idx;
      }
      for (uint32_t uii = 0; uii < vblock_size; ++uii) {
        uintptr_t* cur_write_genovec = &(write_genovec[uii * sample_ctaw2]);
        pgr_plink1_to_plink2_inplace_unsafe(sample_ct, cur_write_genovec);
        zero_trailing_quaters(sample_ct, cur_write_genovec);
        pwc_append_biallelic_genovec(cur_write_genovec, pwcp);
      }
      write_idx += vblock_size;
      read_iter = &(read_iter[kPglQuaterTransposeWords]);
    }
    if ((tidx == calc_thread_ct - 1) || is_last_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

pglerr_t plink1_sample_major_to_pgen(const char* pgenname, uintptr_t variant_ct, uintptr_t sample_ct, uint32_t real_ref_alleles, uint32_t max_thread_ct, FILE* infile) {
  unsigned char* bigstack_mark = g_bigstack_base;
  mt_pgen_writer_t* mpgwp = nullptr;
  pglerr_t reterr = kPglRetSuccess;
  {
    // file size already validated by pgfi_init_phase1()
    LOGPRINTFWW("Sample-major .bed file detected.  Transposing to %s .\n", pgenname);
    fputs("0%", stdout);
    fflush(stdout);
    if ((!variant_ct) || (!sample_ct)) {
      // todo: hardcoded 12-byte write
      logprint("\n");
      logerrprint("Error: Zero-variant/zero-sample .pgen writing is not currently supported.\n");
      reterr = kPglRetNotYetSupported;
      goto plink1_sample_major_to_pgen_ret_1;
    }
    const uint32_t variant_ct4 = QUATERCT_TO_BYTECT(variant_ct);
    unsigned char* raw_loadbuf = nullptr;
    uint32_t raw_load_batch_size = 1;
    if (variant_ct4 < 5120) {
      // assuming 4K block size, fseek won't let us avoid reading many
      // unnecessary disk blocks
      raw_load_batch_size += 131071 / variant_ct4;
      if (bigstack_alloc_uc(raw_load_batch_size * variant_ct4, &raw_loadbuf)) {
        goto plink1_sample_major_to_pgen_ret_NOMEM;
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
    mpgw_init_phase1(nullptr, variant_ct, sample_ct, kfPgenGlobal0, &alloc_base_cacheline_ct, &mpgw_per_thread_cacheline_ct, &vrec_len_byte_ct, &vblock_cacheline_ct);
#ifndef __LP64__
    if ((mpgw_per_thread_cacheline_ct > (0x7fffffff / kCacheline)) || (vblock_cacheline_ct > (0x7fffffff / kCacheline))) {
      goto plink1_sample_major_to_pgen_ret_NOMEM;
    }
#endif

    uint32_t calc_thread_ct = DIV_UP(variant_ct, kPglVblockSize);
    if (calc_thread_ct >= max_thread_ct) {
      calc_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
    }
    mpgwp = (mt_pgen_writer_t*)bigstack_alloc((calc_thread_ct + DIV_UP(sizeof(mt_pgen_writer_t), kBytesPerWord)) * sizeof(intptr_t));
    if (!mpgwp) {
      goto plink1_sample_major_to_pgen_ret_NOMEM;
    }
    mpgwp->pgen_outfile = nullptr;
    pthread_t* threads;
    if (bigstack_alloc_thread(calc_thread_ct, &threads) ||
        bigstack_alloc_vp(calc_thread_ct, &g_thread_vecaligned_bufs) ||
        bigstack_alloc_ulp(calc_thread_ct, &g_thread_write_genovecs)) {
      goto plink1_sample_major_to_pgen_ret_NOMEM;
    }
    g_pwcs = &(mpgwp->pwcs[0]);
    uintptr_t cachelines_avail = bigstack_left() / kCacheline;
    // inner loop transposes kPglQuaterTransposeBatch variants at a time
    const uintptr_t transpose_thread_cacheline_ct = kPglQuaterTransposeBufbytes / kCacheline + QUATERCT_TO_VECCT(sample_ct) * (kPglQuaterTransposeBatch / kVecsPerCacheline);
    if (cachelines_avail < calc_thread_ct * ((uint64_t)transpose_thread_cacheline_ct)) {
      goto plink1_sample_major_to_pgen_ret_NOMEM;
    }
    for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
      g_thread_vecaligned_bufs[tidx] = (vul_t*)bigstack_alloc_raw(kPglQuaterTransposeBufbytes);
      g_thread_write_genovecs[tidx] = (uintptr_t*)bigstack_alloc_raw(QUATERCT_TO_VECCT(sample_ct) * kBytesPerVec * kPglQuaterTransposeBatch);
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
      } else if ((cachelines_avail / 3) >= alloc_base_cacheline_ct + mpgw_per_thread_cacheline_ct) {
        calc_thread_ct = 1;
      } else {
        // possible todo: simple single-threaded fallback
        // report this value since it can plausibly come up
        g_failed_alloc_attempt_size = (alloc_base_cacheline_ct + mpgw_per_thread_cacheline_ct) * kCacheline * 3;
        goto plink1_sample_major_to_pgen_ret_NOMEM;
      }
    }
    // todo: determine appropriate calc_thread_ct limit.  (should not be less
    // than 7-8.)
    unsigned char* mpgw_alloc = bigstack_alloc_raw((alloc_base_cacheline_ct + mpgw_per_thread_cacheline_ct * calc_thread_ct) * kCacheline);
    reterr = mpgw_init_phase2(pgenname, nullptr, nullptr, variant_ct, sample_ct, kfPgenGlobal0, 2 - real_ref_alleles, vrec_len_byte_ct, vblock_cacheline_ct, calc_thread_ct, mpgw_alloc, mpgwp);

    cachelines_avail = bigstack_left() / kCacheline;
    const uint64_t full_load_vecs_req = sample_ct * ((uint64_t)QUATERCT_TO_ALIGNED_WORDCT(variant_ct));
    uintptr_t* plink1_smaj_loadbuf;
    uint32_t load_multiplier;
    uint32_t cur_vidx_ct;
    if (full_load_vecs_req > cachelines_avail * kVecsPerCacheline) {
      // each iteration requires ((kPglVblockSize / 4) * calc_thread_ct *
      //   sample_ct) bytes to be loaded
      load_multiplier = cachelines_avail / ((kPglVblockSize / (4 * kCacheline)) * calc_thread_ct * ((uintptr_t)sample_ct));
      assert(load_multiplier);
      cur_vidx_ct = load_multiplier * calc_thread_ct * kPglVblockSize;
      plink1_smaj_loadbuf = (uintptr_t*)bigstack_alloc_raw_rd((cur_vidx_ct / 4) * ((uintptr_t)sample_ct));
      // bugfix (18 Nov 2017): this may be larger than variant_ct
      if (cur_vidx_ct > variant_ct) {
        cur_vidx_ct = variant_ct;
        load_multiplier = 1 + (cur_vidx_ct - 1) / (kPglVblockSize * calc_thread_ct);
      }
    } else {
      load_multiplier = 1 + ((variant_ct - 1) / (calc_thread_ct * kPglVblockSize));
      cur_vidx_ct = variant_ct;
      plink1_smaj_loadbuf = (uintptr_t*)bigstack_alloc_raw_rd(full_load_vecs_req * kBytesPerVec);
    }
    uint32_t cur_vidx_base = 0;
    uint32_t cur_vidx_ct4 = QUATERCT_TO_BYTECT(cur_vidx_ct);
    uint32_t cur_vidx_ctaw2 = QUATERCT_TO_ALIGNED_WORDCT(cur_vidx_ct);
    uint32_t pass_idx = 0;
    const uint32_t pass_ct = 1 + (variant_ct - 1) / cur_vidx_ct;
    g_sample_ct = sample_ct;
    g_stride = QUATERCT_TO_VECCT(cur_vidx_ct) * kWordsPerVec;
    g_calc_thread_ct = calc_thread_ct;
    while (1) {
      uint32_t raw_load_batch_idx = 0;
      uint32_t cur_raw_load_batch_size = raw_load_batch_size;
      uintptr_t* smaj_loadbuf_iter = plink1_smaj_loadbuf;
      ++pass_idx;
      putc_unlocked('\r', stdout);
      printf("Pass %u/%u: loading... 0%%", pass_idx, pass_ct);
      fflush(stdout);
      uint32_t pct = 0;
      uint32_t next_print_idx = raw_load_batch_ct / 100;
      const uint64_t seek_addl_offset = 3 + cur_vidx_base / 4;
      while (1) {
        if (raw_load_batch_size == 1) {
          if (fseeko(infile, seek_addl_offset + raw_load_batch_idx * ((uint64_t)variant_ct4), SEEK_SET)) {
            goto plink1_sample_major_to_pgen_ret_READ_FAIL;
          }
          if (!fread_unlocked(smaj_loadbuf_iter, cur_vidx_ct4, 1, infile)) {
            goto plink1_sample_major_to_pgen_ret_READ_FAIL;
          }
          smaj_loadbuf_iter = &(smaj_loadbuf_iter[cur_vidx_ctaw2]);
        } else {
          if (!fread_unlocked(raw_loadbuf, cur_raw_load_batch_size * variant_ct4, 1, infile)) {
            goto plink1_sample_major_to_pgen_ret_READ_FAIL;
          }
          unsigned char* raw_loadbuf_iter = &(raw_loadbuf[cur_vidx_base / 4]);
          for (uint32_t uii = 0; uii < cur_raw_load_batch_size; ++uii) {
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
          next_print_idx = (pct * ((uint64_t)raw_load_batch_ct)) / 100;
        }
      }
      const uintptr_t last_tidx = calc_thread_ct - 1;
      uint32_t load_idx = 0;
      g_cur_block_write_ct = calc_thread_ct * kPglVblockSize;
      uint32_t is_last_block;
      putc_unlocked('\r', stdout);
      printf("Pass %u/%u: transposing and compressing... 0%%", pass_idx, pass_ct);
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
          next_print_idx = (pct * ((uint64_t)load_multiplier)) / 100;
        }
        g_plink1_smaj_loadbuf_iter = &(plink1_smaj_loadbuf[load_idx * calc_thread_ct * (kPglVblockSize / kBitsPerWordD2)]);
        is_last_block = (++load_idx == load_multiplier);
        if (is_last_block) {
          g_cur_block_write_ct = cur_vidx_ct - (load_idx - 1) * calc_thread_ct * kPglVblockSize;
        }
        if (last_tidx) {
          if (spawn_threads2z(plink1_smaj_transpose_thread, last_tidx, is_last_block, threads)) {
            goto plink1_sample_major_to_pgen_ret_THREAD_CREATE_FAIL;
          }
        }
        plink1_smaj_transpose_thread((uintptr_t*)last_tidx);
        if (last_tidx) {
          join_threads2z(last_tidx, is_last_block, threads);
        }
        reterr = mpgw_flush(mpgwp);
        if (reterr) {
          if (!is_last_block) {
            g_cur_block_write_ct = 0;
            error_cleanup_threads2z(plink1_smaj_transpose_thread, last_tidx, threads);
          }
          goto plink1_sample_major_to_pgen_ret_WRITE_FAIL;
        }
      } while (!is_last_block);
      cur_vidx_base += cur_vidx_ct;
      if (cur_vidx_base == variant_ct) {
        if (pct > 10) {
          putc_unlocked('\b', stdout);
        }
        break;
      }
      fputs("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b                     ", stdout);
      // assumes pgfi_init_phase1() leaves file pointer at byte 3; otherwise,
      // necessary to put this at top of main loop
      if (fseeko(infile, 3, SEEK_SET)) {
        goto plink1_sample_major_to_pgen_ret_READ_FAIL;
      }
      if (variant_ct - cur_vidx_base <= cur_vidx_ct) {
        cur_vidx_ct = variant_ct - cur_vidx_base;
        cur_vidx_ct4 = QUATERCT_TO_BYTECT(cur_vidx_ct);
        cur_vidx_ctaw2 = QUATERCT_TO_ALIGNED_WORDCT(cur_vidx_ct);
        g_stride = QUATERCT_TO_VECCT(cur_vidx_ct) * kWordsPerVec;
        load_multiplier = 1 + (cur_vidx_ct - 1) / (kPglVblockSize * calc_thread_ct);
      }
    }
    mpgwp = nullptr;
    fputs("\b\bdone.\n", stdout);
    LOGPRINTF("Transpose complete.\n");
  }
  while (0) {
  plink1_sample_major_to_pgen_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  plink1_sample_major_to_pgen_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  plink1_sample_major_to_pgen_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  plink1_sample_major_to_pgen_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 plink1_sample_major_to_pgen_ret_1:
  if (mpgw_cleanup(mpgwp) && (!reterr)) {
    reterr = kPglRetWriteFail;
  }
  bigstack_reset(bigstack_mark);
  return reterr;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
