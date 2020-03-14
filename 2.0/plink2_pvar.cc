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

#include "plink2_pvar.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// this used to employ a backward-growing linked list, with variable-size
// elements, but demultiplexing was relatively expensive.  now we allocate
// size-64k pos[], allele_idxs[], ids[], cms[], etc. blocks, and just memcpy
// those chunks at the end.  (cms[] is lazy-initialized.)
//
// probable todo: try making this more parallel.  could e.g. repeatedly peek at
// end of ReadLineStream buffer, count # of '\n's, and divvy them up between
// worker threads.
CONSTI32(kLoadPvarBlockSize, 65536);
static_assert(!(kLoadPvarBlockSize & (kLoadPvarBlockSize - 1)), "kLoadPvarBlockSize must be a power of 2.");
static_assert(kLoadPvarBlockSize >= (kMaxMediumLine / 8), "kLoadPvarBlockSize cannot be smaller than kMaxMediumLine / 8.");


static_assert(!kChrOffsetX, "ReadChrsetHeaderLine() assumes kChrOffsetX == 0.");
static_assert(kChrOffsetY == 1, "ReadChrsetHeaderLine() assumes kChrOffsetY == 1.");
static_assert(kChrOffsetPAR1 == 4, "ReadChrsetHeaderLine() assumes kChrOffsetPAR1 == 4.");
PglErr ReadChrsetHeaderLine(const char* chrset_iter, const char* file_descrip, MiscFlags misc_flags, uintptr_t line_idx, ChrInfo* cip) {
  // chrset_iter is expected to point to first character after
  // "##chrSet=<".
  PglErr reterr = kPglRetSuccess;
  {
    uint32_t cmdline_autosome_ct = 0;
    uint32_t cmdline_haploid = 0;
    STD_ARRAY_DECL(uint32_t, kChrOffsetCt, cmdline_xymt_codes);
    if (cip->chrset_source == kChrsetSourceCmdline) {
      if (misc_flags & kfMiscChrOverrideCmdline) {
        goto ReadChrsetHeaderLine_ret_1;
      }
      if (!(misc_flags & kfMiscChrOverrideFile)) {
        // save off info we need for consistency check
        cmdline_autosome_ct = cip->autosome_ct;
        cmdline_haploid = cip->haploid_mask[0] & 1;
        STD_ARRAY_COPY(cip->xymt_codes, kChrOffsetCt, cmdline_xymt_codes);
      }
      ZeroWArr(kChrMaskWords, cip->haploid_mask);
    }
    for (uint32_t uii = 0; uii != kChrOffsetCt; ++uii) {
      cip->xymt_codes[uii] = UINT32_MAXM1;
    }
    if (StrStartsWithUnsafe(chrset_iter, "haploidAutosomeCt=")) {
      uint32_t explicit_haploid_ct;
      if (unlikely(ScanPosintCapped(&(chrset_iter[strlen("haploidAutosomeCt=")]), kMaxChrTextnum, &explicit_haploid_ct))) {
        snprintf(g_logbuf, kLogbufSize, "Error: Header line %" PRIuPTR " of %s has an invalid ##chrSet haploid count (max %u).\n", line_idx, file_descrip, kMaxChrTextnum);
        goto ReadChrsetHeaderLine_ret_MALFORMED_INPUT_WW;
      }
      // could verify that X, Y, etc. are not present?
      if (cmdline_autosome_ct) {
        if (unlikely(!cmdline_haploid)) {
          snprintf(g_logbuf, kLogbufSize, "Error: Header line %" PRIuPTR " of %s specifies a haploid genome, while a diploid genome was specified on the command line.\n", line_idx, file_descrip);
          goto ReadChrsetHeaderLine_ret_INCONSISTENT_INPUT_WW;
        }
        if (unlikely(explicit_haploid_ct != cmdline_autosome_ct)) {
          snprintf(g_logbuf, kLogbufSize, "Error: Header line %" PRIuPTR " of %s specifies %u autosome%s, while the command line specified %u.\n", line_idx, file_descrip, explicit_haploid_ct, (explicit_haploid_ct == 1)? "" : "s", cmdline_autosome_ct);
          goto ReadChrsetHeaderLine_ret_INCONSISTENT_INPUT_WW;
        }
      }
      cip->autosome_ct = explicit_haploid_ct;
      SetAllBits(explicit_haploid_ct + 1, cip->haploid_mask);
    } else {
      if (unlikely(!StrStartsWithUnsafe(chrset_iter, "autosomePairCt="))) {
        snprintf(g_logbuf, kLogbufSize, "Error: Header line %" PRIuPTR " of %s does not have expected ##chrSet format.\n", line_idx, file_descrip);
        goto ReadChrsetHeaderLine_ret_MALFORMED_INPUT_WW;
      }
      chrset_iter = &(chrset_iter[strlen("autosomePairCt=")]);
      uint32_t explicit_autosome_ct;
      if (unlikely(ScanmovPosintCapped(kMaxChrTextnum, &chrset_iter, &explicit_autosome_ct))) {
        snprintf(g_logbuf, kLogbufSize, "Error: Header line %" PRIuPTR " of %s has an invalid ##chrSet autosome count (max %u).\n", line_idx, file_descrip, kMaxChrTextnum);
        goto ReadChrsetHeaderLine_ret_MALFORMED_INPUT_WW;
      }
      cip->autosome_ct = explicit_autosome_ct;
      if (*chrset_iter != '>') {
        if (unlikely(*chrset_iter != ',')) {
          snprintf(g_logbuf, kLogbufSize, "Error: Header line %" PRIuPTR " of %s does not have expected ##chrSet format.\n", line_idx, file_descrip);
          goto ReadChrsetHeaderLine_ret_MALFORMED_INPUT_WW;
        }
        // this can theoretically be confused by e.g. a Description="..." field
        // containing commas not followed by spaces
        do {
          ++chrset_iter;
          // uppercase
          uint32_t first_char_ui = ctou32(*chrset_iter) & 0xdf;

          uint32_t second_char_ui = ctou32(chrset_iter[1]);
          // 44 is ',', 62 is '>'
          if ((second_char_ui == 44) || (second_char_ui == 62)) {
            if (first_char_ui == 77) {
              // M
              cip->xymt_codes[kChrOffsetMT] = explicit_autosome_ct + 1 + kChrOffsetMT;
            } else {
              first_char_ui -= 88;  // X = 0, Y = 1, everything else larger
              if (first_char_ui < 2) {
                cip->xymt_codes[first_char_ui] = explicit_autosome_ct + 1 + first_char_ui;
              }
            }
          } else {
            second_char_ui &= 0xdf;
            const uint32_t third_char_ui = ctou32(chrset_iter[2]);
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
              const uint32_t par_idx_m1 = ctou32(chrset_iter[3]) - '1';
              if ((par_idx_m1 < 2) && ((chrset_iter[4] == ',') || (chrset_iter[4] == '>'))) {
                cip->xymt_codes[kChrOffsetPAR1] = explicit_autosome_ct + 1 + kChrOffsetPAR1 + par_idx_m1;
              }
            }
          }
        } while (!strchrnul_n_mov(',', &chrset_iter));
      }
      if (cmdline_autosome_ct) {
        if (unlikely(cmdline_haploid)) {
          snprintf(g_logbuf, kLogbufSize, "Error: Header line %" PRIuPTR " of %s specifies a diploid genome, while a haploid genome was specified on the command line.\n", line_idx, file_descrip);
          goto ReadChrsetHeaderLine_ret_INCONSISTENT_INPUT_WW;
        }
        if (unlikely(explicit_autosome_ct != cmdline_autosome_ct)) {
          snprintf(g_logbuf, kLogbufSize, "Error: Header line %" PRIuPTR " of %s specifies %u autosome%s, while the command line specified %u.\n", line_idx, file_descrip, explicit_autosome_ct, (explicit_autosome_ct == 1)? "" : "s", cmdline_autosome_ct);
          goto ReadChrsetHeaderLine_ret_INCONSISTENT_INPUT_WW;
        }
        for (uint32_t xymt_idx = 0; xymt_idx != kChrOffsetPAR1; ++xymt_idx) {
          // it's okay if the command line doesn't explicitly exclude e.g. chrX
          // while for whatever reason it is excluded from ##chrSet; but the
          // reverse can create problems
          if (unlikely(IsI32Neg(cmdline_xymt_codes[xymt_idx]) && (!IsI32Neg(cip->xymt_codes[xymt_idx])))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Header line %" PRIuPTR " of %s specifies a chromosome set including %s, while the command line excludes it.\n", line_idx, file_descrip, g_xymt_log_names[xymt_idx]);
            goto ReadChrsetHeaderLine_ret_INCONSISTENT_INPUT_WW;
          }
        }
      }
    }
    cip->chrset_source = kChrsetSourceFile;
  }
  while (0) {
  ReadChrsetHeaderLine_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetMalformedInput;
    break;
  ReadChrsetHeaderLine_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetInconsistentInput;
    break;
  }
 ReadChrsetHeaderLine_ret_1:
  return reterr;
}

void VaridTemplateInit(const char* varid_template_str, const char* missing_id_match, char* chr_output_name_buf, uint32_t new_id_max_allele_slen, uint32_t overflow_substitute_blen, VaridTemplate* vtp) {
  // template string was previously validated
  // varid_template is only input, everything else is output values
  const char* varid_template_str_iter = varid_template_str;
  uint32_t template_insert_ct = 0;
  uint32_t template_base_len = 0;
  uint32_t alleles_needed = 0;  // bit 0 = ref, bit 1 = alt, bit 2 = ascii sort
  vtp->chr_output_name_buf = chr_output_name_buf;
  vtp->segs[0] = varid_template_str_iter;
  vtp->chr_slen = 0;
  for (unsigned char ucc = *varid_template_str_iter; ucc; ucc = *(++varid_template_str_iter)) {
    if (ucc > '@') {
      continue;
    }
    uint32_t seg_len;
    uint32_t insert_type;
    if (ucc == '@') {
      seg_len = varid_template_str_iter - vtp->segs[template_insert_ct];
      insert_type = 0;
    } else if (ucc == '#') {
      seg_len = varid_template_str_iter - vtp->segs[template_insert_ct];
      insert_type = 1;
    } else if (ucc != '$') {
      continue;
    } else {
      seg_len = varid_template_str_iter - vtp->segs[template_insert_ct];
      const uint32_t uii = ctou32(*(++varid_template_str_iter));
      if (uii <= '2') {
        alleles_needed += 2;  // this happens twice
        insert_type = uii - 48;  // '1' -> type 2, '2' -> type 3
      } else {
        // 'r' -> type 2, 'a' -> type 3
        insert_type = 1 + ((uii & 0xdf) == 'A');
      }
      alleles_needed += insert_type;
      ++insert_type;
    }
    vtp->seg_lens[template_insert_ct] = seg_len;
    template_base_len += seg_len;
    vtp->insert_types[template_insert_ct++] = insert_type;
    vtp->segs[template_insert_ct] = &(varid_template_str_iter[1]);
  }
  const uint32_t seg_len = varid_template_str_iter - vtp->segs[template_insert_ct];
  vtp->seg_lens[template_insert_ct] = seg_len;
  vtp->insert_ct = template_insert_ct;
  vtp->base_len = template_base_len + seg_len;
  vtp->alleles_needed = alleles_needed;
  vtp->new_id_max_allele_slen = new_id_max_allele_slen;
  vtp->overflow_substitute_blen = overflow_substitute_blen;
  vtp->missing_id_match = missing_id_match;
}

// alt1_end currently allowed to be nullptr in biallelic case
BoolErr VaridTemplateApply(unsigned char* tmp_alloc_base, const VaridTemplate* vtp, const char* ref_start, const char* alt1_start, uint32_t cur_bp, uint32_t ref_token_slen, uint32_t extra_alt_ct, uint32_t alt_token_slen, unsigned char** tmp_alloc_endp, uintptr_t* new_id_allele_len_overflowp, uint32_t* id_slen_ptr) {
  uint32_t insert_slens[4];
  const uint32_t alleles_needed = vtp->alleles_needed;
  const uint32_t new_id_max_allele_slen = vtp->new_id_max_allele_slen;
  uint32_t id_slen = UintSlen(cur_bp);
  insert_slens[0] = vtp->chr_slen;
  insert_slens[1] = id_slen;
  id_slen += vtp->base_len;
  uint32_t ref_slen = 0;
  uint32_t cur_overflow = 0;
  const char* tmp_allele_ptrs[2];
  if (alleles_needed & 1) {
    ref_slen = ref_token_slen;
    if (ref_slen > new_id_max_allele_slen) {
      ref_slen = new_id_max_allele_slen;
      cur_overflow = 1;
    }
    insert_slens[2] = ref_slen;
    id_slen += ref_slen;
    tmp_allele_ptrs[0] = ref_start;
  }
  if (alleles_needed > 1) {
    uint32_t alt1_slen;
    if (!extra_alt_ct) {
      alt1_slen = alt_token_slen;
    } else {
      alt1_slen = AdvToDelim(alt1_start, ',') - alt1_start;
    }
    if (alt1_slen > new_id_max_allele_slen) {
      alt1_slen = new_id_max_allele_slen;
      ++cur_overflow;
    }
    id_slen += alt1_slen;
    if (alleles_needed <= 3) {
    VaridTemplateApply_keep_allele_ascii_order:
      insert_slens[3] = alt1_slen;
      tmp_allele_ptrs[1] = alt1_start;
    } else {
      uint32_t smaller_slen = alt1_slen;
      const int32_t ref_slen_geq = (ref_slen >= alt1_slen);
      if (!ref_slen_geq) {
        smaller_slen = ref_slen;
      }
      int32_t memcmp_result = Memcmp(ref_start, alt1_start, smaller_slen);
      if (!memcmp_result) {
        memcmp_result = ref_slen_geq;
      }
      if (memcmp_result <= 0) {
        goto VaridTemplateApply_keep_allele_ascii_order;
      }
      insert_slens[3] = ref_slen;
      tmp_allele_ptrs[1] = tmp_allele_ptrs[0];
      insert_slens[2] = alt1_slen;
      tmp_allele_ptrs[0] = alt1_start;
    }
  }
  const uint32_t overflow_substitute_blen = vtp->overflow_substitute_blen;
  if (overflow_substitute_blen && cur_overflow) {
    // er, do we really want to make a separate allocation here?  see if we can
    // remove this.  (that would impose more constraints on downstream
    // variant-ID-changing functions, though.)
    if (PtrWSubCk(tmp_alloc_base, overflow_substitute_blen, tmp_alloc_endp)) {
      return 1;
    }
    memcpy(*tmp_alloc_endp, vtp->missing_id_match, overflow_substitute_blen);
    id_slen = 0;
    cur_overflow = 1;
  } else {
    if (PtrWSubCk(tmp_alloc_base, id_slen + 1, tmp_alloc_endp)) {
      return 1;
    }
    char* id_iter = R_CAST(char*, *tmp_alloc_endp);
    const uint32_t insert_ct = vtp->insert_ct;
    char* insert_ptrs[4];

    // maybe-uninitialized warnings
    insert_ptrs[0] = nullptr;
    insert_ptrs[1] = nullptr;
    insert_ptrs[2] = nullptr;

    for (uint32_t insert_idx = 0; insert_idx != insert_ct; ++insert_idx) {
      id_iter = memcpya(id_iter, vtp->segs[insert_idx], vtp->seg_lens[insert_idx]);
      const uint32_t cur_insert_type = vtp->insert_types[insert_idx];
      insert_ptrs[cur_insert_type] = id_iter;
      id_iter = &(id_iter[insert_slens[cur_insert_type]]);
    }
    memcpyx(id_iter, vtp->segs[insert_ct], vtp->seg_lens[insert_ct], '\0');

    memcpy(insert_ptrs[0], vtp->chr_output_name_buf, insert_slens[0]);
    u32toa(cur_bp, insert_ptrs[1]);
    // insert_ct currently guaranteed to be >= 2 since @ and # required
    for (uint32_t insert_type_idx = 2; insert_type_idx != insert_ct; ++insert_type_idx) {
      memcpy(insert_ptrs[insert_type_idx], tmp_allele_ptrs[insert_type_idx - 2], insert_slens[insert_type_idx]);
    }
  }
  *id_slen_ptr = id_slen;
  *new_id_allele_len_overflowp += cur_overflow;
  return 0;
}

// Exported variant of VaridTemplateApply() which appends to a buffer.
// Probable todo: pull out the common parts of the functions.
char* VaridTemplateWrite(const VaridTemplate* vtp, const char* ref_start, const char* alt1_start, uint32_t cur_bp, uint32_t ref_token_slen, uint32_t extra_alt_ct, uint32_t alt_token_slen, char* dst) {
  uint32_t insert_slens[4];
  const uint32_t alleles_needed = vtp->alleles_needed;
  const uint32_t new_id_max_allele_slen = vtp->new_id_max_allele_slen;
  uint32_t id_slen = UintSlen(cur_bp);
  insert_slens[0] = vtp->chr_slen;
  insert_slens[1] = id_slen;
  id_slen += vtp->base_len;
  uint32_t ref_slen = 0;
  uint32_t cur_overflow = 0;
  const char* tmp_allele_ptrs[2];
  if (alleles_needed & 1) {
    ref_slen = ref_token_slen;
    if (ref_slen > new_id_max_allele_slen) {
      ref_slen = new_id_max_allele_slen;
      cur_overflow = 1;
    }
    insert_slens[2] = ref_slen;
    id_slen += ref_slen;
    tmp_allele_ptrs[0] = ref_start;
  }
  if (alleles_needed > 1) {
    uint32_t alt1_slen;
    if (!extra_alt_ct) {
      alt1_slen = alt_token_slen;
    } else {
      alt1_slen = AdvToDelim(alt1_start, ',') - alt1_start;
    }
    if (alt1_slen > new_id_max_allele_slen) {
      alt1_slen = new_id_max_allele_slen;
      cur_overflow = 1;
    }
    id_slen += alt1_slen;
    if (alleles_needed <= 3) {
    VaridTemplateWrite_keep_allele_ascii_order:
      insert_slens[3] = alt1_slen;
      tmp_allele_ptrs[1] = alt1_start;
    } else {
      uint32_t smaller_slen = alt1_slen;
      const int32_t ref_slen_geq = (ref_slen >= alt1_slen);
      if (!ref_slen_geq) {
        smaller_slen = ref_slen;
      }
      int32_t memcmp_result = Memcmp(ref_start, alt1_start, smaller_slen);
      if (!memcmp_result) {
        memcmp_result = ref_slen_geq;
      }
      if (memcmp_result <= 0) {
        goto VaridTemplateWrite_keep_allele_ascii_order;
      }
      insert_slens[3] = ref_slen;
      tmp_allele_ptrs[1] = tmp_allele_ptrs[0];
      insert_slens[2] = alt1_slen;
      tmp_allele_ptrs[0] = alt1_start;
    }
  }
  const uint32_t overflow_substitute_blen = vtp->overflow_substitute_blen;
  if (overflow_substitute_blen && cur_overflow) {
    return memcpya(dst, vtp->missing_id_match, overflow_substitute_blen);
  }
  char* id_iter = dst;
  const uint32_t insert_ct = vtp->insert_ct;
  char* insert_ptrs[4];

  // maybe-uninitialized warnings
  insert_ptrs[0] = nullptr;
  insert_ptrs[1] = nullptr;
  insert_ptrs[2] = nullptr;

  for (uint32_t insert_idx = 0; insert_idx != insert_ct; ++insert_idx) {
    id_iter = memcpya(id_iter, vtp->segs[insert_idx], vtp->seg_lens[insert_idx]);
    const uint32_t cur_insert_type = vtp->insert_types[insert_idx];
    insert_ptrs[cur_insert_type] = id_iter;
    id_iter = &(id_iter[insert_slens[cur_insert_type]]);
  }
  char* id_end = memcpya(id_iter, vtp->segs[insert_ct], vtp->seg_lens[insert_ct]);

  memcpy(insert_ptrs[0], vtp->chr_output_name_buf, insert_slens[0]);
  u32toa(cur_bp, insert_ptrs[1]);
  // insert_ct currently guaranteed to be >= 2 since @ and # required
  for (uint32_t insert_type_idx = 2; insert_type_idx != insert_ct; ++insert_type_idx) {
    memcpy(insert_ptrs[insert_type_idx], tmp_allele_ptrs[insert_type_idx - 2], insert_slens[insert_type_idx]);
  }
  return id_end;
}

uint32_t VaridWorstCaseSlen(const VaridTemplate* vtp, uint32_t max_chr_slen, uint32_t max_allele_slen) {
  // +10 for base-pair coordinate
  return (max_allele_slen * vtp->alleles_needed + vtp->base_len + max_chr_slen + 10);
}

void BackfillChrIdxs(const ChrInfo* cip, uint32_t chrs_encountered_m1, uint32_t offset, uint32_t end_vidx, ChrIdx* chr_idxs) {
  for (uint32_t chr_fo_idx = chrs_encountered_m1; ; --chr_fo_idx) {
    uint32_t start_vidx = cip->chr_fo_vidx_start[chr_fo_idx];
    if (start_vidx < offset) {
      start_vidx = offset;
    }
    ChrIdx* chr_idxs_write_base = &(chr_idxs[start_vidx - offset]);
    const uint32_t vidx_ct = end_vidx - start_vidx;
    const ChrIdx cur_chr_idx = cip->chr_file_order[chr_fo_idx];
    for (uint32_t uii = 0; uii != vidx_ct; ++uii) {
      chr_idxs_write_base[uii] = cur_chr_idx;
    }
    if (start_vidx == offset) {
      return;
    }
    end_vidx = start_vidx;
  }
}

char* PrInInfoToken(uint32_t info_slen, char* info_token) {
  if (memequal_k(info_token, "PR", 2) && ((info_slen == 2) || (info_token[2] == ';'))) {
    return info_token;
  }
  if (memequal_k(&(info_token[S_CAST(int32_t, info_slen) - 3]), ";PR", 3)) {
    return &(info_token[info_slen - 2]);
  }
  info_token[info_slen] = '\0';
  char* first_info_end = strchr(info_token, ';');
  if (!first_info_end) {
    return nullptr;
  }
  char* pr_prestart = strstr(first_info_end, ";PR;");
  // bugfix (27 Sep 2019): had this backward
  return pr_prestart? (&(pr_prestart[1])) : nullptr;
}

typedef struct InfoExistStruct {
  NONCOPYABLE(InfoExistStruct);
  char* prekeys;
  uint32_t key_ct;
  uint32_t key_slens[];
} InfoExist;

PglErr InfoExistInit(const unsigned char* arena_end, const char* require_info_flattened, const char* flagname_p, unsigned char** arena_base_ptr, InfoExist** existpp) {
  const char* read_iter = require_info_flattened;
  uintptr_t prekeys_blen = 0;
  do {
    const char* key_end_or_invalid = strchrnul2(read_iter, ';', '=');
    if (unlikely(*key_end_or_invalid)) {
      if (*key_end_or_invalid == ';') {
        logerrprintfww("Error: Invalid --%s key '%s' (semicolon prohibited).\n", flagname_p, read_iter);
      } else {
        logerrprintfww("Error: Invalid --%s key '%s' ('=' prohibited).\n", flagname_p, read_iter);
      }
      return kPglRetInvalidCmdline;
    }
    ++prekeys_blen;
    read_iter = &(key_end_or_invalid[1]);
  } while (*read_iter);
  prekeys_blen += read_iter - require_info_flattened;
  const uint32_t key_ct = prekeys_blen - S_CAST(uintptr_t, read_iter - require_info_flattened);
  const uintptr_t bytes_used = offsetof(InfoExist, key_slens) + key_ct * sizeof(int32_t) + prekeys_blen;
  const uintptr_t cur_alloc = RoundUpPow2(bytes_used, kCacheline);
  if (unlikely(S_CAST(uintptr_t, arena_end - (*arena_base_ptr)) < cur_alloc)) {
    return kPglRetNomem;
  }
  *existpp = R_CAST(InfoExist*, *arena_base_ptr);
  (*existpp)->prekeys = R_CAST(char*, &((*arena_base_ptr)[bytes_used - prekeys_blen]));
  (*arena_base_ptr) += cur_alloc;
  (*existpp)->key_ct = key_ct;
  read_iter = require_info_flattened;
  char* write_iter = (*existpp)->prekeys;
  for (uint32_t kidx = 0; kidx != key_ct; ++kidx) {
    *write_iter++ = ';';
    const uintptr_t slen = strlen(read_iter);
    (*existpp)->key_slens[kidx] = slen;
    const uintptr_t blen = slen + 1;
    write_iter = memcpya(write_iter, read_iter, blen);
    read_iter = &(read_iter[blen]);
  }
  return kPglRetSuccess;
}

uint32_t InfoExistCheck(const char* info_token, const InfoExist* existp) {
  const uint32_t key_ct = existp->key_ct;
  const char* prekeys_iter = existp->prekeys;
  const uint32_t* key_slens = existp->key_slens;
  for (uint32_t kidx = 0; kidx != key_ct; ++kidx) {
    const uint32_t key_slen = key_slens[kidx];
    const char* possible_hit;
    // similar logic to hardcoded PR, except key can also be followed by '=',
    // and if it is there can't be a lone '.' after the equals
    if (memequal(info_token, &(prekeys_iter[1]), key_slen)) {
      possible_hit = &(info_token[key_slen]);
    } else {
      possible_hit = strstr(info_token, prekeys_iter);
      if (!possible_hit) {
        return 0;
      }
      possible_hit = &(possible_hit[key_slen + 1]);
    }
    while (1) {
      const char cc = *possible_hit;
      if ((!cc) || (cc == ';')) {
        break;
      }
      if (cc == '=') {
        if ((possible_hit[1] != '.') || (possible_hit[2] && (possible_hit[2] != ';'))) {
          break;
        }
        return 0;
      }
      possible_hit = strstr(possible_hit, prekeys_iter);
      if (!possible_hit) {
        return 0;
      }
      possible_hit = &(possible_hit[key_slen + 1]);
    }
    prekeys_iter = &(prekeys_iter[key_slen + 2]);
  }
  return 1;
}

uint32_t InfoNonexistCheck(const char* info_token, const InfoExist* nonexistp) {
  const uint32_t key_ct = nonexistp->key_ct;
  const char* prekeys_iter = nonexistp->prekeys;
  const uint32_t* key_slens = nonexistp->key_slens;
  uint32_t key_slen = 0;
  for (uint32_t kidx = 0; kidx != key_ct; ++kidx, prekeys_iter = &(prekeys_iter[key_slen + 2])) {
    key_slen = key_slens[kidx];
    const char* possible_hit;
    if (memequal(info_token, &(prekeys_iter[1]), key_slen)) {
      possible_hit = &(info_token[key_slen]);
    } else {
      possible_hit = strstr(info_token, prekeys_iter);
      if (!possible_hit) {
        continue;
      }
      possible_hit = &(possible_hit[key_slen + 1]);
    }
    while (1) {
      const char cc = *possible_hit;
      if ((!cc) || (cc == ';')) {
        return 0;
      }
      if (cc == '=') {
        if ((possible_hit[1] != '.') || (possible_hit[2] && (possible_hit[2] != ';'))) {
          return 0;
        }
        break;
      }
      possible_hit = strstr(possible_hit, prekeys_iter);
      if (!possible_hit) {
        break;
      }
      possible_hit = &(possible_hit[key_slen + 1]);
    }
  }
  return 1;
}

typedef struct InfoFilterStruct {
  NONCOPYABLE(InfoFilterStruct);
  char* prekey;
  const char* val_str;
  uint32_t key_slen;
  uint32_t val_slen;
  CmpBinaryOp binary_op;
  double val;
} InfoFilter;

PglErr InfoFilterInit(const unsigned char* arena_end, const CmpExpr* filter_exprp, const char* flagname_p, unsigned char** arena_base_ptr, InfoFilter* filterp) {
  const char* pheno_name = filter_exprp->pheno_name;
  const char* pheno_name_end_or_invalid = strchrnul3(pheno_name, ';', '=', ',');
  if (unlikely(*pheno_name_end_or_invalid)) {
    if (*pheno_name_end_or_invalid == ';') {
      logerrprintfww("Error: Invalid --%s key '%s' (semicolon prohibited).\n", flagname_p, pheno_name);
    } else if (*pheno_name_end_or_invalid == '=') {
      logerrprintfww("Error: Invalid --%s key '%s' ('=' prohibited).\n", flagname_p, pheno_name);
    } else {
      logerrprintfww("Error: Invalid --%s key '%s' (comma prohibited).\n", flagname_p, pheno_name);
    }
    return kPglRetInvalidCmdline;
  }
  uint32_t key_slen = pheno_name_end_or_invalid - pheno_name;
  const uintptr_t cur_alloc = RoundUpPow2(3 + key_slen, kCacheline);
  if (unlikely(S_CAST(uintptr_t, arena_end - (*arena_base_ptr)) < cur_alloc)) {
    return kPglRetNomem;
  }
  filterp->prekey = R_CAST(char*, *arena_base_ptr);
  (*arena_base_ptr) += cur_alloc;
  filterp->prekey[0] = ';';
  char* key_iter = memcpya(&(filterp->prekey[1]), pheno_name, key_slen);
  strcpy_k(key_iter, "=");
  ++key_slen;
  filterp->key_slen = key_slen;
  filterp->val_str = nullptr;
  filterp->val_slen = 0;
  const CmpBinaryOp binary_op = filter_exprp->binary_op;
  filterp->binary_op = binary_op;
  // shouldn't need to initialize val

  const char* val_str = &(pheno_name[key_slen]);
  // bugfix (14 Dec 2017): INFO string constants are not guaranteed to start in
  // a letter.  Only interpret the value as a number if ScanadvDouble()
  // consumes the entire value string.
  const char* val_num_end = ScanadvDouble(val_str, &filterp->val);
  if (val_num_end && (!val_num_end[0])) {
    return kPglRetSuccess;
  }
  if (unlikely((binary_op != kCmpOperatorNoteq) && (binary_op != kCmpOperatorEq))) {
    logerrprintfww("Error: Invalid --%s value '%s' (finite number expected).\n", flagname_p, val_str);
    return kPglRetInvalidCmdline;
  }
  filterp->val_str = val_str;
  const char* val_str_end_or_invalid = strchrnul2(val_str, ';', '=');
  if (unlikely(*val_str_end_or_invalid)) {
    if (*val_str_end_or_invalid == ';') {
      logerrprintfww("Error: Invalid --%s value '%s' (semicolon prohibited).\n", flagname_p, val_str);
    } else {
      logerrprintfww("Error: Invalid --%s value '%s' ('=' prohibited).\n", flagname_p, val_str);
    }
    return kPglRetInvalidCmdline;
  }
  filterp->val_slen = val_str_end_or_invalid - val_str;
  return kPglRetSuccess;
}

uint32_t InfoConditionSatisfied(const char* info_token, const InfoFilter* filterp) {
  const char* prekey = filterp->prekey;
  const uint32_t key_slen = filterp->key_slen;
  const CmpBinaryOp binary_op = filterp->binary_op;
  const char* possible_hit;
  // search for "[key]=" at start or ";[key]=" later; key_slen includes the
  //   trailing =
  if (memequal(info_token, &(prekey[1]), key_slen)) {
    possible_hit = &(info_token[key_slen]);
  } else {
    possible_hit = strstr(info_token, prekey);
    if (!possible_hit) {
      return (binary_op == kCmpOperatorNoteq);
    }
    possible_hit = &(possible_hit[key_slen + 1]);
  }
  if (filterp->val_str) {
    const uint32_t val_slen = filterp->val_slen;
    const uint32_t mismatch = ((!memequal(possible_hit, filterp->val_str, val_slen)) || (possible_hit[val_slen] && (possible_hit[val_slen] != ';')));
    return mismatch ^ (binary_op != kCmpOperatorNoteq);
  }
  // bugfix (3 Jan 2020): semicolon (or \0) terminator expected, can't use
  // ScantokDouble
  double dxx;
  const char* scan_end = ScanadvDouble(possible_hit, &dxx);
  if ((!scan_end) || ((*scan_end != ';') && (*scan_end))) {
    return (binary_op == kCmpOperatorNoteq);
  }
  const double val = filterp->val;
  switch (binary_op) {
  case kCmpOperatorNoteq:
    return (dxx != val);
  case kCmpOperatorLe:
    return (dxx < val);
  case kCmpOperatorLeq:
    return (dxx <= val);
  case kCmpOperatorGe:
    return (dxx > val);
  case kCmpOperatorGeq:
    return (dxx >= val);
  case kCmpOperatorEq:
    return (dxx == val);
  }
  // should be unreachable, but some gcc versions complain
  return 0;
}

PglErr SplitPar(const uint32_t* variant_bps, UnsortedVar vpos_sortstatus, uint32_t splitpar_bound1, uint32_t splitpar_bound2, uintptr_t* variant_include, uintptr_t* loaded_chr_mask, ChrInfo* cip, uint32_t* chrs_encountered_m1_ptr, uint32_t* exclude_ct_ptr) {
  const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
  if (IsI32Neg(x_code) || (!IsSet(loaded_chr_mask, x_code))) {
    logerrputs("Warning: --split-par had no effect (no X chromosome in dataset).\n");
    return kPglRetSuccess;
  }
  const uint32_t par1_code = cip->xymt_codes[kChrOffsetPAR1];
  const uint32_t par2_code = cip->xymt_codes[kChrOffsetPAR2];
  if (unlikely(IsI32Neg(par2_code))) {
    // may want to remove this restriction later
    logerrputs("Error: --split-par cannot currently be used with a custom chromosome set.\n");
    return kPglRetInvalidCmdline;
  }
  if (unlikely(IsSet(loaded_chr_mask, par1_code) || IsSet(loaded_chr_mask, par2_code))) {
    logerrputs("Error: --split-par cannot be used on a dataset which already contains a PAR1 or\nPAR2 region.\n");
    return kPglRetInvalidCmdline;
  }
  if (unlikely(vpos_sortstatus & kfUnsortedVarBp)) {
    logerrputs("Error: --split-par cannot be used with an unsorted .bim/.pvar file.\n");
    return kPglRetInvalidCmdline;
  }
  const uint32_t orig_xchr_fo_idx = cip->chr_idx_to_foidx[x_code];
  const uint32_t orig_x_start = cip->chr_fo_vidx_start[orig_xchr_fo_idx];
  const uint32_t orig_x_end = cip->chr_fo_vidx_start[orig_xchr_fo_idx + 1];
  const uint32_t par1_end = orig_x_start + CountSortedSmallerU32(&(variant_bps[orig_x_start]), orig_x_end - orig_x_start, splitpar_bound1 + 1);
  const uint32_t par2_start = par1_end + CountSortedSmallerU32(&(variant_bps[par1_end]), orig_x_end - par1_end, splitpar_bound2);
  uint32_t tot_codes_changed = (par1_end - orig_x_start) + (orig_x_end - par2_start);
  if (!tot_codes_changed) {
    logerrputs("Warning: --split-par had no effect (no X variants were in the PARs).\n");
    return kPglRetSuccess;
  }
  // one of the PARs, and/or the main chrX body, may be empty; that's not a big
  // deal
  *chrs_encountered_m1_ptr += 2;
  const uint32_t chrs_encountered_m1 = *chrs_encountered_m1_ptr;
  cip->chr_fo_vidx_start[chrs_encountered_m1 + 1] = cip->chr_fo_vidx_start[chrs_encountered_m1 - 1];
  for (uint32_t chr_fo_idx = chrs_encountered_m1 - 2; chr_fo_idx != orig_xchr_fo_idx; --chr_fo_idx) {
    cip->chr_fo_vidx_start[chr_fo_idx + 2] = cip->chr_fo_vidx_start[chr_fo_idx];
    const uint32_t cur_chr_idx = cip->chr_file_order[chr_fo_idx];
    cip->chr_file_order[chr_fo_idx + 2] = cur_chr_idx;
    cip->chr_idx_to_foidx[cur_chr_idx] = chr_fo_idx + 2;
  }
  cip->chr_fo_vidx_start[orig_xchr_fo_idx + 1] = par1_end;
  cip->chr_fo_vidx_start[orig_xchr_fo_idx + 2] = par2_start;
  cip->chr_file_order[orig_xchr_fo_idx] = par1_code;
  cip->chr_file_order[orig_xchr_fo_idx + 1] = x_code;
  cip->chr_file_order[orig_xchr_fo_idx + 2] = par2_code;
  cip->chr_idx_to_foidx[par1_code] = orig_xchr_fo_idx;
  cip->chr_idx_to_foidx[x_code] = orig_xchr_fo_idx + 1;
  cip->chr_idx_to_foidx[par2_code] = orig_xchr_fo_idx + 2;
  uintptr_t* chr_mask = cip->chr_mask;
  if (par1_end > orig_x_start) {
    if (!IsSet(chr_mask, par1_code)) {
      *exclude_ct_ptr += PopcountBitRange(variant_include, orig_x_start, par1_end);
      ClearBitsNz(orig_x_start, par1_end, variant_include);
    } else {
      SetBit(par1_code, loaded_chr_mask);
    }
  }
  if (par1_end == par2_start) {
    ClearBit(x_code, chr_mask);
  } else if (!IsSet(chr_mask, x_code)) {
    ClearBit(x_code, chr_mask);
    *exclude_ct_ptr += PopcountBitRange(variant_include, par1_end, par2_start);
    ClearBitsNz(par1_end, par2_start, variant_include);
  }
  if (par2_start < orig_x_end) {
    if (!IsSet(chr_mask, par2_code)) {
      *exclude_ct_ptr += PopcountBitRange(variant_include, par2_start, orig_x_end);
      ClearBitsNz(par2_start, orig_x_end, variant_include);
    } else {
      SetBit(par2_code, loaded_chr_mask);
    }
  }
  logprintf("--split-par: %u chromosome code%s changed.\n", tot_codes_changed, (tot_codes_changed == 1)? "" : "s");
  return kPglRetSuccess;
}

static_assert((!(kMaxIdSlen % kCacheline)), "LoadPvar() must be updated.");
PglErr LoadPvar(const char* pvarname, const char* var_filter_exceptions_flattened, const char* varid_template_str, const char* varid_multi_template_str, const char* varid_multi_nonsnp_template_str, const char* missing_varid_match, const char* require_info_flattened, const char* require_no_info_flattened, const CmpExpr* extract_if_info_exprp, const CmpExpr* exclude_if_info_exprp, MiscFlags misc_flags, PvarPsamFlags pvar_psam_flags, uint32_t xheader_needed, uint32_t qualfilter_needed, float var_min_qual, uint32_t splitpar_bound1, uint32_t splitpar_bound2, uint32_t new_variant_id_max_allele_slen, uint32_t snps_only, uint32_t split_chr_ok, uint32_t filter_min_allele_ct, uint32_t filter_max_allele_ct, uint32_t max_thread_ct, ChrInfo* cip, uint32_t* max_variant_id_slen_ptr, uint32_t* info_reload_slen_ptr, UnsortedVar* vpos_sortstatus_ptr, char** xheader_ptr, uintptr_t** variant_include_ptr, uint32_t** variant_bps_ptr, char*** variant_ids_ptr, uintptr_t** allele_idx_offsets_ptr, const char*** allele_storage_ptr, uintptr_t** qual_present_ptr, float** quals_ptr, uintptr_t** filter_present_ptr, uintptr_t** filter_npass_ptr, char*** filter_storage_ptr, uintptr_t** nonref_flags_ptr, double** variant_cms_ptr, ChrIdx** chr_idxs_ptr, uint32_t* raw_variant_ct_ptr, uint32_t* variant_ct_ptr, uint32_t* max_allele_ct_ptr, uint32_t* max_allele_slen_ptr, uintptr_t* xheader_blen_ptr, InfoFlags* info_flags_ptr, uint32_t* max_filter_slen_ptr) {
  // chr_info, max_variant_id_slen, and info_reload_slen are in/out; just
  // outparameters after them.  (Due to its large size in some VCFs, INFO is
  // not kept in memory for now.  This has a speed penalty, of course; maybe
  // it's worthwhile to conditionally load it later.)

  // allele_idx_offsets currently assumed to be initialized to nullptr

  // should handle raw_variant_ct == 0 properly

  // possible todo: optionally skip allele code loading
  // probable todo: load INFO:END.  (does this allow the CNV module to be
  //   unified with the rest of the program?)  but this will probably wait
  //   until I need to analyze some sort of CNV data, and that day keeps
  //   getting postponed... for now, the BCF exporter performs its own parsing
  //   of INFO:END so that it can fill each variant's rlen field correctly.
  // possible todo: require FILTER to only contain values declared in header,
  //   and modify its storage accordingly?  (pointless for now, but worthwhile
  //   to keep an eye on what typical VCF files look like.)

  // Workspace is used as follows:
  // |--header, allele_storage->----|--other return arrays---|--linebuf--|-
  //                                                        1/4
  //
  // -temp-->----|----<- filter failures, variant IDs, long alleles--|
  //                                                                end
  // I.e. on successful return, both bigstack_base and bigstack_end will move.
  // This is designed to be called near the start of a program, at a time when
  // no large temporary buffer is needed.
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  uintptr_t line_idx = 0;
  uint32_t max_extra_alt_ct = 0;
  uint32_t max_allele_slen = 1;
  PglErr reterr = kPglRetSuccess;
  TextStream pvar_txs;
  PreinitTextStream(&pvar_txs);
  {
    const uintptr_t quarter_left = RoundDownPow2(bigstack_left() / 4, kCacheline);
    uint32_t max_line_blen;
    if (unlikely(StandardizeMaxLineBlenEx(quarter_left, kLoadPvarBlockSize * 2 * sizeof(intptr_t), &max_line_blen))) {
      goto LoadPvar_ret_NOMEM;
    }
    // bugfix (29 Jun 2018): don't cap allele_storage space at 2 GB, otherwise
    // we're limited to 134M variants
    unsigned char* rlstream_start = &(bigstack_mark[quarter_left]);
    g_bigstack_base = rlstream_start;
    uint32_t decompress_thread_ct = max_thread_ct - 1;
    // 3 still seems best on a heavily multicore Linux test machine
    if (decompress_thread_ct > 3) {
      decompress_thread_ct = 3;
    } else if (!decompress_thread_ct) {
      decompress_thread_ct = 1;
    }
    reterr = InitTextStream(pvarname, max_line_blen, decompress_thread_ct, &pvar_txs);
    if (unlikely(reterr)) {
      goto LoadPvar_ret_TSTREAM_FAIL;
    }
    unsigned char* tmp_alloc_base = g_bigstack_base;
    g_bigstack_base = bigstack_mark;

    char* xheader_end = ((pvar_psam_flags & (kfPvarColXheader | kfPvarColVcfheader)) || xheader_needed)? R_CAST(char*, bigstack_mark) : nullptr;
    uint32_t chrset_present = 0;
    uint32_t info_pr_present = 0;
    uint32_t info_pr_nonflag_present = 0;
    uint32_t info_nonpr_present = 0;
    char* line_iter = TextLineEnd(&pvar_txs);
    char* line_start;
    while (1) {
      ++line_idx;
      if (!TextGetUnsafe2(&pvar_txs, &line_iter)) {
        if (unlikely(TextStreamErrcode2(&pvar_txs, &reterr))) {
          goto LoadPvar_ret_TSTREAM_FAIL;
        }
        line_start = K_CAST(char*, &(g_one_char_strs[0]));
        --line_iter;
        break;
      }
      line_start = line_iter;
      if ((*line_start != '#') || tokequal_k(line_start, "#CHROM")) {
        break;
      }
      if (StrStartsWithUnsafe(line_start, "##INFO=<ID=PR,Number=")) {
        if (unlikely(info_pr_present || info_pr_nonflag_present)) {
          snprintf(g_logbuf, kLogbufSize, "Error: Duplicate INFO:PR header line in %s.\n", pvarname);
          goto LoadPvar_ret_MALFORMED_INPUT_WW;
        }
        info_pr_nonflag_present = !StrStartsWithUnsafe(&(line_start[strlen("##INFO=<ID=PR,Number=")]), "0,Type=Flag,Description=");
        info_pr_present = 1 - info_pr_nonflag_present;
        if (info_pr_nonflag_present) {
          logerrprintfww("Warning: Header line %" PRIuPTR " of %s has an unexpected definition of INFO:PR. This interferes with a few merge and liftover operations.\n", line_idx, pvarname);
        }
      } else if ((!info_nonpr_present) && StrStartsWithUnsafe(line_start, "##INFO=<ID=")) {
        info_nonpr_present = 1;
      }
      if (StrStartsWithUnsafe(line_start, "##chrSet=<")) {
        if (unlikely(chrset_present)) {
          snprintf(g_logbuf, kLogbufSize, "Error: Multiple ##chrSet header lines in %s.\n", pvarname);
          goto LoadPvar_ret_MALFORMED_INPUT_WW;
        }
        chrset_present = 1;
        const uint32_t cmdline_chrset = (cip->chrset_source == kChrsetSourceCmdline) && (!(misc_flags & kfMiscChrOverrideFile));
        reterr = ReadChrsetHeaderLine(&(line_start[strlen("##chrSet=<")]), pvarname, misc_flags, line_idx, cip);
        if (unlikely(reterr)) {
          goto LoadPvar_ret_1;
        }
        if (!cmdline_chrset) {
          const uint32_t autosome_ct = cip->autosome_ct;
          if (cip->haploid_mask[0] & 1) {
            logprintf("chrSet header line: %u autosome%s (haploid).\n", autosome_ct, (autosome_ct == 1)? "" : "s");
          } else {
            logprintf("chrSet header line: %u autosome pair%s.\n", autosome_ct, (autosome_ct == 1)? "" : "s");
          }
        }
      } else if (xheader_end) {
        // if the "pvar file" was actually a VCF, suppress the same lines we'd
        // suppress when importing with --vcf.
        if ((!StrStartsWithUnsafe(line_start, "##fileformat=")) &&
            (!StrStartsWithUnsafe(line_start, "##fileDate=")) &&
            (!StrStartsWithUnsafe(line_start, "##source=")) &&
            (!StrStartsWithUnsafe(line_start, "##FORMAT="))) {
          char* line_end = AdvToDelim(line_start, '\n');
          uint32_t line_slen = line_end - line_start;
          if (line_start[line_slen - 1] == '\r') {
            --line_slen;
          }
          if (unlikely(S_CAST(uintptr_t, R_CAST(char*, rlstream_start) - xheader_end) < line_slen + 2)) {
            goto LoadPvar_ret_NOMEM;
          }
          xheader_end = memcpya(xheader_end, line_start, line_slen);
          AppendBinaryEoln(&xheader_end);
          line_iter = line_end;
        }
      }
      line_iter = AdvPastDelim(line_iter, '\n');
    }
    if (xheader_end) {
      *xheader_ptr = R_CAST(char*, bigstack_mark);
      *xheader_blen_ptr = xheader_end - (*xheader_ptr);
      BigstackBaseSet(xheader_end);
    }
    FinalizeChrset(misc_flags, cip);
    const char** allele_storage = R_CAST(const char**, g_bigstack_base);
    const char** allele_storage_iter = allele_storage;

    uint32_t col_skips[8];
    uint32_t col_types[8];
    uint32_t relevant_postchr_col_ct = 5;
    uint32_t alt_col_idx = 4;
    uint32_t load_qual_col = 0;
    uint32_t load_filter_col = 0;
    uint32_t info_col_present = 0;
    uint32_t cm_col_present = 0;
    if (line_start[0] == '#') {
      *info_flags_ptr = S_CAST(InfoFlags, (info_pr_present * kfInfoPrFlagPresent) | (info_pr_nonflag_present * kfInfoPrNonflagPresent) | (info_nonpr_present * kfInfoNonprPresent));
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
      const char* token_end = &(line_start[6]);
      uint32_t found_header_bitset = 0;
      relevant_postchr_col_ct = 0;
      const char* linebuf_iter;
      for (uint32_t col_idx = 1; ; ++col_idx) {
        linebuf_iter = FirstNonTspace(token_end);
        if (IsEolnKns(*linebuf_iter)) {
          break;
        }
        token_end = CurTokenEnd(linebuf_iter);
        const uint32_t token_slen = token_end - linebuf_iter;
        uint32_t cur_col_type;
        if (token_slen <= 3) {
          if (token_slen == 3) {
            if (memequal_k(linebuf_iter, "POS", 3)) {
              cur_col_type = 0;
            } else if (memequal_k(linebuf_iter, "REF", 3)) {
              cur_col_type = 2;
            } else if (memequal_k(linebuf_iter, "ALT", 3)) {
              cur_col_type = 3;
              alt_col_idx = col_idx;
            } else {
              continue;
            }
          } else if (token_slen == 2) {
            if (memequal_k(linebuf_iter, "ID", 2)) {
              cur_col_type = 1;
            } else if (memequal_k(linebuf_iter, "CM", 2)) {
              cur_col_type = 7;
              cm_col_present = 1;
            } else {
              continue;
            }
          } else {
            continue;
          }
        } else if (strequal_k(linebuf_iter, "QUAL", token_slen)) {
          load_qual_col = 2 * ((pvar_psam_flags & (kfPvarColMaybequal | kfPvarColQual)) || qualfilter_needed) + (var_min_qual != -1);
          if (!load_qual_col) {
            continue;
          }
          cur_col_type = 4;
        } else if (strequal_k(linebuf_iter, "INFO", token_slen)) {
          cur_col_type = 6;
          info_col_present = 1;
        } else if (token_slen == 6) {
          if (memequal_k(linebuf_iter, "FILTER", 6)) {
            load_filter_col = 2 * ((pvar_psam_flags & (kfPvarColMaybefilter | kfPvarColFilter)) || qualfilter_needed) + ((misc_flags / kfMiscExcludePvarFilterFail) & 1);
            if (!load_filter_col) {
              continue;
            }
            cur_col_type = 5;
          } else if (memequal_k(linebuf_iter, "FORMAT", 6)) {
            break;
          } else {
            continue;
          }
        } else {
          continue;
        }
        const uint32_t cur_col_type_shifted = 1 << cur_col_type;
        if (unlikely(found_header_bitset & cur_col_type_shifted)) {
          // known token, so no overflow danger
          char* write_iter = strcpya_k(g_logbuf, "Error: Duplicate column header '");
          write_iter = memcpya(write_iter, linebuf_iter, token_slen);
          write_iter = strcpya_k(write_iter, "' on line ");
          write_iter = wtoa(line_idx, write_iter);
          write_iter = strcpya_k(write_iter, " of ");
          write_iter = strcpya(write_iter, pvarname);
          memcpy_k(write_iter, ".\n\0", 4);
          goto LoadPvar_ret_MALFORMED_INPUT_WW;
        }
        found_header_bitset |= cur_col_type_shifted;
        col_skips[relevant_postchr_col_ct] = col_idx;
        col_types[relevant_postchr_col_ct++] = cur_col_type;
      }
      if (unlikely((found_header_bitset & 0x0f) != 0x0f)) {
        snprintf(g_logbuf, kLogbufSize, "Error: Missing column header(s) on line %" PRIuPTR " of %s. (POS, ID, REF, and ALT are required.)\n", line_idx, pvarname);
        goto LoadPvar_ret_MALFORMED_INPUT_WW;
      }
      if ((var_min_qual != -1) && (!(found_header_bitset & 0x10))) {
        logerrputs("Error: --var-min-qual used on a variant file with no QUAL column.\n");
        goto LoadPvar_ret_INCONSISTENT_INPUT;
      }
      for (uint32_t rpc_col_idx = relevant_postchr_col_ct - 1; rpc_col_idx; --rpc_col_idx) {
        col_skips[rpc_col_idx] -= col_skips[rpc_col_idx - 1];
      }
      // skip this line in main loop
      line_iter = K_CAST(char*, AdvToDelim(linebuf_iter, '\n'));
      line_start = line_iter;
    } else if (line_start[0] != '\n') {
      if (var_min_qual != -1) {
        logerrputs("Error: --var-min-qual used on a variant file with no QUAL column.\n");
        goto LoadPvar_ret_INCONSISTENT_INPUT;
      }
      *info_flags_ptr = kfInfoPrNonrefDefault;
      col_skips[0] = 1;
      col_skips[1] = 1;
      col_skips[2] = 1;
      col_skips[3] = 1;
      col_types[0] = 1;
      // CM column is formally optional in headerless .pvar files (and it was
      // "secretly" optional for the standard plink 1.9 standard .bim loader).
      // If the line has exactly 5 columns, assume CM is omitted.

      // Note that this doesn't break in the line_start = &(g_one_char_strs[0])
      // case.
      const char* linebuf_iter = NextTokenMult(line_start, 4);
      if (unlikely(!linebuf_iter)) {
        goto LoadPvar_ret_MISSING_TOKENS;
      }
      linebuf_iter = NextToken(linebuf_iter);
      if (!linebuf_iter) {
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

    // done with header.  line_start now points to either the beginning of the
    // first real line, or an eoln character.
    uint32_t max_variant_id_slen = *max_variant_id_slen_ptr;
    uint32_t chrs_encountered_m1 = UINT32_MAX;  // intentional overflow
    uint32_t prev_chr_code = UINT32_MAX;  // force initial mismatch
    uint32_t raw_variant_ct = 0;
    uintptr_t* chr_mask = cip->chr_mask;
    const char* missing_allele_str = &(g_one_char_strs[92]);
    double last_cm = -DBL_MAX;
    int32_t last_bp = 0;

    // this way, we only need to check allele_storage_iter against this (i)
    // when processing a multiallelic variant or (ii) at the end of a block
    const char** allele_storage_limit = R_CAST(const char**, &(rlstream_start[kLoadPvarBlockSize * (-2) * sizeof(intptr_t)]));

    uintptr_t* loaded_chr_mask = R_CAST(uintptr_t*, tmp_alloc_base);
    unsigned char* tmp_alloc_end = bigstack_end_mark;
    // guaranteed to succeed since max_line_blen > 128k, etc.
    /*
    if ((uintptr_t)(tmp_alloc_end - tmp_alloc_base) < RoundUpPow2(kChrMaskWords * sizeof(intptr_t), kCacheline)) {
      goto LoadPvar_ret_NOMEM;
    }
    */
    tmp_alloc_base = &(tmp_alloc_base[RoundUpPow2(kChrMaskWords * sizeof(intptr_t), kCacheline)]);
    // bugfix (2 Jun 2017): forgot to zero-initialize loaded_chr_mask
    ZeroWArr(kChrMaskWords, loaded_chr_mask);

    InfoExist* info_existp = nullptr;
    if (require_info_flattened) {
      reterr = InfoExistInit(tmp_alloc_end, require_info_flattened, "require-info", &tmp_alloc_base, &info_existp);
      if (unlikely(reterr)) {
        goto LoadPvar_ret_1;
      }
    }
    InfoExist* info_nonexistp = nullptr;
    if (require_no_info_flattened) {
      reterr = InfoExistInit(tmp_alloc_end, require_no_info_flattened, "require-no-info", &tmp_alloc_base, &info_nonexistp);
      if (unlikely(reterr)) {
        goto LoadPvar_ret_1;
      }
    }
    InfoFilter info_keep;
    info_keep.prekey = nullptr;
    if (extract_if_info_exprp->pheno_name) {
      // todo: also print warning (or optionally error out?) if header line is
      // missing or doesn't match type expectation
      // (same for --require-info)
      reterr = InfoFilterInit(tmp_alloc_end, extract_if_info_exprp, "extract-if-info", &tmp_alloc_base, &info_keep);
      if (unlikely(reterr)) {
        goto LoadPvar_ret_1;
      }
    }
    InfoFilter info_remove;
    info_remove.prekey = nullptr;
    if (exclude_if_info_exprp->pheno_name) {
      reterr = InfoFilterInit(tmp_alloc_end, exclude_if_info_exprp, "exclude-if-info", &tmp_alloc_base, &info_remove);
      if (unlikely(reterr)) {
        goto LoadPvar_ret_1;
      }
    }
    if (!info_col_present) {
      if (unlikely(require_info_flattened)) {
        logerrputs("Error: --require-info used on a variant file with no INFO column.\n");
        goto LoadPvar_ret_INCONSISTENT_INPUT;
      }
      if (unlikely(info_keep.prekey)) {
        logerrputs("Error: --extract-if-info used on a variant file with no INFO column.\n");
        goto LoadPvar_ret_INCONSISTENT_INPUT;
      }
      info_pr_present = 0;
      info_reload_slen = 0;
    } else if ((!info_pr_present) && (!info_reload_slen) && (!info_existp) && (!info_nonexistp) && (!info_keep.prekey) && (!info_remove.prekey)) {
      info_col_present = 0;
    }

    uint32_t fexcept_ct = 0;
    uintptr_t max_fexcept_blen = 2;
    char* sorted_fexcepts = nullptr;
    if (var_filter_exceptions_flattened) {
      if (unlikely(MultistrToStrboxDedupArenaAlloc(tmp_alloc_end, var_filter_exceptions_flattened, &tmp_alloc_base, &sorted_fexcepts, &fexcept_ct, &max_fexcept_blen))) {
        goto LoadPvar_ret_NOMEM;
      }
    }
    const uint32_t new_variant_id_overflow_missing = (misc_flags / kfMiscNewVarIdOverflowMissing) & 1;
    char* chr_output_name_buf = nullptr;
    VaridTemplate* varid_templatep = nullptr;
    VaridTemplate* varid_multi_templatep = nullptr;
    VaridTemplate* varid_multi_nonsnp_templatep = nullptr;
    uint32_t missing_varid_blen = 0;
    uint32_t missing_varid_match_slen = 0;
    if (varid_template_str) {
      if (unlikely(S_CAST(uintptr_t, tmp_alloc_end - tmp_alloc_base) < kMaxIdSlen + 3 * RoundUpPow2(sizeof(VaridTemplate), kCacheline))) {
        goto LoadPvar_ret_NOMEM;
      }
      chr_output_name_buf = R_CAST(char*, tmp_alloc_base);
      tmp_alloc_base = &(tmp_alloc_base[kMaxIdSlen]);
      if (!missing_varid_match) {
        missing_varid_match = &(g_one_char_strs[92]);  // '.'
      }
      missing_varid_blen = strlen(missing_varid_match);
      if (misc_flags & kfMiscSetMissingVarIds) {
        missing_varid_match_slen = missing_varid_blen;
      }
      ++missing_varid_blen;
      varid_templatep = R_CAST(VaridTemplate*, tmp_alloc_base);
      tmp_alloc_base = &(tmp_alloc_base[RoundUpPow2(sizeof(VaridTemplate), kCacheline)]);
      const uint32_t overflow_substitute_blen = new_variant_id_overflow_missing? missing_varid_blen : 0;
      VaridTemplateInit(varid_template_str, missing_varid_match, chr_output_name_buf, new_variant_id_max_allele_slen, overflow_substitute_blen, varid_templatep);
      if (varid_multi_template_str) {
        varid_multi_templatep = R_CAST(VaridTemplate*, tmp_alloc_base);
        tmp_alloc_base = &(tmp_alloc_base[RoundUpPow2(sizeof(VaridTemplate), kCacheline)]);
        VaridTemplateInit(varid_multi_template_str, missing_varid_match, chr_output_name_buf, new_variant_id_max_allele_slen, overflow_substitute_blen, varid_multi_templatep);
      }
      if (varid_multi_nonsnp_template_str) {
        varid_multi_nonsnp_templatep = R_CAST(VaridTemplate*, tmp_alloc_base);
        tmp_alloc_base = &(tmp_alloc_base[RoundUpPow2(sizeof(VaridTemplate), kCacheline)]);
        VaridTemplateInit(varid_multi_nonsnp_template_str, missing_varid_match, chr_output_name_buf, new_variant_id_max_allele_slen, overflow_substitute_blen, varid_multi_nonsnp_templatep);
      }
    }

    // prevent later return-array allocations from overlapping with temporary
    // storage
    g_bigstack_end = tmp_alloc_base;

    // prevent variant_id_htable_find from breaking
    if (R_CAST(const char*, tmp_alloc_end) > (&(g_one_char_strs[512 - kMaxIdSlen]))) {
      tmp_alloc_end = R_CAST(unsigned char*, K_CAST(char*, &(g_one_char_strs[512 - kMaxIdSlen])));
    }
    const uint32_t allow_extra_chrs = (misc_flags / kfMiscAllowExtraChrs) & 1;
    const uint32_t merge_par = ((misc_flags & (kfMiscMergePar | kfMiscMergeX)) != 0);
    const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
    const char input_missing_geno_char = *g_input_missing_geno_ptr;
    uint32_t parx_code = cip->xymt_codes[kChrOffsetPAR1];
    uint32_t par2_code = cip->xymt_codes[kChrOffsetPAR2];
    if (misc_flags & kfMiscMergeX) {
      parx_code = cip->xymt_codes[kChrOffsetXY];
      par2_code = UINT32_MAXM1;
    }
    uint32_t merge_par_ct = 0;

    // Corner case: with --split-par + --not-chr x, we should keep the
    // pseudoautosomal regions.  To facilitate this, we temporarily don't mask
    // out chrX; SplitPar() handles this properly later.
    const uint32_t splitpar_and_exclude_x = splitpar_bound2 && (!IsI32Neg(x_code)) && (!IsSet(cip->chr_mask, x_code));
    if (splitpar_and_exclude_x) {
      SetBit(x_code, cip->chr_mask);
    }

    uint8_t acgtm_bool_table[256] = {
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
    if (snps_only > 1) {
      acgtm_bool_table[ctou32(input_missing_geno_char)] = 1;
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

    // only need this for --set-{missing,all}-var-ids error message
    uint32_t max_chr_slen = 0;

    uint32_t max_filter_slen = 0;
    uint32_t exclude_ct = 0;

    // only allocated when necessary
    // if we want to scale this approach to more fields, we'll need to add a
    // few pointers to the start of each block.  right now, we force cur_cms[]
    // to be allocated before cur_chr_idxs[] when both are present, but this
    // is error-prone.
    uint32_t at_least_one_npass_filter = 0;
    uint32_t at_least_one_nzero_cm = 0;
    uintptr_t new_variant_id_allele_len_overflow = 0;
    double* cur_cms = nullptr;
    uint32_t cms_start_block = UINT32_MAX;

    ChrIdx* cur_chr_idxs = nullptr;
    uint32_t chr_idxs_start_block = UINT32_MAX;
    uint32_t is_split_chr = 0;
    UnsortedVar vpos_sortstatus = kfUnsortedVar0;

    if (IsEolnKns(*line_start)) {
      ++line_iter;
      ++line_idx;
    } else {
      line_iter = line_start;
    }
    for (; TextGetUnsafe2(&pvar_txs, &line_iter); ++line_iter, ++line_idx) {
      if (unlikely(line_iter[0] == '#')) {
        snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of %s starts with a '#'. (This is only permitted before the first nonheader line, and if a #CHROM header line is present it must denote the end of the header block.)\n", line_idx, pvarname);
        goto LoadPvar_ret_MALFORMED_INPUT_WW;
      }
#ifdef __LP64__
      // maximum prime < 2^32 is 4294967291; quadratic hashing guarantee
      // breaks down past that divided by 2.
      if (unlikely(raw_variant_ct == 0x7ffffffd)) {
        logerrputs("Error: " PROG_NAME_STR " does not support more than 2^31 - 3 variants.  We recommend using\nother software for very deep studies of small numbers of genomes.\n");
        goto LoadPvar_ret_MALFORMED_INPUT;
      }
#endif
      const uint32_t variant_idx_lowbits = raw_variant_ct % kLoadPvarBlockSize;
      if (!variant_idx_lowbits) {
        if (unlikely(
                (S_CAST(uintptr_t, tmp_alloc_end - tmp_alloc_base) <=
                  kLoadPvarBlockSize *
                    (sizeof(int32_t) +
                     2 * sizeof(intptr_t) +
                     at_least_one_nzero_cm * sizeof(double)) +
                  is_split_chr * sizeof(ChrIdx) +
                  (1 + info_pr_present) * (kLoadPvarBlockSize / CHAR_BIT) +
                  (load_qual_col? ((kLoadPvarBlockSize / CHAR_BIT) + kLoadPvarBlockSize * sizeof(float)) : 0) +
                  (load_filter_col? (2 * (kLoadPvarBlockSize / CHAR_BIT) + kLoadPvarBlockSize * sizeof(intptr_t)) : 0)) ||
                (allele_storage_iter >= allele_storage_limit))) {
          goto LoadPvar_ret_NOMEM;
        }
        cur_bps = R_CAST(uint32_t*, tmp_alloc_base);
        cur_allele_idxs = R_CAST(uintptr_t*, &(tmp_alloc_base[kLoadPvarBlockSize * sizeof(int32_t)]));
        cur_ids = R_CAST(char**, &(tmp_alloc_base[kLoadPvarBlockSize * (sizeof(int32_t) + sizeof(intptr_t))]));
        cur_include = R_CAST(uintptr_t*, &(tmp_alloc_base[kLoadPvarBlockSize * (sizeof(int32_t) + 2 * sizeof(intptr_t))]));
        SetAllWArr(kLoadPvarBlockSize / kBitsPerWord, cur_include);
        tmp_alloc_base = &(tmp_alloc_base[kLoadPvarBlockSize * (sizeof(int32_t) + 2 * sizeof(intptr_t)) + (kLoadPvarBlockSize / CHAR_BIT)]);
        if (load_qual_col > 1) {
          cur_qual_present = R_CAST(uintptr_t*, tmp_alloc_base);
          ZeroWArr(kLoadPvarBlockSize / kBitsPerWord, cur_qual_present);
          cur_quals = R_CAST(float*, &(tmp_alloc_base[kLoadPvarBlockSize / CHAR_BIT]));
          tmp_alloc_base = &(tmp_alloc_base[kLoadPvarBlockSize * sizeof(float) + (kLoadPvarBlockSize / CHAR_BIT)]);
        }
        if (load_filter_col > 1) {
          cur_filter_present = R_CAST(uintptr_t*, tmp_alloc_base);
          cur_filter_npass = R_CAST(uintptr_t*, &(tmp_alloc_base[kLoadPvarBlockSize / CHAR_BIT]));
          cur_filter_storage = R_CAST(char**, &(tmp_alloc_base[2 * (kLoadPvarBlockSize / CHAR_BIT)]));
          ZeroWArr(kLoadPvarBlockSize / kBitsPerWord, cur_filter_present);
          ZeroWArr(kLoadPvarBlockSize / kBitsPerWord, cur_filter_npass);
          tmp_alloc_base = &(tmp_alloc_base[2 * (kLoadPvarBlockSize / CHAR_BIT) + kLoadPvarBlockSize * sizeof(intptr_t)]);
        }
        if (info_pr_present) {
          cur_nonref_flags = R_CAST(uintptr_t*, tmp_alloc_base);
          ZeroWArr(kLoadPvarBlockSize / kBitsPerWord, cur_nonref_flags);
          tmp_alloc_base = &(tmp_alloc_base[kLoadPvarBlockSize / CHAR_BIT]);
        }
        if (at_least_one_nzero_cm) {
          cur_cms = R_CAST(double*, tmp_alloc_base);
          ZeroDArr(kLoadPvarBlockSize, cur_cms);
          tmp_alloc_base = R_CAST(unsigned char*, &(cur_cms[kLoadPvarBlockSize]));
        }
        if (is_split_chr) {
          cur_chr_idxs = R_CAST(ChrIdx*, tmp_alloc_base);
          tmp_alloc_base = R_CAST(unsigned char*, &(cur_chr_idxs[kLoadPvarBlockSize]));
        }
      }
      char* linebuf_iter = CurTokenEnd(line_iter);
      // #CHROM
      if (unlikely(*linebuf_iter == '\n')) {
        goto LoadPvar_ret_MISSING_TOKENS;
      }
      uint32_t cur_chr_code;
      reterr = GetOrAddChrCodeDestructive(".pvar file", line_idx, allow_extra_chrs, line_iter, linebuf_iter, cip, &cur_chr_code);
      if (unlikely(reterr)) {
        goto LoadPvar_ret_1;
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
      if (cur_chr_code != prev_chr_code) {
        prev_chr_code = cur_chr_code;
        if (!is_split_chr) {
          if (IsSet(loaded_chr_mask, cur_chr_code)) {
            if (unlikely(!split_chr_ok)) {
              snprintf(g_logbuf, kLogbufSize, "Error: %s has a split chromosome. Use --make-pgen + --sort-vars to remedy this.\n", pvarname);
              goto LoadPvar_ret_MALFORMED_INPUT_WW;
            }
            if (unlikely(S_CAST(uintptr_t, tmp_alloc_end - tmp_alloc_base) < kLoadPvarBlockSize * sizeof(ChrIdx))) {
              goto LoadPvar_ret_NOMEM;
            }
            cur_chr_idxs = R_CAST(ChrIdx*, tmp_alloc_base);
            tmp_alloc_base = R_CAST(unsigned char*, &(cur_chr_idxs[kLoadPvarBlockSize]));
            // may want to track the first problem variant index
            // cip->chr_fo_vidx_start[chrs_encountered_m1] = raw_variant_ct;
            BackfillChrIdxs(cip, chrs_encountered_m1, RoundDownPow2(raw_variant_ct, kLoadPvarBlockSize), raw_variant_ct, cur_chr_idxs);
            chr_idxs_start_block = raw_variant_ct / kLoadPvarBlockSize;
            is_split_chr = 1;
            vpos_sortstatus |= kfUnsortedVarBp | kfUnsortedVarCm | kfUnsortedVarSplitChr;
          } else {
            // how much of this do we need in split-chrom case?
            cip->chr_file_order[++chrs_encountered_m1] = cur_chr_code;
            cip->chr_fo_vidx_start[chrs_encountered_m1] = raw_variant_ct;
            cip->chr_idx_to_foidx[cur_chr_code] = chrs_encountered_m1;
            last_cm = -DBL_MAX;
          }
        }
        // always need to set this, if we want to avoid chromosome-length
        // overestimation in --sort-vars + early-variant-filter +
        // pvar-cols=+vcfheader case.
        last_bp = 0;

        SetBit(cur_chr_code, loaded_chr_mask);
        if (chr_output_name_buf) {
          char* chr_name_end = chrtoa(cip, cur_chr_code, chr_output_name_buf);
          const uint32_t chr_slen = chr_name_end - chr_output_name_buf;
          if (chr_slen > max_chr_slen) {
            max_chr_slen = chr_slen;
          }
          const int32_t chr_slen_delta = chr_slen - varid_templatep->chr_slen;
          varid_templatep->chr_slen = chr_slen;
          varid_templatep->base_len += chr_slen_delta;
          if (varid_multi_templatep) {
            varid_multi_templatep->chr_slen = chr_slen;
            varid_multi_templatep->base_len += chr_slen_delta;
          }
          if (varid_multi_nonsnp_templatep) {
            varid_multi_nonsnp_templatep->chr_slen = chr_slen;
            varid_multi_nonsnp_templatep->base_len += chr_slen_delta;
          }
        }
      }
      *linebuf_iter = '\t';

      // could make this store (and cur_allele_idxs[] allocation) conditional
      // on a multiallelic variant being sighted, but unlike the CM column
      // this should become common
      cur_allele_idxs[variant_idx_lowbits] = allele_storage_iter - allele_storage;

      char* token_ptrs[8];
      uint32_t token_slens[8];
      uint32_t extra_alt_ct;
      if (IsSet(chr_mask, cur_chr_code) || info_pr_present) {
        linebuf_iter = TokenLex(linebuf_iter, col_types, col_skips, relevant_postchr_col_ct, token_ptrs, token_slens);
        if (unlikely(!linebuf_iter)) {
          goto LoadPvar_ret_MISSING_TOKENS;
        }

        extra_alt_ct = CountByte(token_ptrs[3], ',', token_slens[3]);
        if (extra_alt_ct > max_extra_alt_ct) {
          max_extra_alt_ct = extra_alt_ct;
        }

        // It is possible for the info_token[info_slen] assignment below to
        // clobber the line terminator, so we advance line_iter to eoln here
        // and never reference it again before the next line.
        line_iter = AdvToDelim(linebuf_iter, '\n');
        if (info_col_present) {
          const uint32_t info_slen = token_slens[6];
          if (info_slen > info_reload_slen) {
            info_reload_slen = info_slen;
          }
          char* info_token = token_ptrs[6];
          info_token[info_slen] = '\0';
          if (info_pr_present) {
            // always load all nonref_flags entries so (i) --ref-from-fa +
            // --make-just-pvar works and (ii) they can be compared against
            // the .pgen.
            if ((memequal_k(info_token, "PR", 2) && ((info_slen == 2) || (info_token[2] == ';'))) || memequal_k(&(info_token[S_CAST(int32_t, info_slen) - 3]), ";PR", 3)) {
              SetBit(variant_idx_lowbits, cur_nonref_flags);
            } else {
              const char* first_info_end = strchr(info_token, ';');
              if (first_info_end && strstr(first_info_end, ";PR;")) {
                SetBit(variant_idx_lowbits, cur_nonref_flags);
              }
            }
            if (!IsSet(chr_mask, cur_chr_code)) {
              goto LoadPvar_skip_variant;
            }
          }
          if (info_existp) {
            if (!InfoExistCheck(info_token, info_existp)) {
              goto LoadPvar_skip_variant;
            }
          }
          if (info_nonexistp) {
            if (!InfoNonexistCheck(info_token, info_nonexistp)) {
              goto LoadPvar_skip_variant;
            }
          }
          if (info_keep.prekey) {
            if (!InfoConditionSatisfied(info_token, &info_keep)) {
              goto LoadPvar_skip_variant;
            }
          }
          if (info_remove.prekey) {
            if (InfoConditionSatisfied(info_token, &info_remove)) {
              goto LoadPvar_skip_variant;
            }
          }
        }
        // POS
        int32_t cur_bp;
        if (unlikely(ScanIntAbsDefcap(token_ptrs[0], &cur_bp))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Invalid bp coordinate on line %" PRIuPTR " of %s.\n", line_idx, pvarname);
          goto LoadPvar_ret_MALFORMED_INPUT_WW;
        }

        if (cur_bp < 0) {
          goto LoadPvar_skip_variant;
        }

        // QUAL
        if (load_qual_col) {
          const char* qual_token = token_ptrs[4];
          if ((qual_token[0] != '.') || (qual_token[1] > ' ')) {
            float cur_qual;
            if (unlikely(ScanFloat(qual_token, &cur_qual))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid QUAL value on line %" PRIuPTR " of %s.\n", line_idx, pvarname);
              goto LoadPvar_ret_MALFORMED_INPUT_WW;
            }
            if ((load_qual_col & 1) && (cur_qual < var_min_qual)) {
              goto LoadPvar_skip_variant;
            }
            if (load_qual_col > 1) {
              SetBit(variant_idx_lowbits, cur_qual_present);
              // possible todo: optimize all-quals-same case
              // possible todo: conditionally allocate, like cur_cms
              cur_quals[variant_idx_lowbits] = cur_qual;
            }
          } else if (load_qual_col & 1) {
            goto LoadPvar_skip_variant;
          }
        }

        // avoid repeating the ALT string split in --set-...-var-ids case
        linebuf_iter = token_ptrs[3];
        uint32_t remaining_alt_char_ct = token_slens[3];
        // handle --snps-only here instead of later, since it reduces the
        // amount of data we need to load
        if (snps_only) {
          if ((token_slens[2] != 1) || (remaining_alt_char_ct != 2 * extra_alt_ct + 1)) {
            goto LoadPvar_skip_variant;
          }
          if (snps_only > 1) {
            // just-acgt
            if (!acgtm_bool_table[ctou32(token_ptrs[2][0])]) {
              goto LoadPvar_skip_variant;
            }
            for (uint32_t uii = 0; uii <= extra_alt_ct; ++uii) {
              if (!acgtm_bool_table[ctou32(linebuf_iter[2 * uii])]) {
                goto LoadPvar_skip_variant;
              }
            }
          }
        }

        // also handle --min-alleles/--max-alleles here
        if (filter_min_allele_ct || (filter_max_allele_ct <= kPglMaxAltAlleleCt)) {
          uint32_t allele_ct = extra_alt_ct + 2;
          if (!extra_alt_ct) {
            // allele_ct == 1 or 2 for filtering purposes, depending on
            // whether ALT allele matches a missing code.
            if ((remaining_alt_char_ct == 1) && ((linebuf_iter[0] == '.') || (linebuf_iter[0] == input_missing_geno_char))) {
              allele_ct = 1;
            }
          }
          if ((allele_ct < filter_min_allele_ct) || (allele_ct > filter_max_allele_ct)) {
            goto LoadPvar_skip_variant;
          }
        }

        // FILTER
        if (load_filter_col) {
          const char* filter_token = token_ptrs[5];
          const uint32_t filter_slen = token_slens[5];
          if ((filter_slen > 1) || (filter_token[0] != '.')) {
            if (!strequal_k(filter_token, "PASS", filter_slen)) {
              if (load_filter_col & 1) {
                if (!fexcept_ct) {
                  goto LoadPvar_skip_variant;
                }
                const char* filter_token_end = &(filter_token[filter_slen]);
                for (const char* filter_token_iter = filter_token; ; ) {
                  const char* cur_filter_name_end = AdvToDelimOrEnd(filter_token_iter, filter_token_end, ';');
                  uint32_t cur_slen = cur_filter_name_end - filter_token_iter;
                  // possible todo: error out on "PASS", since that
                  // shouldn't coexist with other filters
                  // possible todo: maintain a dictionary of FILTER
                  // strings, analogous to what BCF2 does on disk
                  if (bsearch_str(filter_token_iter, sorted_fexcepts, cur_slen, max_fexcept_blen, fexcept_ct) == -1) {
                    goto LoadPvar_skip_variant;
                  }
                  if (cur_filter_name_end == filter_token_end) {
                    break;
                  }
                  filter_token_iter = &(cur_filter_name_end[1]);
                }
              }
              if (load_filter_col > 1) {
                SetBit(variant_idx_lowbits, cur_filter_npass);
                at_least_one_npass_filter = 1;
                // possible todo: detect repeated filter values, store more
                // compactly
                if (filter_slen > max_filter_slen) {
                  max_filter_slen = filter_slen;
                }
                if (StoreStringAtEnd(tmp_alloc_base, filter_token, filter_slen, &tmp_alloc_end, &(cur_filter_storage[variant_idx_lowbits]))) {
                  goto LoadPvar_ret_NOMEM;
                }
              }
            }
            if (load_filter_col > 1) {
              SetBit(variant_idx_lowbits, cur_filter_present);
            }
          }
        }

        if (cur_chr_idxs) {
          cur_chr_idxs[variant_idx_lowbits] = cur_chr_code;
        }
        if (cur_bp < last_bp) {
          vpos_sortstatus |= kfUnsortedVarBp;
        }
        cur_bps[variant_idx_lowbits] = cur_bp;
        last_bp = cur_bp;
        const uint32_t ref_slen = token_slens[2];
        uint32_t id_slen;
        if ((!varid_templatep) || (missing_varid_match_slen && ((token_slens[1] != missing_varid_match_slen) || (!memequal(token_ptrs[1], missing_varid_match, missing_varid_match_slen))))) {
          id_slen = token_slens[1];
          if (PtrWSubCk(tmp_alloc_base, id_slen + 1, &tmp_alloc_end)) {
            goto LoadPvar_ret_NOMEM;
          }
          memcpyx(tmp_alloc_end, token_ptrs[1], id_slen, '\0');
        } else {
          VaridTemplate* cur_varid_templatep = varid_templatep;
          if (extra_alt_ct && (varid_multi_templatep || varid_multi_nonsnp_templatep)) {
            if (varid_multi_templatep) {
              cur_varid_templatep = varid_multi_templatep;
            }
            if (varid_multi_nonsnp_templatep) {
              if ((ref_slen > 1) || (remaining_alt_char_ct != 2 * extra_alt_ct + 1)) {
                cur_varid_templatep = varid_multi_nonsnp_templatep;
              }
            }
          }
          if (unlikely(VaridTemplateApply(tmp_alloc_base, cur_varid_templatep, token_ptrs[2], linebuf_iter, cur_bp, token_slens[2], extra_alt_ct, remaining_alt_char_ct, &tmp_alloc_end, &new_variant_id_allele_len_overflow, &id_slen))) {
            goto LoadPvar_ret_NOMEM;
          }
        }
        if (id_slen > max_variant_id_slen) {
          max_variant_id_slen = id_slen;
        }
        cur_ids[variant_idx_lowbits] = R_CAST(char*, tmp_alloc_end);

        // REF
        const char* ref_allele = token_ptrs[2];
        if (ref_slen == 1) {
          char geno_char = ref_allele[0];
          if (geno_char == input_missing_geno_char) {
            geno_char = '.';
          }
          *allele_storage_iter = &(g_one_char_strs[2 * ctou32(geno_char)]);
        } else {
          if (StoreStringAtEndK(tmp_alloc_base, ref_allele, ref_slen, &tmp_alloc_end, allele_storage_iter)) {
            goto LoadPvar_ret_NOMEM;
          }
          if (ref_slen > max_allele_slen) {
            max_allele_slen = ref_slen;
          }
        }
        ++allele_storage_iter;

        // ALT
        if (extra_alt_ct) {
          if (PtrCheck(allele_storage_limit, allele_storage_iter, extra_alt_ct * sizeof(intptr_t))) {
            goto LoadPvar_ret_NOMEM;
          }
          char* alt_token_end = &(linebuf_iter[remaining_alt_char_ct]);
          for (uint32_t alt_idx = 0; alt_idx != extra_alt_ct; ++alt_idx) {
            char* cur_alt_end = AdvToDelim(linebuf_iter, ',');
            const uint32_t cur_allele_slen = cur_alt_end - linebuf_iter;
            if (cur_allele_slen == 1) {
              char geno_char = linebuf_iter[0];
              if (geno_char == input_missing_geno_char) {
                geno_char = '.';
              }
              *allele_storage_iter = &(g_one_char_strs[2 * ctou32(geno_char)]);
            } else {
              if (unlikely(!cur_allele_slen)) {
                goto LoadPvar_ret_EMPTY_ALLELE_CODE;
              }
              if (StoreStringAtEndK(tmp_alloc_base, linebuf_iter, cur_allele_slen, &tmp_alloc_end, allele_storage_iter)) {
                goto LoadPvar_ret_NOMEM;
              }
              if (cur_allele_slen > max_allele_slen) {
                max_allele_slen = cur_allele_slen;
              }
            }
            ++allele_storage_iter;
            linebuf_iter = &(cur_alt_end[1]);
          }
          remaining_alt_char_ct = alt_token_end - linebuf_iter;
          if (unlikely(!remaining_alt_char_ct)) {
            goto LoadPvar_ret_EMPTY_ALLELE_CODE;
          }
        }
        if (remaining_alt_char_ct == 1) {
          char geno_char = linebuf_iter[0];
          if (geno_char == input_missing_geno_char) {
            geno_char = '.';
          }
          *allele_storage_iter = &(g_one_char_strs[2 * ctou32(geno_char)]);
        } else {
          if (StoreStringAtEndK(tmp_alloc_base, linebuf_iter, remaining_alt_char_ct, &tmp_alloc_end, allele_storage_iter)) {
            goto LoadPvar_ret_NOMEM;
          }
          if (remaining_alt_char_ct > max_allele_slen) {
            max_allele_slen = remaining_alt_char_ct;
          }
        }
        ++allele_storage_iter;

        // CM
        if (cm_col_present) {
          const char* cm_token = token_ptrs[7];
          if ((cm_token[0] != '0') || (cm_token[1] > ' ')) {
            double cur_cm;
            if (unlikely(!ScantokDouble(cm_token, &cur_cm))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid centimorgan position on line %" PRIuPTR " of %s.\n", line_idx, pvarname);
              goto LoadPvar_ret_MALFORMED_INPUT_WW;
            }
            if (cur_cm < last_cm) {
              vpos_sortstatus |= kfUnsortedVarCm;
            } else {
              last_cm = cur_cm;
            }
            if (cur_cm != 0.0) {
              if (!at_least_one_nzero_cm) {
                if (unlikely(S_CAST(uintptr_t, tmp_alloc_end - tmp_alloc_base) < kLoadPvarBlockSize * sizeof(double))) {
                  goto LoadPvar_ret_NOMEM;
                }
                if (cur_chr_idxs) {
                  // reposition cur_chr_idxs[] after cur_cms[]
                  cur_cms = R_CAST(double*, cur_chr_idxs);
                  cur_chr_idxs = R_CAST(ChrIdx*, &(cur_cms[kLoadPvarBlockSize]));
                  memcpy(cur_chr_idxs, cur_cms, kLoadPvarBlockSize * sizeof(ChrIdx));
                  tmp_alloc_base = R_CAST(unsigned char*, &(cur_chr_idxs[kLoadPvarBlockSize]));
                } else {
                  cur_cms = R_CAST(double*, tmp_alloc_base);
                  tmp_alloc_base = R_CAST(unsigned char*, &(cur_cms[kLoadPvarBlockSize]));
                }
                ZeroDArr(kLoadPvarBlockSize, cur_cms);
                cms_start_block = raw_variant_ct / kLoadPvarBlockSize;
                at_least_one_nzero_cm = 1;
              }
              cur_cms[variant_idx_lowbits] = cur_cm;
            }
          }
        }
      } else {
        {
          // linebuf_iter guaranteed to be at '\t' after chromosome code
          char* alt_col_start = NextTokenMult(linebuf_iter, alt_col_idx);
          if (unlikely(!alt_col_start)) {
            goto LoadPvar_ret_MISSING_TOKENS;
          }
          char* alt_col_end = CurTokenEnd(alt_col_start);
          extra_alt_ct = CountByte(alt_col_start, ',', alt_col_end - alt_col_start);
          line_iter = AdvToDelim(alt_col_end, '\n');
        }
      LoadPvar_skip_variant:
        ++exclude_ct;
        ClearBit(variant_idx_lowbits, cur_include);
        cur_bps[variant_idx_lowbits] = last_bp;
        // need to advance allele_storage_iter for later allele_idx_offsets
        // lookups to work properly
        *allele_storage_iter++ = missing_allele_str;
        *allele_storage_iter++ = missing_allele_str;
        if (extra_alt_ct) {
          if (PtrCheck(allele_storage_limit, allele_storage_iter, extra_alt_ct * sizeof(intptr_t))) {
            goto LoadPvar_ret_NOMEM;
          }
          for (uint32_t uii = 0; uii != extra_alt_ct; ++uii) {
            *allele_storage_iter++ = missing_allele_str;
          }
        }
      }
      ++raw_variant_ct;
    }
    if (unlikely(TextStreamErrcode2(&pvar_txs, &reterr))) {
      goto LoadPvar_ret_TSTREAM_FAIL;
    }
    reterr = kPglRetSuccess;
    if (unlikely(max_variant_id_slen > kMaxIdSlen)) {
      logerrputs("Error: Variant names are limited to " MAX_ID_SLEN_STR " characters.\n");
      goto LoadPvar_ret_MALFORMED_INPUT;
    }
    if (new_variant_id_allele_len_overflow) {
      if (new_variant_id_overflow_missing) {
        logerrprintfww("Warning: %" PRIuPTR " variant ID%s %s due to allele code length.\n", new_variant_id_allele_len_overflow, (new_variant_id_allele_len_overflow == 1)? "" : "s", missing_varid_match_slen? "unchanged by --set-missing-var-ids" : "erased by --set-all-var-ids");
        if (max_variant_id_slen < missing_varid_blen - 1) {
          max_variant_id_slen = missing_varid_blen - 1;
        }
      } else if (likely(misc_flags & kfMiscNewVarIdOverflowTruncate)) {
        // this is less likely in practice than the error branch, but let's
        // keep our usage of likely()/unlikely() as easy to interpret as
        // possible.
        logerrprintf("Warning: %" PRIuPTR " allele code%s truncated by --set-%s-var-ids.\n", new_variant_id_allele_len_overflow, (new_variant_id_allele_len_overflow == 1)? "" : "s", missing_varid_match_slen? "missing" : "all");
      } else {
        uint32_t worst_case_id_slen = VaridWorstCaseSlen(varid_templatep, max_chr_slen, max_allele_slen);
        if ((worst_case_id_slen <= kMaxIdSlen) && varid_multi_templatep) {
          const uint32_t tmp_slen = VaridWorstCaseSlen(varid_multi_templatep, max_chr_slen, max_allele_slen);
          if (tmp_slen > worst_case_id_slen) {
            worst_case_id_slen = tmp_slen;
          }
        }
        if ((worst_case_id_slen <= kMaxIdSlen) && varid_multi_nonsnp_templatep) {
          const uint32_t tmp_slen = VaridWorstCaseSlen(varid_multi_nonsnp_templatep, max_chr_slen, max_allele_slen);
          if (tmp_slen > worst_case_id_slen) {
            worst_case_id_slen = tmp_slen;
          }
        }
        logerrprintf("Error: %" PRIuPTR " allele code%s too long for --set-%s-var-ids.\n", new_variant_id_allele_len_overflow, (new_variant_id_allele_len_overflow == 1)? "" : "s", missing_varid_match_slen? "missing" : "all");
        if (worst_case_id_slen <= kMaxIdSlen) {
          logerrprintfww("The longest observed allele code in this dataset has length %u. If you're fine with the corresponding ID length, rerun with \"--new-id-max-allele-len %u\" added to your command line.\n", max_allele_slen, max_allele_slen);
          logerrputs("Otherwise, use \"--new-id-max-allele-len <limit> missing\" to set the IDs of all\nvariants with an allele code longer than the given length-limit to '.' (and\nthen process those variants with another script, if necessary).\n");
        } else {
          logerrprintfww("The longest observed allele code in this dataset has length %u. We recommend deciding on a length-limit, and then adding \"--new-id-max-allele-len <limit> missing\" to your command line to cause variants with longer allele codes to be assigned '.' IDs. (You can then process just those variants with another script, if necessary.)\n", max_allele_slen);
        }
        goto LoadPvar_ret_INCONSISTENT_INPUT;
      }
    }
    *max_variant_id_slen_ptr = max_variant_id_slen;
    *max_allele_ct_ptr = max_extra_alt_ct + 2;
    *max_allele_slen_ptr = max_allele_slen;
    *max_filter_slen_ptr = max_filter_slen;
    *raw_variant_ct_ptr = raw_variant_ct;
    uintptr_t allele_idx_end = allele_storage_iter - allele_storage;
    BigstackFinalizeCp(allele_storage, allele_idx_end);
    // We may clobber this object soon, so close it now (verifying
    // rewindability first, if necessary).
    if (unlikely(CleanupTextStream2(pvarname, &pvar_txs, &reterr))) {
      goto LoadPvar_ret_1;
    }
    uintptr_t* allele_idx_offsets = nullptr;
    const uint32_t full_block_ct = raw_variant_ct / kLoadPvarBlockSize;
    const uintptr_t raw_variant_ct_lowbits = raw_variant_ct % kLoadPvarBlockSize;
    // todo: determine whether we want variant_include to be guaranteed to be
    // terminated by a zero bit
    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    if (unlikely(
            bigstack_alloc_w(raw_variant_ctl, variant_include_ptr) ||
            bigstack_alloc_u32(raw_variant_ct, variant_bps_ptr) ||
            bigstack_alloc_cp(raw_variant_ct, variant_ids_ptr))) {
      goto LoadPvar_ret_NOMEM;
    }
    uintptr_t* qual_present = nullptr;
    float* quals = nullptr;
    if (load_qual_col > 1) {
      if (unlikely(
              bigstack_alloc_w(raw_variant_ctl, qual_present_ptr) ||
              bigstack_alloc_f(raw_variant_ct, quals_ptr))) {
        goto LoadPvar_ret_NOMEM;
      }
      qual_present = *qual_present_ptr;
      quals = *quals_ptr;
    }
    uintptr_t* filter_present = nullptr;
    uintptr_t* filter_npass = nullptr;
    char** filter_storage = nullptr;
    if (load_filter_col > 1) {
      if (unlikely(
              bigstack_alloc_w(raw_variant_ctl, filter_present_ptr) ||
              bigstack_alloc_w(raw_variant_ctl, filter_npass_ptr))) {
        goto LoadPvar_ret_NOMEM;
      }
      filter_present = *filter_present_ptr;
      filter_npass = *filter_npass_ptr;
      if (at_least_one_npass_filter) {
        // possible todo: store this in a sparse manner
        if (unlikely(bigstack_alloc_cp(raw_variant_ct, filter_storage_ptr))) {
          goto LoadPvar_ret_NOMEM;
        }
        filter_storage = *filter_storage_ptr;
      }
    }
    uintptr_t* nonref_flags = nullptr;
    if (info_pr_present) {
      if (unlikely(bigstack_alloc_w(raw_variant_ctl, nonref_flags_ptr))) {
        goto LoadPvar_ret_NOMEM;
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
    //   kLoadPvarBlockSize * sizeof(ChrIdx)
    unsigned char* read_iter = g_bigstack_end;
    uint32_t* variant_bps = *variant_bps_ptr;
    char** variant_ids = *variant_ids_ptr;
    uintptr_t* variant_include = *variant_include_ptr;
    for (uint32_t block_idx = 0; block_idx != full_block_ct; ++block_idx) {
      memcpy(&(variant_bps[block_idx * kLoadPvarBlockSize]), read_iter, kLoadPvarBlockSize * sizeof(int32_t));
      // skip over allele_idx_offsets
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
        read_iter = &(read_iter[kLoadPvarBlockSize * sizeof(ChrIdx)]);
      }
    }
    memcpy(&(variant_bps[full_block_ct * kLoadPvarBlockSize]), read_iter, raw_variant_ct_lowbits * sizeof(int32_t));
    read_iter = &(read_iter[kLoadPvarBlockSize * (sizeof(int32_t) + sizeof(intptr_t))]);
    memcpy(&(variant_ids[full_block_ct * kLoadPvarBlockSize]), read_iter, raw_variant_ct_lowbits * sizeof(intptr_t));
    read_iter = &(read_iter[kLoadPvarBlockSize * sizeof(intptr_t)]);
    const uint32_t last_bitblock_size = DivUp(raw_variant_ct_lowbits, CHAR_BIT);
    memcpy(&(variant_include[full_block_ct * (kLoadPvarBlockSize / kBitsPerWord)]), read_iter, last_bitblock_size);
    ZeroTrailingBits(raw_variant_ct, variant_include);
    read_iter = &(read_iter[kLoadPvarBlockSize / CHAR_BIT]);
    if (qual_present) {
      memcpy(&(qual_present[full_block_ct * (kLoadPvarBlockSize / kBitsPerWord)]), read_iter, last_bitblock_size);
      ZeroTrailingBits(raw_variant_ct, qual_present);
      read_iter = &(read_iter[kLoadPvarBlockSize / CHAR_BIT]);
      memcpy(&(quals[full_block_ct * kLoadPvarBlockSize]), read_iter, raw_variant_ct_lowbits * sizeof(float));
      read_iter = &(read_iter[kLoadPvarBlockSize * sizeof(float)]);
    }
    if (filter_present) {
      memcpy(&(filter_present[full_block_ct * (kLoadPvarBlockSize / kBitsPerWord)]), read_iter, last_bitblock_size);
      ZeroTrailingBits(raw_variant_ct, filter_present);
      read_iter = &(read_iter[kLoadPvarBlockSize / CHAR_BIT]);
      memcpy(&(filter_npass[full_block_ct * (kLoadPvarBlockSize / kBitsPerWord)]), read_iter, last_bitblock_size);
      ZeroTrailingBits(raw_variant_ct, filter_npass);
      read_iter = &(read_iter[kLoadPvarBlockSize / CHAR_BIT]);
      if (filter_storage) {
        memcpy(&(filter_storage[full_block_ct * kLoadPvarBlockSize]), read_iter, raw_variant_ct_lowbits * sizeof(intptr_t));
        read_iter = &(read_iter[kLoadPvarBlockSize * sizeof(intptr_t)]);
      }
    }
    if (info_pr_present) {
      memcpy(&(nonref_flags[full_block_ct * (kLoadPvarBlockSize / kBitsPerWord)]), read_iter, last_bitblock_size);
      ZeroTrailingBits(raw_variant_ct, nonref_flags);
      // read_iter = &(read_iter[kLoadPvarBlockSize / CHAR_BIT]);
    }
    const uintptr_t read_iter_stride_base = kLoadPvarBlockSize * (sizeof(int32_t) + 2 * sizeof(intptr_t)) + (kLoadPvarBlockSize / CHAR_BIT) + (load_qual_col > 1) * ((kLoadPvarBlockSize / CHAR_BIT) + kLoadPvarBlockSize * sizeof(float)) + (load_filter_col > 1) * (2 * (kLoadPvarBlockSize / CHAR_BIT) + kLoadPvarBlockSize * sizeof(intptr_t)) + info_pr_present * (kLoadPvarBlockSize / CHAR_BIT);
    if (allele_idx_end > 2 * S_CAST(uintptr_t, raw_variant_ct)) {
      if (unlikely(bigstack_alloc_w(raw_variant_ct + 1, allele_idx_offsets_ptr))) {
        goto LoadPvar_ret_NOMEM;
      }
      allele_idx_offsets = *allele_idx_offsets_ptr;
      uintptr_t* allele_idx_read_iter = R_CAST(uintptr_t*, &(g_bigstack_end[kLoadPvarBlockSize * sizeof(int32_t)]));
      for (uint32_t block_idx = 0; block_idx != full_block_ct; ++block_idx) {
        memcpy(&(allele_idx_offsets[block_idx * kLoadPvarBlockSize]), allele_idx_read_iter, kLoadPvarBlockSize * sizeof(intptr_t));
        allele_idx_read_iter = R_CAST(uintptr_t*, R_CAST(uintptr_t, allele_idx_read_iter) + read_iter_stride_base + (block_idx >= cms_start_block) * kLoadPvarBlockSize * sizeof(double) + (block_idx >= chr_idxs_start_block) * kLoadPvarBlockSize * sizeof(ChrIdx));
      }
      memcpy(&(allele_idx_offsets[full_block_ct * kLoadPvarBlockSize]), allele_idx_read_iter, raw_variant_ct_lowbits * sizeof(intptr_t));
      allele_idx_offsets[raw_variant_ct] = allele_idx_end;
    }
    if (at_least_one_nzero_cm) {
      if (unlikely(bigstack_alloc_d(raw_variant_ct, variant_cms_ptr))) {
        goto LoadPvar_ret_NOMEM;
      }
      double* variant_cms = *variant_cms_ptr;
      ZeroDArr(cms_start_block * kLoadPvarBlockSize, variant_cms);
      double* cms_read_iter = R_CAST(double*, &(g_bigstack_end[read_iter_stride_base * (cms_start_block + 1)]));
      if (cms_start_block > chr_idxs_start_block) {
        cms_read_iter = R_CAST(double*, R_CAST(uintptr_t, cms_read_iter) + kLoadPvarBlockSize * sizeof(ChrIdx) * (cms_start_block - chr_idxs_start_block));
      }
      for (uint32_t block_idx = cms_start_block; block_idx != full_block_ct; ++block_idx) {
        memcpy(&(variant_cms[block_idx * kLoadPvarBlockSize]), cms_read_iter, kLoadPvarBlockSize * sizeof(double));
        cms_read_iter = R_CAST(double*, R_CAST(uintptr_t, cms_read_iter) + read_iter_stride_base + kLoadPvarBlockSize * sizeof(double) + (block_idx >= chr_idxs_start_block) * kLoadPvarBlockSize * sizeof(ChrIdx));
      }
      memcpy(&(variant_cms[full_block_ct * kLoadPvarBlockSize]), cms_read_iter, raw_variant_ct_lowbits * sizeof(double));
    } else {
      *variant_cms_ptr = nullptr;
    }
    if (!is_split_chr) {
      cip->chr_fo_vidx_start[chrs_encountered_m1 + 1] = raw_variant_ct;
      if (splitpar_bound2) {
        if (splitpar_and_exclude_x) {
          ClearBit(x_code, chr_mask);
        }
        reterr = SplitPar(variant_bps, *vpos_sortstatus_ptr, splitpar_bound1, splitpar_bound2, variant_include, loaded_chr_mask, cip, &chrs_encountered_m1, &exclude_ct);
        if (unlikely(reterr)) {
          goto LoadPvar_ret_1;
        }
      }
      cip->chr_ct = chrs_encountered_m1 + 1;
    } else {
      ChrIdx* chr_idxs;
      if (unlikely(BIGSTACK_ALLOC_X(ChrIdx, raw_variant_ct, &chr_idxs))) {
        goto LoadPvar_ret_NOMEM;
      }
      *chr_idxs_ptr = chr_idxs;
      if (chr_idxs_start_block) {
        const uint32_t end_vidx = chr_idxs_start_block * kLoadPvarBlockSize;
        uint32_t chr_fo_idx = chrs_encountered_m1;
        while (cip->chr_fo_vidx_start[chr_fo_idx] >= end_vidx) {
          --chr_fo_idx;
        }
        BackfillChrIdxs(cip, chr_fo_idx, 0, end_vidx, chr_idxs);
      }
      ChrIdx* chr_idxs_read_iter = R_CAST(ChrIdx*, &(g_bigstack_end[read_iter_stride_base * (chr_idxs_start_block + 1)]));
      if (chr_idxs_start_block >= cms_start_block) {
        chr_idxs_read_iter = R_CAST(ChrIdx*, R_CAST(uintptr_t, chr_idxs_read_iter) + kLoadPvarBlockSize * sizeof(double) * (chr_idxs_start_block + 1 - cms_start_block));
      }
      for (uint32_t block_idx = chr_idxs_start_block; block_idx != full_block_ct; ) {
        memcpy(&(chr_idxs[block_idx * kLoadPvarBlockSize]), chr_idxs_read_iter, kLoadPvarBlockSize * sizeof(ChrIdx));
        ++block_idx;
        chr_idxs_read_iter = R_CAST(ChrIdx*, R_CAST(uintptr_t, chr_idxs_read_iter) + read_iter_stride_base + kLoadPvarBlockSize * sizeof(ChrIdx) + (block_idx >= cms_start_block) * kLoadPvarBlockSize * sizeof(double));
      }
      memcpy(&(chr_idxs[full_block_ct * kLoadPvarBlockSize]), chr_idxs_read_iter, raw_variant_ct_lowbits * sizeof(ChrIdx));
      cip->chr_ct = PopcountWords(loaded_chr_mask, DivUp(cip->max_code + cip->name_ct + 1, kBitsPerWord));
    }
    if (merge_par) {
      const uint32_t is_plink2_par = (misc_flags / kfMiscMergePar) & 1;
      if (merge_par_ct) {
        logprintf("--merge-%s: %u chromosome code%s changed.\n", is_plink2_par? "par" : "x", merge_par_ct, (merge_par_ct == 1)? "" : "s");
      } else if (is_plink2_par) {
        logerrputs("Warning: --merge-par had no effect (no PAR1/PAR2 chromosome codes present).\n");
      } else {
        logerrputs("Warning: --merge-x had no effect (no XY chromosome codes present).\n");
      }
    }
    const uint32_t last_chr_code = cip->max_code + cip->name_ct;
    const uint32_t chr_word_ct = BitCtToWordCt(last_chr_code + 1);
    BitvecAnd(loaded_chr_mask, chr_word_ct, chr_mask);
    BigstackEndSet(tmp_alloc_end);
    *variant_ct_ptr = raw_variant_ct - exclude_ct;
    *vpos_sortstatus_ptr = vpos_sortstatus;
    *allele_storage_ptr = allele_storage;
    // if only INFO:PR flag present, no need to reload
    if (!(info_nonpr_present || info_pr_nonflag_present)) {
      info_reload_slen = 0;
    }
    if (info_reload_slen) {
      // treat open-fail as rewind-fail here
      if (unlikely(ForceNonFifo(pvarname))) {
        logerrprintfww(kErrprintfRewind, pvarname);
        reterr = kPglRetRewindFail;
        goto LoadPvar_ret_1;
      }
    }
    *info_reload_slen_ptr = info_reload_slen;
  }

  while (0) {
  LoadPvar_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  LoadPvar_ret_TSTREAM_FAIL:
    TextStreamErrPrint(pvarname, &pvar_txs);
    break;
  LoadPvar_ret_EMPTY_ALLELE_CODE:
    snprintf(g_logbuf, kLogbufSize, "Error: Empty allele code on line %" PRIuPTR " of %s.\n", line_idx, pvarname);
  LoadPvar_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  LoadPvar_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  LoadPvar_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  LoadPvar_ret_MISSING_TOKENS:
    logerrprintfww("Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, pvarname);
    reterr = kPglRetMalformedInput;
    break;
  }
 LoadPvar_ret_1:
  CleanupTextStream2(pvarname, &pvar_txs, &reterr);
  if (reterr) {
    BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  }
  return reterr;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
