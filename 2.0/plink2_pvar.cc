// This library is part of PLINK 2.0, copyright (C) 2005-2025 Shaun Purcell,
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

#include <limits.h>
#include <stddef.h>
#include <string.h>

#include "include/pgenlib_misc.h"
#include "include/plink2_bits.h"
#include "include/plink2_string.h"
#include "include/plink2_text.h"
#include "plink2_decompress.h"

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
//
// possible todo: https://crates.io/crates/memchr should be faster than strstr
// when the needle is constant, port that to include/plink2_string and use it
// here.
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
    uint32_t prev_autosome_ct = 0;
    uint32_t prev_haploid = 0;
    STD_ARRAY_DECL(uint32_t, kChrOffsetCt, prev_xymt_codes);
    if ((cip->chrset_source == kChrsetSourceCmdline) && (misc_flags & kfMiscChrOverrideCmdline)) {
      goto ReadChrsetHeaderLine_ret_1;
    }
    if (((cip->chrset_source == kChrsetSourceCmdline) && (!(misc_flags & kfMiscChrOverrideFile))) || (cip->chrset_source == kChrsetSourceAnotherFile)) {
      // save off info we need for consistency check
      prev_autosome_ct = cip->autosome_ct;
      prev_haploid = cip->haploid_mask[0] & 1;
      STD_ARRAY_COPY(cip->xymt_codes, kChrOffsetCt, prev_xymt_codes);
    }
    // bugfix (23 Jan 2021): forgot to zero this out in file-only case
    ZeroWArr(kChrMaskWords, cip->haploid_mask);
    for (uint32_t uii = 0; uii != kChrOffsetCt; ++uii) {
      cip->xymt_codes[uii] = UINT32_MAXM1;
    }
    if (StrStartsWithUnsafe(chrset_iter, "ID=")) {
      // Need this field to conform to VCFv4.3 specification.
      // Just ignore the ID value.
      chrset_iter = strchrnul_n(&(chrset_iter[strlen("ID=")]), ',');
      if (*chrset_iter != ',') {
        snprintf(g_logbuf, kLogbufSize, "Error: Header line %" PRIuPTR " of %s does not have expected ##chrSet format.\n", line_idx, file_descrip);
        goto ReadChrsetHeaderLine_ret_MALFORMED_INPUT_WW;
      }
      ++chrset_iter;
    }
    if (StrStartsWithUnsafe(chrset_iter, "haploidAutosomeCt=")) {
      uint32_t explicit_haploid_ct;
      if (unlikely(ScanPosintCapped(&(chrset_iter[strlen("haploidAutosomeCt=")]), kMaxChrTextnum, &explicit_haploid_ct))) {
        snprintf(g_logbuf, kLogbufSize, "Error: Header line %" PRIuPTR " of %s has an unsupported ##chrSet haploid count (max %u).\n", line_idx, file_descrip, kMaxChrTextnum);
        goto ReadChrsetHeaderLine_ret_MALFORMED_INPUT_WW;
      }
      // could verify that X, Y, etc. are not present?
      if (prev_autosome_ct) {
        if (unlikely(!prev_haploid)) {
          snprintf(g_logbuf, kLogbufSize, "Error: Header line %" PRIuPTR " of %s specifies a haploid genome, while a diploid genome was specified %s.\n", line_idx, file_descrip, (cip->chrset_source == kChrsetSourceCmdline)? "on the command line" : "in another .pvar");
          goto ReadChrsetHeaderLine_ret_INCONSISTENT_INPUT_WW;
        }
        if (unlikely(explicit_haploid_ct != prev_autosome_ct)) {
          snprintf(g_logbuf, kLogbufSize, "Error: Header line %" PRIuPTR " of %s specifies %u autosome%s, while %s specified %u.\n", line_idx, file_descrip, explicit_haploid_ct, (explicit_haploid_ct == 1)? "" : "s", (cip->chrset_source == kChrsetSourceCmdline)? "the command line" : "another .pvar", prev_autosome_ct);
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
        snprintf(g_logbuf, kLogbufSize, "Error: Header line %" PRIuPTR " of %s has an unsupported ##chrSet autosome count (max %u).\n", line_idx, file_descrip, kMaxChrTextnum);
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
              SetBit(explicit_autosome_ct + 1 + kChrOffsetMT, cip->haploid_mask);
            } else {
              first_char_ui -= 88;  // X = 0, Y = 1, everything else larger
              if (first_char_ui < 2) {
                cip->xymt_codes[first_char_ui] = explicit_autosome_ct + 1 + first_char_ui;
                SetBit(explicit_autosome_ct + 1 + first_char_ui, cip->haploid_mask);
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
                SetBit(explicit_autosome_ct + 1 + kChrOffsetMT, cip->haploid_mask);
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
      if (prev_autosome_ct) {
        if (unlikely(prev_haploid)) {
          snprintf(g_logbuf, kLogbufSize, "Error: Header line %" PRIuPTR " of %s specifies a diploid genome, while a haploid genome was specified %s.\n", line_idx, file_descrip, (cip->chrset_source == kChrsetSourceCmdline)? "on the command line" : "in another .pvar");
          goto ReadChrsetHeaderLine_ret_INCONSISTENT_INPUT_WW;
        }
        if (unlikely(explicit_autosome_ct != prev_autosome_ct)) {
          snprintf(g_logbuf, kLogbufSize, "Error: Header line %" PRIuPTR " of %s specifies %u autosome%s, while %s specified %u.\n", line_idx, file_descrip, explicit_autosome_ct, (explicit_autosome_ct == 1)? "" : "s", (cip->chrset_source == kChrsetSourceCmdline)? "the command line" : "another .pvar", prev_autosome_ct);
          goto ReadChrsetHeaderLine_ret_INCONSISTENT_INPUT_WW;
        }
        for (uint32_t xymt_idx = 0; xymt_idx != kChrOffsetPAR1; ++xymt_idx) {
          // it's okay if the command line doesn't explicitly exclude e.g. chrX
          // while for whatever reason it is excluded from ##chrSet; but the
          // reverse can create problems
          if (unlikely(IsI32Neg(prev_xymt_codes[xymt_idx]) && (!IsI32Neg(cip->xymt_codes[xymt_idx])))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Header line %" PRIuPTR " of %s specifies a chromosome set including %s, while %s excludes it.\n", line_idx, file_descrip, g_xymt_log_names[xymt_idx], (cip->chrset_source == kChrsetSourceCmdline)? "the command line" : "another .pvar");
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

static_assert(kfVaridTemplateAlleleRefOr1 == 1, "VaridTemplateInit() or VaridTemplateAlleleFlags needs to be updated.");
static_assert(kfVaridTemplateAlleleAltOr2 == 2, "VaridTemplateInit() or VaridTemplateAlleleFlags needs to be updated.");
void VaridTemplateInit(const char* varid_template_str, const char* missing_id_match, char* chr_output_name_buf, uint32_t new_id_max_allele_slen, uint32_t overflow_substitute_blen, VaridTemplate* vtp) {
  // template string was previously validated
  // varid_template is only input, everything else is output values
  const char* varid_template_str_iter = varid_template_str;
  uint32_t template_insert_ct = 0;
  uint32_t template_base_len = 0;
  VaridTemplateAlleleFlags allele_flags = kfVaridTemplateAllele0;
  vtp->chr_output_name_buf = chr_output_name_buf;
  vtp->segs[0] = varid_template_str_iter;
  vtp->chr_slen = 0;
  for (unsigned char ucc = *varid_template_str_iter; ucc; ucc = *(++varid_template_str_iter)) {
    if (ucc > '@') {
      continue;
    }
    // VaridTemplateIsValid() ensures ucc is in {'$', '1', '2', 'A', 'R', 'a',
    // 'r'}.
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
      const uint32_t char_code = ctou32(*(++varid_template_str_iter));
      uint32_t flagbit;  // 1 if '1' or 'R'/'r', 2 if '2' or 'A'/'a'
      if (char_code <= '2') {
        allele_flags |= kfVaridTemplateAlleleAsciiOrder;
        flagbit = char_code - 48;
      } else {
        flagbit = 1 + ((char_code & 0xdf) == 'A');
      }
      allele_flags |= S_CAST(VaridTemplateAlleleFlags, flagbit);
      insert_type = flagbit + 1;
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
  vtp->allele_flags = allele_flags;
  vtp->new_id_max_allele_slen = new_id_max_allele_slen;
  vtp->overflow_substitute_blen = overflow_substitute_blen;
  vtp->missing_id_match = missing_id_match;
}

static_assert((!(kMaxIdSlen % kCacheline)), "VaridInitAll() must be updated.");
BoolErr VaridInitAll(unsigned char* arena_end, const char* varid_template_str, const char* varid_multi_template_str, const char* varid_multi_nonsnp_template_str, MiscFlags misc_flags, uint32_t new_variant_id_max_allele_slen, unsigned char** arena_basep, const char** missing_varid_matchp, char** chr_output_name_bufp, VaridTemplate** varid_templatepp, VaridTemplate** varid_multi_templatepp, VaridTemplate** varid_multi_nonsnp_templatepp, uint32_t* missing_varid_blenp, uint32_t* missing_varid_match_slenp) {
  unsigned char* arena_base = *arena_basep;
  if (unlikely(S_CAST(uintptr_t, arena_end - arena_base) < (chr_output_name_bufp != nullptr) * kMaxIdSlen + 3 * RoundUpPow2(sizeof(VaridTemplate), kCacheline))) {
    return 1;
  }
  char* chr_output_name_buf = nullptr;
  if (chr_output_name_bufp) {
    chr_output_name_buf = R_CAST(char*, arena_base);
    *chr_output_name_bufp = chr_output_name_buf;
    arena_base = &(arena_base[kMaxIdSlen]);
  }
  if (!(*missing_varid_matchp)) {
    *missing_varid_matchp = &(g_one_char_strs[92]);  // '.'
  }
  uint32_t missing_varid_blen = strlen(*missing_varid_matchp);
  if (misc_flags & kfMiscSetMissingVarIds) {
    *missing_varid_match_slenp = missing_varid_blen;
  }
  ++missing_varid_blen;
  if (missing_varid_blenp) {
    *missing_varid_blenp = missing_varid_blen;
  }
  *varid_templatepp = R_CAST(VaridTemplate*, arena_base);
  arena_base = &(arena_base[RoundUpPow2(sizeof(VaridTemplate), kCacheline)]);
  const uint32_t new_variant_id_overflow_missing = (misc_flags / kfMiscNewVarIdOverflowMissing) & 1;
  const uint32_t overflow_substitute_blen = new_variant_id_overflow_missing? missing_varid_blen : 0;
  VaridTemplateInit(varid_template_str, *missing_varid_matchp, chr_output_name_buf, new_variant_id_max_allele_slen, overflow_substitute_blen, *varid_templatepp);
  if (varid_multi_template_str) {
    *varid_multi_templatepp = R_CAST(VaridTemplate*, arena_base);
    arena_base = &(arena_base[RoundUpPow2(sizeof(VaridTemplate), kCacheline)]);
    VaridTemplateInit(varid_multi_template_str, *missing_varid_matchp, chr_output_name_buf, new_variant_id_max_allele_slen, overflow_substitute_blen, *varid_multi_templatepp);
  }
  if (varid_multi_nonsnp_template_str) {
    *varid_multi_nonsnp_templatepp = R_CAST(VaridTemplate*, arena_base);
    arena_base = &(arena_base[RoundUpPow2(sizeof(VaridTemplate), kCacheline)]);
    VaridTemplateInit(varid_multi_nonsnp_template_str, *missing_varid_matchp, chr_output_name_buf, new_variant_id_max_allele_slen, overflow_substitute_blen, *varid_multi_nonsnp_templatepp);
  }
  *arena_basep = arena_base;
  return 0;
}

// alt1_end currently allowed to be nullptr in biallelic case
BoolErr VaridTemplateApply(unsigned char* tmp_alloc_base, const VaridTemplate* vtp, const char* ref_start, const char* alt1_start, uint32_t cur_bp, uint32_t ref_token_slen, uint32_t extra_alt_ct, uint32_t alt_token_slen, unsigned char** tmp_alloc_endp, uintptr_t* new_id_allele_len_overflowp, uint32_t* id_slen_ptr) {
  // (insert_ptrs[x], insert_slens[x]) is the (start, slen) of the
  // insert-segment of type x:
  //   0: chromosome
  //   1: 1-based position
  //   2: $r or $1
  //   3: $a or $2
  // Types 0 and 1 are always present, but 2 and/or 3 can be missing.
  uint32_t insert_slens[4];
  const VaridTemplateAlleleFlags allele_flags = vtp->allele_flags;
  const uint32_t new_id_max_allele_slen = vtp->new_id_max_allele_slen;
  uint32_t id_slen = UintSlen(cur_bp);
  insert_slens[0] = vtp->chr_slen;
  insert_slens[1] = id_slen;
  id_slen += vtp->base_len;
  const uint32_t ref_or_1_exists = (allele_flags & kfVaridTemplateAlleleRefOr1)? 1 : 0;
  uint32_t ref_slen = 0;
  uint32_t cur_overflow = 0;
  const char* tmp_allele_ptrs[2];  // [0] = RefOr1, [1] = AltOr2
  tmp_allele_ptrs[0] = nullptr;
  tmp_allele_ptrs[1] = nullptr;  // maybe-uninitialized warning
  if (ref_or_1_exists) {
    ref_slen = ref_token_slen;
    if (ref_slen > new_id_max_allele_slen) {
      ref_slen = new_id_max_allele_slen;
      cur_overflow = 1;
    }
    insert_slens[2] = ref_slen;
    id_slen += ref_slen;
    tmp_allele_ptrs[0] = ref_start;
  }
  if (allele_flags & (kfVaridTemplateAlleleAltOr2 | kfVaridTemplateAlleleAsciiOrder)) {
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
    if (!(allele_flags & kfVaridTemplateAlleleAsciiOrder)) {
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
    id_slen = overflow_substitute_blen - 1;
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
    insert_ptrs[3] = nullptr;
    // insert_ct currently guaranteed to be >= 2 since @ and # required
    // insert_ct > 2 possibilities:
    //   $r without $a: insert_ct=3, use tmp_allele_ptrs[0] and insert_slens[2]
    //   $a without $r: insert_ct=3, use tmp_allele_ptrs[1] and insert_slens[3]
    //                  (fixed bug on 24 Nov 2024)
    //   $r+$a or $1+$2: insert_ct=4, use (tmp_allele_ptrs[0], insert_slens[2])
    //                   followed by (tmp_allele_ptrs[1], insert_slens[3])
    for (uint32_t insert_idx = 0; insert_idx != insert_ct; ++insert_idx) {
      id_iter = memcpya(id_iter, vtp->segs[insert_idx], vtp->seg_lens[insert_idx]);
      const uint32_t cur_insert_type = vtp->insert_types[insert_idx];
      insert_ptrs[cur_insert_type] = id_iter;
      id_iter = &(id_iter[insert_slens[cur_insert_type]]);
    }
    memcpyx(id_iter, vtp->segs[insert_ct], vtp->seg_lens[insert_ct], '\0');

    memcpy(insert_ptrs[0], vtp->chr_output_name_buf, insert_slens[0]);
    u32toa(cur_bp, insert_ptrs[1]);
    if (insert_ct > 2) {
      for (uint32_t insert_type_idx = 2; insert_type_idx != 4; ++insert_type_idx) {
        char* cur_dst = insert_ptrs[insert_type_idx];
        if (cur_dst == nullptr) {
          // bugfix (24 Nov 2024)
          continue;
        }
        memcpy(cur_dst, tmp_allele_ptrs[insert_type_idx - 2], insert_slens[insert_type_idx]);
      }
    }
  }
  *id_slen_ptr = id_slen;
  *new_id_allele_len_overflowp += cur_overflow;
  return 0;
}

// Exported variant of VaridTemplateApply() which appends to a buffer.
// Probable todo: pull out the common parts of the functions.
char* VaridTemplateWrite(const VaridTemplate* vtp, const char* ref_start, const char* alt1_start, uint32_t cur_bp, uint32_t ref_token_slen, uint32_t extra_alt_ct, uint32_t alt_token_slen, uint32_t* max_overflow_slenp, char* dst) {
  uint32_t insert_slens[4];
  const uint32_t allele_flags = vtp->allele_flags;
  const uint32_t new_id_max_allele_slen = vtp->new_id_max_allele_slen;
  uint32_t id_slen = UintSlen(cur_bp);
  insert_slens[0] = vtp->chr_slen;
  insert_slens[1] = id_slen;
  id_slen += vtp->base_len;
  const uint32_t ref_or_1_exists = (allele_flags & kfVaridTemplateAlleleRefOr1)? 1 : 0;
  uint32_t ref_slen = 0;
  uint32_t cur_max_overflow_slen = 0;
  const char* tmp_allele_ptrs[2];
  tmp_allele_ptrs[0] = nullptr;
  tmp_allele_ptrs[1] = nullptr;
  if (ref_or_1_exists) {
    ref_slen = ref_token_slen;
    if (ref_slen > new_id_max_allele_slen) {
      cur_max_overflow_slen = ref_slen;
      ref_slen = new_id_max_allele_slen;
    }
    insert_slens[2] = ref_slen;
    id_slen += ref_slen;
    tmp_allele_ptrs[0] = ref_start;
  }
  if (allele_flags & (kfVaridTemplateAlleleAltOr2 | kfVaridTemplateAlleleAsciiOrder)) {
    uint32_t alt1_slen;
    if (!extra_alt_ct) {
      alt1_slen = alt_token_slen;
    } else {
      alt1_slen = AdvToDelim(alt1_start, ',') - alt1_start;
    }
    if (alt1_slen > new_id_max_allele_slen) {
      if (alt1_slen > cur_max_overflow_slen) {
        cur_max_overflow_slen = alt1_slen;
      }
      alt1_slen = new_id_max_allele_slen;
    }
    id_slen += alt1_slen;
    if (!(allele_flags & kfVaridTemplateAlleleAsciiOrder)) {
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
  if (cur_max_overflow_slen) {
    if (cur_max_overflow_slen > *max_overflow_slenp) {
      *max_overflow_slenp = cur_max_overflow_slen;
    }
    const uint32_t overflow_substitute_blen = vtp->overflow_substitute_blen;
    if (overflow_substitute_blen) {
      return memcpya(dst, vtp->missing_id_match, overflow_substitute_blen - 1);
    }
  }
  char* id_iter = dst;
  const uint32_t insert_ct = vtp->insert_ct;
  char* insert_ptrs[4];

  // maybe-uninitialized warnings
  insert_ptrs[0] = nullptr;
  insert_ptrs[1] = nullptr;

  insert_ptrs[2] = nullptr;
  insert_ptrs[3] = nullptr;

  for (uint32_t insert_idx = 0; insert_idx != insert_ct; ++insert_idx) {
    id_iter = memcpya(id_iter, vtp->segs[insert_idx], vtp->seg_lens[insert_idx]);
    const uint32_t cur_insert_type = vtp->insert_types[insert_idx];
    insert_ptrs[cur_insert_type] = id_iter;
    id_iter = &(id_iter[insert_slens[cur_insert_type]]);
  }
  char* id_end = memcpya(id_iter, vtp->segs[insert_ct], vtp->seg_lens[insert_ct]);

  memcpy(insert_ptrs[0], vtp->chr_output_name_buf, insert_slens[0]);
  u32toa(cur_bp, insert_ptrs[1]);
  if (insert_ct > 2) {
    for (uint32_t insert_type_idx = 2; insert_type_idx != 4; ++insert_type_idx) {
      char* cur_dst = insert_ptrs[insert_type_idx];
      if (cur_dst == nullptr) {
        continue;
      }
      memcpy(cur_dst, tmp_allele_ptrs[insert_type_idx - 2], insert_slens[insert_type_idx]);
    }
  }
  return id_end;
}

uint32_t VaridWorstCaseSlen(const VaridTemplate* vtp, uint32_t max_chr_slen, uint32_t max_allele_slen) {
  const VaridTemplateAlleleFlags allele_flags = vtp->allele_flags;
  uint32_t allele_ct = (allele_flags & kfVaridTemplateAlleleRefOr1)? 1 : 0;
  allele_ct += (allele_flags & kfVaridTemplateAlleleAltOr2)? 1 : 0;
  // +10 for base-pair coordinate
  return (max_allele_slen * allele_ct + vtp->base_len + max_chr_slen + 10);
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

uint32_t PrInInfo(uint32_t info_slen, char* info_token) {
  if (memequal_sk(info_token, "PR") && ((info_slen == 2) || (info_token[2] == ';'))) {
    return 1;
  }
  if (memequal_sk(&(info_token[S_CAST(int32_t, info_slen) - 3]), ";PR")) {
    return 1;
  }
  info_token[info_slen] = '\0';
  char* first_info_end = strchr(info_token, ';');
  return first_info_end && (strstr(info_token, ";PR;") != nullptr);
}

char* InfoPrStart(uint32_t info_slen, char* info_token) {
  if (memequal_sk(info_token, "PR") && ((info_slen == 2) || (info_token[2] == ';'))) {
    return info_token;
  }
  if (memequal_sk(&(info_token[S_CAST(int32_t, info_slen) - 3]), ";PR")) {
    return &(info_token[info_slen - 2]);
  }
  // bugfix (29 Aug 2023): need this function to be nondestructive
  const char token_end_char = info_token[info_slen];
  info_token[info_slen] = '\0';
  char* first_info_end = strchr(info_token, ';');
  if (!first_info_end) {
    info_token[info_slen] = token_end_char;
    return nullptr;
  }
  char* pr_prestart = strstr(first_info_end, ";PR;");
  info_token[info_slen] = token_end_char;
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

// Returns 1 iff all keys are found.
// (could use InfoFilter for this, but this is simple enough)
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

// Returns 1 iff all keys are absent.
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

struct InfoExprStruct;

typedef struct InfoExprNStruct {
  double value;
} InfoExprN;

typedef struct InfoExprSStruct {
  const char* str_value;
  uint32_t slen;
} InfoExprS;

typedef struct InfoExprJctStruct {
  struct InfoExprStruct* children[2];
} InfoExprJct;

typedef union {
  InfoExprN n;
  InfoExprS s;
  InfoExprJct jct;
} InfoExprArgs;

typedef struct InfoExprStruct {
  CmpExprType etype;  // normalized to remove Noteq, StrNoteq
  uint32_t kidx;
  uint32_t negate;
  InfoExprArgs args;
} InfoExpr;

typedef struct InfoFilterStruct {
  InfoExpr expr;

  const char* const* prekeys;  // ;<key>;
  uint32_t* keyeq_slens;  // includes =, does not include ;

  const char** cur_str_values;  // nullptr = unparsed, ';' = nonexistent
  double* cur_values;  // INFINITY = unparsed, -INFINITY = invalid

  uint32_t key_ct;
} InfoFilter;

PglErr InfoFilterFirstPass(char* arena_end_c, const CmpExpr* cmp_expr, char* prekey_base, char** prekey_alloc_ptr, uint32_t* key_ct_ptr, uintptr_t* val_blen_sum_ptr) {
  PglErr reterr = kPglRetSuccess;
  {
    const CmpExprType etype = cmp_expr->etype;
    if (!CmpExprIsJct(etype)) {
      assert(etype != kCmpExprTypeNull);
      const char* key;
      // in KS case, prepare to make a local copy of the string
      const char* str_value = nullptr;
      if (etype == kCmpExprTypeExists) {
        key = cmp_expr->args.k.key;
      } else if (etype <= kCmpExprTypeGeq) {
        key = cmp_expr->args.kn.key;
      } else {
        key = cmp_expr->args.ks.key;
        str_value = cmp_expr->args.ks.str_value;
        const char* str_value_end_or_invalid = strchrnul2(str_value, ';', '=');
        if (*str_value_end_or_invalid) {
          if (unlikely(*str_value_end_or_invalid == '=')) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --{extract,exclude}-if-info value '%s' ('=' prohibited).\n", str_value);
            goto InfoFilterFirstPass_ret_INVALID_CMDLINE_WW;
          } else if (unlikely(str_value_end_or_invalid[1])) {
            logerrprintfww("Error: Invalid --{extract,exclude}-if-info value '%s' (semicolon prohibited, unless by itself to specify empty-string).\n", str_value);
            goto InfoFilterFirstPass_ret_INVALID_CMDLINE_WW;
          }
          // Empty-string special case.
          str_value = &(str_value[1]);
        }
        // don't bother to deduplicate str_values for now
        *val_blen_sum_ptr += 1 + S_CAST(uintptr_t, str_value_end_or_invalid - str_value);
      }
      const uint32_t key_slen = strlen(key);
      const uint32_t prev_key_ct = *key_ct_ptr;
      char* prekey_iter = prekey_base;
      uint32_t kidx = 0;
      // ok for this to be quadratic
      for (; kidx != prev_key_ct; ++kidx) {
        const uint32_t prekey_slen = strlen(prekey_iter);
        if ((key_slen == prekey_slen - 2) && memequal(key, &(prekey_iter[1]), key_slen)) {
          break;
        }
        prekey_iter = &(prekey_iter[prekey_slen + 1]);
      }
      if (kidx == prev_key_ct) {
        assert(*prekey_alloc_ptr == prekey_iter);
        if (unlikely(S_CAST(uintptr_t, arena_end_c - prekey_iter) < key_slen + 3)) {
          goto InfoFilterFirstPass_ret_NOMEM;
        }
        *prekey_iter++ = ';';
        prekey_iter = memcpya(prekey_iter, key, key_slen);
        prekey_iter = memcpya(prekey_iter, "=", 2);
        *prekey_alloc_ptr = prekey_iter;
        *key_ct_ptr = prev_key_ct + 1;
      }
    } else {
      const uint32_t child_ct = 1 + (etype != kCmpExprTypeNot);
      for (uint32_t child_idx = 0; child_idx != child_ct; ++child_idx) {
        reterr = InfoFilterFirstPass(arena_end_c, cmp_expr->args.jct.children[child_idx], prekey_base, prekey_alloc_ptr, key_ct_ptr, val_blen_sum_ptr);
        if (unlikely(reterr)) {
          goto InfoFilterFirstPass_ret_1;
        }
      }
    }
  }
  while (0) {
  InfoFilterFirstPass_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  InfoFilterFirstPass_ret_INVALID_CMDLINE_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetInvalidCmdline;
    break;
  }
 InfoFilterFirstPass_ret_1:
  return reterr;
}

int32_t bsearch_prekeys(const char* key, const char* const* prekeys, uint32_t key_slen, uint32_t end_idx) {
  uint32_t start_idx = 0;
  while (start_idx < end_idx) {
    const uint32_t mid_idx = (start_idx + end_idx) / 2;
    const char* candidate_key_start = &(prekeys[mid_idx][1]);
    int32_t ii = Memcmp(key, &(prekeys[mid_idx][1]), key_slen);
    if (ii == 0) {
      ii = S_CAST(int32_t, '=') - S_CAST(int32_t, ctou32(candidate_key_start[key_slen]));
      if (ii == 0) {
        return mid_idx;
      }
    }
    if (ii > 0) {
      start_idx = mid_idx + 1;
    } else {
      end_idx = mid_idx;
    }
  }
  return -1;
}

// exprp->negate assumed to be initialized, other fields need not be
PglErr InfoFilterSecondPass(const unsigned char* arena_end, const CmpExpr* cmp_expr, const char* const* prekeys, uint32_t key_ct, unsigned char** arena_base_ptr, char** str_value_alloc_ptr, InfoExpr* exprp) {
  PglErr reterr = kPglRetSuccess;
  {
    CmpExprType etype;
    while (1) {
      etype = cmp_expr->etype;
      if (etype != kCmpExprTypeNot) {
        break;
      }
      exprp->negate = 1 - exprp->negate;
      cmp_expr = cmp_expr->args.jct.children[0];
    }
    exprp->etype = etype;
    if (!CmpExprIsJct(etype)) {
      const char* key;
      if (etype == kCmpExprTypeExists) {
        key = cmp_expr->args.k.key;
      } else if (etype <= kCmpExprTypeGeq) {
        key = cmp_expr->args.kn.key;
      } else {
        key = cmp_expr->args.ks.key;
      }
      int32_t ii = bsearch_prekeys(key, prekeys, strlen(key), key_ct);
      assert(ii != -1);
      exprp->kidx = ii;
      if (etype <= kCmpExprTypeGeq) {
        if (etype != kCmpExprTypeExists) {
          exprp->args.n.value = cmp_expr->args.kn.value;
          if (etype == kCmpExprTypeNoteq) {
            exprp->etype = kCmpExprTypeEq;
            exprp->negate = 1 - exprp->negate;
          }
        }
      } else {
        const char* str_value = cmp_expr->args.ks.str_value;
        if (str_value[0] == ';') {
          // Empty-string special case.  Previously validated.
          str_value = &(str_value[1]);
        }
        const uint32_t str_value_blen = 1 + strlen(str_value);
        memcpy(*str_value_alloc_ptr, str_value, str_value_blen);
        *str_value_alloc_ptr += str_value_blen;
        exprp->args.s.str_value = str_value;
        if (etype == kCmpExprTypeStrNoteq) {
          exprp->etype = kCmpExprTypeStrEq;
          exprp->negate = 1 - exprp->negate;
        }
      }
    } else {
      if (unlikely(S_CAST(uintptr_t, arena_end - (*arena_base_ptr)) < RoundUpPow2(2 * sizeof(InfoExpr), kCacheline))) {
        goto InfoFilterSecondPass_ret_NOMEM;
      }
      InfoExpr* children = R_CAST(InfoExpr*, *arena_base_ptr);
      (*arena_base_ptr) += kCacheline;
      for (uint32_t child_idx = 0; child_idx != 2; ++child_idx) {
        InfoExpr* cur_child = &(children[child_idx]);
        exprp->args.jct.children[child_idx] = cur_child;
        cur_child->negate = 0;
        reterr = InfoFilterSecondPass(arena_end, cmp_expr->args.jct.children[child_idx], prekeys, key_ct, arena_base_ptr, str_value_alloc_ptr, cur_child);
        if (unlikely(reterr)) {
          goto InfoFilterSecondPass_ret_1;
        }
      }
    }
  }
  while (0) {
  InfoFilterSecondPass_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  }
 InfoFilterSecondPass_ret_1:
  return reterr;
}

PglErr InfoFilterInit(unsigned char* arena_end, const CmpExpr* cmp_expr, unsigned char** arena_base_ptr, InfoFilter* filterp) {
  PglErr reterr = kPglRetSuccess;
  {
    // First pass: Copy unique keys to arena (adding ; in front and = after),
    //             count number of distinct keys, validate str_values and
    //             reserve space for them.
    // Intermission: Allocate and initialize filterp->prekeys and
    //               filterp->keyeq_slens.  Allocate cur_values and
    //               cur_str_values.
    // Second pass: Allocate and initialize InfoExpr nodes and str_values.
    char* arena_end_c = R_CAST(char*, arena_end);
    char* prekey_base = R_CAST(char*, *arena_base_ptr);
    char* prekey_and_str_value_alloc = prekey_base;
    uintptr_t val_blen_sum = 0;
    uint32_t key_ct = 0;
    reterr = InfoFilterFirstPass(arena_end_c, cmp_expr, prekey_base, &prekey_and_str_value_alloc, &key_ct, &val_blen_sum);
    if (unlikely(reterr)) {
      goto InfoFilterInit_ret_1;
    }
    if (unlikely(S_CAST(uintptr_t, arena_end_c - prekey_and_str_value_alloc) < val_blen_sum)) {
      goto InfoFilterInit_ret_NOMEM;
    }
    ArenaBaseSet(&(prekey_and_str_value_alloc[val_blen_sum]), arena_base_ptr);
    const char** prekeys_mutable;
    uint32_t* keyeq_slens;
    if (unlikely(arena_alloc_kcp(arena_end, key_ct, arena_base_ptr, &prekeys_mutable) ||
                 arena_alloc_u32(arena_end, key_ct, arena_base_ptr, &keyeq_slens) ||
                 arena_alloc_kcp(arena_end, key_ct, arena_base_ptr, &(filterp->cur_str_values)) ||
                 arena_alloc_d(arena_end, key_ct, arena_base_ptr, &(filterp->cur_values)))) {
      goto InfoFilterInit_ret_NOMEM;
    }
    filterp->keyeq_slens = keyeq_slens;
    char* prekey_iter = prekey_base;
    for (uint32_t kidx = 0; kidx != key_ct; ++kidx) {
      prekeys_mutable[kidx] = prekey_iter;
      const uint32_t prekey_slen = strlen(prekey_iter);
      prekey_iter = &(prekey_iter[prekey_slen + 1]);
    }
    StrptrArrSortOverread(key_ct, prekeys_mutable);
    const char* const* prekeys = prekeys_mutable;
    filterp->prekeys = prekeys;
    for (uint32_t kidx = 0; kidx != key_ct; ++kidx) {
      keyeq_slens[kidx] = strlen(prekeys[kidx]) - 1;
    }
    filterp->key_ct = key_ct;
    filterp->expr.negate = 0;

    reterr = InfoFilterSecondPass(arena_end, cmp_expr, prekeys, key_ct, arena_base_ptr, &prekey_and_str_value_alloc, &(filterp->expr));
    if (unlikely(reterr)) {
      goto InfoFilterInit_ret_1;
    }
  }
  while (0) {
  InfoFilterInit_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  }
 InfoFilterInit_ret_1:
  return reterr;
}

static_assert(kCmpExprTypeOr == kCmpExprTypeAnd + 1, "InfoConditionSatisfiedInternal() must be updated.");
uint32_t InfoConditionSatisfiedInternal(const InfoExpr* exprp, const char* info_token, InfoFilter* filterp) {
  const CmpExprType etype = exprp->etype;
  const uint32_t negate = exprp->negate;
  if (!CmpExprIsJct(etype)) {
    const uint32_t kidx = exprp->kidx;
    const char** str_value_ptr = &(filterp->cur_str_values[kidx]);
    const char* str_value = *str_value_ptr;
    if (!str_value) {
      const char* prekey = filterp->prekeys[kidx];
      const uint32_t keyeq_slen = filterp->keyeq_slens[kidx];
      if (memequal(info_token, &(prekey[1]), keyeq_slen)) {
        str_value = &(info_token[keyeq_slen]);
      } else {
        str_value = strstr(info_token, prekey);
        if (!str_value) {
          // ';'
          str_value = &(g_one_char_strs[118]);
        } else {
          str_value = &(str_value[keyeq_slen + 1]);
        }
      }
      *str_value_ptr = str_value;
    }
    if (*str_value == ';') {
      // key doesn't exist
      return negate;
    }
    if (etype <= kCmpExprTypeGeq) {
      if (etype == kCmpExprTypeExists) {
        return !negate;
      }
      double* cur_value_ptr = &(filterp->cur_values[kidx]);
      double value = *cur_value_ptr;
      if (value == S_CAST(double, INFINITY)) {
        const char* scan_end = ScanadvDouble(str_value, &value);
        if ((!scan_end) || ((*scan_end != ';') && (*scan_end))) {
          value = S_CAST(double, -INFINITY);
        }
        *cur_value_ptr = value;
      }
      if (value == S_CAST(double, -INFINITY)) {
        // value cannot be parsed as a number
        return negate;
      }
      const double cmp_value = exprp->args.n.value;
      if (etype == kCmpExprTypeEq) {
        return (value == cmp_value) != negate;
      } else if (etype < kCmpExprTypeEq) {
        if (etype == kCmpExprTypeLe) {
          return (value < cmp_value) != negate;
        } else {
          return (value <= cmp_value) != negate;
        }
      } else {
        if (etype == kCmpExprTypeGe) {
          return (value > cmp_value) != negate;
        } else {
          return (value >= cmp_value) != negate;
        }
      }
    } else {
      const char* cmp_str = exprp->args.s.str_value;
      const uint32_t cmp_slen = exprp->args.s.slen;
      const uint32_t is_match = memequal(str_value, cmp_str, cmp_slen) && (!(str_value[cmp_slen] && (str_value[cmp_slen] != ';')));
      return is_match != negate;
    }
  }
  const uint32_t first_result = InfoConditionSatisfiedInternal(exprp->args.jct.children[0], info_token, filterp);
  // short-circuit: if kCmpExprTypeAnd, can return 0; if Or, can return 1
  if (first_result == etype - kCmpExprTypeAnd) {
    return first_result != negate;
  }
  return InfoConditionSatisfiedInternal(exprp->args.jct.children[1], info_token, filterp) != negate;
}

uint32_t InfoConditionSatisfied(const char* info_token, InfoFilter* filterp) {
  const uint32_t key_ct = filterp->key_ct;
  ZeroPtrArr(key_ct, filterp->cur_str_values);
  double* cur_values = filterp->cur_values;
  for (uint32_t kidx = 0; kidx != key_ct; ++kidx) {
    cur_values[kidx] = S_CAST(double, INFINITY);
  }
  return InfoConditionSatisfiedInternal(&(filterp->expr), info_token, filterp);
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
  const uint32_t par1_end = LowerBoundConstrainedNonemptyU32(variant_bps, orig_x_start, orig_x_end, splitpar_bound1 + 1);
  const uint32_t par2_start = LowerBoundConstrainedNonemptyU32(variant_bps, par1_end, orig_x_end, splitpar_bound2);
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

PglErr LoadPvar(const char* pvarname, const char* var_filter_exceptions_flattened, const char* varid_template_str, const char* varid_multi_template_str, const char* varid_multi_nonsnp_template_str, const char* missing_varid_match, const char* require_info_flattened, const char* require_no_info_flattened, const CmpExpr* extract_if_info_exprp, const CmpExpr* exclude_if_info_exprp, MiscFlags misc_flags, PvarPsamFlags pvar_psam_flags, LoadFilterLogFlags load_filter_log_flags, uint32_t xheader_needed, uint32_t qualfilter_needed, float var_min_qual, uint32_t splitpar_bound1, uint32_t splitpar_bound2, uint32_t new_variant_id_max_allele_slen, uint32_t snps_only_just_acgt, uint32_t split_chr_ok, uint32_t filter_min_allele_ct, uint32_t filter_max_allele_ct, char input_missing_geno_char, uint32_t max_thread_ct, ChrInfo* cip, uint32_t* max_variant_id_slen_ptr, uint32_t* info_reload_slen_ptr, UnsortedVar* vpos_sortstatus_ptr, char** xheader_ptr, uintptr_t** variant_include_ptr, uint32_t** variant_bps_ptr, char*** variant_ids_ptr, uintptr_t** allele_idx_offsets_ptr, const char*** allele_storage_ptr, uintptr_t** qual_present_ptr, float** quals_ptr, uintptr_t** filter_present_ptr, uintptr_t** filter_npass_ptr, char*** filter_storage_ptr, uintptr_t** nonref_flags_ptr, double** variant_cms_ptr, ChrIdx** chr_idxs_ptr, uint32_t* raw_variant_ct_ptr, uint32_t* variant_ct_ptr, uint32_t* neg_bp_seen_ptr, uint32_t* max_allele_ct_ptr, uint32_t* max_allele_slen_ptr, uintptr_t* xheader_blen_ptr, InfoFlags* info_flags_ptr, uint32_t* max_filter_slen_ptr) {
  // chr_info, max_variant_id_slen, and info_reload_slen are in/out; just
  // outparameters after them.  (Due to its large size in some VCFs, INFO is
  // not kept in memory for now.  This has a speed penalty, of course; maybe
  // it's worthwhile to conditionally load it later.)

  // allele_idx_offsets currently assumed to be initialized to nullptr

  // should handle raw_variant_ct == 0 properly

  // possible todo: optionally skip allele code loading
  // probable todo: conditionally compute the value previously known as
  //   INFO/END.  (does this allow the CNV module to be unified with the rest
  //   of the program?)  but this will probably wait until I need to analyze
  //   some sort of CNV data, and that day keeps getting postponed... for now,
  //   the BCF exporter performs its own parsing of INFO/END to have a decent
  //   shot at filling each variant's rlen field correctly.
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
    uint32_t info_pr_exists = 0;
    uint32_t info_pr_nonflag_exists = 0;
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
      if (StrStartsWithUnsafe(line_start, "##INFO=<")) {
        line_iter = &(line_iter[strlen("##INFO=<")]);
        char* idval;
        uint32_t id_slen;
        if (unlikely(HkvlineId(&line_iter, &idval, &id_slen))) {
          goto LoadPvar_ret_MALFORMED_HEADER_LINE;
        }
        if (strequal_k(idval, "PR", id_slen)) {
          if (unlikely(info_pr_exists || info_pr_nonflag_exists)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Duplicate INFO/PR header line in %s.\n", pvarname);
            goto LoadPvar_ret_MALFORMED_INPUT_WW;
          }
          // Could also enforce numstr == "0".
          char* typestr;
          uint32_t type_slen;
          if (unlikely(HkvlineFind(line_iter, "Type", &typestr, &type_slen))) {
            goto LoadPvar_ret_MALFORMED_HEADER_LINE;
          }
          info_pr_nonflag_exists = !strequal_k(typestr, "Flag", type_slen);
          info_pr_exists = 1 - info_pr_nonflag_exists;
          if (info_pr_nonflag_exists) {
            logerrprintfww("Warning: Header line %" PRIuPTR " of %s has an unexpected definition of INFO/PR. This interferes with a few merge and liftover operations.\n", line_idx, pvarname);
          }
        } else if (!info_nonpr_present) {
          info_nonpr_present = 1;
        }
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
          char* line_lf = AdvToDelim(line_start, '\n');
          uint32_t line_slen = line_lf - line_start;
          if (line_start[line_slen - 1] == '\r') {
            --line_slen;
          }
          if (unlikely(S_CAST(uintptr_t, R_CAST(char*, rlstream_start) - xheader_end) < line_slen + 2)) {
            goto LoadPvar_ret_NOMEM;
          }
          xheader_end = memcpya(xheader_end, line_start, line_slen);
          AppendBinaryEoln(&xheader_end);
          line_iter = line_lf;
        }
      }
      line_iter = AdvPastDelim(line_iter, '\n');
    }
    if (xheader_end) {
      *xheader_ptr = R_CAST(char*, bigstack_mark);
      *xheader_blen_ptr = xheader_end - (*xheader_ptr);
      BigstackBaseSet(xheader_end);
    }
    FinalizeChrset(load_filter_log_flags, cip);
    const char** allele_storage = R_CAST(const char**, g_bigstack_base);
    const char** allele_storage_iter = allele_storage;

    uint32_t col_skips[8];
    uint32_t col_types[8];
    uint32_t relevant_postchr_col_ct = 5;
    uint32_t no_multiallelic_allowed = 0;
    uint32_t alt_col_idx = 4;
    uint32_t load_qual_col = 0;
    uint32_t load_filter_col = 0;
    uint32_t info_col_present = 0;
    uint32_t cm_col_present = 0;
    if (line_start[0] == '#') {
      *info_flags_ptr = S_CAST(InfoFlags, (info_pr_exists * kfInfoPrFlagPresent) | (info_pr_nonflag_exists * kfInfoPrNonflagPresent) | (info_nonpr_present * kfInfoNonprPresent));
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
            if (memequal_sk(linebuf_iter, "POS")) {
              cur_col_type = 0;
            } else if (memequal_sk(linebuf_iter, "REF")) {
              cur_col_type = 2;
            } else if (memequal_sk(linebuf_iter, "ALT")) {
              cur_col_type = 3;
              alt_col_idx = col_idx;
            } else {
              continue;
            }
          } else if (token_slen == 2) {
            if (memequal_sk(linebuf_iter, "ID")) {
              cur_col_type = 1;
            } else if (memequal_sk(linebuf_iter, "CM")) {
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
          if (memequal_sk(linebuf_iter, "FILTER")) {
            load_filter_col = 2 * ((pvar_psam_flags & (kfPvarColMaybefilter | kfPvarColFilter)) || qualfilter_needed) + ((load_filter_log_flags / kfLoadFilterLogVarFilter) & 1);
            if (!load_filter_col) {
              continue;
            }
            cur_col_type = 5;
          } else if (memequal_sk(linebuf_iter, "FORMAT")) {
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
      no_multiallelic_allowed = 1;
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
    const char** allele_storage_limit = R_CAST(const char**, &(rlstream_start[-S_CAST(intptr_t, kLoadPvarBlockSize * 2 * sizeof(intptr_t))]));

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
    InfoFilter info_filter;
    info_filter.expr.etype = kCmpExprTypeNull;
    {
      const CmpExprType etype_incl = extract_if_info_exprp->etype;
      const CmpExprType etype_excl = exclude_if_info_exprp->etype;
      if (etype_incl || etype_excl) {
        // todo: also print warning (or optionally error out?) if header line
        // is missing or doesn't match type expectation
        // (same for --require-info)
        const CmpExpr* cmp_expr;
        CmpExpr not_slot;
        CmpExpr and_slot;
        if (!etype_excl) {
          cmp_expr = extract_if_info_exprp;
        } else {
          not_slot.etype = kCmpExprTypeNot;
          not_slot.args.jct.children[0] = K_CAST(CmpExpr*, exclude_if_info_exprp);
          if (etype_incl) {
            and_slot.etype = kCmpExprTypeAnd;
            and_slot.args.jct.children[0] = K_CAST(CmpExpr*, extract_if_info_exprp);
            and_slot.args.jct.children[1] = &not_slot;
            cmp_expr = &and_slot;
          } else {
            cmp_expr = &not_slot;
          }
        }
        reterr = InfoFilterInit(tmp_alloc_end, cmp_expr, &tmp_alloc_base, &info_filter);
        if (unlikely(reterr)) {
          goto LoadPvar_ret_1;
        }
      }
    }
    if (!info_col_present) {
      if (unlikely(require_info_flattened)) {
        logerrputs("Error: --require-info used on a variant file with no INFO column.\n");
        goto LoadPvar_ret_INCONSISTENT_INPUT;
      }
      if (unlikely(extract_if_info_exprp->etype)) {
        logerrputs("Error: --extract-if-info used on a variant file with no INFO column.\n");
        goto LoadPvar_ret_INCONSISTENT_INPUT;
      }
      info_pr_exists = 0;
      info_reload_slen = 0;
    } else if ((!info_pr_exists) && (!info_reload_slen) && (!info_existp) && (!info_nonexistp) && (!info_filter.expr.etype)) {
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
    char* chr_output_name_buf = nullptr;
    VaridTemplate* varid_templatep = nullptr;
    VaridTemplate* varid_multi_templatep = nullptr;
    VaridTemplate* varid_multi_nonsnp_templatep = nullptr;
    uint32_t missing_varid_blen = 0;
    uint32_t missing_varid_match_slen = 0;
    if (varid_template_str) {
      if (unlikely(VaridInitAll(tmp_alloc_end, varid_template_str, varid_multi_template_str, varid_multi_nonsnp_template_str, misc_flags, new_variant_id_max_allele_slen, &tmp_alloc_base, &missing_varid_match, &chr_output_name_buf, &varid_templatep, &varid_multi_templatep, &varid_multi_nonsnp_templatep, &missing_varid_blen, &missing_varid_match_slen))) {
        goto LoadPvar_ret_NOMEM;
      }
    }

    // prevent later return-array allocations from overlapping with temporary
    // storage
    g_bigstack_end = tmp_alloc_base;

    // prevent variant_id_htable_find from breaking
    if (R_CAST(const char*, tmp_alloc_end) > (&(g_one_char_strs[512 - kMaxIdSlen]))) {
      tmp_alloc_end = R_CAST(unsigned char*, K_CAST(char*, &(g_one_char_strs[512 - kMaxIdSlen])));
    }
    const uint32_t prohibit_extra_chrs = (misc_flags / kfMiscProhibitExtraChr) & 1;
    const uint32_t merge_par = ((misc_flags & (kfMiscMergePar | kfMiscMergeX)) != 0);
    const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
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

    // A/C/G/T/a/c/g/t/missing nonzero
    // A+T = C+G = 7, no other biallelic way to get sum of 7
    uint8_t acgtm_table[256] = {
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 2, 0, 3, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 2, 0, 3, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    acgtm_table[ctou32(input_missing_geno_char)] = 1;

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

    const uint32_t exclude_palindromic_snps = (load_filter_log_flags / kfLoadFilterLogExcludePalindromicSnps) & 1;
    const uint32_t snps_only = (load_filter_log_flags / kfLoadFilterLogSnpsOnly) & 1;
    // only allocated when necessary
    // if we want to scale this approach to more fields, we'll need to add a
    // few pointers to the start of each block.  right now, we force cur_cms[]
    // to be allocated before cur_chr_idxs[] when both are present, but this
    // is error-prone.
    uint32_t at_least_one_npass_filter = 0;
    uint32_t at_least_one_nzero_cm = 0;
    uint32_t neg_bp_seen = 0;
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
      if (unlikely(raw_variant_ct == kPglMaxVariantCt)) {
        logerrputs("Error: " PROG_NAME_STR " does not support more than 2^31 - 3 variants.  We recommend using\nother software for very deep studies of small numbers of genomes.\n");
        goto LoadPvar_ret_MALFORMED_INPUT;
      }
#endif
      const uint32_t variant_idx_lowbits = raw_variant_ct % kLoadPvarBlockSize;
      if (!variant_idx_lowbits) {
        if (unlikely((S_CAST(uintptr_t, tmp_alloc_end - tmp_alloc_base) <=
                     kLoadPvarBlockSize *
                      (sizeof(int32_t) +
                       2 * sizeof(intptr_t) +
                       at_least_one_nzero_cm * sizeof(double)) +
                      is_split_chr * sizeof(ChrIdx) +
                      (1 + info_pr_exists) * (kLoadPvarBlockSize / CHAR_BIT) +
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
        if (info_pr_exists) {
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
      reterr = GetOrAddChrCodeDestructive(".pvar file", line_idx, prohibit_extra_chrs, line_iter, linebuf_iter, cip, &cur_chr_code);
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
              snprintf(g_logbuf, kLogbufSize, "Error: %s has a split chromosome. Use --make-pgen + --sort-vars (without other simultaneous commands) to remedy this.\n", pvarname);
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
      if (IsSet(chr_mask, cur_chr_code) || info_pr_exists) {
        linebuf_iter = TokenLex(linebuf_iter, col_types, col_skips, relevant_postchr_col_ct, token_ptrs, token_slens);
        if (unlikely(!linebuf_iter)) {
          goto LoadPvar_ret_MISSING_TOKENS;
        }

        extra_alt_ct = CountByte(token_ptrs[3], ',', token_slens[3]);
        if (extra_alt_ct > max_extra_alt_ct) {
          if (extra_alt_ct >= kPglMaxAltAlleleCt) {
            logerrprintfww("Error: Too many ALT alleles on line %" PRIuPTR " of %s. (This " PROG_NAME_STR " build is limited to " PGL_MAX_ALT_ALLELE_CT_STR ".)\n", line_idx, pvarname);
            reterr = kPglRetNotYetSupported;
            goto LoadPvar_ret_1;
          }
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
          if (info_pr_exists) {
            // always load all nonref_flags entries so (i) --ref-from-fa +
            // --make-just-pvar works and (ii) they can be compared against
            // the .pgen.
            if ((memequal_sk(info_token, "PR") && ((info_slen == 2) || (info_token[2] == ';'))) || memequal_sk(&(info_token[S_CAST(int32_t, info_slen) - 3]), ";PR")) {
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
          if (info_filter.expr.etype) {
            if (!InfoConditionSatisfied(info_token, &info_filter)) {
              goto LoadPvar_skip_variant;
            }
          }
        }
        // POS
        // could error out on floating-point number?
        int32_t cur_bp;
        if (unlikely(ScanIntAbsDefcap(token_ptrs[0], &cur_bp))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Invalid bp coordinate on line %" PRIuPTR " of %s.\n", line_idx, pvarname);
          goto LoadPvar_ret_MALFORMED_INPUT_WW;
        }

        if (cur_bp < 0) {
          neg_bp_seen = 1;
          goto LoadPvar_skip_variant;
        }

        // QUAL
        if (load_qual_col) {
          const char* qual_token = token_ptrs[4];
          if ((qual_token[0] != '.') || (qual_token[1] > ' ')) {
            float cur_qual;
            if (unlikely(ScanFloatAllowInf(qual_token, &cur_qual))) {
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
          if (snps_only_just_acgt) {
            // just-acgt
            if (!acgtm_table[ctou32(token_ptrs[2][0])]) {
              goto LoadPvar_skip_variant;
            }
            for (uint32_t uii = 0; uii <= extra_alt_ct; ++uii) {
              if (!acgtm_table[ctou32(linebuf_iter[2 * uii])]) {
                goto LoadPvar_skip_variant;
              }
            }
          }
        }
        if (exclude_palindromic_snps) {
          if ((acgtm_table[ctou32(token_ptrs[2][0])] + acgtm_table[ctou32(linebuf_iter[0])] == 7) && (remaining_alt_char_ct == 1) && (token_slens[2] == 1)) {
            goto LoadPvar_skip_variant;
          }
        }

        // also handle --min-alleles/--max-alleles here
        if (filter_min_allele_ct || (filter_max_allele_ct <= kPglMaxAltAlleleCt)) {
          uint32_t allele_ct = extra_alt_ct + 2;
          if (!extra_alt_ct) {
            // allele_ct == 1 or 2 for filtering purposes, depending on
            // whether ALT allele matches a missing code.  (Don't see a reason
            // to subtract 1 for missing REF.)
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
                  if (bsearch_strbox(filter_token_iter, sorted_fexcepts, cur_slen, max_fexcept_blen, fexcept_ct) == -1) {
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
        char ref_allele_first_char = ref_allele[0];
        uint32_t missing_allele_ct = 0;
        if (ref_slen == 1) {
          if (ref_allele_first_char == input_missing_geno_char) {
            ref_allele_first_char = '.';
          }
          missing_allele_ct += (ref_allele_first_char == '.');
          *allele_storage_iter = &(g_one_char_strs[2 * ctou32(ref_allele_first_char)]);
        } else {
          // sanity check: prohibit comma
          // (could also do it in ref_slen == 1 case, but that's much less
          // likely to happen by accident)
          if (unlikely(memchr(ref_allele, ',', ref_slen) != nullptr)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid REF allele on line %" PRIuPTR " of %s.\n", line_idx, pvarname);
            goto LoadPvar_ret_MALFORMED_INPUT_WW;
          }
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
              const char geno_char = linebuf_iter[0];
              if (unlikely((geno_char == '.') || (geno_char == input_missing_geno_char))) {
                goto LoadPvar_ret_MULTIALLELIC_MISSING_ALLELE_CODE;
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
          if (geno_char == '.') {
            ++missing_allele_ct;
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
        if (unlikely(missing_allele_ct && extra_alt_ct)) {
          goto LoadPvar_ret_MULTIALLELIC_MISSING_ALLELE_CODE;
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
      const uint32_t new_variant_id_overflow_missing = (misc_flags / kfMiscNewVarIdOverflowMissing) & 1;
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
    if (unlikely(bigstack_alloc_w(raw_variant_ctl, variant_include_ptr) ||
                 bigstack_alloc_u32(raw_variant_ct, variant_bps_ptr) ||
                 bigstack_alloc_cp(raw_variant_ct, variant_ids_ptr))) {
      goto LoadPvar_ret_NOMEM;
    }
    uintptr_t* qual_present = nullptr;
    float* quals = nullptr;
    if (load_qual_col > 1) {
      if (unlikely(bigstack_alloc_w(raw_variant_ctl, qual_present_ptr) ||
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
      if (unlikely(bigstack_alloc_w(raw_variant_ctl, filter_present_ptr) ||
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
    if (info_pr_exists) {
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
      if (info_pr_exists) {
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
      }
      // bugfix (3 Feb 2021): temporary filter_storage space was allocated
      // regardless of whether any non-PASS values were ultimately found
      read_iter = &(read_iter[kLoadPvarBlockSize * sizeof(intptr_t)]);
    }
    if (info_pr_exists) {
      memcpy(&(nonref_flags[full_block_ct * (kLoadPvarBlockSize / kBitsPerWord)]), read_iter, last_bitblock_size);
      ZeroTrailingBits(raw_variant_ct, nonref_flags);
      // read_iter = &(read_iter[kLoadPvarBlockSize / CHAR_BIT]);
    }
    const uintptr_t read_iter_stride_base = kLoadPvarBlockSize * (sizeof(int32_t) + 2 * sizeof(intptr_t)) + (kLoadPvarBlockSize / CHAR_BIT) + (load_qual_col > 1) * ((kLoadPvarBlockSize / CHAR_BIT) + kLoadPvarBlockSize * sizeof(float)) + (load_filter_col > 1) * (2 * (kLoadPvarBlockSize / CHAR_BIT) + kLoadPvarBlockSize * sizeof(intptr_t)) + info_pr_exists * (kLoadPvarBlockSize / CHAR_BIT);
    if (allele_idx_end > 2 * S_CAST(uintptr_t, raw_variant_ct)) {
      if (no_multiallelic_allowed) {
        snprintf(g_logbuf, kLogbufSize, "Error: %s contains multiallelic variant(s), despite having no #CHROM header line. Add that header line to make it obvious that this isn't a valid .bim.\n", pvarname);
        goto LoadPvar_ret_MALFORMED_INPUT_WW;
      }
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
    const uint32_t name_ct = cip->name_ct;
    if (name_ct) {
      logprintf("Note: %u nonstandard chromosome code%s present.\n", name_ct, (name_ct == 1)? "" : "s");
    }
    const uint32_t last_chr_code = cip->max_code + name_ct;
    const uint32_t chr_word_ct = BitCtToWordCt(last_chr_code + 1);
    BitvecAnd(loaded_chr_mask, chr_word_ct, chr_mask);
    BigstackEndSet(tmp_alloc_end);
    *variant_ct_ptr = raw_variant_ct - exclude_ct;
    *neg_bp_seen_ptr = neg_bp_seen;
    *vpos_sortstatus_ptr = vpos_sortstatus;
    *allele_storage_ptr = allele_storage;
    // if only INFO/PR flag present, no need to reload
    if (!(info_nonpr_present || info_pr_nonflag_exists)) {
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
  LoadPvar_ret_MALFORMED_HEADER_LINE:
    logerrprintf("Error: Header line %" PRIuPTR " of %s is malformed.\n", line_idx);
    reterr = kPglRetMalformedInput;
    break;
  LoadPvar_ret_MULTIALLELIC_MISSING_ALLELE_CODE:
    // Note that we must permit one or even both alleles to be missing in the
    // biallelic case, for .ped data to be usable.
    snprintf(g_logbuf, kLogbufSize, "Error: Missing allele code in multiallelic variant on line %" PRIuPTR " of %s.\n", line_idx, pvarname);
    WordWrapB(0);
    logerrputsb();
  LoadPvar_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  LoadPvar_ret_MISSING_TOKENS:
    logerrprintfww("Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, pvarname);
    reterr = kPglRetMalformedInput;
    break;
    // old NONMISSING_ALT_MATCHES_REF sanity-check replaced with conditional
    // CheckAlleleUniqueness() call.
  }
 LoadPvar_ret_1:
  CleanupTextStream2(pvarname, &pvar_txs, &reterr);
  if (reterr) {
    BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  }
  return reterr;
}

// Useful when we need to read from multiple .pgens (possibly with multiallelic
// variants) at once, and the usual LoadPvar() memory footprint is larger than
// we want.
// allele_idx_offsets is allocated at the bottom of bigstack.  It's assumed to
// be initialized to nullptr.
// max_allele_ct is assumed to be initialized to 2.
// raw_variant_ctp, max_allele_slenp, and/or max_observed_line_blenp can be
// nullptr if the caller is uninterested in them.  If they're all null, this
// can exit immediately on .bim files.
PglErr LoadAlleleIdxOffsetsFromPvar(const char* pvarname, const char* file_descrip, uint32_t max_thread_ct, uint32_t* raw_variant_ctp, uint32_t* max_allele_slenp, uint32_t* max_observed_line_blenp, uintptr_t** allele_idx_offsets_ptr, uint32_t* max_allele_ctp) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  uintptr_t line_idx = 0;
  PglErr reterr = kPglRetSuccess;
  TextStream txs;
  PreinitTextStream(&txs);
  {
    uint32_t max_line_blen;
    if (unlikely(StandardizeMaxLineBlen(bigstack_left() / 4, &max_line_blen))) {
      goto LoadAlleleIdxOffsetsFromPvar_ret_NOMEM;
    }
    reterr = InitTextStreamEx(pvarname, 1, kMaxLongLine, max_line_blen, ClipU32(max_thread_ct - 1, 1, 4), &txs);
    if (unlikely(reterr)) {
      goto LoadAlleleIdxOffsetsFromPvar_ret_TSTREAM_FAIL;
    }
    uint32_t max_observed_line_blen = 0;
    const char* line_start;
    // Skip header.
    do {
      line_start = TextGet(&txs);
      if (unlikely(!line_start)) {
        if (!TextStreamErrcode2(&txs, &reterr)) {
          snprintf(g_logbuf, kLogbufSize, "Error: %s is empty.\n", file_descrip);
          goto LoadAlleleIdxOffsetsFromPvar_ret_MALFORMED_INPUT_WW;
        }
        goto LoadAlleleIdxOffsetsFromPvar_ret_TSTREAM_FAIL;
      }
      ++line_idx;
      const char* line_end = TextLineEnd(&txs);
      const uint32_t cur_line_blen = line_end - line_start;
      if (cur_line_blen > max_observed_line_blen) {
        max_observed_line_blen = cur_line_blen;
      }
    } while ((*line_start == '#') && (!tokequal_k(line_start, "#CHROM")));
    // multiallelics prohibited in .bim files (and .bim lookalikes to reduce
    // confusion)
    const uint32_t no_multiallelic_allowed = (*line_start != '#');
    uint32_t col_skips[2];
    uint32_t col_types[2];
    if (!no_multiallelic_allowed) {
      // [0] = REF
      // [1] = ALT
      const char* token_end = &(line_start[6]);
      uint32_t found_header_bitset = 0;
      uint32_t relevant_col_uidx = 0;
      for (uint32_t col_idx = 1; ; ++col_idx) {
        const char* token_start = FirstNonTspace(token_end);
        if (IsEolnKns(*token_start)) {
          break;
        }
        token_end = CurTokenEnd(token_start);
        const uint32_t token_slen = token_end - token_start;
        uint32_t cur_col_type;
        if (strequal_k(token_start, "REF", token_slen)) {
          cur_col_type = 0;
        } else if (strequal_k(token_start, "ALT", token_slen)) {
          cur_col_type = 1;
        } else {
          continue;
        }
        const uint32_t cur_col_type_shifted = 1 << cur_col_type;
        if (unlikely(found_header_bitset & cur_col_type_shifted)) {
          char* write_iter = strcpya_k(g_logbuf, "Error: Duplicate column header '");
          write_iter = memcpya(write_iter, token_start, token_slen);
          write_iter = strcpya_k(write_iter, "' on line ");
          write_iter = wtoa(line_idx, write_iter);
          write_iter = strcpya_k(write_iter, " of ");
          write_iter = strcpya(write_iter, file_descrip);
          memcpy_k(write_iter, ".\n\0", 4);
          goto LoadAlleleIdxOffsetsFromPvar_ret_MALFORMED_INPUT_WW;
        }
        found_header_bitset |= cur_col_type_shifted;
        col_skips[relevant_col_uidx] = col_idx;
        col_types[relevant_col_uidx++] = cur_col_type;
      }
      if (unlikely(found_header_bitset != 3)) {
        snprintf(g_logbuf, kLogbufSize, "Error: Missing column header(s) on line %" PRIuPTR " of %s. (POS, ID, REF, and ALT are required.)\n", line_idx, file_descrip);
        goto LoadAlleleIdxOffsetsFromPvar_ret_MALFORMED_INPUT_WW;
      }
      col_skips[1] -= col_skips[0];
      line_start = TextGet(&txs);
      if (unlikely(!line_start)) {
        if (!TextStreamErrcode2(&txs, &reterr)) {
          snprintf(g_logbuf, kLogbufSize, "Error: %s is empty.\n", file_descrip);
          goto LoadAlleleIdxOffsetsFromPvar_ret_MALFORMED_INPUT_WW;
        }
        goto LoadAlleleIdxOffsetsFromPvar_ret_TSTREAM_FAIL;
      }
      ++line_idx;
      const char* line_end = TextLineEnd(&txs);
      const uint32_t cur_line_blen = line_end - line_start;
      if (cur_line_blen > max_observed_line_blen) {
        max_observed_line_blen = cur_line_blen;
      }
    } else {
      const char* fifth_column_start = NextTokenMult(line_start, 4);
      if (unlikely(!fifth_column_start)) {
        snprintf(g_logbuf, kLogbufSize, "Error: Fewer columns than expected in %s.\n", file_descrip);
        goto LoadAlleleIdxOffsetsFromPvar_ret_MALFORMED_INPUT_WW;
      }
      if ((!raw_variant_ctp) && (!max_allele_slenp) && (!max_observed_line_blenp)) {
        goto LoadAlleleIdxOffsetsFromPvar_ret_1;
      }
      col_skips[1] = 1;
      col_types[0] = 1;
      col_types[1] = 0;
      const char* sixth_column_start = NextToken(fifth_column_start);
      if (sixth_column_start) {
        // #CHROM ID CM POS ALT REF
        col_skips[0] = 4;
      } else {
        // #CHROM ID POS ALT REF
        col_skips[0] = 3;
      }
    }
    uintptr_t* allele_idx_offsets = R_CAST(uintptr_t*, g_bigstack_base);
    uintptr_t* allele_idx_offsets_limit = R_CAST(uintptr_t*, g_bigstack_end);
    if (allele_idx_offsets_limit == allele_idx_offsets) {
      goto LoadAlleleIdxOffsetsFromPvar_ret_NOMEM;
    }
#ifdef __LP64__
    if (allele_idx_offsets_limit > &(allele_idx_offsets[kPglMaxVariantCt + 1])) {
      allele_idx_offsets_limit = &(allele_idx_offsets[kPglMaxVariantCt + 1]);
    }
#endif
    uintptr_t* allele_idx_offsets_iter = allele_idx_offsets;
    uintptr_t allele_idx = 0;
    *allele_idx_offsets_iter++ = allele_idx;
    uint32_t max_allele_ct = 2;
    uint32_t max_allele_slen = 1;
    while (1) {
      if (allele_idx_offsets_iter == allele_idx_offsets_limit) {
#ifdef __LP64__
        if (allele_idx_offsets_limit == &(allele_idx_offsets[kPglMaxVariantCt + 1])) {
          logerrputs("Error: " PROG_NAME_STR " does not support more than 2^31 - 3 variants.  We recommend using\nother software for very deep studies of small numbers of genomes.\n");
          goto LoadAlleleIdxOffsetsFromPvar_ret_MALFORMED_INPUT;
        }
#endif
        goto LoadAlleleIdxOffsetsFromPvar_ret_NOMEM;
      }
      const char* token_ptrs[2];
      uint32_t token_slens[2];
      const char* last_token_end = TokenLexK(line_start, col_types, col_skips, 2, token_ptrs, token_slens);
      if (unlikely(!last_token_end)) {
        goto LoadAlleleIdxOffsetsFromPvar_ret_MISSING_TOKENS;
      }
      if (token_slens[0] > max_allele_slen) {
        max_allele_slen = token_slens[0];
      }
      const uint32_t alt_col_slen = token_slens[1];
      const char* alt_iter = token_ptrs[1];
      const char* alt_end = &(alt_iter[alt_col_slen]);
      uint32_t extra_alt_ct = 0;
      if (alt_col_slen > 1) {
        extra_alt_ct = CountByte(alt_iter, ',', alt_col_slen);
        for (uint32_t alt_idx = 0; alt_idx != extra_alt_ct; ++alt_idx) {
          const char* cur_alt_end = AdvToDelim(alt_iter, ',');
          const uint32_t cur_allele_slen = cur_alt_end - alt_iter;
          if (cur_allele_slen > max_allele_slen) {
            max_allele_slen = cur_allele_slen;
          }
          alt_iter = &(cur_alt_end[1]);
        }
        const uint32_t cur_allele_slen = alt_end - alt_iter;
        if (cur_allele_slen > max_allele_slen) {
          max_allele_slen = cur_allele_slen;
        }
      }
      const uint32_t cur_allele_ct = 2 + extra_alt_ct;
      if (cur_allele_ct > max_allele_ct) {
        max_allele_ct = cur_allele_ct;
      }
      allele_idx += cur_allele_ct;
      *allele_idx_offsets_iter++ = allele_idx;
      const char* line_end = AdvPastDelim(alt_end, '\n');
      const uint32_t cur_line_blen = line_end - line_start;
      if (cur_line_blen > max_observed_line_blen) {
        max_observed_line_blen = cur_line_blen;
      }
      line_start = line_end;
      if (!TextGetUnsafe2K(&txs, &line_start)) {
        break;
      }
      ++line_idx;
      if (unlikely(*line_start == '#')) {
        snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of %s starts with a '#'. (This is only permitted before the first nonheader line, and if a #CHROM header line is present it must denote the end of the header block.)\n", line_idx, file_descrip);
        goto LoadAlleleIdxOffsetsFromPvar_ret_MALFORMED_INPUT_WW;
      }
    }
    if (unlikely(TextStreamErrcode2(&txs, &reterr))) {
      goto LoadAlleleIdxOffsetsFromPvar_ret_TSTREAM_FAIL;
    }
    const uintptr_t raw_variant_ct_p1 = allele_idx_offsets_iter - allele_idx_offsets;
    const uintptr_t raw_variant_ct = raw_variant_ct_p1 - 1;
    if (allele_idx > 2 * raw_variant_ct) {
      BigstackFinalizeW(allele_idx_offsets, raw_variant_ct_p1);
      *allele_idx_offsets_ptr = allele_idx_offsets;
    }
    *max_allele_ctp = max_allele_ct;
    if (max_allele_slenp) {
      *max_allele_slenp = max_allele_slen;
    }
    if (raw_variant_ctp) {
      *raw_variant_ctp = raw_variant_ct;
    }
    if (max_observed_line_blenp) {
      *max_observed_line_blenp = max_observed_line_blen;
    }
  }
  while (0) {
  LoadAlleleIdxOffsetsFromPvar_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  LoadAlleleIdxOffsetsFromPvar_ret_TSTREAM_FAIL:
    TextStreamErrPrint(file_descrip, &txs);
    break;
  LoadAlleleIdxOffsetsFromPvar_ret_MISSING_TOKENS:
    snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, file_descrip);
  LoadAlleleIdxOffsetsFromPvar_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
#ifdef __LP64__
  LoadAlleleIdxOffsetsFromPvar_ret_MALFORMED_INPUT:
#endif
    reterr = kPglRetMalformedInput;
    break;
  }
 LoadAlleleIdxOffsetsFromPvar_ret_1:
  CleanupTextStream2(file_descrip, &txs, &reterr);
  BigstackEndReset(bigstack_end_mark);
  if (reterr) {
    BigstackReset(bigstack_mark);
  }
  return reterr;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
