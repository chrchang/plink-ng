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

#include "pvar_ffi_support.h"

#include <assert.h>
#include <errno.h>
#include <string.h>

#include "pgenlib_misc.h"
#include "plink2_bits.h"
#include "plink2_htable.h"
#include "plink2_string.h"
#include "plink2_text.h"

namespace plink2 {

RefcountedWptr* CreateRefcountedWptr(uintptr_t size) {
  // "- sizeof(uintptr_t) + size * sizeof(uintptr_t)" because p[1] is intended
  // as a flexible array member
  RefcountedWptr* rwp = S_CAST(RefcountedWptr*, malloc(sizeof(RefcountedWptr)));
  if (!rwp) {
    return nullptr;
  }
  rwp->p = S_CAST(uintptr_t*, malloc(size * sizeof(intptr_t)));
  if (!rwp->p) {
    free(rwp);
    return nullptr;
  }
  rwp->ref_ct = 1;
  return rwp;
}

void CondReleaseRefcountedWptr(RefcountedWptr** rwpp) {
  RefcountedWptr* rwp = *rwpp;
  if (!rwp) {
    return;
  }
  --rwp->ref_ct;
  if (!rwp->ref_ct) {
    free(rwp->p);
    free(rwp);
  }
  *rwpp = nullptr;
}

void PreinitMinimalPvar(MinimalPvar* mpp) {
  mpp->chr_names = nullptr;
  mpp->chr_idxs = nullptr;
  mpp->variant_bps = nullptr;
  mpp->variant_ids = nullptr;
  mpp->allele_idx_offsetsp = nullptr;
  mpp->chr_ct = 0;
  mpp->variant_ct = 0;
  mpp->max_allele_ct = 0;
}

// note that \xxx character constants are interpreted in octal.
// technically no need to represent 0-31, but 64 extra bytes of data is
// probably cheaper than the code to subtract 32 everywhere.
const char g_one_char_strs[] = "\0\0\1\0\2\0\3\0\4\0\5\0\6\0\7\0\10\0\11\0\12\0\13\0\14\0\15\0\16\0\17\0\20\0\21\0\22\0\23\0\24\0\25\0\26\0\27\0\30\0\31\0\32\0\33\0\34\0\35\0\36\0\37\0\40\0\41\0\42\0\43\0\44\0\45\0\46\0\47\0\50\0\51\0\52\0\53\0\54\0\55\0\56\0\57\0\60\0\61\0\62\0\63\0\64\0\65\0\66\0\67\0\70\0\71\0\72\0\73\0\74\0\75\0\76\0\77\0\100\0\101\0\102\0\103\0\104\0\105\0\106\0\107\0\110\0\111\0\112\0\113\0\114\0\115\0\116\0\117\0\120\0\121\0\122\0\123\0\124\0\125\0\126\0\127\0\130\0\131\0\132\0\133\0\134\0\135\0\136\0\137\0\140\0\141\0\142\0\143\0\144\0\145\0\146\0\147\0\150\0\151\0\152\0\153\0\154\0\155\0\156\0\157\0\160\0\161\0\162\0\163\0\164\0\165\0\166\0\167\0\170\0\171\0\172\0\173\0\174\0\175\0\176\0\177\0\200\0\201\0\202\0\203\0\204\0\205\0\206\0\207\0\210\0\211\0\212\0\213\0\214\0\215\0\216\0\217\0\220\0\221\0\222\0\223\0\224\0\225\0\226\0\227\0\230\0\231\0\232\0\233\0\234\0\235\0\236\0\237\0\240\0\241\0\242\0\243\0\244\0\245\0\246\0\247\0\250\0\251\0\252\0\253\0\254\0\255\0\256\0\257\0\260\0\261\0\262\0\263\0\264\0\265\0\266\0\267\0\270\0\271\0\272\0\273\0\274\0\275\0\276\0\277\0\300\0\301\0\302\0\303\0\304\0\305\0\306\0\307\0\310\0\311\0\312\0\313\0\314\0\315\0\316\0\317\0\320\0\321\0\322\0\323\0\324\0\325\0\326\0\327\0\330\0\331\0\332\0\333\0\334\0\335\0\336\0\337\0\340\0\341\0\342\0\343\0\344\0\345\0\346\0\347\0\350\0\351\0\352\0\353\0\354\0\355\0\356\0\357\0\360\0\361\0\362\0\363\0\364\0\365\0\366\0\367\0\370\0\371\0\372\0\373\0\374\0\375\0\376\0\377";

ENUM_U31_DEF_START()
  kPvarColPos = 0,
  kPvarColId,
  kPvarColRef,
  kPvarColAlt,
  kPvarColNull
ENUM_U31_DEF_END(PvarColidx);

FLAGSET_DEF_START()
  kfPvarColset0,
  kfPvarColsetPos = (1 << kPvarColPos),
  kfPvarColsetId = (1 << kPvarColId),
  kfPvarColsetRef = (1 << kPvarColRef),
  kfPvarColsetAlt = (1 << kPvarColAlt)
FLAGSET_DEF_END(PvarColFlags);

PglErr LoadMinimalPvarEx(const char* fname, LoadMinimalPvarFlags flags, MinimalPvar* mpp, char* errstr_buf) {
  // Simple, somewhat inefficient two-pass loader.  Only looks at
  // CHROM/POS/ID/REF/ALT; the first two are optional.
  const char** chr_names_tmp = nullptr;
  uint32_t* chr_names_htable = nullptr;
  uint32_t chr_ct = 0;
  TextStream txs;
  PreinitTextStream(&txs);
  uintptr_t line_idx = 0;
  PglErr reterr = kPglRetSuccess;
  {
    if (unlikely(TextStreamOpen(fname, &txs))) {
      goto LoadMinimalPvar_ret_FILE_FAIL;
    }
    const uint32_t chr_needed = ((flags & kfLoadMinimalPvarOmitChrom) == 0);
    const uint32_t pos_needed = ((flags & kfLoadMinimalPvarOmitPos) == 0);
    // [-1] = #CHROM (must be first column)
    // [0] = POS
    // [1] = ID
    // [2] = REF
    // [3] = ALT
    uint32_t col_skips[kPvarColNull];
    uint32_t col_types[kPvarColNull];
    uint32_t relevant_postchr_col_ct = 4;
    const char* line_iter;
    while (1) {
      ++line_idx;
      line_iter = TextGet(&txs);
      if (unlikely(!line_iter)) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: No variants in %s.\n", fname);
        goto LoadMinimalPvar_ret_MALFORMED_INPUT;
      }
      if (*line_iter != '#') {
        // .bim:
        //   5 columns: #CHROM ID POS ALT REF
        //   6+ columns: #CHROM ID CM POS ALT REF
        const char* fifth_token = NextTokenMult(line_iter, 4);
        if (unlikely(!fifth_token)) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: No variants in %s.\n", fname);
          goto LoadMinimalPvar_ret_MALFORMED_INPUT;
        }
        // bugfix (12 Jan 2025): this was off by 1
        col_skips[0] = 1;
        col_skips[2] = 1;
        col_skips[3] = 1;
        col_types[0] = kPvarColId;
        col_types[1] = kPvarColPos;
        col_types[2] = kPvarColAlt;
        col_types[3] = kPvarColRef;
        if (!NextToken(fifth_token)) {
          col_skips[1] = 1;
        } else {
          col_skips[1] = 2;
        }
        break;
      }
      if (tokequal_k(line_iter, "#CHROM")) {
        const char* token_end = &(line_iter[6]);
        uint32_t found_header_bitset = 0;
        relevant_postchr_col_ct = 0;
        for (uint32_t col_idx = 1; ; ++col_idx) {
          line_iter = FirstNonTspace(token_end);
          if (IsEolnKns(*line_iter)) {
            break;
          }
          token_end = CurTokenEnd(line_iter);
          const uint32_t token_slen = token_end - line_iter;
          uint32_t cur_col_type;
          if (token_slen <= 3) {
            if (token_slen == 3) {
              if (memequal_sk(line_iter, "POS")) {
                cur_col_type = kPvarColPos;
              } else if (memequal_sk(line_iter, "REF")) {
                cur_col_type = kPvarColRef;
              } else if (memequal_sk(line_iter, "ALT")) {
                cur_col_type = kPvarColAlt;
              } else {
                continue;
              }
            } else if (memequal_sk(line_iter, "ID")) {
              cur_col_type = kPvarColId;
            } else {
              continue;
            }
          } else if (strequal_k(line_iter, "FORMAT", token_slen)) {
            break;
          } else {
            continue;
          }
          const uint32_t cur_col_type_shifted = 1 << cur_col_type;
          if (unlikely(found_header_bitset & cur_col_type_shifted)) {
            // known token, so no overflow danger
            char* write_iter = strcpya_k(errstr_buf, "Error: Duplicate column header '");
            write_iter = memcpya(write_iter, line_iter, token_slen);
            write_iter = strcpya_k(write_iter, "' on line ");
            write_iter = wtoa(line_idx, write_iter);
            write_iter = strcpya_k(write_iter, " of ");
            write_iter = strcpya(write_iter, fname);
            memcpy_k(write_iter, ".\n\0", 4);
            goto LoadMinimalPvar_ret_MALFORMED_INPUT;
          }
          found_header_bitset |= cur_col_type_shifted;
          col_skips[relevant_postchr_col_ct] = col_idx;
          col_types[relevant_postchr_col_ct++] = cur_col_type;
        }
        uint32_t required_cols = kfPvarColsetId | kfPvarColsetRef | kfPvarColsetAlt;
        if (pos_needed) {
          required_cols |= kfPvarColsetPos;
        }
        if (unlikely((found_header_bitset & required_cols) != required_cols)) {
          // PRIuPTR does not seem to be working as intended on CRAN
          // win-builder, so we can't have nice things.
          // snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Missing column header(s) on line %" PRIuPTR " of %s. (ID, REF, ALT, and usually POS are required.)\n", line_idx, fname);
          char* write_iter = strcpya_k(errstr_buf, "Error: Missing column header(s) on line ");
          write_iter = wtoa(line_idx, write_iter);
          write_iter = strcpya_k(write_iter, " of ");
          write_iter = strcpya(write_iter, fname);
          strcpy_k(write_iter, ". (ID, REF, ALT, and usually POS are required.)\n");
          goto LoadMinimalPvar_ret_MALFORMED_INPUT;
        }
        for (uint32_t rpc_col_idx = relevant_postchr_col_ct - 1; rpc_col_idx; --rpc_col_idx) {
          col_skips[rpc_col_idx] -= col_skips[rpc_col_idx - 1];
        }
        // skip this line in main loop
        line_iter = &(TextLineEnd(&txs)[-1]);
        break;
      }
    }
    if (unlikely(TextStreamErrcode(&txs))) {
      goto LoadMinimalPvar_ret_FILE_FAIL;
    }

    if (chr_needed) {
      chr_names_tmp = S_CAST(const char**, malloc(kMaxChromosomes * sizeof(intptr_t)));
      chr_names_htable = S_CAST(uint32_t*, malloc(kChrHtableSize * sizeof(int32_t)));
      if (unlikely((!chr_names_tmp) || (!chr_names_htable))) {
        goto LoadMinimalPvar_ret_NOMEM;
      }
      SetAllU32Arr(kChrHtableSize, chr_names_htable);
    }
    uintptr_t variant_ct = 0;
    uintptr_t allele_ct = 0;
    uintptr_t string_byte_ct = 0;
    uint32_t max_extra_alt_ct = 0;
    while (1) {
      if (!IsEolnKns(*line_iter)) {
        if (chr_needed) {
          const char* chr_start = line_iter;
          const char* chr_end = CurTokenEnd(chr_start);
          const uint32_t chr_slen = chr_end - chr_start;
          if (unlikely(chr_slen > kMaxIdSlen)) {
            // snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Line %" PRIuPTR " of %s has an excessively long chromosome/contig name.\n", line_idx, fname);
            char* write_iter = strcpya_k(errstr_buf, "Error: Line ");
            write_iter = wtoa(line_idx, write_iter);
            write_iter = strcpya_k(write_iter, " of ");
            write_iter = strcpya(write_iter, fname);
            strcpy_k(write_iter, " has an excessively long chromosome/contig name.\n");
            goto LoadMinimalPvar_ret_MALFORMED_INPUT;
          }
          const uint32_t chr_idx = IdHtableFindNnt(chr_start, chr_names_tmp, chr_names_htable, chr_slen, kChrHtableSize);
          if (chr_idx == UINT32_MAX) {
            // update chromosome name list.  Wait till next pass to save
            // per-variant info since we haven't allocated chr_idxs yet.
            if (unlikely(chr_ct == kMaxChromosomes)) {
              snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Too many distinct chromosome/contig names in %s (limit %u).\n", fname, kMaxChromosomes);
              goto LoadMinimalPvar_ret_MALFORMED_INPUT;
            }
            // max ~65k, rather than tens or hundreds of millions.  So we don't
            // bother optimizing out most of these allocations.  (But maybe
            // worth compacting the allocations at the end of the first pass?)
            char* chr_name_copy = S_CAST(char*, malloc(chr_slen + 1));
            if (unlikely(!chr_name_copy)) {
              goto LoadMinimalPvar_ret_NOMEM;
            }
            memcpyx(chr_name_copy, chr_start, chr_slen, '\0');
            chr_names_tmp[chr_ct] = chr_name_copy;
            HtableAddNondup(chr_name_copy, chr_slen, kChrHtableSize, chr_ct, chr_names_htable);
            ++chr_ct;
          }
        }
        const char* token_ptrs[kPvarColNull];
        uint32_t token_slens[kPvarColNull];
        line_iter = TokenLexK(line_iter, col_types, col_skips, relevant_postchr_col_ct, token_ptrs, token_slens);
        if (unlikely(!line_iter)) {
          goto LoadMinimalPvar_ret_MISSING_TOKENS;
        }
        // POS always ignored during this pass
        ++variant_ct;
        // could enforce ID length limits
        string_byte_ct += token_slens[kPvarColId]; // +1 added later
        if (token_slens[1] > 1) {
          string_byte_ct += token_slens[kPvarColRef] + 1;
        }
        const uint32_t alt_slen = token_slens[kPvarColAlt];
        if (alt_slen > 1) { // could also error out on isolated ','
          const char* alt_start = token_ptrs[kPvarColAlt];
          const uint32_t extra_alt_ct = CountByte(alt_start, ',', alt_slen);
          if (!extra_alt_ct) {
            string_byte_ct += alt_slen + 1;
          } else {
            if (extra_alt_ct > max_extra_alt_ct) {
              max_extra_alt_ct = extra_alt_ct;
            }
            allele_ct += extra_alt_ct;
            const char* alt_iter = alt_start;
            const char* alt_end = &(alt_start[alt_slen]);
            for (uint32_t alt_idx = 0; alt_idx != extra_alt_ct; ++alt_idx) {
              const char* cur_alt_end = AdvToDelim(alt_iter, ',');
              const uint32_t cur_allele_slen = cur_alt_end - alt_iter;
              if (cur_allele_slen != 1) {
                if (unlikely(!cur_allele_slen)) {
                  goto LoadMinimalPvar_ret_EMPTY_ALLELE_CODE;
                }
                string_byte_ct += cur_allele_slen + 1;
              }
              alt_iter = &(cur_alt_end[1]);
            }
            const uint32_t cur_allele_slen = alt_end - alt_iter;
            if (cur_allele_slen != 1) {
              if (unlikely(!cur_allele_slen)) {
                goto LoadMinimalPvar_ret_EMPTY_ALLELE_CODE;
              }
              string_byte_ct += cur_allele_slen + 1;
            }
          }
        }
      }
      if (TextNextLineLstripK(&txs, &line_iter)) {
        if (unlikely(TextStreamErrcode(&txs))) {
          goto LoadMinimalPvar_ret_FILE_FAIL;
        }
        break;
      }
      ++line_idx;
    }
    if (unlikely(!variant_ct)) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: No variants in %s.\n", fname);
      goto LoadMinimalPvar_ret_MALFORMED_INPUT;
    }
    if (unlikely(variant_ct > kPglMaxVariantCt)) {
      // snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Too many variants in %s (%" PRIuPTR "; max 2^31 - 3).\n", fname, variant_ct);
      char* write_iter = strcpya_k(errstr_buf, "Error: Too many variants in ");
      write_iter = strcpya(write_iter, fname);
      write_iter = strcpya_k(write_iter, " (");
      write_iter = wtoa(variant_ct, write_iter);
      strcpy_k(write_iter, "; max 2^31 - 3).\n");
      goto LoadMinimalPvar_ret_MALFORMED_INPUT;
    }
    ChrIdx* chr_idxs = nullptr;
    if (chr_needed) {
      assert(!mpp->chr_names);
      mpp->chr_names = S_CAST(const char**, realloc(chr_names_tmp, chr_ct * sizeof(intptr_t)));
      if (unlikely(!mpp->chr_names)) {
        goto LoadMinimalPvar_ret_NOMEM;
      }
      mpp->chr_ct = chr_ct;
      chr_names_tmp = nullptr;
      chr_idxs = S_CAST(ChrIdx*, malloc(variant_ct * sizeof(ChrIdx)));
      if (unlikely(!chr_idxs)) {
        goto LoadMinimalPvar_ret_NOMEM;
      }
      assert(!mpp->chr_idxs);
      mpp->chr_idxs = chr_idxs;
    }
    int32_t* variant_bps;
    if (pos_needed) {
      assert(!mpp->variant_bps);
      variant_bps = S_CAST(int32_t*, malloc(variant_ct * sizeof(int32_t)));
      if (unlikely(!variant_bps)) {
        goto LoadMinimalPvar_ret_NOMEM;
      }
      mpp->variant_bps = variant_bps;
    }
    allele_ct += 2 * variant_ct;
    string_byte_ct += variant_ct;
    const char** variant_ids = S_CAST(const char**, malloc((variant_ct + allele_ct) * sizeof(intptr_t)));
    if (unlikely(!variant_ids)) {
      goto LoadMinimalPvar_ret_NOMEM;
    }
    assert(!mpp->variant_ids);
    mpp->variant_ids = variant_ids;
    char* string_space_iter = S_CAST(char*, malloc(string_byte_ct));
    if (unlikely(!string_space_iter)) {
      variant_ids[0] = nullptr;
      goto LoadMinimalPvar_ret_NOMEM;
    }
    // ensure this is freed if we trigger one of the next few error conditions
    variant_ids[0] = string_space_iter;

#ifndef NDEBUG
    char* string_space_end = &(string_space_iter[string_byte_ct]);
#endif
    uintptr_t* allele_idx_offsets = nullptr;
    if (max_extra_alt_ct) {
      if (unlikely(max_extra_alt_ct >= kPglMaxAltAlleleCt)) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Variant in %s has too many ALT alleles (%u; max " PGL_MAX_ALT_ALLELE_CT_STR ").\n", fname, max_extra_alt_ct + 1);
        reterr = kPglRetNotYetSupported;
        goto LoadMinimalPvar_ret_1;
      }
      mpp->allele_idx_offsetsp = CreateRefcountedWptr(variant_ct + 1);
      if (unlikely(!mpp->allele_idx_offsetsp)) {
        goto LoadMinimalPvar_ret_NOMEM;
      }
      allele_idx_offsets = mpp->allele_idx_offsetsp->p;
    }
    if (unlikely(TextRewind(&txs))) {
      goto LoadMinimalPvar_ret_FILE_FAIL;
    }
    const char** allele_storage = &(variant_ids[variant_ct]);
    mpp->allele_storage = allele_storage;
    mpp->variant_ct = variant_ct;
    mpp->max_allele_ct = max_extra_alt_ct + 2;
    uint32_t variant_idx = 0;
    uintptr_t allele_idx = 0;
    while (!TextNextLineLstripK(&txs, &line_iter)) {
      if (*line_iter == '#') {
        continue;
      }
      if (chr_needed) {
        const char* chr_start = line_iter;
        const char* chr_end = CurTokenEnd(chr_start);
        const uint32_t chr_slen = chr_end - chr_start;
        const uint32_t chr_idx = IdHtableFindNnt(chr_start, mpp->chr_names, chr_names_htable, chr_slen, kChrHtableSize);
        if (unlikely(chr_idx == UINT32_MAX)) {
          // should have been loaded during previous pass
          goto LoadMinimalPvar_ret_READ_FAIL;
        }
        chr_idxs[variant_idx] = chr_idx;
      }
      const char* token_ptrs[kPvarColNull];
      uint32_t token_slens[kPvarColNull];
      line_iter = TokenLexK(line_iter, col_types, col_skips, relevant_postchr_col_ct, token_ptrs, token_slens);
      if (unlikely(!line_iter)) {
        // previously validated
        goto LoadMinimalPvar_ret_READ_FAIL;
      }
      if (allele_idx_offsets) {
        allele_idx_offsets[variant_idx] = allele_idx;
      }

      if (pos_needed) {
        if (unlikely(ScanIntAbsDefcap(token_ptrs[kPvarColPos], &(variant_bps[variant_idx])))) {
          // snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid POS on line %" PRIuPTR " of %s.\n", line_idx, fname);
          char* write_iter = strcpya_k(errstr_buf, "Error: Invalid POS on line ");
          write_iter = wtoa(line_idx, write_iter);
          write_iter = strcpya_k(write_iter, " of ");
          write_iter = strcpya(write_iter, fname);
          memcpy_k(write_iter, ".\n\0", 4);
          goto LoadMinimalPvar_ret_MALFORMED_INPUT;
        }
      }

      const uint32_t variant_id_slen = token_slens[kPvarColId];
      memcpyx(string_space_iter, token_ptrs[kPvarColId], variant_id_slen, '\0');
      variant_ids[variant_idx] = string_space_iter;
      string_space_iter = &(string_space_iter[variant_id_slen + 1]);

      const uint32_t ref_slen = token_slens[kPvarColRef];
      if (ref_slen == 1) {
        // don't bother with missing-code standardization for now
        allele_storage[allele_idx] = &(g_one_char_strs[2 * token_ptrs[kPvarColRef][0]]);
      } else {
        memcpyx(string_space_iter, token_ptrs[kPvarColRef], ref_slen, '\0');
        allele_storage[allele_idx] = string_space_iter;
        string_space_iter = &(string_space_iter[ref_slen + 1]);
      }
      ++allele_idx;

      const char* alt_start = token_ptrs[kPvarColAlt];
      const uint32_t alt_slen = token_slens[kPvarColAlt];
      if (alt_slen == 1) {
        allele_storage[allele_idx++] = &(g_one_char_strs[2 * alt_start[0]]);
      } else {
        const uint32_t extra_alt_ct = CountByte(alt_start, ',', alt_slen);
        if (!extra_alt_ct) {
          memcpyx(string_space_iter, alt_start, alt_slen, '\0');
          allele_storage[allele_idx++] = string_space_iter;
          string_space_iter = &(string_space_iter[alt_slen + 1]);
        } else {
          const char* alt_iter = alt_start;
          const char* alt_end = &(alt_start[alt_slen]);
          for (uint32_t alt_idx = 0; alt_idx != extra_alt_ct; ++alt_idx) {
            const char* cur_alt_end = AdvToDelim(alt_iter, ',');
            const uint32_t cur_allele_slen = cur_alt_end - alt_iter;
            if (cur_allele_slen == 1) {
              allele_storage[allele_idx] = &(g_one_char_strs[2 * alt_iter[0]]);
            } else {
              memcpyx(string_space_iter, alt_iter, cur_allele_slen, '\0');
              allele_storage[allele_idx] = string_space_iter;
              string_space_iter = &(string_space_iter[cur_allele_slen + 1]);
            }
            ++allele_idx;
            alt_iter = &(cur_alt_end[1]);
          }
          const uint32_t cur_allele_slen = alt_end - alt_iter;
          if (cur_allele_slen == 1) {
            allele_storage[allele_idx] = &(g_one_char_strs[2 * alt_iter[0]]);
          } else {
            memcpyx(string_space_iter, alt_iter, cur_allele_slen, '\0');
            allele_storage[allele_idx] = string_space_iter;
            string_space_iter = &(string_space_iter[cur_allele_slen + 1]);
          }
          ++allele_idx;
        }
      }
      ++variant_idx;
    }
    if (unlikely(TextStreamErrcode(&txs))) {
      goto LoadMinimalPvar_ret_FILE_FAIL;
    }
    if (unlikely(variant_idx != variant_ct)) {
      goto LoadMinimalPvar_ret_READ_FAIL;
    }
#ifndef NDEBUG
    assert(string_space_iter == string_space_end);
#endif
    if (allele_idx_offsets) {
      allele_idx_offsets[variant_idx] = allele_idx;
    }
  }
  while (0) {
  LoadMinimalPvar_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  LoadMinimalPvar_ret_READ_FAIL:
    errno = 0;
    reterr = kPglRetReadFail;
    break;
  LoadMinimalPvar_ret_FILE_FAIL:
    reterr = TextStreamErrcode(&txs);
    assert(reterr != kPglRetSuccess);
    {
      const char* errmsg = TextStreamError(&txs);
      if (errmsg) {
        memcpy_k(errstr_buf, "Error: ", 7);
        strcpy(&(errstr_buf[7]), errmsg);
      }
    }
    break;
  LoadMinimalPvar_ret_MISSING_TOKENS:
    // snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Fewer tokens than expected on line %" PRIuPTR " of %s.\n", line_idx, fname);
    {
      char* write_iter = strcpya_k(errstr_buf, "Error: Fewer tokens than expected on line ");
      write_iter = wtoa(line_idx, write_iter);
      write_iter = strcpya_k(write_iter, " of ");
      write_iter = strcpya(write_iter, fname);
      memcpy_k(write_iter, ".\n\0", 4);
    }
    // fallthrough
  LoadMinimalPvar_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  LoadMinimalPvar_ret_EMPTY_ALLELE_CODE:
    // snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Empty allele code on line %" PRIuPTR " of %s.\n", line_idx, fname);
    {
      char* write_iter = strcpya_k(errstr_buf, "Error: Empty allele code on line ");
      write_iter = wtoa(line_idx, write_iter);
      write_iter = strcpya_k(write_iter, " of ");
      write_iter = strcpya(write_iter, fname);
      memcpy_k(write_iter, ".\n\0", 4);
    }
    reterr = kPglRetMalformedInput;
    break;
  }
 LoadMinimalPvar_ret_1:
  free_cond(chr_names_htable);
  if (chr_names_tmp) {
    for (uint32_t chr_idx = 0; chr_idx != chr_ct; ++chr_idx) {
      free_const(chr_names_tmp[chr_idx]);
    }
    free(chr_names_tmp);
  }
  CleanupTextStream(&txs, &reterr);
  return reterr;
}

void CleanupMinimalPvar(MinimalPvar* mpp) {
  if (mpp->variant_ids) {
    if (mpp->variant_ids[0]) {
      free_const(mpp->variant_ids[0]);
    }
    CondReleaseRefcountedWptr(&mpp->allele_idx_offsetsp);
    free(mpp->variant_ids);
    mpp->variant_ids = nullptr;
  }
  if (mpp->chr_names) {
    for (uint32_t chr_idx = 0; chr_idx != mpp->chr_ct; ++chr_idx) {
      free_const(mpp->chr_names[chr_idx]);
    }
    free(mpp->chr_names);
    mpp->chr_names = nullptr;
  }
  if (mpp->chr_idxs) {
    free_const(mpp->chr_idxs);
    mpp->chr_idxs = nullptr;
  }
  if (mpp->variant_bps) {
    free_const(mpp->variant_bps);
    mpp->variant_bps = nullptr;
  }
  mpp->chr_ct = 0;
  mpp->variant_ct = 0;
  mpp->max_allele_ct = 0;
}


}  // namespace plink2
