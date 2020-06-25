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

#include <errno.h>
#include <zlib.h>
#include "pvar_ffi_support.h"
#include "include/plink2_text.h"

namespace plink2 {

RefcountedWptr* CreateRefcountedWptr(uintptr_t size) {
  RefcountedWptr* rwp = S_CAST(RefcountedWptr*, malloc(sizeof(RefcountedWptr) + size * sizeof(intptr_t)));
  if (!rwp) {
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
    free(rwp);
  }
  *rwpp = nullptr;
}

void PreinitMinimalPvar(MinimalPvar* mpp) {
  mpp->variant_ids = nullptr;
  mpp->allele_idx_offsetsp = nullptr;
  mpp->variant_ct = 0;
  mpp->max_allele_ct = 0;
}

// note that \xxx character constants are interpreted in octal.
// technically no need to represent 0-31, but 64 extra bytes of data is
// probably cheaper than the code to subtract 32 everywhere.
const char g_one_char_strs[] = "\0\0\1\0\2\0\3\0\4\0\5\0\6\0\7\0\10\0\11\0\12\0\13\0\14\0\15\0\16\0\17\0\20\0\21\0\22\0\23\0\24\0\25\0\26\0\27\0\30\0\31\0\32\0\33\0\34\0\35\0\36\0\37\0\40\0\41\0\42\0\43\0\44\0\45\0\46\0\47\0\50\0\51\0\52\0\53\0\54\0\55\0\56\0\57\0\60\0\61\0\62\0\63\0\64\0\65\0\66\0\67\0\70\0\71\0\72\0\73\0\74\0\75\0\76\0\77\0\100\0\101\0\102\0\103\0\104\0\105\0\106\0\107\0\110\0\111\0\112\0\113\0\114\0\115\0\116\0\117\0\120\0\121\0\122\0\123\0\124\0\125\0\126\0\127\0\130\0\131\0\132\0\133\0\134\0\135\0\136\0\137\0\140\0\141\0\142\0\143\0\144\0\145\0\146\0\147\0\150\0\151\0\152\0\153\0\154\0\155\0\156\0\157\0\160\0\161\0\162\0\163\0\164\0\165\0\166\0\167\0\170\0\171\0\172\0\173\0\174\0\175\0\176\0\177\0\200\0\201\0\202\0\203\0\204\0\205\0\206\0\207\0\210\0\211\0\212\0\213\0\214\0\215\0\216\0\217\0\220\0\221\0\222\0\223\0\224\0\225\0\226\0\227\0\230\0\231\0\232\0\233\0\234\0\235\0\236\0\237\0\240\0\241\0\242\0\243\0\244\0\245\0\246\0\247\0\250\0\251\0\252\0\253\0\254\0\255\0\256\0\257\0\260\0\261\0\262\0\263\0\264\0\265\0\266\0\267\0\270\0\271\0\272\0\273\0\274\0\275\0\276\0\277\0\300\0\301\0\302\0\303\0\304\0\305\0\306\0\307\0\310\0\311\0\312\0\313\0\314\0\315\0\316\0\317\0\320\0\321\0\322\0\323\0\324\0\325\0\326\0\327\0\330\0\331\0\332\0\333\0\334\0\335\0\336\0\337\0\340\0\341\0\342\0\343\0\344\0\345\0\346\0\347\0\350\0\351\0\352\0\353\0\354\0\355\0\356\0\357\0\360\0\361\0\362\0\363\0\364\0\365\0\366\0\367\0\370\0\371\0\372\0\373\0\374\0\375\0\376\0\377";

PglErr LoadMinimalPvar(const char* fname, MinimalPvar* mpp, char* errstr_buf) {
  // Simple, somewhat inefficient two-pass loader.  Ignores everything but the
  // ID, REF, and ALT columns for now (since that's what we need to make
  // multiallelic .pgen files usable).
  TextStream txs;
  PreinitTextStream(&txs);
  uintptr_t line_idx = 0;
  PglErr reterr = kPglRetSuccess;
  {
    if (TextStreamOpen(fname, &txs)) {
      goto LoadMinimalPvar_ret_FILE_FAIL;
    }
    // [0] = ID
    // [1] = REF
    // [2] = ALT
    uint32_t col_skips[3];
    uint32_t col_types[3];
    const char* line_iter;
    while (1) {
      ++line_idx;
      line_iter = TextGet(&txs);
      if (*line_iter != '#') {
        // .bim:
        //   5 columns: #CHROM ID POS ALT REF
        //   6+ columns: #CHROM ID CM POS ALT REF
        const char* fifth_token = NextTokenMult(line_iter, 4);
        if (!fifth_token) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: No variants in %s.\n", fname);
          goto LoadMinimalPvar_ret_MISSING_TOKENS;
        }
        col_skips[0] = 2;
        col_skips[2] = 1;
        col_types[0] = 0;
        col_types[1] = 2;
        col_types[2] = 1;
        if (!NextToken(fifth_token)) {
          col_skips[1] = 2;
        } else {
          col_skips[1] = 3;
        }
        break;
      }
      if (tokequal_k(line_iter, "#CHROM")) {
        const char* token_end = &(line_iter[6]);
        uint32_t found_header_bitset = 0;
        uint32_t relevant_postchr_col_ct = 0;
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
              if (memequal_k(line_iter, "REF", 3)) {
                cur_col_type = 1;
              } else if (memequal_k(line_iter, "ALT", 3)) {
                cur_col_type = 2;
              } else {
                continue;
              }
            } else if (memequal_k(line_iter, "ID", 2)) {
              cur_col_type = 0;
            } else {
              continue;
            }
          } else if (strequal_k(line_iter, "FORMAT", token_slen)) {
            break;
          } else {
            continue;
          }
          const uint32_t cur_col_type_shifted = 1 << cur_col_type;
          if (found_header_bitset & cur_col_type_shifted) {
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
        if (relevant_postchr_col_ct != 3) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Missing column header(s) on line %" PRIuPTR " of %s. (ID, REF, and ALT are required.)\n", line_idx, fname);
          goto LoadMinimalPvar_ret_MALFORMED_INPUT;
        }
        col_skips[2] -= col_skips[1];
        col_skips[1] -= col_skips[0];
        // skip this line in main loop
        line_iter = &(TextLineEnd(&txs)[-1]);
        break;
      }
    }
    if (TextStreamErrcode(&txs)) {
      goto LoadMinimalPvar_ret_FILE_FAIL;
    }
    uintptr_t variant_ct = 0;
    uintptr_t allele_ct = 0;
    uintptr_t string_byte_ct = 0;
    uint32_t max_extra_alt_ct = 0;
    while (1) {
      if (!IsEolnKns(*line_iter)) {
        const char* token_ptrs[3];
        uint32_t token_slens[3];
        line_iter = TokenLexK(line_iter, col_types, col_skips, 3, token_ptrs, token_slens);
        if (!line_iter) {
          goto LoadMinimalPvar_ret_MISSING_TOKENS;
        }
        ++variant_ct;
        string_byte_ct += token_slens[0]; // +1 added later
        if (token_slens[1] > 1) {
          string_byte_ct += token_slens[1] + 1;
        }
        const uint32_t alt_slen = token_slens[2];
        if (alt_slen > 1) { // could also error out on isolated ','
          const char* alt_start = token_ptrs[2];
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
                if (!cur_allele_slen) {
                  goto LoadMinimalPvar_ret_EMPTY_ALLELE_CODE;
                }
                string_byte_ct += cur_allele_slen + 1;
              }
              alt_iter = &(cur_alt_end[1]);
            }
            const uint32_t cur_allele_slen = alt_end - alt_iter;
            if (cur_allele_slen != 1) {
              if (!cur_allele_slen) {
                goto LoadMinimalPvar_ret_EMPTY_ALLELE_CODE;
              }
              string_byte_ct += cur_allele_slen + 1;
            }
          }
        }
      }
      if (TextNextLineLstripK(&txs, &line_iter)) {
        if (TextStreamErrcode(&txs)) {
          goto LoadMinimalPvar_ret_FILE_FAIL;
        }
        break;
      }
      ++line_idx;
    }
    if (!variant_ct) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: No variants in %s.\n", fname);
      goto LoadMinimalPvar_ret_MALFORMED_INPUT;
    }
    if (variant_ct > 0x7ffffffd) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Too many variants in %s (%" PRIuPTR "; max 2^31 - 3).\n", fname, variant_ct);
      goto LoadMinimalPvar_ret_MALFORMED_INPUT;
    }
    allele_ct += 2 * variant_ct;
    string_byte_ct += variant_ct;
    const char** variant_ids = S_CAST(const char**, malloc((variant_ct + allele_ct) * sizeof(intptr_t)));
    if (!variant_ids) {
      goto LoadMinimalPvar_ret_NOMEM;
    }
    mpp->variant_ids = variant_ids;
    char* string_space_iter = S_CAST(char*, malloc(string_byte_ct));
    if (!string_space_iter) {
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
      if (max_extra_alt_ct >= kPglMaxAltAlleleCt) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Variant in %s has too many ALT alleles (%u; max 254).\n", fname, max_extra_alt_ct + 1);
        reterr = kPglRetNotYetSupported;
        goto LoadMinimalPvar_ret_1;
      }
      mpp->allele_idx_offsetsp = CreateRefcountedWptr(variant_ct + 1);
      if (!mpp->allele_idx_offsetsp) {
        goto LoadMinimalPvar_ret_NOMEM;
      }
      allele_idx_offsets = mpp->allele_idx_offsetsp->p;
    }
    if (TextRewind(&txs)) {
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
      const char* token_ptrs[3];
      uint32_t token_slens[3];
      line_iter = TokenLexK(line_iter, col_types, col_skips, 3, token_ptrs, token_slens);
      if (!line_iter) {
        // previously validated
        goto LoadMinimalPvar_ret_READ_FAIL;
      }
      if (allele_idx_offsets) {
        allele_idx_offsets[variant_idx] = allele_idx;
      }

      const uint32_t variant_id_slen = token_slens[0];
      memcpyx(string_space_iter, token_ptrs[0], variant_id_slen, '\0');
      variant_ids[variant_idx] = string_space_iter;
      string_space_iter = &(string_space_iter[variant_id_slen + 1]);

      const uint32_t ref_slen = token_slens[1];
      if (ref_slen == 1) {
        // don't bother with missing-code standardization for now
        allele_storage[allele_idx] = &(g_one_char_strs[2 * token_ptrs[1][0]]);
      } else {
        memcpyx(string_space_iter, token_ptrs[1], ref_slen, '\0');
        allele_storage[allele_idx] = string_space_iter;
        string_space_iter = &(string_space_iter[ref_slen + 1]);
      }
      ++allele_idx;

      const char* alt_start = token_ptrs[2];
      const uint32_t alt_slen = token_slens[2];
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
    if (TextStreamErrcode(&txs)) {
      goto LoadMinimalPvar_ret_FILE_FAIL;
    }
    if (variant_idx != variant_ct) {
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
  LoadMinimalPvar_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  LoadMinimalPvar_ret_MISSING_TOKENS:
    snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Fewer tokens than expected on line %" PRIuPTR " of %s.\n", line_idx, fname);
    reterr = kPglRetMalformedInput;
    break;
  LoadMinimalPvar_ret_EMPTY_ALLELE_CODE:
    snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Empty allele code on line %" PRIuPTR " of %s.\n", line_idx, fname);
    reterr = kPglRetMalformedInput;
    break;
  }
 LoadMinimalPvar_ret_1:
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
  mpp->variant_ct = 0;
  mpp->max_allele_ct = 0;
}


}  // namespace plink2
