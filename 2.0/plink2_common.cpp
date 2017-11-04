// This library is part of PLINK 2.00, copyright (C) 2005-2017 Shaun Purcell,
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

static_assert(kDosageMax == 32768, "print_dosage() needs to be updated.");
char* print_dosage(uint64_t dosage, char* start) {
  // 3 digit precision seems like the best compromise between accuracy and
  // avoidance of rounding ugliness
  // (Rounding ugliness is not actually hidden for e.g. 1000 Genomes phase 1,
  // since there are lots of 0.05 and 0.1 dosages which all get rounded in the
  // same direction; oh well.)

  // +16 since we need to round .99951 up to 1
  const uint64_t dosage_p16 = dosage + 16;
  start = uint32toa(dosage_p16 / kDosageMax, start);
  const uint32_t remainder_p16 = ((uint32_t)dosage_p16) & (kDosageMax - 1);
  if (remainder_p16 < 33) {
    return start;
  }
  // (1000 * remainder + 16384) / 32768
  //   1/16 = .0625 -> print 0.062
  //   3/16 = .1875 -> print 0.188
  //   5/16 = .3125 -> print 0.312
  // const uint32_t three_decimal_places = ((125 * remainder + 2048) / 4096) - ((remainder % 8192) == 2048);
  const uint32_t three_decimal_places = ((125 * remainder_p16 + 48) / 4096) - ((remainder_p16 % 8192) == 4048);
  // three_decimal_places guaranteed to be nonzero here
  *start++ = '.';
  const uint32_t first_decimal_place = three_decimal_places / 100;
  *start++ = '0' + first_decimal_place;
  const uint32_t last_two_digits = three_decimal_places - first_decimal_place * 100;
  if (last_two_digits) {
    memcpy(start, &(kDigitPair[last_two_digits]), 2);
    return &(start[1 + (start[1] != '0')]);
  }
  return start;
}

void populate_dense_dosage(const uintptr_t* genovec, const uintptr_t* dosage_present, const dosage_t* dosage_vals, uint32_t sample_ct, uint32_t dosage_ct, dosage_t* dense_dosage) {
  // see also fill_cur_dosage_ints in plink2_matrix_calc
  const uint32_t sample_ctl2_m1 = (sample_ct - 1) / kBitsPerWordD2;
  dosage_t lookup_table[4];
  lookup_table[0] = 0;
  lookup_table[1] = kDosageMid;
  lookup_table[2] = kDosageMax;
  lookup_table[3] = kDosageMissing;
  uint32_t loop_len = kBitsPerWordD2;
  uint32_t widx = 0;
  dosage_t* dense_dosage_iter = dense_dosage;
  while (1) {
    if (widx >= sample_ctl2_m1) {
      if (widx > sample_ctl2_m1) {
        break;
      }
      loop_len = MOD_NZ(sample_ct, kBitsPerWordD2);
    }
    uintptr_t cur_geno_word = genovec[widx];
    for (uint32_t uii = 0; uii < loop_len; ++uii) {
      const uintptr_t cur_geno = cur_geno_word & 3;
      *dense_dosage_iter++ = lookup_table[cur_geno];
      cur_geno_word >>= 2;
    }
    ++widx;
  }
  // fill trailing bits with missing values so vector operations work
  const uint32_t trailing_entry_ct = (-sample_ct) % kDosagePerVec;
  for (uint32_t uii = 0; uii < trailing_entry_ct; ++uii) {
    dense_dosage_iter[uii] = kDosageMissing;
  }

  uint32_t sample_idx = 0;
  for (uint32_t dosage_idx = 0; dosage_idx < dosage_ct; ++dosage_idx, ++sample_idx) {
    next_set_unsafe_ck(dosage_present, &sample_idx);
    dense_dosage[sample_idx] = dosage_vals[dosage_idx];
  }
}

void set_het_missing(uintptr_t word_ct, uintptr_t* genovec) {
  // 01 -> 11, nothing else changes
#ifdef __LP64__
  const vul_t m1 = VCONST_UL(kMask5555);
  vul_t* geno_vvec_iter = (vul_t*)genovec;
  const uintptr_t full_vec_ct = word_ct / kWordsPerVec;
  if (full_vec_ct & 1) {
    const vul_t cur_geno_vword = *geno_vvec_iter;
    const vul_t cur_geno_vword_low_lshifted = vul_lshift(cur_geno_vword & m1, 1);
    *geno_vvec_iter++ = cur_geno_vword | cur_geno_vword_low_lshifted;
  }
  if (full_vec_ct & 2) {
    vul_t cur_geno_vword = *geno_vvec_iter;
    vul_t cur_geno_vword_low_lshifted = vul_lshift(cur_geno_vword & m1, 1);
    *geno_vvec_iter++ = cur_geno_vword | cur_geno_vword_low_lshifted;
    cur_geno_vword = *geno_vvec_iter;
    cur_geno_vword_low_lshifted = vul_lshift(cur_geno_vword & m1, 1);
    *geno_vvec_iter++ = cur_geno_vword | cur_geno_vword_low_lshifted;
  }
  for (uintptr_t ulii = 3; ulii < full_vec_ct; ulii += 4) {
    vul_t cur_geno_vword = *geno_vvec_iter;
    vul_t cur_geno_vword_low_lshifted = vul_lshift(cur_geno_vword & m1, 1);
    *geno_vvec_iter++ = cur_geno_vword | cur_geno_vword_low_lshifted;
    cur_geno_vword = *geno_vvec_iter;
    cur_geno_vword_low_lshifted = vul_lshift(cur_geno_vword & m1, 1);
    *geno_vvec_iter++ = cur_geno_vword | cur_geno_vword_low_lshifted;
    cur_geno_vword = *geno_vvec_iter;
    cur_geno_vword_low_lshifted = vul_lshift(cur_geno_vword & m1, 1);
    *geno_vvec_iter++ = cur_geno_vword | cur_geno_vword_low_lshifted;
    cur_geno_vword = *geno_vvec_iter;
    cur_geno_vword_low_lshifted = vul_lshift(cur_geno_vword & m1, 1);
    *geno_vvec_iter++ = cur_geno_vword | cur_geno_vword_low_lshifted;
  }
  #ifdef USE_AVX2
  if (word_ct & 2) {
    const uintptr_t base_idx = full_vec_ct * kWordsPerVec;
    uintptr_t geno_word = genovec[base_idx];
    genovec[base_idx] = geno_word | ((geno_word & kMask5555) << 1);
    geno_word = genovec[base_idx + 1];
    genovec[base_idx + 1] = geno_word | ((geno_word & kMask5555) << 1);
  }
  #endif
  if (word_ct & 1) {
    const uintptr_t geno_word = genovec[word_ct - 1];
    genovec[word_ct - 1] = geno_word | ((geno_word & kMask5555) << 1);
  }
#else
  for (uintptr_t widx = 0; widx < word_ct; ++widx) {
    const uintptr_t geno_word = genovec[widx];
    genovec[widx] = geno_word | ((geno_word & kMask5555) << 1);
  }
#endif
}

void genoarr_to_nonmissing(const uintptr_t* genoarr, uint32_t sample_ct, uintptr_t* nonmissing_bitarr) {
  const uint32_t sample_ctl2 = QUATERCT_TO_WORDCT(sample_ct);
  const uintptr_t* genoarr_iter = genoarr;
  halfword_t* nonmissing_bitarr_iter = (halfword_t*)nonmissing_bitarr;
  for (uint32_t widx = 0; widx < sample_ctl2; ++widx) {
    uintptr_t ww = ~(*genoarr_iter++);
    ww = (ww | (ww >> 1)) & kMask5555;
    *nonmissing_bitarr_iter++ = pack_word_to_halfword(ww);
  }
  // zero trailing bits up to word boundary, in a way that doesn't create
  // aliasing issues
  // (if zeroing is needed up to vector boundary, that's the caller's
  // responsibility)
  const uint32_t trail_ct = sample_ct % kBitsPerWordD2;
  if (trail_ct) {
    nonmissing_bitarr_iter[-1] &= (1U << trail_ct) - 1;
  }
  if (sample_ctl2 % 2) {
    *nonmissing_bitarr_iter = 0;
  }
}

uint32_t genoarr_count_missing_notsubset_unsafe(const uintptr_t* genoarr, const uintptr_t* exclude_mask, uint32_t sample_ct) {
  const uint32_t sample_ctl2 = QUATERCT_TO_WORDCT(sample_ct);
  const uintptr_t* genoarr_iter = genoarr;
  const halfword_t* exclude_alias_iter = (halfword_t*)exclude_mask;
  uint32_t missing_ct = 0;
  for (uint32_t widx = 0; widx < sample_ctl2; ++widx) {
    uintptr_t ww = *genoarr_iter++;
    ww = ww & (ww >> 1);
    const uint32_t include_mask = ~(*exclude_alias_iter++);
    missing_ct += popcount01_long(ww & unpack_halfword_to_word(include_mask));
  }
  return missing_ct;
}


// assumes 0b11 == missing
// retired in favor of genoarr_to_nonmissing followed by a for loop?
/*
void copy_when_nonmissing(const uintptr_t* loadbuf, const void* source, uintptr_t elem_size, uintptr_t unfiltered_sample_ct, uintptr_t missing_ct, void* dest) {
  // tried hardcoding elem_size == 4 and 8; that was only ~4% faster, so may
  // as well use a single general-purpose function.
  if (!missing_ct) {
    memcpy(dest, source, unfiltered_sample_ct * elem_size);
    return;
  }
  const uintptr_t* loadbuf_iter = loadbuf;
  const uintptr_t* loadbuf_end = &(loadbuf[QUATERCT_TO_WORDCT(unfiltered_sample_ct)]);
  const unsigned char* source_alias = (const unsigned char*)source;
  char* dest_iter = (char*)dest;
  uintptr_t copy_start_idx = 0;
  uintptr_t sample_idx_offset = 0;
  do {
    uintptr_t cur_word = *loadbuf_iter++;
    cur_word = cur_word & (cur_word >> 1) & kMask5555;
    while (cur_word) {
      const uintptr_t new_missing_idx = sample_idx_offset + (CTZLU(cur_word) / 2);
      if (new_missing_idx != copy_start_idx) {
        const uintptr_t diff = new_missing_idx - copy_start_idx;
        dest_iter = memcpya(dest_iter, &(source_alias[copy_start_idx * elem_size]), diff * elem_size);
      }
      copy_start_idx = new_missing_idx + 1;
      cur_word &= cur_word - 1;
    }
    sample_idx_offset += kBitsPerWordD2;
  } while (loadbuf_iter < loadbuf_end);
  const uintptr_t diff = unfiltered_sample_ct - copy_start_idx;
  if (diff) {
    memcpy(dest_iter, &(source_alias[copy_start_idx * elem_size]), diff * elem_size);
  }
}
*/


uint32_t sid_col_required(const uintptr_t* sample_include, const char* sids, uint32_t sample_ct, uint32_t max_sid_blen, uint32_t maybe_modifier) {
  // note that MAYBESID and SID can both be set
  if (maybe_modifier & 2) {
    return 1;
  }
  if (sids && (maybe_modifier & 1)) {
    uint32_t sample_uidx = 0;
    for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
      next_set_unsafe_ck(sample_include, &sample_uidx);
      if (memcmp(&(sids[sample_uidx * max_sid_blen]), "0", 2)) {
        return 1;
      }
    }
  }
  return 0;
}

// sample_augid_map_ptr == nullptr ok
pglerr_t augid_init_alloc(const uintptr_t* sample_include, const char* sample_ids, const char* sids, uint32_t sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen, uint32_t** sample_augid_map_ptr, char** sample_augids_ptr, uintptr_t* max_sample_augid_blen_ptr) {
  if (!sids) {
    max_sid_blen = 2;
  }
  const uintptr_t max_sample_augid_blen = max_sample_id_blen + max_sid_blen;
  *max_sample_augid_blen_ptr = max_sample_augid_blen;
  uint32_t* sample_augid_map = nullptr;
  if (sample_augid_map_ptr) {
    if (bigstack_alloc_ui(sample_ct, sample_augid_map_ptr)) {
      return kPglRetNomem;
    }
    sample_augid_map = *sample_augid_map_ptr;
  }
  if (bigstack_alloc_c(max_sample_augid_blen * sample_ct, sample_augids_ptr)) {
    return kPglRetNomem;
  }
  char* sample_augids_iter = *sample_augids_ptr;
  uint32_t sample_uidx = 0;
  for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
    next_set_unsafe_ck(sample_include, &sample_uidx);
    char* write_iter = strcpyax(sample_augids_iter, &(sample_ids[sample_uidx * max_sample_id_blen]), '\t');
    if (sids) {
      strcpy(write_iter, &(sids[sample_uidx * max_sid_blen]));
    } else {
      strcpy(write_iter, "0");
    }
    sample_augids_iter = &(sample_augids_iter[max_sample_augid_blen]);
    if (sample_augid_map) {
      sample_augid_map[sample_idx] = sample_uidx;
    }
  }
  return kPglRetSuccess;
}

pglerr_t sorted_xidbox_init_alloc(const uintptr_t* sample_include, const char* sample_ids, const char* sids, uint32_t sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen, xid_mode_t xid_mode, uint32_t use_nsort, char** sorted_xidbox_ptr, uint32_t** xid_map_ptr, uintptr_t* max_xid_blen_ptr) {
  if (!(xid_mode & kfXidModeFlagSid)) {
    // two fields
    *max_xid_blen_ptr = max_sample_id_blen;
    return copy_sort_strbox_subset(sample_include, sample_ids, sample_ct, max_sample_id_blen, 0, 0, use_nsort, sorted_xidbox_ptr, xid_map_ptr);
  }
  // three fields
  if (augid_init_alloc(sample_include, sample_ids, sids, sample_ct, max_sample_id_blen, max_sid_blen, xid_map_ptr, sorted_xidbox_ptr, max_xid_blen_ptr)) {
    return kPglRetNomem;
  }
  if (sort_strbox_indexed(sample_ct, *max_xid_blen_ptr, use_nsort, *sorted_xidbox_ptr, *xid_map_ptr)) {
    return kPglRetNomem;
  }
  char* dup_id = scan_for_duplicate_ids(*sorted_xidbox_ptr, sample_ct, *max_xid_blen_ptr);
  if (dup_id) {
    char* tptr = (char*)rawmemchr(dup_id, '\t');
    *tptr = ' ';
    tptr = (char*)rawmemchr(&(tptr[1]), '\t');
    *tptr = ' ';
    LOGERRPRINTFWW("Error: Duplicate ID '%s'.\n", dup_id);
    return kPglRetMalformedInput;
  }
  return kPglRetSuccess;
}

boolerr_t sorted_xidbox_read_find(const char* __restrict sorted_xidbox, const uint32_t* __restrict xid_map, uintptr_t max_xid_blen, uintptr_t end_idx, uint32_t comma_delim, xid_mode_t xid_mode, char** read_pp, uint32_t* sample_uidx_ptr, char* __restrict idbuf) {
  // idbuf = workspace
  // sorted_xidbox = packed, sorted list of ID strings to search over.
  //
  // input *read_pp must point to beginning of FID; this is a change from plink
  // 1.9.
  //
  // *read_pp is now set to point to the end of the last parsed token instead
  // of the beginning of the next; this is another change from plink 1.9.
  //
  // returns 1 on missing token *or* if the sample ID is not present.  cases
  // can be distinguished by checking whether *read_pp == nullptr.
  char* first_token_start = *read_pp;
  uintptr_t blen_sid = 0;
  char* token_iter;
  char* iid_ptr;
  char* sid_ptr = nullptr;
  uintptr_t slen_fid;
  uintptr_t slen_iid;
  if (comma_delim) {
    token_iter = first_token_start;
    unsigned char ucc = (unsigned char)(*token_iter);
    while (ucc != ',') {
      if (ucc < 32) {
        if (!(xid_mode & kfXidModeFlagOneTokenOk)) {
          *read_pp = nullptr;
          return 1;
        }
        slen_fid = (uintptr_t)(token_iter - first_token_start);
        goto sorted_xidbox_read_find_comma_single_token;
      }
      ucc = (unsigned char)(*(++token_iter));
    }
    slen_fid = (uintptr_t)(token_iter - first_token_start);
    if (xid_mode & kfXidModeFlagNeverFid) {
    sorted_xidbox_read_find_comma_single_token:
      iid_ptr = first_token_start;
      slen_iid = slen_fid;
    } else {
      do {
        ucc = (unsigned char)(*(++token_iter));
      } while ((ucc == ' ') || (ucc == '\t'));
      iid_ptr = token_iter;
      while ((ucc >= 32) && (ucc != ',')) {
        ucc = (unsigned char)(*(++token_iter));
      }
      slen_iid = (uintptr_t)(token_iter - iid_ptr);
    }
    // token_iter now points to comma/eoln at end of IID
    if (xid_mode & kfXidModeFlagSid) {
      if (*token_iter != ',') {
        return 1;
      }
      do {
        ucc = (unsigned char)(*(++token_iter));
      } while ((ucc == ' ') || (ucc == '\t'));
      sid_ptr = token_iter;
      while ((ucc >= 32) && (ucc != ',')) {
        ucc = (unsigned char)(*(++token_iter));
      }
      blen_sid = 1 + (uintptr_t)(token_iter - sid_ptr);
      if (token_iter == sid_ptr) {
        // special case: treat missing SID as '0'
        blen_sid = 2;
        // const_cast, since token_endnn doesn't return const pointer
        // function is too long for me to be comfortable just turning off
        // -Wcast-qual...
        sid_ptr = (char*)((uintptr_t)(&(g_one_char_strs[96])));
      }
    }
  } else {
    assert(!is_eoln_kns(*first_token_start));
    token_iter = token_endnn(first_token_start);
    slen_fid = (uintptr_t)(token_iter - first_token_start);
    if (xid_mode & kfXidModeFlagNeverFid) {
    sorted_xidbox_read_find_space_single_token:
      iid_ptr = first_token_start;
      slen_iid = slen_fid;
    } else {
      token_iter = skip_initial_spaces(token_iter);
      if (is_eoln_kns(*token_iter)) {
        if (!(xid_mode & kfXidModeFlagOneTokenOk)) {
          *read_pp = nullptr;
          return 1;
        }
        // need to backtrack
        token_iter = &(first_token_start[slen_fid]);
        goto sorted_xidbox_read_find_space_single_token;
      }
      iid_ptr = token_iter;
      token_iter = token_endnn(iid_ptr);
      slen_iid = (uintptr_t)(token_iter - iid_ptr);
    }
    // token_iter now points to space/eoln at end of IID
    if (xid_mode & kfXidModeFlagSid) {
      token_iter = skip_initial_spaces(token_iter);
      if (is_eoln_kns(*token_iter)) {
        *read_pp = nullptr;
        return 1;
      }
      sid_ptr = token_iter;
      token_iter = token_endnn(sid_ptr);
      blen_sid = 1 + (uintptr_t)(token_iter - sid_ptr);
    }
  }
  *read_pp = token_iter;
  uintptr_t slen_final = slen_fid + slen_iid + blen_sid + 1;
  if (slen_final >= max_xid_blen) {
    // avoid buffer overflow
    return 1;
  }
  char* idbuf_end = memcpya(memcpyax(idbuf, first_token_start, slen_fid, '\t'), iid_ptr, slen_iid);
  if (blen_sid) {
    *idbuf_end++ = '\t';
    memcpy(idbuf_end, sid_ptr, blen_sid - 1);
  }
  return sorted_idbox_find(idbuf, sorted_xidbox, xid_map, slen_final, max_xid_blen, end_idx, sample_uidx_ptr);
}


pglerr_t load_xid_header(const char* flag_name, sid_detect_mode_t sid_detect_mode, uintptr_t loadbuf_size, char* loadbuf, char** loadbuf_iter_ptr, uintptr_t* line_idx_ptr, char** loadbuf_first_token_ptr, gzFile* gz_infile_ptr, xid_mode_t* xid_mode_ptr) {
  // possible todo: support comma delimiter
  uintptr_t line_idx = *line_idx_ptr;
  uint32_t is_header_line;
  char* loadbuf_first_token;
  do {
    ++line_idx;
    if (!gzgets(*gz_infile_ptr, loadbuf, loadbuf_size)) {
      if (!gzeof(*gz_infile_ptr)) {
        return kPglRetReadFail;
      }
      return kPglRetEmptyFile;
    }
    if (!loadbuf[loadbuf_size - 1]) {
      return kPglRetLongLine;
    }
    loadbuf_first_token = skip_initial_spaces(loadbuf);
    is_header_line = (loadbuf_first_token[0] == '#');
  } while (is_header_line && strcmp_se(&(loadbuf_first_token[1]), "FID", 3) && strcmp_se(&(loadbuf_first_token[1]), "IID", 3));
  xid_mode_t xid_mode = kfXidMode0;
  char* loadbuf_iter;
  if (is_header_line) {
    // the following header leading columns are supported:
    // #FID IID (sid_detect_mode can't be FORCE)
    // #FID IID SID (SID ignored on sid_detect_mode NOT_LOADED)
    // #IID
    // #IID SID
    loadbuf_iter = &(loadbuf_first_token[4]);
    if (loadbuf_first_token[1] == 'I') {
      xid_mode = kfXidModeFlagNeverFid;
    } else {
      loadbuf_iter = skip_initial_spaces(loadbuf_iter);
      if (strcmp_se(loadbuf_iter, "IID", 3)) {
        LOGERRPRINTF("Error: No IID column on line %" PRIuPTR " of --%s file.\n", line_idx, flag_name);
        return kPglRetMalformedInput;
      }
      loadbuf_iter = &(loadbuf_iter[3]);
    }
    loadbuf_iter = skip_initial_spaces(loadbuf_iter);
    if (!strcmp_se(loadbuf_iter, "SID", 3)) {
      if ((uint32_t)sid_detect_mode >= kSidDetectModeLoaded) {
        xid_mode |= kfXidModeFlagSid;
      }
      loadbuf_iter = skip_initial_spaces(&(loadbuf_iter[3]));
    } else if (sid_detect_mode == kSidDetectModeForce) {
      LOGERRPRINTFWW("Error: No SID column on line %" PRIuPTR " of --%s file.\n", line_idx, flag_name);
      return kPglRetMalformedInput;
    }
  } else {
    xid_mode = (sid_detect_mode == kSidDetectModeForce)? kfXidModeFidiidSid : kfXidModeFidiidOrIid;
    loadbuf_iter = loadbuf_first_token;
  }
  if (loadbuf_iter_ptr) {
    *loadbuf_iter_ptr = loadbuf_iter;
  }
  *loadbuf_first_token_ptr = loadbuf_first_token;
  *line_idx_ptr = line_idx;
  *xid_mode_ptr = xid_mode;
  return kPglRetSuccess;
}

pglerr_t open_and_load_xid_header(const char* fname, const char* flag_name, sid_detect_mode_t sid_detect_mode, uintptr_t loadbuf_size, char* loadbuf, char** loadbuf_iter_ptr, uintptr_t* line_idx_ptr, char** loadbuf_first_token_ptr, gzFile* gz_infile_ptr, xid_mode_t* xid_mode_ptr) {
  pglerr_t reterr = gzopen_read_checked(fname, gz_infile_ptr);
  if (reterr) {
    return reterr;
  }
  loadbuf[loadbuf_size - 1] = ' ';
  return load_xid_header(flag_name, sid_detect_mode, loadbuf_size, loadbuf, loadbuf_iter_ptr, line_idx_ptr, loadbuf_first_token_ptr, gz_infile_ptr, xid_mode_ptr);
}


const char g_xymt_log_names[][5] = {"chrX", "chrY", "XY", "chrM", "PAR1", "PAR2"};

static_assert(!(kChrRawEnd % kBytesPerVec), "kChrRawEnd must be a multiple of kBytesPerVec.");
pglerr_t init_chr_info(chr_info_t* cip) {
  // "constructor".  initializes with maximum capacity.  doesn't use bigstack.
  // chr_mask, haploid_mask: bits
  // chr_file_order, chr_idx_to_foidx: int32s
  // chr_fo_vidx_start: int32s, with an extra trailing element
  // nonstd_names: intptr_ts
  // nonstd_id_htable: kChrHtableSize int32s

  // this assumes kChrRawEnd is divisible by kBytesPerVec
  const uintptr_t vecs_required = 2 * BITCT_TO_VECCT(kChrRawEnd) + 3 * (kChrRawEnd / kInt32PerVec) + 1 + (kChrRawEnd / kWordsPerVec) + INT32CT_TO_VECCT(kChrHtableSize);

  // needed for proper cleanup
  cip->name_ct = 0;
  cip->incl_excl_name_stack = nullptr;
  if (vecaligned_malloc(vecs_required * kBytesPerVec, &(cip->chr_mask))) {
    return kPglRetNomem;
  }
  uintptr_t* alloc_iter = &(cip->chr_mask[BITCT_TO_VECCT(kChrRawEnd) * kWordsPerVec]);
  cip->haploid_mask = alloc_iter;
  alloc_iter = &(alloc_iter[BITCT_TO_VECCT(kChrRawEnd) * kWordsPerVec]);
  cip->chr_file_order = (uint32_t*)alloc_iter;
  alloc_iter = &(alloc_iter[(kChrRawEnd / kInt32PerVec) * kWordsPerVec]);
  cip->chr_fo_vidx_start = (uint32_t*)alloc_iter;
  alloc_iter = &(alloc_iter[((kChrRawEnd / kInt32PerVec) + 1) * kWordsPerVec]);
  cip->chr_idx_to_foidx = (uint32_t*)alloc_iter;
  alloc_iter = &(alloc_iter[(kChrRawEnd / kInt32PerVec) * kWordsPerVec]);
  cip->nonstd_names = (char**)alloc_iter;
  alloc_iter = &(alloc_iter[kChrRawEnd]);
  cip->nonstd_id_htable = (uint32_t*)alloc_iter;
  // alloc_iter = &(alloc_iter[((kChrHtableSize + (kInt32PerVec - 1)) / kInt32PerVec) * kWordsPerVec]);
  // fill_uint_one(kChrHtableSize, cip->nonstd_id_htable);

  fill_ulong_zero(kChrMaskWords, cip->chr_mask);
  fill_ulong_zero(kChrExcludeWords, cip->chr_exclude);

  // this is a change from plink 1.x.  MT > M since the former matches Ensembl,
  // while the latter doesn't match any major resource.  no "chr" to reduce
  // file sizes and reduce the impact of this change.
  cip->output_encoding = kfChrOutputMT;

  cip->zero_extra_chrs = 0;
  cip->is_include_stack = 0;
  cip->chrset_source = kChrsetSourceDefault;
  cip->autosome_ct = 22;
  for (uint32_t xymt_idx = 0; xymt_idx < kChrOffsetCt; ++xymt_idx) {
    cip->xymt_codes[xymt_idx] = 23 + xymt_idx;
  }
  cip->haploid_mask[0] = 0x1800000;
  fill_ulong_zero(kChrMaskWords - 1, &(cip->haploid_mask[1]));
  return kPglRetSuccess;
}

// explicit plink 1.07 species (now initialized by command line parser):
// human: 22, X, Y, XY, MT, PAR1, PAR2 (PAR1/PAR2 added, XY deprecated in plink
//   2.0)
// cow: 29, X, Y, MT
// dog: 38, X, Y, XY, MT
// horse: 31, X, Y
// mouse: 19, X, Y
// rice: 12
// sheep: 26, X, Y

// must be safe to call this twice.
void finalize_chrset(misc_flags_t misc_flags, chr_info_t* cip) {
  uint32_t autosome_ct = cip->autosome_ct;
  uint32_t max_code = autosome_ct;
  for (uint32_t xymt_idx_p1 = kChrOffsetCt; xymt_idx_p1; --xymt_idx_p1) {
    if (cip->xymt_codes[xymt_idx_p1 - 1] >= 0) {
      max_code = autosome_ct + xymt_idx_p1;
      break;
    }
  }

  // could initialize haploid_mask bits (after the first) here, instead of
  // earlier...

  cip->max_numeric_code = MINV(max_code, autosome_ct + 4);
  cip->max_code = max_code;
  uintptr_t* chr_mask = cip->chr_mask;
  uintptr_t last_chr_mask_word = chr_mask[kChrMaskWords - 1];
  int32_t* xymt_codes = cip->xymt_codes;
  if (last_chr_mask_word) {
    // avoids repeating some work if this is called twice
    chr_mask[kChrMaskWords - 1] = 0;

    uint32_t xymt_include = last_chr_mask_word >> (kBitsPerWord - kChrOffsetCt);
    do {
      const uint32_t xymt_idx = __builtin_ctz(xymt_include);
      const int32_t cur_chr_code = xymt_codes[xymt_idx];
      if (cur_chr_code >= 0) {
        set_bit(cur_chr_code, chr_mask);
      }
      xymt_include &= xymt_include - 1;
    } while (xymt_include);
  } else if (are_all_words_zero(chr_mask, kChrExcludeWords) && (!cip->is_include_stack)) {
    // init_default_chr_mask()
    fill_all_bits(cip->autosome_ct + 1, chr_mask);
    for (uint32_t xymt_idx = 0; xymt_idx < kChrOffsetCt; ++xymt_idx) {
      const int32_t cur_chr_code = cip->xymt_codes[xymt_idx];
      if (cur_chr_code >= 0) {
        set_bit(cur_chr_code, chr_mask);
      }
    }
  } else if (misc_flags & (kfMiscAutosomePar | kfMiscAutosomeOnly)) {
    fill_bits_nz(1, cip->autosome_ct + 1, chr_mask);
    clear_bits_nz(cip->autosome_ct + 1, kChrExcludeWords * kBitsPerWord, chr_mask);
    if (misc_flags & kfMiscAutosomePar) {
      int32_t par_chr_code = cip->xymt_codes[kChrOffsetXY];
      if (par_chr_code >= 0) {
        set_bit(par_chr_code, chr_mask);
      }
      par_chr_code = cip->xymt_codes[kChrOffsetPAR1];
      if (par_chr_code >= 0) {
        set_bit(par_chr_code, chr_mask);
      }
      par_chr_code = cip->xymt_codes[kChrOffsetPAR2];
      if (par_chr_code >= 0) {
        set_bit(par_chr_code, chr_mask);
      }
    }
  }

  uintptr_t* chr_exclude = cip->chr_exclude;
  uintptr_t last_chr_exclude_word = chr_exclude[kChrExcludeWords - 1];
  uint32_t xymt_exclude = last_chr_exclude_word >> (kBitsPerWord - kChrOffsetCt);
  last_chr_exclude_word &= (k1LU << (kBitsPerWord - kChrOffsetCt)) - k1LU;
  for (uint32_t widx = 0; widx < kChrExcludeWords - 1; ++widx) {
    chr_mask[widx] &= ~chr_exclude[widx];
  }
  chr_mask[kChrExcludeWords - 1] &= ~last_chr_exclude_word;
  if (xymt_exclude) {
    do {
      const uint32_t xymt_idx = __builtin_ctz(xymt_exclude);
      const int32_t cur_chr_code = xymt_codes[xymt_idx];
      if (cur_chr_code >= 0) {
        clear_bit(cur_chr_code, chr_mask);
      }
      xymt_exclude &= xymt_exclude - 1;
    } while (xymt_exclude);
  }
  fill_uint_one(max_code + 1, cip->chr_idx_to_foidx);
}

void forget_extra_chr_names(uint32_t reinitialize, chr_info_t* cip) {
  const uint32_t name_ct = cip->name_ct;
  if (name_ct) {
    char** nonstd_names = cip->nonstd_names;
    const uint32_t chr_idx_end = cip->max_code + 1 + name_ct;
    for (uint32_t chr_idx = cip->max_code + 1; chr_idx < chr_idx_end; ++chr_idx) {
      free(nonstd_names[chr_idx]);
      nonstd_names[chr_idx] = nullptr;
    }
    if (reinitialize) {
      // fill_uint_one(kChrHtableSize, cip->nonstd_id_htable);
      cip->name_ct = 0;
    }
  }
}

// not currently called.  might want to do so in the future.
pglerr_t finalize_chr_info(chr_info_t* cip) {
  const uint32_t chr_ct = cip->chr_ct;
  const uint32_t name_ct = cip->name_ct;
  const uint32_t chr_code_end = cip->max_code + 1 + name_ct;
  const uint32_t chr_code_bitvec_ct = BITCT_TO_VECCT(chr_code_end);
  const uint32_t chr_ct_int32vec_ct = INT32CT_TO_VECCT(chr_ct);
  const uint32_t chr_ct_p1_int32vec_ct = 1 + (chr_ct / kInt32PerVec);
  const uint32_t chr_code_end_int32vec_ct = INT32CT_TO_VECCT(chr_code_end);
  const uint32_t chr_code_end_wordvec_ct = WORDCT_TO_VECCT(chr_code_end);
  uint32_t final_vecs_required = 2 * chr_code_bitvec_ct + chr_ct_int32vec_ct + chr_ct_p1_int32vec_ct + chr_code_end_int32vec_ct;
  if (name_ct) {
    final_vecs_required += chr_code_end_wordvec_ct + INT32CT_TO_VECCT(kChrHtableSize);
  }
  uintptr_t* new_alloc;
  if (vecaligned_malloc(final_vecs_required * kBytesPerVec, &new_alloc)) {
    return kPglRetNomem;
  }
  uintptr_t* old_alloc = cip->chr_mask;
  uintptr_t* new_alloc_iter = new_alloc;

  memcpy(new_alloc_iter, cip->chr_mask, chr_code_bitvec_ct * kBytesPerVec);
  cip->chr_mask = new_alloc_iter;
  new_alloc_iter = &(new_alloc_iter[chr_code_bitvec_ct * kWordsPerVec]);

  memcpy(new_alloc_iter, cip->haploid_mask, chr_code_bitvec_ct * kBytesPerVec);
  cip->haploid_mask = new_alloc_iter;
  new_alloc_iter = &(new_alloc_iter[chr_code_bitvec_ct * kWordsPerVec]);

  memcpy(new_alloc_iter, cip->chr_file_order, chr_ct_int32vec_ct * kBytesPerVec);
  cip->chr_file_order = (uint32_t*)new_alloc_iter;
  new_alloc_iter = &(new_alloc_iter[chr_ct_int32vec_ct * kWordsPerVec]);

  memcpy(new_alloc_iter, cip->chr_fo_vidx_start, chr_ct_p1_int32vec_ct * kBytesPerVec);
  cip->chr_fo_vidx_start = (uint32_t*)new_alloc_iter;
  new_alloc_iter = &(new_alloc_iter[chr_ct_p1_int32vec_ct * kWordsPerVec]);

  memcpy(new_alloc_iter, cip->chr_idx_to_foidx, chr_code_end_int32vec_ct * kBytesPerVec);
  cip->chr_idx_to_foidx = (uint32_t*)new_alloc_iter;

  if (!name_ct) {
    cip->nonstd_names = nullptr;
    cip->nonstd_id_htable = nullptr;
  } else {
    new_alloc_iter = &(new_alloc_iter[chr_code_end_int32vec_ct * kWordsPerVec]);

    memcpy(new_alloc_iter, cip->nonstd_names, chr_code_end_wordvec_ct * kBytesPerVec);
    cip->nonstd_names = (char**)new_alloc_iter;
    new_alloc_iter = &(new_alloc_iter[chr_code_end_wordvec_ct * kWordsPerVec]);

    memcpy(new_alloc_iter, cip->nonstd_id_htable, kChrHtableSize * sizeof(int32_t));
    cip->nonstd_id_htable = (uint32_t*)new_alloc_iter;
  }
  vecaligned_free(old_alloc);
  return kPglRetSuccess;
}

void cleanup_chr_info(chr_info_t* cip) {
  if (cip->chr_mask) {
    forget_extra_chr_names(0, cip);
    vecaligned_free(cip->chr_mask);
    cip->chr_mask = nullptr;
  }
  ll_str_t* llstr_ptr = cip->incl_excl_name_stack;
  while (llstr_ptr) {
    ll_str_t* next_ptr = llstr_ptr->next;
    free(llstr_ptr);
    llstr_ptr = next_ptr;
  }
  cip->incl_excl_name_stack = nullptr;
}

char* chr_name_std(const chr_info_t* cip, uint32_t chr_idx, char* buf) {
  if (chr_idx > cip->max_numeric_code) {
    // this is encoding-independent.  users who require all numbers should use
    // 25 == XY instead.
    // this code will probably need to be changed later if we add more standard
    // nonnumeric codes.
    memcpyl3(buf, "PAR");
    buf[3] = '0' + (((int32_t)chr_idx) - cip->max_numeric_code);
    return &(buf[4]);
  }
  const uint32_t output_encoding = cip->output_encoding;
  if (output_encoding & (kfChrOutputPrefix | kfChrOutput0M)) {
    if (output_encoding == kfChrOutput0M) {
      // force two chars
      if (chr_idx <= cip->autosome_ct) {
        buf = memcpya(buf, &(kDigitPair[chr_idx]), 2);
      } else if ((int32_t)chr_idx == cip->xymt_codes[kChrOffsetY]) {
        buf = strcpya(buf, "XY");
      } else {
        *buf++ = '0';
        if ((int32_t)chr_idx == cip->xymt_codes[kChrOffsetX]) {
          *buf++ = 'X';
        } else {
          // assumes only X/Y/XY/MT defined
          *buf++ = ((int32_t)chr_idx == cip->xymt_codes[kChrOffsetY])? 'Y' : 'M';
        }
      }
      return buf;
    }
    buf = memcpyl3a(buf, "chr");
  }
  if ((!(output_encoding & (kfChrOutputM | kfChrOutputMT))) || (chr_idx <= cip->autosome_ct)) {
    return uint32toa(chr_idx, buf);
  }
  if ((int32_t)chr_idx == cip->xymt_codes[kChrOffsetX]) {
    *buf++ = 'X';
  } else if ((int32_t)chr_idx == cip->xymt_codes[kChrOffsetY]) {
    *buf++ = 'Y';
  } else if ((int32_t)chr_idx == cip->xymt_codes[kChrOffsetXY]) {
    buf = strcpya(buf, "XY");
  } else {
    *buf++ = 'M';
    if (output_encoding & kfChrOutputMT) {
      *buf++ = 'T';
    }
  }
  return buf;
}

char* chr_name_write(const chr_info_t* cip, uint32_t chr_idx, char* buf) {
  // assumes chr_idx is valid
  if (!chr_idx) {
    *buf++ = '0';
    return buf;
  }
  if (chr_idx <= cip->max_code) {
    return chr_name_std(cip, chr_idx, buf);
  }
  if (cip->zero_extra_chrs) {
    *buf++ = '0';
    return buf;
  }
  return strcpya(buf, cip->nonstd_names[chr_idx]);
}

uint32_t get_max_chr_slen(const chr_info_t* cip) {
  // does not include trailing null
  // can be overestimate
  // if more functions start calling this, it should just be built into
  // load_bim() instead
  if (cip->zero_extra_chrs) {
    return 3 + kMaxChrTextnum;
  }
  const uint32_t chr_ct = cip->chr_ct;
  const uint32_t max_code = cip->max_code;
  uint32_t max_chr_slen = 3 + kMaxChrTextnum;
  for (uint32_t chr_fo_idx = 0; chr_fo_idx < chr_ct; ++chr_fo_idx) {
    const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
    if (!is_set(cip->chr_mask, chr_idx)) {
      continue;
    }
    if (chr_idx > max_code) {
      const uint32_t name_slen = strlen(cip->nonstd_names[chr_idx]);
      if (name_slen > max_chr_slen) {
        max_chr_slen = name_slen;
      }
    }
  }
  return max_chr_slen;
}

uint32_t haploid_chr_present(const chr_info_t* cip) {
  const uintptr_t* chr_mask = cip->chr_mask;
  const uintptr_t* haploid_mask = cip->haploid_mask;
  // since we don't load haploid vs. diploid info from ##contig header lines,
  // this is sufficient
  for (uint32_t widx = 0; widx < kChrExcludeWords; ++widx) {
    if (chr_mask[widx] & haploid_mask[widx]) {
      return 1;
    }
  }
  return 0;
}

static inline int32_t single_cap_letter_chrom(uint32_t cap_letter) {
  if (cap_letter == 'X') {
    return kChrRawX;
  }
  if (cap_letter == 'Y') {
    return kChrRawY;
  }
  if (cap_letter == 'M') {
    return kChrRawMT;
  }
  return -1;
}

static_assert(kMaxChrTextnumSlen == 2, "get_chr_code_raw() must be updated.");
int32_t get_chr_code_raw(const char* sptr) {
  // any character <= ' ' is considered a terminator
  // note that char arithmetic tends to be compiled to int32 operations, so we
  // mostly work with ints here
  uint32_t first_char_code = (unsigned char)sptr[0];
  uint32_t first_char_toi;
  if (first_char_code < 58) {
  get_chr_code_raw_digits:
    first_char_toi = first_char_code - '0';
    if (first_char_toi < 10) {
      const uint32_t second_char_code = (unsigned char)sptr[1];
      if (second_char_code <= ' ') {
        return first_char_toi;
      }
      if (((unsigned char)sptr[2]) <= ' ') {
        const uint32_t second_char_toi = second_char_code - '0';
        if (second_char_toi < 10) {
          return first_char_toi * 10 + second_char_toi;
        }
        if (!first_char_toi) {
          // accept '0X', '0Y', '0M' emitted by Oxford software
          return single_cap_letter_chrom(second_char_code & 0xdf);
        }
      }
    }
    return -1;
  }
  first_char_code &= 0xdf;
  uint32_t second_char_code = (unsigned char)sptr[1];
  if (first_char_code == 'P') {
    // chrPAR1 *not* supported; has to be PAR1 by itself.
    // can't do uint16_t compare of multiple characters, since we could be
    // dealing with a length-1 null-terminated string; that IS faster when it's
    // safe, though
    if (((second_char_code & 0xdf) == 'A') && ((((unsigned char)sptr[2]) & 0xdf) == 'R')) {
      const uint32_t par_idx_m1 = ((unsigned char)sptr[3]) - '1';
      if ((par_idx_m1 < 2) && (((unsigned char)sptr[4]) <= ' ')) {
        return kChrRawPAR1 + par_idx_m1;
      }
    }
    return -1;
  }
  if (first_char_code == 'C') {
    if (((second_char_code & 0xdf) != 'H') || ((((unsigned char)sptr[2]) & 0xdf) != 'R')) {
      return -1;
    }
    sptr = &(sptr[3]);
    first_char_code = (unsigned char)sptr[0];
    if (first_char_code < 58) {
      goto get_chr_code_raw_digits;
    }
    first_char_code &= 0xdf;
    second_char_code = (unsigned char)sptr[1];
  }
  if (second_char_code <= ' ') {
    return single_cap_letter_chrom(first_char_code);
  }
  if (((unsigned char)sptr[2]) <= ' ') {
    second_char_code &= 0xdf;
    if ((first_char_code == 'X') && (second_char_code == 'Y')) {
      return kChrRawXY;
    } else if ((first_char_code == 'M') && (second_char_code == 'T')) {
      return kChrRawMT;
    }
  }
  return -1;
}

int32_t get_chr_code(const char* chr_name, const chr_info_t* cip, uint32_t name_slen) {
  // requires chr_name to be null-terminated
  // in practice, name_slen will usually already be known, may as well avoid
  // redundant strlen() calls even though this uglifies the interface
  // does not perform exhaustive error-checking
  // -1 = --allow-extra-chr ok, -2 = total fail
  int32_t chr_code_raw = get_chr_code_raw(chr_name);
  if (((uint32_t)chr_code_raw) <= cip->max_code) {
    return chr_code_raw;
  }
  if (chr_code_raw != -1) {
    if (chr_code_raw >= ((int32_t)kMaxContigs)) {
      return cip->xymt_codes[chr_code_raw - kMaxContigs];
    }
    return -2;
  }
  if (!cip->name_ct) {
    return -1;
  }
  // 0xffffffffU gets casted to -1
  return (int32_t)id_htable_find(chr_name, cip->nonstd_names, cip->nonstd_id_htable, name_slen, kChrHtableSize);
}

int32_t get_chr_code_counted(const chr_info_t* cip, uint32_t name_slen, char* chr_name) {
  // when the chromosome name isn't null-terminated
  char* s_end = &(chr_name[name_slen]);
  const char tmpc = *s_end;
  *s_end = '\0';
  const int32_t chr_code = get_chr_code(chr_name, cip, name_slen);
  *s_end = tmpc;
  return chr_code;
}

void chr_error(const char* chr_name, const char* file_descrip, const chr_info_t* cip, uintptr_t line_idx, int32_t error_code) {
  // assumes chr_name is null-terminated
  const int32_t raw_code = get_chr_code_raw(chr_name);
  logprint("\n");
  if (line_idx) {
    LOGERRPRINTFWW("Error: Invalid chromosome code '%s' on line %" PRIuPTR " of %s.\n", chr_name, line_idx, file_descrip);
  } else {
    LOGERRPRINTFWW("Error: Invalid chromosome code '%s' in %s.\n", chr_name, file_descrip);
  }
  if ((raw_code > ((int32_t)cip->max_code)) && ((raw_code <= (int32_t)(kMaxChrTextnum + kChrOffsetCt)) || (raw_code >= ((int32_t)kMaxContigs)))) {
    if (cip->chrset_source == kChrsetSourceDefault) {
      logerrprint("(This is disallowed for humans.  Check if the problem is with your data, or if\nyou forgot to define a different chromosome set with e.g. --chr-set.).\n");
    } else if (cip->chrset_source == kChrsetSourceCmdline) {
      logerrprint("(This is disallowed by your command-line flags.)\n");
    } else {
      // kChrsetSourceFile
      logerrprint("(This is disallowed by the file's own ##chrSet header line.)\n");
    }
    // maybe want to print message(s) depending on whether chromosome set was
    // defined on the command line or by the input file?
  } else if (error_code == -1) {
    logerrprint("(Use --allow-extra-chr to force it to be accepted.)\n");
  }
}

pglerr_t try_to_add_chr_name(const char* chr_name, const char* file_descrip, uintptr_t line_idx, uint32_t name_slen, uint32_t allow_extra_chrs, int32_t* chr_idx_ptr, chr_info_t* cip) {
  // assumes chr_name is either nonstandard (i.e. not "2", "chr2", "chrX",
  // etc.), or a rejected xymt.
  // requires chr_name to be null-terminated
  // assumes chr_idx currently has the return value of get_chr_code()
  if ((!allow_extra_chrs) || ((*chr_idx_ptr) == -2)) {
    chr_error(chr_name, file_descrip, cip, line_idx, *chr_idx_ptr);
    return kPglRetMalformedInput;
  }

  // quasi-bugfix: remove redundant hash table check

  if (chr_name[0] == '#') {
    // redundant with some of the comment-skipping loaders, but this isn't
    // performance-critical
    logprint("\n");
    logerrprint("Error: Chromosome/contig names may not begin with '#'.\n");
    return kPglRetMalformedInput;
  }
  if (name_slen > kMaxIdSlen) {
    logprint("\n");
    if (line_idx) {
      LOGERRPRINTFWW("Error: Line %" PRIuPTR " of %s has an excessively long chromosome/contig name. (The " PROG_NAME_STR " limit is " MAX_ID_SLEN_STR " characters.)\n", line_idx, file_descrip);
    } else {
      LOGERRPRINTFWW("Error: Excessively long chromosome/contig name in %s. (The " PROG_NAME_STR " limit is " MAX_ID_SLEN_STR " characters.)\n", file_descrip);
    }
    return kPglRetMalformedInput;
  }
  const uint32_t max_code_p1 = cip->max_code + 1;
  const uint32_t name_ct = cip->name_ct;
  const uint32_t chr_code_end = max_code_p1 + name_ct;
  if (chr_code_end == kMaxContigs) {
    logprint("\n");
    logerrprint("Error: Too many distinct nonstandard chromosome/contig names.\n");
    return kPglRetMalformedInput;
  }
  if (!name_ct) {
    // lazy initialization
    fill_uint_one(kChrHtableSize, cip->nonstd_id_htable);
  }
  char** nonstd_names = cip->nonstd_names;
  if (pgl_malloc(name_slen + 1, &(nonstd_names[chr_code_end]))) {
    return kPglRetNomem;
  }
  ll_str_t* name_stack_ptr = cip->incl_excl_name_stack;
  uint32_t in_name_stack = 0;
  while (name_stack_ptr) {
    // there shouldn't be many of these, so sorting is unimportant
    if (!strcmp(chr_name, name_stack_ptr->ss)) {
      in_name_stack = 1;
      break;
    }
    name_stack_ptr = name_stack_ptr->next;
  }
  if ((in_name_stack && cip->is_include_stack) || ((!in_name_stack) && (!cip->is_include_stack))) {
    SET_BIT(chr_code_end, cip->chr_mask);
    if (cip->haploid_mask[0] & 1) {
      SET_BIT(chr_code_end, cip->haploid_mask);
    }
  }
  memcpy(nonstd_names[chr_code_end], chr_name, name_slen + 1);
  *chr_idx_ptr = (int32_t)chr_code_end;
  cip->name_ct = name_ct + 1;
  uint32_t* id_htable = cip->nonstd_id_htable;
  uint32_t hashval = hashceil(chr_name, name_slen, kChrHtableSize);
  while (1) {
    if (id_htable[hashval] == 0xffffffffU) {
      id_htable[hashval] = chr_code_end;
      return kPglRetSuccess;
    }
    if (++hashval == kChrHtableSize) {
      hashval = 0;
    }
  }
}


/*
uintptr_t count_11_vecs(const vul_t* geno_vvec, uintptr_t vec_ct) {
  // Counts number of aligned 11s in vptr[0..(vec_ct-1)].  Assumes vec_ct is a
  // multiple of 6 (0 ok).
  assert(!(vec_ct % 6));
  const vul_t m1 = VCONST_UL(kMask5555);
  const vul_t m2 = VCONST_UL(kMask3333);
  const vul_t m4 = VCONST_UL(kMask0F0F);
  const vul_t m8 = VCONST_UL(kMask00FF);
  const vul_t* geno_vvec_iter = geno_vvec;
  const vul_t* geno_vvec_end = &(geno_vvec[vec_ct]);
  uintptr_t tot = 0;

  while (1) {
    const vul_t* geno_vvec_stop = &(geno_vvec_iter[60]);

    univec_t acc;
    acc.vi = vul_setzero();

    if (geno_vvec_stop > geno_vvec_end) {
      if (geno_vvec_iter == geno_vvec_end) {
        return tot;
      }
      geno_vvec_stop = geno_vvec_end;
    }
    do {
      vul_t cur_geno_vword = *geno_vvec_iter++;
      vul_t count1 = cur_geno_vword & m1;
      count1 = count1 & vul_rshift(cur_geno_vword, 1);

      cur_geno_vword = *geno_vvec_iter++;
      vul_t cur_11 = cur_geno_vword & m1;
      cur_11 = cur_11 & vul_rshift(cur_geno_vword, 1);
      count1 = count1 + cur_11;

      cur_geno_vword = *geno_vvec_iter++;
      cur_11 = cur_geno_vword & m1;
      cur_11 = cur_11 & vul_rshift(cur_geno_vword, 1);
      count1 = count1 + cur_11;
      count1 = (count1 & m2) + (vul_rshift(count1, 2) & m2);

      cur_geno_vword = *geno_vvec_iter++;
      vul_t count2 = cur_geno_vword & m1;
      count2 = count2 & vul_rshift(cur_geno_vword, 1);

      cur_geno_vword = *geno_vvec_iter++;
      vul_t cur_11 = cur_geno_vword & m1;
      cur_11 = cur_11 & vul_rshift(cur_geno_vword, 1);
      count2 = count2 + cur_11;

      cur_geno_vword = *geno_vvec_iter++;
      cur_11 = cur_geno_vword & m1;
      cur_11 = cur_11 & vul_rshift(cur_geno_vword, 1);
      count2 = count2 + cur_11;
      count1 = count1 + (count2 & m2) + (vul_rshift(count2, 2) & m2);

      acc.vi = acc.vi + (count1 & m4) + (vul_rshift(count1, 4) & m4);
    } while (geno_vvec_iter < geno_vvec_stop);
    acc.vi = (acc.vi & m8) + (vul_rshift(acc.vi, 8) & m8);
    tot += univec_hsum_16bit(acc);
  }
}

uintptr_t count_11_longs(const uintptr_t* genovec, uintptr_t word_ct) {
  uintptr_t tot = 0;
  if (word_ct >= (6 * kWordsPerVec)) {
    assert(IS_VEC_ALIGNED(genovec));
    const uintptr_t remainder = word_ct % (6 * kWordsPerVec);
    const uintptr_t main_block_word_ct = word_ct - remainder;
    tot = count_11_vecs((const vul_t*)genovec, main_block_word_ct / kWordsPerVec);
    word_ct = remainder;
    genovec = &(genovec[main_block_word_ct]);
  }
  for (uintptr_t trailing_word_idx = 0; trailing_word_idx < word_ct; ++trailing_word_idx) {
    const uintptr_t cur_geno_word = genovec[trailing_word_idx];
    tot += popcount01_long(cur_geno_word & (cur_geno_word >> 1) & kMask5555);
  }
}
*/

void interleaved_mask_zero(const uintptr_t* __restrict interleaved_mask, uintptr_t vec_ct, uintptr_t* __restrict genovec) {
  const uintptr_t twovec_ct = vec_ct / 2;
#ifdef __LP64__
  const vul_t m1 = VCONST_UL(kMask5555);
  const vul_t* interleaved_mask_iter = (const vul_t*)interleaved_mask;
  vul_t* genovvec_iter = (vul_t*)genovec;
  for (uintptr_t twovec_idx = 0; twovec_idx < twovec_ct; ++twovec_idx) {
    const vul_t mask_vvec = *interleaved_mask_iter++;
    vul_t mask_first = mask_vvec & m1;
    mask_first = mask_first | vul_lshift(mask_first, 1);
    vul_t mask_second = (~m1) & mask_vvec;
    mask_second = mask_second | vul_rshift(mask_second, 1);
    *genovvec_iter = (*genovvec_iter) & mask_first;
    ++genovvec_iter;
    *genovvec_iter = (*genovvec_iter) & mask_second;
    ++genovvec_iter;
  }
  if (vec_ct & 1) {
    vul_t mask_first = *interleaved_mask_iter;
    mask_first = mask_first | vul_lshift(mask_first, 1);
    *genovvec_iter = (*genovvec_iter) & mask_first;
  }
#else
  const uintptr_t* interleaved_mask_iter = interleaved_mask;
  uintptr_t* genovec_iter = genovec;
  for (uintptr_t twovec_idx = 0; twovec_idx < twovec_ct; ++twovec_idx) {
    const uintptr_t mask_word = *interleaved_mask_iter++;
    *genovec_iter &= (mask_word & kMask5555) * 3;
    ++genovec_iter;
    *genovec_iter &= ((mask_word >> 1) & kMask5555) * 3;
    ++genovec_iter;
  }
  if (vec_ct & 1) {
    const uintptr_t mask_word = *interleaved_mask_iter;
    *genovec_iter &= mask_word * 3;
  }
#endif
}

void interleaved_set_missing(const uintptr_t* __restrict interleaved_set, uintptr_t vec_ct, uintptr_t* __restrict genovec) {
  const uintptr_t twovec_ct = vec_ct / 2;
#ifdef __LP64__
  const vul_t m1 = VCONST_UL(kMask5555);
  const vul_t* interleaved_set_iter = (const vul_t*)interleaved_set;
  vul_t* genovvec_iter = (vul_t*)genovec;
  for (uintptr_t twovec_idx = 0; twovec_idx < twovec_ct; ++twovec_idx) {
    const vul_t set_vvec = *interleaved_set_iter++;
    vul_t set_first = set_vvec & m1;
    set_first = set_first | vul_lshift(set_first, 1);
    vul_t set_second = (~m1) & set_vvec;
    set_second = set_second | vul_rshift(set_second, 1);
    *genovvec_iter = (*genovvec_iter) | set_first;
    ++genovvec_iter;
    *genovvec_iter = (*genovvec_iter) | set_second;
    ++genovvec_iter;
  }
  if (vec_ct & 1) {
    vul_t set_first = *interleaved_set_iter;
    set_first = set_first | vul_lshift(set_first, 1);
    *genovvec_iter = (*genovvec_iter) | set_first;
  }
#else
  const uintptr_t* interleaved_set_iter = interleaved_set;
  uintptr_t* genovec_iter = genovec;
  for (uintptr_t twovec_idx = 0; twovec_idx < twovec_ct; ++twovec_idx) {
    const uintptr_t set_word = *interleaved_set_iter++;
    *genovec_iter |= (set_word & kMask5555) * 3;
    ++genovec_iter;
    *genovec_iter |= ((set_word >> 1) & kMask5555) * 3;
    ++genovec_iter;
  }
  if (vec_ct & 1) {
    const uintptr_t set_word = *interleaved_set_iter;
    *genovec_iter |= set_word * 3;
  }
#endif
}

void set_male_het_missing(const uintptr_t* __restrict sex_male_interleaved, uint32_t vec_ct, uintptr_t* __restrict genovec) {
  const uint32_t twovec_ct = vec_ct / 2;
#ifdef __LP64__
  const vul_t m1 = VCONST_UL(kMask5555);
  const vul_t* sex_male_interleaved_iter = (const vul_t*)sex_male_interleaved;
  vul_t* genovvec_iter = (vul_t*)genovec;
  for (uint32_t twovec_idx = 0; twovec_idx < twovec_ct; ++twovec_idx) {
    const vul_t sex_male_vvec = *sex_male_interleaved_iter++;
    // we wish to bitwise-or with (sex_male_quatervec_01 & genovec) << 1
    const vul_t sex_male_first = sex_male_vvec & m1;
    const vul_t sex_male_second_shifted = (~m1) & sex_male_vvec;
    vul_t cur_geno_vword = *genovvec_iter;

    const vul_t missing_male_vword = sex_male_first & cur_geno_vword;

    *genovvec_iter++ = cur_geno_vword | vul_lshift(missing_male_vword, 1);
    cur_geno_vword = *genovvec_iter;
    *genovvec_iter++ = cur_geno_vword | (sex_male_second_shifted & vul_lshift(cur_geno_vword, 1));
  }
  if (vec_ct & 1) {
    const vul_t sex_male_first = (*sex_male_interleaved_iter) & m1;
    const vul_t cur_geno_vword = *genovvec_iter;
    const vul_t missing_male_vword = sex_male_first & cur_geno_vword;
    *genovvec_iter = cur_geno_vword | vul_lshift(missing_male_vword, 1);
  }
#else
  const uintptr_t* sex_male_interleaved_iter = sex_male_interleaved;
  uintptr_t* genovec_iter = genovec;
  for (uint32_t twovec_idx = 0; twovec_idx < twovec_ct; ++twovec_idx) {
    const uintptr_t sex_male_word = *sex_male_interleaved_iter++;
    uintptr_t cur_geno_word = *genovec_iter;
    *genovec_iter++ = cur_geno_word | ((sex_male_word & kMask5555 & cur_geno_word) << 1);
    cur_geno_word = *genovec_iter;
    *genovec_iter++ = cur_geno_word | (sex_male_word & kMaskAAAA & (cur_geno_word << 1));
  }
  if (vec_ct & 1) {
    const uintptr_t sex_male_word = *sex_male_interleaved_iter;
    uintptr_t cur_geno_word = *genovec_iter;
    *genovec_iter = cur_geno_word | ((sex_male_word & kMask5555 & cur_geno_word) << 1);
  }
#endif
}

// Clears each bit in bitarr which doesn't correspond to a genovec het.
// Assumes that either trailing bits of bitarr are already zero, or trailing
// bits of genovec are zero.
//
// Similar to pgr_detect_genovec_hets_unsafe().
void mask_genovec_hets_unsafe(const uintptr_t* __restrict genovec, uint32_t raw_sample_ctl2, uintptr_t* __restrict bitarr) {
  halfword_t* bitarr_alias = (halfword_t*)bitarr;
  for (uint32_t widx = 0; widx < raw_sample_ctl2; ++widx) {
    const uintptr_t cur_word = genovec[widx];
    uintptr_t ww = (~(cur_word >> 1)) & cur_word & kMask5555; // low 1, high 0
    bitarr_alias[widx] &= pack_word_to_halfword(ww);
  }
}

/*
uint32_t chr_window_max(const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_bp, uint32_t chr_fo_idx, uint32_t ct_max, uint32_t bp_max, uint32_t cur_window_max) {
  if (cur_window_max >= ct_max) {
    return ct_max;
  }
  const uint32_t chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
  uint32_t variant_uidx = next_set(variant_include, cip->chr_fo_vidx_start[chr_fo_idx], chr_end);
  const uint32_t variant_ct = popcount_bit_idx(variant_include, variant_uidx, chr_end);
  if (variant_ct <= cur_window_max) {
    return cur_window_max;
  }
  uint32_t window_idx_first = 0;
  uint32_t window_uidx_first = variant_uidx;
  uint32_t window_bp_first = variant_bp[variant_uidx];
  for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_uidx, ++variant_idx) {
    next_set_unsafe_ck(variant_include, &variant_uidx);
    uint32_t variant_bp_thresh = variant_bp[variant_uidx];
    if (variant_bp_thresh < bp_max) {
      variant_bp_thresh = 0;
    } else {
      variant_bp_thresh -= bp_max;
    }
    if (variant_bp_thresh > window_bp_first) {
      do {
        ++window_uidx_first;
        next_set_unsafe_ck(variant_include, &window_uidx_first);
        window_bp_first = variant_bp[window_uidx_first];
        ++window_idx_first;
      } while (variant_bp_thresh > window_bp_first);
    } else if (variant_idx - window_idx_first == cur_window_max) {
      if (++cur_window_max == ct_max) {
        return cur_window_max;
      }
    }
  }
  return cur_window_max;
}
*/

uint32_t not_only_xymt(const uintptr_t* variant_include, const chr_info_t* cip, uint32_t raw_variant_ct, uint32_t xymt_offset) {
  const uint32_t xymt_code = (uint32_t)cip->xymt_codes[xymt_offset];
  const uint32_t cur_chr_fo_idx = cip->chr_idx_to_foidx[xymt_code];
  const uint32_t chr_start = cip->chr_fo_vidx_start[cur_chr_fo_idx];
  if (chr_start) {
    const uint32_t first_uidx = next_set_unsafe(variant_include, 0);
    if (first_uidx < chr_start) {
      return 1;
    }
  }
  const uint32_t chr_end = cip->chr_fo_vidx_start[cur_chr_fo_idx + 1];
  return (chr_end < raw_variant_ct) && (next_set(variant_include, chr_end, raw_variant_ct) != raw_variant_ct);
}

uint32_t count_non_autosomal_variants(const uintptr_t* variant_include, const chr_info_t* cip, uint32_t count_x, uint32_t count_mt) {
  // for backward compatibility, unplaced markers are considered to be
  // autosomal here
  uint32_t ct = 0;
  if (count_x) {
    int32_t x_code;
    if (xymt_exists(cip, kChrOffsetX, &x_code)) {
      ct += count_chr_variants_unsafe(variant_include, cip, x_code);
    }
  }
  int32_t y_code;
  if (xymt_exists(cip, kChrOffsetY, &y_code)) {
    ct += count_chr_variants_unsafe(variant_include, cip, y_code);
  }
  if (count_mt) {
    int32_t mt_code;
    if (xymt_exists(cip, kChrOffsetMT, &mt_code)) {
      ct += count_chr_variants_unsafe(variant_include, cip, mt_code);
    }
  }
  return ct;
}

pglerr_t conditional_allocate_non_autosomal_variants(const chr_info_t* cip, const char* calc_descrip, uint32_t raw_variant_ct, uintptr_t** variant_include_ptr, uint32_t* variant_ct_ptr) {
  const uint32_t non_autosomal_variant_ct = count_non_autosomal_variants(*variant_include_ptr, cip, 1, 1);
  if (!non_autosomal_variant_ct) {
    return kPglRetSuccess;
  }
  LOGPRINTF("Excluding %u variant%s on non-autosomes from %s.\n", non_autosomal_variant_ct, (non_autosomal_variant_ct == 1)? "" : "s", calc_descrip);
  *variant_ct_ptr -= non_autosomal_variant_ct;
  if (!(*variant_ct_ptr)) {
    // this may not always be an error condition
    LOGERRPRINTF("Error: No variants remaining for %s.\n", calc_descrip);
    return kPglRetInconsistentInput;
  }
  const uint32_t raw_variant_ctl = BITCT_TO_WORDCT(raw_variant_ct);
  uintptr_t* working_variant_include;
  if (bigstack_alloc_ul(raw_variant_ctl, &working_variant_include)) {
    return kPglRetNomem;
  }
  memcpy(working_variant_include, *variant_include_ptr, raw_variant_ctl * sizeof(intptr_t));
  int32_t x_code;
  if (xymt_exists(cip, kChrOffsetX, &x_code)) {
    uint32_t chr_fo_idx = cip->chr_idx_to_foidx[(uint32_t)x_code];
    clear_bits_nz(cip->chr_fo_vidx_start[chr_fo_idx], cip->chr_fo_vidx_start[chr_fo_idx + 1], working_variant_include);
  }
  int32_t y_code;
  if (xymt_exists(cip, kChrOffsetX, &y_code)) {
    uint32_t chr_fo_idx = cip->chr_idx_to_foidx[(uint32_t)y_code];
    clear_bits_nz(cip->chr_fo_vidx_start[chr_fo_idx], cip->chr_fo_vidx_start[chr_fo_idx + 1], working_variant_include);
  }
  int32_t mt_code;
  if (xymt_exists(cip, kChrOffsetX, &mt_code)) {
    uint32_t chr_fo_idx = cip->chr_idx_to_foidx[(uint32_t)mt_code];
    clear_bits_nz(cip->chr_fo_vidx_start[chr_fo_idx], cip->chr_fo_vidx_start[chr_fo_idx + 1], working_variant_include);
  }
  *variant_include_ptr = working_variant_include;
  return kPglRetSuccess;
}

void fill_subset_chr_fo_vidx_start(const uintptr_t* variant_include, const chr_info_t* cip, uint32_t* subset_chr_fo_vidx_start) {
  const uint32_t chr_ct = cip->chr_ct;
  subset_chr_fo_vidx_start[0] = 0;
  uint32_t variant_uidx = 0;
  uint32_t variant_idx = 0;
  for (uint32_t chr_fo_idx = 1; chr_fo_idx <= chr_ct; ++chr_fo_idx) {
    const uint32_t chr_end_variant_uidx = cip->chr_fo_vidx_start[chr_fo_idx];
    variant_idx += popcount_bit_idx(variant_include, variant_uidx, chr_end_variant_uidx);
    subset_chr_fo_vidx_start[chr_fo_idx] = variant_idx;
    variant_uidx = chr_end_variant_uidx;
  }
}

boolerr_t allele_set(const char* newval, uint32_t allele_slen, char** allele_ptr) {
  char* newptr;
  if (allele_slen == 1) {
    // const_cast
    newptr = (char*)((uintptr_t)(&(g_one_char_strs[((unsigned char)(*newval)) * 2])));
  } else {
    char* new_alloc;
    if (pgl_malloc(allele_slen + 1, &new_alloc)) {
      return 1;
    }
    memcpyx(new_alloc, newval, allele_slen, '\0');
    newptr = new_alloc;
  }
  *allele_ptr = newptr;
  return 0;
}

boolerr_t allele_reset(const char* newval, uint32_t allele_slen, char** allele_ptr) {
  char* newptr;
  if (allele_slen == 1) {
    // const_cast
    newptr = (char*)((uintptr_t)(&(g_one_char_strs[((unsigned char)(*newval)) * 2])));
  } else {
    char* new_alloc;
    if (pgl_malloc(allele_slen + 1, &new_alloc)) {
      return 1;
    }
    memcpyx(new_alloc, newval, allele_slen, '\0');
    newptr = new_alloc;
  }
  const uintptr_t bigstack_end_addr = (uintptr_t)g_bigstack_end;
  const uintptr_t maxdiff = ((uintptr_t)(&(g_one_char_strs[512]))) - bigstack_end_addr;
  // take advantage of unsigned wraparound
  if ((((uintptr_t)(*allele_ptr)) - bigstack_end_addr) >= maxdiff) {
    free(*allele_ptr);
  }
  *allele_ptr = newptr;
  return 0;
}

void cleanup_allele_storage(uint32_t max_allele_slen, uintptr_t allele_storage_entry_ct, char** allele_storage) {
  // Now doesn't improperly free bigstack allocations (as long as they aren't
  // past g_bigstack_end), and doesn't need to be called at all most of the
  // time.

  // An alternative representation: have a separate bitarray which indicates
  // whether the allele_storage[] element should be interpreted as a heap
  // pointer or an in-place zero-terminated string (i.e. string length can be
  // up to 7 on 64-bit systems).  I expect that to be more efficient for new
  // datasets, but let's get the simple (and 1.9-codebase-compatible)
  // implementation working first, and then benchmark the fancier code later.
  if (allele_storage && (max_allele_slen > 1)) {
    const uintptr_t bigstack_end_addr = (uintptr_t)g_bigstack_end;
    const uintptr_t maxdiff = ((uintptr_t)(&(g_one_char_strs[512]))) - bigstack_end_addr;
    for (uintptr_t idx = 0; idx < allele_storage_entry_ct; ++idx) {
      char* cur_entry = allele_storage[idx];
      assert(cur_entry);
      // take advantage of unsigned wraparound
      if ((((uintptr_t)cur_entry) - bigstack_end_addr) >= maxdiff) {
        free(cur_entry);
      }
    }
  }
}

char g_missing_catname[kMaxMissingPhenostrBlen];
char g_output_missing_pheno[kMaxMissingPhenostrBlen];
char g_legacy_output_missing_pheno[kMaxMissingPhenostrBlen];

void init_pheno() {
  strcpy(g_missing_catname, "NONE");
  strcpy(g_output_missing_pheno, "NA");
  strcpy(g_legacy_output_missing_pheno, "-9");
}

uint32_t is_categorical_phenostr(const char* phenostr) {
  uint32_t first_char_code = (unsigned char)(*phenostr++);
  // allow leading +/-
  if ((first_char_code == 43) || (first_char_code == 45)) {
    first_char_code = (unsigned char)(*phenostr++);
  }
  if (((first_char_code - 48) < 10) || (first_char_code == 44) || (first_char_code < 32)) {
    // the last two conditions are for detecting CSV empty strings
    return 0;
  }
  if (first_char_code == 46) {
    // decimal point.  classify based on whether next character is a digit.
    const uint32_t second_char_code = (unsigned char)phenostr[0];
    return ((second_char_code - 48) >= 10);
  }
  // allow any capitalization of "NA"/"nan", but not "inf"
  if ((first_char_code & 0xdf) != 78) {
    return 1;
  }
  const uint32_t second_char_code = (unsigned char)phenostr[0];
  if ((second_char_code & 0xdf) != 65) {
    return 1;
  }
  const uint32_t third_char_code = (unsigned char)phenostr[1];
  if ((third_char_code & 0xdf) == 78) {
    return (((unsigned char)phenostr[2]) > ' ');
  }
  return (third_char_code > 32);
}

uint32_t is_categorical_phenostr_nocsv(const char* phenostr) {
  uint32_t first_char_code = (unsigned char)(*phenostr++);
  // allow leading +/-
  if ((first_char_code == 43) || (first_char_code == 45)) {
    first_char_code = (unsigned char)(*phenostr++);
  }
  if ((first_char_code - 48) < 10) {
    return 0;
  }
  if (first_char_code == 46) {
    // decimal point.  classify based on whether next character is a digit.
    const uint32_t second_char_code = (unsigned char)phenostr[0];
    return ((second_char_code - 48) >= 10);
  }
  // allow any capitalization of "NA"/"nan", but not "inf"
  if ((first_char_code & 0xdf) != 78) {
    return 1;
  }
  const uint32_t second_char_code = (unsigned char)phenostr[0];
  if ((second_char_code & 0xdf) != 65) {
    return 1;
  }
  const uint32_t third_char_code = (unsigned char)phenostr[1];
  if ((third_char_code & 0xdf) == 78) {
    return (((unsigned char)phenostr[2]) > ' ');
  }
  return (third_char_code > 32);
}

uint32_t first_cc_or_qt_pheno_idx(const pheno_col_t* pheno_cols, uint32_t pheno_ct) {
  for (uint32_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
    if (pheno_cols[pheno_idx].type_code < kPhenoDtypeCat) {
      return pheno_idx;
    }
  }
  return 0xffffffffU;
}

uint32_t is_const_covar(const pheno_col_t* covar_col, const uintptr_t* sample_include, uint32_t sample_ct) {
  if (sample_ct < 2) {
    return 1;
  }
  uint32_t sample_uidx = next_set_unsafe(sample_include, 0);
  if (covar_col->type_code == kPhenoDtypeQt) {
    const double* covar_vals = covar_col->data.qt;
    const double first_covar_val = covar_vals[sample_uidx];
    for (uint32_t sample_idx = 1; sample_idx < sample_ct; ++sample_idx) {
      ++sample_uidx;
      next_set_unsafe_ck(sample_include, &sample_uidx);
      if (covar_vals[sample_uidx] != first_covar_val) {
        return 0;
      }
    }
    return 1;
  }
  assert(covar_col->type_code == kPhenoDtypeCat);
  const uint32_t* covar_vals = covar_col->data.cat;
  const uint32_t first_covar_val = covar_vals[sample_uidx];
  for (uint32_t sample_idx = 1; sample_idx < sample_ct; ++sample_idx) {
    ++sample_uidx;
    next_set_unsafe_ck(sample_include, &sample_uidx);
    if (covar_vals[sample_uidx] != first_covar_val) {
      return 0;
    }
  }
  return 1;
}

uint32_t identify_remaining_cats(const uintptr_t* sample_include, const pheno_col_t* covar_col, uint32_t sample_ct, uintptr_t* observed_cat_bitarr) {
  // assumes covar_col->type_code == kPhenoTypeCat
  const uint32_t nonnull_cat_ct = covar_col->nonnull_category_ct;
  const uint32_t* covar_vals = covar_col->data.cat;
  const uint32_t word_ct = 1 + (nonnull_cat_ct / kBitsPerWord);
  fill_ulong_zero(word_ct, observed_cat_bitarr);
  uint32_t sample_uidx = 0;
  for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
    next_set_unsafe_ck(sample_include, &sample_uidx);
    set_bit(covar_vals[sample_uidx], observed_cat_bitarr);
  }
  return popcount_longs(observed_cat_bitarr, word_ct);
}

uint32_t get_is_cat_include(const uintptr_t* sample_include_base, const pheno_col_t* cat_pheno_col, uint32_t raw_sample_ctl, uint32_t sample_ct, uint32_t cat_uidx, uintptr_t* is_cat_include) {
  fill_ulong_zero(raw_sample_ctl, is_cat_include);
  const uint32_t* cat_vals = cat_pheno_col->data.cat;
  uint32_t sample_uidx = 0;
  for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
    next_set_unsafe_ck(sample_include_base, &sample_uidx);
    if (cat_vals[sample_uidx] == cat_uidx) {
      set_bit(sample_uidx, is_cat_include);
    }
  }
  return popcount_longs(is_cat_include, raw_sample_ctl);
}

void cleanup_pheno_cols(uint32_t pheno_ct, pheno_col_t* pheno_cols) {
  if (pheno_cols) {
    for (uint32_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
      vecaligned_free_cond(pheno_cols[pheno_idx].nonmiss);
    }
    free(pheno_cols);
  }
}

pglerr_t parse_chr_ranges(const char* flagname_p, const char* errstr_append, uint32_t param_ct, uint32_t allow_extra_chrs, uint32_t xymt_subtract, char range_delim, char** argv, chr_info_t* cip, uintptr_t* chr_mask) {
  pglerr_t reterr = kPglRetSuccess;
  {
    char* cur_arg_ptr = argv[1];
    char* range_end = nullptr;
    uint32_t cur_param_idx = 1;
    uint32_t rs_len = 0;
    uint32_t re_len = 0;
    while (1) {
      char* range_start;
      if (parse_next_range(argv, param_ct, range_delim, &cur_param_idx, &cur_arg_ptr, &range_start, &rs_len, &range_end, &re_len)) {
        sprintf(g_logbuf, "Error: Invalid --%s parameter '%s'.\n", flagname_p, argv[cur_param_idx]);
        goto parse_chr_ranges_ret_INVALID_CMDLINE_WWA;
      }
      if (!range_start) {
        break;
      }
      const char cc = range_start[rs_len];
      range_start[rs_len] = '\0';
      int32_t chr_code_start = get_chr_code_raw(range_start);
      if (chr_code_start < 0) {
        if (!allow_extra_chrs) {
          sprintf(g_logbuf, "Error: Invalid --%s chromosome code '%s'.\n", flagname_p, range_start);
          goto parse_chr_ranges_ret_INVALID_CMDLINE_WWA;
        }
        if (range_end) {
          goto parse_chr_ranges_ret_INVALID_CMDLINE_NONSTD;
        }
        if (push_llstr(range_start, &(cip->incl_excl_name_stack))) {
          goto parse_chr_ranges_ret_NOMEM;
        }
      } else {
        if (chr_code_start >= ((int32_t)kMaxContigs)) {
          chr_code_start -= xymt_subtract;
        }
        if (range_end) {
          const char cc2 = range_end[re_len];
          range_end[re_len] = '\0';
          int32_t chr_code_end = get_chr_code_raw(range_end);
          if (chr_code_end < 0) {
            if (!allow_extra_chrs) {
              sprintf(g_logbuf, "Error: Invalid --%s chromosome code '%s'.\n", flagname_p, range_end);
              goto parse_chr_ranges_ret_INVALID_CMDLINE_WWA;
            }
            goto parse_chr_ranges_ret_INVALID_CMDLINE_NONSTD;
          }
          if (chr_code_end >= ((int32_t)kMaxContigs)) {
            // prohibit stuff like "--chr par1-par2", "--chr x-y", "--chr x-26"
            sprintf(g_logbuf, "Error: --%s chromosome code '%s' cannot be the end of a range.\n", flagname_p, range_end);
            goto parse_chr_ranges_ret_INVALID_CMDLINE_WWA;
          }
          if (chr_code_end <= chr_code_start) {
            sprintf(g_logbuf, "Error: --%s chromosome code '%s' is not greater than '%s'.\n", flagname_p, range_end, range_start);
            goto parse_chr_ranges_ret_INVALID_CMDLINE_WWA;
          }
          range_end[re_len] = cc2;
          fill_bits_nz(chr_code_start, chr_code_end + 1, chr_mask);
        } else {
          set_bit(chr_code_start, chr_mask);
        }
      }
      range_start[rs_len] = cc;
    }
    // no compelling reason to prohibit "--not-chr ,"
  }
  while (0) {
  parse_chr_ranges_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  parse_chr_ranges_ret_INVALID_CMDLINE_NONSTD:
    logerrprint("Error: Chromosome ranges cannot include nonstandard names.\n");
    reterr = kPglRetInvalidCmdline;
    break;
  parse_chr_ranges_ret_INVALID_CMDLINE_WWA:
    wordwrapb(0);
    logerrprintb();
    logerrprint(errstr_append);
    reterr = kPglRetInvalidCmdline;
    break;
  }
  return reterr;
}


pglerr_t multithread_load_init(const uintptr_t* variant_include, uint32_t sample_ct, uint32_t variant_ct, uintptr_t pgr_alloc_cacheline_ct, uintptr_t thread_xalloc_cacheline_ct, uintptr_t per_variant_xalloc_byte_ct, pgen_file_info_t* pgfip, uint32_t* calc_thread_ct_ptr, uintptr_t*** genovecs_ptr, uintptr_t*** dosage_present_ptr, dosage_t*** dosage_val_bufs_ptr, uint32_t* read_block_size_ptr, unsigned char** main_loadbufs, pthread_t** threads_ptr, pgen_reader_t*** pgr_pps, uint32_t** read_variant_uidx_starts_ptr) {
  uintptr_t cachelines_avail = bigstack_left() / kCacheline;
  uint32_t read_block_size = kPglVblockSize;
  uint64_t multiread_cacheline_ct;
  while (1) {
    multiread_cacheline_ct = pgfi_multiread_get_cacheline_req(variant_include, pgfip, variant_ct, read_block_size);
    // limit each raw load buffer to 1/4 of remaining workspace
    // if there's an additional per-variant allocation, put it in the same bin
    // as the load buffers
    if ((multiread_cacheline_ct + (((uint64_t)per_variant_xalloc_byte_ct) * read_block_size) / kCacheline) * 4 <= cachelines_avail) {
      break;
    }
    // lots of callers require read_block_size to be either raw_variant_ct or a
    // multiple of kBitsPerVec
#ifdef __LP64__
    if (read_block_size <= kBitsPerVec) {
      return kPglRetNomem;
    }
#else
    if (read_block_size <= kCacheline) {
      return kPglRetNomem;
    }
#endif
    read_block_size /= 2;
  }
#ifndef __LP64__
  if (multiread_cacheline_ct > (kMaxBytesPerIO / kCacheline)) {
    return kPglRetNomem;
  }
#endif
  main_loadbufs[0] = bigstack_alloc_raw(multiread_cacheline_ct * kCacheline);
  main_loadbufs[1] = bigstack_alloc_raw(multiread_cacheline_ct * kCacheline);
  pgfip->block_base = main_loadbufs[0];
  *read_block_size_ptr = read_block_size;
  cachelines_avail -= 2 * (multiread_cacheline_ct + (((uint64_t)per_variant_xalloc_byte_ct) * read_block_size) / kCacheline);
  // reduce calc_thread_ct if necessary
  uint32_t calc_thread_ct = *calc_thread_ct_ptr;
  if (calc_thread_ct > read_block_size) {
    calc_thread_ct = read_block_size;
    *calc_thread_ct_ptr = calc_thread_ct;
  }

  // pgr_pps, threads_ptr, read_variant_uidx_starts_ptr, (*pgr_pps)[tidx],
  //   pgr_alloc; deliberately a slight overestimate
  const uintptr_t pgr_struct_alloc = round_up_pow2(sizeof(pgen_reader_t), kCacheline);
  uintptr_t thread_alloc_cacheline_ct = 1 + 1 + 1 + (pgr_struct_alloc / kCacheline) + pgr_alloc_cacheline_ct + thread_xalloc_cacheline_ct;

  const uint32_t sample_ctcl2 = QUATERCT_TO_CLCT(sample_ct);
  const uint32_t sample_ctcl = BITCT_TO_CLCT(sample_ct);

  // todo: increase in multiallelic case
  const uintptr_t dosage_vals_cl = DIV_UP(sample_ct, (kCacheline / sizeof(dosage_t)));
  if (genovecs_ptr) {
    thread_alloc_cacheline_ct += 1 + sample_ctcl2;
    if (dosage_present_ptr) {
      assert(dosage_val_bufs_ptr);
      thread_alloc_cacheline_ct += 2 + sample_ctcl + dosage_vals_cl;
    }
  }
  if (thread_alloc_cacheline_ct * calc_thread_ct > cachelines_avail) {
    if (thread_alloc_cacheline_ct > cachelines_avail) {
      return kPglRetNomem;
    }
    calc_thread_ct = cachelines_avail / thread_alloc_cacheline_ct;
    *calc_thread_ct_ptr = calc_thread_ct;
  }

  const uint32_t array_of_ptrs_alloc = round_up_pow2(calc_thread_ct * sizeof(intptr_t), kCacheline);
  *pgr_pps = (pgen_reader_t**)bigstack_alloc_raw(array_of_ptrs_alloc);
  *threads_ptr = (pthread_t*)bigstack_alloc_raw(array_of_ptrs_alloc);
  *read_variant_uidx_starts_ptr = (uint32_t*)bigstack_alloc_raw_rd(calc_thread_ct * sizeof(int32_t));
  for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
    (*pgr_pps)[tidx] = (pgen_reader_t*)bigstack_alloc_raw(pgr_struct_alloc);
    // pgr_preinit(g_pgr_ptrs[tidx]);
    unsigned char* pgr_alloc = bigstack_alloc_raw(pgr_alloc_cacheline_ct * kCacheline);

    // shouldn't be possible for this to fail
    pgr_init(nullptr, 0, pgfip, (*pgr_pps)[tidx], pgr_alloc);
  }
  if (genovecs_ptr) {
    *genovecs_ptr = (uintptr_t**)bigstack_alloc_raw(array_of_ptrs_alloc);
    if (dosage_present_ptr) {
      *dosage_present_ptr = (uintptr_t**)bigstack_alloc_raw(array_of_ptrs_alloc);
      *dosage_val_bufs_ptr = (dosage_t**)bigstack_alloc_raw(array_of_ptrs_alloc);
    }
    const uintptr_t genovec_alloc = sample_ctcl2 * kCacheline;
    const uintptr_t dosage_present_alloc = sample_ctcl * kCacheline;
    const uintptr_t dosage_vals_alloc = dosage_vals_cl * kCacheline;
    for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
      (*genovecs_ptr)[tidx] = (uintptr_t*)bigstack_alloc_raw(genovec_alloc);
      if (dosage_present_ptr) {
        (*dosage_present_ptr)[tidx] = (uintptr_t*)bigstack_alloc_raw(dosage_present_alloc);
        (*dosage_val_bufs_ptr)[tidx] = (dosage_t*)bigstack_alloc_raw(dosage_vals_alloc);
      }
    }
  }
  return kPglRetSuccess;
}

pglerr_t write_sample_ids(const uintptr_t* sample_include, const char* sample_ids, const char* sids, const char* outname, uint32_t sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen) {
  FILE* outfile = nullptr;
  pglerr_t reterr = kPglRetSuccess;
  {
    if (fopen_checked(outname, FOPEN_WB, &outfile)) {
      goto write_sample_ids_ret_OPEN_FAIL;
    }
    char* write_iter = g_textbuf;
    char* textbuf_flush = &(write_iter[kMaxMediumLine]);
    uintptr_t sample_uidx = 0;
    for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
      next_set_ul_unsafe_ck(sample_include, &sample_uidx);
      write_iter = strcpya(write_iter, &(sample_ids[sample_uidx * max_sample_id_blen]));
      if (sids) {
        *write_iter++ = '\t';
        write_iter = strcpya(write_iter, &(sids[sample_uidx * max_sid_blen]));
      }
      append_binary_eoln(&write_iter);
      if (fwrite_ck(textbuf_flush, outfile, &write_iter)) {
        goto write_sample_ids_ret_WRITE_FAIL;
      }
    }
    if (fclose_flush_null(textbuf_flush, write_iter, &outfile)) {
      goto write_sample_ids_ret_WRITE_FAIL;
    }
  }
  while (0) {
  write_sample_ids_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  write_sample_ids_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
  fclose_cond(outfile);
  return reterr;
}

#ifdef __cplusplus
} // namespace plink2
#endif
