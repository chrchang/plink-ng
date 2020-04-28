// This file is part of PLINK 1.90, copyright (C) 2005-2020 Shaun Purcell,
// Christopher Chang.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


#include "plink_common.h"

#include "plink_homozyg.h"

// This implementation's use of unfiltered indices is idiotic; it should be
// rewritten to be based on post-filtering variant indices ASAP.

void homozyg_init(Homozyg_info* homozyg_ptr) {
  homozyg_ptr->modifier = 0;
  homozyg_ptr->min_snp = 100;
  homozyg_ptr->min_bases = 1000000;
  homozyg_ptr->max_bases_per_snp = 50000.0 + EPSILON;
  homozyg_ptr->max_hets = 0xffffffffU;
  homozyg_ptr->max_gap = 1000000;
  homozyg_ptr->window_size = 50;
  homozyg_ptr->window_max_hets = 1;
  homozyg_ptr->window_max_missing = 5;
  homozyg_ptr->hit_threshold = 0.05;
  homozyg_ptr->overlap_min = 0.95;
  homozyg_ptr->pool_size_min = 2;
}

void mask_out_homozyg_major(uintptr_t* readbuf_cur, uint32_t sample_ct) {
  // if readbuf_cur were 16-byte aligned, this could be vectorized, but it
  // isn't, and this isn't a limiting step anyway
  uintptr_t* readbuf_cur_end = &(readbuf_cur[QUATERCT_TO_WORDCT(sample_ct)]);
  uintptr_t cur_word;
  do {
    cur_word = *readbuf_cur;
    *readbuf_cur &= ((cur_word ^ (cur_word >> 1)) & FIVEMASK) * (3 * ONELU);
  } while (++readbuf_cur < readbuf_cur_end);
}

void increment_het_missing(uintptr_t* readbuf_cur, uint32_t sample_ct, uint32_t* het_cts, uint32_t* missing_cts, int32_t incr) {
  // assumes homozyg majors masked out
  uint32_t sample_idx_offset;
  uintptr_t cur_word;
  uint32_t last_set_bit;
  for (sample_idx_offset = 0; sample_idx_offset < sample_ct; sample_idx_offset += BITCT2) {
    cur_word = *readbuf_cur++;
    while (cur_word) {
      last_set_bit = CTZLU(cur_word);
      if (last_set_bit & 1) {
        het_cts[sample_idx_offset + (last_set_bit / 2)] += incr;
      } else {
        missing_cts[sample_idx_offset + (last_set_bit / 2)] += incr;
      }
      cur_word &= cur_word - ONELU;
    }
  }
}

void decrement_het_missing_extend(uintptr_t* readbuf_cur, uint32_t sample_ct, uint32_t* het_cts, uint32_t* missing_cts, uint32_t* end_nonhom_uidxs, uint32_t next_uidx) {
  uint32_t sample_idx_offset;
  uint32_t sample_idx;
  uintptr_t cur_word;
  uint32_t last_set_bit;
  for (sample_idx_offset = 0; sample_idx_offset < sample_ct; sample_idx_offset += BITCT2) {
    cur_word = *readbuf_cur++;
    while (cur_word) {
      last_set_bit = CTZLU(cur_word);
      sample_idx = sample_idx_offset + (last_set_bit / 2);
      if (last_set_bit & 1) {
        het_cts[sample_idx] -= 1;
      } else {
        missing_cts[sample_idx] -= 1;
      }
      end_nonhom_uidxs[sample_idx] = next_uidx;
      cur_word &= cur_word - ONELU;
    }
  }
}

void update_end_nonhom(uintptr_t* readbuf_cur, uint32_t sample_ct, uint32_t* end_nonhom_uidxs, uint32_t next_uidx) {
  uint32_t sample_idx_offset;
  uintptr_t cur_word;
  uint32_t last_set_bit;
  for (sample_idx_offset = 0; sample_idx_offset < sample_ct; sample_idx_offset += BITCT2) {
    cur_word = *readbuf_cur++;
    while (cur_word) {
      last_set_bit = CTZLU(cur_word);
      end_nonhom_uidxs[sample_idx_offset + (last_set_bit / 2)] = next_uidx;
      cur_word &= cur_word - ONELU;
    }
  }
}

#ifdef __LP64__
#define ROH_ENTRY_INTS 7
#else
#define ROH_ENTRY_INTS 6
#endif

void roh_extend_forward(uint32_t* marker_pos, uintptr_t* marker_exclude, uint32_t marker_uidx, uint32_t nsnp_incr, uint32_t is_new_lengths, double max_bases_per_snp, uint32_t* roh_list_cur) {
  // marker_uidx should be the index of the marker immediately *after* the last
  // eligible one.
  int32_t subtract_pos = (int32_t)(marker_pos[roh_list_cur[0]] - is_new_lengths);
  uint32_t orig_nsnp = roh_list_cur[2];
  while (nsnp_incr) {
    prev_unset_unsafe_ck(marker_exclude, &marker_uidx);
    if (((int32_t)(orig_nsnp + nsnp_incr)) * max_bases_per_snp >= (double)(marker_pos[marker_uidx] - subtract_pos)) {
      roh_list_cur[1] = marker_uidx;
      roh_list_cur[2] += nsnp_incr;
      roh_list_cur[3] += nsnp_incr;
      return;
    }
    nsnp_incr--;
  }
}

void save_confirmed_roh_extend(uint32_t uidx_first, uint32_t older_uidx, uint32_t cidx_len, uint32_t cur_roh_het_ct, uint32_t cur_roh_missing_ct, uintptr_t* sample_last_roh_idx_ptr, uintptr_t* roh_ct_ptr, uint32_t* roh_list, uint32_t* marker_pos, uintptr_t* marker_exclude, uintptr_t marker_uidx, uint32_t is_new_lengths, double max_bases_per_snp, uint32_t* prev_roh_end_cidx_ptr, uint32_t* cur_roh_earliest_extend_uidx_ptr, uint32_t cidx_cur) {
  // This is straightforward when the 'extend' modifier isn't present, we just
  // append to roh_list.  It's a bit more complicated with 'extend':
  //
  // * When a new ROH is confirmed, we extend it backwards as far as possible.
  //   In a few pathological cases, this actually means we merge it with the
  //   previous ROH (consider:
  //     --homozyg-snp 2
  //     --homozyg-window-het 0
  //     --homozyg-window-threshold 0.0727
  //     chromosome with 99 loci: 24 hets, then 51 homs, then 24 hets), if that
  //   doesn't blow the --homozyg-density or --homozyg-het limit.
  // * We defer the forwards extension until the main loop passes the last
  //   adjacent homozygous call after the ROH (prev_roh_end_cidxs[] is set to
  //   0xffffffffU when this happens).  It is theoretically possible to be in
  //   another ROH candidate at this point (consider what happens to the
  //   example above when --homozyg-window-het is changed to 1, and the 74th
  //   site is changed to het), but w.r.t. the density constraint, it's safe to
  //   first extend the previous ROH forwards and then try to merge if the new
  //   ROH is confirmed.
  //   However, --homozyg-het with a nonzero parameter can lead to ambiguity:
  //     --homozyg-snp 2
  //     --homozyg-window-het 2
  //     --homozyg-het 1
  //     --homozyg-window-threshold 0.0816
  //     chromosome with 96 loci: 24 hets, then 48 homs, then 24 hets
  //   so for now we actually prohibit it from being specified with 'extend' on
  //   the command line (this shouldn't be a big loss since we do permit
  //   --homozyg-het 0).
  //
  // This is biased towards backwards extension if the ROH is near the density
  // bound, but that also shouldn't be a big deal.  (Enough minutiae for
  // now...)
  uintptr_t roh_ct = *roh_ct_ptr;
  uintptr_t sample_last_roh_idx = *sample_last_roh_idx_ptr;
  uint32_t cidx_start = cidx_cur - cidx_len;
  uint32_t* prev_roh_list_entry;
  uint32_t prev_roh_end_cidx;
  uint32_t prev_nsnp;
  uint32_t extended_nsnp;
  uint32_t add_pos;
  uint32_t cur_uidx_first;
  prev_roh_end_cidx = *prev_roh_end_cidx_ptr;
  *prev_roh_end_cidx_ptr = cidx_cur;
  add_pos = marker_pos[older_uidx] + is_new_lengths;
  if (prev_roh_end_cidx == 0xffffffffU) {
    cur_uidx_first = *cur_roh_earliest_extend_uidx_ptr;
  } else {
    // we're still in the same string of homozygous calls that the last ROH
    // ended in.  direct merge?
    prev_roh_list_entry = &(roh_list[sample_last_roh_idx * ROH_ENTRY_INTS]);
    prev_nsnp = prev_roh_list_entry[2];
    extended_nsnp = prev_nsnp + cidx_cur - prev_roh_end_cidx;
    // only thing that can prevent it at this point is --homozyg-density
    // bound violation; every other condition was checked incrementally
    // (including --homozyg-het since we've prevented it from being given any
    // value other than 0 or infinity when 'extend' is present)
    if (((double)((int32_t)extended_nsnp)) * max_bases_per_snp >= ((double)(add_pos - marker_pos[prev_roh_list_entry[0]]))) {
      // technically the merge can unlock even more backward extension, but
      // I'm just gonna a draw a line here and not check for that since we're
      // already inside a corner case that should be irrelevant in practice.
      prev_roh_list_entry[1] = older_uidx;
      prev_roh_list_entry[2] = extended_nsnp;

      // previous hom_ct = prev_nsnp - prev_het_ct - prev_missing_ct
      // -> prev_het_ct + prev_missing_ct = prev_nsnp - prev_hom_ct
      //
      // new hom_ct = extended_nsnp - cur_roh_het_ct - cur_roh_missing_ct -
      //              (prev_het_ct + prev_missing_ct)
      //            = extended_nsnp - cur_roh_het_ct - cur_roh_missing_ct -
      //              (prev_nsnp - prev_hom_ct)
      //            = prev_hom_ct + extended_nsnp - cur_roh_het_ct -
      //              cur_roh_missing_ct - prev_nsnp
      prev_roh_list_entry[3] += extended_nsnp - cur_roh_het_ct - cur_roh_missing_ct - prev_nsnp;
      prev_roh_list_entry[4] += cur_roh_het_ct;
      return;
    }
    // extend previous ROH as far as possible before current ROH's start,
    // then extend current one backward
    roh_extend_forward(marker_pos, marker_exclude, uidx_first, cidx_start - prev_roh_end_cidx, is_new_lengths, max_bases_per_snp, prev_roh_list_entry);
    cur_uidx_first = prev_roh_list_entry[1] + 1;
    next_unset_unsafe_ck(marker_exclude, &cur_uidx_first);
  }
  if (cur_uidx_first < uidx_first) {
    extended_nsnp = cidx_len + uidx_first - cur_uidx_first - popcount_bit_idx(marker_exclude, cur_uidx_first, uidx_first);
    do {
      if (((int32_t)extended_nsnp) * max_bases_per_snp >= (double)(add_pos - marker_pos[cur_uidx_first])) {
	uidx_first = cur_uidx_first;
	cidx_len = extended_nsnp;
	break;
      }
      extended_nsnp--;
      cur_uidx_first++;
      next_unset_unsafe_ck(marker_exclude, &cur_uidx_first);
    } while (cur_uidx_first < uidx_first);
  }
  roh_list = &(roh_list[roh_ct * ROH_ENTRY_INTS]);
  *roh_list++ = uidx_first;
  *roh_list++ = older_uidx;
  *roh_list++ = cidx_len;
  *roh_list++ = cidx_len - cur_roh_het_ct - cur_roh_missing_ct;
  *roh_list++ = cur_roh_het_ct;
#ifdef __LP64__
  *roh_list++ = (uint32_t)sample_last_roh_idx;
  *roh_list++ = (uint32_t)(sample_last_roh_idx >> 32);
#else
  *roh_list++ = (uint32_t)sample_last_roh_idx;
#endif
  *sample_last_roh_idx_ptr = roh_ct++;
  *roh_ct_ptr = roh_ct;
}

uint32_t roh_update(Homozyg_info* hp, uintptr_t* readbuf_cur, uintptr_t* swbuf_cur, uint32_t* het_cts, uint32_t* missing_cts, uint32_t* marker_pos, uintptr_t* cur_sample_male, uint32_t sample_ct, uint32_t swhit_min, uint32_t older_uidx, uint32_t old_uidx, uint32_t marker_cidx, uintptr_t max_roh_ct, uint32_t* swhit_cts, uint32_t* cur_roh_uidx_starts, uint32_t* cur_roh_cidx_starts, uint32_t* cur_roh_het_cts, uint32_t* cur_roh_missing_cts, uintptr_t* sample_to_last_roh, uint32_t* roh_list, uintptr_t* roh_ct_ptr, uintptr_t* marker_exclude, uint32_t* prev_roh_end_cidxs, uint32_t* end_nonhom_uidxs, uint32_t* cur_roh_earliest_extend_uidxs) {
  // Finishes processing current marker, saving all ROH that end at the
  // previous one.  If readbuf_cur is nullptr, the previous marker is assumed to
  // be the last one on the chromosome.
  uint32_t* roh_list_cur = &(roh_list[(*roh_ct_ptr) * ROH_ENTRY_INTS]);
  uint32_t min_snp = hp->min_snp;
  uint32_t min_bases = hp->min_bases;
  double max_bases_per_snp = hp->max_bases_per_snp;
  uint32_t max_hets = hp->max_hets;
  uint32_t max_sw_hets = hp->window_max_hets;
  uint32_t max_sw_missings = hp->window_max_missing;
  uint32_t is_new_lengths = 1 ^ ((hp->modifier / HOMOZYG_OLD_LENGTHS) & 1);
  uint32_t forced_end = (marker_pos[old_uidx] - marker_pos[older_uidx]) > hp->max_gap;
  uint32_t is_cur_hit = 0;
  uintptr_t cur_word = 0;
  uintptr_t cur_call;
  uintptr_t last_roh_idx;
  uint32_t sample_idx;
  uint32_t cidx_len;
  uint32_t uidx_first;
  uint32_t base_len;
  uint32_t cur_het_ct;
  if (!prev_roh_end_cidxs) {
    // ~25% performance penalty if we don't separate this out
    for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
      if (!(sample_idx & (BITCT2 - 1))) {
	if (readbuf_cur) {
	  cur_word = *readbuf_cur++;
	}
      } else {
	cur_word >>= 2;
      }
      if (cur_sample_male && IS_SET(cur_sample_male, sample_idx)) {
	// skip males on Xchr (cur_sample_male should be nullptr when not Xchr)
	continue;
      }
      if (readbuf_cur) {
	if (swbuf_cur) {
	  if ((het_cts[sample_idx] <= max_sw_hets) && (missing_cts[sample_idx] <= max_sw_missings)) {
	    SET_BIT(sample_idx, swbuf_cur);
	    swhit_cts[sample_idx] += 1;
	  }
	}
	is_cur_hit = (swhit_cts[sample_idx] >= swhit_min);
      }
      cur_call = cur_word & (3 * ONELU);
      if (cur_roh_cidx_starts[sample_idx] != 0xffffffffU) {
	if ((!is_cur_hit) || ((cur_call == 2) && (cur_roh_het_cts[sample_idx] == max_hets)) || forced_end) {
	  cidx_len = marker_cidx - cur_roh_cidx_starts[sample_idx];
	  if (cidx_len >= min_snp) {
	    uidx_first = cur_roh_uidx_starts[sample_idx];
	    base_len = marker_pos[older_uidx] + is_new_lengths - marker_pos[uidx_first];
	    if ((base_len >= min_bases) && (((double)((int32_t)cidx_len)) * max_bases_per_snp >= ((double)base_len))) {
	      if (*roh_ct_ptr == max_roh_ct) {
		return 1;
	      }
	      *roh_list_cur++ = uidx_first;
	      *roh_list_cur++ = older_uidx;
	      *roh_list_cur++ = cidx_len;
              cur_het_ct = cur_roh_het_cts[sample_idx];
	      *roh_list_cur++ = cidx_len - cur_het_ct - cur_roh_missing_cts[sample_idx];
	      *roh_list_cur++ = cur_het_ct;
	      last_roh_idx = sample_to_last_roh[sample_idx];
#ifdef __LP64__
	      *roh_list_cur++ = (uint32_t)last_roh_idx;
	      *roh_list_cur++ = (uint32_t)(last_roh_idx >> 32);
#else
	      *roh_list_cur++ = (uint32_t)last_roh_idx;
#endif
              sample_to_last_roh[sample_idx] = *roh_ct_ptr;
	      *roh_ct_ptr += 1;
	    }
	  }
	  cur_roh_cidx_starts[sample_idx] = 0xffffffffU;
	}
      }
      if (is_cur_hit) {
	if (cur_roh_cidx_starts[sample_idx] == 0xffffffffU) {
	  if ((!max_hets) && (cur_call == 2)) {
	    continue;
	  }
	  cur_roh_uidx_starts[sample_idx] = old_uidx;
	  cur_roh_cidx_starts[sample_idx] = marker_cidx;
	  cur_roh_het_cts[sample_idx] = 0;
	  cur_roh_missing_cts[sample_idx] = 0;
	}
	if (cur_call) {
	  if (cur_call == 2) {
	    cur_roh_het_cts[sample_idx] += 1;
	  } else {
	    cur_roh_missing_cts[sample_idx] += 1;
	  }
	}
      }
    }
  } else {
    if (forced_end) {
      for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
        end_nonhom_uidxs[sample_idx] = old_uidx;
      }
    }
    for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
      if (!(sample_idx & (BITCT2 - 1))) {
	if (readbuf_cur) {
	  cur_word = *readbuf_cur++;
	}
      } else {
	cur_word >>= 2;
      }
      if (cur_sample_male && IS_SET(cur_sample_male, sample_idx)) {
	continue;
      }
      if (readbuf_cur) {
	if (swbuf_cur) {
	  if ((het_cts[sample_idx] <= max_sw_hets) && (missing_cts[sample_idx] <= max_sw_missings)) {
	    SET_BIT(sample_idx, swbuf_cur);
	    swhit_cts[sample_idx] += 1;
	  }
	}
	is_cur_hit = (swhit_cts[sample_idx] >= swhit_min);
      }
      cur_call = cur_word & (3 * ONELU);
      if (cur_roh_cidx_starts[sample_idx] != 0xffffffffU) {
	if ((!is_cur_hit) || ((cur_call == 2) && (cur_roh_het_cts[sample_idx] == max_hets)) || forced_end) {
	  cidx_len = marker_cidx - cur_roh_cidx_starts[sample_idx];
	  if (cidx_len >= min_snp) {
	    uidx_first = cur_roh_uidx_starts[sample_idx];
	    base_len = marker_pos[older_uidx] + is_new_lengths - marker_pos[uidx_first];
	    if ((base_len >= min_bases) && (((double)((int32_t)cidx_len)) * max_bases_per_snp >= ((double)base_len))) {
	      if (*roh_ct_ptr == max_roh_ct) {
		return 1;
	      }
	      save_confirmed_roh_extend(uidx_first, older_uidx, cidx_len, cur_roh_het_cts[sample_idx], cur_roh_missing_cts[sample_idx], &(sample_to_last_roh[sample_idx]), roh_ct_ptr, roh_list, marker_pos, marker_exclude, old_uidx, is_new_lengths, max_bases_per_snp, &(prev_roh_end_cidxs[sample_idx]), &(cur_roh_earliest_extend_uidxs[sample_idx]), marker_cidx);
	    }
	  }
	  cur_roh_cidx_starts[sample_idx] = 0xffffffffU;
	  if (cur_call || forced_end) {
	    prev_roh_end_cidxs[sample_idx] = 0xffffffffU;
	  }
	}
      } else if ((cur_call || forced_end) && (prev_roh_end_cidxs[sample_idx] != 0xffffffffU)) {
	roh_extend_forward(marker_pos, marker_exclude, old_uidx, marker_cidx - prev_roh_end_cidxs[sample_idx], is_new_lengths, max_bases_per_snp, &(roh_list[sample_to_last_roh[sample_idx] * ROH_ENTRY_INTS]));
	prev_roh_end_cidxs[sample_idx] = 0xffffffffU;
      }
      if (is_cur_hit) {
	if (cur_roh_cidx_starts[sample_idx] == 0xffffffffU) {
	  if ((!max_hets) && (cur_call == 2)) {
	    continue;
	  }
	  cur_roh_uidx_starts[sample_idx] = old_uidx;
	  cur_roh_cidx_starts[sample_idx] = marker_cidx;
	  cur_roh_het_cts[sample_idx] = 0;
	  cur_roh_missing_cts[sample_idx] = 0;
	  if (cur_call) {
	    cur_roh_earliest_extend_uidxs[sample_idx] = old_uidx;
	    goto roh_update_extend_nonhom;
	  }
	  cur_roh_earliest_extend_uidxs[sample_idx] = end_nonhom_uidxs[sample_idx];
	} else if (cur_call) {
	roh_update_extend_nonhom:
	  if (cur_call == 2) {
	    cur_roh_het_cts[sample_idx] += 1;
	  } else {
	    cur_roh_missing_cts[sample_idx] += 1;
	  }
	}
      }
    }
  }
  return 0;
}

int32_t write_main_roh_reports(char* outname, char* outname_end, uintptr_t* marker_exclude, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, Chrom_info* chrom_info_ptr, uint32_t* marker_pos, uintptr_t sample_ct, uintptr_t* sample_exclude, char* sample_ids, uint32_t plink_maxfid, uint32_t plink_maxiid, uintptr_t max_sample_id_len, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, char* missing_pheno_str, uint32_t omp_is_numeric, uint32_t missing_pheno_len, uint32_t is_new_lengths, uintptr_t roh_ct, uint32_t* roh_list, uintptr_t* roh_list_chrom_starts, uintptr_t* sample_to_last_roh, uint32_t* max_pool_size_ptr, uint32_t* max_roh_len_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  FILE* outfile_indiv = nullptr;
  char* wptr_iid = &(g_textbuf[plink_maxfid + 1]);
  char* wptr_phe = &(g_textbuf[plink_maxfid + plink_maxiid + 2]);
  int32_t* roh_ct_aff_adj = nullptr;
  uintptr_t next_roh_idx = 0;
  uint32_t max_pool_size = 0;
  uint32_t max_roh_len = 0;
  int32_t retval = 0;
  char* cptr;
  char* cptr2;
  char* wptr_chr;
  char* wptr_bp1;
  char* wptr;
  uint32_t* cur_roh;
  int32_t* roh_ct_unaff_adj;
  uintptr_t sample_uidx;
  uintptr_t sample_idx;
  uintptr_t prev_roh_idx;
  uintptr_t cur_roh_idx;
  uintptr_t chrom_roh_start;
  uintptr_t chrom_roh_ct;
  uint32_t cur_roh_ct;
  uint32_t chrom_fo_idx;
  uint32_t marker_uidx1;
  uint32_t marker_uidx2;
  uint32_t chrom_start;
  uint32_t chrom_len;
  double dxx;
  double dyy;
  double kb_tot;
  uint32_t slen;
  uint32_t uii;
  memcpy(outname_end, ".hom", 5);
  if (fopen_checked(outname, "w", &outfile)) {
    goto write_main_roh_reports_ret_OPEN_FAIL;
  }
  sprintf(g_textbuf, "%%%us %%%us      PHE  CHR %%%us %%%us         POS1         POS2         KB     NSNP  DENSITY     PHOM     PHET\n", plink_maxfid, plink_maxiid, plink_maxsnp, plink_maxsnp);
  fprintf(outfile, g_textbuf, "FID", "IID", "SNP1", "SNP2");
  memcpy(&(outname_end[4]), ".indiv", 7);
  if (fopen_checked(outname, "w", &outfile_indiv)) {
    goto write_main_roh_reports_ret_OPEN_FAIL;
  }
  sprintf(g_textbuf, "%%%us %%%us  PHE     NSEG       KB    KBAVG\n", plink_maxfid, plink_maxiid);
  fprintf(outfile_indiv, g_textbuf, "FID", "IID");
  g_textbuf[plink_maxfid] = ' ';
  g_textbuf[plink_maxfid + plink_maxiid + 1] = ' ';
  sample_uidx = 0;
  for (sample_idx = 0; sample_idx < sample_ct; sample_uidx++, sample_idx++) {
    next_unset_ul_unsafe_ck(sample_exclude, &sample_uidx);
    cptr = &(sample_ids[sample_uidx * max_sample_id_len]);
    cptr2 = (char*)memchr(cptr, '\t', max_sample_id_len);
    slen = (uintptr_t)(cptr2 - cptr);
    memcpy(memseta(g_textbuf, 32, plink_maxfid - slen), cptr, slen);
    slen = strlen(++cptr2);
    memcpy(memseta(wptr_iid, 32, plink_maxiid - slen), cptr2, slen);
    if (!IS_SET(pheno_nm, sample_uidx)) {
      wptr_chr = fw_strcpyn(8, missing_pheno_len, missing_pheno_str, wptr_phe);
    } else if (pheno_c) {
      wptr_chr = memseta(wptr_phe, 32, 7);
      *wptr_chr++ = '1' + IS_SET(pheno_c, sample_uidx);
    } else {
      wptr_chr = width_force(8, wptr_phe, dtoa_f_p3(pheno_d[sample_uidx], wptr_phe));
    }
    *wptr_chr++ = ' ';
    // traverse roh_list backwards, reversing the direction of [5], then
    // traverse forwards and reset [5] again to sample_idx
    cur_roh_idx = sample_to_last_roh[sample_idx];
    cur_roh_ct = 0;
    while (cur_roh_idx != ~ZEROLU) {
      cur_roh = &(roh_list[cur_roh_idx * ROH_ENTRY_INTS]);
#ifdef __LP64__
      prev_roh_idx = ((uintptr_t)cur_roh[5]) | (((uintptr_t)cur_roh[6]) << 32);
      cur_roh[5] = (uint32_t)next_roh_idx;
      cur_roh[6] = (uint32_t)(next_roh_idx >> 32);
#else
      prev_roh_idx = (uintptr_t)cur_roh[5];
      cur_roh[5] = (uint32_t)next_roh_idx;
#endif
      next_roh_idx = cur_roh_idx;
      cur_roh_idx = prev_roh_idx;
      cur_roh_ct++;
    }
    cur_roh_idx = next_roh_idx;
    kb_tot = 0;
    for (uii = 0; uii < cur_roh_ct; uii++) {
      cur_roh = &(roh_list[cur_roh_idx * ROH_ENTRY_INTS]);
      marker_uidx1 = cur_roh[0];
      marker_uidx2 = cur_roh[1];
      wptr = width_force(4, wptr_chr, chrom_name_write(chrom_info_ptr, get_variant_chrom(chrom_info_ptr, marker_uidx1), wptr_chr));
      *wptr++ = ' ';
      cptr = &(marker_ids[marker_uidx1 * max_marker_id_len]);
      slen = strlen(cptr);
      wptr = memcpya(memseta(wptr, 32, plink_maxsnp - slen), cptr, slen);
      *wptr++ = ' ';
      cptr = &(marker_ids[marker_uidx2 * max_marker_id_len]);
      slen = strlen(cptr);
      wptr = memcpya(memseta(wptr, 32, plink_maxsnp - slen), cptr, slen);
      wptr = memseta(wptr, 32, 3);
      wptr = uint32toa_w10(marker_pos[marker_uidx1], wptr);
      wptr = memseta(wptr, 32, 3);
      wptr = uint32toa_w10x(marker_pos[marker_uidx2], ' ', wptr);
      dxx = ((double)(marker_pos[marker_uidx2] + is_new_lengths - marker_pos[marker_uidx1])) / (1000.0 - EPSILON);
      kb_tot += dxx;
      wptr = width_force(10, wptr, dtoa_f_p3(dxx, wptr));
      *wptr++ = ' ';
      if (cur_roh[2] > max_roh_len) {
	max_roh_len = cur_roh[2];
      }
      wptr = uint32toa_w8x(cur_roh[2], ' ', wptr);
      dyy = (1.0 + SMALLISH_EPSILON) / ((double)((int32_t)cur_roh[2]));
      wptr = width_force(8, wptr, dtoa_f_p3(dxx * dyy, wptr));
      // next two decimals guaranteed to be length 5
      wptr = memseta(wptr, 32, 4);
      wptr = dtoa_f_p3(((double)((int32_t)cur_roh[3])) * dyy, wptr);
      wptr = memseta(wptr, 32, 4);
      wptr = dtoa_f_p3(((double)((int32_t)cur_roh[4])) * dyy, wptr);
      *wptr++ = '\n';
      if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
	goto write_main_roh_reports_ret_WRITE_FAIL;
      }
#ifdef __LP64__
      cur_roh_idx = ((uintptr_t)cur_roh[5]) | (((uintptr_t)cur_roh[6]) << 32);
#else
      cur_roh_idx = (uintptr_t)cur_roh[5];
#endif
      cur_roh[5] = sample_uidx;
    }
    if (!IS_SET(pheno_nm, sample_uidx)) {
      wptr = fw_strcpyn(4, missing_pheno_len - omp_is_numeric * 4, missing_pheno_str, wptr_phe);
    } else if (pheno_c) {
      wptr = memseta(wptr_phe, 32, 3);
      *wptr++ = '1' + IS_SET(pheno_c, sample_uidx);
    } else {
      wptr = width_force(4, wptr_phe, dtoa_g(pheno_d[sample_uidx], wptr_phe));
    }
    *wptr++ = ' ';
    wptr = uint32toa_w8x(cur_roh_ct, ' ', wptr);
    wptr = width_force(8, wptr, dtoa_g(kb_tot, wptr));
    *wptr++ = ' ';
    if (cur_roh_ct) {
      kb_tot /= (double)((int32_t)cur_roh_ct);
    }
    wptr = width_force(8, wptr, dtoa_g(kb_tot, wptr));
    if (cur_roh_ct) {
      *wptr++ = ' ';
    }
    *wptr++ = '\n';
    if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile_indiv)) {
      goto write_main_roh_reports_ret_WRITE_FAIL;
    }
  }
  if (fclose_null(&outfile)) {
    goto write_main_roh_reports_ret_WRITE_FAIL;
  }
  if (fclose_null(&outfile_indiv)) {
    goto write_main_roh_reports_ret_WRITE_FAIL;
  }
  memcpy(&(outname_end[5]), "summary", 8);
  if (fopen_checked(outname, "w", &outfile)) {
    goto write_main_roh_reports_ret_WRITE_FAIL;
  }
  sprintf(g_textbuf, " CHR %%%us           BP      AFF    UNAFF\n", plink_maxsnp);
  fprintf(outfile, g_textbuf, "SNP");
  for (chrom_fo_idx = 0; chrom_fo_idx < chrom_info_ptr->chrom_ct; chrom_fo_idx++) {
    chrom_roh_start = roh_list_chrom_starts[chrom_fo_idx];
    chrom_roh_ct = roh_list_chrom_starts[chrom_fo_idx + 1] - chrom_roh_start;
    uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
    chrom_start = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx];
    chrom_len = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx + 1] - chrom_start;
    bigstack_reset(bigstack_mark);
    if (bigstack_calloc_i(chrom_len + 1, &roh_ct_unaff_adj)) {
      goto write_main_roh_reports_ret_NOMEM;
    }
    if (pheno_c) {
      if (bigstack_calloc_i(chrom_len + 1, &roh_ct_aff_adj)) {
        goto write_main_roh_reports_ret_NOMEM;
      }
    }
    cur_roh = &(roh_list[chrom_roh_start * ROH_ENTRY_INTS]);
    for (cur_roh_idx = 0; cur_roh_idx < chrom_roh_ct; cur_roh_idx++) {
      sample_uidx = cur_roh[5];
      if ((!pheno_c) || (!IS_SET(pheno_c, sample_uidx))) {
        roh_ct_unaff_adj[cur_roh[0] - chrom_start] += 1;
        roh_ct_unaff_adj[cur_roh[1] + 1 - chrom_start] -= 1;
      } else {
        roh_ct_aff_adj[cur_roh[0] - chrom_start] += 1;
        roh_ct_aff_adj[cur_roh[1] + 1 - chrom_start] -= 1;
      }
      cur_roh = &(cur_roh[ROH_ENTRY_INTS]);
    }
    wptr_chr = width_force(4, g_textbuf, chrom_name_write(chrom_info_ptr, uii, g_textbuf));
    *wptr_chr++ = ' ';
    memset(&(wptr_chr[plink_maxsnp]), 32, 3);
    wptr_bp1 = &(wptr_chr[plink_maxsnp + 3]);
    memset(&(wptr_bp1[10]), 32, 8);
    wptr_bp1[18] = '0';
    wptr_bp1[19] = ' ';
    chrom_len += chrom_start; // now chrom_end
    cur_roh_ct = 0; // unaff ct
    uii = 0; // aff ct
    for (marker_uidx1 = chrom_start; marker_uidx1 < chrom_len; marker_uidx1++) {
      // bugfix: excluded markers can actually be at the end of a ROH.
      if (pheno_c) {
	uii += roh_ct_aff_adj[marker_uidx1 - chrom_start];
      }
      cur_roh_ct += roh_ct_unaff_adj[marker_uidx1 - chrom_start];
      if (IS_SET(marker_exclude, marker_uidx1)) {
	continue;
      }
      cptr = &(marker_ids[marker_uidx1 * max_marker_id_len]);
      slen = strlen(cptr);
      memcpy(memseta(wptr_chr, 32, plink_maxsnp - slen), cptr, slen);
      uint32toa_w10(marker_pos[marker_uidx1], wptr_bp1);
      if (!pheno_c) {
        wptr = &(wptr_bp1[20]);
      } else {
        wptr = uint32toa_w8x(uii, ' ', &(wptr_bp1[11]));
      }
      if (cur_roh_ct + uii > max_pool_size) {
        max_pool_size = cur_roh_ct + uii;
      }
      wptr = uint32toa_w8x(cur_roh_ct, '\n', wptr);
      if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
        goto write_main_roh_reports_ret_WRITE_FAIL;
      }
    }
  }
  if (fclose_null(&outfile)) {
    goto write_main_roh_reports_ret_WRITE_FAIL;
  }
  *max_pool_size_ptr = max_pool_size;
  *max_roh_len_ptr = max_roh_len;
  while (0) {
  write_main_roh_reports_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  write_main_roh_reports_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  write_main_roh_reports_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  bigstack_reset(bigstack_mark);
  fclose_cond(outfile);
  fclose_cond(outfile_indiv);
  return retval;
}

// heap of 64-bit integers with maximum on top
// root at element 1, heap_size offset by 1
// slightly more complicated variant of this in plink_cluster.c

void heapmax64_down(uint32_t cur_pos, uint32_t heap_size, uint64_t* heapmax64) {
  uint64_t cur_val = heapmax64[cur_pos];
  uint32_t child_pos = cur_pos * 2;
  uint64_t tmp_val;
  while (child_pos < heap_size) {
    tmp_val = heapmax64[child_pos];
    if ((child_pos + 1 < heap_size) && (heapmax64[child_pos + 1] > tmp_val)) {
      tmp_val = heapmax64[++child_pos];
    }
    if (cur_val >= tmp_val) {
      break;
    }
    heapmax64[cur_pos] = tmp_val;
    cur_pos = child_pos;
    child_pos *= 2;
  }
  heapmax64[cur_pos] = cur_val;
}

void heapmax64_up_then_down(uint32_t orig_pos, uint64_t* heapmax64, uint32_t heap_size) {
  uint32_t cur_pos = orig_pos;
  uint64_t cur_val = heapmax64[orig_pos];
  uint32_t parent_pos = orig_pos / 2;
  uint64_t tmp_val;
  while (parent_pos) {
    tmp_val = heapmax64[parent_pos];
    if (cur_val <= tmp_val) {
      break;
    }
    heapmax64[cur_pos] = tmp_val;
    cur_pos = parent_pos;
    parent_pos /= 2;
  }
  if (cur_pos != orig_pos) {
    heapmax64[cur_pos] = cur_val;
  }
  heapmax64_down(cur_pos, heap_size, heapmax64);
}

void cur_roh_heap_removemax(uintptr_t* roh_slot_occupied, uint64_t* cur_roh_heap, uint32_t* cur_roh_heap_top_ptr, uint32_t* cur_roh_heap_max_ptr) {
  uint32_t cur_roh_heap_top = *cur_roh_heap_top_ptr;
  uint32_t initial_heap_max = *cur_roh_heap_max_ptr;
  uint32_t new_heap_max;
  do {
    clear_bit((uint32_t)(cur_roh_heap[1]), roh_slot_occupied);
    if ((--cur_roh_heap_top) == 1) {
      new_heap_max = 0;
      break;
    }
    cur_roh_heap[1] = cur_roh_heap[cur_roh_heap_top];
    heapmax64_down(1, cur_roh_heap_top, cur_roh_heap);
    new_heap_max = (uint32_t)(cur_roh_heap[1] >> 32);
  } while (new_heap_max == initial_heap_max);
  *cur_roh_heap_top_ptr = cur_roh_heap_top;
  *cur_roh_heap_max_ptr = new_heap_max;
}

void extract_pool_info(uint32_t pool_size, uintptr_t* roh_list_idxs, uint32_t* roh_list, uint32_t* con_uidx1_ptr, uint32_t* con_uidx2_ptr, uint32_t* union_uidx1_ptr, uint32_t* union_uidx2_ptr) {
  uint32_t uii = 0;
  uint32_t* cur_roh = &(roh_list[roh_list_idxs[0] * ROH_ENTRY_INTS]);
  uint32_t con_uidx1 = cur_roh[0];
  uint32_t con_uidx2 = cur_roh[1];
  uint32_t union_uidx1 = cur_roh[0];
  uint32_t union_uidx2 = cur_roh[1];
  uint32_t cur_uidx;
  for (uii = 1; uii < pool_size; uii++) {
    cur_roh = &(roh_list[roh_list_idxs[uii] * ROH_ENTRY_INTS]);
    cur_uidx = cur_roh[0];
    if (cur_uidx > con_uidx1) {
      con_uidx1 = cur_uidx;
    } else if (cur_uidx < union_uidx1) {
      union_uidx1 = cur_uidx;
    }
    cur_uidx = cur_roh[1];
    if (cur_uidx < con_uidx2) {
      con_uidx2 = cur_uidx;
    } else if (cur_uidx > union_uidx2) {
      union_uidx2 = cur_uidx;
    }
  }
  *con_uidx1_ptr = con_uidx1;
  *con_uidx2_ptr = con_uidx2;
  *union_uidx1_ptr = union_uidx1;
  *union_uidx2_ptr = union_uidx2;
}

void initialize_roh_slot(uint32_t* cur_roh, uint32_t chrom_start, uint32_t* marker_uidx_to_cidx, uintptr_t* roh_slot, uint32_t* roh_slot_cidx_start, uint32_t* roh_slot_cidx_end, uint32_t* roh_slot_end_uidx) {
  uint32_t cidx_first = marker_uidx_to_cidx[cur_roh[0] - chrom_start];
  uint32_t cidx_last = marker_uidx_to_cidx[cur_roh[1] - chrom_start];
#ifdef __LP64__
  uint32_t cidx_first_block = cidx_first & (~63);
  uint32_t cidx_last_block = cidx_last & (~63);
  uint32_t cur_bidx = 2;
#else
  uint32_t cidx_first_block = cidx_first & (~15);
  uint32_t cidx_last_block = cidx_last & (~15);
  uint32_t cur_bidx = 1;
#endif
  uint32_t end_bidx = cur_bidx + (cidx_last_block - cidx_first_block) / BITCT2;
  uint32_t uii;

  // given an interval [a, b], when we can't just use the half-open convention
  // everywhere, "last" means b and "end" means b+1.  (yeah, it wouldn't have
  // been necessary to make this distinction at all if I hadn't used unfiltered
  // indices where they didn't belong.)
  *roh_slot_cidx_start = cidx_first;
  *roh_slot_cidx_end = cidx_last + 1;
  *roh_slot_end_uidx = cur_roh[1] + 1;
  uii = cidx_first & (BITCT2 - 1);
#ifdef __LP64__
  if (cidx_first & 32) {
    roh_slot[0] = FIVEMASK;
    roh_slot[1] = 0x1555555555555555LLU >> (2 * (31 - uii));
  } else {
    roh_slot[0] = 0x1555555555555555LLU >> (2 * (31 - uii));
    roh_slot[1] = 0;
  }
#else
  roh_slot[0] = 0x15555555 >> (2 * (15 - uii));
#endif
  fill_ulong_zero(end_bidx - cur_bidx, &(roh_slot[cur_bidx]));
  uii = cidx_last & (BITCT2 - 1);
#ifdef __LP64__
  if (cidx_last & 32) {
    // |= instead of = in case first_block and last_block are the same
    roh_slot[end_bidx - 1] |= 0x5555555555555554LLU << (2 * uii);
  } else {
    roh_slot[end_bidx - 2] |= 0x5555555555555554LLU << (2 * uii);
    roh_slot[end_bidx - 1] |= FIVEMASK;
  }
#else
  roh_slot[end_bidx - 1] |= 0x55555554 << (2 * uii);
#endif
}

void populate_roh_slot_from_lookahead_nowrap(uintptr_t* lookahead_cur_col, uint32_t read_shift, uintptr_t unfiltered_sample_ctl2, uintptr_t row_ct, uint32_t write_start, uintptr_t* write_ptr) {
  uint32_t write_shift = (write_start & (BITCT2 - 1)) * 2;
  uintptr_t row_idx;
  uintptr_t cur_word;
  write_ptr = &(write_ptr[write_start / BITCT2]);
  cur_word = *write_ptr;
  for (row_idx = 0; row_idx < row_ct; row_idx++) {
    cur_word |= ((lookahead_cur_col[row_idx * unfiltered_sample_ctl2] >> read_shift) & (3 * ONELU)) << write_shift;
    write_shift += 2;
    if (write_shift == BITCT) {
      *write_ptr++ = cur_word;
      cur_word = *write_ptr;
      write_shift = 0;
    }
  }
  if (write_shift) {
    *write_ptr = cur_word;
  }
}

void populate_roh_slots_from_lookahead_buf(uintptr_t* lookahead_buf, uintptr_t max_lookahead, uintptr_t unfiltered_sample_ctl2, uintptr_t cur_lookahead_start, uintptr_t cur_lookahead_size, uint32_t lookahead_first_cidx, uint64_t* roh_slot_map, uintptr_t roh_slot_wsize, uint32_t* roh_slot_cidx_start, uint32_t* roh_slot_cidx_end, uintptr_t* roh_slots) {
  uint32_t cur_lookahead_cidx_end = lookahead_first_cidx + cur_lookahead_size;
  uint32_t cur_lookahead_cidx_wrap = lookahead_first_cidx + max_lookahead - cur_lookahead_start;
  uint32_t sample_uidx = (uint32_t)((*roh_slot_map) >> 32);
  uintptr_t* lookahead_cur_col;
  uintptr_t* write_ptr;
  uintptr_t slot_idx;
  uint32_t read_shift;
  uint32_t cidx_start_block;
  uint32_t cidx_start;
  uint32_t cidx_end;
  do {
    lookahead_cur_col = &(lookahead_buf[sample_uidx / BITCT2]);
    read_shift = 2 * (sample_uidx & (BITCT2 - 1));
    slot_idx = (uintptr_t)((*roh_slot_map) & 0xffffffffU);
    cidx_start = roh_slot_cidx_start[slot_idx];
#ifdef __LP64__
    cidx_start_block = cidx_start & (~63);
#else
    cidx_start_block = cidx_start & (~15);
#endif
    if (lookahead_first_cidx > cidx_start) {
      cidx_start = lookahead_first_cidx;
    }
    cidx_end = roh_slot_cidx_end[slot_idx];
    if (cidx_end > cur_lookahead_cidx_end) {
      cidx_end = cur_lookahead_cidx_end;
    }
    if (cidx_end > cidx_start) {
      write_ptr = &(roh_slots[slot_idx * roh_slot_wsize]);
      if (cidx_end <= cur_lookahead_cidx_wrap) {
	populate_roh_slot_from_lookahead_nowrap(&(lookahead_cur_col[(cur_lookahead_start + cidx_start - lookahead_first_cidx) * unfiltered_sample_ctl2]), read_shift, unfiltered_sample_ctl2, cidx_end - cidx_start, cidx_start - cidx_start_block, write_ptr);
      } else {
	if (cidx_start < cur_lookahead_cidx_wrap) {
	  populate_roh_slot_from_lookahead_nowrap(&(lookahead_cur_col[(cur_lookahead_start + cidx_start - lookahead_first_cidx) * unfiltered_sample_ctl2]), read_shift, unfiltered_sample_ctl2, cur_lookahead_cidx_wrap - cidx_start, cidx_start - cidx_start_block, write_ptr);
	  cidx_start = cur_lookahead_cidx_wrap;
	}
	populate_roh_slot_from_lookahead_nowrap(&(lookahead_cur_col[(cidx_start - cur_lookahead_cidx_wrap) * unfiltered_sample_ctl2]), read_shift, unfiltered_sample_ctl2, cidx_end - cidx_start, cidx_start - cidx_start_block, write_ptr);
      }
    }
    // one-terminated list
    sample_uidx = (uint32_t)((*(++roh_slot_map)) >> 32);
  } while (sample_uidx != 0xffffffffU);
}

int32_t populate_roh_slots_from_disk(FILE* bedfile, uint64_t bed_offset, uintptr_t* rawbuf, uintptr_t* marker_exclude, uint64_t unfiltered_sample_ct4, uint32_t* marker_uidx_to_cidx, uint32_t chrom_start, uint32_t marker_uidx, uint32_t last_uidx, uint64_t* roh_slot_map, uintptr_t roh_slot_wsize, uint32_t* roh_slot_cidx_start, uint32_t* roh_slot_cidx_end, uintptr_t* roh_slots, uint32_t roh_read_slot_ct) {
  uintptr_t roh_write_slot_idx;
  uintptr_t marker_c_bidx;
  uintptr_t start_c_bidx;
  uint32_t roh_read_slot_idx;
  uint32_t marker_cidx;
  uint32_t cidx_start;
  uint32_t write_shift;
  uint32_t sample_uidx;
  if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
    return RET_READ_FAIL;
  }
  marker_uidx--;
  do {
    marker_uidx++;
    if (IS_SET(marker_exclude, marker_uidx)) {
      marker_uidx = next_unset_unsafe(marker_exclude, marker_uidx);
      if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
        return RET_READ_FAIL;
      }
    }
    if (load_raw(unfiltered_sample_ct4, bedfile, rawbuf)) {
      return RET_READ_FAIL;
    }
    marker_cidx = marker_uidx_to_cidx[marker_uidx - chrom_start];
    marker_c_bidx = marker_cidx / BITCT2;
    write_shift = 2 * (marker_cidx % BITCT2);
    for (roh_read_slot_idx = 0; roh_read_slot_idx < roh_read_slot_ct; roh_read_slot_idx++) {
      sample_uidx = (uint32_t)(roh_slot_map[roh_read_slot_idx] >> 32);
      roh_write_slot_idx = (uintptr_t)(roh_slot_map[roh_read_slot_idx] & 0xffffffffU);
      cidx_start = roh_slot_cidx_start[roh_write_slot_idx];
      if ((marker_cidx >= cidx_start) && (marker_cidx < roh_slot_cidx_end[roh_write_slot_idx])) {
#ifdef __LP64__
        start_c_bidx = 2 * (cidx_start / 64);
#else
        start_c_bidx = cidx_start / 16;
#endif
        roh_slots[roh_write_slot_idx * roh_slot_wsize + marker_c_bidx - start_c_bidx] |= ((rawbuf[sample_uidx / BITCT2] >> (2 * (sample_uidx % BITCT2))) & (3 * ONELU)) << write_shift;
      }
    }
  } while (marker_uidx != last_uidx);
  return 0;
}

static inline uint32_t is_allelic_match(double mismatch_max, uintptr_t* roh_slot_idxl, uintptr_t* roh_slot_idxs, uint32_t block_start_idxl, uint32_t block_start_idxs, uint32_t overlap_cidx_start, uint32_t overlap_cidx_end) {
#ifdef __LP64__
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  uint32_t words_left = ((overlap_cidx_end + 31) / 32) - 2 * (overlap_cidx_start / 64);
  uint32_t joint_homozyg_ct = 0;
  uint32_t joint_homozyg_mismatch_ct = 0;
  __m128i loader_l;
  __m128i loader_s;
  __m128i joint_vec;
  __m128i joint_sum1;
  __m128i mismatch_sum1;
  __m128i joint_sum2;
  __m128i mismatch_sum2;
  __univec accj;
  __univec accm;
  __m128i* vptrl;
  __m128i* vptrs;
  __m128i* vptrl_end;
  uintptr_t* roh_idxl_end;
  uintptr_t wloader_l;
  uintptr_t wloader_s;
  uintptr_t joint_word;
  // 1. set joint_vec to 01 at joint-homozygous spots, and 00 elsewhere
  // 2. popcount joint_vec to determine number of jointly homozygous sites
  // 3. then popcount (joint_vec & (roh_vslot_idxl ^ roh_vslot_idxs)) to
  //    determine number of mismatches at these sites
  roh_slot_idxl = &(roh_slot_idxl[2 * ((overlap_cidx_start - block_start_idxl) / 64)]);
  roh_slot_idxs = &(roh_slot_idxs[2 * ((overlap_cidx_start - block_start_idxs) / 64)]);
  roh_idxl_end = &(roh_slot_idxl[words_left]);
  if (words_left >= 11) {
    vptrl = (__m128i*)roh_slot_idxl;
    vptrs = (__m128i*)roh_slot_idxs;
    if (words_left % 12 == 11) {
      words_left++;
    }
    words_left -= words_left % 12;
    while (words_left >= 120) {
      words_left -= 120;
      vptrl_end = &(vptrl[60]);
    is_allelic_match_main_loop:
      accj.vi = _mm_setzero_si128();
      accm.vi = _mm_setzero_si128();
      do {
	loader_l = *vptrl++;
	loader_s = *vptrs++;
	joint_sum1 = _mm_andnot_si128(_mm_or_si128(_mm_xor_si128(loader_l, _mm_srli_epi64(loader_l, 1)), _mm_xor_si128(loader_s, _mm_srli_epi64(loader_s, 1))), m1);
	mismatch_sum1 = _mm_and_si128(joint_sum1, _mm_xor_si128(loader_l, loader_s));

	loader_l = *vptrl++;
	loader_s = *vptrs++;
	joint_vec = _mm_andnot_si128(_mm_or_si128(_mm_xor_si128(loader_l, _mm_srli_epi64(loader_l, 1)), _mm_xor_si128(loader_s, _mm_srli_epi64(loader_s, 1))), m1);
	joint_sum1 = _mm_add_epi64(joint_sum1, joint_vec);
	joint_vec = _mm_and_si128(joint_vec, _mm_xor_si128(loader_l, loader_s));
	mismatch_sum1 = _mm_add_epi64(mismatch_sum1, joint_vec);

	loader_l = *vptrl++;
	loader_s = *vptrs++;
	joint_vec = _mm_andnot_si128(_mm_or_si128(_mm_xor_si128(loader_l, _mm_srli_epi64(loader_l, 1)), _mm_xor_si128(loader_s, _mm_srli_epi64(loader_s, 1))), m1);
	joint_sum1 = _mm_add_epi64(joint_sum1, joint_vec);
	joint_vec = _mm_and_si128(joint_vec, _mm_xor_si128(loader_l, loader_s));
	mismatch_sum1 = _mm_add_epi64(mismatch_sum1, joint_vec);

	joint_sum1 = _mm_add_epi64(_mm_and_si128(joint_sum1, m2), _mm_and_si128(_mm_srli_epi64(joint_sum1, 2), m2));
	mismatch_sum1 = _mm_add_epi64(_mm_and_si128(mismatch_sum1, m2), _mm_and_si128(_mm_srli_epi64(mismatch_sum1, 2), m2));

	loader_l = *vptrl++;
	loader_s = *vptrs++;
	joint_sum2 = _mm_andnot_si128(_mm_or_si128(_mm_xor_si128(loader_l, _mm_srli_epi64(loader_l, 1)), _mm_xor_si128(loader_s, _mm_srli_epi64(loader_s, 1))), m1);
	mismatch_sum2 = _mm_and_si128(joint_sum2, _mm_xor_si128(loader_l, loader_s));

	loader_l = *vptrl++;
	loader_s = *vptrs++;
	joint_vec = _mm_andnot_si128(_mm_or_si128(_mm_xor_si128(loader_l, _mm_srli_epi64(loader_l, 1)), _mm_xor_si128(loader_s, _mm_srli_epi64(loader_s, 1))), m1);
	joint_sum2 = _mm_add_epi64(joint_sum2, joint_vec);
	joint_vec = _mm_and_si128(joint_vec, _mm_xor_si128(loader_l, loader_s));
	mismatch_sum2 = _mm_add_epi64(mismatch_sum2, joint_vec);

	loader_l = *vptrl++;
	loader_s = *vptrs++;
	joint_vec = _mm_andnot_si128(_mm_or_si128(_mm_xor_si128(loader_l, _mm_srli_epi64(loader_l, 1)), _mm_xor_si128(loader_s, _mm_srli_epi64(loader_s, 1))), m1);
	joint_sum2 = _mm_add_epi64(joint_sum2, joint_vec);
	joint_vec = _mm_and_si128(joint_vec, _mm_xor_si128(loader_l, loader_s));
	mismatch_sum2 = _mm_add_epi64(mismatch_sum2, joint_vec);

	joint_sum1 = _mm_add_epi64(joint_sum1, _mm_add_epi64(_mm_and_si128(joint_sum2, m2), _mm_and_si128(_mm_srli_epi64(joint_sum2, 2), m2)));
	mismatch_sum1 = _mm_add_epi64(mismatch_sum1, _mm_add_epi64(_mm_and_si128(mismatch_sum2, m2), _mm_and_si128(_mm_srli_epi64(mismatch_sum2, 2), m2)));

	accj.vi = _mm_add_epi64(accj.vi, _mm_add_epi64(_mm_and_si128(joint_sum1, m4), _mm_and_si128(_mm_srli_epi64(joint_sum1, 4), m4)));
	accm.vi = _mm_add_epi64(accm.vi, _mm_add_epi64(_mm_and_si128(mismatch_sum1, m4), _mm_and_si128(_mm_srli_epi64(mismatch_sum1, 4), m4)));
      } while (vptrl < vptrl_end);
      accj.vi = _mm_add_epi64(_mm_and_si128(accj.vi, m8), _mm_and_si128(_mm_srli_epi64(accj.vi, 8), m8));
      accm.vi = _mm_add_epi64(_mm_and_si128(accm.vi, m8), _mm_and_si128(_mm_srli_epi64(accm.vi, 8), m8));
      joint_homozyg_ct += ((accj.u8[0] + accj.u8[1]) * 0x1000100010001LLU) >> 48;
      joint_homozyg_mismatch_ct += ((accm.u8[0] + accm.u8[1]) * 0x1000100010001LLU) >> 48;
    }
    if (words_left) {
      vptrl_end = &(vptrl[words_left / 2]);
      words_left = 0;
      goto is_allelic_match_main_loop;
    }
    roh_slot_idxl = (uintptr_t*)vptrl;
    roh_slot_idxs = (uintptr_t*)vptrs;
  }
#else
  uint32_t joint_homozyg_ct = 0;
  uint32_t joint_homozyg_mismatch_ct = 0;
  uint32_t words_left = ((overlap_cidx_end + 15) / 16) - (overlap_cidx_start / 16);
  uintptr_t* roh_idxl_end;
  uintptr_t wloader_l;
  uintptr_t wloader_s;
  uintptr_t joint_word;
  uintptr_t joint_sum1;
  uintptr_t mismatch_sum1;
  uintptr_t joint_sum2;
  uintptr_t mismatch_sum2;
  uintptr_t joint_acc;
  uintptr_t mismatch_acc;
  roh_slot_idxl = &(roh_slot_idxl[(overlap_cidx_start - block_start_idxl) / 16]);
  roh_slot_idxs = &(roh_slot_idxs[(overlap_cidx_start - block_start_idxs) / 16]);
  roh_idxl_end = &(roh_slot_idxl[words_left]);
  words_left -= words_left % 12;
  while (words_left) {
    wloader_l = *roh_slot_idxl++;
    wloader_s = *roh_slot_idxs++;
    joint_sum1 = (~((wloader_l ^ (wloader_l >> 1)) | (wloader_s ^ (wloader_s >> 1)))) & FIVEMASK;
    mismatch_sum1 = joint_sum1 & (wloader_l ^ wloader_s);

    wloader_l = *roh_slot_idxl++;
    wloader_s = *roh_slot_idxs++;
    joint_word = (~((wloader_l ^ (wloader_l >> 1)) | (wloader_s ^ (wloader_s >> 1)))) & FIVEMASK;
    joint_sum1 += joint_word;
    mismatch_sum1 += joint_word & (wloader_l ^ wloader_s);

    wloader_l = *roh_slot_idxl++;
    wloader_s = *roh_slot_idxs++;
    joint_word = (~((wloader_l ^ (wloader_l >> 1)) | (wloader_s ^ (wloader_s >> 1)))) & FIVEMASK;
    joint_sum1 += joint_word;
    mismatch_sum1 += joint_word & (wloader_l ^ wloader_s);

    joint_sum1 = (joint_sum1 & 0x33333333) + ((joint_sum1 >> 2) & 0x33333333);
    mismatch_sum1 = (mismatch_sum1 & 0x33333333) + ((mismatch_sum1 >> 2) & 0x33333333);

    wloader_l = *roh_slot_idxl++;
    wloader_s = *roh_slot_idxs++;
    joint_sum2 = (~((wloader_l ^ (wloader_l >> 1)) | (wloader_s ^ (wloader_s >> 1)))) & FIVEMASK;
    mismatch_sum2 = joint_sum2 & (wloader_l ^ wloader_s);

    wloader_l = *roh_slot_idxl++;
    wloader_s = *roh_slot_idxs++;
    joint_word = (~((wloader_l ^ (wloader_l >> 1)) | (wloader_s ^ (wloader_s >> 1)))) & FIVEMASK;
    joint_sum2 += joint_word;
    mismatch_sum2 += joint_word & (wloader_l ^ wloader_s);

    wloader_l = *roh_slot_idxl++;
    wloader_s = *roh_slot_idxs++;
    joint_word = (~((wloader_l ^ (wloader_l >> 1)) | (wloader_s ^ (wloader_s >> 1)))) & FIVEMASK;
    joint_sum2 += joint_word;
    mismatch_sum2 += joint_word & (wloader_l ^ wloader_s);

    joint_sum1 += (joint_sum2 & 0x33333333) + ((joint_sum2 >> 2) & 0x33333333);
    mismatch_sum1 += (mismatch_sum2 & 0x33333333) + ((mismatch_sum2 >> 2) & 0x33333333);
    joint_acc = (joint_sum1 & 0x0f0f0f0f) + ((joint_sum1 >> 4) & 0x0f0f0f0f);
    mismatch_acc = (mismatch_sum1 & 0x0f0f0f0f) + ((mismatch_sum1 >> 4) & 0x0f0f0f0f);

    wloader_l = *roh_slot_idxl++;
    wloader_s = *roh_slot_idxs++;
    joint_sum1 = (~((wloader_l ^ (wloader_l >> 1)) | (wloader_s ^ (wloader_s >> 1)))) & FIVEMASK;
    mismatch_sum1 = joint_sum1 & (wloader_l ^ wloader_s);

    wloader_l = *roh_slot_idxl++;
    wloader_s = *roh_slot_idxs++;
    joint_word = (~((wloader_l ^ (wloader_l >> 1)) | (wloader_s ^ (wloader_s >> 1)))) & FIVEMASK;
    joint_sum1 += joint_word;
    mismatch_sum1 += joint_word & (wloader_l ^ wloader_s);

    wloader_l = *roh_slot_idxl++;
    wloader_s = *roh_slot_idxs++;
    joint_word = (~((wloader_l ^ (wloader_l >> 1)) | (wloader_s ^ (wloader_s >> 1)))) & FIVEMASK;
    joint_sum1 += joint_word;
    mismatch_sum1 += joint_word & (wloader_l ^ wloader_s);

    joint_sum1 = (joint_sum1 & 0x33333333) + ((joint_sum1 >> 2) & 0x33333333);
    mismatch_sum1 = (mismatch_sum1 & 0x33333333) + ((mismatch_sum1 >> 2) & 0x33333333);

    wloader_l = *roh_slot_idxl++;
    wloader_s = *roh_slot_idxs++;
    joint_sum2 = (~((wloader_l ^ (wloader_l >> 1)) | (wloader_s ^ (wloader_s >> 1)))) & FIVEMASK;
    mismatch_sum2 = joint_sum2 & (wloader_l ^ wloader_s);

    wloader_l = *roh_slot_idxl++;
    wloader_s = *roh_slot_idxs++;
    joint_word = (~((wloader_l ^ (wloader_l >> 1)) | (wloader_s ^ (wloader_s >> 1)))) & FIVEMASK;
    joint_sum2 += joint_word;
    mismatch_sum2 += joint_word & (wloader_l ^ wloader_s);

    wloader_l = *roh_slot_idxl++;
    wloader_s = *roh_slot_idxs++;
    joint_word = (~((wloader_l ^ (wloader_l >> 1)) | (wloader_s ^ (wloader_s >> 1)))) & FIVEMASK;
    joint_sum2 += joint_word;
    mismatch_sum2 += joint_word & (wloader_l ^ wloader_s);

    joint_sum1 += (joint_sum2 & 0x33333333) + ((joint_sum2 >> 2) & 0x33333333);
    mismatch_sum1 += (mismatch_sum2 & 0x33333333) + ((mismatch_sum2 >> 2) & 0x33333333);
    joint_acc += (joint_sum1 & 0x0f0f0f0f) + ((joint_sum1 >> 4) & 0x0f0f0f0f);
    mismatch_acc += (mismatch_sum1 & 0x0f0f0f0f) + ((mismatch_sum1 >> 4) & 0x0f0f0f0f);
    joint_homozyg_ct += (joint_acc * 0x01010101) >> 24;
    joint_homozyg_mismatch_ct += (mismatch_sum1 * 0x01010101) >> 24;
    words_left -= 12;
  }
#endif
  while (roh_slot_idxl < roh_idxl_end) {
    wloader_l = *roh_slot_idxl++;
    wloader_s = *roh_slot_idxs++;
    joint_word = (~((wloader_l ^ (wloader_l >> 1)) | (wloader_s ^ (wloader_s >> 1)))) & FIVEMASK;
    joint_homozyg_ct += popcount2_long(joint_word);
    joint_homozyg_mismatch_ct += popcount2_long(joint_word & (wloader_l ^ wloader_s));
  }
  return (((double)((int32_t)joint_homozyg_mismatch_ct)) <= mismatch_max * ((int32_t)joint_homozyg_ct));
}

void compute_allelic_match_matrix(double mismatch_max, uintptr_t roh_slot_wsize, uint32_t pool_size, uintptr_t* roh_slots, uintptr_t* roh_slot_occupied, uintptr_t* roh_slot_uncached, uint32_t* roh_slot_cidx_start, uint32_t* roh_slot_cidx_end, uint64_t* roh_slot_map, uint32_t overlap_cidx_start, uint32_t overlap_cidx_end, uint32_t* allelic_match_cts, uintptr_t* allelic_match_matrix) {
  // consensus_match in effect iff roh_slot_uncached is nullptr
  // may want to make this multithreaded in the future
  uint32_t cidx_end_idxl = 0;
  uint32_t skip_cached = 0;
  uintptr_t* roh_slot_idxl;
  uint32_t map_idxl;
  uint32_t map_idxs;
  uint32_t map_first_unset;
  uint32_t cur_limit;
  uint32_t slot_idxl;
  uint32_t slot_idxs;
  uintptr_t tri_offset_idxl;
  uintptr_t tri_coord;
  uint32_t incr_idxl;
  uint32_t cidx_start_idxl;
  uint32_t cidx_start_idxs;
  uint32_t block_start_idxl;
  uint32_t block_start_idxs;
  uint32_t uii;
  fill_uint_zero(pool_size, allelic_match_cts);
  if (roh_slot_uncached) {
    // count cached results
    map_first_unset = next_unset(roh_slot_uncached, 0, pool_size);
    cur_limit = 1;
    for (map_idxl = map_first_unset + 1; map_idxl < pool_size; map_idxl++) {
      if (IS_SET(roh_slot_uncached, map_idxl)) {
	continue;
      }
      slot_idxl = (uint32_t)(roh_slot_map[map_idxl]);
      tri_offset_idxl = (((uintptr_t)slot_idxl) * (slot_idxl - 1)) / 2;
      incr_idxl = 0;
      map_idxs = map_first_unset;
      uii = 0;
      while (1) {
	slot_idxs = (uint32_t)(roh_slot_map[map_idxs]);
	if (slot_idxs < slot_idxl) {
	  tri_coord = tri_offset_idxl + slot_idxs;
	} else {
	  tri_coord = tri_coord_no_diag(slot_idxl, slot_idxs);
	}
        if (IS_SET(allelic_match_matrix, tri_coord)) {
	  allelic_match_cts[map_idxs] += 1;
	  incr_idxl++;
	}
	if (++uii == cur_limit) {
	  break;
	}
	do {
	  map_idxs++;
	} while (IS_SET(roh_slot_uncached, map_idxs));
      }
      allelic_match_cts[map_idxl] += incr_idxl;
      cur_limit++;
    }
  }
  for (map_idxl = 1; map_idxl < pool_size; map_idxl++) {
    slot_idxl = (uint32_t)(roh_slot_map[map_idxl]);
    tri_offset_idxl = (((uintptr_t)slot_idxl) * (slot_idxl - 1)) / 2;
    incr_idxl = 0;
    roh_slot_idxl = &(roh_slots[slot_idxl * roh_slot_wsize]);
    cidx_start_idxl = roh_slot_cidx_start[slot_idxl];
#ifdef __LP64__
    block_start_idxl = cidx_start_idxl & (~63);
#else
    block_start_idxl = cidx_start_idxl & (~15);
#endif
    if (roh_slot_uncached) {
      cidx_end_idxl = roh_slot_cidx_end[slot_idxl];
      skip_cached = 1 ^ IS_SET(roh_slot_uncached, map_idxl);
    }
    for (map_idxs = 0; map_idxs < map_idxl; map_idxs++) {
      if (skip_cached && (!IS_SET(roh_slot_uncached, map_idxs))) {
	map_idxs = next_set(roh_slot_uncached, map_idxs, map_idxl);
	if (map_idxs == map_idxl) {
	  break;
	}
      }
      slot_idxs = (uint32_t)(roh_slot_map[map_idxs]);
      cidx_start_idxs = roh_slot_cidx_start[slot_idxs];
#ifdef __LP64__
      block_start_idxs = cidx_start_idxs & (~63);
#else
      block_start_idxs = cidx_start_idxs & (~15);
#endif
      if (roh_slot_uncached) {
        overlap_cidx_start = cidx_start_idxs;
        if (overlap_cidx_start < cidx_start_idxl) {
	  overlap_cidx_start = cidx_start_idxl;
	}
	overlap_cidx_end = roh_slot_cidx_end[slot_idxs];
	if (overlap_cidx_end > cidx_end_idxl) {
          overlap_cidx_end = cidx_end_idxl;
	}
      }
      if (is_allelic_match(mismatch_max, roh_slot_idxl, &(roh_slots[slot_idxs * roh_slot_wsize]), block_start_idxl, block_start_idxs, overlap_cidx_start, overlap_cidx_end)) {
	if (slot_idxs < slot_idxl) {
	  tri_coord = tri_offset_idxl + slot_idxs;
	} else {
	  tri_coord = tri_coord_no_diag(slot_idxl, slot_idxs);
	}
        SET_BIT(tri_coord, allelic_match_matrix);
        allelic_match_cts[map_idxs] += 1;
	incr_idxl++;
      }
    }
    allelic_match_cts[map_idxl] += incr_idxl;
  }
}

void assign_allelic_match_groups(uint32_t pool_size, uint32_t* allelic_match_cts, uintptr_t* allelic_match_matrix, uint64_t* roh_slot_map, uintptr_t* cur_pool) {
  uintptr_t group_idx = 1;
  uint32_t nsim_nz_ct = 0;
  uint32_t max_nsim_pidx = 0;
  uintptr_t ulii;
  uintptr_t main_slot_idx;
  uintptr_t tri_offset;
  uintptr_t slot_idx2;
  uintptr_t tri_coord;
  uint32_t pool_idx;
  uint32_t max_nsim;
  uint32_t uii;
  for (pool_idx = 0; pool_idx < pool_size; pool_idx++) {
    ulii = allelic_match_cts[pool_idx];
    if (ulii) {
      nsim_nz_ct++;
    }
#ifdef __LP64__
    cur_pool[pool_idx] = ulii << 32;
#else
    cur_pool[2 * pool_idx + 1] = ulii;
#endif
  }
  while (nsim_nz_ct) {
    max_nsim = 0;
    for (pool_idx = 0; pool_idx < pool_size; pool_idx++) {
      uii = allelic_match_cts[pool_idx];
      if ((uii != 0xffffffffU) && (uii > max_nsim)) {
        max_nsim = uii;
        max_nsim_pidx = pool_idx;
      }
    }
    nsim_nz_ct--;
    main_slot_idx = (uintptr_t)((uint32_t)roh_slot_map[max_nsim_pidx]);
    tri_offset = (main_slot_idx * (main_slot_idx - 1)) / 2;
    allelic_match_cts[max_nsim_pidx] = 0xffffffffU;
    for (pool_idx = 0; pool_idx < pool_size; pool_idx++) {
      if (max_nsim_pidx == pool_idx) {
	continue;
      }
      slot_idx2 = (uintptr_t)((uint32_t)roh_slot_map[pool_idx]);
      if (slot_idx2 < main_slot_idx) {
	tri_coord = tri_offset + slot_idx2;
      } else {
	tri_coord = tri_coord_no_diag(main_slot_idx, slot_idx2);
      }
      if (IS_SET(allelic_match_matrix, tri_coord)) {
	if (allelic_match_cts[pool_idx] != 0xffffffffU) {
	  nsim_nz_ct--;
	  allelic_match_cts[pool_idx] = 0xffffffffU;
	}
#ifdef __LP64__
        cur_pool[pool_idx] = (cur_pool[pool_idx] & 0xffffffff00000000LLU) | group_idx;
#else
        cur_pool[2 * pool_idx] = group_idx;
#endif
      }
    }
#ifdef __LP64__
    cur_pool[max_nsim_pidx] |= 0x80000000LLU | (group_idx++);
#else
    cur_pool[2 * max_nsim_pidx] = 0x80000000U | (group_idx++);
#endif
  }
  for (pool_idx = 0; pool_idx < pool_size; pool_idx++) {
    if (allelic_match_cts[pool_idx] != 0xffffffffU) {
#ifdef __LP64__
      cur_pool[pool_idx] |= 0x80000000LLU | (group_idx++);
#else
      cur_pool[2 * pool_idx] = 0x80000000U | (group_idx++);
#endif
    }
  }
}

char* roh_pool_write_middle(char* wptr, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, uint32_t* marker_pos, uint32_t is_new_lengths, uint32_t marker_uidx1, uint32_t marker_uidx2) {
  *wptr++ = ' ';
  wptr = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx1 * max_marker_id_len]), wptr);
  *wptr++ = ' ';
  wptr = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx2 * max_marker_id_len]), wptr);
  wptr = memseta(wptr, 32, 5);
  wptr = uint32toa_w10(marker_pos[marker_uidx1], wptr);
  wptr = memseta(wptr, 32, 5);
  wptr = uint32toa_w10x(marker_pos[marker_uidx2], ' ', wptr);
  wptr = dtoa_g_wxp8x(((double)(marker_pos[marker_uidx2] + is_new_lengths - marker_pos[marker_uidx1])) / 1000.0, 8, ' ', wptr);
  return wptr;
}

int32_t roh_pool(Homozyg_info* hp, FILE* bedfile, uint64_t bed_offset, char* outname, char* outname_end, uintptr_t* rawbuf, uintptr_t* marker_exclude, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, uintptr_t* marker_reverse, Chrom_info* chrom_info_ptr, uint32_t* marker_pos, uintptr_t sample_ct, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, char* sample_ids, uint32_t plink_maxfid, uint32_t plink_maxiid, uintptr_t max_sample_id_len, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, char* missing_pheno_str, uint32_t missing_pheno_len, uint32_t is_new_lengths, uintptr_t roh_ct, uint32_t* roh_list, uintptr_t* roh_list_chrom_starts, uint32_t max_pool_size, uint32_t max_roh_len) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  uint64_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t unfiltered_sample_ctl2 = QUATERCT_TO_WORDCT(unfiltered_sample_ct);
  double mismatch_max = 1 - (hp->overlap_min * (1 - EPSILON)); // fuzz
  uint32_t is_consensus_match = hp->modifier & HOMOZYG_CONSENSUS_MATCH;
  uint32_t is_verbose = hp->modifier & HOMOZYG_GROUP_VERBOSE;
  uint32_t max_pool_sizel = BITCT_TO_WORDCT(max_pool_size);
  uint32_t pool_size_min = hp->pool_size_min;
  uint32_t pool_size_ct = max_pool_size + 1 - pool_size_min;
  uint32_t marker_uidx2 = 0;
  uint32_t fresh_meat = 0;
  uintptr_t pool_list_size = 0;
  uint32_t pool_ct = 0;
  uint32_t onechar_max = (chrom_info_ptr->max_code > 9)? 9 : chrom_info_ptr->max_code;
  uintptr_t* roh_slot_uncached = nullptr;
  uint64_t* verbose_group_sort_buf = nullptr;
  uint32_t* verbose_uidx_bounds = nullptr;
  uint32_t* verbose_sample_uidx = nullptr;
  char* writebuf = nullptr;
  int32_t retval = 0;
  char* allele_strs[4];

  // Circular lookahead buffer.  This is one of the few analyses which is
  // noticeably inconvenienced by the main data file being SNP-major instead of
  // sample-major; a pool that ends at SNP n could contain a few runs of
  // homozygosity which don't end until SNP n+C for large C, and we'd like to
  // cache just those genotypes.  In most cases, we can simply load everyone's
  // genotype data at the next C SNPs and avoid the need to keep going forwards
  // and backwards in the file.  However, with very large datasets, there may
  // not be sufficient memory to do this; in that case, we load as many full
  // SNPs as we can, and then switch to load-and-forget-and-reload-later for
  // the final SNPs.
  uintptr_t* lookahead_buf;
  uintptr_t* lookahead_row;

  uintptr_t* pool_size_first_plidx;
  uint32_t* marker_uidx_to_cidx;
  uint32_t* chrom_fo_idx_to_pidx; // decreasing order
  uintptr_t* roh_slots;
  uintptr_t* roh_slot_occupied; // bitfield marking roh_slots occupancy

  // round this value down to nearest multiple of 16 (in 32-bit build) or 64
  // (in 64-bit build) to determine where the roh_slot actually begins.  (This
  // enables use of vector popcount.)
  uint32_t* roh_slot_cidx_start;
  uint32_t* roh_slot_cidx_end;

  uint32_t* sample_uidx_sort_buf;
  uint64_t* cur_roh_heap; // high 32 bits = marker_uidx, low 32 bits = slot_idx
  uintptr_t* pool_list;
  uint32_t* cur_roh;
  uintptr_t* cur_pool;
  uint64_t* roh_slot_map; // high 32 bits = sample_uidx, low 32 bits = slot_idx
  uint32_t* roh_slot_end_uidx; // tracks when to flush a roh_slot
  uint32_t* allelic_match_cts;
  uintptr_t* allelic_match_matrix; // pairwise match matrix, potentially huge
  uint32_t* uiptr;
  char* wptr_start;
  char* wptr;
  char* cptr;
  char* cptr2;
  uint64_t ullii;
  uintptr_t cur_lookahead_start;
  uintptr_t cur_lookahead_size;
  uintptr_t max_lookahead;
  uintptr_t old_pool_list_size;
  uintptr_t pool_list_idx;
  uintptr_t roh_slot_wsize;
  uintptr_t max_pool_list_size;
  uintptr_t chrom_roh_start;
  uintptr_t roh_idx;
  uintptr_t ulii;
  uint32_t lookahead_first_cidx;
  uint32_t lookahead_end_uidx;
  uint32_t pool_size;
  uint32_t chrom_fo_idx;
  uint32_t cur_roh_heap_top;
  uint32_t marker_uidx1;
  uint32_t marker_cidx;
  uint32_t slot_idx1;
  uint32_t slot_idx2;
  uint32_t group_slot_end;
  uint32_t chrom_start;
  uint32_t chrom_len;
  uint32_t con_uidx1;
  uint32_t con_uidx2;
  uint32_t union_uidx1;
  uint32_t union_uidx2;
  uint32_t sample_uidx1;
  uint32_t sample_uidx2;
  uint32_t case_ct;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  if (bigstack_alloc_ui(chrom_info_ptr->chrom_ct + 1, &chrom_fo_idx_to_pidx)) {
    goto roh_pool_ret_NOMEM;
  }
  uii = 0; // max chrom len
  for (chrom_fo_idx = 0; chrom_fo_idx < chrom_info_ptr->chrom_ct; chrom_fo_idx++) {
    if (roh_list_chrom_starts[chrom_fo_idx] == roh_list_chrom_starts[chrom_fo_idx + 1]) {
      continue;
    }
    chrom_len = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx + 1] - chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx];
    if (chrom_len > uii) {
      uii = chrom_len;
    }
  }
#ifdef __LP64__
  // want each roh_slots space to be 16-byte aligned, to enable SSE2
  // max_roh_len = 1 -> 1 vec
  // max_roh_len in {2..65} -> 2 vecs
  // max_roh_len in {66..129} -> 3 vecs, etc.
  roh_slot_wsize = 2 * ((max_roh_len + 126) / 64);
#else
  // max_roh_len = 1 -> 1 word
  // max_roh_len in {2..17} -> 2 words, etc.
  roh_slot_wsize = (max_roh_len + 30) / 16;
#endif
  if (bigstack_alloc_ul(pool_size_ct, &pool_size_first_plidx) ||
      bigstack_alloc_ui(uii, &marker_uidx_to_cidx) ||
      bigstack_alloc_ul(max_pool_size * roh_slot_wsize, &roh_slots) ||
      bigstack_calloc_ul(max_pool_sizel, &roh_slot_occupied) ||
      bigstack_alloc_ull(max_pool_size + 1, &roh_slot_map) ||
      bigstack_alloc_ui(max_pool_size, &roh_slot_cidx_start) ||
      bigstack_alloc_ui(max_pool_size, &roh_slot_cidx_end) ||
      bigstack_alloc_ui(max_pool_size, &roh_slot_end_uidx)) {
    goto roh_pool_ret_NOMEM;
  }
  if (!is_consensus_match) {
    if (bigstack_alloc_ul(max_pool_sizel, &roh_slot_uncached)) {
      goto roh_pool_ret_NOMEM;
    }
  }
  if (is_verbose) {
    if (bigstack_alloc_ull(max_pool_size, &verbose_group_sort_buf) ||
        bigstack_alloc_ui(max_pool_size * 2, &verbose_uidx_bounds) ||
        bigstack_alloc_ui(max_pool_size, &verbose_sample_uidx) ||
        bigstack_alloc_c(2 * max_marker_allele_len + 5, &writebuf)) {
      goto roh_pool_ret_NOMEM;
    }
  }
  if (bigstack_alloc_ui(max_pool_size, &allelic_match_cts) ||
      bigstack_alloc_ul((((uintptr_t)max_pool_size) * (max_pool_size - 1)) / 2, &allelic_match_matrix)) {
    goto roh_pool_ret_NOMEM;
  }
  // roh_slot_map / roh_slot_cidx_start... not used at the same time as
  // cur_roh_heap / sample_uidx_sort_buf
  cur_roh_heap = &(roh_slot_map[-1]);
  sample_uidx_sort_buf = roh_slot_cidx_start;

  fill_ulong_one(pool_size_ct, pool_size_first_plidx);

  pool_list = (uintptr_t*)g_bigstack_base;
  max_pool_list_size = bigstack_left() / sizeof(intptr_t);
  // Since our ROH are sorted by *last* SNP, it's easiest to scan for pools
  // from back to front if we wish to painlessly produce sorted lists.
  chrom_fo_idx = chrom_info_ptr->chrom_ct;

  do {
    chrom_fo_idx_to_pidx[chrom_fo_idx--] = pool_ct;
    chrom_roh_start = roh_list_chrom_starts[chrom_fo_idx];
    roh_idx = roh_list_chrom_starts[chrom_fo_idx + 1];
    if (chrom_roh_start == roh_idx) {
      continue;
    }
    cur_roh_heap_top = 1;

    // last marker in next ROH; marker_uidx2 is maximum heap element
    marker_uidx1 = roh_list[1 + (roh_idx - 1) * ROH_ENTRY_INTS];
    do {
      if ((marker_uidx2 <= marker_uidx1) && (roh_idx != chrom_roh_start)) {
	roh_idx--;
	cur_roh = &(roh_list[roh_idx * ROH_ENTRY_INTS]);
        uii = cur_roh[0];
	// check if this ROH doesn't intersect anything
	if ((cur_roh_heap_top > 1) || ((roh_idx != chrom_roh_start) && (cur_roh[1 - ROH_ENTRY_INTS] >= uii))) {
	  slot_idx1 = next_unset_unsafe(roh_slot_occupied, 0);
	  SET_BIT(slot_idx1, roh_slot_occupied);
	  // use roh_slots[0..(max_pool_size - 1)] to store references to
	  // active ROH here
	  roh_slots[slot_idx1] = roh_idx;
          cur_roh_heap[cur_roh_heap_top++] = (((uint64_t)uii) << 32) | ((uint64_t)slot_idx1);
	  heapmax64_up_then_down(cur_roh_heap_top - 1, cur_roh_heap, cur_roh_heap_top);
	  marker_uidx2 = (uint32_t)(cur_roh_heap[1] >> 32);
	  fresh_meat = 1;
	}
	if (roh_idx != chrom_roh_start) {
	  marker_uidx1 = cur_roh[1 - ROH_ENTRY_INTS];
	} else {
	  marker_uidx1 = 0;
	}
      } else {
	if (fresh_meat) {
	  // At every SNP where a ROH begins, we check if another ROH ended
	  // between that SNP (included) and the next SNP where a ROH begins
          // (excluded).  (This is tracked by the fresh_meat flag.)  If, and
	  // only if, that is the case, we are at the beginning of the
	  // consensus region for a maximal pool.
	  pool_size = popcount_longs(roh_slot_occupied, max_pool_sizel);
	  if (pool_size >= pool_size_min) {
	    // pool encoding:
	    // [0]: pool_list index of next pool of same size (~ZEROLU if none)
	    // 64-bit:
            //   [1]: pool size (P) in low 32 bits, 1-based pool idx in high
	    //   [2-(P+1)]: roh indexes, sorted by increasing sample_idx
	    //   [(P+2)-(2P+1)]: allelic-match group assignment (32-bit) with
	    //                   31st bit set if reference member; NSIM count
	    //                   in high 32 bits
	    //   [2P+2]: consensus NSNP in low 32, union NSNP in high 32
	    // 32-bit:
            //   [1]: P
	    //   [2]: 1-based pool idx
	    //   [3-(P+2)]: roh indexes
	    //   [(P+3)-(3P+2)]: allelic-match group assignment (31st bit set
	    //                   if reference), followed by NSIM count,
	    //                   interleaved
	    //   [3P+3]: consensus NSNP
	    //   [3P+4]: union NSNP
	    old_pool_list_size = pool_list_size;
#ifdef __LP64__
	    pool_list_size += 2 * pool_size + 3;
#else
            pool_list_size += 3 * pool_size + 5;
#endif
	    if (pool_list_size > max_pool_list_size) {
	      goto roh_pool_ret_NOMEM;
	    }
	    cur_pool = &(pool_list[old_pool_list_size]);
            *cur_pool++ = pool_size_first_plidx[pool_size - pool_size_min];
            pool_size_first_plidx[pool_size - pool_size_min] = old_pool_list_size;
	    *cur_pool++ = pool_size;
#ifndef __LP64__
	    *cur_pool++ = 0;
#endif
	    uiptr = sample_uidx_sort_buf;
	    slot_idx1 = 0;
	    for (uii = 0; uii < pool_size; slot_idx1++, uii++) {
	      slot_idx1 = next_set_unsafe(roh_slot_occupied, slot_idx1);
	      pool_list_idx = roh_slots[slot_idx1]; // actually a ROH idx
	      *uiptr++ = roh_list[pool_list_idx * ROH_ENTRY_INTS + 5]; // sample_uidx
              *uiptr++ = (uint32_t)pool_list_idx;
#ifdef __LP64__
	      *uiptr++ = (uint32_t)(pool_list_idx >> 32);
#endif
	    }
	    // sort in increasing sample_uidx order, for reproducible results
            qsort(sample_uidx_sort_buf, pool_size, 4 + sizeof(intptr_t), intcmp);
	    for (uii = 0; uii < pool_size; uii++) {
#ifdef __LP64__
              *cur_pool++ = ((uintptr_t)sample_uidx_sort_buf[3 * uii + 1]) | (((uintptr_t)sample_uidx_sort_buf[3 * uii + 2]) << 32);
#else
	      *cur_pool++ = sample_uidx_sort_buf[2 * uii + 1];
#endif
	    }
	    // leave a backward pointer
	    pool_list[pool_list_size - 1] = old_pool_list_size;
	    pool_ct++;
	  }
	  fresh_meat = 0;
	}
	cur_roh_heap_removemax(roh_slot_occupied, cur_roh_heap, &cur_roh_heap_top, &marker_uidx2);
      }
    } while ((roh_idx != chrom_roh_start) || (cur_roh_heap_top != 1));
  } while (chrom_fo_idx);
  chrom_fo_idx_to_pidx[0] = pool_ct;

  wptr = uint32toa(pool_ct, g_logbuf);
  if (pool_size_min > 2) {
    wptr = memcpya(wptr, " size-", 6);
    wptr = uint32toa_x(pool_size_min, '+', wptr);
  }
  wptr = memcpya(wptr, " pool", 5);
  if (pool_ct != 1) {
    *wptr++ = 's';
  }
  strcpy(wptr, " of overlapping ROH present.\n");
  logprintb();

  // Now we know how much memory the pools require, so we can assign the rest
  // to a lookahead buffer.
  bigstack_alloc(pool_list_size * sizeof(intptr_t)); // pool_list
  max_lookahead = bigstack_left() / (unfiltered_sample_ctl2 * sizeof(intptr_t));
  lookahead_buf = (uintptr_t*)g_bigstack_base;

  // Now assign ID numbers.
  // We do not precisely imitate PLINK 1.07 here.  This is because
  // 1. it generates non-maximal pools, includes them during ID assignment, and
  //    then winds up with gaps in its final set of ID numbers; this is ugly.
  // 2. it also sorts by *reverse* physical position.  Since we're already
  //    giving up on using diff for compatibility checking, we may as well
  //    switch this to increasing position.  If there's a script lying around
  //    somewhere which actually depends on the reverse order, we can give them
  //    a 'group-reverse' modifier to keep them happy.
  uii = 0;
  for (pool_size = max_pool_size; pool_size >= pool_size_min; --pool_size) {
    pool_list_idx = pool_size_first_plidx[pool_size - pool_size_min];
    while (pool_list_idx != ~ZEROLU) {
#ifdef __LP64__
      pool_list[pool_list_idx + 1] |= ((uintptr_t)(++uii)) << 32;
#else
      pool_list[pool_list_idx + 2] = ++uii;
#endif
      pool_list_idx = pool_list[pool_list_idx];
    }
  }

  // Now form allelic-match groups in an I/O friendly manner (rescan the file
  // in forward order), i.e. traverse pool_list[] in reverse order.
  memcpy(outname_end, ".hom.overlap.S", 14);
  pool_list_idx = pool_list_size;
  fputs("Determining within-pool allelic similarity... [chromosome   ", stdout);
  for (chrom_fo_idx = 0; chrom_fo_idx < chrom_info_ptr->chrom_ct; chrom_fo_idx++) {
    if (chrom_fo_idx_to_pidx[chrom_fo_idx] == chrom_fo_idx_to_pidx[chrom_fo_idx + 1]) {
      continue;
    }
    if (chrom_info_ptr->chrom_file_order[chrom_fo_idx] <= chrom_info_ptr->max_code) {
      printf("\b\b%u] \b", chrom_info_ptr->chrom_file_order[chrom_fo_idx]);
    } else {
      // nonstandard chromosome name
      fputs("\b\b**] \b", stdout);
    }
    fflush(stdout);

    chrom_start = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx];
    marker_uidx2 = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx + 1];
    marker_cidx = 0;
    for (marker_uidx1 = chrom_start; marker_uidx1 < marker_uidx2; marker_uidx1++) {
      if (IS_SET(marker_exclude, marker_uidx1)) {
	continue;
      }
      marker_uidx_to_cidx[marker_uidx1 - chrom_start] = marker_cidx++;
    }
    fill_ulong_zero(max_pool_sizel, roh_slot_occupied);
    // an extra slot because this is a "1-terminated" list
    fill_ulong_one((max_pool_size + 1) * (sizeof(int64_t) / sizeof(intptr_t)), (uintptr_t*)roh_slot_map);
    fill_uint_zero(max_pool_size, roh_slot_end_uidx);
    fill_ulong_zero((((uintptr_t)max_pool_size) * (max_pool_size - 1)) / 2, allelic_match_matrix);
    lookahead_end_uidx = next_unset_unsafe(marker_exclude, chrom_start);
    cur_lookahead_start = 0;
    cur_lookahead_size = 0;
    lookahead_first_cidx = 0;
    roh_slot_map[0] = ~0LLU;

    for (uii = chrom_fo_idx_to_pidx[chrom_fo_idx + 1]; uii < chrom_fo_idx_to_pidx[chrom_fo_idx]; uii++) {
      pool_list_idx = pool_list[pool_list_idx - 1];
      cur_pool = &(pool_list[pool_list_idx]);
      pool_size = (uint32_t)cur_pool[1];
#ifdef __LP64__
      cur_pool = &(cur_pool[2]);
#else
      cur_pool = &(cur_pool[3]);
#endif
      extract_pool_info(pool_size, cur_pool, roh_list, &con_uidx1, &con_uidx2, &union_uidx1, &union_uidx2);

      // flush all roh_slots with end_uidx entries less than or equal to
      // con_uidx2
      if (is_consensus_match) {
	// do not cache any allelic_match_matrix results if using consensus
	// match
	slot_idx1 = 0;
	while (1) {
          slot_idx1 = next_set(roh_slot_occupied, slot_idx1, max_pool_size);
          if (slot_idx1 == max_pool_size) {
	    break;
	  }
          clear_bits((((uintptr_t)slot_idx1) * (slot_idx1 - 1)) / 2, slot_idx1, allelic_match_matrix);
          slot_idx1++;
	}
      } else {
        fill_ulong_zero(BITCT_TO_WORDCT(pool_size), roh_slot_uncached);
      }
      slot_idx1 = 0;
      while (1) {
        slot_idx1 = next_set(roh_slot_occupied, slot_idx1, max_pool_size);
        if (slot_idx1 == max_pool_size) {
	  break;
	}
	if (roh_slot_end_uidx[slot_idx1] <= con_uidx2) {
          CLEAR_BIT(slot_idx1, roh_slot_occupied);
          if (!is_consensus_match) {
            clear_bits((((uintptr_t)slot_idx1) * (slot_idx1 - 1)) / 2, slot_idx1, allelic_match_matrix);
	    slot_idx2 = slot_idx1;
	    while (1) {
              slot_idx2 = next_set(roh_slot_occupied, slot_idx2 + 1, max_pool_size);
	      if (slot_idx2 == max_pool_size) {
		break;
	      }
	      clear_bit_ul(tri_coord_no_diag(slot_idx1, slot_idx2), allelic_match_matrix);
	    }
	  }
	}
	slot_idx1++;
      }
      // now collapse roh_slot_map[]
      slot_idx1 = 0; // read idx
      slot_idx2 = 0; // write idx
      while (1) {
        ullii = roh_slot_map[slot_idx1];
	if (ullii == ~0LLU) {
	  break;
	}
        if (is_set(roh_slot_occupied, (uint32_t)ullii)) {
	  roh_slot_map[slot_idx2++] = ullii;
	}
	slot_idx1++;
      }
      // now fill it with the rest of the new pool's info, and resort to the
      // right order
      slot_idx1 = 0;
      roh_idx = 0;
      while (roh_idx < slot_idx2) {
	sample_uidx1 = (uint32_t)(roh_slot_map[slot_idx1] >> 32);
	while (1) {
	  cur_roh = &(roh_list[cur_pool[roh_idx++] * ROH_ENTRY_INTS]);
	  sample_uidx2 = cur_roh[5];
          if (sample_uidx2 == sample_uidx1) {
	    break;
	  }
	  ujj = next_unset_unsafe(roh_slot_occupied, 0);
	  SET_BIT(ujj, roh_slot_occupied);
	  if (roh_slot_uncached) {
	    SET_BIT(roh_idx - 1, roh_slot_uncached);
	  }
          roh_slot_map[slot_idx2++] = (((uint64_t)sample_uidx2) << 32) | ((uint64_t)ujj);
	  initialize_roh_slot(cur_roh, chrom_start, marker_uidx_to_cidx, &(roh_slots[ujj * roh_slot_wsize]), &(roh_slot_cidx_start[ujj]), &(roh_slot_cidx_end[ujj]), &(roh_slot_end_uidx[ujj]));
	}
	slot_idx1++;
      }
      while (roh_idx < pool_size) {
	cur_roh = &(roh_list[cur_pool[roh_idx] * ROH_ENTRY_INTS]);
        sample_uidx2 = cur_roh[5];
        ujj = next_unset_unsafe(roh_slot_occupied, 0);
        SET_BIT(ujj, roh_slot_occupied);
	if (roh_slot_uncached) {
	  SET_BIT(roh_idx, roh_slot_uncached);
	}
        roh_slot_map[roh_idx++] = (((uint64_t)sample_uidx2) << 32) | ((uint64_t)ujj);
	initialize_roh_slot(cur_roh, chrom_start, marker_uidx_to_cidx, &(roh_slots[ujj * roh_slot_wsize]), &(roh_slot_cidx_start[ujj]), &(roh_slot_cidx_end[ujj]), &(roh_slot_end_uidx[ujj]));
      }
      roh_slot_map[pool_size] = ~0LLU;

      // Now populate the uncached ROH genotype info.
      // 1. fill genotype info from existing lookahead_buf
      // 2. retire SNPs not needed in the future from the buffer
      // 3. if buffer is now empty, use load-and-forget up to con_uidx2
      // 4. load as many SNPs as possible which will be needed in the future
      //    into lookahead_buf
      // 5. fill more genotype info from lookahead_buf
      // 6. resort to load-and-forget for the last SNPs if lookahead_buf is too
      //    small

      if (cur_lookahead_size) {
	populate_roh_slots_from_lookahead_buf(lookahead_buf, max_lookahead, unfiltered_sample_ctl2, cur_lookahead_start, cur_lookahead_size, lookahead_first_cidx, &(roh_slot_map[slot_idx1]), roh_slot_wsize, roh_slot_cidx_start, roh_slot_cidx_end, roh_slots);
      }
      if (!is_verbose) {
	if (con_uidx2 >= lookahead_end_uidx) {
	  lookahead_first_cidx = marker_uidx_to_cidx[con_uidx2 - chrom_start] + 1;
	  cur_lookahead_start = 0;
	  cur_lookahead_size = 0;
	  retval = populate_roh_slots_from_disk(bedfile, bed_offset, rawbuf, marker_exclude, unfiltered_sample_ct4, marker_uidx_to_cidx, chrom_start, lookahead_end_uidx, con_uidx2, &(roh_slot_map[slot_idx1]), roh_slot_wsize, roh_slot_cidx_start, roh_slot_cidx_end, roh_slots, pool_size - slot_idx1);
	  if (retval) {
	    goto roh_pool_ret_1;
	  }
	  lookahead_end_uidx = con_uidx2 + 1;
	} else {
	  ujj = marker_uidx_to_cidx[con_uidx2 - chrom_start] + 1 - lookahead_first_cidx;
	  lookahead_first_cidx += ujj;
	  cur_lookahead_start += ujj;
	  if (cur_lookahead_start >= max_lookahead) {
	    cur_lookahead_start -= max_lookahead;
	  }
	  cur_lookahead_size -= ujj;
	}
      } else {
	if (union_uidx1 >= lookahead_end_uidx) {
	  lookahead_first_cidx = marker_uidx_to_cidx[union_uidx1 - chrom_start];
	  cur_lookahead_start = 0;
	  cur_lookahead_size = 0;
	  lookahead_end_uidx = union_uidx1;
	} else {
	  ujj = marker_uidx_to_cidx[union_uidx1 - chrom_start] - lookahead_first_cidx;
	  lookahead_first_cidx += ujj;
	  cur_lookahead_start += ujj;
	  if (cur_lookahead_start >= max_lookahead) {
	    cur_lookahead_start -= max_lookahead;
	  }
	  cur_lookahead_size -= ujj;
	}
      }

      if (max_lookahead && (lookahead_end_uidx <= union_uidx2)) {
	if (fseeko(bedfile, bed_offset + ((uint64_t)lookahead_end_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	  goto roh_pool_ret_READ_FAIL;
	}
	ulii = cur_lookahead_start + cur_lookahead_size;
	if (ulii >= max_lookahead) {
	  ulii -= max_lookahead;
	}
	do {
	  if (IS_SET(marker_exclude, lookahead_end_uidx)) {
            lookahead_end_uidx = next_unset_unsafe(marker_exclude, lookahead_end_uidx);
            if (fseeko(bedfile, bed_offset + ((uint64_t)lookahead_end_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	      goto roh_pool_ret_READ_FAIL;
	    }
	  }

	  // last few bytes of each lookahead_buf row may be filled with
	  // garbage, but it doesn't matter
	  if (load_raw(unfiltered_sample_ct4, bedfile, &(lookahead_buf[ulii * unfiltered_sample_ctl2]))) {
	    goto roh_pool_ret_READ_FAIL;
	  }
	  ulii++;
	  cur_lookahead_size++;
	  lookahead_end_uidx++;
	  if (ulii == max_lookahead) {
	    ulii = 0;
	  }
	} while ((cur_lookahead_size < max_lookahead) && (lookahead_end_uidx <= union_uidx2));
      }
      if (is_verbose && (cur_lookahead_size == max_lookahead)) {
	goto roh_pool_ret_NOMEM;
      }
      if (cur_lookahead_size) {
        populate_roh_slots_from_lookahead_buf(lookahead_buf, max_lookahead, unfiltered_sample_ctl2, cur_lookahead_start, cur_lookahead_size, lookahead_first_cidx, &(roh_slot_map[slot_idx1]), roh_slot_wsize, roh_slot_cidx_start, roh_slot_cidx_end, roh_slots);
      }
      if (lookahead_end_uidx < union_uidx2) {
        retval = populate_roh_slots_from_disk(bedfile, bed_offset, rawbuf, marker_exclude, unfiltered_sample_ct4, marker_uidx_to_cidx, chrom_start, lookahead_end_uidx, union_uidx2, &(roh_slot_map[slot_idx1]), roh_slot_wsize, roh_slot_cidx_start, roh_slot_cidx_end, roh_slots, pool_size - slot_idx1);
        if (retval) {
	  goto roh_pool_ret_1;
        }
      }

#ifdef __cplusplus
      std::sort((int64_t*)roh_slot_map, (int64_t*)(&(roh_slot_map[pool_size])));
#else
      qsort((int64_t*)roh_slot_map, pool_size, sizeof(int64_t), llcmp);
#endif

      compute_allelic_match_matrix(mismatch_max, roh_slot_wsize, pool_size, roh_slots, roh_slot_occupied, roh_slot_uncached, roh_slot_cidx_start, roh_slot_cidx_end, roh_slot_map, marker_uidx_to_cidx[con_uidx1 - chrom_start], marker_uidx_to_cidx[con_uidx2 - chrom_start] + 1, allelic_match_cts, allelic_match_matrix);

      assign_allelic_match_groups(pool_size, allelic_match_cts, allelic_match_matrix, roh_slot_map, &(cur_pool[pool_size]));

#ifdef __LP64__
      cur_pool[2 * pool_size] = (((uintptr_t)(marker_uidx_to_cidx[union_uidx2 - chrom_start] + 1 - marker_uidx_to_cidx[union_uidx1 - chrom_start])) << 32) | ((uintptr_t)(marker_uidx_to_cidx[con_uidx2 - chrom_start] + 1 - marker_uidx_to_cidx[con_uidx1 - chrom_start]));
#else
      cur_pool[3 * pool_size] = marker_uidx_to_cidx[con_uidx2 - chrom_start] + 1 - marker_uidx_to_cidx[con_uidx1 - chrom_start];
      cur_pool[3 * pool_size + 1] = marker_uidx_to_cidx[union_uidx2 - chrom_start] + 1 - marker_uidx_to_cidx[union_uidx1 - chrom_start];
#endif

      if (is_verbose) {
#ifdef __LP64__
	wptr = uint32toa((uint32_t)(cur_pool[-1] >> 32), &(outname_end[14]));
#else
	wptr = uint32toa((uint32_t)cur_pool[-1], &(outname_end[14]));
#endif
	memcpy(wptr, ".verbose", 9);
	if (fopen_checked(outname, "w", &outfile)) {
	  goto roh_pool_ret_OPEN_FAIL;
	}

	for (slot_idx1 = 0; slot_idx1 < pool_size; slot_idx1++) {
#ifdef __LP64__
	  verbose_group_sort_buf[slot_idx1] = ((cur_pool[pool_size + slot_idx1] & 0x7fffffffLLU) << 32) | ((uint64_t)slot_idx1);
#else
	  verbose_group_sort_buf[slot_idx1] = (((uint64_t)(cur_pool[pool_size + 2 * slot_idx1] & 0x7fffffff)) << 32) | ((uint64_t)slot_idx1);
#endif
	}
#ifdef __cplusplus
	std::sort((int64_t*)verbose_group_sort_buf, (int64_t*)(&(verbose_group_sort_buf[pool_size])));
#else
	qsort((int64_t*)verbose_group_sort_buf, pool_size, sizeof(int64_t), llcmp);
#endif
        sprintf(g_textbuf, "       %%%us %%%us  GRP \n", plink_maxfid, plink_maxiid);
	fprintf(outfile, g_textbuf, "FID", "IID");

	for (slot_idx1 = 0; slot_idx1 < pool_size; slot_idx1++) {
	  slot_idx2 = (uint32_t)verbose_group_sort_buf[slot_idx1];
	  roh_idx = cur_pool[slot_idx2];
	  cur_roh = &(roh_list[roh_idx * ROH_ENTRY_INTS]);
	  verbose_uidx_bounds[slot_idx1 * 2] = cur_roh[0];
	  verbose_uidx_bounds[slot_idx1 * 2 + 1] = cur_roh[1];
	  verbose_sample_uidx[slot_idx1] = cur_roh[5];
	  sample_uidx1 = cur_roh[5];
	  wptr = uint32toa(slot_idx1 + 1, g_textbuf);
          wptr = width_force(4, g_textbuf, wptr);
	  wptr = memcpyl3a(wptr, ")  ");
	  cptr = &(sample_ids[sample_uidx1 * max_sample_id_len]);
	  cptr2 = (char*)memchr(cptr, '\t', max_sample_id_len);
          wptr = fw_strcpyn(plink_maxfid, cptr2 - cptr, cptr, wptr);
          *wptr++ = ' ';
          wptr = fw_strcpy(plink_maxiid, &(cptr2[1]), wptr);
          wptr = memseta(wptr, 32, 3);
          wptr = uint32toa((uint32_t)(verbose_group_sort_buf[slot_idx1] >> 32), wptr);
	  *wptr++ = '\n';
	  if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
	    goto roh_pool_ret_WRITE_FAIL;
	  }
	}
	putc_unlocked('\n', outfile);
	wptr = memseta(g_textbuf, 32, plink_maxsnp - 3);
	wptr = memcpya(wptr, "SNP ", 4);
        fwrite(g_textbuf, 1, wptr - g_textbuf, outfile);
	for (slot_idx1 = 0; slot_idx1 < pool_size; slot_idx1++) {
	  wptr = uint32toa(slot_idx1 + 1, g_textbuf);
          wptr = width_force(4, g_textbuf, wptr);
	  wptr = memseta(wptr, 32, 2);
	  fwrite(g_textbuf, 1, wptr - g_textbuf, outfile);
	}
	if (fputs_checked("\n\n", outfile)) {
	  goto roh_pool_ret_WRITE_FAIL;
	}
	marker_cidx = marker_uidx_to_cidx[union_uidx1 - chrom_start];
	for (marker_uidx1 = union_uidx1; marker_uidx1 <= union_uidx2; marker_uidx1++) {
	  next_unset_unsafe_ck(marker_exclude, &marker_uidx1);
	  if (marker_uidx1 == con_uidx1) {
	    putc_unlocked('\n', outfile);
	  }
	  ulii = marker_cidx + cur_lookahead_start - lookahead_first_cidx;
	  if (ulii >= max_lookahead) {
	    ulii -= max_lookahead;
	  }
	  lookahead_row = &(lookahead_buf[ulii * unfiltered_sample_ctl2]);
	  wptr = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx1 * max_marker_id_len]), g_textbuf);
	  *wptr++ = ' ';
          fwrite(g_textbuf, 1, wptr - g_textbuf, outfile);
	  allele_strs[2] = marker_allele_ptrs[marker_uidx1 * 2];
	  allele_strs[3] = marker_allele_ptrs[marker_uidx1 * 2 + 1];
	  allele_strs[0] = allele_strs[2 + IS_SET(marker_reverse, marker_uidx1)];
	  allele_strs[1] = allele_strs[3 - IS_SET(marker_reverse, marker_uidx1)];
	  for (slot_idx1 = 0; slot_idx1 < pool_size; slot_idx1++) {
	    sample_uidx1 = verbose_sample_uidx[slot_idx1];
	    ujj = ((marker_uidx1 >= verbose_uidx_bounds[slot_idx1 * 2]) && (marker_uidx1 <= verbose_uidx_bounds[slot_idx1 * 2 + 1]));
	    if (ujj) {
	      writebuf[0] = '[';
	    } else {
	      writebuf[0] = ' ';
	    }
	    ulii = (lookahead_row[(sample_uidx1 / BITCT2)] >> (2 * (sample_uidx1 % BITCT2))) & (3 * ONELU);
	    wptr = &(writebuf[1]);
	    if (ulii == 1) {
	      wptr = memcpyl3a(wptr, "0/0");
	    } else {
	      if (ulii == 3) {
		wptr = strcpyax(wptr, allele_strs[1], '/');
		wptr = strcpya(wptr, allele_strs[1]);
	      } else if (ulii) {
		wptr = strcpyax(wptr, allele_strs[2], '/');
		wptr = strcpya(wptr, allele_strs[3]);
	      } else {
		wptr = strcpyax(wptr, allele_strs[0], '/');
		wptr = strcpya(wptr, allele_strs[0]);
	      }
	    }
	    if (ujj) {
	      *wptr++ = ']';
	    } else {
	      *wptr++ = ' ';
	    }
	    *wptr++ = ' ';
	    fwrite(writebuf, 1, wptr - writebuf, outfile);
	  }
          if (putc_checked('\n', outfile)) {
	    goto roh_pool_ret_WRITE_FAIL;
	  }
	  if (marker_uidx1 == con_uidx2) {
	    putc_unlocked('\n', outfile);
	  }
	  marker_cidx++;
	}

	if (fputs_checked("\n\n", outfile)) {
	  goto roh_pool_ret_WRITE_FAIL;
	}

	slot_idx1 = 0;
	do {
	  group_slot_end = slot_idx1 + 1;
	  ujj = (uint32_t)(verbose_group_sort_buf[slot_idx1] >> 32);
	  while ((group_slot_end < pool_size) && (((uint32_t)(verbose_group_sort_buf[group_slot_end] >> 32)) == ujj)) {
	    group_slot_end++;
	  }
	  wptr = memcpya(g_textbuf, "Group ", 6);
	  wptr = uint32toa(ujj, wptr);
	  wptr = memcpya(wptr, "\n\n", 2);
	  if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
	    goto roh_pool_ret_WRITE_FAIL;
	  }
	  for (slot_idx2 = slot_idx1; slot_idx2 < group_slot_end; slot_idx2++) {
            wptr = uint32toa_w4(slot_idx2 + 1, g_textbuf);
	    wptr = memcpya(wptr, ") ", 2);
	    sample_uidx1 = verbose_sample_uidx[slot_idx2];
	    cptr = &(sample_ids[sample_uidx1 * max_sample_id_len]);
	    cptr2 = (char*)memchr(cptr, '\t', max_sample_id_len);
	    wptr = fw_strcpyn(plink_maxfid, cptr2 - cptr, cptr, wptr);
	    *wptr++ = ' ';
	    wptr = fw_strcpy(plink_maxiid, &(cptr2[1]), wptr);
	    *wptr++ = ' ';
	    if (IS_SET(pheno_nm, sample_uidx1)) {
	      if (pheno_c) {
		wptr = memseta(wptr, 32, 7);
		*wptr++ = '1' + IS_SET(pheno_c, sample_uidx1);
	      } else {
		wptr = dtoa_g_wxp2(pheno_d[sample_uidx1], 8, wptr);
	      }
	    } else {
              wptr = fw_strcpyn(8, missing_pheno_len, missing_pheno_str, wptr);
	    }
	    *wptr++ = '\n';
	    fwrite(g_textbuf, 1, wptr - g_textbuf, outfile);
	  }
	  if (fputs_checked("\n\n", outfile)) {
	    goto roh_pool_ret_WRITE_FAIL;
	  }
	  wptr = memseta(g_textbuf, 32, plink_maxsnp - 3);
	  wptr = memcpya(wptr, "SNP         ", 12);
	  if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
	    goto roh_pool_ret_WRITE_FAIL;
	  }
	  for (slot_idx2 = slot_idx1; slot_idx2 < group_slot_end; slot_idx2++) {
	    wptr = uint32toa_w4(slot_idx2 + 1, g_textbuf);
	    wptr = memseta(wptr, 32, 2);
	    fwrite(g_textbuf, 1, wptr - g_textbuf, outfile);
	  }
	  if (fputs_checked("\n\n", outfile)) {
	    goto roh_pool_ret_WRITE_FAIL;
	  }
	  marker_cidx = marker_uidx_to_cidx[union_uidx1 - chrom_start];
	  for (marker_uidx1 = union_uidx1; marker_uidx1 <= union_uidx2; marker_uidx1++) {
	    next_unset_unsafe_ck(marker_exclude, &marker_uidx1);
	    if (marker_uidx1 == con_uidx1) {
	      putc_unlocked('\n', outfile);
	    }
	    ulii = marker_cidx + cur_lookahead_start - lookahead_first_cidx;
	    if (ulii >= max_lookahead) {
	      ulii -= max_lookahead;
	    }
	    lookahead_row = &(lookahead_buf[ulii * unfiltered_sample_ctl2]);
	    wptr = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx1 * max_marker_id_len]), g_textbuf);
	    *wptr++ = ' ';
	    if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
	      goto roh_pool_ret_WRITE_FAIL;
	    }
	    ujj = 0; // A1 hom ct
	    ukk = 0; // A2 hom ct
	    for (slot_idx2 = slot_idx1; slot_idx2 < group_slot_end; slot_idx2++) {
	      if ((marker_uidx1 >= verbose_uidx_bounds[slot_idx2 * 2]) && (marker_uidx1 <= verbose_uidx_bounds[slot_idx2 * 2 + 1])) {
	        sample_uidx1 = verbose_sample_uidx[slot_idx2];
		ulii = (lookahead_row[sample_uidx1 / BITCT2] >> (2 * (sample_uidx1 % BITCT2))) & (3 * ONELU);
		// no need to actually count hets here
		if (ulii == 3) {
		  ukk++;
		} else if (!ulii) {
		  ujj++;
		}
	      }
	    }
	    allele_strs[2] = marker_allele_ptrs[marker_uidx1 * 2];
	    allele_strs[3] = marker_allele_ptrs[marker_uidx1 * 2 + 1];
	    allele_strs[0] = allele_strs[2 + IS_SET(marker_reverse, marker_uidx1)];
	    allele_strs[1] = allele_strs[3 - IS_SET(marker_reverse, marker_uidx1)];
	    if (ukk > ujj) {
	      wptr = fw_strcpy(2, allele_strs[1], writebuf);
	    } else if (ujj > ukk) {
	      wptr = fw_strcpy(2, allele_strs[0], writebuf);
	    } else {
	      wptr = memcpya(writebuf, " ?", 2);
	    }
	    wptr = memseta(wptr, 32, 6);
	    if (fwrite_checked(writebuf, wptr - writebuf, outfile)) {
	      goto roh_pool_ret_WRITE_FAIL;
	    }
	    for (slot_idx2 = slot_idx1; slot_idx2 < group_slot_end; slot_idx2++) {
	      sample_uidx1 = verbose_sample_uidx[slot_idx2];
	      ujj = ((marker_uidx1 >= verbose_uidx_bounds[slot_idx2 * 2]) && (marker_uidx1 <= verbose_uidx_bounds[slot_idx2 * 2 + 1]));
	      if (ujj) {
		writebuf[0] = '[';
	      } else {
		writebuf[0] = ' ';
	      }
	      ulii = (lookahead_row[(sample_uidx1 / BITCT2)] >> (2 * (sample_uidx1 % BITCT2))) & (3 * ONELU);
	      wptr = &(writebuf[1]);
	      if (ulii == 3) {
		wptr = strcpyax(wptr, allele_strs[1], '/');
		wptr = strcpya(wptr, allele_strs[1]);
	      } else if (ulii == 2) {
		wptr = strcpyax(wptr, allele_strs[2], '/');
		wptr = strcpya(wptr, allele_strs[3]);
	      } else if (ulii) {
		wptr = memcpyl3a(wptr, "0/0");
	      } else {
		wptr = strcpyax(wptr, allele_strs[0], '/');
		wptr = strcpya(wptr, allele_strs[0]);
	      }
	      if (ujj) {
		*wptr++ = ']';
	      } else {
		*wptr++ = ' ';
	      }
	      *wptr++ = ' ';
	      fwrite(writebuf, 1, wptr - writebuf, outfile);
	    }
	    if (putc_checked('\n', outfile)) {
	      goto roh_pool_ret_WRITE_FAIL;
	    }
	    if (marker_uidx1 == con_uidx2) {
	      putc_unlocked('\n', outfile);
	    }
	    marker_cidx++;
	  }
	  if (putc_checked('\n', outfile)) {
	    goto roh_pool_ret_WRITE_FAIL;
	  }
	  slot_idx1 = group_slot_end;
	} while (slot_idx1 < pool_size);

        if (fputs_checked("\n\n", outfile)) {
	  goto roh_pool_ret_WRITE_FAIL;
	}

	marker_cidx = marker_uidx_to_cidx[union_uidx1 - chrom_start];
	for (marker_uidx1 = union_uidx1; marker_uidx1 <= union_uidx2; marker_uidx1++) {
	  next_unset_unsafe_ck(marker_exclude, &marker_uidx1);
	  if (marker_uidx1 == con_uidx1) {
	    putc_unlocked('\n', outfile);
	  }
	  ulii = marker_cidx + cur_lookahead_start - lookahead_first_cidx;
	  if (ulii >= max_lookahead) {
	    ulii -= max_lookahead;
	  }
          lookahead_row = &(lookahead_buf[ulii * unfiltered_sample_ctl2]);
	  wptr = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx1 * max_marker_id_len]), g_textbuf);
	  *wptr++ = ' ';
	  if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
	    goto roh_pool_ret_WRITE_FAIL;
	  }

	  allele_strs[0] = marker_allele_ptrs[marker_uidx1 * 2];
	  allele_strs[1] = marker_allele_ptrs[marker_uidx1 * 2 + 1];
	  slot_idx1 = 0;
	  do {
	    group_slot_end = slot_idx1 + 1;
	    ujj = (uint32_t)(verbose_group_sort_buf[slot_idx1] >> 32);
	    while ((group_slot_end < pool_size) && (((uint32_t)(verbose_group_sort_buf[group_slot_end] >> 32)) == ujj)) {
	      group_slot_end++;
	    }
	    // to conserve memory, recalculate consensus haplotypes
	    ujj = 0;
	    ukk = 0;
	    for (slot_idx2 = slot_idx1; slot_idx2 < group_slot_end; slot_idx2++) {
	      if ((marker_uidx1 >= verbose_uidx_bounds[slot_idx2 * 2]) && (marker_uidx1 <= verbose_uidx_bounds[slot_idx2 * 2 + 1])) {
		sample_uidx1 = verbose_sample_uidx[slot_idx2];
		ulii = (lookahead_row[sample_uidx1 / BITCT2] >> (2 * (sample_uidx1 % BITCT2))) & (3 * ONELU);
		if (ulii == 3) {
		  ukk++;
		} else if (!ulii) {
		  ujj++;
		}
	      }
	    }
	    if (ukk > ujj) {
	      fputs(allele_strs[1], outfile);
	    } else if (ujj > ukk) {
	      fputs(allele_strs[0], outfile);
	    } else {
	      putc_unlocked('?', outfile);
	    }
	    putc_unlocked(' ', outfile);
	    slot_idx1 = group_slot_end;
	  } while (slot_idx1 < pool_size);
          if (putc_checked('\n', outfile)) {
	    goto roh_pool_ret_WRITE_FAIL;
	  }
	  if (marker_uidx1 == con_uidx2) {
	    putc_unlocked('\n', outfile);
	  }
	  marker_cidx++;
	}

	if (fclose_null(&outfile)) {
	  goto roh_pool_ret_WRITE_FAIL;
	}
	LOGPREPRINTFWW("%s written.\n", outname);
        logstr(g_logbuf);
      }
    }
    if (chrom_info_ptr->chrom_file_order[chrom_fo_idx] > onechar_max) {
      putc_unlocked('\b', stdout);
    }
  }
  fputs("\b\b\b\b\b\b\b\b\b\b\b\b\b\b               \b\b\b\b\b\b\b\b\b\b\b\b\b\b\bdone.\n", stdout);

  outname_end[12] = '\0';
  if (fopen_checked(outname, "w", &outfile)) {
    goto roh_pool_ret_OPEN_FAIL;
  }
  sprintf(g_textbuf, " POOL %%%us %%%us      PHE  CHR %%%us %%%us            BP1            BP2       KB     NSNP NSIM    GRP\n", plink_maxfid, plink_maxiid, plink_maxsnp, plink_maxsnp);
  fprintf(outfile, g_textbuf, "FID", "IID", "SNP1", "SNP2");
  uii = 1; // pool ID
  fputs("Writing...", stdout);
  fflush(stdout);
  for (pool_size = max_pool_size; pool_size >= pool_size_min; pool_size--) {
    pool_list_idx = pool_size_first_plidx[pool_size - pool_size_min];
    while (pool_list_idx != ~ZEROLU) {
      cur_pool = &(pool_list[pool_list_idx]);
      pool_list_idx = *cur_pool;
#ifdef __LP64__
      cur_pool = &(cur_pool[2]);
#else
      cur_pool = &(cur_pool[3]);
#endif
      case_ct = 0;
      g_textbuf[0] = 'S';
      wptr_start = width_force(5, g_textbuf, uint32toa(uii, &(g_textbuf[1])));
      *wptr_start++ = ' ';
      cur_roh = &(roh_list[cur_pool[0] * ROH_ENTRY_INTS]);
      con_uidx1 = cur_roh[0];
      union_uidx1 = cur_roh[0];
      con_uidx2 = cur_roh[1];
      union_uidx2 = cur_roh[1];
      chrom_start = get_variant_chrom(chrom_info_ptr, con_uidx1);
      // sort pool members primarily by allelic-match group number, then by
      // internal ID
      for (slot_idx1 = 0; slot_idx1 < pool_size; slot_idx1++) {
#ifdef __LP64__
	roh_slot_map[slot_idx1] = ((cur_pool[pool_size + slot_idx1] & 0x7fffffffLLU) << 32) | ((uint64_t)slot_idx1);
#else
	// would like to just sort 32-bit integers, but if there are >32k
	// allelic-match groups it won't work due to signed integer overflow
        roh_slot_map[slot_idx1] = (((uint64_t)(cur_pool[pool_size + 2 * slot_idx1] & 0x7fffffff)) << 32) | ((uint64_t)slot_idx1);
#endif
      }
#ifdef __cplusplus
      std::sort((int64_t*)roh_slot_map, (int64_t*)(&(roh_slot_map[pool_size])));
#else
      qsort((int64_t*)roh_slot_map, pool_size, sizeof(int64_t), llcmp);
#endif
      for (slot_idx1 = 0; slot_idx1 < pool_size; slot_idx1++) {
	slot_idx2 = (uint32_t)roh_slot_map[slot_idx1];
        roh_idx = cur_pool[slot_idx2];
	cur_roh = &(roh_list[roh_idx * ROH_ENTRY_INTS]);
	sample_uidx1 = cur_roh[5];
	cptr = &(sample_ids[sample_uidx1 * max_sample_id_len]);
	cptr2 = (char*)memchr(cptr, '\t', max_sample_id_len);
	wptr = fw_strcpyn(plink_maxfid, cptr2 - cptr, cptr, wptr_start);
	*wptr++ = ' ';
	wptr = fw_strcpy(plink_maxiid, &(cptr2[1]), wptr);
	*wptr++ = ' ';
	if (IS_SET(pheno_nm, sample_uidx1)) {
	  if (pheno_c) {
	    if (IS_SET(pheno_c, sample_uidx1)) {
	      case_ct++;
              wptr = memcpya(wptr, "       2", 8);
	    } else {
              wptr = memcpya(wptr, "       1", 8);
	    }
	  } else {
	    wptr = dtoa_g_wxp4(pheno_d[sample_uidx1], 8, wptr);
	  }
	} else {
          wptr = fw_strcpyn(8, missing_pheno_len, missing_pheno_str, wptr);
	}
	*wptr++ = ' ';
	wptr = width_force(4, wptr, chrom_name_write(chrom_info_ptr, chrom_start, wptr));
	marker_uidx1 = cur_roh[0];
	marker_uidx2 = cur_roh[1];
	if (marker_uidx1 > con_uidx1) {
	  con_uidx1 = marker_uidx1;
	} else if (marker_uidx1 < union_uidx1) {
	  union_uidx1 = marker_uidx1;
	}
	if (marker_uidx2 < con_uidx2) {
	  con_uidx2 = marker_uidx2;
	} else if (marker_uidx2 > union_uidx2) {
	  union_uidx2 = marker_uidx2;
	}
        wptr = roh_pool_write_middle(wptr, marker_ids, max_marker_id_len, plink_maxsnp, marker_pos, is_new_lengths, marker_uidx1, marker_uidx2);
	wptr = uint32toa_w8x(cur_roh[2], ' ', wptr);
#ifdef __LP64__
	ulii = cur_pool[pool_size + slot_idx2];
        wptr = uint32toa_w4x((uint32_t)(ulii >> 32), ' ', wptr);
        wptr = width_force(5, wptr, uint32toa(ulii & 0x7fffffff, wptr));
        if (ulii & 0x80000000LLU) {
          *wptr++ = '*';
	} else {
	  *wptr++ = ' ';
	}
#else
	ulii = cur_pool[pool_size + 2 * slot_idx2];
        wptr = uint32toa_w4x(cur_pool[pool_size + 2 * slot_idx2 + 1], ' ', wptr);
        wptr = width_force(5, wptr, uint32toa(ulii & 0x7fffffff, wptr));
	if (ulii & 0x80000000U) {
	  *wptr++ = '*';
	} else {
	  *wptr++ = ' ';
	}
#endif
        wptr = memcpya(wptr, " \n", 2);
	if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
	  goto roh_pool_ret_WRITE_FAIL;
	}
      }
      for (ujj = 0; ujj < 2; ujj++) {
	if (!ujj) {
	  wptr = fw_strcpyn(plink_maxfid, 3, "CON", wptr_start);
	  marker_uidx1 = con_uidx1;
	  marker_uidx2 = con_uidx2;
#ifdef __LP64__
	  marker_cidx = (uint32_t)(cur_pool[2 * pool_size]);
#else
	  marker_cidx = cur_pool[3 * pool_size];
#endif
	} else {
	  wptr = fw_strcpyn(plink_maxfid, 5, "UNION", wptr_start);
	  marker_uidx1 = union_uidx1;
	  marker_uidx2 = union_uidx2;
#ifdef __LP64__
	  // NSNP
	  marker_cidx = (uint32_t)(cur_pool[2 * pool_size] >> 32);
#else
	  marker_cidx = cur_pool[3 * pool_size + 1];
#endif
	}
        *wptr++ = ' ';
	wptr = width_force(plink_maxiid, wptr, uint32toa(pool_size, wptr));
        *wptr++ = ' ';
        cptr = uint32toa(case_ct, wptr);
	*cptr++ = ':';
	cptr = uint32toa(pool_size - case_ct, cptr);
        wptr = width_force(8, wptr, cptr);
	*wptr++ = ' ';
	wptr = width_force(4, wptr, chrom_name_write(chrom_info_ptr, chrom_start, wptr));
        wptr = roh_pool_write_middle(wptr, marker_ids, max_marker_id_len, plink_maxsnp, marker_pos, is_new_lengths, marker_uidx1, marker_uidx2);
        wptr = uint32toa_w8(marker_cidx, wptr);
        wptr = memcpya(wptr, "    NA     NA \n", 15);
	if (ujj) {
	  *wptr++ = '\n';
	}
	if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
	  goto roh_pool_ret_WRITE_FAIL;
	}
      }
      uii++;
    }
  }
  if (fclose_null(&outfile)) {
    goto roh_pool_ret_WRITE_FAIL;
  }

  putc_unlocked('\r', stdout);
  LOGPRINTFWW("ROH pool report written to %s .\n", outname);
  if (is_verbose) {
    wptr = strcpya(g_logbuf, "Per-pool report");
    if (pool_ct != 1) {
      *wptr++ = 's';
    }
    wptr = strcpya(wptr, " written to ");
    wptr = strcpya(wptr, outname);
    wptr = memcpya(wptr, ".S", 2);
    if (pool_ct == 1) {
      *wptr++ = '1';
    } else if (pool_ct == 2) {
      wptr = memcpya(wptr, "{1,2}", 5);
    } else {
      wptr = memcpya(wptr, "{1,...,", 7);
      wptr = uint32toa(pool_ct, wptr);
      *wptr++ = '}';
    }
    wptr = memcpya(wptr, ".verbose.\n", 11);
    fputs(g_logbuf, stdout);
  }

  while (0) {
  roh_pool_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  roh_pool_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  roh_pool_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  roh_pool_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
 roh_pool_ret_1:
  bigstack_reset(bigstack_mark);
  fclose_cond(outfile);
  return retval;
}

int32_t calc_homozyg(Homozyg_info* hp, FILE* bedfile, uintptr_t bed_offset, uint32_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, uintptr_t* marker_reverse, Chrom_info* chrom_info_ptr, uint32_t* marker_pos, uintptr_t sample_ct, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, char* sample_ids, uint32_t plink_maxfid, uint32_t plink_maxiid, uintptr_t max_sample_id_len, char* outname, char* outname_end, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, char* output_missing_pheno, uintptr_t* sex_male) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  uint64_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t unfiltered_sample_ctl2 = QUATERCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t sample_ctl = BITCT_TO_WORDCT(sample_ct);
  uintptr_t sample_ctl2 = QUATERCT_TO_WORDCT(sample_ct);
  uintptr_t window_size = hp->window_size;
  double hit_threshold = hp->hit_threshold;
  uint32_t is_new_lengths = 1 ^ ((hp->modifier / HOMOZYG_OLD_LENGTHS) & 1);
  uint32_t chrom_ct = chrom_info_ptr->chrom_ct;
  int32_t x_code = chrom_info_ptr->xymt_codes[X_OFFSET];
  int32_t mt_code = chrom_info_ptr->xymt_codes[MT_OFFSET];
  uintptr_t* haploid_mask = chrom_info_ptr->haploid_mask;
  uintptr_t roh_ct = 0;
  uintptr_t final_mask = get_final_mask(sample_ct);
  uintptr_t* sample_male = nullptr;

  // support for new 'extend' modifier
  uint32_t* prev_roh_end_cidxs = nullptr;
  uint32_t* end_nonhom_uidxs = nullptr;
  uint32_t* cur_roh_earliest_extend_uidxs = nullptr;

  uint32_t swhit_min = 0;
  int32_t retval = 0;
  char missing_pheno_str[35];
  uint32_t missing_pheno_len;
  uintptr_t* roh_list_chrom_starts;
  uintptr_t* rawbuf;
  uintptr_t* readbuf; // circular window of actual genotype data
  uintptr_t* swbuf; // circular window of recent scanning window hits
  uint32_t* het_cts;
  uint32_t* missing_cts;
  uint32_t* swhit_cts;
  uint32_t* cur_roh_uidx_starts;
  // probably want to rewrite this code to mostly just refer to cidxs...
  uint32_t* cur_roh_cidx_starts;
  uint32_t* cur_roh_het_cts;
  uint32_t* cur_roh_missing_cts;
  uintptr_t* sample_to_last_roh;

  // 6 uint32_ts per entry in 32-bit build, 7 on 64-bit
  // [0]: first uidx
  // [1]: last uidx
  // [2]: nsnp
  // [3]: hom ct
  // [4]: het ct
  // [5] (and usually [6]): uintptr_t indicating position of sample's previous
  // ROH, or 0xfff... if there is none; becomes indiv_idx after .hom.indiv file
  // written
  // Note that this is sorted by primarily by LAST uidx, and secondarily by
  // sample_idx.
  uint32_t* roh_list;

  uint32_t* uidx_buf; // circular buffer tracking most recent marker_uidxs
  uintptr_t* readbuf_cur;
  uintptr_t* swbuf_cur;
  uintptr_t* cur_sample_male;
  char* wptr;
  uintptr_t max_roh_ct;
  uintptr_t marker_uidx;
  uintptr_t chrom_end;
  uintptr_t widx;
  uintptr_t widx_next;
  uintptr_t sample_idx;
  uintptr_t ulii;
  double dxx;
  uint32_t chrom_fo_idx;
  uint32_t old_uidx;
  uint32_t older_uidx;
  uint32_t marker_cidx;
  uint32_t marker_cidx_max;
  uint32_t swbuf_full;
  uint32_t max_pool_size;
  uint32_t max_roh_len;
  uint32_t omp_is_numeric;
  uint32_t uii;
  missing_pheno_len = strlen(output_missing_pheno);
  memcpy(missing_pheno_str, output_missing_pheno, missing_pheno_len);
  omp_is_numeric = scan_double(output_missing_pheno, &dxx)? 0 : 1;
  if (omp_is_numeric) {
    wptr = strchr(output_missing_pheno, '.');
    if (!wptr) {
      memcpy(&(missing_pheno_str[missing_pheno_len]), ".000", 4);
      missing_pheno_len += 4;
    } else {
      // enforce 3 decimal places
      uii = (uintptr_t)(wptr - output_missing_pheno);
      if (missing_pheno_len - uii < 4) {
	memset(&(missing_pheno_str[missing_pheno_len]), '0', uii + 4 - missing_pheno_len);
      } else {
	missing_pheno_len = uii + 4;
      }
    }
  }

  if (bigstack_alloc_ul(chrom_ct + 1, &roh_list_chrom_starts) ||
      bigstack_alloc_ul(unfiltered_sample_ctl2, &rawbuf) ||
      bigstack_end_alloc_ul(sample_ctl2 * window_size, &readbuf) ||
      bigstack_end_alloc_ul(sample_ctl * window_size, &swbuf) ||
      bigstack_end_alloc_ui(sample_ct, &het_cts) ||
      bigstack_end_alloc_ui(sample_ct, &missing_cts) ||
      bigstack_end_alloc_ui(sample_ct, &swhit_cts) ||
      bigstack_end_alloc_ui(sample_ct, &cur_roh_uidx_starts) ||
      bigstack_end_alloc_ui(sample_ct, &cur_roh_cidx_starts)) {
    goto calc_homozyg_ret_NOMEM;
  }
  if (hp->modifier & HOMOZYG_EXTEND) {
    if (bigstack_end_alloc_ui(sample_ct, &prev_roh_end_cidxs) ||
        bigstack_end_alloc_ui(sample_ct, &end_nonhom_uidxs) ||
        bigstack_end_alloc_ui(sample_ct, &cur_roh_earliest_extend_uidxs)) {
      goto calc_homozyg_ret_NOMEM;
    }
  }
  if (bigstack_end_alloc_ui(sample_ct, &cur_roh_het_cts) ||
      bigstack_end_alloc_ui(sample_ct, &cur_roh_missing_cts) ||
      bigstack_end_alloc_ul(sample_ct, &sample_to_last_roh) ||
      bigstack_end_alloc_ui(window_size, &uidx_buf)) {
    goto calc_homozyg_ret_NOMEM;
  }
  if ((x_code != -2) && is_set(chrom_info_ptr->chrom_mask, x_code)) {
    if (bigstack_end_alloc_ul(sample_ctl, &sample_male)) {
      goto calc_homozyg_ret_NOMEM;
    }
    copy_bitarr_subset_excl(sex_male, sample_exclude, sample_ct, popcount_longs_exclude(sex_male, sample_exclude, sample_ctl), sample_male);
  }
  // no other workspace allocations during main scan, so we can assign it all
  // to the ROH list
  max_roh_ct = (bigstack_left() & (~(CACHELINE - ONELU))) / (ROH_ENTRY_INTS * sizeof(int32_t));
  roh_list = (uint32_t*)g_bigstack_base;
  ulii = sample_ctl2 - 1;
  rawbuf[unfiltered_sample_ctl2 - 1] = 0;
  for (widx = 0; widx < window_size; widx++) {
    readbuf[widx * sample_ctl2 + ulii] = 0;
  }
  if (fseeko(bedfile, bed_offset, SEEK_SET)) {
    goto calc_homozyg_ret_READ_FAIL;
  }
  fill_ulong_one(sample_ct, sample_to_last_roh);

  for (chrom_fo_idx = 0; chrom_fo_idx < chrom_ct; chrom_fo_idx++) {
    uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
    roh_list_chrom_starts[chrom_fo_idx] = roh_ct;
    chrom_end = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx + 1];
    if ((x_code == -2) || (uii != ((uint32_t)x_code))) {
      if (IS_SET(haploid_mask, uii) || (uii == (uint32_t)mt_code)) {
	marker_uidx = chrom_end;
	if (fseeko(bedfile, bed_offset + (uint64_t)marker_uidx * unfiltered_sample_ct4, SEEK_SET)) {
	  goto calc_homozyg_ret_READ_FAIL;
	}
	continue;
      }
      cur_sample_male = nullptr;
    } else {
      cur_sample_male = sample_male;
    }
    if (uii <= chrom_info_ptr->max_code) {
      printf("\r--homozyg: Scanning chromosome %u. \b", uii);
    } else {
      fputs("\r--homozyg: Scanning chromosome **.\b", stdout);
    }
    fflush(stdout);
    marker_uidx = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx];
    fill_ulong_zero(sample_ctl * window_size, swbuf);
    fill_uint_zero(sample_ct, het_cts);
    fill_uint_zero(sample_ct, missing_cts);
    for (widx = 0; widx < window_size; widx++) {
      if (IS_SET(marker_exclude, marker_uidx)) {
        marker_uidx = next_unset_ul(marker_exclude, marker_uidx, chrom_end);
        if (fseeko(bedfile, bed_offset + (uint64_t)marker_uidx * unfiltered_sample_ct4, SEEK_SET)) {
	  goto calc_homozyg_ret_READ_FAIL;
	}
      }
      if (marker_uidx == chrom_end) {
	break;
      }
      readbuf_cur = &(readbuf[widx * sample_ctl2]);
      if (load_and_collapse(unfiltered_sample_ct, sample_ct, sample_exclude, final_mask, 0, bedfile, rawbuf, readbuf_cur)) {
	goto calc_homozyg_ret_READ_FAIL;
      }
      mask_out_homozyg_major(readbuf_cur, sample_ct);
      increment_het_missing(readbuf_cur, sample_ct, het_cts, missing_cts, 1);
      uidx_buf[widx] = marker_uidx++;
    }
    if (widx == window_size) {
      marker_uidx--;
      fill_ulong_zero(window_size * sample_ctl, swbuf);
      fill_uint_zero(sample_ct, swhit_cts);
      fill_uint_one(sample_ct, cur_roh_cidx_starts);

      widx = 0;
      swbuf_full = 0;
      marker_cidx = 0;
      old_uidx = uidx_buf[0];
      if (prev_roh_end_cidxs) {
        fill_uint_one(sample_ct, prev_roh_end_cidxs);
	for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
          end_nonhom_uidxs[sample_idx] = old_uidx;
	}
      }
      // repurpose widx as next readbuf row
      while (1) {
	older_uidx = old_uidx;
	old_uidx = uidx_buf[widx];
	readbuf_cur = &(readbuf[widx * sample_ctl2]);
	swbuf_cur = &(swbuf[widx * sample_ctl]);
	if (swbuf_full) {
	  vertical_bitct_subtract(swbuf_cur, sample_ct, swhit_cts);
	  fill_ulong_zero(sample_ctl, swbuf_cur);
	} else {
	  swhit_min = (int32_t)(((double)((int32_t)(widx + 1))) * hit_threshold + 1.0 - EPSILON);
	}
	if (roh_update(hp, readbuf_cur, swbuf_cur, het_cts, missing_cts, marker_pos, cur_sample_male, sample_ct, swhit_min, older_uidx, old_uidx, marker_cidx, max_roh_ct, swhit_cts, cur_roh_uidx_starts, cur_roh_cidx_starts, cur_roh_het_cts, cur_roh_missing_cts, sample_to_last_roh, roh_list, &roh_ct, marker_exclude, prev_roh_end_cidxs, end_nonhom_uidxs, cur_roh_earliest_extend_uidxs)) {
	  goto calc_homozyg_ret_NOMEM;
	}
	widx_next = widx + 1;
	if (widx_next == window_size) {
	  widx_next = 0;
	  swbuf_full = 1;
	}
	if (!prev_roh_end_cidxs) {
	  increment_het_missing(readbuf_cur, sample_ct, het_cts, missing_cts, -1);
	} else {
	  decrement_het_missing_extend(readbuf_cur, sample_ct, het_cts, missing_cts, end_nonhom_uidxs, uidx_buf[widx_next]);
	}
	if (++marker_uidx == chrom_end) {
	  break;
	}
	if (IS_SET(marker_exclude, marker_uidx)) {
	  marker_uidx = next_unset_ul(marker_exclude, marker_uidx, chrom_end);
	  if (fseeko(bedfile, bed_offset + (uint64_t)marker_uidx * unfiltered_sample_ct4, SEEK_SET)) {
	    goto calc_homozyg_ret_READ_FAIL;
	  }
	  if (marker_uidx == chrom_end) {
	    break;
	  }
	}
	uidx_buf[widx] = marker_uidx;
	if (load_and_collapse(unfiltered_sample_ct, sample_ct, sample_exclude, final_mask, 0, bedfile, rawbuf, readbuf_cur)) {
	  goto calc_homozyg_ret_READ_FAIL;
	}
	mask_out_homozyg_major(readbuf_cur, sample_ct);
	increment_het_missing(readbuf_cur, sample_ct, het_cts, missing_cts, 1);
	widx = widx_next;
	marker_cidx++;
      }
      // now handle the last (window_size - 1) markers
      marker_cidx_max = marker_cidx + window_size;
      do {
	widx = widx_next;
        marker_cidx++;
	older_uidx = old_uidx;
	if (marker_cidx < marker_cidx_max) {
	  old_uidx = uidx_buf[widx];
	  readbuf_cur = &(readbuf[widx * sample_ctl2]);
	  swbuf_cur = &(swbuf[widx * sample_ctl]);
	  vertical_bitct_subtract(swbuf_cur, sample_ct, swhit_cts);
	  swhit_min = (int32_t)(((double)((int32_t)(marker_cidx_max - marker_cidx))) * hit_threshold + 1.0 - EPSILON);
	} else {
	  readbuf_cur = nullptr;
	}
	if (roh_update(hp, readbuf_cur, nullptr, het_cts, missing_cts, marker_pos, cur_sample_male, sample_ct, swhit_min, older_uidx, old_uidx, marker_cidx, max_roh_ct, swhit_cts, cur_roh_uidx_starts, cur_roh_cidx_starts, cur_roh_het_cts, cur_roh_missing_cts, sample_to_last_roh, roh_list, &roh_ct, marker_exclude, prev_roh_end_cidxs, end_nonhom_uidxs, cur_roh_earliest_extend_uidxs)) {
	  goto calc_homozyg_ret_NOMEM;
	}
	widx_next = widx + 1;
	if (widx_next == window_size) {
          widx_next = 0;
	}
	if (prev_roh_end_cidxs && (marker_cidx + hp->min_snp < marker_cidx_max)) {
	  // no need to update this once it's too late for a new ROH to begin
          update_end_nonhom(readbuf_cur, sample_ct, end_nonhom_uidxs, uidx_buf[widx_next]);
	}
      } while (marker_cidx <= marker_cidx_max);
    }
  }
  putc_unlocked('\r', stdout);
  LOGPRINTF("--homozyg: Scan complete, found %" PRIuPTR " ROH.\n", roh_ct);
  roh_list_chrom_starts[chrom_ct] = roh_ct;
  // "truncate" the completed list so we can start making workspace allocations
  // again
  bigstack_alloc(roh_ct * ROH_ENTRY_INTS * sizeof(int32_t)); // roh_list
  retval = write_main_roh_reports(outname, outname_end, marker_exclude, marker_ids, max_marker_id_len, plink_maxsnp, chrom_info_ptr, marker_pos, sample_ct, sample_exclude, sample_ids, plink_maxfid, plink_maxiid, max_sample_id_len, pheno_nm, pheno_c, pheno_d, missing_pheno_str, omp_is_numeric, missing_pheno_len, is_new_lengths, roh_ct, roh_list, roh_list_chrom_starts, sample_to_last_roh, &max_pool_size, &max_roh_len);
  if (retval) {
    goto calc_homozyg_ret_1;
  }
  *outname_end = '\0';
  LOGPRINTFWW("Results saved to %s.hom + %s.hom.indiv + %s.hom.summary .\n", outname, outname, outname);
  if (hp->modifier & (HOMOZYG_GROUP | HOMOZYG_GROUP_VERBOSE)) {
    if (max_pool_size < hp->pool_size_min) {
      LOGERRPRINTF("Warning: Skipping --homozyg group%s report since there are no pools.\n", (hp->modifier & HOMOZYG_GROUP_VERBOSE)? "-verbose" : "");
#ifndef __LP64__
    } else if (max_pool_size > 65536) {
      logerrprint("Error: 32-bit " PROG_NAME_STR "'s --homozyg group cannot handle a pool of size >65536.\n");
      goto calc_homozyg_ret_NOMEM;
#endif
    } else {
      if (omp_is_numeric) {
	scan_double(output_missing_pheno, &dxx);
	wptr = dtoa_g_wxp4(dxx, 8, missing_pheno_str);
	missing_pheno_len = (uintptr_t)(wptr - missing_pheno_str);
      }
      retval = roh_pool(hp, bedfile, bed_offset, outname, outname_end, rawbuf, marker_exclude, marker_ids, max_marker_id_len, plink_maxsnp, marker_allele_ptrs, max_marker_allele_len, marker_reverse, chrom_info_ptr, marker_pos, sample_ct, unfiltered_sample_ct, sample_exclude, sample_ids, plink_maxfid, plink_maxiid, max_sample_id_len, pheno_nm, pheno_c, pheno_d, missing_pheno_str, missing_pheno_len, is_new_lengths, roh_ct, roh_list, roh_list_chrom_starts, max_pool_size, max_roh_len);
      if (retval) {
	goto calc_homozyg_ret_1;
      }
    }
  }

  while (0) {
  calc_homozyg_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  calc_homozyg_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  }
 calc_homozyg_ret_1:
  bigstack_double_reset(bigstack_mark, bigstack_end_mark);
  return retval;
}
