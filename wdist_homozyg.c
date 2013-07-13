#include "wdist_homozyg.h"

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
  homozyg_ptr->segment_match_snp = 20;
  homozyg_ptr->pool_size_min = 2;
}

void homozyg_update_readbuf(uintptr_t* rawbuf, uint32_t cur_indiv_ct, uintptr_t* cur_indiv_exclude, uintptr_t* readbuf_cur, uint32_t* het_cts, uint32_t* missing_cts) {
  uintptr_t cur_write = 0;
  uint32_t indiv_uidx = 0;
  uint32_t indiv_idx_low = 0;
  uintptr_t cur_read;
  uint32_t indiv_idx;
  for (indiv_idx = 0; indiv_idx < cur_indiv_ct; indiv_idx++) {
    indiv_uidx = next_non_set_unsafe(cur_indiv_exclude, indiv_uidx);
    cur_read = (rawbuf[indiv_uidx / BITCT2] >> (2 * (indiv_uidx % BITCT2))) & (3 * ONELU);
    if (cur_read == 2) {
      het_cts[indiv_idx] += 1;
    } else if (cur_read == 1) {
      missing_cts[indiv_idx] += 1;
    }
    cur_write |= cur_read << (indiv_idx_low * 2);
    if (++indiv_idx_low == BITCT2) {
      *readbuf_cur++ = cur_write;
      cur_write = 0;
    }
    indiv_uidx++;
  }
  if (indiv_idx_low) {
    *readbuf_cur = cur_write;
  }
}

void homozyg_scroll_out(uintptr_t* readbuf_cur, uint32_t cur_indiv_ct, uint32_t* het_cts, uint32_t* missing_cts) {
  uint32_t indiv_idx_offset;
  uintptr_t cur_word;
  uint32_t last_set_bit;
  // oh, why not, may as well make this faster than the original add
  for (indiv_idx_offset = 0; indiv_idx_offset < cur_indiv_ct; indiv_idx_offset += BITCT2) {
    cur_word = *readbuf_cur++;
    // mask out all homozygous majors
    cur_word &= ((cur_word ^ (cur_word >> 1)) & FIVEMASK) * (3 * ONELU);

    while (cur_word) {
      last_set_bit = CTZLU(cur_word);
      if (last_set_bit & 1) {
        het_cts[indiv_idx_offset + (last_set_bit / 2)] -= 1;
      } else {
        missing_cts[indiv_idx_offset + (last_set_bit / 2)] -= 1;
      }
      cur_word &= cur_word - ONELU;
    }
  }
}

#ifdef __LP64__
#define ROH_ENTRY_INTS 7
#else
#define ROH_ENTRY_INTS 6
#endif

uint32_t roh_update(Homozyg_info* hp, uintptr_t* readbuf_cur, uintptr_t* swbuf_cur, uint32_t* het_cts, uint32_t* missing_cts, uint32_t* marker_pos, uint32_t indiv_ct, uint32_t swhit_min, uint32_t older_uidx, uint32_t old_uidx, uint32_t marker_cidx, uintptr_t max_roh_ct, uint32_t* swhit_cts, uint32_t* cur_roh_uidx_starts, uint32_t* cur_roh_cidx_starts, uint32_t* cur_roh_het_cts, uint32_t* cur_roh_missing_cts, uintptr_t* indiv_to_last_roh, uint32_t* roh_list, uintptr_t* roh_ct_ptr) {
  uint32_t min_snp = hp->min_snp;
  uint32_t min_bases = hp->min_bases;
  double max_bases_per_snp = hp->max_bases_per_snp;
  uint32_t max_hets = hp->max_hets;
  uint32_t max_sw_hets = hp->window_max_hets;
  uint32_t max_sw_missings = hp->window_max_missing;
  uint32_t is_new_lengths = 1 ^ ((hp->modifier / HOMOZYG_OLD_LENGTHS) & 1);
  uint32_t forced_end = (marker_pos[old_uidx] - marker_pos[older_uidx]) > hp->max_gap;
  uintptr_t roh_ct = *roh_ct_ptr;
  uint32_t* roh_list_cur = &(roh_list[roh_ct * ROH_ENTRY_INTS]);
  uintptr_t cur_word = *readbuf_cur++;
  uint32_t indiv_idx;
  uint32_t is_cur_hit;
  uint32_t cidx_len;
  uint32_t uidx_first;
  uint32_t base_len;
  uintptr_t indiv_last_roh;
  uintptr_t cur_call;
  for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
    if ((het_cts[indiv_idx] <= max_sw_hets) && (missing_cts[indiv_idx] <= max_sw_missings)) {
      set_bit_noct(swbuf_cur, indiv_idx);
      swhit_cts[indiv_idx] += 1;
    }
    cur_call = cur_word & (3 * ONELU);
    is_cur_hit = (swhit_cts[indiv_idx] >= swhit_min);
    if ((cur_roh_cidx_starts[indiv_idx] != 0xffffffffU) && (forced_end || (!is_cur_hit) || ((cur_call == 2) && (cur_roh_het_cts[indiv_idx] == max_hets)))) {
      cidx_len = marker_cidx - cur_roh_cidx_starts[indiv_idx];
      uidx_first = cur_roh_uidx_starts[indiv_idx];
      base_len = marker_pos[older_uidx] + is_new_lengths - marker_pos[uidx_first];
      if ((cidx_len >= min_snp) && (base_len >= min_bases) && (((double)((int32_t)cidx_len)) * max_bases_per_snp < (double)((int32_t)base_len))) {
	if (roh_ct == max_roh_ct) {
	  return 1;
	}
        *roh_list_cur++ = uidx_first;
	*roh_list_cur++ = older_uidx;
	*roh_list_cur++ = cidx_len;
        *roh_list_cur++ = cidx_len - cur_roh_het_cts[indiv_idx] - cur_roh_missing_cts[indiv_idx];
        *roh_list_cur++ = cur_roh_het_cts[indiv_idx];
	indiv_last_roh = indiv_to_last_roh[indiv_idx];
#ifdef __LP64__
        *roh_list_cur++ = (uint32_t)(indiv_last_roh & 0xffffffffU);
	*roh_list_cur++ = (uint32_t)(indiv_last_roh >> 32);
#else
        *roh_list_cur++ = (uint32_t)indiv_last_roh;
#endif
        indiv_to_last_roh[indiv_idx] = roh_ct++;
      }
      cur_roh_cidx_starts[indiv_idx] = 0xffffffffU;
    }
    if (is_cur_hit) {
      if (cur_roh_cidx_starts[indiv_idx] == 0xffffffffU) {
	cur_roh_uidx_starts[indiv_idx] = old_uidx;
	cur_roh_cidx_starts[indiv_idx] = marker_cidx;
	cur_roh_het_cts[indiv_idx] = 0;
	cur_roh_missing_cts[indiv_idx] = 0;
      }
      if (cur_call == 2) {
	if (!max_hets) {
	  // if the first would-be ROH call is actually a het, and
	  // max_hets == 0, don't start the ROH at all
	  cur_roh_cidx_starts[indiv_idx] = 0xffffffffU;
	} else {
	  cur_roh_het_cts[indiv_idx] += 1;
	}
      } else if (cur_call == 1) {
	cur_roh_missing_cts[indiv_idx] += 1;
      }
    }
    if (indiv_idx & (BITCT2 - 1)) {
      cur_word = *readbuf_cur++;
    } else {
      cur_word >>= 2;
    }
  }
  *roh_ct_ptr = roh_ct;
  return 0;
}

int32_t calc_homozyg(Homozyg_info* hp, FILE* bedfile, uintptr_t bed_offset, uint32_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, char* marker_ids, uintptr_t max_marker_id_len, Chrom_info* chrom_info_ptr, uint32_t* marker_pos, uintptr_t indiv_ct, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, char* person_ids, uint32_t plink_maxfid, uint32_t plink_maxiid, uintptr_t max_person_id_len, char* outname, char* outname_end, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, uintptr_t* sex_male) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* outfile = NULL;
  uint64_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t unfiltered_indiv_ctl2 = (unfiltered_indiv_ct + (BITCT2 - 1)) / BITCT2;
  uintptr_t indiv_ctl = (indiv_ct + (BITCT - 1)) / BITCT;
  uintptr_t indiv_ctl2 = (indiv_ct + (BITCT2 - 1)) / BITCT2;
  uintptr_t window_size = hp->window_size;
  double hit_threshold = hp->hit_threshold;
  // uint32_t is_new_lengths = 1 & ((hp->modifier / HOMOZYG_OLD_LENGTHS) & 1);
  uint32_t species = chrom_info_ptr->species;
  uint32_t chrom_ct = chrom_info_ptr->chrom_ct;
  int32_t x_code = species_x_code[species];
  uint64_t haploid_mask = species_haploid_mask[species];
  uintptr_t topsize = 0;
  uint32_t nonmale_ct = 0;
  uintptr_t roh_ct = 0;
  uintptr_t* indiv_exclude_male = NULL;
  uint32_t swhit_min = 0;
  int32_t retval = 0;
  uintptr_t* rawbuf;
  uintptr_t* readbuf; // circular window of actual genotype data
  uintptr_t* swbuf; // circular window of recent scanning window hits
  uint32_t* het_cts;
  uint32_t* missing_cts;
  uint32_t* swhit_cts;
  uint32_t* cur_roh_uidx_starts;
  uint32_t* cur_roh_cidx_starts;
  uint32_t* cur_roh_het_cts;
  uint32_t* cur_roh_missing_cts;
  uintptr_t* indiv_to_last_roh;

  // 6 uint32_ts per entry in 32-bit build, 7 on 64-bit
  // [0]: first uidx
  // [1]: last uidx
  // [2]: nsnp
  // [3]: hom ct
  // [4]: het ct
  // [5] (and usually [6]): uintptr_t indicating position of individual's
  // previous ROH, or 0xfff... if there is none; becomes indiv_idx after
  // .hom.indiv file written
  uint32_t* roh_list;

  uint32_t* uidx_buf; // circular buffer tracking most recent marker_uidxs
  uintptr_t max_roh_ct;
  uint32_t chrom_fo_idx;
  uintptr_t marker_uidx;
  uint32_t chrom_end;
  uint32_t cur_indiv_ct;
  uintptr_t* cur_indiv_exclude;
  uintptr_t* readbuf_cur;
  uintptr_t* swbuf_cur;
  uintptr_t ulii;
  uintptr_t widx;
  uint32_t old_uidx;
  uint32_t older_uidx;
  uint32_t marker_cidx;
  // uint32_t indiv_idx;
  uint32_t swbuf_full;
  uint32_t uii;
  if (wkspace_alloc_ul_checked(&rawbuf, unfiltered_indiv_ctl2 * sizeof(intptr_t))) {
    goto calc_homozyg_ret_NOMEM;
  }

  readbuf = (uintptr_t*)top_alloc(&topsize, indiv_ctl2 * window_size * sizeof(intptr_t));
  if (!readbuf) {
    goto calc_homozyg_ret_NOMEM;
  }
  swbuf = (uintptr_t*)top_alloc(&topsize, indiv_ctl * window_size * sizeof(intptr_t));
  if (!swbuf) {
    goto calc_homozyg_ret_NOMEM;
  }
  het_cts = (uint32_t*)top_alloc(&topsize, indiv_ct * sizeof(int32_t));
  if (!het_cts) {
    goto calc_homozyg_ret_NOMEM;
  }
  missing_cts = (uint32_t*)top_alloc(&topsize, indiv_ct * sizeof(int32_t));
  if (!missing_cts) {
    goto calc_homozyg_ret_NOMEM;
  }
  swhit_cts = (uint32_t*)top_alloc(&topsize, indiv_ct * sizeof(int32_t));
  if (!swhit_cts) {
    goto calc_homozyg_ret_NOMEM;
  }
  cur_roh_uidx_starts = (uint32_t*)top_alloc(&topsize, indiv_ct * sizeof(int32_t));
  if (!cur_roh_uidx_starts) {
    goto calc_homozyg_ret_NOMEM;
  }
  cur_roh_cidx_starts = (uint32_t*)top_alloc(&topsize, indiv_ct * sizeof(int32_t));
  if (!cur_roh_cidx_starts) {
    goto calc_homozyg_ret_NOMEM;
  }
  cur_roh_het_cts = (uint32_t*)top_alloc(&topsize, indiv_ct * sizeof(int32_t));
  if (!cur_roh_het_cts) {
    goto calc_homozyg_ret_NOMEM;
  }
  cur_roh_missing_cts = (uint32_t*)top_alloc(&topsize, indiv_ct * sizeof(int32_t));
  if (!cur_roh_missing_cts) {
    goto calc_homozyg_ret_NOMEM;
  }
  indiv_to_last_roh = (uintptr_t*)top_alloc(&topsize, indiv_ct * sizeof(intptr_t));
  if (!indiv_to_last_roh) {
    goto calc_homozyg_ret_NOMEM;
  }
  uidx_buf = (uint32_t*)top_alloc(&topsize, window_size * sizeof(int32_t));
  if (!uidx_buf) {
    goto calc_homozyg_ret_NOMEM;
  }
  max_roh_ct = ((wkspace_left - topsize) & (~(CACHELINE - 1))) / (ROH_ENTRY_INTS * sizeof(int32_t));
  roh_list = (uint32_t*)wkspace_alloc(max_roh_ct * ROH_ENTRY_INTS * sizeof(int32_t));
  ulii = indiv_ctl2 - 1;
  for (widx = 0; widx < window_size; widx++) {
    readbuf[widx * indiv_ctl2 + ulii] = 0;
  }
  if ((x_code != -1) && ((chrom_info_ptr->chrom_mask >> ((uint32_t)x_code)) & 1)) {
    indiv_exclude_male = (uintptr_t*)top_alloc(&topsize, indiv_ctl * sizeof(intptr_t));
    if (!indiv_exclude_male) {
      goto calc_homozyg_ret_NOMEM;
    }
    for (ulii = 0; ulii < indiv_ctl; ulii++) {
      indiv_exclude_male[ulii] = indiv_exclude[ulii] | sex_male[ulii];
    }
    nonmale_ct = unfiltered_indiv_ct - popcount_longs(indiv_exclude_male, 0, indiv_ctl);
  }
  if (fseeko(bedfile, bed_offset, SEEK_SET)) {
    goto calc_homozyg_ret_READ_FAIL;
  }
  fill_ulong_one(indiv_to_last_roh, indiv_ct);
  for (chrom_fo_idx = 0; chrom_fo_idx < chrom_ct; chrom_fo_idx++) {
    uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
    chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx + 1];
    if ((x_code == -1) || (uii != ((uint32_t)x_code))) {
      if ((haploid_mask >> uii) & 1) {
	marker_uidx = chrom_end;
	if (fseeko(bedfile, bed_offset + marker_uidx * unfiltered_indiv_ct4, SEEK_SET)) {
	  goto calc_homozyg_ret_READ_FAIL;
	}
	continue;
      }
      cur_indiv_ct = indiv_ct;
      cur_indiv_exclude = indiv_exclude;
    } else {
      cur_indiv_ct = nonmale_ct;
      cur_indiv_exclude = indiv_exclude_male;
    }
    marker_uidx = chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx];
    fill_ulong_zero(swbuf, indiv_ctl * window_size);
    fill_uint_zero(het_cts, cur_indiv_ct);
    fill_uint_zero(missing_cts, cur_indiv_ct);
    for (widx = 0; widx < window_size; widx++) {
      if (is_set(marker_exclude, marker_uidx)) {
        marker_uidx = next_non_set(marker_exclude, marker_uidx + 1, chrom_end);
        if (fseeko(bedfile, bed_offset + marker_uidx * unfiltered_indiv_ct4, SEEK_SET)) {
	  goto calc_homozyg_ret_READ_FAIL;
	}
      }
      if (marker_uidx == chrom_end) {
	break;
      }
      if (fread(rawbuf, 1, unfiltered_indiv_ct4, bedfile)) {
	goto calc_homozyg_ret_READ_FAIL;
      }
      homozyg_update_readbuf(rawbuf, cur_indiv_ct, cur_indiv_exclude, &(readbuf[widx * indiv_ctl2]), het_cts, missing_cts);
      uidx_buf[widx] = marker_uidx++;
    }
    if (widx == window_size) {
      fill_uint_zero(swhit_cts, cur_indiv_ct);
      fill_uint_one(cur_roh_cidx_starts, cur_indiv_ct);
      widx = 0;
      swbuf_full = 0;
      marker_cidx = 0;
      old_uidx = uidx_buf[0];
      // repurpose widx as next readbuf row
      while (1) {
	older_uidx = old_uidx;
	old_uidx = uidx_buf[widx];
	readbuf_cur = &(readbuf[widx * indiv_ctl2]);
	swbuf_cur = &(swbuf[widx * indiv_ctl]);
	if (swbuf_full) {
	  vertical_bitct_subtract(swbuf_cur, cur_indiv_ct, swhit_cts);
	  fill_ulong_zero(swbuf_cur, indiv_ctl);
	} else {
	  swhit_min = (int32_t)(((double)((int32_t)(widx + 1))) * hit_threshold + 1.0 - EPSILON);
	}
	if (roh_update(hp, readbuf_cur, swbuf_cur, het_cts, missing_cts, marker_pos, cur_indiv_ct, swhit_min, older_uidx, old_uidx, marker_cidx, max_roh_ct, swhit_cts, cur_roh_uidx_starts, cur_roh_cidx_starts, cur_roh_het_cts, cur_roh_missing_cts, indiv_to_last_roh, roh_list, &roh_ct)) {
	  goto calc_homozyg_ret_NOMEM;
	}

	uidx_buf[widx] = marker_uidx++;
	if (marker_uidx == chrom_end) {
	  break;
	}
	if (is_set(marker_exclude, marker_uidx)) {
	  marker_uidx = next_non_set(marker_exclude, marker_uidx + 1, chrom_end);
	  if (fseeko(bedfile, bed_offset + marker_uidx * unfiltered_indiv_ct4, SEEK_SET)) {
	    goto calc_homozyg_ret_READ_FAIL;
	  }
	  if (marker_uidx == chrom_end) {
	    break;
	  }
	}
	if (fread(rawbuf, 1, unfiltered_indiv_ct4, bedfile)) {
	  goto calc_homozyg_ret_READ_FAIL;
	}
	homozyg_scroll_out(readbuf_cur, cur_indiv_ct, het_cts, missing_cts);
        homozyg_update_readbuf(rawbuf, cur_indiv_ct, cur_indiv_exclude, readbuf_cur, het_cts, missing_cts);
	widx++;
	marker_cidx++;
	if (widx == window_size) {
	  widx = 0;
	  swbuf_full = 1;
	}
      }
      // now handle the last (window_size - 1) markers
      /*
      for () {
      }
      */
    }
  }
  logprint("Error: --homozyg is currently under development.\n");
  retval = RET_CALC_NOT_YET_SUPPORTED;
  goto calc_homozyg_ret_1;
  while (0) {
  calc_homozyg_ret_NOMEM:
    retval = RET_NOMEM;
    break;
    /*
  calc_homozyg_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
    */
  calc_homozyg_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
    /*
  calc_homozyg_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
    */
  }
 calc_homozyg_ret_1:
  wkspace_reset(wkspace_mark);
  fclose_cond(outfile);
  return retval;
}
