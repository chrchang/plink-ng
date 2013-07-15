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

void mask_out_homozyg_major(uintptr_t* readbuf_cur, uint32_t indiv_ct) {
  // if readbuf_cur were 16-byte aligned, this could be vectorized, but it
  // isn't, and this isn't a limiting step anyway
  uintptr_t* readbuf_cur_end = &(readbuf_cur[(indiv_ct + (BITCT2 - 1)) / BITCT2]);
  uintptr_t cur_word;
  do {
    cur_word = *readbuf_cur;
    *readbuf_cur &= ((cur_word ^ (cur_word >> 1)) & FIVEMASK) * (3 * ONELU);
  } while (++readbuf_cur < readbuf_cur_end);
}

void increment_het_missing(uintptr_t* readbuf_cur, uint32_t indiv_ct, uint32_t* het_cts, uint32_t* missing_cts, int32_t incr) {
  // assumes homozyg majors masked out
  uint32_t indiv_idx_offset;
  uintptr_t cur_word;
  uint32_t last_set_bit;
  for (indiv_idx_offset = 0; indiv_idx_offset < indiv_ct; indiv_idx_offset += BITCT2) {
    cur_word = *readbuf_cur++;
    while (cur_word) {
      last_set_bit = CTZLU(cur_word);
      if (last_set_bit & 1) {
        het_cts[indiv_idx_offset + (last_set_bit / 2)] += incr;
      } else {
        missing_cts[indiv_idx_offset + (last_set_bit / 2)] += incr;
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

uint32_t roh_update(Homozyg_info* hp, uintptr_t* readbuf_cur, uintptr_t* swbuf_cur, uint32_t* het_cts, uint32_t* missing_cts, uint32_t* marker_pos, uintptr_t* cur_indiv_male, uint32_t indiv_ct, uint32_t swhit_min, uint32_t older_uidx, uint32_t old_uidx, uint32_t marker_cidx, uintptr_t max_roh_ct, uint32_t* swhit_cts, uint32_t* cur_roh_uidx_starts, uint32_t* cur_roh_cidx_starts, uint32_t* cur_roh_het_cts, uint32_t* cur_roh_missing_cts, uintptr_t* indiv_to_last_roh, uint32_t* roh_list, uintptr_t* roh_ct_ptr) {
  uint32_t min_snp = hp->min_snp;
  uint32_t min_bases = hp->min_bases;
  double max_bases_per_snp = hp->max_bases_per_snp;
  uint32_t max_hets = hp->max_hets;
  uint32_t max_sw_hets = hp->window_max_hets;
  uint32_t max_sw_missings = hp->window_max_missing;
  uint32_t is_new_lengths = 1 ^ ((hp->modifier / HOMOZYG_OLD_LENGTHS) & 1);
  uint32_t forced_end = (marker_pos[old_uidx] - marker_pos[older_uidx]) > hp->max_gap;
  uintptr_t roh_ct = *roh_ct_ptr;
  uint32_t is_cur_hit = 0;
  uintptr_t cur_word = 0;
  uint32_t* roh_list_cur = &(roh_list[roh_ct * ROH_ENTRY_INTS]);
  uint32_t indiv_idx;
  uintptr_t cur_call;
  uint32_t cidx_len;
  uint32_t uidx_first;
  uint32_t base_len;
  uintptr_t indiv_last_roh;
  for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
    if (!(indiv_idx & (BITCT2 - 1))) {
      if (readbuf_cur) {
        cur_word = *readbuf_cur++;
      }
    } else {
      cur_word >>= 2;
    }
    if ((!cur_indiv_male) || (!is_set(cur_indiv_male, indiv_idx))) {
      // not X-chromosome male
      if (readbuf_cur) {
        if (swbuf_cur) {
	  if ((het_cts[indiv_idx] <= max_sw_hets) && (missing_cts[indiv_idx] <= max_sw_missings)) {
	    set_bit_noct(swbuf_cur, indiv_idx);
	    swhit_cts[indiv_idx] += 1;
	  }
        }
        is_cur_hit = (swhit_cts[indiv_idx] >= swhit_min);
      }
      cur_call = cur_word & (3 * ONELU);
      if ((cur_roh_cidx_starts[indiv_idx] != 0xffffffffU) && ((!is_cur_hit) || ((cur_call == 2) && (cur_roh_het_cts[indiv_idx] == max_hets)) || forced_end)) {
	cidx_len = marker_cidx - cur_roh_cidx_starts[indiv_idx];
	if (cidx_len >= min_snp) {
	  uidx_first = cur_roh_uidx_starts[indiv_idx];
	  base_len = marker_pos[older_uidx] + is_new_lengths - marker_pos[uidx_first];
          if ((base_len >= min_bases) && (((double)((int32_t)cidx_len)) * max_bases_per_snp >= ((double)base_len))) {
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
	    *roh_list_cur++ = (uint32_t)indiv_last_roh;
	    *roh_list_cur++ = (uint32_t)(indiv_last_roh >> 32);
#else
	    *roh_list_cur++ = (uint32_t)indiv_last_roh;
#endif
	    indiv_to_last_roh[indiv_idx] = roh_ct++;
	  }
	}
	cur_roh_cidx_starts[indiv_idx] = 0xffffffffU;
      }
      if (is_cur_hit) {
	if (cur_roh_cidx_starts[indiv_idx] == 0xffffffffU) {
	  if ((!max_hets) && (cur_call == 2)) {
	    continue;
	  }
	  cur_roh_uidx_starts[indiv_idx] = old_uidx;
	  cur_roh_cidx_starts[indiv_idx] = marker_cidx;
	  cur_roh_het_cts[indiv_idx] = 0;
	  cur_roh_missing_cts[indiv_idx] = 0;
	}
	if (cur_call) {
	  if (cur_call == 2) {
	    cur_roh_het_cts[indiv_idx] += 1;
	  } else {
	    cur_roh_missing_cts[indiv_idx] += 1;
	  }
	}
      }
    }
  }
  *roh_ct_ptr = roh_ct;
  return 0;
}

int32_t write_main_roh_reports(char* outname, char* outname_end, uintptr_t* marker_exclude, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, Chrom_info* chrom_info_ptr, uint32_t* marker_pos, uintptr_t indiv_ct, uintptr_t* indiv_exclude, char* person_ids, uint32_t plink_maxfid, uint32_t plink_maxiid, uintptr_t max_person_id_len, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, char* missing_pheno_str, uint32_t missing_pheno_len, uint32_t is_new_lengths, uintptr_t roh_ct, uint32_t* roh_list, uintptr_t* roh_list_chrom_starts, uintptr_t* indiv_to_last_roh, uint32_t* max_pool_size_ptr, uint32_t* max_roh_len_ptr) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* outfile = NULL;
  FILE* outfile_indiv = NULL;
  char* wptr_iid = &(tbuf[plink_maxfid + 1]);
  char* wptr_phe = &(tbuf[plink_maxfid + plink_maxiid + 2]);
  int32_t* roh_ct_aff_adj = NULL;
  uint32_t max_pool_size = 0;
  uint32_t max_roh_len = 0;
  int32_t retval = 0;
  char* cptr;
  char* cptr2;
  char* wptr_chr;
  char* wptr_snp2;
  char* wptr_bp1;
  char* wptr_kb;
  char* wptr;
  uint32_t* cur_roh;
  int32_t* roh_ct_unaff_adj;
  uintptr_t indiv_uidx;
  uintptr_t indiv_idx;
  uintptr_t prev_roh_idx;
  uintptr_t cur_roh_idx;
  uintptr_t next_roh_idx;
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
  if (fopen_checked(&outfile, outname, "w")) {
    goto write_main_roh_reports_ret_OPEN_FAIL;
  }
  sprintf(tbuf, "%%%us %%%us      PHE  CHR %%%us %%%us         POS1         POS2         KB     NSNP  DENSITY     PHOM     PHET\n", plink_maxfid, plink_maxiid, plink_maxsnp, plink_maxsnp);
  if (fprintf(outfile, tbuf, "FID", "IID", "SNP1", "SNP2") < 0) {
    goto write_main_roh_reports_ret_WRITE_FAIL;
  }
  memcpy(&(outname_end[4]), ".indiv", 7);
  if (fopen_checked(&outfile_indiv, outname, "w")) {
    goto write_main_roh_reports_ret_OPEN_FAIL;
  }
  sprintf(tbuf, "%%%us %%%us  PHE     NSEG       KB    KBAVG\n", plink_maxfid, plink_maxiid);
  if (fprintf(outfile_indiv, tbuf, "FID", "IID") < 0) {
    goto write_main_roh_reports_ret_WRITE_FAIL;
  }
  tbuf[plink_maxfid] = ' ';
  tbuf[plink_maxfid + plink_maxiid + 1] = ' ';
  indiv_uidx = 0;
  for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
    indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx);
    cptr = &(person_ids[indiv_uidx * max_person_id_len]);
    cptr2 = (char*)memchr(cptr, '\t', max_person_id_len);
    slen = (uintptr_t)(cptr2 - cptr);
    memcpy(memseta(tbuf, 32, plink_maxfid - slen), cptr, slen);
    slen = strlen(++cptr2);
    memcpy(memseta(wptr_iid, 32, plink_maxiid - slen), cptr2, slen);
    if (!is_set(pheno_nm, indiv_uidx)) {
      wptr_chr = fw_strcpyn(8, missing_pheno_len, missing_pheno_str, wptr_phe);
    } else if (pheno_c) {
      wptr_chr = memseta(wptr_phe, 32, 7);
      *wptr_chr++ = '1' + is_set(pheno_c, indiv_uidx);
    } else {
      wptr_chr = width_force(8, wptr_phe, double_f_writew3(wptr_phe, pheno_d[indiv_uidx]));
    }
    wptr_chr = memseta(wptr_chr, 32, 3);
    wptr_chr[2] = ' ';
    wptr_chr[3 + plink_maxsnp] = ' ';
    wptr_snp2 = &(wptr_chr[4 + plink_maxsnp]);
    wptr_bp1 = memseta(&(wptr_snp2[plink_maxsnp]), 32, 3);
    memset(&(wptr_bp1[10]), 32, 3);
    wptr_bp1[23] = ' ';
    wptr_kb = &(wptr_bp1[24]);
    // traverse roh_list backwards, reversing the direction of [5], then
    // traverse forwards and reset [5] again to indiv_idx
    cur_roh_idx = indiv_to_last_roh[indiv_idx];
    next_roh_idx = ZEROLU;
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
      intprint2(wptr_chr, get_marker_chrom(chrom_info_ptr, marker_uidx1));
      cptr = &(marker_ids[marker_uidx1 * max_marker_id_len]);
      slen = strlen(cptr);
      memcpy(memseta(&(wptr_chr[3]), 32, plink_maxsnp - slen), cptr, slen);
      cptr = &(marker_ids[marker_uidx2 * max_marker_id_len]);
      slen = strlen(cptr);
      memcpy(memseta(wptr_snp2, 32, plink_maxsnp - slen), cptr, slen);
      uint32_writew10(wptr_bp1, marker_pos[marker_uidx1]);
      uint32_writew10(&(wptr_bp1[13]), marker_pos[marker_uidx2]);
      dxx = ((double)(marker_pos[marker_uidx2] + is_new_lengths - marker_pos[marker_uidx1])) / (1000.0 - EPSILON);
      kb_tot += dxx;
      wptr = width_force(10, wptr_kb, double_f_writew3(wptr_kb, dxx));
      *wptr++ = ' ';
      if (cur_roh[2] > max_roh_len) {
	max_roh_len = cur_roh[2];
      }
      wptr = uint32_writew8x(wptr, cur_roh[2], ' ');
      dyy = (1.0 + SMALLISH_EPSILON) / ((double)((int32_t)cur_roh[2]));
      wptr = width_force(8, wptr, double_f_writew3(wptr, dxx * dyy));
      // next two decimals guaranteed to be length 5
      wptr = memseta(wptr, 32, 4);
      wptr = double_f_writew3(wptr, ((double)((int32_t)cur_roh[3])) * dyy);
      wptr = memseta(wptr, 32, 4);
      wptr = double_f_writew3(wptr, ((double)((int32_t)cur_roh[4])) * dyy);
      *wptr++ = '\n';
      if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	goto write_main_roh_reports_ret_WRITE_FAIL;
      }
#ifdef __LP64__
      cur_roh_idx = ((uintptr_t)cur_roh[5]) | (((uintptr_t)cur_roh[6]) << 32);
#else
      cur_roh_idx = (uintptr_t)cur_roh[5];
#endif
      cur_roh[5] = indiv_uidx;
    }
    if (!is_set(pheno_nm, indiv_uidx)) {
      wptr = fw_strcpyn(4, missing_pheno_len - 4, missing_pheno_str, wptr_phe);
    } else if (pheno_c) {
      wptr = memseta(wptr_phe, 32, 3);
      *wptr++ = '1' + is_set(pheno_c, indiv_uidx);
    } else {
      wptr = width_force(4, wptr_phe, double_g_write(wptr_phe, pheno_d[indiv_uidx]));
    }
    *wptr++ = ' ';
    wptr = uint32_writew8x(wptr, cur_roh_ct, ' ');
    wptr = width_force(8, wptr, double_g_write(wptr, kb_tot));
    *wptr++ = ' ';
    if (cur_roh_ct) {
      kb_tot /= (double)((int32_t)cur_roh_ct);
    }
    wptr = width_force(8, wptr, double_g_write(wptr, kb_tot));
    wptr = memcpya(wptr, " \n", 2);
    if (fwrite_checked(tbuf, wptr - tbuf, outfile_indiv)) {
      goto write_main_roh_reports_ret_WRITE_FAIL;
    }
    indiv_uidx++;
  }
  if (fclose_null(&outfile)) {
    goto write_main_roh_reports_ret_WRITE_FAIL;
  }
  if (fclose_null(&outfile_indiv)) {
    goto write_main_roh_reports_ret_WRITE_FAIL;
  }
  memcpy(&(outname_end[5]), "summary", 8);
  if (fopen_checked(&outfile, outname, "w")) {
    goto write_main_roh_reports_ret_WRITE_FAIL;
  }
  sprintf(tbuf, " CHR %%%us           BP      AFF    UNAFF\n", plink_maxsnp);
  if (fprintf(outfile, tbuf, "SNP") < 0) {
    goto write_main_roh_reports_ret_WRITE_FAIL;
  }
  tbuf[1] = ' ';
  wptr_chr = &(tbuf[5]);
  memset(&(wptr_chr[plink_maxsnp]), 32, 3);
  wptr_bp1 = &(wptr_chr[plink_maxsnp + 3]);
  memset(&(wptr_bp1[10]), 32, 8);
  wptr_bp1[18] = '0';
  wptr_bp1[19] = ' ';
  for (chrom_fo_idx = 0; chrom_fo_idx < chrom_info_ptr->chrom_ct; chrom_fo_idx++) {
    chrom_roh_start = roh_list_chrom_starts[chrom_fo_idx];
    chrom_roh_ct = roh_list_chrom_starts[chrom_fo_idx + 1] - chrom_roh_start;
    uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
    chrom_start = chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx];
    chrom_len = chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx + 1] - chrom_start;
    wkspace_reset(wkspace_mark);
    if (wkspace_alloc_i_checked(&roh_ct_unaff_adj, (chrom_len + 1) * sizeof(int32_t))) {
      goto write_main_roh_reports_ret_NOMEM;
    }
    fill_int_zero(roh_ct_unaff_adj, chrom_len);
    if (pheno_c) {
      if (wkspace_alloc_i_checked(&roh_ct_aff_adj, (chrom_len + 1) * sizeof(int32_t))) {
        goto write_main_roh_reports_ret_NOMEM;
      }
      fill_int_zero(roh_ct_aff_adj, chrom_len);
    }
    cur_roh = &(roh_list[chrom_roh_start * ROH_ENTRY_INTS]);
    for (cur_roh_idx = 0; cur_roh_idx < chrom_roh_ct; cur_roh_idx++) {
      indiv_uidx = cur_roh[5];
      if ((!pheno_c) || (!is_set(pheno_c, indiv_uidx))) {
        roh_ct_unaff_adj[cur_roh[0] - chrom_start] += 1;
        roh_ct_unaff_adj[cur_roh[1] + 1 - chrom_start] -= 1;
      } else {
        roh_ct_aff_adj[cur_roh[0] - chrom_start] += 1;
        roh_ct_aff_adj[cur_roh[1] + 1 - chrom_start] -= 1;
      }
      cur_roh = &(cur_roh[ROH_ENTRY_INTS]);
    }
    intprint2(&(tbuf[2]), uii);
    chrom_len += chrom_start; // now chrom_end
    cur_roh_ct = 0; // unaff ct
    uii = 0; // aff ct
    for (marker_uidx1 = chrom_start; marker_uidx1 < chrom_len; marker_uidx1++) {
      if (is_set(marker_exclude, marker_uidx1)) {
	continue;
      }
      cptr = &(marker_ids[marker_uidx1 * max_marker_id_len]);
      slen = strlen(cptr);
      memcpy(memseta(wptr_chr, 32, plink_maxsnp - slen), cptr, slen);
      uint32_writew10(wptr_bp1, marker_pos[marker_uidx1]);
      if (!pheno_c) {
        wptr = &(wptr_bp1[20]);
      } else {
	uii += roh_ct_aff_adj[marker_uidx1 - chrom_start];
        wptr = uint32_writew8x(&(wptr_bp1[11]), uii, ' ');
      }
      cur_roh_ct += roh_ct_unaff_adj[marker_uidx1 - chrom_start];
      if (cur_roh_ct + uii > max_pool_size) {
        max_pool_size = cur_roh_ct + uii;
      }
      wptr = uint32_writew8x(wptr, cur_roh_ct, '\n');
      if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
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
  wkspace_reset(wkspace_mark);
  fclose_cond(outfile);
  fclose_cond(outfile_indiv);
  return retval;
}

// heap of 64-bit integers with maximum on top
// root at element 1, heap_size offset by 1
// slightly more complicated variant of this in wdist_cluster.c

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
    clear_bit_noct(roh_slot_occupied, (uint32_t)(cur_roh_heap[1]));
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

int32_t roh_pool(Homozyg_info* hp, FILE* bedfile, uint64_t bed_offset, char* outname, char* outname_end, uintptr_t* rawbuf, uintptr_t* marker_exclude, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, char* marker_alleles, uintptr_t max_marker_allele_len, Chrom_info* chrom_info_ptr, uint32_t* marker_pos, uintptr_t indiv_ct, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, char* person_ids, uint32_t plink_maxfid, uint32_t plink_maxiid, uintptr_t max_person_id_len, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, char* missing_pheno_str, uint32_t missing_pheno_len, uint32_t is_new_lengths, uintptr_t roh_ct, uint32_t* roh_list, uintptr_t* roh_list_chrom_starts, uint32_t max_pool_size, uint32_t max_roh_len) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* outfile = NULL;
  uint64_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  // uintptr_t unfiltered_indiv_ctl2 = (unfiltered_indiv_ct + (BITCT2 - 1)) / BITCT2;
  // uintptr_t cur_lookahead = 0;
  uint32_t is_verbose = hp->modifier & HOMOZYG_GROUP_VERBOSE;
  uint32_t max_pool_sizel = (max_pool_size + (BITCT - 1)) / BITCT;
  uint32_t pool_size_min = hp->pool_size_min;
  uint32_t pool_size_ct = max_pool_size + 1 - pool_size_min;
  uint32_t marker_uidx2 = 0;
  uint32_t fresh_meat = 0;
  uintptr_t pool_list_size = 0;
  uint32_t pool_ct = 0;
  int32_t retval = 0;
  uint32_t chrom_fo_idx_to_pidx[MAX_POSSIBLE_CHROM + 1]; // decreasing order
  unsigned char* wkspace_mark2;
  // uintptr_t* lookahead_buf;
  uintptr_t* pool_size_first_plidx;
  uint32_t* marker_uidx_to_cidx;
  uintptr_t* roh_slots;
  uintptr_t* roh_slot_occupied; // bitfield marking roh_slots occupancy
  uint32_t* indiv_uidx_sort_buf;
  uint64_t* cur_roh_heap; // high 32 bits = marker_uidx, low 32 bits = slot_idx
  uintptr_t* pool_list;
  uint32_t* cur_roh;
  uintptr_t* cur_pool;
  uint64_t* roh_slot_map; // high 32 bits = indiv_uidx, low 32 bits = slot_idx
  uint32_t* roh_slot_end_uidx; // tracks when to flush a roh_slot
  uint32_t* allelic_match; // counts, then group assignments
  uintptr_t* allelic_match_matrix; // pairwise match matrix, potentially huge
  uint32_t* uiptr;
  char* wptr;
  // uintptr_t max_lookahead;
  uintptr_t old_pool_list_size;
  uintptr_t pool_list_idx;
  uintptr_t roh_slot_wsize;
  uintptr_t max_pool_list_size;
  uintptr_t chrom_roh_start;
  uintptr_t roh_idx;
  // uint32_t lookahead_base_uidx;
  uint32_t pool_size;
  uint32_t chrom_fo_idx;
  uint32_t cur_roh_heap_top;
  uint32_t marker_uidx1;
  uint32_t marker_cidx;
  uint32_t slot_idx;
  uint32_t chrom_start;
  uint32_t chrom_len;
  uint32_t uii;
  logprint("Error: --homozyg group[-verbose] is under development.\n");
  retval = RET_CALC_NOT_YET_SUPPORTED;
  goto roh_pool_ret_1;
  uii = 0; // max chrom len
  for (chrom_fo_idx = 0; chrom_fo_idx < chrom_info_ptr->chrom_ct; chrom_fo_idx++) {
    if (roh_list_chrom_starts[chrom_fo_idx] == roh_list_chrom_starts[chrom_fo_idx + 1]) {
      continue;
    }
    chrom_len = chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx + 1] - chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx];
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
  if (wkspace_alloc_ul_checked(&pool_size_first_plidx, pool_size_ct * sizeof(intptr_t)) ||
      wkspace_alloc_ui_checked(&marker_uidx_to_cidx, uii * sizeof(int32_t)) ||
      wkspace_alloc_ul_checked(&roh_slots, max_pool_size * roh_slot_wsize * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&roh_slot_occupied, max_pool_sizel * sizeof(intptr_t)) ||
      wkspace_alloc_ull_checked(&roh_slot_map, max_pool_size * sizeof(int64_t)) ||
      wkspace_alloc_ui_checked(&roh_slot_end_uidx, max_pool_size * sizeof(int32_t)) ||
      wkspace_alloc_ui_checked(&allelic_match, max_pool_size * sizeof(int32_t)) ||
      wkspace_alloc_ul_checked(&allelic_match_matrix, (((uintptr_t)max_pool_size) * (max_pool_size - 1)) * (sizeof(intptr_t) / 2))) {
    goto roh_pool_ret_NOMEM;
  }
  wkspace_mark2 = wkspace_base;
  // roh_slot_map / roh_slot_end_uidx / allelic_match / allelic_match_matrix
  // not used at the same time as cur_roh_heap / indiv_uidx_sort_buf
  cur_roh_heap = &(roh_slot_map[-1]);
  wkspace_reset((unsigned char*)roh_slot_end_uidx);
  if (wkspace_alloc_ui_checked(&indiv_uidx_sort_buf, 3 * max_pool_size * sizeof(int32_t))) {
      goto roh_pool_ret_NOMEM;
  }
  if (wkspace_base < wkspace_mark2) {
    wkspace_base = wkspace_mark2;
  }

  fill_ulong_one(pool_size_first_plidx, pool_size_ct);
  fill_ulong_zero(roh_slot_occupied, max_pool_sizel);

  pool_list = (uintptr_t*)wkspace_base;
  max_pool_list_size = wkspace_left / sizeof(intptr_t);
  // Since our ROH are sorted by *last* SNP, it's easiest to scan for pools
  // from back to front if we wish to painlessly produce sorted lists.
  chrom_fo_idx = chrom_info_ptr->chrom_ct;

  do {
    chrom_fo_idx_to_pidx[chrom_fo_idx--] = pool_ct;
    chrom_fo_idx--;
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
	  slot_idx = next_non_set_unsafe(roh_slot_occupied, 0);
	  set_bit_noct(roh_slot_occupied, slot_idx);
	  // use roh_slots[0..(max_pool_size - 1)] to store references to
	  // active ROH here
	  roh_slots[slot_idx] = roh_idx;
          cur_roh_heap[cur_roh_heap_top++] = (((uint64_t)uii) << 32) | ((uint64_t)slot_idx);
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
	  pool_size = popcount_longs(roh_slot_occupied, 0, max_pool_sizel);
	  if (pool_size >= pool_size_min) {
	    // pool encoding:
	    // [0]: pool_list index of next pool of same size (~ZEROLU if none)
	    // 64-bit:
            //   [1]: pool size (P) in low 32 bits, 1-based pool idx in high
	    //   [2-(P+1)]: roh indexes, sorted by increasing indiv_idx
	    //   [(P+2)-(2P+1)]: allelic-match group assignment (32-bit) with
	    //                   31st bit set if reference member; NSIM count
	    //                   in high 32 bits
	    // 32-bit:
            //   [1]: P
	    //   [2]: 1-based pool idx
	    //   [3-(P+2)]: roh indexes
	    //   [(P+3)-(2P+2)]: allelic-match group assignment (31st bit set
	    //                   if reference), followed by NSIM count,
	    //                   interleaved
	    old_pool_list_size = pool_list_size;
#ifdef __LP64__
	    pool_list_size += 2 * pool_size + 2;
#else
            pool_list_size += 3 * pool_size + 3;
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
	    slot_idx = 0;
	    uiptr = indiv_uidx_sort_buf;
	    for (uii = 0; uii < pool_size; uii++) {
	      slot_idx = next_set_unsafe(roh_slot_occupied, slot_idx);
	      pool_list_idx = roh_slots[slot_idx]; // actually a ROH idx
	      *uiptr++ = roh_list[pool_list_idx * ROH_ENTRY_INTS + 5]; // indiv_uidx
              *uiptr++ = (uint32_t)pool_list_idx;
#ifdef __LP64__
	      *uiptr++ = (uint32_t)(pool_list_idx >> 32);
#endif
	      slot_idx++;
	    }
	    // sort in increasing indiv_uidx order, for reproducible results
            qsort(indiv_uidx_sort_buf, pool_size, 4 + sizeof(intptr_t), intcmp);
	    for (uii = 0; uii < pool_size; uii++) {
#ifdef __LP64__
              *cur_pool++ = ((uintptr_t)indiv_uidx_sort_buf[3 * uii + 1]) | (((uintptr_t)indiv_uidx_sort_buf[3 * uii + 2]) << 32);
#else
	      *cur_pool++ = indiv_uidx_sort_buf[2 * uii + 1];
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

  // Now we know how much memory the pools require, so we can assign the rest
  // to a lookahead buffer.
  pool_list = (uintptr_t*)wkspace_alloc(pool_list_size * sizeof(intptr_t));
  // max_lookahead = wkspace_left / (unfiltered_indiv_ctl2 * sizeof(intptr_t));
  // lookahead_buf = (uintptr_t*)wkspace_base;

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
  memcpy(&(outname_end[5]), "overlap.S", 9);
  pool_list_idx = pool_list_size;
  for (chrom_fo_idx = 0; chrom_fo_idx < chrom_info_ptr->chrom_ct; chrom_fo_idx++) {
    if (chrom_fo_idx_to_pidx[chrom_fo_idx] == chrom_fo_idx_to_pidx[chrom_fo_idx + 1]) {
      continue;
    }

    chrom_start = chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx];
    if (fseeko(bedfile, bed_offset + ((uint64_t)chrom_start) * unfiltered_indiv_ct4, SEEK_SET)) {
      goto roh_pool_ret_READ_FAIL;
    }
    marker_uidx2 = chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx + 1];
    marker_cidx = 0;
    for (marker_uidx1 = chrom_start; marker_uidx1 < marker_uidx2; marker_uidx1++) {
      if (is_set(marker_exclude, marker_uidx1)) {
	continue;
      }
      marker_uidx_to_cidx[marker_uidx1 - chrom_start] = marker_cidx;
    }
    fill_ulong_zero(roh_slot_occupied, max_pool_sizel);
    fill_ulong_one((uintptr_t*)roh_slot_map, max_pool_size * (sizeof(int64_t) / sizeof(intptr_t)));
    fill_uint_zero(roh_slot_end_uidx, max_pool_size);
    // lookahead_base_uidx = chrom_start;
    // cur_lookahead = 0;

    for (uii = chrom_fo_idx_to_pidx[chrom_fo_idx + 1]; uii < chrom_fo_idx_to_pidx[chrom_fo_idx]; uii++) {
      pool_list_idx = pool_list[pool_list_idx - 1];
      // todo:
      // 1. determine new roh_slots assignment, and which slots need to be
      //    populated
      // 2. populate them from disk and lookahead_buf
      // 3. populate allelic_match_matrix, store NSIM values in allelic_match[]
      // 4. greedily assign allelic match groups, save info
      // 5. handle is_verbose if necessary

      if (is_verbose) {
#ifdef __LP64__
	wptr = uint32_write(&(outname_end[14]), (uint32_t)(cur_pool[1] >> 32));
#else
	wptr = uint32_write(&(outname_end[14]), (uint32_t)cur_pool[2]);
#endif
	memcpy(wptr, ".verbose", 9);
	if (fopen_checked(&outfile, outname, "w")) {
	  goto roh_pool_ret_OPEN_FAIL;
	}

	// todo

	if (fclose_null(&outfile)) {
	  goto roh_pool_ret_WRITE_FAIL;
	}
      }
    }
  }

  outname_end[12] = '\0';
  if (fopen_checked(&outfile, outname, "w")) {
    goto roh_pool_ret_OPEN_FAIL;
  }
  sprintf(tbuf, " POOL %%%us %%%us      PHE  CHR %%%us %%%us            BP1            BP2       KB     NSNP NSIM    GRP\n", plink_maxfid, plink_maxiid, plink_maxsnp, plink_maxsnp);
  if (fprintf(outfile, tbuf, "FID", "IID", "SNP1", "SNP2") < 0) {
    goto roh_pool_ret_WRITE_FAIL;
  }
  for (; max_pool_size >= pool_size_min; max_pool_size--) {
    // todo
  }
  if (fclose_null(&outfile)) {
    goto roh_pool_ret_WRITE_FAIL;
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
  wkspace_reset(wkspace_mark);
  fclose_cond(outfile);
  return retval;
}

int32_t calc_homozyg(Homozyg_info* hp, FILE* bedfile, uintptr_t bed_offset, uint32_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, char* marker_alleles, uintptr_t max_marker_allele_len, Chrom_info* chrom_info_ptr, uint32_t* marker_pos, uintptr_t indiv_ct, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, char* person_ids, uint32_t plink_maxfid, uint32_t plink_maxiid, uintptr_t max_person_id_len, char* outname, char* outname_end, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, int32_t missing_pheno, uintptr_t* sex_male) {
  unsigned char* wkspace_mark = wkspace_base;
  uint64_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t unfiltered_indiv_ctl2 = (unfiltered_indiv_ct + (BITCT2 - 1)) / BITCT2;
  uintptr_t indiv_ctl = (indiv_ct + (BITCT - 1)) / BITCT;
  uintptr_t indiv_ctl2 = (indiv_ct + (BITCT2 - 1)) / BITCT2;
  uintptr_t window_size = hp->window_size;
  double hit_threshold = hp->hit_threshold;
  uint32_t is_new_lengths = 1 ^ ((hp->modifier / HOMOZYG_OLD_LENGTHS) & 1);
  uint32_t species = chrom_info_ptr->species;
  uint32_t chrom_ct = chrom_info_ptr->chrom_ct;
  int32_t x_code = species_x_code[species];
  uint64_t haploid_mask = species_haploid_mask[species];
  uintptr_t topsize = 0;
  uintptr_t roh_ct = 0;
  uintptr_t* indiv_male = NULL;
  uint32_t swhit_min = 0;
  int32_t retval = 0;
  uintptr_t roh_list_chrom_starts[MAX_POSSIBLE_CHROM + 1];
  char missing_pheno_str[15];
  uint32_t missing_pheno_len;
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
  // Note that this is sorted by primarily by LAST uidx, and secondarily by
  // indiv_idx.
  uint32_t* roh_list;

  uint32_t* uidx_buf; // circular buffer tracking most recent marker_uidxs
  uintptr_t max_roh_ct;
  uint32_t chrom_fo_idx;
  uintptr_t marker_uidx;
  uint32_t chrom_end;
  uintptr_t* readbuf_cur;
  uintptr_t* swbuf_cur;
  uintptr_t* cur_indiv_male;
  char* wptr;
  uintptr_t ulii;
  uintptr_t widx;
  uint32_t old_uidx;
  uint32_t older_uidx;
  uint32_t marker_cidx;
  uint32_t marker_cidx_max;
  uint32_t swbuf_full;
  uint32_t uii;
  uint32_t max_pool_size;
  uint32_t max_roh_len;
  wptr = int32_write(missing_pheno_str, missing_pheno);
  wptr = memcpya(wptr, ".000", 4);
  missing_pheno_len = (uintptr_t)(wptr - missing_pheno_str);

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
  // no other workspace allocations during main scan, so we can assign it all
  // to the ROH list
  max_roh_ct = ((wkspace_left - topsize) & (~(CACHELINE - 1))) / (ROH_ENTRY_INTS * sizeof(int32_t));
  roh_list = (uint32_t*)wkspace_base;
  ulii = indiv_ctl2 - 1;
  rawbuf[unfiltered_indiv_ctl2 - 1] = 0;
  for (widx = 0; widx < window_size; widx++) {
    readbuf[widx * indiv_ctl2 + ulii] = 0;
  }
  if ((x_code != -1) && ((chrom_info_ptr->chrom_mask >> ((uint32_t)x_code)) & 1)) {
    indiv_male = (uintptr_t*)top_alloc(&topsize, indiv_ctl * sizeof(intptr_t));
    if (!indiv_male) {
      goto calc_homozyg_ret_NOMEM;
    }
    collapse_copy_bitarr(indiv_ct, sex_male, indiv_exclude, popcount_longs_exclude(sex_male, indiv_exclude, indiv_ctl), indiv_male);
  }
  if (fseeko(bedfile, bed_offset, SEEK_SET)) {
    goto calc_homozyg_ret_READ_FAIL;
  }
  fill_ulong_one(indiv_to_last_roh, indiv_ct);

  for (chrom_fo_idx = 0; chrom_fo_idx < chrom_ct; chrom_fo_idx++) {
    uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
    roh_list_chrom_starts[chrom_fo_idx] = roh_ct;
    chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx + 1];
    if ((x_code == -1) || (uii != ((uint32_t)x_code))) {
      if ((haploid_mask >> uii) & 1) {
	marker_uidx = chrom_end;
	if (fseeko(bedfile, bed_offset + (uint64_t)marker_uidx * unfiltered_indiv_ct4, SEEK_SET)) {
	  goto calc_homozyg_ret_READ_FAIL;
	}
	continue;
      }
      cur_indiv_male = NULL;
    } else {
      cur_indiv_male = indiv_male;
    }
    printf("\r--homozyg: Scanning chromosome %u. \b", uii);
    fflush(stdout);
    marker_uidx = chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx];
    fill_ulong_zero(swbuf, indiv_ctl * window_size);
    fill_uint_zero(het_cts, indiv_ct);
    fill_uint_zero(missing_cts, indiv_ct);
    for (widx = 0; widx < window_size; widx++) {
      if (is_set(marker_exclude, marker_uidx)) {
        marker_uidx = next_non_set(marker_exclude, marker_uidx + 1, chrom_end);
        if (fseeko(bedfile, bed_offset + (uint64_t)marker_uidx * unfiltered_indiv_ct4, SEEK_SET)) {
	  goto calc_homozyg_ret_READ_FAIL;
	}
      }
      if (marker_uidx == chrom_end) {
	break;
      }
      readbuf_cur = &(readbuf[widx * indiv_ctl2]);
      if (load_and_collapse(bedfile, rawbuf, unfiltered_indiv_ct, readbuf_cur, indiv_ct, indiv_exclude)) {
	goto calc_homozyg_ret_READ_FAIL;
      }
      mask_out_homozyg_major(readbuf_cur, indiv_ct);
      increment_het_missing(readbuf_cur, indiv_ct, het_cts, missing_cts, 1);
      uidx_buf[widx] = marker_uidx++;
    }
    if (widx == window_size) {
      marker_uidx--;
      fill_ulong_zero(swbuf, window_size * indiv_ctl);
      fill_uint_zero(swhit_cts, indiv_ct);
      fill_uint_one(cur_roh_cidx_starts, indiv_ct);
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
	  vertical_bitct_subtract(swbuf_cur, indiv_ct, swhit_cts);
	  fill_ulong_zero(swbuf_cur, indiv_ctl);
	} else {
	  swhit_min = (int32_t)(((double)((int32_t)(widx + 1))) * hit_threshold + 1.0 - EPSILON);
	}
	if (roh_update(hp, readbuf_cur, swbuf_cur, het_cts, missing_cts, marker_pos, cur_indiv_male, indiv_ct, swhit_min, older_uidx, old_uidx, marker_cidx, max_roh_ct, swhit_cts, cur_roh_uidx_starts, cur_roh_cidx_starts, cur_roh_het_cts, cur_roh_missing_cts, indiv_to_last_roh, roh_list, &roh_ct)) {
	  goto calc_homozyg_ret_NOMEM;
	}
	increment_het_missing(readbuf_cur, indiv_ct, het_cts, missing_cts, -1);
	if (++marker_uidx == chrom_end) {
	  break;
	}
	if (is_set(marker_exclude, marker_uidx)) {
	  marker_uidx = next_non_set(marker_exclude, marker_uidx + 1, chrom_end);
	  if (fseeko(bedfile, bed_offset + (uint64_t)marker_uidx * unfiltered_indiv_ct4, SEEK_SET)) {
	    goto calc_homozyg_ret_READ_FAIL;
	  }
	  if (marker_uidx == chrom_end) {
	    break;
	  }
	}
	uidx_buf[widx] = marker_uidx;
	if (load_and_collapse(bedfile, rawbuf, unfiltered_indiv_ct, readbuf_cur, indiv_ct, indiv_exclude)) {
	  goto calc_homozyg_ret_READ_FAIL;
	}
	mask_out_homozyg_major(readbuf_cur, indiv_ct);
	increment_het_missing(readbuf_cur, indiv_ct, het_cts, missing_cts, 1);
	widx++;
	marker_cidx++;
	if (widx == window_size) {
	  widx = 0;
	  swbuf_full = 1;
	}
      }
      // now handle the last (window_size - 1) markers
      marker_cidx_max = marker_cidx + window_size;
      do {
	widx++;
	marker_cidx++;
	if (widx == window_size) {
	  widx = 0;
	}
	older_uidx = old_uidx;
	if (marker_cidx < marker_cidx_max) {
	  old_uidx = uidx_buf[widx];
	  readbuf_cur = &(readbuf[widx * indiv_ctl2]);
	  swbuf_cur = &(swbuf[widx * indiv_ctl]);
	  vertical_bitct_subtract(swbuf_cur, indiv_ct, swhit_cts);
	  swhit_min = (int32_t)(((double)((int32_t)(marker_cidx_max - marker_cidx))) * hit_threshold + 1.0 - EPSILON);
	} else {
	  readbuf_cur = NULL;
	}
	if (roh_update(hp, readbuf_cur, NULL, het_cts, missing_cts, marker_pos, cur_indiv_male, indiv_ct, swhit_min, older_uidx, old_uidx, marker_cidx, max_roh_ct, swhit_cts, cur_roh_uidx_starts, cur_roh_cidx_starts, cur_roh_het_cts, cur_roh_missing_cts, indiv_to_last_roh, roh_list, &roh_ct)) {
	  goto calc_homozyg_ret_NOMEM;
	}
      } while (marker_cidx <= marker_cidx_max);
    }
  }
  putchar('\r');
  sprintf(logbuf, "--homozyg: Scan complete, found %" PRIuPTR " ROH.\n", roh_ct);
  logprintb();
  roh_list_chrom_starts[chrom_ct] = roh_ct;
  // "truncate" the completed list so we can start making workspace allocations
  // again
  roh_list = (uint32_t*)wkspace_alloc(roh_ct * ROH_ENTRY_INTS * sizeof(int32_t));
  retval = write_main_roh_reports(outname, outname_end, marker_exclude, marker_ids, max_marker_id_len, plink_maxsnp, chrom_info_ptr, marker_pos, indiv_ct, indiv_exclude, person_ids, plink_maxfid, plink_maxiid, max_person_id_len, pheno_nm, pheno_c, pheno_d, missing_pheno_str, missing_pheno_len, is_new_lengths, roh_ct, roh_list, roh_list_chrom_starts, indiv_to_last_roh, &max_pool_size, &max_roh_len);
  if (retval) {
    goto calc_homozyg_ret_1;
  }
  *outname_end = '\0';
  sprintf(logbuf, "Results saved to %s.hom{,.indiv,.summary}.\n", outname);
  logprintb();
  if (hp->modifier & (HOMOZYG_GROUP | HOMOZYG_GROUP_VERBOSE)) {
    if (max_pool_size < hp->pool_size_min) {
      sprintf(logbuf, "Skipping --homozyg group%s report since there are no pools.\n", (hp->modifier & HOMOZYG_GROUP_VERBOSE)? "-verbose" : "");
      logprintb();
#ifndef __LP64__
    } else if (max_pool_size > 65536) {
      logprint("Error: 32-bit " PROG_NAME_STR "'s --homozyg group cannot handle a pool of size >65536.\n");
      goto calc_homozyg_ret_NOMEM;
#endif
    } else {
      wptr = double_g_writewx4(missing_pheno_str, (double)missing_pheno, 8);
      missing_pheno_len = (uintptr_t)(wptr - missing_pheno_str);
      retval = roh_pool(hp, bedfile, bed_offset, outname, outname_end, rawbuf, marker_exclude, marker_ids, max_marker_id_len, plink_maxsnp, marker_alleles, max_marker_allele_len, chrom_info_ptr, marker_pos, indiv_ct, unfiltered_indiv_ct, indiv_exclude, person_ids, plink_maxfid, plink_maxiid, max_person_id_len, pheno_nm, pheno_c, pheno_d, missing_pheno_str, missing_pheno_len, is_new_lengths, roh_ct, roh_list, roh_list_chrom_starts, max_pool_size, max_roh_len);
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
  wkspace_reset(wkspace_mark);
  return retval;
}
