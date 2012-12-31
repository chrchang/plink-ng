#include "wdist_common.h"

#define PHENO_EPSILON 0.000030517578125

// last prime before 2^23
#define HASHSIZE 8388593
#ifdef __LP64__
#define HASHMEM 67108752
#else
#define HASHMEM 33554384
#endif

typedef union {
  unsigned long long ullv;
  double dv;
} Lld_union;

typedef struct ll_entry_struct {
  Lld_union vu;
  struct ll_entry_struct* next;
  char idstr[0];
} Ll_entry;

static inline int idmatch2(char* idtab, char* id0, unsigned int id0_len_p1, char* id1, unsigned int id1_len_p1) {
  return (!(memcmp(idtab, id0, id0_len_p1) || memcmp(&(idtab[id0_len_p1]), id1, id1_len_p1)));
}

static inline unsigned int hashval(unsigned char* ss) {
  // just interpret as little-endian number and take modulo HASHSIZE
  unsigned int vv = 0;
  unsigned char cc = *ss;
  while (cc) {
    vv = ((vv << 8) + cc) % HASHSIZE;
    cc = *(++ss);
  }
  return vv;
}

static inline unsigned int hashval2(unsigned char* id1, unsigned char* id2) {
  unsigned int vv = 0;
  unsigned char cc = *id1;
  do {
    vv = ((vv << 8) + cc) % HASHSIZE;
    cc = *(++id1);
  } while (cc != 9);
  vv = ((vv << 8) + 9) % HASHSIZE;
  cc = *id2;
  do {
    vv = ((vv << 8) + cc) % HASHSIZE;
    cc = *(++id2);
  } while (cc != 9);
  return vv;
}

static inline Ll_entry* alloc_ll(unsigned long* topsize_ptr, unsigned int size) {
  unsigned long ts = *topsize_ptr + ((size + 15) & (~15LU));
  if (ts > wkspace_left) {
    return NULL;
  } else {
    *topsize_ptr = ts;
    return (Ll_entry*)(&(wkspace_base[wkspace_left - ts]));
  }
}

int merge_fam_id_scan(char* bedname, char* famname, unsigned int* max_person_id_len_ptr, unsigned int* max_person_full_len_ptr, int* is_dichot_pheno_ptr, Ll_entry** htable, unsigned long* topsize_ptr, unsigned long long* tot_indiv_ct_ptr, unsigned long* ped_buflen_ptr, unsigned int* cur_indiv_ct_ptr, int merge_mode) {
  unsigned int max_person_id_len = *max_person_id_len_ptr;
  unsigned int max_person_full_len = *max_person_full_len_ptr;
  int is_dichot_pheno = *is_dichot_pheno_ptr;
  unsigned long topsize = *topsize_ptr;
  unsigned long long tot_indiv_ct = *tot_indiv_ct_ptr;
  unsigned int cur_indiv_ct = 0;
  FILE* infile = NULL;
  unsigned int text_file = 0;
  unsigned int uii;
  unsigned int ujj;
  unsigned int ukk;
  unsigned int umm;
  unsigned int unn;
  unsigned int uoo;
  unsigned long ulii;
  Ll_entry** ll_pptr;
  Ll_entry* ll_ptr;
  char* bufptr;
  char* bufptr2;
  char* bufptr3;
  char* bufptr4;
  char* bufptr5;
  char* bufptr6;
  char* bufptr7;
  int retval;
  double pheno;
  if (!famname) {
    famname = bedname;
    text_file = 1;
  }
  if (fopen_checked(&infile, famname, "r")) {
    goto merge_fam_id_scan_ret_INVALID_FORMAT;
  }
  tbuf[MAXLINELEN - 1] = ' ';
  while (fgets(tbuf, MAXLINELEN, infile)) {
    if ((*tbuf > ' ') && (*tbuf != '#')) {
      bufptr6 = item_endnn(tbuf);
      uii = bufptr6 - tbuf;
      bufptr = skip_initial_spaces(bufptr6);
      bufptr7 = item_endnn(bufptr);
      ujj = bufptr7 - bufptr;
      bufptr2 = skip_initial_spaces(bufptr7);
      bufptr3 = item_endnn(bufptr2);
      unn = bufptr3 - bufptr2;
      bufptr3 = skip_initial_spaces(bufptr3);
      bufptr4 = item_endnn(bufptr3);
      uoo = bufptr4 - bufptr3;
      bufptr4 = skip_initial_spaces(bufptr4);
      bufptr5 = item_endnn(bufptr4);
      umm = bufptr5 - bufptr4;
      if (umm != 1) {
	*bufptr4 = '0';
      }
      bufptr5 = skip_initial_spaces(bufptr5);
      if (no_more_items(bufptr5)) {
	goto merge_fam_id_scan_ret_INVALID_FORMAT;
      }
      *bufptr6 = '\t';
      *bufptr7 = '\t';
      ukk = uii + ujj + 2;
      if (ukk > max_person_id_len) {
	max_person_id_len = ukk;
      }
      ukk += unn + uoo + 4;
      umm = hashval2((unsigned char*)tbuf, (unsigned char*)bufptr);
      ll_pptr = &(htable[umm]);
      ll_ptr = *ll_pptr;
      umm = 1;
      if (is_dichot_pheno) {
	is_dichot_pheno = eval_affection(bufptr5, -9, 2, 0);
      }
      if (sscanf(bufptr5, "%lg", &pheno) != 1) {
	pheno = -9;
      }
      while (ll_ptr) {
	if (idmatch2(ll_ptr->idstr, tbuf, uii + 1, bufptr, ujj + 1)) {
	  // add parental ID/sex check?  (PLINK doesn't have it)
	  umm = 0;
	  if (merge_mode == 1) {
	    if (fabs(pheno - ll_ptr->vu.dv) > PHENO_EPSILON) {
	      ll_ptr->vu.dv = -9;
	    }
	  } else if (merge_mode == 2) {
	    if (ll_ptr->vu.dv == -9) {
	      ll_ptr->vu.dv = pheno;
	    }
	  } else if ((merge_mode == 5) || ((merge_mode == 3) && (pheno != -9))) {
	    ll_ptr->vu.dv = pheno;
	  }
	  break;
	}
        ll_pptr = &(ll_ptr->next);
	ll_ptr = *ll_pptr;
      }
      if (umm) {
	if (ukk > max_person_full_len) {
	  max_person_full_len = ukk;
	}
	ll_ptr = alloc_ll(&topsize, ukk + sizeof(long) + 8);
	ll_ptr->next = NULL;
	ll_ptr->vu.dv = pheno;
	memcpy(ll_ptr->idstr, tbuf, uii);
	ll_ptr->idstr[uii] = '\t';
	memcpy(&(ll_ptr->idstr[uii + 1]), bufptr, ujj);
	ukk = uii + ujj + 1;
	ll_ptr->idstr[ukk++] = '\t';
	memcpy(&(ll_ptr->idstr[ukk]), bufptr2, unn);
        ll_ptr->idstr[ukk + unn] = '\t';
	ukk += unn + 1;
	memcpy(&(ll_ptr->idstr[ukk]), bufptr3, uoo);
	ll_ptr->idstr[ukk + uoo] = '\t';
	ukk += uoo + 1;
	memcpy(&(ll_ptr->idstr[ukk]), bufptr4, 1);
	ll_ptr->idstr[ukk + 1] = '\0';
	*ll_pptr = ll_ptr;
	tot_indiv_ct++;
      }
      cur_indiv_ct++;
    }
    if (!tbuf[MAXLINELEN - 1]) {
      if (!text_file) {
	logprint("Error: .fam line too long.\n");
	goto merge_fam_id_scan_ret_INVALID_FORMAT;
      }
      ulii = 0;
      do {
	tbuf[MAXLINELEN - 1] = ' ';
	if (tbuf[MAXLINELEN - 2] == '\n') {
	  break;
	}
	ulii += MAXLINELEN - 1;
#ifndef __LP64__
	if (ulii > 2147483647) {
	  logprint("Error: .ped line too long (>= 2GB) for 32-bit WDIST.\n");
	  goto merge_fam_id_scan_ret_INVALID_FORMAT;
	}
#endif
        if (!fgets(tbuf, MAXLINELEN, infile)) {
	  goto merge_fam_id_scan_ret_READ_FAIL;
	}
      } while (!tbuf[MAXLINELEN - 1]);
      ulii += strlen(tbuf) + 1;
      if (ulii > (*ped_buflen_ptr)) {
	*ped_buflen_ptr = ulii;
      }
    }
  }
  if (!feof(infile)) {
    goto merge_fam_id_scan_ret_READ_FAIL;
  }
  *max_person_id_len_ptr = max_person_id_len;
  *max_person_full_len_ptr = max_person_full_len;
  *is_dichot_pheno_ptr = is_dichot_pheno;
  *topsize_ptr = topsize;
  *tot_indiv_ct_ptr = tot_indiv_ct;
  *cur_indiv_ct_ptr = cur_indiv_ct;
  retval = 0;
  while (0) {
  merge_fam_id_scan_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  merge_fam_id_scan_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
  }
  fclose_cond(infile);
  return retval;
}

int merge_datasets(char* bedname, char* bimname, char* famname, char* outname, char* outname_end, char* mergename1, char* mergename2, char* mergename3, int calculation_type, int merge_type) {
  FILE* mergelistfile = NULL;
  FILE* bed_outfile = NULL;
  FILE* bim_outfile = NULL;
  FILE* fam_outfile = NULL;
  unsigned char* wkspace_mark = wkspace_base;
  unsigned int max_person_id_len = 0;
  unsigned int max_person_full_len = 0;
  unsigned int max_marker_id_len = 0;
  int is_dichot_pheno = 1;
  int merge_mode = (merge_type & MERGE_MODE_MASK);
  Ll_entry** htable = (Ll_entry**)(&(wkspace_base[wkspace_left - HASHMEM]));
  unsigned long ped_buflen = MAXLINELEN;
  char* pheno_c = NULL;
  double* pheno_d = NULL;
  unsigned long topsize;
  char* person_ids;
  char* person_fids;
  char* marker_ids;
  unsigned long mlpos;
  unsigned long merge_ct;
  char* mergelist_buf;
  char** mergelist_bed;
  char** mergelist_bim;
  char** mergelist_fam;
  unsigned int* cur_indiv_cts;
  unsigned int tot_indiv_ct;
  unsigned int tot_marker_ct;
  unsigned long long ullxx;
  unsigned long ulii;
  unsigned int uii;
  char* bufptr;
  char* bufptr2;
  char* bufptr3;
  char* bufptr4;
  Ll_entry* ll_ptr;
  int retval = 0;
  char cc;

  if (merge_type & MERGE_LIST) {
    if (fopen_checked(&mergelistfile, mergename1, "r")) {
      goto merge_datasets_ret_READ_FAIL;
    }
    merge_ct = 1;
    ullxx = 0;
    // first pass: determine merge_ct, mergelist_buf size, verify no lines have
    // > 3 entries
    while (fgets(tbuf, MAXLINELEN, mergelistfile)) {
      bufptr = skip_initial_spaces(tbuf);
      if (no_more_items(bufptr)) {
	continue;
      }
      bufptr2 = next_item_mult(bufptr, 3);
      if (!no_more_items(bufptr2)) {
	logprint("Error: More than three items on --merge-list line.\n");
        goto merge_datasets_ret_INVALID_FORMAT;
      }
      if (no_more_items(next_item(bufptr))) {
	bufptr2 = item_endnn(bufptr);
	ulii = bufptr2 - bufptr;
	if (ulii > FNAMESIZE - 5) {
	  logprint("Error: Excessively long fileset prefix in --merge-list file.\n");
	  goto merge_datasets_ret_INVALID_FORMAT;
	}
	ullxx += 3 * ulii + 15;
      } else {
	do {
	  bufptr2 = item_endnn(bufptr);
	  ulii = bufptr2 - bufptr;
	  if (ulii > FNAMESIZE - 1) {
	    logprint("Error: Excessively long filename in --merge-list file.\n");
	    goto merge_datasets_ret_INVALID_FORMAT;
	  }
	  ullxx += ulii + 1;
	  bufptr = skip_initial_spaces(bufptr2);
	} while (!no_more_items(bufptr));
      }
      merge_ct++;
    }
    if (!feof(mergelistfile)) {
      goto merge_datasets_ret_READ_FAIL;
    }
    if (!merge_ct) {
      logprint("Error: Empty --merge-list file.\n");
      goto merge_datasets_ret_INVALID_FORMAT;
    }
#ifndef __LP64__
    if (ullxx > 2147483647) {
      goto merge_datasets_ret_NOMEM;
    }
#endif
    mergelist_bed = (char**)wkspace_alloc(merge_ct * sizeof(long));
    mergelist_bim = (char**)wkspace_alloc(merge_ct * sizeof(long));
    mergelist_fam = (char**)wkspace_alloc(merge_ct * sizeof(long));
    if (wkspace_alloc_c_checked(&mergelist_buf, (unsigned long)ullxx) ||
        wkspace_alloc_ui_checked(&cur_indiv_cts, merge_ct * sizeof(int))) {
      goto merge_datasets_ret_NOMEM;
    }
    rewind(mergelistfile);
    bufptr4 = mergelist_buf;
    mlpos = 1;
    while (fgets(tbuf, MAXLINELEN, mergelistfile)) {
      bufptr = skip_initial_spaces(tbuf);
      if (no_more_items(bufptr)) {
	continue;
      }
      bufptr2 = item_endnn(bufptr);
      ulii = (bufptr2 - bufptr);
      bufptr3 = skip_initial_spaces(bufptr2);
      if (no_more_items(bufptr3)) {
        memcpy(bufptr4, bufptr, ulii);
	memcpy(&(bufptr4[ulii]), ".bed", 5);
	mergelist_bed[mlpos] = bufptr4;
	bufptr4 = &(bufptr4[ulii + 5]);
        memcpy(bufptr4, bufptr, ulii);
	memcpy(&(bufptr4[ulii]), ".bim", 5);
	mergelist_bim[mlpos] = bufptr4;
	bufptr4 = &(bufptr4[ulii + 5]);
        memcpy(bufptr4, bufptr, ulii);
	memcpy(&(bufptr4[ulii]), ".fam", 5);
	mergelist_fam[mlpos] = bufptr4;
	bufptr4 = &(bufptr4[ulii + 5]);
      } else {
	memcpy(bufptr4, bufptr, ulii);
	bufptr4[ulii] = '\0';
	mergelist_bed[mlpos] = bufptr4;
	bufptr4 = &(bufptr4[ulii + 1]);
	bufptr2 = item_endnn(bufptr3);
	ulii = bufptr2 - bufptr3;
	bufptr = skip_initial_spaces(bufptr2);
	memcpy(bufptr4, bufptr3, ulii);
	bufptr4[ulii] = '\0';
	mergelist_bim[mlpos] = bufptr4;
	bufptr4 = &(bufptr4[ulii + 1]);
	if (no_more_items(bufptr)) {
	  mergelist_fam[mlpos] = NULL;
	} else {
	  bufptr2 = item_endnn(bufptr);
	  ulii = bufptr2 - bufptr;
	  memcpy(bufptr4, bufptr, ulii);
	  bufptr4[ulii] = '\0';
	  mergelist_fam[mlpos] = bufptr4;
	  bufptr4 = &(bufptr4[ulii + 1]);
	}
      }
      if (++mlpos == merge_ct) {
	break;
      }
    }
    if ((mlpos < merge_ct) && (!feof(mergelistfile))) {
      goto merge_datasets_ret_READ_FAIL;
    }
    fclose_null(&mergelistfile);
  } else {
    merge_ct = 2;
    mergelist_bed = (char**)wkspace_alloc(2 * sizeof(long));
    mergelist_bim = (char**)wkspace_alloc(2 * sizeof(long));
    mergelist_fam = (char**)wkspace_alloc(2 * sizeof(long));
    if (wkspace_alloc_ui_checked(&cur_indiv_cts, 2 * sizeof(int))) {
      goto merge_datasets_ret_NOMEM;
    }
    mergelist_bed[1] = mergename1;
    mergelist_bim[1] = mergename2;
    mergelist_fam[1] = (merge_type & MERGE_BINARY)? mergename3 : NULL;
  }
  mergelist_bed[0] = bedname;
  mergelist_bim[0] = bimname;
  mergelist_fam[0] = (famname[0])? famname : NULL;

  // ID counting/duplicate detection strategy:
  // - We do NOT want to scan through .ped files any more times than absolutely
  // necessary.  So we actually use *gasp* a hash table here.
  // - The hash table is positioned at the FAR end of wkspace, automatically
  // sized to ~64MB (or ~32MB on 32-bit systems).  IDs are then stored
  // backwards from there.  This simplifies copying into a sorted list.
  if (wkspace_left < HASHSIZE * sizeof(long)) {
    goto merge_datasets_ret_NOMEM;
  }
  for (uii = 0; uii < HASHSIZE; uii++) {
    htable[uii] = NULL;
  }
  topsize = HASHMEM;

  ullxx = 0;
  for (mlpos = 0; mlpos < merge_ct; mlpos++) {
    retval = merge_fam_id_scan(mergelist_bed[mlpos], mergelist_fam[mlpos], &max_person_id_len, &max_person_full_len, &is_dichot_pheno, htable, &topsize, &ullxx, &ped_buflen, &(cur_indiv_cts[mlpos]), merge_mode);
    if (retval) {
      goto merge_datasets_ret_1;
    }
  }
#ifdef __LP64__
  if (ullxx > 2147483647) {
    logprint("Error: Too many individuals (max 2147483647).\n");
    goto merge_datasets_ret_INVALID_FORMAT;
  }
#else
  if (ullxx > 268435455) {
    logprint("Error: Too many individuals for 32-bit WDIST (max 268435455).\n");
    goto merge_datasets_ret_INVALID_FORMAT;
  }
#endif
  tot_indiv_ct = ullxx;
  if (wkspace_alloc_c_checked(&person_ids, max_person_id_len * tot_indiv_ct) ||
      wkspace_alloc_c_checked(&person_fids, max_person_full_len * tot_indiv_ct)) {
    goto merge_datasets_ret_NOMEM;
  }
  if (is_dichot_pheno) {
    if (wkspace_left < topsize + tot_indiv_ct) {
      goto merge_datasets_ret_NOMEM;
    }
    wkspace_alloc_c_checked(&pheno_c, tot_indiv_ct);
  } else {
    if (wkspace_left < topsize + tot_indiv_ct * sizeof(double)) {
      goto merge_datasets_ret_NOMEM;
    }
    wkspace_alloc_d_checked(&pheno_d, tot_indiv_ct * sizeof(double));
  }
  ulii = 0;
  for (uii = 0; uii < HASHSIZE; uii++) {
    if (htable[uii]) {
      ll_ptr = htable[uii];
      htable[uii] = NULL;
      do {
	strcpy(&(person_fids[ulii * max_person_full_len]), ll_ptr->idstr);
	if (is_dichot_pheno) {
	  if (ll_ptr->vu.dv == -9) {
            pheno_c[ulii] = -1;
	  } else if (ll_ptr->vu.dv == 1) {
	    pheno_c[ulii] = 0;
	  } else {
	    pheno_c[ulii] = 1;
	  }
	} else {
	  pheno_d[ulii] = ll_ptr->vu.dv;
	}
	ulii++;
	ll_ptr = ll_ptr->next;
      } while (ll_ptr);
    }
  }
  wkspace_left -= topsize;
  if (is_dichot_pheno) {
    if (qsort_ext(person_fids, tot_indiv_ct, max_person_full_len, strcmp_deref, pheno_c, 1)) {
      wkspace_left += topsize;
      goto merge_datasets_ret_NOMEM;
    }
  } else {
    if (qsort_ext(person_fids, tot_indiv_ct, max_person_full_len, strcmp_deref, (char*)pheno_d, sizeof(double))) {
      wkspace_left += topsize;
      goto merge_datasets_ret_NOMEM;
    }
  }
  wkspace_left += topsize;
  memcpy(outname_end, ".fam", 5);
  if (fopen_checked(&fam_outfile, outname, "w")) {
    goto merge_datasets_ret_OPEN_FAIL;
  }
  bufptr = person_fids;
  bufptr3 = person_ids;
  for (ulii = 0; ulii < tot_indiv_ct; ulii++) {
    bufptr2 = next_item_mult(bufptr, 2);
    uii = (bufptr2 - bufptr) - 1;
    memcpy(bufptr3, bufptr, uii);
    bufptr3[uii] = '\0';
    bufptr3 = &(bufptr3[max_person_id_len]);
    uii += strlen(bufptr2) + 1;
    if (fwrite_checked(bufptr, uii, fam_outfile)) {
      goto merge_datasets_ret_WRITE_FAIL;
    }
    bufptr = &(bufptr[max_person_full_len]);
    if (is_dichot_pheno) {
      cc = pheno_c[ulii];
      if (fprintf(fam_outfile, "\t%s\n", cc? ((cc == 1)? "2" : "-9") : "1") < 0) {
	goto merge_datasets_ret_WRITE_FAIL;
      }
    } else {
      if (fprintf(fam_outfile, "\t%g\n", pheno_d[ulii]) < 0) {
	goto merge_datasets_ret_WRITE_FAIL;
      }
    }
  }
  if (fclose_null(&fam_outfile)) {
    goto merge_datasets_ret_WRITE_FAIL;
  }
  wkspace_reset((unsigned char*)person_fids);
  topsize = HASHMEM;

  while (0) {
  merge_datasets_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  merge_datasets_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  merge_datasets_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  merge_datasets_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  merge_datasets_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 merge_datasets_ret_1:
  fclose_cond(mergelistfile);
  fclose_cond(bed_outfile);
  fclose_cond(bim_outfile);
  fclose_cond(fam_outfile);
  wkspace_reset(wkspace_mark);
  return retval;
}
