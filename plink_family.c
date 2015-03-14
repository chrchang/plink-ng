#include "plink_common.h"

#include "plink_assoc.h"
#include "plink_cluster.h"
#include "plink_family.h"
#include "plink_stats.h"

void family_init(Family_info* fam_ip) {
  fam_ip->mendel_max_trio_error = 1.0;
  fam_ip->mendel_max_var_error = 1.0;
  fam_ip->mendel_exclude_one_ratio = 0.0;
  fam_ip->mendel_modifier = 0;
  fam_ip->tdt_modifier = 0;
  fam_ip->tdt_mperm_val = 0;
  fam_ip->dfam_modifier = 0;
  fam_ip->dfam_mperm_val = 0;
  fam_ip->qfam_modifier = 0;
  fam_ip->qfam_mperm_val = 0;
}

// bottom 2 bits of index = child genotype
// middle 2 bits of index = paternal genotype
// top 2 bits of index = maternal genotype

// bits 0-7 = child increment (always 1)
// bits 8-15 = father increment
// bits 16-23 = mother increment
// bits 24-31 = error code
const uint32_t mendel_error_table[] =
{0, 0, 0x1010101, 0x8000001,
 0, 0, 0, 0x7010001,
 0, 0, 0, 0x7010001,
 0x3000101, 0, 0, 0x7010001,
 0, 0, 0, 0x6000101,
 0, 0, 0, 0,
 0, 0, 0, 0,
 0x3000101, 0, 0, 0,
 0, 0, 0, 0x6000101,
 0, 0, 0, 0,
 0, 0, 0, 0,
 0x3000101, 0, 0, 0,
 0x4010001, 0, 0, 0x6000101,
 0x4010001, 0, 0, 0,
 0x4010001, 0, 0, 0,
 0x5000001, 0, 0x2010101, 0};
// necessary to check child gender when dealing with error 9/10
const uint32_t mendel_error_table_x[] =
{0, 0, 0x1010101, 0x8000001,
 0, 0, 0, 0x9010001,
 0, 0, 0, 0x7010001,
 0x3000101, 0, 0, 0x7010001,
 0, 0, 0, 0x6000101,
 0, 0, 0, 0,
 0, 0, 0, 0,
 0x3000101, 0, 0, 0,
 0, 0, 0, 0x6000101,
 0, 0, 0, 0,
 0, 0, 0, 0,
 0x3000101, 0, 0, 0,
 0x4010001, 0, 0, 0x6000101,
 0xa010001, 0, 0, 0,
 0x4010001, 0, 0, 0,
 0x5000001, 0, 0x2010101, 0};

int32_t get_trios_and_families(uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, uintptr_t* founder_info, uintptr_t* sex_nm, uintptr_t* sex_male, char* sample_ids, uintptr_t max_sample_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, char** fids_ptr, uintptr_t* max_fid_len_ptr, char** iids_ptr, uintptr_t* max_iid_len_ptr, uint64_t** family_list_ptr, uint32_t* family_ct_ptr, uint64_t** trio_list_ptr, uintptr_t* trio_ct_ptr, uint32_t** trio_lookup_ptr, uint32_t include_duos, uint32_t toposort) {
  // This mirrors linkRelateds() in genedrop.cpp, and parseTrios() in trio.cpp,
  // in PLINK 1.07.
  //
  // family_list has paternal indices in low 32 bits, maternal indices in high
  // 32, sorted in child ID order.
  // trio_list has child IDs in low 32 bits, family_list indices in high 32
  // bits, in PLINK 1.07 reporting order.
  // If include_duos is set, missing parental IDs are coded as
  // unfiltered_sample_ct.  include_duos must be 0 or 1.
  // If toposort is set, trio_lookup has child IDs in [4n], paternal IDs in
  // [4n+1], maternal IDs in [4n+2], and reporting sequence index in [4n+3].
  // Otherwise, it's [3n]/[3n+1]/[3n+2].  toposort must be 0 or 1.
  //
  // fids is a list of null-terminated FIDs using trio_list indices, and iids
  // is a list of IIDs using regular unfiltered indices.  If include_duos is
  // set, iids has a trailing entry set to '0'.  (fids_ptr, iids_ptr, and the
  // corresponding lengths can be NULL.)
  //
  // PLINK 1.07 enforces <= 1 father and <= 1 mother per sample (and ambiguous
  // sex parents are not permitted), but the IDs CAN be reversed in the .fam
  // with no adverse consequences.  For backward compatibility, we replicate
  // this.  (Possible todo: report a warning exactly once when this happens.)
  // It won't be replicated in PLINK 2.0.
  unsigned char* wkspace_mark = wkspace_base;
  uint64_t* edge_list = NULL;
  uint32_t* toposort_queue = NULL;
  char* fids = NULL;
  char* iids = NULL;
  uintptr_t unfiltered_sample_ctl = (unfiltered_sample_ct + (BITCT - 1)) / BITCT;
  uintptr_t unfiltered_sample_ctp1l = 1 + (unfiltered_sample_ct / BITCT);
  uintptr_t sample_uidx = next_unset_unsafe(sample_exclude, 0);
  // does *not* use populate_id_htable
  uintptr_t htable_size = geqprime(2 * unfiltered_sample_ct + 1);
  uintptr_t topsize = 0;
  uintptr_t max_fid_len = 2;
  uintptr_t max_iid_len = 2;
  uint64_t family_code = 0;
  uint32_t family_ct = 0;
  uint32_t edge_ct = 0;
  uint32_t remaining_edge_ct = 0;
  uint32_t tqueue_start = 0;
  uint32_t tqueue_end = 0;
  int32_t retval = 0;
  uintptr_t* founder_info2; // now decoupled from --make-founders
  uint64_t* family_list;
  uint64_t* trio_list_tmp;

  // Initial "hash value" is twice the lower of the two parental uidxs.  We use
  // an open addressing scheme with prime table size >2n and a form of
  // quadratic probing to ensure performance stays reasonable if someone has
  // waaay too many concubines, etc.  Specifically, given an initial hash
  // collision at position 2j, and larger parent ID k (set to n if second
  // parent is actually missing), next probe is (k+j) positions later, 3rd
  // probe after that is 3(k+j) positions after the second, etc. (all numbers
  // modulo table size).
  uint64_t* family_htable;

  uint64_t* trio_write;
  uint64_t* edge_write;
  uint32_t* family_idxs;
  uint32_t* sample_id_map;
  uint32_t* trio_lookup;
  uint32_t* uiptr;
  char* sorted_sample_ids;
  char* idbuf;
  char* idptr;
  char* iidptr;
  char* pidptr1;
  char* pidptr2;
  uint64_t trio_code;
  uint64_t edge_code;
  uint64_t ullii;
  uintptr_t topsize_bak;
  uintptr_t topsize_bak2;
  uintptr_t family_idx;
  uintptr_t trio_ct;
  uintptr_t trio_idx;
  uintptr_t fidlen;
  uintptr_t sample_idx;
  uintptr_t slen;
  uintptr_t uidx1;
  uintptr_t uidx2;
  uintptr_t next_probe_incr;
  uintptr_t probe_incr_incr;
  uint32_t first_sex;
  uint32_t uii;
  int32_t sorted_idx;
  founder_info2 = (uintptr_t*)top_alloc(&topsize, unfiltered_sample_ctp1l * sizeof(intptr_t));
  if (!founder_info2) {
    goto get_trios_and_families_ret_NOMEM;
  }
  memcpy(founder_info2, founder_info, unfiltered_sample_ctl * sizeof(intptr_t));
  if (unfiltered_sample_ct & (BITCT - 1)) {
    SET_BIT(founder_info2, unfiltered_sample_ct);
  } else {
    founder_info2[unfiltered_sample_ctl] = 1;
  }
  topsize_bak = topsize;
  trio_list_tmp = (uint64_t*)top_alloc(&topsize, sample_ct * sizeof(int64_t));
  if (!trio_list_tmp) {
    goto get_trios_and_families_ret_NOMEM;
  }
  topsize_bak2 = topsize;
  sorted_sample_ids = (char*)top_alloc(&topsize, sample_ct * max_sample_id_len);
  if (!sorted_sample_ids) {
    goto get_trios_and_families_ret_NOMEM;
  }
  sample_id_map = (uint32_t*)top_alloc(&topsize, sample_ct * sizeof(int32_t));
  if (!sample_id_map) {
    goto get_trios_and_families_ret_NOMEM;
  }
  idbuf = (char*)top_alloc(&topsize, max_sample_id_len);
  if (!idbuf) {
    goto get_trios_and_families_ret_NOMEM;
  }
  wkspace_left -= topsize;
  if (sort_item_ids_noalloc(sorted_sample_ids, sample_id_map, unfiltered_sample_ct, sample_exclude, sample_ct, sample_ids, max_sample_id_len, 0, 0, strcmp_deref)) {
    wkspace_left += topsize;
    goto get_trios_and_families_ret_1;
  }
  // over-allocate here, we shrink family_list later when we know how many
  // families there are
  if (wkspace_alloc_ull_checked(&family_list, sample_ct * sizeof(int64_t)) ||
      wkspace_alloc_ull_checked(&family_htable, htable_size * sizeof(int64_t)) ||
      wkspace_alloc_ui_checked(&family_idxs, htable_size * sizeof(int32_t))) {
    goto get_trios_and_families_ret_NOMEM2;
  }
  fill_ull_zero(family_htable, htable_size);
  wkspace_left += topsize;
  // 1. populate family_list (while using family_htable to track duplicates),
  //    determine max_iid_len, count qualifying trios
  trio_write = trio_list_tmp;
  for (sample_idx = 0; sample_idx < sample_ct; sample_uidx++, sample_idx++) {
    next_unset_ul_unsafe_ck(sample_exclude, &sample_uidx);
    idptr = &(sample_ids[sample_uidx * max_sample_id_len]);
    iidptr = (char*)memchr(idptr, '\t', max_sample_id_len);
    fidlen = (uintptr_t)((++iidptr) - idptr);
    if (fidlen > max_fid_len) {
      max_fid_len = fidlen;
    }
    slen = strlen(iidptr);
    if (slen >= max_iid_len) {
      max_iid_len = slen + 1;
    }
    if (IS_SET(founder_info, sample_uidx)) {
      continue;
    }
    memcpy(idbuf, idptr, fidlen);
    pidptr1 = &(paternal_ids[sample_uidx * max_paternal_id_len]);
    pidptr2 = &(maternal_ids[sample_uidx * max_maternal_id_len]);
    slen = strlen(pidptr1);
    if (fidlen + slen < max_sample_id_len) {
      memcpy(&(idbuf[fidlen]), pidptr1, slen + 1);
      sorted_idx = bsearch_str(idbuf, fidlen + slen, sorted_sample_ids, max_sample_id_len, sample_ct);
    } else {
      sorted_idx = -1;
    }
    first_sex = 0;
    if (sorted_idx == -1) {
      if (!include_duos) {
	SET_BIT(founder_info2, sample_uidx);
	continue;
      }
      uidx1 = unfiltered_sample_ct;
    } else {
      uidx1 = sample_id_map[(uint32_t)sorted_idx];
      if (uidx1 == sample_uidx) {
        idbuf[fidlen - 1] = ' ';
	LOGPREPRINTFWW("Error: '%s' is his/her own parent.\n", idbuf);
	goto get_trios_and_families_ret_INVALID_FORMAT_2;
      } else if (!IS_SET(sex_nm, uidx1)) {
        idbuf[fidlen - 1] = ' ';
	LOGPREPRINTFWW("Error: Parent '%s' has unspecified sex.\n", idbuf);
	goto get_trios_and_families_ret_INVALID_FORMAT_2;
      }
      first_sex = 2 - IS_SET(sex_male, uidx1);
    }
    slen = strlen(pidptr2);
    if (fidlen + slen < max_sample_id_len) {
      memcpy(&(idbuf[fidlen]), pidptr2, slen);
      sorted_idx = bsearch_str(idbuf, fidlen + slen, sorted_sample_ids, max_sample_id_len, sample_ct);
    } else {
      sorted_idx = -1;
    }
    if (sorted_idx == -1) {
      if ((!include_duos) || (uidx1 == unfiltered_sample_ct)) {
	SET_BIT(founder_info2, sample_uidx);
	continue;
      }
      if (first_sex == 1) {
	family_code = (((uint64_t)unfiltered_sample_ct) << 32) | ((uint64_t)uidx1);
      } else {
	family_code = (((uint64_t)uidx1) << 32) | ((uint64_t)unfiltered_sample_ct);
      }
      next_probe_incr = unfiltered_sample_ct + uidx1;
    } else {
      uidx2 = sample_id_map[(uint32_t)sorted_idx];
      if (uidx2 == sample_uidx) {
        idbuf[fidlen - 1] = ' ';
	LOGPREPRINTFWW("Error: '%s' is their own parent.\n", idbuf);
	goto get_trios_and_families_ret_INVALID_FORMAT_2;
      } else if (!IS_SET(sex_nm, uidx2)) {
        idbuf[fidlen - 1] = ' ';
	LOGPREPRINTFWW("Error: Parent '%s' has unspecified sex.\n", idbuf);
	goto get_trios_and_families_ret_INVALID_FORMAT_2;
      }
      uii = IS_SET(sex_male, uidx2);
      if (2 - uii == first_sex) {
	idptr[fidlen - 1] = ' ';
	LOGPREPRINTFWW("Error: '%s' has two %sies.\n", idptr, (first_sex == 1)? "dadd" : "momm");
	goto get_trios_and_families_ret_INVALID_FORMAT_2;
      }
      if (uii) {
        family_code = (((uint64_t)uidx1) << 32) | ((uint64_t)uidx2);
      } else {
        family_code = (((uint64_t)uidx2) << 32) | ((uint64_t)uidx1);
      }
      next_probe_incr = uidx1 + uidx2;
      if (uidx1 > uidx2) {
	uidx1 = uidx2;
      }
    }
    probe_incr_incr = next_probe_incr * 2;
    uidx1 *= 2; // now current probe position
    while (1) {
      ullii = family_htable[uidx1];
      if (!ullii) {
	family_idx = family_ct++;
	family_list[family_idx] = family_code;
        family_htable[uidx1] = family_code;
	family_idxs[uidx1] = family_idx;
	break;
      } else if (ullii == family_code) {
	family_idx = family_idxs[uidx1];
	break;
      }
      uidx1 += next_probe_incr;
      if (uidx1 >= htable_size) {
        uidx1 -= htable_size;
      }
      // quadratic probing implementation that avoids 64-bit modulus operation
      // on 32-bit systems
      next_probe_incr += probe_incr_incr;
      if (next_probe_incr >= htable_size) {
	next_probe_incr -= htable_size;
      }
    }
    *trio_write++ = (((uint64_t)family_idx) << 32) | ((uint64_t)sample_uidx);
  }
  trio_ct = (uintptr_t)(trio_write - trio_list_tmp);
  wkspace_reset(wkspace_mark);
  wkspace_alloc(family_ct * sizeof(int64_t)); // family_list
  topsize = topsize_bak2;
  wkspace_left -= topsize;
  if (wkspace_alloc_ull_checked(&trio_write, trio_ct * sizeof(int64_t))) {
    goto get_trios_and_families_ret_NOMEM2;
  }
  wkspace_left += topsize;
  topsize = topsize_bak;
  memcpy(trio_write, trio_list_tmp, trio_ct * sizeof(int64_t));
#ifdef __cplusplus
  std::sort((int64_t*)trio_write, (int64_t*)(&(trio_write[trio_ct])));
#else
  qsort(trio_write, trio_ct, sizeof(int64_t), llcmp);
#endif
  wkspace_left -= topsize;
  if (wkspace_alloc_ui_checked(&trio_lookup, trio_ct * (3 + toposort) * sizeof(int32_t))) {
    goto get_trios_and_families_ret_NOMEM2;
  }
  if (fids_ptr) {
    if (wkspace_alloc_c_checked(&fids, trio_ct * max_fid_len) ||
	wkspace_alloc_c_checked(&iids, (unfiltered_sample_ct + include_duos) * max_iid_len)) {
      goto get_trios_and_families_ret_NOMEM2;
    }
  }
  if (toposort) {
    if (wkspace_alloc_ull_checked(&edge_list, trio_ct * 2 * sizeof(int64_t))) {
      goto get_trios_and_families_ret_NOMEM2;
    }
    // Edge list excludes founder parents; edge codes have parental uidx in
    // high 32 bits, trio idx in low 32.
    edge_write = edge_list;
    for (trio_idx = 0; trio_idx < trio_ct; trio_idx++) {
      family_code = family_list[(uintptr_t)(trio_write[trio_idx] >> 32)];
      if (!is_set(founder_info2, family_code >> 32)) {
	*edge_write++ = (family_code & 0xffffffff00000000LLU) | ((uint64_t)trio_idx);
      }
      if (!is_set(founder_info2, (uint32_t)family_code)) {
        *edge_write++ = ((family_code & 0xffffffffLLU) << 32) | ((uint64_t)trio_idx);
      }
    }
    edge_ct = (uintptr_t)(edge_write - edge_list);
    if (edge_ct) {
#ifdef __cplusplus
      std::sort((int64_t*)edge_list, (int64_t*)(&(edge_list[edge_ct])));
#else
      qsort(edge_list, edge_ct, sizeof(int64_t), llcmp);
#endif
    }
    wkspace_shrink_top(edge_list, edge_ct * sizeof(int64_t));
    if (wkspace_alloc_ui_checked(&toposort_queue, trio_ct * sizeof(int32_t))) {
      goto get_trios_and_families_ret_NOMEM2;
    }
    remaining_edge_ct = edge_ct;
  }
  wkspace_left += topsize;
  *family_list_ptr = family_list;
  *family_ct_ptr = family_ct;
  *trio_list_ptr = trio_write;
  *trio_ct_ptr = trio_ct;
  *trio_lookup_ptr = trio_lookup;
  if (max_fid_len_ptr) {
    *max_fid_len_ptr = max_fid_len;
  }
  if (fids_ptr) {
    *fids_ptr = fids;
    *iids_ptr = iids;
    *max_iid_len_ptr = max_iid_len;
    sample_uidx = next_unset_unsafe(sample_exclude, 0);
    for (sample_idx = 0; sample_idx < sample_ct; sample_uidx++, sample_idx++) {
      next_unset_ul_unsafe_ck(sample_exclude, &sample_uidx);
      idptr = &(sample_ids[sample_uidx * max_sample_id_len]);
      iidptr = (char*)memchr(idptr, '\t', max_fid_len);
      strcpy(&(iids[sample_uidx * max_iid_len]), &(iidptr[1]));
    }
    if (include_duos) {
      memcpy(&(iids[unfiltered_sample_ct * max_iid_len]), "0", 2);
    }
  }
  uiptr = trio_lookup;
  for (trio_idx = 0; trio_idx < trio_ct; trio_idx++) {
    trio_code = trio_write[trio_idx];
    sample_uidx = (uint32_t)trio_code;
    if (fids_ptr) {
      idptr = &(sample_ids[sample_uidx * max_sample_id_len]);
      iidptr = (char*)memchr(idptr, '\t', max_fid_len);
      fidlen = (uintptr_t)(iidptr - idptr);
      memcpyx(&(fids[trio_idx * max_fid_len]), idptr, fidlen, '\0');
    }
    family_code = family_list[(uintptr_t)(trio_code >> 32)];
    if (!toposort) {
      *uiptr++ = sample_uidx;
      *uiptr++ = (uint32_t)family_code;
      *uiptr++ = (uint32_t)(family_code >> 32);
    } else if (is_set(founder_info2, family_code >> 32) && is_set(founder_info2, (uint32_t)family_code)) {
      // Just populate the "no incoming edges" queue here.
      toposort_queue[tqueue_end++] = trio_idx;
    }
  }
  if (toposort) {
    // Kahn A. B. (1962) Topological sorting of large networks.
    // Communications of the ACM, 5.
    // This is a breadth-first implementation.
    while (tqueue_start < tqueue_end) {
      trio_idx = toposort_queue[tqueue_start++];
      trio_code = trio_write[trio_idx];
      sample_uidx = (uint32_t)trio_code;
      family_code = family_list[(uintptr_t)(trio_code >> 32)];
      *uiptr++ = sample_uidx;
      *uiptr++ = (uint32_t)family_code;
      *uiptr++ = (uint32_t)(family_code >> 32);
      *uiptr++ = trio_idx;
      if (remaining_edge_ct) {
        SET_BIT(founder_info2, sample_uidx);
        ullii = ((uint64_t)sample_uidx) << 32;
        uii = uint64arr_greater_than(edge_list, edge_ct, ullii);
        ullii |= 0xffffffffLLU;
        while (uii < edge_ct) {
          edge_code = edge_list[uii];
          if (edge_code > ullii) {
	    break;
	  }
	  remaining_edge_ct--;
	  // look up child, see if other parent is now a founder
	  trio_code = trio_write[(uint32_t)edge_code];
	  family_code = family_list[(uintptr_t)(trio_code >> 32)];
	  if (is_set(founder_info2, family_code >> 32) && is_set(founder_info2, (uint32_t)family_code)) {
	    toposort_queue[tqueue_end++] = (uint32_t)edge_code;
	  }
	  uii++;
	}
      }
    }
    if (remaining_edge_ct) {
      logprint("Error: Pedigree graph is cyclic.  Check for evidence of time travel abuse in\nyour cohort.\n");
      goto get_trios_and_families_ret_INVALID_FORMAT;
    }
    wkspace_reset(edge_list);
  }
  while (0) {
  get_trios_and_families_ret_NOMEM2:
    wkspace_left += topsize;
  get_trios_and_families_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  get_trios_and_families_ret_INVALID_FORMAT_2:
    logprintb();
  get_trios_and_families_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 get_trios_and_families_ret_1:
  if (retval) {
    wkspace_reset(wkspace_mark);
  }
  return retval;
}

uint32_t erase_mendel_errors(uintptr_t unfiltered_sample_ct, uintptr_t* loadbuf, uintptr_t* workbuf, uint32_t* trio_lookup, uint32_t trio_ct, uint32_t multigen) {
  uint32_t* uiptr = trio_lookup;
  uint32_t cur_errors = 0;
  uint32_t trio_idx;
  uint32_t lookup_idx;
  uintptr_t ulii;
  uintptr_t uljj;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  uint32_t unn;
  uint32_t uoo;
  uint32_t upp;
  memcpy(workbuf, loadbuf, (unfiltered_sample_ct + 3) / 4);
  SET_BIT_DBL(workbuf, unfiltered_sample_ct);
  if (!multigen) {
    for (trio_idx = 0; trio_idx < trio_ct; trio_idx++) {
      uii = *uiptr++;
      ujj = *uiptr++;
      ukk = *uiptr++;
      umm = mendel_error_table[((workbuf[uii / BITCT2] >> (2 * (uii % BITCT2))) & 3) | (((workbuf[ujj / BITCT2] >> (2 * (ujj % BITCT2))) & 3) << 2) | (((workbuf[ukk / BITCT2] >> (2 * (ukk % BITCT2))) & 3) << 4)];
      if (umm) {
	ulii = loadbuf[uii / BITCT2];
	uljj = ONELU << (2 * (uii % BITCT2));
        loadbuf[uii / BITCT2] = (ulii & (~(3 * uljj))) | uljj;
	if (umm & 0x100) {
	  ulii = loadbuf[ujj / BITCT2];
	  uljj = ONELU << (2 * (ujj % BITCT2));
	  loadbuf[ujj / BITCT2] = (ulii & (~(3 * uljj))) | uljj;
	}
	if (umm & 0x10000) {
	  ulii = loadbuf[ukk / BITCT2];
	  uljj = ONELU << (2 * (ukk % BITCT2));
	  loadbuf[ukk / BITCT2] = (ulii & (~(3 * uljj))) | uljj;
	}
	cur_errors++;
      }
    }
  } else {
    for (lookup_idx = 0; lookup_idx < trio_ct; lookup_idx++) {
      uii = *uiptr++;
      ujj = *uiptr++;
      ukk = *uiptr++;
      trio_idx = *uiptr++;
      unn = uii / BITCT2;
      uoo = ujj / BITCT2;
      upp = ukk / BITCT2;
      uii = 2 * (uii % BITCT2);
      ujj = 2 * (ujj % BITCT2);
      ukk = 2 * (ukk % BITCT2);
      ulii = (workbuf[unn] >> uii) & 3;
      uljj = ((workbuf[uoo] >> ujj) & 3) | (((workbuf[upp] >> ukk) & 3) << 2);
      if (ulii != 1) {
        umm = mendel_error_table[ulii | (uljj << 2)];
        if (umm) {
	  ulii = loadbuf[unn];
	  uljj = ONELU << uii;
	  loadbuf[unn] = (ulii & (~(3 * uljj))) | uljj;
	  if (umm & 0x100) {
	    ulii = loadbuf[uoo];
	    uljj = ONELU << ujj;
	    loadbuf[uoo] = (ulii & (~(3 * uljj))) | uljj;
	  }
	  if (umm & 0x10000) {
	    ulii = loadbuf[upp];
	    uljj = ONELU << ukk;
	    loadbuf[upp] = (ulii & (~(3 * uljj))) | uljj;
	  }
	  cur_errors++;
	}
      } else if (!uljj) {
	workbuf[unn] &= ~(ONELU << uii);
      } else if (uljj == 15) {
	workbuf[unn] |= (2 * ONELU) << uii;
      }
    }
  }
  return cur_errors;
}

void fill_mendel_errstr(uint32_t error_code, char** allele_ptrs, uint32_t* alens, char* wbuf, uint32_t* len_ptr) {
  char* wptr;
  uint32_t len;
  if (error_code < 10) {
    wbuf[0] = ' ';
    wbuf[1] = error_code + '0';
  } else {
    memcpy(wbuf, "10", 2);
  }
  if (!alens[0]) {
    // lazy fill
    alens[0] = strlen(allele_ptrs[0]);
    alens[1] = strlen(allele_ptrs[1]);
  }
  // PLINK 1.07 may fail to put space here when there are long allele codes; we
  // don't replicate that.  (We actually force two spaces here since the error
  // column itself contains spaces.)
  wptr = memseta(&(wbuf[2]), 32, 2);
  switch (error_code) {
  case 1:
    len = 5 * alens[0] + alens[1] + 10;
    if (len < 20) {
      wptr = memseta(wptr, 32, 20 - len);
    }
    wptr = memcpyax(wptr, allele_ptrs[0], alens[0], '/');
    wptr = memcpya(wptr, allele_ptrs[0], alens[0]);
    wptr = memcpyl3a(wptr, " x ");
    wptr = memcpyax(wptr, allele_ptrs[0], alens[0], '/');
    wptr = memcpya(wptr, allele_ptrs[0], alens[0]);
    wptr = memcpya(wptr, " -> ", 4);
    wptr = memcpyax(wptr, allele_ptrs[0], alens[0], '/');
    wptr = memcpya(wptr, allele_ptrs[1], alens[1]);
    break;
  case 2:
    len = alens[0] + 5 * alens[1] + 10;
    if (len < 20) {
      wptr = memseta(wptr, 32, 20 - len);
    }
    wptr = memcpyax(wptr, allele_ptrs[1], alens[1], '/');
    wptr = memcpya(wptr, allele_ptrs[1], alens[1]);
    wptr = memcpyl3a(wptr, " x ");
    wptr = memcpyax(wptr, allele_ptrs[1], alens[1], '/');
    wptr = memcpya(wptr, allele_ptrs[1], alens[1]);
    wptr = memcpya(wptr, " -> ", 4);
    wptr = memcpyax(wptr, allele_ptrs[0], alens[0], '/');
    wptr = memcpya(wptr, allele_ptrs[1], alens[1]);
    break;
  case 3:
    len = 2 * alens[0] + 2 * alens[1] + 12;
    if (len < 20) {
      wptr = memseta(wptr, 32, 20 - len);
    }
    wptr = memcpyax(wptr, allele_ptrs[1], alens[1], '/');
    wptr = memcpya(wptr, allele_ptrs[1], alens[1]);
    wptr = memcpya(wptr, " x */* -> ", 10);
    wptr = memcpyax(wptr, allele_ptrs[0], alens[0], '/');
    wptr = memcpya(wptr, allele_ptrs[0], alens[0]);
    break;
  case 4:
  case 10:
    len = 2 * alens[0] + 2 * alens[1] + 12;
    if (len < 20) {
      wptr = memseta(wptr, 32, 20 - len);
    }
    wptr = memcpya(wptr, "*/* x ", 6);
    wptr = memcpyax(wptr, allele_ptrs[1], alens[1], '/');
    wptr = memcpya(wptr, allele_ptrs[1], alens[1]);
    wptr = memcpya(wptr, " -> ", 4);
    wptr = memcpyax(wptr, allele_ptrs[0], alens[0], '/');
    wptr = memcpya(wptr, allele_ptrs[0], alens[0]);
    break;
  case 5:
    len = 2 * alens[0] + 4 * alens[1] + 10;
    if (len < 20) {
      wptr = memseta(wptr, 32, 20 - len);
    }
    wptr = memcpyax(wptr, allele_ptrs[1], alens[1], '/');
    wptr = memcpya(wptr, allele_ptrs[1], alens[1]);
    wptr = memcpyl3a(wptr, " x ");
    wptr = memcpyax(wptr, allele_ptrs[1], alens[1], '/');
    wptr = memcpya(wptr, allele_ptrs[1], alens[1]);
    wptr = memcpya(wptr, " -> ", 4);
    wptr = memcpyax(wptr, allele_ptrs[0], alens[0], '/');
    wptr = memcpya(wptr, allele_ptrs[0], alens[0]);
    break;
  case 6:
    len = 2 * alens[0] + 2 * alens[1] + 12;
    if (len < 20) {
      wptr = memseta(wptr, 32, 20 - len);
    }
    wptr = memcpyax(wptr, allele_ptrs[0], alens[0], '/');
    wptr = memcpya(wptr, allele_ptrs[0], alens[0]);
    wptr = memcpya(wptr, " x */* -> ", 10);
    wptr = memcpyax(wptr, allele_ptrs[1], alens[1], '/');
    wptr = memcpya(wptr, allele_ptrs[1], alens[1]);
    break;
  case 7:
  case 9:
    len = 2 * alens[0] + 2 * alens[1] + 12;
    if (len < 20) {
      wptr = memseta(wptr, 32, 20 - len);
    }
    wptr = memcpya(wptr, "*/* x ", 6);
    wptr = memcpyax(wptr, allele_ptrs[0], alens[0], '/');
    wptr = memcpya(wptr, allele_ptrs[0], alens[0]);
    wptr = memcpya(wptr, " -> ", 4);
    wptr = memcpyax(wptr, allele_ptrs[1], alens[1], '/');
    wptr = memcpya(wptr, allele_ptrs[1], alens[1]);
    break;
  case 8:
    len = 4 * alens[0] + 2 * alens[1] + 10;
    if (len < 20) {
      wptr = memseta(wptr, 32, 20 - len);
    }
    wptr = memcpyax(wptr, allele_ptrs[0], alens[0], '/');
    wptr = memcpya(wptr, allele_ptrs[0], alens[0]);
    wptr = memcpyl3a(wptr, " x ");
    wptr = memcpyax(wptr, allele_ptrs[0], alens[0], '/');
    wptr = memcpya(wptr, allele_ptrs[0], alens[0]);
    wptr = memcpya(wptr, " -> ", 4);
    wptr = memcpyax(wptr, allele_ptrs[1], alens[1], '/');
    wptr = memcpya(wptr, allele_ptrs[1], alens[1]);
    break;
  }
  *wptr++ = '\n';
  *len_ptr = (uintptr_t)(wptr - wbuf);
}

int32_t mendel_error_scan(Family_info* fam_ip, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uint32_t plink_maxfid, uint32_t plink_maxiid, uint32_t plink_maxsnp, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_exclude_ct_ptr, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t* sample_exclude_ct_ptr, uintptr_t* founder_info, uintptr_t* sex_nm, uintptr_t* sex_male, char* sample_ids, uintptr_t max_sample_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uint32_t hh_exists, Chrom_info* chrom_info_ptr, uint32_t calc_mendel) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* outfile = NULL;
  FILE* outfile_l = NULL;
  uintptr_t* sample_male_include2 = NULL;
  uintptr_t* error_locs = NULL;
  char* varptr = NULL;
  char* chrom_name_ptr = NULL;
  unsigned char* cur_errors = NULL;
  uint64_t* family_error_cts = NULL;
  uint32_t* child_cts = NULL;
  uintptr_t marker_ct = unfiltered_marker_ct - *marker_exclude_ct_ptr;
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t unfiltered_sample_ctp1l2 = 1 + (unfiltered_sample_ct / BITCT2);
  uintptr_t final_mask = get_final_mask(unfiltered_sample_ct);
  uintptr_t sample_ct = unfiltered_sample_ct - *sample_exclude_ct_ptr;
  uintptr_t marker_uidx = ~ZEROLU;
  uint64_t tot_error_ct = 0;
  uint32_t unfiltered_sample_ctl2m1 = (unfiltered_sample_ct - 1) / BITCT2;
  uint32_t include_duos = (fam_ip->mendel_modifier / MENDEL_DUOS) & 1;
  uint32_t multigen = (fam_ip->mendel_modifier / MENDEL_MULTIGEN) & 1;
  uint32_t var_first = fam_ip->mendel_modifier & MENDEL_FILTER_VAR_FIRST;
  uint32_t varlen = 0;
  uint32_t chrom_name_len = 0;
  uint32_t new_marker_exclude_ct = 0;
  uint32_t error_ct_fill = 0;
  int32_t retval = 0;
  char chrom_name_buf[5];
  char* errstrs[10];
  uint32_t errstr_lens[11];
  uint32_t alens[2];
  const uint32_t* error_table_ptr;
  char* fids;
  char* iids;
  char* wptr;
  char* cptr;
  uint64_t* family_list;
  uint64_t* trio_list;
  uintptr_t* loadbuf;
  uint32_t* trio_lookup;
  uint32_t* error_cts;
  uint32_t* error_cts_tmp;
  uint32_t* error_cts_tmp2;
  uint32_t* uiptr;
#ifdef __LP64__
  __m128i* vptr;
  __m128i* vptr2;
#endif
  uintptr_t trio_ct4;
  uintptr_t max_fid_len;
  uintptr_t max_iid_len;
  uintptr_t trio_ct;
  uintptr_t trio_ctl;
  uintptr_t trio_idx;
  uintptr_t lookup_idx;
  uintptr_t ulii;
  uintptr_t uljj;
  double exclude_one_ratio;
  double dxx;
  double dyy;
  uint64_t trio_code;
  uint64_t family_code;
  uint32_t family_ct;
  uint32_t var_error_max;
  uint32_t cur_error_ct;
  uint32_t chrom_fo_idx;
  uint32_t chrom_idx;
  uint32_t chrom_end;
  uint32_t is_x;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  marker_ct -= count_non_autosomal_markers(chrom_info_ptr, marker_exclude, 0, 1);
  if ((!marker_ct) || is_set(chrom_info_ptr->haploid_mask, 0)) {
    logprint("Warning: Skipping --me/--mendel since there is no autosomal or Xchr data.\n");
    goto mendel_error_scan_ret_1;
  }
  retval = get_trios_and_families(unfiltered_sample_ct, sample_exclude, sample_ct, founder_info, sex_nm, sex_male, sample_ids, max_sample_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, &fids, &max_fid_len, &iids, &max_iid_len, &family_list, &family_ct, &trio_list, &trio_ct, &trio_lookup, include_duos, multigen);
  if (retval) {
    goto mendel_error_scan_ret_1;
  }
  if (!trio_ct) {
    LOGPRINTF("Warning: Skipping --me/--mendel since there are no %strios.\n", include_duos? "duos or " : "");
    goto mendel_error_scan_ret_1;
  }
  trio_ct4 = (trio_ct + 3) / 4;
  trio_ctl = (trio_ct + (BITCT - 1)) / BITCT;
  var_error_max = (int32_t)(fam_ip->mendel_max_var_error * (1 + SMALL_EPSILON) * ((intptr_t)trio_ct));
  if (wkspace_alloc_ul_checked(&loadbuf, unfiltered_sample_ctp1l2 * sizeof(intptr_t)) ||
      wkspace_alloc_ui_checked(&error_cts, trio_ct * 3 * sizeof(int32_t)) ||
      wkspace_alloc_ui_checked(&error_cts_tmp, trio_ct4 * 4 * sizeof(int32_t))) {
    goto mendel_error_scan_ret_NOMEM;
  }
  if (!var_first) {
    error_cts_tmp2 = error_cts_tmp;
  } else {
    if (wkspace_alloc_ui_checked(&error_cts_tmp2, trio_ct4 * 4 * sizeof(int32_t))) {
      goto mendel_error_scan_ret_NOMEM;
    }
    fill_uint_zero(error_cts_tmp2, trio_ct4 * 4);
  }
  loadbuf[unfiltered_sample_ctp1l2 - 1] = 0;
  fill_uint_zero(error_cts, trio_ct * 3);
  fill_uint_zero(error_cts_tmp, trio_ct4 * 4);
  hh_exists &= XMHH_EXISTS;
  if (alloc_raw_haploid_filters(unfiltered_sample_ct, hh_exists, 0, sample_exclude, sex_male, NULL, &sample_male_include2)) {
    goto mendel_error_scan_ret_NOMEM;
  }
  alens[0] = 0;
  alens[1] = 0;
  if (calc_mendel) {
    // ugh, this calculation was totally off before.
    // * max_marker_allele_len includes trailing null (though forgetting this
    //   was harmless)
    // * minimum buffer size was 25, not 21, due to inclusion of error code and
    //   following double-space in the string (forgetting this was NOT
    //   harmless)
    ulii = max_marker_allele_len * 6 + 9;
    if (ulii < 25) {
      ulii = 25;
    }
    if (wkspace_alloc_ull_checked(&family_error_cts, family_ct * 3 * sizeof(int64_t)) ||
        wkspace_alloc_ui_checked(&child_cts, family_ct * sizeof(int32_t)) ||
        wkspace_alloc_c_checked(&(errstrs[0]), ulii * 10)) {
      goto mendel_error_scan_ret_NOMEM;
    }
    for (uii = 1; uii < 10; uii++) {
      errstrs[uii] = &(errstrs[0][uii * ulii]);
    }
    if (multigen) {
      if (wkspace_alloc_ul_checked(&error_locs, trio_ctl * sizeof(intptr_t)) ||
          wkspace_alloc_uc_checked(&cur_errors, trio_ct)) {
	goto mendel_error_scan_ret_NOMEM;
      }
      fill_ulong_zero(error_locs, trio_ctl);
    }
    memcpy(outname_end, ".mendel", 8);
    if (fopen_checked(&outfile, outname, "w")) {
      goto mendel_error_scan_ret_OPEN_FAIL;
    }
    sprintf(tbuf, "%%%us %%%us  CHR %%%us   CODE                 ERROR\n", plink_maxfid, plink_maxiid, plink_maxsnp);
    fprintf(outfile, tbuf, "FID", "KID", "SNP");
    memcpy(outname_end, ".lmendel", 9);
    if (fopen_checked(&outfile_l, outname, "w")) {
      goto mendel_error_scan_ret_OPEN_FAIL;
    }
    // replicate harmless 'N' misalignment bug
    sprintf(tbuf, " CHR %%%us   N\n", plink_maxsnp);
    fprintf(outfile_l, tbuf, "SNP");
  } else {
    // suppress warning
    fill_ulong_zero((uintptr_t*)errstrs, 10);
  }
  for (chrom_fo_idx = 0; chrom_fo_idx < chrom_info_ptr->chrom_ct; chrom_fo_idx++) {
    chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
    is_x = (((uint32_t)chrom_info_ptr->x_code) == chrom_idx);
    if ((IS_SET(chrom_info_ptr->haploid_mask, chrom_idx) && (!is_x)) || (((uint32_t)chrom_info_ptr->mt_code) == chrom_idx)) {
      continue;
    }
    chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx + 1];
    uii = next_unset(marker_exclude, chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx], chrom_end);
    if (uii == chrom_end) {
      continue;
    }
    if (!is_x) {
      error_table_ptr = mendel_error_table;
    } else {
      error_table_ptr = mendel_error_table_x;
    }
    if (calc_mendel) {
      chrom_name_ptr = chrom_name_buf5w4write(chrom_name_buf, chrom_info_ptr, chrom_idx, &chrom_name_len);
    }
    if (uii != marker_uidx) {
      marker_uidx = uii;
      goto mendel_error_scan_seek;
    }
    while (1) {
      if (load_raw2(bedfile, loadbuf, unfiltered_sample_ct4, unfiltered_sample_ctl2m1, final_mask)) {
	goto mendel_error_scan_ret_READ_FAIL;
      }
      if (IS_SET(marker_reverse, marker_uidx)) {
        reverse_loadbuf((unsigned char*)loadbuf, unfiltered_sample_ct);
      }
      if (hh_exists && is_x) {
	hh_reset((unsigned char*)loadbuf, sample_male_include2, unfiltered_sample_ct);
      }
      // missing parents are treated as having uidx equal to
      // unfiltered_sample_ct, and we set the corresponding genotype to always
      // be missing.  This lets us avoid special-casing duos.
      SET_BIT_DBL(loadbuf, unfiltered_sample_ct);
      uiptr = trio_lookup;
      cur_error_ct = 0;
      if (calc_mendel) {
	varptr = &(marker_ids[marker_uidx * max_marker_id_len]);
	varlen = strlen(varptr);
	alens[0] = 0;
	alens[1] = 0;
	fill_uint_zero(errstr_lens, 11);
      }
      if (!multigen) {
	for (trio_idx = 0; trio_idx < trio_ct; trio_idx++) {
	  uii = *uiptr++;
	  ujj = *uiptr++;
	  ukk = *uiptr++;
          umm = error_table_ptr[((loadbuf[uii / BITCT2] >> (2 * (uii % BITCT2))) & 3) | (((loadbuf[ujj / BITCT2] >> (2 * (ujj % BITCT2))) & 3) << 2) | (((loadbuf[ukk / BITCT2] >> (2 * (ukk % BITCT2))) & 3) << 4)];
	  if (umm) {
	    error_cts_tmp2[trio_idx] += umm & 0xffffff;
	    cur_error_ct++;
	    if (calc_mendel) {
	      umm >>= 24;
	      if ((umm > 8) && (!is_set(sex_male, uii))) {
		umm = 34 - 3 * umm; // 9 -> 7, 10 -> 4
	      }
	      wptr = fw_strcpy(plink_maxfid, &(fids[trio_idx * max_fid_len]), tbuf);
	      *wptr++ = ' ';
	      wptr = fw_strcpy(plink_maxiid, &(iids[uii * max_iid_len]), wptr);
	      *wptr++ = ' ';
	      wptr = memcpyax(wptr, chrom_name_ptr, chrom_name_len, ' ');
	      wptr = fw_strcpyn(plink_maxsnp, varlen, varptr, wptr);
	      wptr = memseta(wptr, 32, 5);
	      if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
		goto mendel_error_scan_ret_WRITE_FAIL;
	      }
	      if (!errstr_lens[umm]) {
		fill_mendel_errstr(umm, &(marker_allele_ptrs[2 * marker_uidx]), alens, errstrs[umm - 1], &(errstr_lens[umm]));
	      }
	      if (fwrite_checked(errstrs[umm - 1], errstr_lens[umm], outfile)) {
		goto mendel_error_scan_ret_WRITE_FAIL;
	      }
	    }
	  }
	}
      } else {
	for (lookup_idx = 0; lookup_idx < trio_ct; lookup_idx++) {
	  uii = *uiptr++;
	  ujj = *uiptr++;
	  ukk = *uiptr++;
          trio_idx = *uiptr++;
          uljj = ((loadbuf[ujj / BITCT2] >> (2 * (ujj % BITCT2))) & 3) | (((loadbuf[ukk / BITCT2] >> (2 * (ukk % BITCT2))) & 3) << 2);
	  umm = uii / BITCT2;
	  ujj = 2 * (uii % BITCT2);
	  ulii = (loadbuf[umm] >> ujj) & 3;
	  if (ulii != 1) {
	    umm = error_table_ptr[ulii | (uljj << 2)];
	    if (umm) {
	      error_cts_tmp2[trio_idx] += umm & 0xffffff;
	      cur_error_ct++;
	      if (calc_mendel) {
	        set_bit(error_locs, trio_idx);
		umm >>= 24;
		if ((umm > 8) && (!is_set(sex_male, uii))) {
		  umm = 34 - 3 * umm;
		}
                cur_errors[trio_idx] = (unsigned char)umm;
	      }
	    }
	  } else if (!uljj) {
	    loadbuf[umm] &= ~(ONELU << ujj);
          } else if (uljj == 15) {
	    loadbuf[umm] |= (2 * ONELU) << ujj;
	  }
	}
	if (calc_mendel && cur_error_ct) {
          trio_idx = 0;
	  for (uii = 0; uii < cur_error_ct; trio_idx++, uii++) {
            next_set_ul_unsafe_ck(error_locs, &trio_idx);
	    wptr = fw_strcpy(plink_maxfid, &(fids[trio_idx * max_fid_len]), tbuf);
	    *wptr++ = ' ';
	    wptr = fw_strcpy(plink_maxiid, &(iids[((uint32_t)trio_list[trio_idx]) * max_iid_len]), wptr);
	    *wptr++ = ' ';
	    wptr = memcpyax(wptr, chrom_name_ptr, chrom_name_len, ' ');
	    wptr = fw_strcpyn(plink_maxsnp, varlen, varptr, wptr);
	    wptr = memseta(wptr, 32, 5);
	    if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	      goto mendel_error_scan_ret_WRITE_FAIL;
	    }
	    umm = cur_errors[trio_idx];
	    if (!errstr_lens[umm]) {
	      fill_mendel_errstr(umm, &(marker_allele_ptrs[2 * marker_uidx]), alens, errstrs[umm - 1], &(errstr_lens[umm]));
	    }
	    if (fwrite_checked(errstrs[umm - 1], errstr_lens[umm], outfile)) {
	      goto mendel_error_scan_ret_WRITE_FAIL;
	    }
	  }
          fill_ulong_zero(error_locs, trio_ctl);
	}
      }
      if (calc_mendel) {
	if (fwrite_checked(chrom_name_ptr, chrom_name_len, outfile_l)) {
	  goto mendel_error_scan_ret_WRITE_FAIL;
	}
	tbuf[0] = ' ';
	wptr = fw_strcpyn(plink_maxsnp, varlen, varptr, &(tbuf[1]));
        *wptr++ = ' ';
        wptr = uint32_writew4x(wptr, cur_error_ct, '\n');
	if (fwrite_checked(tbuf, wptr - tbuf, outfile_l)) {
	  goto mendel_error_scan_ret_WRITE_FAIL;
	}
      }
      if (cur_error_ct) {
	if (cur_error_ct > var_error_max) {
	  SET_BIT(marker_exclude, marker_uidx);
	  new_marker_exclude_ct++;
	}
	if ((cur_error_ct <= var_error_max) || (!var_first)) {
	  if (var_first) {
#ifdef __LP64__
	    vptr = (__m128i*)error_cts_tmp;
	    vptr2 = (__m128i*)error_cts_tmp2;
	    for (trio_idx = 0; trio_idx < trio_ct4; trio_idx++) {
	      *vptr = _mm_add_epi64(*vptr, *vptr2++);
	      vptr++;
	    }
#else
            for (trio_idx = 0; trio_idx < trio_ct; trio_idx++) {
	      error_cts_tmp[trio_idx] += error_cts_tmp2[trio_idx];
	    }
#endif
	    fill_uint_zero(error_cts_tmp2, trio_ct);
	  }
	  error_ct_fill++;
	  if (error_ct_fill == 255) {
	    uiptr = error_cts;
            for (trio_idx = 0; trio_idx < trio_ct; trio_idx++) {
	      uii = error_cts_tmp[trio_idx];
	      *uiptr += (unsigned char)uii;
	      uiptr++;
	      *uiptr += (unsigned char)(uii >> 8);
	      uiptr++;
	      *uiptr += uii >> 16;
	      uiptr++;
	    }
	    fill_uint_zero(error_cts_tmp, trio_ct);
	    error_ct_fill = 0;
	  }
	}
      }
      tot_error_ct += cur_error_ct;
      if (++marker_uidx == chrom_end) {
	break;
      }
      if (IS_SET(marker_exclude, marker_uidx)) {
        marker_uidx = next_unset_ul(marker_exclude, marker_uidx, chrom_end);
      mendel_error_scan_seek:
        if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	  goto mendel_error_scan_ret_READ_FAIL;
	}
	if (marker_uidx == chrom_end) {
	  break;
	}
      }
    }
  }
  if (error_ct_fill) {
    uiptr = error_cts;
    for (trio_idx = 0; trio_idx < trio_ct; trio_idx++) {
      uii = error_cts_tmp[trio_idx];
      *uiptr += (unsigned char)uii;
      uiptr++;
      *uiptr += (unsigned char)(uii >> 8);
      uiptr++;
      *uiptr += uii >> 16;
      uiptr++;
    }
  }
  LOGPRINTF("--me/--mendel: %" PRIu64 " Mendel error%s detected.\n", tot_error_ct, (tot_error_ct == 1)? "" : "s");
  if (calc_mendel) {
    if (fclose_null(&outfile)) {
      goto mendel_error_scan_ret_WRITE_FAIL;
    }
    if (fclose_null(&outfile_l)) {
      goto mendel_error_scan_ret_WRITE_FAIL;
    }
    outname_end[1] = 'f';
    if (fopen_checked(&outfile, outname, "w")) {
      goto mendel_error_scan_ret_OPEN_FAIL;
    }
    sprintf(tbuf, "%%%us %%%us %%%us   CHLD    N\n", plink_maxfid, plink_maxiid, plink_maxiid);
    fprintf(outfile, tbuf, "FID", "PAT", "MAT");
    fill_ull_zero(family_error_cts, family_ct * 3);
    fill_uint_zero(child_cts, family_ct);
    for (trio_idx = 0; trio_idx < trio_ct; trio_idx++) {
      uii = (uint32_t)(trio_list[trio_idx] >> 32);
      child_cts[uii] += 1;
      family_error_cts[uii * 3] += error_cts[trio_idx * 3];
      family_error_cts[uii * 3 + 1] += error_cts[trio_idx * 3 + 1];
      family_error_cts[uii * 3 + 2] += error_cts[trio_idx * 3 + 2];
    }
    for (uii = 0; uii < family_ct; uii++) {
      family_code = family_list[uii];
      ujj = (uint32_t)family_code; // paternal uidx
      ukk = (uint32_t)(family_code >> 32); // maternal uidx
      if (ujj < unfiltered_sample_ct) {
	// bleah, fids[] isn't in right order for this lookup
	cptr = &(sample_ids[ujj * max_sample_id_len]);
	wptr = fw_strcpyn(plink_maxfid, (uintptr_t)(((char*)memchr(cptr, '\t', max_sample_id_len)) - cptr), cptr, tbuf);
      } else {
	cptr = &(sample_ids[ukk * max_sample_id_len]);
	wptr = fw_strcpyn(plink_maxfid, (uintptr_t)(((char*)memchr(cptr, '\t', max_sample_id_len)) - cptr), cptr, tbuf);
	// wptr = memseta(tbuf, 32, plink_maxfid - 1);
	// *wptr++ = '0';
      }
      *wptr++ = ' ';
      if (ujj != unfiltered_sample_ct) {
        wptr = fw_strcpy(plink_maxiid, &(iids[ujj * max_iid_len]), wptr);
      } else {
	wptr = memseta(wptr, 32, plink_maxiid - 1);
	*wptr++ = '0';
      }
      *wptr++ = ' ';
      if (ukk != unfiltered_sample_ct) {
        wptr = fw_strcpy(plink_maxiid, &(iids[ukk * max_iid_len]), wptr);
      } else {
	wptr = memseta(wptr, 32, plink_maxiid - 1);
	*wptr++ = '0';
      }
      *wptr++ = ' ';
      wptr = uint32_writew6x(wptr, child_cts[uii], ' ');
      if (family_error_cts[uii * 3] < 10000) {
	wptr = uint32_writew4(wptr, (uint32_t)family_error_cts[uii * 3]);
      } else {
        wptr = int64_write(wptr, family_error_cts[uii * 3]);
      }
      *wptr++ = '\n';
      if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	goto mendel_error_scan_ret_WRITE_FAIL;
      }
    }
    if (fclose_null(&outfile)) {
      goto mendel_error_scan_ret_WRITE_FAIL;
    }
    outname_end[1] = 'i';
    if (fopen_checked(&outfile, outname, "w")) {
      goto mendel_error_scan_ret_OPEN_FAIL;
    }
    sprintf(tbuf, "%%%us %%%us   N\n", plink_maxfid, plink_maxiid);
    fprintf(outfile, tbuf, "FID", "IID");
    uii = 0xffffffffU; // family idx
    for (trio_idx = 0; trio_idx < trio_ct; trio_idx++) {
      trio_code = trio_list[trio_idx];
      ujj = (uint32_t)(trio_code >> 32);
      if (ujj != uii) {
	uii = ujj;
        family_code = family_list[uii];
	wptr = fw_strcpy(plink_maxfid, &(fids[trio_idx * max_fid_len]), tbuf);
	*wptr++ = ' ';
	ujj = (uint32_t)family_code;
	if (ujj != unfiltered_sample_ct) {
	  wptr = fw_strcpy(plink_maxiid, &(iids[ujj * max_iid_len]), wptr);
	  *wptr++ = ' ';
	  if (family_error_cts[3 * uii + 1] < 10000) {
	    wptr = uint32_writew4(wptr, (uint32_t)family_error_cts[3 * uii + 1]);
	  } else {
	    wptr = int64_write(wptr, family_error_cts[3 * uii + 1]);
	  }
	  if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	    goto mendel_error_scan_ret_WRITE_FAIL;
	  }
	}
	ukk = (uint32_t)(family_code >> 32);
	if (ukk != unfiltered_sample_ct) {
	  if (ujj != unfiltered_sample_ct) {
	    putc('\n', outfile);
	  }
	  wptr = fw_strcpy(plink_maxiid, &(iids[ukk * max_iid_len]), &(tbuf[plink_maxfid + 1]));
	  *wptr++ = ' ';
	  if (family_error_cts[3 * uii + 2] < 10000) {
	    wptr = uint32_writew4(wptr, (uint32_t)family_error_cts[3 * uii + 2]);
	  } else {
	    wptr = int64_write(wptr, family_error_cts[3 * uii + 2]);
	  }
	  if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	    goto mendel_error_scan_ret_WRITE_FAIL;
	  }
	}
	putc(' ', outfile); // PLINK 1.07 formatting quirk
	putc('\n', outfile);
      }
      wptr = fw_strcpy(plink_maxiid, &(iids[((uint32_t)trio_code) * max_iid_len]), &(tbuf[plink_maxfid + 1]));
      *wptr++ = ' ';
      wptr = uint32_writew4x(wptr, error_cts[trio_idx * 3], '\n');
      if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	goto mendel_error_scan_ret_WRITE_FAIL;
      }
    }
    *outname_end = '\0';
    LOGPRINTFWW("Reports written to %s.mendel + %s.imendel + %s.fmendel + %s.lmendel .\n", outname, outname, outname, outname);
  }
  if (fam_ip->mendel_modifier & MENDEL_FILTER) {
    *marker_exclude_ct_ptr += new_marker_exclude_ct;
    if (unfiltered_marker_ct == *marker_exclude_ct_ptr) {
      logprint("Error: All variants excluded by --me.\n");
      goto mendel_error_scan_ret_ALL_MARKERS_EXCLUDED;
    }
    if (var_first) {
      marker_ct -= new_marker_exclude_ct;
    }
    uii = (int32_t)(fam_ip->mendel_max_trio_error * (1 + SMALL_EPSILON) * ((intptr_t)marker_ct));
    if (uii < marker_ct) {
      exclude_one_ratio = fam_ip->mendel_exclude_one_ratio;
      for (trio_idx = 0; trio_idx < trio_ct; trio_idx++) {
	if (error_cts[trio_idx * 3] > uii) {
	  trio_code = trio_list[trio_idx];
	  family_code = family_list[(uintptr_t)(trio_code >> 32)];
	  ujj = (uint32_t)family_code;
	  ukk = (uint32_t)(family_code >> 32);
          if (exclude_one_ratio == 0.0) {
	    set_bit(sample_exclude, (uint32_t)trio_code);
	    if (ujj < unfiltered_sample_ct) {
	      set_bit(sample_exclude, ujj);
	    }
	    if (ukk < unfiltered_sample_ct) {
	      set_bit(sample_exclude, ukk);
	    }
	  } else if ((exclude_one_ratio == -1) || (ujj == unfiltered_sample_ct) || (ukk == unfiltered_sample_ct)) {
            set_bit(sample_exclude, (uint32_t)trio_code);
	  } else {
	    dxx = (double)((int32_t)trio_list[trio_idx * 3 + 1]);
	    dyy = (double)((int32_t)trio_list[trio_idx * 3 + 2]);
	    if (dxx > exclude_one_ratio * dyy) {
	      set_bit(sample_exclude, ujj);
	    } else if (dyy > exclude_one_ratio * dxx) {
	      set_bit(sample_exclude, ukk);
	    } else {
	      set_bit(sample_exclude, (uint32_t)trio_code);
	    }
	  }
	}
      }
    }
    ulii = popcount_longs(sample_exclude, (unfiltered_sample_ct + (BITCT - 1)) / BITCT);
    if (unfiltered_sample_ct == ulii) {
      LOGPRINTF("Error: All %s excluded by --me.\n", g_species_plural);
      goto mendel_error_scan_ret_ALL_SAMPLES_EXCLUDED;
    }
    *sample_exclude_ct_ptr = ulii;
    ulii -= unfiltered_sample_ct - sample_ct;
    LOGPRINTF("%u variant%s and %" PRIuPTR " %s excluded.\n", new_marker_exclude_ct, (new_marker_exclude_ct == 1)? "" : "s", ulii, species_str(ulii));
  }
  while (0) {
  mendel_error_scan_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  mendel_error_scan_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  mendel_error_scan_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  mendel_error_scan_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  mendel_error_scan_ret_ALL_MARKERS_EXCLUDED:
    retval = RET_ALL_MARKERS_EXCLUDED;
    break;
  mendel_error_scan_ret_ALL_SAMPLES_EXCLUDED:
    retval = RET_ALL_SAMPLES_EXCLUDED;
    break;
  }
 mendel_error_scan_ret_1:
  wkspace_reset(wkspace_mark);
  fclose_cond(outfile);
  fclose_cond(outfile_l);
  return retval;
}

int32_t populate_pedigree_rel_info(Pedigree_rel_info* pri_ptr, uintptr_t unfiltered_sample_ct, char* sample_ids, uintptr_t max_sample_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* founder_info) {
  // possible todo: if any families have been entirely filtered out, don't
  // construct pedigree for them
  unsigned char* wkspace_mark;
  unsigned char* wkspace_mark2;
  int32_t ii;
  int32_t jj;
  int32_t kk;
  int32_t mm;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t initial_family_blocks;
  uintptr_t ulii;
  uintptr_t sample_uidx;
  uint64_t ullii;
  char* family_ids;
  char* cur_sample_id;
  char* last_family_id = NULL;
  char* cur_family_id;
  char* id_ptr;
  uint32_t* family_sizes;
  uint32_t* uiptr;
  uint32_t* uiptr2 = NULL;
  uint32_t fidx;
  int32_t family_size;
  uint32_t* remaining_sample_idxs;
  int32_t* remaining_sample_parent_idxs; // -1 = no parent (or nonshared)
  uint32_t remaining_sample_ct;
  uint32_t sample_idx_write;
  uintptr_t max_family_id_len = 0;
  char* indiv_ids;
  uint32_t* sample_id_lookup;
  uintptr_t max_indiv_id_len = 0;
  uintptr_t max_pm_id_len;
  uint32_t family_id_ct;
  uint32_t* fis_ptr;
  char* stray_parent_ids;
  intptr_t stray_parent_ct;
  uintptr_t* processed_samples;
  uint32_t founder_ct;
  int32_t max_family_nf = 0;
  uintptr_t unfiltered_sample_ctl = (unfiltered_sample_ct + (BITCT - 1)) / BITCT;
  uintptr_t unfiltered_sample_ctlm = unfiltered_sample_ctl * BITCT;
  uint32_t* complete_sample_idxs;
  uintptr_t complete_sample_idx_ct;
  double* rs_ptr;
  double* rel_writer;
  double dxx;
  double* tmp_rel_space = NULL;
  double* tmp_rel_writer = NULL;

  for (sample_uidx = 0; sample_uidx < unfiltered_sample_ct; sample_uidx++) {
    ujj = strlen_se(&(sample_ids[sample_uidx * max_sample_id_len])) + 1;
    if (ujj > max_family_id_len) {
      max_family_id_len = ujj;
    }
    ujj = strlen(&(sample_ids[sample_uidx * max_sample_id_len + ujj]));
    if (ujj >= max_indiv_id_len) {
      max_indiv_id_len = ujj + 1;
    }
  }
  if (max_paternal_id_len > max_maternal_id_len) {
    max_pm_id_len = max_paternal_id_len;
  } else {
    max_pm_id_len = max_maternal_id_len;
  }
  if (wkspace_alloc_ui_checked(&(pri_ptr->family_info_space), unfiltered_sample_ct * sizeof(int32_t)) ||
      wkspace_alloc_ui_checked(&(pri_ptr->family_rel_nf_idxs), unfiltered_sample_ct * sizeof(int32_t)) ||
      wkspace_alloc_ui_checked(&(pri_ptr->family_idxs), unfiltered_sample_ct * sizeof(int32_t)) ||
      wkspace_alloc_c_checked(&family_ids, unfiltered_sample_ct * max_family_id_len) ||
      wkspace_alloc_ui_checked(&family_sizes, unfiltered_sample_ct * sizeof(int32_t))) {
    return RET_NOMEM;
  }

  // copy all the items over in order, then qsort, then eliminate duplicates
  // and count family sizes.
  cur_family_id = family_ids;
  cur_sample_id = sample_ids;
  uiptr = family_sizes;
  *uiptr = 1;
  jj = strlen_se(cur_sample_id);
  memcpyx(cur_family_id, cur_sample_id, jj, 0);
  for (sample_uidx = 1; sample_uidx < unfiltered_sample_ct; sample_uidx++) {
    cur_sample_id = &(cur_sample_id[max_sample_id_len]);
    mm = strlen_se(cur_sample_id);
    if ((jj != mm) || memcmp(cur_family_id, cur_sample_id, mm)) {
      cur_family_id = &(cur_family_id[max_family_id_len]);
      memcpyx(cur_family_id, cur_sample_id, mm, 0);
      jj = mm;
      *(++uiptr) = 1;
    } else {
      *uiptr += 1;
    }
  }
  initial_family_blocks = 1 + (uint32_t)(uiptr - family_sizes);
  if (qsort_ext(family_ids, initial_family_blocks, max_family_id_len, strcmp_deref, (char*)family_sizes, sizeof(int32_t))) {
    return RET_NOMEM;
  }

  last_family_id = family_ids;
  cur_family_id = &(family_ids[max_family_id_len]);
  family_id_ct = 1;
  uii = 1; // read idx
  if (initial_family_blocks != 1) {
    uiptr = family_sizes;
    while (strcmp(cur_family_id, last_family_id)) {
      family_id_ct++;
      uiptr++;
      last_family_id = cur_family_id;
      cur_family_id = &(cur_family_id[max_family_id_len]);
      uii++;
      if (uii == initial_family_blocks) {
	break;
      }
    }
    if (uii < initial_family_blocks) {
      uiptr2 = uiptr; // family_sizes read pointer
      *uiptr += *(++uiptr2);
      uii++;
      cur_family_id = &(cur_family_id[max_family_id_len]); // read pointer
      while (uii < initial_family_blocks) {
	while (!strcmp(cur_family_id, last_family_id)) {
	  *uiptr += *(++uiptr2);
	  uii++;
	  if (uii == initial_family_blocks) {
	    break;
	  }
	  cur_family_id = &(cur_family_id[max_family_id_len]);
	}
	if (uii < initial_family_blocks) {
	  *(++uiptr) = *(++uiptr2);
	  last_family_id = &(last_family_id[max_family_id_len]);
	  strcpy(last_family_id, cur_family_id);
	  family_id_ct++;
	  uii++;
	  cur_family_id = &(cur_family_id[max_family_id_len]);
	}
      }
    }
  }

  if (family_id_ct < unfiltered_sample_ct) {
    uiptr = family_sizes;
    wkspace_shrink_top(family_ids, family_id_ct * max_family_id_len);
    family_sizes = (uint32_t*)wkspace_alloc(family_id_ct * sizeof(int32_t));
    if (family_sizes < uiptr) {
      // copy back
      for (uii = 0; uii < family_id_ct; uii++) {
	family_sizes[uii] = *uiptr++;
      }
    }
  }
  pri_ptr->family_ids = family_ids;
  pri_ptr->family_id_ct = family_id_ct;
  pri_ptr->max_family_id_len = max_family_id_len;
  pri_ptr->family_sizes = family_sizes;

  if (wkspace_alloc_ui_checked(&(pri_ptr->family_info_offsets), (family_id_ct + 1) * sizeof(int32_t)) ||
      wkspace_alloc_ul_checked(&(pri_ptr->family_rel_space_offsets), (family_id_ct + 1) * sizeof(intptr_t)) ||
      wkspace_alloc_ui_checked(&(pri_ptr->family_founder_cts), family_id_ct * sizeof(int32_t))) {
    return RET_NOMEM;
  }
  fill_int_zero((int32_t*)(pri_ptr->family_founder_cts), family_id_ct);

  ii = 0; // running family_info offset
  for (fidx = 0; fidx < family_id_ct; fidx++) {
    family_size = pri_ptr->family_sizes[fidx];
    pri_ptr->family_info_offsets[fidx] = ii;
    ii += family_size;
  }

  if (wkspace_alloc_ui_checked(&uiptr, family_id_ct * sizeof(int32_t))) {
    return RET_NOMEM;
  }
  fill_uint_zero(uiptr, family_id_ct);

  // Fill family_idxs, family_founder_cts, and founder portion of
  // family_rel_nf_idxs.
  cur_sample_id = sample_ids;
  for (sample_uidx = 0; sample_uidx < unfiltered_sample_ct; sample_uidx++) {
    kk = bsearch_str(cur_sample_id, strlen_se(cur_sample_id), family_ids, max_family_id_len, family_id_ct);
    pri_ptr->family_idxs[sample_uidx] = kk;
    if (IS_SET(founder_info, sample_uidx)) {
      pri_ptr->family_founder_cts[(uint32_t)kk] += 1;
      pri_ptr->family_rel_nf_idxs[sample_uidx] = uiptr[(uint32_t)kk];
      uiptr[kk] += 1;
    }
    cur_sample_id = &(cur_sample_id[max_sample_id_len]);
  }
  wkspace_reset(uiptr);
  ulii = 0; // running rel_space offset
  for (fidx = 0; fidx < family_id_ct; fidx++) {
    family_size = pri_ptr->family_sizes[fidx];
    pri_ptr->family_rel_space_offsets[fidx] = ulii;
    kk = pri_ptr->family_founder_cts[fidx];
    if (family_size - kk > max_family_nf) {
      max_family_nf = family_size - kk;
    }
    // No need to explicitly store the (kk * (kk - 1)) / 2 founder-founder
    // relationships.
    ulii += (((int64_t)family_size) * (family_size - 1) - ((int64_t)kk) * (kk - 1)) / 2;
  }

  // make it safe to determine size of blocks by subtracting from the next
  // offset, even if we're at the last family
  pri_ptr->family_info_offsets[family_id_ct] = unfiltered_sample_ct;
  pri_ptr->family_rel_space_offsets[family_id_ct] = ulii;
  if (wkspace_alloc_d_checked(&(pri_ptr->rel_space), ulii * sizeof(double))) {
    return RET_NOMEM;
  }

  wkspace_mark = wkspace_base;
  if (wkspace_alloc_ui_checked(&uiptr, family_id_ct * sizeof(int32_t))) {
    return RET_NOMEM;
  }
  // populate family_info_space
  for (fidx = 0; fidx < family_id_ct; fidx++) {
    uiptr[fidx] = pri_ptr->family_info_offsets[fidx];
  }
  for (sample_uidx = 0; sample_uidx < unfiltered_sample_ct; sample_uidx++) {
    fidx = pri_ptr->family_idxs[sample_uidx];
    pri_ptr->family_info_space[uiptr[fidx]] = sample_uidx;
    uiptr[fidx] += 1;
  }
  wkspace_reset(wkspace_mark);

  if (wkspace_alloc_ul_checked(&processed_samples, (unfiltered_sample_ctl + (max_family_nf + (BITCT2 - 1)) / BITCT2) * sizeof(intptr_t))) {
    return RET_NOMEM;
  }
  fill_ulong_one(&(processed_samples[unfiltered_sample_ctl]), (max_family_nf + (BITCT2 - 1)) / BITCT2);

  wkspace_mark2 = wkspace_base;
  for (fidx = 0; fidx < family_id_ct; fidx++) {
    family_size = family_sizes[fidx];
    founder_ct = pri_ptr->family_founder_cts[fidx];
    remaining_sample_ct = family_size - founder_ct;
    stray_parent_ct = 0;
    if (remaining_sample_ct) {
      memcpy(processed_samples, founder_info, unfiltered_sample_ctl * sizeof(intptr_t));
      if (wkspace_alloc_ui_checked(&complete_sample_idxs, family_size * sizeof(int32_t)) ||
          wkspace_alloc_ui_checked(&remaining_sample_idxs, remaining_sample_ct * sizeof(int32_t)) ||
          wkspace_alloc_c_checked(&indiv_ids, family_size * max_indiv_id_len) ||
          wkspace_alloc_ui_checked(&sample_id_lookup, family_size * sizeof(int32_t)) ||
          wkspace_alloc_i_checked(&remaining_sample_parent_idxs, remaining_sample_ct * 2 * sizeof(int32_t)) ||
          wkspace_alloc_c_checked(&stray_parent_ids, remaining_sample_ct * 2 * max_pm_id_len)) {
	return RET_NOMEM;
      }
      ii = pri_ptr->family_info_offsets[fidx];
      fis_ptr = &(pri_ptr->family_info_space[ii]);
      rs_ptr = &(pri_ptr->rel_space[pri_ptr->family_rel_space_offsets[fidx]]);
      rel_writer = rs_ptr;
      cur_sample_id = indiv_ids;
      for (ii = 0; ii < family_size; ii++) {
	kk = fis_ptr[(uint32_t)ii];
	jj = strlen_se(&(sample_ids[kk * max_sample_id_len])) + 1;
	strcpy(cur_sample_id, &(sample_ids[kk * max_sample_id_len + jj]));
	cur_sample_id = &(cur_sample_id[max_indiv_id_len]);
	sample_id_lookup[(uint32_t)ii] = ii;
      }

      if (qsort_ext(indiv_ids, family_size, max_indiv_id_len, strcmp_deref, (char*)sample_id_lookup, sizeof(int32_t))) {
	return RET_NOMEM;
      }
      // Compile list of non-founder family member indices, and identify
      // parents who are referred to by sample ID but are NOT in the dataset.
      ii = 0; // family_info_space index
      complete_sample_idx_ct = 0;
      cur_sample_id = stray_parent_ids;
      for (uii = 0; uii < remaining_sample_ct; uii++) {
	while (IS_SET(founder_info, fis_ptr[ii])) {
	  complete_sample_idxs[complete_sample_idx_ct++] = fis_ptr[ii];
	  ii++;
	}
	kk = fis_ptr[ii++];

	// does not track sex for now
	id_ptr = &(paternal_ids[((uint32_t)kk) * max_paternal_id_len]);
	if (memcmp("0", id_ptr, 2)) {
	  ujj = strlen(id_ptr);
	  mm = bsearch_str(id_ptr, ujj, indiv_ids, max_indiv_id_len, family_size);
	  if (mm == -1) {
	    memcpy(cur_sample_id, id_ptr, ujj + 1);
	    cur_sample_id = &(cur_sample_id[max_pm_id_len]);
	    stray_parent_ct++;
	    remaining_sample_parent_idxs[uii * 2] = -2;
	  } else {
            remaining_sample_parent_idxs[uii * 2] = fis_ptr[sample_id_lookup[(uint32_t)mm]];
	  }
	} else {
          remaining_sample_parent_idxs[uii * 2] = -1;
	}
	id_ptr = &(maternal_ids[((uint32_t)kk) * max_maternal_id_len]);
	if (memcmp("0", id_ptr, 2)) {
	  ujj = strlen(id_ptr);
          mm = bsearch_str(id_ptr, ujj, indiv_ids, max_indiv_id_len, family_size);
	  if (mm == -1) {
	    memcpy(cur_sample_id, id_ptr, ujj + 1);
	    cur_sample_id = &(cur_sample_id[max_pm_id_len]);
	    stray_parent_ct++;
	    remaining_sample_parent_idxs[uii * 2 + 1] = -2;
	  } else {
	    remaining_sample_parent_idxs[uii * 2 + 1] = fis_ptr[sample_id_lookup[(uint32_t)mm]];
	  }
	} else {
	  remaining_sample_parent_idxs[uii * 2 + 1] = -1;
	}
        remaining_sample_idxs[uii] = kk;
      }
      while (ii < family_size) {
	complete_sample_idxs[complete_sample_idx_ct++] = fis_ptr[(uint32_t)ii];
	ii++;
      }
      qsort(stray_parent_ids, stray_parent_ct, max_pm_id_len, strcmp_casted);
      cur_sample_id = stray_parent_ids;
      ii = 0; // read idx
      jj = 0; // write idx

      // Now filter out all such parents who aren't referenced at least twice.
      while (ii + 1 < stray_parent_ct) {
        if (strcmp(&(stray_parent_ids[ii * max_pm_id_len]), &(stray_parent_ids[(ii + 1) * max_pm_id_len]))) {
	  ii++;
	  continue;
	}
	ii++;
	strcpy(cur_sample_id, &(stray_parent_ids[ii * max_pm_id_len]));
	do {
	  ii++;
        } while (!(strcmp(cur_sample_id, &(stray_parent_ids[ii * max_pm_id_len])) || (ii > stray_parent_ct)));
        cur_sample_id = &(cur_sample_id[max_pm_id_len]);
	jj++;
      }
      stray_parent_ct = jj;

      // Now allocate temporary relatedness table between nonfounders and
      // stray parents with multiple references.
      if (stray_parent_ct) {
        if (wkspace_alloc_d_checked(&tmp_rel_space, (family_size - founder_ct) * stray_parent_ct * sizeof(double))) {
	  return RET_NOMEM;
        }
	tmp_rel_writer = tmp_rel_space;
      }

      // Now fill in remainder of remaining_sample_parent_idxs.
      for (uii = 0; uii < remaining_sample_ct; uii++) {
	jj = remaining_sample_idxs[uii];
	if (remaining_sample_parent_idxs[uii * 2] == -2) {
	  kk = bsearch_str_nl(&(paternal_ids[((uint32_t)jj) * max_paternal_id_len]), stray_parent_ids, max_pm_id_len, stray_parent_ct);
	  if (kk != -1) {
	    kk += unfiltered_sample_ctlm;
	  }
	  remaining_sample_parent_idxs[uii * 2] = kk;
	}
	if (remaining_sample_parent_idxs[uii * 2 + 1] == -2) {
	  kk = bsearch_str_nl(&(maternal_ids[((uint32_t)jj) * max_maternal_id_len]), stray_parent_ids, max_pm_id_len, stray_parent_ct);
	  if (kk != -1) {
	    kk += unfiltered_sample_ctlm;
	  }
	  remaining_sample_parent_idxs[uii * 2 + 1] = kk;
	}
      }
      ullii = ((uint64_t)founder_ct) * (founder_ct - 1);
      while (remaining_sample_ct) {
	sample_idx_write = 0;
	for (uii = 0; uii < remaining_sample_ct; uii++) {
	  kk = remaining_sample_parent_idxs[uii * 2];
	  mm = remaining_sample_parent_idxs[uii * 2 + 1];
	  jj = remaining_sample_idxs[uii];
	  if (((kk == -1) || is_set(processed_samples, kk)) && ((mm == -1) || is_set(processed_samples, mm))) {
	    for (ujj = 0; ujj < founder_ct; ujj++) {
	      // relationship between kk and ujjth founder
	      if ((kk >= (int32_t)unfiltered_sample_ct) || (kk == -1)) {
		dxx = 0.0;
	      } else if (is_set(founder_info, kk)) {
		if (kk == (int32_t)complete_sample_idxs[ujj]) {
		  dxx = 0.5;
		} else {
		  dxx = 0.0;
		}
	      } else {
		ukk = pri_ptr->family_rel_nf_idxs[(uint32_t)kk];
                dxx = 0.5 * rs_ptr[((uint64_t)ukk * (ukk - 1) - ullii) / 2 + ujj];
	      }
	      if (is_set(founder_info, mm)) {
		if (mm == (int32_t)complete_sample_idxs[ujj]) {
		  dxx += 0.5;
		}
	      } else if ((mm != -1) && (mm < (int32_t)unfiltered_sample_ct)) {
		ukk = pri_ptr->family_rel_nf_idxs[(uint32_t)mm];
		dxx += 0.5 * rs_ptr[((uint64_t)ukk * (ukk - 1) - ullii) / 2 + ujj];
	      }
	      *rel_writer++ = dxx;
	    }
	    for (; ujj < complete_sample_idx_ct; ujj++) {
	      if (kk == -1) {
		dxx = 0.0;
	      } else if (kk >= (int32_t)unfiltered_sample_ct) {
		dxx = 0.5 * tmp_rel_space[(ujj - founder_ct) * stray_parent_ct + kk - unfiltered_sample_ctlm];
	      } else if (is_set(founder_info, kk)) {
                dxx = 0.5 * rs_ptr[((uint64_t)ujj * (ujj - 1) - ullii) / 2 + pri_ptr->family_rel_nf_idxs[kk]];
	      } else {
		ukk = pri_ptr->family_rel_nf_idxs[kk];
		if (ukk == ujj) {
		  dxx = 0.5;
		} else if (ukk < ujj) {
		  dxx = 0.5 * rs_ptr[((uint64_t)ujj * (ujj - 1) - ullii) / 2 + ukk];
		} else {
		  dxx = 0.5 * rs_ptr[((uint64_t)ukk * (ukk - 1) - ullii) / 2 + ujj];
		}
	      }
	      if (mm >= (int32_t)unfiltered_sample_ct) {
		dxx += 0.5 * tmp_rel_space[(ujj - founder_ct) * stray_parent_ct + mm - unfiltered_sample_ctlm];
	      } else if (is_set(founder_info, mm)) {
		dxx += 0.5 * rs_ptr[((uint64_t)ujj * (ujj - 1) - ullii) / 2 + pri_ptr->family_rel_nf_idxs[mm]];
	      } else if (mm != -1) {
		ukk = pri_ptr->family_rel_nf_idxs[mm];
		if (ukk == ujj) {
		  dxx += 0.5;
		} else if (ukk < ujj) {
		  dxx += 0.5 * rs_ptr[((uint64_t)ujj * (ujj - 1) - ullii) / 2 + ukk];
		} else {
		  dxx += 0.5 * rs_ptr[((uint64_t)ukk * (ukk - 1) - ullii) / 2 + ujj];
		}
	      }
	      *rel_writer++ = dxx;
	    }
	    for (ujj = 0; ujj < (uintptr_t)stray_parent_ct; ujj++) {
	      if (kk >= (int32_t)unfiltered_sample_ct) {
		if ((uint32_t)kk == ujj + unfiltered_sample_ctlm) {
		  dxx = 0.5;
		} else {
		  dxx = 0.0;
		}
	      } else if (kk == -1) {
                dxx = 0.0;
	      } else {
		ukk = pri_ptr->family_rel_nf_idxs[kk];
		if (ukk < founder_ct) {
		  dxx = 0.0;
		} else {
                  dxx = 0.5 * tmp_rel_space[(ukk - founder_ct) * stray_parent_ct + ujj];
		}
	      }
	      if (mm >= (int32_t)unfiltered_sample_ct) {
		if ((uint32_t)mm == ujj + unfiltered_sample_ctlm) {
		  dxx += 0.5;
		}
	      } else if (mm != -1) {
		ukk = pri_ptr->family_rel_nf_idxs[mm];
		if (ukk >= founder_ct) {
		  dxx += 0.5 * tmp_rel_space[(ukk - founder_ct) * stray_parent_ct + ujj];
		}
	      }
	      *tmp_rel_writer++ = dxx;
	    }
	    pri_ptr->family_rel_nf_idxs[jj] = complete_sample_idx_ct;
	    complete_sample_idxs[complete_sample_idx_ct++] = jj;
	    set_bit(processed_samples, jj);
	  } else {
            remaining_sample_parent_idxs[sample_idx_write * 2] = kk;
	    remaining_sample_parent_idxs[sample_idx_write * 2 + 1] = mm;
	    remaining_sample_idxs[sample_idx_write++] = jj;
	  }
	}
	if (sample_idx_write == remaining_sample_ct) {
	  logprint("Error: Pedigree graph is cyclic.  Check for evidence of time travel abuse in\nyour cohort.\n");
	  return RET_INVALID_FORMAT;
	}
	remaining_sample_ct = sample_idx_write;
      }
      wkspace_reset(wkspace_mark2);
    }
  }
  wkspace_reset(wkspace_mark);
  return 0;
}

int32_t tdt_poo(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, double output_min_p, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct_ax, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, uintptr_t* marker_reverse, uintptr_t unfiltered_sample_ct, uintptr_t* sample_male_include2, uint32_t* trio_nuclear_lookup, uint32_t family_ct, Aperm_info* apip, uint32_t mperm_save, char* sample_ids, uintptr_t max_sample_id_len, Chrom_info* chrom_info_ptr, uint32_t hh_exists, Family_info* fam_ip, uintptr_t* loadbuf, uintptr_t* workbuf, char* textbuf, double* orig_chisq, uint32_t* trio_error_lookup, uintptr_t trio_ct) {
  FILE* outfile = NULL;
  uint64_t mendel_error_ct = 0;
  double pat_a2transmit_recip = 0.0;
  double mat_a1transmit_recip = 0.0;
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t final_mask = get_final_mask(unfiltered_sample_ct);
  uintptr_t marker_uidx = ~ZEROLU;
  uintptr_t markers_done = 0;
  uintptr_t pct = 1;
  uintptr_t pct_thresh = marker_ct_ax / 100;
  uint32_t unfiltered_sample_ctl2m1 = (unfiltered_sample_ct - 1) / BITCT2;
  uint32_t multigen = (fam_ip->mendel_modifier / MENDEL_MULTIGEN) & 1;
  int32_t retval = 0;
  // index bits 0-1: child genotype
  // index bits 2-3: paternal genotype
  // index bits 4-5: maternal genotype
  // entry bit 1: paternal observation?
  // entry bit 9: maternal observation?
  // entry bit 16 & 24: hhh?
  // entry bit 17: paternal A1 transmitted?
  // entry bit 25: maternal A1 transmitted?
  const uint32_t poo_table[] =
    {0, 0, 0, 0,
     0, 0, 0, 0,
     0x20002, 0, 2, 0,
     0, 0, 0, 0,
     0, 0, 0, 0,
     0, 0, 0, 0,
     0, 0, 0, 0,
     0, 0, 0, 0,
     0x2000200, 0, 0x200, 0,
     0, 0, 0, 0,
     0x2020202, 0, 0x1010202, 0x202,
     0, 0, 0x2000200, 0x200,
     0, 0, 0, 0,
     0, 0, 0, 0,
     0, 0, 0x20002, 2};
  char* wptr_start;
  char* wptr;
  char* wptr2;
  uint32_t* lookup_ptr;
  const uint32_t* poo_table_ptr;
  double chisq;
  double pat_a1transmit;
  double cur_a2transmit;
  double dxx;
  uintptr_t ulii;
  uintptr_t uljj;
  uint32_t chrom_fo_idx;
  uint32_t chrom_idx;
  uint32_t chrom_end;
  uint32_t family_idx;
  uint32_t cur_child_ct;
  uint32_t child_idx;
  uint32_t is_x;

  // use a vertical popcount-like strategy here
  uint32_t poo_acc;
  uint32_t poo_acc_ct;
  uint32_t poo_obs_pat_x2;
  uint32_t poo_obs_mat_x2;
  uint32_t poo_pat_a1transmit_x2;
  uint32_t poo_mat_a1transmit_x2;

  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  memcpy(outname_end, ".tdt.poo", 9);
  if (fopen_checked(&outfile, outname, "w")) {
    goto tdt_poo_ret_OPEN_FAIL;
  }
  sprintf(textbuf, " CHR %%%us  A1:A2      T:U_PAT    CHISQ_PAT        P_PAT      T:U_MAT    CHISQ_MAT        P_MAT        Z_POO        P_POO \n", plink_maxsnp);
  fprintf(outfile, textbuf, "SNP");
  fputs("--tdt poo: 0%", stdout);
  fflush(stdout);
  for (chrom_fo_idx = 0; chrom_fo_idx < chrom_info_ptr->chrom_ct; chrom_fo_idx++) {
    chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
    is_x = ((int32_t)chrom_idx == chrom_info_ptr->x_code);
    if ((IS_SET(chrom_info_ptr->haploid_mask, chrom_idx) && (!is_x)) || (((uint32_t)chrom_info_ptr->mt_code) == chrom_idx)) {
      continue;
    }
    chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx + 1];
    uii = next_unset(marker_exclude, chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx], chrom_end);
    if (uii == chrom_end) {
      continue;
    }
    wptr_start = width_force(4, textbuf, chrom_name_write(textbuf, chrom_info_ptr, chrom_idx));
    *wptr_start++ = ' ';
    if (uii != marker_uidx) {
      marker_uidx = uii;
      goto tdt_poo_scan_seek;
    }
    while (1) {
      if (load_raw2(bedfile, loadbuf, unfiltered_sample_ct4, unfiltered_sample_ctl2m1, final_mask)) {
	goto tdt_poo_ret_READ_FAIL;
      }
      if (IS_SET(marker_reverse, marker_uidx)) {
	reverse_loadbuf((unsigned char*)loadbuf, unfiltered_sample_ct);
      }
      if (hh_exists && is_x) {
        hh_reset((unsigned char*)loadbuf, sample_male_include2, unfiltered_sample_ct);
      }
      mendel_error_ct += erase_mendel_errors(unfiltered_sample_ct, loadbuf, workbuf, trio_error_lookup, trio_ct, multigen);
      lookup_ptr = trio_nuclear_lookup;
      poo_acc = 0;
      poo_acc_ct = 0;
      poo_obs_pat_x2 = 0;
      poo_obs_mat_x2 = 0;
      poo_pat_a1transmit_x2 = 0;
      poo_mat_a1transmit_x2 = 0;
      for (family_idx = 0; family_idx < family_ct; family_idx++) {
        uii = *lookup_ptr++;
        ujj = *lookup_ptr++;
        cur_child_ct = *lookup_ptr++;
        ulii = (loadbuf[uii / BITCT2] >> (2 * (uii % BITCT2))) & 3;
        uljj = (loadbuf[ujj / BITCT2] >> (2 * (ujj % BITCT2))) & 3;
        ukk = ulii | (uljj << 2);
	if ((0x4d04 >> ukk) & 1) {
	  // 1+ het parents, no missing
	  // xor doesn't work here since we need to track fathers and mothers
	  // separately
	  poo_table_ptr = &(poo_table[4 * ukk]);
          for (child_idx = 0; child_idx < cur_child_ct; child_idx++) {
            ukk = *lookup_ptr++;
            poo_acc += poo_table_ptr[(loadbuf[ukk / BITCT2] >> (2 * (ukk % BITCT2))) & 3];
	    if (++poo_acc_ct == 127) {
	      // accumulator about to overflow, unpack it
              poo_obs_pat_x2 += (unsigned char)poo_acc;
              poo_obs_mat_x2 += (unsigned char)(poo_acc >> 8);
              poo_pat_a1transmit_x2 += (unsigned char)(poo_acc >> 16);
              poo_mat_a1transmit_x2 += poo_acc >> 24;
              poo_acc = 0;
              poo_acc_ct = 0;
	    }
	  }
	} else {
          lookup_ptr = &(lookup_ptr[cur_child_ct]);
	}
      }
      if (poo_acc_ct) {
	poo_obs_pat_x2 += (unsigned char)poo_acc;
	poo_obs_mat_x2 += (unsigned char)(poo_acc >> 8);
	poo_pat_a1transmit_x2 += (unsigned char)(poo_acc >> 16);
	poo_mat_a1transmit_x2 += poo_acc >> 24;
      }
      wptr = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), wptr_start);
      *wptr++ = ' ';
      wptr2 = strcpyax(wptr, marker_allele_ptrs[2 * marker_uidx], ':');
      wptr2 = strcpya(wptr2, marker_allele_ptrs[2 * marker_uidx + 1]);
      wptr = width_force(6, wptr, wptr2);
      *wptr++ = ' ';
      pat_a1transmit = 0.5 * ((double)poo_pat_a1transmit_x2);
      cur_a2transmit = 0.5 * ((double)(poo_obs_pat_x2 - poo_pat_a1transmit_x2));
      wptr2 = double_g_writewx4x(wptr, pat_a1transmit, 1, ':');
      wptr2 = double_g_writewx4(wptr2, cur_a2transmit, 1);
      wptr = width_force(12, wptr, wptr2);
      *wptr++ = ' ';
      if (poo_obs_pat_x2) {
	pat_a2transmit_recip = 1.0 / cur_a2transmit;
	dxx = pat_a1transmit - cur_a2transmit;
	chisq = dxx * dxx / (pat_a1transmit + cur_a2transmit);
	wptr = double_g_writewx4x(wptr, chisq, 12, ' ');
	wptr = double_g_writewx4(wptr, chiprob_p(chisq, 1), 12);
      } else {
	wptr = memcpya(wptr, "          NA           NA", 25);
      }
      *wptr++ = ' ';
      dxx = 0.5 * ((double)poo_mat_a1transmit_x2);
      cur_a2transmit = 0.5 * ((double)(poo_obs_mat_x2 - poo_mat_a1transmit_x2));
      wptr2 = double_g_writewx4x(wptr, dxx, 1, ':');
      wptr2 = double_g_writewx4(wptr2, cur_a2transmit, 1);
      wptr = width_force(12, wptr, wptr2);
      *wptr++ = ' ';
      if (poo_obs_mat_x2) {
	mat_a1transmit_recip = 1.0 / dxx;
	chisq = dxx - cur_a2transmit;
	chisq = chisq * chisq / (dxx + cur_a2transmit);
	wptr = double_g_writewx4x(wptr, chisq, 12, ' ');
	wptr = double_g_writewx4(wptr, chiprob_p(chisq, 1), 12);
      } else {
	wptr = memcpya(wptr, "          NA           NA", 25);
      }
      *wptr++ = ' ';
      if (poo_pat_a1transmit_x2 && poo_mat_a1transmit_x2 && (poo_obs_pat_x2 > poo_pat_a1transmit_x2) && (poo_obs_mat_x2 > poo_mat_a1transmit_x2)) {
	// Z-score
	dxx = (log(pat_a1transmit * pat_a2transmit_recip * mat_a1transmit_recip * cur_a2transmit) / sqrt(1.0 / pat_a1transmit + pat_a2transmit_recip + mat_a1transmit_recip + 1.0 / cur_a2transmit));

        wptr = double_g_writewx4x(wptr, dxx, 12, ' ');
	if (orig_chisq) {
	  // todo: --pat/--mat support
	  orig_chisq[markers_done] = dxx * dxx;
	}
	dxx = normdist(-fabs(dxx)) * 2;
	wptr = double_g_writewx4(wptr, MAXV(dxx, output_min_p), 12);
      } else {
	wptr = memcpya(wptr, "          NA           NA", 25);
	if (orig_chisq) {
	  orig_chisq[markers_done] = -1;
	}
      }
      wptr = memcpya(wptr, " \n", 2);
      // 1.07 implementation does not support pfilter; we might want to add it
      if (fwrite_checked(textbuf, wptr - textbuf, outfile)) {
	goto tdt_poo_ret_WRITE_FAIL;
      }
      if (++markers_done >= pct_thresh) {
	if (pct > 10) {
	  putchar('\b');
	}
	pct = (markers_done * 100LLU) / marker_ct_ax;
	if (pct < 100) {
	  printf("\b\b%" PRIuPTR "%%", pct);
	  fflush(stdout);
	  pct_thresh = ((++pct) * ((uint64_t)marker_ct_ax)) / 100;
	}
      }
      if (++marker_uidx == chrom_end) {
	break;
      }
      if (IS_SET(marker_exclude, marker_uidx)) {
	marker_uidx = next_unset_ul(marker_exclude, marker_uidx, chrom_end);
      tdt_poo_scan_seek:
	if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	  goto tdt_poo_ret_READ_FAIL;
	}
	if (marker_uidx == chrom_end) {
	  break;
	}
      }
    }
  }
  if (fclose_null(&outfile)) {
    goto tdt_poo_ret_WRITE_FAIL;
  }
  putchar('\r');
  LOGPRINTF("--tdt poo: Report written to %s .\n", outname);
  while (0) {
  tdt_poo_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  tdt_poo_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  tdt_poo_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  fclose_cond(outfile);
  return retval;
}

int32_t tdt(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, double ci_size, double ci_zt, double pfilter, double output_min_p, uint32_t mtest_adjust, double adjust_lambda, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, uint32_t* marker_pos, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, uintptr_t* marker_reverse, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, Aperm_info* apip, uint32_t mperm_save, uintptr_t* pheno_nm, uintptr_t* pheno_c, uintptr_t* founder_info, uintptr_t* sex_nm, uintptr_t* sex_male, char* sample_ids, uintptr_t max_sample_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, Chrom_info* chrom_info_ptr, uint32_t hh_exists, Family_info* fam_ip) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* outfile = NULL;
  char* textbuf = tbuf;
  double* orig_chisq = NULL; // pval if exact test
  uint64_t last_parents = 0;
  // uint64_t mendel_error_ct = 0;
  double chisq = 0;
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t unfiltered_sample_ctl2 = (unfiltered_sample_ct + (BITCT2 - 1)) / BITCT2;
  uintptr_t unfiltered_sample_ctp1l2 = 1 + (unfiltered_sample_ct / BITCT2);
  uintptr_t final_mask = get_final_mask(unfiltered_sample_ct);
  uintptr_t marker_uidx = ~ZEROLU;
  uintptr_t markers_done = 0;
  uintptr_t pct = 1;
  uint32_t unfiltered_sample_ctl2m1 = (unfiltered_sample_ct - 1) / BITCT2;
  uint32_t multigen = (fam_ip->mendel_modifier / MENDEL_MULTIGEN) & 1;
  uint32_t display_ci = (ci_size > 0);
  uint32_t is_exact = fam_ip->tdt_modifier & TDT_EXACT;
  uint32_t is_midp = fam_ip->tdt_modifier & TDT_MIDP;
  uint32_t poo_test = fam_ip->tdt_modifier & TDT_POO;
  // uint32_t perm_count = fam_ip->tdt_modifier & TDT_PERM_COUNT;
  uint32_t case_trio_ct = 0;
  uint32_t is_discordant = 0;
  uint32_t discord_exists = 0;
  uint32_t cur_child_ct = 0xffffffffU;
  int32_t retval = 0;
  // index bits 0-1: child genotype
  // index bits 2-3: XOR of parental genotypes
  // entry bits 0-15: observation count increment
  // entry bits 16-31: A1 transmission increment
  // het + hom A1 parents, hom A2 child (or vice versa) needs to be in the
  // table because of Xchr; fortunately, that's actually all that's necessary
  // to handle X properly.
  const uint32_t tdt_table[] =
    {0x20002, 0, 0x10002, 2,
     0x10001, 0, 0x10001, 1,
     0x10001, 0, 1, 1};
  // index bits 0-1: case genotype
  // index bits 2-3: ctrl genotype
  // entry bit 0: single-observation?
  // entry bit 8: double-observation?
  // entry bit 16: single-obs = case A2 excess?
  // entry bit 24: double-obs = case A2 excess?
  const uint32_t parentdt_table[] =
    {0, 0, 1, 0x100,
     0, 0, 0, 0,
     0x10001, 0, 0, 1,
     0x1000100, 0, 0x10001, 0};
  char* fids;
  char* iids;
  char* wptr_start;
  char* wptr;
  char* wptr2;
  uint64_t* family_list;
  uint64_t* trio_list;
  uintptr_t* loadbuf;
  uintptr_t* workbuf;
  uintptr_t* sample_male_include2;
  uintptr_t* marker_exclude_tmp;
  uint32_t* trio_error_lookup;

  // sequences of variable-length nuclear family records.
  // [0]: father uidx
  // [1]: mother uidx
  // [2]: number of children (n)
  // [3..(2+n)]: child uidxs
  // [3+n]: father uidx for next record
  // etc.
  uint32_t* trio_nuclear_lookup;
  uint32_t* lookup_ptr;
  uint32_t* marker_idx_to_uidx;
  const uint32_t* tdt_table_ptr;
  uint64_t cur_trio;
  uint64_t cur_parents;
  uintptr_t max_fid_len;
  uintptr_t max_iid_len;
  uintptr_t trio_ct;
  uintptr_t trio_idx;
  uintptr_t pct_thresh;
  uintptr_t ulii;
  uintptr_t uljj;
  double pval;
  double odds_ratio;
  double untransmitted_recip;
  double dxx;
  uint32_t chrom_fo_idx;
  uint32_t chrom_idx;
  uint32_t chrom_end;
  uint32_t is_x;
  uint32_t family_ct;
  uint32_t family_idx;
  uint32_t child_idx;

  // use a vertical popcount-like strategy here
  uint32_t parentdt_acc;
  uint32_t parentdt_acc_ct;
  uint32_t parentdt_obs_ct1;
  uint32_t parentdt_obs_ct2;
  uint32_t parentdt_case_a2_excess1;
  uint32_t parentdt_case_a2_excess2;

  uint32_t tdt_obs_ct;
  uint32_t tdt_a1_trans_ct;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  marker_ct -= count_non_autosomal_markers(chrom_info_ptr, marker_exclude, 0, 1);
  if ((!marker_ct) || is_set(chrom_info_ptr->haploid_mask, 0)) {
    logprint("Warning: Skipping --tdt since there is no autosomal or Xchr data.\n");
    goto tdt_ret_1;
  }
  retval = get_trios_and_families(unfiltered_sample_ct, sample_exclude, sample_ct, founder_info, sex_nm, sex_male, sample_ids, max_sample_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, &fids, &max_fid_len, &iids, &max_iid_len, &family_list, &family_ct, &trio_list, &trio_ct, &trio_error_lookup, 0, multigen);
  if (retval) {
    goto tdt_ret_1;
  }
  if (!trio_ct) {
    logprint("Warning: Skipping --tdt since there are no trios.\n");
    goto tdt_ret_1;
  }
  // now assemble list of nuclear families with at least one case child
  if (wkspace_alloc_ui_checked(&trio_nuclear_lookup, (3 * family_ct + trio_ct) * sizeof(int32_t))) {
    goto tdt_ret_NOMEM;
  }
  lookup_ptr = trio_nuclear_lookup;
  family_ct = 0;
  // note that trio_list is already sorted by family, so we can just traverse
  // it in order.
  for (trio_idx = 0; trio_idx < trio_ct; trio_idx++) {
    cur_trio = trio_list[trio_idx];
    ukk = (uint32_t)cur_trio;
    cur_parents = family_list[(uintptr_t)(cur_trio >> 32)];
    if (cur_parents != last_parents) {
      // cur_child_ct initialized to maxint instead of 0 to prevent this from
      // triggering on the first family
      if (!cur_child_ct) {
        if (!is_discordant) {
	  // don't need to track previous family after all
          lookup_ptr = &(lookup_ptr[-3]);
	} else {
	  family_ct++;
	}
      } else if (cur_child_ct != 0xffffffffU) {
	case_trio_ct += cur_child_ct;
        lookup_ptr[-(1 + (int32_t)cur_child_ct)] = cur_child_ct | (is_discordant << 31);
	family_ct++;
      }
      last_parents = cur_parents;
      uii = (uint32_t)cur_parents;
      ujj = (uint32_t)(cur_parents >> 32);
      if ((uii < unfiltered_sample_ct) && (ujj < unfiltered_sample_ct)) {
	umm = is_set(pheno_c, uii);
	is_discordant = (!poo_test) && is_set(pheno_nm, uii) && is_set(pheno_nm, ujj) && (umm ^ is_set(pheno_c, ujj));
	if (!is_discordant) {
	  *lookup_ptr++ = uii;
	  *lookup_ptr++ = ujj;
	} else {
	  discord_exists = 1;
	  // safe to reverse order of parents so that case parent comes first
	  if (umm) {
	    *lookup_ptr++ = uii;
	    *lookup_ptr++ = ujj;
	  } else {
	    *lookup_ptr++ = ujj;
	    *lookup_ptr++ = uii;
	  }
	  *lookup_ptr = 0x80000000;
	}
	lookup_ptr++;
	cur_child_ct = 0;
	goto tdt_get_case_children;
      } else {
	cur_child_ct = 0xffffffffU;
      }
    } else if (cur_child_ct != 0xffffffffU) {
    tdt_get_case_children:
      if (is_set(pheno_c, ukk)) {
	cur_child_ct++;
	*lookup_ptr++ = ukk;
      }
    }
  }
  if (!cur_child_ct) {
    if (!is_discordant) {
      lookup_ptr = &(lookup_ptr[-3]);
    } else {
      family_ct++;
    }
  } else if (cur_child_ct != 0xffffffffU) {
    case_trio_ct += cur_child_ct;
    lookup_ptr[-(1 + (int32_t)cur_child_ct)] = cur_child_ct | (is_discordant << 31);
    family_ct++;
  }
  if (lookup_ptr == trio_nuclear_lookup) {
    LOGPRINTF("Warning: Skipping --tdt%s since there are no trios with an affected child%s.\n", poo_test? " poo" : "", poo_test? "" : ", and no\ndiscordant parent pairs");
    goto tdt_ret_1;
  }
  wkspace_shrink_top(trio_nuclear_lookup, ((uintptr_t)(lookup_ptr - trio_nuclear_lookup)) * sizeof(int32_t));

  if (mtest_adjust) {
    if (wkspace_alloc_d_checked(&orig_chisq, marker_ct * sizeof(double))) {
      goto tdt_ret_NOMEM;
    }
  }
  if (wkspace_alloc_ul_checked(&loadbuf, unfiltered_sample_ctl2 * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&workbuf, unfiltered_sample_ctp1l2 * sizeof(intptr_t))) {
    goto tdt_ret_NOMEM;
  }
  loadbuf[unfiltered_sample_ctl2 - 1] = 0;
  workbuf[unfiltered_sample_ctp1l2 - 1] = 0;
  hh_exists &= XMHH_EXISTS;
  if (alloc_raw_haploid_filters(unfiltered_sample_ct, hh_exists, 1, sample_exclude, sex_male, NULL, &sample_male_include2)) {
    goto tdt_ret_NOMEM;
  }
  if (fam_ip->tdt_modifier & (TDT_PERM | TDT_MPERM)) {
    logprint("Error: --tdt permutation tests are currently under development.\n");
    retval = RET_CALC_NOT_YET_SUPPORTED;
    goto tdt_ret_1;
  }
  ulii = 2 * max_marker_allele_len + plink_maxsnp + MAX_ID_LEN + 256;
  if (ulii > MAXLINELEN) {
    if (wkspace_alloc_c_checked(&textbuf, ulii)) {
      goto tdt_ret_NOMEM;
    }
  }
  if (poo_test) {
    retval = tdt_poo(threads, bedfile, bed_offset, outname, outname_end, output_min_p, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, plink_maxsnp, marker_allele_ptrs, max_marker_allele_len, marker_reverse, unfiltered_sample_ct, sample_male_include2, trio_nuclear_lookup, family_ct, apip, mperm_save, sample_ids, max_sample_id_len, chrom_info_ptr, hh_exists, fam_ip, loadbuf, workbuf, textbuf, orig_chisq, trio_error_lookup, trio_ct);
    if (retval) {
      goto tdt_ret_1;
    }
    if (mtest_adjust) {
      outname_end[4] = '\0';
      goto tdt_multcomp;
    }
    goto tdt_ret_1;
  }
  pct_thresh = marker_ct / 100;
  memcpy(outname_end, ".tdt", 5);
  if (fopen_checked(&outfile, outname, "w")) {
    goto tdt_ret_OPEN_FAIL;
  }
  sprintf(textbuf, " CHR %%%us           BP  A1  A2      T      U           OR ", plink_maxsnp);
  fprintf(outfile, textbuf, "SNP");
  if (display_ci) {
    uii = (uint32_t)((int32_t)(ci_size * 100));
    if (uii >= 10) {
      fprintf(outfile, "         L%u          U%u ", uii, uii);
    } else {
      fprintf(outfile, "          L%u           U%u ", uii, uii);
    }
  }
  if (!is_exact) {
    fputs("       CHISQ ", outfile);
  }
  fputs("           P ", outfile);
  if (discord_exists) {
    fputs("     A:U_PAR    CHISQ_PAR        P_PAR    CHISQ_COM        P_COM ", outfile);
  }
  if (putc_checked('\n', outfile)) {
    goto tdt_ret_WRITE_FAIL;
  }
  fputs("--tdt: 0%", stdout);
  fflush(stdout);
  for (chrom_fo_idx = 0; chrom_fo_idx < chrom_info_ptr->chrom_ct; chrom_fo_idx++) {
    chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
    is_x = ((int32_t)chrom_idx == chrom_info_ptr->x_code);
    if ((IS_SET(chrom_info_ptr->haploid_mask, chrom_idx) && (!is_x)) || (((uint32_t)chrom_info_ptr->mt_code) == chrom_idx)) {
      continue;
    }
    chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx + 1];
    uii = next_unset(marker_exclude, chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx], chrom_end);
    if (uii == chrom_end) {
      continue;
    }
    wptr_start = width_force(4, textbuf, chrom_name_write(textbuf, chrom_info_ptr, chrom_idx));
    *wptr_start++ = ' ';
    if (uii != marker_uidx) {
      marker_uidx = uii;
      goto tdt_scan_seek;
    }
    while (1) {
      if (load_raw2(bedfile, loadbuf, unfiltered_sample_ct4, unfiltered_sample_ctl2m1, final_mask)) {
	goto tdt_ret_READ_FAIL;
      }
      if (IS_SET(marker_reverse, marker_uidx)) {
	reverse_loadbuf((unsigned char*)loadbuf, unfiltered_sample_ct);
      }
      if (hh_exists && is_x) {
	hh_reset((unsigned char*)loadbuf, sample_male_include2, unfiltered_sample_ct);
      }
      // 1. iterate through all trios, setting Mendel errors to missing
      // 2. iterate through trio_nuclear_lookup:
      //    a. if either parent isn't genotyped, skip the family
      //    b. if parents are discordant, increment parenTDT counts
      //    c. if at least one het parent, iterate through genotyped children,
      //       incrementing regular TDT counts.
      // mendel_error_ct += erase_mendel_errors(unfiltered_sample_ct, loadbuf, workbuf, trio_error_lookup, trio_ct, multigen);
      erase_mendel_errors(unfiltered_sample_ct, loadbuf, workbuf, trio_error_lookup, trio_ct, multigen);
      lookup_ptr = trio_nuclear_lookup;
      parentdt_acc = 0;
      parentdt_acc_ct = 0;
      parentdt_obs_ct1 = 0;
      parentdt_obs_ct2 = 0;
      parentdt_case_a2_excess1 = 0;
      parentdt_case_a2_excess2 = 0;
      tdt_obs_ct = 0;
      tdt_a1_trans_ct = 0;
      for (family_idx = 0; family_idx < family_ct; family_idx++) {
	uii = *lookup_ptr++;
	ujj = *lookup_ptr++;
	cur_child_ct = *lookup_ptr++;
        ulii = (loadbuf[uii / BITCT2] >> (2 * (uii % BITCT2))) & 3;
        uljj = (loadbuf[ujj / BITCT2] >> (2 * (ujj % BITCT2))) & 3;
        ukk = ulii | (uljj << 2);
	if (cur_child_ct & 0x80000000U) {
          // discordant
	  cur_child_ct ^= 0x80000000U;
	  if ((0x22f2 >> ukk) & 1) {
	    // at least one missing parent, skip both parenTDT and regular TDT
	    lookup_ptr = &(lookup_ptr[cur_child_ct]);
	    continue;
	  }
	  parentdt_acc += parentdt_table[ukk];
	  if (++parentdt_acc_ct == 255) {
	    // accumulator about to overflow, unpack it
	    parentdt_obs_ct1 += (unsigned char)parentdt_acc;
	    parentdt_obs_ct2 += (unsigned char)(parentdt_acc >> 8);
	    parentdt_case_a2_excess1 += (unsigned char)(parentdt_acc >> 16);
	    parentdt_case_a2_excess2 += parentdt_acc >> 24;
	    parentdt_acc = 0;
	    parentdt_acc_ct = 0;
	  }
	  if (!cur_child_ct) {
	    continue;
	  }
	}
        if ((0x4d04 >> ukk) & 1) {
	  // 1+ het parents, no missing
          tdt_table_ptr = &(tdt_table[4 * (ulii ^ uljj)]);
          for (child_idx = 0; child_idx < cur_child_ct; child_idx++) {
            ukk = *lookup_ptr++;
	    umm = tdt_table_ptr[(loadbuf[ukk / BITCT2] >> (2 * (ukk % BITCT2))) & 3];
	    tdt_obs_ct += (uint16_t)umm;
            tdt_a1_trans_ct += umm >> 16;
	  }
	} else {
	  lookup_ptr = &(lookup_ptr[cur_child_ct]);
	}
      }
      if (is_exact) {
	pval = binom_2sided(tdt_a1_trans_ct, tdt_obs_ct, is_midp);
	if (mtest_adjust) {
	  orig_chisq[markers_done] = pval;
	}
      } else if (!tdt_obs_ct) {
	pval = -9;
	if (mtest_adjust) {
	  orig_chisq[markers_done] = -9;
	}
      } else {
	dxx = (double)(((int32_t)tdt_obs_ct) - 2 * ((int32_t)tdt_a1_trans_ct));
	chisq = dxx * dxx / ((double)((int32_t)tdt_obs_ct));
	if (mtest_adjust) {
	  orig_chisq[markers_done] = chisq;
	}
	pval = chiprob_p(chisq, 1);
      }
      if ((pfilter == 2.0) || ((pval <= pfilter) && (pval >= 0.0))) {
	wptr = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), wptr_start);
	wptr = memseta(wptr, 32, 3);
	wptr = uint32_writew10x(wptr, marker_pos[marker_uidx], ' ');
	wptr = fw_strcpy(3, marker_allele_ptrs[2 * marker_uidx], wptr);
	*wptr++ = ' ';
	wptr = fw_strcpy(3, marker_allele_ptrs[2 * marker_uidx + 1], wptr);
	*wptr++ = ' ';
	wptr = uint32_writew6x(wptr, tdt_a1_trans_ct, ' ');
	uii = tdt_obs_ct - tdt_a1_trans_ct; // untransmitted
	wptr = uint32_writew6x(wptr, uii, ' ');
	if (uii) {
	  untransmitted_recip = 1.0 / ((double)((int32_t)uii));
	  dxx = (double)((int32_t)tdt_a1_trans_ct);
	  odds_ratio = dxx * untransmitted_recip;
	  wptr = double_g_writewx4x(wptr, odds_ratio, 12, ' ');
	  if (display_ci) {
	    odds_ratio = log(odds_ratio);
	    dxx = ci_zt * sqrt(1.0 / dxx + untransmitted_recip);
	    wptr = double_g_writewx4x(wptr, exp(odds_ratio - dxx), 12, ' ');
	    wptr = double_g_writewx4x(wptr, exp(odds_ratio + dxx), 12, ' ');
	  }
	} else {
	  wptr = memcpya(wptr, "          NA ", 13);
	  if (display_ci) {
	    wptr = memcpya(wptr, "          NA           NA ", 26);
	  }
	}
        if (is_exact) {
	  wptr = double_g_writewx4x(wptr, MAXV(pval, output_min_p), 12, ' ');
	} else {
	  if (pval >= 0) {
	    wptr = double_g_writewx4x(wptr, chisq, 12, ' ');
            wptr = double_g_writewx4x(wptr, MAXV(pval, output_min_p), 12, ' ');
	  } else {
	    wptr = memcpya(wptr, "          NA           NA ", 26);
	  }
	}
	if (discord_exists) {
	  if (parentdt_acc_ct) {
	    parentdt_obs_ct1 += (unsigned char)parentdt_acc;
	    parentdt_obs_ct2 += (unsigned char)(parentdt_acc >> 8);
	    parentdt_case_a2_excess1 += (unsigned char)(parentdt_acc >> 16);
	    parentdt_case_a2_excess2 += parentdt_acc >> 24;
	  }
	  uii = parentdt_case_a2_excess1 + 2 * parentdt_case_a2_excess2;
	  ujj = parentdt_obs_ct1 + 2 * parentdt_obs_ct2;
	  wptr2 = uint32_writex(wptr, uii, ':');
	  wptr2 = uint32_write(wptr2, ujj - uii);
          wptr = width_force(12, wptr, wptr2);
          *wptr++ = ' ';
	  // No exact test for now since we're dealing with a sum of step-1 and
	  // step-2 binomial distributions.  If anyone wants it, though, it can
	  // be implemented as follows:
	  // 1. Preallocate a size-O(n) buffer; it'll actually be useful here.
	  // 2. Compute actual (not just relative) likelihoods for both
	  //    component binomial distributions, and save them in the buffer.
	  //    Note futility thresholds.  (It may be worthwhile to cache the
	  //    most common distributions, especially since that also enables a
	  //    faster binom_2sided() variant.)
	  // 3. The likelihood of sum k is then
	  //      p1(k) * p2(0) + p1(k-2) * p2(1) + ... + p1(0) * p2(k/2)
	  //    if k is even, and
	  //      p1(k) * p2(0) + ... + p1(1) * p2((k-1)/2)
	  //    if k is odd.
	  //    These sums are unimodal, so evaluation should start from the
	  //    center and employ the usual early-termination logic to bring
	  //    complexity down to O(sqrt(n)).  (To check: is the distribution
	  //    of final likelihoods also unimodal?  Poisson binomial
	  //    distribution is still unimodal, but that's variable p, not
	  //    variable step size.)
	  //    Because there's no unknown scaling factor to worry about, we
	  //    can simply guess where the start of the other tail is and
	  //    examine adjacent likelihoods until we've verified where the
	  //    crossover point is.
	  //    Then, we can use a heuristic to decide between summing tail or
	  //    summing central likelihoods for determination of the final
	  //    p-value.  (Due to catastrophic cancellation, we should always
	  //    use the tail approach if our point likelihood is smaller than,
	  //    say, 10^{-6}.  But there are at least a few common cases, e.g.
	  //    pval=1, where the central approach is much better.)
	  //    General-case computational complexity is O(n), and the most
	  //    common cases are O(sqrt(n)).
	  if (!ujj) {
	    wptr = memcpya(wptr, "          NA           NA", 25);
	  } else {
	    dxx = (double)(((int32_t)ujj) - 2 * ((int32_t)uii));
	    chisq = dxx * dxx / ((double)((intptr_t)((uintptr_t)(ujj + 2 * parentdt_obs_ct2))));
	    wptr = double_g_writewx4x(wptr, chisq, 12, ' ');
	    dxx = chiprob_p(chisq, 1);
	    wptr = double_g_writewx4(wptr, MAXV(dxx, output_min_p), 12);
	  }
	  *wptr++ = ' ';
	  uii += tdt_a1_trans_ct;
	  ujj += tdt_obs_ct;
	  if (!ujj) {
	    wptr = memcpya(wptr, "          NA           NA", 25);
	  } else {
	    // okay, these casting choices, which are only relevant if there
	    // are close to a billion samples, are silly, but may as well make
	    // them consistent with each other.
            dxx = (double)((intptr_t)((uintptr_t)ujj) - 2 * ((intptr_t)((uintptr_t)uii)));
	    chisq = dxx * dxx / ((double)((intptr_t)(((uintptr_t)ujj) + 2 * parentdt_obs_ct2)));
	    wptr = double_g_writewx4x(wptr, chisq, 12, ' ');
	    dxx = chiprob_p(chisq, 1);
	    wptr = double_g_writewx4(wptr, MAXV(dxx, output_min_p), 12);
	  }
	}
	wptr = memcpya(wptr, " \n", 2);
        if (fwrite_checked(textbuf, wptr - textbuf, outfile)) {
	  goto tdt_ret_WRITE_FAIL;
	}
      }

      if (++markers_done >= pct_thresh) {
        if (pct > 10) {
	  putchar('\b');
	}
	pct = (markers_done * 100LLU) / marker_ct;
        if (pct < 100) {
	  printf("\b\b%" PRIuPTR "%%", pct);
          fflush(stdout);
          pct_thresh = ((++pct) * ((uint64_t)marker_ct)) / 100;
	}
      }
      if (++marker_uidx == chrom_end) {
	break;
      }
      if (IS_SET(marker_exclude, marker_uidx)) {
        marker_uidx = next_unset_ul(marker_exclude, marker_uidx, chrom_end);
      tdt_scan_seek:
	if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	  goto tdt_ret_READ_FAIL;
	}
	if (marker_uidx == chrom_end) {
	  break;
	}
      }
    }
  }
  putchar('\r');
  LOGPRINTF("--tdt: Report written to %s .\n", outname);
  if (mtest_adjust) {
  tdt_multcomp:
    ulii = (unfiltered_marker_ct + (BITCT - 1)) / BITCT;
    if (wkspace_alloc_ui_checked(&marker_idx_to_uidx, marker_ct * sizeof(int32_t)) ||
        wkspace_alloc_ul_checked(&marker_exclude_tmp, ulii * sizeof(intptr_t))) {
      goto tdt_ret_NOMEM;
    }
    // need a custom marker_exclude that's set at Y/haploid/MT
    memcpy(marker_exclude_tmp, marker_exclude, ulii * sizeof(intptr_t));
    for (chrom_fo_idx = 0; chrom_fo_idx < chrom_info_ptr->chrom_ct; chrom_fo_idx++) {
      chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
      if ((is_set(chrom_info_ptr->haploid_mask, chrom_idx) && ((int32_t)chrom_idx != chrom_info_ptr->x_code)) || ((int32_t)chrom_idx == chrom_info_ptr->mt_code)) {
	uii = chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx];
	fill_bits(marker_exclude_tmp, uii, chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx + 1] - uii);
      }
    }
    fill_idx_to_uidx(marker_exclude_tmp, unfiltered_marker_ct, marker_ct, marker_idx_to_uidx);
    retval = multcomp(outname, outname_end, marker_idx_to_uidx, marker_ct, marker_ids, max_marker_id_len, plink_maxsnp, chrom_info_ptr, is_exact? NULL : orig_chisq, pfilter, output_min_p, mtest_adjust, 0, adjust_lambda, NULL, is_exact? orig_chisq : NULL);
    if (retval) {
      goto tdt_ret_1;
    }
  }
  while (0) {
  tdt_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  tdt_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  tdt_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  tdt_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
 tdt_ret_1:
  wkspace_reset(wkspace_mark);
  fclose_cond(outfile);
  return retval;
}

int32_t get_sibship_info(uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, uintptr_t* pheno_nm, uintptr_t* founder_info, char* sample_ids, uintptr_t max_sample_id_len, uintptr_t max_fid_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uint64_t* family_list, uint64_t* trio_list, uint32_t family_ct, uintptr_t trio_ct, uint32_t test_type, uintptr_t** lm_eligible_ptr, uintptr_t** lm_within2_founder_ptr, uint32_t** fs_starts_ptr, uint32_t** fss_contents_ptr, uint32_t** sample_lm_to_fss_idx_ptr, uint32_t* fs_ct_ptr, uint32_t* lm_ct_ptr, uint32_t* singleton_ct_ptr) {
  // On top of get_trios_and_families()'s return values, we need the following
  // information for the main dfam() and qfam() loops:
  // 1. sample idx -> family/sibship idx array
  // 2. fs_starts[]/fs_contents[] arrays describing family/sibship idx ->
  //    sample idxs mapping.
  // We may as well sort size-1 sibships/singleton founders to the end; this
  // lets us get away with a smaller fs_starts[] array and a faster loop.
  // There is also some qfam-specific initialization here (e.g. a divorcee with
  // children from two different spouses may be excluded from the linear
  // model).  test_type is zero for dfam and nonzero for qfam.
  uintptr_t unfiltered_sample_ctl = (unfiltered_sample_ct + (BITCT - 1)) / BITCT;
  uintptr_t sample_ctl = (sample_ct + (BITCT - 1)) / BITCT;
  uintptr_t max_merged_id_len = max_fid_len + max_paternal_id_len + max_maternal_id_len + sizeof(int32_t);
  uintptr_t trio_idx = 0;
  uintptr_t topsize = 0;
  uintptr_t* tmp_within2_founder = NULL;
  uintptr_t* lm_within2_founder = NULL;
  uintptr_t* lm_eligible = NULL;
  uint32_t is_within2 = (test_type == QFAM_WITHIN2);
  uint32_t family_idx = 0;
  uint32_t fssc_idx = 0;
  int32_t retval = 0;
  char* merged_ids;
  char* bufptr;
  char* bufptr2;
  char* bufptr3;
  uintptr_t* not_in_family;
  uintptr_t* ulptr;
  uintptr_t* ulptr2;
  uint32_t* sample_uidx_to_idx;
  uint32_t* sample_to_fss_idx;
  uint32_t* sample_lm_to_fss_idx;
  uint32_t* fs_starts;
  uint32_t* fss_contents;
  uintptr_t topsize_bak;
  uintptr_t cur_sample_ct;
  uintptr_t sample_uidx;
  uintptr_t sample_idx;
  uintptr_t ulii;
  uint64_t ullii;
  uint32_t slen;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  if (test_type) {
    if (is_within2) {
      if (wkspace_alloc_ul_checked(&lm_within2_founder, sample_ctl * sizeof(intptr_t))) {
	goto get_sibship_info_ret_NOMEM2;
      }
    }
    if (wkspace_alloc_ul_checked(&lm_eligible, sample_ctl * sizeof(intptr_t))) {
      goto get_sibship_info_ret_NOMEM;
    }
  }
  if (test_type) {
    // shrink later
    if (wkspace_alloc_ui_checked(&fss_contents, (sample_ct + 2 * family_ct) * sizeof(int32_t))) {
      goto get_sibship_info_ret_NOMEM;
    }
    // this is the equivalent of PLINK 1.07's family pointers
    sample_to_fss_idx = (uint32_t*)top_alloc(&topsize, sample_ct * sizeof(int32_t));
    if (!sample_to_fss_idx) {
      goto get_sibship_info_ret_NOMEM;
    }
  } else {
    if (wkspace_alloc_ui_checked(&sample_to_fss_idx, sample_ct * sizeof(int32_t))) {
      goto get_sibship_info_ret_NOMEM;
    }
    // shrink later
    if (wkspace_alloc_ui_checked(&fss_contents, (sample_ct + 2 * family_ct) * sizeof(int32_t))) {
      goto get_sibship_info_ret_NOMEM;
    }
  }
  topsize_bak = topsize;
  not_in_family = (uintptr_t*)top_alloc(&topsize, unfiltered_sample_ctl * sizeof(intptr_t));
  if (!not_in_family) {
    goto get_sibship_info_ret_NOMEM;
  }

  // Temporary bitfields used to track which parents are (i) part of multiple
  // families, and (ii) not a child in any of them.  To ensure results are not
  // dependent on the order of samples in the dataset, we now exclude these
  // parents from the permutation.  (todo: compute an average in this case
  // instead?)
  sample_uidx_to_idx = (uint32_t*)top_alloc(&topsize, unfiltered_sample_ct * sizeof(int32_t));
  if (!sample_uidx_to_idx) {
    goto get_sibship_info_ret_NOMEM2;
  }
  ulptr = (uintptr_t*)top_alloc(&topsize, unfiltered_sample_ctl * sizeof(intptr_t));
  if (!ulptr) {
    goto get_sibship_info_ret_NOMEM2;
  }
  ulii = topsize;
  ulptr2 = (uintptr_t*)top_alloc(&topsize, unfiltered_sample_ctl * sizeof(intptr_t));
  if (!ulptr2) {
    goto get_sibship_info_ret_NOMEM2;
  }
  if (is_within2) {
    tmp_within2_founder = (uintptr_t*)top_alloc(&topsize, unfiltered_sample_ctl * sizeof(intptr_t));
    if (!tmp_within2_founder) {
      goto get_sibship_info_ret_NOMEM2;
    }
    fill_ulong_zero(tmp_within2_founder, unfiltered_sample_ctl);
  }

  bitfield_exclude_to_include(sample_exclude, not_in_family, unfiltered_sample_ct);
  fill_uint_one(sample_to_fss_idx, sample_ct);
  fill_uidx_to_idx(sample_exclude, unfiltered_sample_ct, sample_ct, sample_uidx_to_idx);
  fill_ulong_zero(ulptr, unfiltered_sample_ctl); // is a double-parent
  fill_ulong_zero(ulptr2, unfiltered_sample_ctl); // is a child
  if (family_ct) {
    while (1) {
      ullii = family_list[family_idx];
      // uii, ukk = unfiltered idxs of parents
      // ujj, umm = filtered idxs
      uii = (uint32_t)ullii;
      ujj = sample_uidx_to_idx[uii];
      fss_contents[fssc_idx++] = ujj;
      ukk = (uint32_t)(ullii >> 32);
      umm = sample_uidx_to_idx[ukk];
      if (is_within2) {
	if (is_set(pheno_nm, uii) && is_set(pheno_nm, ukk)) {
	  set_bit(tmp_within2_founder, uii);
	  set_bit(tmp_within2_founder, ukk);
	}
      }
      if (is_set(not_in_family, uii)) {
	if (sample_to_fss_idx[ujj] == 0xffffffffU) {
	  sample_to_fss_idx[ujj] = family_idx;
	}
	clear_bit(not_in_family, uii);
      } else {
	set_bit(ulptr, uii);
      }
      fss_contents[fssc_idx++] = umm;
      if (is_set(not_in_family, ukk)) {
	if (sample_to_fss_idx[umm] == 0xffffffffU) {
	  sample_to_fss_idx[umm] = family_idx;
	}
	clear_bit(not_in_family, ukk);
      } else {
	set_bit(ulptr, ukk);
      }

      ullii = trio_list[trio_idx];
      do {
	uii = (uint32_t)ullii;
	ujj = sample_uidx_to_idx[uii];
	fss_contents[fssc_idx++] = ujj;
	sample_to_fss_idx[ujj] = family_idx;
	set_bit(ulptr2, uii);
	if (++trio_idx == trio_ct) {
	  goto get_sibship_info_first_pass_done;
	}
	ullii = trio_list[trio_idx];
      } while ((uint32_t)(ullii >> 32) == family_idx);
      family_idx++;
    }
  }
 get_sibship_info_first_pass_done:
  bitfield_andnot(not_in_family, ulptr2, unfiltered_sample_ctl);
  wkspace_shrink_top(fss_contents, (fssc_idx + popcount_longs(not_in_family, unfiltered_sample_ctl)) * sizeof(int32_t));
  bitfield_andnot(ulptr, ulptr2, unfiltered_sample_ctl);
  if (is_within2) {
    bitfield_andnot(tmp_within2_founder, ulptr, unfiltered_sample_ctl);
    bitfield_and(tmp_within2_founder, founder_info, unfiltered_sample_ctl);
    // now this only consists of founder parents who (i) aren't in multiple
    // families, and (ii) have a different phenotype from their partner.
    collapse_copy_bitarr(unfiltered_sample_ct, tmp_within2_founder, sample_exclude, sample_ct, lm_within2_founder);
  }
  bitfield_andnot_reversed_args(ulptr, pheno_nm, unfiltered_sample_ctl);
  if (test_type) {
    if (test_type == QFAM_WITHIN1) {
      bitfield_andnot(ulptr, founder_info, unfiltered_sample_ctl);
    }
    collapse_copy_bitarr(unfiltered_sample_ct, ulptr, sample_exclude, sample_ct, lm_eligible);
  }
  topsize = ulii;

  memcpy(ulptr, not_in_family, unfiltered_sample_ctl * sizeof(intptr_t));
  bitfield_andnot(ulptr, founder_info, unfiltered_sample_ctl);
  cur_sample_ct = popcount_longs(ulptr, unfiltered_sample_ctl);
  wkspace_left -= topsize;
  if (wkspace_alloc_ui_checked(&fs_starts, (1 + family_ct + (cur_sample_ct / 2)) * sizeof(int32_t))) {
    goto get_sibship_info_ret_NOMEM2;
  }
  wkspace_left += topsize;
  family_idx = 0;
  if (trio_ct) {
    fs_starts[0] = 0;
    ulii = 0; // first trio_idx in family
    for (trio_idx = 1; trio_idx < trio_ct; trio_idx++) {
      uii = (uint32_t)(trio_list[trio_idx] >> 32);
      if (uii > family_idx) {
	family_idx++;
	fs_starts[family_idx] = fs_starts[family_idx - 1] + 2 + trio_idx - ulii;
        ulii = trio_idx;
      }
    }
    family_idx++;
  }
  if (cur_sample_ct > 1) {
    // identify sibships of size >1
    ulii = topsize;
    merged_ids = (char*)top_alloc(&topsize, max_merged_id_len * cur_sample_ct);
    if (!merged_ids) {
      goto get_sibship_info_ret_NOMEM2;
    }
    for (sample_uidx = 0, sample_idx = 0; sample_idx < cur_sample_ct; sample_uidx++, sample_idx++) {
      next_set_ul_unsafe_ck(ulptr, &sample_uidx);
      bufptr = &(sample_ids[sample_uidx * max_sample_id_len]);
      bufptr2 = (char*)memchr(bufptr, '\t', max_sample_id_len);
      bufptr3 = memcpya(&(merged_ids[sample_idx * max_merged_id_len]), bufptr, 1 + ((uintptr_t)(bufptr2 - bufptr)));
      bufptr3 = strcpyax(bufptr3, &(paternal_ids[sample_uidx * max_paternal_id_len]), '\t');
      bufptr3 = strcpyax(bufptr3, &(maternal_ids[sample_uidx * max_maternal_id_len]), '\0');
      memcpy(bufptr3, &sample_uidx, sizeof(int32_t)); // assumes little-endian
    }
    qsort(merged_ids, cur_sample_ct, max_merged_id_len, strcmp_casted);
    bufptr = merged_ids;
    for (sample_idx = 1; sample_idx < cur_sample_ct; sample_idx++) {
      slen = strlen(bufptr) + 1;
      bufptr2 = &(merged_ids[sample_idx * max_merged_id_len]);
      if (!memcmp(bufptr, bufptr2, slen)) {
        fs_starts[family_idx] = fssc_idx;
	uii = *((uint32_t*)(&(bufptr[slen])));
	clear_bit(not_in_family, uii);
	ujj = sample_uidx_to_idx[uii];
        fss_contents[fssc_idx++] = ujj;
        sample_to_fss_idx[ujj] = family_idx;
	do {
	  uii = *((uint32_t*)(&(bufptr2[slen])));
	  clear_bit(not_in_family, uii);
	  ujj = sample_uidx_to_idx[uii];
          sample_to_fss_idx[ujj] = family_idx;
          fss_contents[fssc_idx++] = ujj;
	  if (++sample_idx == cur_sample_ct) {
	    family_idx++;
	    goto get_sibship_info_second_pass_done;
	  }
	  bufptr2 = &(merged_ids[sample_idx * max_merged_id_len]);
	} while (!memcmp(bufptr, bufptr2, slen));
	family_idx++;
      }
      bufptr = bufptr2;
    }
  }
 get_sibship_info_second_pass_done:
  *fs_ct_ptr = family_idx;
  fs_starts[family_idx] = fssc_idx;
  wkspace_shrink_top(fs_starts, (family_idx + 1) * sizeof(int32_t));
  if (test_type) {
    // for qfam, save singletons, and collapse sample_to_fss_idx to
    // sample_lm_to_fss_idx
    ulii = popcount_longs(not_in_family, unfiltered_sample_ctl);
    for (sample_uidx = 0, sample_idx = 0; sample_idx < ulii; sample_uidx++, sample_idx++) {
      next_set_ul_unsafe_ck(not_in_family, &sample_uidx);
      ujj = sample_uidx_to_idx[sample_uidx];
      fss_contents[fssc_idx++] = ujj;
      sample_to_fss_idx[ujj] = family_idx + sample_idx;
    }
    *singleton_ct_ptr = ulii;
    topsize = topsize_bak;
    wkspace_left -= topsize;
    ulii = popcount_longs(lm_eligible, sample_ctl);
    if (wkspace_alloc_ui_checked(&sample_lm_to_fss_idx, ulii * sizeof(int32_t))) {
      goto get_sibship_info_ret_NOMEM2;
    }
    wkspace_left += topsize;
    for (sample_uidx = 0, sample_idx = 0; sample_idx < ulii; sample_uidx++, sample_idx++) {
      next_set_ul_unsafe_ck(lm_eligible, &sample_uidx);
      sample_lm_to_fss_idx[sample_idx] = sample_to_fss_idx[sample_uidx];
    }
    *lm_eligible_ptr = lm_eligible;
    *lm_within2_founder_ptr = lm_within2_founder;
    *sample_lm_to_fss_idx_ptr = sample_lm_to_fss_idx;
    *lm_ct_ptr = ulii;
  } else {
    // return sample_to_fss_idx in place of sample_lm_to_fss_idx
    *sample_lm_to_fss_idx_ptr = sample_to_fss_idx;
  }
  *fs_starts_ptr = fs_starts;
  *fss_contents_ptr = fss_contents;
  // topsize = 0;

  while (0) {
  get_sibship_info_ret_NOMEM2:
    wkspace_left += topsize;
  get_sibship_info_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  }
  return retval;
}

// multithread globals
/*
static double* g_maxt_extreme_stat;
static double* g_maxt_thread_results;
static double* g_mperm_save_all;
static uintptr_t* g_pheno_c;
*/

static uintptr_t* g_loadbuf;
static uintptr_t* g_lm_eligible;
static uintptr_t* g_lm_within2_founder;
static uintptr_t* g_qfam_flip;
static uintptr_t* g_nm_fss;
static uintptr_t* g_nm_lm;
static uint32_t* g_qfam_permute;
static uint32_t* g_permute_edit;
static uint32_t* g_perm_2success_ct;
static uint32_t* g_perm_attempt_ct;
static uint32_t* g_fs_starts;
static uint32_t* g_fss_contents;
static uint32_t* g_sample_lm_to_fss_idx;
static unsigned char* g_perm_adapt_stop;
static uint32_t g_adapt_m_table[MODEL_BLOCKSIZE];
static double* g_orig_stat;
static double* g_pheno_d2;
static double* g_qfam_b;
static double* g_qfam_w;
static double* g_beta_sum;
static double* g_beta_ssq;
static uint32_t* g_beta_fail_cts;
static uintptr_t g_cur_perm_ct;
static double g_qt_sum_all;
static double g_qt_ssq_all;
static uint32_t g_test_type;
static uint32_t g_qfam_thread_ct;
static uint32_t g_fs_ct;
static uint32_t g_singleton_ct;
static uint32_t g_lm_ct;
static uint32_t g_family_ct;
static uint32_t g_block_size;
static uint32_t g_perms_done;
static uint32_t g_first_adapt_check;
static double g_adaptive_intercept;
static double g_adaptive_slope;
static double g_aperm_alpha;
static double g_adaptive_ci_zt;

// tried encoding this in a single 32-bit integer, but that appears to be
// slower.
const uint8_t dfam_allele_ct_table[] =
{0, 0, 3, 0,
 0, 0, 0, 0,
 3, 0, 2, 1,
 0, 0, 1, 0};

void dfam_sibship_calc(uint32_t cur_case_ct, uint32_t case_hom_a1_ct, uint32_t case_het_ct, uint32_t cur_ctrl_ct, uint32_t ctrl_hom_a1_ct, uint32_t ctrl_het_ct, uint32_t* total_a1_count_ptr, double* numer_ptr, double* denom_ptr, double* total_expected_ptr) {
  if (!cur_ctrl_ct) {
    return;
  }
  uint32_t hom_a1_ct = case_hom_a1_ct + ctrl_hom_a1_ct;
  uint32_t het_ct = case_het_ct + ctrl_het_ct;
  uint32_t total_ct = cur_case_ct + cur_ctrl_ct;
  uint32_t case_a1_ct = 2 * case_hom_a1_ct + case_het_ct;
  *total_a1_count_ptr += case_a1_ct;
  if (((!hom_a1_ct) && (!het_ct)) || (het_ct == total_ct) || (hom_a1_ct == total_ct)) {
    *total_expected_ptr += (double)((int32_t)case_a1_ct);
    return;
  }
  double hom_a1_ctd = (double)((int32_t)hom_a1_ct);
  double het_ctd = (double)((int32_t)het_ct);
  double case_ctd = (double)((int32_t)cur_case_ct);
  double ctrl_ctd = (double)((int32_t)cur_ctrl_ct);
  double total_ctd = (double)((int32_t)total_ct);
  double total_ct_recip = 1.0 / total_ctd;
  double case_proportion = case_ctd * total_ct_recip;
  double case_expected_hom_a1 = case_proportion * hom_a1_ctd;
  double case_expected_het = case_proportion * het_ctd;
  double case_ctrl_div_xxxm1 = case_proportion * ctrl_ctd / (total_ctd * (total_ctd - 1));
  double case_var_hom_a1 = case_ctrl_div_xxxm1 * hom_a1_ctd * (total_ctd - hom_a1_ctd);
  double case_var_het = case_ctrl_div_xxxm1 * het_ctd * (total_ctd - het_ctd);
  double case_neg_covar = case_ctrl_div_xxxm1 * hom_a1_ctd * het_ctd;
  double case_expected_a1_ct = 2 * case_expected_hom_a1 + case_expected_het;
  double case_var_a1_ct = 4 * (case_var_hom_a1 + case_neg_covar) + case_var_het;
  *numer_ptr += (double)((int32_t)case_a1_ct) - case_expected_a1_ct;
  *denom_ptr += case_var_a1_ct;
  *total_expected_ptr += case_expected_a1_ct;
}

int32_t dfam(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, double pfilter, double output_min_p, uint32_t mtest_adjust, double adjust_lambda, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude_orig, uintptr_t marker_ct_orig, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, uintptr_t* marker_reverse, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, Aperm_info* apip, uint32_t mperm_save, uintptr_t* pheno_c, uintptr_t* founder_info, uintptr_t* sex_nm, uintptr_t* sex_male, char* sample_ids, uintptr_t max_sample_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, Chrom_info* chrom_info_ptr, uint32_t hh_exists, uint32_t within_cmdflag, uint32_t perm_batch_size, Family_info* fam_ip, Set_info* sip) {
  logprint("Error: --dfam is currently under development.\n");
  return RET_CALC_NOT_YET_SUPPORTED;
  /*
  unsigned char* wkspace_mark = wkspace_base;
  FILE* outfile = NULL;
  FILE* outfile_msa = NULL;
  char* textbuf = tbuf;
  uintptr_t marker_ct_orig_autosomal = marker_ct_orig;
  uintptr_t unfiltered_marker_ctl = (unfiltered_marker_ct + (BITCT - 1)) / BITCT;
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t unfiltered_sample_ctl = (unfiltered_sample_ct + (BITCT - 1)) / BITCT;
  uintptr_t unfiltered_sample_ctl2 = (unfiltered_sample_ct + (BITCT2 - 1)) / BITCT2;
  uintptr_t unfiltered_sample_ctp1l2 = 1 + (unfiltered_sample_ct / BITCT2);
  uintptr_t final_mask = get_final_mask(unfiltered_sample_ct);
  uintptr_t* marker_exclude_orig_autosomal = marker_exclude_orig;
  uintptr_t* founder_pnm = NULL;
  double* orig_chisq = NULL;
  uint32_t unfiltered_sample_ctl2m1 = (unfiltered_sample_ct - 1) / BITCT2;
  uint32_t multigen = (fam_ip->mendel_modifier / MENDEL_MULTIGEN) & 1;
  uint32_t is_set_test = fam_ip->dfam_modifier & DFAM_SET_TEST;
  uint32_t perm_adapt_nst = (fam_ip->dfam_modifier & DFAM_PERM) && (!is_set_test);
  uint32_t perm_maxt_nst = (fam_ip->dfam_modifier & DFAM_MPERM) && (!is_set_test);
  uint32_t do_perms = fam_ip->dfam_modifier & (DFAM_PERM | DFAM_MPERM);
  uint32_t do_perms_nst = do_perms && (!is_set_test);
  uint32_t perm_count = fam_ip->dfam_modifier & DFAM_PERM_COUNT;
  uint32_t fill_orig_chisq = do_perms || mtest_adjust;
  uint32_t no_unrelateds = (fam_ip->dfam_modifier & DFAM_NO_UNRELATEDS) || (within_cmdflag && (!cluster_ct));
  uint32_t family_all_case_children_ct = 0;
  uint32_t family_mixed_ct = 0;
  uint32_t sibship_mixed_ct = 0;
  uint32_t unrelated_cluster_ct = 0;
  uint32_t pct = 0;
  uint32_t max_thread_ct = g_thread_ct;
  uint32_t perm_pass_idx = 0;
  uint32_t perms_total = 0;
  uint32_t perms_done = 0;
  int32_t retval = 0;
  uintptr_t* pheno_nm;
  uintptr_t* dfam_pheno_c;
  uintptr_t* loadbuf_raw;
  uintptr_t* loadbuf_ptr;
  uintptr_t* workbuf;
  uintptr_t* marker_exclude;
  uintptr_t* dfam_sample_exclude;
  double* maxt_extreme_stat = NULL;
  uint32_t mu_table[MODEL_BLOCKSIZE];
  char* outname_end2;
  char* wptr_start;
  char* wptr;
  uint64_t* family_list;
  uint64_t* trio_list;
  uint32_t* trio_error_lookup;
  uint32_t* fs_starts;
  uint32_t* fss_contents;
  uint32_t* sample_to_fss_idx;
  uint32_t* dfam_iteration_order;
  uint32_t* idx_to_uidx;
  uint32_t* uidx_to_idx;
  uint32_t* sample_to_cluster;
  uint32_t* cluster_ctrl_case_cts;
  uint32_t* cluster_write_idxs;
  uint32_t* cur_dfam_ptr;
  uintptr_t marker_ct;
  uintptr_t marker_uidx; // loading
  uintptr_t marker_uidx2; // writing
  uintptr_t trio_ct;
  uintptr_t max_fid_len;
  uintptr_t ulii;
  double numer;
  double denom;
  double total_expected;
  double case_proportion;
  double case_expected_a1_ct;
  double case_var_a1_ct;
  double dxx;
  uint32_t family_ct;
  uint32_t fs_ct;
  uint32_t sample_uidx;
  uint32_t sample_idx;
  uint32_t fs_idx;
  uint32_t fssc_start;
  uint32_t fssc_end;
  uint32_t fssc_idx;
  uint32_t unrelated_cluster_idx;
  uint32_t write_idx;
  uint32_t cur_ctrl_ct;
  uint32_t cur_case_ct;
  uint32_t dfam_sample_ct;
  uint32_t dfam_sample_ctl2;
  uint32_t chrom_fo_idx;
  uint32_t chrom_idx;
  uint32_t chrom_end;
  uint32_t block_size;
  uint32_t block_end;
  uint32_t marker_bidx;
  uint32_t marker_unstopped_ct;
  uint32_t loop_end;
  uint32_t marker_idx;
  uint32_t marker_idx2;
  uint32_t paternal_id;
  uint32_t maternal_id;
  uint32_t paternal_geno;
  uint32_t maternal_geno;
  uint32_t sibling_ct;
  uint32_t parental_a1_ct;
  uint32_t sib_idx;
  uint32_t cur_geno;
  uint32_t case_a1_ct;
  uint32_t quad_denom;
  uint32_t total_count;
  uint32_t twice_total_expected;
  uint32_t case_hom_a1_ct;
  uint32_t case_het_ct;
  uint32_t ctrl_hom_a1_ct;
  uint32_t ctrl_het_ct;
  uint32_t hom_a1_ct;
  uint32_t het_ct;
  uint32_t uii;
  uint32_t ujj;
  int32_t twice_numer;
  uii = count_non_autosomal_markers(chrom_info_ptr, marker_exclude_orig, 1, 1);
  if (uii) {
    LOGPRINTF("Excluding %u X/MT/haploid variant%s from DFAM test.\n", uii, (uii == 1)? "" : "s");
    if (uii == marker_ct_orig_autosomal) {
      logprint("Error: No variants remaining for DFAM analysis.\n");
      goto dfam_ret_INVALID_CMDLINE;
    }
    marker_ct_orig_autosomal -= uii;
    if (wkspace_alloc_ul_checked(&marker_exclude_orig_autosomal, unfiltered_marker_ctl * sizeof(intptr_t))) {
      goto dfam_ret_NOMEM;
    }
    memcpy(marker_exclude_orig_autosomal, marker_exclude_orig, unfiltered_marker_ctl * sizeof(intptr_t));
    for (chrom_fo_idx = 0; chrom_fo_idx < chrom_info_ptr->chrom_ct; chrom_fo_idx++) {
      chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
      if (is_set(chrom_info_ptr->haploid_mask, chrom_idx) || ((int32_t)chrom_idx == chrom_info_ptr->mt_code)) {
	uii = chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx];
	fill_bits(marker_exclude_orig_autosomal, uii, chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx + 1] - uii);
      }
    }
  } else if (is_set(chrom_info_ptr->haploid_mask, 0)) {
    logprint("Error: DFAM test does not support haploid data.\n");
    goto dfam_ret_INVALID_CMDLINE;
  }
  uii = popcount_longs_exclude(pheno_c, sample_exclude, unfiltered_sample_ct);
  if (!uii) {
    logprint("Error: DFAM test requires at least one case.\n");
    goto dfam_ret_INVALID_CMDLINE;
  }
  marker_exclude = marker_exclude_orig_autosomal;
  marker_ct = marker_ct_orig_autosomal;

  // PLINK 1.07 treats missing phenotypes as controls here
  if (wkspace_alloc_ul_checked(&pheno_nm, unfiltered_sample_ctl * sizeof(intptr_t))) {
    goto dfam_ret_NOMEM;
  }
  bitfield_exclude_to_include(sample_exclude, pheno_nm, unfiltered_sample_ct);
  if (is_set_test) {
    if (wkspace_alloc_ul_checked(&founder_pnm, unfiltered_sample_ctl * sizeof(intptr_t))) {
      goto dfam_ret_NOMEM;
    }
    memcpy(founder_pnm, pheno_nm, unfiltered_sample_ctl * sizeof(intptr_t));
    bitfield_and(founder_pnm, founder_info, unfiltered_sample_ctl);
    if (extract_set_union_unfiltered(sip, NULL, unfiltered_marker_ct, marker_exclude_orig_autosomal, &marker_exclude, &marker_ct)) {
      goto dfam_ret_NOMEM;
    }
  }

  // no --mendel-duos support for now
  retval = get_trios_and_families(unfiltered_sample_ct, sample_exclude, sample_ct, founder_info, sex_nm, sex_male, sample_ids, max_sample_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, NULL, &max_fid_len, NULL, NULL, &family_list, &family_ct, &trio_list, &trio_ct, &trio_error_lookup, 0, multigen);
  if (retval) {
    goto dfam_ret_1;
  }
#ifdef __LP64__
  if ((12 * sample_ct + 2 * family_ct) > 0xffffffffLLU) {
    logprint("Error: Too many samples and families for DFAM test.\n");
    goto dfam_ret_INVALID_CMDLINE;
  }
#endif
  if (get_sibship_info(unfiltered_sample_ct, sample_exclude, sample_ct, pheno_nm, founder_info, sample_ids, max_sample_id_len, max_fid_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, family_list, trio_list, family_ct, trio_ct, 0, NULL, NULL, &fs_starts, &fss_contents, &sample_to_fss_idx, &fs_ct, NULL, NULL)) {
    goto dfam_ret_NOMEM;
  }
  // Prepare final family, sibship, and unrelated cluster data structures.
  // * Families with at least one affected child are processed using regular
  //   TDT logic when possible; however, when both parents have homozygous
  //   calls, or they aren't both genotyped, we fall back on sibship logic.
  //   (Families with no affected children are entirely excluded from the
  //   test.)
  // * Only sibships with at least one affected child and one unaffected child
  //   are considered.  (I.e. the sibship fallback never applies to families
  //   with only affected children.)
  // * Only unrelated clusters with at least one affected and one unaffected
  //   member are considered.
  // The data structures are optimized for the permutation test, since the
  // computation is nearly I/O-bound without it.  Phenotypes are permuted
  // within each sibship/unrelated cluster, while transmitted alleles are
  // permuted in case-containing families.
  if (wkspace_alloc_ul_checked(&dfam_sample_exclude, unfiltered_sample_ctl * sizeof(intptr_t)) ||
      // shrink this later
      wkspace_alloc_ui_checked(&dfam_iteration_order, (sample_ct + (sample_ct / 2)) * sizeof(int32_t)) ||
      wkspace_alloc_ui_checked(&idx_to_uidx, sample_ct * sizeof(int32_t))) {
    goto dfam_ret_NOMEM;
  }
  fill_all_bits(dfam_sample_exclude, unfiltered_sample_ct);
  fill_idx_to_uidx(sample_exclude, unfiltered_sample_ct, sample_ct, idx_to_uidx);
  cur_dfam_ptr = dfam_iteration_order;
  for (fs_idx = 0; fs_idx < family_ct; fs_idx++) {
    // Scan for families with only case children.
    fssc_start = fs_starts[fs_idx] + 2;
    fssc_end = fs_starts[fs_idx + 1];
    cur_case_ct = 0;
    for (fssc_idx = fssc_start; fssc_idx < fssc_end; fssc_idx++) {
      cur_case_ct += is_set(pheno_c, idx_to_uidx[fss_contents[fssc_idx]]);
    }
    if (cur_case_ct == fssc_end - fssc_start) {
      family_all_case_children_ct++;
      // Could point to fss_contents, but I assume it's a better idea to
      // optimize the inner loop for data locality and linear access.
      // These family entries are temporarily stored as:
      // [0-1]: parent uidxs
      // [2]: number of children
      // [3...]: child uidxs
      // We collapse the indexes again later.
      sample_uidx = idx_to_uidx[fss_contents[fssc_start - 2]];
      clear_bit(dfam_sample_exclude, sample_uidx);
      *cur_dfam_ptr++ = sample_uidx;

      sample_uidx = idx_to_uidx[fss_contents[fssc_start - 1]];
      clear_bit(dfam_sample_exclude, sample_uidx);
      *cur_dfam_ptr++ = sample_uidx;

      *cur_dfam_ptr++ = cur_case_ct;
      for (fssc_idx = fssc_start; fssc_idx < fssc_end; fssc_idx++) {
	sample_uidx = idx_to_uidx[fss_contents[fssc_idx]];
	clear_bit(dfam_sample_exclude, sample_uidx);
	*cur_dfam_ptr++ = sample_uidx;
      }
    }
  }
  for (fs_idx = 0; fs_idx < family_ct; fs_idx++) {
    // Scan for families with at least one case and one control child.
    fssc_start = fs_starts[fs_idx] + 2;
    fssc_end = fs_starts[fs_idx + 1];
    cur_case_ct = 0;
    for (fssc_idx = fssc_start; fssc_idx < fssc_end; fssc_idx++) {
      cur_case_ct += is_set(pheno_c, idx_to_uidx[fss_contents[fssc_idx]]);
    }
    if (cur_case_ct && (cur_case_ct != fssc_end - fssc_start)) {
      family_mixed_ct++;
      sample_uidx = idx_to_uidx[fss_contents[fssc_start - 2]];
      clear_bit(dfam_sample_exclude, sample_uidx);
      *cur_dfam_ptr++ = sample_uidx;

      sample_uidx = idx_to_uidx[fss_contents[fssc_start - 1]];
      clear_bit(dfam_sample_exclude, sample_uidx);
      *cur_dfam_ptr++ = sample_uidx;

      *cur_dfam_ptr++ = fssc_end - fssc_start;
      for (fssc_idx = fssc_start; fssc_idx < fssc_end; fssc_idx++) {
	sample_uidx = idx_to_uidx[fss_contents[fssc_idx]];
	clear_bit(dfam_sample_exclude, sample_uidx);
	*cur_dfam_ptr++ = sample_uidx;
      }
    }
  }
  for (; fs_idx < fs_ct; fs_idx++) {
    // Scan for sibships with at least one case and one control.
    fssc_start = fs_starts[fs_idx];
    fssc_end = fs_starts[fs_idx + 1];
    cur_case_ct = 0;
    for (fssc_idx = fssc_start; fssc_idx < fssc_end; fssc_idx++) {
      cur_case_ct += is_set(pheno_c, idx_to_uidx[fss_contents[fssc_idx]]);
    }
    if (cur_case_ct && (cur_case_ct != fssc_end - fssc_start)) {
      sibship_mixed_ct++;
      // [0]: sibling ct
      // [1...]: member uidxs
      *cur_dfam_ptr++ = fssc_end - fssc_start;
      for (fssc_idx = fssc_start; fssc_idx < fssc_end; fssc_idx++) {
	sample_uidx = idx_to_uidx[fss_contents[fssc_idx]];
	clear_bit(dfam_sample_exclude, sample_uidx);
	*cur_dfam_ptr++ = sample_uidx;
      }
    }
  }
  if (!no_unrelateds) {
    if (wkspace_alloc_ui_checked(&sample_to_cluster, sample_ct * sizeof(int32_t))) {
      goto dfam_ret_NOMEM;
    }
    // --within on an empty file actually causes --dfam to behave differently
    // than no --within at all in PLINK 1.07.  Replicate this for now.
    if (within_cmdflag) {
      if (fill_sample_to_cluster(unfiltered_sample_ct, sample_exclude, sample_ct, cluster_ct, cluster_map, cluster_starts, sample_to_cluster, NULL)) {
	goto dfam_ret_NOMEM;
      }
    } else {
      // Start everyone in the same cluster.
      fill_uint_zero(sample_to_cluster, sample_ct);
      cluster_ct = 1;
    }
    for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
      // Remove families and sibships.
      if (sample_to_fss_idx[sample_idx] != 0xffffffffU) {
	sample_to_cluster[sample_idx] = 0xffffffffU;
      }
    }

    if (wkspace_alloc_ui_checked(&cluster_ctrl_case_cts, cluster_ct * 2 * sizeof(int32_t)) ||
        wkspace_alloc_ui_checked(&cluster_write_idxs, cluster_ct * sizeof(int32_t))) {
      goto dfam_ret_NOMEM;
    }
    fill_uint_zero(cluster_ctrl_case_cts, 2 * cluster_ct);
    for (sample_uidx = 0, sample_idx = 0; sample_idx < sample_ct; sample_uidx++, sample_idx++) {
      unrelated_cluster_idx = sample_to_cluster[sample_idx];
      if (unrelated_cluster_idx != 0xffffffffU) {
	cluster_ctrl_case_cts[2 * unrelated_cluster_idx + is_set(pheno_c, sample_uidx)] += 1;
      }
    }
    // Construct reduced clusters -> samples map.
    write_idx = 0;
    for (unrelated_cluster_idx = 0; unrelated_cluster_idx < cluster_ct; unrelated_cluster_idx++) {
      cur_ctrl_ct = cluster_ctrl_case_cts[2 * unrelated_cluster_idx];
      cur_case_ct = cluster_ctrl_case_cts[2 * unrelated_cluster_idx + 1];
      if (cur_ctrl_ct && cur_case_ct) {
	unrelated_cluster_ct++;
	cur_dfam_ptr[write_idx++] = cur_ctrl_ct + cur_case_ct;
	cluster_write_idxs[unrelated_cluster_idx] = write_idx;
	write_idx += cur_ctrl_ct + cur_case_ct;
      }
    }
    for (sample_uidx = 0, sample_idx = 0; sample_idx < sample_ct; sample_uidx++, sample_idx++) {
      next_unset_unsafe_ck(sample_exclude, &sample_uidx);
      unrelated_cluster_idx = sample_to_cluster[sample_idx];
      if (unrelated_cluster_idx != 0xffffffffU) {
        cur_ctrl_ct = cluster_ctrl_case_cts[2 * unrelated_cluster_idx];
	cur_case_ct = cluster_ctrl_case_cts[2 * unrelated_cluster_idx + 1];
	if (cur_ctrl_ct && cur_case_ct) {
	  uii = cluster_write_idxs[unrelated_cluster_idx];
	  cur_dfam_ptr[uii] = sample_uidx;
	  clear_bit(dfam_sample_exclude, sample_uidx);
	  cluster_write_idxs[unrelated_cluster_idx] = uii + 1;
	}
      }
    }
    cur_dfam_ptr = &(cur_dfam_ptr[write_idx]);
  }
  wkspace_reset((unsigned char*)idx_to_uidx);
  wkspace_shrink_top(dfam_iteration_order, (cur_dfam_ptr - dfam_iteration_order) * sizeof(int32_t));
  if (do_perms) {
    logprint("Error: --dfam permutation tests are currently under development.\n");
    retval = RET_CALC_NOT_YET_SUPPORTED;
    goto dfam_ret_1;
  }
  if (mtest_adjust || do_perms) {
    if (wkspace_alloc_d_checked(&orig_chisq, marker_ct * sizeof(double))) {
      goto dfam_ret_NOMEM;
    }
  }
  dfam_sample_ct = unfiltered_sample_ct - popcount_longs(dfam_sample_exclude, unfiltered_sample_ctl);
  dfam_sample_ctl2 = (dfam_sample_ct + (BITCT2 - 1)) / BITCT2;
  if (wkspace_alloc_ui_checked(&uidx_to_idx, unfiltered_sample_ct * sizeof(int32_t))) {
    goto dfam_ret_NOMEM;
  }
  fill_uidx_to_idx(dfam_sample_exclude, unfiltered_sample_ct, dfam_sample_ct, uidx_to_idx);
  cur_dfam_ptr = dfam_iteration_order;
  uii = family_all_case_children_ct + family_mixed_ct;
  for (fs_idx = 0; fs_idx < uii; fs_idx++) {
    *cur_dfam_ptr = uidx_to_idx[*cur_dfam_ptr];
    cur_dfam_ptr++;
    *cur_dfam_ptr = uidx_to_idx[*cur_dfam_ptr];
    cur_dfam_ptr++;
    sibling_ct = *cur_dfam_ptr++;
    for (sib_idx = 0; sib_idx < sibling_ct; sib_idx++) {
      *cur_dfam_ptr = uidx_to_idx[*cur_dfam_ptr];
      cur_dfam_ptr++;
    }
  }
  uii = sibship_mixed_ct + unrelated_cluster_ct;
  for (fs_idx = 0; fs_idx < uii; fs_idx++) {
    sibling_ct = *cur_dfam_ptr++;
    for (sib_idx = 0; sib_idx < sibling_ct; sib_idx++) {
      *cur_dfam_ptr = uidx_to_idx[*cur_dfam_ptr];
      cur_dfam_ptr++;
    }
  }
  // DEBUG
  printf("%u %u %u %u\n", family_all_case_children_ct, family_mixed_ct, sibship_mixed_ct, unrelated_cluster_ct);
  wkspace_reset((unsigned char*)uidx_to_idx);
  if (wkspace_alloc_ul_checked(&dfam_pheno_c, dfam_sample_ctl2 * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&loadbuf_raw, unfiltered_sample_ctl2 * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&workbuf, unfiltered_sample_ctp1l2 * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&g_loadbuf, MODEL_BLOCKSIZE * dfam_sample_ctl2 * sizeof(intptr_t))) {
    goto dfam_ret_NOMEM;
  }
  collapse_copy_bitarr(sample_ct, pheno_c, dfam_sample_exclude, dfam_sample_ct, dfam_pheno_c);
  g_pheno_c = dfam_pheno_c;
  loadbuf_raw[unfiltered_sample_ctl2 - 1] = 0;
  workbuf[unfiltered_sample_ctp1l2 - 1] = 0;
  for (ulii = 1; ulii <= MODEL_BLOCKSIZE; ulii++) {
    // defensive
    g_loadbuf[dfam_sample_ctl2 * ulii - 1] = 0;
  }
  // no X/haploid/MT, so no haploid filters

  if (fill_orig_chisq) {
    if (wkspace_alloc_d_checked(&g_orig_stat, marker_ct * sizeof(double))) {
      goto dfam_ret_NOMEM;
    }
  }

  ulii = 2 * max_marker_allele_len + plink_maxsnp + MAX_ID_LEN + 256;
  if (ulii > MAXLINELEN) {
    if (wkspace_alloc_c_checked(&textbuf, ulii)) {
      goto dfam_ret_NOMEM;
    }
  }

  // permutation test boilerplate mostly copied from qassoc() in plink_assoc.c,
  // since it's also restricted to autosomes
  g_mperm_save_all = NULL;
  if (perm_maxt_nst) {
    perms_total = fam_ip->dfam_mperm_val;
    if (wkspace_alloc_d_checked(&maxt_extreme_stat, perms_total * sizeof(double))) {
      goto dfam_ret_NOMEM;
    }
    g_maxt_extreme_stat = maxt_extreme_stat;
    fill_double_zero(maxt_extreme_stat, perms_total);
    if (mperm_save & MPERM_DUMP_ALL) {
      memcpy(outname_end, ".mperm.dump.all", 16);
      if (fopen_checked(&outfile_msa, outname, "w")) {
        goto dfam_ret_OPEN_FAIL;
      }
      if (putc_checked('0', outfile_msa)) {
	goto dfam_ret_WRITE_FAIL;
      }
      LOGPRINTF("Dumping all permutation chi-square values to %s .\n", outname);
    }
  } else {
    mperm_save = 0;
    if (perm_adapt_nst) {
      g_aperm_alpha = apip->alpha;
      perms_total = apip->max;
      if (wkspace_alloc_ui_checked(&g_perm_attempt_ct, marker_ct * sizeof(int32_t)) ||
          wkspace_alloc_uc_checked(&g_perm_adapt_stop, marker_ct)) {
        goto dfam_ret_NOMEM;
      }
      ujj = apip->max;
      for (uii = 0; uii < marker_ct; uii++) {
	g_perm_attempt_ct[uii] = ujj;
      }
      fill_ulong_zero((uintptr_t*)g_perm_adapt_stop, (marker_ct + sizeof(intptr_t) - 1) / sizeof(intptr_t));
      g_adaptive_ci_zt = ltqnorm(1 - apip->beta / (2.0 * ((intptr_t)marker_ct)));
      if (apip->min < apip->init_interval) {
        g_first_adapt_check = (int32_t)(apip->init_interval);
      } else {
	g_first_adapt_check = apip->min;
      }
      g_adaptive_intercept = apip->init_interval;
      g_adaptive_slope = apip->interval_slope;
    }
  }

  outname_end2 = memcpyb(outname_end, ".dfam", 6);
  if (fopen_checked(&outfile, outname, "w")) {
    goto dfam_ret_OPEN_FAIL;
  }
  LOGPRINTFWW5("Writing --dfam results to %s ... ", outname);
  fflush(stdout);
  sprintf(textbuf, " CHR %%%us   A1   A2      OBS      EXP        CHISQ            P \n", plink_maxsnp);
  fprintf(outfile, textbuf, "SNP");
  loop_end = marker_ct / 100;
  marker_unstopped_ct = marker_ct;

  if (do_perms) {
    if (fam_ip->dfam_modifier & DFAM_PERM) {
      if (perm_batch_size > apip->max) {
        perm_batch_size = apip->max;
      }
    } else {
      if (perm_batch_size > fam_ip->dfam_mperm_val) {
        perm_batch_size = fam_ip->dfam_mperm_val;
      }
    }
  }
  
  fputs("0%", stdout);
  fflush(stdout);
  // ----- begin main loop -----
 dfam_more_perms:
  if (do_perms_nst) {
    if (!perm_pass_idx) {
      // ...
    }
  }
  chrom_fo_idx = 0xffffffffU;
  marker_uidx = next_unset_unsafe(marker_exclude, 0);
  if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
    goto dfam_ret_READ_FAIL;
  }
  marker_idx = 0;
  marker_idx2 = 0;
  chrom_end = 0;
  do {
    // since X/haploid/MT is not supported, ignore chromosome boundaries in
    // this loop
    block_size = 0;
    block_end = marker_unstopped_ct - marker_idx;
    if (block_end > MODEL_BLOCKSIZE) {
      block_end = MODEL_BLOCKSIZE;
    }
    do {
      if (perm_adapt_nst && g_perm_adapt_stop[marker_idx2]) {
        do {
	  marker_uidx++;
	  next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
	  marker_idx2++;
	} while (g_perm_adapt_stop[marker_idx2]);
	if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	  goto dfam_ret_READ_FAIL;
	}
      }
      if (load_raw2(bedfile, loadbuf_raw, unfiltered_sample_ct4, unfiltered_sample_ctl2m1, final_mask)) {
	goto dfam_ret_READ_FAIL;
      }
      if (IS_SET(marker_reverse, marker_uidx)) {
	reverse_loadbuf((unsigned char*)loadbuf_raw, unfiltered_sample_ct);
      }
      erase_mendel_errors(unfiltered_sample_ct, loadbuf_raw, workbuf, trio_error_lookup, trio_ct, multigen);
      collapse_copy_2bitarr(loadbuf_raw, &(g_loadbuf[block_size * dfam_sample_ctl2]), unfiltered_sample_ct, dfam_sample_ct, dfam_sample_exclude);
      if (perm_adapt_nst) {
	g_adapt_m_table[block_size] = marker_idx2++;
      }
      mu_table[block_size++] = marker_uidx;
      if (marker_idx + block_size == marker_unstopped_ct) {
	break;
      }
      marker_uidx++;
      if (IS_SET(marker_exclude, marker_uidx)) {
	marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
	if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	  goto dfam_ret_READ_FAIL;
	}
      }
    } while (block_size < block_end);
    if (!perm_pass_idx) {
      // Calculate original chi-square values and write to disk:
      // 1. Iterate through nuclear families with only case children.  If both
      //    parents are not heterozygous, either parent has a missing call, or
      //    all children have missing calls, skip.  Otherwise,
      //      twice_numer          += 2 * [A1 allele count among kids] -
      //                              ([# of kids] * [parental A1 allele ct])
      //      quad_denom           += # of het parents
      //      total_count          += [A1 allele count among kids]
      //      twice_total_expected += [# of kids] * [parental A1 allele ct]
      // 2. Iterate through nuclear families with at least one case and at
      //    least one control child.  If all case children have missing calls,
      //    skip.  Otherwise, if both parents are not heterozygous, or either
      //    parent has a missing call, handle the children as in step 3.
      //    Otherwise,
      //      twice_numer          += 2 * [A1 allele count among case kids] -
      //                              ([# of case kids] * [parental A1 ct])
      //      quad_denom           += # of het parents
      //      total_count          += [A1 allele count among case kids]
      //      twice_total_expected += [# of case kids] * [parental A1 ct]
      // 3. Iterate through sibships.  If all case siblings, or all control
      //    siblings, have missing genotypes, skip.  Otherwise (see lines
      //    420-456 of PLINK 1.07 dfam.cpp),
      //      case_expected_hom_a1 := [case sib ct] * [sib hom A1 ct] /
      //                              [sib ct]
      //      case_expected_het    := [case sib ct] * [sib het ct] / [sib ct]
      //      case_var_hom_a1      := ([case sib ct] * [ctrl sib ct] *
      //                               [sib hom A1 ct] * [sib non-hom-A1]) /
      //                              ([sib ct] * [sib ct] * ([sib ct - 1]))
      //      case_var_het         := ([case sib ct] * [ctrl sib ct] *
      //                               [sib het ct] * [sib non-het]) /
      //                              ([sib ct] * [sib ct] * ([sib ct - 1]))
      //      case_neg_covar       := ([case sib ct] * [ctrl sib ct] *
      //      (between case hom a1     [sib hom A1 ct] * [sib het ct]) /
      //       and case het cts)      ([sib ct] * [sib ct] * ([sib ct] - 1))
      //      case_expected_a1_ct  := 2 * case_expected_hom_a1 +
      //                              case_expected_het
      //      case_var_a1_ct       := 4 * case_var_hom_a1 + case_var_het +
      //                              4 * case_neg_covar
      //      numer          += case_a1_ct - case_expected_a1_ct
      //      denom          += case_var_a1_ct
      //      total_count    += case_a1_ct
      //      total_expected += case_expected_a1_ct
      //    Shortcut when all genotypes are identical (this is common):
      //      total_count    += case_a1_ct
      //      total_expected += case_a1_ct
      //      We could entirely skip this instead, but that would lead to a
      //      different output file than 1.07.
      // 4. Iterate through clusters of unrelateds.  If all genotypes are
      //    missing or identical, skip.  Otherwise (see lines 557-571 of
      //    dfam.cpp),
      //      case_expected_a1_ct := [case ct] * [A1 ct] / [cluster size]
      //      case_var_a1_ct      := ([case ct] * [ctrl ct]
      //                              [A1 ct] * [A2 ct]) /
      //                             (([clst size]^2) * ([clst size] - 1))
      //      numer          += case_a1_ct - case_expected_a1_ct
      //      denom          += case_var_a1_ct
      //      total_count    += case_a1_ct
      //      total_expected += case_expected_a1_ct
      for (marker_bidx = 0; marker_bidx < block_size; marker_bidx++) {
	marker_uidx2 = mu_table[marker_bidx];
	// marker_idx_to_uidx[marker_idx + marker_bidx] = marker_uidx2;
	loadbuf_ptr = &(g_loadbuf[marker_bidx * dfam_sample_ctl2]);
	cur_dfam_ptr = dfam_iteration_order;
	twice_numer = 0;
	quad_denom = 0;
	total_count = 0;
	numer = 0.0;
	denom = 0.0;
	twice_total_expected = 0;
	total_expected = 0;
	for (fs_idx = 0; fs_idx < family_all_case_children_ct; fs_idx++) {
          paternal_id = *cur_dfam_ptr++;
	  maternal_id = *cur_dfam_ptr++;
	  sibling_ct = *cur_dfam_ptr++;
	  paternal_geno = (loadbuf_ptr[paternal_id / BITCT2] >> (2 * (paternal_id % BITCT2))) & 3;
	  maternal_geno = (loadbuf_ptr[maternal_id / BITCT2] >> (2 * (maternal_id % BITCT2))) & 3;
	  parental_a1_ct = dfam_allele_ct_table[paternal_geno * 4 + maternal_geno];
	  if (!parental_a1_ct) {
	    cur_dfam_ptr = &(cur_dfam_ptr[sibling_ct]);
	    continue;
	  }
	  cur_case_ct = 0;
	  case_a1_ct = 0;
          for (sib_idx = 0; sib_idx < sibling_ct; sib_idx++) {
            sample_idx = *cur_dfam_ptr++;
	    cur_geno = (loadbuf_ptr[sample_idx / BITCT2] >> (2 * (sample_idx % BITCT2))) & 3;
	    if (cur_geno == 1) {
	      continue;
	    }
            cur_case_ct++;
	    case_a1_ct += (4 - cur_geno) / 2;
	  }
	  if (cur_case_ct) {
	    twice_numer += (int32_t)(2 * case_a1_ct) - (int32_t)(cur_case_ct * parental_a1_ct);
	    quad_denom += 2 - (parental_a1_ct & 1);
	    total_count += case_a1_ct;
	    twice_total_expected += cur_case_ct * parental_a1_ct;
	  }
	}
	for (fs_idx = 0; fs_idx < family_mixed_ct; fs_idx++) {
          paternal_id = *cur_dfam_ptr++;
	  maternal_id = *cur_dfam_ptr++;
	  sibling_ct = *cur_dfam_ptr++;
	  paternal_geno = (loadbuf_ptr[paternal_id / BITCT2] >> (2 * (paternal_id % BITCT2))) & 3;
	  maternal_geno = (loadbuf_ptr[maternal_id / BITCT2] >> (2 * (maternal_id % BITCT2))) & 3;
	  parental_a1_ct = dfam_allele_ct_table[paternal_geno * 4 + maternal_geno];
	  cur_case_ct = 0;
	  cur_ctrl_ct = 0;
	  case_hom_a1_ct = 0;
	  case_het_ct = 0;
	  ctrl_hom_a1_ct = 0;
	  ctrl_het_ct = 0;
          for (sib_idx = 0; sib_idx < sibling_ct; sib_idx++) {
            sample_idx = *cur_dfam_ptr++;
	    cur_geno = (loadbuf_ptr[sample_idx / BITCT2] >> (2 * (sample_idx % BITCT2))) & 3;
	    if (cur_geno == 1) {
	      continue;
	    }
	    if (IS_SET(dfam_pheno_c, sample_idx)) {
	      cur_case_ct++;
	      if (cur_geno != 3) {
		if (cur_geno == 2) {
		  case_het_ct++;
		} else {
		  case_hom_a1_ct++;
		}
	      }
	    } else {
	      cur_ctrl_ct++;
	      if (cur_geno != 3) {
		if (cur_geno == 2) {
                  ctrl_het_ct++;
		} else {
		  ctrl_hom_a1_ct++;
		}
	      }
	    }
	  }
	  if (!cur_case_ct) {
	    continue;
	  }
          if (!parental_a1_ct) {
	    dfam_sibship_calc(cur_case_ct, case_hom_a1_ct, case_het_ct, cur_ctrl_ct, ctrl_hom_a1_ct, ctrl_het_ct, &total_count, &numer, &denom, &total_expected);
	  } else {
	    case_a1_ct = 2 * case_hom_a1_ct + case_het_ct;
	    twice_numer += (int32_t)(2 * case_a1_ct) - (int32_t)(cur_case_ct * parental_a1_ct);
	    quad_denom += 2 - (parental_a1_ct & 1);
	    total_count += case_a1_ct;
	    twice_total_expected += cur_case_ct * parental_a1_ct;
	  }
	}
	numer += 0.5 * ((double)twice_numer);
	denom += 0.25 * ((double)quad_denom);
	total_expected += 0.5 * ((double)twice_total_expected);
	for (fs_idx = 0; fs_idx < sibship_mixed_ct; fs_idx++) {
	  sibling_ct = *cur_dfam_ptr++;
	  cur_case_ct = 0;
	  cur_ctrl_ct = 0;
	  case_hom_a1_ct = 0;
	  case_het_ct = 0;
	  ctrl_hom_a1_ct = 0;
	  ctrl_het_ct = 0;
          for (sib_idx = 0; sib_idx < sibling_ct; sib_idx++) {
            sample_idx = *cur_dfam_ptr++;
	    cur_geno = (loadbuf_ptr[sample_idx / BITCT2] >> (2 * (sample_idx % BITCT2))) & 3;
	    if (cur_geno == 1) {
	      continue;
	    }
	    if (IS_SET(dfam_pheno_c, sample_idx)) {
	      cur_case_ct++;
	      if (cur_geno != 3) {
		if (cur_geno == 2) {
		  case_het_ct++;
		} else {
		  case_hom_a1_ct++;
		}
	      }
	    } else {
	      cur_ctrl_ct++;
	      if (cur_geno != 3) {
		if (cur_geno == 2) {
                  ctrl_het_ct++;
		} else {
		  ctrl_hom_a1_ct++;
		}
	      }
	    }
	  }
	  if (!cur_case_ct) {
	    continue;
	  }
	  dfam_sibship_calc(cur_case_ct, case_hom_a1_ct, case_het_ct, cur_ctrl_ct, ctrl_hom_a1_ct, ctrl_het_ct, &total_count, &numer, &denom, &total_expected);
	}
	for (unrelated_cluster_idx = 0; unrelated_cluster_idx < unrelated_cluster_ct; unrelated_cluster_idx++) {
	  sibling_ct = *cur_dfam_ptr++; // not actually siblings
	  cur_case_ct = 0;
	  cur_ctrl_ct = 0;
	  case_hom_a1_ct = 0;
	  case_het_ct = 0;
	  ctrl_hom_a1_ct = 0;
	  ctrl_het_ct = 0;
          for (sib_idx = 0; sib_idx < sibling_ct; sib_idx++) {
            sample_idx = *cur_dfam_ptr++;
	    cur_geno = (loadbuf_ptr[sample_idx / BITCT2] >> (2 * (sample_idx % BITCT2))) & 3;
	    if (cur_geno == 1) {
	      continue;
	    }
	    if (IS_SET(dfam_pheno_c, sample_idx)) {
	      cur_case_ct++;
	      if (cur_geno != 3) {
		if (cur_geno == 2) {
		  case_het_ct++;
		} else {
		  case_hom_a1_ct++;
		}
	      }
	    } else {
	      cur_ctrl_ct++;
	      if (cur_geno != 3) {
		if (cur_geno == 2) {
		  ctrl_het_ct++;
		} else {
		  ctrl_hom_a1_ct++;
		}
	      }
	    }
	  }
	  case_a1_ct = 2 * case_hom_a1_ct + case_het_ct;
	  hom_a1_ct = case_hom_a1_ct + ctrl_hom_a1_ct;
	  het_ct = case_het_ct + ctrl_het_ct;
	  uii = cur_case_ct + cur_ctrl_ct;
	  if ((uii <= 1) || ((!hom_a1_ct) && (!het_ct)) || (hom_a1_ct == uii) || (het_ct == uii)) {
	    continue;
	  }
	  total_count += case_a1_ct;
	  if ((!cur_case_ct) || (!cur_ctrl_ct)) {
	    total_expected += (double)((int32_t)case_a1_ct);
	    continue;
	  }
	  dxx = ((double)((int32_t)uii));
	  case_proportion = ((double)((int32_t)cur_case_ct)) / dxx;
	  ujj = 2 * hom_a1_ct + het_ct;
	  case_expected_a1_ct = case_proportion * ((double)((int32_t)ujj));
	  case_var_a1_ct = case_expected_a1_ct * cur_ctrl_ct * ((double)((int32_t)(2 * uii - ujj))) / (dxx * (dxx - 1));
          numer += case_a1_ct - case_expected_a1_ct;
	  denom += case_var_a1_ct;
	  total_expected += case_expected_a1_ct;
	}
	printf("%g %g %u %g\n", numer, denom, total_count, total_expected);
	if (marker_bidx == 2) {
	  exit(1);
	}
      }
    }
    marker_idx += block_size;
    if ((!perm_pass_idx) && (marker_idx >= loop_end)) {
      if (marker_idx < marker_unstopped_ct) {
	if (pct >= 10) {
	  putchar('\b');
	}
        pct = (marker_idx * 100LLU) / marker_unstopped_ct;
	printf("\b\b%u%%", pct);
	fflush(stdout);
	loop_end = (((uint64_t)pct + 1LLU) * marker_unstopped_ct) / 100;
      }
    }
  } while (marker_idx < marker_unstopped_ct);
  if (!perm_pass_idx) {
    if (pct >= 10) {
      putchar('\b');
    }
    fputs("\b\b", stdout);
    logprint("done.\n");
    if (do_perms_nst) {
      // wkspace_reset();
    }
    if (fclose_null(&outfile)) {
      goto dfam_ret_WRITE_FAIL;
    }
    if (!is_set_test) {
      if (wkspace_alloc_ui_checked(&idx_to_uidx, marker_ct * sizeof(int32_t))) {
	goto dfam_ret_NOMEM;
      }
      fill_idx_to_uidx(marker_exclude, unfiltered_marker_ct, marker_ct, idx_to_uidx);
      retval = multcomp(outname, outname_end, idx_to_uidx, marker_ct, marker_ids, max_marker_id_len, plink_maxsnp, chrom_info_ptr, orig_chisq, pfilter, output_min_p, mtest_adjust, 0, adjust_lambda, NULL, NULL);
      if (retval) {
	goto dfam_ret_1;
      }
      wkspace_reset(idx_to_uidx);
      // if (mperm_save & MPERM_DUMP_ALL) { ...
    } else {
      // retval = dfam_set_test(threads, bedfile, bed_offset, outname, outname_end, ...);
      if (retval) {
        goto dfam_ret_1;
      }
    }
  }
  if (do_perms_nst) {
    // if (mperm_save & MPERM_DUMP_ALL) { ...
    // wkspace_reset();
    if (perms_done < perms_total) {
    }
  }
  // ...
  
  while (0) {
  dfam_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  dfam_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  dfam_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  dfam_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  dfam_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  }
 dfam_ret_1:
  wkspace_reset(wkspace_mark);
  fclose_cond(outfile);
  fclose_cond(outfile_msa);
  return retval;
  */
}

void uint32_permute(uint32_t* perm_arr, uint32_t* precomputed_mods, sfmt_t* sfmtp, uint32_t ct) {
  // Sets perm_arr[0..(ct-1)] to a random permutation of 0..(ct-1).  Assumes
  // ct >= 2.
  // Will probably go into plink_common at some point.
  uint32_t write_idx;
  uint32_t lbound;
  uint32_t urand;
  perm_arr[0] = 0;
  for (write_idx = 1; write_idx < ct; write_idx++) {
    lbound = *precomputed_mods++;
    do {
      urand = sfmt_genrand_uint32(sfmtp);
    } while (urand < lbound);
    // integer modulus is slow.  some of the other permutation generator handle
    // many at once, in a manner that allows magic number divide to be employed
    // efficiently.
    urand %= write_idx + 1;
    perm_arr[write_idx] = perm_arr[urand];
    perm_arr[urand] = write_idx;
  }
}

void qfam_compute_bw(uintptr_t* loadbuf, uintptr_t sample_ct, uint32_t* fs_starts, uint32_t* fss_contents, uint32_t* sample_lm_to_fss_idx, uintptr_t* lm_eligible, uintptr_t* lm_within2_founder, uint32_t family_ct, uint32_t fs_ct, uint32_t singleton_ct, uint32_t lm_ct, uintptr_t* nm_fss, uintptr_t* nm_lm, double* pheno_d2, double qt_sum_all, double qt_ssq_all, double* qfam_b, double* qfam_w, double* qt_sum_ptr, double* qt_ssq_ptr) {
  uint32_t* fs_starts_ptr = fs_starts;
  double qt_sum = qt_sum_all;
  double qt_ssq = qt_ssq_all;
  uint32_t fss_ct = fs_ct + singleton_ct;
  uint32_t* fss_ptr;
  uint32_t* fss_end;
  uintptr_t ulii;
  uintptr_t uljj;
  double dxx;
  uint32_t cur_idx;
  uint32_t cur_start;
  uint32_t cur_end;
  uint32_t sib_ct;
  uint32_t sample_uidx;
  uint32_t sample_idx;
  uint32_t fss_idx;
  uint32_t uii;
  fill_all_bits(nm_fss, fss_ct);
  cur_start = *fs_starts_ptr++;
  for (cur_idx = 0; cur_idx < family_ct; cur_idx++) {
    cur_end = *fs_starts_ptr++;
    sample_uidx = fss_contents[cur_start];
    uii = fss_contents[cur_start + 1];
    ulii = (loadbuf[sample_uidx / BITCT2] >> (2 * (sample_uidx % BITCT2))) & 3;
    uljj = (loadbuf[uii / BITCT2] >> (2 * (uii % BITCT2))) & 3;
    if ((ulii != 1) && (uljj != 1)) {
      // both parents nonmissing
      qfam_b[cur_idx] = 0.5 * (double)(4 - ((intptr_t)((ulii + (ulii == 0)) + (uljj + (uljj == 0)))));
    } else {
      // check siblings
      sib_ct = cur_end - cur_start - 2;
      fss_ptr = &(fss_contents[cur_start + 2]);
      fss_end = &(fss_contents[cur_end]);
      uljj = 0;
      do {
        sample_uidx = *fss_ptr++;
        ulii = (loadbuf[sample_uidx / BITCT2] >> (2 * (sample_uidx % BITCT2))) & 3;
        if (ulii != 1) {
          uljj += ulii + (ulii == 0);
	} else {
	  sib_ct--;
	}
      } while (fss_ptr < fss_end);
      if (sib_ct) {
        qfam_b[cur_idx] = ((double)((intptr_t)(2 * (uintptr_t)sib_ct) - ((intptr_t)uljj))) / ((double)((int32_t)sib_ct));
      } else {
	clear_bit(nm_fss, cur_idx);
      }
    }
    cur_start = cur_end;
  }
  fss_ptr = &(fss_contents[cur_start]);
  for (; cur_idx < fs_ct; cur_idx++) {
    cur_end = *fs_starts_ptr++;
    sib_ct = cur_end - cur_start;
    fss_end = &(fss_contents[cur_end]);
    uljj = 0;
    do {
      sample_uidx = *fss_ptr++;
      ulii = (loadbuf[sample_uidx / BITCT2] >> (2 * (sample_uidx % BITCT2))) & 3;
      if (ulii != 1) {
        uljj += ulii + (ulii == 0);
      } else {
        sib_ct--;
      }
    } while (fss_ptr < fss_end);
    if (sib_ct) {
      qfam_b[cur_idx] = ((double)((intptr_t)(2 * (uintptr_t)sib_ct) - ((intptr_t)uljj))) / ((double)((int32_t)sib_ct));
    } else {
      clear_bit(nm_fss, cur_idx);
    }
    cur_start = cur_end;
  }
  for (; cur_idx < fss_ct; cur_idx++) {
    sample_uidx = *fss_ptr++;
    ulii = (loadbuf[sample_uidx / BITCT2] >> (2 * (sample_uidx % BITCT2))) & 3;
    if (ulii != 1) {
      qfam_b[cur_idx] = (double)(2 - (intptr_t)(ulii + (ulii == 0)));
    } else {
      clear_bit(nm_fss, cur_idx);
    }
  }
  fill_all_bits(nm_lm, lm_ct);
  for (sample_uidx = 0, sample_idx = 0; sample_idx < lm_ct; sample_uidx++, sample_idx++) {
    next_set_unsafe_ck(lm_eligible, &sample_uidx);
    ulii = (loadbuf[sample_uidx / BITCT2] >> (2 * (sample_uidx % BITCT2))) & 3;
    if (ulii != 1) {
      fss_idx = sample_lm_to_fss_idx[sample_idx];
      if (!is_set(nm_fss, fss_idx)) {
	goto qfam_compute_bw_skip;
      }
      if (lm_within2_founder && is_set(lm_within2_founder, sample_uidx)) {
	uii = fs_starts[fss_idx];
	if (fss_contents[uii] == sample_uidx) {
	  uii = fss_contents[uii + 1];
	} else {
	  // assert: fss_contents[uii + 1] == sample_uidx
          uii = fss_contents[uii];
	}
        if (((loadbuf[uii / BITCT2] >> (2 * (uii % BITCT2))) & 3) == 1) {
	  goto qfam_compute_bw_skip;
	}
      }
      qfam_w[sample_idx] = ((double)(2 - (intptr_t)(ulii + (ulii == 0)))) - qfam_b[fss_idx];
    } else {
    qfam_compute_bw_skip:
      dxx = pheno_d2[sample_idx];
      qt_sum -= dxx;
      qt_ssq -= dxx * dxx;
      clear_bit(nm_lm, sample_idx);
    }
  }
  // 1.07 also excludes the nonmissing parent when only one out of two parents
  // are genotyped, and no kids are (even when said parent has their own
  // genotyped parents, etc.).  We don't do this.
  *qt_sum_ptr = qt_sum;
  *qt_ssq_ptr = qt_ssq;
}

void flip_precalc(uint32_t lm_ct, double* qfam_w, double* pheno_d2, uintptr_t* nm_lm, double* geno_sum_ptr, double* geno_ssq_ptr, double* qt_g_prod_ptr) {
  // --qfam{-parents} optimizations:
  // * geno_ssq is constant, so we precompute it.
  // * The inner loop can also skip W=0 samples.  So we clear those nm_lm bits.
  // * We can also precompute geno_sum and qt_g_prod under the assumption of no
  //   flips, and then 
  double geno_sum = 0.0;
  double geno_ssq = 0.0;
  double qt_g_prod = 0.0;
  double cur_geno;
  uint32_t sample_idx;
  for (sample_idx = 0; sample_idx < lm_ct; sample_idx++) {
    if (!is_set(nm_lm, sample_idx)) {
      sample_idx = next_set(nm_lm, sample_idx, lm_ct);
      if (sample_idx == lm_ct) {
	break;
      }
    }
    cur_geno = qfam_w[sample_idx];
    if (fabs(cur_geno) < SMALL_EPSILON) {
      clear_bit(nm_lm, sample_idx);
    } else {
      geno_sum += cur_geno;
      geno_ssq += cur_geno * cur_geno;
      qt_g_prod += cur_geno * pheno_d2[sample_idx];
    }
  }
  *geno_sum_ptr = geno_sum * 0.5;
  *geno_ssq_ptr = geno_ssq;
  *qt_g_prod_ptr = qt_g_prod * 0.5;
}

static inline uint32_t qfam_regress(uint32_t test_type, uint32_t nind, uint32_t lm_ct, uint32_t* sample_lm_to_fss_idx, uintptr_t* nm_lm, double* pheno_d2, double* qfam_b, double* qfam_w, uint32_t* qfam_permute, uintptr_t* qfam_flip, double nind_recip, double qt_sum, double qt_ssq, double geno_sum, double geno_ssq, double qt_g_prod, double* beta_ptr, double* tstat_ptr) {
  // returns 0 on success
  if (nind < 3) {
    return 1;
  }
  uint32_t sample_idx = 0;
  uint32_t sample_idx2 = 0;
  uintptr_t* ulptr;
  uintptr_t* ulptr2;
  uintptr_t cur_word;
  double cur_geno;
  double qt_mean;
  double geno_mean;
  double qt_var;
  double geno_var;
  double qt_g_covar;
  double beta;
  double dxx;
  uint32_t fss_idx;
  if (test_type & (QFAM_WITHIN1 | QFAM_WITHIN2)) {
    ulptr = nm_lm;
    ulptr2 = qfam_flip;
    for (; sample_idx < lm_ct; sample_idx += BITCT) {
      cur_word = (*ulptr++) & (*ulptr2++);
      while (cur_word) {
	sample_idx2 = sample_idx + CTZLU(cur_word);
	dxx = -qfam_w[sample_idx2];
	geno_sum += dxx;
	qt_g_prod += dxx * pheno_d2[sample_idx2];
	cur_word &= cur_word - 1;
      }
    }
    geno_sum *= 2;
    qt_g_prod *= 2;
  } else {
    for (; sample_idx2 < nind; sample_idx++) {
      next_set_unsafe_ck(nm_lm, &sample_idx);
      fss_idx = qfam_permute[sample_lm_to_fss_idx[sample_idx]];
      cur_geno = qfam_b[fss_idx];
      if (test_type == QFAM_TOTAL) {
	dxx = qfam_w[sample_idx];
	if (is_set(qfam_flip, fss_idx)) {
	  cur_geno -= dxx;
	} else {
	  cur_geno += dxx;
	}
      }
      geno_sum += cur_geno;
      geno_ssq += cur_geno * cur_geno;
      qt_g_prod += cur_geno * pheno_d2[sample_idx];
      sample_idx2++;
    }
  }
  qt_mean = qt_sum * nind_recip;
  geno_mean = geno_sum * nind_recip;
  dxx = 1.0 / ((double)((int32_t)(nind - 1)));
  qt_var = (qt_ssq - qt_sum * qt_mean) * dxx;
  geno_var = (geno_ssq - geno_sum * geno_mean) * dxx;
  if (geno_var == 0.0) {
    return 1;
  }

  qt_g_covar = (qt_g_prod - qt_sum * geno_mean) * dxx;
  dxx = 1.0 / geno_var;
  beta = qt_g_covar * dxx;
  dxx = qt_var * dxx - beta * beta;
  *beta_ptr = beta;
  *tstat_ptr = beta * sqrt(((double)((int32_t)(nind - 2))) / dxx);
  return 0;
}

THREAD_RET_TYPE qfam_thread(void* arg) {
  uintptr_t tidx = (uintptr_t)arg;
  uint32_t qfam_thread_ct = g_qfam_thread_ct;
  uint32_t fs_ct = g_fs_ct;
  uint32_t lm_ct = g_lm_ct;
  uint32_t singleton_ct = g_singleton_ct;
  uint32_t fss_ct = fs_ct + singleton_ct;
  uint32_t fss_ctl = (fss_ct + (BITCT - 1)) / BITCT;
  uint32_t lm_ctl = (lm_ct + (BITCT - 1)) / BITCT;
  uint32_t test_type = g_test_type;
  uint32_t only_within = (test_type & (QFAM_WITHIN1 | QFAM_WITHIN2))? 1 : 0;
  uintptr_t* lm_eligible = g_lm_eligible;
  uintptr_t* lm_within2_founder = g_lm_within2_founder;
  uintptr_t* qfam_flip = g_qfam_flip;
  uintptr_t* nm_fss = &(g_nm_fss[tidx * CACHEALIGN32_WORD(fss_ctl)]);
  uintptr_t* nm_lm = &(g_nm_lm[tidx * CACHEALIGN32_WORD(lm_ctl)]);
  double* qfam_b = &(g_qfam_b[tidx * CACHEALIGN32_DBL(fss_ct)]);
  double* qfam_w = &(g_qfam_w[tidx * CACHEALIGN32_DBL(lm_ct)]);
  double* pheno_d2 = g_pheno_d2;
  double* beta_sum = g_beta_sum;
  double* beta_ssq = g_beta_ssq;
  uint32_t* qfam_permute = only_within? NULL : g_qfam_permute;
  uint32_t* permute_edit_buf = only_within? NULL : (&(g_permute_edit[tidx * CACHEALIGN32_INT32(fss_ct)]));
  uint32_t* perm_2success_ct = g_perm_2success_ct;
  uint32_t* perm_attempt_ct = g_perm_attempt_ct;
  uint32_t* fs_starts = g_fs_starts;
  uint32_t* fss_contents = g_fss_contents;
  uint32_t* sample_lm_to_fss_idx = g_sample_lm_to_fss_idx;
  uint32_t* perm_ptr = NULL;
  uint32_t* beta_fail_cts = g_beta_fail_cts;
  uintptr_t cur_perm_ct = g_cur_perm_ct;
  uintptr_t sample_ct = g_sample_ct;
  uintptr_t sample_ctl2 = (sample_ct + (BITCT2 - 1)) / BITCT2;
  uintptr_t flip_ctl = only_within? lm_ctl : fss_ctl;
  double adaptive_intercept = g_adaptive_intercept;
  double adaptive_slope = g_adaptive_slope;
  double adaptive_ci_zt = g_adaptive_ci_zt;
  double aperm_alpha = g_aperm_alpha;
  double qt_sum_all = g_qt_sum_all;
  double qt_ssq_all = g_qt_ssq_all;
  double geno_sum = 0.0;
  double geno_ssq = 0.0;
  double qt_g_prod = 0.0;
  uint32_t pidx_offset = g_perms_done;
  uint32_t family_ct = g_family_ct;
  uint32_t first_adapt_check = g_first_adapt_check;
  uintptr_t* __restrict__ loadbuf;
  uintptr_t pidx;
  unsigned char* perm_adapt_stop;
  double* orig_stat;
  double stat_high;
  double stat_low;
  double qt_sum;
  double qt_ssq;
  double nind_recip;
  double beta;
  double cur_beta_sum;
  double cur_beta_ssq;
  double tstat;
  double pval;
  double dxx;
  double dyy;
  double dzz;
  uint32_t marker_bidx;
  uint32_t marker_bceil;
  uint32_t marker_idx;
  uint32_t success_2start;
  uint32_t success_2incr;
  uint32_t next_adapt_check;
  uint32_t cur_beta_fail_cts;
  uint32_t cur_fss_ct;
  uint32_t nind;
  uint32_t orig_fss_idx;
  uint32_t new_fss_idx;
  uint32_t uii;
  uint32_t ujj;
  while (1) {
    if (g_block_size <= qfam_thread_ct) {
      if (g_block_size <= tidx) {
	goto qfam_thread_skip_all;
      }
      marker_bidx = tidx;
      marker_bceil = tidx + 1;
    } else {
      marker_bidx = (((uint64_t)tidx) * g_block_size) / qfam_thread_ct;
      marker_bceil = (((uint64_t)tidx + 1) * g_block_size) / qfam_thread_ct;
    }
    loadbuf = g_loadbuf;
    perm_adapt_stop = g_perm_adapt_stop;
    orig_stat = g_orig_stat;
    for (; marker_bidx < marker_bceil; marker_bidx++) {
      marker_idx = g_adapt_m_table[marker_bidx];
      if (perm_adapt_stop[marker_idx]) {
        continue;
      }
      success_2start = perm_2success_ct[marker_idx];
      success_2incr = 0;
      next_adapt_check = first_adapt_check;
      stat_high = orig_stat[marker_idx] + EPSILON;
      stat_low = orig_stat[marker_idx] - EPSILON;
      qfam_compute_bw(&(loadbuf[sample_ctl2 * marker_bidx]), sample_ct, fs_starts, fss_contents, sample_lm_to_fss_idx, lm_eligible, lm_within2_founder, family_ct, fs_ct, singleton_ct, lm_ct, nm_fss, nm_lm, pheno_d2, qt_sum_all, qt_ssq_all, qfam_b, qfam_w, &qt_sum, &qt_ssq);
      cur_fss_ct = popcount_longs(nm_fss, fss_ctl);
      nind = popcount_longs(nm_lm, lm_ctl);
      nind_recip = 1.0 / ((double)((int32_t)nind));
      if (only_within) {
	flip_precalc(lm_ct, qfam_w, pheno_d2, nm_lm, &geno_sum, &geno_ssq, &qt_g_prod);
      }
      cur_beta_sum = 0.0;
      cur_beta_ssq = 0.0;
      cur_beta_fail_cts = 0;

      for (pidx = 0; pidx < cur_perm_ct;) {
	if (!only_within) {
	  perm_ptr = &(qfam_permute[fss_ct * pidx]);
          if (cur_fss_ct != fss_ct) {
	    memcpy(permute_edit_buf, perm_ptr, fss_ct * sizeof(int32_t));
	    for (orig_fss_idx = 0, uii = 0; uii < cur_fss_ct; orig_fss_idx++, uii++) {
	      // Necessary to edit permutation so that nonmissing families are
	      // mapped to nonmissing families.  See PLINK 1.07 qfam.cpp.
	      next_set_unsafe_ck(nm_fss, &orig_fss_idx);
	      new_fss_idx = permute_edit_buf[orig_fss_idx];
	      if (is_set(nm_fss, new_fss_idx)) {
		continue;
	      }
	      // Walk through permutation cycle, swapping until a nonmissing
	      // family is found.
	      while (1) {
		ujj = permute_edit_buf[new_fss_idx];
		permute_edit_buf[new_fss_idx] = new_fss_idx;
		if (is_set(nm_fss, ujj)) {
		  break;
		}
		new_fss_idx = ujj;
	      }
	      permute_edit_buf[orig_fss_idx] = ujj;
	    }
	    perm_ptr = permute_edit_buf;
	  }
	}
	if (!qfam_regress(test_type, nind, lm_ct, sample_lm_to_fss_idx, nm_lm, pheno_d2, qfam_b, qfam_w, perm_ptr, &(qfam_flip[pidx * flip_ctl]), nind_recip, qt_sum, qt_ssq, geno_sum, geno_ssq, qt_g_prod, &beta, &tstat)) {
	  cur_beta_sum += beta;
	  cur_beta_ssq += beta * beta;
	  tstat = fabs(tstat);
	  if (tstat > stat_high) {
	    success_2incr += 2;
	  } else if (tstat > stat_low) {
	    success_2incr++;
	  }
	} else {
	  // conservative handling of permutation regression failure
	  success_2incr += 2;
	  cur_beta_fail_cts++;
	}
	if (++pidx == next_adapt_check - pidx_offset) {
	  // won't ever get here with fixed number of permutations
          uii = success_2start + success_2incr;
	  if (uii) {
	    pval = ((double)((int32_t)uii + 2)) / ((double)(2 * ((int32_t)next_adapt_check + 1)));
	    dxx = adaptive_ci_zt * sqrt(pval * (1 - pval) / ((int32_t)next_adapt_check));
            dyy = pval - dxx; // lower bound
            dzz = pval + dxx; // upper bound
            if ((dyy > aperm_alpha) || (dzz < aperm_alpha)) {
              perm_adapt_stop[marker_idx] = 1;
              perm_attempt_ct[marker_idx] = next_adapt_check;
              break;
	    }
	  }
	  next_adapt_check += (int32_t)(adaptive_intercept + ((int32_t)next_adapt_check) * adaptive_slope);
	}
      }
      perm_2success_ct[marker_idx] += success_2incr;
      if (beta_sum) {
	beta_sum[marker_idx] += cur_beta_sum;
	beta_ssq[marker_idx] += cur_beta_ssq;
	if (cur_beta_fail_cts) {
	  beta_fail_cts[marker_idx] += cur_beta_fail_cts;
	}
      }
    }
  qfam_thread_skip_all:
    if ((!tidx) || g_is_last_thread_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

int32_t qfam(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, uint32_t* marker_pos, char** marker_allele_ptrs, uintptr_t* marker_reverse, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, Aperm_info* apip, uintptr_t* pheno_nm, double* pheno_d, uintptr_t* founder_info, uintptr_t* sex_nm, uintptr_t* sex_male, char* sample_ids, uintptr_t max_sample_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, Chrom_info* chrom_info_ptr, uint32_t hh_exists, uint32_t perm_batch_size, Family_info* fam_ip) {
  // Fortunately, this can use some of qassoc()'s logic instead of punting to
  // LAPACK, since it doesn't support covariates.
  unsigned char* wkspace_mark = wkspace_base;
  FILE* outfile = NULL;
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t unfiltered_sample_ctl2 = (unfiltered_sample_ct + (BITCT2 - 1)) / BITCT2;
  uintptr_t unfiltered_sample_ctp1l2 = 1 + (unfiltered_sample_ct / BITCT2);
  uintptr_t sample_ctl2 = (sample_ct + (BITCT2 - 1)) / BITCT2;
  uintptr_t final_mask = get_final_mask(unfiltered_sample_ct);
  double qt_sum_all = 0.0;
  double qt_ssq_all = 0.0;
  double geno_sum = 0.0;
  double geno_ssq = 0.0;
  double qt_g_prod = 0.0;
  double* orig_beta = NULL;
  char* chrom_name_ptr = NULL;
  uint32_t unfiltered_sample_ctl2m1 = (unfiltered_sample_ct - 1) / BITCT2;
  uint32_t test_type = fam_ip->qfam_modifier & QFAM_TEST;
  uint32_t perm_adapt = fam_ip->qfam_modifier & QFAM_PERM;
  uint32_t multigen = (fam_ip->mendel_modifier / MENDEL_MULTIGEN) & 1;
  uint32_t only_within = (test_type & (QFAM_WITHIN1 | QFAM_WITHIN2))? 1 : 0;
  uint32_t perm_count = fam_ip->qfam_modifier & QFAM_PERM_COUNT;
  uint32_t emp_se = fam_ip->qfam_modifier & QFAM_EMP_SE;
  uint32_t perms_done = 0;
  uint32_t chrom_idx = 0;
  uint32_t qfam_thread_ct = g_thread_ct;
  uint32_t chrom_name_len = 0;
  uint32_t regress_fail_ct = 0;
  uint32_t pct = 0;
  int32_t mt_code = chrom_info_ptr->mt_code;
  int32_t retval = 0;
  const char qfam_flag_suffixes[][8] = {"within", "parents", "total", "between"};
  const char qfam_test_str[][6] = {"WITH ", " TOT ", " BET "};
  const char* qfam_test_ptr = qfam_test_str[0];
  char chrom_name_buf[3 + MAX_CHROM_TEXTNUM_LEN];
  uint32_t mu_table[MODEL_BLOCKSIZE];
  const char* flag_suffix;
  uintptr_t* lm_within2_founder;
  uintptr_t* nm_fss;
  uintptr_t* nm_lm;
  uintptr_t* lm_eligible;
  uintptr_t* loadbuf_raw;
  uintptr_t* workbuf;
  uintptr_t* dummy_flip;
  uintptr_t* loadbuf_ptr;
  uintptr_t* ulptr;
  char* bufptr;
  unsigned char* perm_adapt_stop;
  double* pheno_d2;
  double* qfam_b;
  double* qfam_w;
  double* orig_stat_ptr;
  uint64_t* family_list;
  uint64_t* trio_list;
  uint32_t* precomputed_mods;
  uint32_t* trio_error_lookup;
  uint32_t* fs_starts;
  uint32_t* fss_contents;
  uint32_t* sample_lm_to_fss_idx;
  uint32_t* dummy_perm;
  uint32_t* uiptr;
  uintptr_t trio_ct;
  uintptr_t max_fid_len;
  uintptr_t fss_ct;
  uintptr_t fss_ctl;
  uintptr_t flip_ctl;
  uintptr_t cur_perm_ct;
  uintptr_t marker_unstopped_ct;
  uintptr_t marker_uidx;
  uintptr_t marker_uidx_cur;
  uintptr_t marker_idx_base;
  uintptr_t marker_idx;
  uintptr_t marker_idx2_base;
  uintptr_t block_idx;
  uintptr_t block_size;
  uintptr_t ulii;
  uintptr_t uljj;
  double nind_recip;
  double qt_sum;
  double qt_ssq;
  double beta;
  double tstat;
  double dxx;
  uint32_t family_ct;
  uint32_t fs_ct;
  uint32_t singleton_ct;
  uint32_t lm_ct;
  uint32_t lm_ctl;
  uint32_t perms_total;
  uint32_t chrom_fo_idx;
  uint32_t chrom_end;
  uint32_t seek_flag;
  uint32_t nind;
  uint32_t loop_end;
  uint32_t is_last_block;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  if (test_type == QFAM_WITHIN1) {
    flag_suffix = qfam_flag_suffixes[0];
  } else if (test_type == QFAM_WITHIN2) {
    flag_suffix = qfam_flag_suffixes[1];
  } else if (test_type == QFAM_TOTAL) {
    flag_suffix = qfam_flag_suffixes[2];
    qfam_test_ptr = qfam_test_str[1];
  } else {
    flag_suffix = qfam_flag_suffixes[3];
    qfam_test_ptr = qfam_test_str[2];
  }

  uii = count_non_autosomal_markers(chrom_info_ptr, marker_exclude, 1, 1);
  if (uii) {
    LOGPRINTF("Excluding %u X/MT/haploid variant%s from QFAM test.\n", uii, (uii == 1)? "" : "s");
    if (uii == marker_ct) {
      logprint("Error: No variants remaining for QFAM analysis.\n");
      goto qfam_ret_INVALID_CMDLINE;
    }
    marker_ct -= uii;
  } else if (is_set(chrom_info_ptr->haploid_mask, 0)) {
    logprint("Error: QFAM test does not currently support haploid data.\n");
    goto qfam_ret_INVALID_CMDLINE;
  }
  // no --mendel-duos support for now
  retval = get_trios_and_families(unfiltered_sample_ct, sample_exclude, sample_ct, founder_info, sex_nm, sex_male, sample_ids, max_sample_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, NULL, &max_fid_len, NULL, NULL, &family_list, &family_ct, &trio_list, &trio_ct, &trio_error_lookup, 0, multigen);
  if (retval) {
    goto qfam_ret_1;
  }
  g_family_ct = family_ct;
#ifdef __LP64__
  // no need to check in 32-bit case since a nomem error would have occurred
  // earlier...
  // (okay, no need to check anyway, but best to document this overflow
  // possibility.)
  if ((sample_ct + 2 * family_ct) > 0xffffffffLLU) {
    logprint("Error: Too many samples and families for QFAM test.\n");
    goto qfam_ret_INVALID_CMDLINE;
  }
#endif
  if (get_sibship_info(unfiltered_sample_ct, sample_exclude, sample_ct, pheno_nm, founder_info, sample_ids, max_sample_id_len, max_fid_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, family_list, trio_list, family_ct, trio_ct, test_type, &lm_eligible, &lm_within2_founder, &fs_starts, &fss_contents, &sample_lm_to_fss_idx, &fs_ct, &lm_ct, &singleton_ct)) {
    goto qfam_ret_NOMEM;
  }
  fss_ct = fs_ct + singleton_ct;
  if (fss_ct < 2) {
    logprint("Error: QFAM test requires at least two families.\n");
    goto qfam_ret_INVALID_CMDLINE;
  } else if (lm_ct < 3) {
    LOGPRINTF("Error: Less than three eligible %ss for QFAM test.\n", (test_type == QFAM_WITHIN1)? "nonfounder" : "sample");
    goto qfam_ret_INVALID_CMDLINE;
  }
  g_fs_starts = fs_starts;
  g_fss_contents = fss_contents;
  g_sample_lm_to_fss_idx = sample_lm_to_fss_idx;
  g_lm_eligible = lm_eligible;
  g_lm_within2_founder = lm_within2_founder;
  g_test_type = test_type;
  g_sample_ct = sample_ct;
  g_fs_ct = fs_ct;
  g_singleton_ct = singleton_ct;
  g_lm_ct = lm_ct;
  g_qfam_thread_ct = qfam_thread_ct;
  fss_ctl = (fss_ct + BITCT - 1) / BITCT;
  lm_ctl = (lm_ct + BITCT - 1) / BITCT;
  flip_ctl = only_within? lm_ctl : fss_ctl;

  if (wkspace_alloc_uc_checked(&perm_adapt_stop, marker_ct) ||
      wkspace_alloc_ui_checked(&g_perm_2success_ct, marker_ct * sizeof(int32_t))) {
    goto qfam_ret_NOMEM;
  }
  fill_ulong_zero((uintptr_t*)perm_adapt_stop, (marker_ct + sizeof(intptr_t) - 1) / sizeof(intptr_t));
  fill_uint_zero(g_perm_2success_ct, marker_ct);
  g_perm_adapt_stop = perm_adapt_stop;

  if (perm_adapt) {
    g_aperm_alpha = apip->alpha;
    perms_total = apip->max;
    if (wkspace_alloc_ui_checked(&g_perm_attempt_ct, marker_ct * sizeof(int32_t))) {
      goto qfam_ret_NOMEM;
    }
    ujj = apip->max;
    for (uii = 0; uii < marker_ct; uii++) {
      g_perm_attempt_ct[uii] = ujj;
    }
    g_adaptive_ci_zt = ltqnorm(1 - apip->beta / (2.0 * ((intptr_t)marker_ct)));
    if (apip->min < apip->init_interval) {
      g_first_adapt_check = (int32_t)(apip->init_interval);
    } else {
      g_first_adapt_check = apip->min;
    }
    g_adaptive_intercept = apip->init_interval;
    g_adaptive_slope = apip->interval_slope;
  } else {
    g_perm_attempt_ct = NULL;
    perms_total = fam_ip->qfam_mperm_val;
    g_first_adapt_check = perms_total + 1;
    g_aperm_alpha = 0.0;
    g_adaptive_ci_zt = 0.0;
    g_adaptive_intercept = 0.0;
    g_adaptive_slope = 0.0;
  }
  outname_end = memcpya(outname_end, ".qfam.", 6);
  outname_end = strcpya(outname_end, flag_suffix);
  *outname_end = '\0';
  if (wkspace_alloc_ul_checked(&loadbuf_raw, unfiltered_sample_ctl2 * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&workbuf, unfiltered_sample_ctp1l2 * sizeof(intptr_t))) {
    goto qfam_ret_NOMEM;
  }
  loadbuf_raw[unfiltered_sample_ctl2 - 1] = 0;
  workbuf[unfiltered_sample_ctp1l2 - 1] = 0;
  if (fopen_checked(&outfile, outname, "w")) {
    goto qfam_ret_OPEN_FAIL;
  }
  if (perms_total < perm_batch_size) {
    perm_batch_size = perms_total;
  }
  ulii = CACHELINE * ((uintptr_t)qfam_thread_ct);
  if (!only_within) {
    if (wkspace_alloc_ui_checked(&g_qfam_permute, perm_batch_size * fss_ct * sizeof(intptr_t)) ||
        wkspace_alloc_ui_checked(&g_permute_edit, ((fss_ct + CACHELINE_INT32 - 1) / CACHELINE_INT32) * ulii)) {
      goto qfam_ret_NOMEM;
    }
  }
  if (emp_se) {
    if (wkspace_alloc_d_checked(&orig_beta, marker_ct * sizeof(double)) ||
        wkspace_alloc_d_checked(&g_beta_sum, marker_ct * sizeof(double)) ||
        wkspace_alloc_d_checked(&g_beta_ssq, marker_ct * sizeof(double)) ||
        wkspace_alloc_ui_checked(&g_beta_fail_cts, marker_ct * sizeof(double))) {
      goto qfam_ret_NOMEM;
    }
    fill_double_zero(g_beta_sum, marker_ct);
    fill_double_zero(g_beta_ssq, marker_ct);
    fill_uint_zero(g_beta_fail_cts, marker_ct);
  } else {
    g_beta_sum = NULL;
    g_beta_ssq = NULL;
    g_beta_fail_cts = NULL;
  }
  if (wkspace_alloc_ul_checked(&g_loadbuf, MODEL_BLOCKSIZE * sample_ctl2 * sizeof(intptr_t)) ||
      wkspace_alloc_d_checked(&g_orig_stat, marker_ct * sizeof(double)) ||
      wkspace_alloc_ul_checked(&g_qfam_flip, perm_batch_size * flip_ctl * sizeof(intptr_t)) ||
      wkspace_alloc_ui_checked(&precomputed_mods, (fss_ct - 1) * sizeof(int32_t)) ||
      wkspace_alloc_ul_checked(&nm_fss, ((fss_ct + (CACHELINE * 8 - 1)) / (CACHELINE * 8)) * ulii) ||
      wkspace_alloc_ul_checked(&nm_lm, ((lm_ct + (CACHELINE * 8 - 1)) / (CACHELINE * 8)) * ulii) ||
      wkspace_alloc_d_checked(&pheno_d2, lm_ct * sizeof(double)) ||
      wkspace_alloc_d_checked(&qfam_b, ((fss_ct + CACHELINE_DBL - 1) / CACHELINE_DBL) * ulii) ||
      wkspace_alloc_d_checked(&qfam_w, ((lm_ct + CACHELINE_DBL - 1) / CACHELINE_DBL) * ulii) ||
      wkspace_alloc_ui_checked(&dummy_perm, fss_ct * sizeof(int32_t)) ||
      wkspace_alloc_ul_checked(&dummy_flip, fss_ctl * sizeof(intptr_t))) {
    goto qfam_ret_NOMEM;
  }
  for (uii = 0, ujj = 0, ukk = 0; ujj < sample_ct; uii++, ujj++) {
    next_unset_unsafe_ck(sample_exclude, &uii);
    if (!is_set(lm_eligible, ujj)) {
      continue;
    }
    dxx = pheno_d[uii];
    pheno_d2[ukk++] = dxx;
    qt_sum_all += dxx;
    qt_ssq_all += dxx * dxx;
  }
  g_nm_fss = nm_fss;
  g_nm_lm = nm_lm;
  g_qfam_b = qfam_b;
  g_qfam_w = qfam_w;
  g_pheno_d2 = pheno_d2;
  g_qt_sum_all = qt_sum_all;
  g_qt_ssq_all = qt_ssq_all;
  precompute_mods(fss_ct, precomputed_mods);
  for (ulii = 1; ulii <= MODEL_BLOCKSIZE; ulii++) {
    g_loadbuf[ulii * sample_ctl2 - 1] = 0;
  }
  for (uii = 0; uii < fss_ct; uii++) {
    dummy_perm[uii] = uii;
  }
  fill_ulong_zero(dummy_flip, fss_ctl);

  LOGPRINTFWW("--qfam-%s: Permuting %" PRIuPTR " families/singletons, and including %u %s in linear regression.\n", flag_suffix, fss_ct, lm_ct, g_species_plural);
  LOGPRINTFWW5("Writing report to %s ... ", outname);
  fputs("0%", stdout);
  fflush(stdout);
  // deliberately rename last field to RAW_P to reduce likelihood of
  // misinterpretation.  --adjust also disabled.
  sprintf(tbuf, " CHR %%%us         BP   A1       TEST     NIND       BETA         STAT        RAW_P\n", plink_maxsnp);
  fprintf(outfile, tbuf, "SNP");
  marker_unstopped_ct = marker_ct;
  loop_end = marker_ct / 100;

  while (1) {
    g_perms_done = perms_done;
    cur_perm_ct = perm_batch_size;
    if (perm_adapt && perms_done) {
      while (g_first_adapt_check <= perms_done) {
	g_first_adapt_check += (int32_t)(apip->init_interval + ((int32_t)g_first_adapt_check) * apip->interval_slope);
      }
    }
    if (cur_perm_ct > perms_total - perms_done) {
      cur_perm_ct = perms_total - perms_done;
    }
    g_cur_perm_ct = cur_perm_ct;

    // PLINK 1.07 actually generates new sets of permutations for each marker;
    // we don't bother with that.  (Maybe we should still use multithreaded
    // permutation generation, since integer divison/modulus sucks that badly?
    // Todo: test with a dataset with ~10k samples.)
    if (only_within) {
      fill_ulong_zero(g_qfam_flip, cur_perm_ct * lm_ctl);
      ujj = fss_ctl * (BITCT / 32);
      ulptr = g_qfam_flip;
      for (ulii = 0; ulii < cur_perm_ct; ulii++) {
        uiptr = (uint32_t*)dummy_flip;
	for (uii = 0; uii < ujj; uii++) {
          uiptr[uii] = sfmt_genrand_uint32(&sfmt);
	}
        for (uii = 0; uii < lm_ct; uii++) {
          if (is_set(dummy_flip, sample_lm_to_fss_idx[uii])) {
	    set_bit(ulptr, uii);
	  }
	}
	ulptr = &(ulptr[lm_ctl]);
      }
      fill_ulong_zero(dummy_flip, fss_ctl);
    } else {
      for (ulii = 0; ulii < cur_perm_ct; ulii++) {
	uint32_permute(&(g_qfam_permute[ulii * fss_ct]), &(precomputed_mods[-1]), &sfmt, fss_ct);
      }
      uiptr = (uint32_t*)g_qfam_flip;
      uljj = cur_perm_ct * fss_ctl * (BITCT / 32);
      for (ulii = 0; ulii < uljj; ulii++) {
	*uiptr++ = sfmt_genrand_uint32(&sfmt);
      }
    }
    marker_uidx = next_unset_unsafe(marker_exclude, 0);
    marker_idx_base = 0; // regular filtered index
    marker_idx2_base = 0; // after also accounting for adaptive stop
    chrom_fo_idx = 0xffffffffU; // deliberate overflow
    chrom_end = 0;
    seek_flag = 1;
    do {
      block_size = marker_unstopped_ct - marker_idx2_base;
      if (block_size > MODEL_BLOCKSIZE) {
	block_size = MODEL_BLOCKSIZE;
      }
      for (block_idx = 0, marker_idx = marker_idx_base; block_idx < block_size; marker_uidx++, marker_idx++) {
	if (IS_SET(marker_exclude, marker_uidx)) {
	  marker_uidx = next_set_ul_unsafe(marker_exclude, marker_uidx);
	  seek_flag = 1;
	}
	if (marker_uidx >= chrom_end) {
	  while (1) {
	    do {
	      chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[(++chrom_fo_idx) + 1U];
	    } while (marker_uidx >= chrom_end);
	    chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	    if ((!IS_SET(chrom_info_ptr->haploid_mask, chrom_idx)) && (chrom_idx != (uint32_t)mt_code)) {
	      break;
	    }
	    seek_flag = 1;
	    marker_uidx = next_unset_unsafe(marker_exclude, chrom_end);
	  }
	}
	if (perm_adapt_stop[marker_idx]) {
	  seek_flag = 1;
	  continue;
	}
	if (seek_flag) {
	  if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	    goto qfam_ret_READ_FAIL;
	  }
	  seek_flag = 0;
	}
	if (load_raw2(bedfile, loadbuf_raw, unfiltered_sample_ct4, unfiltered_sample_ctl2m1, final_mask)) {
	  goto qfam_ret_READ_FAIL;
	}
	if (IS_SET(marker_reverse, marker_uidx)) {
	  reverse_loadbuf((unsigned char*)loadbuf_raw, unfiltered_sample_ct);
	}
	erase_mendel_errors(unfiltered_sample_ct, loadbuf_raw, workbuf, trio_error_lookup, trio_ct, multigen);
	loadbuf_ptr = &(g_loadbuf[block_idx * sample_ctl2]);
	collapse_copy_2bitarr(loadbuf_raw, loadbuf_ptr, unfiltered_sample_ct, sample_ct, sample_exclude);
	g_adapt_m_table[block_idx] = marker_idx;
	mu_table[block_idx++] = marker_uidx;
      }
      if (!perms_done) {
	orig_stat_ptr = &(g_orig_stat[marker_idx_base]);
	chrom_end = 0;
	for (block_idx = 0; block_idx < block_size; block_idx++) {
	  marker_uidx_cur = mu_table[block_idx];
	  if (marker_uidx_cur >= chrom_end) {
	    // note that chrom_fo_idx/chrom_end state actually needs to be
	    // restored at the end of this loop.  fortunately, that
	    // automatically happens.
	    chrom_fo_idx = get_marker_chrom_fo_idx(chrom_info_ptr, marker_uidx_cur);
	    chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	    chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx + 1];
	    chrom_name_ptr = chrom_name_buf5w4write(chrom_name_buf, chrom_info_ptr, chrom_idx, &chrom_name_len);
	  }
	  bufptr = memcpyax(tbuf, chrom_name_ptr, chrom_name_len, ' ');
	  bufptr = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx_cur * max_marker_id_len]), bufptr);
	  *bufptr++ = ' ';
	  bufptr = uint32_writew10x(bufptr, marker_pos[marker_uidx_cur], ' ');
	  if (fwrite_checked(tbuf, bufptr - tbuf, outfile)) {
	    goto qfam_ret_WRITE_FAIL;
	  }
	  fputs_w4(marker_allele_ptrs[marker_uidx_cur * 2], outfile);
	  loadbuf_ptr = &(g_loadbuf[block_idx * sample_ctl2]);
	  qfam_compute_bw(loadbuf_ptr, sample_ct, fs_starts, fss_contents, sample_lm_to_fss_idx, lm_eligible, lm_within2_founder, family_ct, fs_ct, singleton_ct, lm_ct, nm_fss, nm_lm, pheno_d2, qt_sum_all, qt_ssq_all, qfam_b, qfam_w, &qt_sum, &qt_ssq);
	  nind = popcount_longs(nm_lm, lm_ctl);
	  bufptr = memseta(tbuf, 32, 7);
	  bufptr = memcpya(bufptr, qfam_test_ptr, 5);
	  bufptr = uint32_writew8x(bufptr, nind, ' ');
	  nind_recip = 1.0 / ((double)((int32_t)nind));
	  if (only_within) {
	    flip_precalc(lm_ct, qfam_w, pheno_d2, nm_lm, &geno_sum, &geno_ssq, &qt_g_prod);
          }
	  if (!qfam_regress(test_type, nind, lm_ct, sample_lm_to_fss_idx, nm_lm, pheno_d2, qfam_b, qfam_w, dummy_perm, dummy_flip, nind_recip, qt_sum, qt_ssq, geno_sum, geno_ssq, qt_g_prod, &beta, &tstat)) {
	    bufptr = double_g_writewx4x(bufptr, beta, 10, ' ');
	    bufptr = double_g_writewx4x(bufptr, tstat, 12, ' ');
	    // do not apply --output-min-p since only the empirical p-value is
	    // supposed to be postprocessed here, not this one
	    bufptr = double_g_writewx4x(bufptr, calc_tprob(tstat, nind - 2), 12, '\n');
	    if (emp_se) {
	      orig_beta[marker_idx_base + block_idx] = beta;
	    }
	    *orig_stat_ptr++ = fabs(tstat);
	  } else {
	    bufptr = memcpya(bufptr, "        NA           NA           NA\n", 37);
	    perm_adapt_stop[marker_idx_base + block_idx] = 1;
	    *orig_stat_ptr++ = -9;
	    regress_fail_ct++;
	  }
	  if (fwrite_checked(tbuf, bufptr - tbuf, outfile)) {
	    goto qfam_ret_WRITE_FAIL;
	  }
	}
      }
      g_block_size = block_size;
      is_last_block = (marker_idx2_base + block_size == marker_unstopped_ct);
      if (spawn_threads2(threads, &qfam_thread, qfam_thread_ct, is_last_block)) {
	goto qfam_ret_THREAD_CREATE_FAIL;
      }
      ulii = 0;
      qfam_thread((void*)ulii);
      join_threads2(threads, qfam_thread_ct, is_last_block);
      marker_idx_base = marker_idx;
      marker_idx2_base += block_size;
      if ((!perms_done) && (marker_idx2_base >= loop_end) && (marker_idx2_base != marker_unstopped_ct)) {
	if (pct >= 10) {
	  putchar('\b');
	}
	pct = (marker_idx2_base * 100LLU) / marker_unstopped_ct;
	printf("\b\b%u%%", pct);
	fflush(stdout);
	loop_end = (((uint64_t)pct + 1LLU) * marker_unstopped_ct) / 100;
      }
    } while (marker_idx2_base < marker_unstopped_ct);

    if (!perms_done) {
      if (fclose_null(&outfile)) {
	goto qfam_ret_WRITE_FAIL;
      }
      if (pct >= 10) {
	putchar('\b');
      }
      fputs("\b\b", stdout);
      logprint("done.\n");
      if (regress_fail_ct) {
	LOGPRINTF("%u regression failure%s (excluding th%s from permutation test).\n", regress_fail_ct, (regress_fail_ct == 1)? "" : "s", (regress_fail_ct == 1)? "is" : "ese");
      }
    }
    perms_done += cur_perm_ct;
    if (perms_done == perms_total) {
      break;
    }
    printf("\r%u permutations complete.", perms_done);
    fflush(stdout);
    marker_unstopped_ct = marker_ct - popcount_longs((uintptr_t*)perm_adapt_stop, (marker_ct + sizeof(intptr_t) - 1) / sizeof(intptr_t));
  }
  putchar('\r');
  memcpy(outname_end, ".perm", 6);
  if (fopen_checked(&outfile, outname, "w")) {
    goto qfam_ret_OPEN_FAIL;
  }
  sprintf(tbuf, emp_se? " CHR %%%us         BETA     EMP_BETA       EMP_SE         EMP1           NP \n" : " CHR %%%us         EMP1           NP \n", plink_maxsnp);
  fprintf(outfile, tbuf, "SNP");
  chrom_fo_idx = 0xffffffffU;
  chrom_end = 0;
  for (marker_uidx = 0, marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
    next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
    if (marker_uidx >= chrom_end) {
      while (1) {
	do {
	  chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[(++chrom_fo_idx) + 1U];
	} while (marker_uidx >= chrom_end);
	chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	if ((!IS_SET(chrom_info_ptr->haploid_mask, chrom_idx)) && (chrom_idx != (uint32_t)mt_code)) {
	  break;
	}
	marker_uidx = next_unset_unsafe(marker_exclude, chrom_end);
      }
      chrom_name_ptr = chrom_name_buf5w4write(chrom_name_buf, chrom_info_ptr, chrom_idx, &chrom_name_len);
    }
    bufptr = memcpyax(tbuf, chrom_name_ptr, chrom_name_len, ' ');
    bufptr = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), bufptr);
    *bufptr++ = ' ';
    if (g_orig_stat[marker_idx] == -9) {
      if (emp_se) {
	bufptr = memcpya(bufptr, "          NA           NA           NA ", 39);
      }
      bufptr = memcpya(bufptr, "          NA           NA\n", 26);
    } else {
      uii = g_perm_2success_ct[marker_idx];
      if (g_perm_attempt_ct) {
	ujj = g_perm_attempt_ct[marker_idx];
      } else {
	ujj = perms_total;
      }
      if (emp_se) {
	bufptr = double_g_writewx4x(bufptr, orig_beta[marker_idx], 12, ' ');
	ukk = ujj - g_beta_fail_cts[marker_idx];
	if (ukk <= 1) {
          bufptr = memcpya(bufptr, "          NA ", 13);
	} else {
	  dxx = g_beta_sum[marker_idx] / ((double)((int32_t)ukk));
	  bufptr = double_g_writewx4x(bufptr, dxx, 12, ' ');
	  dxx = sqrt((g_beta_ssq[marker_idx] - g_beta_sum[marker_idx] * dxx) / ((double)((int32_t)(ukk - 1))));
          bufptr = double_g_writewx4x(bufptr, dxx, 12, ' ');
	}
      }
      if (!perm_count) {
        dxx = ((double)(uii + 2)) / ((double)(2 * (ujj + 1)));
      } else {
	dxx = ((double)uii) * 0.5;
      }
      bufptr = double_g_writewx4(bufptr, dxx, 12);
      bufptr = memseta(bufptr, 32, 3);
      bufptr = uint32_writew10x(bufptr, ujj, '\n');
    }
    if (fwrite_checked(tbuf, bufptr - tbuf, outfile)) {
      goto qfam_ret_WRITE_FAIL;
    }
  }
  if (fclose_null(&outfile)) {
    goto qfam_ret_WRITE_FAIL;
  }
  LOGPRINTFWW("Permutation test report written to %s .\n", outname);
  while (0) {
  qfam_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  qfam_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  qfam_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  qfam_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  qfam_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  qfam_ret_THREAD_CREATE_FAIL:
    retval = RET_THREAD_CREATE_FAIL;
    break;
  }
 qfam_ret_1:
  wkspace_reset(wkspace_mark);
  fclose_cond(outfile);
  return retval;
}
