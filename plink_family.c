#include "plink_family.h"

void family_init(Family_info* fam_ip) {
  fam_ip->mendel_max_trio_error = 1.0;
  fam_ip->mendel_max_var_error = 1.0;
  fam_ip->mendel_exclude_one_ratio = 0.0;
  fam_ip->mendel_modifier = 0;
  fam_ip->tdt_modifier = 0;
  fam_ip->tdt_mperm_val = 0;
}

uint32_t is_composite6(uintptr_t num) {
  // assumes num is congruent to 1 or 5 mod 6.
  uintptr_t divisor = 5;
  while (divisor * divisor <= num) {
    if (!(num % divisor)) {
      return 1;
    }
    divisor += 2;
    if (!(num % divisor)) {
      return 1;
    }
    divisor += 4;
  }
  return 0;
}

// this is probably going to go back into plink_common later...
uintptr_t geqprime(uintptr_t floor) {
  // assumes floor is odd and greater than 1.  Returns 5 if floor = 3,
  // otherwise returns the first prime >= floor.
  uintptr_t ulii = floor % 3;
  if (!ulii) {
    floor += 2;
  } else if (ulii == 1) {
    goto geqprime_1mod6;
  }
  while (is_composite6(floor)) {
    floor += 2;
  geqprime_1mod6:
    if (!is_composite6(floor)) {
      return floor;
    }
    floor += 4;
  }
  return floor;
}

// bottom 2 bits of index = child genotype
// middle 2 bits of index = paternal genotype
// top 2 bits of index = maternal genotype

// bits 0-7 = child increment (always 1)
// bits 8-15 = father increment
// bits 16-23 = mother increment
// bits 24-31 = error code

// note that \xx is octal, not decimal.
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

int32_t get_trios_and_families(uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_ct, uintptr_t* founder_info, uintptr_t* sex_nm, uintptr_t* sex_male, char* person_ids, uintptr_t max_person_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, char** fids_ptr, uintptr_t* max_fid_len_ptr, char** iids_ptr, uintptr_t* max_iid_len_ptr, uint64_t** family_list_ptr, uint32_t* family_ct_ptr, uint64_t** trio_list_ptr, uintptr_t* trio_ct_ptr, uint32_t** trio_lookup_ptr, uint32_t include_duos, uint32_t toposort) {
  // family_list has paternal indices in low 32 bits, maternal indices in high
  // 32, sorted in child ID order.
  // trio_list has child IDs in low 32 bits, family_list indices in high 32
  // bits, in PLINK 1.07 reporting order.
  // If include_duos is set, missing parental IDs are coded as
  // unfiltered_indiv_ct.  include_duos must be 0 or 1.
  // If toposort is set, trio_lookup has child IDs in [4n], paternal IDs in
  // [4n+1], maternal IDs in [4n+2], and reporting sequence index in [4n+3].
  // Otherwise, it's [3n]/[3n+1]/[3n+2].  toposort must be 0 or 1.
  //
  // fids is a list of null-terminated FIDs using trio_list indices, and iids
  // is a list of IIDs using regular unfiltered indices.  If include_duos is
  // set, iids has a trailing entry set to '0'.  (fids_ptr, iids_ptr, and the
  // corresponding lengths can be NULL.)
  //
  // PLINK 1.07 enforces <= 1 father and <= 1 mother per individual (and
  // ambiguous sex parents are not permitted), but the IDs CAN be reversed in
  // the .fam with no adverse consequences.  For backward compatibility, we
  // replicate this.  (Todo: report a warning exactly once when this happens.)
  //
  // This will probably be extended to enforce some phenotype constraint (for
  // e.g. TDT).
  unsigned char* wkspace_mark = wkspace_base;
  uint64_t* edge_list = NULL;
  uint32_t* toposort_queue = NULL;
  char* fids = NULL;
  char* iids = NULL;
  uintptr_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  uintptr_t unfiltered_indiv_ctp1l = 1 + (unfiltered_indiv_ct / BITCT);
  uintptr_t indiv_uidx = next_unset_unsafe(indiv_exclude, 0);
  uintptr_t htable_size = geqprime(2 * unfiltered_indiv_ct + 1);
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
  uint32_t* person_id_map;
  uint32_t* trio_lookup;
  uint32_t* uiptr;
  char* sorted_person_ids;
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
  uintptr_t indiv_idx;
  uintptr_t slen;
  uintptr_t uidx1;
  uintptr_t uidx2;
  uintptr_t next_probe_incr;
  uintptr_t probe_incr_incr;
  uint32_t first_sex;
  uint32_t uii;
  int32_t sorted_idx;
  founder_info2 = (uintptr_t*)top_alloc(&topsize, unfiltered_indiv_ctp1l * sizeof(intptr_t));
  if (!founder_info2) {
    goto get_trios_and_families_ret_NOMEM;
  }
  memcpy(founder_info2, founder_info, unfiltered_indiv_ctl * sizeof(intptr_t));
  if (unfiltered_indiv_ct & (BITCT - 1)) {
    SET_BIT(founder_info2, unfiltered_indiv_ct);
  } else {
    founder_info2[unfiltered_indiv_ctl] = 1;
  }
  topsize_bak = topsize;
  trio_list_tmp = (uint64_t*)top_alloc(&topsize, indiv_ct * sizeof(int64_t));
  if (!trio_list_tmp) {
    goto get_trios_and_families_ret_NOMEM;
  }
  topsize_bak2 = topsize;
  sorted_person_ids = (char*)top_alloc(&topsize, indiv_ct * max_person_id_len);
  if (!sorted_person_ids) {
    goto get_trios_and_families_ret_NOMEM;
  }
  person_id_map = (uint32_t*)top_alloc(&topsize, indiv_ct * sizeof(int32_t));
  if (!person_id_map) {
    goto get_trios_and_families_ret_NOMEM;
  }
  idbuf = (char*)top_alloc(&topsize, max_person_id_len);
  if (!idbuf) {
    goto get_trios_and_families_ret_NOMEM;
  }
  wkspace_left -= topsize;
  if (sort_item_ids_noalloc(sorted_person_ids, person_id_map, unfiltered_indiv_ct, indiv_exclude, indiv_ct, person_ids, max_person_id_len, 0, 0, strcmp_deref)) {
    wkspace_left += topsize;
    goto get_trios_and_families_ret_1;
  }
  // over-allocate here, we shrink family_list later when we know how many
  // families there are
  if (wkspace_alloc_ull_checked(&family_list, indiv_ct * sizeof(int64_t)) ||
      wkspace_alloc_ull_checked(&family_htable, htable_size * sizeof(int64_t)) ||
      wkspace_alloc_ui_checked(&family_idxs, htable_size * sizeof(int32_t))) {
    goto get_trios_and_families_ret_NOMEM2;
  }
  fill_ull_zero(family_htable, htable_size);
  wkspace_left += topsize;
  // 1. populate family_list (while using family_htable to track duplicates),
  //    determine max_iid_len, count qualifying trios
  trio_write = trio_list_tmp;
  for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_uidx++, indiv_idx++) {
    next_unset_ul_unsafe_ck(indiv_exclude, &indiv_uidx);
    idptr = &(person_ids[indiv_uidx * max_person_id_len]);
    iidptr = (char*)memchr(idptr, '\t', max_person_id_len);
    fidlen = (uintptr_t)((++iidptr) - idptr);
    if (fidlen > max_fid_len) {
      max_fid_len = fidlen;
    }
    slen = strlen(iidptr);
    if (slen >= max_iid_len) {
      max_iid_len = slen + 1;
    }
    if (IS_SET(founder_info, indiv_uidx)) {
      continue;
    }
    memcpy(idbuf, idptr, fidlen);
    pidptr1 = &(paternal_ids[indiv_uidx * max_paternal_id_len]);
    pidptr2 = &(maternal_ids[indiv_uidx * max_maternal_id_len]);
    slen = strlen(pidptr1);
    if (fidlen + slen < max_person_id_len) {
      memcpy(&(idbuf[fidlen]), pidptr1, slen + 1);
      sorted_idx = bsearch_str(idbuf, fidlen + slen, sorted_person_ids, max_person_id_len, indiv_ct);
    } else {
      sorted_idx = -1;
    }
    first_sex = 0;
    if (sorted_idx == -1) {
      if (!include_duos) {
	SET_BIT(founder_info2, indiv_uidx);
	continue;
      }
      uidx1 = unfiltered_indiv_ct;
    } else {
      uidx1 = person_id_map[(uint32_t)sorted_idx];
      if (uidx1 == indiv_uidx) {
        idbuf[fidlen - 1] = ' ';
	LOGPRINTF("Error: %s is his/her own parent.\n", idbuf);
	goto get_trios_and_families_ret_INVALID_FORMAT;
      } else if (!IS_SET(sex_nm, uidx1)) {
        idbuf[fidlen - 1] = ' ';
	LOGPRINTF("Error: Parent %s has unspecified sex.\n", idbuf);
	goto get_trios_and_families_ret_INVALID_FORMAT;
      }
      first_sex = 2 - IS_SET(sex_male, uidx1);
    }
    slen = strlen(pidptr2);
    if (fidlen + slen < max_person_id_len) {
      memcpy(&(idbuf[fidlen]), pidptr2, slen);
      sorted_idx = bsearch_str(idbuf, fidlen + slen, sorted_person_ids, max_person_id_len, indiv_ct);
    } else {
      sorted_idx = -1;
    }
    if (sorted_idx == -1) {
      if ((!include_duos) || (uidx1 == unfiltered_indiv_ct)) {
	SET_BIT(founder_info2, indiv_uidx);
	continue;
      }
      if (first_sex == 1) {
	family_code = (((uint64_t)unfiltered_indiv_ct) << 32) | ((uint64_t)uidx1);
      } else {
	family_code = (((uint64_t)uidx1) << 32) | ((uint64_t)unfiltered_indiv_ct);
      }
      next_probe_incr = unfiltered_indiv_ct + uidx1;
    } else {
      uidx2 = person_id_map[(uint32_t)sorted_idx];
      if (uidx2 == indiv_uidx) {
        idbuf[fidlen - 1] = ' ';
	LOGPRINTF("Error: %s is their own parent.\n", idbuf);
	goto get_trios_and_families_ret_INVALID_FORMAT;
      } else if (!IS_SET(sex_nm, uidx2)) {
        idbuf[fidlen - 1] = ' ';
	LOGPRINTF("Error: Parent %s has unspecified sex.\n", idbuf);
	goto get_trios_and_families_ret_INVALID_FORMAT;
      }
      uii = IS_SET(sex_male, uidx2);
      if (2 - uii == first_sex) {
	idptr[fidlen - 1] = ' ';
	LOGPRINTF("Error: %s has two %sies.\n", idptr, (first_sex == 1)? "dadd" : "momm");
	goto get_trios_and_families_ret_INVALID_FORMAT;
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
    *trio_write++ = (((uint64_t)family_idx) << 32) | ((uint64_t)indiv_uidx);
  }
  trio_ct = (uintptr_t)(trio_write - trio_list_tmp);
  wkspace_reset(wkspace_mark);
  family_list = (uint64_t*)wkspace_alloc(family_ct * sizeof(int64_t));
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
	wkspace_alloc_c_checked(&iids, (unfiltered_indiv_ct + include_duos) * max_iid_len)) {
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
    wkspace_reset((unsigned char*)edge_list);
    edge_list = (uint64_t*)wkspace_alloc(edge_ct * sizeof(int64_t));
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
  if (fids_ptr) {
    *fids_ptr = fids;
    *iids_ptr = iids;
    *max_fid_len_ptr = max_fid_len;
    *max_iid_len_ptr = max_iid_len;
    indiv_uidx = next_unset_unsafe(indiv_exclude, 0);
    for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_uidx++, indiv_idx++) {
      next_unset_ul_unsafe_ck(indiv_exclude, &indiv_uidx);
      idptr = &(person_ids[indiv_uidx * max_person_id_len]);
      iidptr = (char*)memchr(idptr, '\t', max_fid_len);
      strcpy(&(iids[indiv_uidx * max_iid_len]), &(iidptr[1]));
    }
    if (include_duos) {
      memcpy(&(iids[unfiltered_indiv_ct * max_iid_len]), "0", 2);
    }
  }
  uiptr = trio_lookup;
  for (trio_idx = 0; trio_idx < trio_ct; trio_idx++) {
    trio_code = trio_write[trio_idx];
    indiv_uidx = (uint32_t)trio_code;
    if (fids_ptr) {
      idptr = &(person_ids[indiv_uidx * max_person_id_len]);
      iidptr = (char*)memchr(idptr, '\t', max_fid_len);
      fidlen = (uintptr_t)(iidptr - idptr);
      memcpyx(&(fids[trio_idx * max_fid_len]), idptr, fidlen, '\0');
    }
    family_code = family_list[(uintptr_t)(trio_code >> 32)];
    if (!toposort) {
      *uiptr++ = indiv_uidx;
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
      indiv_uidx = (uint32_t)trio_code;
      family_code = family_list[(uintptr_t)(trio_code >> 32)];
      *uiptr++ = indiv_uidx;
      *uiptr++ = (uint32_t)family_code;
      *uiptr++ = (uint32_t)(family_code >> 32);
      *uiptr++ = trio_idx;
      if (remaining_edge_ct) {
        SET_BIT(founder_info2, indiv_uidx);
        ullii = ((uint64_t)indiv_uidx) << 32;
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
    wkspace_reset((unsigned char*)edge_list);
  }
  while (0) {
  get_trios_and_families_ret_NOMEM2:
    wkspace_left += topsize;
  get_trios_and_families_ret_NOMEM:
    retval = RET_NOMEM;
    break;
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

int32_t mendel_error_scan(Family_info* fam_ip, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uint32_t plink_maxfid, uint32_t plink_maxiid, uint32_t plink_maxsnp, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_exclude_ct_ptr, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t* indiv_exclude_ct_ptr, uintptr_t* founder_info, uintptr_t* sex_nm, uintptr_t* sex_male, char* person_ids, uintptr_t max_person_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uint32_t hh_exists, uint32_t zero_extra_chroms, Chrom_info* chrom_info_ptr, uint32_t calc_mendel) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* outfile = NULL;
  FILE* outfile_l = NULL;
  uintptr_t* indiv_male_include2 = NULL;
  uintptr_t* error_locs = NULL;
  char* varptr = NULL;
  char* chrom_name_ptr = NULL;
  unsigned char* cur_errors = NULL;
  uint64_t* family_error_cts = NULL;
  uint32_t* child_cts = NULL;
  uintptr_t marker_ct = unfiltered_marker_ct - *marker_exclude_ct_ptr;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t unfiltered_indiv_ctp1l2 = 1 + (unfiltered_indiv_ct / BITCT2);
  uintptr_t indiv_ct = unfiltered_indiv_ct - *indiv_exclude_ct_ptr;
  uintptr_t marker_uidx = 0xffffffffLU;
  uint64_t tot_error_ct = 0;
  uint32_t do_filter = fam_ip->mendel_modifier & MENDEL_FILTER;
  uint32_t include_duos = (fam_ip->mendel_modifier / MENDEL_DUOS) & 1;
  uint32_t multigen = (fam_ip->mendel_modifier / MENDEL_MULTIGEN) & 1;
  uint32_t var_first = fam_ip->mendel_modifier & MENDEL_FILTER_VAR_FIRST;
  uint32_t varlen = 0;
  uint32_t chrom_name_len = 0;
  uint32_t new_marker_exclude_ct = 0;
  uint32_t error_ct_fill = 0;
  int32_t retval = 0;
  char chrom_name_buf[4];
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
  retval = get_trios_and_families(unfiltered_indiv_ct, indiv_exclude, indiv_ct, founder_info, sex_nm, sex_male, person_ids, max_person_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, &fids, &max_fid_len, &iids, &max_iid_len, &family_list, &family_ct, &trio_list, &trio_ct, &trio_lookup, include_duos, multigen);
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
  if (wkspace_alloc_ul_checked(&loadbuf, unfiltered_indiv_ctp1l2 * sizeof(intptr_t)) ||
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
  loadbuf[unfiltered_indiv_ctp1l2 - 1] = 0;
  fill_uint_zero(error_cts, trio_ct * 3);
  fill_uint_zero(error_cts_tmp, trio_ct4 * 4);
  hh_exists &= XMHH_EXISTS;
  if (alloc_raw_haploid_filters(unfiltered_indiv_ct, hh_exists, 0, indiv_exclude, sex_male, NULL, &indiv_male_include2)) {
    goto mendel_error_scan_ret_NOMEM;
  }
  alens[0] = 0;
  alens[1] = 0;
  if (calc_mendel) {
    ulii = max_marker_allele_len * 6 + 10;
    if (ulii < 21) {
      ulii = 21;
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
    is_x = (((uint32_t)chrom_info_ptr->x_code) == chrom_idx)? 1 : 0;
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
      chrom_name_ptr = chrom_name_buf;
      chrom_name_len = 4;
      if (chrom_idx <= chrom_info_ptr->max_code) {
	chrom_num_write4(chrom_name_buf, chrom_idx);
      } else if (zero_extra_chroms) {
	memcpy(chrom_name_buf, "   0", 4);
      } else {
	ujj = strlen(chrom_info_ptr->nonstd_names[chrom_idx]);
	if (ujj < 4) {
	  fw_strcpyn(4, ujj, chrom_info_ptr->nonstd_names[chrom_idx], chrom_name_buf);
	} else {
	  chrom_name_ptr = chrom_info_ptr->nonstd_names[chrom_idx];
	  chrom_name_len = ujj;
	}
      }
    }
    if (uii != marker_uidx) {
      marker_uidx = uii;
      goto mendel_error_scan_seek;
    }
    while (1) {
      if (fread(loadbuf, 1, unfiltered_indiv_ct4, bedfile) < unfiltered_indiv_ct4) {
	goto mendel_error_scan_ret_READ_FAIL;
      }
      if (IS_SET(marker_reverse, marker_uidx)) {
        reverse_loadbuf((unsigned char*)loadbuf, unfiltered_indiv_ct);
      }
      if (is_x && hh_exists) {
	hh_reset((unsigned char*)loadbuf, indiv_male_include2, unfiltered_indiv_ct);
      }
      // missing parents are treated as having uidx equal to
      // unfiltered_indiv_ct, and we set the corresponding genotype to always
      // be missing.  This lets us avoid special-casing duos.
      SET_BIT_DBL(loadbuf, unfiltered_indiv_ct);
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
	if (marker_uidx == chrom_end) {
	  break;
	}
      mendel_error_scan_seek:
        if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	  goto mendel_error_scan_ret_READ_FAIL;
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
      if (ujj < unfiltered_indiv_ct) {
	// bleah, fids[] isn't in right order for this lookup
	cptr = &(person_ids[ujj * max_person_id_len]);
	wptr = fw_strcpyn(plink_maxfid, (uintptr_t)(((char*)memchr(cptr, '\t', max_person_id_len)) - cptr), cptr, tbuf);
      } else {
	wptr = memseta(tbuf, 32, plink_maxfid - 1);
	*wptr++ = '0';
      }
      *wptr++ = ' ';
      wptr = fw_strcpy(plink_maxiid, &(iids[ujj * max_iid_len]), wptr);
      *wptr++ = ' ';
      wptr = fw_strcpy(plink_maxiid, &(iids[((uintptr_t)(family_code >> 32)) * max_iid_len]), wptr);
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
	if (ujj != unfiltered_indiv_ct) {
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
	if (ukk != unfiltered_indiv_ct) {
	  if (ujj != unfiltered_indiv_ct) {
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
    LOGPRINTF("Report written to %s.{mendel,imendel,fmendel,lmendel}.\n", outname);
  }
  if (do_filter) {
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
	    set_bit(indiv_exclude, (uint32_t)trio_code);
	    if (ujj < unfiltered_indiv_ct) {
	      set_bit(indiv_exclude, ujj);
	    }
	    if (ukk < unfiltered_indiv_ct) {
	      set_bit(indiv_exclude, ukk);
	    }
	  } else if ((exclude_one_ratio == INFINITY) || (ujj == unfiltered_indiv_ct) || (ukk == unfiltered_indiv_ct)) {
            set_bit(indiv_exclude, (uint32_t)trio_code);
	  } else {
	    dxx = (double)((int32_t)trio_list[trio_idx * 3 + 1]);
	    dyy = (double)((int32_t)trio_list[trio_idx * 3 + 2]);
	    if (dxx > exclude_one_ratio * dyy) {
	      set_bit(indiv_exclude, ujj);
	    } else if (dyy > exclude_one_ratio * dxx) {
	      set_bit(indiv_exclude, ukk);
	    } else {
	      set_bit(indiv_exclude, (uint32_t)trio_code);
	    }
	  }
	}
      }
    }
    ulii = popcount_longs(indiv_exclude, (unfiltered_indiv_ct + (BITCT - 1)) / BITCT);
    if (unfiltered_indiv_ct == ulii) {
      LOGPRINTF("Error: All %s excluded by --me.\n", g_species_plural);
      goto mendel_error_scan_ret_ALL_SAMPLES_EXCLUDED;
    }
    *indiv_exclude_ct_ptr = ulii;
    ulii -= unfiltered_indiv_ct - indiv_ct;
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

int32_t populate_pedigree_rel_info(Pedigree_rel_info* pri_ptr, uintptr_t unfiltered_indiv_ct, char* person_ids, uintptr_t max_person_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* founder_info) {
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
  uintptr_t indiv_uidx;
  uint64_t ullii;
  char* family_ids;
  char* cur_person_id;
  char* last_family_id = NULL;
  char* cur_family_id;
  char* id_ptr;
  uint32_t* family_sizes;
  uint32_t* uiptr;
  uint32_t* uiptr2 = NULL;
  uint32_t fidx;
  int32_t family_size;
  uint32_t* remaining_indiv_idxs;
  int32_t* remaining_indiv_parent_idxs; // -1 = no parent (or nonshared)
  uint32_t remaining_indiv_ct;
  uint32_t indiv_idx_write;
  uintptr_t max_family_id_len = 0;
  char* indiv_ids;
  uint32_t* indiv_id_lookup;
  uintptr_t max_indiv_id_len = 0;
  uintptr_t max_pm_id_len;
  uint32_t family_id_ct;
  uint32_t* fis_ptr;
  char* stray_parent_ids;
  intptr_t stray_parent_ct;
  uintptr_t* processed_indivs;
  uint32_t founder_ct;
  int32_t max_family_nf = 0;
  uintptr_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  uintptr_t unfiltered_indiv_ctlm = unfiltered_indiv_ctl * BITCT;
  uint32_t* complete_indiv_idxs;
  uintptr_t complete_indiv_idx_ct;
  double* rs_ptr;
  double* rel_writer;
  double dxx;
  double* tmp_rel_space = NULL;
  double* tmp_rel_writer = NULL;

  for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ct; indiv_uidx++) {
    ujj = strlen_se(&(person_ids[indiv_uidx * max_person_id_len])) + 1;
    if (ujj > max_family_id_len) {
      max_family_id_len = ujj;
    }
    ujj = strlen(&(person_ids[indiv_uidx * max_person_id_len + ujj]));
    if (ujj >= max_indiv_id_len) {
      max_indiv_id_len = ujj + 1;
    }
  }
  if (max_paternal_id_len > max_maternal_id_len) {
    max_pm_id_len = max_paternal_id_len;
  } else {
    max_pm_id_len = max_maternal_id_len;
  }
  if (wkspace_alloc_ui_checked(&(pri_ptr->family_info_space), unfiltered_indiv_ct * sizeof(int32_t)) ||
      wkspace_alloc_ui_checked(&(pri_ptr->family_rel_nf_idxs), unfiltered_indiv_ct * sizeof(int32_t)) ||
      wkspace_alloc_ui_checked(&(pri_ptr->family_idxs), unfiltered_indiv_ct * sizeof(int32_t)) ||
      wkspace_alloc_c_checked(&family_ids, unfiltered_indiv_ct * max_family_id_len) ||
      wkspace_alloc_ui_checked(&family_sizes, unfiltered_indiv_ct * sizeof(int32_t))) {
    return RET_NOMEM;
  }

  // copy all the items over in order, then qsort, then eliminate duplicates
  // and count family sizes.
  cur_family_id = family_ids;
  cur_person_id = person_ids;
  uiptr = family_sizes;
  *uiptr = 1;
  jj = strlen_se(cur_person_id);
  memcpyx(cur_family_id, cur_person_id, jj, 0);
  for (indiv_uidx = 1; indiv_uidx < unfiltered_indiv_ct; indiv_uidx++) {
    cur_person_id = &(cur_person_id[max_person_id_len]);
    mm = strlen_se(cur_person_id);
    if ((jj != mm) || memcmp(cur_family_id, cur_person_id, mm)) {
      cur_family_id = &(cur_family_id[max_family_id_len]);
      memcpyx(cur_family_id, cur_person_id, mm, 0);
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

  if (family_id_ct < unfiltered_indiv_ct) {
    uiptr = family_sizes;
    wkspace_reset(family_ids);
    family_ids = (char*)wkspace_alloc(family_id_ct * max_family_id_len);
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
  cur_person_id = person_ids;
  for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ct; indiv_uidx++) {
    kk = bsearch_str(cur_person_id, strlen_se(cur_person_id), family_ids, max_family_id_len, family_id_ct);
    pri_ptr->family_idxs[indiv_uidx] = kk;
    if (IS_SET(founder_info, indiv_uidx)) {
      pri_ptr->family_founder_cts[(uint32_t)kk] += 1;
      pri_ptr->family_rel_nf_idxs[indiv_uidx] = uiptr[(uint32_t)kk];
      uiptr[kk] += 1;
    }
    cur_person_id = &(cur_person_id[max_person_id_len]);
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
  pri_ptr->family_info_offsets[family_id_ct] = unfiltered_indiv_ct;
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
  for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ct; indiv_uidx++) {
    fidx = pri_ptr->family_idxs[indiv_uidx];
    pri_ptr->family_info_space[uiptr[fidx]] = indiv_uidx;
    uiptr[fidx] += 1;
  }
  wkspace_reset(wkspace_mark);

  if (wkspace_alloc_ul_checked(&processed_indivs, (unfiltered_indiv_ctl + (max_family_nf + (BITCT2 - 1)) / BITCT2) * sizeof(intptr_t))) {
    return RET_NOMEM;
  }
  fill_ulong_one(&(processed_indivs[unfiltered_indiv_ctl]), (max_family_nf + (BITCT2 - 1)) / BITCT2);

  wkspace_mark2 = wkspace_base;
  for (fidx = 0; fidx < family_id_ct; fidx++) {
    family_size = family_sizes[fidx];
    founder_ct = pri_ptr->family_founder_cts[fidx];
    remaining_indiv_ct = family_size - founder_ct;
    stray_parent_ct = 0;
    if (remaining_indiv_ct) {
      memcpy(processed_indivs, founder_info, unfiltered_indiv_ctl * sizeof(intptr_t));
      if (wkspace_alloc_ui_checked(&complete_indiv_idxs, family_size * sizeof(int32_t)) ||
          wkspace_alloc_ui_checked(&remaining_indiv_idxs, remaining_indiv_ct * sizeof(int32_t)) ||
          wkspace_alloc_c_checked(&indiv_ids, family_size * max_indiv_id_len) ||
          wkspace_alloc_ui_checked(&indiv_id_lookup, family_size * sizeof(int32_t)) ||
          wkspace_alloc_i_checked(&remaining_indiv_parent_idxs, remaining_indiv_ct * 2 * sizeof(int32_t)) ||
          wkspace_alloc_c_checked(&stray_parent_ids, remaining_indiv_ct * 2 * max_pm_id_len)) {
	return RET_NOMEM;
      }
      ii = pri_ptr->family_info_offsets[fidx];
      fis_ptr = &(pri_ptr->family_info_space[ii]);
      rs_ptr = &(pri_ptr->rel_space[pri_ptr->family_rel_space_offsets[fidx]]);
      rel_writer = rs_ptr;
      cur_person_id = indiv_ids;
      for (ii = 0; ii < family_size; ii++) {
	kk = fis_ptr[(uint32_t)ii];
	jj = strlen_se(&(person_ids[kk * max_person_id_len])) + 1;
	strcpy(cur_person_id, &(person_ids[kk * max_person_id_len + jj]));
	cur_person_id = &(cur_person_id[max_indiv_id_len]);
	indiv_id_lookup[(uint32_t)ii] = ii;
      }

      if (qsort_ext(indiv_ids, family_size, max_indiv_id_len, strcmp_deref, (char*)indiv_id_lookup, sizeof(int32_t))) {
	return RET_NOMEM;
      }
      // Compile list of non-founder family member indices, and identify
      // parents who are referred to by individual ID but are NOT in the
      // dataset.
      ii = 0; // family_info_space index
      complete_indiv_idx_ct = 0;
      cur_person_id = stray_parent_ids;
      for (uii = 0; uii < remaining_indiv_ct; uii++) {
	while (IS_SET(founder_info, fis_ptr[ii])) {
	  complete_indiv_idxs[complete_indiv_idx_ct++] = fis_ptr[ii];
	  ii++;
	}
	kk = fis_ptr[ii++];

	// does not track sex for now
	id_ptr = &(paternal_ids[((uint32_t)kk) * max_paternal_id_len]);
	if (memcmp("0", id_ptr, 2)) {
	  ujj = strlen(id_ptr);
	  mm = bsearch_str(id_ptr, ujj, indiv_ids, max_indiv_id_len, family_size);
	  if (mm == -1) {
	    memcpy(cur_person_id, id_ptr, ujj + 1);
	    cur_person_id = &(cur_person_id[max_pm_id_len]);
	    stray_parent_ct++;
	    remaining_indiv_parent_idxs[uii * 2] = -2;
	  } else {
            remaining_indiv_parent_idxs[uii * 2] = fis_ptr[indiv_id_lookup[(uint32_t)mm]];
	  }
	} else {
          remaining_indiv_parent_idxs[uii * 2] = -1;
	}
	id_ptr = &(maternal_ids[((uint32_t)kk) * max_maternal_id_len]);
	if (memcmp("0", id_ptr, 2)) {
	  ujj = strlen(id_ptr);
          mm = bsearch_str(id_ptr, ujj, indiv_ids, max_indiv_id_len, family_size);
	  if (mm == -1) {
	    memcpy(cur_person_id, id_ptr, ujj + 1);
	    cur_person_id = &(cur_person_id[max_pm_id_len]);
	    stray_parent_ct++;
	    remaining_indiv_parent_idxs[uii * 2 + 1] = -2;
	  } else {
	    remaining_indiv_parent_idxs[uii * 2 + 1] = fis_ptr[indiv_id_lookup[(uint32_t)mm]];
	  }
	} else {
	  remaining_indiv_parent_idxs[uii * 2 + 1] = -1;
	}
        remaining_indiv_idxs[uii] = kk;
      }
      while (ii < family_size) {
	complete_indiv_idxs[complete_indiv_idx_ct++] = fis_ptr[(uint32_t)ii];
	ii++;
      }
      qsort(stray_parent_ids, stray_parent_ct, max_pm_id_len, strcmp_casted);
      cur_person_id = stray_parent_ids;
      ii = 0; // read idx
      jj = 0; // write idx

      // Now filter out all such parents who aren't referenced at least twice.
      while (ii + 1 < stray_parent_ct) {
        if (strcmp(&(stray_parent_ids[ii * max_pm_id_len]), &(stray_parent_ids[(ii + 1) * max_pm_id_len]))) {
	  ii++;
	  continue;
	}
	ii++;
	strcpy(cur_person_id, &(stray_parent_ids[ii * max_pm_id_len]));
	do {
	  ii++;
        } while (!(strcmp(cur_person_id, &(stray_parent_ids[ii * max_pm_id_len])) || (ii > stray_parent_ct)));
        cur_person_id = &(cur_person_id[max_pm_id_len]);
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

      // Now fill in remainder of remaining_indiv_parent_idxs.
      for (uii = 0; uii < remaining_indiv_ct; uii++) {
	jj = remaining_indiv_idxs[uii];
	if (remaining_indiv_parent_idxs[uii * 2] == -2) {
	  kk = bsearch_str_nl(&(paternal_ids[((uint32_t)jj) * max_paternal_id_len]), stray_parent_ids, max_pm_id_len, stray_parent_ct);
	  if (kk != -1) {
	    kk += unfiltered_indiv_ctlm;
	  }
	  remaining_indiv_parent_idxs[uii * 2] = kk;
	}
	if (remaining_indiv_parent_idxs[uii * 2 + 1] == -2) {
	  kk = bsearch_str_nl(&(maternal_ids[((uint32_t)jj) * max_maternal_id_len]), stray_parent_ids, max_pm_id_len, stray_parent_ct);
	  if (kk != -1) {
	    kk += unfiltered_indiv_ctlm;
	  }
	  remaining_indiv_parent_idxs[uii * 2 + 1] = kk;
	}
      }
      ullii = ((uint64_t)founder_ct) * (founder_ct - 1);
      while (remaining_indiv_ct) {
	indiv_idx_write = 0;
	for (uii = 0; uii < remaining_indiv_ct; uii++) {
	  kk = remaining_indiv_parent_idxs[uii * 2];
	  mm = remaining_indiv_parent_idxs[uii * 2 + 1];
	  jj = remaining_indiv_idxs[uii];
	  if (((kk == -1) || is_set(processed_indivs, kk)) && ((mm == -1) || is_set(processed_indivs, mm))) {
	    for (ujj = 0; ujj < founder_ct; ujj++) {
	      // relationship between kk and ujjth founder
	      if ((kk >= (int32_t)unfiltered_indiv_ct) || (kk == -1)) {
		dxx = 0.0;
	      } else if (is_set(founder_info, kk)) {
		if (kk == (int32_t)complete_indiv_idxs[ujj]) {
		  dxx = 0.5;
		} else {
		  dxx = 0.0;
		}
	      } else {
		ukk = pri_ptr->family_rel_nf_idxs[(uint32_t)kk];
                dxx = 0.5 * rs_ptr[((uint64_t)ukk * (ukk - 1) - ullii) / 2 + ujj];
	      }
	      if (is_set(founder_info, mm)) {
		if (mm == (int32_t)complete_indiv_idxs[ujj]) {
		  dxx += 0.5;
		}
	      } else if ((mm != -1) && (mm < (int32_t)unfiltered_indiv_ct)) {
		ukk = pri_ptr->family_rel_nf_idxs[(uint32_t)mm];
		dxx += 0.5 * rs_ptr[((uint64_t)ukk * (ukk - 1) - ullii) / 2 + ujj];
	      }
	      *rel_writer++ = dxx;
	    }
	    for (; ujj < complete_indiv_idx_ct; ujj++) {
	      if (kk == -1) {
		dxx = 0.0;
	      } else if (kk >= (int32_t)unfiltered_indiv_ct) {
		dxx = 0.5 * tmp_rel_space[(ujj - founder_ct) * stray_parent_ct + kk - unfiltered_indiv_ctlm];
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
	      if (mm >= (int32_t)unfiltered_indiv_ct) {
		dxx += 0.5 * tmp_rel_space[(ujj - founder_ct) * stray_parent_ct + mm - unfiltered_indiv_ctlm];
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
	      if (kk >= (int32_t)unfiltered_indiv_ct) {
		if ((uint32_t)kk == ujj + unfiltered_indiv_ctlm) {
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
	      if (mm >= (int32_t)unfiltered_indiv_ct) {
		if ((uint32_t)mm == ujj + unfiltered_indiv_ctlm) {
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
	    pri_ptr->family_rel_nf_idxs[jj] = complete_indiv_idx_ct;
	    complete_indiv_idxs[complete_indiv_idx_ct++] = jj;
	    set_bit(processed_indivs, jj);
	  } else {
            remaining_indiv_parent_idxs[indiv_idx_write * 2] = kk;
	    remaining_indiv_parent_idxs[indiv_idx_write * 2 + 1] = mm;
	    remaining_indiv_idxs[indiv_idx_write++] = jj;
	  }
	}
	if (indiv_idx_write == remaining_indiv_ct) {
	  logprint("Error: Pedigree graph is cyclic.  Check for evidence of time travel abuse in\nyour cohort.\n");
	  return RET_INVALID_FORMAT;
	}
	remaining_indiv_ct = indiv_idx_write;
      }
      wkspace_reset(wkspace_mark2);
    }
  }
  wkspace_reset(wkspace_mark);
  return 0;
}
