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

#include "plink_assoc.h"
#include "plink_cluster.h"
#include "plink_family.h"
#include "plink_perm.h"
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
  fam_ip->tucc_bed = 0;
}

// bottom 2 bits of index = child genotype
// middle 2 bits of index = paternal genotype
// top 2 bits of index = maternal genotype

// bits 0-7 = child increment (always 1)
// bits 8-15 = father increment
// bits 16-23 = mother increment
// bits 24-31 = female error code
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

// bottom 2 bits of index = child genotype
// next 2 bits of index = maternal genotype
const uint32_t mendel_error_table_male_x[] =
{0, 0, 0, 0x9010001,
 0, 0, 0, 0,
 0, 0, 0, 0,
 0xa010001, 0, 0, 0};

int32_t get_trios_and_families(uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, uintptr_t* founder_info, uintptr_t* sex_nm, uintptr_t* sex_male, char* sample_ids, uintptr_t max_sample_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, char** fids_ptr, uintptr_t* max_fid_len_ptr, char** iids_ptr, uintptr_t* max_iid_len_ptr, uint64_t** family_list_ptr, uint32_t* family_ct_ptr, uint64_t** trio_list_ptr, uintptr_t* trio_ct_ptr, uint32_t** trio_lookup_ptr, uint32_t include_duos, uint32_t toposort) {
  // This mirrors linkRelateds() in genedrop.cpp, and parseTrios() in trio.cpp,
  // in PLINK 1.07.
  //
  // family_list has paternal indices in low 32 bits, maternal indices in high
  // 32, sorted in child ID order.
  // trio_list has child IDs in low 32 bits, family_list indices in high 32
  // bits, in PLINK 1.07 reporting order.  (Actually, reporting order may
  // change, but this isn't a big deal.)
  // If include_duos is set, missing parental IDs are coded as
  // unfiltered_sample_ct.  include_duos must be 0 or 1.
  // If toposort is set, trio_lookup has child IDs in [4n], paternal IDs in
  // [4n+1], maternal IDs in [4n+2], and reporting sequence index in [4n+3].
  // Otherwise, it's [3n]/[3n+1]/[3n+2].  toposort must be 0 or 1.
  //
  // fids is a list of null-terminated FIDs using trio_list indices, and iids
  // is a list of IIDs using regular unfiltered indices.  If include_duos is
  // set, iids has a trailing entry set to '0'.  (fids_ptr, iids_ptr, and the
  // corresponding lengths can be nullptr.)
  //
  // PLINK 1.07 enforces <= 1 father and <= 1 mother per sample (and ambiguous
  // sex parents are not permitted), but the IDs CAN be reversed in the .fam
  // with no adverse consequences.  For backward compatibility, we replicate
  // this.  (Possible todo: report a warning exactly once when this happens.)
  // It won't be replicated in PLINK 2.0.
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  uint64_t* edge_list = nullptr;
  uint32_t* toposort_queue = nullptr;
  char* fids = nullptr;
  char* iids = nullptr;
  uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t unfiltered_sample_ctp1l = 1 + (unfiltered_sample_ct / BITCT);
  uintptr_t sample_uidx = next_unset_unsafe(sample_exclude, 0);
  // does *not* use populate_id_htable
  uintptr_t htable_size = geqprime(2 * unfiltered_sample_ct + 1);
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
  if (bigstack_end_alloc_ul(unfiltered_sample_ctp1l, &founder_info2) ||
      bigstack_end_alloc_ull(sample_ct, &trio_list_tmp) ||
      bigstack_end_alloc_c(sample_ct * max_sample_id_len, &sorted_sample_ids) ||
      bigstack_end_alloc_ui(sample_ct, &sample_id_map) ||
      bigstack_end_alloc_c(max_sample_id_len, &idbuf)) {
    goto get_trios_and_families_ret_NOMEM;
  }
  memcpy(founder_info2, founder_info, unfiltered_sample_ctl * sizeof(intptr_t));
  if (unfiltered_sample_ct & (BITCT - 1)) {
    SET_BIT(unfiltered_sample_ct, founder_info2);
  } else {
    founder_info2[unfiltered_sample_ctl] = 1;
  }
  if (sort_item_ids_noalloc(unfiltered_sample_ct, sample_exclude, sample_ct, sample_ids, max_sample_id_len, 0, 0, strcmp_deref, sorted_sample_ids, sample_id_map)) {
    goto get_trios_and_families_ret_1;
  }
  // over-allocate here, we shrink family_list later when we know how many
  // families there are
  if (bigstack_alloc_ull(sample_ct, &family_list) ||
      bigstack_calloc_ull(htable_size, &family_htable) ||
      bigstack_alloc_ui(htable_size, &family_idxs)) {
    goto get_trios_and_families_ret_NOMEM;
  }
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
	SET_BIT(sample_uidx, founder_info2);
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
	SET_BIT(sample_uidx, founder_info2);
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
  bigstack_reset(bigstack_mark);
  bigstack_alloc(family_ct * sizeof(int64_t)); // family_list
  bigstack_end_reset(trio_list_tmp);
  if (bigstack_alloc_ull(trio_ct, &trio_write)) {
    goto get_trios_and_families_ret_NOMEM;
  }
  memcpy(trio_write, trio_list_tmp, trio_ct * sizeof(int64_t));
#ifdef __cplusplus
  std::sort((int64_t*)trio_write, (int64_t*)(&(trio_write[trio_ct])));
#else
  qsort(trio_write, trio_ct, sizeof(int64_t), llcmp);
#endif
  bigstack_end_reset(founder_info2);
  if (bigstack_alloc_ui(trio_ct * (3 + toposort), &trio_lookup)) {
    goto get_trios_and_families_ret_NOMEM;
  }
  if (fids_ptr) {
    if (bigstack_alloc_c(trio_ct * max_fid_len, &fids) ||
	bigstack_alloc_c((unfiltered_sample_ct + include_duos) * max_iid_len, &iids)) {
      goto get_trios_and_families_ret_NOMEM;
    }
  }
  if (toposort) {
    if (bigstack_alloc_ull(trio_ct * 2, &edge_list)) {
      goto get_trios_and_families_ret_NOMEM;
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
    bigstack_shrink_top(edge_list, edge_ct * sizeof(int64_t));
    if (bigstack_alloc_ui(trio_ct, &toposort_queue)) {
      goto get_trios_and_families_ret_NOMEM;
    }
    remaining_edge_ct = edge_ct;
  }
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
        SET_BIT(sample_uidx, founder_info2);
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
      logerrprint("Error: Pedigree graph is cyclic.  Check for evidence of time travel abuse in\nyour cohort.\n");
      goto get_trios_and_families_ret_INVALID_FORMAT;
    }
    bigstack_reset(edge_list);
  }
  while (0) {
  get_trios_and_families_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  get_trios_and_families_ret_INVALID_FORMAT_2:
    logerrprintb();
  get_trios_and_families_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 get_trios_and_families_ret_1:
  if (retval) {
    bigstack_reset(bigstack_mark);
  }
  bigstack_end_reset(bigstack_end_mark);
  return retval;
}

uint32_t erase_mendel_errors(uintptr_t unfiltered_sample_ct, uintptr_t* loadbuf, uintptr_t* workbuf, uintptr_t* sex_male, uint32_t* trio_lookup, uint32_t trio_ct, uint32_t is_x, uint32_t multigen) {
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
  SET_BIT_DBL(unfiltered_sample_ct, workbuf);
  if (!multigen) {
    for (trio_idx = 0; trio_idx < trio_ct; trio_idx++) {
      uii = *uiptr++;
      ujj = *uiptr++;
      ukk = *uiptr++;
      umm = EXTRACT_2BIT_GENO(workbuf, uii);
      unn = EXTRACT_2BIT_GENO(workbuf, ukk);
      if ((!is_x) || (!is_set(sex_male, uii))) {
        umm = mendel_error_table[umm | (EXTRACT_2BIT_GENO(workbuf, ujj) << 2) | (unn << 4)];
      } else {
	umm = mendel_error_table_male_x[umm | (unn << 2)];
      }
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
      const uint32_t uqq = 2 * (uii % BITCT2);
      ujj = 2 * (ujj % BITCT2);
      ukk = 2 * (ukk % BITCT2);
      ulii = (workbuf[unn] >> uqq) & 3;
      uljj = ((workbuf[uoo] >> ujj) & 3) | (((workbuf[upp] >> ukk) & 3) << 2);
      if (ulii != 1) {
        // bugfix (10 Apr 2018): is_set(sex_male, uii) didn't work since we had
        // replaced the value with twice the low-order bits
	if ((!is_x) || (!is_set(sex_male, uii))) {
          umm = mendel_error_table[ulii | (uljj << 2)];
	} else {
	  umm = mendel_error_table_male_x[ulii | (uljj & 12)];
	}
        if (umm) {
	  ulii = loadbuf[unn];
	  uljj = ONELU << uqq;
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
	// both parents are homozygous for the same allele, so child genotype
	// is "known" for the purpose of checking grandchild genotypes
	workbuf[unn] &= ~(ONELU << uqq);
      } else if (uljj == 15) {
	workbuf[unn] |= (2 * ONELU) << uqq;
      }
      // no need to fill "known" heterozygous genotype, since that's treated
      // the same way as a missing genotype
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

int32_t mendel_error_scan(Family_info* fam_ip, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uint32_t plink_maxfid, uint32_t plink_maxiid, uint32_t plink_maxsnp, uint32_t allow_no_variants, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_exclude_ct_ptr, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t* sample_exclude_ct_ptr, uintptr_t* founder_info, uintptr_t* sex_nm, uintptr_t* sex_male, char* sample_ids, uintptr_t max_sample_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uint32_t hh_exists, Chrom_info* chrom_info_ptr, uint32_t calc_mendel) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  FILE* outfile_l = nullptr;
  uintptr_t* sample_male_include2 = nullptr;
  uintptr_t* error_locs = nullptr;
  char* varptr = nullptr;
  char* chrom_name_ptr = nullptr;
  unsigned char* cur_errors = nullptr;
  uint64_t* family_error_cts = nullptr;
  uint32_t* child_cts = nullptr;
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
  uint32_t full_error_list = calc_mendel && (!(fam_ip->mendel_modifier & MENDEL_SUMMARIES_ONLY));
  uint32_t varlen = 0;
  uint32_t chrom_name_len = 0;
  uint32_t new_marker_exclude_ct = 0;
  uint32_t error_ct_fill = 0;
  int32_t retval = 0;
  char chrom_name_buf[5];
  char* errstrs[10];
  uint32_t errstr_lens[11];
  uint32_t alens[2];
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
  uint32_t unn;
  marker_ct -= count_non_autosomal_markers(chrom_info_ptr, marker_exclude, 0, 1);
  if ((!marker_ct) || is_set(chrom_info_ptr->haploid_mask, 0)) {
    logerrprint("Warning: Skipping --me/--mendel since there is no autosomal or Xchr data.\n");
    goto mendel_error_scan_ret_1;
  }
  retval = get_trios_and_families(unfiltered_sample_ct, sample_exclude, sample_ct, founder_info, sex_nm, sex_male, sample_ids, max_sample_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, &fids, &max_fid_len, &iids, &max_iid_len, &family_list, &family_ct, &trio_list, &trio_ct, &trio_lookup, include_duos, multigen);
  if (retval) {
    goto mendel_error_scan_ret_1;
  }
  if (!trio_ct) {
    LOGERRPRINTF("Warning: Skipping --me/--mendel since there are no %strios.\n", include_duos? "duos or " : "");
    goto mendel_error_scan_ret_1;
  }
  if (family_ct > 0x55555555U) {
    // may as well document this limit
    logerrprint("Error: Too many families for --me/--mendel.\n");
    goto mendel_error_scan_ret_INVALID_CMDLINE;
  }

  trio_ct4 = (trio_ct + 3) / 4;
  trio_ctl = BITCT_TO_WORDCT(trio_ct);
  var_error_max = (int32_t)(fam_ip->mendel_max_var_error * (1 + SMALL_EPSILON) * ((intptr_t)trio_ct));
  if (bigstack_alloc_ul(unfiltered_sample_ctp1l2, &loadbuf) ||
      bigstack_calloc_ui(trio_ct * 3, &error_cts) ||
      bigstack_calloc_ui(trio_ct4 * 4, &error_cts_tmp)) {
    goto mendel_error_scan_ret_NOMEM;
  }
  if (!var_first) {
    error_cts_tmp2 = error_cts_tmp;
  } else {
    if (bigstack_calloc_ui(trio_ct4 * 4, &error_cts_tmp2)) {
      goto mendel_error_scan_ret_NOMEM;
    }
  }
  loadbuf[unfiltered_sample_ctp1l2 - 1] = 0;
  hh_exists &= XMHH_EXISTS;
  if (alloc_raw_haploid_filters(unfiltered_sample_ct, hh_exists, 0, sample_exclude, sex_male, nullptr, &sample_male_include2)) {
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
    if (bigstack_alloc_ull(family_ct * 3, &family_error_cts) ||
        bigstack_alloc_ui(family_ct, &child_cts) ||
        bigstack_alloc_c(ulii * 10, &(errstrs[0]))) {
      goto mendel_error_scan_ret_NOMEM;
    }
    for (uii = 1; uii < 10; uii++) {
      errstrs[uii] = &(errstrs[0][uii * ulii]);
    }
    if (multigen && full_error_list) {
      if (bigstack_calloc_ul(trio_ctl, &error_locs) ||
	  bigstack_alloc_uc(trio_ct, &cur_errors)) {
	goto mendel_error_scan_ret_NOMEM;
      }
    }
    if (full_error_list) {
      memcpy(outname_end, ".mendel", 8);
      if (fopen_checked(outname, "w", &outfile)) {
	goto mendel_error_scan_ret_OPEN_FAIL;
      }
      sprintf(g_textbuf, "%%%us %%%us  CHR %%%us   CODE                 ERROR\n", plink_maxfid, plink_maxiid, plink_maxsnp);
      fprintf(outfile, g_textbuf, "FID", "KID", "SNP");
    }
    memcpy(outname_end, ".lmendel", 9);
    if (fopen_checked(outname, "w", &outfile_l)) {
      goto mendel_error_scan_ret_OPEN_FAIL;
    }
    // replicate harmless 'N' misalignment bug
    sprintf(g_textbuf, " CHR %%%us   N\n", plink_maxsnp);
    fprintf(outfile_l, g_textbuf, "SNP");
  } else {
    // suppress warning
    fill_ulong_zero(10, (uintptr_t*)errstrs);
  }
  for (chrom_fo_idx = 0; chrom_fo_idx < chrom_info_ptr->chrom_ct; chrom_fo_idx++) {
    chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
    is_x = (((uint32_t)chrom_info_ptr->xymt_codes[X_OFFSET]) == chrom_idx);
    if ((IS_SET(chrom_info_ptr->haploid_mask, chrom_idx) && (!is_x)) || (((uint32_t)chrom_info_ptr->xymt_codes[MT_OFFSET]) == chrom_idx)) {
      continue;
    }
    chrom_end = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx + 1];
    uii = next_unset(marker_exclude, chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx], chrom_end);
    if (uii == chrom_end) {
      continue;
    }
    if (calc_mendel) {
      chrom_name_ptr = chrom_name_buf5w4write(chrom_info_ptr, chrom_idx, &chrom_name_len, chrom_name_buf);
    }
    if (uii != marker_uidx) {
      marker_uidx = uii;
      goto mendel_error_scan_seek;
    }
    while (1) {
      if (load_raw2(unfiltered_sample_ct4, unfiltered_sample_ctl2m1, final_mask, bedfile, loadbuf)) {
	goto mendel_error_scan_ret_READ_FAIL;
      }
      if (IS_SET(marker_reverse, marker_uidx)) {
        reverse_loadbuf(unfiltered_sample_ct, (unsigned char*)loadbuf);
      }
      if (hh_exists && is_x) {
	hh_reset((unsigned char*)loadbuf, sample_male_include2, unfiltered_sample_ct);
      }
      // missing parents are treated as having uidx equal to
      // unfiltered_sample_ct, and we set the corresponding genotype to always
      // be missing.  This lets us avoid special-casing duos.
      SET_BIT_DBL(unfiltered_sample_ct, loadbuf);
      uiptr = trio_lookup;
      cur_error_ct = 0;
      if (calc_mendel) {
	varptr = &(marker_ids[marker_uidx * max_marker_id_len]);
	varlen = strlen(varptr);
	alens[0] = 0;
	alens[1] = 0;
	fill_uint_zero(11, errstr_lens);
      }
      if (!multigen) {
	for (trio_idx = 0; trio_idx < trio_ct; trio_idx++) {
	  uii = *uiptr++;
	  ujj = *uiptr++;
	  ukk = *uiptr++;
	  umm = EXTRACT_2BIT_GENO(loadbuf, uii);
	  unn = EXTRACT_2BIT_GENO(loadbuf, ukk);
	  if ((!is_x) || (!is_set(sex_male, uii))) {
            umm = mendel_error_table[umm | (EXTRACT_2BIT_GENO(loadbuf, ujj) << 2) | (unn << 4)];
	  } else {
	    umm = mendel_error_table_male_x[umm | (unn << 2)];
	  }
	  if (umm) {
	    error_cts_tmp2[trio_idx] += umm & 0xffffff;
	    cur_error_ct++;
	    if (full_error_list) {
	      umm >>= 24;
	      wptr = fw_strcpy(plink_maxfid, &(fids[trio_idx * max_fid_len]), g_textbuf);
	      *wptr++ = ' ';
	      wptr = fw_strcpy(plink_maxiid, &(iids[uii * max_iid_len]), wptr);
	      *wptr++ = ' ';
	      wptr = memcpyax(wptr, chrom_name_ptr, chrom_name_len, ' ');
	      wptr = fw_strcpyn(plink_maxsnp, varlen, varptr, wptr);
	      wptr = memseta(wptr, 32, 5);
	      if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
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
          uljj = EXTRACT_2BIT_GENO(loadbuf, ujj) | (EXTRACT_2BIT_GENO(loadbuf, ukk) << 2);
	  umm = uii / BITCT2;
	  ujj = 2 * (uii % BITCT2);
	  ulii = (loadbuf[umm] >> ujj) & 3;
	  if (ulii != 1) {
	    if ((!is_x) || (!is_set(sex_male, uii))) {
	      umm = mendel_error_table[ulii | (uljj << 2)];
	    } else {
	      umm = mendel_error_table_male_x[ulii | (uljj & 12)];
	    }
	    if (umm) {
	      error_cts_tmp2[trio_idx] += umm & 0xffffff;
	      cur_error_ct++;
	      if (full_error_list) {
	        set_bit(trio_idx, error_locs);
		umm >>= 24;
                cur_errors[trio_idx] = (unsigned char)umm;
	      }
	    }
	  } else if (!uljj) {
	    loadbuf[umm] &= ~(ONELU << ujj);
          } else if (uljj == 15) {
	    loadbuf[umm] |= (2 * ONELU) << ujj;
	  }
	}
	if (full_error_list && cur_error_ct) {
          trio_idx = 0;
	  for (uii = 0; uii < cur_error_ct; trio_idx++, uii++) {
            next_set_ul_unsafe_ck(error_locs, &trio_idx);
	    wptr = fw_strcpy(plink_maxfid, &(fids[trio_idx * max_fid_len]), g_textbuf);
	    *wptr++ = ' ';
	    wptr = fw_strcpy(plink_maxiid, &(iids[((uint32_t)trio_list[trio_idx]) * max_iid_len]), wptr);
	    *wptr++ = ' ';
	    wptr = memcpyax(wptr, chrom_name_ptr, chrom_name_len, ' ');
	    wptr = fw_strcpyn(plink_maxsnp, varlen, varptr, wptr);
	    wptr = memseta(wptr, 32, 5);
	    if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
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
          fill_ulong_zero(trio_ctl, error_locs);
	}
      }
      if (calc_mendel) {
	if (fwrite_checked(chrom_name_ptr, chrom_name_len, outfile_l)) {
	  goto mendel_error_scan_ret_WRITE_FAIL;
	}
	g_textbuf[0] = ' ';
	wptr = fw_strcpyn(plink_maxsnp, varlen, varptr, &(g_textbuf[1]));
        *wptr++ = ' ';
        wptr = uint32toa_w4x(cur_error_ct, '\n', wptr);
	if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile_l)) {
	  goto mendel_error_scan_ret_WRITE_FAIL;
	}
      }
      if (cur_error_ct) {
	if (cur_error_ct > var_error_max) {
	  SET_BIT(marker_uidx, marker_exclude);
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
	    fill_uint_zero(trio_ct, error_cts_tmp2);
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
	    fill_uint_zero(trio_ct, error_cts_tmp);
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
    if (full_error_list) {
      if (fclose_null(&outfile)) {
	goto mendel_error_scan_ret_WRITE_FAIL;
      }
    }
    if (fclose_null(&outfile_l)) {
      goto mendel_error_scan_ret_WRITE_FAIL;
    }
    outname_end[1] = 'f';
    if (fopen_checked(outname, "w", &outfile)) {
      goto mendel_error_scan_ret_OPEN_FAIL;
    }
    sprintf(g_textbuf, "%%%us %%%us %%%us   CHLD    N\n", plink_maxfid, plink_maxiid, plink_maxiid);
    fprintf(outfile, g_textbuf, "FID", "PAT", "MAT");
    fill_ull_zero(family_ct * 3, family_error_cts);
    fill_uint_zero(family_ct, child_cts);
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
	wptr = fw_strcpyn(plink_maxfid, (uintptr_t)(((char*)memchr(cptr, '\t', max_sample_id_len)) - cptr), cptr, g_textbuf);
      } else {
	cptr = &(sample_ids[ukk * max_sample_id_len]);
	wptr = fw_strcpyn(plink_maxfid, (uintptr_t)(((char*)memchr(cptr, '\t', max_sample_id_len)) - cptr), cptr, g_textbuf);
	// wptr = memseta(g_textbuf, 32, plink_maxfid - 1);
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
      wptr = uint32toa_w6x(child_cts[uii], ' ', wptr);
      if (family_error_cts[uii * 3] < 10000) {
	wptr = uint32toa_w4((uint32_t)family_error_cts[uii * 3], wptr);
      } else {
        wptr = int64toa(family_error_cts[uii * 3], wptr);
      }
      *wptr++ = '\n';
      if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
	goto mendel_error_scan_ret_WRITE_FAIL;
      }
    }
    if (fclose_null(&outfile)) {
      goto mendel_error_scan_ret_WRITE_FAIL;
    }
    outname_end[1] = 'i';
    if (fopen_checked(outname, "w", &outfile)) {
      goto mendel_error_scan_ret_OPEN_FAIL;
    }
    sprintf(g_textbuf, "%%%us %%%us   N\n", plink_maxfid, plink_maxiid);
    fprintf(outfile, g_textbuf, "FID", "IID");
    uii = 0xffffffffU; // family idx
    for (trio_idx = 0; trio_idx < trio_ct; trio_idx++) {
      trio_code = trio_list[trio_idx];
      ujj = (uint32_t)(trio_code >> 32);
      if (ujj != uii) {
	uii = ujj;
        family_code = family_list[uii];
	wptr = fw_strcpy(plink_maxfid, &(fids[trio_idx * max_fid_len]), g_textbuf);
	*wptr++ = ' ';
	ujj = (uint32_t)family_code;
	if (ujj != unfiltered_sample_ct) {
	  wptr = fw_strcpy(plink_maxiid, &(iids[ujj * max_iid_len]), wptr);
	  *wptr++ = ' ';
	  if (family_error_cts[3 * uii + 1] < 10000) {
	    wptr = uint32toa_w4((uint32_t)family_error_cts[3 * uii + 1], wptr);
	  } else {
	    wptr = int64toa(family_error_cts[3 * uii + 1], wptr);
	  }
	  if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
	    goto mendel_error_scan_ret_WRITE_FAIL;
	  }
	}
	ukk = (uint32_t)(family_code >> 32);
	if (ukk != unfiltered_sample_ct) {
	  if (ujj != unfiltered_sample_ct) {
	    putc_unlocked('\n', outfile);
	  }
	  wptr = fw_strcpy(plink_maxiid, &(iids[ukk * max_iid_len]), &(g_textbuf[plink_maxfid + 1]));
	  *wptr++ = ' ';
	  if (family_error_cts[3 * uii + 2] < 10000) {
	    wptr = uint32toa_w4((uint32_t)family_error_cts[3 * uii + 2], wptr);
	  } else {
	    wptr = int64toa(family_error_cts[3 * uii + 2], wptr);
	  }
	  if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
	    goto mendel_error_scan_ret_WRITE_FAIL;
	  }
	}
	putc_unlocked(' ', outfile); // PLINK 1.07 formatting quirk
	putc_unlocked('\n', outfile);
      }
      wptr = fw_strcpy(plink_maxiid, &(iids[((uint32_t)trio_code) * max_iid_len]), &(g_textbuf[plink_maxfid + 1]));
      *wptr++ = ' ';
      wptr = uint32toa_w4x(error_cts[trio_idx * 3], '\n', wptr);
      if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
	goto mendel_error_scan_ret_WRITE_FAIL;
      }
    }
    *outname_end = '\0';
    LOGPRINTFWW("Reports written to %s%s%s.imendel + %s.fmendel + %s.lmendel .\n", full_error_list? outname : "", full_error_list? ".mendel + " : "", outname, outname, outname);
  }
  if (fam_ip->mendel_modifier & MENDEL_FILTER) {
    *marker_exclude_ct_ptr += new_marker_exclude_ct;
    if ((unfiltered_marker_ct == *marker_exclude_ct_ptr) && (!allow_no_variants)) {
      logerrprint("Error: All variants excluded by --me.\n");
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
	    set_bit((uint32_t)trio_code, sample_exclude);
	    if (ujj < unfiltered_sample_ct) {
	      set_bit(ujj, sample_exclude);
	    }
	    if (ukk < unfiltered_sample_ct) {
	      set_bit(ukk, sample_exclude);
	    }
	  } else if ((exclude_one_ratio == -1) || (ujj == unfiltered_sample_ct) || (ukk == unfiltered_sample_ct)) {
            set_bit((uint32_t)trio_code, sample_exclude);
	  } else {
	    dxx = (double)((int32_t)trio_list[trio_idx * 3 + 1]);
	    dyy = (double)((int32_t)trio_list[trio_idx * 3 + 2]);
	    if (dxx > exclude_one_ratio * dyy) {
	      set_bit(ujj, sample_exclude);
	    } else if (dyy > exclude_one_ratio * dxx) {
	      set_bit(ukk, sample_exclude);
	    } else {
	      set_bit((uint32_t)trio_code, sample_exclude);
	    }
	  }
	}
      }
    }
    ulii = popcount_longs(sample_exclude, BITCT_TO_WORDCT(unfiltered_sample_ct));
    if (unfiltered_sample_ct == ulii) {
      LOGERRPRINTF("Error: All %s excluded by --me.\n", g_species_plural);
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
  mendel_error_scan_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  mendel_error_scan_ret_ALL_MARKERS_EXCLUDED:
    retval = RET_ALL_MARKERS_EXCLUDED;
    break;
  mendel_error_scan_ret_ALL_SAMPLES_EXCLUDED:
    retval = RET_ALL_SAMPLES_EXCLUDED;
    break;
  }
 mendel_error_scan_ret_1:
  bigstack_reset(bigstack_mark);
  fclose_cond(outfile);
  fclose_cond(outfile_l);
  return retval;
}

int32_t populate_pedigree_rel_info(Pedigree_rel_info* pri_ptr, uintptr_t unfiltered_sample_ct, char* sample_ids, uintptr_t max_sample_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* founder_info) {
  // possible todo: if any families have been entirely filtered out, don't
  // construct pedigree for them
  uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t unfiltered_sample_ctlm = unfiltered_sample_ctl * BITCT;
  uintptr_t max_family_id_len = 0;
  uintptr_t max_indiv_id_len = 0;
  uintptr_t max_pm_id_len = MAXV(max_paternal_id_len, max_maternal_id_len);
  char* last_family_id = nullptr;
  double* tmp_rel_space = nullptr;
  double* tmp_rel_writer = nullptr;
  uint32_t* uiptr2 = nullptr;
  int32_t max_family_nf = 0;
  unsigned char* bigstack_mark;
  unsigned char* bigstack_mark2;
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
  char* cur_family_id;
  char* id_ptr;
  uint32_t* family_sizes;
  uint32_t* uiptr;
  uint32_t fidx;
  int32_t family_size;
  uint32_t* remaining_sample_idxs;
  int32_t* remaining_sample_parent_idxs; // -1 = no parent (or nonshared)
  uint32_t remaining_sample_ct;
  uint32_t sample_idx_write;
  char* indiv_ids; // within a single family
  uint32_t* sample_id_lookup;
  uint32_t family_id_ct;
  uint32_t* fis_ptr;
  char* stray_parent_ids;
  intptr_t stray_parent_ct;
  uintptr_t* processed_samples;
  uint32_t founder_ct;
  uint32_t* complete_sample_idxs;
  uintptr_t complete_sample_idx_ct;
  double* rs_ptr;
  double* rel_writer;
  double dxx;

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
  if (bigstack_alloc_ui(unfiltered_sample_ct, &(pri_ptr->family_info_space)) ||
      bigstack_alloc_ui(unfiltered_sample_ct, &(pri_ptr->family_rel_nf_idxs)) ||
      bigstack_alloc_ui(unfiltered_sample_ct, &(pri_ptr->family_idxs)) ||
      bigstack_alloc_c(unfiltered_sample_ct * max_family_id_len, &family_ids) ||
      bigstack_alloc_ui(unfiltered_sample_ct, &family_sizes)) {
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
    bigstack_shrink_top(family_ids, family_id_ct * max_family_id_len);
    bigstack_alloc_ui(family_id_ct, &family_sizes);
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

  if (bigstack_alloc_ui(family_id_ct + 1, &(pri_ptr->family_info_offsets)) ||
      bigstack_alloc_ul(family_id_ct + 1, &(pri_ptr->family_rel_space_offsets)) ||
      bigstack_calloc_ui(family_id_ct, &(pri_ptr->family_founder_cts))) {
    return RET_NOMEM;
  }

  ii = 0; // running family_info offset
  for (fidx = 0; fidx < family_id_ct; fidx++) {
    family_size = pri_ptr->family_sizes[fidx];
    pri_ptr->family_info_offsets[fidx] = ii;
    ii += family_size;
  }

  if (bigstack_calloc_ui(family_id_ct, &uiptr)) {
    return RET_NOMEM;
  }

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
  bigstack_reset(uiptr);
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
  if (bigstack_alloc_d(ulii, &(pri_ptr->rel_space))) {
    return RET_NOMEM;
  }

  bigstack_mark = g_bigstack_base;
  if (bigstack_alloc_ui(family_id_ct, &uiptr)) {
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
  bigstack_reset(bigstack_mark);

  ulii = QUATERCT_TO_WORDCT(max_family_nf);
  if (bigstack_alloc_ul(unfiltered_sample_ctl + ulii, &processed_samples)) {
    return RET_NOMEM;
  }
  fill_ulong_one(ulii, &(processed_samples[unfiltered_sample_ctl]));

  bigstack_mark2 = g_bigstack_base;
  for (fidx = 0; fidx < family_id_ct; fidx++) {
    family_size = family_sizes[fidx];
    founder_ct = pri_ptr->family_founder_cts[fidx];
    remaining_sample_ct = family_size - founder_ct;
    stray_parent_ct = 0;
    if (remaining_sample_ct) {
      memcpy(processed_samples, founder_info, unfiltered_sample_ctl * sizeof(intptr_t));
      if (bigstack_alloc_ui(family_size, &complete_sample_idxs) ||
          bigstack_alloc_ui(remaining_sample_ct, &remaining_sample_idxs) ||
          bigstack_alloc_c(family_size * max_indiv_id_len, &indiv_ids) ||
          bigstack_alloc_ui(family_size, &sample_id_lookup) ||
          bigstack_alloc_i(remaining_sample_ct * 2, &remaining_sample_parent_idxs) ||
          bigstack_alloc_c(remaining_sample_ct * 2 * max_pm_id_len, &stray_parent_ids)) {
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
        if (bigstack_alloc_d((family_size - founder_ct) * stray_parent_ct, &tmp_rel_space)) {
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
	    set_bit(jj, processed_samples);
	  } else {
            remaining_sample_parent_idxs[sample_idx_write * 2] = kk;
	    remaining_sample_parent_idxs[sample_idx_write * 2 + 1] = mm;
	    remaining_sample_idxs[sample_idx_write++] = jj;
	  }
	}
	if (sample_idx_write == remaining_sample_ct) {
	  logerrprint("Error: Pedigree graph is cyclic.  Check for evidence of time travel abuse in\nyour cohort.\n");
	  return RET_INVALID_FORMAT;
	}
	remaining_sample_ct = sample_idx_write;
      }
      bigstack_reset(bigstack_mark2);
    }
  }
  bigstack_reset(bigstack_mark);
  return 0;
}

int32_t tdt_poo(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, double output_min_p, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct_ax, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, uintptr_t* marker_reverse, uintptr_t unfiltered_sample_ct, uintptr_t* sex_male, uintptr_t* sample_male_include2, uint32_t* trio_nuclear_lookup, uint32_t family_ct, Aperm_info* apip, uint32_t mperm_save, char* sample_ids, uintptr_t max_sample_id_len, Chrom_info* chrom_info_ptr, uint32_t hh_exists, Family_info* fam_ip, uintptr_t* loadbuf, uintptr_t* workbuf, char* textbuf, double* orig_chisq, uint32_t* trio_error_lookup, uintptr_t trio_ct) {
  FILE* outfile = nullptr;
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
  // bugfix (10 Oct 2018): missed a few chrX possibilities
  const uint32_t poo_table[] =
    {0, 0, 0, 0,
     0, 0, 0, 0,
     0x20002, 0, 2, 0,
     0, 0, 0, 0,
     0, 0, 0, 0,
     0, 0, 0, 0,
     0, 0, 0, 0,
     0, 0, 0, 0,
     0x2000200, 0, 0x200, 0x200,
     0, 0, 0, 0,
     0x2020202, 0, 0x1010202, 0x202,
     0x2000200, 0, 0x2000200, 0x200,
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
  if (fopen_checked(outname, "w", &outfile)) {
    goto tdt_poo_ret_OPEN_FAIL;
  }
  sprintf(textbuf, " CHR %%%us  A1:A2      T:U_PAT    CHISQ_PAT        P_PAT      T:U_MAT    CHISQ_MAT        P_MAT        Z_POO        P_POO \n", plink_maxsnp);
  fprintf(outfile, textbuf, "SNP");
  fputs("--tdt poo: 0%", stdout);
  fflush(stdout);
  for (chrom_fo_idx = 0; chrom_fo_idx < chrom_info_ptr->chrom_ct; chrom_fo_idx++) {
    chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
    is_x = ((int32_t)chrom_idx == chrom_info_ptr->xymt_codes[X_OFFSET]);
    if ((IS_SET(chrom_info_ptr->haploid_mask, chrom_idx) && (!is_x)) || (((uint32_t)chrom_info_ptr->xymt_codes[MT_OFFSET]) == chrom_idx)) {
      continue;
    }
    chrom_end = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx + 1];
    uii = next_unset(marker_exclude, chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx], chrom_end);
    if (uii == chrom_end) {
      continue;
    }
    wptr_start = width_force(4, textbuf, chrom_name_write(chrom_info_ptr, chrom_idx, textbuf));
    *wptr_start++ = ' ';
    if (uii != marker_uidx) {
      marker_uidx = uii;
      goto tdt_poo_scan_seek;
    }
    while (1) {
      if (load_raw2(unfiltered_sample_ct4, unfiltered_sample_ctl2m1, final_mask, bedfile, loadbuf)) {
	goto tdt_poo_ret_READ_FAIL;
      }
      if (IS_SET(marker_reverse, marker_uidx)) {
	reverse_loadbuf(unfiltered_sample_ct, (unsigned char*)loadbuf);
      }
      if (hh_exists && is_x) {
        hh_reset((unsigned char*)loadbuf, sample_male_include2, unfiltered_sample_ct);
      }
      mendel_error_ct += erase_mendel_errors(unfiltered_sample_ct, loadbuf, workbuf, sex_male, trio_error_lookup, trio_ct, is_x, multigen);
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
        ulii = EXTRACT_2BIT_GENO(loadbuf, uii);
        uljj = EXTRACT_2BIT_GENO(loadbuf, ujj);
        ukk = ulii | (uljj << 2);
	if ((0x4d04 >> ukk) & 1) {
	  // 1+ het parents, no missing
	  // xor doesn't work here since we need to track fathers and mothers
	  // separately
	  poo_table_ptr = &(poo_table[4 * ukk]);
          for (child_idx = 0; child_idx < cur_child_ct; child_idx++) {
            ukk = *lookup_ptr++;
            poo_acc += poo_table_ptr[EXTRACT_2BIT_GENO(loadbuf, ukk)];
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
      wptr2 = dtoa_g_wxp4x(pat_a1transmit, 1, ':', wptr);
      wptr2 = dtoa_g_wxp4(cur_a2transmit, 1, wptr2);
      wptr = width_force(12, wptr, wptr2);
      *wptr++ = ' ';
      if (poo_obs_pat_x2) {
	pat_a2transmit_recip = 1.0 / cur_a2transmit;
	dxx = pat_a1transmit - cur_a2transmit;
	chisq = dxx * dxx / (pat_a1transmit + cur_a2transmit);
	wptr = dtoa_g_wxp4x(chisq, 12, ' ', wptr);
	wptr = dtoa_g_wxp4(chiprob_p(chisq, 1), 12, wptr);
      } else {
	wptr = memcpya(wptr, "          NA           NA", 25);
      }
      *wptr++ = ' ';
      dxx = 0.5 * ((double)poo_mat_a1transmit_x2);
      cur_a2transmit = 0.5 * ((double)(poo_obs_mat_x2 - poo_mat_a1transmit_x2));
      wptr2 = dtoa_g_wxp4x(dxx, 1, ':', wptr);
      wptr2 = dtoa_g_wxp4(cur_a2transmit, 1, wptr2);
      wptr = width_force(12, wptr, wptr2);
      *wptr++ = ' ';
      if (poo_obs_mat_x2) {
	mat_a1transmit_recip = 1.0 / dxx;
	chisq = dxx - cur_a2transmit;
	chisq = chisq * chisq / (dxx + cur_a2transmit);
	wptr = dtoa_g_wxp4x(chisq, 12, ' ', wptr);
	wptr = dtoa_g_wxp4(chiprob_p(chisq, 1), 12, wptr);
      } else {
	wptr = memcpya(wptr, "          NA           NA", 25);
      }
      *wptr++ = ' ';
      if (poo_pat_a1transmit_x2 && poo_mat_a1transmit_x2 && (poo_obs_pat_x2 > poo_pat_a1transmit_x2) && (poo_obs_mat_x2 > poo_mat_a1transmit_x2)) {
	// Z-score
	dxx = (log(pat_a1transmit * pat_a2transmit_recip * mat_a1transmit_recip * cur_a2transmit) / sqrt(1.0 / pat_a1transmit + pat_a2transmit_recip + mat_a1transmit_recip + 1.0 / cur_a2transmit));

        wptr = dtoa_g_wxp4x(dxx, 12, ' ', wptr);
	if (orig_chisq) {
	  // todo: --pat/--mat support
	  orig_chisq[markers_done] = dxx * dxx;
	}
	dxx = normdist(-fabs(dxx)) * 2;
	wptr = dtoa_g_wxp4(MAXV(dxx, output_min_p), 12, wptr);
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
	  putc_unlocked('\b', stdout);
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
  putc_unlocked('\r', stdout);
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
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  char* textbuf = g_textbuf;
  double* orig_chisq = nullptr; // pval if exact test
  uint64_t last_parents = 0;
  // uint64_t mendel_error_ct = 0;
  double chisq = 0;
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t unfiltered_sample_ctl2 = QUATERCT_TO_WORDCT(unfiltered_sample_ct);
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
    logerrprint("Warning: Skipping --tdt since there is no autosomal or Xchr data.\n");
    goto tdt_ret_1;
  }
  retval = get_trios_and_families(unfiltered_sample_ct, sample_exclude, sample_ct, founder_info, sex_nm, sex_male, sample_ids, max_sample_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, &fids, &max_fid_len, &iids, &max_iid_len, &family_list, &family_ct, &trio_list, &trio_ct, &trio_error_lookup, 0, multigen);
  if (retval) {
    goto tdt_ret_1;
  }
  if (!trio_ct) {
    logerrprint("Warning: Skipping --tdt since there are no trios.\n");
    goto tdt_ret_1;
  }
  // now assemble list of nuclear families with at least one case child
  if (bigstack_alloc_ui(3LU * family_ct + trio_ct, &trio_nuclear_lookup)) {
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
	  *lookup_ptr = 0x80000000U;
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
    LOGERRPRINTF("Warning: Skipping --tdt%s since there are no trios with an affected child%s.\n", poo_test? " poo" : "", poo_test? "" : ", and no\ndiscordant parent pairs");
    goto tdt_ret_1;
  }
  bigstack_shrink_top(trio_nuclear_lookup, ((uintptr_t)(lookup_ptr - trio_nuclear_lookup)) * sizeof(int32_t));

  if (mtest_adjust) {
    if (bigstack_alloc_d(marker_ct, &orig_chisq)) {
      goto tdt_ret_NOMEM;
    }
  }
  if (bigstack_alloc_ul(unfiltered_sample_ctl2, &loadbuf) ||
      bigstack_alloc_ul(unfiltered_sample_ctp1l2, &workbuf)) {
    goto tdt_ret_NOMEM;
  }
  loadbuf[unfiltered_sample_ctl2 - 1] = 0;
  workbuf[unfiltered_sample_ctp1l2 - 1] = 0;
  hh_exists &= XMHH_EXISTS;
  if (alloc_raw_haploid_filters(unfiltered_sample_ct, hh_exists, 1, sample_exclude, sex_male, nullptr, &sample_male_include2)) {
    goto tdt_ret_NOMEM;
  }
  if (fam_ip->tdt_modifier & (TDT_PERM | TDT_MPERM)) {
    logerrprint("Error: --tdt permutation tests are currently under development.\n");
    retval = RET_CALC_NOT_YET_SUPPORTED;
    goto tdt_ret_1;
  }
  ulii = 2 * max_marker_allele_len + plink_maxsnp + MAX_ID_SLEN + 256;
  if (ulii > MAXLINELEN) {
    if (bigstack_alloc_c(ulii, &textbuf)) {
      goto tdt_ret_NOMEM;
    }
  }
  if (poo_test) {
    retval = tdt_poo(threads, bedfile, bed_offset, outname, outname_end, output_min_p, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, plink_maxsnp, marker_allele_ptrs, max_marker_allele_len, marker_reverse, unfiltered_sample_ct, sex_male, sample_male_include2, trio_nuclear_lookup, family_ct, apip, mperm_save, sample_ids, max_sample_id_len, chrom_info_ptr, hh_exists, fam_ip, loadbuf, workbuf, textbuf, orig_chisq, trio_error_lookup, trio_ct);
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
  if (fopen_checked(outname, "w", &outfile)) {
    goto tdt_ret_OPEN_FAIL;
  }
  sprintf(textbuf, " CHR %%%us           BP  A1  A2      T      U           OR ", plink_maxsnp);
  fprintf(outfile, textbuf, "SNP");
  if (display_ci) {
    uii = (uint32_t)((int32_t)(ci_size * (100 + EPSILON)));
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
    is_x = ((int32_t)chrom_idx == chrom_info_ptr->xymt_codes[X_OFFSET]);
    if ((IS_SET(chrom_info_ptr->haploid_mask, chrom_idx) && (!is_x)) || (((uint32_t)chrom_info_ptr->xymt_codes[MT_OFFSET]) == chrom_idx)) {
      continue;
    }
    chrom_end = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx + 1];
    uii = next_unset(marker_exclude, chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx], chrom_end);
    if (uii == chrom_end) {
      continue;
    }
    wptr_start = width_force(4, textbuf, chrom_name_write(chrom_info_ptr, chrom_idx, textbuf));
    *wptr_start++ = ' ';
    if (uii != marker_uidx) {
      marker_uidx = uii;
      goto tdt_scan_seek;
    }
    while (1) {
      if (load_raw2(unfiltered_sample_ct4, unfiltered_sample_ctl2m1, final_mask, bedfile, loadbuf)) {
	goto tdt_ret_READ_FAIL;
      }
      if (IS_SET(marker_reverse, marker_uidx)) {
	reverse_loadbuf(unfiltered_sample_ct, (unsigned char*)loadbuf);
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
      erase_mendel_errors(unfiltered_sample_ct, loadbuf, workbuf, sex_male, trio_error_lookup, trio_ct, is_x, multigen);
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
        ulii = EXTRACT_2BIT_GENO(loadbuf, uii);
        uljj = EXTRACT_2BIT_GENO(loadbuf, ujj);
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
	    umm = tdt_table_ptr[EXTRACT_2BIT_GENO(loadbuf, ukk)];
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
	wptr = uint32toa_w10x(marker_pos[marker_uidx], ' ', wptr);
	wptr = fw_strcpy(3, marker_allele_ptrs[2 * marker_uidx], wptr);
	*wptr++ = ' ';
	wptr = fw_strcpy(3, marker_allele_ptrs[2 * marker_uidx + 1], wptr);
	*wptr++ = ' ';
	wptr = uint32toa_w6x(tdt_a1_trans_ct, ' ', wptr);
	uii = tdt_obs_ct - tdt_a1_trans_ct; // untransmitted
	wptr = uint32toa_w6x(uii, ' ', wptr);
	if (uii) {
	  untransmitted_recip = 1.0 / ((double)((int32_t)uii));
	  dxx = (double)((int32_t)tdt_a1_trans_ct);
	  odds_ratio = dxx * untransmitted_recip;
	  wptr = dtoa_g_wxp4x(odds_ratio, 12, ' ', wptr);
	  if (display_ci) {
	    odds_ratio = log(odds_ratio);
	    dxx = ci_zt * sqrt(1.0 / dxx + untransmitted_recip);
	    wptr = dtoa_g_wxp4x(exp(odds_ratio - dxx), 12, ' ', wptr);
	    wptr = dtoa_g_wxp4x(exp(odds_ratio + dxx), 12, ' ', wptr);
	  }
	} else {
	  wptr = memcpya(wptr, "          NA ", 13);
	  if (display_ci) {
	    wptr = memcpya(wptr, "          NA           NA ", 26);
	  }
	}
        if (is_exact) {
	  wptr = dtoa_g_wxp4x(MAXV(pval, output_min_p), 12, ' ', wptr);
	} else {
	  if (pval >= 0) {
	    wptr = dtoa_g_wxp4x(chisq, 12, ' ', wptr);
            wptr = dtoa_g_wxp4x(MAXV(pval, output_min_p), 12, ' ', wptr);
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
	  wptr2 = uint32toa_x(uii, ':', wptr);
	  wptr2 = uint32toa(ujj - uii, wptr2);
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
	    wptr = dtoa_g_wxp4x(chisq, 12, ' ', wptr);
	    dxx = chiprob_p(chisq, 1);
	    wptr = dtoa_g_wxp4(MAXV(dxx, output_min_p), 12, wptr);
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
	    wptr = dtoa_g_wxp4x(chisq, 12, ' ', wptr);
	    dxx = chiprob_p(chisq, 1);
	    wptr = dtoa_g_wxp4(MAXV(dxx, output_min_p), 12, wptr);
	  }
	}
	wptr = memcpya(wptr, " \n", 2);
        if (fwrite_checked(textbuf, wptr - textbuf, outfile)) {
	  goto tdt_ret_WRITE_FAIL;
	}
      }

      if (++markers_done >= pct_thresh) {
        if (pct > 10) {
	  putc_unlocked('\b', stdout);
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
  putc_unlocked('\r', stdout);
  LOGPRINTF("--tdt: Report written to %s .\n", outname);
  if (mtest_adjust) {
  tdt_multcomp:
    ulii = BITCT_TO_WORDCT(unfiltered_marker_ct);
    if (bigstack_alloc_ui(marker_ct, &marker_idx_to_uidx) ||
        bigstack_alloc_ul(ulii, &marker_exclude_tmp)) {
      goto tdt_ret_NOMEM;
    }
    // need a custom marker_exclude that's set at Y/haploid/MT
    memcpy(marker_exclude_tmp, marker_exclude, ulii * sizeof(intptr_t));
    for (chrom_fo_idx = 0; chrom_fo_idx < chrom_info_ptr->chrom_ct; chrom_fo_idx++) {
      chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
      if ((is_set(chrom_info_ptr->haploid_mask, chrom_idx) && ((int32_t)chrom_idx != chrom_info_ptr->xymt_codes[X_OFFSET])) || ((int32_t)chrom_idx == chrom_info_ptr->xymt_codes[MT_OFFSET])) {
	uii = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx];
	fill_bits(uii, chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx + 1] - uii, marker_exclude_tmp);
      }
    }
    fill_idx_to_uidx(marker_exclude_tmp, unfiltered_marker_ct, marker_ct, marker_idx_to_uidx);
    retval = multcomp(outname, outname_end, marker_idx_to_uidx, marker_ct, marker_ids, max_marker_id_len, plink_maxsnp, chrom_info_ptr, is_exact? nullptr : orig_chisq, pfilter, output_min_p, mtest_adjust, 0, adjust_lambda, nullptr, is_exact? orig_chisq : nullptr);
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
  bigstack_reset(bigstack_mark);
  fclose_cond(outfile);
  return retval;
}

int32_t get_sibship_info(uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, uintptr_t* pheno_nm, uintptr_t* founder_info, char* sample_ids, uintptr_t max_sample_id_len, uintptr_t max_fid_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uint64_t* family_list, uint64_t* trio_list, uint32_t family_ct, uintptr_t trio_ct, uint32_t test_type, uintptr_t** size_one_sibships_ptr, uintptr_t** lm_eligible_ptr, uintptr_t** lm_within2_founder_ptr, uint32_t** fs_starts_ptr, uint32_t** fss_contents_ptr, uint32_t** sample_lm_to_fss_idx_ptr, uint32_t* fs_ct_ptr, uint32_t* lm_ct_ptr, uint32_t* singleton_ct_ptr) {
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

  // It's probably appropriate to split this into two functions in the future,
  // one for dfam and one for qfam; the differences make this difficult to
  // maintain.
  unsigned char* bigstack_end_mark = g_bigstack_end;
  uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t sample_ctl = BITCT_TO_WORDCT(sample_ct);
  uintptr_t max_merged_id_len = max_fid_len + max_paternal_id_len + max_maternal_id_len + sizeof(int32_t);
  uintptr_t trio_idx = 0;
  uintptr_t* tmp_within2_founder = nullptr;
  uintptr_t* lm_within2_founder = nullptr;
  uintptr_t* lm_eligible = nullptr;
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
  unsigned char* bigstack_end_mark2;
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
      if (bigstack_alloc_ul(sample_ctl, &lm_within2_founder)) {
	goto get_sibship_info_ret_NOMEM;
      }
    }
    if (bigstack_alloc_ul(sample_ctl, &lm_eligible)) {
      goto get_sibship_info_ret_NOMEM;
    }
  }
  if (test_type) {
    // shrink later
    if (bigstack_alloc_ui(sample_ct + 2 * family_ct, &fss_contents)) {
      goto get_sibship_info_ret_NOMEM;
    }
    // this is the equivalent of PLINK 1.07's family pointers
    if (bigstack_end_alloc_ui(sample_ct, &sample_to_fss_idx)) {
      goto get_sibship_info_ret_NOMEM;
    }
  } else {
    if (bigstack_alloc_ui(sample_ct, &sample_to_fss_idx)) {
      goto get_sibship_info_ret_NOMEM;
    }
    // shrink later
    if (bigstack_alloc_ui(sample_ct + 2 * family_ct, &fss_contents)) {
      goto get_sibship_info_ret_NOMEM;
    }
  }
  bigstack_end_mark2 = g_bigstack_end;

  if (bigstack_end_alloc_ul(unfiltered_sample_ctl, &not_in_family) ||

      // Temporary bitfields used to track which parents are (i) part of
      // multiple families, and (ii) not a child in any of them.  To ensure
      // results are not dependent on the order of samples in the dataset, we
      // now exclude these parents from the QFAM permutation.  (todo: compute
      // an average in this case instead?)
      // ulptr = is a double-parent
      // ulptr2 = is a child
      bigstack_end_alloc_ui(unfiltered_sample_ct, &sample_uidx_to_idx) ||
      bigstack_end_calloc_ul(unfiltered_sample_ctl, &ulptr) ||
      bigstack_end_calloc_ul(unfiltered_sample_ctl, &ulptr2)) {
    goto get_sibship_info_ret_NOMEM;
  }
  if (is_within2) {
    if (bigstack_end_calloc_ul(unfiltered_sample_ctl, &tmp_within2_founder)) {
      goto get_sibship_info_ret_NOMEM;
    }
  }

  bitarr_invert_copy(sample_exclude, unfiltered_sample_ct, not_in_family);
  fill_uint_one(sample_ct, sample_to_fss_idx);
  fill_uidx_to_idx(sample_exclude, unfiltered_sample_ct, sample_ct, sample_uidx_to_idx);
  if (family_ct) {
    // iterate over all parents
    while (1) {
      ullii = family_list[family_idx];
      // uii, ukk = unfiltered idxs of father, mother respectively
      // ujj, umm = filtered idxs
      uii = (uint32_t)ullii;
      ujj = sample_uidx_to_idx[uii];
      fss_contents[fssc_idx++] = ujj;
      ukk = (uint32_t)(ullii >> 32);
      umm = sample_uidx_to_idx[ukk];
      if (is_within2) {
	if (is_set(pheno_nm, uii) && is_set(pheno_nm, ukk)) {
	  set_bit(uii, tmp_within2_founder);
	  set_bit(ukk, tmp_within2_founder);
	}
      }
      if (is_set(not_in_family, uii)) {
	if (sample_to_fss_idx[ujj] == 0xffffffffU) {
	  // missing father
	  sample_to_fss_idx[ujj] = family_idx;
	}
	clear_bit(uii, not_in_family);
      } else {
	set_bit(uii, ulptr);
      }
      fss_contents[fssc_idx++] = umm;
      if (is_set(not_in_family, ukk)) {
	if (sample_to_fss_idx[umm] == 0xffffffffU) {
	  // missing mother
	  sample_to_fss_idx[umm] = family_idx;
	}
	clear_bit(ukk, not_in_family);
      } else {
	set_bit(ukk, ulptr);
      }

      ullii = trio_list[trio_idx];
      // do-while since there's at least one trio for each set of parents
      do {
	uii = (uint32_t)ullii;
	ujj = sample_uidx_to_idx[uii];
	fss_contents[fssc_idx++] = ujj;
	sample_to_fss_idx[ujj] = family_idx;
	set_bit(uii, ulptr2);
	if (++trio_idx == trio_ct) {
	  goto get_sibship_info_first_pass_done;
	}
	ullii = trio_list[trio_idx];
      } while ((uint32_t)(ullii >> 32) == family_idx);
      family_idx++;
    }
  }
 get_sibship_info_first_pass_done:
  bitvec_andnot(ulptr2, unfiltered_sample_ctl, not_in_family);
  bigstack_shrink_top(fss_contents, (fssc_idx + popcount_longs(not_in_family, unfiltered_sample_ctl)) * sizeof(int32_t));
  if (test_type) {
    bitvec_andnot(ulptr2, unfiltered_sample_ctl, ulptr);
  } else {
    bitarr_invert_copy(ulptr2, unfiltered_sample_ct, ulptr);
  }
  // qfam: ulptr = double-parents who aren't also a child of two parents in
  //               immediate dataset
  // dfam: ulptr = everyone who isn't a child of two parents in immediate
  //               dataset

  if (is_within2) {
    bitvec_andnot(ulptr, unfiltered_sample_ctl, tmp_within2_founder);
    bitvec_and(founder_info, unfiltered_sample_ctl, tmp_within2_founder);
    // now this only consists of founder parents who (i) aren't in multiple
    // families, and (ii) have a different phenotype from their partner.
    copy_bitarr_subset_excl(tmp_within2_founder, sample_exclude, unfiltered_sample_ct, sample_ct, lm_within2_founder);
  }
  if (test_type) {
    bitvec_andnot_reversed_args(pheno_nm, unfiltered_sample_ctl, ulptr);
    if (test_type == QFAM_WITHIN1) {
      bitvec_andnot(founder_info, unfiltered_sample_ctl, ulptr);
    }
    copy_bitarr_subset_excl(ulptr, sample_exclude, unfiltered_sample_ct, sample_ct, lm_eligible);
    bitvec_andnot_copy(not_in_family, founder_info, unfiltered_sample_ctl, ulptr);
  } else {
    bitvec_and(pheno_nm, unfiltered_sample_ctl, ulptr);
    bitvec_andnot(founder_info, unfiltered_sample_ctl, ulptr);
  }
  bigstack_end_reset(ulptr);

  // qfam: not a parent or child in a trio, not a founder
  // dfam: not a child in a trio, not a founder; parent ok

  cur_sample_ct = popcount_longs(ulptr, unfiltered_sample_ctl);

  if (bigstack_alloc_ui(1 + family_ct + (cur_sample_ct / 2), &fs_starts)) {
    goto get_sibship_info_ret_NOMEM;
  }
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
    // identify size-2+ sibships
    if (bigstack_end_alloc_c(max_merged_id_len * cur_sample_ct, &merged_ids)) {
      goto get_sibship_info_ret_NOMEM;
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
	clear_bit(uii, not_in_family);
	ujj = sample_uidx_to_idx[uii];
        fss_contents[fssc_idx++] = ujj;
        sample_to_fss_idx[ujj] = family_idx;
	do {
	  uii = *((uint32_t*)(&(bufptr2[slen])));
	  clear_bit(uii, not_in_family);
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
  bigstack_shrink_top(fs_starts, (family_idx + 1) * sizeof(int32_t));
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
    bigstack_end_reset(bigstack_end_mark2);
    ulii = popcount_longs(lm_eligible, sample_ctl);
    if (bigstack_alloc_ui(ulii, &sample_lm_to_fss_idx)) {
      goto get_sibship_info_ret_NOMEM;
    }
    for (sample_uidx = 0, sample_idx = 0; sample_idx < ulii; sample_uidx++, sample_idx++) {
      next_set_ul_unsafe_ck(lm_eligible, &sample_uidx);
      sample_lm_to_fss_idx[sample_idx] = sample_to_fss_idx[sample_uidx];
    }
    *lm_eligible_ptr = lm_eligible;
    *lm_within2_founder_ptr = lm_within2_founder;
    *sample_lm_to_fss_idx_ptr = sample_lm_to_fss_idx;
    *lm_ct_ptr = ulii;
  } else {
    // bugfix: for DFAM, we also need to prevent size-1 sibships from being
    // included in the unrelated cluster
    if (bigstack_alloc_ul(unfiltered_sample_ctl, size_one_sibships_ptr)) {
      goto get_sibship_info_ret_NOMEM;
    }
    memcpy(*size_one_sibships_ptr, not_in_family, unfiltered_sample_ctl * sizeof(intptr_t));
    bitvec_and(ulptr, unfiltered_sample_ctl, *size_one_sibships_ptr);

    // return sample_to_fss_idx in place of sample_lm_to_fss_idx
    *sample_lm_to_fss_idx_ptr = sample_to_fss_idx;
  }
  *fs_starts_ptr = fs_starts;
  *fss_contents_ptr = fss_contents;

  while (0) {
  get_sibship_info_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  }
  bigstack_end_reset(bigstack_end_mark);
  return retval;
}

// multithread globals
static double* g_maxt_extreme_stat;
static double* g_maxt_thread_results;
static double* g_mperm_save_all;
static uintptr_t* g_pheno_c;
static uintptr_t* g_dfam_flipa;
#ifdef __LP64__
static uintptr_t* g_dfam_flipa_shuffled;
#endif
static uintptr_t* g_dfam_perm_vecst; // sample-major, shuffled
static double* g_dfam_numers;
static double* g_dfam_denoms;
static uintptr_t* g_dfam_acc;
static int32_t* g_dfam_twice_numers;
static uint32_t* g_dfam_total_counts;
static uint32_t* g_dfam_iteration_order;
static uint32_t g_dfam_family_all_case_children_ct;
static uint32_t g_dfam_family_mixed_ct;
static uint32_t g_dfam_sibship_mixed_ct;
static uint32_t g_dfam_unrelated_cluster_ct;

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
static uintptr_t g_qfam_sample_ct;
static double g_qt_sum_all;
static double g_qt_ssq_all;
static uint32_t g_test_type;
static uint32_t g_xfam_thread_ct;
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

// note that 4 is encoded as 0 since they're treated the same.
// tried encoding this in a single 32-bit integer, but that appears to be
// slower.
const uint8_t dfam_allele_ct_table[] =
{0, 0, 3, 0,
 0, 0, 0, 0,
 3, 0, 2, 1,
 0, 0, 1, 0};

void dfam_sibship_or_unrelated_perm_calc(uintptr_t* loadbuf_ptr, const uint32_t* cur_dfam_ptr, const uintptr_t* perm_vecst, const uintptr_t* orig_pheno_c, uint32_t sibling_ct, uint32_t is_unrelated_calc, uintptr_t perm_vec_ct,
#ifdef __LP64__
					 __m128i* acc4, __m128i* acc8,
#else
					 uintptr_t* acc4, uintptr_t* acc8,
#endif
					 uint32_t* cur_case_a1_cts, uint32_t* cur_case_missing_cts, int32_t* twice_numers, double* numers, double* denoms, uint32_t* total_counts) {
  // okay, compute array of familial/sibship case_a1_ct values.  Most
  // families/sibships should have 7 or fewer children, so it makes sense
  // to use 4-bit accumulators in the inner loop (similar to calc_git()
  // in plink_assoc.c).
  uintptr_t perm_vec_ct128 = (perm_vec_ct + 127) / 128;
  const uintptr_t perm_vec_wcta = perm_vec_ct128 * (128 / BITCT);
  uint32_t cur_genotype_cts[4];
#ifdef __LP64__
  const __m128i m1x4 = {0x1111111111111111LLU, 0x1111111111111111LLU};
  const __m128i m1x4ls1 = {0x2222222222222222LLU, 0x2222222222222222LLU};
  uintptr_t acc4_word_ct = perm_vec_ct128 * 8;
  // uintptr_t acc8_word_ct = perm_vec_ct128 * 16;
  uintptr_t acc4_vec_ct = perm_vec_ct128 * 4;
  uintptr_t acc8_vec_ct = acc4_word_ct;
  const __m128i* pheno_perm_ptr;
  __m128i* acc4_ptr;
  __m128i loader;
#else
  uintptr_t acc4_word_ct = perm_vec_ct128 * 16;
  uintptr_t acc8_word_ct = perm_vec_ct128 * 32;
  uintptr_t perm_vec_wct = BITCT_TO_WORDCT(perm_vec_ct);
  const uintptr_t* pheno_perm_ptr;
  uintptr_t* acc4_ptr;
  uintptr_t loader;
#endif
  uint32_t case_ct_base = 0;
  uintptr_t perm_idx;
  double total_ctd;
  double total_ct_recip;
  double xxm1_recip;
  double hom_a1_ctd;
  double het_ctd;
  double case_ctd;
  double ctrl_ctd;
  double case_proportion;
  double case_expected_hom_a1;
  double case_expected_het;
  double case_ctrl_div_xxxm1;
  double case_var_hom_a1;
  double case_var_het;
  double case_neg_covar;
  double case_expected_a1_ct;
  double case_var_a1_ct;
  double case_a1_ctd;
  double dbl_total_ctd;
  uint32_t sib_idx;
  uint32_t sample_idx;
  uint32_t cur_geno;
  uint32_t geno_match;
  uint32_t cur_case_ct;
  uint32_t case_missing_ct;
  uint32_t case_a1_ct;
  uint32_t total_ct;
  uint32_t cur_ctrl_ct;
  uint32_t max_incr4;
  uint32_t max_incr8;
  uint32_t uii;
  // first check if all genotypes are identical
  fill_uint_zero(4, cur_genotype_cts);
  for (sib_idx = 0; sib_idx < sibling_ct; sib_idx++) {
    sample_idx = cur_dfam_ptr[sib_idx];
    cur_geno = EXTRACT_2BIT_GENO(loadbuf_ptr, sample_idx);
    cur_genotype_cts[cur_geno] += 1;
    case_ct_base += IS_SET(orig_pheno_c, sample_idx);
  }
  cur_geno = 4;
  for (geno_match = 0; geno_match < 4; geno_match++) {
    if (cur_genotype_cts[geno_match]) {
      if (cur_geno != 4) {
	break;
      }
      cur_geno = geno_match;
    }
  }
  if (geno_match == 4) {
    if ((!is_unrelated_calc) && (!(cur_geno % 2))) {
      if (!cur_geno) {
	uii = cur_genotype_cts[0] * 2;
      } else {
	uii = cur_genotype_cts[0];
      }
      for (perm_idx = 0; perm_idx < perm_vec_ct; perm_idx++) {
	total_counts[perm_idx] += uii;
      }
    }
    return;
  }

#ifdef __LP64__
  fill_vvec_zero(acc4_vec_ct, acc4);
  fill_vvec_zero(acc8_vec_ct, acc8);
#else
  fill_ulong_zero(acc4_word_ct, acc4);
  fill_ulong_zero(acc8_word_ct, acc8);
#endif
  fill_uint_zero(perm_vec_ct, cur_case_a1_cts);
  max_incr4 = 0;
  max_incr8 = 0;
  for (sib_idx = 0; sib_idx < sibling_ct; sib_idx++) {
    sample_idx = cur_dfam_ptr[sib_idx];
    cur_geno = EXTRACT_2BIT_GENO(loadbuf_ptr, sample_idx);
    if (cur_geno & 1) {
      continue;
    }
    uii = 2 - (cur_geno / 2);
#ifdef __LP64__
    if (max_incr4 + uii > 15) {
      unroll_zero_incr_4_8(acc4, acc8, acc4_vec_ct);
      max_incr8 += max_incr4;
      if (max_incr8 > 240) {
	unroll_zero_incr_8_32(acc8, (__m128i*)cur_case_a1_cts, acc8_vec_ct);
	max_incr8 = 0;
      }
      max_incr4 = 0;
    }
    max_incr4 += uii;
    pheno_perm_ptr = (const __m128i*)(&(perm_vecst[sample_idx * perm_vec_wcta]));
    if (cur_geno) {
      unroll_incr_1_4(pheno_perm_ptr, acc4, perm_vec_ct128);
    } else {
      // add 2 whenever this sample x permutation is a case
      acc4_ptr = acc4;
      for (uii = 0; uii < acc4_vec_ct; uii++) {
	loader = *pheno_perm_ptr++;
	acc4_ptr[0] = _mm_add_epi64(acc4_ptr[0], _mm_slli_epi64(_mm_and_si128(loader, m1x4), 1));
	acc4_ptr[1] = _mm_add_epi64(acc4_ptr[1], _mm_and_si128(loader, m1x4ls1));
	loader = _mm_srli_epi64(loader, 1);
	acc4_ptr[2] = _mm_add_epi64(acc4_ptr[2], _mm_and_si128(loader, m1x4ls1));
	loader = _mm_srli_epi64(loader, 1);
	acc4_ptr[3] = _mm_add_epi64(acc4_ptr[3], _mm_and_si128(loader, m1x4ls1));
	acc4_ptr = &(acc4_ptr[4]);
      }
    }
#else
    if (max_incr4 + uii > 15) {
      unroll_zero_incr_4_8(acc4, acc8, acc4_word_ct);
      max_incr8 += max_incr4;
      if (max_incr8 > 240) {
	unroll_zero_incr_8_32(acc8, (uintptr_t*)cur_case_a1_cts, acc8_word_ct);
	max_incr8 = 0;
      }
      max_incr4 = 0;
    }
    max_incr4 += uii;
    pheno_perm_ptr = &(perm_vecst[sample_idx * perm_vec_wcta]);
    if (cur_geno) {
      unroll_incr_1_4(pheno_perm_ptr, acc4, perm_vec_wct);
    } else {
      acc4_ptr = acc4;
      for (uii = 0; uii < perm_vec_wct; uii++) {
	loader = *pheno_perm_ptr++;
	acc4_ptr[0] += (loader & 0x11111111U) << 1;
	acc4_ptr[1] += loader & 0x22222222U;
	acc4_ptr[2] += (loader >> 1) & 0x22222222U;
	acc4_ptr[3] += (loader >> 2) & 0x22222222U;
	acc4_ptr = &(acc4_ptr[4]);
      }
    }
#endif
  }
#ifdef __LP64__
  unroll_incr_4_8(acc4, acc8, acc4_vec_ct);
  unroll_incr_8_32(acc8, (__m128i*)cur_case_a1_cts, acc8_vec_ct);
#else
  unroll_incr_4_8(acc4, acc8, acc4_word_ct);
  unroll_incr_8_32(acc8, (uintptr_t*)cur_case_a1_cts, acc8_word_ct);
#endif

  if (!cur_genotype_cts[1]) {
    // optimize the common no-missing-genotypes case
    total_ctd = (double)((int32_t)sibling_ct);
    total_ct_recip = 1.0 / total_ctd;
    case_ctd = (double)((int32_t)case_ct_base);
    case_proportion = case_ctd * total_ct_recip;
    cur_ctrl_ct = sibling_ct - case_ct_base;
    ctrl_ctd = (double)((int32_t)cur_ctrl_ct);
    if (!is_unrelated_calc) {
      // actually ctrl_ct/(x(x-1)), not 1/(x(x-1))
      xxm1_recip = ctrl_ctd * total_ct_recip / ((double)((int32_t)(sibling_ct - 1)));
      hom_a1_ctd = (double)((int32_t)cur_genotype_cts[0]);
      het_ctd = (double)((int32_t)cur_genotype_cts[2]);
      case_expected_hom_a1 = case_proportion * hom_a1_ctd;
      case_expected_het = case_proportion * het_ctd;
      case_ctrl_div_xxxm1 = case_proportion * xxm1_recip;
      case_var_hom_a1 = case_ctrl_div_xxxm1 * hom_a1_ctd * (total_ctd - hom_a1_ctd);
      case_var_het = case_ctrl_div_xxxm1 * het_ctd * (total_ctd - het_ctd);
      case_neg_covar = case_ctrl_div_xxxm1 * het_ctd;
      case_expected_a1_ct = 2 * case_expected_hom_a1 + case_expected_het;
      case_var_a1_ct = 4 * (case_var_hom_a1 + case_neg_covar) + case_var_het;
      for (perm_idx = 0; perm_idx < perm_vec_ct; perm_idx++) {
	case_a1_ct = cur_case_a1_cts[perm_idx];
	total_counts[perm_idx] += case_a1_ct;
	numers[perm_idx] += (double)((int32_t)case_a1_ct) - case_expected_a1_ct;
	denoms[perm_idx] += case_var_a1_ct;
      }
    } else {
      // actually ctrl_ct/(x(2x-1)), not 1/(x(x-1))
      xxm1_recip = ctrl_ctd * total_ct_recip / ((double)((int32_t)(2 * sibling_ct - 1)));
      dbl_total_ctd = 2 * total_ctd;
      for (perm_idx = 0; perm_idx < perm_vec_ct; perm_idx++) {
	case_a1_ct = cur_case_a1_cts[perm_idx];
	case_a1_ctd = (double)((int32_t)case_a1_ct);
	case_expected_a1_ct = case_proportion * case_a1_ctd;
	case_var_a1_ct = case_expected_a1_ct * (dbl_total_ctd - case_a1_ctd) * xxm1_recip;
	total_counts[perm_idx] += case_a1_ct;
	numers[perm_idx] += case_a1_ctd - case_expected_a1_ct;
	denoms[perm_idx] += case_var_a1_ct;
      }
    }
    return;
  }

#ifdef __LP64__
  fill_vvec_zero(acc4_vec_ct, acc4);
  fill_vvec_zero(acc8_vec_ct, acc8);
#else
  fill_ulong_zero(acc4_word_ct, acc4);
  fill_ulong_zero(acc8_word_ct, acc8);
#endif
  fill_uint_zero(perm_vec_ct, cur_case_missing_cts);
  uii = 0;
  for (sib_idx = 0; sib_idx < sibling_ct; sib_idx++) {
    sample_idx = cur_dfam_ptr[sib_idx];
    cur_geno = EXTRACT_2BIT_GENO(loadbuf_ptr, sample_idx);
    if (cur_geno != geno_match) {
      continue;
    }

#ifdef __LP64__
    pheno_perm_ptr = (const __m128i*)(&(perm_vecst[sample_idx * perm_vec_wcta]));
    unroll_incr_1_4(pheno_perm_ptr, acc4, perm_vec_ct128);
    if (!(uii % 15)) {
      unroll_zero_incr_4_8(acc4, acc8, acc4_vec_ct);
      if (!(uii % 255)) {
	unroll_zero_incr_8_32(acc8, (__m128i*)cur_case_missing_cts, acc8_vec_ct);
      }
    }
#else
    pheno_perm_ptr = &(perm_vecst[sample_idx * perm_vec_wcta]);
    unroll_incr_1_4(pheno_perm_ptr, acc4, perm_vec_wct);
    if (!(uii % 15)) {
      unroll_zero_incr_4_8(acc4, acc8, acc4_word_ct);
      if (!(uii % 255)) {
	unroll_zero_incr_8_32(acc8, (uintptr_t*)cur_case_missing_cts, acc8_word_ct);
      }
    }
#endif
  }
  if (uii % 255) {
#ifdef __LP64__
    if (uii % 15) {
      unroll_incr_4_8(acc4, acc8, acc4_vec_ct);
    }
    unroll_incr_8_32(acc8, (__m128i*)cur_case_missing_cts, acc8_vec_ct);
#else
    if (uii % 15) {
      unroll_incr_4_8(acc4, acc8, acc4_word_ct);
    }
    unroll_incr_8_32(acc8, (uintptr_t*)cur_case_missing_cts, acc8_word_ct);
#endif
  }

  total_ct = sibling_ct - cur_genotype_cts[1];
  total_ctd = (double)((int32_t)total_ct);
  total_ct_recip = 1.0 / total_ctd;
  if (!is_unrelated_calc) {
    xxm1_recip = total_ct_recip / ((double)((int32_t)(total_ct - 1)));
    hom_a1_ctd = (double)((int32_t)cur_genotype_cts[0]);
    het_ctd = (double)((int32_t)cur_genotype_cts[2]);
    for (perm_idx = 0; perm_idx < perm_vec_ct; perm_idx++) {
      case_missing_ct = cur_case_missing_cts[perm_idx];
      cur_case_ct = case_ct_base - case_missing_ct;
      cur_ctrl_ct = total_ct - cur_case_ct;
      if ((!cur_case_ct) || (!cur_ctrl_ct)) {
	continue;
      }
      case_ctd = (double)((int32_t)cur_case_ct);
      ctrl_ctd = (double)((int32_t)cur_ctrl_ct);
      case_a1_ct = cur_case_a1_cts[perm_idx];
      case_proportion = case_ctd * total_ct_recip;
      case_expected_hom_a1 = case_proportion * hom_a1_ctd;
      case_expected_het = case_proportion * het_ctd;
      case_ctrl_div_xxxm1 = case_proportion * ctrl_ctd * xxm1_recip;
      case_var_hom_a1 = case_ctrl_div_xxxm1 * hom_a1_ctd * (total_ctd - hom_a1_ctd);
      case_var_het = case_ctrl_div_xxxm1 * het_ctd * (total_ctd - het_ctd);
      case_neg_covar = case_ctrl_div_xxxm1 * het_ctd;
      case_expected_a1_ct = 2 * case_expected_hom_a1 + case_expected_het;
      case_var_a1_ct = 4 * (case_var_hom_a1 + case_neg_covar) + case_var_het;
      total_counts[perm_idx] += case_a1_ct;
      numers[perm_idx] += (double)((int32_t)case_a1_ct) - case_expected_a1_ct;
      denoms[perm_idx] += case_var_a1_ct;
    }
  } else {
    // actually 1/(x(2x-1)), not 1/(x(x-1))
    xxm1_recip = total_ct_recip / ((double)((int32_t)(2 * total_ct - 1)));
    dbl_total_ctd = 2 * total_ctd;

    for (perm_idx = 0; perm_idx < perm_vec_ct; perm_idx++) {
      case_missing_ct = cur_case_missing_cts[perm_idx];
      cur_case_ct = case_ct_base - case_missing_ct;
      cur_ctrl_ct = total_ct - cur_case_ct;
      if ((!cur_case_ct) || (!cur_ctrl_ct)) {
	continue;
      }
      case_proportion = ((double)((int32_t)cur_case_ct)) * total_ct_recip;
      case_a1_ct = cur_case_a1_cts[perm_idx];
      case_a1_ctd = (double)((int32_t)case_a1_ct);
      case_expected_a1_ct = case_proportion * case_a1_ctd;
      case_var_a1_ct = case_expected_a1_ct * (dbl_total_ctd - case_a1_ctd) * ((double)((int32_t)cur_ctrl_ct)) * xxm1_recip;
      total_counts[perm_idx] += case_a1_ct;
      numers[perm_idx] += case_a1_ctd - case_expected_a1_ct;
      denoms[perm_idx] += case_var_a1_ct;
    }
  }
}

THREAD_RET_TYPE dfam_perm_thread(void* arg) {
  uintptr_t tidx = (uintptr_t)arg;
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uintptr_t perm_vec_ct128 = (perm_vec_ct + 127) / 128;
  uintptr_t perm_vec_cta128 = perm_vec_ct128 * 128;
  uintptr_t perm_vec_ctcl8m = round_up_pow2(perm_vec_ct, CACHELINE_DBL);
  uint32_t dfam_thread_ct = g_xfam_thread_ct;
  uint32_t pidx_offset = g_perms_done;
  uint32_t first_adapt_check = g_first_adapt_check;
  uint32_t family_all_case_children_ct = g_dfam_family_all_case_children_ct;
  uint32_t family_mixed_ct = g_dfam_family_mixed_ct;
  uint32_t sibship_mixed_ct = g_dfam_sibship_mixed_ct;
  uint32_t unrelated_cluster_ct = g_dfam_unrelated_cluster_ct;
  uint32_t dfam_sample_ct = g_perm_pheno_nm_ct;
  uint32_t dfam_sample_ctl2 = QUATERCT_TO_WORDCT(dfam_sample_ct);
  const uintptr_t perm_vec_wcta = perm_vec_ct128 * (128 / BITCT);
  const uintptr_t* flipa = g_dfam_flipa;
  const uintptr_t* perm_vecst = g_dfam_perm_vecst;
  const uintptr_t* orig_pheno_c = g_pheno_c;
  int32_t* __restrict__ twice_numers = &(g_dfam_twice_numers[tidx * perm_vec_cta128]);
  uint32_t* __restrict__ total_counts = &(g_dfam_total_counts[tidx * perm_vec_cta128]);
  uint32_t* __restrict__ perm_attempt_ct = g_perm_attempt_ct;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  double* __restrict__ mperm_save_all = g_mperm_save_all;
  double* msa_ptr = nullptr;
  double* numers = &(g_dfam_numers[tidx * perm_vec_cta128]);
  double* denoms = &(g_dfam_denoms[tidx * perm_vec_cta128]);
  const uint32_t* dfam_iteration_order = g_dfam_iteration_order;
  unsigned char* perm_adapt_stop = nullptr;
  double adaptive_intercept = 0.0;
  double adaptive_slope = 0.0;
  double adaptive_ci_zt = 0.0;
  double aperm_alpha = 0.0;
  double* maxt_results = nullptr;
  uint32_t perm_adapt = g_test_type;
  uint32_t next_adapt_check = 0;
  uint32_t cur_case_a1_ct_flip[2];
#ifdef __LP64__
  const __m128i m1x8 = {0x0101010101010101LLU, 0x0101010101010101LLU};
  const __m128i m1x4 = {0x1111111111111111LLU, 0x1111111111111111LLU};
  const __m128i m1x4ls1 = {0x2222222222222222LLU, 0x2222222222222222LLU};
  __m128i diff_vec;
  __m128i incr8;
  __m128i loader;
  // acc8 (8-bit accumulator) requires (perm_vec_ct + 7) / 8 words; this is
  //   16-byte aligned when perm_vec_ct is divisible by 16
  // acc4 requires (perm_vec_ct + 15) / 16 words
  // sum reduces to (perm_vec_ct128 * 248) since we have 3 acc8s and 2 acc32s
  const uintptr_t acc_thread_offset = perm_vec_ct128 * 184;
  const uintptr_t acc4_word_ct = perm_vec_ct128 * 8;
  const uintptr_t acc8_word_ct = perm_vec_ct128 * 16;
  const uintptr_t acc4_vec_ct = perm_vec_ct128 * 4;
  const uintptr_t acc8_vec_ct = acc4_word_ct;
  __m128i* acc4 = (__m128i*)(&(g_dfam_acc[tidx * acc_thread_offset]));
  __m128i* acc8 = (__m128i*)(&(g_dfam_acc[tidx * acc_thread_offset + acc4_word_ct]));
  __m128i* case_a1_ct_acc8 = (__m128i*)(&(g_dfam_acc[tidx * acc_thread_offset + acc4_word_ct + acc8_word_ct]));
  // __m128i* cur_case_ct_acc8 = (__m128i*)(&(g_dfam_acc[tidx * acc_thread_offset + acc4_word_ct + 2 * acc8_word_ct]));

  uint32_t* cur_case_a1_cts = (uint32_t*)(&(g_dfam_acc[tidx * acc_thread_offset + acc4_word_ct + 3 * acc8_word_ct]));
  uint32_t* cur_case_missing_cts = (uint32_t*)(&(g_dfam_acc[tidx * acc_thread_offset + acc4_word_ct + 7 * acc8_word_ct]));

  const uintptr_t* flipa_shuffled = g_dfam_flipa_shuffled;
  const __m128i* pheno_perm_ptr;
  const __m128i* flipa_perm_ptr;
  __m128i* acc4_ptr;
  __m128i* acc8_ptr;
  uintptr_t vidx;
#else
  const uintptr_t perm_vec_wct = BITCT_TO_WORDCT(perm_vec_ct);
  // acc8 requires (perm_vec_ct + 3) / 4 words
  // acc4 requires (perm_vec_ct + 7) / 8 words
  // sum reduces to perm_vec_ct128 * 304 since we also have 2 acc32s
  const uintptr_t acc_thread_offset = perm_vec_ct128 * 304;
  const uintptr_t acc4_word_ct = perm_vec_ct128 * 16;
  const uintptr_t acc8_word_ct = perm_vec_ct128 * 32;
  uintptr_t* acc4 = &(g_dfam_acc[tidx * acc_thread_offset]);
  uintptr_t* acc8 = &(g_dfam_acc[tidx * acc_thread_offset + acc4_word_ct]);
  uint32_t* cur_case_a1_cts = (uint32_t*)(&(g_dfam_acc[tidx * acc_thread_offset + acc4_word_ct + acc8_word_ct]));
  uint32_t* cur_case_missing_cts = (uint32_t*)(&(g_dfam_acc[tidx * acc_thread_offset + acc4_word_ct + 5 * acc8_word_ct]));
  const uintptr_t* pheno_perm_ptr;
  uintptr_t* acc4_ptr;
  uintptr_t loader;
  uintptr_t widx;
#endif
  uintptr_t perm_idx;
  const uintptr_t* cur_flipa;
  double* orig_chisq;
  const uint32_t* cur_dfam_ptr;
  uintptr_t* loadbuf_ptr;
  double chisq_high;
  double chisq_low;
  double chisq;
  double pval;
  double dxx;
  double dyy;
  double dzz;
  uint32_t marker_bidx;
  uint32_t marker_bceil;
  uint32_t marker_idx;
  uint32_t sample_idx;
  uint32_t fs_idx;
  uint32_t unrelated_cluster_idx;
  uint32_t quad_denom;
  uint32_t twice_numer_subtract;
  uint32_t paternal_id;
  uint32_t maternal_id;
  uint32_t sibling_ct;
  uint32_t paternal_geno;
  uint32_t maternal_geno;
  uint32_t parental_a1_ct;
  uint32_t sib_idx;
  uint32_t nonmissing_sib_ct;
  uint32_t cur_geno;
  uint32_t is_flipped;
  uint32_t max_incr4;
  uint32_t max_incr8;
  uint32_t cur_max_incr;
  uint32_t orig_case_ct;
  uint32_t success_2start;
  uint32_t success_2incr;
  uint32_t uii;
  uint32_t ujj;
  if (perm_adapt) {
    perm_adapt_stop = g_perm_adapt_stop;
    adaptive_intercept = g_adaptive_intercept;
    adaptive_slope = g_adaptive_slope;
    adaptive_ci_zt = g_adaptive_ci_zt;
    aperm_alpha = g_aperm_alpha;
  } else {
    maxt_results = &(g_maxt_thread_results[perm_vec_ctcl8m * tidx]);
  }
  while (1) {
    if (g_block_size <= dfam_thread_ct) {
      if (g_block_size <= tidx) {
	goto dfam_perm_thread_skip_all;
      }
      marker_bidx = tidx;
      marker_bceil = tidx + 1;
    } else {
      marker_bidx = (((uint64_t)tidx) * g_block_size) / dfam_thread_ct;
      marker_bceil = (((uint64_t)tidx + 1) * g_block_size) / dfam_thread_ct;
    }
    orig_chisq = g_orig_stat;
    for (; marker_bidx < marker_bceil; marker_bidx++) {
      marker_idx = g_adapt_m_table[marker_bidx];
      loadbuf_ptr = &(g_loadbuf[marker_bidx * dfam_sample_ctl2]);
      if (perm_adapt) {
	next_adapt_check = first_adapt_check;
      } else if (mperm_save_all) {
	msa_ptr = &(mperm_save_all[marker_idx * perm_vec_ct]);
      }
      quad_denom = 0;
      twice_numer_subtract = 0;
      fill_uint_zero(perm_vec_ct, total_counts);
      fill_double_zero(perm_vec_ct, numers);
      fill_double_zero(perm_vec_ct, denoms);

      cur_dfam_ptr = dfam_iteration_order;
      success_2start = perm_2success_ct[marker_idx];
      success_2incr = 0;
      chisq_high = orig_chisq[marker_idx] + EPSILON;
      chisq_low = orig_chisq[marker_idx] - EPSILON;
#ifdef __LP64__
      fill_vvec_zero(acc8_vec_ct, case_a1_ct_acc8);
      max_incr4 = 0;
      max_incr8 = 0;
#endif
      for (fs_idx = 0; fs_idx < family_all_case_children_ct; fs_idx++, cur_dfam_ptr = &(cur_dfam_ptr[sibling_ct])) {
	paternal_id = *cur_dfam_ptr++;
	maternal_id = *cur_dfam_ptr++;
	sibling_ct = *cur_dfam_ptr++;
	paternal_geno = EXTRACT_2BIT_GENO(loadbuf_ptr, paternal_id);
	maternal_geno = EXTRACT_2BIT_GENO(loadbuf_ptr, maternal_id);
	parental_a1_ct = dfam_allele_ct_table[paternal_geno * 4 + maternal_geno];
	// skip if parent has missing genotype, or neither parent is het
	if (!parental_a1_ct) {
	  continue;
	}

	for (sib_idx = 0, nonmissing_sib_ct = 0; sib_idx < sibling_ct; sib_idx++) {
	  sample_idx = cur_dfam_ptr[sib_idx];
	  cur_geno = EXTRACT_2BIT_GENO(loadbuf_ptr, sample_idx);
	  nonmissing_sib_ct += (cur_geno != 1);
	}
	// skip if all children have missing genotypes
	if (!nonmissing_sib_ct) {
	  continue;
	}

	cur_case_a1_ct_flip[0] = 0;
	cur_case_a1_ct_flip[1] = 0;
	for (sib_idx = 0; sib_idx < sibling_ct; sib_idx++) {
	  sample_idx = cur_dfam_ptr[sib_idx];
	  cur_geno = EXTRACT_2BIT_GENO(loadbuf_ptr, sample_idx);
	  if (cur_geno == 1) {
	    continue;
	  }
	  cur_case_a1_ct_flip[0] += (4 - cur_geno) / 2;
	}
	cur_case_a1_ct_flip[1] = nonmissing_sib_ct * parental_a1_ct - cur_case_a1_ct_flip[0];
	quad_denom += (2 - (parental_a1_ct & 1)) * nonmissing_sib_ct;
	twice_numer_subtract += nonmissing_sib_ct * parental_a1_ct;

#ifdef __LP64__
	cur_max_incr = MAXV(cur_case_a1_ct_flip[0], cur_case_a1_ct_flip[1]);
	max_incr8 += cur_max_incr;
	// also tried 16-bit accumulators, but that has ~50% greater runtime on
	// typical datasets
	if (max_incr8 >= 256) {
	  if (max_incr4) {
	    loader = _mm_set1_epi8(max_incr4);
	    acc8_ptr = case_a1_ct_acc8;
	    for (vidx = 0; vidx < acc8_vec_ct; vidx++) {
	      *acc8_ptr = _mm_add_epi8(*acc8_ptr, loader);
	      acc8_ptr++;
	    }
	  }
	  unroll_zero_incr_8_32(case_a1_ct_acc8, (__m128i*)total_counts, acc8_vec_ct);
	  max_incr8 = cur_max_incr;
	  max_incr4 = 0;
	}
	if (cur_max_incr < 256) {
          max_incr4 += cur_case_a1_ct_flip[0];
	  diff_vec = _mm_set1_epi8((uint8_t)(cur_case_a1_ct_flip[1] - cur_case_a1_ct_flip[0]));
	  acc8_ptr = case_a1_ct_acc8;
	  flipa_perm_ptr = (__m128i*)(&(flipa_shuffled[fs_idx * perm_vec_wcta]));
	  for (vidx = 0; vidx < perm_vec_ct128; vidx++) {
	    loader = *flipa_perm_ptr++;
	    for (uii = 0; uii < 8; uii++) {
	      // set incr8 to (cur_case_a1_ct_flip[1] - cur_case_a1_ct_flip[0])
	      // where (specially permuted) flipA is set, zero when it is not
	      incr8 = _mm_and_si128(_mm_sub_epi8(_mm_setzero_si128(), _mm_and_si128(loader, m1x8)), diff_vec);
	      *acc8_ptr = _mm_add_epi8(*acc8_ptr, incr8);
	      acc8_ptr++;
	      loader = _mm_srli_epi64(loader, 1);
	    }
	  }
	} else {
	  cur_flipa = &(flipa[fs_idx * perm_vec_wcta]);
	  for (uii = 0; uii < perm_vec_ct; uii++) {
	    is_flipped = IS_SET(cur_flipa, uii);
	    total_counts[uii] += cur_case_a1_ct_flip[is_flipped];
	  }
	  max_incr8 = 0;
	}
#else
	cur_flipa = &(flipa[fs_idx * perm_vec_wcta]);
	for (uii = 0; uii < perm_vec_ct; uii++) {
	  is_flipped = IS_SET(cur_flipa, uii);
	  total_counts[uii] += cur_case_a1_ct_flip[is_flipped];
	}
#endif
      }
      for (perm_idx = 0; perm_idx < perm_vec_ct; perm_idx++) {
	twice_numers[perm_idx] = 2 * total_counts[perm_idx] - twice_numer_subtract;
      }
      for (fs_idx = 0; fs_idx < family_mixed_ct; fs_idx++, cur_dfam_ptr = &(cur_dfam_ptr[sibling_ct])) {
	paternal_id = *cur_dfam_ptr++;
	maternal_id = *cur_dfam_ptr++;
	sibling_ct = *cur_dfam_ptr++;
	paternal_geno = EXTRACT_2BIT_GENO(loadbuf_ptr, paternal_id);
	maternal_geno = EXTRACT_2BIT_GENO(loadbuf_ptr, maternal_id);
	parental_a1_ct = dfam_allele_ct_table[paternal_geno * 4 + maternal_geno];
	if (!parental_a1_ct) {
	  dfam_sibship_or_unrelated_perm_calc(loadbuf_ptr, cur_dfam_ptr, perm_vecst, orig_pheno_c, sibling_ct, 0, perm_vec_ct, acc4, acc8, cur_case_a1_cts, cur_case_missing_cts, twice_numers, numers, denoms, total_counts);
	} else {
	  for (sib_idx = 0, nonmissing_sib_ct = 0; sib_idx < sibling_ct; sib_idx++) {
	    sample_idx = cur_dfam_ptr[sib_idx];
	    cur_geno = EXTRACT_2BIT_GENO(loadbuf_ptr, sample_idx);
	    nonmissing_sib_ct += (cur_geno != 1);
	  }
	  // skip if all children have missing genotypes
	  if (!nonmissing_sib_ct) {
	    continue;
	  }

	  quad_denom += (2 - (parental_a1_ct & 1)) * nonmissing_sib_ct;
	  cur_flipa = &(flipa[fs_idx * perm_vec_wcta]);
	  fill_uint_zero(perm_vec_ct, cur_case_a1_cts);
#ifdef __LP64__
	  fill_vvec_zero(acc4_vec_ct, acc4);
	  fill_vvec_zero(acc8_vec_ct, acc8);
#else
	  fill_ulong_zero(acc4_word_ct, acc4);
	  fill_ulong_zero(acc8_word_ct, acc8);
#endif
	  // compute (unflipped) case_a1_ct for each permutation
	  max_incr4 = 0; // maximum possible value in acc4
	  max_incr8 = 0; // maximum possible value in acc8
	  orig_case_ct = 0;
	  for (sib_idx = 0; sib_idx < sibling_ct; sib_idx++) {
	    sample_idx = cur_dfam_ptr[sib_idx];
	    orig_case_ct += IS_SET(orig_pheno_c, sample_idx);
	    cur_geno = EXTRACT_2BIT_GENO(loadbuf_ptr, sample_idx);
	    // nothing to do here when cur_geno == 3, since a1_ct is zero
	    if (cur_geno & 1) {
	      continue;
	    }
	    cur_max_incr = (4 - cur_geno) / 2;
#ifdef __LP64__
	    if (max_incr4 + cur_max_incr > 15) {
	      unroll_zero_incr_4_8(acc4, acc8, acc4_vec_ct);
	      max_incr8 += max_incr4;
	      if (max_incr8 > 240) {
	        unroll_zero_incr_8_32(acc8, (__m128i*)cur_case_a1_cts, acc8_vec_ct);
		max_incr8 = 0;
	      }
	      max_incr4 = 0;
	    }
	    max_incr4 += cur_max_incr;

	    pheno_perm_ptr = (const __m128i*)(&(perm_vecst[sample_idx * perm_vec_wcta]));
	    if (cur_max_incr == 1) {
	      unroll_incr_1_4(pheno_perm_ptr, acc4, perm_vec_ct128);
	    } else {
	      // add 2 whenever this sample is a case
	      acc4_ptr = acc4;
	      for (vidx = 0; vidx < acc4_vec_ct; vidx++) {
		loader = *pheno_perm_ptr++;
		acc4_ptr[0] = _mm_add_epi64(acc4_ptr[0], _mm_slli_epi64(_mm_and_si128(loader, m1x4), 1));
		acc4_ptr[1] = _mm_add_epi64(acc4_ptr[1], _mm_and_si128(loader, m1x4ls1));
		loader = _mm_srli_epi64(loader, 1);
		acc4_ptr[2] = _mm_add_epi64(acc4_ptr[2], _mm_and_si128(loader, m1x4ls1));
		loader = _mm_srli_epi64(loader, 1);
		acc4_ptr[3] = _mm_add_epi64(acc4_ptr[3], _mm_and_si128(loader, m1x4ls1));
		acc4_ptr = &(acc4_ptr[4]);
	      }
	    }
#else
	    if (max_incr4 + cur_max_incr > 15) {
	      unroll_zero_incr_4_8(acc4, acc8, acc4_word_ct);
	      max_incr8 += max_incr4;
	      if (max_incr8 > 240) {
		unroll_zero_incr_8_32(acc8, (uintptr_t*)cur_case_a1_cts, acc8_word_ct);
		max_incr8 = 0;
	      }
	      max_incr4 = 0;
	    }
	    max_incr4 += cur_max_incr;

	    pheno_perm_ptr = &(perm_vecst[sample_idx * perm_vec_wcta]);
	    if (cur_max_incr == 1) {
	      unroll_incr_1_4(pheno_perm_ptr, acc4, perm_vec_wct);
	    } else {
	      acc4_ptr = acc4;
	      for (widx = 0; widx < perm_vec_wct; widx++) {
		loader = *pheno_perm_ptr++;
		acc4_ptr[0] += (loader & 0x11111111U) << 1;
		acc4_ptr[1] += loader & 0x22222222U;
		acc4_ptr[2] += (loader >> 1) & 0x22222222U;
		acc4_ptr[3] += (loader >> 2) & 0x22222222U;
		acc4_ptr = &(acc4_ptr[4]);
	      }
	    }
#endif
	  }
#ifdef __LP64__
	  // max_incr4 guaranteed to be nonzero unless no child had any A1
	  // alleles
	  if (max_incr4) {
	    unroll_incr_4_8(acc4, acc8, acc4_vec_ct);
	    unroll_incr_8_32(acc8, (__m128i*)cur_case_a1_cts, acc8_vec_ct);
	  }
#else
	  if (max_incr4) {
	    unroll_incr_4_8(acc4, acc8, acc4_word_ct);
	    unroll_incr_8_32(acc8, (uintptr_t*)cur_case_a1_cts, acc8_word_ct);
	  }
#endif
	  if (nonmissing_sib_ct == sibling_ct) {
	    cur_flipa = &(flipa[fs_idx * perm_vec_wcta]);
	    cur_max_incr = orig_case_ct * parental_a1_ct;
	    for (perm_idx = 0; perm_idx < perm_vec_ct; perm_idx++) {
	      uii = cur_case_a1_cts[perm_idx];
	      if (IS_SET(cur_flipa, perm_idx)) {
		uii = cur_max_incr - uii;
	      }
	      total_counts[perm_idx] += uii;
	      twice_numers[perm_idx] += (int32_t)(2 * uii) - (int32_t)cur_max_incr;
	    }
	  } else {
	    // cur_case_ct also varies; need to compute case_missing_ct for
	    // each permutation, and twice_numers/total_counts updates are more
	    // complex.
	    // (technically could separate out >50% missingness as a special
	    // case, but we focus our attention on the far more common sparse
	    // missingness scenario.)
	    cur_max_incr = 0;
	    fill_uint_zero(perm_vec_ct, cur_case_missing_cts);
#ifdef __LP64__
	    fill_vvec_zero(acc4_vec_ct, acc4);
	    fill_vvec_zero(acc8_vec_ct, acc8);
#else
	    fill_ulong_zero(acc4_word_ct, acc4);
	    fill_ulong_zero(acc8_word_ct, acc8);
#endif
	    for (sib_idx = 0; sib_idx < sibling_ct; sib_idx++) {
	      sample_idx = cur_dfam_ptr[sib_idx];
	      cur_geno = EXTRACT_2BIT_GENO(loadbuf_ptr, sample_idx);
	      if (cur_geno != 1) {
		continue;
	      }
#ifdef __LP64__
	      if (!(cur_max_incr % 15)) {
		unroll_zero_incr_4_8(acc4, acc8, acc4_vec_ct);
		if (!(cur_max_incr % 255)) {
		  unroll_zero_incr_8_32(acc8, (__m128i*)cur_case_missing_cts, acc8_vec_ct);
		}
	      }
	      unroll_incr_1_4((__m128i*)(&(perm_vecst[sample_idx * perm_vec_wcta])), acc4, perm_vec_ct128);
#else
	      if (!(cur_max_incr % 15)) {
		unroll_zero_incr_4_8(acc4, acc8, acc4_word_ct);
		if (!(cur_max_incr % 255)) {
		  unroll_zero_incr_8_32(acc8, (uintptr_t*)cur_case_missing_cts, acc8_word_ct);
		}
	      }
	      unroll_incr_1_4(&(perm_vecst[sample_idx * perm_vec_wcta]), acc4, perm_vec_wct);
#endif
	      cur_max_incr++;
	    }
#ifdef __LP64__
	    unroll_incr_4_8(acc4, acc8, acc4_vec_ct);
	    unroll_incr_8_32(acc8, (__m128i*)cur_case_missing_cts, acc8_vec_ct);
#else
	    unroll_incr_4_8(acc4, acc8, acc4_word_ct);
	    unroll_incr_8_32(acc8, (uintptr_t*)cur_case_missing_cts, acc8_word_ct);
#endif
	    for (perm_idx = 0; perm_idx < perm_vec_ct; perm_idx++) {
	      uii = cur_case_a1_cts[perm_idx];
	      ujj = (orig_case_ct - cur_case_missing_cts[perm_idx]) * parental_a1_ct;
	      twice_numers[perm_idx] += (int32_t)(2 * uii) - ((int32_t)ujj);
	      total_counts[perm_idx] += uii;
	    }
	  }
	}
      }
      for (fs_idx = 0; fs_idx < sibship_mixed_ct; fs_idx++, cur_dfam_ptr = &(cur_dfam_ptr[sibling_ct])) {
	sibling_ct = *cur_dfam_ptr++;
	dfam_sibship_or_unrelated_perm_calc(loadbuf_ptr, cur_dfam_ptr, perm_vecst, orig_pheno_c, sibling_ct, 0, perm_vec_ct, acc4, acc8, cur_case_a1_cts, cur_case_missing_cts, twice_numers, numers, denoms, total_counts);
      }
      for (unrelated_cluster_idx = 0; unrelated_cluster_idx < unrelated_cluster_ct; unrelated_cluster_idx++, cur_dfam_ptr = &(cur_dfam_ptr[sibling_ct])) {
	sibling_ct = *cur_dfam_ptr++;
	// call sibling permutation routine with unrelated bool set (most of
	// the code should be identical so this should be one function)
	dfam_sibship_or_unrelated_perm_calc(loadbuf_ptr, cur_dfam_ptr, perm_vecst, orig_pheno_c, sibling_ct, 1, perm_vec_ct, acc4, acc8, cur_case_a1_cts, cur_case_missing_cts, twice_numers, numers, denoms, total_counts);
      }
      if (perm_adapt) {
	for (perm_idx = 0; perm_idx < perm_vec_ct;) {
	  // now harvest the chi-square values, check adaptive termination
	  // condition, etc.
	  dxx = numers[perm_idx] + ((double)((int32_t)twice_numers[perm_idx])) * 0.5;
	  dyy = denoms[perm_idx] + ((double)((int32_t)quad_denom)) * 0.25;
	  chisq = dxx * dxx / dyy;
	  if (chisq > chisq_high) {
	    success_2incr += 2;
	  } else if (chisq > chisq_low) {
	    success_2incr++;
	  }
	  if (++perm_idx == next_adapt_check - pidx_offset) {
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
      } else {
	for (perm_idx = 0; perm_idx < perm_vec_ct; ++perm_idx) {
	  dxx = numers[perm_idx] + ((double)((int32_t)twice_numers[perm_idx])) * 0.5;
	  dyy = denoms[perm_idx] + ((double)((int32_t)quad_denom)) * 0.25;
	  chisq = dxx * dxx / dyy;
	  if (chisq > chisq_high) {
	    success_2incr += 2;
	  } else if (chisq > chisq_low) {
	    success_2incr++;
	  }
	  if (maxt_results[perm_idx] < chisq) {
	    maxt_results[perm_idx] = chisq;
	  }
	  if (msa_ptr) {
	    *msa_ptr++ = chisq;
	  }
	}
      }
      perm_2success_ct[marker_idx] += success_2incr;
    }
  dfam_perm_thread_skip_all:
    if ((!tidx) || g_is_last_thread_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

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

#ifdef __LP64__
void dfam_flipa_shuffle(uintptr_t* perms, uintptr_t* shuffled_perms, uint32_t perm_ct) {
  // 0 16 32 48 64 80 96 112 4 20 36 52 68 84 100 116 8 24 40 56 72 88 104 120 12 28 44 60 76 92 108 124
  // 1 17 ...
  uint32_t vct = BITCT_TO_VECCT(perm_ct);
  uint32_t vidx;
  uint32_t offset1;
  uint32_t offset8;
  uint32_t read_offset;
  uint32_t write_offset;
  for (vidx = 0; vidx < vct; ++vidx) {
    shuffled_perms[0] = 0;
    shuffled_perms[1] = 0;
    for (offset1 = 0; offset1 < 8; offset1++) {
      for (offset8 = 0; offset8 < 4; offset8++) {
	read_offset = offset1 * 16 + offset8 * 4;
	write_offset = offset1 + offset8 * 8;
	shuffled_perms[0] |= IS_SET(perms, read_offset) << write_offset;
	shuffled_perms[0] |= IS_SET(perms, read_offset + 1) << (write_offset + 32);
	shuffled_perms[1] |= IS_SET(perms, read_offset + 2) << write_offset;
	shuffled_perms[1] |= IS_SET(perms, read_offset + 3) << (write_offset + 32);
      }
    }
    perms = &(perms[2]);
    shuffled_perms = &(shuffled_perms[2]);
  }
}
#endif

int32_t dfam(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, double pfilter, double output_min_p, uint32_t mtest_adjust, double adjust_lambda, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude_orig, uintptr_t marker_ct_orig, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, uintptr_t* marker_reverse, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, Aperm_info* apip, uint32_t mperm_save, uintptr_t* pheno_c, uintptr_t* founder_info, uintptr_t* sex_nm, uintptr_t* sex_male, char* sample_ids, uintptr_t max_sample_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, Chrom_info* chrom_info_ptr, uint32_t hh_exists, uint32_t within_cmdflag, uint32_t perm_batch_size, Family_info* fam_ip, Set_info* sip) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  FILE* outfile_msa = nullptr;
  char* textbuf = g_textbuf;
  uintptr_t marker_ct_orig_autosomal = marker_ct_orig;
  uintptr_t unfiltered_marker_ctl = BITCT_TO_WORDCT(unfiltered_marker_ct);
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t unfiltered_sample_ctl2 = QUATERCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t unfiltered_sample_ctp1l2 = 1 + (unfiltered_sample_ct / BITCT2);
  uintptr_t final_mask = get_final_mask(unfiltered_sample_ct);
  uintptr_t perm_vec_ct128 = 0;
  uintptr_t perm_vec_cta128 = 0;
  uintptr_t perm_vec_wct = 0;
  uintptr_t perm_vec_wcta = 0;
  uintptr_t perm_vec_ctcl8m = 0;
  uintptr_t* marker_exclude_orig_autosomal = marker_exclude_orig;
  uintptr_t* founder_pnm = nullptr;
  double* orig_chisq = nullptr;
  double* maxt_extreme_stat = nullptr;
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
  uint32_t max_thread_ct = MINV(g_thread_ct, MODEL_BLOCKSIZE);
  uint32_t perm_pass_idx = 0;
  uint32_t perms_total = 0;
  uint32_t dfam_cluster_map_size = 0;
  int32_t retval = 0;
  uintptr_t* pheno_nm;
  uintptr_t* dfam_pheno_c;
  uintptr_t* loadbuf_raw;
  uintptr_t* loadbuf_ptr;
  uintptr_t* workbuf;
  uintptr_t* marker_exclude;
  uintptr_t* dfam_sample_exclude;
  uintptr_t* size_one_sibships;
  uint32_t mu_table[MODEL_BLOCKSIZE];
  char* outname_end2;
  char* wptr;
  char* wptr_start;
  uint64_t* family_list;
  uint64_t* trio_list;
  uint32_t* trio_error_lookup;
  uint32_t* fs_starts;
  uint32_t* fss_contents;
  uint32_t* sample_to_fss_idx;
  uint32_t* dfam_iteration_order;
  uint32_t* idx_to_uidx;
  uint32_t* sample_uidx_to_idx;
  uint32_t* sample_to_cluster;
  uint32_t* cluster_ctrl_case_cts;
  uint32_t* cluster_write_idxs;
  uint32_t* cur_dfam_ptr;
  uint32_t* dfam_mixed_start;
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
  double chisq;
  double pval;
  double dxx;
  double dyy;
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
  uint32_t dfam_sample_ctl;
  uint32_t dfam_sample_ctl2;
  uint32_t dfam_sample_ctv;
  uint32_t chrom_fo_idx;
  uint32_t chrom_end;
  uint32_t chrom_idx;
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
      logerrprint("Error: No variants remaining for DFAM analysis.\n");
      goto dfam_ret_INVALID_CMDLINE;
    }
    marker_ct_orig_autosomal -= uii;
    if (bigstack_alloc_ul(unfiltered_marker_ctl, &marker_exclude_orig_autosomal)) {
      goto dfam_ret_NOMEM;
    }
    memcpy(marker_exclude_orig_autosomal, marker_exclude_orig, unfiltered_marker_ctl * sizeof(intptr_t));
    for (chrom_fo_idx = 0; chrom_fo_idx < chrom_info_ptr->chrom_ct; chrom_fo_idx++) {
      chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
      if (is_set(chrom_info_ptr->haploid_mask, chrom_idx) || ((int32_t)chrom_idx == chrom_info_ptr->xymt_codes[MT_OFFSET])) {
	uii = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx];
	fill_bits(uii, chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx + 1] - uii, marker_exclude_orig_autosomal);
      }
    }
  } else if (is_set(chrom_info_ptr->haploid_mask, 0)) {
    logerrprint("Error: DFAM test does not support haploid data.\n");
    goto dfam_ret_INVALID_CMDLINE;
  }
  // bugfix (3 Jan 2020): last argument needs to be a word, not bit, count
  uii = popcount_longs_exclude(pheno_c, sample_exclude, unfiltered_sample_ctl);
  if (!uii) {
    logerrprint("Error: DFAM test requires at least one case.\n");
    goto dfam_ret_INVALID_CMDLINE;
  }
  marker_exclude = marker_exclude_orig_autosomal;
  marker_ct = marker_ct_orig_autosomal;

  // PLINK 1.07 treats missing phenotypes as controls here
  if (bigstack_alloc_ul(unfiltered_sample_ctl, &pheno_nm)) {
    goto dfam_ret_NOMEM;
  }
  bitarr_invert_copy(sample_exclude, unfiltered_sample_ct, pheno_nm);
  if (is_set_test) {
    if (bigstack_alloc_ul(unfiltered_sample_ctl, &founder_pnm)) {
      goto dfam_ret_NOMEM;
    }
    memcpy(founder_pnm, pheno_nm, unfiltered_sample_ctl * sizeof(intptr_t));
    bitvec_and(founder_info, unfiltered_sample_ctl, founder_pnm);
    if (extract_set_union_unfiltered(sip, nullptr, unfiltered_marker_ct, marker_exclude_orig_autosomal, &marker_exclude, &marker_ct)) {
      goto dfam_ret_NOMEM;
    }
  }

  // no --mendel-duos support for now
  retval = get_trios_and_families(unfiltered_sample_ct, sample_exclude, sample_ct, founder_info, sex_nm, sex_male, sample_ids, max_sample_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, nullptr, &max_fid_len, nullptr, nullptr, &family_list, &family_ct, &trio_list, &trio_ct, &trio_error_lookup, 0, multigen);
  if (retval) {
    goto dfam_ret_1;
  }
#ifdef __LP64__
  if ((12 * sample_ct + 2 * family_ct) > 0xffffffffLLU) {
    logerrprint("Error: Too many samples and families for DFAM test.\n");
    goto dfam_ret_INVALID_CMDLINE;
  }
#endif
  if (get_sibship_info(unfiltered_sample_ct, sample_exclude, sample_ct, pheno_nm, founder_info, sample_ids, max_sample_id_len, max_fid_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, family_list, trio_list, family_ct, trio_ct, 0, &size_one_sibships, nullptr, nullptr, &fs_starts, &fss_contents, &sample_to_fss_idx, &fs_ct, nullptr, nullptr)) {
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
  if (bigstack_alloc_ul(unfiltered_sample_ctl, &dfam_sample_exclude) ||
      // shrink this later
      bigstack_alloc_ui(sample_ct + (sample_ct / 2), &dfam_iteration_order) ||
      bigstack_alloc_ui(sample_ct, &idx_to_uidx)) {
    goto dfam_ret_NOMEM;
  }
  fill_all_bits(unfiltered_sample_ct, dfam_sample_exclude);
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
      clear_bit(sample_uidx, dfam_sample_exclude);
      *cur_dfam_ptr++ = sample_uidx;

      sample_uidx = idx_to_uidx[fss_contents[fssc_start - 1]];
      clear_bit(sample_uidx, dfam_sample_exclude);
      *cur_dfam_ptr++ = sample_uidx;

      *cur_dfam_ptr++ = cur_case_ct;
      for (fssc_idx = fssc_start; fssc_idx < fssc_end; fssc_idx++) {
	sample_uidx = idx_to_uidx[fss_contents[fssc_idx]];
	clear_bit(sample_uidx, dfam_sample_exclude);
	*cur_dfam_ptr++ = sample_uidx;
      }
    }
  }
  dfam_mixed_start = cur_dfam_ptr;
  for (fs_idx = 0; fs_idx < family_ct; fs_idx++) {
    // Scan for families with at least one case and one control child.
    fssc_start = fs_starts[fs_idx] + 2;
    fssc_end = fs_starts[fs_idx + 1];
    cur_case_ct = 0;
    for (fssc_idx = fssc_start; fssc_idx < fssc_end; fssc_idx++) {
      cur_case_ct += is_set(pheno_c, idx_to_uidx[fss_contents[fssc_idx]]);
    }
    sibling_ct = fssc_end - fssc_start;
    if (cur_case_ct && (cur_case_ct != sibling_ct)) {
      family_mixed_ct++;
      sample_uidx = idx_to_uidx[fss_contents[fssc_start - 2]];
      clear_bit(sample_uidx, dfam_sample_exclude);
      *cur_dfam_ptr++ = sample_uidx;

      sample_uidx = idx_to_uidx[fss_contents[fssc_start - 1]];
      clear_bit(sample_uidx, dfam_sample_exclude);
      *cur_dfam_ptr++ = sample_uidx;

      dfam_cluster_map_size += sibling_ct;
      *cur_dfam_ptr++ = sibling_ct;
      for (fssc_idx = fssc_start; fssc_idx < fssc_end; fssc_idx++) {
	sample_uidx = idx_to_uidx[fss_contents[fssc_idx]];
	clear_bit(sample_uidx, dfam_sample_exclude);
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
    sibling_ct = fssc_end - fssc_start;
    if (cur_case_ct && (cur_case_ct != sibling_ct)) {
      sibship_mixed_ct++;
      // [0]: sibling ct
      // [1...]: member uidxs
      dfam_cluster_map_size += sibling_ct;
      *cur_dfam_ptr++ = sibling_ct;
      for (fssc_idx = fssc_start; fssc_idx < fssc_end; fssc_idx++) {
	sample_uidx = idx_to_uidx[fss_contents[fssc_idx]];
	clear_bit(sample_uidx, dfam_sample_exclude);
	*cur_dfam_ptr++ = sample_uidx;
      }
    }
  }
  dfam_cluster_map_size = ((uintptr_t)(cur_dfam_ptr - dfam_mixed_start)) - 3 * family_mixed_ct - sibship_mixed_ct;
  if (!no_unrelateds) {
    if (bigstack_alloc_ui(sample_ct, &sample_to_cluster)) {
      goto dfam_ret_NOMEM;
    }
    // --within on an empty file actually causes --dfam to behave differently
    // (no unrelated cluster at all) than no --within at all (one big unrelated
    // cluster) in PLINK 1.07.  Replicate this for now.
    if (within_cmdflag) {
      if (fill_sample_to_cluster(unfiltered_sample_ct, sample_exclude, sample_ct, cluster_ct, cluster_map, cluster_starts, sample_to_cluster, nullptr)) {
	goto dfam_ret_NOMEM;
      }
    } else {
      fill_uint_zero(sample_ct, sample_to_cluster);
      cluster_ct = 1;
    }
    for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
      // Remove families and size-2+ sibships.
      // bugfix: also remove size-1 sibships.
      if ((sample_to_fss_idx[sample_idx] != 0xffffffffU) || IS_SET(size_one_sibships, idx_to_uidx[sample_idx])) {
	sample_to_cluster[sample_idx] = 0xffffffffU;
      }
    }

    if (bigstack_calloc_ui(cluster_ct * 2, &cluster_ctrl_case_cts) ||
        bigstack_alloc_ui(cluster_ct, &cluster_write_idxs)) {
      goto dfam_ret_NOMEM;
    }
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
	uii = cur_ctrl_ct + cur_case_ct;
	cur_dfam_ptr[write_idx++] = uii;
	cluster_write_idxs[unrelated_cluster_idx] = write_idx;
	write_idx += uii;
      }
    }
    dfam_cluster_map_size += write_idx - unrelated_cluster_ct;
    for (sample_uidx = 0, sample_idx = 0; sample_idx < sample_ct; sample_uidx++, sample_idx++) {
      next_unset_unsafe_ck(sample_exclude, &sample_uidx);
      unrelated_cluster_idx = sample_to_cluster[sample_idx];
      if (unrelated_cluster_idx != 0xffffffffU) {
        cur_ctrl_ct = cluster_ctrl_case_cts[2 * unrelated_cluster_idx];
	cur_case_ct = cluster_ctrl_case_cts[2 * unrelated_cluster_idx + 1];
	if (cur_ctrl_ct && cur_case_ct) {
	  uii = cluster_write_idxs[unrelated_cluster_idx];
	  cur_dfam_ptr[uii] = sample_uidx;
	  clear_bit(sample_uidx, dfam_sample_exclude);
	  cluster_write_idxs[unrelated_cluster_idx] = uii + 1;
	}
      }
    }
    cur_dfam_ptr = &(cur_dfam_ptr[write_idx]);
  }
  bigstack_reset((unsigned char*)idx_to_uidx);
  bigstack_shrink_top(dfam_iteration_order, (cur_dfam_ptr - dfam_iteration_order) * sizeof(int32_t));
  dfam_sample_ct = unfiltered_sample_ct - popcount_longs(dfam_sample_exclude, unfiltered_sample_ctl);
  dfam_sample_ctl = BITCT_TO_WORDCT(dfam_sample_ct);
  dfam_sample_ctl2 = QUATERCT_TO_WORDCT(dfam_sample_ct);
  dfam_sample_ctv = BITCT_TO_ALIGNED_WORDCT(dfam_sample_ct);
  if (bigstack_alloc_ui(unfiltered_sample_ct, &sample_uidx_to_idx)) {
    goto dfam_ret_NOMEM;
  }
  fill_uidx_to_idx(dfam_sample_exclude, unfiltered_sample_ct, dfam_sample_ct, sample_uidx_to_idx);
  cur_dfam_ptr = dfam_iteration_order;
  uii = family_all_case_children_ct + family_mixed_ct;
  for (fs_idx = 0; fs_idx < uii; fs_idx++) {
    *cur_dfam_ptr = sample_uidx_to_idx[*cur_dfam_ptr];
    cur_dfam_ptr++;
    *cur_dfam_ptr = sample_uidx_to_idx[*cur_dfam_ptr];
    cur_dfam_ptr++;
    sibling_ct = *cur_dfam_ptr++;
    for (sib_idx = 0; sib_idx < sibling_ct; sib_idx++) {
      *cur_dfam_ptr = sample_uidx_to_idx[*cur_dfam_ptr];
      cur_dfam_ptr++;
    }
  }
  uii = sibship_mixed_ct + unrelated_cluster_ct;
  for (fs_idx = 0; fs_idx < uii; fs_idx++) {
    sibling_ct = *cur_dfam_ptr++;
    for (sib_idx = 0; sib_idx < sibling_ct; sib_idx++) {
      *cur_dfam_ptr = sample_uidx_to_idx[*cur_dfam_ptr];
      cur_dfam_ptr++;
    }
  }
  bigstack_reset((unsigned char*)sample_uidx_to_idx);
  if (bigstack_alloc_ul(dfam_sample_ctl2, &dfam_pheno_c) ||
      bigstack_alloc_ul(unfiltered_sample_ctl2, &loadbuf_raw) ||
      bigstack_alloc_ul(unfiltered_sample_ctp1l2, &workbuf) ||
      bigstack_alloc_ul(MODEL_BLOCKSIZE * ((uintptr_t)dfam_sample_ctl2), &g_loadbuf)) {
    goto dfam_ret_NOMEM;
  }
  copy_bitarr_subset_excl(pheno_c, dfam_sample_exclude, sample_ct, dfam_sample_ct, dfam_pheno_c);
  g_pheno_c = dfam_pheno_c;
  g_dfam_iteration_order = dfam_iteration_order;
  g_dfam_family_all_case_children_ct = family_all_case_children_ct;
  g_dfam_family_mixed_ct = family_mixed_ct;
  g_perm_pheno_nm_ct = dfam_sample_ct;
  g_dfam_sibship_mixed_ct = sibship_mixed_ct;
  g_dfam_unrelated_cluster_ct = unrelated_cluster_ct;
  g_test_type = perm_adapt_nst;
  loadbuf_raw[unfiltered_sample_ctl2 - 1] = 0;
  workbuf[unfiltered_sample_ctp1l2 - 1] = 0;
  for (ulii = 1; ulii <= MODEL_BLOCKSIZE; ulii++) {
    // defensive
    g_loadbuf[dfam_sample_ctl2 * ulii - 1] = 0;
  }
  // no X/haploid/MT, so no haploid filters

  if (fill_orig_chisq) {
    if (bigstack_alloc_d(marker_ct, &orig_chisq)) {
      goto dfam_ret_NOMEM;
    }
    g_orig_stat = orig_chisq;
  }

  g_perm_cluster_ct = family_mixed_ct + sibship_mixed_ct + unrelated_cluster_ct;
  if (do_perms_nst) {
    logerrprint("Error: --dfam permutation tests are currently under development.\n");
    retval = RET_CALC_NOT_YET_SUPPORTED;
    goto dfam_ret_1;
    if (bigstack_alloc_ui(dfam_cluster_map_size, &g_perm_cluster_map) ||
        bigstack_alloc_ui(g_perm_cluster_ct + 1, &g_perm_cluster_starts) ||
        bigstack_alloc_ui(g_perm_cluster_ct, &g_perm_cluster_case_cts) ||
        bigstack_calloc_ul(dfam_sample_ctl, &g_perm_cluster_cc_preimage)) {
      goto dfam_ret_NOMEM;
    }
    cur_dfam_ptr = dfam_mixed_start;
    write_idx = 0;
    for (uii = 0; uii < family_mixed_ct; uii++) {
      g_perm_cluster_starts[uii] = write_idx;
      cur_dfam_ptr = &(cur_dfam_ptr[2]);
      sibling_ct = *cur_dfam_ptr++;
      cur_case_ct = 0;
      for (sib_idx = 0; sib_idx < sibling_ct; sib_idx++) {
	sample_idx = cur_dfam_ptr[sib_idx];
	g_perm_cluster_map[write_idx++] = sample_idx;
	cur_case_ct += IS_SET(dfam_pheno_c, sample_idx);
      }
      if (cur_case_ct * 2 >= sibling_ct) {
	for (sib_idx = 0; sib_idx < sibling_ct; sib_idx++) {
	  SET_BIT(cur_dfam_ptr[sib_idx], g_perm_cluster_cc_preimage);
	}
      }
      cur_dfam_ptr = &(cur_dfam_ptr[sibling_ct]);
      g_perm_cluster_case_cts[uii] = cur_case_ct;
    }
    for (; uii < g_perm_cluster_ct; uii++) {
      g_perm_cluster_starts[uii] = write_idx;
      sibling_ct = *cur_dfam_ptr++;
      cur_case_ct = 0;
      for (sib_idx = 0; sib_idx < sibling_ct; sib_idx++) {
	sample_idx = cur_dfam_ptr[sib_idx];
	g_perm_cluster_map[write_idx++] = sample_idx;
	cur_case_ct += IS_SET(dfam_pheno_c, sample_idx);
      }
      if (cur_case_ct * 2 >= sibling_ct) {
	for (sib_idx = 0; sib_idx < sibling_ct; sib_idx++) {
	  SET_BIT(cur_dfam_ptr[sib_idx], g_perm_cluster_cc_preimage);
	}
      }
      cur_dfam_ptr = &(cur_dfam_ptr[sibling_ct]);
      g_perm_cluster_case_cts[uii] = cur_case_ct;
    }
    if (write_idx != dfam_cluster_map_size) {
      logerrprint("assert failure: write_idx != dfam_cluster_map_size\n");
      exit(1);
    }
    g_perm_cluster_starts[g_perm_cluster_ct] = write_idx;

    retval = cluster_alloc_and_populate_magic_nums(g_perm_cluster_ct, g_perm_cluster_map, g_perm_cluster_starts, &g_perm_tot_quotients, &g_perm_totq_magics, &g_perm_totq_preshifts, &g_perm_totq_postshifts, &g_perm_totq_incrs);
    if (retval) {
      goto dfam_ret_1;
    }
  }

  ulii = 2 * max_marker_allele_len + plink_maxsnp + MAX_ID_SLEN + 256;
  if (ulii > MAXLINELEN) {
    if (bigstack_alloc_c(ulii, &textbuf)) {
      goto dfam_ret_NOMEM;
    }
  }

  // permutation test boilerplate mostly copied from qassoc() in plink_assoc.c,
  // since it's also restricted to autosomes
  g_perms_done = 0;
  g_mperm_save_all = nullptr;
  g_perm_vecs = nullptr;
  if (perm_maxt_nst) {
    perms_total = fam_ip->dfam_mperm_val;
    if (bigstack_calloc_d(perms_total, &maxt_extreme_stat)) {
      goto dfam_ret_NOMEM;
    }
    g_maxt_extreme_stat = maxt_extreme_stat;
    if (mperm_save & MPERM_DUMP_ALL) {
      memcpy(outname_end, ".mperm.dump.all", 16);
      if (fopen_checked(outname, "w", &outfile_msa)) {
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
      if (bigstack_alloc_ui(marker_ct, &g_perm_attempt_ct) ||
          bigstack_calloc_uc(round_up_pow2(marker_ct, BYTECT), &g_perm_adapt_stop)) {
        goto dfam_ret_NOMEM;
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
    }
  }

  outname_end2 = memcpyb(outname_end, ".dfam", 6);
  if (fopen_checked(outname, "w", &outfile)) {
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
    if (bigstack_init_sfmtp(max_thread_ct)) {
      goto dfam_ret_NOMEM;
    }
    g_perm_is_1bit = 1;
  }

  fputs("0%", stdout);
  fflush(stdout);
  // ----- begin main loop -----
 dfam_more_perms:
  if (do_perms_nst) {
    if (perm_adapt_nst && perm_pass_idx) {
      while (g_first_adapt_check <= g_perms_done) {
	// APERM_MAX prevents infinite loop here
	g_first_adapt_check += (int32_t)(apip->init_interval + ((int32_t)g_first_adapt_check) * apip->interval_slope);
      }
    }
    // todo: check whether larger batches make sense
    g_perm_vec_ct = perm_batch_size;
    if (g_perm_vec_ct > perms_total - g_perms_done) {
      g_perm_vec_ct = perms_total - g_perms_done;
    }
    perm_vec_ct128 = (g_perm_vec_ct + 127) / 128;
    perm_vec_cta128 = perm_vec_ct128 * 128;
    perm_vec_wct = BITCT_TO_WORDCT(g_perm_vec_ct);
    perm_vec_wcta = perm_vec_ct128 * (128 / BITCT);
    perm_vec_ctcl8m = round_up_pow2(g_perm_vec_ct, CACHELINE_DBL);

    if (bigstack_alloc_ul(dfam_sample_ct * perm_vec_wcta, &g_dfam_perm_vecst) ||
        bigstack_alloc_ul(g_perm_vec_ct * dfam_sample_ctv, &g_perm_vecs)) {
      goto dfam_ret_NOMEM;
    }
    // initialize phenotype permutations.
    g_perm_generation_thread_ct = MINV(max_thread_ct, g_perm_vec_ct);
    if (spawn_threads(threads, &generate_cc_cluster_perms_thread, g_perm_generation_thread_ct)) {
      goto dfam_ret_THREAD_CREATE_FAIL;
    }
    ulii = 0;
    generate_cc_cluster_perms_thread((void*)ulii);
    join_threads(threads, g_perm_generation_thread_ct);

    transpose_perm1s(g_perm_vecs, g_perm_vec_ct, sample_ct, (uint32_t*)g_dfam_perm_vecst);
    bigstack_reset(g_perm_vecs);

    if (bigstack_alloc_ul(family_ct * perm_vec_wct, &g_dfam_flipa) ||
#ifdef __LP64__
        bigstack_alloc_ul(family_all_case_children_ct * perm_vec_wcta, &g_dfam_flipa_shuffled) ||
#endif
	bigstack_alloc_i(max_thread_ct * perm_vec_cta128, &g_dfam_twice_numers) ||
	bigstack_alloc_ui(max_thread_ct * perm_vec_cta128, &g_dfam_total_counts) ||
	bigstack_alloc_d(max_thread_ct * perm_vec_cta128, &g_dfam_numers) ||
	bigstack_alloc_d(max_thread_ct * perm_vec_cta128, &g_dfam_denoms)) {
      goto dfam_ret_NOMEM;
    }
    // initialize flipa permutations.
    ;;;

    /*
    for () {
    }
    */

#ifdef __LP64__
    for (fs_idx = 0; fs_idx < family_all_case_children_ct; fs_idx++) {
      dfam_flipa_shuffle(&(g_dfam_flipa[fs_idx * perm_vec_wcta]), &(g_dfam_flipa_shuffled[fs_idx * perm_vec_wcta]), g_perm_vec_ct);
    }
#endif
    if (perm_maxt_nst) {
      if (bigstack_alloc_d(max_thread_ct * perm_vec_ctcl8m, &g_maxt_thread_results)) {
	goto dfam_ret_NOMEM;
      }
      if (mperm_save & MPERM_DUMP_ALL) {
	if (bigstack_alloc_d(marker_ct * g_perm_vec_ct, &g_mperm_save_all)) {
	  goto dfam_ret_NOMEM;
	}
      }
    }
  }
  chrom_fo_idx = 0xffffffffU;
  marker_uidx = next_unset_unsafe(marker_exclude, 0);
  if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
    goto dfam_ret_READ_FAIL;
  }
  marker_idx = 0;
  marker_idx2 = 0;
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
      if (load_raw2(unfiltered_sample_ct4, unfiltered_sample_ctl2m1, final_mask, bedfile, loadbuf_raw)) {
	goto dfam_ret_READ_FAIL;
      }
      if (IS_SET(marker_reverse, marker_uidx)) {
	reverse_loadbuf(unfiltered_sample_ct, (unsigned char*)loadbuf_raw);
      }
      erase_mendel_errors(unfiltered_sample_ct, loadbuf_raw, workbuf, sex_male, trio_error_lookup, trio_ct, 0, multigen);
      copy_quaterarr_nonempty_subset_excl(loadbuf_raw, dfam_sample_exclude, unfiltered_sample_ct, dfam_sample_ct, &(g_loadbuf[block_size * dfam_sample_ctl2]));
      if (do_perms_nst) {
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
      // 3. Iterate through mixed sibships.  If all case siblings, or all
      //    control siblings, have missing genotypes, skip.  Otherwise (see
      //    lines 420-456 of PLINK 1.07 dfam.cpp),
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
      //                             (([clst size]^2) * (2 * [clst size] - 1))
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
	  paternal_geno = EXTRACT_2BIT_GENO(loadbuf_ptr, paternal_id);
	  maternal_geno = EXTRACT_2BIT_GENO(loadbuf_ptr, maternal_id);
	  parental_a1_ct = dfam_allele_ct_table[paternal_geno * 4 + maternal_geno];
	  if (!parental_a1_ct) {
	    cur_dfam_ptr = &(cur_dfam_ptr[sibling_ct]);
	    continue;
	  }
	  cur_case_ct = 0;
	  case_a1_ct = 0;
          for (sib_idx = 0; sib_idx < sibling_ct; sib_idx++) {
            sample_idx = *cur_dfam_ptr++;
	    cur_geno = EXTRACT_2BIT_GENO(loadbuf_ptr, sample_idx);
	    if (cur_geno == 1) {
	      continue;
	    }
            cur_case_ct++;
	    case_a1_ct += (4 - cur_geno) / 2;
	  }
	  if (cur_case_ct) {
	    twice_numer += (int32_t)(2 * case_a1_ct) - (int32_t)(cur_case_ct * parental_a1_ct);
	    quad_denom += (2 - (parental_a1_ct & 1)) * cur_case_ct;
	    total_count += case_a1_ct;
	    twice_total_expected += cur_case_ct * parental_a1_ct;
	  }
	}
	for (fs_idx = 0; fs_idx < family_mixed_ct; fs_idx++) {
          paternal_id = *cur_dfam_ptr++;
	  maternal_id = *cur_dfam_ptr++;
	  sibling_ct = *cur_dfam_ptr++;
	  paternal_geno = EXTRACT_2BIT_GENO(loadbuf_ptr, paternal_id);
	  maternal_geno = EXTRACT_2BIT_GENO(loadbuf_ptr, maternal_id);
	  parental_a1_ct = dfam_allele_ct_table[paternal_geno * 4 + maternal_geno];
	  cur_case_ct = 0;
	  cur_ctrl_ct = 0;
	  case_hom_a1_ct = 0;
	  case_het_ct = 0;
	  ctrl_hom_a1_ct = 0;
	  ctrl_het_ct = 0;
          for (sib_idx = 0; sib_idx < sibling_ct; sib_idx++) {
            sample_idx = *cur_dfam_ptr++;
	    cur_geno = EXTRACT_2BIT_GENO(loadbuf_ptr, sample_idx);
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
	    quad_denom += (2 - (parental_a1_ct & 1)) * (cur_case_ct + cur_ctrl_ct);
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
	    cur_geno = EXTRACT_2BIT_GENO(loadbuf_ptr, sample_idx);
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
	    cur_geno = EXTRACT_2BIT_GENO(loadbuf_ptr, sample_idx);
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
	  case_var_a1_ct = case_expected_a1_ct * ((double)((int32_t)(2 * uii - ujj))) * ((double)((int32_t)cur_ctrl_ct)) / (dxx * (2 * dxx - 1));
          numer += case_a1_ct - case_expected_a1_ct;
	  denom += case_var_a1_ct;
	  total_expected += case_expected_a1_ct;
	}
	chisq = numer * numer / denom;
	if (fill_orig_chisq) {
	  orig_chisq[marker_idx + marker_bidx] = chisq;
	}
	pval = chiprob_p(chisq, 1);
	if ((pfilter == 2.0) || ((pval <= pfilter) && (pval >= 0.0))) {
	  wptr = width_force(4, textbuf, chrom_name_write(chrom_info_ptr, get_variant_chrom(chrom_info_ptr, marker_uidx2), textbuf));
	  *wptr++ = ' ';
	  wptr = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx2 * max_marker_id_len]), wptr);
	  *wptr++ = ' ';
	  wptr = fw_strcpy(4, marker_allele_ptrs[2 * marker_uidx2], wptr);
	  *wptr++ = ' ';
	  wptr = fw_strcpy(4, marker_allele_ptrs[2 * marker_uidx2 + 1], wptr);
	  *wptr++ = ' ';
	  wptr = uint32toa_w8x(total_count, ' ', wptr);
	  wptr = dtoa_g_wxp4x(total_expected, 8, ' ', wptr);
	  if (denom != 0.0) {
	    wptr = dtoa_g_wxp4x(chisq, 12, ' ', wptr);
	    wptr = dtoa_g_wxp4(pval, 12, wptr);
	  } else {
	    wptr = memcpya(wptr, "          NA           NA", 25);
	  }
	  wptr = memcpya(wptr, " \n", 2);
	  if (fwrite_checked(textbuf, wptr - textbuf, outfile)) {
	    goto dfam_ret_WRITE_FAIL;
	  }
	}
      }
    }
    if (do_perms_nst) {
      // g_xfam_thread_ct = ;;; // f(block size)
      // ...
      g_perms_done += g_perm_vec_ct;
    }
    marker_idx += block_size;
    if ((!perm_pass_idx) && (marker_idx >= loop_end)) {
      if (marker_idx < marker_unstopped_ct) {
	if (pct >= 10) {
	  putc_unlocked('\b', stdout);
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
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    logprint("done.\n");
    if (do_perms_nst) {
      // bigstack_reset();
    }
    if (fclose_null(&outfile)) {
      goto dfam_ret_WRITE_FAIL;
    }
    if (!is_set_test) {
      if (do_perms_nst) {
	bigstack_reset(g_dfam_perm_vecst);
      }
      if (mtest_adjust) {
	if (bigstack_alloc_ui(marker_ct, &idx_to_uidx)) {
	  goto dfam_ret_NOMEM;
	}
	fill_idx_to_uidx(marker_exclude, unfiltered_marker_ct, marker_ct, idx_to_uidx);
	retval = multcomp(outname, outname_end, idx_to_uidx, marker_ct, marker_ids, max_marker_id_len, plink_maxsnp, chrom_info_ptr, orig_chisq, pfilter, output_min_p, mtest_adjust, 0, adjust_lambda, nullptr, nullptr);
	if (retval) {
	  goto dfam_ret_1;
	}
	bigstack_reset(idx_to_uidx);
      }
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
    bigstack_reset(g_dfam_perm_vecst);
    if (g_perms_done < perms_total) {
      if (perm_adapt_nst) {
	marker_unstopped_ct = marker_ct - popcount_longs((uintptr_t*)g_perm_adapt_stop, (marker_ct + sizeof(intptr_t) - 1) / sizeof(intptr_t));
	if (!marker_unstopped_ct) {
	  goto dfam_adapt_perm_count;
	}
      }
      printf("\r%u permutation%s complete.", g_perms_done, (g_perms_done != 1)? "s" : "");
      fflush(stdout);
      perm_pass_idx++;
      goto dfam_more_perms;
    }
    if (perm_adapt_nst) {
    dfam_adapt_perm_count:
      g_perms_done = 0;
      for (uii = 0; uii < marker_ct; uii++) {
	if (g_perm_attempt_ct[uii] > g_perms_done) {
	  g_perms_done = g_perm_attempt_ct[uii];
	  if (g_perms_done == perms_total) {
	    break;
	  }
	}
      }
    }
    putc_unlocked('\r', stdout);
    LOGPRINTF("%u %s permutation%s complete.\n", g_perms_done, perm_maxt_nst? "max(T)" : "adaptive", (g_perms_done != 1)? "s" : "");
    if (perm_adapt_nst) {
      memcpy(outname_end2, ".perm", 6);
    } else {
      if (mperm_save & MPERM_DUMP_BEST) {
	memcpy(outname_end, ".mperm.dump.best", 17);
	// ...
	memcpy(outname_end, ".qassoc", 7);
      }
      memcpy(outname_end2, ".mperm", 7);
    }
    if (fopen_checked(outname, "w", &outfile)) {
      goto dfam_ret_OPEN_FAIL;
    }
    if (perm_adapt_nst) {
      sprintf(g_textbuf, " CHR %%%us    CHISQ_TDT         EMP1           NP \n", plink_maxsnp);
    } else {
      sprintf(g_textbuf, " CHR %%%us    CHISQ_TDT         EMP1         EMP2 \n", plink_maxsnp);
#ifdef __cplusplus
      std::sort(g_maxt_extreme_stat, &(g_maxt_extreme_stat[perms_total]));
#else
      qsort(g_maxt_extreme_stat, perms_total, sizeof(double), double_cmp);
#endif
    }
    fprintf(outfile, g_textbuf, "SNP");
    chrom_fo_idx = 0xffffffffU;
    marker_uidx = next_unset_unsafe(marker_exclude, 0);
    marker_idx = 0;
    dyy = 1.0 / ((double)((int32_t)perms_total + 1));
    dxx = 0.5 * dyy;
    while (1) {
      do {
	chrom_end = chrom_info_ptr->chrom_fo_vidx_start[(++chrom_fo_idx) + 1U];
      } while (marker_uidx >= chrom_end);
      uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
      wptr_start = width_force(4, g_textbuf, chrom_name_write(chrom_info_ptr, uii, g_textbuf));
      *wptr_start++ = ' ';
      wptr_start[plink_maxsnp] = ' ';
      for (; marker_uidx < chrom_end;) {
	if (perm_adapt_nst) {
	  pval = ((double)(g_perm_2success_ct[marker_idx] + 2)) / ((double)(2 * (g_perm_attempt_ct[marker_idx + 1])));
	} else {
	  pval = ((double)(g_perm_2success_ct[marker_idx + 2])) * dxx;
	}
	if (pval <= pfilter) {
	  fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), wptr_start);
	  wptr = &(wptr_start[1 + plink_maxsnp]);
	  if (perm_adapt_nst && (!g_perm_attempt_ct[marker_idx])) {
	    // invalid
	    wptr = memcpya(wptr, "          NA           NA           NA", 38);
	  } else {
	    wptr = dtoa_g_wxp4x(orig_chisq[marker_idx], 12, ' ', wptr);
	    if (!perm_count) {
	      wptr = dtoa_g_wxp4(pval, 12, wptr);
	    } else {
	      wptr = dtoa_g_wxp4(((double)g_perm_2success_ct[marker_idx]) * 0.5, 12, wptr);
	    }
	    *wptr++ = ' ';
	    if (perm_adapt_nst) {
	      wptr = memseta(wptr, 32, 2);
	      wptr = uint32toa_w10(g_perm_attempt_ct[marker_idx], wptr);
	    } else {
	      // ...
	      if (!perm_count) {
	      } else {
	      }
	    }
	    *wptr++ = '\n';
	    if (fwrite_checked(g_textbuf, wptr - g_textbuf, outfile)) {
	      goto dfam_ret_WRITE_FAIL;
	    }
	  }
	  if (++marker_idx == marker_ct) {
	    goto dfam_loop_end;
	  }
	  marker_uidx++;
	  next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
	}
      }
    }
  dfam_loop_end:
    if (fclose_null(&outfile)) {
      goto dfam_ret_WRITE_FAIL;
    }
    LOGPRINTFWW("Permutation test report written to %s .\n", outname);
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
  dfam_ret_THREAD_CREATE_FAIL:
    retval = RET_THREAD_CREATE_FAIL;
    break;
  }
 dfam_ret_1:
  bigstack_reset(bigstack_mark);
  fclose_cond(outfile);
  fclose_cond(outfile_msa);
  return retval;
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
    // integer modulus is slow.  some of the other permutation generators
    // handle many at once, in a manner that allows magic number divide to be
    // employed efficiently.
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
  fill_all_bits(fss_ct, nm_fss);
  cur_start = *fs_starts_ptr++;
  for (cur_idx = 0; cur_idx < family_ct; cur_idx++) {
    cur_end = *fs_starts_ptr++;
    sample_uidx = fss_contents[cur_start];
    uii = fss_contents[cur_start + 1];
    ulii = EXTRACT_2BIT_GENO(loadbuf, sample_uidx);
    uljj = EXTRACT_2BIT_GENO(loadbuf, uii);
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
        ulii = EXTRACT_2BIT_GENO(loadbuf, sample_uidx);
        if (ulii != 1) {
          uljj += ulii + (ulii == 0);
	} else {
	  sib_ct--;
	}
      } while (fss_ptr < fss_end);
      if (sib_ct) {
        qfam_b[cur_idx] = ((double)((intptr_t)(2 * (uintptr_t)sib_ct) - ((intptr_t)uljj))) / ((double)((int32_t)sib_ct));
      } else {
	clear_bit(cur_idx, nm_fss);
      }
    }
    cur_start = cur_end;
  }
  fss_ptr = &(fss_contents[cur_start]);
  for (; cur_idx < fs_ct; cur_idx++) {
    // sibships
    cur_end = *fs_starts_ptr++;
    sib_ct = cur_end - cur_start;
    fss_end = &(fss_contents[cur_end]);
    uljj = 0;
    do {
      sample_uidx = *fss_ptr++;
      ulii = EXTRACT_2BIT_GENO(loadbuf, sample_uidx);
      if (ulii != 1) {
        uljj += ulii + (ulii == 0);
      } else {
        sib_ct--;
      }
    } while (fss_ptr < fss_end);
    if (sib_ct) {
      qfam_b[cur_idx] = ((double)((intptr_t)(2 * (uintptr_t)sib_ct) - ((intptr_t)uljj))) / ((double)((int32_t)sib_ct));
    } else {
      clear_bit(cur_idx, nm_fss);
    }
    cur_start = cur_end;
  }
  for (; cur_idx < fss_ct; cur_idx++) {
    // singletons
    sample_uidx = *fss_ptr++;
    ulii = EXTRACT_2BIT_GENO(loadbuf, sample_uidx);
    if (ulii != 1) {
      qfam_b[cur_idx] = (double)(2 - (intptr_t)(ulii + (ulii == 0)));
    } else {
      clear_bit(cur_idx, nm_fss);
    }
  }
  fill_all_bits(lm_ct, nm_lm);
  for (sample_uidx = 0, sample_idx = 0; sample_idx < lm_ct; sample_uidx++, sample_idx++) {
    next_set_unsafe_ck(lm_eligible, &sample_uidx);
    ulii = EXTRACT_2BIT_GENO(loadbuf, sample_uidx);
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
        if (EXTRACT_2BIT_GENO(loadbuf, uii) == 1) {
	  goto qfam_compute_bw_skip;
	}
      }
      qfam_w[sample_idx] = ((double)(2 - (intptr_t)(ulii + (ulii == 0)))) - qfam_b[fss_idx];
    } else {
    qfam_compute_bw_skip:
      dxx = pheno_d2[sample_idx];
      qt_sum -= dxx;
      qt_ssq -= dxx * dxx;
      clear_bit(sample_idx, nm_lm);
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
  //   flips, and then patch them for each flip
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
      clear_bit(sample_idx, nm_lm);
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
  uint32_t qfam_thread_ct = g_xfam_thread_ct;
  uint32_t fs_ct = g_fs_ct;
  uint32_t lm_ct = g_lm_ct;
  uint32_t singleton_ct = g_singleton_ct;
  uint32_t fss_ct = fs_ct + singleton_ct;
  uint32_t fss_ctl = BITCT_TO_WORDCT(fss_ct);
  uint32_t lm_ctl = BITCT_TO_WORDCT(lm_ct);
  uint32_t test_type = g_test_type;
  uint32_t only_within = (test_type & (QFAM_WITHIN1 | QFAM_WITHIN2))? 1 : 0;
  uintptr_t* lm_eligible = g_lm_eligible;
  uintptr_t* lm_within2_founder = g_lm_within2_founder;
  uintptr_t* qfam_flip = g_qfam_flip;
  uintptr_t* nm_fss = &(g_nm_fss[tidx * round_up_pow2(fss_ctl, CACHELINE_WORD)]);
  uintptr_t* nm_lm = &(g_nm_lm[tidx * round_up_pow2(lm_ctl, CACHELINE_WORD)]);
  double* qfam_b = &(g_qfam_b[tidx * round_up_pow2(fss_ct, CACHELINE_DBL)]);
  double* qfam_w = &(g_qfam_w[tidx * round_up_pow2(lm_ct, CACHELINE_DBL)]);
  double* pheno_d2 = g_pheno_d2;
  double* beta_sum = g_beta_sum;
  double* beta_ssq = g_beta_ssq;
  uint32_t* qfam_permute = only_within? nullptr : g_qfam_permute;
  uint32_t* permute_edit_buf = only_within? nullptr : (&(g_permute_edit[tidx * round_up_pow2(fss_ct, CACHELINE_INT32)]));
  uint32_t* perm_2success_ct = g_perm_2success_ct;
  uint32_t* perm_attempt_ct = g_perm_attempt_ct;
  uint32_t* fs_starts = g_fs_starts;
  uint32_t* fss_contents = g_fss_contents;
  uint32_t* sample_lm_to_fss_idx = g_sample_lm_to_fss_idx;
  uint32_t* perm_ptr = nullptr;
  uint32_t* beta_fail_cts = g_beta_fail_cts;
  uintptr_t cur_perm_ct = g_cur_perm_ct;
  uintptr_t sample_ct = g_qfam_sample_ct;
  uintptr_t sample_ctl2 = QUATERCT_TO_WORDCT(sample_ct);
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
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uintptr_t unfiltered_sample_ctl2 = QUATERCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t unfiltered_sample_ctp1l2 = 1 + (unfiltered_sample_ct / BITCT2);
  uintptr_t sample_ctl2 = QUATERCT_TO_WORDCT(sample_ct);
  uintptr_t final_mask = get_final_mask(unfiltered_sample_ct);
  double qt_sum_all = 0.0;
  double qt_ssq_all = 0.0;
  double geno_sum = 0.0;
  double geno_ssq = 0.0;
  double qt_g_prod = 0.0;
  double* orig_beta = nullptr;
  char* chrom_name_ptr = nullptr;
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
  int32_t mt_code = chrom_info_ptr->xymt_codes[MT_OFFSET];
  int32_t retval = 0;
  const char qfam_flag_suffixes[][8] = {"within", "parents", "total", "between"};
  const char qfam_test_str[][6] = {"WITH ", " TOT ", " BET "};
  const char* qfam_test_ptr = qfam_test_str[0];
  char chrom_name_buf[3 + MAX_CHROM_TEXTNUM_SLEN];
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
      logerrprint("Error: No variants remaining for QFAM analysis.\n");
      goto qfam_ret_INVALID_CMDLINE;
    }
    marker_ct -= uii;
  } else if (is_set(chrom_info_ptr->haploid_mask, 0)) {
    logerrprint("Error: QFAM test does not currently support haploid data.\n");
    goto qfam_ret_INVALID_CMDLINE;
  }
  // no --mendel-duos support for now
  retval = get_trios_and_families(unfiltered_sample_ct, sample_exclude, sample_ct, founder_info, sex_nm, sex_male, sample_ids, max_sample_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, nullptr, &max_fid_len, nullptr, nullptr, &family_list, &family_ct, &trio_list, &trio_ct, &trio_error_lookup, 0, multigen);
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
    logerrprint("Error: Too many samples and families for QFAM test.\n");
    goto qfam_ret_INVALID_CMDLINE;
  }
#endif
  if (get_sibship_info(unfiltered_sample_ct, sample_exclude, sample_ct, pheno_nm, founder_info, sample_ids, max_sample_id_len, max_fid_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, family_list, trio_list, family_ct, trio_ct, test_type, nullptr, &lm_eligible, &lm_within2_founder, &fs_starts, &fss_contents, &sample_lm_to_fss_idx, &fs_ct, &lm_ct, &singleton_ct)) {
    goto qfam_ret_NOMEM;
  }
  fss_ct = fs_ct + singleton_ct;
  if (fss_ct < 2) {
    logerrprint("Error: QFAM test requires at least two families.\n");
    goto qfam_ret_INVALID_CMDLINE;
  } else if (lm_ct < 3) {
    LOGERRPRINTF("Error: Less than three eligible %ss for QFAM test.\n", (test_type == QFAM_WITHIN1)? "nonfounder" : "sample");
    goto qfam_ret_INVALID_CMDLINE;
  }
  g_fs_starts = fs_starts;
  g_fss_contents = fss_contents;
  g_sample_lm_to_fss_idx = sample_lm_to_fss_idx;
  g_lm_eligible = lm_eligible;
  g_lm_within2_founder = lm_within2_founder;
  g_test_type = test_type;
  g_qfam_sample_ct = sample_ct;
  g_fs_ct = fs_ct;
  g_singleton_ct = singleton_ct;
  g_lm_ct = lm_ct;
  g_xfam_thread_ct = qfam_thread_ct;
  fss_ctl = BITCT_TO_WORDCT(fss_ct);
  lm_ctl = BITCT_TO_WORDCT(lm_ct);
  flip_ctl = only_within? lm_ctl : fss_ctl;

  if (bigstack_calloc_uc(round_up_pow2(marker_ct, BYTECT), &perm_adapt_stop) ||
      bigstack_calloc_ui(marker_ct, &g_perm_2success_ct)) {
    goto qfam_ret_NOMEM;
  }
  g_perm_adapt_stop = perm_adapt_stop;

  if (perm_adapt) {
    g_aperm_alpha = apip->alpha;
    perms_total = apip->max;
    if (bigstack_alloc_ui(marker_ct, &g_perm_attempt_ct)) {
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
    g_perm_attempt_ct = nullptr;
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
  if (bigstack_alloc_ul(unfiltered_sample_ctl2, &loadbuf_raw) ||
      bigstack_alloc_ul(unfiltered_sample_ctp1l2, &workbuf)) {
    goto qfam_ret_NOMEM;
  }
  loadbuf_raw[unfiltered_sample_ctl2 - 1] = 0;
  workbuf[unfiltered_sample_ctp1l2 - 1] = 0;
  if (fopen_checked(outname, "w", &outfile)) {
    goto qfam_ret_OPEN_FAIL;
  }
  if (perms_total < perm_batch_size) {
    perm_batch_size = perms_total;
  }
  if (!only_within) {
    if (bigstack_alloc_ui(perm_batch_size * fss_ct, &g_qfam_permute) ||
        bigstack_alloc_ui(round_up_pow2(fss_ct, CACHELINE_INT32) * qfam_thread_ct, &g_permute_edit)) {
      goto qfam_ret_NOMEM;
    }
  }
  if (emp_se) {
    if (bigstack_alloc_d(marker_ct, &orig_beta) ||
        bigstack_calloc_d(marker_ct, &g_beta_sum) ||
        bigstack_calloc_d(marker_ct, &g_beta_ssq) ||
        bigstack_calloc_ui(marker_ct, &g_beta_fail_cts)) {
      goto qfam_ret_NOMEM;
    }
  } else {
    g_beta_sum = nullptr;
    g_beta_ssq = nullptr;
    g_beta_fail_cts = nullptr;
  }
  if (bigstack_alloc_ul(MODEL_BLOCKSIZE * sample_ctl2, &g_loadbuf) ||
      bigstack_alloc_d(marker_ct, &g_orig_stat) ||
      bigstack_alloc_ul(perm_batch_size * flip_ctl, &g_qfam_flip) ||
      bigstack_alloc_ui(fss_ct - 1, &precomputed_mods) ||
      bigstack_alloc_ul(round_up_pow2(fss_ct, CACHELINE_BIT) * qfam_thread_ct, &nm_fss) ||
      bigstack_alloc_ul(round_up_pow2(lm_ct, CACHELINE_BIT) * qfam_thread_ct, &nm_lm) ||
      bigstack_alloc_d(lm_ct, &pheno_d2) ||
      bigstack_alloc_d(round_up_pow2(fss_ct, CACHELINE_DBL) * qfam_thread_ct, &qfam_b) ||
      bigstack_alloc_d(round_up_pow2(lm_ct, CACHELINE_DBL) * qfam_thread_ct, &qfam_w) ||
      bigstack_alloc_ui(fss_ct, &dummy_perm) ||
      bigstack_alloc_ul(flip_ctl, &dummy_flip)) {
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
  fill_ulong_zero(flip_ctl, dummy_flip);

  LOGPRINTFWW("--qfam-%s: Permuting %" PRIuPTR " families/singletons, and including %u %s in linear regression.\n", flag_suffix, fss_ct, lm_ct, g_species_plural);
  LOGPRINTFWW5("Writing report to %s ... ", outname);
  fputs("0%", stdout);
  fflush(stdout);
  // deliberately rename last field to RAW_P to reduce likelihood of
  // misinterpretation.  --adjust also disabled.
  sprintf(g_textbuf, " CHR %%%us         BP   A1       TEST     NIND       BETA         STAT        RAW_P\n", plink_maxsnp);
  fprintf(outfile, g_textbuf, "SNP");
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
      fill_ulong_zero(cur_perm_ct * lm_ctl, g_qfam_flip);
      ujj = fss_ctl * (BITCT / 32);
      ulptr = g_qfam_flip;
      for (ulii = 0; ulii < cur_perm_ct; ulii++) {
        uiptr = (uint32_t*)dummy_flip;
	for (uii = 0; uii < ujj; uii++) {
          uiptr[uii] = sfmt_genrand_uint32(&g_sfmt);
	}
        for (uii = 0; uii < lm_ct; uii++) {
          if (is_set(dummy_flip, sample_lm_to_fss_idx[uii])) {
	    set_bit(uii, ulptr);
	  }
	}
	ulptr = &(ulptr[lm_ctl]);
      }
      fill_ulong_zero(lm_ctl, dummy_flip);
    } else {
      for (ulii = 0; ulii < cur_perm_ct; ulii++) {
	uint32_permute(&(g_qfam_permute[ulii * fss_ct]), &(precomputed_mods[-1]), &g_sfmt, fss_ct);
      }
      uiptr = (uint32_t*)g_qfam_flip;
      uljj = cur_perm_ct * fss_ctl * (BITCT / 32);
      for (ulii = 0; ulii < uljj; ulii++) {
	*uiptr++ = sfmt_genrand_uint32(&g_sfmt);
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
	      chrom_end = chrom_info_ptr->chrom_fo_vidx_start[(++chrom_fo_idx) + 1U];
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
	if (load_raw2(unfiltered_sample_ct4, unfiltered_sample_ctl2m1, final_mask, bedfile, loadbuf_raw)) {
	  goto qfam_ret_READ_FAIL;
	}
	if (IS_SET(marker_reverse, marker_uidx)) {
	  reverse_loadbuf(unfiltered_sample_ct, (unsigned char*)loadbuf_raw);
	}
	erase_mendel_errors(unfiltered_sample_ct, loadbuf_raw, workbuf, sex_male, trio_error_lookup, trio_ct, 0, multigen);
	loadbuf_ptr = &(g_loadbuf[block_idx * sample_ctl2]);
	copy_quaterarr_nonempty_subset_excl(loadbuf_raw, sample_exclude, unfiltered_sample_ct, sample_ct, loadbuf_ptr);
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
	    chrom_fo_idx = get_variant_chrom_fo_idx(chrom_info_ptr, marker_uidx_cur);
	    chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	    chrom_end = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx + 1];
	    chrom_name_ptr = chrom_name_buf5w4write(chrom_info_ptr, chrom_idx, &chrom_name_len, chrom_name_buf);
	  }
	  bufptr = memcpyax(g_textbuf, chrom_name_ptr, chrom_name_len, ' ');
	  bufptr = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx_cur * max_marker_id_len]), bufptr);
	  *bufptr++ = ' ';
	  bufptr = uint32toa_w10x(marker_pos[marker_uidx_cur], ' ', bufptr);
	  if (fwrite_checked(g_textbuf, bufptr - g_textbuf, outfile)) {
	    goto qfam_ret_WRITE_FAIL;
	  }
	  fputs_w4(marker_allele_ptrs[marker_uidx_cur * 2], outfile);
	  loadbuf_ptr = &(g_loadbuf[block_idx * sample_ctl2]);
	  qfam_compute_bw(loadbuf_ptr, sample_ct, fs_starts, fss_contents, sample_lm_to_fss_idx, lm_eligible, lm_within2_founder, family_ct, fs_ct, singleton_ct, lm_ct, nm_fss, nm_lm, pheno_d2, qt_sum_all, qt_ssq_all, qfam_b, qfam_w, &qt_sum, &qt_ssq);
	  nind = popcount_longs(nm_lm, lm_ctl);
	  bufptr = memseta(g_textbuf, 32, 7);
	  bufptr = memcpya(bufptr, qfam_test_ptr, 5);
	  bufptr = uint32toa_w8x(nind, ' ', bufptr);
	  nind_recip = 1.0 / ((double)((int32_t)nind));
	  if (only_within) {
	    flip_precalc(lm_ct, qfam_w, pheno_d2, nm_lm, &geno_sum, &geno_ssq, &qt_g_prod);
          }
	  if (!qfam_regress(test_type, nind, lm_ct, sample_lm_to_fss_idx, nm_lm, pheno_d2, qfam_b, qfam_w, dummy_perm, dummy_flip, nind_recip, qt_sum, qt_ssq, geno_sum, geno_ssq, qt_g_prod, &beta, &tstat)) {
	    bufptr = dtoa_g_wxp4x(beta, 10, ' ', bufptr);
	    bufptr = dtoa_g_wxp4x(tstat, 12, ' ', bufptr);
	    // do not apply --output-min-p since only the empirical p-value is
	    // supposed to be postprocessed here, not this one
	    bufptr = dtoa_g_wxp4x(calc_tprob(tstat, nind - 2), 12, '\n', bufptr);
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
	  if (fwrite_checked(g_textbuf, bufptr - g_textbuf, outfile)) {
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
	  putc_unlocked('\b', stdout);
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
	putc_unlocked('\b', stdout);
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
  putc_unlocked('\r', stdout);
  memcpy(outname_end, ".perm", 6);
  if (fopen_checked(outname, "w", &outfile)) {
    goto qfam_ret_OPEN_FAIL;
  }
  sprintf(g_textbuf, emp_se? " CHR %%%us         BETA     EMP_BETA       EMP_SE         EMP1           NP \n" : " CHR %%%us         EMP1           NP \n", plink_maxsnp);
  fprintf(outfile, g_textbuf, "SNP");
  chrom_fo_idx = 0xffffffffU;
  chrom_end = 0;
  for (marker_uidx = 0, marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
    next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
    if (marker_uidx >= chrom_end) {
      while (1) {
	do {
	  chrom_end = chrom_info_ptr->chrom_fo_vidx_start[(++chrom_fo_idx) + 1U];
	} while (marker_uidx >= chrom_end);
	chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	if ((!IS_SET(chrom_info_ptr->haploid_mask, chrom_idx)) && (chrom_idx != (uint32_t)mt_code)) {
	  break;
	}
	marker_uidx = next_unset_unsafe(marker_exclude, chrom_end);
      }
      chrom_name_ptr = chrom_name_buf5w4write(chrom_info_ptr, chrom_idx, &chrom_name_len, chrom_name_buf);
    }
    bufptr = memcpyax(g_textbuf, chrom_name_ptr, chrom_name_len, ' ');
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
	bufptr = dtoa_g_wxp4x(orig_beta[marker_idx], 12, ' ', bufptr);
	ukk = ujj - g_beta_fail_cts[marker_idx];
	if (ukk <= 1) {
          bufptr = memcpya(bufptr, "          NA ", 13);
	} else {
	  dxx = g_beta_sum[marker_idx] / ((double)((int32_t)ukk));
	  bufptr = dtoa_g_wxp4x(dxx, 12, ' ', bufptr);
	  dxx = sqrt((g_beta_ssq[marker_idx] - g_beta_sum[marker_idx] * dxx) / ((double)((int32_t)(ukk - 1))));
          bufptr = dtoa_g_wxp4x(dxx, 12, ' ', bufptr);
	}
      }
      if (!perm_count) {
        dxx = ((double)(uii + 2)) / ((double)(2 * (ujj + 1)));
      } else {
	dxx = ((double)uii) * 0.5;
      }
      bufptr = dtoa_g_wxp4(dxx, 12, bufptr);
      bufptr = memseta(bufptr, 32, 3);
      bufptr = uint32toa_w10x(ujj, '\n', bufptr);
    }
    if (fwrite_checked(g_textbuf, bufptr - g_textbuf, outfile)) {
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
  bigstack_reset(bigstack_mark);
  fclose_cond(outfile);
  return retval;
}

// bottom 2 bits of index = child genotype
// middle 2 bits of index = paternal genotype
// top 2 bits of index = maternal genotype

// bottom 2 bits of result = child genotype
// top 2 bits of result = untransmitted genotype
const uintptr_t tucc_table[] =
{0, 5, 5, 5,
 5, 5, 5, 5,
 8, 5, 2, 5,
 5, 5, 10, 5,
 5, 5, 5, 5,
 5, 5, 5, 5,
 5, 5, 5, 5,
 5, 5, 5, 5,
 8, 5, 2, 5,
 5, 5, 5, 5,
 12, 5, 10, 3,
 5, 5, 14, 11,
 5, 5, 10, 5,
 5, 5, 5, 5,
 5, 5, 14, 11,
 5, 5, 5, 15};

int32_t make_pseudocontrols(FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_blen, double* marker_cms, uint32_t* marker_pos, char** marker_allele_ptrs, uintptr_t max_marker_allele_blen, uintptr_t* marker_reverse, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, uintptr_t* founder_info, uintptr_t* sex_nm, uintptr_t* sex_male, char* sample_ids, uintptr_t max_sample_id_blen, char* paternal_ids, uintptr_t max_paternal_id_blen, char* maternal_ids, uintptr_t max_maternal_id_blen, Chrom_info* chrom_info_ptr, Family_info* fam_ip) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = NULL;
  int32_t retval = 0;
  {
    const uint32_t unfiltered_marker_ctl = BITCT_TO_WORDCT(unfiltered_marker_ct);
    const uint32_t sex_mt_variant_ct = count_non_autosomal_markers(chrom_info_ptr, marker_exclude, 1, 1);
    if (sex_mt_variant_ct) {
      LOGPRINTF("Excluding %u X/MT/haploid variant%s from --tucc dataset.\n", sex_mt_variant_ct, (sex_mt_variant_ct == 1)? "" : "s");
      if (sex_mt_variant_ct == marker_ct) {
	logerrprint("Error: No variants remaining for --tucc.\n");
	goto make_pseudocontrols_ret_INVALID_CMDLINE;
      }
      marker_ct -= sex_mt_variant_ct;
      uintptr_t* marker_exclude_new;
      if (bigstack_alloc_ul(unfiltered_marker_ctl, &marker_exclude_new)) {
	goto make_pseudocontrols_ret_NOMEM;
      }
      memcpy(marker_exclude_new, marker_exclude, unfiltered_marker_ctl * sizeof(intptr_t));
      for (uint32_t chrom_fo_idx = 0; chrom_fo_idx < chrom_info_ptr->chrom_ct; ++chrom_fo_idx) {
	const uint32_t chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	if (is_set(chrom_info_ptr->haploid_mask, chrom_idx) || ((int32_t)chrom_idx == chrom_info_ptr->xymt_codes[MT_OFFSET])) {
	  const uint32_t variant_uidx_start = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx];
	  fill_bits(variant_uidx_start, chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx + 1] - variant_uidx_start, marker_exclude_new);
	}
      }
      marker_exclude = marker_exclude_new;
    } else if (is_set(chrom_info_ptr->haploid_mask, 0)) {
      logerrprint("Error: --tucc test does not support haploid data.\n");
      goto make_pseudocontrols_ret_INVALID_CMDLINE;
    }
    const uint32_t multigen = (fam_ip->mendel_modifier / MENDEL_MULTIGEN) & 1;
    uint64_t* family_list;
    uint64_t* trio_list;
    uint32_t* trio_error_lookup;
    uint32_t family_ct;
    uintptr_t trio_ct;
    retval = get_trios_and_families(unfiltered_sample_ct, sample_exclude, sample_ct, founder_info, sex_nm, sex_male, sample_ids, max_sample_id_blen, paternal_ids, max_paternal_id_blen, maternal_ids, max_maternal_id_blen, nullptr, nullptr, nullptr, nullptr, &family_list, &family_ct, &trio_list, &trio_ct, &trio_error_lookup, 0, multigen);
    if (retval) {
      goto make_pseudocontrols_ret_1;
    }
    if (!trio_ct) {
      LOGERRPRINTF("Warning: Skipping --tucc since there are no trios.\n");
      goto make_pseudocontrols_ret_1;
    } else if (trio_ct > 0x3fffffff) {
      logerrprint("Error: Too many trios for --tucc.\n");
      goto make_pseudocontrols_ret_INVALID_CMDLINE;
    }
    // main write buffer size requirements:
    //   write-bed:
    //     MAXLINELEN + max(max_sample_id_blen + 16,
    //                      max_chrom_slen + max_marker_id_blen +
    //                        2 * max_marker_allele_blen + 48,
    //                      trio_ct2 + sizeof(intptr_t) - 1)
    //   ped:
    //     MAXLINELEN + max(max_sample_id_blen + 16,
    //                      2 * max_marker_allele_blen)
    const uint32_t write_bed = fam_ip->tucc_bed;
    const uintptr_t trio_ct2 = (trio_ct + 1) / 2;
    char* chrom_buf = nullptr;
    uintptr_t writebuf_blen = max_sample_id_blen + 16;
    if (write_bed) {
      if (writebuf_blen < trio_ct2 + BYTECT - 1) {
	writebuf_blen = trio_ct2 + BYTECT - 1;
      }
      const uint32_t max_chrom_slen = get_max_chrom_slen(chrom_info_ptr);
      if (bigstack_alloc_c(max_chrom_slen + 1, &chrom_buf)) {
	goto make_pseudocontrols_ret_NOMEM;
      }
      if (writebuf_blen < max_chrom_slen + max_marker_id_blen + 2 * max_marker_id_blen + 48) {
	writebuf_blen = max_chrom_slen + max_marker_id_blen + 2 * max_marker_id_blen + 48;
      }
    } else {
      if (writebuf_blen < 2 * max_marker_allele_blen) {
	writebuf_blen = 2 * max_marker_allele_blen;
      }
    }
    writebuf_blen += MAXLINELEN;
    const uintptr_t unfiltered_sample_ctl2m1 = (unfiltered_sample_ct - 1) / BITCT2;
    const uintptr_t unfiltered_sample_ctp1l2 = 1 + (unfiltered_sample_ct / BITCT2);
    uintptr_t* loadbuf;
    uintptr_t* workbuf;
    char* writebuf;
    if (bigstack_alloc_ul(unfiltered_sample_ctl2m1 + 1, &loadbuf) ||
	bigstack_alloc_ul(unfiltered_sample_ctp1l2, &workbuf) ||
	bigstack_alloc_c(writebuf_blen, &writebuf)) {
      goto make_pseudocontrols_ret_NOMEM;
    }
    loadbuf[unfiltered_sample_ctl2m1] = 0;
    workbuf[unfiltered_sample_ctp1l2 - 1] = 0;
    char* writebuf_flush = &(writebuf[MAXLINELEN]);
    const char* output_missing_geno_ptr = g_output_missing_geno_ptr;
    unsigned char* new_bed_contents = nullptr;
    unsigned char* uwrite_iter;
    if (write_bed) {
      strcpy(outname_end, ".tucc.fam");
      if (fopen_checked(outname, "wb", &outfile)) {
	goto make_pseudocontrols_ret_OPEN_FAIL;
      }
      char* write_iter = writebuf;
      for (uintptr_t trio_idx = 0; trio_idx < trio_ct; ++trio_idx) {
	const uintptr_t child_uidx = (uint32_t)trio_list[trio_idx];
	const char* child_sample_id = (const char*)(&(sample_ids[child_uidx * max_sample_id_blen]));
	const uint32_t child_sample_id_slen = strlen(child_sample_id);
	const char child_sex_char = sexchar(sex_nm, sex_male, child_uidx);
	for (uint32_t is_pseudocontrol = 0; is_pseudocontrol < 2; ++is_pseudocontrol) {
	  write_iter = memcpyax(write_iter, child_sample_id, child_sample_id_slen, '_');
	  *write_iter++ = 'T' + is_pseudocontrol;
	  write_iter = strcpya(write_iter, "\t0\t0\t");
	  *write_iter++ = child_sex_char;
	  *write_iter++ = '\t';
	  *write_iter++ = '2' - is_pseudocontrol;
#ifdef _WIN32
	  *write_iter++ = '\r';
#endif
	  *write_iter++ = '\n';
	  if (write_iter >= writebuf_flush) {
	    if (fwrite_checked(writebuf, write_iter - writebuf, outfile)) {
	      goto make_pseudocontrols_ret_WRITE_FAIL;
	    }
	    write_iter = writebuf;
	  }
	}
      }
      if (write_iter != writebuf) {
	if (fwrite_checked(writebuf, write_iter - writebuf, outfile)) {
	  goto make_pseudocontrols_ret_WRITE_FAIL;
	}
      }
      if (fclose_null(&outfile)) {
	goto make_pseudocontrols_ret_WRITE_FAIL;
      }

      memcpy(&(outname_end[6]), "bi", 2);
      if (fopen_checked(outname, "wb", &outfile)) {
	goto make_pseudocontrols_ret_OPEN_FAIL;
      }
      write_iter = writebuf;
      const char* missing_geno_ptr = g_missing_geno_ptr;
      uint32_t variant_uidx = 0;
      uint32_t chrom_fo_idx = 0xffffffffU;
      uint32_t chrom_end = 0;
      uint32_t chrom_name_blen = 0;
      for (uint32_t variant_idx = 0; variant_idx < marker_ct; ++variant_idx, ++variant_uidx) {
	next_unset_unsafe_ck(marker_exclude, &variant_uidx);
	if (variant_uidx >= chrom_end) {
	  do {
	    ++chrom_fo_idx;
	    chrom_end = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx + 1];
	  } while (variant_uidx >= chrom_end);
	  const uint32_t chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	  char* chrom_name_end = chrom_name_write(chrom_info_ptr, chrom_idx, chrom_buf);
	  *chrom_name_end = '\t';
	  chrom_name_blen = 1 + (uintptr_t)(chrom_name_end - chrom_buf);
	}
	write_iter = memcpya(write_iter, chrom_buf, chrom_name_blen);
	write_iter = strcpyax(write_iter, &(marker_ids[variant_uidx * max_marker_id_blen]), '\t');
	if (!marker_cms) {
	  *write_iter++ = '0';
	} else {
	  write_iter = dtoa_g_wxp8(marker_cms[variant_uidx], 1, write_iter);
	}
	*write_iter++ = '\t';
	write_iter = uint32toa_x(marker_pos[variant_uidx], '\t', write_iter);
	write_iter = strcpyax(write_iter, cond_replace(marker_allele_ptrs[2 * variant_uidx], missing_geno_ptr, output_missing_geno_ptr), '\t');
	write_iter = strcpya(write_iter, cond_replace(marker_allele_ptrs[2 * variant_uidx + 1], missing_geno_ptr, output_missing_geno_ptr));
#ifdef _WIN32
	*write_iter++ = '\r';
#endif
	*write_iter++ = '\n';
	if (write_iter >= writebuf_flush) {
	  if (fwrite_checked(writebuf, write_iter - writebuf, outfile)) {
	    goto make_pseudocontrols_ret_WRITE_FAIL;
	  }
	  write_iter = writebuf;
	}
      }
      if (write_iter != writebuf) {
	if (fwrite_checked(writebuf, write_iter - writebuf, outfile)) {
	  goto make_pseudocontrols_ret_WRITE_FAIL;
	}
      }
      if (fclose_null(&outfile)) {
	goto make_pseudocontrols_ret_WRITE_FAIL;
      }

      memcpy(&(outname_end[7]), "ed", 2);
      if (fopen_checked(outname, "wb", &outfile)) {
	goto make_pseudocontrols_ret_OPEN_FAIL;
      }
      uwrite_iter = (unsigned char*)memcpyl3a(writebuf, "l\x1b\x01");
    } else {
      if (bigstack_alloc_uc(trio_ct2 * marker_ct + BYTECT, &new_bed_contents)) {
	goto make_pseudocontrols_ret_NOMEM;
      }
      uwrite_iter = new_bed_contents;

      // never flush
      writebuf_flush = (char*)(&(new_bed_contents[trio_ct2 * marker_ct + BYTECT]));
    }
    if (fseeko(bedfile, bed_offset, SEEK_SET)) {
      goto make_pseudocontrols_ret_READ_FAIL;
    }
    const uintptr_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
    const uintptr_t final_mask = get_final_mask(unfiltered_sample_ct);
    const uint32_t write_word_ct_m1 = (trio_ct - 1) / (BITCT / 4);
    uint32_t variant_uidx = 0;
    for (uint32_t variant_idx = 0; variant_idx < marker_ct; ++variant_idx, ++variant_uidx) {
      if (!is_set(marker_exclude, variant_uidx)) {
	variant_uidx = next_unset_unsafe(marker_exclude, variant_uidx);
	if (fseeko(bedfile, bed_offset + ((uint64_t)variant_uidx) * unfiltered_sample_ct4, SEEK_SET)) {
	  goto make_pseudocontrols_ret_READ_FAIL;
	}
      }
      if (load_raw2(unfiltered_sample_ct4, unfiltered_sample_ctl2m1, final_mask, bedfile, loadbuf)) {
	goto make_pseudocontrols_ret_READ_FAIL;
      }
      if (IS_SET(marker_reverse, variant_uidx)) {
	reverse_loadbuf(unfiltered_sample_ct, (unsigned char*)loadbuf);
      }
      // 1. iterate through all trios, setting Mendel errors to missing
      // 2. iterate through trio_list:
      //    a. if child or either parent isn't genotyped, save missing genos
      //    b. otherwise, save appropriate 4-bit genotype pair
      erase_mendel_errors(unfiltered_sample_ct, loadbuf, workbuf, sex_male, trio_error_lookup, trio_ct, 0, multigen);
      uint32_t loop_len = BITCT / 4;
      const uint64_t* trio_list_iter = trio_list;
      uintptr_t* uwrite_alias = (uintptr_t*)uwrite_iter;
      uint32_t widx = 0;
      while (1) {
	if (widx >= write_word_ct_m1) {
	  if (widx > write_word_ct_m1) {
	    break;
	  }
	  loop_len = 1 + ((trio_ct - 1) % (BITCT / 4));
	}
	uintptr_t cur_write_word = 0;
	for (uint32_t trio_idx_lowbits = 0; trio_idx_lowbits < loop_len; ++trio_idx_lowbits) {
	  const uint64_t trio_code = *trio_list_iter++;
	  const uint64_t family_code = family_list[trio_code >> 32];
	  const uint32_t table_index = EXTRACT_2BIT_GENO(loadbuf, ((uint32_t)trio_code)) + 4 * EXTRACT_2BIT_GENO(loadbuf, ((uint32_t)family_code)) + 16 * EXTRACT_2BIT_GENO(loadbuf, (family_code >> 32));
	  cur_write_word |= tucc_table[table_index] << (trio_idx_lowbits * 4);
	}
        uwrite_alias[widx++] = cur_write_word;
      }
      uwrite_iter = &(uwrite_iter[trio_ct2]);
      if (uwrite_iter >= ((unsigned char*)writebuf_flush)) {
	if (fwrite_checked(writebuf, uwrite_iter - ((unsigned char*)writebuf), outfile)) {
	  goto make_pseudocontrols_ret_WRITE_FAIL;
	}
	uwrite_iter = (unsigned char*)writebuf;
      }
    }
    if (write_bed) {
      if (uwrite_iter != (unsigned char*)writebuf) {
	if (fwrite_checked(writebuf, uwrite_iter - ((unsigned char*)writebuf), outfile)) {
	  goto make_pseudocontrols_ret_WRITE_FAIL;
	}
      }
      if (fclose_null(&outfile)) {
	goto make_pseudocontrols_ret_WRITE_FAIL;
      }
      outname_end[6] = '\0';
      LOGPRINTFWW("--tucc write-bed: Pseudo cases/controls written to %sbed + %sbim + %sfam .\n", outname, outname, outname);
    } else {
      strcpy(outname_end, ".tucc.ped");
      if (fopen_checked(outname, "wb", &outfile)) {
	goto make_pseudocontrols_ret_OPEN_FAIL;
      }
      writebuf_flush = &(writebuf[MAXLINELEN]);
      char missing4[4];
      missing4[0] = ' ';
      missing4[1] = *output_missing_geno_ptr;
      missing4[2] = ' ';
      missing4[3] = *output_missing_geno_ptr;
      char* write_iter = writebuf;
      for (uintptr_t trio_idx = 0; trio_idx < trio_ct; ++trio_idx) {
	const uintptr_t child_uidx = (uint32_t)trio_list[trio_idx];
	const char* child_fid = (const char*)(&(sample_ids[child_uidx * max_sample_id_blen]));
	const char* iid_start = (const char*)memchr(child_fid, '\t', max_sample_id_blen);
	const uint32_t fid_slen = (uintptr_t)(iid_start - child_fid);
	++iid_start;
	const uint32_t iid_slen = strlen(iid_start);
	const char child_sex_char = sexchar(sex_nm, sex_male, child_uidx);
	for (uint32_t is_pseudocontrol = 0; is_pseudocontrol < 2; ++is_pseudocontrol) {
	  write_iter = memcpyax(write_iter, child_fid, fid_slen, ' ');
	  write_iter = memcpyax(write_iter, iid_start, iid_slen, '_');
	  *write_iter++ = 'T' + is_pseudocontrol;
	  write_iter = strcpya(write_iter, " 0 0 ");
	  *write_iter++ = child_sex_char;
	  *write_iter++ = ' ';
	  *write_iter++ = '2' - is_pseudocontrol;
	  *write_iter++ = ' ';
	  if (write_iter >= writebuf_flush) {
	    if (fwrite_checked(writebuf, write_iter - writebuf, outfile)) {
	      goto make_pseudocontrols_ret_WRITE_FAIL;
	    }
	    write_iter = writebuf;
	  }
	  variant_uidx = 0;
	  const unsigned char* cur_geno_base = &(new_bed_contents[trio_idx / 2]);
	  const uint32_t cur_shift = ((trio_idx % 2) * 2 + is_pseudocontrol) * 2;
	  for (uint32_t variant_idx = 0; variant_idx < marker_ct; ++variant_idx, ++variant_uidx) {
	    next_unset_unsafe_ck(marker_exclude, &variant_uidx);
	    const uint32_t cur_geno = (cur_geno_base[variant_uidx * trio_ct2] >> cur_shift) & 3;
	    if (cur_geno == 1) {
	      write_iter = memcpya(write_iter, missing4, 4);
	    } else {
	      *write_iter++ = ' ';
	      const char* a1_allele = marker_allele_ptrs[2 * variant_uidx];
	      const char* a2_allele = marker_allele_ptrs[2 * variant_uidx + 1];
	      if (cur_geno == 3) {
		write_iter = strcpya(write_iter, a2_allele);
	      } else {
		write_iter = strcpya(write_iter, a1_allele);
	      }
	      *write_iter++ = ' ';
	      if (cur_geno == 0) {
		write_iter = strcpya(write_iter, a1_allele);
	      } else {
		write_iter = strcpya(write_iter, a2_allele);
	      }
	    }
	    if (write_iter >= writebuf_flush) {
	      if (fwrite_checked(writebuf, write_iter - writebuf, outfile)) {
		goto make_pseudocontrols_ret_WRITE_FAIL;
	      }
	      write_iter = writebuf;
	    }
	  }
#ifdef _WIN32
	  *write_iter++ = '\r';
#endif
	  *write_iter++ = '\n';
	}
      }
      if (write_iter != writebuf_flush) {
	if (fwrite_checked(writebuf, write_iter - writebuf, outfile)) {
	  goto make_pseudocontrols_ret_WRITE_FAIL;
	}
      }
      if (fclose_null(&outfile)) {
	goto make_pseudocontrols_ret_WRITE_FAIL;
      }
      LOGPRINTFWW("--tucc: Pseudo cases/controls written to %s .\n", outname);
    }
  }
  while (0) {
  make_pseudocontrols_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  make_pseudocontrols_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  make_pseudocontrols_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  make_pseudocontrols_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  make_pseudocontrols_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  }
 make_pseudocontrols_ret_1:
  bigstack_reset(bigstack_mark);
  fclose_cond(outfile);
  return retval;
}
