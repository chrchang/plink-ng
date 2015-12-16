#include "plink_common.h"

void generate_cc_perm_vec(uint32_t tot_ct, uint32_t set_ct, uint32_t tot_quotient, uint64_t totq_magic, uint32_t totq_preshift, uint32_t totq_postshift, uint32_t totq_incr, uintptr_t* perm_vec, sfmt_t* sfmtp) {
  // Assumes tot_quotient is 2^32 / tot_ct, and
  // totq_magic/totq_preshift/totq_postshift/totq_incr have been precomputed
  // from magic_num().
  uint32_t num_set = 0;
  uint32_t upper_bound = tot_ct * tot_quotient - 1;
  uintptr_t widx;
  uintptr_t wcomp;
  uintptr_t pv_val;
  uint32_t urand;
  uint32_t uii;
  if (set_ct * 2 < tot_ct) {
    fill_ulong_zero(perm_vec, 2 * ((tot_ct + (BITCT - 1)) / BITCT));
    for (; num_set < set_ct; num_set++) {
      do {
	do {
	  urand = sfmt_genrand_uint32(sfmtp);
	} while (urand > upper_bound);
	uii = (totq_magic * ((urand >> totq_preshift) + totq_incr)) >> totq_postshift;
        widx = uii / BITCT2;
	wcomp = ONELU << (2 * (uii % BITCT2));
	pv_val = perm_vec[widx];
      } while (pv_val & wcomp);
      perm_vec[widx] = pv_val | wcomp;
    }
  } else {
    fill_vec_55(perm_vec, tot_ct);
    set_ct = tot_ct - set_ct;
    for (; num_set < set_ct; num_set++) {
      do {
	do {
	  urand = sfmt_genrand_uint32(sfmtp);
	} while (urand > upper_bound);
	uii = (totq_magic * ((urand >> totq_preshift) + totq_incr)) >> totq_postshift;
        widx = uii / BITCT2;
	wcomp = ONELU << (2 * (uii % BITCT2));
	pv_val = perm_vec[widx];
      } while (!(pv_val & wcomp));
      perm_vec[widx] = pv_val - wcomp;
    }
  }
}

void generate_cc_perm1(uint32_t tot_ct, uint32_t set_ct, uint32_t tot_quotient, uint64_t totq_magic, uint32_t totq_preshift, uint32_t totq_postshift, uint32_t totq_incr, uintptr_t* perm_vec, sfmt_t* sfmtp) {
  // generate_cc_perm_vec() variant which uses 1-bit packing instead of 2.
  uint32_t num_set = 0;
  uint32_t upper_bound = tot_ct * tot_quotient - 1;
  uintptr_t widx;
  uintptr_t wcomp;
  uintptr_t pv_val;
  uint32_t urand;
  uint32_t uii;
  if (set_ct * 2 < tot_ct) {
    fill_ulong_zero(perm_vec, (tot_ct + (BITCT - 1)) / BITCT);
    for (; num_set < set_ct; num_set++) {
      do {
	do {
	  urand = sfmt_genrand_uint32(sfmtp);
	} while (urand > upper_bound);
	uii = (totq_magic * ((urand >> totq_preshift) + totq_incr)) >> totq_postshift;
        widx = uii / BITCT;
	wcomp = ONELU << (uii % BITCT);
	pv_val = perm_vec[widx];
      } while (pv_val & wcomp);
      perm_vec[widx] = pv_val | wcomp;
    }
  } else {
    fill_all_bits(perm_vec, tot_ct);
    set_ct = tot_ct - set_ct;
    for (; num_set < set_ct; num_set++) {
      do {
	do {
	  urand = sfmt_genrand_uint32(sfmtp);
	} while (urand > upper_bound);
	uii = (totq_magic * ((urand >> totq_preshift) + totq_incr)) >> totq_postshift;
        widx = uii / BITCT;
	wcomp = ONELU << (uii % BITCT);
	pv_val = perm_vec[widx];
      } while (!(pv_val & wcomp));
      perm_vec[widx] = pv_val - wcomp;
    }
  }
}

void generate_cc_cluster_perm_vec(uint32_t tot_ct, uintptr_t* preimage, uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, uint32_t* cluster_case_cts, uint32_t* tot_quotients, uint64_t* totq_magics, uint32_t* totq_preshifts, uint32_t* totq_postshifts, uint32_t* totq_incrs, uintptr_t* perm_vec, sfmt_t* sfmtp) {
  uint32_t tot_ctl2 = 2 * ((tot_ct + (BITCT - 1)) / BITCT);
  uint32_t cluster_idx;
  uint32_t target_ct;
  uint32_t cluster_end;
  uint32_t* map_ptr;
  uint32_t num_swapped;
  uint32_t cluster_size;
  uint32_t upper_bound;
  uint64_t totq_magic;
  uint32_t totq_preshift;
  uint32_t totq_postshift;
  uint32_t totq_incr;
  uintptr_t widx;
  uintptr_t wcomp;
  uintptr_t pv_val;
  uint32_t urand;
  uint32_t uii;
  memcpy(perm_vec, preimage, tot_ctl2 * sizeof(intptr_t));
  for (cluster_idx = 0; cluster_idx < cluster_ct; cluster_idx++) {
    target_ct = cluster_case_cts[cluster_idx];
    cluster_end = cluster_starts[cluster_idx + 1];
    cluster_size = cluster_end - cluster_starts[cluster_idx];
    if (target_ct && (target_ct != cluster_size)) {
      upper_bound = cluster_size * tot_quotients[cluster_idx] - 1;
      totq_magic = totq_magics[cluster_idx];
      totq_preshift = totq_preshifts[cluster_idx];
      totq_postshift = totq_postshifts[cluster_idx];
      totq_incr = totq_incrs[cluster_idx];
      map_ptr = &(cluster_map[cluster_starts[cluster_idx]]);
      if (target_ct * 2 < cluster_size) {
	for (num_swapped = 0; num_swapped < target_ct; num_swapped++) {
	  do {
	    do {
	      urand = sfmt_genrand_uint32(sfmtp);
	    } while (urand > upper_bound);
	    uii = map_ptr[(uint32_t)((totq_magic * ((urand >> totq_preshift) + totq_incr)) >> totq_postshift)];
	    widx = uii / BITCT2;
	    wcomp = ONELU << (2 * (uii % BITCT2));
	    pv_val = perm_vec[widx];
	  } while (pv_val & wcomp);
	  perm_vec[widx] = pv_val | wcomp;
	}
      } else {
	target_ct = cluster_size - target_ct;
	for (num_swapped = 0; num_swapped < target_ct; num_swapped++) {
	  do {
	    do {
	      urand = sfmt_genrand_uint32(sfmtp);
	    } while (urand > upper_bound);
	    uii = map_ptr[(uint32_t)((totq_magic * ((urand >> totq_preshift) + totq_incr)) >> totq_postshift)];
	    widx = uii / BITCT2;
	    wcomp = ONELU << (2 * (uii % BITCT2));
	    pv_val = perm_vec[widx];
	  } while (!(pv_val & wcomp));
	  perm_vec[widx] = pv_val - wcomp;
	}
      }
    }
  }
}

void generate_cc_cluster_perm1(uint32_t tot_ct, uintptr_t* preimage, uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, uint32_t* cluster_case_cts, uint32_t* tot_quotients, uint64_t* totq_magics, uint32_t* totq_preshifts, uint32_t* totq_postshifts, uint32_t* totq_incrs, uintptr_t* perm_vec, sfmt_t* sfmtp) {
  uint32_t tot_ctl = (tot_ct + (BITCT - 1)) / BITCT;
  uint32_t cluster_idx;
  uint32_t target_ct;
  uint32_t cluster_end;
  uint32_t cluster_size;
  uint32_t* map_ptr;
  uint32_t num_swapped;
  uint32_t upper_bound;
  uint64_t totq_magic;
  uint32_t totq_preshift;
  uint32_t totq_postshift;
  uint32_t totq_incr;
  uintptr_t widx;
  uintptr_t wcomp;
  uintptr_t pv_val;
  uint32_t urand;
  uint32_t uii;
  memcpy(perm_vec, preimage, tot_ctl * sizeof(intptr_t));
  for (cluster_idx = 0; cluster_idx < cluster_ct; cluster_idx++) {
    target_ct = cluster_case_cts[cluster_idx];
    cluster_end = cluster_starts[cluster_idx + 1];
    cluster_size = cluster_end - cluster_starts[cluster_idx];
    if (target_ct && (target_ct != cluster_size)) {
      upper_bound = cluster_size * tot_quotients[cluster_idx] - 1;
      totq_magic = totq_magics[cluster_idx];
      totq_preshift = totq_preshifts[cluster_idx];
      totq_postshift = totq_postshifts[cluster_idx];
      totq_incr = totq_incrs[cluster_idx];
      map_ptr = &(cluster_map[cluster_starts[cluster_idx]]);
      if (target_ct * 2 < cluster_size) {
	for (num_swapped = 0; num_swapped < target_ct; num_swapped++) {
	  do {
	    do {
	      urand = sfmt_genrand_uint32(sfmtp);
	    } while (urand > upper_bound);
	    uii = map_ptr[(uint32_t)((totq_magic * ((urand >> totq_preshift) + totq_incr)) >> totq_postshift)];
	    widx = uii / BITCT;
	    wcomp = ONELU << (uii % BITCT);
	    pv_val = perm_vec[widx];
	  } while (pv_val & wcomp);
	  perm_vec[widx] = pv_val | wcomp;
	}
      } else {
	target_ct = cluster_size - target_ct;
	for (num_swapped = 0; num_swapped < target_ct; num_swapped++) {
	  do {
	    do {
	      urand = sfmt_genrand_uint32(sfmtp);
	    } while (urand > upper_bound);
	    uii = map_ptr[(uint32_t)((totq_magic * ((urand >> totq_preshift) + totq_incr)) >> totq_postshift)];
	    widx = uii / BITCT;
	    wcomp = ONELU << (uii % BITCT);
	    pv_val = perm_vec[widx];
	  } while (!(pv_val & wcomp));
	  perm_vec[widx] = pv_val - wcomp;
	}
      }
    }
  }
}

void transpose_perms(uintptr_t* perm_vecs, uint32_t perm_vec_ct, uint32_t pheno_nm_ct, uint32_t* perm_vecst) {
  // Transpose permutations so PRESTO/PERMORY-style genotype indexing can work.
  //
  // We used a 32-ply interleaved format, to allow counts up to the uint32_t
  // limit without giving up highly parallel adds in the calc_git() inner loop
  // (performed with a combination of unroll_incr_1_4, unroll_incr_4_8, and
  // unroll_incr_8_32).  The index order is:
  // 64-bit build:
  //   first 16 bytes: 0 32 64 96 16 48 80 112 4 36 68 100 20 52 84 116
  //     8 40 72 104 24 56 88 120 12 44 76 108 28 60 92 124 1...
  //   next 16 bytes: 128 160 192...
  //
  // 32-bit build:
  //   first 4 bytes: 0 8 16 24 4 12 20 28 1 9 17 25 5 13 21 29 2 10 18...
  //   next 4 bytes: 32 40 48...
  uintptr_t sample_idx = 0;
  uintptr_t pheno_nm_ctl2 = 2 * ((pheno_nm_ct + (BITCT - 1)) / BITCT);
#ifdef __LP64__
  uint32_t wbuf[4];
  uint32_t* wbptr;
#else
  uint32_t wval;
#endif
  uint32_t rshift;
  uint32_t wshift;
  uintptr_t* pvptr;
  uintptr_t perm_idx;
  for (; sample_idx < pheno_nm_ct; sample_idx++) {
    perm_idx = 0;
    pvptr = &(perm_vecs[sample_idx / BITCT2]);
    rshift = 2 * (sample_idx % BITCT2);
    goto transpose_perms_loop_start;
#ifdef __LP64__
    do {
      if (!(perm_idx % 4)) {
	if (perm_idx % 128) {
	  wshift = ((perm_idx & 96) >> 5) | ((perm_idx & 16) >> 2) | ((perm_idx & 12) << 1);
	} else {
	  memcpy(perm_vecst, wbuf, 16);
	  perm_vecst = &(perm_vecst[4]);
	transpose_perms_loop_start:
	  fill_uint_zero(wbuf, 4);
	  wshift = 0;
	}
	wbptr = wbuf;
      }
      *wbptr |= ((pvptr[perm_idx * pheno_nm_ctl2] >> rshift) & 1) << wshift;
      wbptr++;
    } while (++perm_idx < perm_vec_ct);
    memcpy(perm_vecst, wbuf, 16);
    perm_vecst = &(perm_vecst[4]);
#else
    do {
      if (perm_idx % 32) {
	wshift = ((perm_idx & 24) >> 3) | (perm_idx & 4) | ((perm_idx & 3) << 3);
      } else {
	*perm_vecst++ = wval;
      transpose_perms_loop_start:
	wval = 0;
	wshift = 0;
      }
      wval |= ((pvptr[perm_idx * pheno_nm_ctl2] >> rshift) & 1) << wshift;
    } while (++perm_idx < perm_vec_ct);
    *perm_vecst++ = wval;
#endif
  }
}

void transpose_perm1s(uintptr_t* perm_vecs, uint32_t perm_vec_ct, uint32_t pheno_nm_ct, uint32_t* perm_vecst) {
  uintptr_t sample_idx = 0;
  uintptr_t pheno_nm_ctl = (pheno_nm_ct + (BITCT - 1)) / BITCT;
#ifdef __LP64__
  uint32_t wbuf[4];
  uint32_t* wbptr;
#else
  uint32_t wval;
#endif
  uint32_t rshift;
  uint32_t wshift;
  uintptr_t* pvptr;
  uintptr_t perm_idx;
  for (; sample_idx < pheno_nm_ct; sample_idx++) {
    perm_idx = 0;
    pvptr = &(perm_vecs[sample_idx / BITCT]);
    rshift = sample_idx % BITCT;
    goto transpose_perm1s_loop_start;
#ifdef __LP64__
    do {
      if (!(perm_idx % 4)) {
	if (perm_idx % 128) {
	  wshift = ((perm_idx & 96) >> 5) | ((perm_idx & 16) >> 2) | ((perm_idx & 12) << 1);
	} else {
	  memcpy(perm_vecst, wbuf, 16);
	  perm_vecst = &(perm_vecst[4]);
	transpose_perm1s_loop_start:
	  fill_uint_zero(wbuf, 2);
	  wshift = 0;
	}
	wbptr = wbuf;
      }
      *wbptr |= ((pvptr[perm_idx * pheno_nm_ctl] >> rshift) & 1) << wshift;
      wbptr++;
    } while (++perm_idx < perm_vec_ct);
    memcpy(perm_vecst, wbuf, 16);
    perm_vecst = &(perm_vecst[4]);
#else
    do {
      if (perm_idx % 32) {
	wshift = ((perm_idx & 24) >> 3) | (perm_idx & 4) | ((perm_idx & 3) << 3);
      } else {
	*perm_vecst++ = wval;
      transpose_perm1s_loop_start:
	wval = 0;
	wshift = 0;
      }
      wval |= ((pvptr[perm_idx * pheno_nm_ctl] >> rshift) & 1) << wshift;
    } while (++perm_idx < perm_vec_ct);
    *perm_vecst++ = wval;
#endif
  }
}

// todo: add multithread globals with extern linkage
