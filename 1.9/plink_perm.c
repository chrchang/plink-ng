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

#include "plink_cluster.h"

// Inputs/outputs for multithreaded permutation generators.
uint32_t g_perm_pheno_nm_ct;
uint32_t g_perm_case_ct;
uint32_t g_perm_tot_quotient;
uint64_t g_perm_totq_magic;
uint32_t g_perm_totq_preshift;
uint32_t g_perm_totq_postshift;
uint32_t g_perm_totq_incr;
uint32_t g_perm_is_1bit;
uint32_t g_perm_generation_thread_ct;
uintptr_t g_perm_vec_ct;

uint32_t g_perm_cluster_ct;
uint32_t* g_perm_cluster_map;
uint32_t* g_perm_cluster_starts;
uint32_t* g_perm_cluster_case_cts;
uintptr_t* g_perm_cluster_cc_preimage;
uint32_t* g_perm_tot_quotients;
uint64_t* g_perm_totq_magics;
uint32_t* g_perm_totq_preshifts;
uint32_t* g_perm_totq_postshifts;
uint32_t* g_perm_totq_incrs;

uintptr_t* g_perm_vecs;

// always use genotype indexing for QT --assoc
double* g_perm_vecstd;
double* g_perm_pheno_d2;
uint32_t* g_perm_sample_to_cluster;
uint32_t* g_perm_qt_cluster_thread_wkspace;

// permutation-major instead of sample-major order for --linear (PERMORY
// speedups do not apply)
double* g_perm_pmajor;
uint32_t* g_perm_precomputed_mods; // [n] = 2^32 mod (n-2)

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
    fill_ulong_zero(QUATERCT_TO_ALIGNED_WORDCT(tot_ct), perm_vec);
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
    fill_quatervec_55(tot_ct, perm_vec);
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
    fill_ulong_zero(BITCT_TO_WORDCT(tot_ct), perm_vec);
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
    fill_all_bits(tot_ct, perm_vec);
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
  memcpy(perm_vec, preimage, QUATERCT_TO_ALIGNED_WORDCT(tot_ct) * sizeof(intptr_t));
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
  uint32_t tot_ctl = BITCT_TO_WORDCT(tot_ct);
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

THREAD_RET_TYPE generate_cc_perms_thread(void* arg) {
  intptr_t tidx = (intptr_t)arg;
  uint32_t pheno_nm_ct = g_perm_pheno_nm_ct;
  uint32_t case_ct = g_perm_case_ct;
  uint32_t tot_quotient = g_perm_tot_quotient;
  uint64_t totq_magic = g_perm_totq_magic;
  uint32_t totq_preshift = g_perm_totq_preshift;
  uint32_t totq_postshift = g_perm_totq_postshift;
  uint32_t totq_incr = g_perm_totq_incr;
  uintptr_t* __restrict__ perm_vecs = g_perm_vecs;
  sfmt_t* __restrict__ sfmtp = g_sfmtp_arr[tidx];
  uintptr_t pheno_nm_ctv = BITCT_TO_WORDCT(pheno_nm_ct);
  uint32_t pidx = (((uint64_t)tidx) * g_perm_vec_ct) / g_perm_generation_thread_ct;
  uint32_t pmax = (((uint64_t)tidx + 1) * g_perm_vec_ct) / g_perm_generation_thread_ct;
  if (!g_perm_is_1bit) {
    pheno_nm_ctv *= 2;
    for (; pidx < pmax; pidx++) {
      generate_cc_perm_vec(pheno_nm_ct, case_ct, tot_quotient, totq_magic, totq_preshift, totq_postshift, totq_incr, &(perm_vecs[pidx * pheno_nm_ctv]), sfmtp);
    }
  } else {
    // 16-byte alignment currently isn't needed; but it might be useful in the
    // future, and the cost is low enough that I won't bother with writing the
    // tiny-bit-more-efficient-half-the-time 8-byte alignment version for now.
    pheno_nm_ctv = round_up_pow2(pheno_nm_ctv, 2);
    for (; pidx < pmax; pidx++) {
      generate_cc_perm1(pheno_nm_ct, case_ct, tot_quotient, totq_magic, totq_preshift, totq_postshift, totq_incr, &(perm_vecs[pidx * pheno_nm_ctv]), sfmtp);
    }
  }
  THREAD_RETURN;
}

THREAD_RET_TYPE generate_cc_cluster_perms_thread(void* arg) {
  intptr_t tidx = (intptr_t)arg;
  uint32_t pheno_nm_ct = g_perm_pheno_nm_ct;
  uintptr_t* __restrict__ perm_vecs = g_perm_vecs;
  sfmt_t* __restrict__ sfmtp = g_sfmtp_arr[tidx];
  uintptr_t pheno_nm_ctv = BITCT_TO_WORDCT(pheno_nm_ct);
  uint32_t pidx = (((uint64_t)tidx) * g_perm_vec_ct) / g_perm_generation_thread_ct;
  uint32_t pmax = (((uint64_t)tidx + 1) * g_perm_vec_ct) / g_perm_generation_thread_ct;
  uint32_t cluster_ct = g_perm_cluster_ct;
  uint32_t* cluster_map = g_perm_cluster_map;
  uint32_t* cluster_starts = g_perm_cluster_starts;
  uint32_t* cluster_case_cts = g_perm_cluster_case_cts;
  uintptr_t* perm_cluster_cc_preimage = g_perm_cluster_cc_preimage;
  uint32_t* tot_quotients = g_perm_tot_quotients;
  uint64_t* totq_magics = g_perm_totq_magics;
  uint32_t* totq_preshifts = g_perm_totq_preshifts;
  uint32_t* totq_postshifts = g_perm_totq_postshifts;
  uint32_t* totq_incrs = g_perm_totq_incrs;
  if (!g_perm_is_1bit) {
    pheno_nm_ctv *= 2;
    for (; pidx < pmax; pidx++) {
      generate_cc_cluster_perm_vec(pheno_nm_ct, perm_cluster_cc_preimage, cluster_ct, cluster_map, cluster_starts, cluster_case_cts, tot_quotients, totq_magics, totq_preshifts, totq_postshifts, totq_incrs, &(perm_vecs[pidx * pheno_nm_ctv]), sfmtp);
    }
  } else {
    pheno_nm_ctv = round_up_pow2(pheno_nm_ctv, 2);
    for (; pidx < pmax; pidx++) {
      generate_cc_cluster_perm1(pheno_nm_ct, perm_cluster_cc_preimage, cluster_ct, cluster_map, cluster_starts, cluster_case_cts, tot_quotients, totq_magics, totq_preshifts, totq_postshifts, totq_incrs, &(perm_vecs[pidx * pheno_nm_ctv]), sfmtp);
    }
  }
  THREAD_RETURN;
}

THREAD_RET_TYPE generate_qt_perms_smajor_thread(void* arg) {
  // Used by QT --assoc and --make-perm-pheno.
  //
  // Takes an array of phenotype values in g_perm_pheno_d2 of length
  // g_perm_pheno_nm_ct, and populates g_perm_vecstd[] with permutations of
  // those values.  Also requires g_sfmtp_arr[] and
  // g_perm_generation_thread_ct to be initialized.
  //
  // g_perm_vecstd is sample-major.  The nth permutation is stored across
  //   g_perm_vecstd[n]
  //   g_perm_vecstd[n + perm_vec_ctcl8m]
  //   g_perm_vecstd[n + 2 * perm_vec_ctcl8m]
  //   ...
  uintptr_t tidx = (uintptr_t)arg;
  uint32_t pheno_nm_ct = g_perm_pheno_nm_ct;
  uintptr_t perm_vec_ctcl8 = (g_perm_vec_ct + (CACHELINE_DBL - 1)) / CACHELINE_DBL;
  uintptr_t perm_vec_ctcl8m = perm_vec_ctcl8 * CACHELINE_DBL;
  double* pheno_d2 = g_perm_pheno_d2;
  sfmt_t* sfmtp = g_sfmtp_arr[tidx];
  uint32_t pmin = CACHELINE_DBL * ((((uint64_t)tidx) * perm_vec_ctcl8) / g_perm_generation_thread_ct);
  uint32_t pmax = CACHELINE_DBL * ((((uint64_t)tidx + 1) * perm_vec_ctcl8) / g_perm_generation_thread_ct);
  double* perm_vecstd = &(g_perm_vecstd[pmin]);
  uint32_t poffset = 0;
  uint32_t sample_idx = 1;
  uint32_t pdiff;
  uint32_t tot_quotient;
  uint32_t upper_bound;
  uint64_t totq_magic;
  uint32_t totq_preshift;
  uint32_t totq_postshift;
  uint32_t totq_incr;
  uint32_t urand;
  uint32_t uii;
  double* wptr;
  double* wptr2;
  double* wptr3;
  double cur_source;
  if (tidx + 1 == g_perm_generation_thread_ct) {
    pmax = g_perm_vec_ct;
  }
  pdiff = pmax - pmin;
  cur_source = *pheno_d2++;
  wptr = perm_vecstd;
  for (; poffset < pdiff; poffset++) {
    *wptr++ = cur_source;
  }
  for (; sample_idx < pheno_nm_ct; sample_idx++) {
    tot_quotient = 0x100000000LLU / (sample_idx + 1);
    upper_bound = (sample_idx + 1) * tot_quotient - 1;
    magic_num(tot_quotient, &totq_magic, &totq_preshift, &totq_postshift, &totq_incr);
    cur_source = *pheno_d2++;
    wptr = &(perm_vecstd[sample_idx * perm_vec_ctcl8m]);
    wptr2 = perm_vecstd;
    for (poffset = 0; poffset < pdiff; poffset++) {
      do {
	urand = sfmt_genrand_uint32(sfmtp);
      } while (urand > upper_bound);
      uii = (totq_magic * ((urand >> totq_preshift) + totq_incr)) >> totq_postshift;
      wptr3 = &(wptr2[uii * perm_vec_ctcl8m]);
      *wptr++ = *wptr3;
      *wptr3 = cur_source;
      wptr2++;
    }
  }
  THREAD_RETURN;
}

THREAD_RET_TYPE generate_qt_cluster_perms_smajor_thread(void* arg) {
  // Variant of generate_qt_perms_smajor_thread() which restricts permutations
  // to be within-cluster.
  // On top of the generate_qt_perms_smajor_thread requirements, this also
  // needs g_perm_cluster_ct, g_perm_cluster_map, g_perm_cluster_starts,
  // g_perm_qt_cluster_thread_wkspace, and g_perm_sample_to_cluster to be
  // initialized.
  uintptr_t tidx = (uintptr_t)arg;
  uint32_t pheno_nm_ct = g_perm_pheno_nm_ct;
  uintptr_t perm_vec_ctcl8 = (g_perm_vec_ct + (CACHELINE_DBL - 1)) / CACHELINE_DBL;
  uintptr_t perm_vec_ctcl8m = perm_vec_ctcl8 * CACHELINE_DBL;
  double* pheno_d2 = g_perm_pheno_d2;
  sfmt_t* sfmtp = g_sfmtp_arr[tidx];
  uint32_t pmin = CACHELINE_DBL * ((((uint64_t)tidx) * perm_vec_ctcl8) / g_perm_generation_thread_ct);
  uint32_t pmax = CACHELINE_DBL * ((((uint64_t)tidx + 1) * perm_vec_ctcl8) / g_perm_generation_thread_ct);
  double* perm_vecstd = &(g_perm_vecstd[pmin]);
  uint32_t cluster_ct = g_perm_cluster_ct;
  uint32_t cluster_ctcl = (cluster_ct + (CACHELINE_INT32 - 1)) / CACHELINE_INT32;
  uint32_t* cluster_map = g_perm_cluster_map;
  uint32_t* cluster_starts = g_perm_cluster_starts;
  uint32_t* in_cluster_positions = &(g_perm_qt_cluster_thread_wkspace[tidx * cluster_ctcl * CACHELINE_INT32]);
  uint32_t* sample_to_cluster = g_perm_sample_to_cluster;
  uint32_t poffset = 0;
  uint32_t sample_idx = 0;
  uint32_t* cur_map_start;
  uint32_t pdiff;
  uint32_t cluster_idx;
  uint32_t cur_in_cluster_pos;
  uint32_t tot_quotient;
  uint32_t upper_bound;
  uint64_t totq_magic;
  uint32_t totq_preshift;
  uint32_t totq_postshift;
  uint32_t totq_incr;
  uint32_t urand;
  uint32_t uii;
  double* wptr;
  double* wptr2;
  double* wptr3;
  double cur_source;
  if (tidx + 1 == g_perm_generation_thread_ct) {
    pmax = g_perm_vec_ct;
  }
  pdiff = pmax - pmin;
  fill_uint_zero(cluster_ct, in_cluster_positions);
  for (; sample_idx < pheno_nm_ct; sample_idx++) {
    cur_source = *pheno_d2++;
    cluster_idx = sample_to_cluster[sample_idx];
    if (cluster_idx == 0xffffffffU) {
      cur_in_cluster_pos = 0;
    } else {
      cur_in_cluster_pos = in_cluster_positions[cluster_idx];
      in_cluster_positions[cluster_idx] += 1;
    }
    wptr = &(perm_vecstd[sample_idx * perm_vec_ctcl8m]);
    if (!cur_in_cluster_pos) {
      for (poffset = 0; poffset < pdiff; poffset++) {
        *wptr++ = cur_source;
      }
    } else {
      cur_map_start = &(cluster_map[cluster_starts[cluster_idx]]);
      tot_quotient = 0x100000000LLU / (cur_in_cluster_pos + 1);
      upper_bound = (cur_in_cluster_pos + 1) * tot_quotient - 1;
      magic_num(tot_quotient, &totq_magic, &totq_preshift, &totq_postshift, &totq_incr);
      wptr2 = perm_vecstd;
      for (poffset = 0; poffset < pdiff; poffset++) {
	do {
	  urand = sfmt_genrand_uint32(sfmtp);
	} while (urand > upper_bound);
	uii = (totq_magic * ((urand >> totq_preshift) + totq_incr)) >> totq_postshift;
	wptr3 = &(wptr2[cur_map_start[uii] * perm_vec_ctcl8m]);
	*wptr++ = *wptr3;
	*wptr3 = cur_source;
	wptr2++;
      }
    }
  }
  THREAD_RETURN;
}

THREAD_RET_TYPE generate_qt_perms_pmajor_thread(void* arg) {
  // Used by --linear.  Requires g_perm_pheno_nm_ct, g_perm_pheno_d2,
  // g_sfmtp_arr, g_perm_generation_thread_ct, and g_perm_vec_ct to be
  // initialized, and space must be allocated for g_perm_pmajor.  The nth
  // permutation (0-based) is stored in g_perm_pmajor indices
  //   [n * sample_valid_ct] to [(n + 1) * sample_valid_ct - 1]
  // inclusive.
  uintptr_t tidx = (uintptr_t)arg;
  uint32_t sample_valid_ct = g_perm_pheno_nm_ct;
  uintptr_t perm_vec_ctcl = (g_perm_vec_ct + (CACHELINE_INT32 - 1)) / CACHELINE_INT32;
  sfmt_t* sfmtp = g_sfmtp_arr[tidx];
  uintptr_t pmin = CACHELINE_INT32 * ((((uint64_t)tidx) * perm_vec_ctcl) / g_perm_generation_thread_ct);
  uintptr_t pmax = CACHELINE_INT32 * ((((uint64_t)tidx + 1) * perm_vec_ctcl) / g_perm_generation_thread_ct);
  double* perm_pmajor = &(g_perm_pmajor[pmin * sample_valid_ct]);
  double* pheno_d2 = g_perm_pheno_d2;
  uint32_t* precomputed_mods = g_perm_precomputed_mods;
  uint32_t* lbound_ptr;
  double* pheno_ptr;
  uint32_t poffset;
  uint32_t pdiff;
  uint32_t sample_idx;
  uint32_t urand;
  uint32_t lbound;
  if (tidx + 1 == g_perm_generation_thread_ct) {
    pmax = g_perm_vec_ct;
  }
  pdiff = pmax - pmin;
  for (poffset = 0; poffset < pdiff; poffset++) {
    lbound_ptr = precomputed_mods;
    pheno_ptr = pheno_d2;
    perm_pmajor[0] = *pheno_ptr++;
    for (sample_idx = 1; sample_idx < sample_valid_ct; sample_idx++) {
      lbound = *lbound_ptr++;
      do {
        urand = sfmt_genrand_uint32(sfmtp);
      } while (urand < lbound);
      // er, this modulus operation is slow.  but doesn't seem to be worthwhile
      // to use magic numbers here.
      urand %= sample_idx + 1;
      perm_pmajor[sample_idx] = perm_pmajor[urand];
      perm_pmajor[urand] = *pheno_ptr++;
    }
    perm_pmajor = &(perm_pmajor[sample_valid_ct]);
  }
  THREAD_RETURN;
}

THREAD_RET_TYPE generate_qt_cluster_perms_pmajor_thread(void* arg) {
  // On top of the linear_gen_perms_thread requirements, this also needs
  // g_perm_cluster_ct, g_perm_cluster_map, g_perm_cluster_starts,
  // g_perm_qt_cluster_thread_wkspace, and g_perm_sample_to_cluster to be
  // initialized.
  uintptr_t tidx = (uintptr_t)arg;
  uint32_t sample_valid_ct = g_perm_pheno_nm_ct;
  uintptr_t perm_vec_ctcl = (g_perm_vec_ct + (CACHELINE_INT32 - 1)) / CACHELINE_INT32;
  sfmt_t* sfmtp = g_sfmtp_arr[tidx];
  uintptr_t pmin = CACHELINE_INT32 * ((((uint64_t)tidx) * perm_vec_ctcl) / g_perm_generation_thread_ct);
  uintptr_t pmax = CACHELINE_INT32 * ((((uint64_t)tidx + 1) * perm_vec_ctcl) / g_perm_generation_thread_ct);
  double* perm_pmajor = &(g_perm_pmajor[pmin * sample_valid_ct]);
  double* pheno_d2 = g_perm_pheno_d2;
  uint32_t* precomputed_mods = &(g_perm_precomputed_mods[-1]);
  uint32_t cluster_ct = g_perm_cluster_ct;
  uint32_t cluster_ctcl = (cluster_ct + (CACHELINE_INT32 - 1)) / CACHELINE_INT32;
  uint32_t* cluster_map = g_perm_cluster_map;
  uint32_t* cluster_starts = g_perm_cluster_starts;
  uint32_t* in_cluster_positions = &(g_perm_qt_cluster_thread_wkspace[tidx * cluster_ctcl * CACHELINE_INT32]);
  uint32_t* sample_to_cluster = g_perm_sample_to_cluster;
  double* pheno_ptr;
  uint32_t poffset;
  uint32_t pdiff;
  uint32_t cluster_idx;
  uint32_t cur_in_cluster_pos;
  uint32_t sample_idx;
  uint32_t urand;
  uint32_t lbound;
  uint32_t uii;
  if (tidx + 1 == g_perm_generation_thread_ct) {
    pmax = g_perm_vec_ct;
  }
  pdiff = pmax - pmin;
  for (poffset = 0; poffset < pdiff; poffset++) {
    fill_uint_zero(cluster_ct, in_cluster_positions);
    pheno_ptr = pheno_d2;
    for (sample_idx = 0; sample_idx < sample_valid_ct; sample_idx++) {
      cluster_idx = sample_to_cluster[sample_idx];
      if (cluster_idx == 0xffffffffU) {
	cur_in_cluster_pos = 0;
      } else {
	cur_in_cluster_pos = in_cluster_positions[cluster_idx];
	in_cluster_positions[cluster_idx] += 1;
      }
      if (!cur_in_cluster_pos) {
        perm_pmajor[sample_idx] = *pheno_ptr++;
      } else {
        lbound = precomputed_mods[cur_in_cluster_pos];
        do {
	  urand = sfmt_genrand_uint32(sfmtp);
	} while (urand < lbound);
	urand %= (cur_in_cluster_pos + 1);
	uii = cluster_map[cluster_starts[cluster_idx] + urand];
        perm_pmajor[sample_idx] = perm_pmajor[uii];
	perm_pmajor[uii] = *pheno_ptr++;
      }
    }
    perm_pmajor = &(perm_pmajor[sample_valid_ct]);
  }
  THREAD_RETURN;
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
  uintptr_t pheno_nm_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(pheno_nm_ct);
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
	  fill_uint_zero(4, wbuf);
	  wshift = 0;
	}
	wbptr = wbuf;
      }
      *wbptr |= ((pvptr[perm_idx * pheno_nm_ctv2] >> rshift) & 1) << wshift;
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
      wval |= ((pvptr[perm_idx * pheno_nm_ctv2] >> rshift) & 1) << wshift;
    } while (++perm_idx < perm_vec_ct);
    *perm_vecst++ = wval;
#endif
  }
}

void transpose_perm1s(uintptr_t* perm_vecs, uint32_t perm_vec_ct, uint32_t pheno_nm_ct, uint32_t* perm_vecst) {
  uintptr_t sample_idx = 0;
  uintptr_t pheno_nm_ctv = BITCT_TO_ALIGNED_WORDCT(pheno_nm_ct);
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
	  fill_uint_zero(2, wbuf);
	  wshift = 0;
	}
	wbptr = wbuf;
      }
      *wbptr |= ((pvptr[perm_idx * pheno_nm_ctv] >> rshift) & 1) << wshift;
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
      wval |= ((pvptr[perm_idx * pheno_nm_ctv] >> rshift) & 1) << wshift;
    } while (++perm_idx < perm_vec_ct);
    *perm_vecst++ = wval;
#endif
  }
}


int32_t make_perm_pheno(pthread_t* threads, char* outname, char* outname_end, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, char* sample_ids, uintptr_t max_sample_id_len, uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, uint32_t pheno_nm_ct, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, char* output_missing_pheno, uint32_t permphe_ct) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t pheno_nm_ctl = BITCT_TO_WORDCT(pheno_nm_ct);
  uintptr_t pheno_nm_ctv = round_up_pow2(pheno_nm_ctl, VEC_WORDS);
  uintptr_t perm_vec_ctcl8m = 0;
  char* writebuf = nullptr;
  int32_t retval = 0;
  uintptr_t* ulptr;
  double* dptr;
  char* wptr;
  uintptr_t sample_uidx;
  uintptr_t sample_idx;
  uintptr_t perm_idx;
  uintptr_t ulii;
  uint32_t sample_nmidx;
  uint32_t rshift;
  if (!pheno_nm_ct) {
    logerrprint("Error: --make-perm-pheno requires phenotype data.\n");
    goto make_perm_pheno_ret_INVALID_CMDLINE;
  }
  g_perm_generation_thread_ct = MINV(g_thread_ct, permphe_ct);
  if (bigstack_init_sfmtp(g_perm_generation_thread_ct)) {
    goto make_perm_pheno_ret_NOMEM;
  }
  g_perm_pheno_nm_ct = pheno_nm_ct;
  g_perm_vec_ct = permphe_ct;
  ulii = 0;
  if (pheno_c) {
    g_perm_is_1bit = 1;
    g_perm_case_ct = popcount_longs(pheno_c, unfiltered_sample_ctl);
    // could seamlessly support multipass by using different permutation logic,
    // but pointless in practice; better to just generate multiple files
    if (bigstack_alloc_ul(permphe_ct * pheno_nm_ctv, &g_perm_vecs)) {
      goto make_perm_pheno_ret_NOMEM;
    }
    if (cluster_starts) {
      // most similar to testmiss()
      retval = cluster_include_and_reindex(unfiltered_sample_ct, pheno_nm, 1, pheno_c, pheno_nm_ct, 1, cluster_ct, cluster_map, cluster_starts, &g_perm_cluster_ct, &g_perm_cluster_map, &g_perm_cluster_starts, &g_perm_cluster_case_cts, &g_perm_cluster_cc_preimage);
      if (retval) {
	goto make_perm_pheno_ret_1;
      }
      if (!g_perm_cluster_ct) {
        logerrprint("Error: Degenerate --make-perm-pheno invocation (no size 2+ clusters).\n");
        goto make_perm_pheno_ret_INVALID_CMDLINE;
      }
      retval = cluster_alloc_and_populate_magic_nums(g_perm_cluster_ct, g_perm_cluster_map, g_perm_cluster_starts, &g_perm_tot_quotients, &g_perm_totq_magics, &g_perm_totq_preshifts, &g_perm_totq_postshifts, &g_perm_totq_incrs);
      if (retval) {
        goto make_perm_pheno_ret_1;
      }
      // not actually much of a point to multithreading since this is I/O
      // bound, but what the hell, the permutation generators already support
      // it
      if (spawn_threads(threads, &generate_cc_cluster_perms_thread, g_perm_generation_thread_ct)) {
	goto make_perm_pheno_ret_THREAD_CREATE_FAIL;
      }
      generate_cc_cluster_perms_thread((void*)ulii);
    } else {
      g_perm_cluster_starts = nullptr;
      g_perm_tot_quotient = 0x100000000LLU / pheno_nm_ct;
      magic_num(g_perm_tot_quotient, &g_perm_totq_magic, &g_perm_totq_preshift, &g_perm_totq_postshift, &g_perm_totq_incr);
      if (spawn_threads(threads, &generate_cc_perms_thread, g_perm_generation_thread_ct)) {
	goto make_perm_pheno_ret_THREAD_CREATE_FAIL;
      }
      generate_cc_perms_thread((void*)ulii);
    }
  } else {
    g_perm_pheno_d2 = (double*)alloc_and_init_collapsed_arr_incl((char*)pheno_d, sizeof(double), unfiltered_sample_ct, pheno_nm, pheno_nm_ct, 1);
    if (!g_perm_pheno_d2) {
      goto make_perm_pheno_ret_NOMEM;
    }
    perm_vec_ctcl8m = round_up_pow2(permphe_ct, CACHELINE_DBL);
    if (bigstack_alloc_d(perm_vec_ctcl8m * pheno_nm_ct, &g_perm_vecstd)) {
      goto make_perm_pheno_ret_NOMEM;
    }
    if (cluster_starts) {
      retval = cluster_include_and_reindex(unfiltered_sample_ct, pheno_nm, 1, nullptr, pheno_nm_ct, 0, cluster_ct, cluster_map, cluster_starts, &g_perm_cluster_ct, &g_perm_cluster_map, &g_perm_cluster_starts, nullptr, nullptr);
      if (retval) {
	goto make_perm_pheno_ret_1;
      }
      if (!g_perm_cluster_ct) {
        logerrprint("Error: Degenerate --make-perm-pheno invocation (no size 2+ clusters).\n");
        goto make_perm_pheno_ret_INVALID_CMDLINE;
      }
      if (bigstack_alloc_ui(pheno_nm_ct, &g_perm_sample_to_cluster) ||
          bigstack_alloc_ui(g_perm_generation_thread_ct * round_up_pow2(g_perm_cluster_ct, CACHELINE_INT32), &g_perm_qt_cluster_thread_wkspace)) {
	goto make_perm_pheno_ret_NOMEM;
      }
      fill_unfiltered_sample_to_cluster(pheno_nm_ct, g_perm_cluster_ct, g_perm_cluster_map, g_perm_cluster_starts, g_perm_sample_to_cluster);
      if (spawn_threads(threads, &generate_qt_cluster_perms_smajor_thread, g_perm_generation_thread_ct)) {
	goto make_perm_pheno_ret_THREAD_CREATE_FAIL;
      }
      generate_qt_cluster_perms_smajor_thread((void*)ulii);
    } else {
      if (spawn_threads(threads, &generate_qt_perms_smajor_thread, g_perm_generation_thread_ct)) {
	goto make_perm_pheno_ret_THREAD_CREATE_FAIL;
      }
      generate_qt_perms_smajor_thread((void*)ulii);
    }
    if (bigstack_alloc_c(permphe_ct * 16LU, &writebuf)) {
      goto make_perm_pheno_ret_NOMEM;
    }
  }
  join_threads(threads, g_perm_generation_thread_ct);
  memcpy(outname_end, ".pphe", 6);
  if (fopen_checked(outname, "w", &outfile)) {
    goto make_perm_pheno_ret_OPEN_FAIL;
  }
  sample_nmidx = 0;
  for (sample_uidx = 0, sample_idx = 0; sample_idx < sample_ct; sample_uidx++, sample_idx++) {
    next_unset_ul_unsafe_ck(sample_exclude, &sample_uidx);
    fputs(&(sample_ids[sample_uidx * max_sample_id_len]), outfile);
    if (!IS_SET(pheno_nm, sample_uidx)) {
      for (perm_idx = 0; perm_idx < permphe_ct; perm_idx++) {
	putc_unlocked('\t', outfile);
	fputs(output_missing_pheno, outfile);
      }
    } else if (pheno_c) {
      ulptr = &(g_perm_vecs[sample_nmidx / BITCT]);
      rshift = sample_nmidx % BITCT;
      for (perm_idx = 0; perm_idx < permphe_ct; perm_idx++) {
	putc_unlocked('\t', outfile);
        putc_unlocked('1' + ((ulptr[perm_idx * pheno_nm_ctv] >> rshift) & 1), outfile);
      }
      sample_nmidx++;
    } else {
      wptr = writebuf;
      dptr = &(g_perm_vecstd[sample_nmidx * perm_vec_ctcl8m]);
      for (perm_idx = 0; perm_idx < permphe_ct; perm_idx++) {
	*wptr++ = '\t';
        wptr = dtoa_g(*dptr++, wptr);
      }
      if (fwrite_checked(writebuf, wptr - writebuf, outfile)) {
	goto make_perm_pheno_ret_WRITE_FAIL;
      }
      sample_nmidx++;
    }
    if (putc_checked('\n', outfile)) {
      goto make_perm_pheno_ret_WRITE_FAIL;
    }
  }
  if (fclose_null(&outfile)) {
    goto make_perm_pheno_ret_WRITE_FAIL;
  }
  LOGPRINTFWW("--make-perm-pheno: Permuted phenotypes written to %s .\n", outname);
  while (0) {
  make_perm_pheno_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  make_perm_pheno_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  make_perm_pheno_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  make_perm_pheno_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  make_perm_pheno_ret_THREAD_CREATE_FAIL:
    retval = RET_THREAD_CREATE_FAIL;
    break;
  }
 make_perm_pheno_ret_1:
  bigstack_reset(bigstack_mark);
  fclose_cond(outfile);
  return retval;
}
