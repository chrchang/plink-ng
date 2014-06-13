#include "plink_common.h"

#include "plink_assoc.h"
#include "plink_cluster.h"
#include "plink_matrix.h"
#include "plink_set.h"
#include "plink_stats.h"

// currently assumed to be no larger than MODEL_BLOCKSIZE
#define GLM_BLOCKSIZE 512

// multithread globals
static double* g_orig_chisq;
static double* g_mperm_save_all;

// A separated-low-and-high-bit format was tried, and found to not really be
// any better than the usual PLINK 2-bit format.
static uintptr_t* g_loadbuf;

static uintptr_t* g_perm_vecs;

static double* g_pheno_d2;
#ifndef NOLAPACK
static double g_pheno_sum;
static double g_pheno_ssq;
#endif

// permutation-major instead of individual-major order for --linear (PERMORY
// speedups do not apply)
static double* g_perm_pmajor;
static uint32_t* g_precomputed_mods; // g_precomputed_mods[n] = 2^32 mod (n-2)

static uint32_t* g_nm_cts;

// This is *twice* the number of successes, because PLINK 1.07 counts tie as
// 0.5.  (Actually, it randomizes instead of deterministically adding 0.5; this
// randomization just adds noise so we don't replicate it.)
static uint32_t* g_perm_2success_ct;
static uint32_t* g_perm_attempt_ct;
static double* g_maxt_extreme_stat;
static double* g_maxt_thread_results;

// to avoid pathological multithreading issues, this is not a bitset
static unsigned char* g_perm_adapt_stop;

static uint32_t g_adapt_m_table[MODEL_BLOCKSIZE];
static uint32_t g_assoc_thread_ct;
static uintptr_t g_perm_vec_ct;
static uint32_t g_block_diff;
static uint32_t g_perms_done;
static uint32_t g_first_adapt_check;
static uint32_t g_pheno_nm_ct;
static uint32_t g_case_ct;
static double g_adaptive_intercept;
static double g_adaptive_slope;
static double g_aperm_alpha;
static double g_adaptive_ci_zt;
static uint32_t g_is_x;
static uint32_t g_is_y;
static uint32_t g_is_haploid;

static uint32_t g_tot_quotient;
static uint64_t g_totq_magic;
static uint32_t g_totq_preshift;
static uint32_t g_totq_postshift;
static uint32_t g_totq_incr;

static uint32_t g_cluster_ct;
static uint32_t* g_cluster_map;
static uint32_t* g_cluster_starts;
static uint32_t* g_cluster_case_cts;

// per-cluster magic number sets
static uintptr_t* g_cluster_cc_perm_preimage;
static uint32_t* g_tot_quotients;
static uint64_t* g_totq_magics;
static uint32_t* g_totq_preshifts;
static uint32_t* g_totq_postshifts;
static uint32_t* g_totq_incrs;

static uint32_t* g_indiv_to_cluster;
static uint32_t* g_qassoc_cluster_thread_wkspace;

THREAD_RET_TYPE logistic_gen_perms_thread(void* arg) {
  // just a clone of model_assoc_gen_perms_thread()
  intptr_t tidx = (intptr_t)arg;
  uint32_t pheno_nm_ct = g_pheno_nm_ct;
  uint32_t case_ct = g_case_ct;
  uint32_t tot_quotient = g_tot_quotient;
  uint64_t totq_magic = g_totq_magic;
  uint32_t totq_preshift = g_totq_preshift;
  uint32_t totq_postshift = g_totq_postshift;
  uint32_t totq_incr = g_totq_incr;
  uintptr_t* __restrict__ perm_vecs = g_perm_vecs;
  sfmt_t* __restrict__ sfmtp = g_sfmtp_arr[tidx];
  uintptr_t pheno_nm_ctl2 = 2 * ((pheno_nm_ct + (BITCT - 1)) / BITCT);
  uint32_t pidx = (((uint64_t)tidx) * g_perm_vec_ct) / g_assoc_thread_ct;
  uint32_t pmax = (((uint64_t)tidx + 1) * g_perm_vec_ct) / g_assoc_thread_ct;
  for (; pidx < pmax; pidx++) {
    generate_cc_perm_vec(pheno_nm_ct, case_ct, tot_quotient, totq_magic, totq_preshift, totq_postshift, totq_incr, &(perm_vecs[pidx * pheno_nm_ctl2]), sfmtp);
  }
  THREAD_RETURN;
}

THREAD_RET_TYPE logistic_gen_cluster_perms_thread(void* arg) {
  intptr_t tidx = (intptr_t)arg;
  uint32_t pheno_nm_ct = g_pheno_nm_ct;
  uintptr_t* __restrict__ perm_vecs = g_perm_vecs;
  sfmt_t* __restrict__ sfmtp = g_sfmtp_arr[tidx];
  uintptr_t pheno_nm_ctl2 = 2 * ((pheno_nm_ct + (BITCT - 1)) / BITCT);
  uint32_t pidx = (((uint64_t)tidx) * g_perm_vec_ct) / g_assoc_thread_ct;
  uint32_t pmax = (((uint64_t)tidx + 1) * g_perm_vec_ct) / g_assoc_thread_ct;
  uint32_t cluster_ct = g_cluster_ct;
  uint32_t* cluster_map = g_cluster_map;
  uint32_t* cluster_starts = g_cluster_starts;
  uint32_t* cluster_case_cts = g_cluster_case_cts;
  uintptr_t* cluster_cc_perm_preimage = g_cluster_cc_perm_preimage;
  uint32_t* tot_quotients = g_tot_quotients;
  uint64_t* totq_magics = g_totq_magics;
  uint32_t* totq_preshifts = g_totq_preshifts;
  uint32_t* totq_postshifts = g_totq_postshifts;
  uint32_t* totq_incrs = g_totq_incrs;
  for (; pidx < pmax; pidx++) {
    generate_cc_cluster_perm_vec(pheno_nm_ct, cluster_cc_perm_preimage, cluster_ct, cluster_map, cluster_starts, cluster_case_cts, tot_quotients, totq_magics, totq_preshifts, totq_postshifts, totq_incrs, &(perm_vecs[pidx * pheno_nm_ctl2]), sfmtp);
  }
  THREAD_RETURN;
}

THREAD_RET_TYPE linear_gen_perms_thread(void* arg) {
  // Used by --linear.  Requires g_pheno_nm_ct, g_pheno_d2, g_sfmtp_arr,
  // g_assoc_thread_ct, and g_perm_vec_ct to be initialized, and space must be
  // allocated for g_perm_pmajor.  The nth permutation (0-based) is stored in
  // g_perm_pmajor indices
  //   [n * indiv_valid_ct] to [(n + 1) * indiv_valid_ct - 1]
  // inclusive.
  intptr_t tidx = (intptr_t)arg;
  uint32_t indiv_valid_ct = g_pheno_nm_ct;
  uintptr_t perm_vec_ctcl = (g_perm_vec_ct + (CACHELINE_INT32 - 1)) / CACHELINE_INT32;
  sfmt_t* sfmtp = g_sfmtp_arr[tidx];
  uintptr_t pmin = CACHELINE_INT32 * ((((uint64_t)tidx) * perm_vec_ctcl) / g_assoc_thread_ct);
  uintptr_t pmax = CACHELINE_INT32 * ((((uint64_t)tidx + 1) * perm_vec_ctcl) / g_assoc_thread_ct);
  double* perm_pmajor = &(g_perm_pmajor[pmin * indiv_valid_ct]);
  double* pheno_d2 = g_pheno_d2;
  uint32_t* precomputed_mods = g_precomputed_mods;
  uint32_t* lbound_ptr;
  double* pheno_ptr;
  uint32_t poffset;
  uint32_t pdiff;
  uint32_t indiv_idx;
  uint32_t urand;
  uint32_t lbound;
  if (((uintptr_t)tidx) + 1 == g_assoc_thread_ct) {
    pmax = g_perm_vec_ct;
  }
  pdiff = pmax - pmin;
  for (poffset = 0; poffset < pdiff; poffset++) {
    lbound_ptr = precomputed_mods;
    pheno_ptr = pheno_d2;
    perm_pmajor[0] = *pheno_ptr++;
    for (indiv_idx = 1; indiv_idx < indiv_valid_ct; indiv_idx++) {
      lbound = *lbound_ptr++;
      do {
        urand = sfmt_genrand_uint32(sfmtp);
      } while (urand < lbound);
      // er, this modulus operation is slow.  but doesn't seem to be worthwhile
      // to use magic numbers here.
      urand %= indiv_idx + 1;
      perm_pmajor[indiv_idx] = perm_pmajor[urand];
      perm_pmajor[urand] = *pheno_ptr++;
    }
    perm_pmajor = &(perm_pmajor[indiv_valid_ct]);
  }
  THREAD_RETURN;
}

THREAD_RET_TYPE linear_gen_cluster_perms_thread(void* arg) {
  // On top of the linear_gen_perms_thread requirements, this also needs
  // g_cluster_ct, g_cluster_map, g_cluster_starts,
  // g_qassoc_cluster_thread_wkspace, and g_indiv_to_cluster to be initialized.
  intptr_t tidx = (intptr_t)arg;
  uint32_t indiv_valid_ct = g_pheno_nm_ct;
  uintptr_t perm_vec_ctcl = (g_perm_vec_ct + (CACHELINE_INT32 - 1)) / CACHELINE_INT32;
  sfmt_t* sfmtp = g_sfmtp_arr[tidx];
  uintptr_t pmin = CACHELINE_INT32 * ((((uint64_t)tidx) * perm_vec_ctcl) / g_assoc_thread_ct);
  uintptr_t pmax = CACHELINE_INT32 * ((((uint64_t)tidx + 1) * perm_vec_ctcl) / g_assoc_thread_ct);
  double* perm_pmajor = &(g_perm_pmajor[pmin * indiv_valid_ct]);
  double* pheno_d2 = g_pheno_d2;
  uint32_t* precomputed_mods = &(g_precomputed_mods[-1]);
  uint32_t cluster_ct = g_cluster_ct;
  uint32_t cluster_ctcl = (cluster_ct + (CACHELINE_INT32 - 1)) / CACHELINE_INT32;
  uint32_t* cluster_map = g_cluster_map;
  uint32_t* cluster_starts = g_cluster_starts;
  uint32_t* in_cluster_positions = &(g_qassoc_cluster_thread_wkspace[tidx * cluster_ctcl * CACHELINE_INT32]);
  uint32_t* indiv_to_cluster = g_indiv_to_cluster;
  double* pheno_ptr;
  uint32_t poffset;
  uint32_t pdiff;
  uint32_t cluster_idx;
  uint32_t cur_in_cluster_pos;
  uint32_t indiv_idx;
  uint32_t urand;
  uint32_t lbound;
  uint32_t uii;
  if (((uintptr_t)tidx) + 1 == g_assoc_thread_ct) {
    pmax = g_perm_vec_ct;
  }
  pdiff = pmax - pmin;
  for (poffset = 0; poffset < pdiff; poffset++) {
    fill_uint_zero(in_cluster_positions, cluster_ct);
    pheno_ptr = pheno_d2;
    for (indiv_idx = 0; indiv_idx < indiv_valid_ct; indiv_idx++) {
      cluster_idx = indiv_to_cluster[indiv_idx];
      if (cluster_idx == 0xffffffffU) {
	cur_in_cluster_pos = 0;
      } else {
	cur_in_cluster_pos = in_cluster_positions[cluster_idx];
	in_cluster_positions[cluster_idx] += 1;
      }
      if (!cur_in_cluster_pos) {
        perm_pmajor[indiv_idx] = *pheno_ptr++;
      } else {
        lbound = precomputed_mods[cur_in_cluster_pos];
        do {
	  urand = sfmt_genrand_uint32(sfmtp);
	} while (urand < lbound);
	urand %= (cur_in_cluster_pos + 1);
	uii = cluster_map[cluster_starts[cluster_idx] + urand];
        perm_pmajor[indiv_idx] = perm_pmajor[uii];
	perm_pmajor[uii] = *pheno_ptr++;
      }
    }
    perm_pmajor = &(perm_pmajor[indiv_valid_ct]);
  }
  THREAD_RETURN;
}

uint32_t glm_init_load_mask(uintptr_t* indiv_exclude, uintptr_t* pheno_nm, uintptr_t* covar_nm, uint32_t indiv_ct, uintptr_t unfiltered_indiv_ctv2, uintptr_t** load_mask_ptr) {
  uint32_t indiv_uidx = 0;
  uintptr_t* load_mask;
  uint32_t indiv_idx;
  if (wkspace_alloc_ul_checked(load_mask_ptr, unfiltered_indiv_ctv2 * (sizeof(intptr_t) / 2))) {
    return 1;
  }
  load_mask = *load_mask_ptr;
  fill_ulong_zero(load_mask, unfiltered_indiv_ctv2 / 2);
  if (covar_nm) {
    for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_uidx++, indiv_idx++) {
      next_unset_unsafe_ck(indiv_exclude, &indiv_uidx);
      if (IS_SET(pheno_nm, indiv_uidx) & IS_SET(covar_nm, indiv_idx)) {
	SET_BIT(load_mask, indiv_uidx);
      }
    }
  } else {
    memcpy(load_mask, pheno_nm, unfiltered_indiv_ctv2 * (sizeof(intptr_t) / 2));
  }
  return 0;
}

int32_t glm_scan_conditions(char* condition_mname, char* condition_fname, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, Chrom_info* chrom_info_ptr, uint32_t hh_exists, uintptr_t* loadbuf_raw, FILE* bedfile, uintptr_t bed_offset, uintptr_t unfiltered_indiv_ct, uintptr_t* sex_male, uintptr_t* load_mask, uintptr_t* indiv_valid_ct_ptr, uintptr_t* condition_ct_ptr, uint32_t** condition_uidxs_ptr, uintptr_t* indiv_raw_include2, uintptr_t* indiv_raw_male_include2) {
  // side effects: load_mask and indiv_valid_ct potentially updated,
  //   condition_ct should be changed, condition_uidxs should be malloc'd
  unsigned char* wkspace_mark = wkspace_base;
  FILE* condition_file = NULL;
  uint32_t* condition_uidxs = NULL;
  uintptr_t marker_ctl = (marker_ct + (BITCT - 1)) / BITCT;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t unfiltered_indiv_ctv2 = 2 * ((unfiltered_indiv_ct + (BITCT - 1)) / BITCT);
  uintptr_t indiv_valid_ct = *indiv_valid_ct_ptr;
  uintptr_t miss_ct = 0;
  uintptr_t condition_ct = 0;
  uintptr_t line_idx = 0;
  int32_t retval = 0;
#ifdef __LP64__
  __m128i* loadbuf_vptr;
  __m128i* loadbuf_mask_vptr;
  __m128i* loadbuf_vend;
  __m128i vii;
#else
  uintptr_t unfiltered_indiv_ctl2 = (unfiltered_indiv_ct + (BITCT2 - 1)) / BITCT2;
  uintptr_t* loadbuf_end;
#endif
  uintptr_t* loadbuf_ptr;
  uintptr_t* loadbuf_mask_ptr;
  char* sorted_ids;
  uint32_t* id_map;
  uintptr_t* already_seen;
  uintptr_t* loadbuf_mask_orig;
  uintptr_t* loadbuf_mask;
  uint32_t* marker_idx_to_uidx;
  uint32_t* condition_uidxs_tmp;
  char* bufptr;
  char* bufptr2;
  uintptr_t condition_ct_max;
  uintptr_t condition_idx;
  uintptr_t indiv_uidx_offset;
  uintptr_t ulii;
  uintptr_t uljj;
  uint32_t marker_uidx;
  uint32_t chrom_idx;
  uint32_t is_x;
  uint32_t is_y;
  int32_t ii;
  if (condition_mname) {
    ii = get_uidx_from_unsorted(condition_mname, marker_exclude, marker_ct, marker_ids, max_marker_id_len);
    if (ii == -1) {
      LOGPRINTFWW("Warning: --condition variant ID '%s' not found.\n", condition_mname);
      return 0;
    }
    condition_ct = 1;
    condition_uidxs = (uint32_t*)malloc(sizeof(int32_t));
    condition_uidxs[0] = (uint32_t)ii;
  } else {
    if (wkspace_alloc_c_checked(&sorted_ids, marker_ct * max_marker_id_len) ||
        wkspace_alloc_ui_checked(&id_map, marker_ct * sizeof(int32_t)) ||
        wkspace_alloc_ul_checked(&already_seen, marker_ctl * sizeof(intptr_t)) ||
        wkspace_alloc_ui_checked(&marker_idx_to_uidx, marker_ct * sizeof(int32_t))) {
      goto glm_scan_conditions_ret_NOMEM;
    }
    fill_idx_to_uidx(marker_exclude, unfiltered_marker_ct, marker_ct, marker_idx_to_uidx);
    fill_ulong_zero(already_seen, marker_ctl);
    retval = sort_item_ids_noalloc(sorted_ids, id_map, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, 0, 1, strcmp_deref);
    if (retval) {
      goto glm_scan_conditions_ret_1;
    }
    condition_uidxs_tmp = (uint32_t*)wkspace_base;
    if (wkspace_left > marker_ct * sizeof(int32_t)) {
      condition_ct_max = marker_ct;
    } else {
      condition_ct_max = wkspace_left / sizeof(int32_t);
    }
    if (fopen_checked(&condition_file, condition_fname, "r")) {
      goto glm_scan_conditions_ret_OPEN_FAIL;
    }
    tbuf[MAXLINELEN - 1] = ' ';
    while (fgets(tbuf, MAXLINELEN, condition_file)) {
      line_idx++;
      if (!tbuf[MAXLINELEN - 1]) {
        sprintf(logbuf, "Error: Line %" PRIuPTR " of --condition-list file is pathologically long.\n", line_idx);
        goto glm_scan_conditions_ret_INVALID_FORMAT_2;
      }
      bufptr = skip_initial_spaces(tbuf);
      while (!is_eoln_kns(*bufptr)) {
        bufptr2 = token_endnn(bufptr);
	ii = bsearch_str(bufptr, (uintptr_t)(bufptr2 - bufptr), sorted_ids, max_marker_id_len, marker_ct);
	if (ii == -1) {
	  miss_ct++;
	} else {
	  if (is_set(already_seen, ii)) {
	    LOGPRINTFWW("Error: Duplicate variant '%s' in --condition-list file.\n", bufptr);
	    goto glm_scan_conditions_ret_INVALID_FORMAT;
	  }
	  if (condition_ct == condition_ct_max) {
	    goto glm_scan_conditions_ret_NOMEM;
	  }
	  set_bit(already_seen, ii);
	  condition_uidxs_tmp[condition_ct++] = marker_idx_to_uidx[id_map[(uint32_t)ii]];
	}
        bufptr = skip_initial_spaces(bufptr2);
      }
    }
    if (!feof(condition_file)) {
      goto glm_scan_conditions_ret_READ_FAIL;
    }
    if (condition_ct) {
      condition_uidxs = (uint32_t*)malloc(condition_ct * sizeof(int32_t));
      memcpy(condition_uidxs, condition_uidxs_tmp, condition_ct * sizeof(int32_t));
    } else if (!miss_ct) {
      logprint("Warning: --condition-list file is empty.\n");
      goto glm_scan_conditions_ret_1;
    }
    if (miss_ct) {
      sprintf(logbuf, "--condition-list: %" PRIuPTR " of %" PRIuPTR " variant ID%s loaded from %s.\n", condition_ct, condition_ct + miss_ct, (condition_ct + miss_ct == 1)? "" : "s", condition_fname);
    } else {
      sprintf(logbuf, "--condition-list: %" PRIuPTR " variant ID%s loaded from %s.\n", condition_ct, (condition_ct == 1)? "" : "s", condition_fname);
    }
    logprintb();
  }
  if (condition_ct) {
    if (wkspace_alloc_ul_checked(&loadbuf_mask_orig, unfiltered_indiv_ctv2 * sizeof(intptr_t)) ||
        wkspace_alloc_ul_checked(&loadbuf_mask, unfiltered_indiv_ctv2 * sizeof(intptr_t))) {
      goto glm_scan_conditions_ret_NOMEM;
    }
    vec_include_init(unfiltered_indiv_ct, loadbuf_mask_orig, load_mask);
    memcpy(loadbuf_mask, loadbuf_mask_orig, unfiltered_indiv_ctv2 * sizeof(intptr_t));
#ifdef __LP64__
    loadbuf_vend = (__m128i*)(&(loadbuf_raw[unfiltered_indiv_ctv2]));
#else
    loadbuf_end = &(loadbuf_raw[unfiltered_indiv_ctl2]);
#endif
    for (condition_idx = 0; condition_idx < condition_ct; condition_idx++) {
      // scan for missing values to update indiv_valid_ct
      marker_uidx = condition_uidxs[condition_idx];
      if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	goto glm_scan_conditions_ret_READ_FAIL;
      }
      // don't use load_and_collapse since collapse bitmask not finalized
      if (fread(loadbuf_raw, 1, unfiltered_indiv_ct4, bedfile) < unfiltered_indiv_ct4) {
	goto glm_scan_conditions_ret_READ_FAIL;
      }
      chrom_idx = get_marker_chrom(chrom_info_ptr, marker_uidx);
      if (IS_SET(chrom_info_ptr->haploid_mask, chrom_idx)) {
	is_x = ((int32_t)chrom_idx == chrom_info_ptr->x_code);
	is_y = ((int32_t)chrom_idx == chrom_info_ptr->y_code);
	haploid_fix(hh_exists, indiv_raw_include2, indiv_raw_male_include2, unfiltered_indiv_ct, is_x, is_y, (unsigned char*)loadbuf_raw);
      }
      // clear loadbuf_mask bits where loadbuf is 01.
#ifdef __LP64__
      loadbuf_vptr = (__m128i*)loadbuf_raw;
      loadbuf_mask_vptr = (__m128i*)loadbuf_mask;
      do {
        vii = *loadbuf_vptr++;
        vii = _mm_andnot_si128(_mm_srli_epi64(vii, 1), vii);
        *loadbuf_mask_vptr = _mm_andnot_si128(vii, *loadbuf_mask_vptr);
        loadbuf_mask_vptr++;
      } while (loadbuf_vptr < loadbuf_vend);
#else
      loadbuf_ptr = loadbuf_raw;
      loadbuf_mask_ptr = loadbuf_mask;
      do {
        ulii = *loadbuf_ptr++;
        ulii = ((~ulii) >> 1) & ulii;
        *loadbuf_mask_ptr &= ~ulii;
	loadbuf_mask_ptr++;
      } while (loadbuf_ptr < loadbuf_end);
#endif
    }
    loadbuf_ptr = loadbuf_mask_orig;
    loadbuf_mask_ptr = loadbuf_mask;
    for (indiv_uidx_offset = 0; indiv_uidx_offset < unfiltered_indiv_ct; indiv_uidx_offset += BITCT2) {
      ulii = (*loadbuf_ptr++) & (~(*loadbuf_mask_ptr++));
      while (ulii) {
        uljj = CTZLU(ulii);
        clear_bit_ul(load_mask, indiv_uidx_offset + (uljj / 2));
        indiv_valid_ct--;
        ulii &= ulii - 1;
      }
    }
    *condition_uidxs_ptr = condition_uidxs;
  }
  *condition_ct_ptr = condition_ct;
  *indiv_valid_ct_ptr = indiv_valid_ct;
  while (0) {
  glm_scan_conditions_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  glm_scan_conditions_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  glm_scan_conditions_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  glm_scan_conditions_ret_INVALID_FORMAT_2:
    logprintb();
  glm_scan_conditions_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 glm_scan_conditions_ret_1:
  wkspace_reset(wkspace_mark);
  fclose_cond(condition_file);
  return retval;
}

uint32_t glm_loadbuf_to_doubles(uintptr_t* loadbuf_collapsed, uint32_t indiv_valid_ct, double* covar_row, double* geno_map, uintptr_t* cur_missing) {
  // ok for cur_missing to be NULL if there can't possibly be any missing
  uintptr_t* ulptr_end = &(loadbuf_collapsed[indiv_valid_ct / BITCT2]);
  uint32_t cur_missing_ct = 0;
  uint32_t indiv_idx = 0;
  uint32_t indiv_idx_stop = BITCT2;
  uintptr_t cur_word;
  uintptr_t cur_genotype;
  while (1) {
    while (loadbuf_collapsed < ulptr_end) {
      cur_word = *loadbuf_collapsed++;
      for (; indiv_idx < indiv_idx_stop; indiv_idx++, covar_row++, cur_word >>= 2) {
        cur_genotype = cur_word & 3;
        if (cur_genotype != 1) {
          *covar_row = geno_map[cur_genotype];
	} else {
          SET_BIT(cur_missing, indiv_idx);
	}
      }
      indiv_idx_stop += BITCT2;
    }
    if (indiv_idx == indiv_valid_ct) {
      return cur_missing_ct;
    }
    ulptr_end++;
    indiv_idx_stop = indiv_valid_ct;
  }
}

uint32_t glm_loadbuf_to_doubles_x(uintptr_t* loadbuf_collapsed, uintptr_t* sex_male_collapsed, uint32_t indiv_valid_ct, double* covar_row, double* geno_map, uintptr_t* cur_missing) {
  uintptr_t* ulptr_end = &(loadbuf_collapsed[indiv_valid_ct / BITCT2]);
  uint32_t cur_missing_ct = 0;
  uint32_t indiv_idx = 0;
  uint32_t indiv_idx_stop = BITCT2;
  uintptr_t cur_word;
  uintptr_t cur_genotype;
  while (1) {
    while (loadbuf_collapsed < ulptr_end) {
      cur_word = *loadbuf_collapsed++;
      for (; indiv_idx < indiv_idx_stop; indiv_idx++, covar_row++, cur_word >>= 2) {
        cur_genotype = cur_word & 3;
        if (cur_genotype != 1) {
          *covar_row = geno_map[cur_genotype + 4 * IS_SET(sex_male_collapsed, indiv_idx)];
	} else {
          SET_BIT(cur_missing, indiv_idx);
	}
      }
      indiv_idx_stop += BITCT2;
    }
    if (indiv_idx == indiv_valid_ct) {
      return cur_missing_ct;
    }
    ulptr_end++;
    indiv_idx_stop = indiv_valid_ct;
  }
}

uint32_t glm_loadbuf_to_floats(uintptr_t* loadbuf_collapsed, uint32_t indiv_valid_ct, float* covar_row, float* geno_map, uintptr_t* cur_missing) {
  // ok for cur_missing to be NULL if there can't possibly be any missing
  uintptr_t* ulptr_end = &(loadbuf_collapsed[indiv_valid_ct / BITCT2]);
  uint32_t cur_missing_ct = 0;
  uint32_t indiv_idx = 0;
  uint32_t indiv_idx_stop = BITCT2;
  uintptr_t cur_word;
  uintptr_t cur_genotype;
  while (1) {
    while (loadbuf_collapsed < ulptr_end) {
      cur_word = *loadbuf_collapsed++;
      for (; indiv_idx < indiv_idx_stop; indiv_idx++, covar_row++, cur_word >>= 2) {
        cur_genotype = cur_word & 3;
        if (cur_genotype != 1) {
          *covar_row = geno_map[cur_genotype];
	} else {
          SET_BIT(cur_missing, indiv_idx);
	}
      }
      indiv_idx_stop += BITCT2;
    }
    if (indiv_idx == indiv_valid_ct) {
      return cur_missing_ct;
    }
    ulptr_end++;
    indiv_idx_stop = indiv_valid_ct;
  }
}

uint32_t glm_loadbuf_to_floats_x(uintptr_t* loadbuf_collapsed, uintptr_t* sex_male_collapsed, uint32_t indiv_valid_ct, float* covar_row, float* geno_map, uintptr_t* cur_missing) {
  uintptr_t* ulptr_end = &(loadbuf_collapsed[indiv_valid_ct / BITCT2]);
  uint32_t cur_missing_ct = 0;
  uint32_t indiv_idx = 0;
  uint32_t indiv_idx_stop = BITCT2;
  uintptr_t cur_word;
  uintptr_t cur_genotype;
  while (1) {
    while (loadbuf_collapsed < ulptr_end) {
      cur_word = *loadbuf_collapsed++;
      for (; indiv_idx < indiv_idx_stop; indiv_idx++, covar_row++, cur_word >>= 2) {
        cur_genotype = cur_word & 3;
        if (cur_genotype != 1) {
          *covar_row = geno_map[cur_genotype + 4 * IS_SET(sex_male_collapsed, indiv_idx)];
	} else {
          SET_BIT(cur_missing, indiv_idx);
	}
      }
      indiv_idx_stop += BITCT2;
    }
    if (indiv_idx == indiv_valid_ct) {
      return cur_missing_ct;
    }
    ulptr_end++;
    indiv_idx_stop = indiv_valid_ct;
  }
}

int32_t glm_check_vif(double vif_thresh, uintptr_t param_ct, uintptr_t indiv_valid_ct, double* covars_collapsed, double* param_2d_buf, MATRIX_INVERT_BUF1_TYPE* mi_buf, double* param_2d_buf2) {
  __CLPK_integer dim = ((uint32_t)param_ct) - 1;
  uintptr_t param_ct_m1 = param_ct - 1;
  double indiv_ct_d = (double)((intptr_t)indiv_valid_ct);
  double indiv_ct_recip = 1.0 / indiv_ct_d;
  double indiv_ct_m1_recip = 1.0 / ((double)((intptr_t)(indiv_valid_ct - 1)));
  double* dptr;
  double* dptr2;
  double* dptr3;
  uintptr_t param_idx;
  uintptr_t param_idx2;
  uintptr_t indiv_idx;
  double dxx;
  double dyy;
  int32_t ii;
  for (param_idx = 1; param_idx < param_ct; param_idx++) {
    dyy = 0; // sum
    dptr = &(covars_collapsed[param_idx * indiv_valid_ct]);
    for (indiv_idx = 0; indiv_idx < indiv_valid_ct; indiv_idx++) {
      dxx = *dptr++;
      dyy += dxx;
    }
    param_2d_buf2[param_idx] = dyy * indiv_ct_recip;
  }
  for (param_idx = 1; param_idx < param_ct; param_idx++) {
    dptr = &(param_2d_buf[(param_idx - 1) * param_ct]);
    dyy = param_2d_buf2[param_idx];
    for (param_idx2 = param_idx; param_idx2 < param_ct; param_idx2++) {
      dxx = 0;
      dptr2 = &(covars_collapsed[param_idx * indiv_valid_ct]);
      dptr3 = &(covars_collapsed[param_idx2 * indiv_valid_ct]);
      for (indiv_idx = 0; indiv_idx < indiv_valid_ct; indiv_idx++) {
	dxx += (*dptr2++) * (*dptr3++);
      }
      dxx -= dyy * param_2d_buf2[param_idx2] * indiv_ct_d;
      *dptr++ = dxx * indiv_ct_m1_recip;
    }
  }
  for (param_idx = 0; param_idx < param_ct_m1; param_idx++) {
    dxx = param_2d_buf[param_idx * param_ct]; // diagonal element
    if ((dxx == 0) || (!realnum(dxx))) {
      return 1;
    }
    param_2d_buf2[param_idx] = 1.0 / sqrt(dxx); // now inverse sqrt of diagonal
  }
  for (param_idx = 1; param_idx < param_ct_m1; param_idx++) {
    dyy = param_2d_buf2[param_idx];
    dptr = &(param_2d_buf[param_idx * param_ct_m1]);
    dptr3 = param_2d_buf2;
    for (param_idx2 = 0; param_idx2 < param_idx; param_idx2++) {
      dptr2 = &(param_2d_buf[param_idx2 * param_ct_m1 + param_idx]);
      dxx = (*dptr2) * dyy * (*dptr3++);
      if (dxx > 0.999) {
	return 1;
      }
      *dptr2 = dxx;
      *dptr++ = dxx;
    }
  }
  for (param_idx = 0; param_idx < param_ct_m1; param_idx++) {
    param_2d_buf[param_idx * param_ct] = 1;
  }
  
  ii = invert_matrix(dim, param_2d_buf, mi_buf, param_2d_buf2);
  if (ii) {
    return ii;
  }
  for (param_idx = 0; param_idx < param_ct_m1; param_idx++) {
    if (param_2d_buf[param_idx * param_ct] > vif_thresh) {
      return 1;
    }
  }
  return 0;
}

uint32_t glm_linear_robust_cluster_covar(uintptr_t cur_batch_size, uintptr_t param_ct, uintptr_t indiv_valid_ct, uint32_t missing_ct, uintptr_t* loadbuf, uint32_t standard_beta, double pheno_sum_base, double pheno_ssq_base, double* covars_cov_major, double* covars_indiv_major, double* perm_pmajor, double* coef, double* param_2d_buf, MATRIX_INVERT_BUF1_TYPE* mi_buf, double* param_2d_buf2, uint32_t cluster_ct1, uint32_t* indiv_to_cluster1, double* cluster_param_buf, double* cluster_param_buf2, double* indiv_1d_buf, double* linear_results, uintptr_t constraint_ct, double* constraints_con_major, double* param_df_buf, double* param_df_buf2, double* df_df_buf, double* df_buf, uint32_t* perm_fail_ct_ptr, uintptr_t* perm_fails) {
  // See the second half of PLINK 1.07 linear.cpp fitLM(), and
  // validParameters().
  // Diagonals of the final covariance matrices (not including the intercept
  // element) are saved to linear_results[(perm_idx * (param_ct - 1))..
  // ((perm_idx + 1) * (param_ct - 1) - 1)].
  // If not all permutations yield valid results, bits in perm_fails[] are set.
  // A return value of 1 reports that ALL permutations failed.  In this case,
  // perm_fails is not necessarily updated.
  uintptr_t param_ct_p1 = param_ct + 1; // diagonals of param * param matrix
  uintptr_t param_ct_m1 = param_ct - 1;
  uintptr_t joint_test_requested = (constraints_con_major? 1 : 0);
  uintptr_t param_ctx = param_ct + joint_test_requested;
  uintptr_t param_ctx_m1 = param_ctx - 1;
  uintptr_t cur_word = 0;
  uint32_t cluster_ct1_p1 = cluster_ct1 + 1;
  uint32_t perm_fail_ct = 0;
  double* dptr;
  double* dptr2;
  double* dptr3;
  uintptr_t* ulptr;
  double* pheno_ptr;
  double partial;
  double min_sigma;
  double pheno_sum;
  double pheno_ssq;
  uintptr_t indiv_uidx;
  uintptr_t indiv_idx;
  uintptr_t param_idx;
  uintptr_t param_idx2;
  uintptr_t perm_idx;

  // +1 since "no cluster" is represented as -1
  // 32-bit since that -1 is actually 0xffffffffU
  uint32_t cluster_idx_p1;

  double dxx;
  double dyy;
  double dzz;
  fill_ulong_zero(perm_fails, ((cur_batch_size + (BITCT - 1)) / BITCT) * sizeof(intptr_t));
  col_major_matrix_multiply((uint32_t)param_ct, (uint32_t)param_ct, (uint32_t)indiv_valid_ct, covars_indiv_major, covars_cov_major, param_2d_buf);
  if (invert_matrix((uint32_t)param_ct, param_2d_buf, mi_buf, param_2d_buf2)) {
    return 1;
  }

  if (!cluster_ct1) {
    // only need to perform S[i][j] / sqrt(S[i][i] * S[j][j]) check once, since
    // it's independent of sigma

    // may as well cache square roots
    for (param_idx = 0; param_idx < param_ct; param_idx++) {
      param_2d_buf2[param_idx] = sqrt(param_2d_buf[param_idx * param_ct_p1]);
    }
    for (param_idx = 1; param_idx < param_ct; param_idx++) {
      dxx = 0.99999 * param_2d_buf2[param_idx];
      dptr = &(param_2d_buf[param_idx * param_ct]); // S[i][j]
      dptr2 = param_2d_buf2; // sqrt(S[j][j])
      for (param_idx2 = 0; param_idx2 < param_idx; param_idx2++) {
	// S[i][j] / sqrt(S[i][i] * S[j][j]) > 0.99999
	// -> S[i][j] > 0.99999 * sqrt(S[i][i]) * sqrt(S[j][j])
	if ((*dptr++) > dxx * (*dptr2++)) {
	  return 1;
	}
      }
    }
    // now determine min(S[i][i]), i >= 1
    min_sigma = param_2d_buf[param_ct_p1];
    for (param_idx = 2; param_idx < param_ct; param_idx++) {
      dxx = param_2d_buf[param_idx * param_ct_p1];
      if (min_sigma > dxx) {
	min_sigma = dxx;
      }
    }
    if (min_sigma <= 0) {
      return 1;
    }
    // S[i][i] * sigma < 1e-20 iff sigma < 1e-20 / S[i][i]
    min_sigma = 1e-20 / min_sigma;

    // now temporarily save sigmas in linear_results[perm_idx * param_ctx_m1]
    pheno_ptr = perm_pmajor;
    dyy = 1;
    dzz = 0;
    for (perm_idx = 0; perm_idx < cur_batch_size; perm_idx++) {
      // not coef[perm_idx * param_ct] due to how dgels() works
      dptr = &(coef[perm_idx * indiv_valid_ct]);

      dptr2 = covars_indiv_major;
      dxx = 0;
      if (!missing_ct) {
	for (indiv_idx = 0; indiv_idx < indiv_valid_ct + missing_ct; indiv_idx++) {
	  partial = 0;
	  dptr3 = dptr;
	  for (param_idx = 0; param_idx < param_ct; param_idx++) {
	    partial += (*dptr2++) * (*dptr3++);
	  }
	  partial -= *pheno_ptr++;
	  dxx += partial * partial;
	}
      } else {
	if (standard_beta) {
	  pheno_sum = pheno_sum_base;
	  pheno_ssq = pheno_ssq_base;
	  ulptr = loadbuf;
	  for (indiv_uidx = 0; indiv_uidx < indiv_valid_ct; indiv_uidx += BITCT2) {
	    cur_word = *ulptr++;
	    cur_word = cur_word & (~(cur_word >> 1)) & FIVEMASK;
	    while (cur_word) {
	      dyy = pheno_ptr[indiv_uidx + (CTZLU(cur_word) / 2)];
	      pheno_sum -= dyy;
	      pheno_ssq -= dyy * dyy;
	      cur_word &= cur_word - 1;
	    }
	  }
	  dzz = pheno_sum / ((double)((intptr_t)indiv_valid_ct));
          dyy = sqrt(((double)((intptr_t)(indiv_valid_ct - 1))) / (pheno_ssq - pheno_sum * dzz));
	}
	for (indiv_uidx = 0, indiv_idx = 0; indiv_idx < indiv_valid_ct; indiv_uidx++, indiv_idx++) {
	  partial = 0;
	  dptr3 = dptr;
	  for (param_idx = 0; param_idx < param_ct; param_idx++) {
	    partial += (*dptr2++) * (*dptr3++);
	  }
	  while (1) {
	    if (indiv_uidx % BITCT2) {
	      cur_word >>= 2;
	    } else {
	      cur_word = loadbuf[indiv_uidx / BITCT2];
	      cur_word = cur_word & (~(cur_word >> 1)) & FIVEMASK;
	    }
	    if (!(cur_word & 1)) {
	      break;
	    }
	    indiv_uidx++;
	  }
	  partial -= (pheno_ptr[indiv_uidx] - dzz) * dyy;
	  dxx += partial * partial;
	}
	pheno_ptr = &(pheno_ptr[indiv_valid_ct + missing_ct]);
      }
      linear_results[perm_idx * param_ctx_m1] = dxx;
    }
    dxx = 1.0 / ((double)((intptr_t)(indiv_valid_ct - param_ct)));
    if (!constraint_ct) {
      // only need diagonal of S
      dptr2 = param_2d_buf2;
      for (param_idx = 1; param_idx < param_ct; param_idx++) {
	*dptr2++ = param_2d_buf[param_idx * param_ct_p1]; // S0[i][i]
      }
    }
    dptr = linear_results;
    for (perm_idx = 0; perm_idx < cur_batch_size; perm_idx++) {
      dyy = dxx * (*dptr); // sigma = (previous sigma) / (nind-np)
      if (dyy < min_sigma) {
	dyy = 0;
	perm_fail_ct++;
	dptr = &(dptr[param_ctx_m1]);
	SET_BIT(perm_fails, perm_idx);
      } else {
	dptr2 = param_2d_buf2;
	if (!joint_test_requested) {
	  for (param_idx = 1; param_idx < param_ct; param_idx++) {
	    *dptr++ = (*dptr2++) * dyy; // S[i][i] = S0[i][i] * sigma
	  }
	} else if (constraint_ct) {
	  // need all of S for linear_hypothesis_chisq()
	  dptr3 = param_2d_buf;
	  param_idx2 = param_ct * param_ct;
	  for (param_idx = 0; param_idx < param_idx2; param_idx++) {
	    *dptr2++ = (*dptr3++) * dyy;
	  }
	  for (param_idx = 1; param_idx < param_ct; param_idx++) {
	    *dptr++ = param_2d_buf2[param_idx * param_ct_p1];
	  }
	  if (!linear_hypothesis_chisq(constraint_ct, param_ct, constraints_con_major, &(coef[perm_idx * indiv_valid_ct]), param_2d_buf2, param_df_buf, param_df_buf2, df_df_buf, mi_buf, df_buf, &dyy)) {
	    *dptr++ = dyy;
	  } else {
	    // joint test failure does not make the entire permutation count
	    // as failed
	    *dptr++ = -9;
	  }
	} else {
	  *dptr++ = -9;
	}
      }
    }
  } else {
    pheno_ptr = perm_pmajor; // Y[i]
    for (perm_idx = 0; perm_idx < cur_batch_size; perm_idx++) {
      fill_double_zero(cluster_param_buf, cluster_ct1_p1 * param_ct);
      dptr2 = covars_indiv_major; // X[i][j]
      if (!missing_ct) {
	for (indiv_idx = 0; indiv_idx < indiv_valid_ct; indiv_idx++) {
	  dptr3 = &(coef[perm_idx * indiv_valid_ct]);
	  partial = 0;
	  for (param_idx = 0; param_idx < param_ct; param_idx++) {
	    partial += dptr2[param_idx] * (*dptr3++);
	  }
	  partial -= *pheno_ptr++;
	  cluster_idx_p1 = indiv_to_cluster1[indiv_idx] + 1;
	  dptr3 = &(cluster_param_buf[cluster_idx_p1 * param_ct]);
	  for (param_idx = 0; param_idx < param_ct; param_idx++) {
	    // sc[clst[i]][j] += partial * X[i][j]
	    *dptr3 += partial * (*dptr2++);
	    dptr3++;
	  }
	}
      } else {
	dyy = 1;
	dzz = 0;
	if (standard_beta) {
	  pheno_sum = pheno_sum_base;
	  pheno_ssq = pheno_ssq_base;
	  ulptr = loadbuf;
	  for (indiv_uidx = 0; indiv_uidx < indiv_valid_ct + missing_ct; indiv_uidx += BITCT2) {
	    cur_word = *ulptr++;
	    cur_word = cur_word & (~(cur_word >> 1)) & FIVEMASK;
	    while (cur_word) {
	      dyy = pheno_ptr[indiv_uidx + (CTZLU(cur_word) / 2)];
	      pheno_sum -= dyy;
	      pheno_ssq -= dyy * dyy;
	      cur_word &= cur_word - 1;
	    }
	  }
	  dzz = pheno_sum / ((double)((intptr_t)indiv_valid_ct));
          dyy = sqrt(((double)((intptr_t)(indiv_valid_ct - 1))) / (pheno_ssq - pheno_sum * dzz));
	}
	for (indiv_uidx = 0, indiv_idx = 0; indiv_idx < indiv_valid_ct; indiv_uidx++, indiv_idx++) {
	  dptr3 = &(coef[perm_idx * indiv_valid_ct]);
	  partial = 0;
	  for (param_idx = 0; param_idx < param_ct; param_idx++) {
	    partial += dptr2[param_idx] * (*dptr3++);
	  }
	  while (1) {
	    if (indiv_uidx % BITCT2) {
	      cur_word >>= 2;
	    } else {
	      cur_word = loadbuf[indiv_uidx / BITCT2];
	      cur_word = cur_word & (~(cur_word >> 1)) & FIVEMASK;
	    }
	    if (!(cur_word & 1)) {
	      break;
	    }
	    indiv_uidx++;
	  }
	  partial -= (pheno_ptr[indiv_uidx] - dzz) * dyy;
	  cluster_idx_p1 = indiv_to_cluster1[indiv_idx] + 1;
	  dptr3 = &(cluster_param_buf[cluster_idx_p1 * param_ct]);
	  for (param_idx = 0; param_idx < param_ct; param_idx++) {
	    *dptr3 += partial * (*dptr2++);
	    dptr3++;
	  }
	}
	pheno_ptr = &(pheno_ptr[indiv_valid_ct + missing_ct]);
      }
      transpose_copy(cluster_ct1_p1, param_ct, cluster_param_buf, cluster_param_buf2);
      // can't overwrite param_2d_buf (= S0), everything else is fair game
      col_major_matrix_multiply(param_ct, param_ct, cluster_ct1_p1, cluster_param_buf, cluster_param_buf2, param_2d_buf2);
      col_major_matrix_multiply(param_ct, param_ct, param_ct, param_2d_buf, param_2d_buf2, cluster_param_buf); // multMatrix (S0, meat, tmp1)
      col_major_matrix_multiply(param_ct, param_ct, param_ct, cluster_param_buf, param_2d_buf, param_2d_buf2); // multMatrix (tmp1, S0, S)

      // now do validParameters() check.  validate and cache sqrt(S[i][i])...
      for (param_idx = 1; param_idx < param_ct; param_idx++) {
	dxx = param_2d_buf2[param_idx * param_ct_p1];
        if ((dxx < 1e-20) || (!realnum(dxx))) {
	  break;
	}
        cluster_param_buf[param_idx] = sqrt(dxx);
      }
      if (param_idx == param_ct) {
	cluster_param_buf[0] = sqrt(param_2d_buf2[0]);
        for (param_idx = 1; param_idx < param_ct; param_idx++) {
	  dxx = 0.99999 * cluster_param_buf[param_idx];
	  dptr = &(param_2d_buf2[param_idx * param_ct]); // S[i][j]
	  dptr2 = cluster_param_buf;
	  for (param_idx2 = 0; param_idx2 < param_idx; param_idx2++) {
	    if ((*dptr++) > dxx * (*dptr2++)) {
	      goto glm_linear_robust_cluster_covar_multicollinear;
	    }
	  }
	}
	dptr = &(linear_results[perm_idx * param_ctx_m1]);
        for (param_idx = 1; param_idx < param_ct; param_idx++) {
	  *dptr++ = param_2d_buf2[param_idx * param_ct_p1];
	}
	if (joint_test_requested) {
	  if (constraint_ct && (!linear_hypothesis_chisq(constraint_ct, param_ct, constraints_con_major, &(coef[perm_idx * indiv_valid_ct]), param_2d_buf2, param_df_buf, param_df_buf2, df_df_buf, mi_buf, df_buf, &dxx))) {
	    *dptr++ = dxx;
	  } else {
            *dptr++ = -9;
	  }
	}
      } else {
      glm_linear_robust_cluster_covar_multicollinear:
	// technically may not need to fill with zeroes/-9
        fill_double_zero(&(linear_results[perm_idx * param_ctx_m1]), param_ct_m1);
	SET_BIT(perm_fails, perm_idx);
	perm_fail_ct++;
	if (joint_test_requested) {
          linear_results[perm_idx * param_ctx_m1 + param_ct_m1] = -9;
	}
      }
    }
  }
  *perm_fail_ct_ptr = perm_fail_ct;
  return 0;
}

// #####
// The following code is based on the winning submissino of Pascal Pons in the
// "GWASSpeedup" contest run in April 2013 by Babbage Analytics & Innovation
// and TopCoder, who have donated the results to be used in PLINK.  The contest
// was designed by Po-Ru Loh; subsequent analysis and code preparation were
// performed by Andrew Hill, Ragu Bharadwaj, and Scott Jelinsky.  A manuscript
// is in preparation by these authors and Iain Kilty, Kevin Boudreau, Karim
// Lakhani and Eva Guinan.
// #####

#ifdef __LP64__
// exp_ps is a C port of Shigeo Mitsunari's fast math library posted at
// http://homepage1.nifty.com/herumi/ .  License is
// http://opensource.org/licenses/BSD-3-Clause .

// programmatically generated by:
// typedef union {
//   float f4;
//   uint32_t u4;
// } __uni4;
//
// __uni4 u4;
// int32_t ii;
// for (ii = 0; ii < 1024; ii++) {
//   u4.f4 = pow(2.0f, ((float)ii) / 1024.0);
//   printf("0x%08x", u4.u4 & 0x7fffff);
//   if (ii % 4 != 3) {
//     printf(", ");
//   } else {
//     printf(",\n");
//   }
// }
const uint32_t float_exp_lookup_int[] __attribute__((aligned(16))) = {
0x00000000, 0x00001630, 0x00002c64, 0x0000429c,
0x000058d8, 0x00006f17, 0x0000855b, 0x00009ba2,
0x0000b1ed, 0x0000c83c, 0x0000de8f, 0x0000f4e6,
0x00010b41, 0x0001219f, 0x00013802, 0x00014e68,
0x000164d2, 0x00017b40, 0x000191b2, 0x0001a828,
0x0001bea1, 0x0001d51f, 0x0001eba1, 0x00020226,
0x000218af, 0x00022f3c, 0x000245ce, 0x00025c63,
0x000272fc, 0x00028998, 0x0002a039, 0x0002b6de,
0x0002cd87, 0x0002e433, 0x0002fae4, 0x00031198,
0x00032850, 0x00033f0d, 0x000355cd, 0x00036c91,
0x00038359, 0x00039a25, 0x0003b0f5, 0x0003c7c9,
0x0003dea1, 0x0003f57d, 0x00040c5d, 0x00042341,
0x00043a29, 0x00045115, 0x00046804, 0x00047ef8,
0x000495f0, 0x0004aceb, 0x0004c3eb, 0x0004daef,
0x0004f1f6, 0x00050902, 0x00052012, 0x00053725,
0x00054e3d, 0x00056558, 0x00057c78, 0x0005939c,
0x0005aac3, 0x0005c1ef, 0x0005d91f, 0x0005f052,
0x0006078a, 0x00061ec6, 0x00063606, 0x00064d4a,
0x00066491, 0x00067bdd, 0x0006932d, 0x0006aa81,
0x0006c1d9, 0x0006d935, 0x0006f095, 0x000707f9,
0x00071f62, 0x000736ce, 0x00074e3e, 0x000765b3,
0x00077d2b, 0x000794a8, 0x0007ac28, 0x0007c3ad,
0x0007db35, 0x0007f2c2, 0x00080a53, 0x000821e8,
0x00083981, 0x0008511e, 0x000868c0, 0x00088065,
0x0008980f, 0x0008afbc, 0x0008c76e, 0x0008df23,
0x0008f6dd, 0x00090e9b, 0x0009265d, 0x00093e24,
0x000955ee, 0x00096dbc, 0x0009858f, 0x00099d66,
0x0009b541, 0x0009cd20, 0x0009e503, 0x0009fcea,
0x000a14d5, 0x000a2cc5, 0x000a44b9, 0x000a5cb1,
0x000a74ad, 0x000a8cad, 0x000aa4b1, 0x000abcba,
0x000ad4c6, 0x000aecd7, 0x000b04ec, 0x000b1d05,
0x000b3523, 0x000b4d44, 0x000b656a, 0x000b7d94,
0x000b95c2, 0x000badf4, 0x000bc62b, 0x000bde65,
0x000bf6a4, 0x000c0ee7, 0x000c272f, 0x000c3f7a,
0x000c57ca, 0x000c701e, 0x000c8876, 0x000ca0d2,
0x000cb933, 0x000cd198, 0x000cea01, 0x000d026e,
0x000d1adf, 0x000d3355, 0x000d4bcf, 0x000d644d,
0x000d7cd0, 0x000d9556, 0x000dade1, 0x000dc671,
0x000ddf04, 0x000df79c, 0x000e1038, 0x000e28d8,
0x000e417d, 0x000e5a25, 0x000e72d3, 0x000e8b84,
0x000ea43a, 0x000ebcf3, 0x000ed5b2, 0x000eee74,
0x000f073b, 0x000f2006, 0x000f38d5, 0x000f51a9,
0x000f6a81, 0x000f835d, 0x000f9c3e, 0x000fb523,
0x000fce0c, 0x000fe6fa, 0x000fffec, 0x001018e2,
0x001031dc, 0x00104adb, 0x001063de, 0x00107ce6,
0x001095f2, 0x0010af02, 0x0010c816, 0x0010e12f,
0x0010fa4d, 0x0011136e, 0x00112c94, 0x001145be,
0x00115eed, 0x00117820, 0x00119158, 0x0011aa93,
0x0011c3d3, 0x0011dd18, 0x0011f661, 0x00120fae,
0x00122900, 0x00124256, 0x00125bb0, 0x0012750f,
0x00128e72, 0x0012a7da, 0x0012c146, 0x0012dab7,
0x0012f42c, 0x00130da5, 0x00132723, 0x001340a5,
0x00135a2b, 0x001373b6, 0x00138d46, 0x0013a6d9,
0x0013c072, 0x0013da0e, 0x0013f3af, 0x00140d55,
0x001426ff, 0x001440ae, 0x00145a60, 0x00147418,
0x00148dd4, 0x0014a794, 0x0014c159, 0x0014db22,
0x0014f4f0, 0x00150ec2, 0x00152898, 0x00154274,
0x00155c53, 0x00157637, 0x00159020, 0x0015aa0d,
0x0015c3ff, 0x0015ddf5, 0x0015f7ef, 0x001611ee,
0x00162bf2, 0x001645fa, 0x00166006, 0x00167a18,
0x0016942d, 0x0016ae47, 0x0016c866, 0x0016e289,
0x0016fcb1, 0x001716dd, 0x0017310e, 0x00174b43,
0x0017657d, 0x00177fbc, 0x001799ff, 0x0017b446,
0x0017ce92, 0x0017e8e3, 0x00180338, 0x00181d92,
0x001837f0, 0x00185253, 0x00186cbb, 0x00188727,
0x0018a197, 0x0018bc0d, 0x0018d686, 0x0018f105,
0x00190b88, 0x0019260f, 0x0019409c, 0x00195b2c,
0x001975c2, 0x0019905c, 0x0019aafa, 0x0019c59e,
0x0019e046, 0x0019faf2, 0x001a15a3, 0x001a3059,
0x001a4b13, 0x001a65d2, 0x001a8096, 0x001a9b5e,
0x001ab62b, 0x001ad0fd, 0x001aebd3, 0x001b06ae,
0x001b218d, 0x001b3c71, 0x001b575a, 0x001b7248,
0x001b8d3a, 0x001ba831, 0x001bc32c, 0x001bde2c,
0x001bf931, 0x001c143b, 0x001c2f49, 0x001c4a5c,
0x001c6573, 0x001c8090, 0x001c9bb1, 0x001cb6d6,
0x001cd201, 0x001ced30, 0x001d0864, 0x001d239c,
0x001d3eda, 0x001d5a1c, 0x001d7562, 0x001d90ae,
0x001dabfe, 0x001dc753, 0x001de2ad, 0x001dfe0b,
0x001e196e, 0x001e34d6, 0x001e5043, 0x001e6bb4,
0x001e872a, 0x001ea2a5, 0x001ebe25, 0x001ed9a9,
0x001ef532, 0x001f10c0, 0x001f2c53, 0x001f47eb,
0x001f6387, 0x001f7f28, 0x001f9ace, 0x001fb679,
0x001fd228, 0x001feddc, 0x00200996, 0x00202553,
0x00204116, 0x00205cde, 0x002078aa, 0x0020947b,
0x0020b051, 0x0020cc2c, 0x0020e80b, 0x002103f0,
0x00211fd9, 0x00213bc7, 0x002157ba, 0x002173b2,
0x00218faf, 0x0021abb0, 0x0021c7b7, 0x0021e3c2,
0x0021ffd2, 0x00221be7, 0x00223801, 0x0022541f,
0x00227043, 0x00228c6b, 0x0022a899, 0x0022c4cb,
0x0022e102, 0x0022fd3e, 0x0023197f, 0x002335c5,
0x0023520f, 0x00236e5f, 0x00238ab3, 0x0023a70d,
0x0023c36b, 0x0023dfce, 0x0023fc37, 0x002418a4,
0x00243516, 0x0024518d, 0x00246e08, 0x00248a89,
0x0024a70f, 0x0024c39a, 0x0024e029, 0x0024fcbe,
0x00251958, 0x002535f6, 0x00255299, 0x00256f42,
0x00258bef, 0x0025a8a2, 0x0025c559, 0x0025e215,
0x0025fed7, 0x00261b9d, 0x00263868, 0x00265538,
0x0026720e, 0x00268ee8, 0x0026abc7, 0x0026c8ac,
0x0026e595, 0x00270283, 0x00271f76, 0x00273c6f,
0x0027596c, 0x0027766e, 0x00279376, 0x0027b082,
0x0027cd94, 0x0027eaaa, 0x002807c6, 0x002824e6,
0x0028420c, 0x00285f37, 0x00287c66, 0x0028999b,
0x0028b6d5, 0x0028d414, 0x0028f158, 0x00290ea1,
0x00292bef, 0x00294942, 0x0029669b, 0x002983f8,
0x0029a15b, 0x0029bec2, 0x0029dc2f, 0x0029f9a1,
0x002a1718, 0x002a3494, 0x002a5215, 0x002a6f9b,
0x002a8d26, 0x002aaab7, 0x002ac84c, 0x002ae5e7,
0x002b0387, 0x002b212c, 0x002b3ed6, 0x002b5c85,
0x002b7a3a, 0x002b97f3, 0x002bb5b2, 0x002bd376,
0x002bf13f, 0x002c0f0d, 0x002c2ce0, 0x002c4ab9,
0x002c6897, 0x002c867a, 0x002ca462, 0x002cc24f,
0x002ce041, 0x002cfe39, 0x002d1c36, 0x002d3a38,
0x002d583f, 0x002d764b, 0x002d945d, 0x002db274,
0x002dd090, 0x002deeb1, 0x002e0cd8, 0x002e2b03,
0x002e4934, 0x002e676b, 0x002e85a6, 0x002ea3e7,
0x002ec22d, 0x002ee078, 0x002efec8, 0x002f1d1e,
0x002f3b79, 0x002f59d9, 0x002f783e, 0x002f96a9,
0x002fb519, 0x002fd38e, 0x002ff209, 0x00301089,
0x00302f0e, 0x00304d98, 0x00306c28, 0x00308abd,
0x0030a957, 0x0030c7f7, 0x0030e69c, 0x00310546,
0x003123f6, 0x003142aa, 0x00316165, 0x00318024,
0x00319ee9, 0x0031bdb3, 0x0031dc83, 0x0031fb57,
0x00321a32, 0x00323911, 0x003257f6, 0x003276e0,
0x003295d0, 0x0032b4c5, 0x0032d3bf, 0x0032f2bf,
0x003311c4, 0x003330cf, 0x00334fde, 0x00336ef4,
0x00338e0e, 0x0033ad2e, 0x0033cc54, 0x0033eb7e,
0x00340aaf, 0x003429e4, 0x0034491f, 0x00346860,
0x003487a6, 0x0034a6f1, 0x0034c642, 0x0034e598,
0x003504f3, 0x00352454, 0x003543bb, 0x00356327,
0x00358298, 0x0035a20f, 0x0035c18b, 0x0035e10d,
0x00360094, 0x00362020, 0x00363fb2, 0x00365f4a,
0x00367ee7, 0x00369e89, 0x0036be31, 0x0036dddf,
0x0036fd92, 0x00371d4a, 0x00373d08, 0x00375ccc,
0x00377c95, 0x00379c63, 0x0037bc37, 0x0037dc11,
0x0037fbf0, 0x00381bd4, 0x00383bbe, 0x00385bae,
0x00387ba3, 0x00389b9e, 0x0038bb9e, 0x0038dba4,
0x0038fbaf, 0x00391bc0, 0x00393bd7, 0x00395bf3,
0x00397c14, 0x00399c3b, 0x0039bc68, 0x0039dc9a,
0x0039fcd2, 0x003a1d10, 0x003a3d53, 0x003a5d9b,
0x003a7dea, 0x003a9e3e, 0x003abe97, 0x003adef6,
0x003aff5b, 0x003b1fc5, 0x003b4035, 0x003b60aa,
0x003b8126, 0x003ba1a6, 0x003bc22d, 0x003be2b9,
0x003c034a, 0x003c23e2, 0x003c447f, 0x003c6521,
0x003c85ca, 0x003ca678, 0x003cc72b, 0x003ce7e5,
0x003d08a4, 0x003d2968, 0x003d4a33, 0x003d6b03,
0x003d8bd8, 0x003dacb4, 0x003dcd95, 0x003dee7c,
0x003e0f68, 0x003e305a, 0x003e5152, 0x003e7250,
0x003e9353, 0x003eb45c, 0x003ed56b, 0x003ef67f,
0x003f179a, 0x003f38ba, 0x003f59df, 0x003f7b0b,
0x003f9c3c, 0x003fbd73, 0x003fdeb0, 0x003ffff2,
0x0040213b, 0x00404289, 0x004063dc, 0x00408536,
0x0040a695, 0x0040c7fb, 0x0040e966, 0x00410ad6,
0x00412c4d, 0x00414dc9, 0x00416f4b, 0x004190d3,
0x0041b261, 0x0041d3f5, 0x0041f58e, 0x0042172d,
0x004238d2, 0x00425a7d, 0x00427c2e, 0x00429de4,
0x0042bfa1, 0x0042e163, 0x0043032b, 0x004324f9,
0x004346cd, 0x004368a7, 0x00438a86, 0x0043ac6b,
0x0043ce57, 0x0043f048, 0x0044123f, 0x0044343c,
0x0044563f, 0x00447848, 0x00449a56, 0x0044bc6b,
0x0044de85, 0x004500a5, 0x004522cc, 0x004544f8,
0x0045672a, 0x00458962, 0x0045aba0, 0x0045cde4,
0x0045f02e, 0x0046127e, 0x004634d3, 0x0046572f,
0x00467991, 0x00469bf8, 0x0046be66, 0x0046e0d9,
0x00470353, 0x004725d2, 0x00474858, 0x00476ae3,
0x00478d75, 0x0047b00c, 0x0047d2aa, 0x0047f54d,
0x004817f7, 0x00483aa6, 0x00485d5b, 0x00488017,
0x0048a2d8, 0x0048c5a0, 0x0048e86d, 0x00490b41,
0x00492e1b, 0x004950fa, 0x004973e0, 0x004996cc,
0x0049b9be, 0x0049dcb5, 0x0049ffb3, 0x004a22b7,
0x004a45c1, 0x004a68d1, 0x004a8be8, 0x004aaf04,
0x004ad226, 0x004af54f, 0x004b187d, 0x004b3bb2,
0x004b5eed, 0x004b822e, 0x004ba575, 0x004bc8c2,
0x004bec15, 0x004c0f6e, 0x004c32ce, 0x004c5633,
0x004c799f, 0x004c9d11, 0x004cc089, 0x004ce407,
0x004d078c, 0x004d2b16, 0x004d4ea7, 0x004d723d,
0x004d95da, 0x004db97e, 0x004ddd27, 0x004e00d6,
0x004e248c, 0x004e4848, 0x004e6c0a, 0x004e8fd2,
0x004eb3a1, 0x004ed775, 0x004efb50, 0x004f1f31,
0x004f4319, 0x004f6706, 0x004f8afa, 0x004faef4,
0x004fd2f4, 0x004ff6fb, 0x00501b08, 0x00503f1b,
0x00506334, 0x00508753, 0x0050ab79, 0x0050cfa5,
0x0050f3d7, 0x00511810, 0x00513c4f, 0x00516094,
0x005184df, 0x0051a931, 0x0051cd89, 0x0051f1e7,
0x0052164c, 0x00523ab7, 0x00525f28, 0x005283a0,
0x0052a81e, 0x0052cca2, 0x0052f12c, 0x005315bd,
0x00533a54, 0x00535ef2, 0x00538396, 0x0053a840,
0x0053ccf1, 0x0053f1a8, 0x00541665, 0x00543b29,
0x00545ff3, 0x005484c3, 0x0054a99a, 0x0054ce77,
0x0054f35b, 0x00551845, 0x00553d35, 0x0055622c,
0x00558729, 0x0055ac2d, 0x0055d137, 0x0055f647,
0x00561b5e, 0x0056407b, 0x0056659f, 0x00568ac9,
0x0056affa, 0x0056d531, 0x0056fa6e, 0x00571fb2,
0x005744fd, 0x00576a4e, 0x00578fa5, 0x0057b503,
0x0057da67, 0x0057ffd2, 0x00582543, 0x00584abb,
0x00587039, 0x005895be, 0x0058bb49, 0x0058e0db,
0x00590673, 0x00592c12, 0x005951b8, 0x00597763,
0x00599d16, 0x0059c2cf, 0x0059e88e, 0x005a0e54,
0x005a3421, 0x005a59f4, 0x005a7fcd, 0x005aa5ae,
0x005acb94, 0x005af182, 0x005b1776, 0x005b3d70,
0x005b6371, 0x005b8979, 0x005baf87, 0x005bd59c,
0x005bfbb8, 0x005c21da, 0x005c4802, 0x005c6e32,
0x005c9468, 0x005cbaa4, 0x005ce0e7, 0x005d0731,
0x005d2d82, 0x005d53d9, 0x005d7a36, 0x005da09b,
0x005dc706, 0x005ded77, 0x005e13f0, 0x005e3a6f,
0x005e60f5, 0x005e8781, 0x005eae14, 0x005ed4ae,
0x005efb4e, 0x005f21f5, 0x005f48a3, 0x005f6f58,
0x005f9613, 0x005fbcd5, 0x005fe39e, 0x00600a6d,
0x00603143, 0x00605820, 0x00607f03, 0x0060a5ee,
0x0060ccdf, 0x0060f3d7, 0x00611ad5, 0x006141db,
0x006168e7, 0x00618ffa, 0x0061b713, 0x0061de34,
0x0062055b, 0x00622c89, 0x006253be, 0x00627af9,
0x0062a23c, 0x0062c985, 0x0062f0d5, 0x0063182c,
0x00633f89, 0x006366ee, 0x00638e59, 0x0063b5cb,
0x0063dd44, 0x006404c4, 0x00642c4b, 0x006453d8,
0x00647b6d, 0x0064a308, 0x0064caaa, 0x0064f253,
0x00651a03, 0x006541b9, 0x00656977, 0x0065913c,
0x0065b907, 0x0065e0d9, 0x006608b2, 0x00663092,
0x00665879, 0x00668067, 0x0066a85c, 0x0066d058,
0x0066f85b, 0x00672064, 0x00674875, 0x0067708c,
0x006798ab, 0x0067c0d0, 0x0067e8fd, 0x00681130,
0x0068396a, 0x006861ac, 0x006889f4, 0x0068b243,
0x0068da99, 0x006902f7, 0x00692b5b, 0x006953c6,
0x00697c38, 0x0069a4b1, 0x0069cd32, 0x0069f5b9,
0x006a1e47, 0x006a46dd, 0x006a6f79, 0x006a981c,
0x006ac0c7, 0x006ae978, 0x006b1231, 0x006b3af1,
0x006b63b7, 0x006b8c85, 0x006bb55a, 0x006bde36,
0x006c0719, 0x006c3003, 0x006c58f4, 0x006c81ec,
0x006caaec, 0x006cd3f2, 0x006cfd00, 0x006d2614,
0x006d4f30, 0x006d7853, 0x006da17d, 0x006dcaae,
0x006df3e7, 0x006e1d26, 0x006e466d, 0x006e6fbb,
0x006e9910, 0x006ec26c, 0x006eebcf, 0x006f1539,
0x006f3eab, 0x006f6824, 0x006f91a4, 0x006fbb2b,
0x006fe4ba, 0x00700e4f, 0x007037ec, 0x00706190,
0x00708b3b, 0x0070b4ee, 0x0070dea8, 0x00710868,
0x00713231, 0x00715c00, 0x007185d7, 0x0071afb5,
0x0071d99a, 0x00720386, 0x00722d7a, 0x00725775,
0x00728177, 0x0072ab81, 0x0072d592, 0x0072ffaa,
0x007329c9, 0x007353f0, 0x00737e1e, 0x0073a853,
0x0073d290, 0x0073fcd4, 0x0074271f, 0x00745172,
0x00747bcc, 0x0074a62d, 0x0074d096, 0x0074fb06,
0x0075257d, 0x00754ffc, 0x00757a82, 0x0075a50f,
0x0075cfa4, 0x0075fa40, 0x007624e4, 0x00764f8f,
0x00767a41, 0x0076a4fb, 0x0076cfbc, 0x0076fa85,
0x00772555, 0x0077502d, 0x00777b0b, 0x0077a5f2,
0x0077d0df, 0x0077fbd5, 0x007826d1, 0x007851d5,
0x00787ce1, 0x0078a7f4, 0x0078d30e, 0x0078fe30,
0x0079295a, 0x0079548b, 0x00797fc3, 0x0079ab03,
0x0079d64a, 0x007a0199, 0x007a2cf0, 0x007a584d,
0x007a83b3, 0x007aaf20, 0x007ada94, 0x007b0610,
0x007b3194, 0x007b5d1f, 0x007b88b2, 0x007bb44c,
0x007bdfed, 0x007c0b97, 0x007c3748, 0x007c6300,
0x007c8ec0, 0x007cba88, 0x007ce657, 0x007d122e,
0x007d3e0c, 0x007d69f2, 0x007d95e0, 0x007dc1d5,
0x007dedd2, 0x007e19d6, 0x007e45e2, 0x007e71f6,
0x007e9e11, 0x007eca34, 0x007ef65f, 0x007f2291,
0x007f4ecb, 0x007f7b0d, 0x007fa756, 0x007fd3a7
};

const float* const float_exp_lookup = (const float*)float_exp_lookup_int;

static inline __m128 fmath_exp_ps(__m128 xx) {
  const __m128i mask7ff = {0x7fffffff7fffffffLLU, 0x7fffffff7fffffffLLU};

  // 88
  const __m128i max_x = {0x42b0000042b00000LLU, 0x42b0000042b00000LLU};
  // -88
  // more sensible 0xc2b00000... not used here due to "narrowing conversion"
  // warning
  const __m128i min_x = {-0x3d4fffff3d500000LL, -0x3d4fffff3d500000LL};
  // 2^10 / log(2)
  const __m128i const_aa = {0x44b8aa3b44b8aa3bLLU, 0x44b8aa3b44b8aa3bLLU};
  // log(2) / 2^10
  const __m128i const_bb = {0x3a3172183a317218LLU, 0x3a3172183a317218LLU};

  const __m128i f1 = {0x3f8000003f800000LLU, 0x3f8000003f800000LLU};
  const __m128i mask_s = {0x3ff000003ffLLU, 0x3ff000003ffLLU};
  const __m128i i127s = {0x1fc000001fc00LLU, 0x1fc000001fc00LLU};
  __m128i limit = _mm_castps_si128(_mm_and_ps(xx, (__m128)mask7ff));
  int32_t over = _mm_movemask_epi8(_mm_cmpgt_epi32(limit, max_x));
  if (over) {
    xx = _mm_min_ps(xx, (__m128)max_x);
    xx = _mm_max_ps(xx, (__m128)min_x);
  }
  __m128i rr = _mm_cvtps_epi32(_mm_mul_ps(xx, (__m128)const_aa));
  __m128 tt = _mm_sub_ps(xx, _mm_mul_ps(_mm_cvtepi32_ps(rr), (__m128)const_bb));
  tt = _mm_add_ps(tt, (__m128)f1);
  __m128i v4 = _mm_and_si128(rr, mask_s);
  __m128i u4 = _mm_add_epi32(rr, i127s);
  u4 = _mm_srli_epi32(u4, 10);
  u4 = _mm_slli_epi32(u4, 23);
  uint32_t v0 = _mm_cvtsi128_si32(v4);
  uint32_t v1 = ((int32_t)(uint16_t)__builtin_ia32_vec_ext_v8hi((__v8hi)(__m128i)(v4), (int32_t)(2)));
  uint32_t v2 = ((int32_t)(uint16_t)__builtin_ia32_vec_ext_v8hi((__v8hi)(__m128i)(v4), (int32_t)(4)));
  uint32_t v3 = ((int32_t)(uint16_t)__builtin_ia32_vec_ext_v8hi((__v8hi)(__m128i)(v4), (int32_t)(6)));
  __m128 t0 = _mm_set_ss(float_exp_lookup[v0]);
  __m128 t1 = _mm_set_ss(float_exp_lookup[v1]);
  __m128 t2 = _mm_set_ss(float_exp_lookup[v2]);
  __m128 t3 = _mm_set_ss(float_exp_lookup[v3]);
  t1 = _mm_movelh_ps(t1, t3);
  t1 = _mm_castsi128_ps(_mm_slli_epi64(_mm_castps_si128(t1), 32));
  t0 = _mm_movelh_ps(t0, t2);
  t0 = _mm_or_ps(t0, t1);
  t0 = _mm_or_ps(t0, _mm_castsi128_ps(u4));
  tt = _mm_mul_ps(tt, t0);
  return tt;
}

// For equivalent "normal" C/C++ code, see the non-__LP64__ versions of these
// functions.
static inline void logistic_sse(float* vect, uint32_t nn) {
  __m128 zero = _mm_setzero_ps();
  __m128 one = _mm_set1_ps(1.0);
  uint32_t uii;
  for (uii = 0; uii < nn; uii += 4) {
    __m128 aa = _mm_load_ps(&(vect[uii]));
    aa = _mm_sub_ps(zero, aa);
    aa = fmath_exp_ps(aa);
    aa = _mm_add_ps(aa, one);
    aa = _mm_div_ps(one, aa);
    _mm_store_ps(&(vect[uii]), aa);
  }
}

static inline void compute_v_and_p_minus_y(float* pp, float* vv, const float* yy, uint32_t nn) {
  __m128 one = _mm_set1_ps(1.0);
  uint32_t uii;
  for (uii = 0; uii < nn; uii += 4) {
    __m128 ptmp = _mm_load_ps(&(pp[uii]));
    __m128 one_minus_ptmp = _mm_sub_ps(one, ptmp);
    _mm_store_ps(&(vv[uii]), _mm_mul_ps(ptmp, one_minus_ptmp));
    __m128 ytmp = _mm_load_ps(&(yy[uii]));
    _mm_store_ps(&(pp[uii]), _mm_sub_ps(ptmp, ytmp));
  }
}

static inline void mult_tmatrix_nxd_vect_d(const float* tm, const float* vect, float* dest, uint32_t col_ct, uint32_t row_ct) {
  // tm is row-major, cols are packed to 16-byte alignment
  // "col_cta4" = col_ct, aligned to multiple of 4.  Since 16-byte blocks
  // contain 4 floats each, this is the actual length (in floats) of each tm
  // row.  (Yes, I need to standardize a zillion other variable names of this
  // sort...)
  __m128 w1;
  __m128 w2;
  __m128 w3;
  __m128 w4;
  __m128 r1;
  __m128 r2;
  __m128 r3;
  __m128 r4;
  uintptr_t col_cta4 = (col_ct + 3) & (~3);
  uint32_t row_idx = 0;
  uint32_t row_ctm3;
  uint32_t col_idx;
  if (row_ct < 4) {
    memset(dest, 0, col_ct * sizeof(float));
  } else {
    w1 = _mm_load1_ps(vect);
    w2 = _mm_load1_ps(&(vect[1]));
    w3 = _mm_load1_ps(&(vect[2]));
    w4 = _mm_load1_ps(&(vect[3]));
    for (col_idx = 0; col_idx < col_ct; col_idx += 4) {
      r1 = _mm_load_ps(&(tm[col_idx]));
      r2 = _mm_load_ps(&(tm[col_idx + col_cta4]));
      r3 = _mm_load_ps(&(tm[col_idx + 2 * col_cta4]));
      r4 = _mm_load_ps(&(tm[col_idx + 3 * col_cta4]));
      r1 = _mm_mul_ps(r1, w1);
      r2 = _mm_mul_ps(r2, w2);
      r3 = _mm_mul_ps(r3, w3);
      r4 = _mm_mul_ps(r4, w4);
      r1 = _mm_add_ps(r1, r2);
      r3 = _mm_add_ps(r3, r4);
      r1 = _mm_add_ps(r1, r3);
      _mm_store_ps(&(dest[col_idx]), r1);
    }
    row_ctm3 = row_ct - 3;
    for (row_idx = 4; row_idx < row_ctm3; row_idx += 4) {
      w1 = _mm_load1_ps(&(vect[row_idx]));
      w2 = _mm_load1_ps(&(vect[row_idx + 1]));
      w3 = _mm_load1_ps(&(vect[row_idx + 2]));
      w4 = _mm_load1_ps(&(vect[row_idx + 3]));
      for (col_idx = 0; col_idx < col_ct; col_idx += 4) {
        r1 = _mm_load_ps(&(tm[col_idx + row_idx * col_cta4]));
        r2 = _mm_load_ps(&(tm[col_idx + (row_idx + 1) * col_cta4]));
        r3 = _mm_load_ps(&(tm[col_idx + (row_idx + 2) * col_cta4]));
        r4 = _mm_load_ps(&(tm[col_idx + (row_idx + 3) * col_cta4]));
        r1 = _mm_mul_ps(r1, w1);
        r2 = _mm_mul_ps(r2, w2);
        r3 = _mm_mul_ps(r3, w3);
        r4 = _mm_mul_ps(r4, w4);
        r1 = _mm_add_ps(r1, r2);
        r3 = _mm_add_ps(r3, r4);
        r1 = _mm_add_ps(r1, r3);
	r1 = _mm_add_ps(r1, _mm_load_ps(&(dest[col_idx])));
	_mm_store_ps(&(dest[col_idx]), r1);
      }
    }
  }
  switch(row_ct % 4) {
  case 3:
    w1 = _mm_load1_ps(&(vect[row_idx]));
    w2 = _mm_load1_ps(&(vect[row_idx + 1]));
    w3 = _mm_load1_ps(&(vect[row_idx + 2]));
    for (col_idx = 0; col_idx < col_ct; col_idx += 4) {
      r1 = _mm_load_ps(&(tm[col_idx + row_idx * col_cta4]));
      r2 = _mm_load_ps(&(tm[col_idx + (row_idx + 1) * col_cta4]));
      r3 = _mm_load_ps(&(tm[col_idx + (row_idx + 2) * col_cta4]));
      r1 = _mm_mul_ps(r1, w1);
      r2 = _mm_mul_ps(r2, w2);
      r3 = _mm_mul_ps(r3, w3);
      r1 = _mm_add_ps(r1, r2);
      r3 = _mm_add_ps(r3, _mm_load_ps(&(dest[col_idx])));
      r1 = _mm_add_ps(r1, r3);
      _mm_store_ps(&(dest[col_idx]), r1);
    }
    break;
  case 2:
    w1 = _mm_load1_ps(&(vect[row_idx]));
    w2 = _mm_load1_ps(&(vect[row_idx + 1]));
    for (col_idx = 0; col_idx < col_ct; col_idx += 4) {
      r1 = _mm_load_ps(&(tm[col_idx + row_idx * col_cta4]));
      r2 = _mm_load_ps(&(tm[col_idx + (row_idx + 1) * col_cta4]));
      r1 = _mm_mul_ps(r1, w1);
      r2 = _mm_mul_ps(r2, w2);
      r1 = _mm_add_ps(r1, r2);
      r1 = _mm_add_ps(r1, _mm_load_ps(&(dest[col_idx])));
      _mm_store_ps(&(dest[col_idx]), r1);
    }
    break;
  case 1:
    w1 = _mm_load1_ps(&(vect[row_idx]));
    for (col_idx = 0; col_idx < col_ct; col_idx += 4) {
      r1 = _mm_load_ps(&(tm[col_idx + row_idx * col_cta4]));
      r1 = _mm_mul_ps(r1, w1);
      r1 = _mm_add_ps(r1, _mm_load_ps(&(dest[col_idx])));
      _mm_store_ps(&(dest[col_idx]), r1);
    }
  }
}

static inline void mult_matrix_dxn_vect_n(const float* mm, const float* vect, float* dest, uint32_t col_ct, uint32_t row_ct) {
  uintptr_t col_cta4 = (col_ct + 3) & (~3);
  uint32_t row_idx = 0;
  const float* mm_ptr;
  __m128 s1;
  __m128 s2;
  __m128 s3;
  __m128 s4;
  __m128 vv;
  __m128 a1;
  __m128 a2;
  __m128 a3;
  __m128 a4;
  __uni16 u16;
  uint32_t row_ctm3;
  uint32_t col_idx;
  if (row_ct > 3) {
    row_ctm3 = row_ct - 3;
    for (; row_idx < row_ctm3; row_idx += 4) {
      s1 = _mm_setzero_ps();
      s2 = _mm_setzero_ps();
      s3 = _mm_setzero_ps();
      s4 = _mm_setzero_ps();
      for (col_idx = 0; col_idx < col_ct; col_idx += 4) {
	mm_ptr = &(mm[row_idx * col_cta4 + col_idx]);
        vv = _mm_load_ps(&(vect[col_idx]));
        a1 = _mm_load_ps(mm_ptr);
        a2 = _mm_load_ps(&(mm_ptr[col_cta4]));
        a3 = _mm_load_ps(&(mm_ptr[2 * col_cta4]));
        a4 = _mm_load_ps(&(mm_ptr[3 * col_cta4]));
	// want to switch this to fused multiply-add...
        a1 = _mm_mul_ps(a1, vv);
        a2 = _mm_mul_ps(a2, vv);
        a3 = _mm_mul_ps(a3, vv);
        a4 = _mm_mul_ps(a4, vv);
        s1 = _mm_add_ps(s1, a1);
        s2 = _mm_add_ps(s2, a2);
        s3 = _mm_add_ps(s3, a3);
        s4 = _mm_add_ps(s4, a4);
      }
      // refrain from using SSE3 _mm_hadd_ps() for now
      u16.vf = s1;
      *dest++ = u16.f4[0] + u16.f4[1] + u16.f4[2] + u16.f4[3];
      u16.vf = s2;
      *dest++ = u16.f4[0] + u16.f4[1] + u16.f4[2] + u16.f4[3];
      u16.vf = s3;
      *dest++ = u16.f4[0] + u16.f4[1] + u16.f4[2] + u16.f4[3];
      u16.vf = s4;
      *dest++ = u16.f4[0] + u16.f4[1] + u16.f4[2] + u16.f4[3];
    }
  }
  s1 = _mm_setzero_ps();
  s2 = _mm_setzero_ps();
  s3 = _mm_setzero_ps();
  switch (row_ct % 4) {
  case 3:
    for (col_idx = 0; col_idx < col_ct; col_idx += 4) {
      mm_ptr = &(mm[row_idx * col_cta4 + col_idx]);
      vv = _mm_load_ps(&(vect[col_idx]));
      a1 = _mm_load_ps(mm_ptr);
      a2 = _mm_load_ps(&(mm_ptr[col_cta4]));
      a3 = _mm_load_ps(&(mm_ptr[2 * col_cta4]));
      a1 = _mm_mul_ps(a1, vv);
      a2 = _mm_mul_ps(a2, vv);
      a3 = _mm_mul_ps(a3, vv);
      s1 = _mm_add_ps(s1, a1);
      s2 = _mm_add_ps(s2, a2);
      s3 = _mm_add_ps(s3, a3);
    }
    u16.vf = s1;
    *dest++ = u16.f4[0] + u16.f4[1] + u16.f4[2] + u16.f4[3];
    u16.vf = s2;
    *dest++ = u16.f4[0] + u16.f4[1] + u16.f4[2] + u16.f4[3];
    u16.vf = s3;
    *dest = u16.f4[0] + u16.f4[1] + u16.f4[2] + u16.f4[3];
    break;
  case 2:
    for (col_idx = 0; col_idx < col_ct; col_idx += 4) {
      mm_ptr = &(mm[row_idx * col_cta4 + col_idx]);
      vv = _mm_load_ps(&(vect[col_idx]));
      a1 = _mm_load_ps(mm_ptr);
      a2 = _mm_load_ps(&(mm_ptr[col_cta4]));
      a1 = _mm_mul_ps(a1, vv);
      a2 = _mm_mul_ps(a2, vv);
      s1 = _mm_add_ps(s1, a1);
      s2 = _mm_add_ps(s2, a2);
    }
    u16.vf = s1;
    *dest++ = u16.f4[0] + u16.f4[1] + u16.f4[2] + u16.f4[3];
    u16.vf = s2;
    *dest = u16.f4[0] + u16.f4[1] + u16.f4[2] + u16.f4[3];
    break;
  case 1:
    for (col_idx = 0; col_idx < col_ct; col_idx += 4) {
      vv = _mm_load_ps(&(vect[col_idx]));
      a1 = _mm_load_ps(&(mm[row_idx * col_cta4 + col_idx]));
      a1 = _mm_mul_ps(a1, vv);
      s1 = _mm_add_ps(s1, a1);
    }
    u16.vf = s1;
    *dest = u16.f4[0] + u16.f4[1] + u16.f4[2] + u16.f4[3];
    break;
  }
}

static inline float triple_product(const float* v1, const float* v2, const float* v3, uint32_t nn) {
  __m128 sum = _mm_setzero_ps();
  __m128 aa;
  __m128 bb;
  __m128 cc;
  __uni16 u16;
  uint32_t uii;
  for (uii = 0; uii < nn; uii += 4) {
    aa = _mm_load_ps(&(v1[uii]));
    bb = _mm_load_ps(&(v2[uii]));
    cc = _mm_load_ps(&(v3[uii]));
    sum = _mm_add_ps(sum, _mm_mul_ps(_mm_mul_ps(aa, bb), cc));
  }
  u16.vf = sum;
  return u16.f4[0] + u16.f4[1] + u16.f4[2] + u16.f4[3];
}

static inline void compute_two_diag_triple_product(const float* aa, const float* bb, const float* vv, float* raa_ptr, float* rab_ptr, float* rbb_ptr, uint32_t nn) {
  __m128 saa = _mm_setzero_ps();
  __m128 sab = _mm_setzero_ps();
  __m128 sbb = _mm_setzero_ps();
  __m128 vtmp;
  __m128 atmp;
  __m128 btmp;
  __m128 av;
  __m128 bv;
  __uni16 u16;
  uint32_t uii;
  for (uii = 0; uii < nn; uii += 4) {
    vtmp = _mm_load_ps(&(vv[uii]));
    atmp = _mm_load_ps(&(aa[uii]));
    btmp = _mm_load_ps(&(bb[uii]));
    av = _mm_mul_ps(atmp, vtmp);
    bv = _mm_mul_ps(btmp, vtmp);
    saa = _mm_add_ps(saa, _mm_mul_ps(atmp, av));
    sab = _mm_add_ps(sab, _mm_mul_ps(atmp, bv));
    sbb = _mm_add_ps(sbb, _mm_mul_ps(btmp, bv));
  }
  u16.vf = saa;
  *raa_ptr = u16.f4[0] + u16.f4[1] + u16.f4[2] + u16.f4[3];
  u16.vf = sab;
  *rab_ptr = u16.f4[0] + u16.f4[1] + u16.f4[2] + u16.f4[3];
  u16.vf = sbb;
  *rbb_ptr = u16.f4[0] + u16.f4[1] + u16.f4[2] + u16.f4[3];
}

static inline void compute_three_triple_product(const float* bb, const float* a1, const float* a2, const float* a3, const float* vv, float* r1_ptr, float* r2_ptr, float* r3_ptr, uint32_t nn) {
  __m128 s1 = _mm_setzero_ps();
  __m128 s2 = _mm_setzero_ps();
  __m128 s3 = _mm_setzero_ps();
  __m128 a1tmp;
  __m128 a2tmp;
  __m128 a3tmp;
  __m128 vtmp;
  __m128 btmp;
  __uni16 u16;
  uint32_t uii;
  for (uii = 0; uii < nn; uii += 4) {
    a1tmp = _mm_load_ps(&(a1[uii]));
    a2tmp = _mm_load_ps(&(a2[uii]));
    a3tmp = _mm_load_ps(&(a3[uii]));
    vtmp = _mm_load_ps(&(vv[uii]));
    btmp = _mm_load_ps(&(bb[uii]));
    btmp = _mm_mul_ps(btmp, vtmp);
    s1 = _mm_add_ps(s1, _mm_mul_ps(a1tmp, btmp));
    s2 = _mm_add_ps(s2, _mm_mul_ps(a2tmp, btmp));
    s3 = _mm_add_ps(s3, _mm_mul_ps(a3tmp, btmp));
  }
  u16.vf = s1;
  *r1_ptr = u16.f4[0] + u16.f4[1] + u16.f4[2] + u16.f4[3];
  u16.vf = s2;
  *r2_ptr = u16.f4[0] + u16.f4[1] + u16.f4[2] + u16.f4[3];
  u16.vf = s3;
  *r3_ptr = u16.f4[0] + u16.f4[1] + u16.f4[2] + u16.f4[3];
}

static inline void compute_two_plus_one_triple_product(const float* bb, const float* a1, const float* a2, const float* vv, float* r1_ptr, float* r2_ptr, float* r3_ptr, uint32_t nn) {
  __m128 s1 = _mm_setzero_ps();
  __m128 s2 = _mm_setzero_ps();
  __m128 s3 = _mm_setzero_ps();
  __m128 a1tmp;
  __m128 a2tmp;
  __m128 btmp;
  __m128 vtmp;
  __m128 bv;
  __uni16 u16;
  uint32_t uii;
  for (uii = 0; uii < nn; uii += 4) {
    a1tmp = _mm_load_ps(&(a1[uii]));
    a2tmp = _mm_load_ps(&(a2[uii]));
    btmp = _mm_load_ps(&(bb[uii]));
    vtmp = _mm_load_ps(&(vv[uii]));
    bv = _mm_mul_ps(btmp, vtmp);
    s1 = _mm_add_ps(s1, _mm_mul_ps(btmp, bv));
    s2 = _mm_add_ps(s2, _mm_mul_ps(a1tmp, bv));
    s3 = _mm_add_ps(s3, _mm_mul_ps(a2tmp, bv));
  }
  u16.vf = s1;
  *r1_ptr = u16.f4[0] + u16.f4[1] + u16.f4[2] + u16.f4[3];
  u16.vf = s2;
  *r2_ptr = u16.f4[0] + u16.f4[1] + u16.f4[2] + u16.f4[3];
  u16.vf = s3;
  *r3_ptr = u16.f4[0] + u16.f4[1] + u16.f4[2] + u16.f4[3];
}
#else // no __LP64__ (and hence, unsafe to assume presence of SSE2)
static inline void logistic_sse(float* vect, uint32_t nn) {
  uint32_t uii;
  for (uii = 0; uii < nn; uii++) {
    vect[uii] = 1.0 / (1 + exp(-vect[uii]));
  }
}

static inline void compute_v_and_p_minus_y(float* pp, float* vv, const float* yy, uint32_t nn) {
  uint32_t uii;
  for (uii = 0; uii < nn; uii++) {
    vv[uii] = pp[uii] * (1.0 - pp[uii]);
    pp[uii] -= yy[uii];
  }
}

static inline void mult_tmatrix_nxd_vect_d(const float* tm, const float* vect, float* dest, uint32_t col_ct, uint32_t row_ct) {
  uintptr_t col_cta4 = (col_ct + 3) & (~3);
  const float* tm_ptr;
  float vect_val;
  uint32_t col_idx;
  uint32_t row_idx;
  fill_float_zero(dest, col_ct);
  for (row_idx = 0; row_idx < row_ct; row_idx++) {
    vect_val = vect[row_idx];
    tm_ptr = &(tm[row_idx * col_cta4]);
    for (col_idx = 0; col_idx < col_ct; col_idx++) {
      dest[col_idx] += (*tm_ptr++) * vect_val;
    }
  }
}

static inline void mult_matrix_dxn_vect_n(const float* mm, const float* vect, float* dest, uint32_t col_ct, uint32_t row_ct) {
  uintptr_t col_cta4 = (col_ct + 3) & (~3);
  const float* mm_ptr;
  const float* vect_ptr;
  uint32_t row_idx;
  uint32_t col_idx;
  float fxx;
  for (row_idx = 0; row_idx < row_ct; row_idx++) {
    fxx = 0.0;
    vect_ptr = vect;
    mm_ptr = &(mm[row_idx * col_cta4]);
    for (col_idx = 0; col_idx < col_ct; col_idx++) {
      fxx += (*mm_ptr++) * (*vect_ptr++);
    }
    *dest++ = fxx;
  }
}

static inline float triple_product(const float* v1, const float* v2, const float* v3, uint32_t nn) {
  float fxx = 0.0;
  uint32_t uii;
  for (uii = 0; uii < nn; uii++) {
    fxx += (*v1++) * (*v2++) * (*v3++);
  }
  return fxx;
}

static inline void compute_two_diag_triple_product(const float* aa, const float* bb, const float* vv, float* raa_ptr, float* rab_ptr, float* rbb_ptr, uint32_t nn) {
  float raa = 0.0;
  float rab = 0.0;
  float rbb = 0.0;
  uint32_t uii;
  float fxx;
  float fyy;
  float fzz;
  for (uii = 0; uii < nn; uii++) {
    fxx = (*aa++);
    fyy = (*bb++);
    fzz = (*vv++);
    raa += fxx * fxx * fzz;
    fzz *= fyy;
    rab += fxx * fzz;
    rbb += fyy * fzz;
  }
  *raa_ptr = raa;
  *rab_ptr = rab;
  *rbb_ptr = rbb;
}

static inline void compute_three_triple_product(const float* bb, const float* a1, const float* a2, const float* a3, const float* vv, float* r1_ptr, float* r2_ptr, float* r3_ptr, uint32_t nn) {
  float r1 = 0.0;
  float r2 = 0.0;
  float r3 = 0.0;
  uint32_t uii;
  float fxx;
  for (uii = 0; uii < nn; uii++) {
    fxx = (*bb++) * (*vv++);
    r1 += (*a1++) * fxx;
    r2 += (*a2++) * fxx;
    r3 += (*a3++) * fxx;
  }
  *r1_ptr = r1;
  *r2_ptr = r2;
  *r3_ptr = r3;
}

static inline void compute_two_plus_one_triple_product(const float* bb, const float* a1, const float* a2, const float* vv, float* r1_ptr, float* r2_ptr, float* r3_ptr, uint32_t nn) {
  float r1 = 0.0;
  float r2 = 0.0;
  float r3 = 0.0;
  uint32_t uii;
  float fxx;
  float fyy;
  for (uii = 0; uii < nn; uii++) {
    fxx = (*bb++);
    fyy = fxx * (*vv++);
    r1 += fxx * fyy;
    r2 += (*a1++) * fyy;
    r3 += (*a2++) * fyy;
  }
  *r1_ptr = r1;
  *r2_ptr = r2;
  *r3_ptr = r3;
}
#endif

static inline void compute_hessian(const float* mm, const float* vv, float* dest, uint32_t col_ct, uint32_t row_ct) {
  uintptr_t col_cta4 = (col_ct + 3) & (~3);
  uintptr_t row_cta4 = (row_ct + 3) & (~3);
  uintptr_t row_cta4p1 = row_cta4 + 1;
  const float* mm_cur;
  uint32_t row_ctm3;
  uint32_t row_idx;
  uint32_t row_idx2;
  if (row_ct >= 3) {
    row_ctm3 = row_ct - 3;
    for (row_idx = 0; row_idx < row_ctm3; row_idx += 3) {
      mm_cur = &(mm[row_idx * col_cta4]);
      compute_two_diag_triple_product(mm_cur, &(mm_cur[col_cta4]), vv, &(dest[row_idx * row_cta4p1]), &(dest[(row_idx + 1) * row_cta4p1 - 1]), &(dest[(row_idx + 1) * row_cta4p1]), col_ct);
      compute_two_plus_one_triple_product(&(mm_cur[2 * col_cta4]), &(mm_cur[col_cta4]), mm_cur, vv, &(dest[(row_idx + 2) * row_cta4p1]), &(dest[(row_idx + 2) * row_cta4p1 - 1]), &(dest[(row_idx + 2) * row_cta4p1 - 2]), col_ct);
      for (row_idx2 = row_idx + 3; row_idx2 < row_ct; row_idx2++) {
        compute_three_triple_product(&(mm[row_idx2 * col_cta4]), mm_cur, &(mm_cur[col_cta4]), &(mm_cur[2 * col_cta4]), vv, &(dest[row_idx2 * row_cta4 + row_idx]), &(dest[row_idx2 * row_cta4 + row_idx + 1]), &(dest[row_idx2 * row_cta4 + row_idx + 2]), col_ct);
      }
    }
  }
  switch (row_ct % 3) {
  case 0:
    compute_two_plus_one_triple_product(&(mm[(row_ct - 3) * col_cta4]), &(mm[(row_ct - 2) * col_cta4]), &(mm[(row_ct - 1) * col_cta4]), vv, &(dest[(row_ct - 3) * row_cta4p1]), &(dest[(row_ct - 2) * row_cta4p1 - 1]), &(dest[(row_ct - 1) * row_cta4p1 - 2]), col_ct);
    // fall through
  case 2:
    compute_two_diag_triple_product(&(mm[(row_ct - 2) * col_cta4]), &(mm[(row_ct - 1) * col_cta4]), vv, &(dest[(row_ct - 2) * row_cta4p1]), &(dest[(row_ct - 1) * row_cta4p1 - 1]), &(dest[(row_ct - 1) * row_cta4p1]), col_ct);
    break;
  case 1:
    dest[(row_ct - 1) * row_cta4p1] = triple_product(&(mm[(row_ct - 1) * col_cta4]), &(mm[(row_ct - 1) * col_cta4]), vv, col_ct);
  }
}

void solve_linear_system(const float* ll, const float* yy, float* xx, uint32_t dd) {
  // if we're ever able to produce 32-bit Linux builds with statically linked
  // LAPACK, we might want to use it in place of this hardcoded stuff
  uintptr_t dim_cta4 = (dd + 3) & (~3);
  const float* ll_ptr;
  float* xx_ptr;
  uint32_t row_idx;
  uint32_t col_idx;
  float fxx;
  for (row_idx = 0; row_idx < dd; row_idx++) {
    fxx = yy[row_idx];
    ll_ptr = &(ll[row_idx * dim_cta4]);
    xx_ptr = xx;
    for (col_idx = 0; col_idx < row_idx; col_idx++) {
      fxx -= (*ll_ptr++) * (*xx_ptr++);
    }
    *xx_ptr = fxx / (*ll_ptr);
  }
  for (col_idx = dd; col_idx;) {
    fxx = xx[--col_idx];
    xx_ptr = &(xx[dd - 1]);
    for (row_idx = dd - 1; row_idx > col_idx; row_idx--) {
      fxx -= ll[row_idx * dim_cta4 + col_idx] * (*xx_ptr--);
    }
    *xx_ptr = fxx / ll[row_idx * (dim_cta4 + 1)];
  }
}

float compute_wald(const float* ll, uint32_t dd, float* xbuf) {
  uintptr_t dim_cta4 = (dd + 3) & (~3);
  uint32_t row_idx = 0;
  const float* ll_ptr;
  float* xbuf_ptr;
  uint32_t col_idx;
  float fxx;
  float fyy;
  while (1) {
    fxx = 0.0;
    ll_ptr = &(ll[row_idx * dim_cta4]);
    xbuf_ptr = xbuf;
    for (col_idx = 0; col_idx < row_idx; col_idx++) {
      fxx -= (*ll_ptr++) * (*xbuf_ptr++);
    }
    if (++row_idx == dd) {
      fyy = (*ll_ptr);
      return (fxx + 1.0) / (fyy * fyy);
    }
    *xbuf_ptr = fxx / (*ll_ptr);
  }
}

void cholesky_decomposition(const float* aa, float* ll, uint32_t dd) {
  uintptr_t dim_cta4 = (dd + 3) & (~3);
  uintptr_t dim_cta4p1 = dim_cta4 + 1;
  float* ll_ptr;
  float* ll_ptr2;
  uint32_t row_idx;
  uint32_t row_idx2;
  uint32_t col_idx;
  float fxx;
  float fyy;
  for (row_idx = 0; row_idx < dd; row_idx++) {
    fxx = aa[row_idx * dim_cta4p1];
    ll_ptr = &(ll[row_idx * dim_cta4]);
    for (col_idx = 0; col_idx < row_idx; col_idx++) {
      fyy = (*ll_ptr++);
      fxx -= fyy * fyy;
    }
    if (fxx >= 0.0) {
      fyy = sqrt(fxx);
    } else {
      fyy = 1e-6;
    }
    ll[row_idx * dim_cta4p1] = fyy;
    fyy = 1.0 / fyy; // now 1.0 / L[j][j]
    for (row_idx2 = row_idx + 1; row_idx2 < dd; row_idx2++) {
      fxx = aa[row_idx2 * dim_cta4 + row_idx];
      ll_ptr = &(ll[row_idx * dim_cta4]);
      ll_ptr2 = &(ll[row_idx2 * dim_cta4]);
      for (col_idx = 0; col_idx < row_idx; col_idx++) {
        fxx -= (*ll_ptr++) * (*ll_ptr2++);
      }
      ll[row_idx2 * dim_cta4 + row_idx] = fxx * fyy;
    }
  }
}

uint32_t logistic_regression(uint32_t indiv_ct, uint32_t param_ct, float* vv, float* hh, float* grad, float* ll, float* dcoef, const float* xx, const float* yy, float* coef, float* pp) {
  // Similar to first part of logistic.cpp fitLM(), but incorporates changes
  // from Pascal Pons et al.'s TopCoder code.
  //
  // Preallocated buffers (initial contents irrelevant):
  // vv    = sample variance buffer
  // hh    = hessian matrix buffer, param_ct x param_ct, rows 16-byte aligned
  // grad  = gradient buffer (length param_ct)
  // ll    = linear system buffer, param_ct x param_ct, rows 16-byte aligned
  // dcoef = current coefficient change buffer (length param_ct)
  // 
  // Inputs:
  // xx    = covariate (and usually genotype) matrix, covariate-major, rows are
  //         16-byte aligned, trailing row elements must be zeroed out
  // yy    = case/control phenotype
  //
  // Input/output:
  // coef  = starting point, overwritten with logistic regression result.  Must
  //         be 16-byte aligned.
  //
  // Output:
  // pp    = final likelihoods minus Y[]
  //
  // Returns 0 on success, 1 on convergence failure.
  uintptr_t param_cta4 = (param_ct + 3) & (~3);
  uint32_t iteration = 0;
  float min_delta_coef = 1e9;
  float delta_coef;
  float fxx;
  uint32_t param_idx;

  fill_float_zero(ll, param_ct * param_cta4);
  while (1) {
    iteration++;

    // P[i] = \sum_j coef[j] * X[i][j];
    mult_tmatrix_nxd_vect_d(xx, coef, pp, indiv_ct, param_ct);

    // P[i] = 1 / (1 + exp(-P[i]));
    logistic_sse(pp, indiv_ct);

    // V[i] = P[i] * (1 - P[i]);
    // P[i] -= Y[i];
    compute_v_and_p_minus_y(pp, vv, yy, indiv_ct);

    compute_hessian(xx, vv, hh, indiv_ct, param_ct);

    mult_matrix_dxn_vect_n(xx, pp, grad, indiv_ct, param_ct);

    cholesky_decomposition(hh, ll, param_ct);

    fill_float_zero(dcoef, param_ct);
    solve_linear_system(ll, grad, dcoef, param_ct);

    delta_coef = 0.0;
    for (param_idx = 0; param_idx < param_ct; param_idx++) {
      fxx = dcoef[param_idx];
      delta_coef += fabsf(fxx);
      coef[param_idx] -= fxx;
    }
    if (delta_coef < min_delta_coef) {
      min_delta_coef = delta_coef;
    }
    if (delta_coef != delta_coef) {
      return 1;
    }
    if (iteration > 4) {
      if (((delta_coef > 20.0) && (delta_coef > 2 * min_delta_coef)) || ((iteration >= 8) && fabsf(1.0 - delta_coef) < 1e-3)) {
	return 1;
      }
      if (iteration >= 15) {
	return 0;
      }
    }
    // Pons reported that 1.1e-3 was dangerous, so I agree with the decision to
    // tighten this threshold from 1e-3 to 1e-4.
    if (delta_coef < 1e-4) {
      return 0;
    }
  }
}

uint32_t glm_logistic_robust_cluster_covar(uintptr_t cur_batch_size, uintptr_t param_ct, uintptr_t indiv_valid_ct, uint32_t missing_ct, uintptr_t* loadbuf, float* covars_cov_major, float* covars_indiv_major, uintptr_t* perm_vecs, float* coef, float* pp, float* indiv_1d_buf, float* pheno_buf, float* param_1d_buf, float* param_1d_buf2, float* param_2d_buf, float* param_2d_buf2, float* logistic_results, uint32_t cluster_ct1, uint32_t* indiv_to_cluster1, float* cluster_param_buf, float* cluster_param_buf2, uintptr_t constraint_ct, double* constraints_con_major, double* param_1d_dbuf, double* param_2d_dbuf, double* param_2d_dbuf2, double* param_df_dbuf, double* df_df_dbuf, MATRIX_INVERT_BUF1_TYPE* mi_buf, double* df_dbuf, uintptr_t* perm_fails) {
  // Similar to logistic.cpp fitLM(), but incorporates changes from the
  // postprocessed TopCoder contest code.
  // * coef is now assumed to be initialized with a good starting point for
  //   each permutation.  The initializer must remove columns corresponding to
  //   missing genotypes.  (Todo: figure out when/whether we actually want to
  //   take advantage of this with a nonzero starting point.  If we do, we'll
  //   also need to add restart logic.)
  // * covars_cov_major must now have 16-byte aligned rows.
  // Returns number of regression failures.
  uintptr_t param_cta4 = (param_ct + 3) & (~3);
  uintptr_t param_ct_p1 = param_ct + 1;
  uintptr_t param_ct_m1 = param_ct - 1;
  uintptr_t joint_test_requested = (constraints_con_major? 1 : 0);
  uintptr_t param_ctx = param_ct + joint_test_requested;
  uintptr_t param_ctx_m1 = param_ctx - 1;
  uintptr_t indiv_validx_ctv2 = 2 * ((indiv_valid_ct + missing_ct + (BITCT - 1)) / BITCT);
  uintptr_t perm_fail_ct = 0;
  uintptr_t cur_word = 0;
  uint32_t cluster_ct1_p1 = cluster_ct1 + 1;
  uintptr_t perm_idx;
  uintptr_t indiv_uidx;
  uintptr_t indiv_idx;
  uintptr_t param_idx;
  uintptr_t param_idx2;
  float* fptr;
  float* fptr2;

  double dxx;

  // +1 since "no cluster" is represented as -1
  // 32-bit since that -1 is actually 0xffffffffU
  uint32_t cluster_idx_p1;
  float fxx;

  fill_ulong_zero(perm_fails, ((cur_batch_size + (BITCT - 1)) / BITCT) * sizeof(intptr_t));
  for (perm_idx = 0; perm_idx < cur_batch_size; perm_idx++) {
    fptr = pheno_buf;
    if (!missing_ct) {
      for (indiv_idx = 0; indiv_idx < indiv_valid_ct; indiv_idx++) {
	// strictly speaking, we can use 1 bit per permutation instead of 2
	// bits here, but the gain is probably too small to justify even adding
	// a parameter to generate_cc_[cluster_]perm_vec.
	*fptr++ = (float)((int32_t)is_set_ul(perm_vecs, indiv_idx * 2));
      }
    } else {
      for (indiv_uidx = 0, indiv_idx = 0; indiv_idx < indiv_valid_ct; indiv_uidx++, indiv_idx++) {
        while (1) {
	  if (indiv_uidx % BITCT2) {
            cur_word >>= 2;
	  } else {
	    cur_word = loadbuf[indiv_uidx / BITCT2];
	    cur_word = cur_word & (~(cur_word >> 1)) & FIVEMASK;
	  }
	  if (!(cur_word & 1)) {
	    break;
	  }
	  indiv_uidx++;
	}
        *fptr++ = (float)((int32_t)is_set_ul(perm_vecs, indiv_uidx * 2));
      }
    }
    if (logistic_regression(indiv_valid_ct, param_ct, indiv_1d_buf, param_2d_buf, param_1d_buf, param_2d_buf2, param_1d_buf2, covars_cov_major, pheno_buf, coef, pp)) {
      goto glm_logistic_robust_cluster_covar_fail;
    }

    // compute S
    // param_2d_buf = S, param_1d_buf = y, param_1d_buf2 = x
    for (param_idx = 0; param_idx < param_ct; param_idx++) {
      fill_float_zero(param_1d_buf, param_ct);
      param_1d_buf[param_idx] = 1.0;
      solve_linear_system(param_2d_buf2, param_1d_buf, param_1d_buf2, param_ct);
      // S does *not* currently have 16-byte aligned rows
      memcpy(&(param_2d_buf[param_idx * param_ct]), param_1d_buf2, param_ct * sizeof(float));
    }

    if (cluster_ct1) {
      // HuberWhite()
      fill_float_zero(cluster_param_buf, cluster_ct1_p1 * param_ct);
      if (!missing_ct) {
	for (indiv_idx = 0; indiv_idx < indiv_valid_ct; indiv_idx++) {
          cluster_idx_p1 = indiv_to_cluster1[indiv_idx] + 1;
          fxx = -pp[indiv_idx];
	  fptr = &(cluster_param_buf[cluster_idx_p1 * param_ct]);
	  fptr2 = &(covars_indiv_major[indiv_idx * param_ct]);
	  for (param_idx = 0; param_idx < param_ct; param_idx++) {
	    *fptr += fxx * (*fptr2++);
	    fptr++;
	  }
	}
      } else {
	for (indiv_uidx = 0, indiv_idx = 0; indiv_idx < indiv_valid_ct; indiv_uidx++, indiv_idx++) {
          cluster_idx_p1 = indiv_to_cluster1[indiv_idx] + 1;
	  while (1) {
	    if (indiv_uidx % BITCT2) {
              cur_word >>= 2;
	    } else {
	      cur_word = loadbuf[indiv_uidx / BITCT2];
              cur_word = cur_word & (~(cur_word >> 1)) & FIVEMASK;
	    }
	    if (!(cur_word & 1)) {
	      break;
	    }
	    indiv_uidx++;
	  }
          fxx = -pp[indiv_idx];
          fptr = &(cluster_param_buf[cluster_idx_p1 * param_ct]);
	  fptr2 = &(covars_indiv_major[indiv_idx * param_ct]);
	  for (param_idx = 0; param_idx < param_ct; param_idx++) {
            *fptr += fxx * (*fptr2++);
            fptr++;
	  }
	}
      }
      transpose_copy_float(cluster_ct1_p1, param_ct, param_ct, cluster_param_buf, cluster_param_buf2);
      col_major_fmatrix_multiply(param_ct, param_ct, cluster_ct1_p1, cluster_param_buf, cluster_param_buf2, param_2d_buf2);
      col_major_fmatrix_multiply(param_ct, param_ct, param_ct, param_2d_buf, param_2d_buf2, cluster_param_buf);
      col_major_fmatrix_multiply(param_ct, param_ct, param_ct, cluster_param_buf, param_2d_buf, param_2d_buf2);
      memcpy(param_2d_buf, param_2d_buf2, param_ct * param_ct * sizeof(float));
    }
    // validParameters() check
    for (param_idx = 1; param_idx < param_ct; param_idx++) {
      fxx = param_2d_buf[param_idx * param_ct_p1];
      if ((fxx < 1e-20) || (!realnum(fxx))) {
	goto glm_logistic_robust_cluster_covar_fail;
      }
      param_2d_buf2[param_idx] = sqrtf(fxx);
    }
    param_2d_buf2[0] = sqrtf(param_2d_buf[0]);
    for (param_idx = 1; param_idx < param_ct; param_idx++) {
      fxx = 0.99999 * param_2d_buf2[param_idx];
      fptr = &(param_2d_buf[param_idx * param_ct]);
      fptr2 = param_2d_buf2;
      for (param_idx2 = 0; param_idx2 < param_idx; param_idx2++) {
        if ((*fptr++) > fxx * (*fptr2++)) {
	  goto glm_logistic_robust_cluster_covar_fail;
	}
      }
    }
    fptr = &(logistic_results[perm_idx * param_ctx_m1]);
    for (param_idx = 1; param_idx < param_ct; param_idx++) {
      *fptr++ = param_2d_buf[param_idx * param_ct_p1];
    }
    if (joint_test_requested) {
      // ugly to do this with doubles while doing everything else with floats,
      // but unless I need to write a float matrix inversion later for some
      // other reason, it's probably still the least-bad approach
      if (!constraint_ct) {
	*fptr++ = -9;
      } else {
	for (param_idx = 0; param_idx < param_ct; param_idx++) {
	  param_1d_dbuf[param_idx] = (double)coef[param_idx];
	}
	param_idx2 = param_ct * param_ct;
	for (param_idx = 0; param_idx < param_idx2; param_idx++) {
	  param_2d_dbuf[param_idx] = (double)(param_2d_buf[param_idx]);
	}
	if (!linear_hypothesis_chisq(constraint_ct, param_ct, constraints_con_major, param_1d_dbuf, param_2d_dbuf, param_2d_dbuf2, param_df_dbuf, df_df_dbuf, mi_buf, df_dbuf, &dxx)) {
	  *fptr++ = (float)dxx;
	} else {
	  *fptr++ = -9;
	}
      }
    }
    if (0) {
    glm_logistic_robust_cluster_covar_fail:
      fill_float_zero(&(logistic_results[perm_idx * param_ctx_m1]), param_ct_m1);
      SET_BIT(perm_fails, perm_idx);
      perm_fail_ct++;
      if (joint_test_requested) {
        logistic_results[perm_idx * param_ctx_m1 + param_ct_m1] = -9;
      }
    }
    coef = &(coef[param_cta4]);
    perm_vecs = &(perm_vecs[indiv_validx_ctv2]);
  }
  return perm_fail_ct;
}

uint32_t glm_fill_design(uintptr_t* loadbuf_collapsed, double* fixed_covars_cov_major, uintptr_t indiv_valid_ct, uint32_t* indiv_to_cluster1, uint32_t cur_param_ct, uint32_t coding_flags, uint32_t xchr_model, uintptr_t condition_list_start_idx, uintptr_t interaction_start_idx, uintptr_t sex_start_idx, uintptr_t* active_params, uintptr_t* haploid_params, uint32_t include_sex, uint32_t male_x_01, uintptr_t* sex_male_collapsed, uint32_t is_nonx_haploid, double* cur_covars_cov_major, double* cur_covars_indiv_major, uint32_t* cur_indiv_to_cluster1_buf, uint32_t** cur_indiv_to_cluster1_ptr, uint32_t standard_beta) {
  double* dptr = cur_covars_cov_major;
  uintptr_t* ulptr_end_init = &(loadbuf_collapsed[indiv_valid_ct / BITCT2]);
  uintptr_t fixed_covar_nonsex_ct = interaction_start_idx - condition_list_start_idx;
  uintptr_t interactions_present = sex_start_idx - interaction_start_idx;
  uint32_t genotypic_or_hethom = condition_list_start_idx - 2;
  uint32_t missing_ct = 0;
  double* dptr2;
  uintptr_t* ulptr;
  uintptr_t* ulptr_end;
  uintptr_t fixed_covar_idx;
  uintptr_t cur_indiv_valid_ct;
  uintptr_t indiv_idx;
  uintptr_t indiv_idx_stop;
  uintptr_t cur_word;
  uintptr_t cur_genotype;
  uintptr_t param_idx;
  double dxx;
  double dyy;
  double dzz;
  // don't need to recompute this during permutations, but it's so cheap that
  // it hardly matters
  missing_ct = count_01(loadbuf_collapsed, (indiv_valid_ct + BITCT2 - 1) / BITCT2);
  if (missing_ct >= indiv_valid_ct - 1) {
    // regression will be skipped in this case
    return missing_ct;
  }
  cur_indiv_valid_ct = indiv_valid_ct - missing_ct;
  for (indiv_idx = 0; indiv_idx < cur_indiv_valid_ct; indiv_idx++) {
    *dptr++ = 1;
  }
  if (IS_SET(active_params, 1)) {
    ulptr = loadbuf_collapsed;
    ulptr_end = ulptr_end_init;
    indiv_idx = 0;
    indiv_idx_stop = BITCT2;
    if (coding_flags & GLM_DOMINANT) {
      while (1) {
	while (ulptr < ulptr_end) {
	  cur_word = *ulptr++;
	  for (; indiv_idx < indiv_idx_stop; indiv_idx++, cur_word >>= 2) {
	    cur_genotype = cur_word & 3;
	    if (cur_genotype != 1) {
	      // 0/1/1, i.e.
	      //   3 = hom A2 -> 0
	      //   2 = het -> 1
	      //   0 = hom A1 -> 1

	      // Alexandrescu states that non-branching comparisons are
	      // actually *lower*-strength than arithmetic operations, so I'll
	      // try this...
	      *dptr++ = (double)(cur_genotype != 3);
	    }
	  }
	  indiv_idx_stop += BITCT2;
	}
	if (indiv_idx == indiv_valid_ct) {
	  break;
	}
	ulptr_end++;
	indiv_idx_stop = indiv_valid_ct;
      }      
    } else if ((!(coding_flags & (GLM_HETHOM | GLM_RECESSIVE))) && (!is_nonx_haploid)) {
      if (!male_x_01) {
	while (1) {
	  while (ulptr < ulptr_end) {
	    cur_word = *ulptr++;
	    for (; indiv_idx < indiv_idx_stop; indiv_idx++, cur_word >>= 2) {
	      cur_genotype = cur_word & 3;
	      if (cur_genotype != 1) {
		// 0/1/2, i.e.
		//   3 = hom A2 -> 0
		//   2 = het -> 1
		//   0 = hom A1 -> 2
		*dptr++ = (double)((intptr_t)(2 + (cur_genotype / 2) - cur_genotype));
	      }
	    }
	    indiv_idx_stop += BITCT2;
	  }
	  if (indiv_idx == indiv_valid_ct) {
	    break;
	  }
	  ulptr_end++;
	  indiv_idx_stop = indiv_valid_ct;
	}
      } else {
	while (1) {
	  while (ulptr < ulptr_end) {
	    cur_word = *ulptr++;
	    for (; indiv_idx < indiv_idx_stop; indiv_idx++, cur_word >>= 2) {
	      cur_genotype = cur_word & 3;
	      if (cur_genotype != 1) {
		// 0/1/2, but downshifted for males
		*dptr++ = (double)((intptr_t)((2 + (cur_genotype / 2) - cur_genotype) >> IS_SET(sex_male_collapsed, indiv_idx)));
	      }
	    }
	    indiv_idx_stop += BITCT2;
	  }
	  if (indiv_idx == indiv_valid_ct) {
	    break;
	  }
	  ulptr_end++;
	  indiv_idx_stop = indiv_valid_ct;
	}
      }
    } else {
      while (1) {
	while (ulptr < ulptr_end) {
	  cur_word = *ulptr++;
	  for (; indiv_idx < indiv_idx_stop; indiv_idx++, cur_word >>= 2) {
	    cur_genotype = cur_word & 3;
	    if (cur_genotype != 1) {
	      // 0/0/1
	      *dptr++ = (double)(cur_genotype == 0);
	    }
	  }
          indiv_idx_stop += BITCT2;
	}
        if (indiv_idx == indiv_valid_ct) {
          break;
	}
        ulptr_end++;
	indiv_idx_stop = indiv_valid_ct;
      }
    }
  }
  if (genotypic_or_hethom && (!is_nonx_haploid) && IS_SET(active_params, 2)) {
    ulptr = loadbuf_collapsed;
    ulptr_end = ulptr_end_init;
    indiv_idx = 0;
    indiv_idx_stop = BITCT2;
    while (1) {
      while (ulptr < ulptr_end) {
	cur_word = *ulptr++;
	for (; indiv_idx < indiv_idx_stop; indiv_idx++, cur_word >>= 2) {
	  cur_genotype = cur_word & 3;
	  if (cur_genotype != 1) {
	    // 0/1/0
            *dptr++ = (double)(cur_genotype == 2);
	  }
	}
	indiv_idx_stop += BITCT2;
      }
      if (indiv_idx == indiv_valid_ct) {
	break;
      }
      ulptr_end++;
      indiv_idx_stop = indiv_valid_ct;
    }
  }
  for (fixed_covar_idx = 0; fixed_covar_idx < fixed_covar_nonsex_ct; fixed_covar_idx++) {
    if (!is_set(active_params, fixed_covar_idx + condition_list_start_idx)) {
      continue;
    }
    copy_when_nonmissing(loadbuf_collapsed, (char*)(&(fixed_covars_cov_major[fixed_covar_idx * indiv_valid_ct])), sizeof(double), indiv_valid_ct, missing_ct, (char*)dptr);
    dptr = &(dptr[cur_indiv_valid_ct]);
  }
  if (interactions_present) {
    param_idx = interaction_start_idx;
    for (fixed_covar_idx = 0; fixed_covar_idx < fixed_covar_nonsex_ct; fixed_covar_idx++, param_idx++) {
      if (IS_SET(active_params, param_idx)) {
	ulptr = loadbuf_collapsed;
	ulptr_end = ulptr_end_init;
	indiv_idx = 0;
	indiv_idx_stop = BITCT2;
	dptr2 = &(fixed_covars_cov_major[fixed_covar_idx * indiv_valid_ct]);
        if (coding_flags & GLM_DOMINANT) {
	  while (1) {
	    while (ulptr < ulptr_end) {
	      cur_word = *ulptr++;
	      for (; indiv_idx < indiv_idx_stop; indiv_idx++, cur_word >>= 2) {
		cur_genotype = cur_word & 3;
		if (cur_genotype != 1) {
		  // 0/1/1
		  *dptr++ = ((double)(cur_genotype != 3)) * (*dptr2);
		}
		dptr2++;
	      }
	      indiv_idx_stop += BITCT2;
	    }
	    if (indiv_idx == indiv_valid_ct) {
	      break;
	    }
	    ulptr_end++;
	    indiv_idx_stop = indiv_valid_ct;
	  }
	} else if (!(coding_flags & (GLM_HETHOM | GLM_RECESSIVE))) {
	  if (!male_x_01) {
	    while (1) {
	      while (ulptr < ulptr_end) {
		cur_word = *ulptr++;
		for (; indiv_idx < indiv_idx_stop; indiv_idx++, cur_word >>= 2) {
		  cur_genotype = cur_word & 3;
		  if (cur_genotype != 1) {
		    // 0/1/2
		    *dptr++ = ((double)((intptr_t)(2 + (cur_genotype / 2) - cur_genotype))) * (*dptr2);
		  }
		  dptr2++;
		}
		indiv_idx_stop += BITCT2;
	      }
	      if (indiv_idx == indiv_valid_ct) {
		break;
	      }
	      ulptr_end++;
	      indiv_idx_stop = indiv_valid_ct;
	    }
	  } else {
	    while (1) {
	      while (ulptr < ulptr_end) {
		cur_word = *ulptr++;
		for (; indiv_idx < indiv_idx_stop; indiv_idx++, cur_word >>= 2) {
		  cur_genotype = cur_word & 3;
		  if (cur_genotype != 1) {
		    // 0/1/2
		    *dptr++ = ((double)((intptr_t)((2 + (cur_genotype / 2) - cur_genotype) >> IS_SET(sex_male_collapsed, indiv_idx)))) * (*dptr2);
		  }
		  dptr2++;
		}
		indiv_idx_stop += BITCT2;
	      }
	      if (indiv_idx == indiv_valid_ct) {
		break;
	      }
	      ulptr_end++;
	      indiv_idx_stop = indiv_valid_ct;
	    }
	  }
	} else {
	  while (1) {
	    while (ulptr < ulptr_end) {
	      cur_word = *ulptr++;
	      for (; indiv_idx < indiv_idx_stop; indiv_idx++, cur_word >>= 2) {
		cur_genotype = cur_word & 3;
		if (cur_genotype != 1) {
		  // 0/0/1
	          *dptr++ = ((double)(cur_genotype == 0)) * (*dptr2);
		}
		dptr2++;
	      }
	      indiv_idx_stop += BITCT2;
	    }
	    if (indiv_idx == indiv_valid_ct) {
	      break;
	    }
	    ulptr_end++;
	    indiv_idx_stop = indiv_valid_ct;
	  }
	}
      }
      if (genotypic_or_hethom) {
	param_idx++;
	if ((!is_nonx_haploid) && IS_SET(active_params, param_idx)) {
	  ulptr = loadbuf_collapsed;
	  ulptr_end = ulptr_end_init;
	  indiv_idx = 0;
	  indiv_idx_stop = BITCT2;
	  dptr2 = &(fixed_covars_cov_major[fixed_covar_idx * indiv_valid_ct]);
	  while (1) {
	    while (ulptr < ulptr_end) {
	      cur_word = *ulptr++;
	      for (; indiv_idx < indiv_idx_stop; indiv_idx++, cur_word >>= 2) {
		cur_genotype = cur_word & 3;
		if (cur_genotype != 1) {
		  // 0/1/0
                  *dptr++ = ((double)(cur_genotype == 2)) * (*dptr2);
		}
		dptr2++;
	      }
	      indiv_idx_stop += BITCT2;
	    }
	    if (indiv_idx == indiv_valid_ct) {
	      break;
	    }
	    ulptr_end++;
	    indiv_idx_stop = indiv_valid_ct;
	  }
	}
      }
    }
  }
  if (include_sex) {
    if (IS_SET(active_params, sex_start_idx)) {
      copy_when_nonmissing(loadbuf_collapsed, (char*)(&(fixed_covars_cov_major[fixed_covar_nonsex_ct * indiv_valid_ct])), sizeof(double), indiv_valid_ct, missing_ct, (char*)dptr);
      dptr = &(dptr[cur_indiv_valid_ct]);
    }
    if (interactions_present) {
      if (is_set(active_params, sex_start_idx + 1)) {
	ulptr = loadbuf_collapsed;
	ulptr_end = ulptr_end_init;
	indiv_idx = 0;
	indiv_idx_stop = BITCT2;
	dptr2 = &(fixed_covars_cov_major[fixed_covar_nonsex_ct * indiv_valid_ct]);
	if (!(coding_flags & GLM_HETHOM)) {
	  if (!male_x_01) {
	    while (1) {
	      while (ulptr < ulptr_end) {
		cur_word = *ulptr++;
		for (; indiv_idx < indiv_idx_stop; indiv_idx++, cur_word >>= 2) {
		  cur_genotype = cur_word & 3;
		  if (cur_genotype != 1) {
		    // 0/1/2
		    *dptr++ = ((double)((intptr_t)(2 + (cur_genotype / 2) - cur_genotype))) * (*dptr2);
		  }
		  dptr2++;
		}
		indiv_idx_stop += BITCT2;
	      }
	      if (indiv_idx == indiv_valid_ct) {
		break;
	      }
	      ulptr_end++;
	      indiv_idx_stop = indiv_valid_ct;
	    }
	  } else {
	    while (1) {
	      while (ulptr < ulptr_end) {
		cur_word = *ulptr++;
		for (; indiv_idx < indiv_idx_stop; indiv_idx++, cur_word >>= 2) {
		  cur_genotype = cur_word & 3;
		  if (cur_genotype != 1) {
		    // 0/1/2
		    *dptr++ = ((double)((intptr_t)((2 + (cur_genotype / 2) - cur_genotype) >> IS_SET(sex_male_collapsed, indiv_idx)))) * (*dptr2);
		  }
		  dptr2++;
		}
		indiv_idx_stop += BITCT2;
	      }
	      if (indiv_idx == indiv_valid_ct) {
		break;
	      }
	      ulptr_end++;
	      indiv_idx_stop = indiv_valid_ct;
	    }
	  }
	} else {
	  while (1) {
	    while (ulptr < ulptr_end) {
	      cur_word = *ulptr++;
	      for (; indiv_idx < indiv_idx_stop; indiv_idx++, cur_word >>= 2) {
		cur_genotype = cur_word & 3;
		if (cur_genotype != 1) {
		  // 0/0/1
	          *dptr++ = ((double)(cur_genotype == 0)) * (*dptr2);
		}
		dptr2++;
	      }
	      indiv_idx_stop += BITCT2;
	    }
	    if (indiv_idx == indiv_valid_ct) {
	      break;
	    }
	    ulptr_end++;
	    indiv_idx_stop = indiv_valid_ct;
	  }
	}
      }
      if (genotypic_or_hethom && (!is_nonx_haploid) && is_set(active_params, sex_start_idx + 2)) {
	ulptr = loadbuf_collapsed;
	ulptr_end = ulptr_end_init;
	indiv_idx = 0;
	indiv_idx_stop = BITCT2;
	dptr2 = &(fixed_covars_cov_major[fixed_covar_nonsex_ct * indiv_valid_ct]);
	while (1) {
	  while (ulptr < ulptr_end) {
	    cur_word = *ulptr++;
	    for (; indiv_idx < indiv_idx_stop; indiv_idx++, cur_word >>= 2) {
	      cur_genotype = cur_word & 3;
	      if (cur_genotype != 1) {
		// 0/1/0
		*dptr++ = ((double)(cur_genotype == 2)) * (*dptr2);
	      }
	      dptr2++;
	    }
	    indiv_idx_stop += BITCT2;
	  }
	  if (indiv_idx == indiv_valid_ct) {
	    break;
	  }
	  ulptr_end++;
	  indiv_idx_stop = indiv_valid_ct;
	}
      }
    }
  }
  // if (dptr != &(cur_covars_cov_major[cur_param_ct * cur_indiv_valid_ct])) {
  //   printf("assert failure:\n  cur_param_ct = %u\n  cur_indiv_valid_ct = %" PRIuPTR "\n  dptr - cur_covars_cov_major = %" PRIuPTR "\n", cur_param_ct, cur_indiv_valid_ct, (uintptr_t)(dptr - cur_covars_cov_major));
  //   exit(1);
  // }
  if (indiv_to_cluster1) {
    if (!missing_ct) {
      *cur_indiv_to_cluster1_ptr = indiv_to_cluster1;
    } else {
      copy_when_nonmissing(loadbuf_collapsed, (char*)indiv_to_cluster1, sizeof(int32_t), indiv_valid_ct, missing_ct, (char*)cur_indiv_to_cluster1_buf);
      *cur_indiv_to_cluster1_ptr = cur_indiv_to_cluster1_buf;
    }
  }
  if (standard_beta) {
    for (param_idx = 1; param_idx < cur_param_ct; param_idx++) {
      dxx = 0; // sum
      dyy = 0; // ssq
      dptr = &(cur_covars_cov_major[param_idx * cur_indiv_valid_ct]);
      for (indiv_idx = 0; indiv_idx < cur_indiv_valid_ct; indiv_idx++) {
	dzz = *dptr++;
	dxx += dzz;
	dyy += dzz * dzz;
      }
      dptr = &(cur_covars_cov_major[param_idx * cur_indiv_valid_ct]);
      dzz = dxx / ((double)((intptr_t)cur_indiv_valid_ct));
      dyy = sqrt((dyy - dxx * dzz) / ((double)((intptr_t)(cur_indiv_valid_ct - 1))));
      if (dyy == 0) {
	fill_double_zero(dptr, cur_indiv_valid_ct);
      } else {
	dyy = 1.0 / dyy;
	for (indiv_idx = 0; indiv_idx < cur_indiv_valid_ct; indiv_idx++) {
	  *dptr = ((*dptr) - dzz) * dyy;
	  dptr++;
	}
      }
    }
  }

  transpose_copy(cur_param_ct, cur_indiv_valid_ct, cur_covars_cov_major, cur_covars_indiv_major);
  return missing_ct;
}

uint32_t glm_fill_design_float(uintptr_t* loadbuf_collapsed, float* fixed_covars_cov_major, uintptr_t indiv_valid_ct, uint32_t* indiv_to_cluster1, uint32_t cur_param_ct, uint32_t coding_flags, uint32_t xchr_model, uintptr_t condition_list_start_idx, uintptr_t interaction_start_idx, uintptr_t sex_start_idx, uintptr_t* active_params, uintptr_t* haploid_params, uint32_t include_sex, uint32_t male_x_01, uintptr_t* sex_male_collapsed, uint32_t is_nonx_haploid, float* cur_covars_cov_major, float* cur_covars_indiv_major, uint32_t* cur_indiv_to_cluster1_buf, uint32_t** cur_indiv_to_cluster1_ptr) {
  // rows of cur_covars_cov_major must be 16-byte aligned
  float* fptr = cur_covars_cov_major;
  uintptr_t* ulptr_end_init = &(loadbuf_collapsed[indiv_valid_ct / BITCT2]);
  uintptr_t fixed_covar_nonsex_ct = interaction_start_idx - condition_list_start_idx;
  uintptr_t interactions_present = sex_start_idx - interaction_start_idx;
  uint32_t genotypic_or_hethom = condition_list_start_idx - 2;
  uint32_t missing_ct = 0;
  float* fptr2;
  uintptr_t* ulptr;
  uintptr_t* ulptr_end;
  uintptr_t fixed_covar_idx;
  uintptr_t cur_indiv_valid_ct;
  uintptr_t cur_indiv_valid_cta4;
  uintptr_t indiv_idx;
  uintptr_t indiv_idx_stop;
  uintptr_t cur_word;
  uintptr_t cur_genotype;
  uintptr_t param_idx;
  uint32_t align_skip;
  // don't need to recompute this during permutations, but it's so cheap that
  // it hardly matters
  missing_ct = count_01(loadbuf_collapsed, (indiv_valid_ct + BITCT2 - 1) / BITCT2);
  if (missing_ct >= indiv_valid_ct - 1) {
    // regression will be skipped in this case
    return missing_ct;
  }
  cur_indiv_valid_ct = indiv_valid_ct - missing_ct;
  cur_indiv_valid_cta4 = (cur_indiv_valid_ct + 3) & (~3);
  align_skip = cur_indiv_valid_cta4 - cur_indiv_valid_ct;
  for (indiv_idx = 0; indiv_idx < cur_indiv_valid_ct; indiv_idx++) {
    *fptr++ = 1;
  }
  fill_float_zero(fptr, align_skip);
  fptr = &(fptr[align_skip]);
  if (IS_SET(active_params, 1)) {
    ulptr = loadbuf_collapsed;
    ulptr_end = ulptr_end_init;
    indiv_idx = 0;
    indiv_idx_stop = BITCT2;
    if (coding_flags & GLM_DOMINANT) {
      while (1) {
	while (ulptr < ulptr_end) {
	  cur_word = *ulptr++;
	  for (; indiv_idx < indiv_idx_stop; indiv_idx++, cur_word >>= 2) {
	    cur_genotype = cur_word & 3;
	    if (cur_genotype != 1) {
	      // 0/1/1, i.e.
	      //   3 = hom A2 -> 0
	      //   2 = het -> 1
	      //   0 = hom A1 -> 1
	      *fptr++ = (float)(cur_genotype != 3);
	    }
	  }
	  indiv_idx_stop += BITCT2;
	}
	if (indiv_idx == indiv_valid_ct) {
	  break;
	}
	ulptr_end++;
	indiv_idx_stop = indiv_valid_ct;
      }      
    } else if ((!(coding_flags & (GLM_HETHOM | GLM_RECESSIVE))) && (!is_nonx_haploid)) {
      if (!male_x_01) {
	while (1) {
	  while (ulptr < ulptr_end) {
	    cur_word = *ulptr++;
	    for (; indiv_idx < indiv_idx_stop; indiv_idx++, cur_word >>= 2) {
	      cur_genotype = cur_word & 3;
	      if (cur_genotype != 1) {
		// 0/1/2, i.e.
		//   3 = hom A2 -> 0
		//   2 = het -> 1
		//   0 = hom A1 -> 2
		*fptr++ = (float)((intptr_t)(2 + (cur_genotype / 2) - cur_genotype));
	      }
	    }
	    indiv_idx_stop += BITCT2;
	  }
	  if (indiv_idx == indiv_valid_ct) {
	    break;
	  }
	  ulptr_end++;
	  indiv_idx_stop = indiv_valid_ct;
	}
      } else {
	while (1) {
	  while (ulptr < ulptr_end) {
	    cur_word = *ulptr++;
	    for (; indiv_idx < indiv_idx_stop; indiv_idx++, cur_word >>= 2) {
	      cur_genotype = cur_word & 3;
	      if (cur_genotype != 1) {
		// 0/1/2, but downshifted for males
		*fptr++ = (float)((intptr_t)((2 + (cur_genotype / 2) - cur_genotype) >> IS_SET(sex_male_collapsed, indiv_idx)));
	      }
	    }
	    indiv_idx_stop += BITCT2;
	  }
	  if (indiv_idx == indiv_valid_ct) {
	    break;
	  }
	  ulptr_end++;
	  indiv_idx_stop = indiv_valid_ct;
	}
      }
    } else {
      while (1) {
	while (ulptr < ulptr_end) {
	  cur_word = *ulptr++;
	  for (; indiv_idx < indiv_idx_stop; indiv_idx++, cur_word >>= 2) {
	    cur_genotype = cur_word & 3;
	    if (cur_genotype != 1) {
	      // 0/0/1
	      *fptr++ = (float)(cur_genotype == 0);
	    }
	  }
          indiv_idx_stop += BITCT2;
	}
        if (indiv_idx == indiv_valid_ct) {
          break;
	}
        ulptr_end++;
	indiv_idx_stop = indiv_valid_ct;
      }
    }
    fill_float_zero(fptr, align_skip);
    fptr = &(fptr[align_skip]);
  }
  if (genotypic_or_hethom && (!is_nonx_haploid) && IS_SET(active_params, 2)) {
    ulptr = loadbuf_collapsed;
    ulptr_end = ulptr_end_init;
    indiv_idx = 0;
    indiv_idx_stop = BITCT2;
    while (1) {
      while (ulptr < ulptr_end) {
	cur_word = *ulptr++;
	for (; indiv_idx < indiv_idx_stop; indiv_idx++, cur_word >>= 2) {
	  cur_genotype = cur_word & 3;
	  if (cur_genotype != 1) {
	    // 0/1/0
            *fptr++ = (float)(cur_genotype == 2);
	  }
	}
	indiv_idx_stop += BITCT2;
      }
      if (indiv_idx == indiv_valid_ct) {
	break;
      }
      ulptr_end++;
      indiv_idx_stop = indiv_valid_ct;
    }
    fill_float_zero(fptr, align_skip);
    fptr = &(fptr[align_skip]);
  }
  for (fixed_covar_idx = 0; fixed_covar_idx < fixed_covar_nonsex_ct; fixed_covar_idx++) {
    if (!is_set(active_params, fixed_covar_idx + condition_list_start_idx)) {
      continue;
    }
    copy_when_nonmissing(loadbuf_collapsed, (char*)(&(fixed_covars_cov_major[fixed_covar_idx * indiv_valid_ct])), sizeof(float), indiv_valid_ct, missing_ct, (char*)fptr);
    fill_float_zero(&(fptr[cur_indiv_valid_ct]), align_skip);
    fptr = &(fptr[cur_indiv_valid_cta4]);
  }
  if (interactions_present) {
    param_idx = interaction_start_idx;
    for (fixed_covar_idx = 0; fixed_covar_idx < fixed_covar_nonsex_ct; fixed_covar_idx++, param_idx++) {
      if (IS_SET(active_params, param_idx)) {
	ulptr = loadbuf_collapsed;
	ulptr_end = ulptr_end_init;
	indiv_idx = 0;
	indiv_idx_stop = BITCT2;
	fptr2 = &(fixed_covars_cov_major[fixed_covar_idx * indiv_valid_ct]);
	if (coding_flags & GLM_DOMINANT) {
	  while (1) {
	    while (ulptr < ulptr_end) {
	      cur_word = *ulptr++;
	      for (; indiv_idx < indiv_idx_stop; indiv_idx++, cur_word >>= 2) {
		cur_genotype = cur_word & 3;
		if (cur_genotype != 1) {
		  // 0/1/1
		  *fptr++ = ((float)(cur_genotype != 3)) * (*fptr2);
		}
		fptr2++;
	      }
	      indiv_idx_stop += BITCT2;
	    }
	    if (indiv_idx == indiv_valid_ct) {
	      break;
	    }
	    ulptr_end++;
	    indiv_idx_stop = indiv_valid_ct;
	  }
	} else if (!(coding_flags & (GLM_HETHOM | GLM_RECESSIVE))) {
	  if (!male_x_01) {
	    while (1) {
	      while (ulptr < ulptr_end) {
		cur_word = *ulptr++;
		for (; indiv_idx < indiv_idx_stop; indiv_idx++, cur_word >>= 2) {
		  cur_genotype = cur_word & 3;
		  if (cur_genotype != 1) {
		    // 0/1/2
		    *fptr++ = ((float)((intptr_t)(2 + (cur_genotype / 2) - cur_genotype))) * (*fptr2);
		  }
		  fptr2++;
		}
		indiv_idx_stop += BITCT2;
	      }
	      if (indiv_idx == indiv_valid_ct) {
		break;
	      }
	      ulptr_end++;
	      indiv_idx_stop = indiv_valid_ct;
	    }
	  } else {
	    while (1) {
	      while (ulptr < ulptr_end) {
		cur_word = *ulptr++;
		for (; indiv_idx < indiv_idx_stop; indiv_idx++, cur_word >>= 2) {
		  cur_genotype = cur_word & 3;
		  if (cur_genotype != 1) {
		    // 0/1/2
		    *fptr++ = ((float)((intptr_t)((2 + (cur_genotype / 2) - cur_genotype) >> IS_SET(sex_male_collapsed, indiv_idx)))) * (*fptr2);
		  }
		  fptr2++;
		}
		indiv_idx_stop += BITCT2;
	      }
	      if (indiv_idx == indiv_valid_ct) {
		break;
	      }
	      ulptr_end++;
	      indiv_idx_stop = indiv_valid_ct;
	    }
	  }
	} else {
	  while (1) {
	    while (ulptr < ulptr_end) {
	      cur_word = *ulptr++;
	      for (; indiv_idx < indiv_idx_stop; indiv_idx++, cur_word >>= 2) {
		cur_genotype = cur_word & 3;
		if (cur_genotype != 1) {
		  // 0/0/1
	          *fptr++ = ((float)(cur_genotype == 0)) * (*fptr2);
		}
		fptr2++;
	      }
	      indiv_idx_stop += BITCT2;
	    }
	    if (indiv_idx == indiv_valid_ct) {
	      break;
	    }
	    ulptr_end++;
	    indiv_idx_stop = indiv_valid_ct;
	  }
	}
	fill_float_zero(fptr, align_skip);
	fptr = &(fptr[align_skip]);
      }
      if (genotypic_or_hethom) {
	param_idx++;
	if ((!is_nonx_haploid) && IS_SET(active_params, param_idx)) {
	  ulptr = loadbuf_collapsed;
	  ulptr_end = ulptr_end_init;
	  indiv_idx = 0;
	  indiv_idx_stop = BITCT2;
	  fptr2 = &(fixed_covars_cov_major[fixed_covar_idx * indiv_valid_ct]);
	  while (1) {
	    while (ulptr < ulptr_end) {
	      cur_word = *ulptr++;
	      for (; indiv_idx < indiv_idx_stop; indiv_idx++, cur_word >>= 2) {
		cur_genotype = cur_word & 3;
		if (cur_genotype != 1) {
		  // 0/1/0
                  *fptr++ = ((float)(cur_genotype == 2)) * (*fptr2);
		}
		fptr2++;
	      }
	      indiv_idx_stop += BITCT2;
	    }
	    if (indiv_idx == indiv_valid_ct) {
	      break;
	    }
	    ulptr_end++;
	    indiv_idx_stop = indiv_valid_ct;
	  }
	  fill_float_zero(fptr, align_skip);
	  fptr = &(fptr[align_skip]);
	}
      }
    }
  }
  if (include_sex) {
    if (IS_SET(active_params, sex_start_idx)) {
      copy_when_nonmissing(loadbuf_collapsed, (char*)(&(fixed_covars_cov_major[fixed_covar_nonsex_ct * indiv_valid_ct])), sizeof(float), indiv_valid_ct, missing_ct, (char*)fptr);
      fill_float_zero(&(fptr[cur_indiv_valid_ct]), align_skip);
      fptr = &(fptr[cur_indiv_valid_cta4]);
    }
    if (interactions_present) {
      if (is_set(active_params, sex_start_idx + 1)) {
	ulptr = loadbuf_collapsed;
	ulptr_end = ulptr_end_init;
	indiv_idx = 0;
	indiv_idx_stop = BITCT2;
	fptr2 = &(fixed_covars_cov_major[fixed_covar_nonsex_ct * indiv_valid_ct]);
	if (coding_flags & GLM_DOMINANT) {
	  while (1) {
	    while (ulptr < ulptr_end) {
	      cur_word = *ulptr++;
	      for (; indiv_idx < indiv_idx_stop; indiv_idx++, cur_word >>= 2) {
		cur_genotype = cur_word & 3;
		if (cur_genotype != 1) {
		  // 0/1/1
		  *fptr++ = ((float)(cur_genotype != 3)) * (*fptr2);
		}
		fptr2++;
	      }
	      indiv_idx_stop += BITCT2;
	    }
	    if (indiv_idx == indiv_valid_ct) {
	      break;
	    }
	    ulptr_end++;
	    indiv_idx_stop = indiv_valid_ct;
	  }
	} else if (!(coding_flags & (GLM_HETHOM | GLM_RECESSIVE))) {
	  if (!male_x_01) {
	    while (1) {
	      while (ulptr < ulptr_end) {
		cur_word = *ulptr++;
		for (; indiv_idx < indiv_idx_stop; indiv_idx++, cur_word >>= 2) {
		  cur_genotype = cur_word & 3;
		  if (cur_genotype != 1) {
		    // 0/1/2
		    *fptr++ = ((float)((intptr_t)(2 + (cur_genotype / 2) - cur_genotype))) * (*fptr2);
		  }
		  fptr2++;
		}
		indiv_idx_stop += BITCT2;
	      }
	      if (indiv_idx == indiv_valid_ct) {
		break;
	      }
	      ulptr_end++;
	      indiv_idx_stop = indiv_valid_ct;
	    }
	  } else {
	    while (1) {
	      while (ulptr < ulptr_end) {
		cur_word = *ulptr++;
		for (; indiv_idx < indiv_idx_stop; indiv_idx++, cur_word >>= 2) {
		  cur_genotype = cur_word & 3;
		  if (cur_genotype != 1) {
		    // 0/1/2
		    *fptr++ = ((float)((intptr_t)((2 + (cur_genotype / 2) - cur_genotype) >> IS_SET(sex_male_collapsed, indiv_idx)))) * (*fptr2);
		  }
		  fptr2++;
		}
		indiv_idx_stop += BITCT2;
	      }
	      if (indiv_idx == indiv_valid_ct) {
		break;
	      }
	      ulptr_end++;
	      indiv_idx_stop = indiv_valid_ct;
	    }
	  }
	} else {
	  while (1) {
	    while (ulptr < ulptr_end) {
	      cur_word = *ulptr++;
	      for (; indiv_idx < indiv_idx_stop; indiv_idx++, cur_word >>= 2) {
		cur_genotype = cur_word & 3;
		if (cur_genotype != 1) {
		  // 0/0/1
	          *fptr++ = ((float)(cur_genotype == 0)) * (*fptr2);
		}
		fptr2++;
	      }
	      indiv_idx_stop += BITCT2;
	    }
	    if (indiv_idx == indiv_valid_ct) {
	      break;
	    }
	    ulptr_end++;
	    indiv_idx_stop = indiv_valid_ct;
	  }
	}
	fill_float_zero(fptr, align_skip);
	fptr = &(fptr[align_skip]);
      }
      if (genotypic_or_hethom && (!is_nonx_haploid) && is_set(active_params, sex_start_idx + 2)) {
	ulptr = loadbuf_collapsed;
	ulptr_end = ulptr_end_init;
	indiv_idx = 0;
	indiv_idx_stop = BITCT2;
	fptr2 = &(fixed_covars_cov_major[fixed_covar_nonsex_ct * indiv_valid_ct]);
	while (1) {
	  while (ulptr < ulptr_end) {
	    cur_word = *ulptr++;
	    for (; indiv_idx < indiv_idx_stop; indiv_idx++, cur_word >>= 2) {
	      cur_genotype = cur_word & 3;
	      if (cur_genotype != 1) {
		// 0/1/0
		*fptr++ = ((float)(cur_genotype == 2)) * (*fptr2);
	      }
	      fptr2++;
	    }
	    indiv_idx_stop += BITCT2;
	  }
	  if (indiv_idx == indiv_valid_ct) {
	    break;
	  }
	  ulptr_end++;
	  indiv_idx_stop = indiv_valid_ct;
	}
	fill_float_zero(fptr, align_skip);
	fptr = &(fptr[align_skip]);
      }
    }
  }
#ifndef STABLE_BUILD
  if (fptr != &(cur_covars_cov_major[cur_param_ct * cur_indiv_valid_cta4])) {
    printf("assert failure:\n  cur_param_ct = %u\n  cur_indiv_valid_cta4 = %" PRIuPTR "\n  fptr - cur_covars_cov_major = %" PRIuPTR "\n", cur_param_ct, cur_indiv_valid_cta4, (uintptr_t)(fptr - cur_covars_cov_major));
    exit(1);
  }
#endif
  if (indiv_to_cluster1) {
    if (!missing_ct) {
      *cur_indiv_to_cluster1_ptr = indiv_to_cluster1;
    } else {
      copy_when_nonmissing(loadbuf_collapsed, (char*)indiv_to_cluster1, sizeof(int32_t), indiv_valid_ct, missing_ct, (char*)cur_indiv_to_cluster1_buf);
      *cur_indiv_to_cluster1_ptr = cur_indiv_to_cluster1_buf;
    }
    transpose_copy_float(cur_param_ct, cur_indiv_valid_cta4, cur_indiv_valid_ct, cur_covars_cov_major, cur_covars_indiv_major);
  }

  return missing_ct;
}

// glm main loop-specific multithread globals

#ifndef NOLAPACK
typedef struct linear_multithread_struct {
  double* indiv_1d_buf;
  double* param_2d_buf;
  double* param_2d_buf2;
  double* cluster_param_buf; // guaranteed to be size >= param_ct^2
  double* cluster_param_buf2; // not guaranteed
  MATRIX_INVERT_BUF1_TYPE* mi_buf;
  double* df_df_buf;
  double* df_buf;
  double* cur_covars_cov_major;
  double* cur_covars_indiv_major;
  uint32_t* cur_indiv_to_cluster1_buf;
  uintptr_t* perm_fails;
  double* dgels_a;
  double* dgels_b;
  double* dgels_work;
  double* param_df_buf;
  double* param_df_buf2;
  double* regression_results;
} Linear_multithread;

static Linear_multithread* g_linear_mt;
#endif

typedef struct logistic_multithread_struct {
  float* cur_covars_cov_major;
  float* cur_covars_indiv_major;
  float* coef;
  float* pp;
  float* indiv_1d_buf;
  float* pheno_buf;
  float* param_1d_buf;
  float* param_1d_buf2;
  float* param_2d_buf;
  float* param_2d_buf2;
  float* regression_results;
  uint32_t* cur_indiv_to_cluster1_buf;
  float* cluster_param_buf;
  float* cluster_param_buf2;
  double* param_1d_dbuf;
  double* param_2d_dbuf;
  double* param_2d_dbuf2;
  double* param_df_dbuf;
  double* df_df_dbuf;
  MATRIX_INVERT_BUF1_TYPE* mi_buf;
  double* df_dbuf;
  uintptr_t* perm_fails;
} Logistic_multithread;

static Logistic_multithread* g_logistic_mt;
static uintptr_t* g_joint_test_params;
static double* g_orig_stats;
#ifndef NOLAPACK
static double* g_fixed_covars_cov_major;
#endif
static uint32_t* g_indiv_to_cluster1;
static uintptr_t g_cur_param_ct;
static uintptr_t g_cur_constraint_ct;
static uint32_t g_standard_beta;
static uint32_t g_coding_flags;
static uint32_t g_glm_xchr_model;
static uintptr_t g_condition_list_start_idx;
static uintptr_t g_interaction_start_idx;
static uintptr_t g_sex_start_idx;
static uintptr_t* g_active_params;
static uintptr_t* g_haploid_params;
static uint32_t g_include_sex;
static uint32_t g_male_x_01;
static uintptr_t* g_sex_male_collapsed;
#ifndef NOLAPACK
static __CLPK_integer g_dgels_lwork;
#endif
static uint32_t g_cluster_ct1;
static double* g_constraints_con_major;
static uint32_t g_perm_batch_max;
static float* g_fixed_covars_cov_major_f;

const char glm_main_effects[] = "REC\0DOM\0HOM\0ADD";

#ifndef NOLAPACK
THREAD_RET_TYPE glm_linear_adapt_thread(void* arg) {
  intptr_t tidx = (intptr_t)arg;
  uintptr_t indiv_valid_ct = g_pheno_nm_ct;
  uintptr_t indiv_valid_ctv2 = 2 * ((indiv_valid_ct + (BITCT - 1)) / BITCT);
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  // unlike the other permutation loops, g_perms_done is not preincremented
  // here
  uint32_t pidx_offset = g_perms_done;
  uintptr_t marker_blocks = g_block_diff / CACHELINE_INT32;
  uint32_t marker_bidx = CACHELINE_INT32 * ((((uint64_t)tidx) * marker_blocks) / g_assoc_thread_ct);
  uint32_t marker_bceil = CACHELINE_INT32 * ((((uint64_t)(tidx + 1)) * marker_blocks) / g_assoc_thread_ct);
  uint32_t first_adapt_check = g_first_adapt_check;
  uintptr_t* loadbuf = g_loadbuf;
  uint32_t* adapt_m_table = &(g_adapt_m_table[marker_bidx]);
  double* perm_pmajor = g_perm_pmajor;
  unsigned char* __restrict__ perm_adapt_stop = g_perm_adapt_stop;
  uint32_t* __restrict__ perm_attempt_ct = g_perm_attempt_ct;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  uintptr_t* joint_test_params = g_joint_test_params;
  double* __restrict__ orig_stats = g_orig_stats;
  double adaptive_intercept = g_adaptive_intercept;
  double adaptive_slope = g_adaptive_slope;
  double adaptive_ci_zt = g_adaptive_ci_zt;
  double aperm_alpha = g_aperm_alpha;
  double pheno_sum_base = g_pheno_sum;
  double pheno_ssq_base = g_pheno_ssq;
  uintptr_t cur_param_ct = g_cur_param_ct;
  uintptr_t cur_constraint_ct = g_cur_constraint_ct;
  char dgels_trans = 'N';
  __CLPK_integer dgels_n = (int32_t)((uint32_t)cur_param_ct);
  __CLPK_integer dgels_lwork = g_dgels_lwork;
  uint32_t standard_beta = g_standard_beta;
  uint32_t coding_flags = g_coding_flags;
  uint32_t glm_xchr_model = g_glm_xchr_model;
  uintptr_t condition_list_start_idx = g_condition_list_start_idx;
  uintptr_t interaction_start_idx = g_interaction_start_idx;
  uintptr_t sex_start_idx = g_sex_start_idx;
  uintptr_t* active_params = g_active_params;
  uintptr_t* haploid_params = g_haploid_params;
  uint32_t include_sex = g_include_sex;
  uint32_t male_x_01 = g_male_x_01;
  uint32_t cluster_ct1 = g_cluster_ct1;
  uintptr_t* sex_male_collapsed = g_sex_male_collapsed;
  uint32_t is_nonx_haploid = g_is_haploid && (!g_is_x);
  double* fixed_covars_cov_major = g_fixed_covars_cov_major;
  uint32_t* indiv_to_cluster1 = g_indiv_to_cluster1;
  double* constraints_con_major = g_constraints_con_major;
  double* indiv_1d_buf = g_linear_mt[tidx].indiv_1d_buf;
  double* param_2d_buf = g_linear_mt[tidx].param_2d_buf;
  double* param_2d_buf2 = g_linear_mt[tidx].param_2d_buf2;
  double* cluster_param_buf = g_linear_mt[tidx].cluster_param_buf;
  double* cluster_param_buf2 = g_linear_mt[tidx].cluster_param_buf2;
  MATRIX_INVERT_BUF1_TYPE* mi_buf = g_linear_mt[tidx].mi_buf;
  double* df_df_buf = g_linear_mt[tidx].df_df_buf;
  double* df_buf = g_linear_mt[tidx].df_buf;
  double* cur_covars_cov_major = g_linear_mt[tidx].cur_covars_cov_major;
  double* cur_covars_indiv_major = g_linear_mt[tidx].cur_covars_indiv_major;
  uint32_t* cur_indiv_to_cluster1_buf = g_linear_mt[tidx].cur_indiv_to_cluster1_buf;
  uintptr_t* perm_fails = g_linear_mt[tidx].perm_fails;
  double* dgels_a = g_linear_mt[tidx].dgels_a;
  double* dgels_b = g_linear_mt[tidx].dgels_b;
  double* dgels_work = g_linear_mt[tidx].dgels_work;
  double* param_df_buf = g_linear_mt[tidx].param_df_buf;
  double* param_df_buf2 = g_linear_mt[tidx].param_df_buf2;
  double* regression_results = g_linear_mt[tidx].regression_results;
  double* dptr;
  uintptr_t* loadbuf_ptr;
  uint32_t* cur_indiv_to_cluster1;
  uintptr_t cur_missing_ct;
  uintptr_t cur_indiv_valid_ct;
  uintptr_t pidx;
  uintptr_t indiv_idx;
  uintptr_t param_ctx_m1;
  uint32_t marker_idx;
  uint32_t next_adapt_check;
  uint32_t success_2start;
  uint32_t success_2incr;
  uint32_t attempts;
  uint32_t cur_attempts;
  uint32_t perm_fail_ct;
  uint32_t cur_fail_ct;
  __CLPK_integer dgels_m;
  __CLPK_integer dgels_nrhs;
  __CLPK_integer dgels_ldb;
  __CLPK_integer dgels_info;
  double stat_high;
  double stat_low;
  double pval;
  double dxx;
  double dyy;
  double dzz;
  if (cur_constraint_ct) {
    param_ctx_m1 = cur_param_ct;
  } else {
    param_ctx_m1 = cur_param_ct - 1;
  }
  if ((uintptr_t)tidx + 1 == g_assoc_thread_ct) {
    marker_bceil = g_block_diff;
  }
  for (; marker_bidx < marker_bceil; marker_bidx++) {
    marker_idx = *adapt_m_table++;
    if (perm_adapt_stop[marker_idx]) {
      continue;
    }
    success_2start = perm_2success_ct[marker_idx];
    attempts = perm_attempt_ct[marker_idx];
    next_adapt_check = first_adapt_check;
    stat_high = orig_stats[marker_idx] + EPSILON;
    stat_low = orig_stats[marker_idx] - EPSILON;
    loadbuf_ptr = &(loadbuf[marker_bidx * indiv_valid_ctv2]);
    cur_missing_ct = glm_fill_design(loadbuf_ptr, fixed_covars_cov_major, indiv_valid_ct, indiv_to_cluster1, cur_param_ct, coding_flags, glm_xchr_model, condition_list_start_idx, interaction_start_idx, sex_start_idx, active_params, haploid_params, include_sex, male_x_01, sex_male_collapsed, is_nonx_haploid, cur_covars_cov_major, cur_covars_indiv_major, cur_indiv_to_cluster1_buf, &cur_indiv_to_cluster1, standard_beta);
    cur_indiv_valid_ct = indiv_valid_ct - cur_missing_ct;
    dgels_m = (int32_t)((uint32_t)cur_indiv_valid_ct);
    dgels_ldb = dgels_m;
    dptr = dgels_b;
    success_2incr = 0;
    cur_fail_ct = 0;
    memcpy(dgels_a, cur_covars_cov_major, cur_param_ct * cur_indiv_valid_ct * sizeof(double));
    if (!cur_missing_ct) {
      memcpy(dgels_b, perm_pmajor, perm_vec_ct * indiv_valid_ct * sizeof(double));
    } else {
      for (pidx = 0; pidx < perm_vec_ct; pidx++) {
	copy_when_nonmissing(loadbuf_ptr, (char*)(&(perm_pmajor[pidx * indiv_valid_ct])), sizeof(double), indiv_valid_ct, cur_missing_ct, (char*)(&(dgels_b[pidx * cur_indiv_valid_ct])));
	if (standard_beta) {
	  dxx = 0;
	  dyy = 0;
	  dptr = &(dgels_b[pidx * cur_indiv_valid_ct]);
	  for (indiv_idx = 0; indiv_idx < cur_indiv_valid_ct; indiv_idx++) {
	    dzz = *dptr++;
	    dxx += dzz;
	    dyy += dzz * dzz;
	  }
	  dptr = &(dgels_b[pidx * cur_indiv_valid_ct]);
	  dzz = dxx / ((double)((intptr_t)cur_indiv_valid_ct));
	  dyy = sqrt(((double)((intptr_t)(cur_indiv_valid_ct - 1))) / (dyy - dxx * dzz));
	  for (indiv_idx = 0; indiv_idx < cur_indiv_valid_ct; indiv_idx++) {
	    *dptr = ((*dptr) - dzz) * dyy;
	    dptr++;
	  }
	}
      }
    }
    dgels_nrhs = (int32_t)((uint32_t)perm_vec_ct);
    dgels_(&dgels_trans, &dgels_m, &dgels_n, &dgels_nrhs, dgels_a, &dgels_m, dgels_b, &dgels_ldb, dgels_work, &dgels_lwork, &dgels_info);
    glm_linear_robust_cluster_covar(perm_vec_ct, cur_param_ct, cur_indiv_valid_ct, cur_missing_ct, loadbuf_ptr, standard_beta, pheno_sum_base, pheno_ssq_base, cur_covars_cov_major, cur_covars_indiv_major, perm_pmajor, dgels_b, param_2d_buf, mi_buf, param_2d_buf2, cluster_ct1, cur_indiv_to_cluster1, cluster_param_buf, cluster_param_buf2, indiv_1d_buf, regression_results, cur_constraint_ct, constraints_con_major, param_df_buf, param_df_buf2, df_df_buf, df_buf, &perm_fail_ct, perm_fails);
    for (pidx = 0; pidx < perm_vec_ct;) {
      if (!IS_SET(perm_fails, pidx)) {
	if (!joint_test_params) {
	  dxx = dgels_b[pidx * cur_indiv_valid_ct + 1];
	  dyy = sqrt(regression_results[pidx * param_ctx_m1]);
	  dxx = fabs(dxx / dyy);
	  if (dxx > stat_high) {
	    success_2incr += 2;
	  } else if (dxx > stat_low) {
	    success_2incr++;
	  }
	} else {
	  dxx = regression_results[(pidx + 1) * param_ctx_m1 - 1];
	  if (dxx > stat_high) {
	    success_2incr += 2;
	  } else if (dxx > stat_low) {
	    success_2incr++;
	  } else if (dxx == -9) {
	    cur_fail_ct++;
	  }
	}
      } else {
	cur_fail_ct++;
      }
      if (++pidx == next_adapt_check - pidx_offset) {
	if (success_2start + success_2incr) {
	  cur_attempts = attempts + pidx - cur_fail_ct;
	  pval = ((double)((int32_t)(success_2start + success_2incr + 2))) / ((double)(2 * ((int32_t)(cur_attempts + 1))));
	  dxx = adaptive_ci_zt * sqrt(pval * (1 - pval) / ((int32_t)cur_attempts));
	  dyy = pval - dxx; // lower bound
	  dzz = pval + dxx; // upper bound
	  if ((dyy > aperm_alpha) || (dzz < aperm_alpha)) {
	    perm_adapt_stop[marker_idx] = 1;
	    goto glm_linear_adapt_thread_done;
	  }
	}
	next_adapt_check += (int32_t)(adaptive_intercept + ((int32_t)next_adapt_check) * adaptive_slope);
      }
    }
    cur_attempts = attempts + perm_vec_ct - cur_fail_ct;
  glm_linear_adapt_thread_done:
    perm_2success_ct[marker_idx] = success_2start + success_2incr;
    perm_attempt_ct[marker_idx] = cur_attempts;
  }
  THREAD_RETURN;
}
#endif

THREAD_RET_TYPE glm_logistic_adapt_thread(void* arg) {
  intptr_t tidx = (intptr_t)arg;
  uintptr_t indiv_valid_ct = g_pheno_nm_ct;
  uintptr_t indiv_valid_ctv2 = 2 * ((indiv_valid_ct + (BITCT - 1)) / BITCT);
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uint32_t pidx_offset = g_perms_done;
  uintptr_t marker_blocks = g_block_diff / CACHELINE_INT32;
  uint32_t marker_bidx = CACHELINE_INT32 * ((((uint64_t)tidx) * marker_blocks) / g_assoc_thread_ct);
  uint32_t marker_bceil = CACHELINE_INT32 * ((((uint64_t)(tidx + 1)) * marker_blocks) / g_assoc_thread_ct);
  uint32_t first_adapt_check = g_first_adapt_check;
  uintptr_t* loadbuf = g_loadbuf;
  uint32_t* adapt_m_table = &(g_adapt_m_table[marker_bidx]);
  uintptr_t* perm_vecs = g_perm_vecs;
  unsigned char* __restrict__ perm_adapt_stop = g_perm_adapt_stop;
  uint32_t* __restrict__ perm_attempt_ct = g_perm_attempt_ct;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  uintptr_t* joint_test_params = g_joint_test_params;
  double* __restrict__ orig_stats = g_orig_stats;
  double adaptive_intercept = g_adaptive_intercept;
  double adaptive_slope = g_adaptive_slope;
  double adaptive_ci_zt = g_adaptive_ci_zt;
  double aperm_alpha = g_aperm_alpha;
  uintptr_t cur_param_ct = g_cur_param_ct;
  uintptr_t cur_constraint_ct = g_cur_constraint_ct;
  uint32_t coding_flags = g_coding_flags;
  uint32_t glm_xchr_model = g_glm_xchr_model;
  uintptr_t condition_list_start_idx = g_condition_list_start_idx;
  uintptr_t interaction_start_idx = g_interaction_start_idx;
  uintptr_t sex_start_idx = g_sex_start_idx;
  uintptr_t* active_params = g_active_params;
  uintptr_t* haploid_params = g_haploid_params;
  uint32_t include_sex = g_include_sex;
  uint32_t male_x_01 = g_male_x_01;
  uint32_t cluster_ct1 = g_cluster_ct1;
  uintptr_t* sex_male_collapsed = g_sex_male_collapsed;
  uint32_t is_nonx_haploid = g_is_haploid && (!g_is_x);
  float* fixed_covars_cov_major = g_fixed_covars_cov_major_f;
  uint32_t* indiv_to_cluster1 = g_indiv_to_cluster1;
  double* constraints_con_major = g_constraints_con_major;
  float* cur_covars_cov_major = g_logistic_mt[tidx].cur_covars_cov_major;
  float* cur_covars_indiv_major = g_logistic_mt[tidx].cur_covars_indiv_major;
  float* coef = g_logistic_mt[tidx].coef;
  float* pp = g_logistic_mt[tidx].pp;
  float* indiv_1d_buf = g_logistic_mt[tidx].indiv_1d_buf;
  float* pheno_buf = g_logistic_mt[tidx].pheno_buf;
  float* param_1d_buf = g_logistic_mt[tidx].param_1d_buf;
  float* param_1d_buf2 = g_logistic_mt[tidx].param_1d_buf2;
  float* param_2d_buf = g_logistic_mt[tidx].param_2d_buf;
  float* param_2d_buf2 = g_logistic_mt[tidx].param_2d_buf2;
  float* regression_results = g_logistic_mt[tidx].regression_results;
  uint32_t* cur_indiv_to_cluster1_buf = g_logistic_mt[tidx].cur_indiv_to_cluster1_buf;
  float* cluster_param_buf = g_logistic_mt[tidx].cluster_param_buf;
  float* cluster_param_buf2 = g_logistic_mt[tidx].cluster_param_buf2;
  double* param_1d_dbuf = g_logistic_mt[tidx].param_1d_dbuf;
  double* param_2d_dbuf = g_logistic_mt[tidx].param_2d_dbuf;
  double* param_2d_dbuf2 = g_logistic_mt[tidx].param_2d_dbuf2;
  double* param_df_dbuf = g_logistic_mt[tidx].param_df_dbuf;
  MATRIX_INVERT_BUF1_TYPE* mi_buf = g_logistic_mt[tidx].mi_buf;
  double* df_df_dbuf = g_logistic_mt[tidx].df_df_dbuf;
  double* df_dbuf = g_logistic_mt[tidx].df_dbuf;
  uintptr_t* perm_fails = g_logistic_mt[tidx].perm_fails;
  uintptr_t* loadbuf_ptr;
  uint32_t* cur_indiv_to_cluster1;
  uintptr_t cur_missing_ct;
  uintptr_t cur_indiv_valid_ct;
  uintptr_t pidx;
  uintptr_t param_ctx_m1;
  uint32_t marker_idx;
  uint32_t next_adapt_check;
  uint32_t success_2start;
  uint32_t success_2incr;
  uint32_t attempts;
  uint32_t cur_attempts;
  uint32_t cur_fail_ct;
  double stat_high;
  double stat_low;
  double pval;
  double dxx;
  double dyy;
  double dzz;
  if (cur_constraint_ct) {
    param_ctx_m1 = cur_param_ct;
  } else {
    param_ctx_m1 = cur_param_ct - 1;
  }
  if ((uintptr_t)tidx + 1 == g_assoc_thread_ct) {
    marker_bceil = g_block_diff;
  }
  for (; marker_bidx < marker_bceil; marker_bidx++) {
    marker_idx = *adapt_m_table++;
    if (perm_adapt_stop[marker_idx]) {
      continue;
    }
    success_2start = perm_2success_ct[marker_idx];
    attempts = perm_attempt_ct[marker_idx];
    next_adapt_check = first_adapt_check;
    stat_high = orig_stats[marker_idx] + EPSILON;
    stat_low = orig_stats[marker_idx] - EPSILON;
    loadbuf_ptr = &(loadbuf[marker_bidx * indiv_valid_ctv2]);
    cur_missing_ct = glm_fill_design_float(loadbuf_ptr, fixed_covars_cov_major, indiv_valid_ct, indiv_to_cluster1, cur_param_ct, coding_flags, glm_xchr_model, condition_list_start_idx, interaction_start_idx, sex_start_idx, active_params, haploid_params, include_sex, male_x_01, sex_male_collapsed, is_nonx_haploid, cur_covars_cov_major, cur_covars_indiv_major, cur_indiv_to_cluster1_buf, &cur_indiv_to_cluster1);
    cur_indiv_valid_ct = indiv_valid_ct - cur_missing_ct;
    success_2incr = 0;
    cur_fail_ct = 0;
    // todo: try better starting position
    fill_float_zero(coef, ((cur_param_ct + 3) & (~3)) * perm_vec_ct);
    glm_logistic_robust_cluster_covar(perm_vec_ct, cur_param_ct, cur_indiv_valid_ct, cur_missing_ct, loadbuf_ptr, cur_covars_cov_major, cur_covars_indiv_major, perm_vecs, coef, pp, indiv_1d_buf, pheno_buf, param_1d_buf, param_1d_buf2, param_2d_buf, param_2d_buf2, regression_results, cluster_ct1, cur_indiv_to_cluster1, cluster_param_buf, cluster_param_buf2, cur_constraint_ct, constraints_con_major, param_1d_dbuf, param_2d_dbuf, param_2d_dbuf2, param_df_dbuf, df_df_dbuf, mi_buf, df_dbuf, perm_fails);
    for (pidx = 0; pidx < perm_vec_ct;) {
      if (!IS_SET(perm_fails, pidx)) {
	if (!joint_test_params) {
	  dxx = (double)coef[pidx * cur_param_ct + 1];
	  dyy = sqrt((double)regression_results[pidx * param_ctx_m1]);
	  dxx /= dyy;
	  dxx *= dxx;
	  if (dxx > stat_high) {
	    success_2incr += 2;
	  } else if (dxx > stat_low) {
	    success_2incr++;
	  }
	} else {
	  dxx = (double)regression_results[(pidx + 1) * param_ctx_m1 - 1];
	  if (dxx > stat_high) {
	    success_2incr += 2;
	  } else if (dxx > stat_low) {
	    success_2incr++;
	  } else if (dxx == -9) {
	    cur_fail_ct++;
	  }
	}
      } else {
	cur_fail_ct++;
      }
      if (++pidx == next_adapt_check - pidx_offset) {
	if (success_2start + success_2incr) {
	  cur_attempts = attempts + pidx - cur_fail_ct;
	  pval = ((double)((int32_t)(success_2start + success_2incr + 2))) / ((double)(2 * ((int32_t)(cur_attempts + 1))));
	  dxx = adaptive_ci_zt * sqrt(pval * (1 - pval) / ((int32_t)cur_attempts));
	  dyy = pval - dxx; // lower bound
	  dzz = pval + dxx; // upper bound
	  if ((dyy > aperm_alpha) || (dzz < aperm_alpha)) {
	    perm_adapt_stop[marker_idx] = 1;
	    goto glm_logistic_adapt_thread_done;
	  }
	}
	next_adapt_check += (int32_t)(adaptive_intercept + ((int32_t)next_adapt_check) * adaptive_slope);
      }
    }
    cur_attempts = attempts + perm_vec_ct - cur_fail_ct;
  glm_logistic_adapt_thread_done:
    perm_2success_ct[marker_idx] = success_2start + success_2incr;
    perm_attempt_ct[marker_idx] = cur_attempts;
  }
  THREAD_RETURN;
}

#ifndef NOLAPACK
THREAD_RET_TYPE glm_linear_maxt_thread(void* arg) {
  intptr_t tidx = (intptr_t)arg;
  uintptr_t indiv_valid_ct = g_pheno_nm_ct;
  uintptr_t indiv_valid_ctv2 = 2 * ((indiv_valid_ct + (BITCT - 1)) / BITCT);
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uint32_t pidx_offset = g_perms_done;
  uintptr_t marker_blocks = g_block_diff / CACHELINE_INT32;
  uint32_t marker_bidx = CACHELINE_INT32 * ((((uint64_t)tidx) * marker_blocks) / g_assoc_thread_ct);
  uint32_t marker_bceil = CACHELINE_INT32 * ((((uint64_t)(tidx + 1)) * marker_blocks) / g_assoc_thread_ct);
  uintptr_t* loadbuf = g_loadbuf;
  uint32_t* adapt_m_table = &(g_adapt_m_table[marker_bidx]);
  double* perm_pmajor = g_perm_pmajor;
  uintptr_t perm_vec_ctcl8m = (perm_vec_ct + (CACHELINE_DBL - 1)) & (~(CACHELINE_DBL - 1));
  double* __restrict__ results = &(g_maxt_thread_results[perm_vec_ctcl8m * tidx]);
  unsigned char* __restrict__ perm_adapt_stop = g_perm_adapt_stop;
  uint32_t* __restrict__ perm_fail_cts = g_perm_attempt_ct;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  uintptr_t* joint_test_params = g_joint_test_params;
  double* __restrict__ orig_stats = g_orig_stats;
  double pheno_sum_base = g_pheno_sum;
  double pheno_ssq_base = g_pheno_ssq;
  uintptr_t cur_param_ct = g_cur_param_ct;
  uintptr_t cur_constraint_ct = g_cur_constraint_ct;
  char dgels_trans = 'N';
  __CLPK_integer dgels_n = (int32_t)((uint32_t)cur_param_ct);
  __CLPK_integer dgels_lwork = g_dgels_lwork;
  uint32_t standard_beta = g_standard_beta;
  uint32_t coding_flags = g_coding_flags;
  uint32_t glm_xchr_model = g_glm_xchr_model;
  uintptr_t condition_list_start_idx = g_condition_list_start_idx;
  uintptr_t interaction_start_idx = g_interaction_start_idx;
  uintptr_t sex_start_idx = g_sex_start_idx;
  uintptr_t* active_params = g_active_params;
  uintptr_t* haploid_params = g_haploid_params;
  uint32_t include_sex = g_include_sex;
  uint32_t male_x_01 = g_male_x_01;
  uint32_t cluster_ct1 = g_cluster_ct1;
  uintptr_t perm_batch_max = g_perm_batch_max;
  uintptr_t* sex_male_collapsed = g_sex_male_collapsed;
  uint32_t is_nonx_haploid = g_is_haploid && (!g_is_x);
  double* fixed_covars_cov_major = g_fixed_covars_cov_major;
  uint32_t* indiv_to_cluster1 = g_indiv_to_cluster1;
  double* constraints_con_major = g_constraints_con_major;
  double* msa_ptr = NULL;
  double* indiv_1d_buf = g_linear_mt[tidx].indiv_1d_buf;
  double* param_2d_buf = g_linear_mt[tidx].param_2d_buf;
  double* param_2d_buf2 = g_linear_mt[tidx].param_2d_buf2;
  double* cluster_param_buf = g_linear_mt[tidx].cluster_param_buf;
  double* cluster_param_buf2 = g_linear_mt[tidx].cluster_param_buf2;
  MATRIX_INVERT_BUF1_TYPE* mi_buf = g_linear_mt[tidx].mi_buf;
  double* df_df_buf = g_linear_mt[tidx].df_df_buf;
  double* df_buf = g_linear_mt[tidx].df_buf;
  double* cur_covars_cov_major = g_linear_mt[tidx].cur_covars_cov_major;
  double* cur_covars_indiv_major = g_linear_mt[tidx].cur_covars_indiv_major;
  uint32_t* cur_indiv_to_cluster1_buf = g_linear_mt[tidx].cur_indiv_to_cluster1_buf;
  uintptr_t* perm_fails = g_linear_mt[tidx].perm_fails;
  double* dgels_a = g_linear_mt[tidx].dgels_a;
  double* dgels_b = g_linear_mt[tidx].dgels_b;
  double* dgels_work = g_linear_mt[tidx].dgels_work;
  double* param_df_buf = g_linear_mt[tidx].param_df_buf;
  double* param_df_buf2 = g_linear_mt[tidx].param_df_buf2;
  double* regression_results = g_linear_mt[tidx].regression_results;
  double* dptr;
  uintptr_t* loadbuf_ptr;
  uint32_t* cur_indiv_to_cluster1;
  uintptr_t cur_missing_ct;
  uintptr_t cur_indiv_valid_ct;
  uintptr_t pidx;
  uintptr_t indiv_idx;
  uintptr_t param_ctx_m1;
  uint32_t marker_idx;
  uint32_t success_2incr;
  uint32_t perm_fail_ct;
  __CLPK_integer dgels_m;
  __CLPK_integer dgels_nrhs;
  __CLPK_integer dgels_ldb;
  __CLPK_integer dgels_info;
  double stat_high;
  double stat_low;
  double dxx;
  double dyy;
  double dzz;
  if (cur_constraint_ct) {
    param_ctx_m1 = cur_param_ct;
  } else {
    param_ctx_m1 = cur_param_ct - 1;
  }
  if ((uintptr_t)tidx + 1 == g_assoc_thread_ct) {
    marker_bceil = g_block_diff;
  }
  memcpy(results, &(g_maxt_extreme_stat[pidx_offset]), perm_vec_ct * sizeof(double));
  for (; marker_bidx < marker_bceil; marker_bidx++) {
    marker_idx = *adapt_m_table++;
    if (perm_adapt_stop[marker_idx]) {
      if (g_mperm_save_all && (perm_batch_max != perm_vec_ct)) {
	msa_ptr = &(g_mperm_save_all[marker_idx * perm_vec_ct]);
	for (pidx = 0; pidx < perm_vec_ct; pidx++) {
	  *msa_ptr++ = -9;
	}
      }
      continue;
    }
    if (g_mperm_save_all) {
      msa_ptr = &(g_mperm_save_all[marker_idx * perm_vec_ct]);
    }
    stat_high = orig_stats[marker_idx] + EPSILON;
    stat_low = orig_stats[marker_idx] - EPSILON;
    loadbuf_ptr = &(loadbuf[marker_bidx * indiv_valid_ctv2]);
    cur_missing_ct = glm_fill_design(loadbuf_ptr, fixed_covars_cov_major, indiv_valid_ct, indiv_to_cluster1, cur_param_ct, coding_flags, glm_xchr_model, condition_list_start_idx, interaction_start_idx, sex_start_idx, active_params, haploid_params, include_sex, male_x_01, sex_male_collapsed, is_nonx_haploid, cur_covars_cov_major, cur_covars_indiv_major, cur_indiv_to_cluster1_buf, &cur_indiv_to_cluster1, standard_beta);
    cur_indiv_valid_ct = indiv_valid_ct - cur_missing_ct;
    dgels_m = (int32_t)((uint32_t)cur_indiv_valid_ct);
    dgels_ldb = dgels_m;
    dptr = dgels_b;
    success_2incr = 0;
    memcpy(dgels_a, cur_covars_cov_major, cur_param_ct * cur_indiv_valid_ct * sizeof(double));
    if (!cur_missing_ct) {
      memcpy(dgels_b, perm_pmajor, perm_vec_ct * indiv_valid_ct * sizeof(double));
    } else {
      for (pidx = 0; pidx < perm_vec_ct; pidx++) {
	copy_when_nonmissing(loadbuf_ptr, (char*)(&(perm_pmajor[pidx * indiv_valid_ct])), sizeof(double), indiv_valid_ct, cur_missing_ct, (char*)(&(dgels_b[pidx * cur_indiv_valid_ct])));
	if (standard_beta) {
	  dxx = 0;
	  dyy = 0;
	  dptr = &(dgels_b[pidx * cur_indiv_valid_ct]);
	  for (indiv_idx = 0; indiv_idx < cur_indiv_valid_ct; indiv_idx++) {
	    dzz = *dptr++;
	    dxx += dzz;
	    dyy += dzz * dzz;
	  }
	  dptr = &(dgels_b[pidx * cur_indiv_valid_ct]);
	  dzz = dxx / ((double)((intptr_t)cur_indiv_valid_ct));
	  dyy = sqrt(((double)((intptr_t)(cur_indiv_valid_ct - 1))) / (dyy - dxx * dzz));
	  for (indiv_idx = 0; indiv_idx < cur_indiv_valid_ct; indiv_idx++) {
	    *dptr = ((*dptr) - dzz) * dyy;
	    dptr++;
	  }
	}
      }
    }
    dgels_nrhs = (int32_t)((uint32_t)perm_vec_ct);
    dgels_(&dgels_trans, &dgels_m, &dgels_n, &dgels_nrhs, dgels_a, &dgels_m, dgels_b, &dgels_ldb, dgels_work, &dgels_lwork, &dgels_info);
    glm_linear_robust_cluster_covar(perm_vec_ct, cur_param_ct, cur_indiv_valid_ct, cur_missing_ct, loadbuf_ptr, standard_beta, pheno_sum_base, pheno_ssq_base, cur_covars_cov_major, cur_covars_indiv_major, perm_pmajor, dgels_b, param_2d_buf, mi_buf, param_2d_buf2, cluster_ct1, cur_indiv_to_cluster1, cluster_param_buf, cluster_param_buf2, indiv_1d_buf, regression_results, cur_constraint_ct, constraints_con_major, param_df_buf, param_df_buf2, df_df_buf, df_buf, &perm_fail_ct, perm_fails);
    for (pidx = 0; pidx < perm_vec_ct; pidx++) {
      if (!IS_SET(perm_fails, pidx)) {
	if (!joint_test_params) {
	  dxx = dgels_b[pidx * cur_indiv_valid_ct + 1];
	  dyy = sqrt(regression_results[pidx * param_ctx_m1]);
	  dxx = fabs(dxx / dyy);
	  if (dxx > stat_high) {
	    success_2incr += 2;
	  } else if (dxx > stat_low) {
	    success_2incr++;
	  }
	} else {
	  dxx = regression_results[(pidx + 1) * param_ctx_m1 - 1];
	  if (dxx > stat_high) {
	    success_2incr += 2;
	  } else if (dxx > stat_low) {
	    success_2incr++;
	  } else if (dxx == -9) {
	    perm_fail_ct++;
	  }
	}
	if (results[pidx] < dxx) {
	  results[pidx] = dxx;
	}
	if (msa_ptr) {
	  *msa_ptr++ = dxx;
	}
      } else if (msa_ptr) {
	*msa_ptr++ = -9;
      }
    }
    perm_2success_ct[marker_idx] += success_2incr;
    if (perm_fail_ct) {
      perm_fail_cts[marker_idx] += perm_fail_ct;
    }
  }
  THREAD_RETURN;
}
#endif

THREAD_RET_TYPE glm_logistic_maxt_thread(void* arg) {
  intptr_t tidx = (intptr_t)arg;
  uintptr_t indiv_valid_ct = g_pheno_nm_ct;
  uintptr_t indiv_valid_ctv2 = 2 * ((indiv_valid_ct + (BITCT - 1)) / BITCT);
  uintptr_t perm_vec_ct = g_perm_vec_ct;
  uint32_t pidx_offset = g_perms_done;
  uintptr_t marker_blocks = g_block_diff / CACHELINE_INT32;
  uint32_t marker_bidx = CACHELINE_INT32 * ((((uint64_t)tidx) * marker_blocks) / g_assoc_thread_ct);
  uint32_t marker_bceil = CACHELINE_INT32 * ((((uint64_t)(tidx + 1)) * marker_blocks) / g_assoc_thread_ct);
  uintptr_t* loadbuf = g_loadbuf;
  uint32_t* adapt_m_table = &(g_adapt_m_table[marker_bidx]);
  uintptr_t* perm_vecs = g_perm_vecs;
  uintptr_t perm_vec_ctcl8m = (perm_vec_ct + (CACHELINE_DBL - 1)) & (~(CACHELINE_DBL - 1));
  double* __restrict__ results = &(g_maxt_thread_results[perm_vec_ctcl8m * tidx]);
  unsigned char* __restrict__ perm_adapt_stop = g_perm_adapt_stop;
  uint32_t* __restrict__ perm_fail_cts = g_perm_attempt_ct;
  uint32_t* __restrict__ perm_2success_ct = g_perm_2success_ct;
  uintptr_t* joint_test_params = g_joint_test_params;
  double* __restrict__ orig_stats = g_orig_stats;
  uintptr_t cur_param_ct = g_cur_param_ct;
  uintptr_t cur_param_cta4 = (cur_param_ct + 3) & (~3);
  uintptr_t cur_constraint_ct = g_cur_constraint_ct;
  uint32_t coding_flags = g_coding_flags;
  uint32_t glm_xchr_model = g_glm_xchr_model;
  uintptr_t condition_list_start_idx = g_condition_list_start_idx;
  uintptr_t interaction_start_idx = g_interaction_start_idx;
  uintptr_t sex_start_idx = g_sex_start_idx;
  uintptr_t* active_params = g_active_params;
  uintptr_t* haploid_params = g_haploid_params;
  uint32_t include_sex = g_include_sex;
  uint32_t male_x_01 = g_male_x_01;
  uint32_t cluster_ct1 = g_cluster_ct1;
  uintptr_t perm_batch_max = g_perm_batch_max;
  uintptr_t* sex_male_collapsed = g_sex_male_collapsed;
  uint32_t is_nonx_haploid = g_is_haploid && (!g_is_x);
  float* fixed_covars_cov_major = g_fixed_covars_cov_major_f;
  uint32_t* indiv_to_cluster1 = g_indiv_to_cluster1;
  double* constraints_con_major = g_constraints_con_major;
  double* msa_ptr = NULL;
  float* cur_covars_cov_major = g_logistic_mt[tidx].cur_covars_cov_major;
  float* cur_covars_indiv_major = g_logistic_mt[tidx].cur_covars_indiv_major;
  float* coef = g_logistic_mt[tidx].coef;
  float* pp = g_logistic_mt[tidx].pp;
  float* indiv_1d_buf = g_logistic_mt[tidx].indiv_1d_buf;
  float* pheno_buf = g_logistic_mt[tidx].pheno_buf;
  float* param_1d_buf = g_logistic_mt[tidx].param_1d_buf;
  float* param_1d_buf2 = g_logistic_mt[tidx].param_1d_buf2;
  float* param_2d_buf = g_logistic_mt[tidx].param_2d_buf;
  float* param_2d_buf2 = g_logistic_mt[tidx].param_2d_buf2;
  float* regression_results = g_logistic_mt[tidx].regression_results;
  uint32_t* cur_indiv_to_cluster1_buf = g_logistic_mt[tidx].cur_indiv_to_cluster1_buf;
  float* cluster_param_buf = g_logistic_mt[tidx].cluster_param_buf;
  float* cluster_param_buf2 = g_logistic_mt[tidx].cluster_param_buf2;
  double* param_1d_dbuf = g_logistic_mt[tidx].param_1d_dbuf;
  double* param_2d_dbuf = g_logistic_mt[tidx].param_2d_dbuf;
  double* param_2d_dbuf2 = g_logistic_mt[tidx].param_2d_dbuf2;
  double* param_df_dbuf = g_logistic_mt[tidx].param_df_dbuf;
  MATRIX_INVERT_BUF1_TYPE* mi_buf = g_logistic_mt[tidx].mi_buf;
  double* df_df_dbuf = g_logistic_mt[tidx].df_df_dbuf;
  double* df_dbuf = g_logistic_mt[tidx].df_dbuf;
  uintptr_t* perm_fails = g_logistic_mt[tidx].perm_fails;
  uintptr_t* loadbuf_ptr;
  uint32_t* cur_indiv_to_cluster1;
  uintptr_t cur_missing_ct;
  uintptr_t cur_indiv_valid_ct;
  uintptr_t pidx;
  uintptr_t param_ctx_m1;
  uint32_t marker_idx;
  uint32_t success_2incr;
  uint32_t perm_fail_ct;
  double stat_high;
  double stat_low;
  double dxx;
  double dyy;
  if (cur_constraint_ct) {
    param_ctx_m1 = cur_param_ct;
  } else {
    param_ctx_m1 = cur_param_ct - 1;
  }
  if ((uintptr_t)tidx + 1 == g_assoc_thread_ct) {
    marker_bceil = g_block_diff;
  }
  memcpy(results, &(g_maxt_extreme_stat[pidx_offset]), perm_vec_ct * sizeof(double));
  for (; marker_bidx < marker_bceil; marker_bidx++) {
    marker_idx = *adapt_m_table++;
    if (perm_adapt_stop[marker_idx]) {
      if (g_mperm_save_all && (perm_batch_max != perm_vec_ct)) {
	msa_ptr = &(g_mperm_save_all[marker_idx * perm_vec_ct]);
	for (pidx = 0; pidx < perm_vec_ct; pidx++) {
	  *msa_ptr++ = -9;
	}
      }
      continue;
    }
    if (g_mperm_save_all) {
      msa_ptr = &(g_mperm_save_all[marker_idx * perm_vec_ct]);
    }
    stat_high = orig_stats[marker_idx] + EPSILON;
    stat_low = orig_stats[marker_idx] - EPSILON;
    loadbuf_ptr = &(loadbuf[marker_bidx * indiv_valid_ctv2]);
    cur_missing_ct = glm_fill_design_float(loadbuf_ptr, fixed_covars_cov_major, indiv_valid_ct, indiv_to_cluster1, cur_param_ct, coding_flags, glm_xchr_model, condition_list_start_idx, interaction_start_idx, sex_start_idx, active_params, haploid_params, include_sex, male_x_01, sex_male_collapsed, is_nonx_haploid, cur_covars_cov_major, cur_covars_indiv_major, cur_indiv_to_cluster1_buf, &cur_indiv_to_cluster1);
    cur_indiv_valid_ct = indiv_valid_ct - cur_missing_ct;
    success_2incr = 0;
    // todo: try better starting position
    fill_float_zero(coef, cur_param_cta4 * perm_vec_ct);
    perm_fail_ct = glm_logistic_robust_cluster_covar(perm_vec_ct, cur_param_ct, cur_indiv_valid_ct, cur_missing_ct, loadbuf_ptr, cur_covars_cov_major, cur_covars_indiv_major, perm_vecs, coef, pp, indiv_1d_buf, pheno_buf, param_1d_buf, param_1d_buf2, param_2d_buf, param_2d_buf2, regression_results, cluster_ct1, cur_indiv_to_cluster1, cluster_param_buf, cluster_param_buf2, cur_constraint_ct, constraints_con_major, param_1d_dbuf, param_2d_dbuf, param_2d_dbuf2, param_df_dbuf, df_df_dbuf, mi_buf, df_dbuf, perm_fails);
    for (pidx = 0; pidx < perm_vec_ct; pidx++) {
      if (!IS_SET(perm_fails, pidx)) {
	if (!joint_test_params) {
	  dxx = (double)coef[pidx * cur_param_cta4 + 1];
	  dyy = sqrt((double)regression_results[pidx * param_ctx_m1]);
	  dxx /= dyy;
	  dxx *= dxx;
	  if (dxx > stat_high) {
	    success_2incr += 2;
	  } else if (dxx > stat_low) {
	    success_2incr++;
	  }
	} else {
	  dxx = (double)regression_results[(pidx + 1) * param_ctx_m1 - 1];
	  if (dxx > stat_high) {
	    success_2incr += 2;
	  } else if (dxx > stat_low) {
	    success_2incr++;
	  } else if (dxx == -9) {
	    perm_fail_ct++;
	  }
	}
	if (results[pidx] < dxx) {
	  results[pidx] = dxx;
	}
	if (msa_ptr) {
	  *msa_ptr++ = dxx;
	}
      } else if (msa_ptr) {
	*msa_ptr++ = -9;
      }
    }
    perm_2success_ct[marker_idx] += success_2incr;
    if (perm_fail_ct) {
      perm_fail_cts[marker_idx] += perm_fail_ct;
    }
  }
  THREAD_RETURN;
}

int32_t glm_common_init(FILE* bedfile, uintptr_t bed_offset, uint32_t glm_modifier, uint32_t standard_beta, uint32_t glm_xchr_model, Range_list* parameters_range_list_ptr, Range_list* tests_range_list_ptr, uint32_t mtest_adjust, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t max_marker_allele_len, uintptr_t* marker_reverse, char* condition_mname, char* condition_fname, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_indiv_ct, uintptr_t indiv_ct, uintptr_t* indiv_exclude, uint32_t mperm_save, uint32_t pheno_nm_ct, uintptr_t* pheno_nm, uintptr_t covar_ct, char* covar_names, uintptr_t max_covar_name_len, uintptr_t* covar_nm, double* covar_d, uintptr_t* sex_nm, uintptr_t* sex_male, uint32_t hh_exists, uint32_t do_perms, uint32_t* perm_batch_size_ptr, uintptr_t* max_param_name_len_ptr, uintptr_t* np_diploid_ptr, uintptr_t* condition_ct_ptr, uintptr_t* np_sex_ptr, uint32_t* marker_initial_ct_ptr, uint32_t* sex_covar_everywhere_ptr, uint32_t* variation_in_sex_ptr, uint32_t* x_sex_interaction_ptr, uintptr_t* constraint_ct_max_ptr, uintptr_t* param_idx_end_ptr, uintptr_t** loadbuf_raw_ptr, uintptr_t** load_mask_ptr, uintptr_t** sex_male_collapsed_ptr, uintptr_t** indiv_include2_ptr, uintptr_t** indiv_male_include2_ptr, uintptr_t** active_params_ptr, uintptr_t** haploid_params_ptr, uint32_t** condition_uidxs_ptr, char** writebuf_ptr, uintptr_t* indiv_valid_ct_ptr, uintptr_t* param_ctx_max_ptr, uintptr_t* param_ctl_max_ptr, uintptr_t* condition_list_start_idx_ptr, uintptr_t* covar_start_idx_ptr, uintptr_t* interaction_start_idx_ptr, uintptr_t* sex_start_idx_ptr, uintptr_t* np_base_ptr, uintptr_t* param_ct_max_ptr) {
  uintptr_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + BITCT - 1) / BITCT;
  uintptr_t unfiltered_indiv_ctv2 = 2 * unfiltered_indiv_ctl;
  uintptr_t indiv_uidx = 0;
  uintptr_t topsize = 0;
  uintptr_t max_param_name_len = 2;
  uintptr_t np_base_raw = 2; // intercept, additive effect
  uintptr_t np_diploid_raw = 0; // genotypic, hethom
  uintptr_t np_sex_raw = 0; // X-only sex covariate, sex-dosage interaction
  uintptr_t np_diploid = 0;
  uintptr_t constraint_ct_max = 0;
  uintptr_t condition_ct = 0;
  uintptr_t np_sex = 0;
  uintptr_t param_idx_end = 2;
  uintptr_t* load_mask = NULL;
  uintptr_t* indiv_include2 = NULL;
  uintptr_t* indiv_male_include2 = NULL;
  uint32_t* condition_uidxs = NULL;
  uint32_t covar_interactions = (glm_modifier / GLM_INTERACTION) & 1;
  uint32_t genotypic = glm_modifier & GLM_GENOTYPIC;
  uint32_t genotypic_or_hethom = (glm_modifier & (GLM_GENOTYPIC | GLM_HETHOM))? 1 : 0;
  uint32_t marker_initial_ct = marker_ct;
  uint32_t slen_add = 0;
  uint32_t sex_covar_everywhere = glm_modifier & GLM_SEX;
  uint32_t x_sex_interaction = (glm_xchr_model == 3);
  uint32_t x_present = (chrom_info_ptr->x_code != -1) && is_set(chrom_info_ptr->chrom_mask, chrom_info_ptr->x_code);
  uint32_t hide_covar = glm_modifier & GLM_HIDE_COVAR;
  uint32_t variation_in_sex = 0; // zero if no-x-sex specified
  int32_t retval = 0;
  uintptr_t* loadbuf_raw;
  uintptr_t* sex_male_collapsed;
  uintptr_t* active_params;
  uintptr_t* haploid_params;
  uintptr_t indiv_valid_ct;
  uintptr_t indiv_valid_ctv2;
  uintptr_t indiv_uidx_stop;
  uintptr_t indiv_idx;
  uintptr_t param_raw_ct_max;
  uintptr_t param_raw_ctl;
  uintptr_t param_ctx_max;
  uintptr_t param_ctl_max;
  uintptr_t condition_list_start_idx;
  uintptr_t covar_start_idx;
  uintptr_t interaction_start_idx;
  uintptr_t sex_start_idx;
  uintptr_t np_base;
  uintptr_t param_ct_max;
  uintptr_t param_idx;
  uint32_t slen;
  uint32_t uii;
  uint32_t ujj;
  g_joint_test_params = NULL;
  g_indiv_to_cluster1 = NULL;
  if (max_marker_allele_len > MAXLINELEN) {
    if (wkspace_alloc_c_checked(writebuf_ptr, max_marker_allele_len + MAXLINELEN)) {
      goto glm_common_init_ret_NOMEM;
    }
  } else {
    *writebuf_ptr = tbuf;
  }
  g_standard_beta = standard_beta;
  g_coding_flags = glm_modifier & (GLM_HETHOM | GLM_DOMINANT | GLM_RECESSIVE);
  g_glm_xchr_model = glm_xchr_model;
  g_include_sex = 0;
  g_male_x_01 = 0;
  if (!glm_xchr_model) {
    if (is_set(chrom_info_ptr->haploid_mask, 0)) {
      logprint("Error: --xchr-model 0 cannot be used with haploid genomes.\n");
      goto glm_common_init_ret_INVALID_CMDLINE;
    }
    uii = count_non_autosomal_markers(chrom_info_ptr, marker_exclude, 1, 1);
    if (uii) {
      LOGPRINTF("Excluding %u nonautosomal variant%s from --linear/--logistic analysis\n(--xchr-model 0).\n", uii, (uii == 1)? "" : "s");
      marker_initial_ct -= uii;
      if (!marker_initial_ct) {
	logprint("Error: No variants remaining for --linear/--logistic analysis.\n");
	goto glm_common_init_ret_INVALID_CMDLINE;
      }
    }
  }
  if (glm_init_load_mask(indiv_exclude, pheno_nm, covar_nm, indiv_ct, unfiltered_indiv_ctv2, &load_mask)) {
    goto glm_common_init_ret_NOMEM;
  }
  if (wkspace_alloc_ul_checked(&loadbuf_raw, unfiltered_indiv_ctv2 * sizeof(intptr_t))) {
    goto glm_common_init_ret_NOMEM;
  }
  loadbuf_raw[unfiltered_indiv_ctv2 - 2] = 0;
  loadbuf_raw[unfiltered_indiv_ctv2 - 1] = 0;
  indiv_valid_ct = popcount_longs(load_mask, unfiltered_indiv_ctl);
  if (condition_mname || condition_fname) {
    // temporary allocation of unfiltered indiv_include2 and
    // indiv_male_include2 for glm_scan_conditions()
    if (hh_exists & (Y_FIX_NEEDED | NXMHH_EXISTS)) {
      indiv_include2 = (uintptr_t*)top_alloc(&topsize, unfiltered_indiv_ctv2 * sizeof(intptr_t));
      if (!indiv_include2) {
	goto glm_common_init_ret_NOMEM;
      }
      fill_vec_55(indiv_include2, unfiltered_indiv_ct);
    }
    if (hh_exists & (XMHH_EXISTS | Y_FIX_NEEDED)) {
      indiv_male_include2 = (uintptr_t*)top_alloc(&topsize, unfiltered_indiv_ctv2 * sizeof(intptr_t));
      if (!indiv_male_include2) {
        goto glm_common_init_ret_NOMEM;
      }
      fill_ulong_zero(indiv_male_include2, unfiltered_indiv_ctv2);
      vec_include_init(unfiltered_indiv_ct, indiv_male_include2, sex_male);
    }
    wkspace_left -= topsize;
    retval = glm_scan_conditions(condition_mname, condition_fname, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, chrom_info_ptr, hh_exists, loadbuf_raw, bedfile, bed_offset, unfiltered_indiv_ct, sex_male, load_mask, &indiv_valid_ct, &condition_ct, &condition_uidxs, indiv_include2, indiv_male_include2);
    wkspace_left += topsize;
    if (retval) {
      goto glm_common_init_ret_1;
    }
    // topsize = 0;

    // need to set to null for next alloc_collapsed_haploid_filters() call to
    // work properly
    indiv_include2 = NULL;
    indiv_male_include2 = NULL;
    np_base_raw += condition_ct * (1 + covar_interactions);
  }
  np_base_raw += covar_ct * (1 + covar_interactions);
  if (!((glm_modifier & GLM_NO_X_SEX) || ((!sex_covar_everywhere) && (!x_present)))) {
    indiv_uidx = 0;
    indiv_idx = 0;
    uii = 2;
    // slight change from PLINK 1.07: missing sex no longer counts as female
    // when both --allow-no-sex and --sex are active (--linear/--logistic
    // normally throws out all individuals missing a covariate, so we shouldn't
    // make an exception for sex; --allow-no-sex may still have an effect on
    // other analyses in the current run)
    //
    // todo: implement this change for --allow-no-sex without --sex as well.
    do {
      indiv_uidx = next_set_ul_unsafe(load_mask, indiv_uidx);
      indiv_uidx_stop = next_unset_ul(load_mask, indiv_uidx, unfiltered_indiv_ct);
      indiv_idx += indiv_uidx_stop - indiv_uidx;
      do {
	if (IS_SET(sex_nm, indiv_uidx)) {
	  ujj = is_set(sex_male, indiv_uidx);
	  if (uii == ujj) {
	    variation_in_sex = 1;
	    indiv_idx = indiv_valid_ct;
	    break;
	  }
          uii = 1 - ujj;
	}
      } while (++indiv_uidx < indiv_uidx_stop);
    } while (indiv_idx < indiv_valid_ct);
    if (variation_in_sex) {
      if (sex_covar_everywhere) {
	np_base_raw += 1 + covar_interactions;
	if (covar_interactions) {
	  np_diploid_raw = genotypic_or_hethom;
	}
        bitfield_and(load_mask, sex_nm, unfiltered_indiv_ctl);
        indiv_valid_ct = popcount_longs(load_mask, unfiltered_indiv_ctl);
      } else {
	np_sex_raw = 1;
	if (covar_interactions) {
	  np_sex_raw += 1 + genotypic_or_hethom;
	}
      }
      if (x_sex_interaction) {
	if (covar_interactions) {
	  logprint("Note: Ignoring --xchr-model 3 since it's redundant with --linear\n'interaction' modifier.\n");
	  x_sex_interaction = 0;
	} else {
	  np_sex_raw++;
	}
      }
    } else {
      if (sex_covar_everywhere) {
        LOGPRINTFWW("Warning: Ignoring --linear/--logistic 'sex' modifier%s since sex is%cinvariant.\n", x_sex_interaction? " and --xchr-model 3" : "", x_sex_interaction? '\n' : ' ');
        sex_covar_everywhere = 0;
      } else if (x_sex_interaction) {
        logprint("Warning: Ignoring --xchr-model 3 since sex is invariant.\n");
      }
      x_sex_interaction = 0;
    }
  }
  if (genotypic_or_hethom) {
    np_diploid_raw++;
    if (covar_interactions) {
      np_diploid_raw += condition_ct + covar_ct;
    }
  }

  indiv_valid_ctv2 = 2 * ((indiv_valid_ct + BITCT - 1) / BITCT);
  if (alloc_collapsed_haploid_filters(unfiltered_indiv_ct, indiv_valid_ct, hh_exists, 1, load_mask, sex_male, &indiv_include2, &indiv_male_include2)) {
    goto glm_common_init_ret_NOMEM;
  }
  if (wkspace_alloc_ul_checked(&g_loadbuf, GLM_BLOCKSIZE * indiv_valid_ctv2 * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&sex_male_collapsed, indiv_valid_ctv2 * sizeof(intptr_t))) {
    goto glm_common_init_ret_NOMEM;
  }
  for (uii = 1; uii <= GLM_BLOCKSIZE; uii++) {
    g_loadbuf[uii * indiv_valid_ctv2 - 2] = 0;
    g_loadbuf[uii * indiv_valid_ctv2 - 1] = 0;
  }
  collapse_copy_bitarr_incl(unfiltered_indiv_ct, sex_male, load_mask, indiv_valid_ct, sex_male_collapsed);
  param_raw_ct_max = np_base_raw + np_diploid_raw + np_sex_raw;
  param_raw_ctl = (param_raw_ct_max + BITCT - 1) / BITCT;
  if (wkspace_alloc_ul_checked(&active_params, param_raw_ctl * sizeof(intptr_t))) {
    goto glm_common_init_ret_NOMEM;
  }
  condition_list_start_idx = 2 + genotypic_or_hethom;
  covar_start_idx = condition_list_start_idx + condition_ct;
  interaction_start_idx = covar_start_idx + covar_ct;
  if (!covar_interactions) {
    sex_start_idx = interaction_start_idx;
  } else {
    if (genotypic_or_hethom) {
      sex_start_idx = interaction_start_idx * 3 - 6;
    } else {
      sex_start_idx = interaction_start_idx * 2 - 2;
    }
  }
  g_condition_list_start_idx = condition_list_start_idx;
  g_interaction_start_idx = interaction_start_idx;
  g_sex_start_idx = sex_start_idx;
  g_active_params = active_params;
  g_sex_male_collapsed = sex_male_collapsed;

  if (parameters_range_list_ptr->name_ct) {
    fill_ulong_zero(active_params, param_raw_ctl);
    active_params[0] = 1;
    numeric_range_list_to_bitfield(parameters_range_list_ptr, param_raw_ct_max, active_params, 0, 1);
    if ((!(active_params[0] & 2)) && ((!np_diploid_raw) || (active_params[0] & 4)) && ((!covar_interactions) || ((!popcount_bit_idx(active_params, interaction_start_idx, sex_start_idx)) && ((!variation_in_sex) || (!popcount_bit_idx(active_params, sex_start_idx + 1, param_raw_ct_max)))))) {
      // force the user to explicitly use no-snp if that's their intention
      logprint("Error: --parameters must retain at least one dosage-dependent variable.  To\nperform one-off regression(s), use the --linear 'no-snp' modifier instead.\n");
      goto glm_common_init_ret_INVALID_CMDLINE;
    }
    param_ct_max = popcount_longs(active_params, param_raw_ctl);
    if (np_diploid_raw) {
      np_diploid = IS_SET(active_params, 2);
      if (covar_interactions) {
	for (param_idx = interaction_start_idx + 1; param_idx < sex_start_idx; param_idx += 2) {
          np_diploid += IS_SET(active_params, param_idx);
	}
      }
    }
    if (!sex_covar_everywhere) {
      np_sex = popcount_bit_idx(active_params, sex_start_idx, param_raw_ct_max);
    }
    np_base = param_ct_max - np_diploid - np_sex;
    if (!np_sex) {
      variation_in_sex = 0;
    }
  } else {
    fill_all_bits(active_params, param_raw_ct_max);
    param_ct_max = param_raw_ct_max;
    np_base = np_base_raw;
    np_diploid = np_diploid_raw;
    np_sex = np_sex_raw;
  }
  if (indiv_valid_ct <= param_ct_max) {
    logprint("Warning: Skipping --linear since # variables >= # samples.\n");
    if (pheno_nm_ct > param_ct_max) {
      logprint("(Check your covariates--all samples with at least one missing covariate are\nexcluded from this analysis.)\n");
    }
    goto glm_common_init_ret_1;
  }
  // parameter sequence:
  // 1. intercept
  // 2. allelic dosage
  // 3. dominance deviation
  // 4. --condition-list
  // 5. --covar
  // 6. interactions (if DOMDEV is present, both interactions with first
  //    covariate occur before either interaction with the second covariate,
  //    etc.
  // 7. sex
  // 8. sex interactions
  if (IS_SET(active_params, 1)) {
    max_param_name_len = 4;
  } else if (mtest_adjust) {
    logprint("Error: --adjust cannot be used when --parameters excludes the main effect.\n");
    goto glm_common_init_ret_INVALID_CMDLINE;
  }
  if (hide_covar) {
    if (genotypic_or_hethom && IS_SET(active_params, 2)) {
      if (IS_SET(active_params, 1)) {
        param_idx_end = 3;
      }
    } else if (!IS_SET(active_params, 1)) {
      if (tests_range_list_ptr->name_ct || (glm_modifier & GLM_TEST_ALL)) {
        param_idx_end = 1;
      } else {
        logprint("Error: 'hide-covar' modifier suppresses all output due to --parameters setting.\n");
        goto glm_common_init_ret_INVALID_CMDLINE;
      }
    }
  }

  if (genotypic_or_hethom) {
    if (genotypic) {
      slen_add = 8; // "DOMDEVx"...
    } else {
      slen_add = 5; // "HETx"...
    }
    if (IS_SET(active_params, 2)) {
      max_param_name_len = slen_add - 1;
    }
  }
  for (uii = 0; uii < condition_ct; uii++) {
    if (is_set(active_params, uii + condition_list_start_idx)) {
      slen = strlen(&(marker_ids[condition_uidxs[uii] * max_marker_id_len]));
      if (max_param_name_len <= slen) {
	max_param_name_len = slen + 1;
      }
    }
    if (covar_interactions) {
      // ugh, special case
      // "CSNP"+int2str(c+1)
      slen = 4 + intlen(uii + 1);
      if (!genotypic_or_hethom) {
	if (is_set(active_params, interaction_start_idx + uii)) {
	  if (max_param_name_len < slen + 5) {
	    max_param_name_len = slen + 5;
	  }
	}
      } else {
	if (is_set(active_params, interaction_start_idx + 2 * uii)) {
	  if (max_param_name_len < slen + 5) {
	    max_param_name_len = slen + 5;
	  }
	}
	if (is_set(active_params, interaction_start_idx + 2 * uii + 1)) {
	  if (max_param_name_len < slen + slen_add) {
	    max_param_name_len = slen + slen_add;
	  }
	}
      }
    }
  }
  for (uii = 0; uii < covar_ct; uii++) {
    slen = strlen(&(covar_names[uii * max_covar_name_len]));
    if (is_set(active_params, uii + covar_start_idx)) {
      if (max_param_name_len <= slen) {
	max_param_name_len = slen + 1;
      }
    }
    if (covar_interactions) {
      if (!genotypic_or_hethom) {
	if (is_set(active_params, interaction_start_idx + uii + condition_ct)) {
	  if (max_param_name_len < slen + 5) {
	    max_param_name_len = slen + 5;
	  }
	}
      } else {
	if (is_set(active_params, interaction_start_idx + 2 * (uii + condition_ct))) {
	  if (max_param_name_len < slen + 5) {
	    max_param_name_len = slen + 5;
	  }
	}
	if (is_set(active_params, interaction_start_idx + 2 * (uii + condition_ct) + 1)) {
	  if (max_param_name_len < slen + slen_add) {
	    max_param_name_len = slen + slen_add;
	  }
	}
      }
    }
  }
  if (variation_in_sex) {
    if (is_set(active_params, sex_start_idx)) {
      if (max_param_name_len < 4) {
	max_param_name_len = 4;
      }
    }
    if (covar_interactions) {
      if (is_set(active_params, sex_start_idx + 1)) {
        if (max_param_name_len < 8) {
	  max_param_name_len = 8; // "ADDxSEX", etc.
	}
      }
      if (genotypic_or_hethom) {
	if (is_set(active_params, sex_start_idx + 2)) {
          if (max_param_name_len < 3 + slen_add) {
	    max_param_name_len = 3 + slen_add; // "DOMDEVxSEX", "HETxSEX"
	  }
	}
      }
    } else if (x_sex_interaction) {
      // --xchr-model 3, use "XxSEX" name instead
      if (is_set(active_params, sex_start_idx + 1)) {
        if (max_param_name_len < 6) {
	  max_param_name_len = 6;
	}
      } else {
	x_sex_interaction = 0;
      }
    }
  }
  param_ctx_max = param_ct_max;
  param_ctl_max = (param_ct_max + BITCT - 1) / BITCT;
  if (wkspace_alloc_ul_checked(&haploid_params, param_ctl_max * sizeof(intptr_t))) {
    goto glm_common_init_ret_NOMEM;
  }
  g_haploid_params = haploid_params;
  if (genotypic_or_hethom) {
    fill_ulong_zero(haploid_params, param_ctl_max);
    ujj = np_base;
    for (uii = 0, param_idx = 0; uii < ujj; uii++, param_idx++) {
      next_set_unsafe_ck(active_params, &uii);
      if ((uii != 2) && ((uii < interaction_start_idx) || (!((uii - interaction_start_idx) & 1)))) {
        SET_BIT(haploid_params, param_idx);
      }
    }
  } else {
    fill_all_bits(haploid_params, param_ct_max - np_sex);
  }
  uii = 0;
  if ((genotypic_or_hethom && ((active_params[0] & 6) == 6)) || tests_range_list_ptr->name_ct || (glm_modifier & GLM_TEST_ALL)) {
    if (wkspace_alloc_ul_checked(&g_joint_test_params, param_ctl_max * sizeof(intptr_t))) {
      goto glm_common_init_ret_NOMEM;
    }
    fill_ulong_zero(g_joint_test_params, param_ctl_max);
    if (tests_range_list_ptr->name_ct) {
      numeric_range_list_to_bitfield(tests_range_list_ptr, param_ct_max - 1, g_joint_test_params, 1, 1);
      constraint_ct_max = popcount_longs(g_joint_test_params, param_ctl_max);
    } else if (glm_modifier & GLM_TEST_ALL) {
      constraint_ct_max = param_ct_max - 1;
      fill_bits(g_joint_test_params, 0, constraint_ct_max);
    } else {
      // genotypic/hethom, neither of first two terms excluded by --parameters,
      // no --tests
      constraint_ct_max = 2;
      g_joint_test_params[0] = 3;
      uii = 1; // flag this to prevent 'interaction' permutation test
    }
    if ((constraint_ct_max > 1) || ((constraint_ct_max == 1) && (hide_covar || (covar_interactions && do_perms)))) {
      // "USER_2DF", etc.
      // Allow 1df if:
      // - 'hide-covar' is in effect (so only this covariate is reported on)
      // - and/or 'interaction' is being used with a permutation test
      // In the latter case, we don't stop the user from shooting themselves in
      // the foot with an inappropriate permutation test target; forcing them
      // to use --tests should prevent enough mistakes.
      slen = 8 + intlen(constraint_ct_max);
      if (max_param_name_len < slen) {
	max_param_name_len = slen;
      }
      param_ctx_max++;
    } else {
      wkspace_reset(g_joint_test_params);
      g_joint_test_params = NULL;
      constraint_ct_max = 0;
      logprint("Warning: Ignoring --tests since too few parameter indices are in range.\n");
    }
  }
  if (covar_interactions && do_perms && ((!constraint_ct_max) || uii)) {
    logprint("Error: --linear/--logistic 'interaction' modifier cannot be used with\npermutation except with --tests.\n");
    goto glm_common_init_ret_INVALID_CMDLINE;
  }
  if (do_perms && (!IS_SET(active_params, 1)) && (!constraint_ct_max)) {
    logprint("Error: --linear/--logistic permutation test cannot occur when --parameters\nexcludes the main effect and no joint test is active.\n");
    goto glm_common_init_ret_INVALID_CMDLINE;
  }

  g_cluster_ct = 0;
  g_pheno_nm_ct = indiv_valid_ct;
  g_perms_done = 0;
  g_mperm_save_all = NULL;
  if (!do_perms) {
    *perm_batch_size_ptr = 1;
  }
  *loadbuf_raw_ptr = loadbuf_raw;
  *load_mask_ptr = load_mask;
  *sex_male_collapsed_ptr = sex_male_collapsed;
  *constraint_ct_max_ptr = constraint_ct_max;
  *max_param_name_len_ptr = max_param_name_len;
  *np_diploid_ptr = np_diploid;
  *condition_ct_ptr = condition_ct;
  *np_sex_ptr = np_sex;
  *param_idx_end_ptr = param_idx_end;
  *marker_initial_ct_ptr = marker_initial_ct;
  *sex_covar_everywhere_ptr = sex_covar_everywhere;
  *variation_in_sex_ptr = variation_in_sex;
  *x_sex_interaction_ptr = x_sex_interaction;
  *indiv_valid_ct_ptr = indiv_valid_ct;
  *indiv_include2_ptr = indiv_include2;
  *indiv_male_include2_ptr = indiv_male_include2;
  *active_params_ptr = active_params;
  *haploid_params_ptr = haploid_params;
  *condition_uidxs_ptr = condition_uidxs;
  *param_ctx_max_ptr = param_ctx_max;
  *param_ctl_max_ptr = param_ctl_max;
  *condition_list_start_idx_ptr = condition_list_start_idx;
  *covar_start_idx_ptr = covar_start_idx;
  *interaction_start_idx_ptr = interaction_start_idx;
  *sex_start_idx_ptr = sex_start_idx;
  *np_base_ptr = np_base;
  *param_ct_max_ptr = param_ct_max;
  while (0) {
  glm_common_init_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  glm_common_init_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  }
 glm_common_init_ret_1:
  return retval;
}

#ifndef NOLAPACK
int32_t glm_linear_assoc(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uint32_t glm_modifier, double glm_vif_thresh, uint32_t glm_xchr_model, uint32_t glm_mperm_val, Range_list* parameters_range_list_ptr, Range_list* tests_range_list_ptr, double ci_size, double ci_zt, double pfilter, uint32_t mtest_adjust, double adjust_lambda, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, uint32_t* marker_pos, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, uintptr_t* marker_reverse, uint32_t zero_extra_chroms, char* condition_mname, char* condition_fname, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_indiv_ct, uintptr_t indiv_ct, uintptr_t* indiv_exclude, uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, Aperm_info* apip, uint32_t mperm_save, uint32_t pheno_nm_ct, uintptr_t* pheno_nm, double* pheno_d, uintptr_t covar_ct, char* covar_names, uintptr_t max_covar_name_len, uintptr_t* covar_nm, double* covar_d, uintptr_t* sex_nm, uintptr_t* sex_male, uint32_t hh_exists, uint32_t perm_batch_size, Set_info* sip) {
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  FILE* outfile = NULL;
  FILE* outfile_msa = NULL;
  uintptr_t indiv_uidx = 0;
  uintptr_t cur_constraint_ct = 0;
  uintptr_t cur_param_ct = 0;
  uintptr_t cur_param_ctx = 0;
  uint32_t perms_total = 0;
  uint32_t cluster_ct1 = 0;
  uint32_t perm_adapt = glm_modifier & GLM_PERM;
  uint32_t perm_maxt = glm_modifier & GLM_MPERM;
  uint32_t covar_interactions = (glm_modifier / GLM_INTERACTION) & 1;
  uint32_t hethom = glm_modifier & GLM_HETHOM;
  uint32_t standard_beta = glm_modifier & GLM_STANDARD_BETA;
  uint32_t genotypic_or_hethom = (glm_modifier & (GLM_GENOTYPIC | GLM_HETHOM))? 1 : 0;
  uint32_t do_perms = perm_adapt | perm_maxt;
  uint32_t perm_count = glm_modifier & GLM_PERM_COUNT;
  uint32_t hide_covar = glm_modifier & GLM_HIDE_COVAR;
  uint32_t mperm_save_all = mperm_save & MPERM_DUMP_ALL;
  uint32_t fill_orig_chiabs = do_perms || mtest_adjust;
  uint32_t display_ci = (ci_size > 0);
  uint32_t perm_pass_idx = 0;
  uint32_t perm_fail_ct = 0;
  uint32_t pct = 0;
  uint32_t max_thread_ct = g_thread_ct;
  int32_t retval = 0;
  uint32_t linear_intercept = glm_modifier & GLM_INTERCEPT;
  char dgels_trans = 'N';
  __CLPK_integer dgels_m = 0;
  __CLPK_integer dgels_n = 0;
  __CLPK_integer dgels_nrhs = 0;
  __CLPK_integer dgels_ldb = 0;
  uintptr_t* ulptr;
  double* dptr3;
  uintptr_t cur_word;
  __CLPK_integer dgels_info;
  double* constraints_con_major = NULL;
  uint32_t* condition_uidxs = NULL;
  uint32_t* cluster_map1 = NULL;
  uint32_t* cluster_starts1 = NULL;
  uint32_t* tcnt = NULL;
  uint32_t* marker_idx_to_uidx = NULL;
  uint32_t* cur_indiv_to_cluster1 = NULL;
  char* cur_param_names = NULL;
  char* haploid_param_names = NULL;
  char* wptr_start = NULL;
  double geno_map[12];
  uint32_t mu_table[GLM_BLOCKSIZE];
  double* geno_map_ptr;
  char* writebuf;
  char* param_names;
  const char* main_effect;

  // linear: pval = calc_tprob(g_orig_stats[i], indiv_valid_ct - param_ct)
  double* orig_stats_ptr;

  char* outname_end2;
  char* wptr;
  char* wptr_start2;
  char* cptr;
  double* dptr;
  double* dptr2;
  uintptr_t* loadbuf_raw;
  uintptr_t* load_mask;
  uintptr_t* sex_male_collapsed;
  uintptr_t* indiv_include2;
  uintptr_t* indiv_male_include2;
  uintptr_t* active_params;
  uintptr_t* haploid_params;
  uintptr_t* loadbuf_ptr;
  uintptr_t max_param_name_len;
  uintptr_t np_diploid;
  uintptr_t condition_ct;
  uintptr_t constraint_ct_max;
  uintptr_t np_sex;
  uintptr_t param_idx_end;
  uintptr_t indiv_valid_ct;
  uintptr_t indiv_valid_ctv2;
  uintptr_t indiv_uidx_stop;
  uintptr_t indiv_idx;
  uintptr_t param_ctx_max;
  uintptr_t param_ctl_max;
  uintptr_t param_ctx_max_m1;
  uintptr_t condition_list_start_idx;
  uintptr_t covar_start_idx;
  uintptr_t interaction_start_idx;
  uintptr_t sex_start_idx;
  uintptr_t np_base;
  uintptr_t param_ct_max;
  uintptr_t cur_missing_ct;
  uintptr_t cur_indiv_valid_ct;
  uintptr_t param_idx;
  uintptr_t param_idx_fixed;
  uintptr_t constraint_idx;
  uintptr_t ulii;
  uintptr_t uljj;
  double* msa_ptr;
  double se;
  double zval;
  double pval;
  double dxx;
  double dyy;
  double dzz;
  uint32_t marker_initial_ct;
  uint32_t sex_covar_everywhere;
  uint32_t variation_in_sex;
  uint32_t x_sex_interaction;
  uint32_t marker_unstopped_ct;
  uint32_t chrom_fo_idx;
  uint32_t chrom_end;
  uint32_t block_size;
  uint32_t block_end;
  uint32_t tidx;
  uint32_t perm_idx;
  uint32_t marker_uidx;
  uint32_t marker_uidx2;
  uint32_t marker_idx;
  uint32_t marker_idx2;
  uint32_t marker_idx3;
  uint32_t marker_bidx;
  uint32_t chrom_idx;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t loop_end;
  uint32_t regression_fail;

  retval = glm_common_init(bedfile, bed_offset, glm_modifier, standard_beta, glm_xchr_model, parameters_range_list_ptr, tests_range_list_ptr, mtest_adjust, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, max_marker_allele_len, marker_reverse, condition_mname, condition_fname, chrom_info_ptr, unfiltered_indiv_ct, indiv_ct, indiv_exclude, mperm_save, pheno_nm_ct, pheno_nm, covar_ct, covar_names, max_covar_name_len, covar_nm, covar_d, sex_nm, sex_male, hh_exists, do_perms, &perm_batch_size, &max_param_name_len, &np_diploid, &condition_ct, &np_sex, &marker_initial_ct, &sex_covar_everywhere, &variation_in_sex, &x_sex_interaction, &constraint_ct_max, &param_idx_end, &loadbuf_raw, &load_mask, &sex_male_collapsed, &indiv_include2, &indiv_male_include2, &active_params, &haploid_params, &condition_uidxs, &writebuf, &indiv_valid_ct, &param_ctx_max, &param_ctl_max, &condition_list_start_idx, &covar_start_idx, &interaction_start_idx, &sex_start_idx, &np_base, &param_ct_max);
  if (retval) {
    goto glm_linear_assoc_ret_1;
  }
  indiv_valid_ctv2 = 2 * ((indiv_valid_ct + BITCT - 1) / BITCT);
  param_ctx_max_m1 = param_ctx_max - 1;
  if (wkspace_alloc_d_checked(&g_orig_stats, marker_initial_ct * sizeof(double)) ||
      wkspace_alloc_c_checked(&param_names, param_ctx_max * max_param_name_len) ||
      wkspace_alloc_d_checked(&g_fixed_covars_cov_major, (variation_in_sex + interaction_start_idx - condition_list_start_idx) * indiv_valid_ct * sizeof(double)) ||
      wkspace_alloc_uc_checked(&g_perm_adapt_stop, marker_initial_ct) ||
      wkspace_alloc_ui_checked(&g_nm_cts, marker_initial_ct * sizeof(int32_t))) {
    goto glm_linear_assoc_ret_NOMEM;
  }
  // use this array to track regression failures even in max(T) case
  fill_ulong_zero((uintptr_t*)g_perm_adapt_stop, (marker_initial_ct + sizeof(intptr_t) - 1) / sizeof(intptr_t));

  param_idx = 1;
  if (glm_modifier & GLM_RECESSIVE) {
    main_effect = glm_main_effects;
  } else if (glm_modifier & GLM_DOMINANT) {
    main_effect = &(glm_main_effects[4]);
  } else if (hethom) {
    main_effect = &(glm_main_effects[8]);
  } else {
    main_effect = &(glm_main_effects[12]);
  }
  if (IS_SET(active_params, 1)) {
    memcpy(&(param_names[max_param_name_len]), main_effect, 4);
    param_idx++;
  }
  if (genotypic_or_hethom && IS_SET(active_params, 2)) {
    if (hethom) {
      memcpy(&(param_names[param_idx * max_param_name_len]), "HET", 4);
    } else {
      memcpy(&(param_names[param_idx * max_param_name_len]), "DOMDEV", 7);
    }
    param_idx++;
  }
  // 0..3: diploid chromosomes, X chromosome female
  // 4..7: X chromosome male
  // 8..11: haploid
  fill_double_zero(geno_map, 12);
  geno_map[0] = 1;
  geno_map[2] = 1;
  geno_map[4] = 1;
  geno_map[8] = 1;
  if (glm_modifier & GLM_CONDITION_RECESSIVE) {
    geno_map[2] = 0;
  } else if (!(glm_modifier & GLM_CONDITION_DOMINANT)) {
    geno_map[0] = 2;
    if (glm_xchr_model == 2) {
      geno_map[4] = 2;
    }
  }
  for (param_idx_fixed = 0; param_idx_fixed < condition_ct; param_idx_fixed++) {
    marker_uidx = condition_uidxs[param_idx_fixed];
    if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
      goto glm_linear_assoc_ret_READ_FAIL;
    }
    if (load_and_collapse_incl(bedfile, loadbuf_raw, unfiltered_indiv_ct, g_loadbuf, indiv_valid_ct, load_mask, IS_SET(marker_reverse, marker_uidx))) {
      goto glm_linear_assoc_ret_READ_FAIL;
    }
    chrom_idx = get_marker_chrom(chrom_info_ptr, marker_uidx);
    geno_map_ptr = geno_map;
    if (IS_SET(chrom_info_ptr->haploid_mask, chrom_idx)) {
      g_is_x = ((int32_t)chrom_idx == chrom_info_ptr->x_code);
      g_is_y = ((int32_t)chrom_idx == chrom_info_ptr->y_code);
      if (hh_exists) {
	haploid_fix(hh_exists, indiv_include2, indiv_male_include2, indiv_valid_ct, g_is_x, g_is_y, (unsigned char*)g_loadbuf);
      }
      if (!g_is_x) {
	geno_map_ptr = &(geno_map[8]);
      }
    } else {
      g_is_x = 0;
    }
    if (!g_is_x) {
      glm_loadbuf_to_doubles(g_loadbuf, indiv_valid_ct, &(g_fixed_covars_cov_major[param_idx_fixed * indiv_valid_ct]), geno_map_ptr, NULL);
    } else {
      glm_loadbuf_to_doubles_x(g_loadbuf, sex_male_collapsed, indiv_valid_ct, &(g_fixed_covars_cov_major[param_idx_fixed * indiv_valid_ct]), geno_map_ptr, NULL);
    }
    if (is_set(active_params, param_idx_fixed + condition_list_start_idx)) {
      strcpy(&(param_names[param_idx * max_param_name_len]), &(marker_ids[marker_uidx * max_marker_id_len]));
      param_idx++;
    }
  }
  for (uii = 0; uii < covar_ct; param_idx_fixed++, uii++) {
    // indiv-major to covariate-major
    indiv_uidx = 0;
    dptr = &(g_fixed_covars_cov_major[param_idx_fixed * indiv_valid_ct]);
    dptr2 = &(covar_d[uii]);
    for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_uidx++, indiv_idx++) {
      next_unset_ul_unsafe_ck(indiv_exclude, &indiv_uidx);
      if (IS_SET(load_mask, indiv_uidx)) {
        *dptr++ = dptr2[indiv_idx * covar_ct];
      }
    }
  }
  for (uii = 0; uii < covar_ct; uii++) {
    if (is_set(active_params, uii + covar_start_idx)) {
      strcpy(&(param_names[param_idx * max_param_name_len]), &(covar_names[uii * max_covar_name_len]));
      param_idx++;
    }
  }
  if (covar_interactions) {
    ujj = interaction_start_idx;
    for (uii = 0; uii < condition_ct; uii++) {
      if (is_set(active_params, ujj++)) {
        wptr = memcpyl3a(&(param_names[param_idx * max_param_name_len]), main_effect);
	wptr = memcpya(wptr, "xCSNP", 5);
        uint32_writex(wptr, uii + 1, '\0');
	param_idx++;
      }
      if (genotypic_or_hethom) {
	if (is_set(active_params, ujj++)) {
	  wptr = &(param_names[param_idx * max_param_name_len]);
	  if (hethom) {
            wptr = memcpya(wptr, "HETxCSNP", 8);
	  } else {
	    wptr = memcpya(wptr, "DOMDEVxCSNP", 11);
	  }
	  uint32_writex(wptr, uii + 1, '\0');
	  param_idx++;
	}
      }
    }
    for (uii = 0; uii < covar_ct; uii++) {
      if (is_set(active_params, ujj++)) {
        wptr = memcpya(&(param_names[param_idx * max_param_name_len]), main_effect, 3);
	*wptr++ = 'x';
	strcpy(wptr, &(covar_names[uii * max_covar_name_len]));
        param_idx++;
      }
      if (genotypic_or_hethom) {
	if (is_set(active_params, ujj++)) {
	  wptr = &(param_names[param_idx * max_param_name_len]);
	  if (hethom) {
	    wptr = memcpya(wptr, "HETx", 4);
	  } else {
	    wptr = memcpya(wptr, "DOMDEVx", 7);
	  }
	  strcpy(wptr, &(covar_names[uii * max_covar_name_len]));
	  param_idx++;
	}
      }
    }
  }
  if (variation_in_sex) {
    dptr = &(g_fixed_covars_cov_major[param_idx_fixed * indiv_valid_ct]);
    fill_double_zero(dptr, indiv_valid_ct);
    indiv_idx = 0;
    while (1) {
      next_set_ul_ck(sex_male_collapsed, &indiv_idx, indiv_valid_ct);
      if (indiv_idx == indiv_valid_ct) {
	break;
      }
      dptr[indiv_idx++] = 1;
    }
    if (IS_SET(active_params, sex_start_idx)) {
      memcpy(&(param_names[param_idx * max_param_name_len]), "SEX", 4);
      param_idx++;
    }
    if (covar_interactions) {
      if (is_set(active_params, sex_start_idx + 1)) {
        wptr = memcpyl3a(&(param_names[param_idx * max_param_name_len]), main_effect);
        memcpy(wptr, "xSEX", 5);
	param_idx++;
      }
      if (genotypic_or_hethom && is_set(active_params, sex_start_idx + 2)) {
	wptr = &(param_names[param_idx * max_param_name_len]);
	if (hethom) {
	  memcpy(wptr, "HETxSEX", 8);
	} else {
	  memcpy(wptr, "DOMDEVxSEX", 11);
	}
	param_idx++;
      }
    } else if (x_sex_interaction) {
      // active_params bit must be set
      memcpy(&(param_names[param_idx * max_param_name_len]), "XxSEX", 6);
      param_idx++;
    }
  }
  if (genotypic_or_hethom) {
    if (wkspace_alloc_c_checked(&haploid_param_names, np_base * max_param_name_len)) {
      goto glm_linear_assoc_ret_NOMEM;
    }
    uii = 1;
    for (param_idx = 1; param_idx < np_base; uii++, param_idx++) {
      next_set_unsafe_ck(haploid_params, &uii);
      strcpy(&(haploid_param_names[param_idx * max_param_name_len]), &(param_names[uii * max_param_name_len]));
    }
  }

  if (constraint_ct_max) {
    if (wkspace_alloc_d_checked(&constraints_con_major, constraint_ct_max * param_ct_max * sizeof(double))) {
      goto glm_linear_assoc_ret_NOMEM;
    }
    // special case: df may vary between chromosomes, so refill suffix at
    // beginning of each chromosome
    wptr = &(param_names[param_ct_max * max_param_name_len]);
    if (tests_range_list_ptr->name_ct) {
      memcpy(wptr, "USER_", 6);
    } else if (glm_modifier & GLM_TEST_ALL) {
      memcpy(wptr, "FULL_", 6);
    } else {
      memcpy(wptr, "GENO_", 6);
    }
  }
  g_constraints_con_major = constraints_con_major;
  if (!perm_maxt) {
    g_aperm_alpha = apip->alpha;
    mperm_save = 0;
  }
  if (fill_orig_chiabs) {
    if (mtest_adjust) {
      if (wkspace_alloc_d_checked(&g_orig_chisq, marker_initial_ct * sizeof(double)) ||
          wkspace_alloc_ui_checked(&tcnt, marker_initial_ct * sizeof(int32_t)) ||
          wkspace_alloc_ui_checked(&marker_idx_to_uidx, marker_initial_ct * sizeof(int32_t))) {
	goto glm_linear_assoc_ret_NOMEM;
      }
    }
    if (do_perms) {
      if (!perm_batch_size) {
	perm_batch_size = 512;
      }
      if (wkspace_alloc_ui_checked(&g_perm_2success_ct, marker_initial_ct * sizeof(int32_t)) ||
	  wkspace_alloc_ui_checked(&g_perm_attempt_ct, marker_initial_ct * sizeof(int32_t))) {
	goto glm_linear_assoc_ret_NOMEM;
      }
      // need this for max(T) now since we need to track permutation failures
      // bugfix: g_perm_2success_ct was uninitialized.  add a
      // wkspace_calloc_...() idiom to reduce the frequency of that mistake?
      fill_uint_zero(g_perm_2success_ct, marker_initial_ct);
      fill_uint_zero(g_perm_attempt_ct, marker_initial_ct);
      if (perm_adapt) {
	perms_total = apip->max;
      } else {
	perms_total = glm_mperm_val;
      }
      if (perms_total < perm_batch_size) {
	perm_batch_size = perms_total;
      }
      uii = perm_batch_size / CACHELINE_INT32;
      if (!uii) {
	uii = 1;
      } else if (uii > GLM_BLOCKSIZE / CACHELINE_INT32) {
	uii = GLM_BLOCKSIZE / CACHELINE_INT32;
      }
      if (max_thread_ct > uii) {
	max_thread_ct = uii;
      }
      if (!perm_adapt) {
	ulii = (perm_batch_size + (CACHELINE_DBL - 1)) & (~(CACHELINE_DBL - 1));
        if (wkspace_alloc_d_checked(&g_maxt_thread_results, ulii * max_thread_ct * sizeof(double)),
            wkspace_alloc_d_checked(&g_maxt_extreme_stat, perms_total * sizeof(double))) {
          goto glm_linear_assoc_ret_NOMEM;
	}
	fill_double_zero(g_maxt_extreme_stat, perms_total);
	if (mperm_save_all) {
	  if (wkspace_alloc_d_checked(&g_mperm_save_all, marker_initial_ct * perm_batch_size * sizeof(double))) {
	    goto glm_linear_assoc_ret_NOMEM;
	  }
	  memcpy(outname_end, ".mperm.dump.all", 16);
	  if (fopen_checked(&outfile_msa, outname, "w")) {
	    goto glm_linear_assoc_ret_OPEN_FAIL;
	  }
	  if (putc_checked('0', outfile_msa)) {
	    goto glm_linear_assoc_ret_WRITE_FAIL;
	  }
	  LOGPRINTF("Dumping all permutation %s to %s .\n", (!constraint_ct_max)? "absolute t-stats" : "chi-square values", outname);
	}
      }
      g_perm_batch_max = perm_batch_size;
    }
  }
  if (cluster_starts) {
    // If there are any size-1 clusters, we actually want two cluster indexes:
    // - one for permutation, which excludes the size-1 clusters, and
    // - one for use by the robust cluster variance estimator, which does not.
    //
    // cluster_ct1, etc. includes size-1 clusters, while the same variables
    // without the "1" at the end exclude them.
    retval = cluster_include_and_reindex(unfiltered_indiv_ct, load_mask, 0, NULL, indiv_valid_ct, 0, cluster_ct, cluster_map, cluster_starts, &cluster_ct1, &cluster_map1, &cluster_starts1, NULL, NULL);
    if (retval) {
      goto glm_linear_assoc_ret_1;
    }
    g_cluster_ct1 = cluster_ct1;
    if (cluster_ct1) {
      if (wkspace_alloc_ui_checked(&g_indiv_to_cluster1, indiv_valid_ct * sizeof(int32_t))) {
	goto glm_linear_assoc_ret_NOMEM;
      }
      fill_unfiltered_indiv_to_cluster(indiv_valid_ct, cluster_ct1, cluster_map1, cluster_starts1, g_indiv_to_cluster1);
      if (do_perms) {
	retval = cluster_include_and_reindex(unfiltered_indiv_ct, load_mask, 1, NULL, indiv_valid_ct, 0, cluster_ct, cluster_map, cluster_starts, &g_cluster_ct, &g_cluster_map, &g_cluster_starts, NULL, NULL);
	if (retval) {
	  goto glm_linear_assoc_ret_1;
	}
	if (!g_cluster_ct) {
	  goto glm_linear_assoc_ret_NO_PERMUTATION_CLUSTERS;
	}
	if (wkspace_alloc_ui_checked(&g_qassoc_cluster_thread_wkspace, max_thread_ct * ((g_cluster_ct + (CACHELINE_INT32 - 1)) / CACHELINE_INT32) * CACHELINE)) {
	  goto glm_linear_assoc_ret_NOMEM;
	}
	if (wkspace_alloc_ui_checked(&g_indiv_to_cluster, indiv_valid_ct * sizeof(int32_t))) {
	  goto glm_linear_assoc_ret_NOMEM;
	}
	fill_unfiltered_indiv_to_cluster(indiv_valid_ct, g_cluster_ct, g_cluster_map, g_cluster_starts, g_indiv_to_cluster);
      }
    }
  }
  if (do_perms) {
    if (wkspace_init_sfmtp(max_thread_ct)) {
      goto glm_linear_assoc_ret_NOMEM;
    }
  }
  if (wkspace_alloc_d_checked(&g_pheno_d2, indiv_valid_ct * sizeof(double))) {
    goto glm_linear_assoc_ret_NOMEM;
  }
  g_linear_mt = (Linear_multithread*)wkspace_alloc(max_thread_ct * sizeof(Linear_multithread));
  ulii = (perm_batch_size + (BITCT - 1)) / BITCT;
  for (tidx = 0; tidx < max_thread_ct; tidx++) {
    if (wkspace_alloc_d_checked(&(g_linear_mt[tidx].cur_covars_cov_major), param_ct_max * indiv_valid_ct * sizeof(double)) ||
        wkspace_alloc_d_checked(&(g_linear_mt[tidx].cur_covars_indiv_major), param_ct_max * indiv_valid_ct * sizeof(double)) ||
        wkspace_alloc_ul_checked(&(g_linear_mt[tidx].perm_fails), ulii * sizeof(intptr_t)) ||
        wkspace_alloc_d_checked(&(g_linear_mt[tidx].indiv_1d_buf), indiv_valid_ct * sizeof(double)) ||
        wkspace_alloc_d_checked(&(g_linear_mt[tidx].param_2d_buf), param_ct_max * param_ct_max * sizeof(double)) ||
        wkspace_alloc_d_checked(&(g_linear_mt[tidx].param_2d_buf2), param_ct_max * param_ct_max * sizeof(double)) ||
        wkspace_alloc_d_checked(&(g_linear_mt[tidx].regression_results), perm_batch_size * param_ctx_max_m1 * sizeof(double))) {
      goto glm_linear_assoc_ret_NOMEM;
    }

    g_linear_mt[tidx].mi_buf = (MATRIX_INVERT_BUF1_TYPE*)wkspace_alloc(param_ct_max * sizeof(MATRIX_INVERT_BUF1_TYPE));
    if (!(g_linear_mt[tidx].mi_buf)) {
      goto glm_linear_assoc_ret_NOMEM;
    }
    if (cluster_ct1) {
      ulii = MAXV(cluster_ct1 + 1, param_ct_max);
      if (wkspace_alloc_ui_checked(&(g_linear_mt[tidx].cur_indiv_to_cluster1_buf), indiv_valid_ct * sizeof(int32_t)) ||
          wkspace_alloc_d_checked(&(g_linear_mt[tidx].cluster_param_buf), ulii * param_ct_max * sizeof(double)) ||
          wkspace_alloc_d_checked(&(g_linear_mt[tidx].cluster_param_buf2), (cluster_ct1 + 1) * param_ct_max * sizeof(double))) {
	goto glm_linear_assoc_ret_NOMEM;
      }
    }
    if (constraint_ct_max) {
      if (wkspace_alloc_d_checked(&(g_linear_mt[tidx].df_df_buf), constraint_ct_max * constraint_ct_max * sizeof(double)) ||
	  wkspace_alloc_d_checked(&(g_linear_mt[tidx].df_buf), constraint_ct_max * sizeof(double))) {
	goto glm_linear_assoc_ret_NOMEM;
      }
    }
    if (wkspace_alloc_d_checked(&(g_linear_mt[tidx].param_df_buf), constraint_ct_max * param_ct_max * sizeof(double)) ||
	wkspace_alloc_d_checked(&(g_linear_mt[tidx].param_df_buf2), constraint_ct_max * param_ct_max * sizeof(double)) ||
	wkspace_alloc_d_checked(&(g_linear_mt[tidx].dgels_a), param_ct_max * indiv_valid_ct * sizeof(double)) ||
	wkspace_alloc_d_checked(&(g_linear_mt[tidx].dgels_b), perm_batch_size * indiv_valid_ct * sizeof(double))) {
      goto glm_linear_assoc_ret_NOMEM;
    }
    if (!tidx) {
      dgels_m = (int32_t)((uint32_t)indiv_valid_ct);
      dgels_n = (int32_t)((uint32_t)param_ct_max);
      dgels_nrhs = perm_batch_size;
      dgels_ldb = dgels_m;
      g_dgels_lwork = -1;
      // no other parameters are needed for workspace query
      dgels_(&dgels_trans, &dgels_m, &dgels_n, &dgels_nrhs, g_linear_mt[0].dgels_a, &dgels_m, g_linear_mt[0].dgels_b, &dgels_ldb, &dxx, &g_dgels_lwork, &dgels_info);
      if (dxx > 2147483647.0) {
	logprint("Error: Multiple linear regression problem too large for current LAPACK version.\n");
	retval = RET_CALC_NOT_YET_SUPPORTED;
	goto glm_linear_assoc_ret_1;
      }
      g_dgels_lwork = (int32_t)dxx;
    }
    if (wkspace_alloc_d_checked(&(g_linear_mt[tidx].dgels_work), g_dgels_lwork * sizeof(double))) {
      goto glm_linear_assoc_ret_NOMEM;
    }
  }

  dptr = g_pheno_d2;
  g_pheno_sum = 0;
  g_pheno_ssq = 0;
  indiv_uidx = 0;
  indiv_idx = 0;
  do {
    indiv_uidx = next_set_ul_unsafe(load_mask, indiv_uidx);
    indiv_uidx_stop = next_unset_ul(load_mask, indiv_uidx, unfiltered_indiv_ct);
    indiv_idx += indiv_uidx_stop - indiv_uidx;
    dptr2 = &(pheno_d[indiv_uidx]);
    indiv_uidx = indiv_uidx_stop;
    dptr3 = &(pheno_d[indiv_uidx_stop]);
    do {
      dxx = *dptr2++;
      *dptr++ = dxx;
      g_pheno_sum += dxx;
      g_pheno_ssq += dxx * dxx;
    } while (dptr2 < dptr3);
  } while (indiv_idx < indiv_valid_ct);
  if (g_pheno_ssq * ((double)((intptr_t)indiv_valid_ct)) == g_pheno_sum * g_pheno_sum) {
    goto glm_linear_assoc_ret_PHENO_CONSTANT;
  }
  if (standard_beta) {
    dxx = g_pheno_sum / ((double)((intptr_t)indiv_valid_ct));
    dyy = sqrt(((double)((intptr_t)(indiv_valid_ct - 1))) / (g_pheno_ssq - g_pheno_sum * dxx));
    dptr = g_pheno_d2;
    for (indiv_idx = 0; indiv_idx < indiv_valid_ct; indiv_idx++) {
      *dptr = ((*dptr) - dxx) * dyy;
      dptr++;
    }
    g_pheno_sum = 0;
    g_pheno_ssq = (double)((intptr_t)(indiv_valid_ct - 1));
  }
  if (do_perms) {
    if (wkspace_alloc_ui_checked(&g_precomputed_mods, (indiv_valid_ct - 1) * sizeof(int32_t)) ||
	wkspace_alloc_d_checked(&g_perm_pmajor, perm_batch_size * indiv_valid_ct * sizeof(double))) {
      goto glm_linear_assoc_ret_NOMEM;
    }
    precompute_mods(indiv_valid_ct, g_precomputed_mods);
  }

  outname_end2 = memcpyb(outname_end, ".assoc.linear", 14);
  if (fopen_checked(&outfile, outname, "w")) {
    goto glm_linear_assoc_ret_OPEN_FAIL;
  }
  LOGPRINTFWW5("Writing linear model association results to %s ... ", outname);
  fflush(stdout);
  sprintf(tbuf, " CHR %%%us         BP   A1       TEST    NMISS       BETA ", plink_maxsnp);
  fprintf(outfile, tbuf, "SNP");
  if (display_ci) {
    uii = (uint32_t)((int32_t)(ci_size * 100));
    if (uii >= 10) {
      fprintf(outfile, "      SE      L%u      U%u ", uii, uii);
    } else {
      fprintf(outfile, "      SE       L%u       U%u ", uii, uii);
    }
  }
  fputs("        STAT            P \n", outfile);
  loop_end = marker_initial_ct / 100;
  marker_unstopped_ct = marker_initial_ct;
  g_adaptive_ci_zt = ltqnorm(1 - apip->beta / (2.0 * ((int32_t)marker_initial_ct)));

  // ----- begin main loop -----
 glm_linear_assoc_more_perms:
  if (do_perms) {
    if (perm_adapt) {
      if (perm_pass_idx) {
	while (g_first_adapt_check <= g_perms_done) {
	  g_first_adapt_check += (int32_t)(apip->init_interval + ((int32_t)g_first_adapt_check) * apip->interval_slope);
	}
      } else {
	if (apip->min < apip->init_interval) {
	  g_first_adapt_check = (int32_t)(apip->init_interval);
	} else {
	  g_first_adapt_check = apip->min;
	}
	g_adaptive_intercept = apip->init_interval;
        g_adaptive_slope = apip->interval_slope;
      }
    }
    g_perm_vec_ct = perm_batch_size;
    if (perm_adapt && (g_perms_done < perm_batch_size)) {
      // special case: split first batch to reduce adaptive overshoot
      ulii = perm_batch_size;
      uljj = (intptr_t)(apip->init_interval);
      uljj = MAXV(uljj, apip->min);
      uljj *= 2;
      uljj = MAXV(64, uljj);
      while (ulii >= (uljj << perm_pass_idx)) {
	ulii >>= 1;
      }
      g_perm_vec_ct = ulii - g_perms_done;
    }
    if (g_perm_vec_ct > perms_total - g_perms_done) {
      g_perm_vec_ct = perms_total - g_perms_done;
    }
    ulii = 0;
    if (g_perm_vec_ct >= CACHELINE_INT32 * max_thread_ct) {
      g_assoc_thread_ct = max_thread_ct;
    } else {
      g_assoc_thread_ct = g_perm_vec_ct / CACHELINE_INT32;
      if (!g_assoc_thread_ct) {
	g_assoc_thread_ct = 1;
      }
    }
    if (!g_cluster_ct) {
      if (spawn_threads(threads, &linear_gen_perms_thread, g_assoc_thread_ct)) {
	goto glm_linear_assoc_ret_THREAD_CREATE_FAIL;
      }
      linear_gen_perms_thread((void*)ulii);
    } else {
      if (spawn_threads(threads, &linear_gen_cluster_perms_thread, g_assoc_thread_ct)) {
	goto glm_linear_assoc_ret_THREAD_CREATE_FAIL;
      }
      linear_gen_cluster_perms_thread((void*)ulii);
    }
    join_threads(threads, g_assoc_thread_ct);
  }
  chrom_fo_idx = 0xffffffffU;
  marker_uidx = next_unset_unsafe(marker_exclude, 0);
  if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
    goto glm_linear_assoc_ret_READ_FAIL;
  }
  if (!perm_pass_idx) {
    fputs("0%", stdout);
    fflush(stdout);
  }
  marker_idx = 0;
  marker_idx2 = 0;
  chrom_end = 0;
  do {
    if (marker_uidx >= chrom_end) {
      // exploit overflow
      do {
        chrom_fo_idx++;
        refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &g_is_x, &g_is_y, &uii, &g_is_haploid);
      } while ((!glm_xchr_model) && (g_is_haploid || uii));
      uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
      wptr_start = width_force(4, writebuf, chrom_name_write(writebuf, chrom_info_ptr, uii, zero_extra_chroms));
      *wptr_start++ = ' ';
      fill_double_zero(constraints_con_major, constraint_ct_max * param_ct_max);
      g_male_x_01 = 0;
      if (g_is_x) {
        cur_param_ct = param_ct_max;
	cur_constraint_ct = constraint_ct_max;
	cur_param_names = param_names;
	if (glm_xchr_model != 2) {
	  g_male_x_01 = 1;
	}
      } else if ((!g_is_haploid) || (!genotypic_or_hethom)) {
	cur_param_ct = np_base + np_diploid;
	if (constraint_ct_max) {
          cur_constraint_ct = popcount_bit_idx(g_joint_test_params, 0, constraint_ct_max - np_sex);
	} else {
	  cur_constraint_ct = 0;
	}
	cur_param_names = param_names;
      } else {
	cur_param_ct = np_base;
	cur_constraint_ct = 0;
	if (constraint_ct_max) {
	  for (ulii = 0; ulii < param_ctl_max; ulii++) {
	    uljj = g_joint_test_params[ulii] & haploid_params[ulii];
	    while (uljj) {
	      uii = CTZLU(uljj) + ulii * BITCT;
              constraints_con_major[cur_constraint_ct * cur_param_ct + uii + 1] = 1;
	      cur_constraint_ct++;
	      uljj &= uljj - 1;
	    }
	  }
	}
	cur_param_names = haploid_param_names;
      }
      g_cur_param_ct = cur_param_ct;
      g_cur_constraint_ct = cur_constraint_ct;
      if (!hide_covar) {
	param_idx_end = cur_param_ct;
      }
      cur_param_ctx = cur_param_ct;
      g_include_sex = sex_covar_everywhere || (g_is_x && np_sex);
      if (cur_constraint_ct) {
	cur_param_ctx++;
	if (g_is_x || (!g_is_haploid)) {
	  ulii = 0;
	  for (constraint_idx = 0; constraint_idx < cur_constraint_ct; ulii++, constraint_idx++) {
            next_set_ul_unsafe_ck(g_joint_test_params, &ulii);
	    constraints_con_major[constraint_idx * cur_param_ct + ulii + 1] = 1;
	  }
	}
        wptr = uint32_write(&(param_names[param_ct_max * max_param_name_len + 5]), cur_constraint_ct);
	memcpy(wptr, "DF", 3);
      }
    }
    block_size = 0;
    block_end = marker_unstopped_ct - marker_idx;
    if (block_end > GLM_BLOCKSIZE) {
      block_end = GLM_BLOCKSIZE;
    }
    do {
      if (g_perm_adapt_stop[marker_idx2]) {
	do {
	  marker_uidx++;
	  next_unset_unsafe_ck(marker_exclude, &marker_uidx);
	  marker_idx2++;
	} while ((marker_uidx < chrom_end) && g_perm_adapt_stop[marker_idx2]);
	if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	  goto glm_linear_assoc_ret_READ_FAIL;
	}
	if (marker_uidx >= chrom_end) {
	  break;
	}
      }
      loadbuf_ptr = &(g_loadbuf[block_size * indiv_valid_ctv2]);
      if (load_and_collapse_incl(bedfile, loadbuf_raw, unfiltered_indiv_ct, loadbuf_ptr, indiv_valid_ct, load_mask, IS_SET(marker_reverse, marker_uidx))) {
	goto glm_linear_assoc_ret_READ_FAIL;
      }
      if (hh_exists) {
	haploid_fix(hh_exists, indiv_include2, indiv_male_include2, indiv_valid_ct, g_is_x, g_is_y, (unsigned char*)loadbuf_ptr);
      }
      g_adapt_m_table[block_size] = marker_idx2++;
      mu_table[block_size++] = marker_uidx;
      if (marker_idx + block_size == marker_unstopped_ct) {
	break;
      }
      marker_uidx++;
      if (IS_SET(marker_exclude, marker_uidx)) {
	marker_uidx = next_unset_unsafe(marker_exclude, marker_uidx);
	if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	  goto glm_linear_assoc_ret_READ_FAIL;
	}
      }
    } while ((block_size < block_end) && (marker_uidx < chrom_end));
    if (!block_size) {
      continue;
    }
    g_block_diff = block_size;
    if (!perm_pass_idx) {
#ifndef NOLAPACK
      dgels_nrhs = 1;
#endif
      for (marker_bidx = 0; marker_bidx < block_size; marker_bidx++) {
	// 1. fill design matrix and make transposed copy; fill g_nm_cts[] in
	//    the process
	// 2. 'standard-beta' adjustment if necessary
	// 3. linear/logistic regression, set g_perm_adapt_stop byte on failure
	// 4. write basic report to disk
	marker_uidx2 = mu_table[marker_bidx];
        loadbuf_ptr = &(g_loadbuf[marker_bidx * indiv_valid_ctv2]);
	marker_idx3 = g_adapt_m_table[marker_bidx];
	if (marker_idx_to_uidx) {
	  marker_idx_to_uidx[marker_idx3] = marker_uidx2;
	}
        cur_missing_ct = glm_fill_design(loadbuf_ptr, g_fixed_covars_cov_major, indiv_valid_ct, g_indiv_to_cluster1, cur_param_ct, glm_modifier & (GLM_HETHOM | GLM_DOMINANT | GLM_RECESSIVE), glm_xchr_model, condition_list_start_idx, interaction_start_idx, sex_start_idx, active_params, haploid_params, g_include_sex, g_male_x_01, sex_male_collapsed, g_is_haploid && (!g_is_x), g_linear_mt[0].cur_covars_cov_major, g_linear_mt[0].cur_covars_indiv_major, g_linear_mt[0].cur_indiv_to_cluster1_buf, &cur_indiv_to_cluster1, standard_beta);
	cur_indiv_valid_ct = indiv_valid_ct - cur_missing_ct;
	g_nm_cts[marker_idx3] = cur_indiv_valid_ct;
	if ((cur_indiv_valid_ct > cur_param_ct) && (!glm_check_vif(glm_vif_thresh, cur_param_ct, cur_indiv_valid_ct, g_linear_mt[0].cur_covars_cov_major, g_linear_mt[0].param_2d_buf, g_linear_mt[0].mi_buf, g_linear_mt[0].param_2d_buf2))) {
	  regression_fail = 0;
	  memcpy(g_linear_mt[0].dgels_a, g_linear_mt[0].cur_covars_cov_major, cur_param_ct * cur_indiv_valid_ct * sizeof(double));
	  copy_when_nonmissing(loadbuf_ptr, (char*)g_pheno_d2, sizeof(double), indiv_valid_ct, cur_missing_ct, (char*)(g_linear_mt[0].dgels_b));
	  if (standard_beta && cur_missing_ct) {
	    dxx = g_pheno_sum;
	    dyy = g_pheno_ssq;
	    ulptr = loadbuf_ptr;
	    for (indiv_idx = 0; indiv_idx < indiv_valid_ct; indiv_idx += BITCT2) {
	      cur_word = *ulptr++;
	      cur_word = cur_word & (~(cur_word >> 1)) & FIVEMASK;
	      while (cur_word) {
		dzz = g_pheno_d2[indiv_idx + (CTZLU(cur_word) / 2)];
		dxx -= dzz;
		dyy -= dzz * dzz;
		cur_word &= cur_word - 1;
	      }
	    }
	    dzz = dxx / ((double)((intptr_t)cur_indiv_valid_ct));
	    dyy = sqrt(((double)((intptr_t)(cur_indiv_valid_ct - 1))) / (dyy - dxx * dzz));
	    dptr = g_linear_mt[0].dgels_b;
	    for (indiv_idx = 0; indiv_idx < cur_indiv_valid_ct; indiv_idx++) {
	      *dptr = ((*dptr) - dzz) * dyy;
	      dptr++;
	    }
	  }
	  dgels_m = (int32_t)((uint32_t)cur_indiv_valid_ct);
	  dgels_n = (int32_t)((uint32_t)cur_param_ct);
	  dgels_ldb = dgels_m;

	  dgels_(&dgels_trans, &dgels_m, &dgels_n, &dgels_nrhs, g_linear_mt[0].dgels_a, &dgels_m, g_linear_mt[0].dgels_b, &dgels_ldb, g_linear_mt[0].dgels_work, &g_dgels_lwork, &dgels_info);
	  if (glm_linear_robust_cluster_covar(1, cur_param_ct, cur_indiv_valid_ct, cur_missing_ct, loadbuf_ptr, standard_beta, g_pheno_sum, g_pheno_ssq, g_linear_mt[0].cur_covars_cov_major, g_linear_mt[0].cur_covars_indiv_major, g_pheno_d2, g_linear_mt[0].dgels_b, g_linear_mt[0].param_2d_buf, g_linear_mt[0].mi_buf, g_linear_mt[0].param_2d_buf2, cluster_ct1, cur_indiv_to_cluster1, g_linear_mt[0].cluster_param_buf, g_linear_mt[0].cluster_param_buf2, g_linear_mt[0].indiv_1d_buf, g_linear_mt[0].regression_results, cur_constraint_ct, constraints_con_major, g_linear_mt[0].param_df_buf, g_linear_mt[0].param_df_buf2, g_linear_mt[0].df_df_buf, g_linear_mt[0].df_buf, &perm_fail_ct, g_linear_mt[0].perm_fails) || perm_fail_ct) {
	    regression_fail = 1;
	  }
	} else {
	  regression_fail = 1;
	}
	wptr_start2 = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx2 * max_marker_id_len]), wptr_start);
	*wptr_start2++ = ' ';
	wptr_start2 = uint32_writew10(wptr_start2, marker_pos[marker_uidx2]);
	*wptr_start2++ = ' ';
	wptr_start2 = fw_strcpy(4, marker_allele_ptrs[marker_uidx2 * 2], wptr_start2);
	*wptr_start2++ = ' ';
	orig_stats_ptr = &(g_orig_stats[marker_idx3]);
	if (!regression_fail) {
	  for (param_idx = 1; param_idx < cur_param_ct; param_idx++) {
	    dxx = g_linear_mt[0].dgels_b[param_idx]; // coef[p]
	    se = sqrt(g_linear_mt[0].regression_results[param_idx - 1]);
	    zval = dxx / se;
	    if (param_idx == 1) {
	      if (mtest_adjust) {
		g_orig_chisq[marker_idx3] = fabs(zval);
		tcnt[marker_idx3] = cur_indiv_valid_ct;
	      }
	      if (mperm_save_all) {
		fprintf(outfile_msa, " %g", fabs(zval));
	      }
	      if (!cur_constraint_ct) {
		*orig_stats_ptr = fabs(zval);
	      }
	    }
	    pval = calc_tprob(zval, cur_indiv_valid_ct - cur_param_ct);
	    if ((param_idx < param_idx_end) && (pval <= pfilter)) {
	      wptr = fw_strcpy(10, &(cur_param_names[param_idx * max_param_name_len]), wptr_start2);
	      *wptr++ = ' ';
	      wptr = uint32_writew8x(wptr, (uint32_t)cur_indiv_valid_ct, ' ');
	      wptr = double_g_writewx4x(wptr, dxx, 10, ' ');
	      if (display_ci) {
		dyy = ci_zt * se;
		wptr = double_g_writewx4x(wptr, se, 8, ' ');
		wptr = double_g_writewx4x(wptr, dxx - dyy, 8, ' ');
		wptr = double_g_writewx4x(wptr, dxx + dyy, 8, ' ');
	      }
	      wptr = double_g_writewx4x(wptr, zval, 12, ' ');
	      wptr = double_g_writewx4x(wptr, pval, 12, '\n');
	      if (fwrite_checked(writebuf, wptr - writebuf, outfile)) {
		goto glm_linear_assoc_ret_WRITE_FAIL;
	      }
	    }
	  }
	  if (linear_intercept) {
	    dxx = g_linear_mt[0].dgels_b[0];
	    wptr = memcpya(wptr_start2, " INTERCEPT ", 11);
	    wptr = uint32_writew8x(wptr, (uint32_t)cur_indiv_valid_ct, ' ');
	    wptr = double_g_writewx4x(wptr, dxx, 10, ' ');
	    if (display_ci) {
	      // okay, this should be made more maintainable...
	      se = sqrt(g_linear_mt[0].param_2d_buf2[0]);
	      zval = dxx / se;
	      dyy = ci_zt * se;
	      wptr = double_g_writewx4x(wptr, se, 8, ' ');
	      wptr = double_g_writewx4x(wptr, dxx - dyy, 8, ' ');
	      wptr = double_g_writewx4x(wptr, dxx + dyy, 8, ' ');
	    }
	    wptr = memcpya(wptr, "          NA           NA\n", 26);
	    if (fwrite_checked(writebuf, wptr - writebuf, outfile)) {
	      goto glm_linear_assoc_ret_WRITE_FAIL;
	    }
	  }
	  if (cur_constraint_ct) {
	    dxx = g_linear_mt[0].regression_results[cur_param_ct - 1];
	    *orig_stats_ptr = dxx;
	    pval = chiprob_p(dxx, cur_constraint_ct);
	    if (pval <= pfilter) {
	      wptr = fw_strcpy(10, &(param_names[param_ct_max * max_param_name_len]), wptr_start2);
              *wptr++ = ' ';
              wptr = uint32_writew8(wptr, (uint32_t)cur_indiv_valid_ct);
              wptr = memcpya(wptr, "         NA ", 12);
              if (display_ci) {
		wptr = memcpya(wptr, "      NA       NA       NA ", 27);
	      }
              wptr = double_g_writewx4x(wptr, dxx, 12, ' ');
              wptr = double_g_writewx4x(wptr, pval, 12, '\n');
              if (fwrite_checked(writebuf, wptr - writebuf, outfile)) {
		goto glm_linear_assoc_ret_WRITE_FAIL;
	      }
	    }
	  }
	} else {
	  g_perm_adapt_stop[marker_idx3] = 1;
	  if (mtest_adjust && (param_idx == 1)) {
	    g_orig_chisq[marker_idx3] = -9;
	    if (tcnt) {
	      tcnt[marker_idx3] = 0;
	    }
	  }
	  *orig_stats_ptr = -9;
	  if (mperm_save_all) {
	    fputs(" NA", outfile_msa);
	    msa_ptr = &(g_mperm_save_all[marker_idx3 * perm_batch_size]);
	    for (perm_idx = 0; perm_idx < perm_batch_size; perm_idx++) {
              *msa_ptr++ = -9;
	    }
	  }
	  for (param_idx = 1; param_idx < cur_param_ctx; param_idx++) {
	    if ((param_idx < param_idx_end) || (param_idx == cur_param_ct)) {
	      if (!(param_idx == cur_param_ct)) {
	        wptr = fw_strcpy(10, &(cur_param_names[param_idx * max_param_name_len]), wptr_start2);
	      } else {
		wptr = fw_strcpy(10, &(param_names[param_ct_max * max_param_name_len]), wptr_start2);
	      }
	      *wptr++ = ' ';
	      wptr = uint32_writew8(wptr, (uint32_t)cur_indiv_valid_ct);
	      wptr = memcpya(wptr, "         NA ", 12);
	      if (display_ci) {
		wptr = memcpya(wptr, "      NA       NA       NA ", 27);
	      }
	      wptr = memcpya(wptr, "          NA           NA\n", 26);
	      if (fwrite_checked(writebuf, wptr - writebuf, outfile)) {
		goto glm_linear_assoc_ret_WRITE_FAIL;
	      }
	    }
	  }
	}
      }
    }
    if (do_perms) {
      // split loaded markers evenly between threads.  each thread does the
      // following for every marker in its set:
      // 1. fill design matrix and make transposed copy, do this more
      //    efficiently if g_nm_cts[] indicates there are no missing genotypes
      // 2. 'standard-beta' adjustment if necessary (these first two steps
      //    should go into another function since they're mostly duplicated)
      // 3. linear/logistic regression on permutations, handle adaptive logic,
      //    --mperm-save[-all], etc.
      g_assoc_thread_ct = block_size / CACHELINE_INT32;
      if (g_assoc_thread_ct > max_thread_ct) {
	g_assoc_thread_ct = max_thread_ct;
      } else if (!g_assoc_thread_ct) {
	g_assoc_thread_ct = 1;
      }
      ulii = 0;
      if (perm_adapt) {
	if (spawn_threads(threads, &glm_linear_adapt_thread, g_assoc_thread_ct)) {
	  goto glm_linear_assoc_ret_THREAD_CREATE_FAIL;
	}
	glm_linear_adapt_thread((void*)ulii);
	join_threads(threads, g_assoc_thread_ct);
      } else {
	if (spawn_threads(threads, &glm_linear_maxt_thread, g_assoc_thread_ct)) {
	  goto glm_linear_assoc_ret_THREAD_CREATE_FAIL;
	}
	glm_linear_maxt_thread((void*)ulii);
	join_threads(threads, g_assoc_thread_ct);
        ulii = (g_perm_vec_ct + (CACHELINE_DBL - 1)) & (~(CACHELINE_DBL - 1));
	ukk = g_perms_done + g_perm_vec_ct;
        for (uii = 0; uii < g_assoc_thread_ct; uii++) {
          dptr = &(g_maxt_thread_results[uii * ulii]);
          for (ujj = g_perms_done; ujj < ukk; ujj++) {
            dxx = *dptr++;
            if (dxx > g_maxt_extreme_stat[ujj]) {
	      g_maxt_extreme_stat[ujj] = dxx;
	    }
	  }
	}
      }
    }
    marker_idx += block_size;
    if ((!perm_pass_idx) && (marker_idx >= loop_end)) {
      if (marker_idx < marker_initial_ct) {
	if (pct >= 10) {
	  putchar('\b');
	}
        pct = (marker_idx * 100LLU) / marker_initial_ct;
        printf("\b\b%u%%", pct);
        fflush(stdout);
        loop_end = (((uint64_t)pct + 1LLU) * marker_initial_ct) / 100;
      }
    }
  } while (marker_idx < marker_unstopped_ct);
  // if more permutations, reevaluate marker_unstopped_ct, etc.
  if (!perm_pass_idx) {
    if (pct >= 10) {
      putchar('\b');
    }
    fputs("\b\b", stdout);
    logprint("done.\n");
    if (fclose_null(&outfile)) {
      goto glm_linear_assoc_ret_WRITE_FAIL;
    }
    if (mtest_adjust) {
      retval = multcomp(outname, outname_end, marker_idx_to_uidx, marker_initial_ct, marker_ids, max_marker_id_len, plink_maxsnp, zero_extra_chroms, chrom_info_ptr, g_orig_chisq, pfilter, mtest_adjust, adjust_lambda, 1, tcnt);
      if (retval) {
	goto glm_linear_assoc_ret_1;
      }
    }
    if (mperm_save_all) {
      if (putc_checked('\n', outfile_msa)) {
	goto glm_linear_assoc_ret_WRITE_FAIL;
      }
    }
  }
  if (do_perms) {
    if (mperm_save_all) {
      if (perm_pass_idx) {
	putchar(' ');
      }
      fputs("[dumping stats]", stdout);
      fflush(stdout);
      ulii = g_perm_vec_ct;
      ujj = 1 + g_perms_done;
      wptr = tbuf;
      cptr = &(tbuf[MAXLINELEN]);
      for (uii = 0; uii < ulii; uii++) {
	wptr = uint32_write(wptr, uii + ujj);
	dptr = &(g_mperm_save_all[uii]);
	for (ukk = 0; ukk < marker_ct; ukk++) {
	  *wptr++ = ' ';
	  dxx = dptr[ukk * ulii];
	  if (dxx >= 0) {
	    wptr = double_g_write(wptr, dxx);
	  } else {
	    wptr = memcpya(wptr, "NA", 2);
	  }
	  if (wptr >= cptr) {
	    if (fwrite_checked(tbuf, (uintptr_t)(wptr - tbuf), outfile_msa)) {
	      goto glm_linear_assoc_ret_WRITE_FAIL;
	    }
	    wptr = tbuf;
	  }
	}
	*wptr++ = '\n';
      }
      if (fwrite_checked(tbuf, (uintptr_t)(wptr - tbuf), outfile_msa)) {
	goto glm_linear_assoc_ret_WRITE_FAIL;
      }
      fputs("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b               ", stdout);
    }
    g_perms_done += g_perm_vec_ct;
    if (g_perms_done < perms_total) {
      if (perm_adapt || (!perm_pass_idx)) {
        marker_unstopped_ct = marker_initial_ct - popcount_longs((uintptr_t*)g_perm_adapt_stop, (marker_initial_ct + sizeof(intptr_t) - 1) / sizeof(intptr_t));
        if (!marker_unstopped_ct) {
          goto glm_linear_assoc_perm_count;
	}
      }
      printf("\r%u permutation%s complete.", g_perms_done, (g_perms_done != 1)? "s" : "");
      fflush(stdout);
      perm_pass_idx++;
      goto glm_linear_assoc_more_perms;
    }
  glm_linear_assoc_perm_count:
    if (perm_adapt) {
      g_perms_done = 0;
      for (uii = 0; uii < marker_initial_ct; uii++) {
        if (g_perm_attempt_ct[uii] > g_perms_done) {
	  g_perms_done = g_perm_attempt_ct[uii];
	  if (g_perms_done == perms_total) {
	    break;
	  }
	}
      }
    }
    putchar('\r');
    LOGPRINTF("%u %s permutation%s complete.\n", g_perms_done, perm_maxt? "max(T)" : "(adaptive)", (g_perms_done != 1)? "s" : "");

    if (perm_adapt) {
      memcpy(outname_end2, ".perm", 6);
    } else {
      if (mperm_save & MPERM_DUMP_BEST) {
	memcpy(outname_end, ".mperm.dump.best", 17);
	LOGPRINTFWW("Dumping best permutation %s to %s .\n", (!constraint_ct_max)? "absolute t-stats" : "chi-square values", outname);
	if (fopen_checked(&outfile, outname, "w")) {
	  goto glm_linear_assoc_ret_OPEN_FAIL;
	}
	dxx = 0;
	for (marker_idx = 0; marker_idx < marker_initial_ct; marker_idx++) {
          if (g_orig_stats[marker_idx] > dxx) {
	    dxx = g_orig_stats[marker_idx];
	  }
	}
	memcpy(tbuf, "0 ", 2);
	wptr = double_g_writex(&(tbuf[2]), dxx, '\n');
        if (fwrite_checked(tbuf, (uintptr_t)(wptr - tbuf), outfile)) {
	  goto glm_linear_assoc_ret_WRITE_FAIL;
	}
        for (uii = 0; uii < perms_total; uii++) {
          wptr = uint32_writex(tbuf, uii + 1, ' ');
          wptr = double_g_writex(wptr, g_maxt_extreme_stat[uii], '\n');
          if (fwrite_checked(tbuf, (uintptr_t)(wptr - tbuf), outfile)) {
	    goto glm_linear_assoc_ret_WRITE_FAIL;
	  }
	}
        if (fclose_null(&outfile)) {
	  goto glm_linear_assoc_ret_WRITE_FAIL;
	}
	memcpy(outname_end, ".assoc.linear", 13);
      }
      memcpy(outname_end2, ".mperm", 7);
    }
    if (fopen_checked(&outfile, outname, "w")) {
      goto glm_linear_assoc_ret_OPEN_FAIL;
    }
    if (perm_adapt) {
      sprintf(tbuf, " CHR %%%us         EMP1           NP \n", plink_maxsnp);
    } else {
      sprintf(tbuf, " CHR %%%us         EMP1         EMP2 \n", plink_maxsnp);
#ifdef __cplusplus
      std::sort(g_maxt_extreme_stat, &(g_maxt_extreme_stat[perms_total]));
#else
      qsort(g_maxt_extreme_stat, perms_total, sizeof(double), double_cmp);
#endif
      // (debugging)
      // printf("extreme stats: %g %g %g\n", g_maxt_extreme_stat[0], g_maxt_extreme_stat[(perms_total - 1) / 2], g_maxt_extreme_stat[perms_total - 1]);

      for (marker_idx = 0; marker_idx < marker_initial_ct; marker_idx++) {
	g_perm_attempt_ct[marker_idx] = perms_total - g_perm_attempt_ct[marker_idx];
      }
    }
    fprintf(outfile, tbuf, "SNP");
    chrom_fo_idx = 0xffffffffU;
    marker_uidx = next_unset_unsafe(marker_exclude, 0);
    marker_idx = 0;
    while (1) {
      do {
	chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[(++chrom_fo_idx) + 1];
      } while (marker_uidx >= chrom_end);
      uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
      wptr_start = width_force(4, tbuf, chrom_name_write(tbuf, chrom_info_ptr, uii, zero_extra_chroms));
      *wptr_start++ = ' ';
      wptr_start[plink_maxsnp] = ' ';
      for (; marker_uidx < chrom_end;) {
	pval = ((double)(g_perm_2success_ct[marker_idx] + 2)) / ((double)(2 * (g_perm_attempt_ct[marker_idx] + 1)));
        if (pval <= pfilter) {
          fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), wptr_start);
          wptr = &(wptr_start[1 + plink_maxsnp]);
	  if (g_orig_stats[marker_idx] == -9) {
            wptr = memcpya(wptr, "          NA           NA", 25);
	  } else {
	    if (!perm_count) {
	      wptr = double_g_writewx4x(wptr, pval, 12, ' ');
	    } else {
	      wptr = double_g_writewx4x(wptr, ((double)g_perm_2success_ct[marker_idx]) / 2.0, 12, ' ');
	    }
	    if (perm_adapt) {
	      wptr = memseta(wptr, 32, 2);
	      wptr = uint32_writew10(wptr, g_perm_attempt_ct[marker_idx]);
	    } else {
	      dzz = (int32_t)(perms_total - doublearr_greater_than(g_maxt_extreme_stat, perms_total, g_orig_stats[marker_idx] - EPSILON) + 1);
              if (!perm_count) {
		wptr = double_g_writewx4(wptr, dzz / ((double)((int32_t)perms_total + 1)), 12);
	      } else {
                wptr = double_g_writewx4(wptr, dzz, 12);
	      }
	    }
	  }
	  wptr = memcpya(wptr, " \n", 2);
	  if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	    goto glm_linear_assoc_ret_WRITE_FAIL;
	  }
	}
	if (++marker_idx == marker_ct) {
	  goto glm_linear_assoc_loop_end;
	}
        marker_uidx++;
        next_unset_unsafe_ck(marker_exclude, &marker_uidx);
      }
    }
  glm_linear_assoc_loop_end:
    if (fclose_null(&outfile)) {
      goto glm_linear_assoc_ret_WRITE_FAIL;
    }
    LOGPRINTFWW("Permutation test report written to %s .\n", outname);
  }

  while (0) {
  glm_linear_assoc_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  glm_linear_assoc_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  glm_linear_assoc_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  glm_linear_assoc_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  glm_linear_assoc_ret_NO_PERMUTATION_CLUSTERS:
    logprint("Error: No size 2+ clusters for permutation test.\n");
    retval = RET_INVALID_CMDLINE;
    break;
  glm_linear_assoc_ret_PHENO_CONSTANT:
    logprint("Warning: Skipping --linear since phenotype is constant.\n");
    break;
  glm_linear_assoc_ret_THREAD_CREATE_FAIL:
    logprint(errstr_thread_create);
    retval = RET_THREAD_CREATE_FAIL;
    break;
  }
 glm_linear_assoc_ret_1:
  wkspace_reset(wkspace_mark);
  fclose_cond(outfile);
  fclose_cond(outfile_msa);
  free_cond(condition_uidxs);
  return retval;
}
#endif

int32_t glm_logistic_assoc(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uint32_t glm_modifier, double glm_vif_thresh, uint32_t glm_xchr_model, uint32_t glm_mperm_val, Range_list* parameters_range_list_ptr, Range_list* tests_range_list_ptr, double ci_size, double ci_zt, double pfilter, uint32_t mtest_adjust, double adjust_lambda, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, uint32_t* marker_pos, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, uintptr_t* marker_reverse, uint32_t zero_extra_chroms, char* condition_mname, char* condition_fname, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_indiv_ct, uintptr_t indiv_ct, uintptr_t* indiv_exclude, uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, Aperm_info* apip, uint32_t mperm_save, uint32_t pheno_nm_ct, uintptr_t* pheno_nm, uintptr_t* pheno_c, uintptr_t covar_ct, char* covar_names, uintptr_t max_covar_name_len, uintptr_t* covar_nm, double* covar_d, uintptr_t* sex_nm, uintptr_t* sex_male, uint32_t hh_exists, uint32_t perm_batch_size, Set_info* sip) {
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  FILE* outfile = NULL;
  FILE* outfile_msa = NULL;
  uintptr_t indiv_uidx = 0;
  uintptr_t cur_constraint_ct = 0;
  uintptr_t cur_param_ct = 0;
  uintptr_t cur_param_ctx = 0;
  uint32_t perms_total = 0;
  uint32_t cluster_ct1 = 0;
  uint32_t perm_adapt = glm_modifier & GLM_PERM;
  uint32_t perm_maxt = glm_modifier & GLM_MPERM;
  uint32_t covar_interactions = (glm_modifier / GLM_INTERACTION) & 1;
  uint32_t hethom = glm_modifier & GLM_HETHOM;
  uint32_t genotypic_or_hethom = (glm_modifier & (GLM_GENOTYPIC | GLM_HETHOM))? 1 : 0;
  uint32_t do_perms = perm_adapt | perm_maxt;
  uint32_t perm_count = glm_modifier & GLM_PERM_COUNT;
  uint32_t hide_covar = glm_modifier & GLM_HIDE_COVAR;
  uint32_t report_odds = !(glm_modifier & GLM_BETA);
  uint32_t mperm_save_all = mperm_save & MPERM_DUMP_ALL;
  uint32_t fill_orig_chiabs = do_perms || mtest_adjust;
  uint32_t display_ci = (ci_size > 0);
  uint32_t perm_pass_idx = 0;
  uint32_t pct = 0;
  uint32_t max_thread_ct = g_thread_ct;
  int32_t retval = 0;
  double* constraints_con_major = NULL;
  uintptr_t* pheno_c_collapsed = NULL;
  uint32_t* condition_uidxs = NULL;
  uint32_t* cluster_map1 = NULL;
  uint32_t* cluster_starts1 = NULL;
  uint32_t* tcnt = NULL;
  uint32_t* marker_idx_to_uidx = NULL;
  uint32_t* cur_indiv_to_cluster1 = NULL;
  char* cur_param_names = NULL;
  char* haploid_param_names = NULL;
  char* wptr_start = NULL;
  float geno_map[12];
  uint32_t mu_table[GLM_BLOCKSIZE];
  float* geno_map_ptr;
  char* writebuf;
  char* param_names;
  const char* main_effect;

  // logistic: pval = chiprob_p(g_orig_stats[i] * g_orig_stats[i], 1);
  double* orig_stats_ptr;

  char* outname_end2;
  char* wptr;
  char* wptr_start2;
  char* cptr;
  float* fptr;
  double* dptr;
  uintptr_t* loadbuf_raw;
  uintptr_t* load_mask;
  uintptr_t* sex_male_collapsed;
  uintptr_t* indiv_include2;
  uintptr_t* indiv_male_include2;
  uintptr_t* active_params;
  uintptr_t* haploid_params;
  uintptr_t* loadbuf_ptr;
  uintptr_t max_param_name_len;
  uintptr_t np_diploid;
  uintptr_t condition_ct;
  uintptr_t constraint_ct_max;
  uintptr_t np_sex;
  uintptr_t param_idx_end;
  uintptr_t indiv_valid_ct;
  uintptr_t indiv_valid_cta4;
  uintptr_t indiv_valid_ctv2;
  uintptr_t indiv_idx;
  uintptr_t param_ctx_max;
  uintptr_t param_ctl_max;
  uintptr_t param_ctx_max_m1;
  uintptr_t condition_list_start_idx;
  uintptr_t covar_start_idx;
  uintptr_t interaction_start_idx;
  uintptr_t sex_start_idx;
  uintptr_t np_base;
  uintptr_t param_ct_max;
  uintptr_t param_ct_maxa4;
  uintptr_t cur_missing_ct;
  uintptr_t cur_indiv_valid_ct;
  uintptr_t param_idx;
  uintptr_t param_idx_fixed;
  uintptr_t constraint_idx;
  uintptr_t ulii;
  uintptr_t uljj;
  double* msa_ptr;
  double se;
  double zval;
  double pval;
  double dxx;
  double dyy;
  double dzz;
  uint32_t marker_initial_ct;
  uint32_t sex_covar_everywhere;
  uint32_t variation_in_sex;
  uint32_t x_sex_interaction;
  uint32_t marker_unstopped_ct;
  uint32_t chrom_fo_idx;
  uint32_t chrom_end;
  uint32_t block_size;
  uint32_t block_end;
  uint32_t tidx;
  uint32_t perm_idx;
  uint32_t marker_uidx;
  uint32_t marker_uidx2;
  uint32_t marker_idx;
  uint32_t marker_idx2;
  uint32_t marker_idx3;
  uint32_t marker_bidx;
  uint32_t chrom_idx;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t loop_end;
  uint32_t regression_fail;

  retval = glm_common_init(bedfile, bed_offset, glm_modifier, 0, glm_xchr_model, parameters_range_list_ptr, tests_range_list_ptr, mtest_adjust, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, max_marker_allele_len, marker_reverse, condition_mname, condition_fname, chrom_info_ptr, unfiltered_indiv_ct, indiv_ct, indiv_exclude, mperm_save, pheno_nm_ct, pheno_nm, covar_ct, covar_names, max_covar_name_len, covar_nm, covar_d, sex_nm, sex_male, hh_exists, do_perms, &perm_batch_size, &max_param_name_len, &np_diploid, &condition_ct, &np_sex, &marker_initial_ct, &sex_covar_everywhere, &variation_in_sex, &x_sex_interaction, &constraint_ct_max, &param_idx_end, &loadbuf_raw, &load_mask, &sex_male_collapsed, &indiv_include2, &indiv_male_include2, &active_params, &haploid_params, &condition_uidxs, &writebuf, &indiv_valid_ct, &param_ctx_max, &param_ctl_max, &condition_list_start_idx, &covar_start_idx, &interaction_start_idx, &sex_start_idx, &np_base, &param_ct_max);
  if (retval) {
    goto glm_logistic_assoc_ret_1;
  }
  indiv_valid_cta4 = (indiv_valid_ct + 3) & (~3);
  indiv_valid_ctv2 = 2 * ((indiv_valid_ct + BITCT - 1) / BITCT);
  param_ctx_max_m1 = param_ctx_max - 1;
  param_ct_maxa4 = (param_ct_max + 3) & (~3);
  if (wkspace_alloc_d_checked(&g_orig_stats, marker_initial_ct * sizeof(double)) ||
      wkspace_alloc_c_checked(&param_names, param_ctx_max * max_param_name_len) ||
      wkspace_alloc_f_checked(&g_fixed_covars_cov_major_f, (variation_in_sex + interaction_start_idx - condition_list_start_idx) * indiv_valid_ct * sizeof(float)) ||
      wkspace_alloc_uc_checked(&g_perm_adapt_stop, marker_initial_ct) ||
      wkspace_alloc_ui_checked(&g_nm_cts, marker_initial_ct * sizeof(int32_t))) {
    goto glm_logistic_assoc_ret_NOMEM;
  }
  // use this array to track regression failures even in max(T) case
  fill_ulong_zero((uintptr_t*)g_perm_adapt_stop, (marker_initial_ct + sizeof(intptr_t) - 1) / sizeof(intptr_t));

  param_idx = 1;
  if (glm_modifier & GLM_RECESSIVE) {
    main_effect = glm_main_effects;
  } else if (glm_modifier & GLM_DOMINANT) {
    main_effect = &(glm_main_effects[4]);
  } else if (hethom) {
    main_effect = &(glm_main_effects[8]);
  } else {
    main_effect = &(glm_main_effects[12]);
  }
  if (IS_SET(active_params, 1)) {
    memcpy(&(param_names[max_param_name_len]), main_effect, 4);
    param_idx++;
  }
  if (genotypic_or_hethom && IS_SET(active_params, 2)) {
    if (hethom) {
      memcpy(&(param_names[param_idx * max_param_name_len]), "HET", 4);
    } else {
      memcpy(&(param_names[param_idx * max_param_name_len]), "DOMDEV", 7);
    }
    param_idx++;
  }
  // 0..3: diploid chromosomes, X chromosome female
  // 4..7: X chromosome male
  // 8..11: haploid
  fill_float_zero(geno_map, 12);
  geno_map[0] = 1;
  geno_map[2] = 1;
  geno_map[4] = 1;
  geno_map[8] = 1;
  if (glm_modifier & GLM_CONDITION_RECESSIVE) {
    geno_map[2] = 0;
  } else if (!(glm_modifier & GLM_CONDITION_DOMINANT)) {
    geno_map[0] = 2;
    if (glm_xchr_model == 2) {
      geno_map[4] = 2;
    }
  }
  for (param_idx_fixed = 0; param_idx_fixed < condition_ct; param_idx_fixed++) {
    marker_uidx = condition_uidxs[param_idx_fixed];
    if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
      goto glm_logistic_assoc_ret_READ_FAIL;
    }
    if (load_and_collapse_incl(bedfile, loadbuf_raw, unfiltered_indiv_ct, g_loadbuf, indiv_valid_ct, load_mask, IS_SET(marker_reverse, marker_uidx))) {
      goto glm_logistic_assoc_ret_READ_FAIL;
    }
    chrom_idx = get_marker_chrom(chrom_info_ptr, marker_uidx);
    geno_map_ptr = geno_map;
    if (IS_SET(chrom_info_ptr->haploid_mask, chrom_idx)) {
      g_is_x = ((int32_t)chrom_idx == chrom_info_ptr->x_code);
      g_is_y = ((int32_t)chrom_idx == chrom_info_ptr->y_code);
      if (hh_exists) {
	haploid_fix(hh_exists, indiv_include2, indiv_male_include2, indiv_valid_ct, g_is_x, g_is_y, (unsigned char*)g_loadbuf);
      }
      if (!g_is_x) {
	geno_map_ptr = &(geno_map[8]);
      }
    } else {
      g_is_x = 0;
    }
    if (!g_is_x) {
      glm_loadbuf_to_floats(g_loadbuf, indiv_valid_ct, &(g_fixed_covars_cov_major_f[param_idx_fixed * indiv_valid_ct]), geno_map_ptr, NULL);
    } else {
      glm_loadbuf_to_floats_x(g_loadbuf, sex_male_collapsed, indiv_valid_ct, &(g_fixed_covars_cov_major_f[param_idx_fixed * indiv_valid_ct]), geno_map_ptr, NULL);
    }
    if (is_set(active_params, param_idx_fixed + condition_list_start_idx)) {
      strcpy(&(param_names[param_idx * max_param_name_len]), &(marker_ids[marker_uidx * max_marker_id_len]));
      param_idx++;
    }
  }
  for (uii = 0; uii < covar_ct; param_idx_fixed++, uii++) {
    // indiv-major to covariate-major
    fptr = &(g_fixed_covars_cov_major_f[param_idx_fixed * indiv_valid_ct]);
    dptr = &(covar_d[uii]);
    for (indiv_uidx = 0, indiv_idx = 0; indiv_idx < indiv_ct; indiv_uidx++, indiv_idx++) {
      next_unset_ul_unsafe_ck(indiv_exclude, &indiv_uidx);
      if (IS_SET(load_mask, indiv_uidx)) {
        *fptr++ = (float)dptr[indiv_idx * covar_ct];
      }
    }
  }
  for (uii = 0; uii < covar_ct; uii++) {
    if (is_set(active_params, uii + covar_start_idx)) {
      strcpy(&(param_names[param_idx * max_param_name_len]), &(covar_names[uii * max_covar_name_len]));
      param_idx++;
    }
  }
  if (covar_interactions) {
    ujj = interaction_start_idx;
    for (uii = 0; uii < condition_ct; uii++) {
      if (is_set(active_params, ujj++)) {
        wptr = memcpyl3a(&(param_names[param_idx * max_param_name_len]), main_effect);
	wptr = memcpya(wptr, "xCSNP", 5);
        uint32_writex(wptr, uii + 1, '\0');
	param_idx++;
      }
      if (genotypic_or_hethom) {
	if (is_set(active_params, ujj++)) {
	  wptr = &(param_names[param_idx * max_param_name_len]);
	  if (hethom) {
            wptr = memcpya(wptr, "HETxCSNP", 8);
	  } else {
	    wptr = memcpya(wptr, "DOMDEVxCSNP", 11);
	  }
	  uint32_writex(wptr, uii + 1, '\0');
	  param_idx++;
	}
      }
    }
    for (uii = 0; uii < covar_ct; uii++) {
      if (is_set(active_params, ujj++)) {
        wptr = memcpyl3a(&(param_names[param_idx * max_param_name_len]), main_effect);
	*wptr++ = 'x';
	strcpy(wptr, &(covar_names[uii * max_covar_name_len]));
        param_idx++;
      }
      if (genotypic_or_hethom) {
	if (is_set(active_params, ujj++)) {
	  wptr = &(param_names[param_idx * max_param_name_len]);
	  if (hethom) {
	    wptr = memcpya(wptr, "HETx", 4);
	  } else {
	    wptr = memcpya(wptr, "DOMDEVx", 7);
	  }
	  strcpy(wptr, &(covar_names[uii * max_covar_name_len]));
	  param_idx++;
	}
      }
    }
  }
  if (variation_in_sex) {
    fptr = &(g_fixed_covars_cov_major_f[param_idx_fixed * indiv_valid_ct]);
    fill_float_zero(fptr, indiv_valid_ct);
    indiv_idx = 0;
    while (1) {
      next_set_ul_ck(sex_male_collapsed, &indiv_idx, indiv_valid_ct);
      if (indiv_idx == indiv_valid_ct) {
	break;
      }
      fptr[indiv_idx++] = 1;
    }
    if (IS_SET(active_params, sex_start_idx)) {
      memcpy(&(param_names[param_idx * max_param_name_len]), "SEX", 4);
      param_idx++;
    }
    if (covar_interactions) {
      if (is_set(active_params, sex_start_idx + 1)) {
        wptr = memcpyl3a(&(param_names[param_idx * max_param_name_len]), main_effect);
        memcpy(wptr, "xSEX", 5);
	param_idx++;
      }
      if (genotypic_or_hethom && is_set(active_params, sex_start_idx + 2)) {
	wptr = &(param_names[param_idx * max_param_name_len]);
	if (hethom) {
	  memcpy(wptr, "HETxSEX", 8);
	} else {
	  memcpy(wptr, "DOMDEVxSEX", 11);
	}
	param_idx++;
      }
    } else if (x_sex_interaction) {
      // active_params bit must be set
      memcpy(&(param_names[param_idx * max_param_name_len]), "XxSEX", 6);
      param_idx++;
    }
  }
  if (genotypic_or_hethom) {
    if (wkspace_alloc_c_checked(&haploid_param_names, np_base * max_param_name_len)) {
      goto glm_logistic_assoc_ret_NOMEM;
    }
    for (uii = 1, param_idx = 1; param_idx < np_base; uii++, param_idx++) {
      next_set_unsafe_ck(haploid_params, &uii);
      strcpy(&(haploid_param_names[param_idx * max_param_name_len]), &(param_names[uii * max_param_name_len]));
    }
  }

  if (constraint_ct_max) {
    if (wkspace_alloc_d_checked(&constraints_con_major, constraint_ct_max * param_ct_max * sizeof(double))) {
      goto glm_logistic_assoc_ret_NOMEM;
    }
    // special case: df may vary between chromosomes, so refill suffix at
    // beginning of each chromosome
    wptr = &(param_names[param_ct_max * max_param_name_len]);
    if (tests_range_list_ptr->name_ct) {
      memcpy(wptr, "USER_", 6);
    } else if (glm_modifier & GLM_TEST_ALL) {
      memcpy(wptr, "FULL_", 6);
    } else {
      memcpy(wptr, "GENO_", 6);
    }
  }
  g_constraints_con_major = constraints_con_major;
  if (!perm_maxt) {
    g_aperm_alpha = apip->alpha;
    mperm_save = 0;
  }
  if (fill_orig_chiabs) {
    if (mtest_adjust) {
      if (wkspace_alloc_d_checked(&g_orig_chisq, marker_initial_ct * sizeof(double)) ||
          wkspace_alloc_ui_checked(&tcnt, marker_initial_ct * sizeof(int32_t)) ||
          wkspace_alloc_ui_checked(&marker_idx_to_uidx, marker_initial_ct * sizeof(int32_t))) {
	goto glm_logistic_assoc_ret_NOMEM;
      }
    }
    if (do_perms) {
      if (!perm_batch_size) {
	perm_batch_size = 512;
      }
      if (wkspace_alloc_ui_checked(&g_perm_2success_ct, marker_initial_ct * sizeof(int32_t)) ||
	  wkspace_alloc_ui_checked(&g_perm_attempt_ct, marker_initial_ct * sizeof(int32_t))) {
	goto glm_logistic_assoc_ret_NOMEM;
      }
      // need this for max(T) now since we need to track permutation failures
      fill_uint_zero(g_perm_2success_ct, marker_initial_ct);
      fill_uint_zero(g_perm_attempt_ct, marker_initial_ct);
      if (perm_adapt) {
	perms_total = apip->max;
	if (perms_total < perm_batch_size) {
	  perm_batch_size = apip->max;
	}
      } else {
	perms_total = glm_mperm_val;
	if (perms_total < perm_batch_size) {
	  perm_batch_size = perms_total;
	}
      }
      uii = perm_batch_size;
      if (uii > GLM_BLOCKSIZE / CACHELINE_INT32) {
	uii = GLM_BLOCKSIZE / CACHELINE_INT32;
      }
      if (max_thread_ct > uii) {
	max_thread_ct = uii;
      }
      if (!perm_adapt) {
	ulii = (perm_batch_size + (CACHELINE_DBL - 1)) & (~(CACHELINE_DBL - 1));
        if (wkspace_alloc_d_checked(&g_maxt_thread_results, ulii * max_thread_ct * sizeof(double)),
            wkspace_alloc_d_checked(&g_maxt_extreme_stat, perms_total * sizeof(double))) {
          goto glm_logistic_assoc_ret_NOMEM;
	}
	fill_double_zero(g_maxt_extreme_stat, perms_total);
	if (mperm_save_all) {
	  if (wkspace_alloc_d_checked(&g_mperm_save_all, marker_initial_ct * perm_batch_size * sizeof(double))) {
	    goto glm_logistic_assoc_ret_NOMEM;
	  }
	  memcpy(outname_end, ".mperm.dump.all", 16);
	  if (fopen_checked(&outfile_msa, outname, "w")) {
	    goto glm_logistic_assoc_ret_OPEN_FAIL;
	  }
	  if (putc_checked('0', outfile_msa)) {
	    goto glm_logistic_assoc_ret_WRITE_FAIL;
	  }
	  LOGPRINTFWW("Dumping all permutation chi-square values to %s .\n", outname);
	}
      }
      g_perm_batch_max = perm_batch_size;
    }
  }
  if (cluster_starts) {
    // If there are any size-1 clusters, we actually want two cluster indexes:
    // - one for permutation, which excludes the size-1 clusters, and
    // - one for use by the robust cluster variance estimator, which does not.
    //
    // cluster_ct1, etc. includes size-1 clusters, while the same variables
    // without the "1" at the end exclude them.
    retval = cluster_include_and_reindex(unfiltered_indiv_ct, load_mask, 0, NULL, indiv_valid_ct, 0, cluster_ct, cluster_map, cluster_starts, &cluster_ct1, &cluster_map1, &cluster_starts1, NULL, NULL);
    if (retval) {
      goto glm_logistic_assoc_ret_1;
    }
    g_cluster_ct1 = cluster_ct1;
    if (cluster_ct1) {
      if (wkspace_alloc_ui_checked(&g_indiv_to_cluster1, indiv_valid_ct * sizeof(int32_t))) {
	goto glm_logistic_assoc_ret_NOMEM;
      }
      fill_unfiltered_indiv_to_cluster(indiv_valid_ct, cluster_ct1, cluster_map1, cluster_starts1, g_indiv_to_cluster1);
      if (do_perms) {
	retval = cluster_include_and_reindex(unfiltered_indiv_ct, load_mask, 1, pheno_c, indiv_valid_ct, 0, cluster_ct, cluster_map, cluster_starts, &g_cluster_ct, &g_cluster_map, &g_cluster_starts, &g_cluster_case_cts, &g_cluster_cc_perm_preimage);
	if (retval) {
	  goto glm_logistic_assoc_ret_1;
	}
	if (!g_cluster_ct) {
	  goto glm_logistic_assoc_ret_NO_PERMUTATION_CLUSTERS;
	}
	if (cluster_alloc_and_populate_magic_nums(g_cluster_ct, g_cluster_map, g_cluster_starts, &g_tot_quotients, &g_totq_magics, &g_totq_preshifts, &g_totq_postshifts, &g_totq_incrs)) {
	  goto glm_logistic_assoc_ret_NOMEM;
	}
	if (wkspace_alloc_ui_checked(&g_indiv_to_cluster, indiv_valid_ct * sizeof(int32_t))) {
	  goto glm_logistic_assoc_ret_NOMEM;
	}
	fill_unfiltered_indiv_to_cluster(indiv_valid_ct, g_cluster_ct, g_cluster_map, g_cluster_starts, g_indiv_to_cluster);
      }
    }
  }
  if (do_perms) {
    g_tot_quotient = 0x100000000LLU / indiv_valid_ct;
    magic_num(g_tot_quotient, &g_totq_magic, &g_totq_preshift, &g_totq_postshift, &g_totq_incr);
    if (wkspace_init_sfmtp(max_thread_ct)) {
      goto glm_logistic_assoc_ret_NOMEM;
    }
  }
  g_logistic_mt = (Logistic_multithread*)wkspace_alloc(max_thread_ct * sizeof(Logistic_multithread));
  ulii = (perm_batch_size + (BITCT - 1)) / BITCT;
  for (tidx = 0; tidx < max_thread_ct; tidx++) {
    // covars_cov_major, param_2d_buf, param_2d_buf2 matrices must have 16-byte
    // aligned rows
    // (no need to worry about 1D 16-byte alignment requirements since
    // wkspace_alloc actually forces 64-byte alignment, and allocation sizes
    // are automatically rounded up)
    if (wkspace_alloc_f_checked(&(g_logistic_mt[tidx].cur_covars_cov_major), param_ct_max * indiv_valid_cta4 * sizeof(float)) ||
	wkspace_alloc_f_checked(&(g_logistic_mt[tidx].coef), param_ct_maxa4 * perm_batch_size * sizeof(float)) ||
	wkspace_alloc_f_checked(&(g_logistic_mt[tidx].pp), indiv_valid_cta4 * sizeof(float)) ||
        wkspace_alloc_f_checked(&(g_logistic_mt[tidx].indiv_1d_buf), indiv_valid_ct * sizeof(float)) ||
        wkspace_alloc_f_checked(&(g_logistic_mt[tidx].pheno_buf), indiv_valid_ct * sizeof(float)) ||
        wkspace_alloc_f_checked(&(g_logistic_mt[tidx].param_1d_buf), param_ct_max * sizeof(float)) ||
        wkspace_alloc_f_checked(&(g_logistic_mt[tidx].param_1d_buf2), param_ct_max * sizeof(float)) ||
        wkspace_alloc_f_checked(&(g_logistic_mt[tidx].param_2d_buf), param_ct_max * param_ct_maxa4 * sizeof(float)) ||
        wkspace_alloc_f_checked(&(g_logistic_mt[tidx].param_2d_buf2), param_ct_max * param_ct_maxa4 * sizeof(float)) ||
        wkspace_alloc_f_checked(&(g_logistic_mt[tidx].regression_results), perm_batch_size * param_ctx_max_m1 * sizeof(float)) ||
        wkspace_alloc_ul_checked(&(g_logistic_mt[tidx].perm_fails), ulii * sizeof(intptr_t))) {
      goto glm_logistic_assoc_ret_NOMEM;
    }
    if (cluster_ct1) {
      ulii = MAXV(cluster_ct1 + 1, param_ct_max);
      if (wkspace_alloc_f_checked(&(g_logistic_mt[tidx].cur_covars_indiv_major), param_ct_max * indiv_valid_ct * sizeof(float)) ||
          wkspace_alloc_ui_checked(&(g_logistic_mt[tidx].cur_indiv_to_cluster1_buf), indiv_valid_ct * sizeof(int32_t)) ||
          wkspace_alloc_f_checked(&(g_logistic_mt[tidx].cluster_param_buf), ulii * param_ct_max * sizeof(float)) ||
          wkspace_alloc_f_checked(&(g_logistic_mt[tidx].cluster_param_buf2), (cluster_ct1 + 1) * param_ct_max * sizeof(float))) {
	goto glm_logistic_assoc_ret_NOMEM;
      }
    } else {
      // defensive
      g_logistic_mt[tidx].cur_covars_indiv_major = NULL;
      g_logistic_mt[tidx].cur_indiv_to_cluster1_buf = NULL;
      g_logistic_mt[tidx].cluster_param_buf = NULL;
      g_logistic_mt[tidx].cluster_param_buf2 = NULL;
    }
    if (constraint_ct_max) {
      g_logistic_mt[tidx].mi_buf = (MATRIX_INVERT_BUF1_TYPE*)wkspace_alloc(param_ct_max * sizeof(MATRIX_INVERT_BUF1_TYPE));
      if (!(g_logistic_mt[tidx].mi_buf)) {
	goto glm_logistic_assoc_ret_NOMEM;
      }
      if (wkspace_alloc_d_checked(&(g_logistic_mt[tidx].param_1d_dbuf), param_ct_max * sizeof(double)) ||
          wkspace_alloc_d_checked(&(g_logistic_mt[tidx].param_2d_dbuf), param_ct_max * param_ct_max * sizeof(double)) ||
          wkspace_alloc_d_checked(&(g_logistic_mt[tidx].param_2d_dbuf2), param_ct_max * param_ct_max * sizeof(double)) ||
          wkspace_alloc_d_checked(&(g_logistic_mt[tidx].param_df_dbuf), param_ct_max * constraint_ct_max * sizeof(double)) ||
          wkspace_alloc_d_checked(&(g_logistic_mt[tidx].df_df_dbuf), constraint_ct_max * constraint_ct_max * sizeof(double)) ||
	  wkspace_alloc_d_checked(&(g_logistic_mt[tidx].df_dbuf), constraint_ct_max * sizeof(double))) {
	goto glm_logistic_assoc_ret_NOMEM;
      }
    } else {
      g_logistic_mt[tidx].mi_buf = NULL;
      g_logistic_mt[tidx].param_1d_dbuf = NULL;
      g_logistic_mt[tidx].param_2d_dbuf = NULL;
      g_logistic_mt[tidx].param_2d_dbuf2 = NULL;
      g_logistic_mt[tidx].param_df_dbuf = NULL;
      g_logistic_mt[tidx].df_df_dbuf = NULL;
      g_logistic_mt[tidx].df_dbuf = NULL;
    }
  }

  if (wkspace_alloc_ul_checked(&pheno_c_collapsed, indiv_valid_ctv2 * sizeof(intptr_t))) {
    goto glm_logistic_assoc_ret_NOMEM;
  }
  vec_collapse_init(pheno_c, unfiltered_indiv_ct, load_mask, indiv_valid_ct, pheno_c_collapsed);
  g_case_ct = popcount_longs(pheno_c_collapsed, indiv_valid_ctv2);
  if ((!g_case_ct) || (g_case_ct == indiv_valid_ct)) {
    goto glm_logistic_assoc_ret_PHENO_CONSTANT;
  }
  if (do_perms) {
    if (wkspace_alloc_ul_checked(&g_perm_vecs, perm_batch_size * indiv_valid_ctv2 * sizeof(intptr_t))) {
      goto glm_logistic_assoc_ret_NOMEM;
    }
  }

  outname_end2 = memcpyb(outname_end, ".assoc.logistic", 16);
  if (fopen_checked(&outfile, outname, "w")) {
    goto glm_logistic_assoc_ret_OPEN_FAIL;
  }
  LOGPRINTFWW5("Writing logistic model association results to %s ... ", outname);
  fflush(stdout);
  sprintf(tbuf, " CHR %%%us         BP   A1       TEST    NMISS       %s ", plink_maxsnp, report_odds? "  OR" : "BETA");
  fprintf(outfile, tbuf, "SNP");
  if (display_ci) {
    uii = (uint32_t)((int32_t)(ci_size * 100));
    if (uii >= 10) {
      fprintf(outfile, "      SE      L%u      U%u ", uii, uii);
    } else {
      fprintf(outfile, "      SE       L%u       U%u ", uii, uii);
    }
  }
  fputs("        STAT            P \n", outfile);
  loop_end = marker_initial_ct / 100;
  marker_unstopped_ct = marker_initial_ct;
  g_adaptive_ci_zt = ltqnorm(1 - apip->beta / (2.0 * ((int32_t)marker_initial_ct)));

  // ----- begin main loop -----
 glm_logistic_assoc_more_perms:
  if (do_perms) {
    if (perm_adapt) {
      if (perm_pass_idx) {
	while (g_first_adapt_check <= g_perms_done) {
	  g_first_adapt_check += (int32_t)(apip->init_interval + ((int32_t)g_first_adapt_check) * apip->interval_slope);
	}
      } else {
	if (apip->min < apip->init_interval) {
	  g_first_adapt_check = (int32_t)(apip->init_interval);
	} else {
	  g_first_adapt_check = apip->min;
	}
	g_adaptive_intercept = apip->init_interval;
        g_adaptive_slope = apip->interval_slope;
      }
    }
    g_perm_vec_ct = perm_batch_size;
    if (perm_adapt && (g_perms_done < perm_batch_size)) {
      // special case: split first batch to reduce adaptive overshoot
      ulii = perm_batch_size;
      uljj = (intptr_t)(apip->init_interval);
      uljj = MAXV(uljj, apip->min);
      uljj *= 2;
      uljj = MAXV(64, uljj);
      while (ulii >= (uljj << perm_pass_idx)) {
	ulii >>= 1;
      }
      g_perm_vec_ct = ulii - g_perms_done;
    }
    if (g_perm_vec_ct > perms_total - g_perms_done) {
      g_perm_vec_ct = perms_total - g_perms_done;
    }
    ulii = 0;
    if (g_perm_vec_ct > max_thread_ct) {
      g_assoc_thread_ct = max_thread_ct;
    } else {
      g_assoc_thread_ct = g_perm_vec_ct;
    }
    if (!g_cluster_ct) {
      if (spawn_threads(threads, &logistic_gen_perms_thread, g_assoc_thread_ct)) {
	goto glm_logistic_assoc_ret_THREAD_CREATE_FAIL;
      }
      logistic_gen_perms_thread((void*)ulii);
    } else {
      if (spawn_threads(threads, &logistic_gen_cluster_perms_thread, g_assoc_thread_ct)) {
	goto glm_logistic_assoc_ret_THREAD_CREATE_FAIL;
      }
      logistic_gen_cluster_perms_thread((void*)ulii);
    }
    join_threads(threads, g_assoc_thread_ct);
  }
  chrom_fo_idx = 0xffffffffU;
  marker_uidx = next_unset_unsafe(marker_exclude, 0);
  if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
    goto glm_logistic_assoc_ret_READ_FAIL;
  }
  if (!perm_pass_idx) {
    fputs("0%", stdout);
    fflush(stdout);
  }
  marker_idx = 0;
  marker_idx2 = 0;
  chrom_end = 0;
  do {
    if (marker_uidx >= chrom_end) {
      // exploit overflow
      do {
        chrom_fo_idx++;
        refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &g_is_x, &g_is_y, &uii, &g_is_haploid);
      } while ((!glm_xchr_model) && (g_is_haploid || uii));
      uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
      wptr_start = width_force(4, writebuf, chrom_name_write(writebuf, chrom_info_ptr, uii, zero_extra_chroms));
      *wptr_start++ = ' ';
      fill_double_zero(constraints_con_major, constraint_ct_max * param_ct_max);
      g_male_x_01 = 0;
      if (g_is_x) {
        cur_param_ct = param_ct_max;
	cur_constraint_ct = constraint_ct_max;
	cur_param_names = param_names;
	if (glm_xchr_model != 2) {
	  g_male_x_01 = 1;
	}
      } else if ((!g_is_haploid) || (!genotypic_or_hethom)) {
	cur_param_ct = np_base + np_diploid;
	if (constraint_ct_max) {
          cur_constraint_ct = popcount_bit_idx(g_joint_test_params, 0, cur_param_ct);
	} else {
	  cur_constraint_ct = 0;
	}
	cur_param_names = param_names;
      } else {
	cur_param_ct = np_base;
	cur_constraint_ct = 0;
	if (constraint_ct_max) {
	  for (ulii = 0; ulii < param_ctl_max; ulii++) {
	    uljj = g_joint_test_params[ulii] & haploid_params[ulii];
	    while (uljj) {
	      uii = CTZLU(uljj) + ulii * BITCT;
              constraints_con_major[cur_constraint_ct * cur_param_ct + uii + 1] = 1;
	      cur_constraint_ct++;
	      uljj &= uljj - 1;
	    }
	  }
	}
	cur_param_names = haploid_param_names;
      }
      g_cur_param_ct = cur_param_ct;
      g_cur_constraint_ct = cur_constraint_ct;
      if (!hide_covar) {
	param_idx_end = cur_param_ct;
      }
      cur_param_ctx = cur_param_ct;
      g_include_sex = sex_covar_everywhere || (g_is_x && np_sex);
      if (cur_constraint_ct) {
	cur_param_ctx++;
	if (g_is_x || (!g_is_haploid)) {
	  ulii = 0;
	  for (constraint_idx = 0; constraint_idx < cur_constraint_ct; ulii++, constraint_idx++) {
            next_set_ul_unsafe_ck(g_joint_test_params, &ulii);
	    constraints_con_major[constraint_idx * cur_param_ct + ulii + 1] = 1;
	  }
	}
        wptr = uint32_write(&(param_names[param_ct_max * max_param_name_len + 5]), cur_constraint_ct);
	memcpy(wptr, "DF", 3);
      }
    }
    block_size = 0;
    block_end = marker_unstopped_ct - marker_idx;
    if (block_end > GLM_BLOCKSIZE) {
      block_end = GLM_BLOCKSIZE;
    }
    do {
      if (g_perm_adapt_stop[marker_idx2]) {
	do {
	  marker_uidx++;
	  next_unset_unsafe_ck(marker_exclude, &marker_uidx);
	  marker_idx2++;
	} while ((marker_uidx < chrom_end) && g_perm_adapt_stop[marker_idx2]);
	if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	  goto glm_logistic_assoc_ret_READ_FAIL;
	}
	if (marker_uidx >= chrom_end) {
	  break;
	}
      }
      loadbuf_ptr = &(g_loadbuf[block_size * indiv_valid_ctv2]);
      if (load_and_collapse_incl(bedfile, loadbuf_raw, unfiltered_indiv_ct, loadbuf_ptr, indiv_valid_ct, load_mask, IS_SET(marker_reverse, marker_uidx))) {
	goto glm_logistic_assoc_ret_READ_FAIL;
      }
      if (hh_exists) {
	haploid_fix(hh_exists, indiv_include2, indiv_male_include2, indiv_valid_ct, g_is_x, g_is_y, (unsigned char*)loadbuf_ptr);
      }
      g_adapt_m_table[block_size] = marker_idx2++;
      mu_table[block_size++] = marker_uidx;
      if (marker_idx + block_size == marker_unstopped_ct) {
	break;
      }
      marker_uidx++;
      if (IS_SET(marker_exclude, marker_uidx)) {
	marker_uidx = next_unset_unsafe(marker_exclude, marker_uidx);
	if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	  goto glm_logistic_assoc_ret_READ_FAIL;
	}
      }
    } while ((block_size < block_end) && (marker_uidx < chrom_end));
    if (!block_size) {
      continue;
    }
    g_block_diff = block_size;
    if (!perm_pass_idx) {
      for (marker_bidx = 0; marker_bidx < block_size; marker_bidx++) {
	// 1. fill design matrix and make transposed copy; fill g_nm_cts[] in
	//    the process
	// 2. logistic regression, set g_perm_adapt_stop byte on failure
	// 3. write basic report to disk
	marker_uidx2 = mu_table[marker_bidx];
        loadbuf_ptr = &(g_loadbuf[marker_bidx * indiv_valid_ctv2]);
	marker_idx3 = g_adapt_m_table[marker_bidx];
	if (marker_idx_to_uidx) {
	  marker_idx_to_uidx[marker_idx3] = marker_uidx2;
	}
	cur_missing_ct = glm_fill_design_float(loadbuf_ptr, g_fixed_covars_cov_major_f, indiv_valid_ct, g_indiv_to_cluster1, cur_param_ct, glm_modifier & (GLM_HETHOM | GLM_DOMINANT | GLM_RECESSIVE), glm_xchr_model, condition_list_start_idx, interaction_start_idx, sex_start_idx, active_params, haploid_params, g_include_sex, g_male_x_01, sex_male_collapsed, g_is_haploid && (!g_is_x), g_logistic_mt[0].cur_covars_cov_major, g_logistic_mt[0].cur_covars_indiv_major, g_logistic_mt[0].cur_indiv_to_cluster1_buf, &cur_indiv_to_cluster1);
	cur_indiv_valid_ct = indiv_valid_ct - cur_missing_ct;
	g_nm_cts[marker_idx3] = cur_indiv_valid_ct;
	if (cur_indiv_valid_ct > cur_param_ct) {
	  // todo: try better starting position
	  fill_float_zero(g_logistic_mt[0].coef, (cur_param_ct + 3) & (~3));
	  regression_fail = glm_logistic_robust_cluster_covar(1, cur_param_ct, cur_indiv_valid_ct, cur_missing_ct, loadbuf_ptr, g_logistic_mt[0].cur_covars_cov_major, g_logistic_mt[0].cur_covars_indiv_major, pheno_c_collapsed, g_logistic_mt[0].coef, g_logistic_mt[0].pp, g_logistic_mt[0].indiv_1d_buf, g_logistic_mt[0].pheno_buf, g_logistic_mt[0].param_1d_buf, g_logistic_mt[0].param_1d_buf2, g_logistic_mt[0].param_2d_buf, g_logistic_mt[0].param_2d_buf2, g_logistic_mt[0].regression_results, cluster_ct1, cur_indiv_to_cluster1, g_logistic_mt[0].cluster_param_buf, g_logistic_mt[0].cluster_param_buf2, cur_constraint_ct, constraints_con_major, g_logistic_mt[0].param_1d_dbuf, g_logistic_mt[0].param_2d_dbuf, g_logistic_mt[0].param_2d_dbuf2, g_logistic_mt[0].param_df_dbuf, g_logistic_mt[0].df_df_dbuf, g_logistic_mt[0].mi_buf, g_logistic_mt[0].df_dbuf, g_logistic_mt[0].perm_fails);
	} else {
	  regression_fail = 1;
	}
	wptr_start2 = fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx2 * max_marker_id_len]), wptr_start);
	*wptr_start2++ = ' ';
	wptr_start2 = uint32_writew10(wptr_start2, marker_pos[marker_uidx2]);
	*wptr_start2++ = ' ';
	wptr_start2 = fw_strcpy(4, marker_allele_ptrs[marker_uidx2 * 2], wptr_start2);
	*wptr_start2++ = ' ';
	orig_stats_ptr = &(g_orig_stats[marker_idx3]);
	if (!regression_fail) {
	  for (param_idx = 1; param_idx < cur_param_ct; param_idx++) {
	    dxx = (double)g_logistic_mt[0].coef[param_idx];
	    se = sqrt((double)g_logistic_mt[0].regression_results[param_idx - 1]);
	    zval = dxx / se;
	    if (param_idx == 1) {
	      if (mtest_adjust) {
		g_orig_chisq[marker_idx3] = zval * zval;
	      }
	      if (mperm_save_all) {
		fprintf(outfile_msa, " %g", zval * zval);
	      }
	      if (!cur_constraint_ct) {
		*orig_stats_ptr = zval * zval;
	      }
	    }
	    pval = chiprob_p(zval * zval, 1);
	    if ((param_idx < param_idx_end) && (pval <= pfilter)) {
	      wptr = fw_strcpy(10, &(cur_param_names[param_idx * max_param_name_len]), wptr_start2);
	      *wptr++ = ' ';
	      wptr = uint32_writew8x(wptr, (uint32_t)cur_indiv_valid_ct, ' ');
	      wptr = double_g_writewx4x(wptr, report_odds? exp(dxx) : dxx, 10, ' ');
	      if (display_ci) {
		dyy = ci_zt * se;
		wptr = double_g_writewx4x(wptr, se, 8, ' ');
		if (report_odds) {
		  wptr = double_g_writewx4x(wptr, exp(dxx - dyy), 8, ' ');
		  wptr = double_g_writewx4x(wptr, exp(dxx + dyy), 8, ' ');
		} else {
		  wptr = double_g_writewx4x(wptr, dxx - dyy, 8, ' ');
		  wptr = double_g_writewx4x(wptr, dxx + dyy, 8, ' ');
		}
	      }
	      wptr = double_g_writewx4x(wptr, zval, 12, ' ');
	      wptr = double_g_writewx4x(wptr, pval, 12, '\n');
	      if (fwrite_checked(writebuf, wptr - writebuf, outfile)) {
		goto glm_logistic_assoc_ret_WRITE_FAIL;
	      }
	    }
	  }
	  if (cur_constraint_ct) {
	    dxx = (double)g_logistic_mt[0].regression_results[cur_param_ct - 1];
	    *orig_stats_ptr = dxx;
	    pval = chiprob_p(dxx, cur_constraint_ct);
	    if (pval <= pfilter) {
	      wptr = fw_strcpy(10, &(param_names[param_ct_max * max_param_name_len]), wptr_start2);
              *wptr++ = ' ';
              wptr = uint32_writew8(wptr, (uint32_t)cur_indiv_valid_ct);
              wptr = memcpya(wptr, "         NA ", 12);
              if (display_ci) {
		wptr = memcpya(wptr, "      NA       NA       NA ", 27);
	      }
              wptr = double_g_writewx4x(wptr, dxx, 12, ' ');
              wptr = double_g_writewx4x(wptr, pval, 12, '\n');
              if (fwrite_checked(writebuf, wptr - writebuf, outfile)) {
		goto glm_logistic_assoc_ret_WRITE_FAIL;
	      }
	    }
	  }
	} else {
	  g_perm_adapt_stop[marker_idx3] = 1;
	  if (mtest_adjust && (param_idx == 1)) {
	    g_orig_chisq[marker_idx3] = -9;
	    if (tcnt) {
	      tcnt[marker_idx3] = 0;
	    }
	  }
	  *orig_stats_ptr = -9;
	  if (mperm_save_all) {
	    fputs(" NA", outfile_msa);
	    msa_ptr = &(g_mperm_save_all[marker_idx3 * perm_batch_size]);
	    for (perm_idx = 0; perm_idx < perm_batch_size; perm_idx++) {
              *msa_ptr++ = -9;
	    }
	  }
	  for (param_idx = 1; param_idx < cur_param_ctx; param_idx++) {
	    if ((param_idx < param_idx_end) || (param_idx == cur_param_ct)) {
	      if (!(param_idx == cur_param_ct)) {
	        wptr = fw_strcpy(10, &(cur_param_names[param_idx * max_param_name_len]), wptr_start2);
	      } else {
		wptr = fw_strcpy(10, &(param_names[param_ct_max * max_param_name_len]), wptr_start2);
	      }
	      *wptr++ = ' ';
	      wptr = uint32_writew8(wptr, (uint32_t)cur_indiv_valid_ct);
	      wptr = memcpya(wptr, "         NA ", 12);
	      if (display_ci) {
		wptr = memcpya(wptr, "      NA       NA       NA ", 27);
	      }
	      wptr = memcpya(wptr, "          NA           NA\n", 26);
	      if (fwrite_checked(writebuf, wptr - writebuf, outfile)) {
		goto glm_logistic_assoc_ret_WRITE_FAIL;
	      }
	    }
	  }
	}
      }
    }
    if (do_perms) {
      // split loaded markers evenly between threads.  each thread does the
      // following for every marker in its set:
      // 1. fill design matrix and make transposed copy, do this more
      //    efficiently if g_nm_cts[] indicates there are no missing genotypes
      // 2. logistic regression on permutations, handle adaptive logic,
      //    --mperm-save[-all], etc.
      g_assoc_thread_ct = block_size / CACHELINE_INT32;
      if (g_assoc_thread_ct > max_thread_ct) {
	g_assoc_thread_ct = max_thread_ct;
      } else if (!g_assoc_thread_ct) {
	g_assoc_thread_ct = 1;
      }
      ulii = 0;
      if (perm_adapt) {
	if (spawn_threads(threads, &glm_logistic_adapt_thread, g_assoc_thread_ct)) {
	  goto glm_logistic_assoc_ret_THREAD_CREATE_FAIL;
	}
	glm_logistic_adapt_thread((void*)ulii);
	join_threads(threads, g_assoc_thread_ct);
      } else {
	if (spawn_threads(threads, &glm_logistic_maxt_thread, g_assoc_thread_ct)) {
	  goto glm_logistic_assoc_ret_THREAD_CREATE_FAIL;
	}
	glm_logistic_maxt_thread((void*)ulii);
	join_threads(threads, g_assoc_thread_ct);
        ulii = (g_perm_vec_ct + (CACHELINE_DBL - 1)) & (~(CACHELINE_DBL - 1));
	ukk = g_perms_done + g_perm_vec_ct;
        for (uii = 0; uii < g_assoc_thread_ct; uii++) {
          dptr = &(g_maxt_thread_results[uii * ulii]);
          for (ujj = g_perms_done; ujj < ukk; ujj++) {
            dxx = *dptr++;
            if (dxx > g_maxt_extreme_stat[ujj]) {
	      g_maxt_extreme_stat[ujj] = dxx;
	    }
	  }
	}
      }
    }
    marker_idx += block_size;
    if ((!perm_pass_idx) && (marker_idx >= loop_end)) {
      if (marker_idx < marker_initial_ct) {
	if (pct >= 10) {
	  putchar('\b');
	}
        pct = (marker_idx * 100LLU) / marker_initial_ct;
        printf("\b\b%u%%", pct);
        fflush(stdout);
        loop_end = (((uint64_t)pct + 1LLU) * marker_initial_ct) / 100;
      }
    }
  } while (marker_idx < marker_unstopped_ct);
  // if more permutations, reevaluate marker_unstopped_ct, etc.
  if (!perm_pass_idx) {
    if (pct >= 10) {
      putchar('\b');
    }
    fputs("\b\b", stdout);
    logprint("done.\n");
    if (fclose_null(&outfile)) {
      goto glm_logistic_assoc_ret_WRITE_FAIL;
    }
    if (mtest_adjust) {
      retval = multcomp(outname, outname_end, marker_idx_to_uidx, marker_initial_ct, marker_ids, max_marker_id_len, plink_maxsnp, zero_extra_chroms, chrom_info_ptr, g_orig_chisq, pfilter, mtest_adjust, adjust_lambda, 0, NULL);
      if (retval) {
	goto glm_logistic_assoc_ret_1;
      }
    }
    if (mperm_save_all) {
      if (putc_checked('\n', outfile_msa)) {
	goto glm_logistic_assoc_ret_WRITE_FAIL;
      }
    }
  }
  if (do_perms) {
    if (mperm_save_all) {
      if (perm_pass_idx) {
	putchar(' ');
      }
      fputs("[dumping stats]", stdout);
      fflush(stdout);
      ulii = g_perm_vec_ct;
      ujj = 1 + g_perms_done;
      wptr = tbuf;
      cptr = &(tbuf[MAXLINELEN]);
      for (uii = 0; uii < ulii; uii++) {
	wptr = uint32_write(wptr, uii + ujj);
	dptr = &(g_mperm_save_all[uii]);
	for (ukk = 0; ukk < marker_ct; ukk++) {
	  *wptr++ = ' ';
	  dxx = dptr[ukk * ulii];
	  if (dxx >= 0) {
	    wptr = double_g_write(wptr, dxx);
	  } else {
	    wptr = memcpya(wptr, "NA", 2);
	  }
	  if (wptr >= cptr) {
	    if (fwrite_checked(tbuf, (uintptr_t)(wptr - tbuf), outfile_msa)) {
	      goto glm_logistic_assoc_ret_WRITE_FAIL;
	    }
	    wptr = tbuf;
	  }
	}
	*wptr++ = '\n';
      }
      if (fwrite_checked(tbuf, (uintptr_t)(wptr - tbuf), outfile_msa)) {
	goto glm_logistic_assoc_ret_WRITE_FAIL;
      }
      fputs("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b               ", stdout);
    }
    g_perms_done += g_perm_vec_ct;
    if (g_perms_done < perms_total) {
      if (perm_adapt || (!perm_pass_idx)) {
        marker_unstopped_ct = marker_initial_ct - popcount_longs((uintptr_t*)g_perm_adapt_stop, (marker_initial_ct + sizeof(intptr_t) - 1) / sizeof(intptr_t));
        if (!marker_unstopped_ct) {
          goto glm_logistic_assoc_perm_count;
	}
      }
      printf("\r%u permutation%s complete.", g_perms_done, (g_perms_done != 1)? "s" : "");
      fflush(stdout);
      perm_pass_idx++;
      goto glm_logistic_assoc_more_perms;
    }
  glm_logistic_assoc_perm_count:
    if (perm_adapt) {
      g_perms_done = 0;
      for (uii = 0; uii < marker_initial_ct; uii++) {
        if (g_perm_attempt_ct[uii] > g_perms_done) {
	  g_perms_done = g_perm_attempt_ct[uii];
	  if (g_perms_done == perms_total) {
	    break;
	  }
	}
      }
    }
    putchar('\r');
    LOGPRINTF("%u %s permutation%s complete.\n", g_perms_done, perm_maxt? "max(T)" : "(adaptive)", (g_perms_done != 1)? "s" : "");

    if (perm_adapt) {
      memcpy(outname_end2, ".perm", 6);
    } else {
      if (mperm_save & MPERM_DUMP_BEST) {
	memcpy(outname_end, ".mperm.dump.best", 17);
	LOGPRINTF("Dumping best permutation chi-square values to %s .\n", outname);
	if (fopen_checked(&outfile, outname, "w")) {
	  goto glm_logistic_assoc_ret_OPEN_FAIL;
	}
	dxx = 0;
	for (marker_idx = 0; marker_idx < marker_initial_ct; marker_idx++) {
          if (g_orig_stats[marker_idx] > dxx) {
	    dxx = g_orig_stats[marker_idx];
	  }
	}
	memcpy(tbuf, "0 ", 2);
	wptr = double_g_writex(&(tbuf[2]), dxx, '\n');
        if (fwrite_checked(tbuf, (uintptr_t)(wptr - tbuf), outfile)) {
	  goto glm_logistic_assoc_ret_WRITE_FAIL;
	}
        for (uii = 0; uii < perms_total; uii++) {
          wptr = uint32_writex(tbuf, uii + 1, ' ');
          wptr = double_g_writex(wptr, g_maxt_extreme_stat[uii], '\n');
          if (fwrite_checked(tbuf, (uintptr_t)(wptr - tbuf), outfile)) {
	    goto glm_logistic_assoc_ret_WRITE_FAIL;
	  }
	}
        if (fclose_null(&outfile)) {
	  goto glm_logistic_assoc_ret_WRITE_FAIL;
	}
	memcpy(outname_end, ".assoc.logistic", 15);
      }
      memcpy(outname_end2, ".mperm", 7);
    }
    if (fopen_checked(&outfile, outname, "w")) {
      goto glm_logistic_assoc_ret_OPEN_FAIL;
    }
    if (perm_adapt) {
      sprintf(tbuf, " CHR %%%us         EMP1           NP \n", plink_maxsnp);
    } else {
      sprintf(tbuf, " CHR %%%us         EMP1         EMP2 \n", plink_maxsnp);
#ifdef __cplusplus
      std::sort(g_maxt_extreme_stat, &(g_maxt_extreme_stat[perms_total]));
#else
      qsort(g_maxt_extreme_stat, perms_total, sizeof(double), double_cmp);
#endif
      // (debugging)
      // printf("extreme stats: %g %g %g\n", g_maxt_extreme_stat[0], g_maxt_extreme_stat[(perms_total - 1) / 2], g_maxt_extreme_stat[perms_total - 1]);

      for (marker_idx = 0; marker_idx < marker_initial_ct; marker_idx++) {
	g_perm_attempt_ct[marker_idx] = perms_total - g_perm_attempt_ct[marker_idx];
      }
    }
    fprintf(outfile, tbuf, "SNP");
    chrom_fo_idx = 0xffffffffU;
    marker_uidx = next_unset_unsafe(marker_exclude, 0);
    marker_idx = 0;
    while (1) {
      do {
	chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[(++chrom_fo_idx) + 1];
      } while (marker_uidx >= chrom_end);
      uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
      wptr_start = width_force(4, tbuf, chrom_name_write(tbuf, chrom_info_ptr, uii, zero_extra_chroms));
      *wptr_start++ = ' ';
      wptr_start[plink_maxsnp] = ' ';
      for (; marker_uidx < chrom_end;) {
	pval = ((double)(g_perm_2success_ct[marker_idx] + 2)) / ((double)(2 * (g_perm_attempt_ct[marker_idx] + 1)));
        if (pval <= pfilter) {
          fw_strcpy(plink_maxsnp, &(marker_ids[marker_uidx * max_marker_id_len]), wptr_start);
          wptr = &(wptr_start[1 + plink_maxsnp]);
	  if (g_orig_stats[marker_idx] == -9) {
            wptr = memcpya(wptr, "          NA           NA", 25);
	  } else {
	    if (!perm_count) {
	      wptr = double_g_writewx4x(wptr, pval, 12, ' ');
	    } else {
	      wptr = double_g_writewx4x(wptr, ((double)g_perm_2success_ct[marker_idx]) / 2.0, 12, ' ');
	    }
	    if (perm_adapt) {
	      wptr = memseta(wptr, 32, 2);
	      wptr = uint32_writew10(wptr, g_perm_attempt_ct[marker_idx]);
	    } else {
	      dzz = (int32_t)(perms_total - doublearr_greater_than(g_maxt_extreme_stat, perms_total, g_orig_stats[marker_idx] - EPSILON) + 1);
              if (!perm_count) {
		wptr = double_g_writewx4(wptr, dzz / ((double)((int32_t)perms_total + 1)), 12);
	      } else {
                wptr = double_g_writewx4(wptr, dzz, 12);
	      }
	    }
	  }
	  wptr = memcpya(wptr, " \n", 2);
	  if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	    goto glm_logistic_assoc_ret_WRITE_FAIL;
	  }
	}
	if (++marker_idx == marker_ct) {
	  goto glm_logistic_assoc_loop_end;
	}
        marker_uidx++;
        next_unset_unsafe_ck(marker_exclude, &marker_uidx);
      }
    }
  glm_logistic_assoc_loop_end:
    if (fclose_null(&outfile)) {
      goto glm_logistic_assoc_ret_WRITE_FAIL;
    }
    LOGPRINTF("Permutation test report written to %s .\n", outname);
  }

  while (0) {
  glm_logistic_assoc_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  glm_logistic_assoc_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  glm_logistic_assoc_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  glm_logistic_assoc_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  glm_logistic_assoc_ret_NO_PERMUTATION_CLUSTERS:
    logprint("Error: No size 2+ clusters for permutation test.\n");
    retval = RET_INVALID_CMDLINE;
    break;
  glm_logistic_assoc_ret_PHENO_CONSTANT:
    logprint("Warning: Skipping --linear/--logistic since phenotype is constant.\n");
    break;
  glm_logistic_assoc_ret_THREAD_CREATE_FAIL:
    logprint(errstr_thread_create);
    retval = RET_THREAD_CREATE_FAIL;
    break;
  }
 glm_logistic_assoc_ret_1:
  wkspace_reset(wkspace_mark);
  fclose_cond(outfile);
  fclose_cond(outfile_msa);
  free_cond(condition_uidxs);
  return retval;
}

#ifndef NOLAPACK
int32_t glm_linear_nosnp(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uint32_t glm_modifier, double glm_vif_thresh, uint32_t glm_xchr_model, uint32_t glm_mperm_val, Range_list* parameters_range_list_ptr, Range_list* tests_range_list_ptr, double ci_size, double ci_zt, double pfilter, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t* marker_reverse, char* condition_mname, char* condition_fname, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_indiv_ct, uintptr_t indiv_ct, uintptr_t* indiv_exclude, uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, uint32_t mperm_save, uint32_t pheno_nm_ct, uintptr_t* pheno_nm, double* pheno_d, uintptr_t covar_ct, char* covar_names, uintptr_t max_covar_name_len, uintptr_t* covar_nm, double* covar_d, uintptr_t* sex_nm, uintptr_t* sex_male, uint32_t hh_exists, uint32_t perm_batch_size, Set_info* sip) {
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + BITCT - 1) / BITCT;
  uintptr_t unfiltered_indiv_ctv2 = 2 * unfiltered_indiv_ctl;
  FILE* outfile = NULL;
  uintptr_t indiv_uidx = 0;
  uintptr_t topsize = 0;
  uintptr_t max_param_name_len = 2;
  uintptr_t param_raw_ct = 1;
  uintptr_t condition_ct = 0;
  uintptr_t constraint_ct = 0;
  uintptr_t ulii = 0;
  uint32_t cur_batch_size = 1;
  uint32_t cluster_ct1 = 0;
  uint32_t do_perms = (glm_modifier / GLM_MPERM) & 1;
  uint32_t perm_count = glm_modifier & GLM_PERM_COUNT;
  uint32_t hide_covar = glm_modifier & GLM_HIDE_COVAR;
  uint32_t display_ci = (ci_size > 0);
  uint32_t variation_in_sex = 0; // no need to initialize if no-x-sex specified
  uint32_t perms_done = 0;
  uint32_t perm_fail_total = 0;
  uint32_t joint_perm_fail_extra = 0;
  uint32_t perm_fail_ct = 0;
  uint32_t max_thread_ct = g_thread_ct;
  int32_t retval = 0;
  uint32_t linear_intercept = glm_modifier & GLM_INTERCEPT;
  char dgels_trans = 'N';
  __CLPK_integer dgels_m = 0;
  __CLPK_integer dgels_n = 0;
  __CLPK_integer dgels_nrhs = 0;
  double* dgels_a = NULL;
  double* dgels_b = NULL;
  double* dgels_work = NULL;
  double* param_df_buf = NULL;
  double* param_df_buf2 = NULL;
  __CLPK_integer dgels_ldb = 0;
  __CLPK_integer dgels_lwork = -1;
  __CLPK_integer dgels_info;
  double dzz;
  double* regression_results = NULL;
  double* cluster_param_buf = NULL;
  double* cluster_param_buf2 = NULL;
  double* indiv_1d_buf = NULL;
  double* mperm_save_stats = NULL;
  double* constraints_con_major = NULL;
  double* df_df_buf = NULL;
  double* df_buf = NULL;
  uintptr_t* loadbuf_raw = NULL;
  uintptr_t* loadbuf_collapsed = NULL;
  uintptr_t* load_mask = NULL;
  uintptr_t* sex_male_collapsed = NULL;
  uintptr_t* indiv_include2 = NULL;
  uintptr_t* indiv_male_include2 = NULL;
  uintptr_t* active_params = NULL;
  uintptr_t* joint_test_params = NULL;
  uintptr_t* perm_fails = NULL;
  uint32_t* perm_2success_ct = NULL;
  uint32_t* condition_uidxs = NULL;
  uint32_t* cluster_map1 = NULL;
  uint32_t* cluster_starts1 = NULL;
  uint32_t* indiv_to_cluster1 = NULL;
  double geno_map[12];
  double* geno_map_ptr;
  double* param_2d_buf;
  double* param_2d_buf2;
  MATRIX_INVERT_BUF1_TYPE* mi_buf;
  char* param_names;
  double* covars_cov_major;
  double* covars_indiv_major;
  double* orig_stats;
  char* outname_end2;
  char* wptr;
  char* wptr2;
  double* dptr;
  double* dptr2;
  uintptr_t indiv_valid_ct;
  uintptr_t indiv_valid_ctv2;
  uintptr_t indiv_uidx_stop;
  uintptr_t indiv_idx;
  uintptr_t param_raw_ctl;
  uintptr_t param_ct;
  uintptr_t param_ctx; // param_ct + 1 if joint test needed, param_ct otherwise
  uintptr_t param_idx;
  uintptr_t constraint_idx;
  uintptr_t uljj;
  double* msa_ptr;
  double* dptr3;
  double se;
  double zval;
  double pval;
  double dxx;
  double dyy;
  uint32_t perm_idx;
  uint32_t marker_uidx;
  uint32_t chrom_idx;
  uint32_t is_x;
  uint32_t is_y;
  uint32_t uii;
  uint32_t ujj;
  uint32_t slen;
  int32_t ii;
  if (hide_covar && (!tests_range_list_ptr->name_ct) && (!(glm_modifier & GLM_TEST_ALL))) {
    logprint("Error: --linear hide-covar no-snp produces no output.\n");
    goto glm_linear_nosnp_ret_INVALID_CMDLINE;
  }
  if (glm_init_load_mask(indiv_exclude, pheno_nm, covar_nm, indiv_ct, unfiltered_indiv_ctv2, &load_mask)) {
    goto glm_linear_nosnp_ret_NOMEM;
  }
  indiv_valid_ct = popcount_longs(load_mask, unfiltered_indiv_ctl);
  if (condition_mname || condition_fname) {
    loadbuf_raw = (uintptr_t*)top_alloc(&topsize, unfiltered_indiv_ctv2 * sizeof(intptr_t));
    if (!loadbuf_raw) {
      goto glm_linear_nosnp_ret_NOMEM;
    }
    loadbuf_raw[unfiltered_indiv_ctv2 - 2] = 0;
    loadbuf_raw[unfiltered_indiv_ctv2 - 1] = 0;
    ulii = topsize;

    if (hh_exists & (Y_FIX_NEEDED | NXMHH_EXISTS)) {
      indiv_include2 = (uintptr_t*)top_alloc(&topsize, unfiltered_indiv_ctv2 * sizeof(intptr_t));
      if (!indiv_include2) {
        goto glm_linear_nosnp_ret_NOMEM;
      }
      fill_vec_55(indiv_include2, unfiltered_indiv_ct); // harmless
    }
    if (hh_exists & (XMHH_EXISTS | Y_FIX_NEEDED)) {
      indiv_male_include2 = (uintptr_t*)top_alloc(&topsize, unfiltered_indiv_ctv2 * sizeof(intptr_t));
      if (!indiv_male_include2) {
	goto glm_linear_nosnp_ret_NOMEM;
      }
      fill_ulong_zero(indiv_male_include2, unfiltered_indiv_ctv2);
      vec_include_init(unfiltered_indiv_ct, indiv_male_include2, sex_male);
    }
    wkspace_left -= topsize;
    retval = glm_scan_conditions(condition_mname, condition_fname, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, chrom_info_ptr, hh_exists, loadbuf_raw, bedfile, bed_offset, unfiltered_indiv_ct, sex_male, load_mask, &indiv_valid_ct, &condition_ct, &condition_uidxs, indiv_include2, indiv_male_include2);
    wkspace_left += topsize;
    if (retval) {
      goto glm_linear_nosnp_ret_1;
    }
    topsize = ulii; // deallocate temporary indiv[_male]_include2
    param_raw_ct += condition_ct;
  }
  param_raw_ct += covar_ct;
  if (glm_modifier & GLM_SEX) {
    indiv_uidx = 0;
    indiv_idx = 0;
    uii = 2;
    do {
      indiv_uidx = next_set_ul_unsafe(load_mask, indiv_uidx);
      indiv_uidx_stop = next_unset_ul(load_mask, indiv_uidx, unfiltered_indiv_ct);
      indiv_idx += indiv_uidx_stop - indiv_uidx;
      do {
	if (IS_SET(sex_nm, indiv_uidx)) {
	  ujj = is_set(sex_male, indiv_uidx);
	  if (uii == ujj) {
	    variation_in_sex = 1;
	    indiv_idx = indiv_valid_ct;
	    break;
	  }
          uii = 1 - ujj;
	}
      } while (++indiv_uidx < indiv_uidx_stop);
    } while (indiv_idx < indiv_valid_ct);
    if (variation_in_sex) {
      param_raw_ct++;
      bitfield_and(load_mask, sex_nm, unfiltered_indiv_ctl);
      indiv_valid_ct = popcount_longs(load_mask, unfiltered_indiv_ctl);
    } else {
      logprint("Warning: Ignoring --linear 'sex' modifier since sex is invariant.\n");
    }
  }
  indiv_valid_ctv2 = 2 * ((indiv_valid_ct + BITCT - 1) / BITCT);

  if (condition_mname || condition_fname) {
    // now that we've determined which individuals will be in the regression,
    // initialize collapsed indiv_include2, indiv_male_include2, sex_male
    if (hh_exists & (Y_FIX_NEEDED | NXMHH_EXISTS)) {
      indiv_include2 = (uintptr_t*)top_alloc(&topsize, indiv_valid_ctv2 * sizeof(intptr_t));
      fill_vec_55(indiv_include2, indiv_valid_ct);
    }
    if (hh_exists & (XMHH_EXISTS | Y_FIX_NEEDED)) {
      indiv_male_include2 = (uintptr_t*)top_alloc(&topsize, indiv_valid_ctv2 * sizeof(intptr_t));
      alloc_collapsed_haploid_filters(unfiltered_indiv_ct, indiv_valid_ct, hh_exists, 1, load_mask, sex_male, &indiv_include2, &indiv_male_include2);
    }
    loadbuf_collapsed = (uintptr_t*)top_alloc(&topsize, indiv_valid_ctv2 * sizeof(intptr_t));
    if (!loadbuf_collapsed) {
      goto glm_linear_nosnp_ret_NOMEM;
    }
    loadbuf_collapsed[indiv_valid_ctv2 - 1] = 0;
    sex_male_collapsed = (uintptr_t*)top_alloc(&topsize, indiv_valid_ctv2 * (sizeof(intptr_t) / 2));
    if (!sex_male_collapsed) {
      goto glm_linear_nosnp_ret_NOMEM;
    }
    collapse_copy_bitarr_incl(unfiltered_indiv_ct, sex_male, load_mask, indiv_valid_ct, sex_male_collapsed);
  }
  param_raw_ctl = (param_raw_ct + BITCT - 1) / BITCT;
  if (aligned_malloc(&active_params, param_raw_ctl * sizeof(intptr_t))) {
    goto glm_linear_nosnp_ret_NOMEM;
  }
  if (parameters_range_list_ptr->name_ct) {
    fill_ulong_zero(active_params, param_raw_ctl);
    active_params[0] = 1;
    numeric_range_list_to_bitfield(parameters_range_list_ptr, param_raw_ct, active_params, 0, 1);
    param_ct = popcount_longs(active_params, param_raw_ctl);
  } else {
    fill_all_bits(active_params, param_raw_ct);
    param_ct = param_raw_ct;
  }
  if (param_ct == 1) {
    logprint("Warning: Skipping --linear since the intercept is the only variable.\n");
    goto glm_linear_nosnp_ret_1;
  } else if (indiv_valid_ct <= param_ct) {
    logprint("Warning: Skipping --linear since # variables >= # samples.\n");
    if (pheno_nm_ct > param_ct) {
      logprint("(Check your covariates--all samples with at least one missing covariate are\nexcluded from this analysis.)\n");
    }
    goto glm_linear_nosnp_ret_1;
  }
  // parameter sequence:
  // 1. intercept
  // 2. --condition-list (1 to condition_ct in raw list)
  // 3. --covar (condition_ct + 1 to condition_ct + covar_ct in raw list)
  // 4. sex (condition_ct + covar_ct + 1 in raw list)
  for (uii = 0; uii < condition_ct; uii++) {
    if (is_set(active_params, uii + 1)) {
      slen = strlen(&(marker_ids[condition_uidxs[uii] * max_marker_id_len]));
      if (max_param_name_len <= slen) {
        max_param_name_len = slen + 1;
      }
    }
  }
  uii = condition_ct + 1;
  for (ujj = 0; ujj < covar_ct; ujj++) {
    if (is_set(active_params, ujj + uii)) {
      slen = strlen(&(covar_names[ujj * max_covar_name_len]));
      if (max_param_name_len <= slen) {
	max_param_name_len = slen + 1;
      }
    }
  }
  param_ctx = param_ct;
  if (tests_range_list_ptr->name_ct || (glm_modifier & GLM_TEST_ALL)) {
    ulii = (param_ct + (BITCT - 1)) / BITCT;
    if (aligned_malloc(&joint_test_params, ulii * sizeof(intptr_t))) {
      goto glm_linear_nosnp_ret_NOMEM;
    }
    fill_ulong_zero(joint_test_params, ulii);
    if (tests_range_list_ptr->name_ct) {
      numeric_range_list_to_bitfield(tests_range_list_ptr, param_ct - 1, joint_test_params, 1, 1);
      constraint_ct = popcount_longs(joint_test_params, ulii);
    } else {
      constraint_ct = param_ct - 1;
      fill_bits(joint_test_params, 0, constraint_ct);
    }
    if ((constraint_ct > 1) || (hide_covar && (constraint_ct == 1))) {
      // permit hide-covar + single --tests parameter combination
      slen = 8 + intlen((uint32_t)constraint_ct);
      if (max_param_name_len < slen) {
	max_param_name_len = slen;
      }
      param_ctx++;
    } else {
      constraint_ct = 0;
      aligned_free_null(&joint_test_params);
      logprint("Warning: Ignoring --tests since fewer than two parameter indices are in range.\n");
    }
  }

  ulii = condition_ct + covar_ct + 1;
  if (variation_in_sex && IS_SET(active_params, ulii)) {
    if (max_param_name_len < 4) {
      max_param_name_len = 4;
    }
  }
  wkspace_left -= topsize;
  if (wkspace_alloc_c_checked(&param_names, param_ctx * max_param_name_len) ||
      wkspace_alloc_d_checked(&covars_cov_major, param_ct * indiv_valid_ct * sizeof(double))) {
    goto glm_linear_nosnp_ret_NOMEM2;
  }
  wkspace_left += topsize;
  for (indiv_idx = 0; indiv_idx < indiv_valid_ct; indiv_idx++) {
    covars_cov_major[indiv_idx] = 1;
  }
  param_idx = 1;
  fill_double_zero(geno_map, 12);
  geno_map[0] = 1;
  geno_map[2] = 1;
  geno_map[4] = 1;
  geno_map[8] = 1;
  if (glm_modifier & GLM_CONDITION_RECESSIVE) {
    geno_map[2] = 0;
  } else if (!(glm_modifier & GLM_CONDITION_DOMINANT)) {
    geno_map[0] = 2;
    if (glm_xchr_model == 2) {
      geno_map[4] = 2;
    }
  }
  for (uii = 0; uii < condition_ct; uii++) {
    if (is_set(active_params, uii + 1)) {
      marker_uidx = condition_uidxs[uii];
      if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	goto glm_linear_nosnp_ret_READ_FAIL;
      }
      if (load_and_collapse_incl(bedfile, loadbuf_raw, unfiltered_indiv_ct, loadbuf_collapsed, indiv_valid_ct, load_mask, IS_SET(marker_reverse, marker_uidx))) {
	goto glm_linear_nosnp_ret_READ_FAIL;
      }
      chrom_idx = get_marker_chrom(chrom_info_ptr, marker_uidx);
      geno_map_ptr = geno_map;
      if (IS_SET(chrom_info_ptr->haploid_mask, chrom_idx)) {
	is_x = ((int32_t)chrom_idx == chrom_info_ptr->x_code);
	is_y = ((int32_t)chrom_idx == chrom_info_ptr->y_code);
	if (hh_exists) {
	  haploid_fix(hh_exists, indiv_include2, indiv_male_include2, indiv_valid_ct, is_x, is_y, (unsigned char*)loadbuf_collapsed);
	}
	if (!is_x) {
	  geno_map_ptr = &(geno_map[8]);
	}
      } else {
	is_x = 0;
      }
      if (!is_x) {
        glm_loadbuf_to_doubles(loadbuf_collapsed, indiv_valid_ct, &(covars_cov_major[param_idx * indiv_valid_ct]), geno_map_ptr, NULL);
      } else {
	glm_loadbuf_to_doubles_x(loadbuf_collapsed, sex_male_collapsed, indiv_valid_ct, &(covars_cov_major[param_idx * indiv_valid_ct]), geno_map_ptr, NULL);
      }
      strcpy(&(param_names[param_idx * max_param_name_len]), &(marker_ids[marker_uidx * max_marker_id_len]));
      param_idx++;
    }
  }
  // topsize = 0;
  if (constraint_ct) {
    if (wkspace_alloc_d_checked(&constraints_con_major, constraint_ct * param_ct * sizeof(double)) ||
        wkspace_alloc_d_checked(&df_df_buf, constraint_ct * constraint_ct * sizeof(double)) ||
        wkspace_alloc_d_checked(&df_buf, constraint_ct * sizeof(double))) {
      goto glm_linear_nosnp_ret_NOMEM;
    }
    if (wkspace_alloc_d_checked(&param_df_buf, constraint_ct * param_ct * sizeof(double)) ||
	wkspace_alloc_d_checked(&param_df_buf2, constraint_ct * param_ct * sizeof(double))) {
      goto glm_linear_nosnp_ret_NOMEM;
    }
    fill_double_zero(constraints_con_major, constraint_ct * param_ct);
    uljj = 0;
    for (constraint_idx = 0; constraint_idx < constraint_ct; uljj++, constraint_idx++) {
      next_set_ul_unsafe_ck(joint_test_params, &uljj);
      constraints_con_major[constraint_idx * param_ct + uljj + 1] = 1;
    }
    wptr = memcpya(&(param_names[param_ct * max_param_name_len]), (glm_modifier & GLM_TEST_ALL)? "FULL" : "USER", 4);
    *wptr++ = '_';
    wptr = uint32_write(wptr, constraint_ct);
    memcpy(wptr, "DF", 3);
  }
  if (wkspace_alloc_d_checked(&covars_indiv_major, param_ct * indiv_valid_ct * sizeof(double)) ||
      wkspace_alloc_ui_checked(&perm_2success_ct, (param_ctx - 1) * sizeof(int32_t)) ||
      wkspace_alloc_d_checked(&orig_stats, (param_ctx - 1) * sizeof(double)) ||
      wkspace_alloc_d_checked(&param_2d_buf, param_ct * param_ct * sizeof(double)) ||
      wkspace_alloc_d_checked(&param_2d_buf2, param_ct * param_ct * sizeof(double))) {
    goto glm_linear_nosnp_ret_NOMEM;
  }
  mi_buf = (MATRIX_INVERT_BUF1_TYPE*)wkspace_alloc(param_ct * sizeof(MATRIX_INVERT_BUF1_TYPE));
  if (!mi_buf) {
    goto glm_linear_nosnp_ret_NOMEM;
  }
  fill_uint_zero(perm_2success_ct, param_ctx - 1);
  indiv_uidx = 0;
  indiv_idx = 0;
  if (wkspace_alloc_d_checked(&g_pheno_d2, indiv_valid_ct * sizeof(double))) {
    goto glm_linear_nosnp_ret_NOMEM;
  }
  dptr = g_pheno_d2;
  g_pheno_sum = 0;
  g_pheno_ssq = 0;
  do {
    indiv_uidx = next_set_ul_unsafe(load_mask, indiv_uidx);
    indiv_uidx_stop = next_unset_ul(load_mask, indiv_uidx, unfiltered_indiv_ct);
    indiv_idx += indiv_uidx_stop - indiv_uidx;
    dptr2 = &(pheno_d[indiv_uidx]);
    indiv_uidx = indiv_uidx_stop;
    dptr3 = &(pheno_d[indiv_uidx_stop]);
    do {
      dxx = *dptr2++;
      *dptr++ = dxx;
      g_pheno_sum += dxx;
      g_pheno_ssq += dxx * dxx;
    } while (dptr2 < dptr3);
  } while (indiv_idx < indiv_valid_ct);
  if (g_pheno_ssq * ((double)((intptr_t)indiv_valid_ct)) == g_pheno_sum * g_pheno_sum) {
    goto glm_linear_nosnp_ret_PHENO_CONSTANT;
  }
  uii = condition_ct + 1;
  for (ujj = 0; ujj < covar_ct; ujj++) {
    if (is_set(active_params, ujj + uii)) {
      // indiv-major to covariate-major
      indiv_uidx = 0;
      dptr = &(covars_cov_major[param_idx * indiv_valid_ct]);
      dptr2 = &(covar_d[ujj]);
      for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_uidx++, indiv_idx++) {
	next_unset_ul_unsafe_ck(indiv_exclude, &indiv_uidx);
	if (IS_SET(load_mask, indiv_uidx)) {
	  *dptr++ = dptr2[indiv_idx * covar_ct];
	}
      }
      strcpy(&(param_names[param_idx * max_param_name_len]), &(covar_names[ujj * max_covar_name_len]));
      param_idx++;
    }
  }
  if (variation_in_sex && IS_SET(active_params, ulii)) {
    indiv_uidx = 0;
    dptr = &(covars_cov_major[param_idx * indiv_valid_ct]);
    dptr2 = &(dptr[indiv_valid_ct]);
    do {
      indiv_uidx = next_set_ul_unsafe(load_mask, indiv_uidx);
      indiv_uidx_stop = next_unset_ul(load_mask, indiv_uidx, unfiltered_indiv_ct);
      do {
	*dptr++ = (double)((intptr_t)is_set_ul(sex_male, indiv_uidx));
      } while (++indiv_uidx < indiv_uidx_stop);
    } while (dptr < dptr2);
    memcpy(&(param_names[param_idx * max_param_name_len]), "SEX", 4);
  }
  if (glm_modifier & GLM_STANDARD_BETA) {
    // with no SNPs, only need to do this once
    dzz = g_pheno_sum / ((double)((intptr_t)indiv_valid_ct)); // mean
    dyy = sqrt(((double)((intptr_t)(indiv_valid_ct - 1))) / (g_pheno_ssq - g_pheno_sum * dzz)); // 1/stdev
    dptr = g_pheno_d2;
    for (indiv_idx = 0; indiv_idx < indiv_valid_ct; indiv_idx++) {
      *dptr = ((*dptr) - dzz) * dyy;
      dptr++;
    }
    for (param_idx = 1; param_idx < param_ct; param_idx++) {
      dxx = 0; // sum
      dyy = 0; // ssq
      dptr = &(covars_cov_major[param_idx * indiv_valid_ct]);
      for (indiv_idx = 0; indiv_idx < indiv_valid_ct; indiv_idx++) {
        dzz = *dptr++;
        dxx += dzz;
        dyy += dzz * dzz;
      }
      dptr = &(covars_cov_major[param_idx * indiv_valid_ct]);
      dzz = dxx / ((double)((intptr_t)indiv_valid_ct));
      dyy = sqrt((dyy - dxx * dzz) / ((double)((intptr_t)(indiv_valid_ct - 1))));
      if (dyy == 0) {
	fill_double_zero(dptr, indiv_valid_ct);
      } else {
	dyy = 1.0 / dyy;
        for (indiv_idx = 0; indiv_idx < indiv_valid_ct; indiv_idx++) {
          *dptr = ((*dptr) - dzz) * dyy;
          dptr++;
        }
      }
    }
  }
  ii = glm_check_vif(glm_vif_thresh, param_ct, indiv_valid_ct, covars_cov_major, param_2d_buf, mi_buf, param_2d_buf2);
  if (ii == -1) {
    goto glm_linear_nosnp_ret_NOMEM;
  } else if (ii == 1) {
    logprint("Warning: Skipping --linear no-snp since VIF check failed.\n");
    goto glm_linear_nosnp_ret_1;
  }

  // required for multithreaded permutation generation
  g_cluster_ct = 0;
  g_pheno_nm_ct = indiv_valid_ct;

  perms_done = 0;

  if (do_perms) {
    if (!perm_batch_size) {
      // er, maybe this should be initialized in main()
      perm_batch_size = 512;
    }
    // not actually max(T), just fixed permutation count.
    if (perm_batch_size > glm_mperm_val) {
      perm_batch_size = glm_mperm_val;
    }
  } else {
    perm_batch_size = 1;
    mperm_save = 0;
  }
  ulii = (perm_batch_size + (BITCT - 1)) / BITCT;
  if (wkspace_alloc_ul_checked(&perm_fails, ulii * sizeof(intptr_t)) ||
      wkspace_alloc_d_checked(&regression_results, perm_batch_size * (param_ctx - 1) * sizeof(double)) ||
      wkspace_alloc_d_checked(&indiv_1d_buf, indiv_valid_ct * sizeof(double))) {
    goto glm_linear_nosnp_ret_NOMEM;
  }

  if (do_perms) {
    if (wkspace_alloc_ui_checked(&g_precomputed_mods, (indiv_valid_ct - 1) * sizeof(int32_t))) {
      goto glm_linear_nosnp_ret_NOMEM;
    }
    precompute_mods(indiv_valid_ct, g_precomputed_mods);
  }
  // may want to put a multiple linear regression wrapper function in
  // plink_matrix, perhaps with the PLINK 1.07 svdcmp/svbksb no-LAPACK
  // fallback

  // multiple linear regression-specific allocations and preprocessing
  dgels_m = (int32_t)((uint32_t)indiv_valid_ct);
  dgels_n = (int32_t)((uint32_t)param_ct);
  dgels_nrhs = perm_batch_size;
  dgels_ldb = dgels_m;
  if (wkspace_alloc_d_checked(&g_perm_pmajor, perm_batch_size * indiv_valid_ct * sizeof(double)) ||
      wkspace_alloc_d_checked(&dgels_a, param_ct * indiv_valid_ct * sizeof(double)) ||
      wkspace_alloc_d_checked(&dgels_b, perm_batch_size * indiv_valid_ct * sizeof(double))) {
    goto glm_linear_nosnp_ret_NOMEM;
  }
  fill_double_zero(regression_results, param_ctx - 1);
  memcpy(dgels_a, covars_cov_major, param_ct * indiv_valid_ct * sizeof(double));
  memcpy(dgels_b, g_pheno_d2, indiv_valid_ct * sizeof(double));
  dgels_(&dgels_trans, &dgels_m, &dgels_n, &dgels_nrhs, dgels_a, &dgels_m, dgels_b, &dgels_ldb, &dxx, &dgels_lwork, &dgels_info);
  if (dxx > 2147483647.0) {
    // maybe this can't actually happen, but just in case...
    // (todo: update this, and all the other matrix logic, once LAPACK can
    // handle 64-bit integers; currently that's on their wishlist)
    logprint("Error: Multiple linear regression problem too large for current LAPACK version.\n");
    retval = RET_CALC_NOT_YET_SUPPORTED;
    goto glm_linear_nosnp_ret_1;
  }
  dgels_lwork = (int32_t)dxx;
  if (wkspace_alloc_d_checked(&dgels_work, dgels_lwork * sizeof(double))) {
    goto glm_linear_nosnp_ret_NOMEM;
  }
  dgels_nrhs = 1;

  dgels_(&dgels_trans, &dgels_m, &dgels_n, &dgels_nrhs, dgels_a, &dgels_m, dgels_b, &dgels_ldb, dgels_work, &dgels_lwork, &dgels_info);
  if (dgels_info) {
    logprint("Warning: Skipping --linear no-snp since regression failed.\n");
    goto glm_linear_nosnp_ret_1;
  }

  uii = perm_batch_size / CACHELINE_INT32;
  if (!uii) {
    uii = 1;
  }
  if (max_thread_ct > uii) {
    max_thread_ct = uii;
  }
  if (cluster_starts) {
    retval = cluster_include_and_reindex(unfiltered_indiv_ct, load_mask, 0, NULL, indiv_valid_ct, 0, cluster_ct, cluster_map, cluster_starts, &cluster_ct1, &cluster_map1, &cluster_starts1, NULL, NULL);
    if (retval) {
      goto glm_linear_nosnp_ret_1;
    }
    if (cluster_ct1) {
      ulii = MAXV(cluster_ct1 + 1, param_ct);
      if (wkspace_alloc_d_checked(&cluster_param_buf, ulii * param_ct * sizeof(double)) ||
          wkspace_alloc_d_checked(&cluster_param_buf2, (cluster_ct1 + 1) * param_ct * sizeof(double)) ||
          wkspace_alloc_ui_checked(&indiv_to_cluster1, indiv_valid_ct * sizeof(int32_t))) {
	goto glm_linear_nosnp_ret_NOMEM;
      }
      fill_unfiltered_indiv_to_cluster(indiv_valid_ct, cluster_ct1, cluster_map1, cluster_starts1, indiv_to_cluster1);
      if (do_perms) {
	retval = cluster_include_and_reindex(unfiltered_indiv_ct, load_mask, 1, NULL, indiv_valid_ct, 0, cluster_ct, cluster_map, cluster_starts, &g_cluster_ct, &g_cluster_map, &g_cluster_starts, NULL, NULL);
	if (retval) {
	  goto glm_linear_nosnp_ret_1;
	}
	if (!g_cluster_ct) {
	  goto glm_linear_nosnp_ret_NO_PERMUTATION_CLUSTERS;
	}
	if (wkspace_alloc_ui_checked(&g_qassoc_cluster_thread_wkspace, max_thread_ct * ((g_cluster_ct + (CACHELINE_INT32 - 1)) / CACHELINE_INT32) * CACHELINE)) {
	  goto glm_linear_nosnp_ret_NOMEM;
	}
	if (wkspace_alloc_ui_checked(&g_indiv_to_cluster, indiv_valid_ct * sizeof(int32_t))) {
	  goto glm_linear_nosnp_ret_NOMEM;
	}
	fill_unfiltered_indiv_to_cluster(indiv_valid_ct, g_cluster_ct, g_cluster_map, g_cluster_starts, g_indiv_to_cluster);
      }
    }
  }
  if (do_perms) {
    // Note that, for now, the main nosnp regression loop is not multithreaded;
    // only the permutation generation process is.
    if (wkspace_init_sfmtp(max_thread_ct)) {
      goto glm_linear_nosnp_ret_NOMEM;
    }
  }

  transpose_copy(param_ct, indiv_valid_ct, covars_cov_major, covars_indiv_major);
  if (glm_linear_robust_cluster_covar(1, param_ct, indiv_valid_ct, 0, NULL, 0, 0, 0, covars_cov_major, covars_indiv_major, g_pheno_d2, dgels_b, param_2d_buf, mi_buf, param_2d_buf2, cluster_ct1, indiv_to_cluster1, cluster_param_buf, cluster_param_buf2, indiv_1d_buf, regression_results, constraint_ct, constraints_con_major, param_df_buf, param_df_buf2, df_df_buf, df_buf, &perm_fail_ct, perm_fails) || perm_fail_ct) {
    logprint("Warning: Skipping --linear no-snp due to multicollinearity.\n");
    goto glm_linear_nosnp_ret_1;
  }
  if (constraint_ct && (regression_results[param_ct - 1] == -9)) {
    logprint("Warning: Ignoring --tests due to regression failure.\n");
    constraint_ct = 0;
  }

  if (mperm_save) {
    // --mperm-save prevented during command-line parsing, so must be
    // --mperm-save-all
    if (wkspace_alloc_d_checked(&mperm_save_stats, glm_mperm_val * (param_ctx - 1) * sizeof(double))) {
      goto glm_linear_nosnp_ret_NOMEM;
    }
    *outname_end = '\0';
    LOGPREPRINTFWW(logbuf, "Dumping all permutation absolute t-stats to %s.[testID].mperm.dump.all.\n", outname);
    fputs(logbuf, stdout);
    if (constraint_ct) {
      logprint("(exception: chi-square values will be dumped for joint test)\n");
    }
  }
  outname_end2 = memcpyb(outname_end, ".assoc.linear", 14);
  if (fopen_checked(&outfile, outname, "w")) {
    goto glm_linear_nosnp_ret_OPEN_FAIL;
  }
  LOGPRINTFWW5("Writing linear model association results to %s ... ", outname);
  fflush(stdout);
  fputs("      TEST    NMISS       BETA ", outfile);
  if (display_ci) {
    uii = (uint32_t)((int32_t)(ci_size * 100));
    if (uii >= 10) {
      fprintf(outfile, "      SE      L%u      U%u ", uii, uii);
    } else {
      fprintf(outfile, "      SE       L%u       U%u ", uii, uii);
    }
  }
  if (fputs_checked("        STAT            P \n", outfile)) {
    goto glm_linear_nosnp_ret_WRITE_FAIL;
  }

  for (param_idx = 1; param_idx < param_ct; param_idx++) {
    dxx = dgels_b[param_idx]; // coef[p]
    se = sqrt(regression_results[param_idx - 1]);
    zval = dxx / se;
    orig_stats[param_idx - 1] = fabs(zval);
    pval = calc_tprob(zval, indiv_valid_ct - param_ct);
    if ((!hide_covar) && (pval <= pfilter)) {
      wptr = fw_strcpy(10, &(param_names[param_idx * max_param_name_len]), tbuf);
      *wptr++ = ' ';
      wptr = uint32_writew8x(wptr, (uint32_t)indiv_valid_ct, ' ');
      wptr = double_g_writewx4x(wptr, dxx, 10, ' ');
      if (display_ci) {
	dyy = ci_zt * se;
	wptr = double_g_writewx4x(wptr, se, 8, ' ');
	wptr = double_g_writewx4x(wptr, dxx - dyy, 8, ' ');
	wptr = double_g_writewx4x(wptr, dxx + dyy, 8, ' ');
      }
      wptr = double_g_writewx4x(wptr, zval, 12, ' ');
      wptr = double_g_writewx4x(wptr, pval, 12, '\n');
      if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	goto glm_linear_nosnp_ret_WRITE_FAIL;
      }
    }      
  }
  if (linear_intercept) {
    wptr = memcpya(tbuf, " INTERCEPT ", 11);
    wptr = uint32_writew8x(wptr, (uint32_t)indiv_valid_ct, ' ');
    wptr = double_g_writewx4x(wptr, dgels_b[0], 10, ' ');
    if (display_ci) {
      se = sqrt(param_2d_buf2[0]);
      zval = dgels_b[0] / se;
      dyy = ci_zt * se;
      wptr = double_g_writewx4x(wptr, se, 8, ' ');
      wptr = double_g_writewx4x(wptr, dgels_b[0] - dyy, 8, ' ');
      wptr = double_g_writewx4x(wptr, dgels_b[0] + dyy, 8, ' ');
    }
    wptr = memcpya(wptr, "          NA           NA\n", 26);
    if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
      goto glm_linear_nosnp_ret_WRITE_FAIL;
    }
  }
  if (constraint_ct) {
    dxx = regression_results[param_ct - 1];
    orig_stats[param_ct - 1] = dxx;
    pval = chiprob_p(dxx, constraint_ct);
    if (pval <= pfilter) {
      wptr = fw_strcpy(10, &(param_names[param_ct * max_param_name_len]), tbuf);
      *wptr++ = ' ';
      wptr = uint32_writew8(wptr, (uint32_t)indiv_valid_ct);
      wptr = memcpya(wptr, "         NA ", 12);
      if (display_ci) {
	wptr = memcpya(wptr, "      NA       NA       NA ", 27);
      }
      wptr = double_g_writewx4x(wptr, dxx, 12, ' ');
      wptr = double_g_writewx4x(wptr, pval, 12, '\n');
      if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	goto glm_linear_nosnp_ret_WRITE_FAIL;
      }
    }
  }
  if (fclose_null(&outfile)) {
    goto glm_linear_nosnp_ret_WRITE_FAIL;
  }
  logprint("done.\n");

  msa_ptr = mperm_save_stats;
  while (perms_done < glm_mperm_val) {
    cur_batch_size = perm_batch_size;
    if (cur_batch_size > glm_mperm_val - perms_done) {
      cur_batch_size = glm_mperm_val - perms_done;
    }
    g_perm_vec_ct = cur_batch_size;
    ulii = 0;

    if (cur_batch_size >= CACHELINE_INT32 * max_thread_ct) {
      g_assoc_thread_ct = max_thread_ct;
    } else {
      g_assoc_thread_ct = cur_batch_size / CACHELINE_INT32;
      if (!g_assoc_thread_ct) {
	g_assoc_thread_ct = 1;
      }
    }
    if (!g_cluster_ct) {
      if (spawn_threads(threads, &linear_gen_perms_thread, g_assoc_thread_ct)) {
	goto glm_linear_nosnp_ret_THREAD_CREATE_FAIL;
      }
      linear_gen_perms_thread((void*)ulii);
    } else {
      if (spawn_threads(threads, &linear_gen_cluster_perms_thread, g_assoc_thread_ct)) {
	goto glm_linear_nosnp_ret_THREAD_CREATE_FAIL;
      }
      linear_gen_cluster_perms_thread((void*)ulii);
    }
    join_threads(threads, g_assoc_thread_ct);
    dgels_nrhs = cur_batch_size;
    fill_double_zero(regression_results, (param_ctx - 1) * cur_batch_size);
    memcpy(dgels_a, covars_cov_major, param_ct * indiv_valid_ct * sizeof(double));
    memcpy(dgels_b, g_perm_pmajor, cur_batch_size * indiv_valid_ct * sizeof(double));
    dgels_(&dgels_trans, &dgels_m, &dgels_n, &dgels_nrhs, dgels_a, &dgels_m, dgels_b, &dgels_ldb, dgels_work, &dgels_lwork, &dgels_info);
    if (glm_linear_robust_cluster_covar(cur_batch_size, param_ct, indiv_valid_ct, 0, NULL, 0, 0, 0, covars_cov_major, covars_indiv_major, g_perm_pmajor, dgels_b, param_2d_buf, mi_buf, param_2d_buf2, cluster_ct1, indiv_to_cluster1, cluster_param_buf, cluster_param_buf2, indiv_1d_buf, regression_results, constraint_ct, constraints_con_major, param_df_buf, param_df_buf2, df_df_buf, df_buf, &perm_fail_ct, perm_fails)) {
      perm_fail_ct = cur_batch_size;
      fill_bits(perm_fails, 0, cur_batch_size);
    }
    perm_fail_total += perm_fail_ct;
    ulii = param_ct - 1;
    uljj = param_ctx - 1;
    for (perm_idx = 0; perm_idx < cur_batch_size; perm_idx++) {
      if (IS_SET(perm_fails, perm_idx)) {
	if (mperm_save) {
	  for (param_idx = 0; param_idx < uljj; param_idx++) {
	    *msa_ptr++ = -9;
	  }
	}
	continue;
      }
      // permutation-major regression coefficients
      dptr = &(dgels_b[perm_idx * indiv_valid_ct + 1]);
      // permutation-major variances
      dptr2 = &(regression_results[perm_idx * uljj]);
      dptr3 = orig_stats;
      for (param_idx = 0; param_idx < ulii; param_idx++) {
	dxx = *dptr++;
	se = sqrt(*dptr2++);
	zval = fabs(dxx / se);
	dyy = *dptr3++;
	if (zval > dyy + EPSILON) {
	  perm_2success_ct[param_idx] += 2;
	} else if (zval > dyy - EPSILON) {
	  perm_2success_ct[param_idx] += 1;
	}
	if (mperm_save) {
	  *msa_ptr++ = zval;
	}
      }
      if (joint_test_params) {
	if (constraint_ct) {
	  dxx = *dptr2;
	  dyy = *dptr3;
	  if (dxx > dyy + EPSILON) {
	    perm_2success_ct[ulii] += 2;
	  } else if (dxx > dyy - EPSILON) {
	    perm_2success_ct[ulii] += 1;
	  }
	  if (mperm_save) {
	    *msa_ptr++ = dxx;
	  }
	  if (dxx == -9) {
	    joint_perm_fail_extra++;
	  }
	} else if (mperm_save) {
	  msa_ptr++;
	}
      }
    }
    perms_done += cur_batch_size;
    putchar('\r');
    LOGPRINTF("%u permutation%s complete.", perms_done, (perms_done != 1)? "s" : "");
    fflush(stdout);
  }
  if (do_perms) {
    putchar('\n');
    memcpy(outname_end2, ".mperm", 7);
    if (fopen_checked(&outfile, outname, "w")) {
      goto glm_linear_nosnp_ret_OPEN_FAIL;
    }
    if (fputs_checked("      TEST         EMP1           NP \n", outfile)) {
      goto glm_linear_nosnp_ret_WRITE_FAIL;
    }
    dxx = 0.5 / ((double)((int32_t)(glm_mperm_val - perm_fail_total) + 1));
    for (param_idx = 1; param_idx < param_ct; param_idx++) {
      wptr = fw_strcpy(10, &(param_names[param_idx * max_param_name_len]), tbuf);
      *wptr++ = ' ';
      pval = ((double)(perm_2success_ct[param_idx - 1] + 2)) * dxx;
      if (pval <= pfilter) {
	if (!perm_count) {
	  wptr = double_g_writewx4(wptr, pval, 12);
	} else {
          wptr = double_g_writewx4(wptr, ((double)perm_2success_ct[param_idx - 1]) / 2.0, 12);
	}
        wptr = memseta(wptr, 32, 3);
        wptr = uint32_writew10(wptr, glm_mperm_val - perm_fail_total);
        wptr = memcpya(wptr, " \n", 2);
	if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	  goto glm_linear_nosnp_ret_WRITE_FAIL;
	}
      }
    }
    if (constraint_ct) {
      wptr = fw_strcpy(10, &(param_names[param_ct * max_param_name_len]), tbuf);
      *wptr++ = ' ';
      pval = ((double)(perm_2success_ct[param_ct - 1] + 2)) * 0.5 / ((double)((int32_t)(glm_mperm_val - perm_fail_total - joint_perm_fail_extra) + 1));
      if (pval <= pfilter) {
	if (!perm_count) {
	  wptr = double_g_writewx4(wptr, pval, 12);
	} else {
          wptr = double_g_writewx4(wptr, ((double)perm_2success_ct[param_ct - 1]) / 2.0, 12);
	}
        wptr = memseta(wptr, 32, 3);
        wptr = uint32_writew10(wptr, glm_mperm_val - perm_fail_total - joint_perm_fail_extra);
        wptr = memcpya(wptr, " \n", 2);
	if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	  goto glm_linear_nosnp_ret_WRITE_FAIL;
	}
      }
    }
    if (fclose_null(&outfile)) {
      goto glm_linear_nosnp_ret_WRITE_FAIL;
    }
    LOGPRINTFWW("Permutation test report written to %s .\n", outname);
    if (mperm_save) {
      *outname_end = '.';
      if ((param_ct != param_ctx) && constraint_ct) {
	// need to distinguish between --tests being cancelled due to inversion
	// failure, vs. normal completion
	param_ct++;
      }
      for (param_idx = 1; param_idx < param_ct; param_idx++) {
	wptr = strcpya(&(outname_end[1]), &(param_names[param_idx * max_param_name_len]));
	memcpy(wptr, ".mperm.dump.all", 17);
	if (fopen_checked(&outfile, outname, "w")) {
	  goto glm_linear_nosnp_ret_OPEN_FAIL;
	}
	ulii = param_ctx - 1;
	wptr = memcpya(tbuf, "0 ", 2);
	wptr = double_g_writex(wptr, orig_stats[param_idx - 1], '\n');
	wptr2 = &(tbuf[MAXLINELEN]);
	dptr = &(mperm_save_stats[param_idx - 1]);
	for (perm_idx = 0; perm_idx < glm_mperm_val; perm_idx++) {
	  wptr = uint32_writex(wptr, perm_idx + 1, ' ');
	  dxx = dptr[perm_idx * ulii];
	  if (dxx >= 0) {
	    wptr = double_g_writex(wptr, dxx, '\n');
	  } else {
	    wptr = memcpyl3a(wptr, "NA\n");
	  }
	  if (wptr >= wptr2) {
	    if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	      goto glm_linear_nosnp_ret_WRITE_FAIL;
	    }
	    wptr = tbuf;
	  }
	}
	if (fwrite_checkedz(tbuf, wptr - tbuf, outfile)) {
	  goto glm_linear_nosnp_ret_WRITE_FAIL;
	}
	if (fclose_null(&outfile)) {
	  goto glm_linear_nosnp_ret_WRITE_FAIL;
	}
	LOGPREPRINTFWW("%s written.\n", outname);
	logstr(logbuf);
      }
    }
  }
  while (0) {
  glm_linear_nosnp_ret_NOMEM2:
    wkspace_left += topsize;
  glm_linear_nosnp_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  glm_linear_nosnp_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  glm_linear_nosnp_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  glm_linear_nosnp_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  glm_linear_nosnp_ret_THREAD_CREATE_FAIL:
    logprint(errstr_thread_create);
    retval = RET_THREAD_CREATE_FAIL;
    break;
  glm_linear_nosnp_ret_NO_PERMUTATION_CLUSTERS:
    logprint("Error: No size 2+ clusters for permutation test.\n");
  glm_linear_nosnp_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  glm_linear_nosnp_ret_PHENO_CONSTANT:
    logprint("Warning: Skipping --linear since phenotype is constant.\n");
    break;
  }
 glm_linear_nosnp_ret_1:
  wkspace_reset(wkspace_mark);
  fclose_cond(outfile);
  aligned_free_cond(active_params);
  aligned_free_cond(joint_test_params);
  free_cond(condition_uidxs);
  return retval;
}
#endif

int32_t glm_logistic_nosnp(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uint32_t glm_modifier, double glm_vif_thresh, uint32_t glm_xchr_model, uint32_t glm_mperm_val, Range_list* parameters_range_list_ptr, Range_list* tests_range_list_ptr, double ci_size, double ci_zt, double pfilter, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t* marker_reverse, char* condition_mname, char* condition_fname, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_indiv_ct, uintptr_t indiv_ct, uintptr_t* indiv_exclude, uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, uint32_t mperm_save, uint32_t pheno_nm_ct, uintptr_t* pheno_nm, uintptr_t* pheno_c, uintptr_t covar_ct, char* covar_names, uintptr_t max_covar_name_len, uintptr_t* covar_nm, double* covar_d, uintptr_t* sex_nm, uintptr_t* sex_male, uint32_t hh_exists, uint32_t perm_batch_size, Set_info* sip) {
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + BITCT - 1) / BITCT;
  uintptr_t unfiltered_indiv_ctv2 = 2 * unfiltered_indiv_ctl;
  FILE* outfile = NULL;
  uintptr_t indiv_uidx = 0;
  uintptr_t topsize = 0;
  uintptr_t max_param_name_len = 2;
  uintptr_t param_raw_ct = 1;
  uintptr_t condition_ct = 0;
  uintptr_t constraint_ct = 0;
  uintptr_t ulii = 0;
  uint32_t cur_batch_size = 1;
  uint32_t cluster_ct1 = 0;
  uint32_t do_perms = (glm_modifier / GLM_MPERM) & 1;
  uint32_t perm_count = glm_modifier & GLM_PERM_COUNT;
  uint32_t hide_covar = glm_modifier & GLM_HIDE_COVAR;
  uint32_t report_odds = !(glm_modifier & GLM_BETA);
  uint32_t display_ci = (ci_size > 0);
  uint32_t variation_in_sex = 0; // no need to initialize if no-x-sex specified
  uint32_t perms_done = 0;
  uint32_t perm_fail_total = 0;
  uint32_t joint_perm_fail_extra = 0;
  uint32_t max_thread_ct = g_thread_ct;
  int32_t retval = 0;
  uintptr_t* loadbuf_raw = NULL;
  uintptr_t* loadbuf_collapsed = NULL;
  uintptr_t* load_mask = NULL;
  uintptr_t* sex_male_collapsed = NULL;
  uintptr_t* indiv_include2 = NULL;
  uintptr_t* indiv_male_include2 = NULL;
  uintptr_t* active_params = NULL;
  uintptr_t* joint_test_params = NULL;
  float* covars_indiv_major = NULL;
  float* cluster_param_buf = NULL;
  float* cluster_param_buf2 = NULL;
  double* mperm_save_stats = NULL;
  double* constraints_con_major = NULL;
  double* param_1d_dbuf = NULL;
  double* param_2d_dbuf = NULL;
  double* param_2d_dbuf2 = NULL;
  double* param_df_dbuf = NULL;
  double* df_df_dbuf = NULL;
  MATRIX_INVERT_BUF1_TYPE* mi_buf = NULL;
  double* df_dbuf = NULL;
  uintptr_t* perm_fails = NULL;
  uint32_t* perm_2success_ct = NULL;
  uint32_t* condition_uidxs = NULL;
  uint32_t* cluster_map1 = NULL;
  uint32_t* cluster_starts1 = NULL;
  uint32_t* indiv_to_cluster1 = NULL;
  float geno_map[12];
  float* geno_map_ptr;
  char* param_names;
  float* covars_cov_major;
  float* coef;
  float* pp;
  float* indiv_1d_buf;
  float* pheno_buf;
  float* param_1d_buf;
  float* param_1d_buf2;
  float* param_2d_buf;
  float* param_2d_buf2;
  float* regression_results;
  double* orig_stats;
  char* outname_end2;
  char* wptr;
  char* wptr2;
  float* fptr;
  float* fptr2;
  double* dptr;
  uintptr_t indiv_valid_ct;
  uintptr_t indiv_valid_cta4;
  uintptr_t indiv_valid_ctv2;
  uintptr_t indiv_uidx_stop;
  uintptr_t indiv_idx;
  uintptr_t param_raw_ctl;
  uintptr_t param_ct;
  uintptr_t param_cta4;
  uintptr_t param_ctx; // param_ct + 1 if joint test needed, param_ct otherwise
  uintptr_t param_idx;
  uintptr_t constraint_idx;
  uintptr_t uljj;
  double* msa_ptr;
  double se;
  double zval;
  double pval;
  double dxx;
  double dyy;
  uint32_t perm_idx;
  uint32_t marker_uidx;
  uint32_t chrom_idx;
  uint32_t is_x;
  uint32_t is_y;
  uint32_t uii;
  uint32_t ujj;
  uint32_t slen;
  if (hide_covar && (!tests_range_list_ptr->name_ct) && (!(glm_modifier & GLM_TEST_ALL))) {
    logprint("Error: --logistic hide-covar no-snp produces no output.\n");
    goto glm_logistic_nosnp_ret_INVALID_CMDLINE;
  }
  if (glm_init_load_mask(indiv_exclude, pheno_nm, covar_nm, indiv_ct, unfiltered_indiv_ctv2, &load_mask)) {
    goto glm_logistic_nosnp_ret_NOMEM;
  }
  indiv_valid_ct = popcount_longs(load_mask, unfiltered_indiv_ctl);
  if (condition_mname || condition_fname) {
    loadbuf_raw = (uintptr_t*)top_alloc(&topsize, unfiltered_indiv_ctv2 * sizeof(intptr_t));
    if (!loadbuf_raw) {
      goto glm_logistic_nosnp_ret_NOMEM;
    }
    loadbuf_raw[unfiltered_indiv_ctv2 - 2] = 0;
    loadbuf_raw[unfiltered_indiv_ctv2 - 1] = 0;
    ulii = topsize;

    if (hh_exists & (Y_FIX_NEEDED | NXMHH_EXISTS)) {
      indiv_include2 = (uintptr_t*)top_alloc(&topsize, unfiltered_indiv_ctv2 * sizeof(intptr_t));
      if (!indiv_include2) {
        goto glm_logistic_nosnp_ret_NOMEM;
      }
      fill_vec_55(indiv_include2, unfiltered_indiv_ct); // harmless
    }
    if (hh_exists & (XMHH_EXISTS | Y_FIX_NEEDED)) {
      indiv_male_include2 = (uintptr_t*)top_alloc(&topsize, unfiltered_indiv_ctv2 * sizeof(intptr_t));
      if (!indiv_male_include2) {
	goto glm_logistic_nosnp_ret_NOMEM;
      }
      fill_ulong_zero(indiv_male_include2, unfiltered_indiv_ctv2);
      vec_include_init(unfiltered_indiv_ct, indiv_male_include2, sex_male);
    }
    wkspace_left -= topsize;
    retval = glm_scan_conditions(condition_mname, condition_fname, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, chrom_info_ptr, hh_exists, loadbuf_raw, bedfile, bed_offset, unfiltered_indiv_ct, sex_male, load_mask, &indiv_valid_ct, &condition_ct, &condition_uidxs, indiv_include2, indiv_male_include2);
    wkspace_left += topsize;
    if (retval) {
      goto glm_logistic_nosnp_ret_1;
    }
    topsize = ulii; // deallocate temporary indiv[_male]_include2
    param_raw_ct += condition_ct;
  }
  param_raw_ct += covar_ct;
  if (glm_modifier & GLM_SEX) {
    indiv_uidx = 0;
    indiv_idx = 0;
    uii = 2;
    do {
      indiv_uidx = next_set_ul_unsafe(load_mask, indiv_uidx);
      indiv_uidx_stop = next_unset_ul(load_mask, indiv_uidx, unfiltered_indiv_ct);
      indiv_idx += indiv_uidx_stop - indiv_uidx;
      do {
	if (IS_SET(sex_nm, indiv_uidx)) {
	  ujj = is_set(sex_male, indiv_uidx);
	  if (uii == ujj) {
	    variation_in_sex = 1;
	    indiv_idx = indiv_valid_ct;
	    break;
	  }
          uii = 1 - ujj;
	}
      } while (++indiv_uidx < indiv_uidx_stop);
    } while (indiv_idx < indiv_valid_ct);
    if (variation_in_sex) {
      param_raw_ct++;
      bitfield_and(load_mask, sex_nm, unfiltered_indiv_ctl);
      indiv_valid_ct = popcount_longs(load_mask, unfiltered_indiv_ctl);
    } else {
      logprint("Warning: Ignoring --logistic 'sex' modifier since sex is invariant.\n");
    }
  }
  indiv_valid_cta4 = (indiv_valid_ct + 3) & (~3);
  indiv_valid_ctv2 = 2 * ((indiv_valid_ct + BITCT - 1) / BITCT);

  if (condition_mname || condition_fname) {
    // now that we've determined which individuals will be in the regression,
    // initialize collapsed indiv_include2, indiv_male_include2, sex_male
    if (hh_exists & (Y_FIX_NEEDED | NXMHH_EXISTS)) {
      indiv_include2 = (uintptr_t*)top_alloc(&topsize, indiv_valid_ctv2 * sizeof(intptr_t));
      fill_vec_55(indiv_include2, indiv_valid_ct);
    }
    if (hh_exists & (XMHH_EXISTS | Y_FIX_NEEDED)) {
      indiv_male_include2 = (uintptr_t*)top_alloc(&topsize, indiv_valid_ctv2 * sizeof(intptr_t));
      alloc_collapsed_haploid_filters(unfiltered_indiv_ct, indiv_valid_ct, hh_exists, 1, load_mask, sex_male, &indiv_include2, &indiv_male_include2);
    }
    loadbuf_collapsed = (uintptr_t*)top_alloc(&topsize, indiv_valid_ctv2 * sizeof(intptr_t));
    if (!loadbuf_collapsed) {
      goto glm_logistic_nosnp_ret_NOMEM;
    }
    loadbuf_collapsed[indiv_valid_ctv2 - 1] = 0;
    sex_male_collapsed = (uintptr_t*)top_alloc(&topsize, indiv_valid_ctv2 * (sizeof(intptr_t) / 2));
    if (!sex_male_collapsed) {
      goto glm_logistic_nosnp_ret_NOMEM;
    }
    collapse_copy_bitarr_incl(unfiltered_indiv_ct, sex_male, load_mask, indiv_valid_ct, sex_male_collapsed);
  }
  param_raw_ctl = (param_raw_ct + BITCT - 1) / BITCT;
  if (aligned_malloc(&active_params, param_raw_ctl * sizeof(intptr_t))) {
    goto glm_logistic_nosnp_ret_NOMEM;
  }
  if (parameters_range_list_ptr->name_ct) {
    fill_ulong_zero(active_params, param_raw_ctl);
    active_params[0] = 1;
    numeric_range_list_to_bitfield(parameters_range_list_ptr, param_raw_ct, active_params, 0, 1);
    param_ct = popcount_longs(active_params, param_raw_ctl);
  } else {
    fill_all_bits(active_params, param_raw_ct);
    param_ct = param_raw_ct;
  }
  param_cta4 = (param_ct + 3) & (~3);
  if (param_ct == 1) {
    logprint("Warning: Skipping --logistic since the intercept is the only variable.\n");
    goto glm_logistic_nosnp_ret_1;
  } else if (indiv_valid_ct <= param_ct) {
    logprint("Warning: Skipping --logistic since # variables >= # samples.\n");
    if (pheno_nm_ct > param_ct) {
      logprint("(Check your covariates--all samples with at least one missing covariate are\nexcluded from this analysis.)\n");
    }
    goto glm_logistic_nosnp_ret_1;
  }
  // parameter sequence:
  // 1. intercept
  // 2. --condition-list (1 to condition_ct in raw list)
  // 3. --covar (condition_ct + 1 to condition_ct + covar_ct in raw list)
  // 4. sex (condition_ct + covar_ct + 1 in raw list)
  for (uii = 0; uii < condition_ct; uii++) {
    if (is_set(active_params, uii + 1)) {
      slen = strlen(&(marker_ids[condition_uidxs[uii] * max_marker_id_len]));
      if (max_param_name_len <= slen) {
        max_param_name_len = slen + 1;
      }
    }
  }
  uii = condition_ct + 1;
  for (ujj = 0; ujj < covar_ct; ujj++) {
    if (is_set(active_params, ujj + uii)) {
      slen = strlen(&(covar_names[ujj * max_covar_name_len]));
      if (max_param_name_len <= slen) {
	max_param_name_len = slen + 1;
      }
    }
  }
  param_ctx = param_ct;
  if (tests_range_list_ptr->name_ct || (glm_modifier & GLM_TEST_ALL)) {
    ulii = (param_ct + (BITCT - 1)) / BITCT;
    if (aligned_malloc(&joint_test_params, ulii * sizeof(intptr_t))) {
      goto glm_logistic_nosnp_ret_NOMEM;
    }
    fill_ulong_zero(joint_test_params, ulii);
    if (tests_range_list_ptr->name_ct) {
      numeric_range_list_to_bitfield(tests_range_list_ptr, param_ct - 1, joint_test_params, 1, 1);
      constraint_ct = popcount_longs(joint_test_params, ulii);
    } else {
      constraint_ct = param_ct - 1;
      fill_bits(joint_test_params, 0, constraint_ct);
    }
    if ((constraint_ct > 1) || (hide_covar && (constraint_ct == 1))) {
      // permit hide-covar + single --tests parameter combination
      slen = 8 + intlen((uint32_t)constraint_ct);
      if (max_param_name_len < slen) {
	max_param_name_len = slen;
      }
      param_ctx++;
    } else {
      constraint_ct = 0;
      aligned_free_null(&joint_test_params);
      logprint("Warning: Ignoring --tests since fewer than two parameter indices are in range.\n");
    }
  }

  ulii = condition_ct + covar_ct + 1;
  if (variation_in_sex && IS_SET(active_params, ulii)) {
    if (max_param_name_len < 4) {
      max_param_name_len = 4;
    }
  }
  wkspace_left -= topsize;
  if (wkspace_alloc_c_checked(&param_names, param_ctx * max_param_name_len) ||
      wkspace_alloc_f_checked(&covars_cov_major, param_ct * indiv_valid_cta4 * sizeof(float))) {
    goto glm_logistic_nosnp_ret_NOMEM2;
  }
  wkspace_left += topsize;
  for (indiv_idx = 0; indiv_idx < indiv_valid_ct; indiv_idx++) {
    covars_cov_major[indiv_idx] = 1;
  }
  if (indiv_valid_ct % 4) {
    // need to zero out end of every row
    uii = indiv_valid_cta4 - indiv_valid_ct;
    for (param_idx = 0; param_idx < param_ct; param_idx++) {
      fill_float_zero(&(covars_cov_major[param_idx * indiv_valid_cta4 + indiv_valid_ct]), uii);
    }
  }
  param_idx = 1;
  fill_float_zero(geno_map, 12);
  geno_map[0] = 1;
  geno_map[2] = 1;
  geno_map[4] = 1;
  geno_map[8] = 1;
  if (glm_modifier & GLM_CONDITION_RECESSIVE) {
    geno_map[2] = 0;
  } else if (!(glm_modifier & GLM_CONDITION_DOMINANT)) {
    geno_map[0] = 2;
    if (glm_xchr_model == 2) {
      geno_map[4] = 2;
    }
  }
  for (uii = 0; uii < condition_ct; uii++) {
    if (is_set(active_params, uii + 1)) {
      marker_uidx = condition_uidxs[uii];
      if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	goto glm_logistic_nosnp_ret_READ_FAIL;
      }
      if (load_and_collapse_incl(bedfile, loadbuf_raw, unfiltered_indiv_ct, loadbuf_collapsed, indiv_valid_ct, load_mask, IS_SET(marker_reverse, marker_uidx))) {
	goto glm_logistic_nosnp_ret_READ_FAIL;
      }
      chrom_idx = get_marker_chrom(chrom_info_ptr, marker_uidx);
      geno_map_ptr = geno_map;
      if (IS_SET(chrom_info_ptr->haploid_mask, chrom_idx)) {
	is_x = ((int32_t)chrom_idx == chrom_info_ptr->x_code);
	is_y = ((int32_t)chrom_idx == chrom_info_ptr->y_code);
	if (hh_exists) {
	  haploid_fix(hh_exists, indiv_include2, indiv_male_include2, indiv_valid_ct, is_x, is_y, (unsigned char*)loadbuf_collapsed);
	}
	if (!is_x) {
	  geno_map_ptr = &(geno_map[8]);
	}
      } else {
	is_x = 0;
      }
      if (!is_x) {
        glm_loadbuf_to_floats(loadbuf_collapsed, indiv_valid_ct, &(covars_cov_major[param_idx * indiv_valid_cta4]), geno_map_ptr, NULL);
      } else {
	glm_loadbuf_to_floats_x(loadbuf_collapsed, sex_male_collapsed, indiv_valid_ct, &(covars_cov_major[param_idx * indiv_valid_cta4]), geno_map_ptr, NULL);
      }
      strcpy(&(param_names[param_idx * max_param_name_len]), &(marker_ids[marker_uidx * max_marker_id_len]));
      param_idx++;
    }
  }
  // topsize = 0;
  if (wkspace_alloc_ui_checked(&perm_2success_ct, (param_ctx - 1) * sizeof(int32_t)) ||
      wkspace_alloc_d_checked(&orig_stats, (param_ctx - 1) * sizeof(double))) {
    goto glm_logistic_nosnp_ret_NOMEM;
  }

  if (constraint_ct) {
    if (wkspace_alloc_d_checked(&constraints_con_major, constraint_ct * param_ct * sizeof(double))) {
      goto glm_logistic_nosnp_ret_NOMEM;
    }
    fill_double_zero(constraints_con_major, constraint_ct * param_ct);
    uljj = 0;
    for (constraint_idx = 0; constraint_idx < constraint_ct; uljj++, constraint_idx++) {
      next_set_ul_unsafe_ck(joint_test_params, &uljj);
      constraints_con_major[constraint_idx * param_ct + uljj + 1] = 1;
    }
    wptr = memcpya(&(param_names[param_ct * max_param_name_len]), (glm_modifier & GLM_TEST_ALL)? "FULL" : "USER", 4);
    *wptr++ = '_';
    wptr = uint32_write(wptr, constraint_ct);
    memcpy(wptr, "DF", 3);
  }
  fill_uint_zero(perm_2success_ct, param_ctx - 1);
  indiv_uidx = 0;
  indiv_idx = 0;
  uii = condition_ct + 1;
  for (ujj = 0; ujj < covar_ct; ujj++) {
    if (is_set(active_params, ujj + uii)) {
      // indiv-major to covariate-major
      indiv_uidx = 0;
      fptr = &(covars_cov_major[param_idx * indiv_valid_cta4]);
      dptr = &(covar_d[ujj]);
      for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_uidx++, indiv_idx++) {
	next_unset_ul_unsafe_ck(indiv_exclude, &indiv_uidx);
	if (IS_SET(load_mask, indiv_uidx)) {
	  *fptr++ = (float)dptr[indiv_idx * covar_ct];
	}
      }
      strcpy(&(param_names[param_idx * max_param_name_len]), &(covar_names[ujj * max_covar_name_len]));
      param_idx++;
    }
  }
  if (variation_in_sex && IS_SET(active_params, ulii)) {
    indiv_uidx = 0;
    fptr = &(covars_cov_major[param_idx * indiv_valid_cta4]);
    fptr2 = &(fptr[indiv_valid_ct]);
    do {
      indiv_uidx = next_set_ul_unsafe(load_mask, indiv_uidx);
      indiv_uidx_stop = next_unset_ul(load_mask, indiv_uidx, unfiltered_indiv_ct);
      do {
	*fptr++ = (float)((intptr_t)is_set_ul(sex_male, indiv_uidx));
      } while (++indiv_uidx < indiv_uidx_stop);
    } while (fptr < fptr2);
    memcpy(&(param_names[param_idx * max_param_name_len]), "SEX", 4);
  }
  // no more VIF check

  // required for multithreaded permutation generation
  g_cluster_ct = 0;
  g_pheno_nm_ct = indiv_valid_ct;

  perms_done = 0;

  if (do_perms) {
    if (!perm_batch_size) {
      // er, maybe this should be initialized in main()
      perm_batch_size = 512;
    }
    // not actually max(T), just fixed permutation count.
    if (perm_batch_size > glm_mperm_val) {
      perm_batch_size = glm_mperm_val;
    }
  } else {
    perm_batch_size = 1;
    mperm_save = 0;
  }
  ulii = (perm_batch_size + (BITCT - 1)) / BITCT;
  if (wkspace_alloc_f_checked(&coef, param_cta4 * perm_batch_size * sizeof(float)) ||
      wkspace_alloc_f_checked(&pp, indiv_valid_cta4 * sizeof(float)) ||
      wkspace_alloc_f_checked(&indiv_1d_buf, indiv_valid_ct * sizeof(float)) ||
      wkspace_alloc_f_checked(&pheno_buf, indiv_valid_ct * sizeof(float)) ||
      wkspace_alloc_f_checked(&param_1d_buf, param_ct * sizeof(float)) ||
      wkspace_alloc_f_checked(&param_1d_buf2, param_ct * sizeof(float)) ||
      wkspace_alloc_f_checked(&param_2d_buf, param_ct * param_cta4 * sizeof(float)) ||
      wkspace_alloc_f_checked(&param_2d_buf2, param_ct * param_cta4 * sizeof(float)) ||
      wkspace_alloc_f_checked(&regression_results, perm_batch_size * (param_ctx - 1) * sizeof(float)) ||
      wkspace_alloc_ul_checked(&perm_fails, ulii * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&g_perm_vecs, perm_batch_size * indiv_valid_ctv2 * sizeof(intptr_t))) {
    goto glm_logistic_nosnp_ret_NOMEM;
  }
  vec_collapse_init(pheno_c, unfiltered_indiv_ct, load_mask, indiv_valid_ct, g_perm_vecs);
  g_case_ct = popcount01_longs(g_perm_vecs, indiv_valid_ctv2);
  if ((!g_case_ct) || (g_case_ct == indiv_valid_ct)) {
    goto glm_logistic_nosnp_ret_PHENO_CONSTANT;
  }

  uii = perm_batch_size;
  if (max_thread_ct > uii) {
    max_thread_ct = uii;
  }
  if (cluster_starts) {
    retval = cluster_include_and_reindex(unfiltered_indiv_ct, load_mask, 0, NULL, indiv_valid_ct, 0, cluster_ct, cluster_map, cluster_starts, &cluster_ct1, &cluster_map1, &cluster_starts1, NULL, NULL);
    if (retval) {
      goto glm_logistic_nosnp_ret_1;
    }
    if (cluster_ct1) {
      ulii = MAXV(cluster_ct1 + 1, param_ct);
      if (wkspace_alloc_f_checked(&covars_indiv_major, param_ct * indiv_valid_ct * sizeof(float)) ||
          wkspace_alloc_ui_checked(&indiv_to_cluster1, indiv_valid_ct * sizeof(int32_t)) ||
          wkspace_alloc_f_checked(&cluster_param_buf, ulii * param_ct * sizeof(float)) ||
          wkspace_alloc_f_checked(&cluster_param_buf2, (cluster_ct1 + 1) * param_ct * sizeof(float))) {
	goto glm_logistic_nosnp_ret_NOMEM;
      }
      fill_unfiltered_indiv_to_cluster(indiv_valid_ct, cluster_ct1, cluster_map1, cluster_starts1, indiv_to_cluster1);
      if (do_perms) {
	retval = cluster_include_and_reindex(unfiltered_indiv_ct, load_mask, 1, pheno_c, indiv_valid_ct, 0, cluster_ct, cluster_map, cluster_starts, &g_cluster_ct, &g_cluster_map, &g_cluster_starts, &g_cluster_case_cts, &g_cluster_cc_perm_preimage);
	if (retval) {
	  goto glm_logistic_nosnp_ret_1;
	}
	if (!g_cluster_ct) {
	  goto glm_logistic_nosnp_ret_NO_PERMUTATION_CLUSTERS;
	}
	if (cluster_alloc_and_populate_magic_nums(g_cluster_ct, g_cluster_map, g_cluster_starts, &g_tot_quotients, &g_totq_magics, &g_totq_preshifts, &g_totq_postshifts, &g_totq_incrs)) {
	  goto glm_logistic_nosnp_ret_NOMEM;
	}
	if (wkspace_alloc_ui_checked(&g_indiv_to_cluster, indiv_valid_ct * sizeof(int32_t))) {
	  goto glm_logistic_nosnp_ret_NOMEM;
	}
	fill_unfiltered_indiv_to_cluster(indiv_valid_ct, g_cluster_ct, g_cluster_map, g_cluster_starts, g_indiv_to_cluster);
      }
    }
  }
  if (constraint_ct) {
    mi_buf = (MATRIX_INVERT_BUF1_TYPE*)wkspace_alloc(param_ct * sizeof(MATRIX_INVERT_BUF1_TYPE));
    if (!mi_buf) {
      goto glm_logistic_nosnp_ret_NOMEM;
    }
    if (wkspace_alloc_d_checked(&param_1d_dbuf, param_ct * sizeof(double)) ||
        wkspace_alloc_d_checked(&param_2d_dbuf, param_ct * param_ct * sizeof(double)) ||
        wkspace_alloc_d_checked(&param_2d_dbuf2, param_ct * param_ct * sizeof(double)) ||
        wkspace_alloc_d_checked(&param_df_dbuf, param_ct * constraint_ct * sizeof(double)) ||
        wkspace_alloc_d_checked(&df_df_dbuf, constraint_ct * constraint_ct * sizeof(double)) ||
        wkspace_alloc_d_checked(&df_dbuf, constraint_ct * sizeof(double))) {
      goto glm_logistic_nosnp_ret_NOMEM;
    }
  }
  if (do_perms) {
    g_tot_quotient = 0x100000000LLU / indiv_valid_ct;
    magic_num(g_tot_quotient, &g_totq_magic, &g_totq_preshift, &g_totq_postshift, &g_totq_incr);
    // Note that, for now, the main nosnp regression loop is not multithreaded;
    // only the permutation generation process is.
    if (wkspace_init_sfmtp(max_thread_ct)) {
      goto glm_logistic_nosnp_ret_NOMEM;
    }
  }

  if (covars_indiv_major) {
    transpose_copy_float(param_ct, indiv_valid_cta4, indiv_valid_ct, covars_cov_major, covars_indiv_major);
  }
  fill_float_zero(coef, param_cta4);
  if (glm_logistic_robust_cluster_covar(1, param_ct, indiv_valid_ct, 0, NULL, covars_cov_major, covars_indiv_major, g_perm_vecs, coef, pp, indiv_1d_buf, pheno_buf, param_1d_buf, param_1d_buf2, param_2d_buf, param_2d_buf2, regression_results, cluster_ct1, indiv_to_cluster1, cluster_param_buf, cluster_param_buf2, constraint_ct, constraints_con_major, param_1d_dbuf, param_2d_dbuf, param_2d_dbuf2, param_df_dbuf, df_df_dbuf, mi_buf, df_dbuf, perm_fails)) {
    logprint("Warning: Skipping --logistic no-snp due to multicollinearity.\n");
    goto glm_logistic_nosnp_ret_1;
  }
  if (constraint_ct && (regression_results[param_ct - 1] == -9)) {
    logprint("Warning: Ignoring --tests due to regression failure.\n");
    constraint_ct = 0;
  }

  if (mperm_save) {
    // --mperm-save prevented during command-line parsing, so must be
    // --mperm-save-all
    if (wkspace_alloc_d_checked(&mperm_save_stats, glm_mperm_val * (param_ctx - 1) * sizeof(double))) {
      goto glm_logistic_nosnp_ret_NOMEM;
    }
    *outname_end = '\0';
    LOGPREPRINTFWW("Dumping all permutation chi-square values to %s.[testID].mperm.dump.all.\n", outname);
    fputs(logbuf, stdout);
  }
  outname_end2 = memcpyb(outname_end, ".assoc.logistic", 16);
  if (fopen_checked(&outfile, outname, "w")) {
    goto glm_logistic_nosnp_ret_OPEN_FAIL;
  }
  LOGPRINTFWW5("Writing logistic model association results to %s ... ", outname);
  fflush(stdout);
  fprintf(outfile, "      TEST    NMISS       %s ", report_odds? "  OR" : "BETA");
  if (display_ci) {
    uii = (uint32_t)((int32_t)(ci_size * 100));
    if (uii >= 10) {
      fprintf(outfile, "      SE      L%u      U%u ", uii, uii);
    } else {
      fprintf(outfile, "      SE       L%u       U%u ", uii, uii);
    }
  }
  if (fputs_checked("        STAT            P \n", outfile)) {
    goto glm_logistic_nosnp_ret_WRITE_FAIL;
  }

  for (param_idx = 1; param_idx < param_ct; param_idx++) {
    dxx = (double)coef[param_idx];
    se = sqrt((double)regression_results[param_idx - 1]);
    zval = dxx / se;
    orig_stats[param_idx - 1] = zval * zval;
    pval = chiprob_p(zval * zval, 1);
    if ((!hide_covar) && (pval <= pfilter)) {
      wptr = fw_strcpy(10, &(param_names[param_idx * max_param_name_len]), tbuf);
      *wptr++ = ' ';
      wptr = uint32_writew8x(wptr, (uint32_t)indiv_valid_ct, ' ');
      wptr = double_g_writewx4x(wptr, report_odds? exp(dxx) : dxx, 10, ' ');
      if (display_ci) {
	dyy = ci_zt * se;
	wptr = double_g_writewx4x(wptr, se, 8, ' ');
	if (report_odds) {
	  wptr = double_g_writewx4x(wptr, exp(dxx - dyy), 8, ' ');
	  wptr = double_g_writewx4x(wptr, exp(dxx + dyy), 8, ' ');
	} else {
	  wptr = double_g_writewx4x(wptr, dxx - dyy, 8, ' ');
	  wptr = double_g_writewx4x(wptr, dxx + dyy, 8, ' ');
	}
      }
      wptr = double_g_writewx4x(wptr, zval, 12, ' ');
      wptr = double_g_writewx4x(wptr, pval, 12, '\n');
      if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	goto glm_logistic_nosnp_ret_WRITE_FAIL;
      }
    }
  }
  if (constraint_ct) {
    dxx = (double)regression_results[param_ct - 1];
    orig_stats[param_ct - 1] = dxx;
    pval = chiprob_p(dxx, constraint_ct);
    if (pval <= pfilter) {
      wptr = fw_strcpy(10, &(param_names[param_ct * max_param_name_len]), tbuf);
      *wptr++ = ' ';
      wptr = uint32_writew8(wptr, (uint32_t)indiv_valid_ct);
      wptr = memcpya(wptr, "         NA ", 12);
      if (display_ci) {
	wptr = memcpya(wptr, "      NA       NA       NA ", 27);
      }
      wptr = double_g_writewx4x(wptr, dxx, 12, ' ');
      wptr = double_g_writewx4x(wptr, pval, 12, '\n');
      if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	goto glm_logistic_nosnp_ret_WRITE_FAIL;
      }
    }
  }
  if (fclose_null(&outfile)) {
    goto glm_logistic_nosnp_ret_WRITE_FAIL;
  }
  logprint("done.\n");

  msa_ptr = mperm_save_stats;
  while (perms_done < glm_mperm_val) {
    cur_batch_size = perm_batch_size;
    if (cur_batch_size > glm_mperm_val - perms_done) {
      cur_batch_size = glm_mperm_val - perms_done;
    }
    g_perm_vec_ct = cur_batch_size;
    ulii = 0;

    if (cur_batch_size > max_thread_ct) {
      g_assoc_thread_ct = max_thread_ct;
    } else {
      g_assoc_thread_ct = g_perm_vec_ct;
    }
    if (!g_cluster_ct) {
      if (spawn_threads(threads, &logistic_gen_perms_thread, g_assoc_thread_ct)) {
	goto glm_logistic_nosnp_ret_THREAD_CREATE_FAIL;
      }
      logistic_gen_perms_thread((void*)ulii);
    } else {
      if (spawn_threads(threads, &logistic_gen_cluster_perms_thread, g_assoc_thread_ct)) {
	goto glm_logistic_nosnp_ret_THREAD_CREATE_FAIL;
      }
      logistic_gen_cluster_perms_thread((void*)ulii);
    }
    join_threads(threads, g_assoc_thread_ct);
    fill_float_zero(coef, cur_batch_size * param_cta4);
    perm_fail_total += glm_logistic_robust_cluster_covar(cur_batch_size, param_ct, indiv_valid_ct, 0, NULL, covars_cov_major, covars_indiv_major, g_perm_vecs, coef, pp, indiv_1d_buf, pheno_buf, param_1d_buf, param_1d_buf2, param_2d_buf, param_2d_buf2, regression_results, cluster_ct1, indiv_to_cluster1, cluster_param_buf, cluster_param_buf2, constraint_ct, constraints_con_major, param_1d_dbuf, param_2d_dbuf, param_2d_dbuf2, param_df_dbuf, df_df_dbuf, mi_buf, df_dbuf, perm_fails);
    ulii = param_ct - 1;
    uljj = param_ctx - 1;
    for (perm_idx = 0; perm_idx < cur_batch_size; perm_idx++) {
      if (IS_SET(perm_fails, perm_idx)) {
	if (mperm_save) {
	  for (param_idx = 0; param_idx < uljj; param_idx++) {
	    *msa_ptr++ = -9;
	  }
	}
	continue;
      }
      fptr = &(coef[perm_idx * param_ct + 1]);
      fptr2 = &(regression_results[perm_idx * uljj]);
      dptr = orig_stats;
      for (param_idx = 0; param_idx < ulii; param_idx++) {
	dxx = (double)(*fptr++);
	se = sqrt((double)(*fptr2++));
	zval = fabs(dxx / se);
	dyy = *dptr++;
	zval *= zval;
	if (zval > dyy + EPSILON) {
	  perm_2success_ct[param_idx] += 2;
	} else if (zval > dyy - EPSILON) {
	  perm_2success_ct[param_idx] += 1;
	}
	if (mperm_save) {
	  *msa_ptr++ = zval;
	}
      }
      if (joint_test_params) {
	if (constraint_ct) {
	  dxx = (double)(*fptr2);
	  dyy = *dptr;
	  if (dxx > dyy + EPSILON) {
	    perm_2success_ct[ulii] += 2;
	  } else if (dxx > dyy - EPSILON) {
	    perm_2success_ct[ulii] += 1;
	  }
	  if (mperm_save) {
	    *msa_ptr++ = dxx;
	  }
	  if (dxx == -9) {
	    joint_perm_fail_extra++;
	  }
	} else if (mperm_save) {
	  msa_ptr++;
	}
      }
    }
    perms_done += cur_batch_size;
    putchar('\r');
    LOGPRINTF("%u permutation%s complete.", perms_done, (perms_done != 1)? "s" : "");
    fflush(stdout);
  }
  if (do_perms) {
    putchar('\n');
    memcpy(outname_end2, ".mperm", 7);
    if (fopen_checked(&outfile, outname, "w")) {
      goto glm_logistic_nosnp_ret_OPEN_FAIL;
    }
    if (fputs_checked("      TEST         EMP1           NP \n", outfile)) {
      goto glm_logistic_nosnp_ret_WRITE_FAIL;
    }
    dxx = 0.5 / ((double)((int32_t)(glm_mperm_val - perm_fail_total) + 1));
    for (param_idx = 1; param_idx < param_ct; param_idx++) {
      wptr = fw_strcpy(10, &(param_names[param_idx * max_param_name_len]), tbuf);
      *wptr++ = ' ';
      pval = ((double)(perm_2success_ct[param_idx - 1] + 2)) * dxx;
      if (pval <= pfilter) {
	if (!perm_count) {
	  wptr = double_g_writewx4(wptr, pval, 12);
	} else {
          wptr = double_g_writewx4(wptr, ((double)perm_2success_ct[param_idx - 1]) / 2.0, 12);
	}
        wptr = memseta(wptr, 32, 3);
        wptr = uint32_writew10(wptr, glm_mperm_val - perm_fail_total);
        wptr = memcpya(wptr, " \n", 2);
	if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	  goto glm_logistic_nosnp_ret_WRITE_FAIL;
	}
      }
    }
    if (constraint_ct) {
      wptr = fw_strcpy(10, &(param_names[param_ct * max_param_name_len]), tbuf);
      *wptr++ = ' ';
      pval = ((double)(perm_2success_ct[param_ct - 1] + 2)) * 0.5 / ((double)((int32_t)(glm_mperm_val - perm_fail_total - joint_perm_fail_extra) + 1));
      if (pval <= pfilter) {
	if (!perm_count) {
	  wptr = double_g_writewx4(wptr, pval, 12);
	} else {
          wptr = double_g_writewx4(wptr, ((double)perm_2success_ct[param_ct - 1]) / 2.0, 12);
	}
        wptr = memseta(wptr, 32, 3);
        wptr = uint32_writew10(wptr, glm_mperm_val - perm_fail_total - joint_perm_fail_extra);
        wptr = memcpya(wptr, " \n", 2);
	if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	  goto glm_logistic_nosnp_ret_WRITE_FAIL;
	}
      }
    }
    if (fclose_null(&outfile)) {
      goto glm_logistic_nosnp_ret_WRITE_FAIL;
    }
    LOGPRINTFWW("Permutation test report written to %s .\n", outname);
    if (mperm_save) {
      *outname_end = '.';
      if ((param_ct != param_ctx) && constraint_ct) {
	// need to distinguish between --tests being cancelled due to inversion
	// failure, vs. normal completion
	param_ct++;
      }
      for (param_idx = 1; param_idx < param_ct; param_idx++) {
	wptr = strcpya(&(outname_end[1]), &(param_names[param_idx * max_param_name_len]));
	memcpy(wptr, ".mperm.dump.all", 17);
	if (fopen_checked(&outfile, outname, "w")) {
	  goto glm_logistic_nosnp_ret_OPEN_FAIL;
	}
	ulii = param_ctx - 1;
	wptr = memcpya(tbuf, "0 ", 2);
	wptr = double_g_writex(wptr, orig_stats[param_idx - 1], '\n');
	wptr2 = &(tbuf[MAXLINELEN]);
	dptr = &(mperm_save_stats[param_idx - 1]);
	for (perm_idx = 0; perm_idx < glm_mperm_val; perm_idx++) {
	  wptr = uint32_writex(wptr, perm_idx + 1, ' ');
	  dxx = dptr[perm_idx * ulii];
	  if (dxx >= 0) {
	    wptr = double_g_writex(wptr, dxx, '\n');
	  } else {
	    wptr = memcpyl3a(wptr, "NA\n");
	  }
	  if (wptr >= wptr2) {
	    if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
	      goto glm_logistic_nosnp_ret_WRITE_FAIL;
	    }
	    wptr = tbuf;
	  }
	}
	if (fwrite_checkedz(tbuf, wptr - tbuf, outfile)) {
	  goto glm_logistic_nosnp_ret_WRITE_FAIL;
	}
	if (fclose_null(&outfile)) {
	  goto glm_logistic_nosnp_ret_WRITE_FAIL;
	}
	LOGPREPRINTFWW("%s written.\n", outname);
	logstr(logbuf);
      }
    }
  }
  while (0) {
  glm_logistic_nosnp_ret_NOMEM2:
    wkspace_left += topsize;
  glm_logistic_nosnp_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  glm_logistic_nosnp_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  glm_logistic_nosnp_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  glm_logistic_nosnp_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  glm_logistic_nosnp_ret_THREAD_CREATE_FAIL:
    logprint(errstr_thread_create);
    retval = RET_THREAD_CREATE_FAIL;
    break;
  glm_logistic_nosnp_ret_NO_PERMUTATION_CLUSTERS:
    logprint("Error: No size 2+ clusters for permutation test.\n");
  glm_logistic_nosnp_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  glm_logistic_nosnp_ret_PHENO_CONSTANT:
    logprint("Warning: Skipping --logistic since phenotype is constant.\n");
    break;
  }
 glm_logistic_nosnp_ret_1:
  wkspace_reset(wkspace_mark);
  fclose_cond(outfile);
  aligned_free_cond(active_params);
  aligned_free_cond(joint_test_params);
  free_cond(condition_uidxs);
  return retval;
}

#ifndef NOLAPACK
uint32_t glm_linear_dosage(uintptr_t indiv_ct, uintptr_t* cur_indivs, uintptr_t indiv_valid_ct, uintptr_t* pheno_nm, double* pheno_d, uintptr_t* perm_fails, uintptr_t covar_ct, double* covar_d, double* cur_dosages, double* pheno_d2, double* covars_cov_major, double* covars_indiv_major, double* param_2d_buf, MATRIX_INVERT_BUF1_TYPE* mi_buf, double* param_2d_buf2, uint32_t cluster_ct1, uint32_t* indiv_to_cluster1, double* cluster_param_buf, double* cluster_param_buf2, double* regression_results, double* indiv_1d_buf, double* dgels_a, double* dgels_b, double* dgels_work, __CLPK_integer dgels_lwork, uint32_t standard_beta, double vif_thresh, double* beta_ptr, double* se_ptr, double* pval_ptr) {
  uintptr_t param_ct = covar_ct + 2;
  if (indiv_valid_ct <= param_ct) {
    return 0;
  }
  double* dptr = pheno_d2;
  __CLPK_integer dgels_m = (int32_t)((uint32_t)indiv_valid_ct);
  __CLPK_integer dgels_n = (int32_t)((uint32_t)param_ct);
  __CLPK_integer dgels_nrhs = 1;
  __CLPK_integer dgels_ldb = dgels_m;
  char dgels_trans = 'N';
  double* dptr2;
  uintptr_t indiv_uidx;
  uintptr_t indiv_idx;
  uintptr_t covar_idx;
  double dxx;
  double dyy;
  double dzz;
  __CLPK_integer dgels_info;
  uint32_t perm_fail_ct;
  int32_t ii;
  if (!standard_beta) {
    dptr2 = covars_indiv_major;
    for (indiv_uidx = 0, indiv_idx = 0; indiv_idx < indiv_valid_ct; indiv_uidx++, indiv_idx++) {
      next_set_ul_unsafe_ck(cur_indivs, &indiv_uidx);
      *dptr++ = pheno_d[indiv_uidx];
      *dptr2++ = 1.0;
      *dptr2++ = cur_dosages[indiv_uidx];
      if (covar_ct) {
	memcpy(dptr2, &(covar_d[indiv_uidx * covar_ct]), covar_ct * sizeof(double));
	dptr2 = &(dptr2[covar_ct]);
      }
    }
    transpose_copy(indiv_valid_ct, param_ct, covars_indiv_major, covars_cov_major);
  } else {
    // kind of silly that this flag makes us switch from indiv-major ->
    // cov-major to cov-major -> indiv-major initialization, but what can you
    // do.
    dxx = 0; // sum
    dyy = 0; // ssq
    for (indiv_uidx = 0, indiv_idx = 0; indiv_idx < indiv_valid_ct; indiv_uidx++, indiv_idx++) {
      next_set_ul_unsafe_ck(cur_indivs, &indiv_uidx);
      dzz = pheno_d[indiv_uidx];
      *dptr++ = dzz;
      dxx += dzz;
      dyy += dzz * dzz;
    }
    dzz = dxx / ((double)((intptr_t)indiv_valid_ct)); // mean
    dyy = sqrt(((double)((intptr_t)(indiv_valid_ct - 1))) / (dyy - dxx * dzz)); // 1/stdev
    dptr = pheno_d2;
    for (indiv_idx = 0; indiv_idx < indiv_valid_ct; indiv_idx++) {
      *dptr = ((*dptr) - dzz) * dyy;
      dptr++;
    }
    dptr = covars_cov_major;
    for (indiv_idx = 0; indiv_idx < indiv_valid_ct; indiv_idx++) {
      *dptr++ = 1;
    }
    for (indiv_uidx = 0, indiv_idx = 0; indiv_idx < indiv_valid_ct; indiv_uidx++, indiv_idx++) {
      next_set_ul_unsafe_ck(cur_indivs, &indiv_uidx);
      *dptr++ = cur_dosages[indiv_uidx];
    }
    for (covar_idx = 0; covar_idx < covar_ct; covar_idx++) {
      dxx = 0;
      dyy = 0;
      dptr2 = dptr;
      for (indiv_uidx = 0, indiv_idx = 0; indiv_idx < indiv_valid_ct; indiv_uidx++, indiv_idx++) {
        next_set_ul_unsafe_ck(cur_indivs, &indiv_uidx);
        dzz = covar_d[indiv_uidx * covar_ct];
	*dptr2++ = dzz;
	dxx += dzz;
	dyy += dzz * dzz;
      }
      dzz = dxx / ((double)((intptr_t)indiv_valid_ct));
      dyy = sqrt(((double)((intptr_t)(indiv_valid_ct - 1))) / (dyy - dxx * dzz));
      for (indiv_idx = 0; indiv_idx < indiv_valid_ct; indiv_idx++) {
	*dptr = ((*dptr) - dzz) * dyy;
	dptr++;
      }
      // "destructive" populate, we don't need covar_d again in this function
      // call
      covar_d++;
    }
    transpose_copy(param_ct, indiv_valid_ct, covars_cov_major, covars_indiv_major);
  }
  ii = glm_check_vif(vif_thresh, param_ct, indiv_valid_ct, covars_cov_major, param_2d_buf, mi_buf, param_2d_buf2);
  if (ii == -1) {
    // kludge, might be better to calculate maximum VIF check memory usage
    // and verify that much workspace is free ahead of time
    return 2;
  } else if (ii == 1) {
    return 0;
  }
  transpose_copy(param_ct, indiv_valid_ct, covars_cov_major, covars_indiv_major);
  fill_double_zero(regression_results, param_ct - 1);
  memcpy(dgels_a, covars_cov_major, param_ct * indiv_valid_ct * sizeof(double));
  memcpy(dgels_b, pheno_d2, indiv_valid_ct * sizeof(double));
  dgels_(&dgels_trans, &dgels_m, &dgels_n, &dgels_nrhs, dgels_a, &dgels_m, dgels_b, &dgels_ldb, dgels_work, &dgels_lwork, &dgels_info);
  glm_linear_robust_cluster_covar(1, param_ct, indiv_valid_ct, 0, NULL, 0, 0, 0, covars_cov_major, covars_indiv_major, pheno_d2, dgels_b, param_2d_buf, mi_buf, param_2d_buf2, cluster_ct1, indiv_to_cluster1, cluster_param_buf, cluster_param_buf2, indiv_1d_buf, regression_results, 0, NULL, NULL, NULL, NULL, NULL, &perm_fail_ct, perm_fails);
  if (perm_fail_ct) {
    return 0;
  }
  dxx = dgels_b[1];
  dyy = sqrt(regression_results[0]);
  *beta_ptr = dxx;
  *se_ptr = dyy;
  *pval_ptr = calc_tprob(dxx / dyy, indiv_valid_ct - param_ct);
  return 1;
}
#endif

uint32_t glm_logistic_dosage(uintptr_t indiv_ct, uintptr_t* cur_indivs, uintptr_t indiv_valid_ct, uintptr_t* pheno_nm, uintptr_t* pheno_c, uintptr_t* perm_vec, uintptr_t* perm_fails, uintptr_t covar_ct, float* covar_f, double* cur_dosages, float* coef, float* pp, float* pheno_buf, float* covars_cov_major, float* covars_indiv_major, float* param_1d_buf, float* param_1d_buf2, float* param_2d_buf, float* param_2d_buf2, uint32_t cluster_ct1, uint32_t* indiv_to_cluster1, float* cluster_param_buf, float* cluster_param_buf2, float* regression_results, float* indiv_1d_buf, double* beta_ptr, double* se_ptr, double* pval_ptr) {
  uintptr_t param_ct = covar_ct + 2;
  if (indiv_valid_ct <= param_ct) {
    return 0;
  }
  uintptr_t indiv_valid_cta4 = (indiv_valid_ct + 3) & (~3);
  uintptr_t indiv_valid_ctv2 = 2 * ((indiv_valid_ct + BITCT - 1) / BITCT);
  uintptr_t param_cta4 = (param_ct + 3) & (~3);
  float* fptr = covars_cov_major;
  uintptr_t case_ct;
  uintptr_t indiv_uidx;
  uintptr_t indiv_idx;
  uintptr_t covar_idx;
  double dxx;
  double dyy;
  vec_collapse_init(pheno_c, indiv_ct, cur_indivs, indiv_valid_ct, perm_vec);
  case_ct = popcount01_longs(perm_vec, indiv_valid_ctv2);
  if ((!case_ct) || (case_ct == indiv_valid_ct)) {
    return 0;
  }
  for (indiv_idx = 0; indiv_idx < indiv_valid_ct; indiv_idx++) {
    *fptr++ = 1;
  }
  for (indiv_uidx = 0, indiv_idx = 0; indiv_idx < indiv_valid_ct; indiv_uidx++, indiv_idx++) {
    next_set_ul_unsafe_ck(cur_indivs, &indiv_uidx);
    *fptr++ = (float)cur_dosages[indiv_uidx];
  }
  for (covar_idx = 0; covar_idx < covar_ct; covar_idx++) {
    for (indiv_uidx = 0, indiv_idx = 0; indiv_idx < indiv_valid_ct; indiv_uidx++, indiv_idx++) {
      next_set_ul_unsafe_ck(cur_indivs, &indiv_uidx);
      *fptr++ = covar_f[indiv_uidx * covar_ct];
    }
    covar_f++;
  }
  if (covars_indiv_major) {
    transpose_copy_float(param_ct, indiv_valid_cta4, indiv_valid_ct, covars_cov_major, covars_indiv_major);
  }
  fill_float_zero(coef, param_cta4);
  if (glm_logistic_robust_cluster_covar(1, param_ct, indiv_valid_ct, 0, NULL, covars_cov_major, covars_indiv_major, perm_vec, coef, pp, indiv_1d_buf, pheno_buf, param_1d_buf, param_1d_buf2, param_2d_buf, param_2d_buf2, regression_results, cluster_ct1, indiv_to_cluster1, cluster_param_buf, cluster_param_buf2, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, perm_fails) || perm_fails[0]) {
    return 0;
  }
  dxx = (double)coef[1];
  dyy = sqrt((double)regression_results[0]);
  *beta_ptr = dxx;
  *se_ptr = dyy;
  *pval_ptr = calc_tprob(dxx / dyy, indiv_valid_ct - param_ct);
  return 1;
}
