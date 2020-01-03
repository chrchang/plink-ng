#ifndef __PLINK_LD_H__
#define __PLINK_LD_H__

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


#include "plink_set.h"

#define LD_MATRIX_SQ 1
#define LD_MATRIX_SQ0 2
#define LD_MATRIX_TRI 3
#define LD_MATRIX_SHAPEMASK 3
#define LD_MATRIX_BIN 4
#define LD_MATRIX_BIN4 8
#define LD_MATRIX_SPACES 0x10
#define LD_R2 0x20
#define LD_INTER_CHR 0x40
#define LD_REPORT_GZ 0x80
#define LD_INPHASE 0x100
#define LD_D 0x200
#define LD_DPRIME 0x400
#define LD_DPRIME_SIGNED 0x800
#define LD_DX (LD_D | LD_DPRIME | LD_DPRIME_SIGNED)
#define LD_WITH_FREQS 0x1000
#define LD_YES_REALLY 0x2000
#define LD_PRUNE_PAIRWISE 0x4000
#define LD_PRUNE_PAIRPHASE 0x8000
#define LD_PRUNE_KB_WINDOW 0x10000
#define LD_IGNORE_X 0x20000
#define LD_WEIGHTED_X 0x40000
#define LD_SNP_LIST_FILE 0x80000
#define LD_BLOCKS_NO_PHENO_REQ 0x100000
#define LD_BLOCKS_NO_SMALL_MAX_SPAN 0x200000
#define LD_FLIPSCAN_VERBOSE 0x400000
#define LD_SHOW_TAGS_LIST_ALL 0x800000
#define LD_SHOW_TAGS_MODE2 0x1000000

typedef struct {
  double prune_last_param; // VIF or r^2 threshold
  double window_cm;
  double window_r2;
  double blocks_min_maf;
  double blocks_inform_frac;
  double flipscan_thresh;
  double show_tags_r2;
  char* snpstr;
  char* show_tags_fname;
  Range_list snps_rl;
  uint32_t modifier;
  uint32_t prune_window_size;
  uint32_t prune_window_incr;
  uint32_t window_size;
  uint32_t window_bp;
  uint32_t blocks_max_bp;
  // need two values here to replicate > vs. >= inconsistency in Haploview
  uint32_t blocks_strong_lowci_outer;
  uint32_t blocks_strong_lowci;
  uint32_t blocks_strong_highci;
  uint32_t blocks_recomb_highci;
  uint32_t flipscan_window_size;
  uint32_t flipscan_window_bp;
  uint32_t show_tags_bp;
} Ld_info;

// fast epistasis test is really similar to LD scan so we put it in the same
// place
#define EPI_FAST 1
#define EPI_FAST_CASE_ONLY 2
#define EPI_FAST_NO_UEKI 4
#define EPI_FAST_BOOST 8
#define EPI_FAST_JOINT_EFFECTS 0x10
#define EPI_FAST_NO_P_VALUE 0x20
#define EPI_REG 0x40
#define EPI_SET_BY_SET 0x80
#define EPI_SET_BY_ALL 0x100
#define EPI_HWE_MIDP 0x200

typedef struct {
  uint32_t modifier;
  uint32_t case_only_gap;
  double epi1;
  double epi2;
  uint32_t je_cellmin;
  // twolocus() handles --ld too
  char* ld_mkr1;
  char* ld_mkr2;
  char* twolocus_mkr1;
  char* twolocus_mkr2;
  char* summary_merge_prefix;
  uint32_t summary_merge_ct;
} Epi_info;

#define CLUMP_ALLOW_OVERLAP 1
#define CLUMP_VERBOSE 2
#define CLUMP_INDEX_FIRST 4
#define CLUMP_REPLICATE 8
#define CLUMP_BEST 0x10

typedef struct {
  uint32_t modifier;
  uint32_t fname_ct;
  uint32_t bp_radius; // distance must be less than or equal to this
  uint32_t range_border;
  char* fnames_flattened;
  char* annotate_flattened;
  char* snpfield_search_order;
  char* pfield_search_order;
  char* range_fname;
  double p1;
  double p2;
  double r2;
} Clump_info;

void ld_epi_init(Ld_info* ldip, Epi_info* epi_ip, Clump_info* clump_ip);

void ld_epi_cleanup(Ld_info* ldip, Epi_info* epi_ip, Clump_info* clump_ip);

int32_t ld_prune(Ld_info* ldip, FILE* bedfile, uintptr_t bed_offset, uintptr_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, Chrom_info* chrom_info_ptr, double* set_allele_freqs, uint32_t* marker_pos, uintptr_t unfiltered_sample_ct, uintptr_t* founder_info, uintptr_t* sex_male, char* outname, char* outname_end, uint32_t hh_exists);

int32_t flipscan(Ld_info* ldip, FILE* bedfile, uintptr_t bed_offset, uintptr_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, Chrom_info* chrom_info_ptr, double* set_allele_freqs, uint32_t* marker_pos, uintptr_t unfiltered_sample_ct, uintptr_t* pheno_nm, uintptr_t* pheno_c, uintptr_t* founder_info, uintptr_t* sex_male, char* outname, char* outname_end, uint32_t hh_exists);

int32_t ld_report(pthread_t* threads, Ld_info* ldip, FILE* bedfile, uintptr_t bed_offset, uintptr_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, double* set_allele_freqs, Chrom_info* chrom_info_ptr, uint32_t* marker_pos, double* marker_cms, uintptr_t unfiltered_sample_ct, uintptr_t* founder_info, uint32_t parallel_idx, uint32_t parallel_tot, uintptr_t* sex_male, char* outname, char* outname_end, uint32_t hh_exists);

int32_t show_tags(Ld_info* ldip, FILE* bedfile, uintptr_t bed_offset, uintptr_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, uint32_t* marker_pos, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_sample_ct, uintptr_t* founder_info, uintptr_t* sex_male, char* outname, char* outname_end, uint32_t hh_exists);

int32_t haploview_blocks(Ld_info* ldip, FILE* bedfile, uintptr_t bed_offset, uintptr_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, char* marker_ids, uintptr_t max_marker_id_len, uint32_t* marker_pos, Chrom_info* chrom_info_ptr, double* set_allele_freqs, uintptr_t unfiltered_sample_ct, uintptr_t* founder_info, uintptr_t* pheno_nm, uintptr_t* sex_male, char* outname, char* outname_end, uint32_t hh_exists);

int32_t twolocus(Epi_info* epi_ip, FILE* bedfile, uintptr_t bed_offset, uintptr_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, char** marker_allele_ptrs, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, uintptr_t* pheno_nm, uint32_t pheno_nm_ct, uint32_t pheno_ctrl_ct, uintptr_t* pheno_c, uintptr_t* sex_male, char* outname, char* outname_end, uint32_t hh_exists);

int32_t epistasis_report(pthread_t* threads, Epi_info* epi_ip, FILE* bedfile, uintptr_t bed_offset, uintptr_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, uint32_t* marker_pos, uint32_t plink_maxsnp, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_sample_ct, uintptr_t* pheno_nm, uint32_t pheno_nm_ct, uint32_t ctrl_ct, uintptr_t* pheno_c, double* pheno_d, uint32_t parallel_idx, uint32_t parallel_tot, char* outname, char* outname_end, double output_min_p, double glm_vif_thresh, Set_info* sip);

int32_t indep_pairphase(Ld_info* ldip, FILE* bedfile, uintptr_t bed_offset, uintptr_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, Chrom_info* chrom_info_ptr, double* set_allele_freqs, uint32_t* marker_pos, uintptr_t unfiltered_sample_ct, uintptr_t* founder_info, uintptr_t* sex_male, char* outname, char* outname_end, uint32_t hh_exists);

int32_t epi_summary_merge(Epi_info* epi_ip, char* outname, char* outname_end);

int32_t test_mishap(FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, double output_min_p, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_reverse, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, char** marker_allele_ptrs, double min_maf, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct);

void set_test_score(uintptr_t marker_ct, double chisq_threshold, uint32_t set_max, double* chisq_arr, uint32_t** ld_map, uint32_t* cur_setdef, double* sorted_chisq_buf, uint32_t* sorted_marker_idx_buf, uint32_t* proxy_arr, uint32_t* raw_sig_ct_ptr, uint32_t* final_sig_ct_ptr, double* set_score_ptr);

int32_t set_test_common_init(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude_orig, uintptr_t marker_ct_orig, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t* marker_reverse, double* orig_chisq, Set_info* sip, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_sample_ct, uintptr_t* sex_male, uintptr_t* founder_pnm, uint32_t ld_ignore_x, uint32_t hh_exists, const char* flag_descrip, uintptr_t* marker_ct_ptr, uintptr_t** marker_exclude_ptr, uintptr_t** set_incl_ptr, uint32_t** marker_idx_to_uidx_ptr, uint32_t*** setdefs_ptr, uintptr_t* set_ct_ptr, uint32_t* max_sigset_size_ptr, uint32_t*** ld_map_ptr, double* chisq_threshold_ptr, double** orig_set_scores_ptr, double** sorted_chisq_buf_ptr, uint32_t** sorted_marker_idx_buf_ptr, uint32_t** proxy_arr_ptr, uintptr_t** perm_adapt_set_unstopped_ptr, uint32_t** perm_2success_ct_ptr, uint32_t** perm_attempt_ct_ptr, uintptr_t** unstopped_markers_ptr);

void compute_set_scores(uintptr_t marker_ct, uintptr_t perm_vec_ct, uintptr_t set_ct, double* chisq_matrix, double* orig_set_scores, double* sorted_chisq_buf, uint32_t* sorted_marker_idx_buf, uint32_t* proxy_arr, uint32_t** setdefs, uint32_t** ld_map, Aperm_info* apip, double chisq_threshold, double adaptive_ci_zt, uint32_t first_adapt_check, uint32_t perms_done, uint32_t set_max, uintptr_t* perm_adapt_set_unstopped, uint32_t* perm_2success_ct, uint32_t* perm_attempt_ct);

int32_t write_set_test_results(char* outname, char* outname_end2, Set_info* sip, uint32_t** ld_map, uint32_t** setdefs, uintptr_t* set_incl, uintptr_t set_ct, uintptr_t marker_ct_orig, uintptr_t marker_ct, uint32_t* marker_idx_to_uidx, char* marker_ids, uintptr_t max_marker_id_len, uint32_t* perm_2success_ct, uint32_t* perm_attempt_ct, uint32_t mtest_adjust, uint32_t perm_count, double pfilter, double output_min_p, double chisq_threshold, double* orig_stats, double* sorted_chisq_buf, uint32_t* sorted_marker_idx_buf, uint32_t* proxy_arr);

int32_t clump_reports(FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, uint32_t* marker_pos, char** marker_allele_ptrs, uintptr_t* marker_reverse, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_sample_ct, uintptr_t* founder_info, Clump_info* clump_ip, uintptr_t* sex_male, uint32_t hh_exists);

#endif // __PLINK_LD_H__
