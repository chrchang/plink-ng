#ifndef __WDIST_CNV_H__
#define __WDIST_CNV_H__

#include "wdist_common.h"

int32_t wdist_cnv(char* outname, char* outname_end, char* cnvname, char* mapname, char* famname, char* phenoname, char* extractname, char* excludename, char* keepname, char* removename, char* filtername, uint64_t misc_flags, Two_col_params* update_chr, Two_col_params* update_cm, Two_col_params* update_map, Two_col_params* update_name, char* update_ids_fname, char* update_parents_fname, char* update_sex_fname, char* filterval, uint32_t filter_binary, uint32_t cnv_calc_type, uint32_t min_seglen, uint32_t max_seglen, double min_score, double max_score, uint32_t min_sites, uint32_t max_sites, uint32_t intersect_filter_type, char* intersect_filter_fname, char* subset_fname, uint32_t overlap_type, double overlap_val, uint32_t freq_type, uint32_t freq_val, double freq_val2, uint32_t test_window, uint32_t segment_modifier, char* segment_spanning_fname, uint32_t indiv_mperms, uint32_t test_mperms, uint32_t test_region_mperms, uint32_t enrichment_test_mperms, uint32_t snp_window_size, char* markername_from, char* markername_to, char* markername_snp, int32_t marker_pos_start, int32_t marker_pos_end, Range_list* snps_range_list_ptr, Chrom_info* chrom_info_ptr);

int32_t wdist_gvar(char* outname, char* outname_end, char* gvarname, char* mapname, char* famname);

#endif // __WDIST_CNV_H__
