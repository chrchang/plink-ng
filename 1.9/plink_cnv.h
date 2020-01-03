#ifndef __PLINK_CNV_H__
#define __PLINK_CNV_H__

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


int32_t plink_cnv(char* outname, char* outname_end, char* cnvname, char* mapname, char* famname, char* phenoname, char* keepname, char* removename, char* filtername, uint64_t misc_flags, Two_col_params* update_chr, Two_col_params* update_cm, Two_col_params* update_map, Two_col_params* update_name, char* update_ids_fname, char* update_parents_fname, char* update_sex_fname, char* filtervals_flattened, uint64_t filter_flags, uint32_t cnv_calc_type, uint32_t min_seglen, uint32_t max_seglen, double min_score, double max_score, uint32_t min_sites, uint32_t max_sites, uint32_t intersect_filter_type, char* intersect_filter_fname, char* subset_fname, uint32_t overlap_type, double overlap_val, uint32_t freq_type, uint32_t freq_val, double freq_val2, uint32_t test_window, uint32_t segment_modifier, char* segment_spanning_fname, uint32_t sample_mperms, uint32_t test_mperms, uint32_t test_region_mperms, uint32_t enrichment_test_mperms, int32_t marker_pos_start, int32_t marker_pos_end, Chrom_info* chrom_info_ptr);

int32_t plink_gvar(char* outname, char* outname_end, char* gvarname, char* mapname, char* famname);

#endif // __PLINK_CNV_H__
