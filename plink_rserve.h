#ifndef __PLINK_RSERVE_H__
#define __PLINK_RSERVE_H__

int32_t rserve_call(char* rplugin_fname, uint32_t rplugin_port, FILE* bedfile, uintptr_t bed_offset, uintptr_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, uint32_t* marker_pos, uint32_t plink_maxsnp, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_sample_ct, uintptr_t* pheno_nm, uint32_t pheno_nm_ct, uintptr_t* pheno_c, double* pheno_d, char* outname, char* outname_end);

#endif // __PLINK_RSERVE_H__
