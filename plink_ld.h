#ifndef __PLINK_LD_H__
#define __PLINK_LD_H__

#include "plink_common.h"

#define LD_MATRIX_SQ 1
#define LD_MATRIX_SQ0 2
#define LD_MATRIX_TRI 3
#define LD_MATRIX_SHAPEMASK 3
#define LD_MATRIX_GZ 4
#define LD_MATRIX_BIN 8
#define LD_MATRIX_SPACES 0x10
#define LD_R2 0x20
#define LD_INTER_CHR 0x40
#define LD_SINGLE_PREC 0x80
#define LD_YES_REALLY 0x100
#define LD_PRUNE_PAIRWISE 0x200
#define LD_IGNORE_X 0x400
#define LD_WEIGHTED_X 0x800

typedef struct {
  uint32_t modifier;
  uint32_t prune_window_size;
  uint32_t prune_window_incr;
  uint32_t prune_window_kb;
  double prune_last_param; // VIF or r^2 threshold
  uint32_t window_size;
  uint32_t window_bp;
  double window_r2;
  char* snpstr;
  Range_list snps_rl;
} Ld_info;

void ld_init(Ld_info* ldip);

void ld_cleanup(Ld_info* ldip);

int32_t ld_prune(Ld_info* ldip, FILE* bedfile, uintptr_t bed_offset, uint32_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, Chrom_info* chrom_info_ptr, double* set_allele_freqs, uint32_t* marker_pos, uintptr_t unfiltered_indiv_ct, uintptr_t* founder_info, uintptr_t* sex_male, char* outname, char* outname_end, uint32_t hh_exists);

#endif // __PLINK_LD_H__
