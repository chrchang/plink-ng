#ifndef __WDIST_CALC_H__
#define __WDIST_CALC_H__

void update_rel_ibc(double* rel_ibc, uintptr_t* geno, double* set_allele_freqs, int ibc_type, uint32_t indiv_ct);

void update_rel_f_ibc(float* rel_ibc, uintptr_t* geno, float* set_allele_freqs, int ibc_type, uint32_t indiv_ct);

#endif // __WDIST_CALC_H__
