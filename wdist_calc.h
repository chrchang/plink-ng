#ifndef __WDIST_CALC_H__
#define __WDIST_CALC_H__

void update_rel_ibc(double* rel_ibc, unsigned long* geno, double* set_allele_freqs, int ibc_type, unsigned int indiv_ct);

void update_rel_f_ibc(float* rel_ibc, unsigned long* geno, float* set_allele_freqs, int ibc_type, unsigned int indiv_ct);

#endif // __WDIST_CALC_H__x
