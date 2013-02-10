#include "wdist_common.h"

double fisher22(uint32_t* cts) {
  // Basic 2x2 Fisher p-value calculation.
  uint32_t ma = cts[0];
  uint32_t mb = cts[1];
  uint32_t mc = cts[2];
  uint32_t md = cts[3];
  double tprob = 1;
  double cur_prob = 1;
  double cprob = 0;
  double tail_stop;
  uint32_t uii;
  int64_t cur_a;
  int64_t cur_b;
  int64_t cur_c;
  int64_t cur_d;
  // Input:
  //   [ ma  mb ]
  //   [ mc  md ]
  //
  // Ensure we are left of the distribution center, ma <= md, and mb <= mc.
  if (mb > mc) {
    uii = mb;
    mb = mc;
    mc = uii;
  }
  if (ma > md) {
    uii = ma;
    ma = md;
    md = uii;
  }
  if ((((int64_t)ma) * md) > (((int64_t)mb) * mc)) {
    uii = ma;
    ma = mb;
    mb = uii;
    uii = mc;
    mc = md;
    md = uii;
  }
  cur_a = ma;
  cur_b = mb;
  cur_c = mc;
  cur_d = md;
  while (cur_b) {
    cur_prob *= ((double)((cur_b--) * (cur_c--))) / ((double)((++cur_a) * (++cur_d)));
    if (cur_prob == INFINITY) {
      return 0;
    }
    if (cur_prob < 1 + SMALL_EPSILON) {
      tprob += cur_prob;
      break;
    }
    cprob += cur_prob;
  }
  if (cprob == 0) {
    return 1;
  }
  if (cur_b) {
    tail_stop = tprob * DOUBLE_PREC_LIMIT;
    do {
      cur_prob *= ((double)((cur_b--) * (cur_c--))) / ((double)((++cur_a) * (++cur_d)));
      if (cur_prob < tail_stop) {
	break;
      }
      tprob += cur_prob;
    } while (cur_b);
  }
  if (ma) {
    cur_a = ma;
    cur_b = mb;
    cur_c = mc;
    cur_d = md;
    cur_prob = 1;
    tail_stop = tprob * DOUBLE_PREC_LIMIT;
    do {
      cur_prob *= ((double)((cur_a--) * (cur_d--))) / ((double)((++cur_b) * (++cur_c)));
      if (cur_prob < tail_stop) {
        return tprob / (cprob + tprob);
      }
      tprob += cur_prob;
    } while (cur_a);
  }
  return tprob / (cprob + tprob);
}

int32_t model_assoc(pthread_t* threads, FILE* bedfile, int32_t bed_offset, char* outname, char* outname_end, int32_t calculation_type, uint32_t model_modifier, uint32_t model_cell_ct, uint32_t model_mperm_val, double ci_size, double pfilter, uint32_t mtest_adjust, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uint32_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, char* marker_alleles, uintptr_t max_marker_allele_len, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude) {
  unsigned char* wkspace_mark = wkspace_base;
  int32_t retval = 0;
  FILE* outfile = NULL;

  logprint("Error: --assoc and --model are currently under development.\n");
  retval = RET_CALC_NOT_YET_SUPPORTED;
  goto model_assoc_ret_1;

  while (0) {
  }
 model_assoc_ret_1:
  wkspace_reset(wkspace_mark);
  fclose_cond(outfile);
  return retval;
}
