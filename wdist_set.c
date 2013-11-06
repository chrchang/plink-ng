#include "wdist_set.h"

void set_init(Set_info* sip) {
  sip->fname = NULL;
  sip->subset_fname = NULL;
  sip->merged_set_name = NULL;
  sip->genekeep_flattened = NULL;
}

void set_cleanup(Set_info* sip) {
  free_cond(sip->fname);
  free_cond(sip->subset_fname);
  free_cond(sip->merged_set_name);
  free_cond(sip->genekeep_flattened);
}

int32_t define_sets(Set_info* sip, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_exclude_ct_ptr, char* marker_ids, uintptr_t max_marker_id_len) {
  FILE* infile = NULL;
  uintptr_t topsize = 0;
  char* sorted_marker_ids = NULL;
  //char* ;
  uintptr_t marker_exclude_ct = *marker_exclude_ct_ptr;
  uintptr_t set_ct = 0;
  uint32_t make_set = sip->modifier & SET_MAKE_FROM_RANGES;
  uint32_t complement_sets = sip->modifier & SET_COMPLEMENTS;
  uint32_t c_prefix = sip->modifier & SET_C_PREFIX;
  uint32_t gene_all = sip->modifier & SET_GENE_ALL;
  int32_t retval = 0;
  uintptr_t genekeep_ct = 0;
  uintptr_t max_genekeep_len = 0;
  uintptr_t max_set_name_len = 0; // includes trailing null
  char* set_names;
  uint32_t* set_range_ptrs;
  uint32_t* set_bounds;
  uintptr_t* set_bitfield_ptrs;
  uintptr_t* include_out_of_bounds;
  if (fopen_checked(&infile, sip->fname, "r")) {
    goto define_sets_ret_NOMEM;
  }
  /*
  // 1. if --gene or --gene-all is present, perform an extra scan to pre-filter
  //    variants.
  if (gene_all || sip->genekeep_flattened) {
    if (sip->genekeep_flattened) {
      bufptr = sip->genekeep_flattened;
      if (sip->merged_set_name) {
	// degenerate case: --gene-all if merged set name present, fail
	// otherwise
        while (1) {
          ;;;
          if (!(*bufptr)) {
            goto define_sets_ret_ALL_MARKERS_EXCLUDED;
	  }
	}
      } else {
	do {
	  slen = strlen(bufptr);
	  if ((!c_prefix) || (!memcmp(bufptr, "C_", 2))) {
	    if (slen >= max_genekeep_len) {
	      max_genekeep_len = slen + 1;
	    }
	    genekeep_ct++;
	  }
	  bufptr = &(bufptr[slen + 1]);
	} while (*bufptr);
	if (!genekeep_ct) {
	  goto define_sets_ret_ALL_MARKERS_EXCLUDED;
	}
      }
    }
    if (make_set) {
    } else {
      ;;;
    }
    if (marker_exclude_ct == unfiltered_marker_ct) {
      goto define_sets_ret_ALL_MARKERS_EXCLUDED;
    }
    *marker_exclude_ct_ptr = ;
    rewind(infile);
  }
  */
  // 1. initial scan.  count number of relevant sets if --set-collapse-all is
  //    not present, determine max_name_len, determine
  while (0) {
  define_sets_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  define_sets_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  define_sets_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  define_sets_ret_ALL_MARKERS_EXCLUDED:
    logprint("Error: All variants excluded by --gene/--gene-all.\n");
    retval = RET_ALL_MARKERS_EXCLUDED;
    break;
  define_sets_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
  if (topsize) {
    wkspace_left += topsize;
  }
  return retval;
}
