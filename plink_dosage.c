#include "plink_common.h"

#include "plink_cluster.h"
#include "plink_data.h"
#include "plink_dosage.h"
#include "plink_filter.h"
#include "plink_misc.h"

void dosage_init(Dosage_info* doip) {
  doip->fname = NULL;
  doip->modifier = 0;
  doip->skip0 = 0;
  doip->skip1 = 0;
  doip->skip2 = 0;
  doip->format = 2;
}

void dosage_cleanup(Dosage_info* doip) {
  free_cond(doip->fname);
}

typedef struct ll_ctstr_entry_struct {
  struct ll_ctstr_entry_struct* next;
  uint32_t ct;
  char ss[];
} Ll_ctstr_entry;

int32_t plink1_dosage(Dosage_info* doip, char* famname, char* mapname, char* outname, char* outname_end, char* phenoname, char* extractname, char* excludename, char* keepname, char* removename, char* keepfamname, char* removefamname, char* filtername, char* makepheno_str, char* phenoname_str, char* covar_fname, Two_col_params* qual_filter, Two_col_params* update_map, Two_col_params* update_name, char* update_ids_fname, char* update_parents_fname, char* update_sex_fname, char* filtervals_flattened, char* filter_attrib_fname, char* filter_attrib_liststr, char* filter_attrib_indiv_fname, char* filter_attrib_indiv_liststr, double qual_min_thresh, double qual_max_thresh, double thin_keep_prob, uint32_t thin_keep_ct, uint32_t min_bp_space, uint32_t mfilter_col, uint32_t fam_cols, int32_t missing_pheno, uint32_t mpheno_col, uint32_t pheno_modifier, Chrom_info* chrom_info_ptr, double tail_bottom, double tail_top, uint64_t misc_flags, uint64_t filter_flags, uint32_t sex_missing_pheno, uint32_t update_sex_col, Cluster_info* cluster_ptr, int32_t marker_pos_start, int32_t marker_pos_end, int32_t snp_window_size, char* markername_from, char* markername_to, char* markername_snp, Range_list* snps_range_list_ptr, uint32_t covar_modifier, Range_list* covar_range_list_ptr, uint32_t mwithin_col, uint32_t glm_modifier, double glm_vif_thresh, uint32_t glm_xchr_model, Range_list* parameters_range_list_ptr) {
  logprint("Error: --dosage is currently under development.\n");
  return RET_CALC_NOT_YET_SUPPORTED;
  /*
  // sucks to duplicate so much, but this code will be thrown out later so
  // there's no long-term maintenance problem
  FILE* phenofile = NULL;
  FILE* infile = NULL;
  FILE* outfile = NULL;
  gzFile* gz_infiles = NULL;
  gzFile gz_outfile = NULL;
  char* marker_ids = NULL;
  char* person_ids = NULL;
  char* paternal_ids = NULL;
  char* maternal_ids = NULL;
  char* cluster_ids = NULL;
  char* covar_names = NULL;
  char* sorted_marker_ids = NULL;
  char* sorted_person_ids = NULL;
  char* sep_fnames = NULL;
  char* cur_marker_id_buf = NULL;
  char* cur_a1 = NULL;
  char* cur_a2 = NULL;
  uintptr_t* marker_exclude = NULL;
  uintptr_t* indiv_exclude = NULL;
  uintptr_t* sex_nm = NULL;
  uintptr_t* sex_male = NULL;
  uintptr_t* pheno_nm = NULL;
  uintptr_t* pheno_c = NULL;
  uintptr_t* founder_info = NULL;
  uintptr_t* covar_nm = NULL;
  double* pheno_d = NULL;
  double* covar_d = NULL;
  Ll_ctstr_entry** htable = NULL;
  uint32_t* marker_pos = NULL;
  uint32_t* cluster_map = NULL;
  uint32_t* cluster_starts = NULL;
  uint32_t* marker_id_map = NULL;
  uint32_t* person_id_map = NULL;
  uint32_t* batch_sizes = NULL;
  uint32_t* uiptr = NULL;
  uint32_t* uiptr2 = NULL;
  uint32_t* uiptr3 = NULL;
  uintptr_t topsize = 0;
  uintptr_t unfiltered_marker_ct = 0;
  uintptr_t marker_exclude_ct = 0;
  uintptr_t max_marker_id_len = 0;
  uintptr_t unfiltered_indiv_ct = 0;
  uintptr_t indiv_exclude_ct = 0;
  uintptr_t max_person_id_len = 4;
  uintptr_t max_paternal_id_len = 2;
  uintptr_t max_maternal_id_len = 2;
  uintptr_t cluster_ct = 0;
  uintptr_t max_cluster_id_len = 2;
  uintptr_t covar_ct = 0;
  uintptr_t max_covar_name_len = 0;
  uintptr_t max_fn_len = 0;
  uintptr_t max_sepheader_len = 0;
  uintptr_t file_idx_start = 0;
  uintptr_t distinct_id_ct = 0; // occur mode
  uintptr_t max_occur_id_len = 0;
  uintptr_t marker_idx = 0;
  uintptr_t ulii = 0;
  double missing_phenod = (double)missing_pheno;
  uint32_t load_map = (mapname[0] != '\0');
  uint32_t do_glm = (doip->modifier / DOSAGE_GLM) & 1;
  uint32_t count_occur = doip->modifier & DOSAGE_OCCUR;
  uint32_t sepheader = (doip->modifier / DOSAGE_SEPHEADER) & 1;
  uint32_t noheader = doip->modifier & DOSAGE_NOHEADER;
  uint32_t output_gz = doip->modifier & DOSAGE_ZOUT;
  uint32_t dose1 = doip->modifier & DOSAGE_DOSE1;
  uint32_t skip0 = doip->skip0;
  uint32_t skip1p1 = doip->skip1 + 1;
  uint32_t skip2 = doip->skip2;
  uint32_t format_val = doip->format;
  uint32_t map_cols = 3;
  uint32_t map_is_unsorted = 0;
  uint32_t affection = 0;
  uint32_t plink_maxsnp = 0;
  uint32_t infile_ct = 0;
  uint32_t pheno_ctrl_ct = 0;
  uint32_t batch_ct = 1;
  uint32_t max_batch_size = 1;
  uint32_t cur_marker_id_len = 0;
  int32_t retval = 0;
  char* missing_mid_templates[2];
  unsigned char* wkspace_mark;
  char* fnames;
  char* loadbuf;
  char* a1_ptr;
  char* a2_ptr;
  char* bufptr;
  char* bufptr2;
  char* bufptr3;
  char* bufptr4;
  char* bufptr5;
  char* bufptr6;
  uintptr_t* line_idx_arr;
  uintptr_t* cur_indivs;
  double* cur_dosages;
  Ll_ctstr_entry** ll_pptr;
  Ll_ctstr_entry* ll_ptr;
  uint32_t* file_icts;
  uint32_t* read_idx_to_indiv_idx;
  uint32_t* skip_vals;
  uintptr_t marker_ct;
  uintptr_t unfiltered_indiv_ctl;
  uintptr_t indiv_ct;
  uintptr_t indiv_ctl;
  uintptr_t line_idx;
  uintptr_t file_idx;
  uintptr_t indiv_idx;
  uintptr_t cur_batch_size;
  uintptr_t loadbuf_size;
  uintptr_t uljj;
  double dxx;
  double dyy;
  uint32_t gender_unk_ct;
  uint32_t pheno_nm_ct;
  uint32_t batch_idx;
  uint32_t read_idx_start;
  uint32_t read_idx;
  uint32_t slen;
  uint32_t uii;
  uint32_t ujj;
  int32_t ii;
  missing_mid_templates[0] = NULL;
  missing_mid_templates[1] = NULL;
  if (load_map) {
    retval = load_bim(mapname, &map_cols, &unfiltered_marker_ct, &marker_exclude_ct, &max_marker_id_len, &marker_exclude, NULL, NULL, &ulii, &marker_ids, missing_mid_templates, NULL, chrom_info_ptr, NULL, &marker_pos, misc_flags, filter_flags, marker_pos_start, marker_pos_end, snp_window_size, markername_from, markername_to, markername_snp, snps_range_list_ptr, &map_is_unsorted, do_glm || min_bp_space || (misc_flags & (MISC_EXTRACT_RANGE | MISC_EXCLUDE_RANGE)), 0, 0, NULL, ".map file");
    if (retval) {
      goto plink1_dosage_ret_1;
    }
    if (map_is_unsorted & UNSORTED_SPLIT_CHROM) {
      logprint("Error: .map file has a split chromosome.\n");
      goto plink1_dosage_ret_INVALID_FORMAT;
    }
  }
  uii = fam_cols & FAM_COL_6;
  if (uii && phenoname) {
    uii = (pheno_modifier & PHENO_MERGE) && (!makepheno_str);
  }
  if (!uii) {
    pheno_modifier &= ~PHENO_MERGE;
  }
  if (update_ids_fname) {
    ulii = 0;
    retval = scan_max_fam_indiv_strlen(update_ids_fname, 3, &max_person_id_len);
    if (retval) {
      goto plink1_dosage_ret_1;
    }
  } else if (update_parents_fname) {
    retval = scan_max_strlen(update_parents_fname, 3, 4, 0, '\0', &max_paternal_id_len, &max_maternal_id_len);
    if (retval) {
      goto plink1_dosage_ret_1;
    }
  }
  retval = load_fam(famname, fam_cols, uii, missing_pheno, (misc_flags / MISC_AFFECTION_01) & 1, &unfiltered_indiv_ct, &person_ids, &max_person_id_len, &paternal_ids, &max_paternal_id_len, &maternal_ids, &max_maternal_id_len, &sex_nm, &sex_male, &affection, &pheno_nm, &pheno_c, &pheno_d, &founder_info, &indiv_exclude);
  if (retval) {
    goto plink1_dosage_ret_1;
  }
  unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  if (misc_flags & MISC_MAKE_FOUNDERS_FIRST) {
    if (make_founders(unfiltered_indiv_ct, unfiltered_indiv_ct, person_ids, max_person_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, (misc_flags / MISC_MAKE_FOUNDERS_REQUIRE_2_MISSING) & 1, indiv_exclude, founder_info)) {
      goto plink1_dosage_ret_NOMEM;
    }
  }
  count_genders(sex_nm, sex_male, unfiltered_indiv_ct, indiv_exclude, &uii, &ujj, &gender_unk_ct);
  marker_ct = unfiltered_marker_ct - marker_exclude_ct;
  if (gender_unk_ct) {
    LOGPRINTF("%" PRIuPTR " %s (%u male%s, %u female%s, %u ambiguous) loaded from .fam.\n", unfiltered_indiv_ct, species_str(unfiltered_indiv_ct), uii, (uii == 1)? "" : "s", ujj, (ujj == 1)? "" : "s", gender_unk_ct);
    retval = write_nosex(outname, outname_end, unfiltered_indiv_ct, indiv_exclude, sex_nm, gender_unk_ct, person_ids, max_person_id_len);
    if (retval) {
      goto plink1_dosage_ret_1;
    }
  } else {
    LOGPRINTF("%" PRIuPTR " %s (%d male%s, %d female%s) loaded from .fam.\n", unfiltered_indiv_ct, species_str(unfiltered_indiv_ct), uii, (uii == 1)? "" : "s", ujj, (ujj == 1)? "" : "s");
  }
  uii = popcount_longs(pheno_nm, unfiltered_indiv_ctl);
  if (uii) {
    LOGPRINTF("%u phenotype value%s loaded from .fam.\n", uii, (uii == 1)? "" : "s");
  }
  if (phenoname && fopen_checked(&phenofile, phenoname, "r")) {
    goto plink1_dosage_ret_OPEN_FAIL;
  }
  if (phenofile || update_ids_fname || update_parents_fname || update_sex_fname || (filter_flags & FILTER_TAIL_PHENO)) {
    wkspace_mark = wkspace_base;
    retval = sort_item_ids(&sorted_person_ids, &person_id_map, unfiltered_indiv_ct, indiv_exclude, 0, person_ids, max_person_id_len, 0, 0, strcmp_deref);
    if (retval) {
      goto plink1_dosage_ret_1;
    }

    if (makepheno_str) {
      retval = makepheno_load(phenofile, makepheno_str, unfiltered_indiv_ct, sorted_person_ids, max_person_id_len, person_id_map, pheno_nm, &pheno_c);
      if (retval) {
        goto plink1_dosage_ret_1;
      }
    } else if (phenofile) {
      retval = load_pheno(phenofile, unfiltered_indiv_ct, 0, sorted_person_ids, max_person_id_len, person_id_map, missing_pheno, (misc_flags / MISC_AFFECTION_01) & 1, mpheno_col, phenoname_str, pheno_nm, &pheno_c, &pheno_d, NULL, 0);
      if (retval) {
	if (retval == LOAD_PHENO_LAST_COL) {
	  logprintb();
	  retval = RET_INVALID_FORMAT;
	  wkspace_reset(wkspace_mark);
	}
        goto plink1_dosage_ret_1;
      }
    }
    if (filter_flags & FILTER_TAIL_PHENO) {
      retval = convert_tail_pheno(unfiltered_indiv_ct, pheno_nm, &pheno_c, &pheno_d, tail_bottom, tail_top, missing_phenod);
      if (retval) {
        goto plink1_dosage_ret_1;
      }
    }
    wkspace_reset(wkspace_mark);
  }
  if (load_map) {
    uii = update_map || update_name || filter_attrib_fname || qual_filter;
    if (uii || extractname || excludename) {
      wkspace_mark = wkspace_base;
      retval = sort_item_ids(&sorted_marker_ids, &marker_id_map, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_ids, max_marker_id_len, !uii, 0, strcmp_deref);
      if (retval) {
	goto plink1_dosage_ret_1;
      }
      ulii = unfiltered_marker_ct - marker_exclude_ct;

      if (update_map) {
	retval = update_marker_pos(update_map, sorted_marker_ids, ulii, max_marker_id_len, marker_id_map, marker_exclude, &marker_exclude_ct, marker_pos, &map_is_unsorted, chrom_info_ptr);
      } else if (update_name) {
	retval = update_marker_names(update_name, sorted_marker_ids, ulii, max_marker_id_len, marker_id_map, marker_ids);
	if (retval) {
	  goto plink1_dosage_ret_1;
	}
	if (extractname || excludename) {
	  wkspace_reset(wkspace_mark);
	  retval = sort_item_ids(&sorted_marker_ids, &marker_id_map, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_ids, max_marker_id_len, 0, 0, strcmp_deref);
	  if (retval) {
	    goto plink1_dosage_ret_1;
	  }
	  ulii = unfiltered_marker_ct - marker_exclude_ct;
	}
      }
      if (extractname) {
	if (!(misc_flags & MISC_EXTRACT_RANGE)) {
	  retval = include_or_exclude(extractname, sorted_marker_ids, ulii, max_marker_id_len, marker_id_map, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, 0);
	  if (retval) {
	    goto plink1_dosage_ret_1;
	  }
	} else {
	  if (map_is_unsorted & UNSORTED_BP) {
	    logprint("Error: '--extract range' requires a sorted .bim.  Retry this command after\nusing --make-bed to sort your data.\n");
	    goto plink1_dosage_ret_INVALID_CMDLINE;
	  }
	  retval = extract_exclude_range(extractname, marker_pos, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, 0, chrom_info_ptr);
	  if (retval) {
	    goto plink1_dosage_ret_1;
	  }
	  uljj = unfiltered_marker_ct - marker_exclude_ct;
	  LOGPRINTF("--extract range: %" PRIuPTR " variant%s remaining.\n", uljj, (uljj == 1)? "" : "s");
	}
      }
      if (excludename) {
	if (!(misc_flags & MISC_EXCLUDE_RANGE)) {
	  retval = include_or_exclude(excludename, sorted_marker_ids, ulii, max_marker_id_len, marker_id_map, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, 1);
	  if (retval) {
	    goto plink1_dosage_ret_1;
	  }
	} else {
	  if (map_is_unsorted & UNSORTED_BP) {
	    logprint("Error: '--exclude range' requires a sorted .bim.  Retry this command after\nusing --make-bed to sort your data.\n");
	    goto plink1_dosage_ret_INVALID_CMDLINE;
	  }
	  retval = extract_exclude_range(excludename, marker_pos, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, 1, chrom_info_ptr);
	  if (retval) {
	    goto plink1_dosage_ret_1;
	  }
	  uljj = unfiltered_marker_ct - marker_exclude_ct;
	  LOGPRINTF("--exclude range: %" PRIuPTR " variant%s remaining.\n", uljj, (uljj == 1)? "" : "s");
	}
      }
      if (filter_attrib_fname) {
	retval = filter_attrib(filter_attrib_fname, filter_attrib_liststr, sorted_marker_ids, ulii, max_marker_id_len, marker_id_map, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, 0);
	if (retval) {
	  goto plink1_dosage_ret_1;
	}
      }
      if (qual_filter) {
	retval = filter_qual_scores(qual_filter, qual_min_thresh, qual_max_thresh, sorted_marker_ids, ulii, max_marker_id_len, marker_id_map, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct);
	if (retval) {
	  goto plink1_dosage_ret_1;
	}
      }
      wkspace_reset(wkspace_mark);
    }
    if (thin_keep_prob != 1.0) {
      if (random_thin_markers(thin_keep_prob, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct)) {
	goto plink1_dosage_ret_ALL_MARKERS_EXCLUDED;
      }
    } else if (thin_keep_ct) {
      retval = random_thin_markers_ct(thin_keep_ct, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct);
      if (retval) {
	goto plink1_dosage_ret_1;
      }
    }
  }
  if (update_ids_fname || update_parents_fname || update_sex_fname || keepname || keepfamname || removename || removefamname || filter_attrib_indiv_fname || filtername) {
    wkspace_mark = wkspace_base;
    retval = sort_item_ids(&sorted_person_ids, &person_id_map, unfiltered_indiv_ct, indiv_exclude, indiv_exclude_ct, person_ids, max_person_id_len, 0, 0, strcmp_deref);
    if (retval) {
      goto plink1_dosage_ret_1;
    }
    ulii = unfiltered_indiv_ct - indiv_exclude_ct;
    if (update_ids_fname) {
      retval = update_indiv_ids(update_ids_fname, sorted_person_ids, ulii, max_person_id_len, person_id_map, person_ids);
      if (retval) {
	goto plink1_dosage_ret_1;
      }
    } else {
      if (update_parents_fname) {
	retval = update_indiv_parents(update_parents_fname, sorted_person_ids, ulii, max_person_id_len, person_id_map, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, founder_info);
	if (retval) {
	  goto plink1_dosage_ret_1;
	}
      }
      if (update_sex_fname) {
        retval = update_indiv_sexes(update_sex_fname, update_sex_col, sorted_person_ids, ulii, max_person_id_len, person_id_map, sex_nm, sex_male);
	if (retval) {
	  goto plink1_dosage_ret_1;
	}
      }
    }
    if (keepfamname) {
      retval = include_or_exclude(keepfamname, sorted_person_ids, ulii, max_person_id_len, person_id_map, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, 6);
      if (retval) {
	goto plink1_dosage_ret_1;
      }
    }
    if (keepname) {
      retval = include_or_exclude(keepname, sorted_person_ids, ulii, max_person_id_len, person_id_map, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, 2);
      if (retval) {
	goto plink1_dosage_ret_1;
      }
    }
    if (removefamname) {
      retval = include_or_exclude(removefamname, sorted_person_ids, ulii, max_person_id_len, person_id_map, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, 7);
      if (retval) {
	goto plink1_dosage_ret_1;
      }
    }
    if (removename) {
      retval = include_or_exclude(removename, sorted_person_ids, ulii, max_person_id_len, person_id_map, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, 3);
      if (retval) {
	goto plink1_dosage_ret_1;
      }
    }
    if (filter_attrib_indiv_fname) {
      retval = filter_attrib(filter_attrib_indiv_fname, filter_attrib_indiv_liststr, sorted_person_ids, ulii, max_person_id_len, person_id_map, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, 1);
      if (retval) {
	goto plink1_dosage_ret_1;
      }
    }
    if (filtername) {
      if (!mfilter_col) {
	mfilter_col = 1;
      }
      retval = filter_indivs_file(filtername, sorted_person_ids, ulii, max_person_id_len, person_id_map, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, filtervals_flattened, mfilter_col);
      if (retval) {
	goto plink1_dosage_ret_1;
      }
    }
    wkspace_reset(wkspace_mark);
  }
  if (gender_unk_ct && (!(sex_missing_pheno & ALLOW_NO_SEX))) {
    uii = popcount_longs_exclude(pheno_nm, sex_nm, unfiltered_indiv_ctl);
    if (uii) {
      bitfield_and(pheno_nm, sex_nm, unfiltered_indiv_ctl);
      logprint("Warning: Ignoring phenotypes of missing-sex samples.  If you don't want those\nphenotypes to be ignored, use the --allow-no-sex flag.\n");
    }
  }
  if (filter_flags & FILTER_PRUNE) {
    bitfield_ornot(indiv_exclude, pheno_nm, unfiltered_indiv_ctl);
    zero_trailing_bits(indiv_exclude, unfiltered_indiv_ct);
    indiv_exclude_ct = popcount_longs(indiv_exclude, unfiltered_indiv_ctl);
    if (indiv_exclude_ct == unfiltered_indiv_ct) {
      LOGPRINTF("Error: All %s removed by --prune.\n", g_species_plural);
      goto plink1_dosage_ret_ALL_SAMPLES_EXCLUDED;
    }
    LOGPRINTF("--prune: %" PRIuPTR " %s remaining.\n", unfiltered_indiv_ct - indiv_exclude_ct, species_str(unfiltered_indiv_ct == indiv_exclude_ct + 1));
  }

  if (filter_flags & (FILTER_BINARY_CASES | FILTER_BINARY_CONTROLS)) {
    if (!pheno_c) {
      logprint("Error: --filter-cases/--filter-controls requires a case/control phenotype.\n");
      goto plink1_dosage_ret_INVALID_CMDLINE;
    }
    ii = indiv_exclude_ct;
    filter_indivs_bitfields(unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, pheno_c, (filter_flags / FILTER_BINARY_CASES) & 1, pheno_nm);
    if (indiv_exclude_ct == unfiltered_indiv_ct) {
      LOGPRINTF("Error: All %s removed due to case/control status (--filter-%s).\n", g_species_plural, (filter_flags & FILTER_BINARY_CASES)? "cases" : "controls");
      goto plink1_dosage_ret_ALL_SAMPLES_EXCLUDED;
    }
    ii = indiv_exclude_ct - ii;
    LOGPRINTF("%d %s removed due to case/control status (--filter-%s).\n", ii, species_str(ii), (filter_flags & FILTER_BINARY_CASES)? "cases" : "controls");
  }
  if (filter_flags & (FILTER_BINARY_FEMALES | FILTER_BINARY_MALES)) {
    ii = indiv_exclude_ct;
    filter_indivs_bitfields(unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, sex_male, (filter_flags / FILTER_BINARY_MALES) & 1, sex_nm);
    if (indiv_exclude_ct == unfiltered_indiv_ct) {
      LOGPRINTF("Error: All %s removed due to gender filter (--filter-%s).\n", g_species_plural, (filter_flags & FILTER_BINARY_MALES)? "males" : "females");
      goto plink1_dosage_ret_ALL_SAMPLES_EXCLUDED;
    }
    ii = indiv_exclude_ct - ii;
    LOGPRINTF("%d %s removed due to gender filter (--filter-%s).\n", ii, species_str(ii), (filter_flags & FILTER_BINARY_MALES)? "males" : "females");
  }
  if (filter_flags & (FILTER_BINARY_FOUNDERS | FILTER_BINARY_NONFOUNDERS)) {
    ii = indiv_exclude_ct;
    filter_indivs_bitfields(unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, founder_info, (filter_flags / FILTER_BINARY_FOUNDERS) & 1, NULL);
    if (indiv_exclude_ct == unfiltered_indiv_ct) {
      LOGPRINTF("Error: All %s removed due to founder status (--filter-%s).\n", g_species_plural, (filter_flags & FILTER_BINARY_FOUNDERS)? "founders" : "nonfounders");
      goto plink1_dosage_ret_ALL_SAMPLES_EXCLUDED;
    }
    ii = indiv_exclude_ct - ii;
    LOGPRINTF("%d %s removed due to founder status (--filter-%s).\n", ii, species_str(ii), (filter_flags & FILTER_BINARY_FOUNDERS)? "founders" : "nonfounders");
  }
  if (cluster_ptr->fname || (misc_flags & MISC_FAMILY_CLUSTERS)) {
    retval = load_clusters(cluster_ptr->fname, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, person_ids, max_person_id_len, mwithin_col, (misc_flags / MISC_LOAD_CLUSTER_KEEP_NA) & 1, &cluster_ct, &cluster_map, &cluster_starts, &cluster_ids, &max_cluster_id_len, cluster_ptr->keep_fname, cluster_ptr->keep_flattened, cluster_ptr->remove_fname, cluster_ptr->remove_flattened);
    if (retval) {
      goto plink1_dosage_ret_1;
    }
  }
  indiv_ct = unfiltered_indiv_ct - indiv_exclude_ct;
  if (!indiv_ct) {
    LOGPRINTF("Error: No %s pass QC.\n", g_species_plural);
    goto plink1_dosage_ret_ALL_SAMPLES_EXCLUDED;
  }
  indiv_ctl = (indiv_ct + (BITCT - 1)) / BITCT;
  if (g_thread_ct > 1) {
    logprint("Using 1 thread (no multithreaded calculations invoked).\n");
  } else {
    logprint("Using 1 thread.\n");
  }
  if (load_map) {
    plink_maxsnp = calc_plink_maxsnp(unfiltered_marker_ct, marker_exclude, unfiltered_marker_ct - marker_exclude_ct, marker_ids, max_marker_id_len);
  }
  if ((filter_flags & FILTER_MAKE_FOUNDERS) && (!(misc_flags & MISC_MAKE_FOUNDERS_FIRST))) {
    if (make_founders(unfiltered_indiv_ct, indiv_ct, person_ids, max_person_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, (misc_flags / MISC_MAKE_FOUNDERS_REQUIRE_2_MISSING) & 1, indiv_exclude, founder_info)) {
      goto plink1_dosage_ret_NOMEM;
    }
  }
  if (covar_fname) {
    // update this as more covariate-referencing commands are added
    if (!do_glm) {
      logprint("Warning: Ignoring --covar since no commands reference the covariates.\n");
    } else {
      retval = load_covars(covar_fname, unfiltered_indiv_ct, indiv_exclude, indiv_ct, person_ids, max_person_id_len, missing_phenod, covar_modifier, covar_range_list_ptr, 0, &covar_ct, &covar_names, &max_covar_name_len, pheno_nm, &covar_nm, &covar_d, NULL, NULL);
      if (retval) {
	goto plink1_dosage_ret_1;
      }
    }
  }
  bitfield_andnot(pheno_nm, indiv_exclude, unfiltered_indiv_ctl);
  if (pheno_c) {
    bitfield_and(pheno_c, pheno_nm, unfiltered_indiv_ctl);
  }
  bitfield_andnot(founder_info, indiv_exclude, unfiltered_indiv_ctl);
  bitfield_andnot(sex_nm, indiv_exclude, unfiltered_indiv_ctl);
  if (gender_unk_ct) {
    gender_unk_ct = indiv_ct - popcount_longs(sex_nm, unfiltered_indiv_ctl);
  }
  bitfield_and(sex_male, sex_nm, unfiltered_indiv_ctl);

  pheno_nm_ct = popcount_longs(pheno_nm, unfiltered_indiv_ctl);

  if (load_map) {
    if (unfiltered_marker_ct == marker_exclude_ct) {
      logprint("Error: No variants remaining.\n");
      goto plink1_dosage_ret_ALL_MARKERS_EXCLUDED;
    }
    if (min_bp_space) {
      if (map_is_unsorted & UNSORTED_BP) {
	logprint("Error: --bp-space requires a sorted .bim file.  Retry this command after using\n--make-bed to sort your data.\n");
	goto plink1_dosage_ret_INVALID_FORMAT;
      }
      enforce_min_bp_space(min_bp_space, unfiltered_marker_ct, marker_exclude, marker_pos, &marker_exclude_ct, chrom_info_ptr);
    }
    marker_ct = unfiltered_marker_ct - marker_exclude_ct;
    retval = sort_item_ids(&sorted_marker_ids, &marker_id_map, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_ids, max_marker_id_len, 0, 0, strcmp_deref);
    if (retval) {
      goto plink1_dosage_ret_1;
    }
    LOGPRINTF("%" PRIuPTR " variant%s and %" PRIuPTR " %s pass filters and QC.\n", marker_ct, (marker_ct == 1)? "" : "s", indiv_ct, species_str(indiv_ct));
  } else {
    LOGPRINTF("%" PRIuPTR " %s pass filters and QC.\n", indiv_ct, species_str(indiv_ct));
  }
  if (!pheno_nm_ct) {
    if (do_glm) {
      logprint("Error: No phenotypes present.\n");
      goto plink1_dosage_ret_ALL_SAMPLES_EXCLUDED;
    }
    logprint("Note: No phenotypes present.\n");
  } else if (pheno_c) {
    pheno_ctrl_ct = pheno_nm_ct - popcount_longs(pheno_c, unfiltered_indiv_ctl);
    if (pheno_nm_ct != indiv_ct) {
      sprintf(logbuf, "Among remaining phenotypes, %u %s and %u %s.  (%" PRIuPTR " phenotype%s missing.)\n", pheno_nm_ct - pheno_ctrl_ct, (pheno_nm_ct - pheno_ctrl_ct == 1)? "is a case" : "are cases", pheno_ctrl_ct, (pheno_ctrl_ct == 1)? "is a control" : "are controls", indiv_ct - pheno_nm_ct, (indiv_ct - pheno_nm_ct == 1)? " is" : "s are");
    } else {
      sprintf(logbuf, "Among remaining phenotypes, %u %s and %u %s.\n", pheno_nm_ct - pheno_ctrl_ct, (pheno_nm_ct - pheno_ctrl_ct == 1)? "is a case" : "are cases", pheno_ctrl_ct, (pheno_ctrl_ct == 1)? "is a control" : "are controls");
    }
    wordwrap(logbuf, 0);
    logprintb();
  } else {
    logprint("Phenotype data is quantitative.\n");
  }

  // now the actual --dosage logic begins
  // 1. either load single file, or
  //    a. determine whether batch numbers are present; if yes, determine set
  //       of batch numbers and frequencies, sort
  //    b. determine max file path length(s) on first pass as well
  //    c. actually load on second pass
  // 2. initialize output writer/compressor, header line if necessary
  // 3. loop through batches
  //    a. partial memory reset, loop through file headers (determine sample
  //       counts in each file, skip counts, and ID mappings, error out on FID
  //       duplication)
  //    b. if association analysis, collapse phenotype/covariate data to
  //       account for missing samples in current batch
  //    c. process one variant at a time, all files in batch must have them in
  //       same order
  //    d. to enable parallel gzip, we can append to a write buffer and flush
  //       it when it's within range of overflowing on the next line.  this
  //       requires stable write buffer sizing, though, so we'll just use
  //       serial compression until the rest of the code is finished and we
  //       know how large all our other memory allocations are.
  // 4. final write loop if necessary
  if (doip->modifier & DOSAGE_LIST) {
    retval = open_and_load_to_first_token(&infile, doip->fname, MAXLINELEN, '\0', "--dosage list file", tbuf, &bufptr, &line_idx);
    if (retval) {
      goto plink1_dosage_ret_1;
    }
    batch_ct = count_tokens(bufptr) - sepheader - 1; // underflow ok
    if (batch_ct) {
      if (batch_ct != 1) {
        logprint("Error: Unexpected number of columns in --dosage list file.\n");
	goto plink1_dosage_ret_INVALID_FORMAT;
      }
      batch_sizes = (uint32_t*)wkspace_base;
      uiptr = batch_sizes;
      uiptr2 = (uint32_t*)(&(wkspace_base[wkspace_left / 2]));
    } else {
      if (wkspace_alloc_ui_checked(&batch_sizes, sizeof(int32_t))) {
	goto plink1_dosage_ret_NOMEM;
      }
    }
    // now batch_ct = 0 indicates no batch column (but one actual batch)
    while (1) {
      if (batch_ct) {
	if (uiptr == uiptr2) {
	  goto plink1_dosage_ret_NOMEM;
	}
	bufptr2 = token_endnn(bufptr);
        if (scan_int32(bufptr, (int32_t*)uiptr)) {
	  sprintf(logbuf, "Error: Invalid batch number on line %" PRIuPTR " of --dosage list file.\n", line_idx);
          goto plink1_dosage_ret_INVALID_FORMAT_2;
	}
	uiptr++;
	bufptr = skip_initial_spaces(bufptr2);
	if (is_eoln_kns(*bufptr)) {
          sprintf(logbuf, "Error: Line %" PRIuPTR " of --dosage list file has fewer tokens than expected.\n", line_idx);
	  goto plink1_dosage_ret_INVALID_FORMAT_2;
	}
      }
      bufptr2 = token_endnn(bufptr);
      ulii = (uintptr_t)(bufptr2 - bufptr);
      if (ulii >= max_fn_len) {
	max_fn_len = ulii + 1;
      }
      if (sepheader) {
	bufptr = skip_initial_spaces(bufptr2);
	if (!bufptr) {
          sprintf(logbuf, "Error: Line %" PRIuPTR " of --dosage list file has fewer tokens than expected.\n", line_idx);
	  goto plink1_dosage_ret_INVALID_FORMAT_2;
	}
	bufptr2 = token_endnn(bufptr);
	ulii = (uintptr_t)(bufptr2 - bufptr);
	if (ulii >= max_sepheader_len) {
          max_sepheader_len = ulii + 1;
	}
      }
      bufptr = skip_initial_spaces(bufptr2);
      if (!is_eoln_kns(*bufptr)) {
	sprintf(logbuf, "Error: Line %" PRIuPTR " of --dosage list file has more tokens than expected.\n", line_idx);
	goto plink1_dosage_ret_INVALID_FORMAT_2;
      }
      infile_ct++;
    plink1_dosage_next_list_line:
      if (!fgets(tbuf, MAXLINELEN, infile)) {
	if (ferror(infile)) {
	  goto plink1_dosage_ret_READ_FAIL;
	}
	break;
      }
      line_idx++;
      if (!tbuf[MAXLINELEN - 1]) {
        sprintf(logbuf, "Error: Line %" PRIuPTR " of --dosage list file is pathologically long.\n", line_idx);
        goto plink1_dosage_ret_INVALID_FORMAT_2;
      }
      bufptr = skip_initial_spaces(tbuf);
      if (is_eoln_kns(*bufptr)) {
	goto plink1_dosage_next_list_line;
      }
    }
    // infile_ct must be at least 1 due to use of open_and_load_to_first_token
    if (batch_ct) {
#ifdef __cplusplus
      std::sort((int32_t*)batch_sizes, (int32_t*)uiptr);
#else
      qsort(batch_sizes, infile_ct, sizeof(int32_t), intcmp2);
#endif
      batch_sizes = (uint32_t*)wkspace_alloc(infile_ct * sizeof(int32_t));
      // temporary batch size buffer
      uiptr3 = (uint32_t*)top_alloc(&topsize, infile_ct * sizeof(int32_t));

      uii = batch_sizes[0];
      uiptr2 = &(batch_sizes[1]);
      uiptr3[0] = 1;
      while (uiptr2 < uiptr) {
	ujj = *uiptr2++;
	if (ujj != uii) {
          uiptr3[batch_ct++] = 1;
	  uii = ujj;
	} else {
	  uiptr3[batch_ct - 1] += 1;
	  while (uiptr2 < uiptr) {
	    ujj = *uiptr2++;
	    if (ujj != uii) {
	      uiptr3[batch_ct] = 1;
	      batch_sizes[batch_ct++] = ujj;
	      uii = ujj;
	    } else {
	      uiptr3[batch_ct - 1] += 1;
	    }
	  }
	  break;
	}
      }

      // batch numbers
      uiptr = (uint32_t*)top_alloc(&topsize, batch_ct * sizeof(int32_t));
      if (!uiptr) {
	goto plink1_dosage_ret_NOMEM;
      }
      memcpy(uiptr, batch_sizes, batch_ct * sizeof(int32_t));
      wkspace_shrink_top(batch_sizes, batch_ct * sizeof(int32_t));
      memcpy(batch_sizes, uiptr3, batch_ct * sizeof(int32_t));
      // convert uiptr3 to write offset array
      uii = uiptr3[0];
      for (ujj = 1; ujj < batch_ct; ujj++) {
        uii += uiptr3[ujj];
	uiptr3[ujj] = uii;
      }
      wkspace_left -= topsize;
    } else {
      uiptr3 = batch_sizes;
      uiptr3[0] = 0;
    }
    if (wkspace_alloc_c_checked(&fnames, infile_ct * max_fn_len)) {
      goto plink1_dosage_ret_NOMEM2;
    }
    if (sepheader) {
      if (wkspace_alloc_c_checked(&sep_fnames, infile_ct * max_sepheader_len)) {
	goto plink1_dosage_ret_NOMEM2;
      }
    }
    wkspace_left += topsize;
    topsize = 0;
    rewind(infile);
    while (fgets(tbuf, MAXLINELEN, infile)) {
      bufptr = skip_initial_spaces(tbuf);
      if (is_eoln_kns(*bufptr)) {
	continue;
      }
      if (batch_ct) {
	scan_int32(bufptr, &ii);
	uii = int32arr_greater_than((int32_t*)uiptr, batch_ct, ii);
	file_idx = uiptr3[ujj];
	uiptr3[uii] += 1;
	bufptr = next_token(bufptr);
      } else {
        file_idx = uiptr3[0];
	uiptr3[0] += 1;
      }
      bufptr2 = token_endnn(bufptr);
      memcpyx(&(fnames[file_idx * max_fn_len]), bufptr, (uintptr_t)(bufptr2 - bufptr), '\0');
      if (sepheader) {
	bufptr = skip_initial_spaces(bufptr2);
	bufptr2 = token_endnn(bufptr);
	memcpyx(&(sep_fnames[file_idx * max_sepheader_len]), bufptr, (uintptr_t)(bufptr2 - bufptr), '\0');
      }
    }
    if (fclose_null(&infile)) {
      goto plink1_dosage_ret_READ_FAIL;
    }
    if (batch_ct) {
      for (uii = 0; uii < batch_ct; uii++) {
	if (batch_sizes[uii] > max_batch_size) {
	  max_batch_size = batch_sizes[uii];
	}
      }
    } else {
      batch_ct = 1;
      max_batch_size = infile_ct;
    }
  } else {
    uii = strlen(doip->fname) + 1;
    if (wkspace_alloc_ui_checked(&batch_sizes, sizeof(int32_t)) ||
        wkspace_alloc_c_checked(&fnames, uii + 1)) {
      goto plink1_dosage_ret_NOMEM;
    }
    batch_sizes[0] = 1;
    memcpy(fnames, doip->fname, uii);
    infile_ct = 1;
  }
  if (wkspace_alloc_ui_checked(&file_icts, max_batch_size * sizeof(int32_t)) ||
      wkspace_alloc_ul_checked(&line_idx_arr, max_batch_size * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&cur_indivs, indiv_ctl * sizeof(intptr_t)) ||
      wkspace_alloc_ui_checked(&read_idx_to_indiv_idx, indiv_ct * sizeof(int32_t)) ||
      wkspace_alloc_ui_checked(&skip_vals, indiv_ct * sizeof(int32_t)) ||
      wkspace_alloc_d_checked(&cur_dosages, indiv_ct * sizeof(double))) {
    goto plink1_dosage_ret_NOMEM;
  }
  if (max_batch_size > 1) {
    if (wkspace_alloc_c_checked(&cur_marker_id_buf, MAX_ID_LEN)) {
      goto plink1_dosage_ret_NOMEM;
    }
  }
  gz_infiles = (gzFile*)wkspace_alloc(infile_ct * sizeof(gzFile));
  if (!gz_infiles) {
    infile_ct = 0;
    goto plink1_dosage_ret_NOMEM;
  }
  for (uii = 0; uii < infile_ct; uii++) {
    gz_infiles[uii] = NULL;
  }
  if (noheader) {
    if (infile_ct != 1) {
      logprint("Error: --dosage 'noheader' modifier cannot be used with multiple input files.\n");
      goto plink1_dosage_ret_INVALID_CMDLINE;
    }
    // sorted_person_ids = NULL;
  } else {
    retval = sort_item_ids(&sorted_person_ids, &person_id_map, unfiltered_indiv_ct, indiv_exclude, indiv_exclude_ct, person_ids, max_person_id_len, 0, 1, strcmp_deref);
    if (retval) {
      goto plink1_dosage_ret_1;
    }
  }
  if (do_glm) {
    if (load_map) {
      bufptr = memcpya(tbuf, " CHR         SNP          BP", 28);
    } else {
      bufptr = memcpya(tbuf, "         SNP", 12);
    }
    bufptr = memcpya(tbuf, "  A1  A2     FRQ    INFO    ", 56);
    bufptr = memcpya(bufptr, pheno_c? "  OR" : "BETA", 4);
    bufptr = memcpya(bufptr, "      SE       P\n", 16);
    bufptr2 = memcpyb(outname_end, ".assoc.dosage", 14);
  } else if (count_occur) {
    // could just use a uint32_t array if .map provided
    htable = (Ll_ctstr_entry**)wkspace_alloc(HASHSIZE * sizeof(intptr_t));
    if (!htable) {
      goto plink1_dosage_ret_NOMEM;
    }
    for (uii = 0; uii < HASHSIZE; uii++) {
      htable[uii] = NULL;
    }
    bufptr2 = memcpyb(outname_end, ".occur.dosage", 14);
  } else {
    bufptr2 = memcpyb(outname_end, ".out.dosage", 12);
  }
  if (output_gz) {
    memcpy(bufptr2, ".gz", 4);
    if (gzopen_checked(&gz_outfile, outname, "wb")) {
      goto plink1_dosage_ret_OPEN_FAIL;
    }
    if (do_glm) {
      if (!gzwrite(gz_outfile, tbuf, bufptr - tbuf)) {
	goto plink1_dosage_ret_WRITE_FAIL;
      }
    } else {
      if (gzputs(gz_outfile, "SNP A1 A2 ") == -1) {
	goto plink1_dosage_ret_WRITE_FAIL;
      }
      for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
	bufptr = &(person_ids[indiv_idx * max_person_id_len]);
	bufptr2 = strchr(bufptr, '\t');
	*bufptr2 = ' ';
        if (gzputs(gz_outfile, bufptr) == -1) {
	  goto plink1_dosage_ret_WRITE_FAIL;
	}
	if (gzputc(gz_outfile, ' ') == -1) {
	  goto plink1_dosage_ret_WRITE_FAIL;
	}
	*bufptr2 = '\t';
      }
      if (gzputc(gz_outfile, '\n') == -1) {
	goto plink1_dosage_ret_WRITE_FAIL;
      }
    }
  } else {
    if (fopen_checked(&outfile, outname, "w")) {
      goto plink1_dosage_ret_OPEN_FAIL;
    }
    if (do_glm) {
      if (fwrite_checked(tbuf, bufptr - tbuf, outfile)) {
	goto plink1_dosage_ret_WRITE_FAIL;
      }
    } else {
      fputs("SNP A1 A2 ", outfile);
      for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
	bufptr = &(person_ids[indiv_idx * max_person_id_len]);
	bufptr2 = strchr(bufptr, '\t');
	*bufptr2 = ' ';
        fputs(bufptr, outfile);
	putc(' ', outfile);
	*bufptr2 = '\t';
      }
      if (putc_checked('\n', outfile)) {
	goto plink1_dosage_ret_WRITE_FAIL;
      }
    }
  }
  wkspace_mark = wkspace_base;
  for (batch_idx = 0; batch_idx < batch_ct; batch_idx++, file_idx_start += cur_batch_size) {
    cur_batch_size = batch_sizes[batch_idx];
    read_idx = 0;
    loadbuf_size = wkspace_left;
    if (loadbuf_size > MAXLINEBUFLEN) {
      loadbuf_size = MAXLINEBUFLEN;
    } else if (loadbuf_size <= MAXLINELEN) {
      goto plink1_dosage_ret_NOMEM;
    }
    loadbuf = (char*)wkspace_base;
    loadbuf[loadbuf_size - 1] = ' ';
    fill_ulong_zero(cur_indivs, indiv_ctl);
    for (file_idx = 0; file_idx < cur_batch_size; file_idx++) {
      read_idx_start = read_idx;
      if (sepheader) {
	if (gzopen_checked(&(gz_infiles[file_idx]), &(sep_fnames[(file_idx + file_idx_start) * max_sepheader_len]), "rb")) {
	  goto plink1_dosage_ret_OPEN_FAIL;
	}
	line_idx = 0;
	uii = 1; // current skip value
	while (gzgets(gz_infiles[file_idx], tbuf, MAXLINELEN)) {
	  line_idx++;
	  if (!tbuf[MAXLINELEN - 1]) {
	    sprintf(logbuf, "Error: Line %" PRIuPTR " of %s is pathologically long.\n", line_idx, &(sep_fnames[(file_idx + file_idx_start) * max_sepheader_len]));
	    goto plink1_dosage_ret_INVALID_FORMAT_WW;
	  }
          bufptr = skip_initial_spaces(tbuf);
          if (is_eoln_kns(*bufptr)) {
	    continue;
	  }
          if (bsearch_read_fam_indiv(&(tbuf[MAXLINELEN]), sorted_person_ids, max_person_id_len, indiv_ct, bufptr, &bufptr2, &ii)) {
            sprintf(logbuf, "Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, &(sep_fnames[(file_idx + file_idx_start) * max_sepheader_len]));
	    goto plink1_dosage_ret_INVALID_FORMAT_WW;
	  }
	  if (ii == -1) {
	    uii += format_val;
	  } else {
	    ii = person_id_map[(uint32_t)ii];
	    if (is_set(cur_indivs, ii)) {
	      bufptr = &(sorted_person_ids[((uint32_t)ii) * max_person_id_len]);
	      *strchr(bufptr, '\t') = ' ';
	      sprintf(logbuf, "Error: '%s' appears multiple times.\n", bufptr);
	      goto plink1_dosage_ret_INVALID_FORMAT_WW;
	    }
	    set_bit(cur_indivs, ii);
	    read_idx_to_indiv_idx[read_idx] = (uint32_t)ii;
            skip_vals[read_idx++] = uii;
	    uii = 1;
	  }
	}
	if (gzclose(gz_infiles[file_idx]) != Z_OK) {
	  gz_infiles[file_idx] = NULL;
	  goto plink1_dosage_ret_READ_FAIL;
	}
	gz_infiles[file_idx] = NULL;
	if (read_idx_start == read_idx) {
          sprintf(logbuf, "Error: %s is empty.\n", &(sep_fnames[(file_idx + file_idx_start) * max_sepheader_len]));
          goto plink1_dosage_ret_INVALID_FORMAT_WW;
	}
      }
      if (gzopen_checked(&(gz_infiles[file_idx]), &(fnames[(file_idx + file_idx_start) * max_fn_len]), "rb")) {
	goto plink1_dosage_ret_OPEN_FAIL;
      }
      line_idx = 0;
      if (noheader) {
	for (read_idx = 0; read_idx < indiv_ct; read_idx++) {
	  read_idx_to_indiv_idx[read_idx] = read_idx;
	}
	for (read_idx = 0; read_idx < indiv_ct; read_idx++) {
	  skip_vals[read_idx] = 1;
	}
	fill_all_bits(cur_indivs, indiv_ct);
      } else if (!sepheader) {
	do {
	  if (!gzgets(gz_infiles[file_idx], loadbuf, loadbuf_size)) {
            sprintf(logbuf, "Error: %s is empty.\n", &(fnames[(file_idx + file_idx_start) * max_fn_len]));
	    goto plink1_dosage_ret_INVALID_FORMAT_WW;
	  }
	  line_idx++;
	  if (!loadbuf[loadbuf_size - 1]) {
	    goto plink1_dosage_ret_LONG_LINE;
	  }
	  bufptr = skip_initial_spaces(loadbuf);
	} while (is_eoln_kns(*bufptr));
	bufptr = next_token_multz(bufptr, skip0);
	bufptr2 = next_token_mult(bufptr, skip1p1);
	if (no_more_tokens(bufptr2)) {
	  goto plink1_dosage_ret_MISSING_TOKENS;
	}
	if (strcmp_se(bufptr, "SNP", 3)) {
	  sprintf(logbuf, "Error: Column %u of %s's header isn't 'SNP'.\n", skip0 + 1, &(fnames[(file_idx + file_idx_start) * max_fn_len]));
	  goto plink1_dosage_ret_INVALID_FORMAT_WW;
	} else if (strcmp_se(bufptr2, "A1", 2)) {
	  sprintf(logbuf, "Error: Column %u of %s's header isn't 'A1'.\n", skip0 + skip1p1 + 1, &(fnames[(file_idx + file_idx_start) * max_fn_len]));
	  goto plink1_dosage_ret_INVALID_FORMAT_WW;
	}
	bufptr = next_token(bufptr2);
	bufptr2 = next_token_multz(bufptr, skip2);
	if (no_more_tokens(bufptr2)) {
	  goto plink1_dosage_ret_MISSING_TOKENS;
	}
	if (strcmp_se(bufptr, "A2", 2)) {
	  sprintf(logbuf, "Error: Column %u of %s's header isn't 'A2'.\n", skip0 + skip1p1 + 2, &(fnames[(file_idx + file_idx_start) * max_fn_len]));
	  goto plink1_dosage_ret_INVALID_FORMAT_WW;
	}
	uii = 1;
	bufptr = skip_initial_spaces(token_endnn(bufptr2));
	while (!is_eoln_kns(*bufptr)) {
          if (bsearch_read_fam_indiv(tbuf, sorted_person_ids, max_person_id_len, indiv_ct, bufptr, &bufptr2, &ii)) {
	    sprintf(logbuf, "Error: Header of %s has an odd number of tokens in the FID/IID section.\n", &(fnames[(file_idx + file_idx_start) * max_fn_len]));
	    goto plink1_dosage_ret_INVALID_FORMAT_WW;
	  }
	  if (ii == -1) {
	    uii += format_val;
	  } else {
	    ii = person_id_map[(uint32_t)ii];
	    if (is_set(cur_indivs, ii)) {
	      bufptr = &(sorted_person_ids[((uint32_t)ii) * max_person_id_len]);
	      *strchr(bufptr, '\t') = ' ';
	      sprintf(logbuf, "Error: '%s' appears multiple times.\n", bufptr);
	      goto plink1_dosage_ret_INVALID_FORMAT_WW;
	    }
	    set_bit(cur_indivs, ii);
	    read_idx_to_indiv_idx[read_idx] = (uint32_t)ii;
            skip_vals[read_idx++] = uii;
	    uii = 1;
	  }
	  bufptr = bufptr2;
	}
	if (read_idx_start == read_idx) {
	  sprintf(logbuf, "Error: Header of %s has no tokens in the FID/IID section.\n", &(fnames[(file_idx + file_idx_start) * max_fn_len]));
	  goto plink1_dosage_ret_INVALID_FORMAT_WW;
	}
      }
      file_icts[file_idx] = read_idx - read_idx_start;
      line_idx_arr[file_idx] = line_idx;
    }

    while (1) {
      read_idx_start = 0;
      for (file_idx = 0; file_idx < cur_batch_size; file_idx++) {
	line_idx = line_idx_arr[file_idx];
	do {
	  if (!gzgets(gz_infiles[file_idx], loadbuf, loadbuf_size)) {
	    if (file_idx) {
	      logprint("Error: Misaligned dosage data files.\n");
	      goto plink1_dosage_ret_INVALID_FORMAT;
	    }
	    goto plink1_dosage_end_loop;
	  }
	  line_idx++;
	  if (!loadbuf[loadbuf_size - 1]) {
	    goto plink1_dosage_ret_LONG_LINE;
	  }
          bufptr = skip_initial_spaces(loadbuf);
	} while (is_eoln_kns(*bufptr));
	line_idx_arr[file_idx] = line_idx;
	bufptr = next_token_multz(bufptr, skip0);
	bufptr3 = next_token_mult(bufptr, skip1p1);
	bufptr5 = next_token(bufptr3);
	if (no_more_tokens(bufptr5)) {
	  goto plink1_dosage_ret_MISSING_TOKENS;
	}
	bufptr2 = token_endnn(bufptr);
	bufptr4 = token_endnn(bufptr3);
	bufptr6 = token_endnn(bufptr5);
        slen = (uintptr_t)(bufptr2 - bufptr);
	if (slen > MAX_ID_LEN) {
	  sprintf(logbuf, "Error: Line %" PRIuPTR " of %s has an excessively long variant ID.\n", line_idx, &(fnames[(file_idx + file_idx_start) * max_fn_len]));
	  goto plink1_dosage_ret_INVALID_FORMAT_WW;
	}
	if (!file_idx) {
	  memcpyx(cur_marker_id_buf, bufptr, slen, '\0');
	  cur_marker_id_len = slen;
	  if (sorted_marker_ids) {
	    ii = bsearch_str(bufptr, slen, sorted_marker_ids, max_marker_id_len, marker_ct);
	    if (ii == -1) {
	      marker_idx = ~ZEROLU;
	      continue;
	    }
	    marker_idx = marker_id_map[(uint32_t)ii];
	  }
	  if (allele_set(&a1_ptr, bufptr3, (uintptr_t)(bufptr4 - bufptr3))) {
	    goto plink1_dosage_ret_NOMEM;
	  }
	  if (allele_set(&a2_ptr, bufptr5, (uintptr_t)(bufptr6 - bufptr5))) {
	    goto plink1_dosage_ret_NOMEM;
	  }
	} else {
	  if ((slen != cur_marker_id_len) || memcmp(bufptr, cur_marker_id_buf, slen)) {
	    sprintf(logbuf, "Error: Variant ID mismatch between line %" PRIuPTR " of %s and line %" PRIuPTR " of %s.\n", line_idx_arr[0], &(fnames[file_idx_start * max_fn_len]), line_idx, &(fnames[(file_idx + file_idx_start) * max_fn_len]));
	    goto plink1_dosage_ret_INVALID_FORMAT_WW;
	  }
	  *bufptr4 = '\0';
	  *bufptr6 = '\0';
	  if (strcmp(bufptr3, cur_a1) || strcmp(bufptr5, cur_a2)) {
	    sprintf(logbuf, "Error: Allele code mismatch between line %" PRIuPTR " of %s and line %" PRIuPTR " of %s.\n", line_idx_arr[0], &(fnames[file_idx_start * max_fn_len]), line_idx, &(fnames[(file_idx + file_idx_start) * max_fn_len]));
	    goto plink1_dosage_ret_INVALID_FORMAT_WW;
	  }
	  if (marker_idx == ~ZEROLU) {
	    continue;
	  }
	}

	if (count_occur) {
	  *bufptr2 = '\0';
          uii = hashval2(bufptr, slen++);
	  ll_pptr = &(htable[uii]);
	  while (1) {
	    ll_ptr = *ll_pptr;
	    if (!ll_ptr) {
	      distinct_id_ct++;
	      topsize += ((slen + sizeof(Ll_ctstr_entry) + 15) & (~(15 * ONELU)));
	      loadbuf_size = wkspace_left - topsize;
              ll_ptr = (Ll_ctstr_entry*)(&(wkspace_base[loadbuf_size]));
	      ll_ptr->next = NULL;
	      memcpy(ll_ptr->ss, bufptr, slen);
	      if (slen > max_occur_id_len) {
		max_occur_id_len = slen;
	      }
	      if (loadbuf_size >= MAXLINEBUFLEN) {
		loadbuf_size = MAXLINEBUFLEN;
	      } else if (loadbuf_size > MAXLINELEN) {
                loadbuf[loadbuf_size - 1] = ' ';
	      } else {
		goto plink1_dosage_ret_NOMEM2;
	      }
	      break;
	    }
            if (!strcmp(ll_ptr->ss, bufptr)) {
	      ll_ptr->ct += 1;
	      break;
	    }
	    ll_pptr = &(ll_ptr->next);
	  }
	} else {
	  *bufptr6 = ' ';
	  bufptr = bufptr5;
	  read_idx = read_idx_start + file_icts[file_idx];
	  if (format_val == 1) {
	    for (; read_idx_start < read_idx; read_idx_start++) {
	      bufptr = next_token_mult(bufptr, skip_vals[read_idx_start]);
	      if (!bufptr) {
		goto plink1_dosage_ret_MISSING_TOKENS;
	      }
	      if (scan_double(bufptr, &dxx)) {
		goto plink1_dosage_ret_INVALID_DOSAGE;
	      }
	      if (!dose1) {
		dxx *= 0.5;
	      }
	      if ((dxx > 1.0) || (dxx < 0.0)) {
		goto plink1_dosage_ret_INVALID_DOSAGE;
	      }
	      cur_dosages[read_idx_to_indiv_idx[read_idx_start]] = dxx;
	    }
	  } else {
	    for (; read_idx_start < read_idx; read_idx_start++) {
	      bufptr = next_token_mult();
	      bufptr2 = next_token(bufptr);
	      if (no_more_tokens(bufptr2)) {
		goto plink1_dosage_ret_MISSING_TOKENS;
	      }
	      if (scan_double(bufptr, &dxx) || scan_double(bufptr2, &dyy)) {
		goto plink1_dosage_ret_INVALID_DOSAGE;
	      }
	      dxx += dyy * 0.5;
	      if ((dxx > 1.0) || (dxx < 0.0)) {
		goto plink1_dosage_ret_INVALID_DOSAGE;
	      }
	      cur_dosages[read_idx_to_indiv_idx[read_idx_start]] = dxx;
	      bufptr = next_token_mult(bufptr2, formatm1);
	    }
	  }
	}
      }
      if (a1_ptr) {
	if (a1_ptr[1]) {
	  free(a1_ptr);
	}
	a1_ptr = NULL;
	if (a2_ptr[1]) {
	  free(a2_ptr);
	}
	a2_ptr = NULL;
      }
    }
  plink1_dosage_end_loop:
    for (file_idx = 0; file_idx < cur_batch_size; file_idx++) {
      if (gzclose(gz_infiles[file_idx]) != Z_OK) {
        gz_infiles[file_idx] = NULL;
        goto plink1_dosage_ret_READ_FAIL;
      }
      gz_infiles[file_idx] = NULL;
    }
    wkspace_reset(wkspace_mark);
  }

  while (0) {
  plink1_dosage_ret_NOMEM2:
    wkspace_left += topsize;
  plink1_dosage_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  plink1_dosage_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  plink1_dosage_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  plink1_dosage_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  plink1_dosage_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  plink1_dosage_ret_LONG_LINE:
    if (loadbuf_size == MAXLINEBUFLEN) {
      LOGPRINTFWW("Error: Line %" PRIuPTR " of %s is pathologically long.\n", line_idx, &(fnames[(file_idx + file_idx_start) * max_fn_len]));
      retval = RET_INVALID_FORMAT;
      break;
    }
    retval = RET_NOMEM;
    break;
  plink1_dosage_ret_INVALID_DOSAGE:
    LOGPRINTFWW("Error: Invalid dosage value on line %" PRIuPTR " of %s.\n", line_idx, &(fnames[(file_idx + file_idx_start) * max_fn_len]));
    retval = RET_INVALID_FORMAT;
    break;
  plink1_dosage_ret_MISSING_TOKENS:
    sprintf(logbuf, "Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, &(fnames[(file_idx + file_idx_start) * max_fn_len]));
  plink1_dosage_ret_INVALID_FORMAT_WW:
    wordwrap(logbuf, 0);
  plink1_dosage_ret_INVALID_FORMAT_2:
    logprintb();
  plink1_dosage_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  plink1_dosage_ret_ALL_MARKERS_EXCLUDED:
    retval = RET_ALL_MARKERS_EXCLUDED;
    break;
  plink1_dosage_ret_ALL_SAMPLES_EXCLUDED:
    retval = RET_ALL_SAMPLES_EXCLUDED;
    break;
  }
 plink1_dosage_ret_1:
  aligned_free_cond(pheno_c);
  free_cond(pheno_d);
  fclose_cond(phenofile);
  fclose_cond(infile);
  fclose_cond(outfile);
  gzclose_cond(gz_outfile);
  if (a1_ptr && a1_ptr[1]) {
    free(a1_ptr);
  }
  if (a2_ptr && a2_ptr[1]) {
    free(a2_ptr);
  }
  if (infile_ct) {
    for (uii = 0; uii < infile_ct; uii++) {
      gzclose_cond(gz_infiles[uii]);
    }
  }

  return retval;
  */
}
