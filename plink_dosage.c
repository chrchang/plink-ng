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

int32_t plink1_dosage(Dosage_info* doip, char* famname, char* mapname, char* outname, char* outname_end, char* phenoname, char* extractname, char* excludename, char* keepname, char* removename, char* keepfamname, char* removefamname, char* filtername, char* makepheno_str, char* phenoname_str, char* covar_fname, Two_col_params* qual_filter, Two_col_params* update_map, Two_col_params* update_name, char* update_ids_fname, char* update_parents_fname, char* update_sex_fname, char* filtervals_flattened, char* filter_attrib_fname, char* filter_attrib_liststr, char* filter_attrib_indiv_fname, char* filter_attrib_indiv_liststr, double qual_min_thresh, double qual_max_thresh, double thin_keep_prob, uint32_t thin_keep_ct, uint32_t min_bp_space, uint32_t mfilter_col, uint32_t fam_cols, int32_t missing_pheno, uint32_t mpheno_col, uint32_t pheno_modifier, Chrom_info* chrom_info_ptr, double tail_bottom, double tail_top, uint64_t misc_flags, uint64_t filter_flags, uint32_t sex_missing_pheno, uint32_t update_sex_col, Cluster_info* cluster_ptr, int32_t marker_pos_start, int32_t marker_pos_end, int32_t snp_window_size, char* markername_from, char* markername_to, char* markername_snp, Range_list* snps_range_list_ptr, uint32_t covar_modifier, Range_list* covar_range_list_ptr, uint32_t mwithin_col, uint32_t glm_modifier, double glm_vif_thresh, uint32_t glm_xchr_model, Range_list* parameters_range_list_ptr) {
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
  uint32_t* marker_pos = NULL;
  uint32_t* cluster_map = NULL;
  uint32_t* cluster_starts = NULL;
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
  uintptr_t ulii = 0;
  double missing_phenod = (double)missing_pheno;
  uint32_t load_map = (mapname[0] != '\0');
  uint32_t do_glm = (doip->modifier / DOSAGE_GLM) & 1;
  // uint32_t output_gz = doip->modifier & DOSAGE_ZOUT;
  uint32_t map_cols = 3;
  uint32_t map_is_unsorted = 0;
  uint32_t affection = 0;
  uint32_t plink_maxsnp = 0;
  uint32_t infile_ct = 0;
  int32_t retval = 0;
  char* missing_mid_templates[2];
  unsigned char* wkspace_mark;
  char* cptr;
  uint32_t* uiptr;
  uintptr_t marker_ct;
  uintptr_t unfiltered_indiv_ctl;
  uintptr_t indiv_ct;
  uintptr_t uljj;
  uint32_t gender_unk_ct;
  uint32_t pheno_nm_ct;
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
    retval = sort_item_ids(&cptr, &uiptr, unfiltered_indiv_ct, indiv_exclude, 0, person_ids, max_person_id_len, 0, 0, strcmp_deref);
    if (retval) {
      goto plink1_dosage_ret_1;
    }

    if (makepheno_str) {
      retval = makepheno_load(phenofile, makepheno_str, unfiltered_indiv_ct, cptr, max_person_id_len, uiptr, pheno_nm, &pheno_c);
      if (retval) {
        goto plink1_dosage_ret_1;
      }
    } else if (phenofile) {
      retval = load_pheno(phenofile, unfiltered_indiv_ct, 0, cptr, max_person_id_len, uiptr, missing_pheno, (misc_flags / MISC_AFFECTION_01) & 1, mpheno_col, phenoname_str, pheno_nm, &pheno_c, &pheno_d, NULL, 0);
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
      retval = sort_item_ids(&cptr, &uiptr, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_ids, max_marker_id_len, !uii, 0, strcmp_deref);
      if (retval) {
	goto plink1_dosage_ret_1;
      }
      ulii = unfiltered_marker_ct - marker_exclude_ct;

      if (update_map) {
	retval = update_marker_pos(update_map, cptr, ulii, max_marker_id_len, uiptr, marker_exclude, &marker_exclude_ct, marker_pos, &map_is_unsorted, chrom_info_ptr);
      } else if (update_name) {
	retval = update_marker_names(update_name, cptr, ulii, max_marker_id_len, uiptr, marker_ids);
	if (retval) {
	  goto plink1_dosage_ret_1;
	}
	if (extractname || excludename) {
	  wkspace_reset(wkspace_mark);
	  retval = sort_item_ids(&cptr, &uiptr, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_ids, max_marker_id_len, 0, 0, strcmp_deref);
	  if (retval) {
	    goto plink1_dosage_ret_1;
	  }
	  ulii = unfiltered_marker_ct - marker_exclude_ct;
	}
      }
      if (extractname) {
	if (!(misc_flags & MISC_EXTRACT_RANGE)) {
	  retval = include_or_exclude(extractname, cptr, ulii, max_marker_id_len, uiptr, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, 0);
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
	  retval = include_or_exclude(excludename, cptr, ulii, max_marker_id_len, uiptr, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, 1);
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
	retval = filter_attrib(filter_attrib_fname, filter_attrib_liststr, cptr, ulii, max_marker_id_len, uiptr, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, 0);
	if (retval) {
	  goto plink1_dosage_ret_1;
	}
      }
      if (qual_filter) {
	retval = filter_qual_scores(qual_filter, qual_min_thresh, qual_max_thresh, cptr, ulii, max_marker_id_len, uiptr, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct);
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
    retval = sort_item_ids(&cptr, &uiptr, unfiltered_indiv_ct, indiv_exclude, indiv_exclude_ct, person_ids, max_person_id_len, 0, 0, strcmp_deref);
    if (retval) {
      goto plink1_dosage_ret_1;
    }
    ulii = unfiltered_indiv_ct - indiv_exclude_ct;
    if (update_ids_fname) {
      retval = update_indiv_ids(update_ids_fname, cptr, ulii, max_person_id_len, uiptr, person_ids);
      if (retval) {
	goto plink1_dosage_ret_1;
      }
    } else {
      if (update_parents_fname) {
	retval = update_indiv_parents(update_parents_fname, cptr, ulii, max_person_id_len, uiptr, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, founder_info);
	if (retval) {
	  goto plink1_dosage_ret_1;
	}
      }
      if (update_sex_fname) {
        retval = update_indiv_sexes(update_sex_fname, update_sex_col, cptr, ulii, max_person_id_len, uiptr, sex_nm, sex_male);
	if (retval) {
	  goto plink1_dosage_ret_1;
	}
      }
    }
    if (keepfamname) {
      retval = include_or_exclude(keepfamname, cptr, ulii, max_person_id_len, uiptr, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, 6);
      if (retval) {
	goto plink1_dosage_ret_1;
      }
    }
    if (keepname) {
      retval = include_or_exclude(keepname, cptr, ulii, max_person_id_len, uiptr, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, 2);
      if (retval) {
	goto plink1_dosage_ret_1;
      }
    }
    if (removefamname) {
      retval = include_or_exclude(removefamname, cptr, ulii, max_person_id_len, uiptr, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, 7);
      if (retval) {
	goto plink1_dosage_ret_1;
      }
    }
    if (removename) {
      retval = include_or_exclude(removename, cptr, ulii, max_person_id_len, uiptr, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, 3);
      if (retval) {
	goto plink1_dosage_ret_1;
      }
    }
    if (filter_attrib_indiv_fname) {
      retval = filter_attrib(filter_attrib_indiv_fname, filter_attrib_indiv_liststr, cptr, ulii, max_person_id_len, uiptr, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, 1);
      if (retval) {
	goto plink1_dosage_ret_1;
      }
    }
    if (filtername) {
      if (!mfilter_col) {
	mfilter_col = 1;
      }
      retval = filter_indivs_file(filtername, cptr, ulii, max_person_id_len, uiptr, unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, filtervals_flattened, mfilter_col);
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
  }
  logprint("Error: --dosage is currently under development.\n");
  retval = RET_CALC_NOT_YET_SUPPORTED;

  // now the actual --dosage logic begins
  // 1. initialize some memory buffers
  // 2. either load single file, or
  //    a. determine whether batch numbers are present; if yes, determine set
  //       of batch numbers and frequencies, sort
  //    b. determine max file path length(s) on first pass as well
  //    c. actually load on second pass
  // 3. initialize output writer/compressor, header line if necessary
  // 4. loop through batches
  //    a. partial memory reset, loop through file headers (determine sample
  //       counts in each file, skip counts, and ID mappings, error out on FID
  //       duplication)
  //    b. if association analysis, collapse phenotype/covariate data to
  //       account for missing samples in current batch
  //    c. process one variant at a time, all files in batch must have them in
  //       same order
  //    d. argh, this doesn't play well with parallel gzip.  oh well, just do
  //       ordinary write
  // 5. final write loop if necessary
  if (doip->modifier & DOSAGE_LIST) {
    if (fopen_checked(&infile, doip->fname, "r")) {
      goto plink1_dosage_ret_OPEN_FAIL;
    }
  } else {
  }

  while (0) {
  plink1_dosage_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  plink1_dosage_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  plink1_dosage_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
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
  if (infile_ct) {
    for (uii = 0; uii < infile_ct; uii++) {
      gzclose_cond(gz_infiles[uii]);
    }
  }
  return retval;
}

