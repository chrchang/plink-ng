#include <fcntl.h>
#include <time.h>

#ifndef _WIN32
#include <sys/mman.h>
#include <unistd.h>
#endif

#include <sys/stat.h>
#include <sys/types.h>
#include "plink_common.h"

#define PHENO_EPSILON 0.000030517578125

// last prime before 2^19
// (maybe want to pick the second-to-last prime instead since this is only 1
// less?  though 19-character offsets are unlikely to matter in practice...)
#define HASHSIZE 524287
#define HASHSIZE_S 524287

#ifdef __LP64__
#define HASHMEM 4194304
#define HASHMEM_S 4194304
#else
#define HASHMEM 2097152
#define HASHMEM_S 2097152
#endif

int32_t load_pheno(FILE* phenofile, uintptr_t unfiltered_indiv_ct, uintptr_t indiv_exclude_ct, char* sorted_person_ids, uintptr_t max_person_id_len, uint32_t* id_map, int32_t missing_pheno, uint32_t missing_pheno_len, uint32_t affection_01, uint32_t mpheno_col, char* phenoname_str, uintptr_t* pheno_nm, uintptr_t** pheno_c_ptr, double** pheno_d_ptr) {
  uint32_t affection = 1;
  uintptr_t* pheno_c = *pheno_c_ptr;
  double* pheno_d = *pheno_d_ptr;
  int32_t header_processed = 0;
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  uintptr_t indiv_ct = unfiltered_indiv_ct - indiv_exclude_ct;
  char case_char = affection_01? '1' : '2';
  uintptr_t* isz = NULL;
  double missing_phenod = (double)missing_pheno;
  int32_t retval = 0;
  char* loadbuf;
  uint32_t loadbuf_size;
  int32_t person_idx;
  char* bufptr0;
  char* bufptr;
  uint32_t tmp_len;
  uint32_t tmp_len2;
  uint32_t uii;
  double dxx;
  double dyy;
  if (pheno_d) {
    affection = 0;
  } else {
    if (wkspace_alloc_ul_checked(&isz, unfiltered_indiv_ctl * sizeof(intptr_t))) {
      goto load_pheno_ret_NOMEM;
    }
    fill_ulong_zero(isz, unfiltered_indiv_ctl);
    if (!pheno_c) {
      pheno_c = (uintptr_t*)malloc(unfiltered_indiv_ctl * sizeof(intptr_t));
      if (!pheno_c) {
	goto load_pheno_ret_NOMEM;
      }
      fill_ulong_zero(pheno_c, unfiltered_indiv_ctl);
      *pheno_c_ptr = pheno_c;
    }
  }
  // ----- phenotype file load -----
  // worthwhile to support very long lines here...
  if (wkspace_left > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  } else if (wkspace_left > MAXLINELEN) {
    loadbuf_size = wkspace_left;
  } else {
    goto load_pheno_ret_NOMEM;
  }
  loadbuf = (char*)wkspace_base;
  loadbuf[loadbuf_size - 1] = ' ';
  while (fgets(loadbuf, loadbuf_size, phenofile) != NULL) {
    if (!loadbuf[loadbuf_size - 1]) {
      if (loadbuf_size == MAXLINEBUFLEN) {
	goto load_pheno_ret_INVALID_FORMAT;
      } else {
	goto load_pheno_ret_NOMEM;
      }
    }
    bufptr0 = skip_initial_spaces(loadbuf);
    if (is_eoln_kns(*bufptr0)) {
      continue;
    }
    tmp_len = strlen_se(bufptr0);
    bufptr = next_item(bufptr0);
    if (no_more_items_kns(bufptr)) {
      logprint(errstr_phenotype_format);
      logprint("At least two tokens expected in line.\n");
      sprintf(logbuf, "Original line: %s", bufptr0);
      logprintb();
      goto load_pheno_ret_INVALID_FORMAT;
    }
    tmp_len2 = strlen_se(bufptr);
    if (!header_processed) {
      if (phenoname_str || ((tmp_len == 3) && (tmp_len2 == 3) && (!memcmp("FID", bufptr0, 3)) && (!memcmp("IID", bufptr, 3)))) {
	if (phenoname_str) {
	  tmp_len = strlen(phenoname_str);
	  do {
	    bufptr = next_item(bufptr);
	    if (no_more_items_kns(bufptr)) {
	      logprint("Error: --pheno-name column not found.\n");
              goto load_pheno_ret_INVALID_FORMAT;
	    }
	    mpheno_col++;
	    tmp_len2 = strlen_se(bufptr);
	  } while ((tmp_len2 != tmp_len) || memcmp(bufptr, phenoname_str, tmp_len));
	}
      } else {
	header_processed = 1;
        if (!mpheno_col) {
	  mpheno_col = 1;
        }
      }
    }
    if (!header_processed) {
      header_processed = 1;
    } else {
      person_idx = bsearch_fam_indiv(tbuf, sorted_person_ids, max_person_id_len, indiv_ct, bufptr0, bufptr);
      if (person_idx != -1) {
	person_idx = id_map[(uint32_t)person_idx];
	bufptr = next_item_mult(bufptr, mpheno_col);
	if (no_more_items_kns(bufptr)) {
	  // not always an error
	  return LOAD_PHENO_LAST_COL;
	}
	if (affection) {
	  if (eval_affection(bufptr, missing_pheno, missing_pheno_len, affection_01)) {
	    if (is_missing_pheno(bufptr, missing_pheno, missing_pheno_len, affection_01)) {
	      // Since we're only making one pass through the file, we don't
	      // have the luxury of knowing in advance whether the phenotype is
	      // binary or scalar.  If there is a '0' entry that occurs before
	      // we know the phenotype is scalar, we need to not set the
	      // phenotype to zero during the binary -> scalar conversion step.
	      if (*bufptr == '0') {
		set_bit(isz, person_idx);
	      }
	      clear_bit(pheno_nm, person_idx);
	      clear_bit(pheno_c, person_idx);
	    } else {
	      if (*bufptr == case_char) {
		set_bit(pheno_c, person_idx);
	      } else {
		clear_bit(pheno_c, person_idx);
	      }
	      set_bit(pheno_nm, person_idx);
	    }
	  } else {
	    pheno_d = (double*)malloc(unfiltered_indiv_ct * sizeof(double));
	    if (!pheno_d) {
	      goto load_pheno_ret_NOMEM;
	    }
	    *pheno_d_ptr = pheno_d;
	    if (affection_01) {
	      dxx = 0.0;
	      dyy = 1.0;
	    } else {
	      dxx = 1.0;
	      dyy = 2.0;
	    }
	    for (uii = 0; uii < unfiltered_indiv_ct; uii++) {
	      if (IS_SET(isz, uii)) {
		pheno_d[uii] = 0.0;
		SET_BIT(pheno_nm, uii);
	      } else if (IS_SET(pheno_nm, uii)) {
		pheno_d[uii] = IS_SET(pheno_c, uii)? dyy : dxx;
	      }
	    }
	    free(pheno_c);
	    *pheno_c_ptr = NULL;
	    affection = 0;
	  }
	}
	if (!affection) {
	  if ((!scan_double(bufptr, &dxx)) && (dxx != missing_phenod)) {
	    pheno_d[(uint32_t)person_idx] = dxx;
	    set_bit(pheno_nm, person_idx);
	  }
	}
      }
    }
  }
  if (!feof(phenofile)) {
    goto load_pheno_ret_READ_FAIL;
  }
  while (0) {
  load_pheno_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  load_pheno_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  load_pheno_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
  wkspace_reset(wkspace_mark);
  return retval;
}

int32_t sort_item_ids_nx(char** sorted_ids_ptr, uint32_t** id_map_ptr, uintptr_t item_ct, char* item_ids, uintptr_t max_id_len) {
  // Version of sort_item_ids() with no exclusion.
  uintptr_t ulii;
  char* sorted_ids;
  char* dup_id;
  char* tptr;
  sorted_ids = (char*)wkspace_alloc(item_ct * max_id_len);
  if (!sorted_ids) {
    return RET_NOMEM;
  }
  *sorted_ids_ptr = sorted_ids;
  *id_map_ptr = (uint32_t*)wkspace_alloc(item_ct * sizeof(int32_t));
  if (!(*id_map_ptr)) {
    return RET_NOMEM;
  }
  if (!item_ct) {
    return 0;
  }
  for (ulii = 0; ulii < item_ct; ulii++) {
    memcpy(&(sorted_ids[ulii * max_id_len]), &(item_ids[ulii * max_id_len]), max_id_len);
    (*id_map_ptr)[ulii] = ulii;
  }
  if (qsort_ext(sorted_ids, item_ct, max_id_len, strcmp_deref, (char*)(*id_map_ptr), sizeof(int32_t))) {
    return RET_NOMEM;
  }
  dup_id = scan_for_duplicate_ids(sorted_ids, item_ct, max_id_len);
  if (dup_id) {
    tptr = strchr(dup_id, '\t');
    if (tptr) {
      *tptr = ' ';
    }
    sprintf(logbuf, "Error: Duplicate ID %s.\n", dup_id);
    logprintb();
    return RET_INVALID_FORMAT;
  }
  return 0;
}

#ifdef _WIN32
int32_t indiv_major_to_snp_major(char* indiv_major_fname, char* outname, uintptr_t unfiltered_marker_ct) {
  logprint("Error: Win32 " PROG_NAME_CAPS " does not yet support transposition of individual-major .bed\nfiles.  Contact the developers if you need this.\n");
  return RET_CALC_NOT_YET_SUPPORTED;
}
#else
int32_t indiv_major_to_snp_major(char* indiv_major_fname, char* outname, uintptr_t unfiltered_marker_ct) {
  // This implementation only handles large files on 64-bit Unix systems; a
  // more portable version needs to be written for Windows and 32-bit Unix.
  FILE* outfile = NULL;
  int32_t in_fd = open(indiv_major_fname, O_RDONLY);
  unsigned char* in_contents = (unsigned char*)MAP_FAILED;
  uintptr_t unfiltered_marker_ct4 = (unfiltered_marker_ct + 3) / 4;
  struct stat sb;
  unsigned char* icoff;
  int32_t retval;
  uintptr_t indiv_ct;
  uintptr_t indiv_ct4;
  uintptr_t indiv_ct4l;
  uintptr_t max_4blocks_in_mem;
  uintptr_t superblock_offset;
  uint32_t block_last_marker;
  uint32_t uii;
  uint32_t add_val;
  uint32_t rshift_val;
  uint32_t ujj;
  unsigned char* write_ptr;
  if (in_fd == -1) {
    sprintf(logbuf, errstr_fopen, indiv_major_fname);
    logprintb();
    return RET_OPEN_FAIL;
  }
  // obtain file size, see example in OS X mmap() documentation
  if (fstat(in_fd, &sb) == -1) {
    goto indiv_major_to_snp_major_ret_READ_FAIL;
  }
  in_contents = (unsigned char*)mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE, in_fd, 0);
  if (in_contents == MAP_FAILED) {
    goto indiv_major_to_snp_major_ret_READ_FAIL;
  }
  if (fopen_checked(&outfile, outname, "wb")) {
    goto indiv_major_to_snp_major_ret_OPEN_FAIL;
  }
  if ((in_contents[0] == 'l') && (in_contents[1] == '\x1b') && (!in_contents[2])) {
    uii = 3;
  } else if ((!in_contents[0]) && (!((sb.st_size - 1) % unfiltered_marker_ct4))) {
    uii = 1;
  } else {
    uii = 0;
  }
  icoff = &(in_contents[uii]);
  if ((sb.st_size - uii) % unfiltered_marker_ct4) {
    goto indiv_major_to_snp_major_ret_INVALID_FORMAT;
  }
  if (fwrite_checked("l\x1b\x01", 3, outfile)) {
    goto indiv_major_to_snp_major_ret_WRITE_FAIL;
  }
  indiv_ct = sb.st_size / unfiltered_marker_ct4;
  indiv_ct4l = indiv_ct / 4;
  indiv_ct4 = (indiv_ct + 3) / 4;
  // 4 * indiv_ct4 bytes needed per 4-marker block
  max_4blocks_in_mem = wkspace_left / (4 * indiv_ct4);
  superblock_offset = 0;
  while (superblock_offset < unfiltered_marker_ct4) {
    block_last_marker = unfiltered_marker_ct - (superblock_offset * 4);
    if (block_last_marker > (max_4blocks_in_mem * 4)) {
      block_last_marker = max_4blocks_in_mem * 4;
    }
    write_ptr = wkspace_base;
    for (uii = 0; uii < block_last_marker; uii++) {
      rshift_val = (uii % 4) * 2;
      add_val = uii / 4;
      for (ujj = 0; ujj < indiv_ct4l; ujj++) {
        *write_ptr++ = ((icoff[4 * ujj * unfiltered_marker_ct4 + add_val] >> rshift_val) & 3) + (((icoff[(4 * ujj + 1) * unfiltered_marker_ct4 + add_val] >> rshift_val) & 3) << 2) + (((icoff[(4 * ujj + 2) * unfiltered_marker_ct4 + add_val] >> rshift_val) & 3) << 4) + (((icoff[(4 * ujj + 3) * unfiltered_marker_ct4 + add_val] >> rshift_val) & 3) << 6);
      }
      if (indiv_ct % 4) {
	*write_ptr = 0;
	for (ujj = 0; ujj < (indiv_ct % 4); ujj++) {
	  *write_ptr |= ((icoff[(ujj + indiv_ct) * unfiltered_marker_ct4 + add_val] >> rshift_val) & 3) << (ujj * 2);
	}
	write_ptr++;
      }
    }
    if (fwrite_checked(wkspace_base, ((int64_t)block_last_marker) * indiv_ct4, outfile)) {
      goto indiv_major_to_snp_major_ret_WRITE_FAIL;
    }
    superblock_offset += max_4blocks_in_mem;
  }
  if (fclose_null(&outfile)) {
    goto indiv_major_to_snp_major_ret_WRITE_FAIL;
  }
  retval = 0;
  while (0) {
  indiv_major_to_snp_major_ret_INVALID_FORMAT:
    sprintf(logbuf, "Error: %s's file size is inconsistent with the variant count.\n", indiv_major_fname);
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  indiv_major_to_snp_major_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  indiv_major_to_snp_major_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  indiv_major_to_snp_major_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  }
  fclose_cond(outfile);
  if (in_contents != MAP_FAILED) {
    munmap(in_contents, sb.st_size);
  }
  close(in_fd);
  return retval;
}
#endif

const char errstr_map_format[] = "Error: Improperly formatted .map file.\n";

int32_t load_map(FILE** mapfile_ptr, char* mapname, uint32_t* map_cols_ptr, uintptr_t* unfiltered_marker_ct_ptr, uintptr_t* marker_exclude_ct_ptr, uintptr_t* max_marker_id_len_ptr, uintptr_t** marker_exclude_ptr, char** marker_ids_ptr, Chrom_info* chrom_info_ptr, uint32_t** marker_pos_ptr, uint32_t* map_is_unsorted_ptr) {
  // todo: some cleanup
  uintptr_t marker_exclude_ct = *marker_exclude_ct_ptr;
  uintptr_t max_marker_id_len = 0;
  uintptr_t unfiltered_marker_ct = 0;
  int32_t last_chrom = -1;
  int32_t last_pos = 0;
  int32_t marker_pos_needed = 0;
  int32_t chroms_encountered_m1 = -1;
  int32_t retval = 0;
  uintptr_t loaded_chrom_mask[CHROM_MASK_WORDS];
  uintptr_t* marker_exclude;
  uintptr_t unfiltered_marker_ctl;
  char* bufptr;
  uintptr_t marker_uidx;
  uintptr_t ulii;
  int32_t jj;
  int32_t cur_pos;
  fill_ulong_zero(loaded_chrom_mask, CHROM_MASK_WORDS);
  if (fopen_checked(mapfile_ptr, mapname, "r")) {
    goto load_map_ret_OPEN_FAIL;
  }
  // first pass: count columns, determine raw marker count, determine maximum
  // marker ID length if necessary.
  tbuf[MAXLINELEN - 6] = ' ';
  while (fgets(tbuf, MAXLINELEN - 5, *mapfile_ptr) != NULL) {
    if (!tbuf[MAXLINELEN - 6]) {
      sprintf(logbuf, "Error: Excessively long line in .map file (max %u chars).\n", MAXLINELEN - 8);
      logprintb();
      goto load_map_ret_INVALID_FORMAT;
    }
    bufptr = skip_initial_spaces(tbuf);
    if (is_eoln_or_comment(*bufptr)) {
      continue;
    }
    bufptr = next_item(bufptr);
    if (no_more_items_kns(bufptr)) {
      logprint(errstr_map_format);
      goto load_map_ret_INVALID_FORMAT;
    }
    ulii = strlen_se(bufptr) + 1;
    if (ulii > max_marker_id_len) {
      max_marker_id_len = ulii;
    }
    if (!unfiltered_marker_ct) {
      bufptr = next_item_mult(bufptr, 2);
      if (!bufptr) {
	logprint(errstr_map_format);
        goto load_map_ret_INVALID_FORMAT;
      }
      if (*bufptr > ' ') {
	*map_cols_ptr = 4;
      }
    }
    unfiltered_marker_ct++;
  }
  if (!feof(*mapfile_ptr)) {
    goto load_map_ret_READ_FAIL;
  }
  if (!unfiltered_marker_ct) {
    logprint("Error: No variants in .map file.");
    goto load_map_ret_INVALID_FORMAT;
  }
  *unfiltered_marker_ct_ptr = unfiltered_marker_ct;
  *max_marker_id_len_ptr = max_marker_id_len;
  rewind(*mapfile_ptr);
  unfiltered_marker_ctl = (unfiltered_marker_ct + (BITCT - 1)) / BITCT;

  // unfiltered_marker_ct can be very large, so use wkspace for all allocations
  // that are a multiple of it

  // permanent stack allocation #1: marker_exclude
  if (wkspace_alloc_ul_checked(marker_exclude_ptr, unfiltered_marker_ctl * sizeof(intptr_t))) {
    goto load_map_ret_NOMEM;
  }
  marker_exclude = *marker_exclude_ptr;
  fill_ulong_zero(marker_exclude, unfiltered_marker_ctl);
  fill_uint_zero(chrom_info_ptr->chrom_file_order, MAX_POSSIBLE_CHROM);
  fill_uint_zero(chrom_info_ptr->chrom_file_order_marker_idx, MAX_POSSIBLE_CHROM + 1);
  fill_uint_zero(chrom_info_ptr->chrom_start, MAX_POSSIBLE_CHROM);
  fill_uint_zero(chrom_info_ptr->chrom_end, MAX_POSSIBLE_CHROM);
  // permanent stack allocation #3, if needed: marker_pos
  if (marker_pos_needed) {
    *marker_pos_ptr = (uint32_t*)wkspace_alloc(unfiltered_marker_ct * sizeof(int32_t));
    if (!(*marker_pos_ptr)) {
      goto load_map_ret_NOMEM;
    }
  }
  *marker_ids_ptr = (char*)wkspace_alloc(unfiltered_marker_ct * max_marker_id_len);
  if (!(*marker_ids_ptr)) {
    goto load_map_ret_NOMEM;
  }

  // second pass: actually load stuff
  for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
    if (get_next_noncomment(*mapfile_ptr, &bufptr)) {
      goto load_map_ret_READ_FAIL;
    }
    jj = get_chrom_code(chrom_info_ptr, bufptr);
    if (jj == -1) {
      logprint("Error: Invalid chromosome code in .map file.\n");
      goto load_map_ret_INVALID_FORMAT;
    }
    if (jj != last_chrom) {
      if (last_chrom != -1) {
	chrom_info_ptr->chrom_end[last_chrom] = marker_uidx;
      }
      if (jj < last_chrom) {
	*map_is_unsorted_ptr |= UNSORTED_CHROM;
      }
      last_chrom = jj;
      if (is_set(loaded_chrom_mask, jj)) {
	*map_is_unsorted_ptr |= UNSORTED_SPLIT_CHROM;
      } else {
	set_bit(loaded_chrom_mask, jj);
	chrom_info_ptr->chrom_start[(uint32_t)jj] = marker_uidx;
	chrom_info_ptr->chrom_file_order[++chroms_encountered_m1] = jj;
	chrom_info_ptr->chrom_file_order_marker_idx[chroms_encountered_m1] = marker_uidx;
      }
      last_pos = 0;
    }

    if (!is_set(chrom_info_ptr->chrom_mask, jj)) {
      SET_BIT(marker_exclude, marker_uidx);
      marker_exclude_ct++;
    } else {
      bufptr = next_item(bufptr);
      if (no_more_items_kns(bufptr)) {
	logprint(errstr_map_format);
	goto load_map_ret_INVALID_FORMAT;
      }
      read_next_terminate(&((*marker_ids_ptr)[marker_uidx * max_marker_id_len]), bufptr);
      bufptr = next_item_mult(bufptr, *map_cols_ptr - 2);
      if (no_more_items_kns(bufptr)) {
        logprint(errstr_map_format);
	goto load_map_ret_INVALID_FORMAT;
      }
      if (*bufptr == '-') {
	// negative marker positions now have same effect in .bim as .map
	SET_BIT(marker_exclude, marker_uidx);
	marker_exclude_ct++;
      } else {
	cur_pos = atoi(bufptr);
	if (cur_pos < last_pos) {
	  *map_is_unsorted_ptr |= UNSORTED_BP;
	} else {
	  last_pos = cur_pos;
	}
	if (marker_pos_needed && jj) {
	  (*marker_pos_ptr)[marker_uidx] = cur_pos;
	}
      }
    }
  }
  chrom_info_ptr->chrom_end[last_chrom] = marker_uidx;
  chrom_info_ptr->chrom_ct = ++chroms_encountered_m1;
  chrom_info_ptr->chrom_file_order_marker_idx[chroms_encountered_m1] = marker_uidx;
  *marker_exclude_ct_ptr = marker_exclude_ct;
  if (*marker_exclude_ct_ptr == unfiltered_marker_ct) {
    logprint("Error: All variants excluded from .map file.\n");
    goto load_map_ret_ALL_MARKERS_EXCLUDED;
  }
  while (0) {
  load_map_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  load_map_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  load_map_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  load_map_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  load_map_ret_ALL_MARKERS_EXCLUDED:
    retval = RET_ALL_MARKERS_EXCLUDED;
    break;
  }
  return retval;
}

uint32_t load_bim_scan_position(uint32_t* pos_ptr, char* bufptr, uint32_t mcm2) {
  int32_t ii;
  bufptr = next_item_mult(bufptr, mcm2);
  if (no_more_items_kns(bufptr)) {
    return -1;
  }
  ii = atoi(bufptr);
  if ((ii < 1) && ((bufptr[0] != '0') || ((bufptr[1] != ' ') && (bufptr[1] != '\t')))) {
    return -1;
  }
  *pos_ptr = ii;
  return 0;
}

void load_bim_sf_insert(uint32_t chrom_idx, uint32_t pos_start, uint32_t pos_end, uint32_t* start_idxs, uint32_t* llbuf, uint32_t* lltop_ptr, uint32_t* entry_ct_ptr) {
  uint32_t lltop = *lltop_ptr;
  uint32_t entry_ct = *entry_ct_ptr;
  uint32_t llidx;
  uint32_t new_start;
  uint32_t new_end;
  uint32_t new_llidx;
  uint32_t old_llidx;
  if (start_idxs[chrom_idx] == 1) {
    start_idxs[chrom_idx] = lltop;
    llbuf[lltop++] = pos_start;
    llbuf[lltop++] = pos_end;
    llbuf[lltop++] = 1;
    entry_ct++;
  } else {
    llidx = start_idxs[chrom_idx];
    while (1) {
      if (llbuf[llidx] > pos_end) {
	if (llbuf[llidx] == pos_end + 1) {
	  llbuf[llidx] = pos_start;
	} else {
	  new_llidx = llidx;
	  do {
	    llidx = new_llidx;
	    new_start = llbuf[llidx];
	    llbuf[llidx] = pos_start;
	    pos_start = new_start;
	    new_end = llbuf[llidx + 1];
	    llbuf[llidx + 1] = pos_end;
	    pos_end = new_end;
	    new_llidx = llbuf[llidx + 2];
	  } while (new_llidx != 1);
	  llbuf[llidx + 2] = lltop;
	  llbuf[lltop++] = pos_start;
	  llbuf[lltop++] = pos_end;
	  llbuf[lltop++] = 1;
	  entry_ct++;
	}
	break;
      } else if (llbuf[llidx + 1] + 1 >= pos_start) {
	// mergeable
	if (llbuf[llidx] > pos_start) {
	  llbuf[llidx] = pos_start;
	}
	if (llbuf[llidx + 1] < pos_end) {
	  // scan forward, attempt to collapse entries
	  old_llidx = llidx;
          new_llidx = llbuf[llidx + 2];
	  while (new_llidx != 1) {
	    llidx = new_llidx;
	    if (llbuf[llidx] > pos_end + 1) {
	      break;
	    }
	    entry_ct--;
	    new_llidx = llbuf[llidx + 2];
	    llbuf[old_llidx + 2] = new_llidx;
	    if (llbuf[llidx + 1] >= pos_end) {
	      llbuf[old_llidx + 1] = llbuf[llidx + 1];
	      break;
	    }
	  }
	}
	break;
      }
      new_llidx = llbuf[llidx + 2];
      if (new_llidx == 1) {
	llbuf[llidx + 2] = lltop;
	llbuf[lltop++] = pos_start;
	llbuf[lltop++] = pos_end;
	llbuf[lltop++] = 1;
	entry_ct++;
	break;
      }
      llidx = new_llidx;
    }
  }
  *lltop_ptr = lltop;
  *entry_ct_ptr = entry_ct;
}

static inline uint32_t sf_out_of_range(uint32_t cur_pos, uint32_t chrom_idx, uint32_t* sf_start_idxs, uint32_t* sf_pos) {
  uint32_t cur_idx = sf_start_idxs[chrom_idx];
  uint32_t end_idx = sf_start_idxs[chrom_idx + 1];
  while (cur_idx < end_idx) {
    if ((cur_pos >= sf_pos[cur_idx]) && (cur_pos <= sf_pos[cur_idx + 1])) {
      return 0;
    }
    cur_idx += 2;
  }
  return 1;
}

int32_t load_bim(char* bimname, uint32_t* map_cols_ptr, uintptr_t* unfiltered_marker_ct_ptr, uintptr_t* marker_exclude_ct_ptr, uintptr_t* max_marker_id_len_ptr, uintptr_t** marker_exclude_ptr, double** set_allele_freqs_ptr, char*** marker_allele_pp, uintptr_t* max_marker_allele_len_ptr, char** marker_ids_ptr, uint32_t allow_extra_chroms, Chrom_info* chrom_info_ptr, double** marker_cms_ptr, uint32_t** marker_pos_ptr, char* freqname, uint64_t calculation_type, uint32_t recode_modifier, int32_t marker_pos_start, int32_t marker_pos_end, uint32_t snp_window_size, char* markername_from, char* markername_to, char* markername_snp, uint32_t exclude_snp, Range_list* sf_range_list_ptr, uint32_t* map_is_unsorted_ptr, uint32_t marker_pos_needed, uint32_t marker_cms_needed, uint32_t marker_alleles_needed, const char* extension, const char* split_chrom_cmd) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* bimfile = NULL;
  uintptr_t unfiltered_marker_ct = 0;
  uintptr_t marker_exclude_ct = *marker_exclude_ct_ptr;
  uintptr_t max_marker_id_len = *max_marker_id_len_ptr;
  uintptr_t max_marker_allele_len = *max_marker_allele_len_ptr;
  int32_t prev_chrom = -1;
  uint32_t last_pos = 0;
  uint32_t from_slen = markername_from? strlen(markername_from) : 0;
  uint32_t to_slen = markername_to? strlen(markername_to) : 0;
  uint32_t snp_slen = markername_snp? strlen(markername_snp) : 0;
  uint32_t sf_ct = sf_range_list_ptr->name_ct;
  uint32_t sf_max_len = sf_range_list_ptr->name_max_len;
  uint32_t slen_check = from_slen || to_slen || snp_slen || sf_ct;
  uint32_t from_chrom = MAX_POSSIBLE_CHROM;
  uint32_t to_chrom = MAX_POSSIBLE_CHROM;
  uint32_t snp_chrom = MAX_POSSIBLE_CHROM;
  uint32_t chroms_encountered_m1 = 0xffffffffU; // intentional overflow
  uint32_t snp_pos = 0;
  uint32_t mcm2 = 1;
  uint32_t split_chrom = 0;
  uint32_t max_loadbuf = 0;
  int32_t exclude_window_start = 0;
  int32_t exclude_window_end = -1;
  char* missing_geno_ptr = (char*)g_missing_geno_ptr;
  uint32_t* sf_pos = NULL;
  uint32_t* sf_str_chroms = NULL;
  uint32_t* sf_str_pos = NULL;
  uint32_t* sf_str_lens = NULL;
  uint32_t* sf_llbuf = NULL;
  char* loadbuf2 = NULL; // on heap, second pass
  char** marker_allele_ptrs = NULL;
  uintptr_t loaded_chrom_mask[CHROM_MASK_WORDS];
  uintptr_t sf_mask[CHROM_MASK_WORDS];
  uint32_t sf_start_idxs[MAX_POSSIBLE_CHROM + 1];
  uintptr_t* marker_exclude;
  char* loadbuf; // on stack, first pass
  uintptr_t loadbuf_size;
  uintptr_t unfiltered_marker_ctl;
  uint32_t sf_entry_ct;
  uint32_t sf_lltop;
  char* bufptr;
  char* bufptr2;
  char* bufptr3;
  uintptr_t ulii;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  uint32_t unn;
  uint32_t uoo;
  int32_t jj;
  uint32_t cur_pos;
  uintptr_t marker_uidx;
  int32_t retval;
  fill_ulong_zero(loaded_chrom_mask, CHROM_MASK_WORDS);
  if (sf_ct) {
    if (wkspace_alloc_ui_checked(&sf_str_chroms, sf_ct * sizeof(int32_t)) ||
	wkspace_alloc_ui_checked(&sf_str_pos, sf_ct * sizeof(int32_t)) ||
        wkspace_alloc_ui_checked(&sf_str_lens, sf_ct * sizeof(int32_t)) ||
        wkspace_alloc_ui_checked(&sf_llbuf, 3 * (MAX_POSSIBLE_CHROM + sf_ct) * sizeof(int32_t))) {
      goto load_bim_ret_NOMEM;
    }
    for (uii = 0; uii < sf_ct; uii++) {
      sf_str_chroms[uii] = MAX_POSSIBLE_CHROM;
      sf_str_lens[uii] = strlen(&(sf_range_list_ptr->names[uii * sf_max_len]));
    }
  }
  if (fopen_checked(&bimfile, bimname, "r")) {
    goto load_bim_ret_OPEN_FAIL;
  }
  // first pass: count columns, determine raw marker count, determine maximum
  // marker ID length and/or marker allele length if necessary, save
  // nonstandard chromosome names.
  loadbuf_size = wkspace_left;
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  } else if (loadbuf_size <= MAXLINELEN) {
    goto load_bim_ret_NOMEM;
  }
  loadbuf = (char*)wkspace_base;
  loadbuf[loadbuf_size - 1] = ' ';
  while (fgets(loadbuf, loadbuf_size, bimfile) != NULL) {
    if (!loadbuf[loadbuf_size - 1]) {
      if (loadbuf_size == MAXLINEBUFLEN) {
        logprint("Error: Pathologically long line in .bim file.\n");
        goto load_bim_ret_INVALID_FORMAT;
      } else {
	goto load_bim_ret_NOMEM;
      }
    }
    uii = strlen(loadbuf);
    if (uii >= max_loadbuf) {
      max_loadbuf = uii + 1;
    }
    // bufptr = col 1 start
    bufptr = skip_initial_spaces(loadbuf);
    if (is_eoln_or_comment(*bufptr)) {
      continue;
    }
    jj = get_chrom_code(chrom_info_ptr, bufptr);
    if (jj == -1) {
      if (!allow_extra_chroms) {
	uii = strlen_se(bufptr);
	bufptr[uii] = '\0';
	sprintf(logbuf, "Error: Invalid chromosome code '%s' in .bim file.\n(Use --allow-extra-chr to force it to be accepted.)\n", bufptr);
	goto load_bim_ret_INVALID_FORMAT_2;
      }
      retval = resolve_or_add_chrom_name(chrom_info_ptr, bufptr, &jj);
      if (retval) {
	goto load_bim_ret_1;
      }
    }

    // bufptr = col 2 start
    bufptr = next_item(bufptr);
    if (no_more_items_kns(bufptr)) {
      goto load_bim_ret_INVALID_FORMAT_5;
    }
    ulii = strlen_se(bufptr) + 1;
    if (ulii > max_marker_id_len) {
      max_marker_id_len = ulii;
    }
    if (!unfiltered_marker_ct) {
      if (!strcmp(extension, "cnv.map")) {
	bufptr2 = next_item(bufptr);
      } else {
	// bufptr2 = col 5 start
	bufptr2 = next_item_mult(bufptr, 3);
      }
      if (no_more_items_kns(bufptr2)) {
	goto load_bim_ret_INVALID_FORMAT_5;
      }
      // check if col 6 exists
      if (*(skip_initial_spaces(item_endnn(bufptr2))) > ' ') {
	*map_cols_ptr = 4;
	mcm2 = 2;
      }
    }
    if (marker_alleles_needed) {
      bufptr2 = next_item_mult(bufptr, mcm2 + 1);
      bufptr3 = next_item(bufptr2);
      if (no_more_items_kns(bufptr3)) {
	goto load_bim_ret_INVALID_FORMAT_5;
      }
      uii = strlen_se(bufptr2);
      if (uii >= max_marker_allele_len) {
	max_marker_allele_len = uii + 1;
      }
      uii = strlen_se(bufptr3);
      if (uii >= max_marker_allele_len) {
	max_marker_allele_len = uii + 1;
      }
    }
    if (slen_check) {
      ulii--;
      if (sf_ct) {
	uii = 0;
	do {
	  if ((ulii == sf_str_lens[uii]) && (!memcmp(bufptr, &(sf_range_list_ptr->names[uii * sf_max_len]), ulii))) {
	    if (sf_str_chroms[uii] != MAX_POSSIBLE_CHROM) {
	      goto load_bim_ret_INVALID_FORMAT_3;
	    }
	    sf_str_chroms[uii] = jj;
	    if (load_bim_scan_position(&(sf_str_pos[uii]), bufptr, mcm2)) {
	      goto load_bim_ret_INVALID_FORMAT;
	    }
	    break;
	  }
	} while (++uii < sf_ct);
      } else {
	if ((ulii == from_slen) && (!memcmp(bufptr, markername_from, ulii))) {
	  if (from_chrom != MAX_POSSIBLE_CHROM) {
	    goto load_bim_ret_INVALID_FORMAT_3;
	  }
	  if (load_bim_scan_position(&uii, bufptr, mcm2)) {
	    goto load_bim_ret_INVALID_FORMAT;
	  }
	  marker_pos_start = uii;
	  from_chrom = jj;
	  if (to_chrom != MAX_POSSIBLE_CHROM) {
	    if (from_chrom != to_chrom) {
	      goto load_bim_ret_INVALID_FORMAT_4;
	    }
	  }
	  fill_ulong_zero(chrom_info_ptr->chrom_mask, CHROM_MASK_WORDS);
	  SET_BIT(chrom_info_ptr->chrom_mask, from_chrom);
	}
	if ((ulii == to_slen) && (!memcmp(bufptr, markername_to, ulii))) {
	  if (to_chrom != MAX_POSSIBLE_CHROM) {
	    goto load_bim_ret_INVALID_FORMAT_3;
	  }
	  if (load_bim_scan_position(&uii, bufptr, mcm2)) {
	    goto load_bim_ret_INVALID_FORMAT;
	  }
	  marker_pos_end = uii;
	  to_chrom = jj;
	  if (from_chrom != MAX_POSSIBLE_CHROM) {
	    if (to_chrom != from_chrom) {
	      goto load_bim_ret_INVALID_FORMAT_4;
	    }
	  }
	  fill_ulong_zero(chrom_info_ptr->chrom_mask, CHROM_MASK_WORDS);
	  SET_BIT(chrom_info_ptr->chrom_mask, to_chrom);
	}
	if ((ulii == snp_slen) && (!memcmp(bufptr, markername_snp, ulii))) {
	  if (snp_chrom != MAX_POSSIBLE_CHROM) {
	    goto load_bim_ret_INVALID_FORMAT_3;
	  }
	  if (load_bim_scan_position(&snp_pos, bufptr, mcm2)) {
	    goto load_bim_ret_INVALID_FORMAT;
	  }
	  snp_chrom = jj;
	  if (!exclude_snp) {
	    fill_ulong_zero(chrom_info_ptr->chrom_mask, CHROM_MASK_WORDS);
	    SET_BIT(chrom_info_ptr->chrom_mask, snp_chrom);
	  }
	}
      }
    }
    unfiltered_marker_ct++;
  }
  if (sf_ct) {
    for (uii = 0; uii < sf_ct; uii++) {
      if (sf_str_chroms[uii] == MAX_POSSIBLE_CHROM) {
	sprintf(logbuf, "Error: Variant '%s' not found in .%s file.\n", &(sf_range_list_ptr->names[uii * sf_max_len]), extension);
	goto load_bim_ret_INVALID_FORMAT_2;
      }
    }
    // effectively build out one linked list per chromosome
    memcpy(sf_mask, chrom_info_ptr->chrom_mask, CHROM_MASK_WORDS * sizeof(intptr_t));
    sf_entry_ct = 0;
    sf_lltop = 0;
    ujj = chrom_info_ptr->max_code + chrom_info_ptr->name_ct;
    for (uii = 0; uii <= ujj; uii++) {
      sf_start_idxs[uii] = 1; // impossible (multiples of 3)
    }
    uii = 0;
    do {
      ujj = sf_str_chroms[uii];
      ukk = sf_str_pos[uii];
      if (sf_range_list_ptr->starts_range[uii]) {
	umm = sf_str_chroms[uii + 1];
	unn = sf_str_pos[uii + 1];
	if (ujj != umm) {
	  if (ujj > umm) {
	    uoo = ujj;
	    ujj = umm;
	    umm = uoo;
	    uoo = ukk;
	    ukk = unn;
	    unn = uoo;
	  }
	  if (IS_SET(sf_mask, ujj)) {
	    load_bim_sf_insert(ujj, ukk, 0x7fffffff, sf_start_idxs, sf_llbuf, &sf_lltop, &sf_entry_ct);
	  }
	  for (uoo = ujj + 1; uoo < umm; uoo++) {
	    if (IS_SET(sf_mask, uoo)) {
	      load_bim_sf_insert(uoo, 0, 0x7fffffff, sf_start_idxs, sf_llbuf, &sf_lltop, &sf_entry_ct);
	    }
	  }
	  if (IS_SET(sf_mask, umm)) {
	    load_bim_sf_insert(umm, 0, unn, sf_start_idxs, sf_llbuf, &sf_lltop, &sf_entry_ct);
	  }
	} else {
	  if (ukk > unn) {
            umm = ukk;
	    ukk = unn;
	    unn = umm;
	  }
	  if (IS_SET(sf_mask, ujj)) {
	    load_bim_sf_insert(ujj, ukk, unn, sf_start_idxs, sf_llbuf, &sf_lltop, &sf_entry_ct);
	  }
	}
	uii += 2;
      } else {
	if (IS_SET(sf_mask, ujj)) {
	  load_bim_sf_insert(ujj, ukk, ukk, sf_start_idxs, sf_llbuf, &sf_lltop, &sf_entry_ct);
	}
	uii++;
      }
    } while (uii < sf_ct);
    // now compactify
    sf_pos = (uint32_t*)malloc(sf_entry_ct * 2 * sizeof(int32_t));
    if (!sf_pos) {
      goto load_bim_ret_NOMEM;
    }
    ujj = chrom_info_ptr->max_code + chrom_info_ptr->name_ct;
    ukk = 0;
    for (uii = 0; uii <= ujj; uii++) {
      if (sf_start_idxs[uii] == 1) {
	CLEAR_BIT(sf_mask, uii);
	sf_start_idxs[uii] = ukk;
	continue;
      }
      umm = sf_start_idxs[uii];
      sf_start_idxs[uii] = ukk;
      do {
	sf_pos[ukk++] = sf_llbuf[umm];
	sf_pos[ukk++] = sf_llbuf[umm + 1];
	umm = sf_llbuf[umm + 2];
      } while (umm != 1);
    }
    sf_start_idxs[ujj + 1] = ukk;
    if (!exclude_snp) {
      memcpy(chrom_info_ptr->chrom_mask, sf_mask, CHROM_MASK_WORDS * sizeof(intptr_t));
    }
    wkspace_reset(wkspace_mark);
  }
  if (!feof(bimfile)) {
    goto load_bim_ret_READ_FAIL;
  }
  if (!unfiltered_marker_ct) {
    sprintf(logbuf, "Error: No variants in .%s file.", extension);
    goto load_bim_ret_INVALID_FORMAT_2;
  }
  if (from_slen || to_slen) {
    if (from_slen && (from_chrom == MAX_POSSIBLE_CHROM)) {
      sprintf(logbuf, "Error: --from variant '%s' not found.\n", markername_from);
      goto load_bim_ret_INVALID_FORMAT_2;
    }
    if (to_slen && (to_chrom == MAX_POSSIBLE_CHROM)) {
      sprintf(logbuf, "Error: --to variant '%s' not found.\n", markername_to);
      goto load_bim_ret_INVALID_FORMAT_2;
    }
    if (marker_pos_start == -1) {
      marker_pos_start = 0;
    }
    if (marker_pos_end == -1) {
      marker_pos_end = 0x7fffffff;
    }
    if (marker_pos_start > marker_pos_end) {
      jj = marker_pos_start;
      marker_pos_start = marker_pos_end;
      marker_pos_end = jj;
    }
  }
  if (snp_slen) {
    if (snp_chrom == MAX_POSSIBLE_CHROM) {
      sprintf(logbuf, "Error: --%ssnp variant '%s' not found.\n", exclude_snp? "exclude-" : "", markername_snp);
      goto load_bim_ret_INVALID_FORMAT_2;
    }
    if (!exclude_snp) {
      if (snp_window_size > snp_pos) {
	marker_pos_start = 0;
      } else {
	marker_pos_start = snp_pos - snp_window_size;
      }
      if (snp_window_size > (0x7fffffff - snp_pos)) {
	marker_pos_end = 0x7fffffff;
      } else {
	marker_pos_end = snp_pos + snp_window_size;
      }
    } else {
      if (snp_window_size <= snp_pos) {
        exclude_window_start = snp_pos - snp_window_size;
      }
      if (snp_window_size > (0x7fffffff - snp_pos)) {
	exclude_window_end = 0x7fffffff;
      } else {
        exclude_window_end = snp_pos + snp_window_size;
      }
    }
  }

  *unfiltered_marker_ct_ptr = unfiltered_marker_ct;
  *max_marker_id_len_ptr = max_marker_id_len;
  rewind(bimfile);
  unfiltered_marker_ctl = (unfiltered_marker_ct + (BITCT - 1)) / BITCT;

  // unfiltered_marker_ct can be very large, so use wkspace for all allocations
  // that are a multiple of it

  // permanent stack allocation #1: marker_exclude
  // permanent stack allocation #2: set_allele_freqs
  if (wkspace_alloc_ul_checked(marker_exclude_ptr, unfiltered_marker_ctl * sizeof(intptr_t)) ||
      wkspace_alloc_d_checked(set_allele_freqs_ptr, unfiltered_marker_ct * sizeof(double))) {
    goto load_bim_ret_NOMEM;
  }
  marker_exclude = *marker_exclude_ptr;
  fill_ulong_zero(marker_exclude, unfiltered_marker_ctl);
  for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
    (*set_allele_freqs_ptr)[marker_uidx] = -1.0;
  }
  fill_uint_zero(chrom_info_ptr->chrom_file_order, MAX_POSSIBLE_CHROM);
  fill_uint_zero(chrom_info_ptr->chrom_file_order_marker_idx, MAX_POSSIBLE_CHROM + 1);
  fill_uint_zero(chrom_info_ptr->chrom_start, MAX_POSSIBLE_CHROM);
  fill_uint_zero(chrom_info_ptr->chrom_end, MAX_POSSIBLE_CHROM);
  // permanent stack allocation #3, if needed: marker_pos
  if (marker_pos_needed) {
    if (wkspace_alloc_ui_checked(marker_pos_ptr, unfiltered_marker_ct * sizeof(int32_t))) {
      goto load_bim_ret_NOMEM;
    }
  }
  if (marker_alleles_needed) {
    *max_marker_allele_len_ptr = max_marker_allele_len;
    marker_allele_ptrs = (char**)wkspace_alloc(unfiltered_marker_ct * 2 * sizeof(intptr_t));
    if (!marker_allele_ptrs) {
      goto load_bim_ret_NOMEM;
    }
    *marker_allele_pp = marker_allele_ptrs;
    ujj = unfiltered_marker_ct * 2;
    for (uii = 0; uii < ujj; uii++) {
      marker_allele_ptrs[uii] = missing_geno_ptr;
    }
  }
  if (wkspace_alloc_c_checked(marker_ids_ptr, unfiltered_marker_ct * max_marker_id_len)) {
    goto load_bim_ret_NOMEM;
  }
  // todo: check whether marker_cms can be unloaded before
  // marker_ids/marker_alleles, or vice versa
  if (marker_cms_needed & MARKER_CMS_FORCED) {
    if (wkspace_alloc_d_checked(marker_cms_ptr, unfiltered_marker_ct * sizeof(double))) {
      goto load_bim_ret_NOMEM;
    }
    fill_double_zero(*marker_cms_ptr, unfiltered_marker_ct);
  }

  // second pass: actually load stuff
  loadbuf2 = (char*)malloc(max_loadbuf);
  if (!loadbuf2) {
    goto load_bim_ret_NOMEM;
  }
  for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
    do {
      if (!fgets(loadbuf2, max_loadbuf, bimfile)) {
	goto load_bim_ret_READ_FAIL;
      }
      bufptr = skip_initial_spaces(loadbuf2);
    } while (is_eoln_or_comment(*bufptr));
    jj = get_chrom_code(chrom_info_ptr, bufptr);
    if (jj != prev_chrom) {
      if (!split_chrom) {
	if (prev_chrom != -1) {
	  chrom_info_ptr->chrom_end[(uint32_t)prev_chrom] = marker_uidx;
	}
	if (jj < prev_chrom) {
	  *map_is_unsorted_ptr |= UNSORTED_CHROM;
	}
	prev_chrom = jj;
	if (is_set(loaded_chrom_mask, jj)) {
	  if (split_chrom_cmd) {
	    sprintf(logbuf, "Error: .%s file has a split chromosome.  Use --%s by itself to\nremedy this.\n", extension, split_chrom_cmd);
	    goto load_bim_ret_INVALID_FORMAT;
	  }
	  split_chrom = 1;
	  *map_is_unsorted_ptr = UNSORTED_CHROM | UNSORTED_BP | UNSORTED_SPLIT_CHROM;
	} else {
	  chrom_info_ptr->chrom_start[(uint32_t)jj] = marker_uidx;
	  chrom_info_ptr->chrom_file_order[++chroms_encountered_m1] = jj;
	  chrom_info_ptr->chrom_file_order_marker_idx[chroms_encountered_m1] = marker_uidx;
	}
        last_pos = 0;
      }
      set_bit(loaded_chrom_mask, jj);
    }

    if (is_set(chrom_info_ptr->chrom_mask, jj)) {
      bufptr = next_item(bufptr);
      if (no_more_items_kns(bufptr)) {
	goto load_bim_ret_INVALID_FORMAT_5;
      }
      read_next_terminate(&((*marker_ids_ptr)[marker_uidx * max_marker_id_len]), bufptr);
      if (marker_cms_needed) {
	bufptr = next_item(bufptr);
	if (no_more_items_kns(bufptr)) {
	  goto load_bim_ret_INVALID_FORMAT_5;
	}
	if ((*bufptr != '0') || (bufptr[1] > ' ')) {
	  if (!(*marker_cms_ptr)) {
	    if (wkspace_alloc_d_checked(marker_cms_ptr, unfiltered_marker_ct * sizeof(double))) {
	      goto load_bim_ret_NOMEM;
	    }
	    fill_double_zero(*marker_cms_ptr, unfiltered_marker_ct);
	  }
	  if (scan_double(bufptr, &((*marker_cms_ptr)[marker_uidx]))) {
	    sprintf(logbuf, "Error: Invalid centimorgan position in .%s file.\n", extension);
	    goto load_bim_ret_INVALID_FORMAT_2;
	  }
	}
	bufptr = next_item(bufptr);
      } else {
        bufptr = next_item_mult(bufptr, mcm2);
      }
      if (no_more_items_kns(bufptr)) {
	goto load_bim_ret_INVALID_FORMAT_5;
      }
      if (*bufptr == '-') {
	// negative marker positions now have same effect in .bim as .map
	goto load_bim_skip_marker;
      }
      cur_pos = atoi(bufptr);
      if (cur_pos < last_pos) {
	*map_is_unsorted_ptr |= UNSORTED_BP;
      } else {
	last_pos = cur_pos;
      }
      if ((sf_ct && (exclude_snp ^ sf_out_of_range(cur_pos, (uint32_t)jj, sf_start_idxs, sf_pos))) || ((marker_pos_start != -1) && ((((int32_t)cur_pos) < marker_pos_start) || (((int32_t)cur_pos) > marker_pos_end))) || (exclude_snp && (((int32_t)cur_pos) <= exclude_window_end) && (((int32_t)cur_pos) >= exclude_window_start) && ((uint32_t)jj == snp_chrom))) {
	goto load_bim_skip_marker;
      }
      if (marker_pos_needed) {
	(*marker_pos_ptr)[marker_uidx] = cur_pos;
      }
      if (marker_alleles_needed) {
	bufptr = next_item(bufptr);
	bufptr2 = next_item(bufptr);
	uii = strlen_se(bufptr);
	ulii = marker_uidx * 2;
        if (allele_set(&(marker_allele_ptrs[ulii]), bufptr, uii)) {
	  goto load_bim_ret_NOMEM;
	}
	uii = strlen_se(bufptr2);
	ulii++;
	if (allele_set(&(marker_allele_ptrs[ulii]), bufptr2, uii)) {
	  goto load_bim_ret_NOMEM;
	}
      }
    } else {
    load_bim_skip_marker:
      SET_BIT(marker_exclude, marker_uidx);
      marker_exclude_ct++;
      if (marker_pos_needed) {
        // support unfiltered marker_pos search
        (*marker_pos_ptr)[marker_uidx] = last_pos;
      }
    }
  }
  if (unfiltered_marker_ct == marker_exclude_ct) {
    logprint("Error: All variants excluded.\n");
    goto load_bim_ret_ALL_MARKERS_EXCLUDED;
  }
  for (uii = 0; uii < CHROM_MASK_WORDS; uii++) {
    chrom_info_ptr->chrom_mask[uii] &= loaded_chrom_mask[uii];
  }
  chrom_info_ptr->chrom_end[prev_chrom] = marker_uidx;
  chrom_info_ptr->chrom_ct = ++chroms_encountered_m1;
  chrom_info_ptr->chrom_file_order_marker_idx[chroms_encountered_m1] = marker_uidx;
  *marker_exclude_ct_ptr = marker_exclude_ct;
  retval = 0;
  while (0) {
  load_bim_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  load_bim_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  load_bim_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  load_bim_ret_INVALID_FORMAT_5:
    sprintf(logbuf, "Error: Improperly formatted .%s file.\n", extension);
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  load_bim_ret_INVALID_FORMAT_4:
    logprint("Error: --from and --to variants are not on the same chromosome.\n");
    retval = RET_INVALID_FORMAT;
    break;
  load_bim_ret_INVALID_FORMAT_3:
    uii = strlen_se(bufptr);
    bufptr[uii] = '\0';
    sprintf(logbuf, "Error: Duplicate variant ID '%s' in .%s file.\n", bufptr, extension);
  load_bim_ret_INVALID_FORMAT_2:
    logprintb();
  load_bim_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  load_bim_ret_ALL_MARKERS_EXCLUDED:
    retval = RET_ALL_MARKERS_EXCLUDED;
    break;
  }
 load_bim_ret_1:
  fclose_cond(bimfile);
  free_cond(sf_pos);
  free_cond(loadbuf2);
  return retval;
}

int32_t update_marker_cms(Two_col_params* update_cm, char* sorted_marker_ids, uintptr_t marker_ct, uintptr_t max_marker_id_len, uint32_t* marker_id_map, double* marker_cms) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* infile = NULL;
  char skipchar = update_cm->skipchar;
  uint32_t colid_first = (update_cm->colid < update_cm->colx);
  uintptr_t marker_ctl = (marker_ct + (BITCT - 1)) / BITCT;
  uintptr_t hit_ct = 0;
  uintptr_t miss_ct = 0;
  uintptr_t* already_seen;
  char* loadbuf;
  uintptr_t loadbuf_size;
  uint32_t colmin;
  uint32_t coldiff;
  char* colid_ptr;
  char* colx_ptr;
  int32_t sorted_idx;
  uint32_t marker_uidx;
  char cc;
  int32_t retval;
  if (wkspace_alloc_ul_checked(&already_seen, marker_ctl * sizeof(intptr_t))) {
    goto update_marker_cms_ret_NOMEM;
  }
  fill_ulong_zero(already_seen, marker_ctl);

  loadbuf = (char*)wkspace_base;
  loadbuf_size = wkspace_left;
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  }
  if (loadbuf_size <= MAXLINELEN) {
    goto update_marker_cms_ret_NOMEM;
  }
  retval = open_and_skip_first_lines(&infile, update_cm->fname, loadbuf, loadbuf_size, update_cm->skip);
  if (retval) {
    goto update_marker_cms_ret_1;
  }
  if (colid_first) {
    colmin = update_cm->colid - 1;
    coldiff = update_cm->colx - update_cm->colid;
  } else {
    colmin = update_cm->colx - 1;
    coldiff = update_cm->colid - update_cm->colx;
  }
  while (fgets(loadbuf, loadbuf_size, infile)) {
    if (!(loadbuf[loadbuf_size - 1])) {
      goto update_marker_cms_ret_NOMEM;
    }
    colid_ptr = skip_initial_spaces(loadbuf);
    cc = *colid_ptr;
    if (is_eoln_kns(cc) || (cc == skipchar)) {
      continue;
    }
    if (colid_first) {
      if (colmin) {
        colid_ptr = next_item_mult(colid_ptr, colmin);
      }
      colx_ptr = next_item_mult(colid_ptr, coldiff);
      if (no_more_items_kns(colx_ptr)) {
	goto update_marker_cms_ret_INVALID_FORMAT_2;
      }
    } else {
      colx_ptr = colid_ptr;
      if (colmin) {
	colx_ptr = next_item_mult(colx_ptr, colmin);
      }
      colid_ptr = next_item_mult(colx_ptr, coldiff);
      if (no_more_items_kns(colid_ptr)) {
	goto update_marker_cms_ret_INVALID_FORMAT_2;
      }
    }
    colid_ptr[strlen_se(colid_ptr)] = '\0';
    sorted_idx = bsearch_str(colid_ptr, sorted_marker_ids, max_marker_id_len, 0, marker_ct - 1);
    if (sorted_idx == -1) {
      miss_ct++;
      continue;
    }
    if (is_set(already_seen, sorted_idx)) {
      sprintf(logbuf, "Error: Duplicate variant %s in --update-cm file.\n", colid_ptr);
      logprintb();
      goto update_marker_cms_ret_INVALID_FORMAT;
    }
    set_bit(already_seen, sorted_idx);
    marker_uidx = marker_id_map[(uint32_t)sorted_idx];
    if (scan_double(colx_ptr, &(marker_cms[marker_uidx]))) {
      logprint("Error: Invalid centimorgan value in --update-cm file.\n");
      goto update_marker_cms_ret_INVALID_FORMAT;
    }
    hit_ct++;
  }
  if (!feof(infile)) {
    goto update_marker_cms_ret_READ_FAIL;
  }
  if (miss_ct) {
    sprintf(logbuf, "--update-cm: %" PRIuPTR " value%s changed, %" PRIuPTR " variant ID%s not present.\n", hit_ct, (hit_ct == 1)? "" : "s", miss_ct, (miss_ct == 1)? "" : "s");
  } else {
    sprintf(logbuf, "--update-cm: %" PRIuPTR " value%s changed.\n", hit_ct, (hit_ct == 1)? "" : "s");
  }
  logprintb();
  while (0) {
  update_marker_cms_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  update_marker_cms_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  update_marker_cms_ret_INVALID_FORMAT_2:
    logprint("Error: Fewer tokens than expected in --update-cm file line.\n");
  update_marker_cms_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 update_marker_cms_ret_1:
  fclose_cond(infile);
  wkspace_reset(wkspace_mark);
  return retval;
}

int32_t update_marker_pos(Two_col_params* update_map, char* sorted_marker_ids, uintptr_t marker_ct, uintptr_t max_marker_id_len, uint32_t* marker_id_map, uintptr_t* marker_exclude, uintptr_t* marker_exclude_ct_ptr, uint32_t* marker_pos, uint32_t* map_is_unsorted_ptr, Chrom_info* chrom_info_ptr) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* infile = NULL;
  char skipchar = update_map->skipchar;
  uint32_t colid_first = (update_map->colid < update_map->colx);
  uintptr_t marker_ctl = (marker_ct + (BITCT - 1)) / BITCT;
  uintptr_t hit_ct = 0;
  uintptr_t miss_ct = 0;
  uint32_t map_is_unsorted = ((*map_is_unsorted_ptr) & UNSORTED_CHROM);
  uintptr_t orig_exclude_ct = *marker_exclude_ct_ptr;
  uint32_t chrom_fo_idx_p1 = 0;
  uint32_t chrom_end = 0;
  uint32_t last_pos = 0;
  uintptr_t* already_seen;
  char* loadbuf;
  uintptr_t loadbuf_size;
  uint32_t colmin;
  uint32_t coldiff;
  char* colid_ptr;
  char* colx_ptr;
  int32_t sorted_idx;
  uint32_t marker_uidx;
  uint32_t marker_idx;
  char cc;
  int32_t retval;
  if (wkspace_alloc_ul_checked(&already_seen, marker_ctl * sizeof(intptr_t))) {
    goto update_marker_pos_ret_NOMEM;
  }
  fill_ulong_zero(already_seen, marker_ctl);

  loadbuf = (char*)wkspace_base;
  loadbuf_size = wkspace_left;
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  }
  if (loadbuf_size <= MAXLINELEN) {
    goto update_marker_pos_ret_NOMEM;
  }
  retval = open_and_skip_first_lines(&infile, update_map->fname, loadbuf, loadbuf_size, update_map->skip);
  if (retval) {
    goto update_marker_pos_ret_1;
  }
  if (colid_first) {
    colmin = update_map->colid - 1;
    coldiff = update_map->colx - update_map->colid;
  } else {
    colmin = update_map->colx - 1;
    coldiff = update_map->colid - update_map->colx;
  }
  while (fgets(loadbuf, loadbuf_size, infile)) {
    if (!(loadbuf[loadbuf_size - 1])) {
      goto update_marker_pos_ret_NOMEM;
    }
    colid_ptr = skip_initial_spaces(loadbuf);
    cc = *colid_ptr;
    if (is_eoln_kns(cc) || (cc == skipchar)) {
      continue;
    }
    if (colid_first) {
      if (colmin) {
        colid_ptr = next_item_mult(colid_ptr, colmin);
      }
      colx_ptr = next_item_mult(colid_ptr, coldiff);
      if (no_more_items_kns(colx_ptr)) {
	goto update_marker_pos_ret_INVALID_FORMAT_2;
      }
    } else {
      colx_ptr = colid_ptr;
      if (colmin) {
	colx_ptr = next_item_mult(colx_ptr, colmin);
      }
      colid_ptr = next_item_mult(colx_ptr, coldiff);
      if (no_more_items_kns(colid_ptr)) {
	goto update_marker_pos_ret_INVALID_FORMAT_2;
      }
    }
    colid_ptr[strlen_se(colid_ptr)] = '\0';
    sorted_idx = bsearch_str(colid_ptr, sorted_marker_ids, max_marker_id_len, 0, marker_ct - 1);
    if (sorted_idx == -1) {
      miss_ct++;
      continue;
    }
    if (is_set(already_seen, sorted_idx)) {
      sprintf(logbuf, "Error: Duplicate variant %s in --update-map file.\n", colid_ptr);
      logprintb();
      goto update_marker_pos_ret_INVALID_FORMAT;
    }
    set_bit(already_seen, sorted_idx);
    marker_uidx = marker_id_map[(uint32_t)sorted_idx];
    sorted_idx = atoi(colx_ptr);
    if (!(sorted_idx || ((colx_ptr[0] == '0') && (colx_ptr[1] == '\0')))) {
      logprint("Error: Invalid bp position value in --update-map file.\n");
      goto update_marker_pos_ret_INVALID_FORMAT;
    }
    if (sorted_idx < 0) {
      SET_BIT(marker_exclude, marker_uidx);
      *marker_exclude_ct_ptr += 1;
    } else {
      marker_pos[marker_uidx] = sorted_idx;
    }
    hit_ct++;
  }
  if (!feof(infile)) {
    goto update_marker_pos_ret_READ_FAIL;
  }
  if (miss_ct) {
    sprintf(logbuf, "--update-map: %" PRIuPTR " value%s updated, %" PRIuPTR " variant ID%s not present.\n", hit_ct, (hit_ct == 1)? "" : "s", miss_ct, (miss_ct == 1)? "" : "s");
  } else {
    sprintf(logbuf, "--update-map: %" PRIuPTR " value%s updated.\n", hit_ct, (hit_ct == 1)? "" : "s");
  }
  logprintb();
  marker_uidx = 0;
  marker_ct -= (*marker_exclude_ct_ptr) - orig_exclude_ct;
  if (!marker_ct) {
    logprint("Error: All variants excluded by --update-map (due to negative marker\npositions).\n");
    goto update_marker_pos_ret_ALL_MARKERS_EXCLUDED;
  }
  for (marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
    next_unset_unsafe_ck(marker_exclude, &marker_uidx);
    while (marker_uidx >= chrom_end) {
      chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[++chrom_fo_idx_p1];
      last_pos = 0;
    }
    if (last_pos > marker_pos[marker_uidx]) {
      map_is_unsorted |= UNSORTED_BP;
      if (!((*map_is_unsorted_ptr) & UNSORTED_BP)) {
	logprint("Warning: Base-pair positions are now unsorted!\n");
      }
      break;
    }
    last_pos = marker_pos[marker_uidx];
  }
  if (((*map_is_unsorted_ptr) & UNSORTED_BP) && (!(map_is_unsorted & UNSORTED_BP))) {
    logprint("Base-pair positions are now sorted.\n");
  }
  *map_is_unsorted_ptr = map_is_unsorted;
  while (0) {
  update_marker_pos_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  update_marker_pos_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  update_marker_pos_ret_INVALID_FORMAT_2:
    logprint("Error: Fewer tokens than expected in --update-map file line.\n");
  update_marker_pos_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  update_marker_pos_ret_ALL_MARKERS_EXCLUDED:
    retval = RET_ALL_MARKERS_EXCLUDED;
    break;
  }
 update_marker_pos_ret_1:
  fclose_cond(infile);
  wkspace_reset(wkspace_mark);
  return retval;
}

int32_t update_marker_names(Two_col_params* update_name, char* sorted_marker_ids, uintptr_t marker_ct, uintptr_t max_marker_id_len, uint32_t* marker_id_map, char* true_marker_ids) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* infile = NULL;
  char skipchar = update_name->skipchar;
  uint32_t colold_first = (update_name->colid < update_name->colx);
  uintptr_t marker_ctl = (marker_ct + (BITCT - 1)) / BITCT;
  uintptr_t hit_ct = 0;
  uintptr_t miss_ct = 0;
  uintptr_t* already_seen;
  char* loadbuf;
  uintptr_t loadbuf_size;
  uint32_t colmin;
  uint32_t coldiff;
  char* colold_ptr;
  char* colnew_ptr;
  int32_t sorted_idx;
  uintptr_t marker_uidx;
  uint32_t slen;
  char cc;
  int32_t retval;
  if (wkspace_alloc_ul_checked(&already_seen, marker_ctl * sizeof(intptr_t))) {
    goto update_marker_names_ret_NOMEM;
  }
  fill_ulong_zero(already_seen, marker_ctl);
  loadbuf = (char*)wkspace_base;
  loadbuf_size = wkspace_left;
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  }
  if (loadbuf_size <= MAXLINELEN) {
    goto update_marker_names_ret_NOMEM;
  }
  retval = open_and_skip_first_lines(&infile, update_name->fname, loadbuf, loadbuf_size, update_name->skip);
  if (retval) {
    goto update_marker_names_ret_1;
  }
  if (colold_first) {
    colmin = update_name->colid - 1;
    coldiff = update_name->colx - update_name->colid;
  } else {
    colmin = update_name->colx - 1;
    coldiff = update_name->colid - update_name->colx;
  }
  while (fgets(loadbuf, loadbuf_size, infile)) {
    if (!(loadbuf[loadbuf_size - 1])) {
      goto update_marker_names_ret_NOMEM;
    }
    colold_ptr = skip_initial_spaces(loadbuf);
    cc = *colold_ptr;
    if (is_eoln_kns(cc) || (cc == skipchar)) {
      continue;
    }
    if (colold_first) {
      if (colmin) {
        colold_ptr = next_item_mult(colold_ptr, colmin);
      }
      colnew_ptr = next_item_mult(colold_ptr, coldiff);
    } else {
      colnew_ptr = colold_ptr;
      if (colmin) {
	colnew_ptr = next_item_mult(colnew_ptr, colmin);
      }
      colold_ptr = next_item_mult(colnew_ptr, coldiff);
    }
    colold_ptr[strlen_se(colold_ptr)] = '\0';
    sorted_idx = bsearch_str(colold_ptr, sorted_marker_ids, max_marker_id_len, 0, marker_ct - 1);
    if (sorted_idx == -1) {
      miss_ct++;
      continue;
    }
    if (is_set(already_seen, sorted_idx)) {
      sprintf(logbuf, "Error: Duplicate variant %s in --update-name file.\n", colold_ptr);
      logprintb();
      goto update_marker_names_ret_INVALID_FORMAT;
    }
    set_bit(already_seen, sorted_idx);
    marker_uidx = marker_id_map[((uint32_t)sorted_idx)];
    slen = strlen_se(colnew_ptr);
    colnew_ptr[slen] = '\0';
    memcpy(&(true_marker_ids[marker_uidx * max_marker_id_len]), colnew_ptr, slen + 1);
    hit_ct++;
  }
  if (!feof(infile)) {
    goto update_marker_names_ret_READ_FAIL;
  }
  if (miss_ct) {
    sprintf(logbuf, "--update-name: %" PRIuPTR " value%s updated, %" PRIuPTR " variant ID%s not present.\n", hit_ct, (hit_ct == 1)? "" : "s", miss_ct, (miss_ct == 1)? "" : "s");
  } else {
    sprintf(logbuf, "--update-name: %" PRIuPTR " value%s updated.\n", hit_ct, (hit_ct == 1)? "" : "s");
  }
  while (0) {
  update_marker_names_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  update_marker_names_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  update_marker_names_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 update_marker_names_ret_1:
  fclose_cond(infile);
  wkspace_reset(wkspace_mark);
  return retval;
}

int32_t update_marker_alleles(char* update_alleles_fname, char* sorted_marker_ids, uintptr_t marker_ct, uintptr_t max_marker_id_len, uint32_t* marker_id_map, char** marker_allele_ptrs, uintptr_t* max_marker_allele_len_ptr, char* outname, char* outname_end) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* infile = NULL;
  FILE* errfile = NULL;
  int32_t retval = 0;
  uintptr_t marker_ctl = (marker_ct + (BITCT - 1)) / BITCT;
  uintptr_t max_marker_allele_len = *max_marker_allele_len_ptr;
  uintptr_t hit_ct = 0;
  uintptr_t miss_ct = 0;
  uintptr_t err_ct = 0;
  uintptr_t* already_seen;
  char* loadbuf;
  char* bufptr;
  char* bufptr2;
  char* bufptr3;
  uintptr_t loadbuf_size;
  uint32_t len2;
  uint32_t len;
  int32_t sorted_idx;
  uintptr_t marker_uidx;
  if (wkspace_alloc_ul_checked(&already_seen, marker_ctl * sizeof(intptr_t))) {
    goto update_marker_alleles_ret_NOMEM;
  }
  fill_ulong_zero(already_seen, marker_ctl);
  if (fopen_checked(&infile, update_alleles_fname, "r")) {
    goto update_marker_alleles_ret_OPEN_FAIL;
  }
  loadbuf_size = wkspace_left;
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  } else if (loadbuf_size <= MAXLINELEN) {
    goto update_marker_alleles_ret_NOMEM;
  }
  loadbuf = (char*)wkspace_alloc(loadbuf_size);
  loadbuf[loadbuf_size - 1] = ' ';
  while (fgets(loadbuf, loadbuf_size, infile)) {
    if (!loadbuf[loadbuf_size - 1]) {
      if (loadbuf_size == MAXLINEBUFLEN) {
	logprint("Error: Pathologically long line in --update-alleles file.\n");
	goto update_marker_alleles_ret_INVALID_FORMAT;
      } else {
	goto update_marker_alleles_ret_NOMEM;
      }
    }
    bufptr3 = skip_initial_spaces(loadbuf);
    if (is_eoln_kns(*bufptr3)) {
      continue;
    }
    bufptr2 = item_endnn(bufptr3);
    *bufptr2 = '\0';
    bufptr2 = skip_initial_spaces(&(bufptr2[1]));
    sorted_idx = bsearch_str(bufptr3, sorted_marker_ids, max_marker_id_len, 0, marker_ct - 1);
    if (sorted_idx == -1) {
      miss_ct++;
      continue;
    }
    if (is_set(already_seen, sorted_idx)) {
      sprintf(logbuf, "Error: Duplicate variant %s in --update-alleles file.\n", bufptr3);
      logprintb();
      goto update_marker_alleles_ret_INVALID_FORMAT;
    }
    set_bit(already_seen, sorted_idx);
    marker_uidx = marker_id_map[((uint32_t)sorted_idx)];
    bufptr = item_endnn(bufptr2);
    len2 = (uintptr_t)(bufptr - bufptr2);
    *bufptr = '\0';
    bufptr = skip_initial_spaces(&(bufptr[1]));
    len = strlen_se(bufptr);
    bufptr[len] = '\0';
    if ((!strcmp(bufptr2, marker_allele_ptrs[2 * marker_uidx])) && (!strcmp(bufptr, marker_allele_ptrs[2 * marker_uidx + 1]))) {
      bufptr2 = skip_initial_spaces(&(bufptr[len + 1]));
      bufptr = next_item(bufptr2);
      goto update_marker_alleles_match;
    } else if ((!strcmp(bufptr, marker_allele_ptrs[2 * marker_uidx])) && (!strcmp(bufptr2, marker_allele_ptrs[2 * marker_uidx + 1]))) {
      bufptr = skip_initial_spaces(&(bufptr[len + 1]));
      bufptr2 = next_item(bufptr);
    update_marker_alleles_match:
      len = strlen_se(bufptr);
      len2 = strlen_se(bufptr2);
      if (len >= max_marker_allele_len) {
	max_marker_allele_len = len + 1;
      }
      if (len2 >= max_marker_allele_len) {
	max_marker_allele_len = len2 + 1;
      }
      allele_reset(&(marker_allele_ptrs[2 * marker_uidx]), bufptr2, len2);
      allele_reset(&(marker_allele_ptrs[2 * marker_uidx + 1]), bufptr, len);
      hit_ct++;
    } else {
      if (!err_ct) {
	memcpy(outname_end, ".allele.no.snp", 15);
	if (fopen_checked(&errfile, outname, "w")) {
	  goto update_marker_alleles_ret_OPEN_FAIL;
	}
      }
      fputs(bufptr3, errfile);
      putc('\t', errfile);
      fputs(bufptr2, errfile);
      putc('\t', errfile);
      fputs(bufptr, errfile);
      if (putc_checked('\n', errfile)) {
	goto update_marker_alleles_ret_WRITE_FAIL;
      }
      err_ct++;
    }
  }
  if (!feof(infile)) {
    goto update_marker_alleles_ret_READ_FAIL;
  }
  *max_marker_allele_len_ptr = max_marker_allele_len;
  if (miss_ct) {
    sprintf(logbuf, "--update-alleles: %" PRIuPTR " variant%s updated, %" PRIuPTR " ID%s not present.\n", hit_ct, (hit_ct == 1)? "" : "s", miss_ct, (miss_ct == 1)? "" : "s");
  } else {
    sprintf(logbuf, "--update-alleles: %" PRIuPTR " variant%s updated.\n", hit_ct, (hit_ct == 1)? "" : "s");
  }
  logprintb();
  if (err_ct) {
    sprintf(logbuf, "%" PRIuPTR " update failure%s logged to %s.\n", err_ct, (err_ct == 1)? "" : "s", outname);
    logprintb();
  }

  while (0) {
  update_marker_alleles_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  update_marker_alleles_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  update_marker_alleles_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  update_marker_alleles_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  update_marker_alleles_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
  fclose_cond(infile);
  fclose_cond(errfile);
  wkspace_reset(wkspace_mark);
  return retval;
}

int32_t update_indiv_ids(char* update_ids_fname, char* sorted_person_ids, uintptr_t indiv_ct, uintptr_t max_person_id_len, uint32_t* indiv_id_map, char* person_ids) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* infile = NULL;
  int32_t retval = 0;
  uintptr_t indiv_ctl = (indiv_ct + (BITCT - 1)) / BITCT;
  uintptr_t hit_ct = 0;
  uintptr_t miss_ct = 0;
  char* idbuf;
  uintptr_t* already_seen;
  char* bufptr;
  char* bufptr2;
  char* wptr;
  uint32_t len;
  uint32_t len2;
  int32_t sorted_idx;
  uintptr_t indiv_uidx;
  if (wkspace_alloc_c_checked(&idbuf, max_person_id_len) ||
      wkspace_alloc_ul_checked(&already_seen, indiv_ctl * sizeof(intptr_t))) {
    goto update_indiv_ids_ret_NOMEM;
  }
  fill_ulong_zero(already_seen, indiv_ctl);
  if (fopen_checked(&infile, update_ids_fname, "r")) {
    goto update_indiv_ids_ret_OPEN_FAIL;
  }
  tbuf[MAXLINELEN - 1] = ' ';
  while (fgets(tbuf, MAXLINELEN, infile)) {
    if (!tbuf[MAXLINELEN - 1]) {
      logprint("Error: Pathologically long line in --update-ids file.\n");
      goto update_indiv_ids_ret_INVALID_FORMAT;
    }
    bufptr = skip_initial_spaces(tbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    len = strlen_se(bufptr);
    bufptr2 = skip_initial_spaces(&(bufptr[len + 1]));
    len2 = strlen_se(bufptr2);
    if (len + len2 + 2 > max_person_id_len) {
      miss_ct++;
      continue;
    }
    memcpy(idbuf, bufptr, len);
    idbuf[len] = '\t';
    memcpy(&(idbuf[len + 1]), bufptr2, len2);
    idbuf[len + len2 + 1] = '\0';
    sorted_idx = bsearch_str(idbuf, sorted_person_ids, max_person_id_len, 0, indiv_ct - 1);
    if (sorted_idx == -1) {
      miss_ct++;
      continue;
    }
    if (is_set(already_seen, sorted_idx)) {
      idbuf[len] = ' ';
      sprintf(logbuf, "Error: Duplicate individual %s in --update-ids file.\n", idbuf);
      logprintb();
      goto update_indiv_ids_ret_INVALID_FORMAT;
    }
    set_bit(already_seen, sorted_idx);
    indiv_uidx = indiv_id_map[((uint32_t)sorted_idx)];
    wptr = &(person_ids[indiv_uidx * max_person_id_len]);
    bufptr = skip_initial_spaces(&(bufptr2[len2 + 1]));
    bufptr2 = item_endnn(bufptr);
    wptr = memcpyax(wptr, bufptr, (uintptr_t)(bufptr2 - bufptr), '\t');
    bufptr = skip_initial_spaces(&(bufptr2[1]));
    memcpyx(wptr, bufptr, strlen_se(bufptr), '\0');
    hit_ct++;
  }
  if (!feof(infile)) {
    goto update_indiv_ids_ret_READ_FAIL;
  }
  if (miss_ct) {
    sprintf(logbuf, "--update-ids: %" PRIuPTR " %s updated, %" PRIuPTR " ID%s not present.\n", hit_ct, species_str(hit_ct), miss_ct, (miss_ct == 1)? "" : "s");
  } else {
    sprintf(logbuf, "--update-ids: %" PRIuPTR " %s updated.\n", hit_ct, species_str(hit_ct));
  }
  logprintb();

  while (0) {
  update_indiv_ids_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  update_indiv_ids_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  update_indiv_ids_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  update_indiv_ids_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
  fclose_cond(infile);
  wkspace_reset(wkspace_mark);
  return retval;
}

int32_t update_indiv_parents(char* update_parents_fname, char* sorted_person_ids, uintptr_t indiv_ct, uintptr_t max_person_id_len, uint32_t* indiv_id_map, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* founder_info) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* infile = NULL;
  int32_t retval = 0;
  uintptr_t indiv_ctl = (indiv_ct + (BITCT - 1)) / BITCT;
  uintptr_t hit_ct = 0;
  uintptr_t miss_ct = 0;
  char* idbuf;
  uintptr_t* already_seen;
  char* bufptr;
  char* bufptr2;
  char* bufptr3;
  char* wptr;
  uint32_t len;
  uint32_t len2;
  int32_t sorted_idx;
  uintptr_t indiv_uidx;
  if (wkspace_alloc_c_checked(&idbuf, max_person_id_len) ||
      wkspace_alloc_ul_checked(&already_seen, indiv_ctl * sizeof(intptr_t))) {
    goto update_indiv_parents_ret_NOMEM;
  }
  fill_ulong_zero(already_seen, indiv_ctl);
  if (fopen_checked(&infile, update_parents_fname, "r")) {
    goto update_indiv_parents_ret_OPEN_FAIL;
  }
  tbuf[MAXLINELEN - 1] = ' ';
  while (fgets(tbuf, MAXLINELEN, infile)) {
    if (!tbuf[MAXLINELEN - 1]) {
      logprint("Error: Pathologically long line in --update-parents file.\n");
      goto update_indiv_parents_ret_INVALID_FORMAT;
    }
    bufptr = skip_initial_spaces(tbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    len = strlen_se(bufptr);
    bufptr2 = skip_initial_spaces(&(bufptr[len + 1]));
    len2 = strlen_se(bufptr2);
    if (len + len2 + 2 > max_person_id_len) {
      miss_ct++;
      continue;
    }
    memcpy(idbuf, bufptr, len);
    idbuf[len] = '\t';
    memcpy(&(idbuf[len + 1]), bufptr2, len2);
    idbuf[len + len2 + 1] = '\0';
    sorted_idx = bsearch_str(idbuf, sorted_person_ids, max_person_id_len, 0, indiv_ct - 1);
    if (sorted_idx == -1) {
      miss_ct++;
      continue;
    }
    if (is_set(already_seen, sorted_idx)) {
      idbuf[len] = ' ';
      sprintf(logbuf, "Error: Duplicate individual %s in --update-parents file.\n", idbuf);
      logprintb();
      goto update_indiv_parents_ret_INVALID_FORMAT;
    }
    set_bit(already_seen, sorted_idx);
    indiv_uidx = indiv_id_map[((uint32_t)sorted_idx)];
    wptr = &(paternal_ids[indiv_uidx * max_paternal_id_len]);
    bufptr = skip_initial_spaces(&(bufptr2[len2 + 1]));
    bufptr2 = item_endnn(bufptr);
    len = (uintptr_t)(bufptr2 - bufptr);
    memcpyx(wptr, bufptr, len, '\0');
    wptr = &(maternal_ids[indiv_uidx * max_maternal_id_len]);
    bufptr3 = skip_initial_spaces(&(bufptr2[1]));
    len2 = strlen_se(bufptr3);
    memcpyx(wptr, bufptr3, len2, '\0');
    if ((len == 1) && (*bufptr == '0') && (len2 == 1) && (*bufptr3 == '0')) {
      SET_BIT(founder_info, indiv_uidx);
    } else {
      CLEAR_BIT(founder_info, indiv_uidx);
    }
    hit_ct++;
  }
  if (!feof(infile)) {
    goto update_indiv_parents_ret_READ_FAIL;
  }
  if (miss_ct) {
    sprintf(logbuf, "--update-parents: %" PRIuPTR " %s updated, %" PRIuPTR " ID%s not present.\n", hit_ct, species_str(hit_ct), miss_ct, (miss_ct == 1)? "" : "s");
  } else {
    sprintf(logbuf, "--update-parents: %" PRIuPTR " %s updated.\n", hit_ct, species_str(hit_ct));
  }
  logprintb();

  while (0) {
  update_indiv_parents_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  update_indiv_parents_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  update_indiv_parents_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  update_indiv_parents_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
  fclose_cond(infile);
  wkspace_reset(wkspace_mark);
  return retval;
}

int32_t update_indiv_sexes(char* update_sex_fname, char* sorted_person_ids, uintptr_t indiv_ct, uintptr_t max_person_id_len, uint32_t* indiv_id_map, uintptr_t* sex_nm, uintptr_t* sex_male) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* infile = NULL;
  int32_t retval = 0;
  uintptr_t indiv_ctl = (indiv_ct + (BITCT - 1)) / BITCT;
  uintptr_t hit_ct = 0;
  uintptr_t miss_ct = 0;
  char* idbuf;
  uintptr_t* already_seen;
  char* bufptr;
  char* bufptr2;
  uint32_t len;
  uint32_t len2;
  int32_t sorted_idx;
  uint32_t indiv_uidx;
  char cc;
  unsigned char ucc;
  if (wkspace_alloc_c_checked(&idbuf, max_person_id_len) ||
      wkspace_alloc_ul_checked(&already_seen, indiv_ctl * sizeof(intptr_t))) {
    goto update_indiv_sexes_ret_NOMEM;
  }
  fill_ulong_zero(already_seen, indiv_ctl);
  if (fopen_checked(&infile, update_sex_fname, "r")) {
    goto update_indiv_sexes_ret_OPEN_FAIL;
  }
  tbuf[MAXLINELEN - 1] = ' ';
  while (fgets(tbuf, MAXLINELEN, infile)) {
    if (!tbuf[MAXLINELEN - 1]) {
      logprint("Error: Pathologically long line in --update-sex file.\n");
      goto update_indiv_sexes_ret_INVALID_FORMAT;
    }
    bufptr = skip_initial_spaces(tbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    len = strlen_se(bufptr);
    bufptr2 = skip_initial_spaces(&(bufptr[len + 1]));
    len2 = strlen_se(bufptr2);
    if (len + len2 + 2 > max_person_id_len) {
      miss_ct++;
      continue;
    }
    memcpy(idbuf, bufptr, len);
    idbuf[len] = '\t';
    memcpy(&(idbuf[len + 1]), bufptr2, len2);
    idbuf[len + len2 + 1] = '\0';
    sorted_idx = bsearch_str(idbuf, sorted_person_ids, max_person_id_len, 0, indiv_ct - 1);
    if (sorted_idx == -1) {
      miss_ct++;
      continue;
    }
    if (is_set(already_seen, sorted_idx)) {
      idbuf[len] = ' ';
      sprintf(logbuf, "Error: Duplicate individual %s in --update-sex file.\n", idbuf);
      logprintb();
      goto update_indiv_sexes_ret_INVALID_FORMAT;
    }
    set_bit(already_seen, sorted_idx);
    indiv_uidx = indiv_id_map[((uint32_t)sorted_idx)];
    bufptr = skip_initial_spaces(&(bufptr2[len2 + 1]));
    if (no_more_items_kns(bufptr)) {
      logprint("Error: Fewer tokens than expected in --update-sex line.\n");
      goto update_indiv_sexes_ret_INVALID_FORMAT;
    }
    cc = *bufptr;
    ucc = ((unsigned char)cc) & 0xdfU;
    if ((cc < '0') || ((cc > '2') && (ucc != 'M') && (ucc != 'F')) || (bufptr[1] > ' ')) {
      logprint("Error: Invalid sex value in --update-sex line.  (Acceptable values: 1/M = male,\n2/F = female, 0 = missing.)\n");
      goto update_indiv_sexes_ret_INVALID_FORMAT;
    }
    if (cc == '0') {
      CLEAR_BIT(sex_nm, indiv_uidx);
      CLEAR_BIT(sex_male, indiv_uidx);
    } else {
      SET_BIT(sex_nm, indiv_uidx);
      if ((cc == '1') || (ucc == 'M')) {
	SET_BIT(sex_male, indiv_uidx);
      } else {
	CLEAR_BIT(sex_male, indiv_uidx);
      }
    }
    hit_ct++;
  }
  if (!feof(infile)) {
    goto update_indiv_sexes_ret_READ_FAIL;
  }
  if (miss_ct) {
    sprintf(logbuf, "--update-sex: %" PRIuPTR " %s updated, %" PRIuPTR " ID%s not present.\n", hit_ct, species_str(hit_ct), miss_ct, (miss_ct == 1)? "" : "s");
  } else {
    sprintf(logbuf, "--update-sex: %" PRIuPTR " %s updated.\n", hit_ct, species_str(hit_ct));
  }
  logprintb();

  while (0) {
  update_indiv_sexes_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  update_indiv_sexes_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  update_indiv_sexes_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  update_indiv_sexes_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
  fclose_cond(infile);
  wkspace_reset(wkspace_mark);
  return retval;
}

uint32_t flip_str(char** allele_str_ptr) {
  char* allele_str = *allele_str_ptr;
  if (allele_str[1] != '\0') {
    return 1;
  }
  if (allele_str == &(g_one_char_strs[130])) { // "A"
    *allele_str_ptr = (char*)(&(g_one_char_strs[168])); // "T"
  } else if (allele_str == &(g_one_char_strs[134])) {
    *allele_str_ptr = (char*)(&(g_one_char_strs[142]));
  } else if (allele_str == &(g_one_char_strs[142])) {
    *allele_str_ptr = (char*)(&(g_one_char_strs[134]));
  } else if (allele_str == &(g_one_char_strs[168])) {
    *allele_str_ptr = (char*)(&(g_one_char_strs[130]));
  } else if (allele_str == g_missing_geno_ptr) {
    return 0;
  } else {
    return 1;
  }
  return 2;
}

int32_t flip_strand(char* flip_fname, char* sorted_marker_ids, uintptr_t marker_ct, uintptr_t max_marker_id_len, uint32_t* marker_id_map, char** marker_allele_ptrs) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* flipfile = NULL;
  uint32_t non_acgt_ct = 0;
  int32_t retval = 0;
  uintptr_t marker_ctl = (marker_ct + (BITCT - 1)) / BITCT;
  uintptr_t hit_ct = 0;
  uintptr_t miss_ct = 0;
  uintptr_t* already_seen;
  char* bufptr;
  uint32_t slen;
  uint32_t cur_non_acgt0;
  int32_t sorted_idx;
  uintptr_t marker_uidx;
  if (wkspace_alloc_ul_checked(&already_seen, marker_ctl * sizeof(intptr_t))) {
    goto flip_strand_ret_NOMEM;
  }
  fill_ulong_zero(already_seen, marker_ctl);
  if (fopen_checked(&flipfile, flip_fname, "r")) {
    goto flip_strand_ret_OPEN_FAIL;
  }
  tbuf[MAXLINELEN - 1] = ' ';
  while (fgets(tbuf, MAXLINELEN, flipfile)) {
    if (!tbuf[MAXLINELEN - 1]) {
      logprint("Error: Pathologically long line in --flip file.\n");
      goto flip_strand_ret_INVALID_FORMAT;
    }
    bufptr = skip_initial_spaces(tbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    slen = strlen_se(bufptr);
    bufptr[slen] = '\0';
    sorted_idx = bsearch_str(bufptr, sorted_marker_ids, max_marker_id_len, 0, marker_ct - 1);
    if (sorted_idx == -1) {
      continue;
    }
    if (is_set(already_seen, sorted_idx)) {
      sprintf(logbuf, "Error: Duplicate SNP %s in --flip file.\n", bufptr);
      logprintb();
      goto flip_strand_ret_INVALID_FORMAT;
    }
    set_bit(already_seen, sorted_idx);
    marker_uidx = marker_id_map[((uint32_t)sorted_idx)];
    cur_non_acgt0 = 0;
    cur_non_acgt0 |= flip_str(&(marker_allele_ptrs[2 * marker_uidx]));
    cur_non_acgt0 |= flip_str(&(marker_allele_ptrs[2 * marker_uidx + 1]));
    non_acgt_ct += (cur_non_acgt0 & 1);
    hit_ct += (cur_non_acgt0 >> 1);
  }
  if (!feof(flipfile)) {
    goto flip_strand_ret_READ_FAIL;
  }
  if (miss_ct) {
    sprintf(logbuf, "--flip: %" PRIuPTR " SNP%s flipped, %" PRIuPTR " SNP ID%s not present.\n", hit_ct, (hit_ct == 1)? "" : "s", miss_ct, (miss_ct == 1)? "" : "s");
  } else {
    sprintf(logbuf, "--flip: %" PRIuPTR " SNP%s flipped.\n", hit_ct, (hit_ct == 1)? "" : "s");
  }
  logprintb();
  if (non_acgt_ct) {
    sprintf(logbuf, "Warning: %u variant%s had at least one non-A/C/G/T allele name.\n", non_acgt_ct, (non_acgt_ct == 1)? "" : "s");
    logprintb();
  }
  while (0) {
  flip_strand_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  flip_strand_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  flip_strand_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  flip_strand_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
  fclose_cond(flipfile);
  wkspace_reset(wkspace_mark);
  return retval;
}

const char keep_str[] = "keep";
const char keep_fam_str[] = "keep-fam";
const char remove_str[] = "remove";
const char remove_fam_str[] = "remove-fam";

const char* include_or_exclude_flag_str(uint32_t flags) {
  switch (flags) {
  case 2:
    return keep_str;
  case 3:
    return remove_str;
  case 6:
    return keep_fam_str;
  case 7:
    return remove_fam_str;
  }
  return NULL;
}

int32_t include_or_exclude(char* fname, char* sorted_ids, uintptr_t sorted_ids_ct, uintptr_t max_id_len, uint32_t* id_map, uintptr_t unfiltered_ct, uintptr_t* exclude_arr, uintptr_t* exclude_ct_ptr, uint32_t flags) {
  FILE* infile = NULL;
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t* exclude_arr_new = NULL;
  uintptr_t include_ct = 0;
  uintptr_t unfiltered_ctl = (unfiltered_ct + (BITCT - 1)) / BITCT;
  uint32_t do_exclude = flags & 1;
  // flags & 2 = indivs
  uint32_t families_only = flags & 4;
  int32_t retval = 0;
  char* id_buf;
  char* bufptr0;
  char* bufptr;
  uintptr_t ulii;
  uintptr_t uljj;
  int32_t ii;
  uint32_t unfiltered_idx;
  uint32_t cur_idx;
  uint32_t last_idx;

  if (!do_exclude) {
    if (wkspace_alloc_ul_checked(&exclude_arr_new, unfiltered_ctl * sizeof(intptr_t))) {
      goto include_or_exclude_ret_NOMEM;
    }
    fill_all_bits(exclude_arr_new, unfiltered_ct);
  }
  if (fopen_checked(&infile, fname, "r")) {
    goto include_or_exclude_ret_OPEN_FAIL;
  }
  tbuf[MAXLINELEN - 1] = ' ';
  if (flags & 2) {
    if (wkspace_alloc_c_checked(&id_buf, max_id_len)) {
      goto include_or_exclude_ret_NOMEM;
    }
    while (fgets(tbuf, MAXLINELEN, infile) != NULL) {
      if (!tbuf[MAXLINELEN - 1]) {
	sprintf(logbuf, "Error: Excessively long line in --%s file (max %d chars).\n", include_or_exclude_flag_str(flags), MAXLINELEN - 3);
	goto include_or_exclude_ret_INVALID_FORMAT;
      }
      bufptr0 = skip_initial_spaces(tbuf);
      if (is_eoln_kns(*bufptr0)) {
	continue;
      }
      if (!families_only) {
	bufptr = next_item(bufptr0);
	if (no_more_items_kns(bufptr)) {
	  sprintf(logbuf, "Error: Fewer tokens than expected in --%s line.\n", include_or_exclude_flag_str(flags));
	  goto include_or_exclude_ret_INVALID_FORMAT;
	}
	ii = bsearch_fam_indiv(id_buf, sorted_ids, max_id_len, sorted_ids_ct, bufptr0, bufptr);
	if (ii != -1) {
	  unfiltered_idx = id_map[(uint32_t)ii];
	  if (do_exclude) {
	    if (IS_SET(exclude_arr, unfiltered_idx)) {
	      goto include_or_exclude_ret_INVALID_FORMAT_2;
	    }
	    SET_BIT(exclude_arr, unfiltered_idx);
	    *exclude_ct_ptr += 1;
	  } else if (!IS_SET(exclude_arr, unfiltered_idx)) {
	    if (!IS_SET(exclude_arr_new, unfiltered_idx)) {
	      goto include_or_exclude_ret_INVALID_FORMAT_2;
	    }
	    CLEAR_BIT(exclude_arr_new, unfiltered_idx);
	    include_ct++;
	  }
	}
      } else {
	bsearch_fam(id_buf, sorted_ids, max_id_len, sorted_ids_ct, bufptr0, &cur_idx, &last_idx);
	while (cur_idx < last_idx) {
	  unfiltered_idx = id_map[cur_idx++];
	  if (do_exclude) {
	    if (IS_SET(exclude_arr, unfiltered_idx)) {
	      goto include_or_exclude_ret_INVALID_FORMAT_2;
	    }
	    SET_BIT(exclude_arr, unfiltered_idx);
	    *exclude_ct_ptr += 1;
	  } else {
	    if (!IS_SET(exclude_arr_new, unfiltered_idx)) {
	      goto include_or_exclude_ret_INVALID_FORMAT_2;
	    }
	    CLEAR_BIT(exclude_arr_new, unfiltered_idx);
	    include_ct++;
	  }
	}
      }
    }
  } else {
    while (fgets(tbuf, MAXLINELEN, infile) != NULL) {
      if (!tbuf[MAXLINELEN - 1]) {
	sprintf(logbuf, "Error: Excessively long line in --%s file (max %d chars).\n", flags? "exclude" : "extract", MAXLINELEN - 3);
	goto include_or_exclude_ret_INVALID_FORMAT;
      }
      bufptr0 = skip_initial_spaces(tbuf);
      if (is_eoln_kns(*bufptr0)) {
	continue;
      }
      bufptr = item_endnn(bufptr0);
      *bufptr = '\0';
      ulii = bsearch_str_lb(sorted_ids, sorted_ids_ct, max_id_len, bufptr0);
      if (ulii != sorted_ids_ct) {
        *bufptr = ' ';
        bufptr[1] = '\0';
        uljj = bsearch_str_lb(sorted_ids, sorted_ids_ct, max_id_len, bufptr0);
        if (uljj > ulii) {
          do {
	    unfiltered_idx = id_map[ulii];
	    if (do_exclude) {
	      if (!IS_SET(exclude_arr, unfiltered_idx)) {
	        SET_BIT(exclude_arr, unfiltered_idx);
		*exclude_ct_ptr += 1;
	      }
	    } else if (!IS_SET(exclude_arr, unfiltered_idx)) {
	      if (IS_SET(exclude_arr_new, unfiltered_idx)) {
		CLEAR_BIT(exclude_arr_new, unfiltered_idx);
                include_ct++;
	      }
	    }
	  } while (++ulii < uljj);
	}
      }
    }
  }
  if (!feof(infile)) {
    goto include_or_exclude_ret_READ_FAIL;
  }
  if (!do_exclude) {
    memcpy(exclude_arr, exclude_arr_new, unfiltered_ctl * sizeof(intptr_t));
    *exclude_ct_ptr = unfiltered_ct - include_ct;
  }
  if (*exclude_ct_ptr == unfiltered_ct) {
    sprintf(logbuf, "Error: No %s remaining after --%s.\n", (flags & 2)? g_species_plural : "variants", include_or_exclude_flag_str(flags));
    logprintb();
    goto include_or_exclude_ret_ALL_SAMPLES_EXCLUDED;
  }
  while (0) {
  include_or_exclude_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  include_or_exclude_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  include_or_exclude_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  include_or_exclude_ret_INVALID_FORMAT_2:
    sprintf(logbuf, "Error: Duplicate ID in --%s file.\n", include_or_exclude_flag_str(flags));
  include_or_exclude_ret_INVALID_FORMAT:
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  include_or_exclude_ret_ALL_SAMPLES_EXCLUDED:
    retval = RET_ALL_SAMPLES_EXCLUDED;
    break;
  }
  wkspace_reset(wkspace_mark);
  fclose_cond(infile);
  return retval;
}

int32_t filter_attrib(char* fname, char* condition_str, char* sorted_ids, uintptr_t sorted_ids_ct, uintptr_t max_id_len, uint32_t* id_map, uintptr_t unfiltered_ct, uintptr_t* exclude_arr, uintptr_t* exclude_ct_ptr, uint32_t is_indiv) {
  FILE* infile = NULL;
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t include_ct = 0;
  uintptr_t unfiltered_ctl = (unfiltered_ct + (BITCT - 1)) / BITCT;
  uintptr_t* cur_neg_matches = NULL;
  char* sorted_pos_match = NULL;
  char* sorted_neg_match = NULL;
  char* id_buf = NULL;
  uint32_t pos_match_ct = 0;
  uint32_t neg_match_ct = 0;
  uint32_t neg_match_ctl = 0;
  uintptr_t max_pos_match_len = 0;
  uintptr_t max_neg_match_len = 0;
  uint32_t is_neg = 0;
  int32_t retval = 0;
  uintptr_t* exclude_arr_new;
  uintptr_t* already_seen;
  char* loadbuf;
  char* cond_ptr;
  char* bufptr;
  char* bufptr2;
  uintptr_t loadbuf_size;
  uintptr_t pos_match_idx;
  uintptr_t neg_match_idx;
  uintptr_t ulii;
  uint32_t unfiltered_idx;
  uint32_t pos_match_needed;
  uint32_t cur_neg_match_ct;
  int32_t sorted_idx;
  
  if (wkspace_alloc_ul_checked(&exclude_arr_new, unfiltered_ctl * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(&already_seen, unfiltered_ctl * sizeof(intptr_t)) ||
      (is_indiv && wkspace_alloc_c_checked(&id_buf, max_id_len))) { 
    goto filter_attrib_ret_NOMEM;
  }
  fill_all_bits(exclude_arr_new, unfiltered_ct);
  fill_ulong_zero(already_seen, unfiltered_ctl);
  if (condition_str) {
    // allow NULL condition_str; this means all indivs/variants named in the
    // file are included
    cond_ptr = condition_str;
    while (1) {
      while (*cond_ptr == ',') {
	cond_ptr++;
      }
      if (*cond_ptr == '-') {
	cond_ptr++;
	if (*cond_ptr == ',') {
	  continue;
	} else if (*cond_ptr == '-') {
	  sprintf(logbuf, "Error: --filter-attrib%s condition cannot contain consecutive dashes.\n", is_indiv? "-indiv" : "");
	  goto filter_attrib_ret_INVALID_CMDLINE_2;
	}
	is_neg = 1;
      } else if (!(*cond_ptr)) {
        // command-line parser actually comma-terminates the string for us, so
        // no need to check for null elsewhere
	break;
      }
      bufptr = strchr(cond_ptr, ',');
      ulii = (uintptr_t)(bufptr - cond_ptr);
      if (is_neg) {
	neg_match_ct++;
	if (ulii >= max_neg_match_len) {
	  max_neg_match_len = ulii + 1;
	}
      } else {
        pos_match_ct++;
	if (ulii >= max_pos_match_len) {
          max_pos_match_len = ulii + 1;
	}
      }
      cond_ptr = bufptr;
      is_neg = 0;
    }
    if (pos_match_ct) {
      if (wkspace_alloc_c_checked(&sorted_pos_match, max_pos_match_len * pos_match_ct)) {
	goto filter_attrib_ret_NOMEM;
      }
    }
    if (neg_match_ct) {
      neg_match_ctl = (neg_match_ct + (BITCT - 1)) / BITCT;
      if (wkspace_alloc_c_checked(&sorted_neg_match, max_neg_match_len * neg_match_ct) ||
	  wkspace_alloc_ul_checked(&cur_neg_matches, neg_match_ctl * sizeof(intptr_t))) {
        goto filter_attrib_ret_NOMEM;
      }
    }
    pos_match_idx = 0;
    neg_match_idx = 0;
    cond_ptr = condition_str;
    is_neg = 0;
    while (1) {
      while (*cond_ptr == ',') {
	cond_ptr++;
      }
      if (*cond_ptr == '-') {
	cond_ptr++;
	if (*cond_ptr == ',') {
	  continue;
	}
	is_neg = 1;
      } else if (!(*cond_ptr)) {
	break;
      }
      bufptr = strchr(cond_ptr, ',');
      ulii = (uintptr_t)(bufptr - cond_ptr);
      if (is_neg) {
	memcpyx(&(sorted_neg_match[neg_match_idx * max_neg_match_len]), cond_ptr, ulii, '\0');
	neg_match_idx++;
      } else {
	memcpyx(&(sorted_pos_match[pos_match_idx * max_pos_match_len]), cond_ptr, ulii, '\0');
        pos_match_idx++;
      }
      cond_ptr = bufptr;
      is_neg = 0;
    }
    if (pos_match_ct) {
      qsort(sorted_pos_match, pos_match_ct, max_pos_match_len, strcmp_casted);
      bufptr = scan_for_duplicate_ids(sorted_pos_match, pos_match_ct, max_pos_match_len);
      if (bufptr) {
	sprintf(logbuf, "Error: Duplicate attribute '%s' in --filter-attrib%s argument.\n", bufptr, is_indiv? "-indiv" : "");
	goto filter_attrib_ret_INVALID_CMDLINE_2;
      }
    }
    if (neg_match_ct) {
      qsort(sorted_neg_match, neg_match_ct, max_neg_match_len, strcmp_casted);
      bufptr = scan_for_duplicate_ids(sorted_neg_match, neg_match_ct, max_neg_match_len);
      if (bufptr) {
	sprintf(logbuf, "Error: Duplicate attribute '%s' in --filter-attrib%s argument.\n", bufptr, is_indiv? "-indiv" : "");
	goto filter_attrib_ret_INVALID_CMDLINE_2;
      }
      // actually may make sense to have same attribute as a positive and
      // negative condition, so we don't check for that
    }
  }
  loadbuf_size = wkspace_left;
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  } else if (loadbuf_size <= MAXLINELEN) {
    goto filter_attrib_ret_NOMEM;
  }
  if (fopen_checked(&infile, fname, "r")) {
    goto filter_attrib_ret_OPEN_FAIL;
  }
  loadbuf = (char*)wkspace_base;
  loadbuf[loadbuf_size - 1] = ' ';
  while (fgets(loadbuf, loadbuf_size, infile)) {
    if (!loadbuf[loadbuf_size - 1]) {
      if (loadbuf_size == MAXLINEBUFLEN) {
	sprintf(logbuf, "Error: Pathologically long line in --filter-attrib%s file.\n", is_indiv? "-indiv" : "");
        goto filter_attrib_ret_INVALID_FORMAT_2;
      }
      goto filter_attrib_ret_NOMEM;
    }
    bufptr = skip_initial_spaces(loadbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    if (is_indiv) {
      bufptr2 = next_item(bufptr);
      if (!bufptr2) {
        logprint("Error: --filter-attrib-indiv line has only one token.\n");
        goto filter_attrib_ret_INVALID_FORMAT;
      }
      sorted_idx = bsearch_fam_indiv(id_buf, sorted_ids, max_id_len, sorted_ids_ct, bufptr, bufptr2);
      cond_ptr = skip_initial_spaces(item_endnn(bufptr2));
    } else {
      bufptr2 = item_endnn(bufptr);
      cond_ptr = skip_initial_spaces(bufptr2);
      *bufptr2 = '\0';
      sorted_idx = bsearch_str(bufptr, sorted_ids, max_id_len, 0, sorted_ids_ct - 1);
    }
    if (sorted_idx == -1) {
      continue;
    }
    if (is_set(already_seen, sorted_idx)) {
      if (is_indiv) {
        bufptr = (char*)memchr(id_buf, '\t', max_id_len);
	*bufptr = ' ';
        sprintf(logbuf, "Error: Duplicate individual %s in --filter-attrib-indiv file.\n", id_buf);
      } else {
	sprintf(logbuf, "Error: Duplicate variant %s in --filter-attrib file.\n", bufptr);
      }
      goto filter_attrib_ret_INVALID_FORMAT_2;
    }
    set_bit(already_seen, sorted_idx);
    unfiltered_idx = id_map[(uint32_t)sorted_idx];
    pos_match_needed = pos_match_ct;
    cur_neg_match_ct = 0;
    if (neg_match_ct) {
      fill_ulong_zero(cur_neg_matches, neg_match_ctl);
    }
    while (!is_eoln_kns(*cond_ptr)) {
      bufptr2 = cond_ptr;
      bufptr = item_endnn(cond_ptr);
      ulii = (uintptr_t)(bufptr - bufptr2);
      cond_ptr = skip_initial_spaces(bufptr);
      *bufptr = '\0';
      if (pos_match_needed) {
	if (max_pos_match_len > ulii) {
	  sorted_idx = bsearch_str(bufptr2, sorted_pos_match, max_pos_match_len, 0, pos_match_ct - 1);
          if (sorted_idx != -1) {
            pos_match_needed = 0;
	  }
	}
      }
      if (cur_neg_match_ct < neg_match_ct) {
        if (max_neg_match_len > ulii) {
          sorted_idx = bsearch_str(bufptr2, sorted_neg_match, max_neg_match_len, 0, neg_match_ct - 1);
	  if ((sorted_idx != -1) && (!is_set(cur_neg_matches, sorted_idx))) {
            cur_neg_match_ct++;
	    if (cur_neg_match_ct == neg_match_ct) {
	      // fail
	      pos_match_needed = 1;
              break;
	    }
            set_bit(cur_neg_matches, sorted_idx);
	  }
	}
      }
    }
    if (pos_match_needed) {
      continue;
    }
    // full negative match causes pos_match_needed to be set, so no further
    // check required
    clear_bit(exclude_arr_new, unfiltered_idx);
    include_ct++;
  }
  if (!feof(infile)) {
    goto filter_attrib_ret_READ_FAIL;
  }
  if (!include_ct) {
    sprintf(logbuf, "Error: No %s remaining after --filter-attrib%s.\n", is_indiv? g_species_plural : "variants", is_indiv? "-indiv" : "");
    if (is_indiv) {
      retval = RET_ALL_SAMPLES_EXCLUDED;
    } else {
      retval = RET_ALL_MARKERS_EXCLUDED;
    }
    goto filter_attrib_ret_1;
  }
  sprintf(logbuf, "--filter-attrib%s: %" PRIuPTR " %s remaining.\n", is_indiv? "-indiv" : "", include_ct, is_indiv? species_str(include_ct) : "variants");
  logprintb();
  memcpy(exclude_arr, exclude_arr_new, unfiltered_ctl * sizeof(intptr_t));
  *exclude_ct_ptr = unfiltered_ct - include_ct;
  while (0) {
  filter_attrib_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  filter_attrib_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  filter_attrib_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  filter_attrib_ret_INVALID_CMDLINE_2:
    logprintb();
    retval = RET_INVALID_CMDLINE;
    break;
  filter_attrib_ret_INVALID_FORMAT_2:
    logprintb();
  filter_attrib_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 filter_attrib_ret_1:
  wkspace_reset(wkspace_mark);
  fclose_cond(infile);
  return retval;
}

int32_t read_dists(char* dist_fname, char* id_fname, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_ct, char* person_ids, uintptr_t max_person_id_len, uintptr_t cluster_ct, uint32_t* cluster_starts, uint32_t* indiv_to_cluster, uint32_t for_cluster_flag, uint32_t is_max_dist, double* dists, uint32_t neighbor_n2, double* neighbor_quantiles, uint32_t* neighbor_qindices) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* dist_file = NULL;
  FILE* id_file = NULL;
  uintptr_t id_entry_ct = indiv_ct;
  uintptr_t matching_entry_ct = indiv_ct;
  char* id_buf = &(tbuf[MAXLINELEN]);
  uint64_t* fidx_to_memidx = NULL; // high 32 bits = fidx, low 32 = memidx
  uint32_t is_presorted = cluster_ct? 0 : 1;
  int32_t retval = 0;
  char* sorted_ids;
  uint32_t* id_map;
  char* fam_id;
  char* indiv_id;
  double* dptr;
  uint64_t fpos;
  uint64_t fpos2;
  uintptr_t memidx1;
  uintptr_t memidx2;
  uintptr_t fidx1;
  uint64_t fidx2;
  uintptr_t trimem;
  uint64_t trif;
  uintptr_t clidx1;
  uintptr_t clidx2;
  uint64_t ullii;
  uintptr_t ulii;
  uintptr_t uljj;
  double cur_ibs;
  uint32_t uii;
  int32_t ii;
  if (fopen_checked(&dist_file, dist_fname, "rb")) {
    goto read_dists_ret_OPEN_FAIL;
  }
  if (fseeko(dist_file, 0, SEEK_END)) {
    goto read_dists_ret_READ_FAIL;
  }
  if (id_fname) {
    if (wkspace_alloc_ull_checked(&fidx_to_memidx, indiv_ct * sizeof(int64_t))) {
      goto read_dists_ret_NOMEM;
    }
    for (ulii = 0; ulii < indiv_ct; ulii++) {
      fidx_to_memidx[ulii] = 0xffffffffffffffffLLU;
    }
    if (fopen_checked(&id_file, id_fname, "r")) {
      goto read_dists_ret_OPEN_FAIL;
    }
    retval = sort_item_ids(&sorted_ids, &id_map, unfiltered_indiv_ct, indiv_exclude, unfiltered_indiv_ct - indiv_ct, person_ids, max_person_id_len, 0, 1, strcmp_deref);
    if (retval) {
      goto read_dists_ret_1;
    }
    id_entry_ct = 0;
    matching_entry_ct = 0;
    tbuf[MAXLINELEN - 1] = ' ';
    while (fgets(tbuf, MAXLINELEN, id_file)) {
      if (!tbuf[MAXLINELEN - 1]) {
        sprintf(logbuf, "Error: Pathologically long line in %s.\n", id_fname);
        logprintb();
        goto read_dists_ret_INVALID_FORMAT;
      }
      fam_id = skip_initial_spaces(tbuf);
      if (is_eoln_kns(*fam_id)) {
        continue;
      }
      indiv_id = next_item(fam_id);
      if (no_more_items(indiv_id)) {
        logprint("Error: Improperly formatted --read-dists ID file.\n");
        goto read_dists_ret_INVALID_FORMAT;
      }
      ii = bsearch_fam_indiv(id_buf, sorted_ids, max_person_id_len, indiv_ct, fam_id, indiv_id);
      if (ii == -1) {
        is_presorted = 0;
        id_entry_ct++;
        continue;
      }
      uii = id_map[(uint32_t)ii];
      if (uii != id_entry_ct) {
        is_presorted = 0;
      }
      if (fidx_to_memidx[uii] != 0xffffffffffffffffLLU) {
        logprint("Error: Duplicate ID in --read-dists ID file.\n");
        goto read_dists_ret_INVALID_FORMAT;
      }
      if (cluster_ct && (!neighbor_n2)) {
	// if cluster_ct && neighbor_n2, best to
	// postpone indiv_to_cluster dereference
        fidx_to_memidx[uii] = ((uint64_t)(indiv_to_cluster[uii])) | (((uint64_t)id_entry_ct) << 32);
      } else {
        fidx_to_memidx[uii] = ((uint64_t)uii) | (((uint64_t)id_entry_ct) << 32);
      }
      id_entry_ct++;
      matching_entry_ct++;
    }
    if (!feof(id_file)) {
      goto read_dists_ret_READ_FAIL;
    }
    fclose_null(&id_file);
    if (matching_entry_ct < indiv_ct) {
      logprint("Error: --read-dists ID file does not contain all individuals in current run.\n");
      goto read_dists_ret_INVALID_FORMAT;
    }
  } else if (cluster_ct) {
    if (wkspace_alloc_ull_checked(&fidx_to_memidx, indiv_ct * sizeof(int64_t))) {
      goto read_dists_ret_NOMEM;
    }
    if (neighbor_n2) {
      for (id_entry_ct = 0; id_entry_ct < indiv_ct; id_entry_ct++) {
	fidx_to_memidx[id_entry_ct] = ((uint64_t)id_entry_ct) | (((uint64_t)id_entry_ct) << 32);
      }
    } else {
      for (id_entry_ct = 0; id_entry_ct < indiv_ct; id_entry_ct++) {
	fidx_to_memidx[id_entry_ct] = ((uint64_t)(indiv_to_cluster[id_entry_ct])) | (((uint64_t)id_entry_ct) << 32);
      }
    }
  }
  if (!is_presorted) {
#ifdef __cplusplus
    // std::sort is faster than qsort for basic types.  See e.g. Anders
    // Kaseorg's answer to
    // http://www.quora.com/Software-Engineering/Generally-how-much-faster-is-C-compared-to-C++
    std::sort((int64_t*)fidx_to_memidx, (int64_t*)(&(fidx_to_memidx[indiv_ct])));
#else
    qsort(fidx_to_memidx, indiv_ct, sizeof(int64_t), llcmp);
#endif
  }
  fpos = (((uint64_t)id_entry_ct) * (id_entry_ct - 1)) * (sizeof(double) / 2);
  if (ftello(dist_file) != (int64_t)fpos) {
    sprintf(logbuf, "Error: Invalid --read-dists filesize (%" PRIu64 " bytes expected).\n", fpos);
    logprintb();
    goto read_dists_ret_INVALID_FORMAT;
  }
  rewind(dist_file);
  if (is_presorted) {
    if (fread(dists, 1, fpos, dist_file) < fpos) {
      goto read_dists_ret_READ_FAIL;
    }
    if (neighbor_n2) {
      dptr = dists;
      for (memidx1 = 1; memidx1 < indiv_ct; memidx1++) {
	for (memidx2 = 0; memidx2 < memidx1; memidx2++) {
          update_neighbor(indiv_ct, neighbor_n2, memidx1, memidx2, *dptr++, neighbor_quantiles, neighbor_qindices);
	}
      }
    }
  } else {
    fpos = 0;
    if ((!cluster_ct) || (!neighbor_n2)) {
      for (ulii = 1; ulii < indiv_ct; ulii++) {
	ullii = fidx_to_memidx[ulii];
	memidx1 = (uintptr_t)(ullii & (ONELU * 0xffffffffU));
	fidx1 = (uintptr_t)(ullii >> 32);
	trimem = (memidx1 * (memidx1 - 1)) / 2;
	trif = (((uint64_t)fidx1) * (fidx1 - 1)) * (sizeof(double) / 2);
	for (uljj = 0; uljj < ulii; uljj++) {
	  ullii = fidx_to_memidx[uljj];
	  memidx2 = (uintptr_t)(ullii & (ONELU * 0xffffffffU));
	  fidx2 = (uint64_t)(ullii >> 32);
	  fpos2 = trif + (fidx2 * sizeof(double));
	  if (fpos2 > fpos) {
	    fpos = fpos2;
	    if (fseeko(dist_file, fpos, SEEK_SET)) {
	      goto read_dists_ret_READ_FAIL;
	    }
	  }
	  if (fread(&cur_ibs, 1, sizeof(double), dist_file) < sizeof(double)) {
	    goto read_dists_ret_READ_FAIL;
	  }
	  if (memidx2 < memidx1) {
	    dists[trimem + memidx2] += cur_ibs;
	  } else if (memidx2 > memidx1) {
	    dists[((memidx2 * (memidx2 - 1)) / 2) + memidx1] += cur_ibs;
	  }
	  if (neighbor_n2) {
	    update_neighbor(indiv_ct, neighbor_n2, memidx1, memidx2, cur_ibs, neighbor_quantiles, neighbor_qindices);
	  }
	  fpos += sizeof(double);
	}
      }
    } else {
      for (ulii = 1; ulii < indiv_ct; ulii++) {
	ullii = fidx_to_memidx[ulii];
	memidx1 = (uintptr_t)(ullii & (ONELU * 0xffffffffU));
	fidx1 = (uintptr_t)(ullii >> 32);
        clidx1 = indiv_to_cluster[memidx1];
	trimem = (clidx1 * (clidx1 - 1)) / 2;
	trif = (((uint64_t)fidx1) * (fidx1 - 1)) * (sizeof(double) / 2);
	for (uljj = 0; uljj < ulii; uljj++) {
	  ullii = fidx_to_memidx[uljj];
	  memidx2 = (uintptr_t)(ullii & (ONELU * 0xffffffffU));
	  fidx2 = (uint64_t)(ullii >> 32);
	  clidx2 = indiv_to_cluster[memidx2];
	  fpos2 = trif + (fidx2 * sizeof(double));
	  if (fpos2 > fpos) {
	    fpos = fpos2;
	    if (fseeko(dist_file, fpos, SEEK_SET)) {
	      goto read_dists_ret_READ_FAIL;
	    }
	  }
	  if (fread(&cur_ibs, 1, sizeof(double), dist_file) < sizeof(double)) {
	    goto read_dists_ret_READ_FAIL;
	  }
	  if (!is_max_dist) {
	    if (clidx2 < clidx1) {
	      dists[trimem + clidx2] += cur_ibs;
	    } else if (clidx2 > clidx1) {
	      dists[((clidx2 * (clidx2 - 1)) / 2) + clidx1] += cur_ibs;
	    }
	  } else {
	    if (clidx1 != clidx2) {
	      if (clidx2 < clidx1) {
	        dptr = &(dists[trimem + clidx2]);
	      } else {
		dptr = &(dists[((clidx2 * (clidx2 - 1)) / 2) + clidx1]);
	      }
	      if (cur_ibs < (*dptr)) {
		*dptr = cur_ibs;
	      }
	    }
	  }
	  if (neighbor_n2) {
	    update_neighbor(indiv_ct, neighbor_n2, memidx1, memidx2, cur_ibs, neighbor_quantiles, neighbor_qindices);
	  }
	  fpos += sizeof(double);
	}
      }
    }
  }
  if (for_cluster_flag != 2) {
    sprintf(logbuf, "--read-dists: %" PRIuPTR " values loaded%s.\n", (indiv_ct * (indiv_ct - 1)) / 2, for_cluster_flag? " for --cluster/--neighbor" : "");
    logprintb();
  }
  while (0) {
  read_dists_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  read_dists_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  read_dists_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  read_dists_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 read_dists_ret_1:
  wkspace_reset(wkspace_mark);
  fclose_cond(dist_file);
  fclose_cond(id_file);
  return retval;
}

int32_t load_covars(char* covar_fname, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_ct, char* person_ids, uintptr_t max_person_id_len, double missing_phenod, uint32_t covar_modifier, Range_list* covar_range_list_ptr, uint32_t gxe_mcovar, uintptr_t* covar_ct_ptr, char** covar_names_ptr, uintptr_t* max_covar_name_len_ptr, uintptr_t* pheno_nm, uint32_t* pheno_nm_ct_ptr, uintptr_t** covar_nm_ptr, double** covar_d_ptr, uintptr_t** gxe_covar_nm_ptr, uintptr_t** gxe_covar_c_ptr) {
  // similar to load_clusters() in plink_cluster.c
  unsigned char* wkspace_mark = wkspace_base;
  unsigned char* wkspace_mark2 = NULL;
  FILE* covar_file = NULL;
  uintptr_t indiv_ctl = (indiv_ct + (BITCT - 1)) / BITCT;
  uintptr_t topsize = 0;
  uintptr_t covar_raw_ct = 0;
  uintptr_t loaded_indiv_ct = 0;
  uintptr_t missing_cov_ct = 0;
  uint32_t* indiv_idx_to_uidx = NULL;
  char* sorted_covar_name_flag_ids = NULL;
  uint32_t* covar_name_flag_id_map = NULL;
  int32_t* covar_name_flag_seen_idxs = NULL;
  char* covar_names = NULL;
  uintptr_t* covar_nm = NULL;
  double* covar_d = NULL;
  uintptr_t* gxe_covar_nm = NULL;
  uintptr_t* gxe_covar_c = NULL;
  double* dptr = NULL;
  char* bufptr = NULL;
  uintptr_t max_covar_name_len = 1;
  double dxx = 0.0;
  uint32_t keep_pheno_on_missing_cov = covar_modifier & COVAR_KEEP_PHENO_ON_MISSING_COV;
  int32_t retval = 0;
  uintptr_t covar_raw_ctl;
  uintptr_t covar_ct;
  uintptr_t* covars_active;
  uintptr_t* already_seen;
  char* sorted_ids;
  uint32_t* id_map;
  char* loadbuf;
  uintptr_t loadbuf_size;
  char* bufptr2;
  uintptr_t covar_uidx;
  uintptr_t covar_idx;
  uint32_t header_absent;
  uint32_t min_covar_col_ct;
  uint32_t uii;
  uint32_t indiv_idx;
  uint32_t indiv_uidx;
  uint32_t covar_missing;
  int32_t ii;

  if ((!keep_pheno_on_missing_cov) || gxe_mcovar) {
    indiv_idx_to_uidx = (uint32_t*)top_alloc(&topsize, indiv_ct * sizeof(int32_t));
    if (!indiv_idx_to_uidx) {
      goto load_covars_ret_NOMEM;
    }
    fill_idx_to_uidx(indiv_exclude, unfiltered_indiv_ct, indiv_ct, indiv_idx_to_uidx);
  }
  sorted_ids = (char*)top_alloc(&topsize, indiv_ct * max_person_id_len);
  if (!sorted_ids) {
    goto load_covars_ret_NOMEM;
  }
  id_map = (uint32_t*)top_alloc(&topsize, indiv_ct * sizeof(int32_t));
  if (!id_map) {
    goto load_covars_ret_NOMEM;
  }
  already_seen = (uintptr_t*)top_alloc(&topsize, indiv_ctl * sizeof(intptr_t));
  if (!already_seen) {
    goto load_covars_ret_NOMEM;
  }
  fill_ulong_zero(already_seen, indiv_ctl);
  if (covar_modifier & COVAR_NAME) {
    sorted_covar_name_flag_ids = (char*)top_alloc(&topsize, covar_range_list_ptr->name_ct * covar_range_list_ptr->name_max_len);
    if (!sorted_covar_name_flag_ids) {
      goto load_covars_ret_NOMEM;
    }
    covar_name_flag_id_map = (uint32_t*)top_alloc(&topsize, covar_range_list_ptr->name_ct * sizeof(int32_t));
    if (!covar_name_flag_id_map) {
      goto load_covars_ret_NOMEM;
    }
    covar_name_flag_seen_idxs = (int32_t*)top_alloc(&topsize, covar_range_list_ptr->name_ct * sizeof(int32_t));
    if (!covar_name_flag_seen_idxs) {
      goto load_covars_ret_NOMEM;
    }

    wkspace_left -= topsize;
    // kludge to use sort_item_ids_noalloc()
    fill_ulong_zero((uintptr_t*)covar_name_flag_seen_idxs, (covar_range_list_ptr->name_ct + (BITCT - 1)) / BITCT);
    retval = sort_item_ids_noalloc(sorted_covar_name_flag_ids, covar_name_flag_id_map, covar_range_list_ptr->name_ct, (uintptr_t*)covar_name_flag_seen_idxs, covar_range_list_ptr->name_ct, covar_range_list_ptr->names, covar_range_list_ptr->name_max_len, 0, 0, strcmp_deref);
    if (retval) {
      wkspace_left += topsize;
      if (retval == RET_INVALID_FORMAT) {
	logprint("(in --covar-name parameter sequence)\n");
	retval = RET_INVALID_CMDLINE;
      }
      goto load_covars_ret_1;
    }
    fill_int_one(covar_name_flag_seen_idxs, covar_range_list_ptr->name_ct);
  } else {
    wkspace_left -= topsize;
  }
  retval = sort_item_ids_noalloc(sorted_ids, id_map, unfiltered_indiv_ct, indiv_exclude, indiv_ct, person_ids, max_person_id_len, 0, 1, strcmp_deref);
  wkspace_left += topsize;
  if (retval) {
    goto load_covars_ret_1;
  }

  // To simplify loading sequence, guarantee enough space for covars_active[]
  // bitfield on first pass.  Each covariate corresponds to at least 16 bits in
  // the first nonempty line (a value and a space = 2 bytes), so reserving the
  // first 1/17 (rounded up) always works.
  loadbuf_size = ((wkspace_left - topsize) / 68) * 64;
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  } else if (loadbuf_size <= MAXLINELEN) {
    goto load_covars_ret_NOMEM;
  }
  loadbuf = (char*)(&(wkspace_base[wkspace_left - topsize - loadbuf_size]));
  retval = open_and_load_to_first_token(&covar_file, covar_fname, loadbuf_size, '\0', "--covar file", loadbuf, &bufptr);
  if (retval) {
    goto load_covars_ret_1;
  }
  covar_raw_ct = count_tokens(bufptr);
  if (covar_raw_ct < 3) {
    goto load_covars_ret_INVALID_FORMAT_2;
  } else if (covar_raw_ct < 2 + gxe_mcovar) {
    logprint("Error: Fewer tokens in --covar file than required by --gxe parameter.\n");
    goto load_covars_ret_INVALID_CMDLINE;
  }
  covar_raw_ct -= 2;
  covar_raw_ctl = (covar_raw_ct + (BITCT - 1)) / BITCT;
  covars_active = (uintptr_t*)wkspace_alloc(covar_raw_ctl * sizeof(intptr_t));

  // no header line present?
  bufptr2 = next_item(bufptr);
  header_absent = (strcmp_se(bufptr, "FID", 3) || strcmp_se(bufptr2, "IID", 3));
  bufptr = next_item(bufptr2);

  if (covar_modifier & (COVAR_NAME | COVAR_NUMBER)) {
    fill_ulong_zero(covars_active, covar_raw_ctl);
    if (covar_modifier & COVAR_NUMBER) {
      if (numeric_range_list_to_bitfield(covar_range_list_ptr, covar_raw_ct, covars_active, 1, 0)) {
	logprint("Error: Fewer tokens in --covar file than required by --covar-number parameters.\n");
	goto load_covars_ret_INVALID_CMDLINE;
      }
    } else {
      if (header_absent) {
	logprint("Error: --covar file doesn't have a header line for --covar-name.\n");
	goto load_covars_ret_INVALID_CMDLINE;
      }
      retval = string_range_list_to_bitfield(bufptr, covar_raw_ct, covar_range_list_ptr, sorted_covar_name_flag_ids, covar_name_flag_id_map, covar_name_flag_seen_idxs, "covar-name", "--covar file header line", covars_active);
      if (retval) {
	goto load_covars_ret_1;
      }
      // deallocate --covar-name support
      topsize -= (uintptr_t)(((unsigned char*)already_seen) - ((unsigned char*)covar_name_flag_seen_idxs));
    }
    covar_ct = popcount_longs(covars_active, 0, covar_raw_ctl);
  } else if (covar_range_list_ptr) {
    fill_all_bits(covars_active, covar_raw_ct);
    covar_ct = covar_raw_ct;
  } else {
    fill_ulong_zero(covars_active, covar_raw_ctl);
    covar_ct = 0;
  }
  min_covar_col_ct = last_set_bit(covars_active, covar_raw_ctl) + 1;
  if (min_covar_col_ct < gxe_mcovar) {
    min_covar_col_ct = gxe_mcovar;
  }
  if (header_absent) {
    max_covar_name_len = 4 + intlen(min_covar_col_ct);
  } else {
    uii = 0;
    while (1) {
      bufptr2 = item_endnn(bufptr);
      if (IS_SET(covars_active, uii)) {
        if (max_covar_name_len <= (uintptr_t)(bufptr2 - bufptr)) {
	  max_covar_name_len = 1 + (uintptr_t)(bufptr2 - bufptr);
	}
      }
      if (++uii == min_covar_col_ct) {
	break;
      }
      bufptr = skip_initial_spaces(&(bufptr2[1]));
    }
  }

  wkspace_left -= topsize;
  // * covar_nm does NOT have a separate entry per covariate; instead,
  //   if a single covariate is missing for a person, that person's covar_nm
  //   bit is zero.
  // * covar_d is in individual-major order (i.e. the value of covariate m for
  //   individual n is covar_d[n * covar_ct + m], where m and n are
  //   zero-based).  It does track when some covariates are missing and others
  //   aren't (missing covariates are represented as the --missing-phenotype
  //   value).
  if (covar_range_list_ptr) {
    // not only --gxe
    *covar_ct_ptr = covar_ct;
    *max_covar_name_len_ptr = max_covar_name_len;
    if (wkspace_alloc_c_checked(covar_names_ptr, covar_ct * max_covar_name_len) ||
        wkspace_alloc_ul_checked(covar_nm_ptr, indiv_ctl * sizeof(intptr_t)) ||
        wkspace_alloc_d_checked(covar_d_ptr, covar_ct * indiv_ct * sizeof(double))) {
      goto load_covars_ret_NOMEM2;
    }
    covar_names = *covar_names_ptr;
    covar_nm = *covar_nm_ptr;
    covar_d = *covar_d_ptr;
    fill_ulong_zero(covar_nm, indiv_ctl);
  }
  if (gxe_mcovar) {
    if (wkspace_alloc_ul_checked(gxe_covar_nm_ptr, indiv_ctl * sizeof(intptr_t)) ||
        wkspace_alloc_ul_checked(gxe_covar_c_ptr, indiv_ctl * sizeof(intptr_t))) {
      goto load_covars_ret_NOMEM2;
    }
    gxe_covar_nm = *gxe_covar_nm_ptr;
    gxe_covar_c = *gxe_covar_c_ptr;
    fill_ulong_zero(gxe_covar_nm, indiv_ctl);
    fill_ulong_zero(gxe_covar_c, indiv_ctl);
  }
  if (wkspace_left <= MAXLINELEN) {
    goto load_covars_ret_NOMEM2;
  }
  wkspace_mark2 = wkspace_base;
  loadbuf = (char*)wkspace_base;
  loadbuf_size = wkspace_left;
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  }
  wkspace_left += topsize;
  loadbuf[loadbuf_size - 1] = ' ';

  rewind(covar_file);
  if (header_absent) {
    if (covar_range_list_ptr) {
      covar_uidx = 0;
      for (covar_idx = 0; covar_idx < covar_ct; covar_idx++) {
	covar_uidx = next_set_ul_unsafe(covars_active, covar_uidx);
	uint32_writex(memcpyl3a(&(covar_names[covar_idx * max_covar_name_len]), "COV"), ++covar_uidx, '\0');
      }
    }
  } else {
    covar_idx = 0;
    retval = load_to_first_token(covar_file, loadbuf_size, '\0', "--covar file", loadbuf, &bufptr);
    if (retval) {
      goto load_covars_ret_1;
    }
    if (covar_range_list_ptr) {
      bufptr2 = item_endnn(next_item(bufptr));
      for (uii = 0; uii < min_covar_col_ct; uii++) {
	bufptr = skip_initial_spaces(bufptr2);
	bufptr2 = item_endnn(bufptr);
	if (IS_SET(covars_active, uii)) {
	  memcpyx(&(covar_names[covar_idx * max_covar_name_len]), bufptr, (uintptr_t)(bufptr2 - bufptr), '\0');
	  covar_idx++;
	}
      }
    }
  }
  while (fgets(loadbuf, loadbuf_size, covar_file)) {
    if (!loadbuf[loadbuf_size - 1]) {
      if (loadbuf_size == MAXLINEBUFLEN) {
	logprint("Error: Pathologically long line in --covar file.\n");
	goto load_covars_ret_INVALID_FORMAT;
      } else {
	goto load_covars_ret_NOMEM;
      }
    }
    bufptr = skip_initial_spaces(loadbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    bufptr2 = next_item(bufptr);
    if (!bufptr2) {
      goto load_covars_ret_INVALID_FORMAT_2;
    }
    ii = bsearch_fam_indiv(tbuf, sorted_ids, max_person_id_len, indiv_ct, bufptr, bufptr2);
    if (ii == -1) {
      continue;
    }
    if (is_set(already_seen, ii)) {
      bufptr[strlen_se(bufptr)] = '\0';
      bufptr2[strlen_se(bufptr2)] = '\0';
      sprintf(logbuf, "Error: %s %s appears multiple times in --covar file.\n", bufptr, bufptr2);
      logprintb();
      goto load_covars_ret_INVALID_FORMAT;
    }
    set_bit(already_seen, ii);
    indiv_idx = id_map[(uint32_t)ii];
    if (!next_item_mult(bufptr2, min_covar_col_ct)) {
      goto load_covars_ret_INVALID_FORMAT_2;
    }
    bufptr2 = item_endnn(bufptr2);
    if (covar_range_list_ptr) {
      dptr = &(covar_d[indiv_idx * covar_ct]);
    }
    covar_missing = 0;
    for (uii = 0; uii < min_covar_col_ct; uii++) {
      bufptr = skip_initial_spaces(bufptr2);
      bufptr2 = item_endnn(bufptr);
      if (IS_SET(covars_active, uii)) {
        if (scan_double(bufptr, &dxx)) {
	  covar_missing = 1;
	  dxx = missing_phenod;
	} else if (dxx == missing_phenod) {
	  covar_missing = 1;
	}
	*dptr++ = dxx;
      }
      if (uii + 1 == gxe_mcovar) {
        indiv_uidx = indiv_idx_to_uidx[indiv_idx];
	if (IS_SET(pheno_nm, indiv_uidx)) {
	  if (scan_double(bufptr, &dxx) || (dxx == missing_phenod)) {
	    // PLINK 1.07 quirk: --gxe treats 0 the same as -9/nonnumeric, but
	    // --write-covar does not, so handle 0 separately for backward
	    // compatibility
	    if (!keep_pheno_on_missing_cov) {
	      CLEAR_BIT(pheno_nm, indiv_uidx);
	      *pheno_nm_ct_ptr -= 1;
	    }
	  } else if (dxx != 0.0) {
	    SET_BIT(gxe_covar_nm, indiv_idx);
	    if (dxx == 2.0) {
	      SET_BIT(gxe_covar_c, indiv_idx);
	    }
	  }
	}
      }
    }
    if (covar_range_list_ptr) {
      if (!covar_missing) {
	SET_BIT(covar_nm, indiv_idx);
      } else {
	missing_cov_ct++;
	if (!keep_pheno_on_missing_cov) {
	  indiv_uidx = indiv_idx_to_uidx[indiv_idx];
	  if (IS_SET(pheno_nm, indiv_uidx)) {
	    CLEAR_BIT(pheno_nm, indiv_uidx);
	    *pheno_nm_ct_ptr -= 1;
	  }
	}
      }
    }
    loaded_indiv_ct++;
  }
  if (!feof(covar_file)) {
    goto load_covars_ret_READ_FAIL;
  }
  if (loaded_indiv_ct == missing_cov_ct) {
    logprint("Error: No --covar values loaded.\n");
    goto load_covars_ret_INVALID_FORMAT;
  }
  if (covar_range_list_ptr) {
    if ((covar_ct < covar_raw_ct - 1) || ((covar_ct == covar_raw_ct - 1) && ((!gxe_mcovar) || is_set(covars_active, gxe_mcovar - 1)))) {
      if (gxe_mcovar && (!is_set(covars_active, gxe_mcovar - 1))) {
        sprintf(logbuf, "--covar: 1 C/C cov. loaded for --gxe, %" PRIuPTR "/%" PRIuPTR " for other operations.\n", covar_ct, covar_raw_ct);
      } else {
        sprintf(logbuf, "--covar: %" PRIuPTR " out of %" PRIuPTR " covariates loaded.\n", covar_ct, covar_raw_ct);
      }
    } else {
      sprintf(logbuf, "--covar: %" PRIuPTR " covariate%s loaded.\n", covar_ct, (covar_ct == 1)? "" : "s");
    }
    logprintb();
  } else {
    logprint("--covar: 1 case/control covariate loaded for --gxe.\n");
  }
  if (missing_cov_ct) {
    sprintf(logbuf, "%" PRIuPTR " %s had missing value(s).\n", missing_cov_ct, species_str(missing_cov_ct));
    logprintb();
  }

  wkspace_reset(wkspace_mark2);
  while (0) {
  load_covars_ret_NOMEM2:
    wkspace_left += topsize;
  load_covars_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  load_covars_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  load_covars_ret_INVALID_FORMAT_2:
    logprint("Error: Fewer tokens than expected in --covar file line.\n");
  load_covars_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  load_covars_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  }
 load_covars_ret_1:
  if (retval) {
    wkspace_reset(wkspace_mark);
  }
  fclose_cond(covar_file);
  return retval;
}

int32_t write_covars(char* outname, char* outname_end, uint32_t write_covar_modifier, uint32_t write_covar_dummy_max_categories, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_ct, char* person_ids, uintptr_t max_person_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* sex_nm, uintptr_t* sex_male, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, double missing_phenod, char* output_missing_pheno, uintptr_t covar_ct, char* covar_names, uintptr_t max_covar_name_len, uintptr_t* covar_nm, double* covar_d) {
  FILE* outfile = NULL;
  uint32_t write_pheno = write_covar_modifier & WRITE_COVAR_PHENO;
  uint32_t exclude_parents = write_covar_modifier & WRITE_COVAR_NO_PARENTS;
  uint32_t exclude_sex = write_covar_modifier & WRITE_COVAR_NO_SEX;
  uint32_t female_2 = write_covar_modifier & WRITE_COVAR_FEMALE_2;
  uintptr_t indiv_uidx = 0;
  uint32_t* downcoding_level = NULL;
  uint32_t* downcoding_values = NULL;
  char* zbuf = NULL;
  char* out_missing_buf = NULL;
  uintptr_t omplen_p1 = strlen(output_missing_pheno) + 1;
  uint32_t do_downcoding = (write_covar_modifier & WRITE_COVAR_DUMMY) && (indiv_ct > 2);
  uint32_t downcoding_no_round = (write_covar_modifier & WRITE_COVAR_DUMMY_NO_ROUND);
  uintptr_t downcoding_covar_ct = 0;
  uintptr_t covar_nm_ct = 0;
  uintptr_t topsize = 0;
  int32_t retval = 0;
  uint32_t* downcoding_buf_idxs;
  int64_t* sorted_downcoding_intbuf;
  int64_t* category_idx_sort_buf;
  uint32_t* category_remap;
  uint32_t* uiptr;
  double* sorted_downcoding_buf;
  double* dptr;
  char* downcoding_string_buf;
  char* wptr_start;
  char* wptr;
  char* bufptr;
  uintptr_t covar_idx;
  uintptr_t indiv_idx;
  uintptr_t indiv_idx2;
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t slen;
  uint64_t ullii;
  uint64_t ulljj;
  uintptr_t downcode_category_ct;
  uintptr_t downcode_idx;
  double dxx;
  uint32_t uii;
  uint32_t ujj;
  memcpy(outname_end, ".cov", 5);
  if (fopen_checked(&outfile, outname, "w")) {
    goto write_covars_ret_OPEN_FAIL;
  }
  if (fputs_checked("FID IID ", outfile)) {
    goto write_covars_ret_WRITE_FAIL;
  }
  if (write_pheno) {
    if (!exclude_parents) {
      fputs("PAT MAT ", outfile);
    }
    if (!exclude_sex) {
      fputs("SEX PHENOTYPE ", outfile);
    } else {
      fputs("PHENOTYPE ", outfile);
    }
  }
  if (do_downcoding) {
    // could make downcoding_values allocation incremental (top_alloc() calls
    // have been arranged to make this a simple change; would just need to
    // wrap the qsort_ext() calls)
    if (wkspace_alloc_ui_checked(&downcoding_level, covar_ct * sizeof(int32_t)) ||
        wkspace_alloc_ui_checked(&downcoding_values, covar_ct * indiv_ct * sizeof(int32_t))) {
      goto write_covars_ret_NOMEM;
    }
    if (write_covar_dummy_max_categories > indiv_ct) {
      write_covar_dummy_max_categories = indiv_ct;
    }
    downcoding_string_buf = (char*)top_alloc(&topsize, 16 * (write_covar_dummy_max_categories + 1));
    if (!downcoding_string_buf) {
      goto write_covars_ret_NOMEM;
    }
    category_idx_sort_buf = (int64_t*)top_alloc(&topsize, write_covar_dummy_max_categories * sizeof(int64_t));
    if (!category_idx_sort_buf) {
      goto write_covars_ret_NOMEM;
    }
    category_remap = (uint32_t*)top_alloc(&topsize, write_covar_dummy_max_categories * sizeof(int32_t));
    if (!category_idx_sort_buf) {
      goto write_covars_ret_NOMEM;
    }
    uiptr = downcoding_values;
    if (!downcoding_no_round) {
      sorted_downcoding_intbuf = (int64_t*)top_alloc(&topsize, indiv_ct * sizeof(int64_t));
      if (!sorted_downcoding_intbuf) {
        goto write_covars_ret_NOMEM;
      }
      for (covar_idx = 0; covar_idx < covar_ct; covar_idx++) {
        dptr = &(covar_d[covar_idx]);
	indiv_idx2 = 0;
        for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
	  dxx = dptr[indiv_idx * covar_ct];
	  if (dxx != missing_phenod) {
            // bits 0-30: original position
            // bits 31-62: 2^31 + (int32_t)phenotype
            // top bit zero to defend against >= 2^31 differences
            sorted_downcoding_intbuf[indiv_idx2++] = (int64_t)((((uint64_t)(((uint32_t)((int32_t)dxx)) ^ 0x80000000U)) << 31) | ((uint64_t)indiv_idx));
	  }
	}
	downcode_category_ct = 0;
	if (indiv_idx2 > 2) {
	  covar_nm_ct = indiv_idx2;
	  fill_uint_one(uiptr, indiv_ct);
#ifdef __cplusplus
	  std::sort(sorted_downcoding_intbuf, &(sorted_downcoding_intbuf[covar_nm_ct]));
#else
          qsort(sorted_downcoding_intbuf, covar_nm_ct, sizeof(int64_t), llcmp);
#endif

          ullii = ((uint64_t)sorted_downcoding_intbuf[0]);
	  ulii = ((uintptr_t)ullii) & (ONELU * 0x7fffffff);
	  ullii >>= 31;
          uiptr[ulii] = 0;
          for (indiv_idx2 = 1; indiv_idx2 < covar_nm_ct; indiv_idx2++) {
	    ulljj = (uint64_t)sorted_downcoding_intbuf[indiv_idx2];
	    uljj = ((uintptr_t)ulljj) & (ONELU * 0x7fffffff);
	    ulljj >>= 31;
	    if (ulljj != ullii) {
	      if (downcode_category_ct == write_covar_dummy_max_categories - 1) {
		downcode_category_ct = 0;
		break;
	      }
	      // save phenotype string
	      int32_writex(&(downcoding_string_buf[16 * downcode_category_ct]), (int32_t)(((uint32_t)ullii) ^ 0x80000000U), '\0');

	      // bits 0-31: initial category assignment
	      // bits 32-63: smallest indiv_idx2
	      category_idx_sort_buf[downcode_category_ct] = (int64_t)((((uint64_t)ulii) << 32) | ((uint64_t)downcode_category_ct));
	      downcode_category_ct++;
	      ulii = uljj;
	      ullii = ulljj;
	    } else {
	      if (uljj < ulii) {
		ulii = uljj;
	      }
	    }
	    uiptr[uljj] = downcode_category_ct;
	  }
	}

	// probably want to make this part its own function since it's
	// practically identical when downcoding_no_round is set
        if (downcode_category_ct > 1) {
	  int32_writex(&(downcoding_string_buf[16 * downcode_category_ct]), (int32_t)(((uint32_t)ullii) ^ 0x80000000U), '\0');
	  category_idx_sort_buf[downcode_category_ct] = (int64_t)((((uint64_t)ulii) << 32) | ((uint64_t)downcode_category_ct));
	  downcode_category_ct++;
          // now recover PLINK 1.07 category order
#ifdef __cplusplus
	  std::sort(category_idx_sort_buf, &(category_idx_sort_buf[downcode_category_ct]));
#else
          qsort(category_idx_sort_buf, downcode_category_ct, sizeof(int64_t), llcmp);
#endif

	  downcoding_level[covar_idx] = downcode_category_ct;
          wptr_start = strcpyax(tbuf, &(covar_names[covar_idx * max_covar_name_len]), '_');
	  for (downcode_idx = 0; downcode_idx < downcode_category_ct; downcode_idx++) {
	    uii = (uint32_t)(category_idx_sort_buf[downcode_idx]);
	    if (downcode_idx) {
	      wptr = strcpyax(wptr_start, &(downcoding_string_buf[16 * uii]), ' ');
	      fwrite(tbuf, 1, wptr - tbuf, outfile);
	    }
            category_remap[uii] = downcode_idx;
	  }
          for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
	    if (uiptr[indiv_idx] != 0xffffffff) {
              uiptr[indiv_idx] = category_remap[uiptr[indiv_idx]];
	    }
	  }
	  downcoding_covar_ct++;
	  uiptr = &(uiptr[indiv_ct]);
	} else {
          downcoding_level[covar_idx] = 0;
          fputs(&(covar_names[covar_idx * max_covar_name_len]), outfile);
	  putc(' ', outfile);
	}
      }
    } else {
      downcoding_buf_idxs = (uint32_t*)top_alloc(&topsize, indiv_ct * sizeof(int32_t));
      if (!downcoding_buf_idxs) {
	goto write_covars_ret_NOMEM;
      }
      sorted_downcoding_buf = (double*)top_alloc(&topsize, indiv_ct * sizeof(double));
      if (!sorted_downcoding_buf) {
	goto write_covars_ret_NOMEM;
      }
      for (covar_idx = 0; covar_idx < covar_ct; covar_idx++) {
	dptr = &(covar_d[covar_idx]);
	indiv_idx2 = 0;
	for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
	  dxx = dptr[indiv_idx * covar_ct];
	  if (dxx != missing_phenod) {
	    sorted_downcoding_buf[indiv_idx2] = dxx;
	    downcoding_buf_idxs[indiv_idx2++] = indiv_idx;
	  }
	}
	downcode_category_ct = 0;
	if (indiv_idx2 > 2) {
	  covar_nm_ct = indiv_idx2;
	  fill_uint_one(uiptr, indiv_ct);
	  if (qsort_ext((char*)sorted_downcoding_buf, covar_nm_ct, sizeof(double), double_cmp_deref, (char*)downcoding_buf_idxs, sizeof(int32_t))) {
	    goto write_covars_ret_NOMEM;
	  }
	  wptr_start = double_g_write(downcoding_string_buf, sorted_downcoding_buf[0]);
	  slen = (uintptr_t)(wptr_start - downcoding_string_buf);
	  *wptr_start = '\0';
	  wptr_start = downcoding_string_buf;
	  uii = downcoding_buf_idxs[0];
	  uiptr[uii] = 0;
	  bufptr = &(downcoding_string_buf[16]);
	  for (indiv_idx2 = 1; indiv_idx2 < covar_nm_ct; indiv_idx2++) {
	    // a bit inefficient, but this is a safe way to achieve "two
	    // doubles are equal if they yield the same printf %g output"
	    // behavior
	    wptr = double_g_write(bufptr, sorted_downcoding_buf[indiv_idx2]);
	    ujj = downcoding_buf_idxs[indiv_idx2];
	    if (((uintptr_t)(wptr - bufptr) != slen) || memcmp(wptr_start, bufptr, slen)) {
	      *wptr = '\0';
	      if (downcode_category_ct == write_covar_dummy_max_categories - 1) {
		downcode_category_ct = 0;
		break;
	      }
	      category_idx_sort_buf[downcode_category_ct] = (int64_t)((((uint64_t)uii) << 32) | ((uint64_t)downcode_category_ct));
	      downcode_category_ct++;
	      slen = (uintptr_t)(wptr - bufptr);
	      *wptr = '\0';
	      wptr_start = bufptr;
              uii = ujj;
	      bufptr = &(bufptr[16]);
	    } else {
	      if (ujj < uii) {
		uii = ujj;
	      }
	    }
	    uiptr[ujj] = downcode_category_ct;
	  }
	}
	if (downcode_category_ct > 1) {
	  *wptr = '\0';
	  category_idx_sort_buf[downcode_category_ct] = (int64_t)((((uint64_t)uii) << 32) | ((uint64_t)downcode_category_ct));
	  downcode_category_ct++;

#ifdef __cplusplus
	  std::sort(category_idx_sort_buf, &(category_idx_sort_buf[downcode_category_ct]));
#else
          qsort(category_idx_sort_buf, downcode_category_ct, sizeof(int64_t), llcmp);
#endif

	  downcoding_level[covar_idx] = downcode_category_ct;
          wptr_start = strcpyax(tbuf, &(covar_names[covar_idx * max_covar_name_len]), '_');
          for (downcode_idx = 0; downcode_idx < downcode_category_ct; downcode_idx++) {
	    uii = (uint32_t)(category_idx_sort_buf[downcode_idx]);
	    if (downcode_idx) {
	      wptr = strcpyax(wptr_start, &(downcoding_string_buf[16 * uii]), ' ');
	      fwrite(tbuf, 1, wptr - tbuf, outfile);
	    }
	    category_remap[uii] = downcode_idx;
          }
	  for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
	    if (uiptr[indiv_idx] != 0xffffffffU) {
	      uiptr[indiv_idx] = category_remap[uiptr[indiv_idx]];
	    }
	  }
	  downcoding_covar_ct++;
	  uiptr = &(uiptr[indiv_ct]);
	} else {
	  downcoding_level[covar_idx] = 0;
          fputs(&(covar_names[covar_idx * max_covar_name_len]), outfile);
	  putc(' ', outfile);
	}
      }
    }
    wkspace_reset((unsigned char*)downcoding_values);
    downcoding_values = (uint32_t*)wkspace_alloc(downcoding_covar_ct * indiv_ct * sizeof(int32_t));
    // topsize = 0;

    // (write_covar_dummy_max_categories - 1) columns, then divide by two
    // rounding up; the -1 and +1 cancel
    ujj = write_covar_dummy_max_categories / 2;
    if (wkspace_alloc_c_checked(&zbuf, ujj * sizeof(int32_t)) ||
        wkspace_alloc_c_checked(&out_missing_buf, (write_covar_dummy_max_categories - 1) * omplen_p1)) {
      goto write_covars_ret_NOMEM;
    }
    uiptr = (uint32_t*)zbuf;
    for (uii = 0; uii < ujj; uii++) {
      *uiptr++ = 0x20302030; // "0 0 "
    }
    bufptr = memcpyax(out_missing_buf, output_missing_pheno, omplen_p1 - 1, ' ');
    for (uii = 2; uii < write_covar_dummy_max_categories; uii++) {
      bufptr = memcpya(bufptr, out_missing_buf, omplen_p1);
    }
  } else {
    for (covar_idx = 0; covar_idx < covar_ct; covar_idx++) {
      fputs(&(covar_names[covar_idx * max_covar_name_len]), outfile);
      putc(' ', outfile);
    }
  }
  if (putc_checked('\n', outfile)) {
    goto write_covars_ret_WRITE_FAIL;
  }
  for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_uidx++, indiv_idx++) {
    indiv_uidx = next_unset_ul_unsafe(indiv_exclude, indiv_uidx);
    bufptr = &(person_ids[indiv_uidx * max_person_id_len]);
    wptr = (char*)memchr(bufptr, '\t', max_person_id_len);
    *wptr = ' ';
    if (fputs_checked(bufptr, outfile)) {
      goto write_covars_ret_WRITE_FAIL;
    }
    *wptr = '\t';
    putc(' ', outfile);
    if (write_pheno) {
      if (!exclude_parents) {
	fputs(&(paternal_ids[indiv_uidx * max_paternal_id_len]), outfile);
	putc(' ', outfile);
	fputs(&(maternal_ids[indiv_uidx * max_maternal_id_len]), outfile);
	putc(' ', outfile);
      }
      if (!exclude_sex) {
	if (!female_2) {
          putc((uint32_t)(48 + IS_SET(sex_male, indiv_uidx)), outfile);
	} else {
          if (IS_SET(sex_nm, indiv_uidx)) {
	    putc((uint32_t)(50 - IS_SET(sex_male, indiv_uidx)), outfile);
	  } else {
	    putc('0', outfile);
	  }
	}
        putc(' ', outfile);
      }
      if (IS_SET(pheno_nm, indiv_uidx)) {
        wptr = double_g_write(tbuf, pheno_d[indiv_uidx]);
	fwrite(tbuf, 1, wptr - tbuf, outfile);
      } else {
        fputs(output_missing_pheno, outfile);
      }
      putc(' ', outfile);
    }
    if (do_downcoding) {
      dptr = &(covar_d[indiv_idx * covar_ct]);
      downcode_idx = 0;
      uiptr = &(downcoding_values[indiv_idx]);
      for (covar_idx = 0; covar_idx < covar_ct; covar_idx++) {
	ujj = downcoding_level[covar_idx];
	dxx = *dptr++;
	if (dxx != missing_phenod) {
	  if (ujj) {
	    uii = uiptr[downcode_idx * indiv_ct];
	    if (uii) {
	      if (uii > 1) {
		fwrite(zbuf, 1, 2 * (uii - 1), outfile);
	      }
	      fputs("1 ", outfile);
	      if (uii < ujj - 1) {
		fwrite(zbuf, 1, 2 * (ujj - 1 - uii), outfile);
	      }
	    } else {
	      fwrite(zbuf, 1, 2 * (ujj - 1), outfile);
	    }
	    downcode_idx++;
	  } else {
	    wptr = double_g_writex(tbuf, dxx, ' ');
	    fwrite(tbuf, 1, wptr - tbuf, outfile);
	  }
	} else {
	  if (ujj) {
            fwrite(out_missing_buf, 1, omplen_p1 * (ujj - 1), outfile);
	  } else {
            fputs(output_missing_pheno, outfile);
            putc(' ', outfile);
	  }
	}
      }
    } else {
      dptr = &(covar_d[indiv_idx * covar_ct]);
      for (covar_idx = 0; covar_idx < covar_ct; covar_idx++) {
	wptr = double_g_writex(tbuf, dptr[covar_idx], ' ');
        fwrite(tbuf, 1, wptr - tbuf, outfile);
      }
    }
    if (putc_checked('\n', outfile)) {
      goto write_covars_ret_WRITE_FAIL;
    }
  }
  sprintf(logbuf, "Covariates written to %s.\n", outname);
  logprintb();
  while (0) {
  write_covars_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  write_covars_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  write_covars_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  fclose_cond(outfile);
  return retval;
}

int32_t write_fam(char* outname, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_ct, char* person_ids, uintptr_t max_person_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* sex_nm, uintptr_t* sex_male, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, char* output_missing_pheno, char delim, uint32_t* indiv_sort_map) {
  FILE* outfile = NULL;
  uintptr_t indiv_uidx = 0;
  uintptr_t indiv_uidx2 = 0;
  uintptr_t omplen = strlen(output_missing_pheno);
  int32_t retval = 0;
  char* cptr;
  char* bufptr;
  uintptr_t indiv_idx;
  uintptr_t clen;
  if (fopen_checked(&outfile, outname, "w")) {
    goto write_fam_ret_OPEN_FAIL;
  }
  for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
    if (indiv_sort_map) {
      do {
	indiv_uidx = indiv_sort_map[indiv_uidx2++];
      } while (IS_SET(indiv_exclude, indiv_uidx));
    } else {
      next_unset_ul_unsafe_ck(indiv_exclude, &indiv_uidx);
    }
    cptr = &(person_ids[indiv_uidx * max_person_id_len]);
    clen = strlen_se(cptr);
    bufptr = memcpyax(tbuf, cptr, clen, delim);
    bufptr = strcpyax(bufptr, &(cptr[clen + 1]), delim);
    if (paternal_ids) {
      bufptr = strcpya(bufptr, &(paternal_ids[indiv_uidx * max_paternal_id_len]));
    } else {
      *bufptr++ = '0';
    }
    *bufptr++ = delim;
    if (maternal_ids) {
      bufptr = strcpya(bufptr, &(maternal_ids[indiv_uidx * max_maternal_id_len]));
    } else {
      *bufptr++ = '0';
    }
    *bufptr++ = delim;
    *bufptr++ = sexchar(sex_nm, sex_male, indiv_uidx);
    *bufptr++ = delim;
    if (!IS_SET(pheno_nm, indiv_uidx)) {
      bufptr = memcpya(bufptr, output_missing_pheno, omplen);
    } else if (pheno_c) {
      *bufptr++ = '1' + IS_SET(pheno_c, indiv_uidx);
    } else {
      bufptr = double_g_write(bufptr, pheno_d[indiv_uidx]);
    }
    *bufptr++ = '\n';
    if (fwrite_checked(tbuf, bufptr - tbuf, outfile)) {
      goto write_fam_ret_WRITE_FAIL;
    }
    if (!indiv_sort_map) {
      indiv_uidx++;
    }
  }
  if (fclose_null(&outfile)) {
    goto write_fam_ret_WRITE_FAIL;
  }
  while (0) {
  write_fam_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  write_fam_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  return retval;
}

int32_t write_map_or_bim(char* outname, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, double* marker_cms, uint32_t* marker_pos, char** marker_allele_ptrs, char delim, uint32_t zero_extra_chroms, Chrom_info* chrom_info_ptr) {
  // write a .map if marker_allele_ptrs is NULL, .bim otherwise
  FILE* outfile = NULL;
  uintptr_t marker_uidx = 0;
  int32_t retval = 0;
  uint32_t chrom_end = 0;
  uint32_t chrom_fo_idx = 0xffffffffU;
  uint32_t chrom_idx = 0;
  char* buf_start = NULL;
  uintptr_t marker_idx;
  char* bufptr;
  if (fopen_checked(&outfile, outname, "w")) {
    goto write_map_or_bim_ret_OPEN_FAIL;
  }
  for (marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
    next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
    while (marker_uidx >= chrom_end) {
      chrom_idx = chrom_info_ptr->chrom_file_order[++chrom_fo_idx];
      chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx + 1];
      buf_start = chrom_name_write(tbuf, chrom_info_ptr, chrom_idx, zero_extra_chroms);
      *buf_start++ = delim;
    }
    bufptr = strcpyax(buf_start, &(marker_ids[marker_uidx * max_marker_id_len]), delim);
    if (!marker_cms) {
      *bufptr++ = '0';
    } else {
      bufptr = double_g_write(bufptr, marker_cms[marker_uidx]);
    }
    *bufptr++ = delim;
    bufptr = uint32_write(bufptr, marker_pos[marker_uidx]);
    if (fwrite_checked(tbuf, bufptr - tbuf, outfile)) {
      goto write_map_or_bim_ret_WRITE_FAIL;
    }
    if (marker_allele_ptrs) {
      putc(delim, outfile);
      fputs(marker_allele_ptrs[2 * marker_uidx], outfile);
      putc(delim, outfile);
      fputs(marker_allele_ptrs[2 * marker_uidx + 1], outfile);
    }
    putc('\n', outfile);
  }
  if (fclose_null(&outfile)) {
    goto write_map_or_bim_ret_WRITE_FAIL;
  }
  while (0) {
  write_map_or_bim_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  write_map_or_bim_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  fclose_cond(outfile);
  return retval;
}

int32_t load_bim_split_chrom(char* bimname, uintptr_t* marker_exclude, uintptr_t marker_ct, Chrom_info* chrom_info_ptr, int64_t* ll_buf) {
  FILE* infile = NULL;
  uint32_t marker_uidx = 0xffffffffU; // deliberate overflow
  uintptr_t marker_idx = 0;
  int32_t retval = 0;
  char* bufptr;
  uint64_t chrom_idx;
  if (fopen_checked(&infile, bimname, "r")) {
    goto load_bim_split_chrom_ret_OPEN_FAIL;
  }
  do {
  load_bim_split_chrom_reread:
    if (!fgets(tbuf, MAXLINELEN, infile)) {
      goto load_bim_split_chrom_ret_READ_FAIL;
    }
    bufptr = skip_initial_spaces(tbuf);
    if (is_eoln_or_comment(*bufptr)) {
      goto load_bim_split_chrom_reread;
    }
    marker_uidx++;
    if (IS_SET(marker_exclude, marker_uidx)) {
      goto load_bim_split_chrom_reread;
    }
    // already validated
    chrom_idx = ((uint32_t)get_chrom_code(chrom_info_ptr, bufptr));
    ll_buf[marker_idx] = (int64_t)((chrom_idx << 32) | ((uint64_t)marker_idx));
  } while ((++marker_idx) < marker_ct);
  while (0) {
  load_bim_split_chrom_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  load_bim_split_chrom_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  }
  fclose_cond(infile);
  return retval;
}

void fill_ll_buf(uintptr_t* marker_exclude, uintptr_t marker_ct, Chrom_info* chrom_info_ptr, int64_t* ll_buf) {
  uint32_t marker_uidx = 0;
  uint32_t chrom_idx_p1 = 0;
  uint32_t chrom_end = 0;
  uint64_t chrom_idx_shifted = 0;
  uintptr_t marker_idx;
  for (marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
    next_unset_unsafe_ck(marker_exclude, &marker_uidx);
    if (marker_uidx >= chrom_end) {
      do {
	chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[++chrom_idx_p1];
      } while (marker_uidx >= chrom_end);
      chrom_idx_shifted = ((uint64_t)(chrom_info_ptr->chrom_file_order[chrom_idx_p1 - 1])) << 32;
    }
    ll_buf[marker_idx] = (int64_t)(chrom_idx_shifted | ((uint64_t)marker_idx));
  }
}

int32_t update_marker_chroms(Two_col_params* update_chr, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t allow_extra_chroms, Chrom_info* chrom_info_ptr, int64_t* ll_buf) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* infile = NULL;
  char skipchar = update_chr->skipchar;
  uint32_t colid_first = (update_chr->colid < update_chr->colx);
  uint32_t marker_ctl = (marker_ct + (BITCT - 1)) / BITCT;
  uintptr_t hit_ct = 0;
  uintptr_t miss_ct = 0;
  uintptr_t marker_uidx = 0;
  char* sorted_marker_ids;
  uint32_t* marker_id_map;
  uintptr_t* already_seen;
  char* loadbuf;
  uintptr_t loadbuf_size;
  uint32_t colmin;
  uint32_t coldiff;
  char* colid_ptr;
  char* colx_ptr;
  uintptr_t marker_idx;
  uintptr_t delta;
  int32_t sorted_idx;
  char cc;
  int32_t retval;
  if (wkspace_alloc_ul_checked(&already_seen, marker_ctl * sizeof(intptr_t)) ||
      wkspace_alloc_ui_checked(&marker_id_map, marker_ct * sizeof(int32_t)) ||
      wkspace_alloc_c_checked(&sorted_marker_ids, marker_ct * max_marker_id_len)) {
    goto update_marker_chroms_ret_NOMEM;
  }
  for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
    marker_id_map[marker_idx] = marker_idx;
  }
  marker_uidx = 0;
  marker_idx = 0;
  do {
    marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
    delta = next_set_ul(marker_exclude, marker_uidx, unfiltered_marker_ct) - marker_uidx;
    memcpy(&(sorted_marker_ids[marker_idx * max_marker_id_len]), &(marker_ids[marker_uidx * max_marker_id_len]), delta * max_marker_id_len);
    marker_idx += delta;
    marker_uidx += delta;
  } while (marker_idx < marker_ct);
  if (qsort_ext(sorted_marker_ids, marker_ct, max_marker_id_len, strcmp_deref, (char*)marker_id_map, sizeof(int32_t))) {
    goto update_marker_chroms_ret_NOMEM;
  }
  colid_ptr = scan_for_duplicate_ids(sorted_marker_ids, marker_ct, max_marker_id_len);
  if (colid_ptr) {
    sprintf(logbuf, "Error: Duplicate variant ID %s.\n", colid_ptr);
    logprintb();
    goto update_marker_chroms_ret_INVALID_FORMAT;
  }

  fill_ulong_zero(already_seen, marker_ctl);
  loadbuf = (char*)wkspace_base;
  loadbuf_size = wkspace_left;
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  }
  if (loadbuf_size <= MAXLINELEN) {
    goto update_marker_chroms_ret_NOMEM;
  }
  retval = open_and_skip_first_lines(&infile, update_chr->fname, loadbuf, loadbuf_size, update_chr->skip);
  if (retval) {
    goto update_marker_chroms_ret_1;
  }
  if (colid_first) {
    colmin = update_chr->colid - 1;
    coldiff = update_chr->colx - update_chr->colid;
  } else {
    colmin = update_chr->colx - 1;
    coldiff = update_chr->colid - update_chr->colx;
  }
  while (fgets(loadbuf, loadbuf_size, infile)) {
    if (!(loadbuf[loadbuf_size - 1])) {
      goto update_marker_chroms_ret_NOMEM;
    }
    colid_ptr = skip_initial_spaces(loadbuf);
    cc = *colid_ptr;
    if (is_eoln_kns(cc) || (cc == skipchar)) {
      continue;
    }
    if (colid_first) {
      if (colmin) {
        colid_ptr = next_item_mult(colid_ptr, colmin);
      }
      colx_ptr = next_item_mult(colid_ptr, coldiff);
      if (no_more_items_kns(colx_ptr)) {
	goto update_marker_chroms_ret_INVALID_FORMAT_2;
      }
    } else {
      colx_ptr = colid_ptr;
      if (colmin) {
	colx_ptr = next_item_mult(colx_ptr, colmin);
      }
      colid_ptr = next_item_mult(colx_ptr, coldiff);
      if (no_more_items_kns(colid_ptr)) {
	goto update_marker_chroms_ret_INVALID_FORMAT_2;
      }
    }
    colid_ptr[strlen_se(colid_ptr)] = '\0';
    sorted_idx = bsearch_str(colid_ptr, sorted_marker_ids, max_marker_id_len, 0, marker_ct - 1);
    if (sorted_idx == -1) {
      miss_ct++;
      continue;
    }
    if (is_set(already_seen, sorted_idx)) {
      sprintf(logbuf, "Error: Duplicate variant %s in --update-chr file.\n", colid_ptr);
      logprintb();
      goto update_marker_chroms_ret_INVALID_FORMAT;
    }
    set_bit(already_seen, sorted_idx);
    marker_idx = marker_id_map[((uint32_t)sorted_idx)];
    sorted_idx = get_chrom_code(chrom_info_ptr, colx_ptr);
    if (sorted_idx == -1) {
      if (!allow_extra_chroms) {
	logprint("Error: Invalid chromosome code in --update-chr file.\n");
	goto update_marker_chroms_ret_INVALID_FORMAT;
      }
      retval = resolve_or_add_chrom_name(chrom_info_ptr, colx_ptr, &sorted_idx);
      if (retval) {
	goto update_marker_chroms_ret_1;
      }
    }
    ll_buf[marker_idx] = (int64_t)((((uint64_t)((uint32_t)sorted_idx)) << 32) | (((uint64_t)ll_buf[marker_idx]) & 0xffffffffLLU));
    hit_ct++;
  }
  if (!feof(infile)) {
    goto update_marker_chroms_ret_READ_FAIL;
  }
  if (miss_ct) {
    sprintf(logbuf, "--update-chr: %" PRIuPTR " value%s updated, %" PRIuPTR " variant ID%s not present.\n", hit_ct, (hit_ct == 1)? "" : "s", miss_ct, (miss_ct == 1)? "" : "s");
  } else {
    sprintf(logbuf, "--update-chr: %" PRIuPTR " value%s updated.\n", hit_ct, (hit_ct == 1)? "" : "s");
  }
  logprintb();
  while (0) {
  update_marker_chroms_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  update_marker_chroms_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  update_marker_chroms_ret_INVALID_FORMAT_2:
    logprint("Error: Fewer tokens than expected in --update-chr file line.\n");
  update_marker_chroms_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 update_marker_chroms_ret_1:
  fclose_cond(infile);
  wkspace_reset(wkspace_mark);
  return retval;
}

void sort_marker_chrom_pos(int64_t* ll_buf, uintptr_t marker_ct, uint32_t* pos_buf, uint32_t* chrom_start, uint32_t* chrom_id, uint32_t* chrom_ct_ptr) {
  // assumes ll_buf is filled with chromosome idxs in high 32 bits, and
  // internal marker indices in low 32 bits.  pos_buf is expected to have
  // base-pair positions.
  uintptr_t marker_idx;
  uint32_t uii;
  uint32_t cur_chrom;
  uint32_t chrom_ct;
#ifdef __cplusplus
  std::sort(ll_buf, &(ll_buf[marker_ct]));
#else
  qsort(ll_buf, marker_ct, sizeof(int64_t), llcmp);
#endif
  cur_chrom = ll_buf[0] >> 32;
  chrom_ct = 0;
  chrom_start[0] = 0;
  chrom_id[0] = cur_chrom;
  uii = (uint32_t)ll_buf[0];
  ll_buf[0] = ((uint64_t)uii) | (((uint64_t)pos_buf[uii]) << 32);
  for (marker_idx = 1; marker_idx < marker_ct; marker_idx++) {
    if ((ll_buf[marker_idx] >> 32) != cur_chrom) {
      cur_chrom = ll_buf[marker_idx] >> 32;
      chrom_start[++chrom_ct] = marker_idx;
      chrom_id[chrom_ct] = cur_chrom;
    }
    uii = (uint32_t)ll_buf[marker_idx];
    ll_buf[marker_idx] = ((uint64_t)uii) | (((uint64_t)pos_buf[uii]) << 32);
  }
  chrom_start[++chrom_ct] = marker_ct;
  for (uii = 0; uii < chrom_ct; uii++) {
#ifdef __cplusplus
    std::sort(&(ll_buf[chrom_start[uii]]), &(ll_buf[chrom_start[uii + 1]]));
#else
    qsort(&(ll_buf[chrom_start[uii]]), chrom_start[uii + 1] - chrom_start[uii], sizeof(int64_t), llcmp);
#endif
  }
  *chrom_ct_ptr = chrom_ct;
}

int32_t sort_and_write_bim(uint32_t* map_reverse, uint32_t map_cols, char* outname, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, double* marker_cms, uint32_t* marker_pos, char** marker_allele_ptrs, int64_t* ll_buf, uint32_t zero_extra_chroms, Chrom_info* chrom_info_ptr) {
  FILE* outfile = NULL;
  uint32_t max_code = chrom_info_ptr->max_code;
  int32_t retval = 0;
  const char* missing_geno_ptr = g_missing_geno_ptr;
  const char* output_missing_geno_ptr = g_output_missing_geno_ptr;
  char wbuf[16];
  uint32_t chrom_start[MAX_POSSIBLE_CHROM + 1];
  uint32_t chrom_id[MAX_POSSIBLE_CHROM];
  uint32_t* unpack_map;
  uintptr_t marker_uidx;
  uintptr_t marker_idx;
  uintptr_t marker_idx2;
  char* chrom_ptr;
  char* bufptr;
  uint32_t uii;
  uint32_t ujj;
  uint32_t cur_chrom;
  uint32_t chrom_ct;
  // There can be a LOT of markers (some 1000 Genomes files we've been offered
  // have ~40 million), so speeding up the sorting step over just calling
  // qsort_ext() may not be a complete waste of effort.
  // Strategy:
  // 1. fill ll_buf with chromosome idx in high-order bits, original position
  // in low-order.
  // 2. std::sort() ll_buf, read off chromosome boundaries
  // 3. then replace high-order bits in ll_buf with marker positions, and
  // std::sort() each chromosome separately.
  // Would be even faster if this was performed in a single sort, in the
  // super-common case where all three numbers can be squeezed together in 64
  // bits.  But we care most about performance when this can't be done, so I
  // haven't bothered with that optimization.
  if (wkspace_alloc_ui_checked(&unpack_map, marker_ct * sizeof(int32_t))) {
    goto sort_and_write_bim_ret_NOMEM;
  }
  fill_idx_to_uidx(marker_exclude, unfiltered_marker_ct, marker_ct, unpack_map);
  sort_marker_chrom_pos(ll_buf, marker_ct, marker_pos, chrom_start, chrom_id, &chrom_ct);
  if (fopen_checked(&outfile, outname, "w")) {
    goto sort_and_write_bim_ret_OPEN_FAIL;
  }

  marker_idx = 0;
  for (uii = 0; uii < chrom_ct; uii++) {
    cur_chrom = chrom_id[uii];
    ujj = chrom_start[uii + 1];
    if (cur_chrom <= max_code) {
      chrom_ptr = uint32_write(wbuf, cur_chrom);
      *chrom_ptr = '\0';
    } else if (zero_extra_chroms) {
      chrom_ptr = (char*)(&(g_one_char_strs[96]));
    } else {
      chrom_ptr = chrom_info_ptr->nonstd_names[cur_chrom];
    }
    for (; marker_idx < ujj; marker_idx++) {
      marker_idx2 = (uint32_t)ll_buf[marker_idx];
      marker_uidx = unpack_map[marker_idx2];
      fputs(chrom_ptr, outfile);
      putc('\t', outfile);
      fputs(&(marker_ids[marker_idx2 * max_marker_id_len]), outfile);
      wbuf[0] = '\t';
      if (!marker_cms) {
        bufptr = memcpya(&(wbuf[1]), "0\t", 2); 
      } else {
        bufptr = double_g_writex(&(wbuf[1]), marker_cms[marker_idx2], '\t');
      }
      fwrite(wbuf, 1, bufptr - wbuf, outfile);
      bufptr = uint32_write(wbuf, (uint32_t)(ll_buf[marker_idx] >> 32));
      fwrite(wbuf, 1, bufptr - wbuf, outfile);
      putc('\t', outfile);
      if (marker_allele_ptrs[2 * marker_uidx] != missing_geno_ptr) {
        fputs(marker_allele_ptrs[2 * marker_uidx], outfile);
      } else {
        fputs(output_missing_geno_ptr, outfile);
      }
      putc('\t', outfile);
      if (marker_allele_ptrs[2 * marker_uidx] != missing_geno_ptr) {
        fputs(marker_allele_ptrs[2 * marker_uidx + 1], outfile);
      } else {
	fputs(output_missing_geno_ptr, outfile);
      }
      if (putc_checked('\n', outfile)) {
	goto sort_and_write_bim_ret_WRITE_FAIL;
      }
      map_reverse[marker_uidx] = marker_idx;
    }
  }
  while (0) {
  sort_and_write_bim_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  sort_and_write_bim_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  sort_and_write_bim_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
  }
  fclose_cond(outfile);
  return retval;
}

int32_t load_sort_and_write_map(uint32_t** map_reverse_ptr, FILE* mapfile, uint32_t map_cols, char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, uintptr_t max_marker_id_len, int32_t compact_map_reverse, Chrom_info* chrom_info_ptr) {
  // no need to support allow_extra_chroms on .map files
  FILE* map_outfile = NULL;
  int64_t* ll_buf = NULL;
  int32_t retval = 0;
  char zstr[2];
  char wbuf[16];
  char* marker_ids;
  double* marker_cms;
  uint32_t* pos_buf;
  uint32_t* unpack_map;
  uintptr_t marker_uidx;
  uintptr_t marker_idx;
  uintptr_t marker_idx2;
  char* bufptr;
  uint32_t uii;
  uint32_t ujj;
  uint32_t chrom_start[MAX_POSSIBLE_CHROM + 1];
  uint32_t chrom_id[MAX_POSSIBLE_CHROM];
  uint32_t cur_chrom;
  uint32_t chrom_ct;
  // See sort_and_write_bim() for discussion.
  if (wkspace_alloc_ui_checked(map_reverse_ptr, (compact_map_reverse? marker_ct : unfiltered_marker_ct) * sizeof(int32_t)) ||
      wkspace_alloc_ll_checked(&ll_buf, marker_ct * sizeof(int64_t)) ||
      wkspace_alloc_c_checked(&marker_ids, marker_ct * max_marker_id_len) ||
      wkspace_alloc_d_checked(&marker_cms, marker_ct * sizeof(double)) ||
      wkspace_alloc_ui_checked(&pos_buf, marker_ct * sizeof(int32_t)) ||
      wkspace_alloc_ui_checked(&unpack_map, marker_ct * sizeof(int32_t))) {
    goto load_sort_and_write_map_ret_NOMEM;
  }
  rewind(mapfile);
  marker_idx = 0;
  for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
    if (get_next_noncomment(mapfile, &bufptr)) {
      goto load_sort_and_write_map_ret_READ_FAIL;
    }
    if (IS_SET(marker_exclude, marker_uidx)) {
      continue;
    }
    ll_buf[marker_idx] = (((int64_t)get_chrom_code(chrom_info_ptr, bufptr)) << 32) + marker_idx;
    bufptr = next_item(bufptr);
    uii = strlen_se(bufptr);
    memcpyx(&(marker_ids[marker_idx * max_marker_id_len]), bufptr, uii, 0);
    bufptr = next_item(bufptr);
    if (map_cols == 4) {
      marker_cms[marker_idx] = atof(bufptr);
      bufptr = next_item(bufptr);
    } else {
      marker_cms[marker_idx] = 0.0;
    }
    unpack_map[marker_idx] = marker_uidx;
    pos_buf[marker_idx++] = atoi(bufptr);
  }
  sort_marker_chrom_pos(ll_buf, marker_ct, pos_buf, chrom_start, chrom_id, &chrom_ct);

  strcpy(outname_end, ".map.tmp");
  if (fopen_checked(&map_outfile, outname, "w")) {
    goto load_sort_and_write_map_ret_OPEN_FAIL;
  }

  marker_idx = 0;
  *zstr = '0';
  zstr[1] = '\0';
  for (uii = 0; uii < chrom_ct; uii++) {
    cur_chrom = chrom_id[uii];
    ujj = chrom_start[uii + 1];
    for (; marker_idx < ujj; marker_idx++) {
      marker_idx2 = (uint32_t)ll_buf[marker_idx];
      marker_uidx = unpack_map[marker_idx2];
      bufptr = uint32_writex(wbuf, cur_chrom, '\t');
      fwrite(wbuf, 1, bufptr - wbuf, map_outfile);
      fputs(&(marker_ids[marker_idx2 * max_marker_id_len]), map_outfile);
      wbuf[0] = '\t';
      bufptr = double_g_writex(&(wbuf[1]), marker_cms[marker_idx2], '\t');
      fwrite(wbuf, 1, bufptr - wbuf, map_outfile);
      bufptr = uint32_write(wbuf, (uint32_t)(ll_buf[marker_idx] >> 32));
      fwrite(wbuf, 1, bufptr - wbuf, map_outfile);
      if (putc_checked('\n', map_outfile)) {
	goto load_sort_and_write_map_ret_WRITE_FAIL;
      }
      (*map_reverse_ptr)[compact_map_reverse? marker_idx2 : marker_uidx] = marker_idx;
    }
  }
  if (fclose_null(&map_outfile)) {
    goto load_sort_and_write_map_ret_WRITE_FAIL;
  }
  while (0) {
  load_sort_and_write_map_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  load_sort_and_write_map_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  load_sort_and_write_map_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  load_sort_and_write_map_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
  }
  if (ll_buf) {
    wkspace_reset((unsigned char*)ll_buf);
  }
  return retval;
}

int32_t make_bed_one_marker(FILE* bedfile, uintptr_t* loadbuf, uint32_t unfiltered_indiv_ct, uintptr_t unfiltered_indiv_ct4, uintptr_t* indiv_exclude, uint32_t indiv_ct, uint32_t* indiv_sort_map, uint32_t is_reverse, uintptr_t* writebuf) {
  uintptr_t* writeptr = writebuf;
  uintptr_t cur_word = 0;
  uint32_t indiv_uidx = 0;
  uint32_t ii_rem = 0;
  uint32_t indiv_idx = 0;
  uint32_t indiv_uidx2;
  if (indiv_sort_map) {
    if (fread(loadbuf, 1, unfiltered_indiv_ct4, bedfile) < unfiltered_indiv_ct4) {
      return RET_READ_FAIL;
    }
    for (; indiv_idx < indiv_ct; indiv_idx++) {
      do {
	indiv_uidx2 = indiv_sort_map[indiv_uidx++];
      } while (IS_SET(indiv_exclude, indiv_uidx2));
      cur_word |= (((loadbuf[indiv_uidx2 / BITCT2] >> ((indiv_uidx2 % BITCT2) * 2)) & 3) << (ii_rem * 2));
      if (++ii_rem == BITCT2) {
	*writeptr++ = cur_word;
	cur_word = 0;
	ii_rem = 0;
      }
    }
    if (ii_rem) {
      *writeptr = cur_word;
    }
    if (is_reverse) {
      reverse_loadbuf((unsigned char*)writebuf, indiv_ct);
    }
  } else {
    if (load_and_collapse(bedfile, loadbuf, unfiltered_indiv_ct, writeptr, indiv_ct, indiv_exclude, is_reverse)) {
      return RET_READ_FAIL;
    }
  }
  return 0;
}

int32_t make_bed(FILE* bedfile, uintptr_t bed_offset, char* bimname, uint32_t map_cols, char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, double* marker_cms, uint32_t* marker_pos, char** marker_allele_ptrs, uintptr_t* marker_reverse, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_ct, char* person_ids, uintptr_t max_person_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* sex_nm, uintptr_t* sex_male, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, char* output_missing_pheno, uint32_t map_is_unsorted, uint32_t* indiv_sort_map, uint64_t misc_flags, Two_col_params* update_chr, char* flip_subset_fname, uint32_t hh_exists, Chrom_info* chrom_info_ptr) {
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t unfiltered_indiv_ctl2 = (unfiltered_indiv_ct + BITCT2 - 1) / BITCT2;
  uintptr_t indiv_ct4 = (indiv_ct + 3) / 4;
  uintptr_t indiv_ctv2 = 2 * ((indiv_ct + (BITCT - 1)) / BITCT);
  uintptr_t marker_uidx = 0;
  uintptr_t marker_idx = 0;
  FILE* bedoutfile = NULL;
  int64_t* ll_buf = NULL;
  uintptr_t* indiv_include2 = NULL;
  uintptr_t* indiv_male_include2 = NULL;
  uint32_t set_hh_missing = (misc_flags / MISC_SET_HH_MISSING) & 1;
  uint32_t allow_extra_chroms = (misc_flags / MISC_ALLOW_EXTRA_CHROMS) & 1;
  uint32_t zero_extra_chroms = (misc_flags / MISC_ZERO_EXTRA_CHROMS) & 1;
  uint32_t chrom_end = 0;
  uint32_t chrom_fo_idx = 0xffffffffU; // deliberate overflow
  int32_t retval = 0;
  uintptr_t* loadbuf;
  uintptr_t* writebuf;
  uint32_t* map_reverse;
  uint32_t is_haploid;
  uint32_t is_x;
  uint32_t is_y;
  uint32_t pct;
  uint32_t loop_end;
  if (flip_subset_fname) {
    logprint("Error: --flip-subset is not implemented yet.\n");
    retval = RET_CALC_NOT_YET_SUPPORTED;
    goto make_bed_ret_1;
  }
  if (wkspace_alloc_ul_checked(&loadbuf, unfiltered_indiv_ctl2 * sizeof(intptr_t))) {
    goto make_bed_ret_NOMEM;
  }
  loadbuf[unfiltered_indiv_ctl2 - 1] = 0;

  memcpy(outname_end, ".bed", 5);
  if (fopen_checked(&bedoutfile, outname, "wb")) {
    goto make_bed_ret_OPEN_FAIL;
  }

  if (fwrite_checked("l\x1b\x01", 3, bedoutfile)) {
    goto make_bed_ret_WRITE_FAIL;
  }
  fflush(stdout);
  if (fseeko(bedfile, bed_offset, SEEK_SET)) {
    goto make_bed_ret_READ_FAIL;
  }
  if (map_is_unsorted || update_chr) {
    if (set_hh_missing) {
      logprint("\nError: --set-hh-missing cannot currently be used on an unsorted .bim file.  Use\n--make-bed without --set-hh-missing to sort by position first; then run\n--make-bed + --set-hh-missing on the new fileset.\n");
      retval = RET_CALC_NOT_YET_SUPPORTED;
      goto make_bed_ret_1;
    }
    if (wkspace_alloc_ui_checked(&map_reverse, unfiltered_marker_ct * sizeof(int32_t)) ||
        wkspace_alloc_ll_checked(&ll_buf, marker_ct * sizeof(int64_t))) {
      goto make_bed_ret_NOMEM;
    }
    if ((map_is_unsorted & UNSORTED_SPLIT_CHROM) || update_chr) {
      if (map_is_unsorted & UNSORTED_SPLIT_CHROM) {
        retval = load_bim_split_chrom(bimname, marker_exclude, marker_ct, chrom_info_ptr, ll_buf);
        if (retval) {
	  goto make_bed_ret_1;
        }
      } else {
	fill_ll_buf(marker_exclude, marker_ct, chrom_info_ptr, ll_buf);
      }
      if (update_chr) {
        retval = update_marker_chroms(update_chr, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, allow_extra_chroms, chrom_info_ptr, ll_buf);
	if (retval) {
	  goto make_bed_ret_1;
	}
      }
    } else {
      fill_ll_buf(marker_exclude, marker_ct, chrom_info_ptr, ll_buf);
    }
    memcpy(outname_end, ".bim", 5);
    retval = sort_and_write_bim(map_reverse, map_cols, outname, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, marker_cms, marker_pos, marker_allele_ptrs, ll_buf, zero_extra_chroms, chrom_info_ptr);
    if (retval) {
      goto make_bed_ret_1;
    }
    wkspace_reset((unsigned char*)ll_buf);

    if (wkspace_alloc_ul_checked(&writebuf, marker_ct * indiv_ctv2)) {
      // todo: implement multipass.  should now be a minimal change
      logprint("\nError: Insufficient memory for current --make-bed implementation.  Try raising\nthe --memory value for now.\n");
      retval = RET_CALC_NOT_YET_SUPPORTED;
      goto make_bed_ret_1;
    }
    sprintf(logbuf, "--make-bed to %s + .bim + .fam... ", outname);
    logprintb();
    fputs("0%", stdout);
    for (pct = 1; pct <= 100; pct++) {
      loop_end = (pct * ((uint64_t)marker_ct)) / 100;
      for (; marker_idx < loop_end; marker_uidx++, marker_idx++) {
	if (IS_SET(marker_exclude, marker_uidx)) {
	  marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
	  if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	    goto make_bed_ret_READ_FAIL;
	  }
	}
        retval = make_bed_one_marker(bedfile, loadbuf, unfiltered_indiv_ct, unfiltered_indiv_ct4, indiv_exclude, indiv_ct, indiv_sort_map, IS_SET(marker_reverse, marker_uidx), &(writebuf[indiv_ctv2 * map_reverse[marker_uidx]]));
	if (retval) {
	  goto make_bed_ret_1;
	}
      }
      if (pct < 100) {
	if (pct > 10) {
	  putchar('\b');
	}
	printf("\b\b%u%%", pct);
	fflush(stdout);
      }
    }
    for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
      if (fwrite_checked(writebuf, indiv_ct4, bedoutfile)) {
	goto make_bed_ret_WRITE_FAIL;
      }
      writebuf = &(writebuf[indiv_ctv2]);
    }
  } else {
    if (!hh_exists) {
      set_hh_missing = 0;
    } else if (!set_hh_missing) {
      hh_exists = 0;
    }
    if (alloc_collapsed_haploid_filters(unfiltered_indiv_ct, indiv_ct, hh_exists, 0, indiv_exclude, sex_male, &indiv_include2, &indiv_male_include2)) {
      goto make_bed_ret_NOMEM;
    }

    if (wkspace_alloc_ul_checked(&writebuf, indiv_ctv2)) {
      goto make_bed_ret_NOMEM;
    }
    sprintf(logbuf, "--make-bed to %s + .bim + .fam... ", outname);
    logprintb();
    fputs("0%", stdout);
    marker_uidx = 0;
    for (pct = 1; pct <= 100; pct++) {
      loop_end = (pct * ((uint64_t)marker_ct)) / 100;
      for (; marker_idx < loop_end; marker_uidx++, marker_idx++) {
        if (IS_SET(marker_exclude, marker_uidx)) {
          marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
	  if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
            goto make_bed_ret_READ_FAIL;
	  }
	}
	if (marker_uidx >= chrom_end) {
          chrom_fo_idx++;
	  refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_haploid);
	}
        retval = make_bed_one_marker(bedfile, loadbuf, unfiltered_indiv_ct, unfiltered_indiv_ct4, indiv_exclude, indiv_ct, indiv_sort_map, IS_SET(marker_reverse, marker_uidx), writebuf);
	if (retval) {
	  goto make_bed_ret_1;
	}
	if (is_haploid && set_hh_missing) {
	  haploid_fix(hh_exists, indiv_include2, indiv_male_include2, indiv_ct, is_x, is_y, (unsigned char*)writebuf);
	}
	if (fwrite_checked(writebuf, indiv_ct4, bedoutfile)) {
	  goto make_bed_ret_WRITE_FAIL;
	}
      }
      if (pct < 100) {
	if (pct > 10) {
	  putchar('\b');
	}
	printf("\b\b%u%%", pct);
	fflush(stdout);
      }
    }
  }
  if (fclose_null(&bedoutfile)) {
    goto make_bed_ret_WRITE_FAIL;
  }

  memcpy(outname_end, ".fam", 5);
  retval = write_fam(outname, unfiltered_indiv_ct, indiv_exclude, indiv_ct, person_ids, max_person_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, sex_nm, sex_male, pheno_nm, pheno_c, pheno_d, output_missing_pheno, ' ', indiv_sort_map);
  if (retval) {
    goto make_bed_ret_1;
  }

  if (!map_is_unsorted) {
    memcpy(outname_end, ".bim", 5);
    retval = write_map_or_bim(outname, marker_exclude, marker_ct, marker_ids, max_marker_id_len, marker_cms, marker_pos, marker_allele_ptrs, '\t', zero_extra_chroms, chrom_info_ptr);
    if (retval) {
      goto make_bed_ret_1;
    }
  }

  fputs("\b\b\b", stdout);
  logprint("done.\n");
  while (0) {
  make_bed_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  make_bed_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  make_bed_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  make_bed_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
 make_bed_ret_1:
  fclose_cond(bedoutfile);
  wkspace_reset(wkspace_mark);
  return retval;
}

const char errstr_fam_format[] = "Error: Improperly formatted .fam file.\n";

int32_t load_fam(FILE* famfile, uint32_t buflen, uint32_t fam_cols, uint32_t tmp_fam_col_6, int32_t missing_pheno, uint32_t missing_pheno_len, uint32_t affection_01, uintptr_t* unfiltered_indiv_ct_ptr, char** person_ids_ptr, uintptr_t* max_person_id_len_ptr, char** paternal_ids_ptr, uintptr_t* max_paternal_id_len_ptr, char** maternal_ids_ptr, uintptr_t* max_maternal_id_len_ptr, uintptr_t** sex_nm_ptr, uintptr_t** sex_male_ptr, uint32_t* affection_ptr, uintptr_t** pheno_nm_ptr, uintptr_t** pheno_c_ptr, double** pheno_d_ptr, uintptr_t** founder_info_ptr, uintptr_t** indiv_exclude_ptr) {
  uintptr_t unfiltered_indiv_ct = 0;
  uintptr_t max_person_id_len = *max_person_id_len_ptr;
  uintptr_t max_paternal_id_len = *max_paternal_id_len_ptr;
  uintptr_t max_maternal_id_len = *max_maternal_id_len_ptr;
  uint32_t affection = 1;
  uint64_t last_tell = 0;
  uint32_t new_buflen = 0;
  unsigned char* wkspace_mark = wkspace_base;
  double missing_phenod = (double)missing_pheno;
  char case_char = affection_01? '1' : '2';
  char* paternal_ids = NULL;
  char* maternal_ids = NULL;
  uintptr_t* pheno_c = NULL;
  double* pheno_d = NULL;
  uintptr_t* sex_nm;
  uintptr_t* sex_male;
  uintptr_t* pheno_nm;
  uintptr_t* founder_info;
  char* bufptr0;
  char* bufptr;
  char* linebuf;
  char* person_ids;
  char cc;
  uint64_t* line_locs;
  uint64_t* tmp_ullp;
  uint32_t max_people;
  uintptr_t tmp_len;
  uintptr_t tmp_len2;
  uintptr_t indiv_uidx;
  uintptr_t unfiltered_indiv_ctl;
  double dxx;
  int32_t ii;
  char* fgets_return;

  if (wkspace_alloc_c_checked(&linebuf, buflen)) {
    return RET_NOMEM;
  }
  linebuf[buflen - 1] = ' ';
  line_locs = (uint64_t*)wkspace_base;
  max_people = wkspace_left / sizeof(int64_t);
  // ----- .fam read, first pass -----
  // count number of people, determine maximum person/father/mother ID lengths,
  // affection status, verify all floating point phenotype values are valid
  while (fgets(linebuf, buflen, famfile) != NULL) {
    bufptr0 = skip_initial_spaces(linebuf);
    if (!is_eoln_kns(*bufptr0)) {
      if (fam_cols & FAM_COL_1) {
	bufptr = next_item(bufptr0);
      } else {
	bufptr = bufptr0;
      }
      tmp_len = strlen_se(bufptr0) + strlen_se(bufptr) + 2;
      if (tmp_len > max_person_id_len) {
	max_person_id_len = tmp_len;
      }
      if (fam_cols & FAM_COL_34) {
	bufptr = next_item(bufptr);
	tmp_len = strlen_se(bufptr) + 1;
	if (tmp_len > max_paternal_id_len) {
	  max_paternal_id_len = tmp_len;
	}
	bufptr = next_item(bufptr);
	tmp_len = strlen_se(bufptr) + 1;
	if (tmp_len > max_maternal_id_len) {
	  max_maternal_id_len = tmp_len;
	}
      }
      if (fam_cols & FAM_COL_5) {
	bufptr = next_item(bufptr);
      }
      if (tmp_fam_col_6) {
	bufptr = next_item(bufptr);
	if (no_more_items_kns(bufptr)) {
	  logprint(errstr_fam_format);
	  return RET_INVALID_FORMAT;
	}
	if (affection) {
	  affection = eval_affection(bufptr, missing_pheno, missing_pheno_len, affection_01);
	}
      }
      if (unfiltered_indiv_ct == max_people) {
	return RET_NOMEM;
      }
      line_locs[unfiltered_indiv_ct++] = last_tell;
    }
    if (!linebuf[buflen - 1]) {
      // determine extended read buffer length needed to handle unexpectedly
      // long line
      linebuf[buflen - 1] = ' ';
      if (linebuf[buflen - 2] != '\n') {
	tmp_len = 0;
	do {
	  tmp_len += buflen - 1;
	  linebuf[buflen - 1] = ' ';
	  fgets_return = fgets(linebuf, buflen, famfile);
	} while (fgets_return && (!linebuf[buflen - 1]) && (linebuf[buflen - 2] != '\n'));
	tmp_len += strlen(linebuf) + 1;
	linebuf[buflen - 1] = ' ';
	if (tmp_len > new_buflen) {
	  new_buflen = tmp_len;
	}
      }
    }
    last_tell = ftello(famfile);
  }
  if (ferror(famfile)) {
    return RET_READ_FAIL;
  }
  if (!unfiltered_indiv_ct) {
    logprint("Error: Nobody in .fam file.\n");
    return RET_INVALID_FORMAT;
  }
  wkspace_reset(wkspace_mark);
  unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  // could make paternal_ids/maternal_ids conditional, but memory footprint is
  // typically negligible
  if (wkspace_alloc_c_checked(person_ids_ptr, unfiltered_indiv_ct * max_person_id_len) ||
      wkspace_alloc_c_checked(paternal_ids_ptr, unfiltered_indiv_ct * max_paternal_id_len) ||
      wkspace_alloc_c_checked(maternal_ids_ptr, unfiltered_indiv_ct * max_maternal_id_len) ||
      wkspace_alloc_ul_checked(sex_nm_ptr, unfiltered_indiv_ctl * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(sex_male_ptr, unfiltered_indiv_ctl * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(founder_info_ptr, unfiltered_indiv_ctl * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(indiv_exclude_ptr, unfiltered_indiv_ctl * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(pheno_nm_ptr, unfiltered_indiv_ctl * sizeof(intptr_t))) {
    return RET_NOMEM;
  }
  // currently NOT safe to initialize any of the above until line_locs has been
  // copied to tmp_ullp

  if (tmp_fam_col_6) {
    if (affection) {
      pheno_c = (uintptr_t*)malloc(unfiltered_indiv_ctl * sizeof(intptr_t));
      if (!pheno_c) {
	return RET_NOMEM;
      }
      fill_ulong_zero(pheno_c, unfiltered_indiv_ctl);
      *pheno_c_ptr = pheno_c;
    } else {
      pheno_d = (double*)malloc(unfiltered_indiv_ct * sizeof(double));
      if (!pheno_d) {
	return RET_NOMEM;
      }
      fill_double_zero(pheno_d, unfiltered_indiv_ct);
      *pheno_d_ptr = pheno_d;
    }
  }
  wkspace_mark = wkspace_base;
  if (new_buflen) {
    buflen = new_buflen;
  }
  if (wkspace_alloc_c_checked(&linebuf, buflen) ||
      wkspace_alloc_ull_checked(&tmp_ullp, unfiltered_indiv_ct * sizeof(int64_t))) {
    return RET_NOMEM;
  }
  for (ii = unfiltered_indiv_ct - 1; ii >= 0; ii--) {
    tmp_ullp[(uint32_t)ii] = line_locs[(uint32_t)ii];
  }
  line_locs = tmp_ullp;
  person_ids = *person_ids_ptr;
  paternal_ids = *paternal_ids_ptr;
  maternal_ids = *maternal_ids_ptr;
  if (!(fam_cols & FAM_COL_34)) {
    for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ct; indiv_uidx++) {
      memcpy(&(paternal_ids[indiv_uidx * max_paternal_id_len]), "0", 2);
    }
    for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ct; indiv_uidx++) {
      memcpy(&(maternal_ids[indiv_uidx * max_maternal_id_len]), "0", 2);
    }
  }
  sex_nm = *sex_nm_ptr;
  sex_male = *sex_male_ptr;
  pheno_nm = *pheno_nm_ptr;
  founder_info = *founder_info_ptr;
  if (fam_cols & FAM_COL_34) {
    fill_ulong_zero(founder_info, unfiltered_indiv_ctl);
  } else {
    fill_all_bits(founder_info, unfiltered_indiv_ct);
  }
  fill_ulong_zero(sex_nm, unfiltered_indiv_ctl);
  fill_ulong_zero(sex_male, unfiltered_indiv_ctl);
  fill_ulong_zero(*indiv_exclude_ptr, unfiltered_indiv_ctl);
  fill_ulong_zero(pheno_nm, unfiltered_indiv_ctl);

  // ----- .fam read, second pass -----
  for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ct; indiv_uidx++) {
    if (fseeko(famfile, line_locs[indiv_uidx], SEEK_SET)) {
      return RET_READ_FAIL;
    }
    if (fgets(linebuf, buflen, famfile) == NULL) {
      return RET_READ_FAIL;
    }
    if (fam_cols & FAM_COL_1) {
      bufptr = next_item(linebuf);
    } else {
      bufptr = linebuf;
    }
    tmp_len = strlen_se(linebuf);
    memcpyx(memcpyax(&(person_ids[indiv_uidx * max_person_id_len]), linebuf, tmp_len, '\t'), bufptr, strlen_se(bufptr), '\0');
    if (fam_cols & FAM_COL_34) {
      bufptr = next_item(bufptr);
      cc = *bufptr;
      tmp_len = strlen_se(bufptr);
      memcpyx(&(paternal_ids[indiv_uidx * max_paternal_id_len]), bufptr, tmp_len, '\0');
      bufptr = next_item(bufptr);
      tmp_len2 = strlen_se(bufptr);
      memcpyx(&(maternal_ids[indiv_uidx * max_maternal_id_len]), bufptr, tmp_len2, '\0');
      if ((tmp_len == 1) && (tmp_len2 == 1) && (cc == '0') && (*bufptr == '0')) {
	SET_BIT(founder_info, indiv_uidx);
      }
    }
    if (fam_cols & FAM_COL_5) {
      bufptr = next_item(bufptr);
      if (strlen_se(bufptr) == 1) {
	if (*bufptr == '1') {
	  SET_BIT(sex_nm, indiv_uidx);
	  SET_BIT(sex_male, indiv_uidx);
	} else if (*bufptr == '2') {
	  SET_BIT(sex_nm, indiv_uidx);
	}
      }
    }
    if (tmp_fam_col_6) {
      bufptr = next_item(bufptr);
      if (affection) {
	if (!is_missing_pheno(bufptr, missing_pheno, missing_pheno_len, affection_01)) {
	  SET_BIT(pheno_nm, indiv_uidx);
	  if (*bufptr == case_char) {
	    SET_BIT(pheno_c, indiv_uidx);
	  }
	}
      } else {
	if ((!scan_double(bufptr, &dxx)) && (dxx != missing_phenod)) {
	  pheno_d[indiv_uidx] = dxx;
	  SET_BIT(pheno_nm, indiv_uidx);
	}
      }
    }
  }

  *unfiltered_indiv_ct_ptr = unfiltered_indiv_ct;
  *max_person_id_len_ptr = max_person_id_len;
  *max_paternal_id_len_ptr = max_paternal_id_len;
  *max_maternal_id_len_ptr = max_maternal_id_len;
  *affection_ptr = affection;
  wkspace_reset(wkspace_mark);
  return 0;
}

// side effect: initializes tbuf to first nonempty line of .map/.bim
int32_t check_cm_col(FILE* bimfile, char* tbuf, uint32_t is_binary, uint32_t bufsize, uint32_t* gd_col_ptr) {
  char* bufptr;
  while (fgets(tbuf, bufsize, bimfile)) {
    bufptr = skip_initial_spaces(tbuf);
    if (is_eoln_or_comment(*bufptr)) {
      continue;
    }
    bufptr = next_item_mult(bufptr, 2 + 2 * is_binary);
    if (no_more_items_kns(bufptr)) {
      return -1;
    }
    if (no_more_items_kns(next_item(bufptr))) {
      *gd_col_ptr = 0;
    }
    *gd_col_ptr = 1;
    return 0;
  }
  return -1;
}

int32_t incr_text_allele0(char cc, char* marker_alleles, uint32_t* marker_allele_cts) {
  uint32_t uii;
  for (uii = 0; uii < 4; uii++) {
    if (marker_alleles[uii] == '\0') {
      marker_alleles[uii] = cc;
      marker_allele_cts[uii] = 1;
      return 0;
    } else if (marker_alleles[uii] == cc) {
      marker_allele_cts[uii] += 1;
      return 0;
    }
  }
  return -1;
}

typedef struct ll_str_fixed_struct {
  struct ll_str_struct* next;
#ifdef __LP64__
  char ss[8];
#else
  char ss[12];
#endif
} Ll_str_fixed;

int32_t incr_text_allele_str(uintptr_t* topsize_ptr, char* allele_name, uint32_t an_len, Ll_str* allele_list_start, uint32_t* marker_allele_cts) {
  // Start with preallocated array of 16-byte Ll_strs.
  // Ll_str.ss is a null-terminated sequence of ordered, tab-delimited allele
  // names.  If the starting 8 (or 12 bytes, on 32-bit systems) is adequate,
  // Ll_str.next is NULL.  Otherwise, Ll_str.ss stores the first few (or
  // possibly 0, if the very first allele name is too long) allele names, and
  // Ll_str.next is a pointer to a linked list entry storing the next 1+ allele
  // names.  Worst case, the linked list is of length 4 (beyond that we error
  // out).  Most of the time, it'll be length 0.
  // This type of function will become more important when it's time to parse
  // .vcf files, etc.
  uint32_t allele_num = 0;
  char* cur_allele_name_start = allele_list_start->ss;
  Ll_str* llptr;
  uint32_t slen;
  uintptr_t chars_left;
  if (!(*cur_allele_name_start)) {
    if (!(allele_list_start->next)) {
      if (an_len >= (16 - sizeof(intptr_t))) {
	llptr = top_alloc_llstr(topsize_ptr, an_len + 1);
	if (!llptr) {
	  return RET_NOMEM;
	}
        allele_list_start->next = llptr;
	llptr->next = NULL;
	cur_allele_name_start = llptr->ss;
      }
      memcpyx(cur_allele_name_start, allele_name, an_len, '\0');
      marker_allele_cts[0] = 1;
      return 0;
    }
    allele_list_start = allele_list_start->next;
    cur_allele_name_start = allele_list_start->ss;
  }
  while (allele_num < 4) {
    slen = strlen_se(cur_allele_name_start);
    if ((slen == an_len) && (!memcmp(cur_allele_name_start, allele_name, an_len))) {
      marker_allele_cts[allele_num] += 1;
      return 0;
    } else {
      allele_num++;
      if (cur_allele_name_start[slen] == '\t') {
	cur_allele_name_start = &(cur_allele_name_start[slen + 1]);
      } else if (allele_list_start->next) {
	allele_list_start = allele_list_start->next;
	cur_allele_name_start = allele_list_start->ss;
      } else {
	chars_left = ((0x7ffffff0 - sizeof(intptr_t)) - ((uintptr_t)(&(cur_allele_name_start[slen + 1]) - allele_list_start->ss))) & 15;
	if (chars_left > an_len) {
	  cur_allele_name_start[slen] = '\t';
	  memcpyx(&(cur_allele_name_start[slen + 1]), allele_name, an_len, '\0');
	} else {
	  llptr = top_alloc_llstr(topsize_ptr, an_len + 1);
	  if (!llptr) {
	    return RET_NOMEM;
	  }
	  allele_list_start->next = llptr;
	  llptr->next = NULL;
	  cur_allele_name_start = llptr->ss;
	  memcpyx(cur_allele_name_start, allele_name, an_len, '\0');
	}
	marker_allele_cts[allele_num] = 1;
	return 0;
      }
    }
  }
  return RET_INVALID_FORMAT;
}

char* get_llstr(Ll_str* llptr, uint32_t allele_idx) {
  char* cptr = llptr->ss;
  if (*cptr == '\0') {
    llptr = llptr->next;
    if (!llptr) {
      return NULL;
    }
    cptr = llptr->ss;
  }
  while (allele_idx) {
    cptr = item_endnn(cptr);
    allele_idx--;
    if (*cptr) {
      cptr++;
    } else {
      llptr = llptr->next;
      if (!llptr) {
	return NULL;
      }
      cptr = llptr->ss;
    }
  }
  return cptr;
}

int32_t ped_to_bed_multichar_allele(FILE** pedfile_ptr, FILE** outfile_ptr, char* outname, char* outname_end, FILE** mapfile_ptr, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_alleles_f, uint32_t map_is_unsorted, uint32_t fam_cols, uint32_t ped_col_skip, uint32_t gd_col, uint32_t* map_reverse, int64_t ped_size) {
  // maintain allele counts and linked lists of observed alleles at FAR end of
  // wkspace.
  int32_t retval = 0;
  uintptr_t topsize = marker_ct * (4LU * sizeof(int32_t) + 16);
  uint32_t ped_buflen = 0;
  uintptr_t indiv_ct = 0;
  uint32_t pct = 1;
  int64_t ped_next_thresh = ped_size / 100;
  uint32_t last_pass = 0;
  int64_t* line_starts = NULL;
  char* missing_geno_ptr = (char*)g_missing_geno_ptr;
  char** marker_allele_ptrs = NULL;
  FILE* outfile;
  uint32_t pass_ct;
  uintptr_t indiv_ct4;
  char* loadbuf;
  char* col1_ptr;
  char* col2_ptr;
  char* bufptr;
  char* bufptr2;
  uintptr_t loadbuf_size;
  uintptr_t cur_slen;
  uintptr_t cur_slen_rdup;
  Ll_str_fixed* marker_alleles_tmp;
  uint32_t* marker_allele_cts;
  uint32_t markers_per_pass;
  uint32_t marker_start;
  uint32_t marker_end;
  uintptr_t marker_uidx;
  uintptr_t marker_idx;
  uint32_t loop_end;
  uintptr_t indiv_idx;
  uintptr_t ulii;
  uintptr_t uljj;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  uint32_t alen1;
  uint32_t alen2;
  char* aptr1;
  char* aptr2;
  unsigned char ucc;
  uint32_t ii_shift;
  unsigned char* writebuf;
  unsigned char* wbufptr;
  wkspace_reset((unsigned char*)marker_alleles_f);
  if ((wkspace_left / (4LU * sizeof(int32_t) + 16)) <= marker_ct) {
    goto ped_to_bed_multichar_allele_ret_NOMEM;
  }
  loadbuf = (char*)wkspace_base;
  marker_allele_cts = (uint32_t*)(&(wkspace_base[wkspace_left - marker_ct * 4LU * sizeof(int32_t)]));
  marker_alleles_tmp = (Ll_str_fixed*)(&(wkspace_base[wkspace_left - marker_ct * (4LU * sizeof(int32_t) + 16)]));
  memset(marker_alleles_tmp, 0, marker_ct * (4LU * sizeof(int32_t) + 16));

  if (fclose_null(outfile_ptr)) {
    goto ped_to_bed_multichar_allele_ret_WRITE_FAIL;
  }
  memcpy(outname_end, ".fam", 5);
  if (fopen_checked(outfile_ptr, outname, "w")) {
    goto ped_to_bed_multichar_allele_ret_OPEN_FAIL;
  }
  outfile = *outfile_ptr;
  rewind(*pedfile_ptr);
  fputs("Rescanning .ped file... 0%", stdout);
  fflush(stdout);
  while (1) {
    loadbuf_size = wkspace_left - topsize;
    if (loadbuf_size > MAXLINEBUFLEN) {
      loadbuf_size = MAXLINEBUFLEN;
    }
    loadbuf[loadbuf_size - 1] = ' ';
  ped_to_bed_multichar_allele_loop_1_start:
    if (!fgets(loadbuf, loadbuf_size, *pedfile_ptr)) {
      break;
    }
    if (!loadbuf[loadbuf_size - 1]) {
      putchar('\n');
      if (loadbuf_size == MAXLINEBUFLEN) {
        sprintf(logbuf, "Error: " PROG_NAME_CAPS " does not currently support .ped lines longer than %u\nbytes.\n", MAXLINEBUFLEN);
	goto ped_to_bed_multichar_allele_ret_INVALID_FORMAT_2;
      } else {
        goto ped_to_bed_multichar_allele_ret_NOMEM;
      }
    }
    cur_slen = strlen(loadbuf);
    ulii = cur_slen + 1;
    if (ulii > ped_buflen) {
      ped_buflen = ulii;
    }
    col1_ptr = skip_initial_spaces(loadbuf);
    if (is_eoln_or_comment(*col1_ptr)) {
      goto ped_to_bed_multichar_allele_loop_1_start;
    }
    // check for top-of-stack allocations colliding with load buffer
    cur_slen_rdup = (cur_slen + CACHELINE) & (~(CACHELINE - 1));
    if (fam_cols & FAM_COL_1) {
      col2_ptr = next_item(col1_ptr);
    } else {
      col2_ptr = col1_ptr;
    }
    bufptr = next_item_mult(col2_ptr, ped_col_skip - 1);
    if (no_more_items_kns(bufptr)) {
      putchar('\n');
      sprintf(logbuf, "Error: Missing token(s) in .ped line: %s\n", col1_ptr);
      goto ped_to_bed_multichar_allele_ret_INVALID_FORMAT_2;
    }
    if (fwrite_checked(col1_ptr, strlen_se(col1_ptr), outfile)) {
      goto ped_to_bed_multichar_allele_ret_WRITE_FAIL;
    }
    putc('\t', outfile);
    bufptr2 = write_item(col2_ptr, outfile);
    if (fam_cols & FAM_COL_34) {
      bufptr2 = write_item(bufptr2, outfile);
      bufptr2 = write_item(bufptr2, outfile);
    } else {
      fputs("0\t0\t", outfile);
    }
    if (fam_cols & FAM_COL_5) {
      bufptr2 = write_item_nt(bufptr2, outfile);
    } else {
      putc('0', outfile);
    }
    putc('\t', outfile);
    if (fam_cols & FAM_COL_6) {
      uii = strlen_se(bufptr2);
      fwrite(bufptr2, 1, uii, outfile);
      if (putc_checked('\n', outfile)) {
        goto ped_to_bed_multichar_allele_ret_WRITE_FAIL;
      }
    } else {
      if (fputs_checked("-9\n", outfile)) {
        goto ped_to_bed_multichar_allele_ret_WRITE_FAIL;
      }
    }
    wkspace_base += cur_slen_rdup;
    wkspace_left -= cur_slen_rdup;
    for (marker_uidx = 0, marker_idx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
      alen1 = strlen_se(bufptr);
      aptr1 = bufptr;
      bufptr = skip_initial_spaces(&(bufptr[alen1]));
      alen2 = strlen_se(bufptr);
      if (!alen2) {
	wkspace_base -= cur_slen_rdup;
	wkspace_left += cur_slen_rdup;
	goto ped_to_bed_multichar_allele_ret_INVALID_FORMAT_3;
      }
      aptr2 = bufptr;
      bufptr = skip_initial_spaces(&(bufptr[alen2]));
      if (IS_SET(marker_exclude, marker_uidx)) {
	continue;
      }
      if ((*aptr1 == '0') && (alen1 == 1)) {
	if ((alen2 != 1) || (*aptr2 != '0')) {
          goto ped_to_bed_multichar_allele_ret_INVALID_FORMAT_4;
	}
	marker_idx++;
	continue;
      } else if ((*aptr2 == '0') && (alen2 == 1)) {
	goto ped_to_bed_multichar_allele_ret_INVALID_FORMAT_4;
      }
      uii = map_is_unsorted? map_reverse[marker_idx] : marker_idx;
      retval = incr_text_allele_str(&topsize, aptr1, alen1, (Ll_str*)(&(marker_alleles_tmp[uii])), &(marker_allele_cts[4 * uii]));
      if (retval) {
	goto ped_to_bed_multichar_allele_ret_INVALID_FORMAT_6;
      }
      retval = incr_text_allele_str(&topsize, aptr2, alen2, (Ll_str*)(&(marker_alleles_tmp[uii])), &(marker_allele_cts[4 * uii]));
      if (retval) {
	goto ped_to_bed_multichar_allele_ret_INVALID_FORMAT_6;
      }
      marker_idx++;
    }
    wkspace_base -= cur_slen_rdup;
    wkspace_left += cur_slen_rdup;
    if (!is_eoln_kns(*bufptr)) {
      putchar('\n');
      sprintf(logbuf, "Error: Too many tokens in .ped line for indiv %" PRIuPTR ".\n", indiv_ct);
      goto ped_to_bed_multichar_allele_ret_INVALID_FORMAT_2;
    }
    indiv_ct++;
    if (ftello(*pedfile_ptr) >= ped_next_thresh) {
      uii = (ftello(*pedfile_ptr) * 100) / ped_size;
      if (pct >= 10) {
	putchar('\b');
      }
      printf("\b\b%u%%", uii);
      fflush(stdout);
      pct = uii;
    }
  }
  if (!feof(*pedfile_ptr)) {
    goto ped_to_bed_multichar_allele_ret_READ_FAIL;
  }
  putchar('\r');
  logprint(".ped scan complete (for binary autoconversion).\n");
  if (!indiv_ct) {
    sprintf(logbuf, "Error: No %s in .ped file.\n", g_species_plural);
    goto ped_to_bed_multichar_allele_ret_INVALID_FORMAT_2;
  }
  if (fclose_null(outfile_ptr)) {
    goto ped_to_bed_multichar_allele_ret_WRITE_FAIL;
  }
  if (marker_ct * 2 * sizeof(intptr_t) + topsize > wkspace_left) {
    goto ped_to_bed_multichar_allele_ret_NOMEM;
  }
  marker_allele_ptrs = (char**)wkspace_alloc(marker_ct * 2 * sizeof(intptr_t));
  memcpy(outname_end, ".bim", 5);
  if (fopen_checked(outfile_ptr, outname, "w")) {
    goto ped_to_bed_multichar_allele_ret_OPEN_FAIL;
  }
  outfile = *outfile_ptr;
  if (map_is_unsorted) {
    memcpy(outname_end, ".map.tmp", 9);
    if (fopen_checked(mapfile_ptr, outname, "r")) {
      goto ped_to_bed_multichar_allele_ret_OPEN_FAIL;
    }
  } else {
    rewind(*mapfile_ptr);
  }
  marker_uidx = 0;
  for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
    if (map_is_unsorted) {
      if (!fgets(tbuf, MAXLINELEN, *mapfile_ptr)) {
	goto ped_to_bed_multichar_allele_ret_READ_FAIL;
      }
    } else {
      if (get_next_noncomment_excl(*mapfile_ptr, &bufptr, marker_exclude, &marker_uidx)) {
	goto ped_to_bed_multichar_allele_ret_READ_FAIL;
      }
    }
    if (marker_allele_cts[4 * marker_idx + 2]) {
      uii = marker_allele_cts[4 * marker_idx + 3];
      if (map_is_unsorted) {
        sprintf(logbuf, "Warning: Variant %u (post-sort/filter) %sallelic; setting rarest missing.\n", map_reverse[marker_idx] + 1, (uii? "quad" : "tri"));
      } else {
        sprintf(logbuf, "Warning: Variant %" PRIuPTR " %sallelic; setting rarest alleles missing.\n", marker_idx + 1, (uii? "quad" : "tri"));
      }
      get_top_two(&(marker_allele_cts[4 * marker_idx]), uii? 4 : 3, &ulii, &uljj);
      uii = map_reverse[marker_idx];
    } else {
      ulii = (marker_allele_cts[4 * marker_idx] < marker_allele_cts[4 * marker_idx + 1])? 1 : 0;
      uljj = ulii ^ 1;
      uii = marker_idx;
    }

    aptr1 = get_llstr((Ll_str*)(&(marker_alleles_tmp[uii])), uljj);
    aptr2 = get_llstr((Ll_str*)(&(marker_alleles_tmp[uii])), ulii);
    if (!aptr1) {
      aptr1 = missing_geno_ptr;
    } else {
      alen1 = strlen_se(aptr1);
      if (aptr1[alen1] == '\t') {
	aptr1[alen1] = '\0';
      }
    }
    if (!aptr2) {
      aptr2 = missing_geno_ptr;
    } else {
      alen2 = strlen_se(aptr2);
      if (aptr2[alen2] == '\t') {
	aptr2[alen2] = '\0';
      }
    }
    marker_allele_ptrs[2 * marker_idx] = aptr1;
    marker_allele_ptrs[2 * marker_idx + 1] = aptr2;
    if (map_is_unsorted) {
      bufptr = (char*)memchr(tbuf, '\n', MAXLINELEN);
      if (fwrite_checked(tbuf, bufptr - tbuf, outfile)) {
        goto ped_to_bed_multichar_allele_ret_WRITE_FAIL;
      }
    } else {
      uii = strlen_se(bufptr);
      if (fwrite_checked(bufptr, uii, outfile)) {
        goto ped_to_bed_multichar_allele_ret_WRITE_FAIL;
      }
      putc('\t', outfile);
      bufptr = skip_initial_spaces(&(bufptr[uii + 1]));
      bufptr = write_item(bufptr, outfile);
      if (gd_col) {
	bufptr = write_item_nt(bufptr, outfile);
      } else {
	putc('0', outfile);
      }
      putc('\t', outfile);
      uii = strlen_se(bufptr);
      fwrite(bufptr, 1, uii, outfile);
    }
    putc('\t', outfile);
    fputs(aptr1, outfile);
    putc('\t', outfile);
    fputs(aptr2, outfile);
    if (putc_checked('\n', outfile)) {
      goto ped_to_bed_multichar_allele_ret_WRITE_FAIL;
    }
    marker_uidx++;
  }
  indiv_ct4 = (indiv_ct + 3) / 4;
  fclose_null(mapfile_ptr);
  if (map_is_unsorted) {
    unlink(outname);
  }
  fclose_null(outfile_ptr);
  if (wkspace_alloc_c_checked(&loadbuf, ped_buflen)) {
    goto ped_to_bed_multichar_allele_ret_NOMEM;
  }
  if (wkspace_left >= marker_ct * indiv_ct4) {
    markers_per_pass = marker_ct;
    sprintf(logbuf, "Performing single-pass .bed write (%" PRIuPTR " variant%s, %" PRIuPTR " %s).\n", marker_ct, (marker_ct == 1)? "" : "s", indiv_ct, species_str(indiv_ct));
    pass_ct = 1;
  } else {
    if (!map_is_unsorted) {
      if (wkspace_alloc_ll_checked(&line_starts, indiv_ct * sizeof(int64_t))) {
	goto ped_to_bed_multichar_allele_ret_NOMEM;
      }
    }
    markers_per_pass = wkspace_left / indiv_ct4;
    if (!markers_per_pass) {
      goto ped_to_bed_multichar_allele_ret_NOMEM;
    }
    pass_ct = (marker_ct + markers_per_pass - 1) / markers_per_pass;
    sprintf(logbuf, "Performing %u-pass .bed write (%u/%" PRIuPTR " variant%s/pass, %" PRIuPTR " %s).\n", pass_ct, markers_per_pass, marker_ct, (markers_per_pass == 1)? "" : "s", indiv_ct, species_str(indiv_ct));
  }
  logprintb();
  writebuf = wkspace_base;
  memcpy(outname_end, ".bed", 5);
  if (fopen_checked(outfile_ptr, outname, "wb")) {
    goto ped_to_bed_multichar_allele_ret_OPEN_FAIL;
  }
  if (fwrite_checked("l\x1b\x01", 3, *outfile_ptr)) {
    goto ped_to_bed_multichar_allele_ret_WRITE_FAIL;
  }
  rewind(*pedfile_ptr);
  umm = 0;
  for (uii = 0; uii < pass_ct; uii++) {
    marker_start = uii * markers_per_pass;
    if (uii + 1 == pass_ct) {
      ujj = marker_ct - marker_start;
      last_pass = 1;
    } else {
      ujj = markers_per_pass;
    }
    memset(writebuf, 0, ujj * indiv_ct4);
    marker_end = marker_start + ujj;
    fputs("0%", stdout);
    indiv_idx = 0;
    // 94 instead of 100 due to big fwrite at the end
    for (pct = 1; pct <= 94; pct++) {
      loop_end = (((uint64_t)pct) * indiv_ct) / 94LLU;
      for (; indiv_idx < loop_end; indiv_idx++) {
	if ((!uii) || map_is_unsorted) {
	  do {
	    if (!last_pass) {
	      ped_next_thresh = ftello(*pedfile_ptr);
	    }
	    if (!fgets(loadbuf, ped_buflen, *pedfile_ptr)) {
	      goto ped_to_bed_multichar_allele_ret_READ_FAIL_2;
	    }
	    col1_ptr = skip_initial_spaces(loadbuf);
	  } while (is_eoln_or_comment(*col1_ptr));
	  bufptr = next_item_mult(col1_ptr, ped_col_skip);
	} else {
	  ped_next_thresh = line_starts[indiv_idx];
	  if (fseeko(*pedfile_ptr, line_starts[indiv_idx], SEEK_SET)) {
	    goto ped_to_bed_multichar_allele_ret_READ_FAIL_2;
	  }
	  if (!fgets(loadbuf, ped_buflen, *pedfile_ptr)) {
	    goto ped_to_bed_multichar_allele_ret_READ_FAIL_2;
	  }
	  bufptr = loadbuf;
	}
	marker_idx = uii * markers_per_pass;
	ii_shift = (indiv_idx % 4) * 2;
	wbufptr = &(writebuf[indiv_idx / 4]);
	if (map_is_unsorted) {
	  umm = 0;
	  for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
	    aptr1 = bufptr;
	    bufptr = item_endnn(bufptr);
	    alen1 = (uintptr_t)(bufptr - aptr1);
	    bufptr = skip_initial_spaces(bufptr);
	    aptr1[alen1] = '\0';
	    aptr2 = bufptr;
	    bufptr = item_endnn(bufptr);
	    alen2 = (uintptr_t)(bufptr - aptr2);
	    bufptr = skip_initial_spaces(bufptr);
	    aptr2[alen2] = '\0';
	    if (IS_SET(marker_exclude, marker_uidx)) {
	      continue;
	    }
	    ukk = map_reverse[umm++];
	    if ((ukk >= marker_start) && (ukk < marker_end)) {
	      ucc = 1;
	      if (!strcmp(aptr1, marker_allele_ptrs[2 * ukk + 1])) {
		if ((alen1 == alen2) && (!memcmp(aptr1, aptr2, alen1))) {
		  ucc = 3;
		} else if (!strcmp(aptr2, marker_allele_ptrs[2 * ukk])) {
		  ucc = 2;
		}
	      } else if (!strcmp(aptr1, marker_allele_ptrs[2 * ukk])) {
		if ((alen1 == alen2) && (!memcmp(aptr1, aptr2, alen1))) {
		  ucc = 0;
		} else if (!strcmp(aptr2, marker_allele_ptrs[2 * ukk + 1])) {
		  ucc = 2;
		}
	      }
	      wbufptr[(ukk - marker_start) * indiv_ct4] |= ucc << ii_shift;
	      marker_idx++;
	    }
	  }
	} else {
	  for (marker_uidx = umm; marker_idx < marker_end; marker_uidx++) {
	    aptr1 = bufptr;
	    bufptr = item_endnn(bufptr);
	    alen1 = (uintptr_t)(bufptr - aptr1);
	    bufptr = skip_initial_spaces(bufptr);
	    aptr1[alen1] = '\0';
	    aptr2 = bufptr;
	    bufptr = item_endnn(bufptr);
	    alen2 = (uintptr_t)(bufptr - aptr2);
	    bufptr = skip_initial_spaces(bufptr);
	    aptr2[alen2] = '\0';
	    if (is_set(marker_exclude, marker_uidx)) {
	      continue;
	    }
	    ucc = 1;
	    if (!strcmp(aptr1, marker_allele_ptrs[2 * marker_idx + 1])) {
	      if ((alen1 == alen2) && (!memcmp(aptr1, aptr2, alen1))) {
		ucc = 3;
	      } else if (!strcmp(aptr2, marker_allele_ptrs[2 * marker_idx])) {
		ucc = 2;
	      }
	    } else if (!strcmp(aptr1, marker_allele_ptrs[2 * marker_idx])) {
	      if ((alen1 == alen2) && (!memcmp(aptr1, aptr2, alen1))) {
		ucc = 0;
	      } else if (!strcmp(aptr2, marker_allele_ptrs[2 * marker_idx + 1])) {
		ucc = 2;
	      }
	    }
	    *wbufptr |= ucc << ii_shift;
	    wbufptr = &(wbufptr[indiv_ct4]);
	    marker_idx++;
	  }
	  if (!last_pass) {
	    line_starts[indiv_idx] = ped_next_thresh + (uintptr_t)(bufptr - loadbuf);
	  }
	}
      }
      if (pct > 10) {
	putchar('\b');
      }
      printf("\b\b%u%%", pct);
      fflush(stdout);
    }
    if (fwrite_checked(writebuf, ujj * indiv_ct4, *outfile_ptr)) {
      goto ped_to_bed_multichar_allele_ret_WRITE_FAIL_2;
    }
    if (!last_pass) {
      printf("\rPass %u:    \b\b\b", uii + 2);
      fflush(stdout);
      if (map_is_unsorted) {
	rewind(*pedfile_ptr);
      } else {
	umm = marker_uidx;
      }
    }
  }

  while (0) {
  ped_to_bed_multichar_allele_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  ped_to_bed_multichar_allele_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  ped_to_bed_multichar_allele_ret_READ_FAIL_2:
    putchar('\n');
  ped_to_bed_multichar_allele_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  ped_to_bed_multichar_allele_ret_WRITE_FAIL_2:
    putchar('\n');
  ped_to_bed_multichar_allele_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  ped_to_bed_multichar_allele_ret_INVALID_FORMAT_6:
    wkspace_base -= cur_slen_rdup;
    wkspace_left += cur_slen_rdup;
    putchar('\n');
    if (retval != RET_NOMEM) {
      sprintf(logbuf, "Error: More than 4 different alleles at variant %u%s.\n", uii + 1, map_is_unsorted? " (post-sort/filter)" : "");
      logprintb();
    }
    break;
  ped_to_bed_multichar_allele_ret_INVALID_FORMAT_4:
    wkspace_base -= cur_slen_rdup;
    wkspace_left += cur_slen_rdup;
    putchar('\n');
    sprintf(logbuf, "Error: Half-missing call in .ped file at variant %" PRIuPTR ", indiv %" PRIuPTR ".\n", marker_uidx + 1, indiv_ct + 1);
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  ped_to_bed_multichar_allele_ret_INVALID_FORMAT_3:
    putchar('\n');
    sprintf(logbuf, "Error: Not enough variants in .ped line %" PRIuPTR ".\n", indiv_ct + 1);
  ped_to_bed_multichar_allele_ret_INVALID_FORMAT_2:
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
  // no marker_allele_ptrs free since all strings were allocated on top of
  // stack
  return retval;
}

int32_t ped_to_bed(char* pedname, char* mapname, char* outname, char* outname_end, uint32_t fam_cols, uint32_t affection_01, int32_t missing_pheno, Chrom_info* chrom_info_ptr) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* mapfile = NULL;
  FILE* pedfile = NULL;
  FILE* outfile = NULL;
  uintptr_t* marker_exclude;
  uintptr_t max_marker_id_len = 0;
  uintptr_t unfiltered_marker_ct = 0;
  uintptr_t marker_exclude_ct = 0;
  uintptr_t marker_ct = 0;
  uint32_t map_is_unsorted = 0;
  int32_t last_chrom = 0;
  uint32_t last_mpos = 0;
  uintptr_t indiv_ct = 0;
  uint32_t ped_buflen = 1;
  int32_t retval = 0;
  uint32_t ped_col_skip = 1 + ((fam_cols & FAM_COL_1) / FAM_COL_1) + 2 * ((fam_cols & FAM_COL_34) / FAM_COL_34) + ((fam_cols & FAM_COL_5) / FAM_COL_5) + ((fam_cols & FAM_COL_6) / FAM_COL_6);
  uint32_t last_pass = 0;
  int64_t* line_starts = NULL;
  uintptr_t max_marker_allele_len = 1;
  uint32_t pass_ct;
  uintptr_t indiv_ct4;
  uint32_t pct;
  char* marker_alleles_f;
  char* marker_alleles;
  uint32_t* marker_allele_cts;
  uint32_t* map_reverse;
  uint32_t cm_col;
  uint32_t markers_per_pass;
  uint32_t marker_start;
  uint32_t marker_end;
  uintptr_t marker_uidx;
  uintptr_t marker_idx;
  uint32_t loop_end;
  uintptr_t indiv_idx;
  uintptr_t ulii;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  int32_t ii;
  char* loadbuf;
  uintptr_t loadbuf_size;
  char* col1_ptr;
  char* col2_ptr;
  char* bufptr;
  char* bufptr2;
  char cc;
  char cc2;
  unsigned char ucc;
  uint32_t ii_shift;
  unsigned char* writebuf;
  unsigned char* wbufptr;
  int64_t ped_size;
  int64_t ped_next_thresh;
  marker_exclude = (uintptr_t*)wkspace_base;
  marker_exclude[0] = 0;
  if (fopen_checked(&mapfile, mapname, "r")) {
    goto ped_to_bed_ret_OPEN_FAIL;
  }
  tbuf[MAXLINELEN - 6] = ' ';
  if (check_cm_col(mapfile, tbuf, 0, MAXLINELEN - 5, &cm_col)) {
    sprintf(logbuf, "Error: Missing token(s) in .map line: %s", skip_initial_spaces(tbuf));
    goto ped_to_bed_ret_INVALID_FORMAT_2;
  }
  do {
    if (!tbuf[MAXLINELEN - 6]) {
      sprintf(logbuf, "Error: Excessively long line in .map file (max %u chars).\n", MAXLINELEN - 8);
      goto ped_to_bed_ret_INVALID_FORMAT_2;
    }
    col1_ptr = skip_initial_spaces(tbuf);
    if (is_eoln_or_comment(*col1_ptr)) {
      continue;
    }
    col2_ptr = next_item(col1_ptr);
    bufptr = next_item_mult(col2_ptr, 1 + cm_col);
    if (no_more_items_kns(bufptr)) {
      sprintf(logbuf, "Error: Missing token(s) in .map line: %s", col1_ptr);
      goto ped_to_bed_ret_INVALID_FORMAT_2;
    }
    ii = get_chrom_code(chrom_info_ptr, col1_ptr);
    if (ii == -1) {
      // no need to support allow_extra_chroms in .map files
      logprint("Error: Invalid chromosome index in .map file.\n");
      goto ped_to_bed_ret_INVALID_FORMAT;
    }
    if ((*bufptr == '-') || (!is_set(chrom_info_ptr->chrom_mask, ii))) {
      SET_BIT(marker_exclude, unfiltered_marker_ct);
      marker_exclude_ct++;
    } else {
      if ((*bufptr < '0') || (*bufptr > '9')) {
	logprint("Error: Non-numeric variant position in .map file.\n");
	goto ped_to_bed_ret_INVALID_FORMAT;
      }
      if (!map_is_unsorted) {
	uii = atoi(bufptr);
	if ((ii < last_chrom) || ((ii == last_chrom) && (uii < last_mpos))) {
	  map_is_unsorted = 1;
	}
	last_chrom = ii;
	last_mpos = uii;
      }
      uii = strlen_se(col2_ptr) + 1;
      if (uii > max_marker_id_len) {
	max_marker_id_len = uii;
      }
    }
    unfiltered_marker_ct++;
    if (unfiltered_marker_ct > 0x7fffffff) {
      logprint("Error: Too many variants in .map file (max 2147483647).\n");
      goto ped_to_bed_ret_INVALID_FORMAT;
    }
    if (!(unfiltered_marker_ct & (BITCT - 1))) {
      if ((unfiltered_marker_ct / 8) == wkspace_left) {
	goto ped_to_bed_ret_NOMEM;
      }
      marker_exclude[unfiltered_marker_ct / BITCT] = 0;
    }
  } while (fgets(tbuf, MAXLINELEN - 5, mapfile));
  if (!feof(mapfile)) {
    goto ped_to_bed_ret_READ_FAIL;
  }
  marker_ct = unfiltered_marker_ct - marker_exclude_ct;
  if (!marker_ct) {
    logprint("Error: No variants in current analysis.\n");
    goto ped_to_bed_ret_ALL_MARKERS_EXCLUDED;
  }
  marker_exclude = (uintptr_t*)wkspace_alloc(((unfiltered_marker_ct + (BITCT - 1)) / BITCT) * sizeof(intptr_t));

  // provisionally assume max_marker_allele_len == 1
  if (wkspace_alloc_c_checked(&marker_alleles_f, marker_ct * 2)) {
    goto ped_to_bed_ret_NOMEM;
  }
  if (map_is_unsorted) {
    retval = load_sort_and_write_map(&map_reverse, mapfile, 3 + cm_col, outname, outname_end, unfiltered_marker_ct, marker_exclude, marker_ct, max_marker_id_len, 1, chrom_info_ptr);
    if (retval) {
      goto ped_to_bed_ret_1;
    }
    cm_col = 1;
    fclose_null(&mapfile);
  }
  if (wkspace_alloc_c_checked(&marker_alleles, marker_ct * 4) ||
      wkspace_alloc_ui_checked(&marker_allele_cts, marker_ct * 4 * sizeof(int32_t))) {
    goto ped_to_bed_ret_NOMEM;
  }
  memset(marker_alleles, 0, marker_ct * 4);

  // first .ped scan: count indivs, write .fam, note alleles at each locus
  if (fopen_checked(&pedfile, pedname, "rb")) {
    goto ped_to_bed_ret_OPEN_FAIL;
  }
  memcpy(outname_end, ".fam", 5);
  if (fopen_checked(&outfile, outname, "w")) {
    goto ped_to_bed_ret_OPEN_FAIL;
  }
  loadbuf = (char*)wkspace_base;
  loadbuf_size = wkspace_left;
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  } else if (loadbuf_size <= MAXLINELEN) {
    goto ped_to_bed_ret_NOMEM;
  }
  if (fseeko(pedfile, 0, SEEK_END)) {
    goto ped_to_bed_ret_READ_FAIL;
  }
  ped_size = ftello(pedfile);
  rewind(pedfile);
  logprint("Scanning .ped file...");
  fputs(" 0%", stdout);
  fflush(stdout);
  ped_next_thresh = ped_size / 100;
  loadbuf[loadbuf_size - 1] = ' ';
  pct = 0;
  while (fgets(loadbuf, loadbuf_size, pedfile)) {
    if (!loadbuf[loadbuf_size - 1]) {
      if (loadbuf_size == MAXLINEBUFLEN) {
	sprintf(logbuf, "Error: " PROG_NAME_CAPS " does not currently support .ped lines longer than %u\nbytes.\n", MAXLINEBUFLEN);
	goto ped_to_bed_ret_INVALID_FORMAT_2;
      } else {
        goto ped_to_bed_ret_NOMEM;
      }
    }
    col1_ptr = skip_initial_spaces(loadbuf);
    if (is_eoln_or_comment(*col1_ptr)) {
      ulii = strlen(loadbuf) + 1;
      if (ulii > ped_buflen) {
	ped_buflen = ulii;
      }
      continue;
    }
    if (fam_cols & FAM_COL_1) {
      col2_ptr = next_item(col1_ptr);
    } else {
      col2_ptr = col1_ptr;
    }
    bufptr = next_item_mult(col2_ptr, ped_col_skip - 1);
    if (no_more_items_kns(bufptr)) {
      sprintf(logbuf, "\nError: Missing token(s) in .ped line: %s\n", col1_ptr);
      goto ped_to_bed_ret_INVALID_FORMAT_2;
    }
    if ((bufptr - col1_ptr) > (MAXLINELEN / 2) - 4) {
      logprint("\nError: Pathologically long token(s) in .ped file.\n");
      goto ped_to_bed_ret_INVALID_FORMAT;
    }
    if (fwrite_checked(col1_ptr, strlen_se(col1_ptr), outfile)) {
      goto ped_to_bed_ret_WRITE_FAIL;
    }
    putc('\t', outfile);
    bufptr2 = write_item(col2_ptr, outfile);
    if (fam_cols & FAM_COL_34) {
      bufptr2 = write_item(bufptr2, outfile);
      bufptr2 = write_item(bufptr2, outfile);
    } else {
      fwrite("0\t0\t", 1, 4, outfile);
    }
    if (fam_cols & FAM_COL_5) {
      bufptr2 = write_item_nt(bufptr2, outfile);
    } else {
      putc('0', outfile);
    }
    putc('\t', outfile);
    if (fam_cols & FAM_COL_6) {
      fwrite(bufptr2, 1, strlen_se(bufptr2), outfile);
    } else {
      putc('-', outfile);
      putc('9', outfile);
    }
    if (putc_checked('\n', outfile)) {
      goto ped_to_bed_ret_WRITE_FAIL;
    }
    marker_idx = 0;
    for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
      cc = *bufptr++;
      if (!cc) {
        goto ped_to_bed_ret_INVALID_FORMAT_3;
      }
      bufptr = skip_initial_spaces(bufptr);
      cc2 = *bufptr++;
      if (!cc2) {
	goto ped_to_bed_ret_INVALID_FORMAT_3;
      }
      bufptr = skip_initial_spaces(bufptr);
      if (IS_SET(marker_exclude, marker_uidx)) {
	continue;
      }
      if (cc == '0') {
	if (cc2 != '0') {
	  max_marker_allele_len = 2;
	  break;
	}
	marker_idx++;
	continue;
      } else if (cc2 == '0') {
	max_marker_allele_len = 2;
	break;
      }
      uii = 4 * (map_is_unsorted? map_reverse[marker_idx] : marker_idx);
      if (incr_text_allele0(cc, &(marker_alleles[uii]), &(marker_allele_cts[uii])) ||
	  incr_text_allele0(cc2, &(marker_alleles[uii]), &(marker_allele_cts[uii]))) {
	max_marker_allele_len = 2;
	break;
      }
      marker_idx++;
    }
    if ((max_marker_allele_len == 2) || (!is_eoln_kns(*bufptr))) {
      // either multi-character alleles, or invalid format.  Restart scan.
      putchar('\r');
      logstr("\n");
      logprint("Possibly irregular .ped line.  Restarting scan, assuming multichar alleles.\n");
      max_marker_allele_len = 2;
      break;
    }
    ulii = strlen(bufptr) + (uintptr_t)(bufptr - loadbuf) + 1;
    if (ulii > ped_buflen) {
      ped_buflen = ulii;
    }
    indiv_ct++;
    if (ftello(pedfile) >= ped_next_thresh) {
      uii = (ftello(pedfile) * 100) / ped_size;
      if (pct >= 10) {
	putchar('\b');
      }
      printf("\b\b%u%%", uii);
      fflush(stdout);
      pct = uii;
    }
  }
  if (max_marker_allele_len == 1) {
    if (!feof(pedfile)) {
      goto ped_to_bed_ret_READ_FAIL;
    }
    if (!indiv_ct) {
      sprintf(logbuf, "\nError: No %s in .ped file.\n", g_species_plural);
      goto ped_to_bed_ret_INVALID_FORMAT_2;
    }
    if (fclose_null(&outfile)) {
      goto ped_to_bed_ret_WRITE_FAIL;
    }
    memcpy(outname_end, ".bim", 5);
    if (fopen_checked(&outfile, outname, "w")) {
      goto ped_to_bed_ret_OPEN_FAIL;
    }
    if (map_is_unsorted) {
      memcpy(outname_end, ".map.tmp", 9);
      if (fopen_checked(&mapfile, outname, "r")) {
	goto ped_to_bed_ret_OPEN_FAIL;
      }
    } else {
      rewind(mapfile);
    }
    logstr(" done.\n");
    fputs("\r.ped scan complete (for binary autoconversion).\n", stdout);
    marker_uidx = 0;
    for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
      if (map_is_unsorted) {
	if (!fgets(tbuf, MAXLINELEN, mapfile)) {
	  goto ped_to_bed_ret_READ_FAIL;
	}
      } else {
	if (get_next_noncomment_excl(mapfile, &bufptr, marker_exclude, &marker_uidx)) {
	  goto ped_to_bed_ret_READ_FAIL;
	}
      }
      if (marker_alleles[marker_idx * 4 + 2]) {
	cc = marker_alleles[marker_idx * 4 + 3];
	if (map_is_unsorted) {
	  sprintf(logbuf, "Warning: Variant %u (post-sort/filter) %sallelic; setting rarest missing.\n", map_reverse[marker_idx] + 1, (cc? "quad" : "tri"));
	} else {
	  sprintf(logbuf, "Warning: Variant %" PRIuPTR " %sallelic; setting rarest alleles missing.\n", marker_idx + 1, (cc? "quad" : "tri"));
	}
	logprintb();
	ujj = (cc? 4 : 3);
	// insertion sort
	for (uii = 1; uii < ujj; uii++) {
	  ukk = marker_allele_cts[4 * marker_idx + uii];
	  if (marker_allele_cts[4 * marker_idx + uii - 1] < ukk) {
	    cc = marker_alleles[4 * marker_idx + uii];
	    umm = uii;
	    do {
	      umm--;
	      marker_alleles[4 * marker_idx + umm + 1] = marker_alleles[4 * marker_idx + umm];
	      marker_allele_cts[4 * marker_idx + umm + 1] = marker_allele_cts[4 * marker_idx + umm];
	    } while (umm && (marker_allele_cts[4 * marker_idx + umm - 1] < ukk));
	    marker_alleles[4 * marker_idx + umm] = cc;
	    marker_allele_cts[4 * marker_idx + umm] = ukk;
	  }
	}
	cc = marker_alleles[marker_idx * 4 + 1];
	cc2 = marker_alleles[marker_idx * 4];
      } else {
	if (marker_allele_cts[marker_idx * 4] >= marker_allele_cts[marker_idx * 4 + 1]) {
	  cc = marker_alleles[marker_idx * 4 + 1];
	  cc2 = marker_alleles[marker_idx * 4];
	} else {
	  cc = marker_alleles[marker_idx * 4];
	  cc2 = marker_alleles[marker_idx * 4 + 1];
	}
      }
      marker_alleles_f[marker_idx * 2] = cc;
      marker_alleles_f[marker_idx * 2 + 1] = cc2;
      if (!cc) {
	cc = '0';
      }
      if (!cc2) {
	cc2 = '0';
      }
      if (map_is_unsorted) {
	bufptr = (char*)memchr(tbuf, '\n', MAXLINELEN);
	if (fwrite_checked(tbuf, bufptr - tbuf, outfile)) {
	  goto ped_to_bed_ret_WRITE_FAIL;
	}
      } else {
	bufptr = write_item(bufptr, outfile);
	bufptr = write_item(bufptr, outfile);
	if (cm_col) {
	  bufptr = write_item_nt(bufptr, outfile);
	} else {
	  putc('0', outfile);
	}
	putc('\t', outfile);
	fwrite(bufptr, 1, strlen_se(bufptr), outfile);
      }
      putc('\t', outfile);
      putc(cc, outfile);
      putc('\t', outfile);
      putc(cc2, outfile);
      if (putc_checked('\n', outfile)) {
	goto ped_to_bed_ret_WRITE_FAIL;
      }
      marker_uidx++;
    }
    indiv_ct4 = (indiv_ct + 3) / 4;
    wkspace_reset((unsigned char*)marker_alleles);
    fclose_null(&mapfile);
    if (map_is_unsorted) {
      unlink(outname);
    }
    fclose_null(&outfile);
    if (wkspace_alloc_c_checked(&loadbuf, ped_buflen)) {
      goto ped_to_bed_ret_NOMEM;
    }
    if (wkspace_left >= marker_ct * indiv_ct4) {
      markers_per_pass = marker_ct;
      sprintf(logbuf, "Performing single-pass .bed write (%" PRIuPTR " variant%s, %" PRIuPTR " %s).\n", marker_ct, (marker_ct == 1)? "" : "s", indiv_ct, species_str(indiv_ct));
      pass_ct = 1;
    } else {
      if (!map_is_unsorted) {
	if (wkspace_alloc_ll_checked(&line_starts, indiv_ct * sizeof(int64_t))) {
	  goto ped_to_bed_ret_NOMEM;
	}
      }
      markers_per_pass = wkspace_left / indiv_ct4;
      if (!markers_per_pass) {
	goto ped_to_bed_ret_NOMEM;
      }
      pass_ct = (marker_ct + markers_per_pass - 1) / markers_per_pass;
      sprintf(logbuf, "Performing %u-pass .bed write (%u/%" PRIuPTR " variant%s/pass, %" PRIuPTR " %s).\n", pass_ct, markers_per_pass, marker_ct, (markers_per_pass == 1)? "" : "s", indiv_ct, species_str(indiv_ct));
    }
    logprintb();
    writebuf = wkspace_base;
    memcpy(outname_end, ".bed", 5);
    if (fopen_checked(&outfile, outname, "wb")) {
      goto ped_to_bed_ret_OPEN_FAIL;
    }
    if (fwrite_checked("l\x1b\x01", 3, outfile)) {
      goto ped_to_bed_ret_WRITE_FAIL;
    }
    rewind(pedfile);
    umm = 0;
    for (uii = 0; uii < pass_ct; uii++) {
      marker_start = uii * markers_per_pass;
      if (uii + 1 == pass_ct) {
	ujj = marker_ct - marker_start;
	last_pass = 1;
      } else {
	ujj = markers_per_pass;
      }
      memset(writebuf, 0, ujj * indiv_ct4);
      marker_end = marker_start + ujj;
      fputs("0%", stdout);
      indiv_idx = 0;
      // 94 instead of 100 due to big fwrite at the end
      for (pct = 1; pct <= 94; pct++) {
	loop_end = (((uint64_t)pct) * indiv_ct) / 94LLU;
	for (; indiv_idx < loop_end; indiv_idx++) {
	  if ((!uii) || map_is_unsorted) {
	    do {
	      if (!last_pass) {
		ped_next_thresh = ftello(pedfile);
	      }
	      if (!fgets(loadbuf, ped_buflen, pedfile)) {
		goto ped_to_bed_ret_READ_FAIL_2;
	      }
	      col1_ptr = skip_initial_spaces(loadbuf);
	    } while (is_eoln_or_comment(*col1_ptr));
	    bufptr = next_item_mult(col1_ptr, ped_col_skip);
	  } else {
	    ped_next_thresh = line_starts[indiv_idx];
	    if (fseeko(pedfile, line_starts[indiv_idx], SEEK_SET)) {
	      goto ped_to_bed_ret_READ_FAIL_2;
	    }
	    if (!fgets(loadbuf, ped_buflen, pedfile)) {
	      goto ped_to_bed_ret_READ_FAIL_2;
	    }
	    bufptr = loadbuf;
	  }
	  marker_idx = uii * markers_per_pass;
	  ii_shift = (indiv_idx % 4) * 2;
	  wbufptr = &(writebuf[indiv_idx / 4]);
	  if (map_is_unsorted) {
	    // multipass optimizations are possible, but we won't bother,
	    // especially since the .map should rarely be unsorted in the first
	    // place...
	    umm = 0;
	    for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
	      cc = *bufptr++;
	      bufptr = skip_initial_spaces(bufptr);
	      cc2 = *bufptr++;
	      bufptr = skip_initial_spaces(bufptr);
	      if (IS_SET(marker_exclude, marker_uidx)) {
		continue;
	      }
	      ukk = map_reverse[umm++];
	      if ((ukk >= marker_start) && (ukk < marker_end)) {
		ucc = 1;
		if (cc == marker_alleles_f[2 * ukk + 1]) {
		  if (cc2 == cc) {
		    ucc = 3;
		  } else if (cc2 == marker_alleles_f[2 * ukk]) {
		    ucc = 2;
		  }
		} else if (cc == marker_alleles_f[2 * ukk]) {
		  if (cc2 == cc) {
		    ucc = 0;
		  } else if (cc2 == marker_alleles_f[2 * ukk + 1]) {
		    ucc = 2;
		  }
		}
		wbufptr[(ukk - marker_start) * indiv_ct4] |= ucc << ii_shift;
		marker_idx++;
	      }
	    }
	  } else {
	    for (marker_uidx = umm; marker_idx < marker_end; marker_uidx++) {
	      cc = *bufptr++;
	      bufptr = skip_initial_spaces(bufptr);
	      cc2 = *bufptr++;
	      bufptr = skip_initial_spaces(bufptr);
	      if (IS_SET(marker_exclude, marker_uidx)) {
		continue;
	      }
	      ucc = 1;
	      if (cc == marker_alleles_f[2 * marker_idx + 1]) {
		if (cc2 == cc) {
		  ucc = 3;
		} else if (cc2 == marker_alleles_f[2 * marker_idx]) {
		  ucc = 2;
		}
	      } else if (cc == marker_alleles_f[2 * marker_idx]) {
		if (cc2 == cc) {
		  ucc = 0;
		} else if (cc2 == marker_alleles_f[2 * marker_idx + 1]) {
		  ucc = 2;
		}
	      }
	      *wbufptr |= ucc << ii_shift;
	      wbufptr = &(wbufptr[indiv_ct4]);
	      marker_idx++;
	    }
	    if (!last_pass) {
	      line_starts[indiv_idx] = ped_next_thresh + (uintptr_t)(bufptr - loadbuf);
	    }
	  }
	}
	if (pct > 10) {
	  putchar('\b');
	}
	printf("\b\b%u%%", pct);
	fflush(stdout);
      }
      if (fwrite_checked(writebuf, ujj * indiv_ct4, outfile)) {
	goto ped_to_bed_ret_WRITE_FAIL_2;
      }
      if (!last_pass) {
	printf("\rPass %u:    \b\b\b", uii + 2);
	fflush(stdout);
	if (map_is_unsorted) {
	  rewind(pedfile);
	} else {
	  umm = marker_uidx;
	}
      }
    }
  } else {
    retval = ped_to_bed_multichar_allele(&pedfile, &outfile, outname, outname_end, &mapfile, unfiltered_marker_ct, marker_exclude, marker_ct, marker_alleles_f, map_is_unsorted, fam_cols, ped_col_skip, cm_col, map_reverse, ped_size);
    if (retval) {
      goto ped_to_bed_ret_1;
    }
  }

  if (fclose_null(&outfile)) {
    goto ped_to_bed_ret_WRITE_FAIL_2;
  }
  putchar('\r');
  sprintf(logbuf, "--file: binary fileset written to %s + .bim + .fam.\n", outname);
  logprintb();

  while (0) {
  ped_to_bed_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  ped_to_bed_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  ped_to_bed_ret_READ_FAIL_2:
    putchar('\n');
  ped_to_bed_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  ped_to_bed_ret_WRITE_FAIL_2:
    putchar('\n');
  ped_to_bed_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  ped_to_bed_ret_INVALID_FORMAT_3:
    sprintf(logbuf, "Error: Not enough variants in .ped line %" PRIuPTR ".\n", indiv_ct + 1);
  ped_to_bed_ret_INVALID_FORMAT_2:
    logprintb();
  ped_to_bed_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  ped_to_bed_ret_ALL_MARKERS_EXCLUDED:
    retval = RET_ALL_MARKERS_EXCLUDED;
    break;
  }
 ped_to_bed_ret_1:
  fclose_cond(pedfile);
  fclose_cond(mapfile);
  fclose_cond(outfile);
  wkspace_reset(wkspace_mark);
  return retval;
}

void fill_bmap_short(unsigned char* bmap_short, uint32_t ct_mod4) {
  unsigned char imap[4] = {3, 1, 2, 0};
  unsigned char* bmap = bmap_short;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  uint32_t unn;
  uint32_t extra_ct;
  for (uii = 0; uii < 4; uii++) {
    umm = imap[uii] * 64;
    for (ujj = 0; ujj < 4; ujj++) {
      unn = umm + imap[ujj] * 16;
      for (ukk = 0; ukk < 4; ukk++) {
	extra_ct = unn + imap[ukk] * 4;
	*bmap++ = extra_ct + 3;
	*bmap++ = extra_ct + 1;
	*bmap++ = extra_ct + 2;
	*bmap++ = extra_ct;
      }
    }
  }
  if (ct_mod4) {
    switch (ct_mod4) {
    case 1:
      bmap[0] = 3;
      bmap[1] = 1;
      bmap[2] = 2;
      bmap[3] = 0;
      return;
    case 2:
      for (uii = 0; uii < 4; uii++) {
	ukk = imap[uii] * 4;
	*bmap++ = ukk + 3;
	*bmap++ = ukk + 1;
	*bmap++ = ukk + 2;
	*bmap++ = ukk;
      }
      return;
    case 3:
      for (uii = 0; uii < 4; uii++) {
	ukk = imap[uii] * 16;
	for (ujj = 0; ujj < 4; ujj++) {
	  umm = ukk + imap[ujj] * 4;
	  *bmap++ = umm + 3;
	  *bmap++ = umm + 1;
	  *bmap++ = umm + 2;
	  *bmap++ = umm;
	}
      }
    }
  }
}

void fill_bmap_raw(unsigned char* bmap_raw, uint32_t ct_mod4) {
  // possibilities:
  // 0. A1 -> A1, A2 -> A2 [0..255]
  // 1. A1 -> A2, A2 -> A1 [256..511]
  // 1 last char. [512..575]
  uint32_t uii = 0;
  do {
    *bmap_raw++ = uii++;
  } while (uii < 256);
  fill_bmap_short(bmap_raw, ct_mod4);
}

int32_t lgen_to_bed(char* lgen_namebuf, char* outname, char* outname_end, int32_t missing_pheno, uint64_t misc_flags, uint32_t lgen_modifier, char* lgen_reference_fname, Chrom_info* chrom_info_ptr) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* infile = NULL;
  FILE* outfile = NULL;
  char* name_end = (char*)memchr(lgen_namebuf, 0, FNAMESIZE);
  uint32_t lgen_allele_count = lgen_modifier & LGEN_ALLELE_COUNT;
  uint32_t affection_01 = (misc_flags / MISC_AFFECTION_01) & 1;
  uint32_t map_cols = 3;
  uintptr_t* marker_exclude = NULL;
  uintptr_t marker_exclude_ct = 0;
  uintptr_t max_marker_id_len = 0;
  uintptr_t unfiltered_marker_ct = 0;
  uintptr_t max_person_id_len = 4;
  uintptr_t marker_ct = 0;
  char* missing_geno_ptr = (char*)g_missing_geno_ptr;
  char** marker_allele_ptrs = NULL;
  char* marker_ids = NULL;
  uint32_t* marker_pos = NULL;
  uintptr_t indiv_ct = 0;
  char* person_ids = NULL;
  char* paternal_ids = NULL;
  uintptr_t max_paternal_id_len = 2;
  char* maternal_ids = NULL;
  uintptr_t max_maternal_id_len = 2;
  uintptr_t* sex_nm = NULL;
  uintptr_t* sex_male = NULL;
  uint32_t affection = 0;
  uintptr_t* founder_info = NULL;
  uintptr_t* indiv_exclude = NULL;
  uint32_t map_is_unsorted = 0;
  uint32_t compound_genotypes = 1; // 0 = no, 1 = unresolved, 2 = yes
  char missing_geno = *missing_geno_ptr;
  uintptr_t* pheno_nm = NULL;
  uintptr_t* pheno_c = NULL;
  double* pheno_d = NULL;
  char* sorted_marker_ids;
  uint32_t* marker_id_map;
  uint32_t* map_reverse;
  char* sorted_indiv_ids;
  uint32_t* indiv_id_map;
  unsigned char* writebuf;
  uintptr_t indiv_ct4;
  uintptr_t marker_idx;
  unsigned char ucc;
  unsigned char* ucptr;
  unsigned char bmap_short[320];
  unsigned char* bmap2;
  char* loadbuf;
  char* id_buf;
  char* cptr;
  char* cptr2;
  char* cptr3;
  char* cptr4;
  char* a1ptr;
  char* a2ptr;
  char* sptr;
  char* sptr2;
  char** ma_end;
  uint32_t a1len;
  uint32_t a2len;
  uintptr_t loadbuf_size;
  uintptr_t indiv_idx;
  uintptr_t ulii;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  int64_t lgen_size;
  uint32_t pct;
  int64_t lgen_next_thresh;
  int32_t ii;
  int32_t retval;
  if (lgen_modifier == LGEN_ALLELE_COUNT) {
    logprint("Error: --allele-count must be used with --reference.\n");
    goto lgen_to_bed_ret_INVALID_CMDLINE;
  }

  memcpy(name_end, ".map", 5);
  retval = load_map(&infile, lgen_namebuf, &map_cols, &unfiltered_marker_ct, &marker_exclude_ct, &max_marker_id_len, &marker_exclude, &marker_ids, chrom_info_ptr, &marker_pos, &map_is_unsorted);
  if (retval) {
    goto lgen_to_bed_ret_1;
  }
  marker_ct = unfiltered_marker_ct - marker_exclude_ct;
  retval = sort_item_ids(&sorted_marker_ids, &marker_id_map, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_ids, max_marker_id_len, 0, 0, strcmp_deref);
  if (retval) {
    goto lgen_to_bed_ret_1;
  }
  if (map_is_unsorted) {
    // Writes a temporary .map which is read later, and then deleted.
    retval = load_sort_and_write_map(&map_reverse, infile, map_cols, outname, outname_end, unfiltered_marker_ct, marker_exclude, marker_ct, max_marker_id_len, 0, chrom_info_ptr);
    if (retval) {
      goto lgen_to_bed_ret_1;
    }
    for (uii = 0; uii < marker_ct; uii++) {
      marker_id_map[uii] = map_reverse[marker_id_map[uii]];
    }
  }
  // collapse
  if (wkspace_alloc_ui_checked(&indiv_id_map, unfiltered_marker_ct * sizeof(int32_t))) {
    goto lgen_to_bed_ret_NOMEM;
  }
  fill_uidx_to_idx(marker_exclude, unfiltered_marker_ct, marker_ct, indiv_id_map);
  for (uii = 0; uii < marker_ct; uii++) {
    marker_id_map[uii] = indiv_id_map[marker_id_map[uii]];
  }
  fclose_null(&infile);
  memcpy(marker_ids, sorted_marker_ids, marker_ct * max_marker_id_len);
  wkspace_reset((unsigned char*)sorted_marker_ids);

  memcpy(name_end, ".fam", 5);
  if (fopen_checked(&infile, lgen_namebuf, "r")) {
    goto lgen_to_bed_ret_OPEN_FAIL;
  }
  retval = load_fam(infile, MAXLINELEN, FAM_COL_13456, 1, missing_pheno, intlen(missing_pheno), affection_01, &indiv_ct, &person_ids, &max_person_id_len, &paternal_ids, &max_paternal_id_len, &maternal_ids, &max_maternal_id_len, &sex_nm, &sex_male, &affection, &pheno_nm, &pheno_c, &pheno_d, &founder_info, &indiv_exclude);
  if (retval) {
    goto lgen_to_bed_ret_1;
  }
  fclose_null(&infile);
  retval = sort_item_ids_nx(&sorted_indiv_ids, &indiv_id_map, indiv_ct, person_ids, max_person_id_len);
  if (retval) {
    goto lgen_to_bed_ret_1;
  }
  if (wkspace_alloc_c_checked(&id_buf, MAXV(max_marker_id_len, max_person_id_len))) {
    goto lgen_to_bed_ret_NOMEM;
  }
  marker_allele_ptrs = (char**)wkspace_alloc(2 * marker_ct * sizeof(char*));
  if (!marker_allele_ptrs) {
    goto lgen_to_bed_ret_NOMEM;
  }
  memset(marker_allele_ptrs, 0, 2 * marker_ct * sizeof(char*));
  indiv_ct4 = (indiv_ct + 3) / 4;
  if (wkspace_alloc_uc_checked(&writebuf, ((uintptr_t)marker_ct) * indiv_ct4)) {
    logprint("Error: Multipass .lgen -> .bed autoconversions are not yet supported.  Try\nusing --chr and/or --memory (perhaps with a better machine).\n");
    goto lgen_to_bed_ret_CALC_NOT_YET_SUPPORTED;
  }
  if (indiv_ct % 4) {
    ucc = 0x15 >> (6 - 2 * (indiv_ct % 4));
    for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
      memset(&(writebuf[marker_idx * indiv_ct4]), 0x55, indiv_ct4 - 1);
      writebuf[(marker_idx + 1) * indiv_ct4 - 1] = ucc;
    }
  } else {
    memset(writebuf, 0x55, marker_ct * indiv_ct4);
  }
  loadbuf_size = wkspace_left;
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  } else if (loadbuf_size <= MAXLINELEN) {
    goto lgen_to_bed_ret_NOMEM;
  }
  loadbuf = (char*)wkspace_base;
  loadbuf[loadbuf_size - 1] = ' ';
  if (lgen_modifier & LGEN_REFERENCE) {
    if (fopen_checked(&infile, lgen_reference_fname, "r")) {
      goto lgen_to_bed_ret_OPEN_FAIL;
    }
    while (fgets(loadbuf, loadbuf_size, infile)) {
      if (!loadbuf[loadbuf_size - 1]) {
	if (loadbuf_size == MAXLINEBUFLEN) {
	  logprint("Error: Pathologically long line in .ref file.\n");
	  goto lgen_to_bed_ret_INVALID_FORMAT_6;
	}
	goto lgen_to_bed_ret_NOMEM;
      }
      cptr = skip_initial_spaces(loadbuf);
      if (is_eoln_kns(*cptr)) {
	continue;
      }
      cptr2 = item_end(cptr);
      a1ptr = skip_initial_spaces(cptr2);
      if (no_more_items_kns(a1ptr)) {
	goto lgen_to_bed_ret_INVALID_FORMAT_4;
      }
      a1len = strlen_se(cptr);
      if (a1len < max_marker_id_len) {
	memcpyx(id_buf, cptr, a1len, '\0');
	ii = bsearch_str(id_buf, marker_ids, max_marker_id_len, 0, marker_ct - 1);
	if (ii != -1) {
	  marker_idx = marker_id_map[ii];
	  if (marker_allele_ptrs[2 * marker_idx + 1]) {
	    goto lgen_to_bed_ret_INVALID_FORMAT_4;
	  }
	  sptr = item_end(a1ptr);
	  a2ptr = skip_initial_spaces(sptr);
	  a1len = (uintptr_t)(sptr - a1ptr);
	  a1ptr[a1len] = '\0';
	  if (allele_set(&(marker_allele_ptrs[2 * marker_idx + 1]), a1ptr, a1len)) {
	    goto lgen_to_bed_ret_NOMEM;
	  }
	  if (no_more_items_kns(a2ptr)) {
	    if (lgen_allele_count) {
	      a1ptr[a1len++] = 'v';
	      a1ptr[a1len] = '\0';
	      if (allele_set(&(marker_allele_ptrs[2 * marker_idx]), a1ptr, a1len)) {
		goto lgen_to_bed_ret_NOMEM;
	      }
	    }
	  } else {
	    a2len = strlen_se(a2ptr);
	    a2ptr[a2len] = '\0';
	    if (allele_set(&(marker_allele_ptrs[2 * marker_idx]), a2ptr, a2len)) {
	      goto lgen_to_bed_ret_NOMEM;
	    }
	  }
	  memset(&(writebuf[marker_idx * indiv_ct4]), 0xff, indiv_ct / 4);
	  if (indiv_ct % 4) {
	    writebuf[(marker_idx + 1) * indiv_ct4 - 1] = 0x3f >> (6 - 2 * (indiv_ct % 4));
	  }
	}
      }
    }
    if (!feof(infile)) {
      goto lgen_to_bed_ret_READ_FAIL;
    }
    fclose_null(&infile);
  }
  // PLINK 1.07 reports an error whenever there are 3+ alleles at one locus, so
  // backwards compatibility does not mandate that we worry about that case.
  // Thus we just use the obvious one-pass load, and save proper handling of
  // triallelic sites, etc. for the future .pgen engine.
  memcpy(outname_end, ".bed", 5);
  if (fopen_checked(&outfile, outname, "wb")) {
    goto lgen_to_bed_ret_OPEN_FAIL;
  }
  if (fwrite_checked("l\x1b\x01", 3, outfile)) {
    goto lgen_to_bed_ret_WRITE_FAIL;
  }
  memcpy(name_end, ".lgen", 6);
  if (fopen_checked(&infile, lgen_namebuf, "r")) {
    goto lgen_to_bed_ret_OPEN_FAIL;
  }
  if (fseeko(infile, 0, SEEK_END)) {
    goto lgen_to_bed_ret_READ_FAIL;
  }
  lgen_size = ftello(infile);
  rewind(infile);
  logprint("Processing .lgen file... ");
  fputs("0%", stdout);
  fflush(stdout);
  lgen_next_thresh = lgen_size / 100;
  pct = 0;
  if (!lgen_allele_count) {
    while (fgets(loadbuf, loadbuf_size, infile)) {
      if (!loadbuf[loadbuf_size - 1]) {
	if (loadbuf_size == MAXLINEBUFLEN) {
	  logprint("Error: Pathologically long line in .lgen file.\n");
	  goto lgen_to_bed_ret_INVALID_FORMAT_6;
	}
	goto lgen_to_bed_ret_NOMEM;
      }
      cptr = skip_initial_spaces(loadbuf);
      if (is_eoln_kns(*cptr)) {
	continue;
      }
      cptr2 = next_item(cptr);
      cptr3 = next_item(cptr2);
      cptr4 = item_end(cptr3);
      if (!cptr4) {
	goto lgen_to_bed_ret_INVALID_FORMAT;
      }
      a1ptr = skip_initial_spaces(cptr4);
      sptr = item_end(a1ptr);
      a1len = (uintptr_t)(sptr - a1ptr);
      a2ptr = next_item(sptr);
      if (compound_genotypes == 1) {
	if (no_more_items_kns(a2ptr)) {
	  compound_genotypes = 2;
	} else {
	  compound_genotypes = 0;
	}
      }
      if (!compound_genotypes) {
	if (no_more_items_kns(a2ptr)) {
	  printf("fail2\n");
	  goto lgen_to_bed_ret_INVALID_FORMAT;
	}
	a2len = strlen_se(a2ptr);
      } else {
	if (a1len != 2) {
	  goto lgen_to_bed_ret_INVALID_FORMAT;
	}
	a1len = 1;
	a2ptr = &(a1ptr[2]);
	*a2ptr = a1ptr[1];
	a2len = 1;
      }
      ii = bsearch_fam_indiv(id_buf, sorted_indiv_ids, max_person_id_len, indiv_ct, cptr, cptr2);
      if (ii == -1) {
	goto lgen_to_bed_ret_INVALID_FORMAT_2;
      }
      indiv_idx = indiv_id_map[(uint32_t)ii];
      ulii = (uintptr_t)(cptr4 - cptr3);
      memcpyx(id_buf, cptr3, ulii, '\0');
      ii = bsearch_str(id_buf, marker_ids, max_marker_id_len, 0, marker_ct - 1);
      if (ii != -1) {
	marker_idx = marker_id_map[ii];
	sptr = marker_allele_ptrs[2 * marker_idx + 1]; // existing A2
	a1ptr[a1len] = '\0';
	a2ptr[a2len] = '\0';
	if ((*a1ptr == missing_geno) && (a1len == 1)) {
	  if ((*a2ptr == missing_geno) && (a2len == 1)) {
	    uii = 1;
	  } else {
	    goto lgen_to_bed_ret_INVALID_FORMAT_5;
	  }
	} else if ((*a2ptr == missing_geno) && (a2len == 1)) {
	  goto lgen_to_bed_ret_INVALID_FORMAT_5;
        } else {
          if (!sptr) {
	    if (allele_set(&(marker_allele_ptrs[2 * marker_idx + 1]), a1ptr, a1len)) {
	      goto lgen_to_bed_ret_NOMEM;
	    }
	    if (!strcmp(a1ptr, a2ptr)) {
	      uii = 2;
	    } else {
	      uii = 1;
	      if (allele_set(&(marker_allele_ptrs[2 * marker_idx]), a2ptr, a2len)) {
		goto lgen_to_bed_ret_NOMEM;
	      }
	    }
	  } else {
	    sptr2 = marker_allele_ptrs[2 * marker_idx];
	    if (!sptr2) {
	      if (!strcmp(a1ptr, sptr)) {
		if (!strcmp(a2ptr, sptr)) {
		  uii = 2;
		} else {
		  uii = 1;
		  if (allele_set(&(marker_allele_ptrs[2 * marker_idx]), a2ptr, a2len)) {
		    goto lgen_to_bed_ret_NOMEM;
		  }
		}
	      } else {
		if (allele_set(&(marker_allele_ptrs[2 * marker_idx]), a1ptr, a1len)) {
		  goto lgen_to_bed_ret_NOMEM;
		}
		if (!strcmp(a2ptr, sptr)) {
		  uii = 1;
		} else if (!strcmp(a2ptr, a1ptr)) {
		  uii = 0;
		} else {
		  goto lgen_to_bed_ret_INVALID_FORMAT_3;
		}
	      }
	    } else {
	      if (!strcmp(a1ptr, sptr)) {
		uii = 1;
	      } else if (!strcmp(a1ptr, sptr2)) {
		uii = 0;
	      } else {
		goto lgen_to_bed_ret_INVALID_FORMAT_3;
	      }
	      if (!strcmp(a2ptr, sptr)) {
		uii++;
	      } else if (strcmp(a2ptr, sptr2)) {
		goto lgen_to_bed_ret_INVALID_FORMAT_3;
	      }
	    }
	  }
	  if (uii) {
	    uii++;
	  }
	}
	ulii = marker_idx * indiv_ct4 + (indiv_idx / 4);
	ujj = (indiv_idx % 4) * 2;
	writebuf[ulii] = (writebuf[ulii] & (~(3 << ujj))) | (uii << ujj);
      }
      if (ftello(infile) >= lgen_next_thresh) {
	uii = (ftello(infile) * 100) / lgen_size;
	if (pct >= 10) {
	  putchar('\b');
	}
	printf("\b\b%u%%", uii);
	fflush(stdout);
	pct = uii;
	lgen_next_thresh = ((pct + 1) * lgen_size) / 100;
      }
    }
  } else {
    while (fgets(loadbuf, loadbuf_size, infile)) {
      if (!loadbuf[loadbuf_size - 1]) {
	if (loadbuf_size == MAXLINEBUFLEN) {
	  logprint("Error: Pathologically long line in .lgen file.\n");
	  goto lgen_to_bed_ret_INVALID_FORMAT_6;
	}
	goto lgen_to_bed_ret_NOMEM;
      }
      cptr = skip_initial_spaces(loadbuf);
      if (is_eoln_kns(*cptr)) {
	continue;
      }
      cptr2 = next_item(cptr);
      cptr3 = next_item(cptr2);
      cptr4 = item_end(cptr3);
      if (!cptr4) {
	goto lgen_to_bed_ret_INVALID_FORMAT;
      }
      a1ptr = skip_initial_spaces(cptr4);
      if (no_more_items_kns(a1ptr)) {
	goto lgen_to_bed_ret_INVALID_FORMAT;
      }
      ii = bsearch_fam_indiv(id_buf, sorted_indiv_ids, max_person_id_len, indiv_ct, cptr, cptr2);
      if (ii == -1) {
	goto lgen_to_bed_ret_INVALID_FORMAT_2;
      }
      indiv_idx = indiv_id_map[(uint32_t)ii];
      ulii = (uintptr_t)(cptr4 - cptr3);
      memcpyx(id_buf, cptr3, ulii, '\0');
      ii = bsearch_str(id_buf, marker_ids, max_marker_id_len, 0, marker_ct - 1);
      if (ii != -1) {
	marker_idx = marker_id_map[ii];
	a1len = strlen_se(a1ptr);
	ucc = (unsigned char)(*a1ptr);
	if ((a1len != 1) || (ucc < 48) || (ucc > 50)) {
	  uii = 1;
	} else {
	  uii = ucc - 48;
	  if (uii) {
	    uii++;
	  }
	}
	ulii = marker_idx * indiv_ct4 + (indiv_idx / 4);
	ujj = (indiv_idx % 4) * 2;
	writebuf[ulii] = (writebuf[ulii] & (~(3 << ujj))) | (uii << ujj);
      }
      if (ftello(infile) >= lgen_next_thresh) {
	uii = (ftello(infile) * 100) / lgen_size;
	if (pct >= 10) {
	  putchar('\b');
	}
	printf("\b\b%u%%", uii);
	fflush(stdout);
	pct = uii;
	lgen_next_thresh = ((pct + 1) * lgen_size) / 100;
      }
    }
  }
  if (!feof(infile)) {
    goto lgen_to_bed_ret_READ_FAIL;
  }
  fclose_null(&infile);
  if (pct < 10) {
    fputs("\b\b", stdout);
  } else if (pct < 100) {
    fputs("\b\b\b", stdout);
  } else {
    fputs("\b\b\b\b", stdout);
  }
  logprint("done.\n");
  ukk = indiv_ct / 4;
  umm = indiv_ct % 4;
  fill_bmap_short(bmap_short, umm);
  bmap2 = &(bmap_short[256]);
  for (uii = 0; uii < marker_ct; uii++) {
    if (popcount_chars((uintptr_t*)writebuf, uii * indiv_ct4, (uii + 1) * indiv_ct4) < indiv_ct) {
      ucptr = &(writebuf[uii * indiv_ct4]);
      for (ujj = 0; ujj < ukk; ujj++) {
	*ucptr = bmap_short[*ucptr];
	ucptr++;
      }
      if (umm) {
        *ucptr = bmap2[*ucptr];
	ucptr++;
      }
      cptr = marker_allele_ptrs[uii * 2];
      marker_allele_ptrs[uii * 2] = marker_allele_ptrs[uii * 2 + 1];
      marker_allele_ptrs[uii * 2 + 1] = cptr;
    }
  }
  if (fwrite_checked(writebuf, ((uintptr_t)marker_ct) * indiv_ct4, outfile)) {
    goto lgen_to_bed_ret_WRITE_FAIL;
  }
  if (fclose_null(&outfile)) {
    goto lgen_to_bed_ret_WRITE_FAIL;
  }
  if (map_is_unsorted) {
    memcpy(outname_end, ".map.tmp", 9);
    if (fopen_checked(&infile, outname, "r")) {
      goto lgen_to_bed_ret_OPEN_FAIL;
    }
  } else {
    memcpy(name_end, ".map", 5);
    if (fopen_checked(&infile, lgen_namebuf, "r")) {
      goto lgen_to_bed_ret_OPEN_FAIL;
    }
  }
  memcpy(outname_end, ".bim", 5);
  if (fopen_checked(&outfile, outname, "w")) {
    goto lgen_to_bed_ret_OPEN_FAIL;
  }
  uii = 2 * marker_ct;
  for (ujj = 0; ujj < uii; ujj++) {
    if (!marker_allele_ptrs[ujj]) {
      marker_allele_ptrs[ujj] = missing_geno_ptr;
    }
  }
  uii = 0;
  marker_idx = 0;
  while (fgets(tbuf, MAXLINELEN, infile)) {
    if (is_eoln_or_comment(*(skip_initial_spaces(tbuf)))) {
      continue;
    }
    if (IS_SET(marker_exclude, uii)) {
      uii++;
      continue;
    }
    cptr = (char*)memchr(tbuf, 0, MAXLINELEN);
    if (cptr[-1] == '\n') {
      cptr--;
      if (cptr[-1] == '\r') {
	cptr--;
      }
    }
    *cptr++ = '\t';
    fwrite(tbuf, 1, cptr - tbuf, outfile);
    fputs(marker_allele_ptrs[marker_idx * 2], outfile);
    putc('\t', outfile);
    fputs(marker_allele_ptrs[marker_idx * 2 + 1], outfile);
    if (putc_checked('\n', outfile)) {
      goto lgen_to_bed_ret_WRITE_FAIL;
    }
    uii++;
    marker_idx++;
  }
  if (!feof(infile)) {
    goto lgen_to_bed_ret_READ_FAIL;
  }
  fclose_null(&infile);
  if (map_is_unsorted) {
    memcpy(outname_end, ".map.tmp", 9);
    unlink(outname);
  }
  if (fclose_null(&outfile)) {
    goto lgen_to_bed_ret_WRITE_FAIL;
  }
  memcpy(name_end, ".fam", 5);
  memcpy(outname_end, ".fam", 5);
#if _WIN32
  uii = GetFullPathName(lgen_namebuf, FNAMESIZE, tbuf, NULL);
  if ((!uii) || (uii > FNAMESIZE))
#else
  if (!realpath(lgen_namebuf, tbuf))
#endif
  {
    sprintf(logbuf, "Error: Failed to open %s.\n", outname);
    logprintb();
    goto lgen_to_bed_ret_OPEN_FAIL;
  }
#if _WIN32
  uii = GetFullPathName(outname, FNAMESIZE, &(tbuf[FNAMESIZE + 64]), NULL);
  if (!(uii && (uii <= FNAMESIZE) && (!strcmp(tbuf, &(tbuf[FNAMESIZE + 64])))))
#else
  cptr = realpath(outname, &(tbuf[FNAMESIZE + 64]));
  if (!(cptr && (!strcmp(tbuf, &(tbuf[FNAMESIZE + 64])))))
#endif
  {
    if (fopen_checked(&infile, lgen_namebuf, "r")) {
      goto lgen_to_bed_ret_OPEN_FAIL;
    }
    if (fopen_checked(&outfile, outname, "w")) {
      goto lgen_to_bed_ret_OPEN_FAIL;
    }
    while (fgets(tbuf, MAXLINELEN, infile)) {
      cptr = skip_initial_spaces(tbuf);
      if (is_eoln_kns(*cptr)) {
	continue;
      }
      ulii = strlen(cptr);
      if (cptr[ulii - 1] != '\n') {
	cptr[ulii++] = '\n';
      }
      if (fwrite_checked(cptr, ulii, outfile)) {
	goto lgen_to_bed_ret_WRITE_FAIL;
      }
    }
    if (!feof(infile)) {
      goto lgen_to_bed_ret_READ_FAIL;
    }
  }
  memcpy(outname_end, ".bed", 5);
  sprintf(logbuf, "Binary fileset written to %s + .bim + .fam.\n", outname);
  logprintb();

  while (0) {
  lgen_to_bed_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  lgen_to_bed_ret_CALC_NOT_YET_SUPPORTED:
    retval = RET_CALC_NOT_YET_SUPPORTED;
    break;
  lgen_to_bed_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  lgen_to_bed_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  lgen_to_bed_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  lgen_to_bed_ret_INVALID_FORMAT_2:
    cptr[strlen_se(cptr)] = '\0';
    cptr2[strlen_se(cptr2)] = '\0';
    sprintf(logbuf, "Error: Person %s %s in .lgen file but missing from .fam.\n", cptr, cptr2);
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  lgen_to_bed_ret_INVALID_FORMAT_5:
    logprint("Error: Half-missing call in .lgen file.\n");
    retval = RET_INVALID_FORMAT;
    break;
  lgen_to_bed_ret_INVALID_FORMAT_4:
    logprint("Error: Improperly formatted --reference file.\n");
    retval = RET_INVALID_FORMAT;
    break;
  lgen_to_bed_ret_INVALID_FORMAT:
    logprint("Error: Improperly formatted .lgen file.\n");
  lgen_to_bed_ret_INVALID_FORMAT_6:
    retval = RET_INVALID_FORMAT;
    break;
  lgen_to_bed_ret_INVALID_FORMAT_3:
    sprintf(logbuf, "Error: Variant %s in .lgen file has 3+ different alleles.\n", id_buf);
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  lgen_to_bed_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  }
 lgen_to_bed_ret_1:
  if (marker_allele_ptrs) {
    ma_end = &(marker_allele_ptrs[2 * marker_ct]);
    while (marker_allele_ptrs < ma_end) {
      sptr = *marker_allele_ptrs++;
      if (sptr && ((sptr < g_one_char_strs) || (sptr >= (&(g_one_char_strs[512]))))) {
	free(sptr);
      }
    }
  }
  wkspace_reset(wkspace_mark);
  if (infile) {
    fclose(infile);
  }
  if (outfile) {
    fclose(outfile);
  }
  return retval;
}

static inline uint32_t update_tped_alleles_and_cts(uint32_t* allele_tot_ptr, char** alleles, uint32_t* alens, uint32_t* allele_cts, char* ss, uint32_t slen, uint32_t* allele_idx_ptr) {
  uint32_t allele_idx;
  for (allele_idx = 0; allele_idx < (*allele_tot_ptr); allele_idx++) {
    if ((slen == alens[allele_idx]) && (!memcmp(alleles[allele_idx], ss, slen))) {
      allele_cts[allele_idx] += 1;
      *allele_idx_ptr = allele_idx;
      return 0;
    }
  }
  *allele_idx_ptr = allele_idx;
  if (allele_idx < 4) {
    alens[allele_idx] = slen;
    allele_cts[allele_idx] = 1;
    *allele_tot_ptr = allele_idx + 1;
    if (slen < 2) {
      alleles[allele_idx] = (char*)(&(g_one_char_strs[((unsigned char)(*ss)) * 2]));
    } else {
      alleles[allele_idx] = (char*)malloc(slen + 1);
      if (!alleles[allele_idx]) {
        return RET_NOMEM;
      }
      memcpyx(alleles[allele_idx], ss, slen, '\0');
    }
    return 0;
  } else {
    return RET_INVALID_FORMAT;
  }
}

void transposed_to_bed_print_pct(uint32_t pct) {
  printf("Processing .tped file... %u%%", pct);
  fflush(stdout);
}

int32_t transposed_to_bed(char* tpedname, char* tfamname, char* outname, char* outname_end, uint64_t misc_flags, Chrom_info* chrom_info_ptr) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* infile = NULL;
  FILE* bimfile = NULL;
  FILE* outfile = NULL;
  char** marker_allele_ptrs = NULL;
  uintptr_t topsize = 0;
  uintptr_t indiv_ct = 0;
  uint32_t no_extra_cols = 1;
  int32_t retval = 0;
  uint32_t pct = 0;
  uint32_t map_is_unsorted = 0;
  int64_t last_mapval = 0;
  uint32_t allele_tot = 0;
  uintptr_t marker_ct = 0;
  uintptr_t max_marker_id_len = 0;
  uintptr_t max_marker_allele_len = 2; // for .bim.tmp reloading
  const char* missing_geno_ptr = g_missing_geno_ptr;
  char missing_geno = *missing_geno_ptr;

  // We include handling of nonstandard chromosome names here since our own
  // pipeline requires this.  (To set them to zero, --make-bed must be used
  // in addition to this autoconverter.)
  uint32_t allow_extra_chroms = (misc_flags / MISC_ALLOW_EXTRA_CHROMS) & 1;

  // It's unnecessary for now, but we choose to include nice handling of
  // triallelic and quadallelic sites in this routine.
  char* alleles[4];
  char* salleles[4];
  uint32_t alens[4];
  uint32_t allele_cts[4];
  unsigned char writemap[17];
  uint32_t chrom_start[MAX_POSSIBLE_CHROM + 1];
  uint32_t chrom_id[MAX_POSSIBLE_CHROM];
  uintptr_t max_markers;
  uintptr_t indiv_ct4;
  uintptr_t indiv_idx;
  uintptr_t ulii;
  uint32_t uii;
  uint32_t ujj;
  int32_t ii;
  unsigned char* writebuf;
  unsigned char* prewritebuf;
  unsigned char ucc;
  uint32_t loadbuf_size;
  char* allele_buf;
  char* loadbuf;
  char* cptr;
  char* cptr2;
  char* cptr3;
  char* cptr4;
  char* axptr;
  uint32_t axlen;
  unsigned char* ucptr;
  unsigned char* ucptr2;
  int64_t tped_size;
  int64_t tped_next_thresh;
  int64_t cur_mapval;
  int64_t* mapvals;
  uintptr_t marker_idx;
  uintptr_t marker_uidx;
  uint32_t* map_reverse;
  int64_t* ll_buf;
  uint32_t* pos_buf;
  char* marker_ids;
  uint32_t cur_chrom;
  uint32_t chrom_ct;
  double* marker_cms;

  logstr("Processing .tped file.\n");
  transposed_to_bed_print_pct(0);
  fflush(stdout);
  if (fopen_checked(&infile, tfamname, "r")) {
    goto transposed_to_bed_ret_OPEN_FAIL;
  }
  memcpy(outname_end, ".fam", 5);
  if (fopen_checked(&outfile, outname, "w")) {
    goto transposed_to_bed_ret_OPEN_FAIL;
  }
  tbuf[MAXLINELEN - 1] = ' ';
  while (fgets(tbuf, MAXLINELEN, infile)) {
    if (!tbuf[MAXLINELEN - 1]) {
      logprint("Error: Pathologically long line in .tfam file.\n");
      goto transposed_to_bed_ret_INVALID_FORMAT_2;
    }
    cptr = skip_initial_spaces(tbuf);
    if (is_eoln_kns(*cptr)) {
      continue;
    }
    ulii = strlen(cptr);
    if (cptr[ulii - 1] != '\n') {
      cptr[ulii++] = '\n';
    }
    if (fwrite_checked(cptr, ulii, outfile)) {
      goto transposed_to_bed_ret_WRITE_FAIL;
    }
    indiv_ct++;
  }
  if (!feof(infile)) {
    goto transposed_to_bed_ret_READ_FAIL;
  }
  if (!indiv_ct) {
    sprintf(logbuf, "Error: No %s in .tfam file.\n", g_species_plural);
    goto transposed_to_bed_ret_INVALID_FORMAT_5;
  }
  indiv_ct4 = (indiv_ct + 3) / 4;
  fclose_null(&infile);
  fclose_null(&outfile);

  memcpy(outname_end, ".bim.tmp", 9);
  if (fopen_checked(&bimfile, outname, "w")) {
    goto transposed_to_bed_ret_OPEN_FAIL;
  }
  memcpy(outname_end, ".bed.tmp", 9);
  if (fopen_checked(&outfile, outname, "wb")) {
    goto transposed_to_bed_ret_OPEN_FAIL;
  }
  if (wkspace_alloc_uc_checked(&writebuf, indiv_ct4) ||
      wkspace_alloc_uc_checked(&prewritebuf, indiv_ct)) {
    goto transposed_to_bed_ret_NOMEM;
  }
  // long allele names are allocated outside workspace anyway, so it makes
  // sense for max allele length to be related to reserved non-workspace memory
  allele_buf = (char*)top_alloc(&topsize, NON_WKSPACE_MIN);
  if (!allele_buf) {
    goto transposed_to_bed_ret_NOMEM;
  }
  max_markers = (wkspace_left - topsize) / sizeof(int64_t);
  mapvals = (int64_t*)wkspace_base;
  writemap[16] = 1;
  if (fwrite_checked("l\x1b\x01", 3, outfile)) {
    goto transposed_to_bed_ret_WRITE_FAIL;
  }

  // given e.g. 6MB indels in real datasets, there's legitimate reason for a
  // .tped line to be even longer than 2GB, so we use ftoken_...() over
  // fgets().
  if (fopen_checked(&infile, tpedname, "r")) {
    goto transposed_to_bed_ret_OPEN_FAIL;
  }
  if (fseeko(infile, 0, SEEK_END)) {
    goto transposed_to_bed_ret_READ_FAIL;
  }
  tped_size = ftello(infile);
  rewind(infile);
  tped_next_thresh = tped_size / 100;

  while (1) {
    tbuf[MAXLINELEN - 1] = ' ';
    if (!fgets(tbuf, MAXLINELEN, infile)) {
      break;
    }
    // assume first four fields are within MAXLINELEN characters, but after
    // that, anything goes
    cptr = skip_initial_spaces(tbuf);
    if (is_eoln_kns(*cptr)) {
      if (!tbuf[MAXLINELEN - 1]) {
        goto transposed_to_bed_ret_INVALID_FORMAT;
      }
      continue;
    }
    cptr2 = next_item(cptr);
    cptr3 = next_item_mult(cptr2, 2);
    cptr4 = next_item(cptr3);
    if (no_more_items_kns(cptr4)) {
      goto transposed_to_bed_ret_INVALID_FORMAT;
    }
    if (ftello(infile) >= tped_next_thresh) {
      uii = (ftello(infile) * 100) / tped_size;
      if (pct >= 10) {
	putchar('\b');
      }
      printf("\b\b%u%%", uii);
      fflush(stdout);
      pct = uii;
      tped_next_thresh = ((pct + 1) * tped_size) / 100;
    }
    ii = get_chrom_code(chrom_info_ptr, cptr);
    if (ii == -1) {
      if (!allow_extra_chroms) {
	logprint("Error: Unrecognized chromosome code in .tped file.  (Did you forget\n--allow-extra-chroms?).\n");
	goto transposed_to_bed_ret_INVALID_FORMAT_2;
      }
      retval = resolve_or_add_chrom_name(chrom_info_ptr, cptr, &ii);
      if (retval) {
	goto transposed_to_bed_ret_1;
      }
    }

    if ((*cptr3 == '-') || (!is_set(chrom_info_ptr->chrom_mask, ii))) {
      cptr2 = cptr4;
      goto transposed_to_bed_nextline;
    }
    uii = strlen_se(cptr2);
    if (uii >= max_marker_id_len) {
      max_marker_id_len = uii + 1;
    }
    cur_mapval = (int64_t)((((uint64_t)((uint32_t)ii)) << 32) | ((uint32_t)atoi(cptr3)));
    if (marker_ct == max_markers) {
      goto transposed_to_bed_ret_NOMEM;
    }
    mapvals[marker_ct++] = cur_mapval;
    if (last_mapval > cur_mapval) {
      map_is_unsorted = 1;
    } else {
      last_mapval = cur_mapval;
    }
    for (uii = 0; uii < 3; uii++) {
      cptr2 = item_endnn(cptr);
      *cptr2++ = '\t';
      fwrite(cptr, 1, cptr2 - cptr, bimfile);
      cptr = skip_initial_spaces(cptr2);
    }
    cptr2 = item_endnn(cptr);
    *cptr2++ = '\t';
    if (fwrite_checked(cptr, cptr2 - cptr, bimfile)) {
      goto transposed_to_bed_ret_WRITE_FAIL;
    }
    cptr2 = cptr4;
    alleles[0] = NULL;
    alleles[1] = NULL;
    alleles[2] = NULL;
    alleles[3] = NULL;
    allele_cts[0] = 0;
    allele_cts[1] = 0;
    allele_cts[2] = 0;
    allele_cts[3] = 0;
    for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
      cptr2 = skip_initial_spaces(cptr2);
      while (cptr2 == &(tbuf[MAXLINELEN - 1])) {
	if (cptr2[-1] == '\n') {
	  goto transposed_to_bed_ret_INVALID_FORMAT;
	}
        if (!fgets(tbuf, MAXLINELEN, infile)) {
          if (ferror(infile)) {
	    goto transposed_to_bed_ret_READ_FAIL;
	  }
	  goto transposed_to_bed_ret_INVALID_FORMAT;
	}
	cptr2 = skip_initial_spaces(tbuf);
      }
      axptr = cptr2;
      axlen = strlen_se(cptr2);
      cptr2 = &(axptr[axlen]);
      // only way for this to happen if it isn't at end of buffer is if we're
      // at EOF, which is an error anyway
      if (!(*cptr2)) {
	if (!axlen) {
	  goto transposed_to_bed_ret_INVALID_FORMAT;
	}
	cptr3 = memcpya(allele_buf, axptr, axlen);
        axptr = allele_buf;
	do {
	  if (!fgets(tbuf, MAXLINELEN, infile)) {
	    if (ferror(infile)) {
	      goto transposed_to_bed_ret_READ_FAIL;
	    }
            goto transposed_to_bed_ret_INVALID_FORMAT;
	  }
	  cptr2 = tbuf;
          if (!is_space_or_eoln(*cptr2)) {
	    cptr2 = item_endnn(cptr2);
	  }
	  if ((((uintptr_t)(cptr3 - allele_buf)) + ((uintptr_t)(cptr2 - tbuf))) >= NON_WKSPACE_MIN) {
	    goto transposed_to_bed_ret_NOMEM;
	  }
	  cptr3 = memcpya(cptr3, tbuf, cptr2 - tbuf);
	} while (!(*cptr2));
	axlen = (uintptr_t)(cptr3 - allele_buf);
      }
      if ((*axptr != missing_geno) || (axlen != 1)) {
	retval = update_tped_alleles_and_cts(&allele_tot, alleles, alens, allele_cts, axptr, axlen, &uii);
	if (retval) {
	  if (retval == RET_INVALID_FORMAT) {
	    goto transposed_to_bed_ret_INVALID_FORMAT_6;
	  }
	  goto transposed_to_bed_ret_NOMEM;
	}
      } else {
	uii = 4;
      }
      cptr2 = skip_initial_spaces(cptr2);
      while (cptr2 == &(tbuf[MAXLINELEN - 1])) {
	if (cptr2[-1] == '\n') {
	  goto transposed_to_bed_ret_INVALID_FORMAT;
	}
        if (!fgets(tbuf, MAXLINELEN, infile)) {
          if (ferror(infile)) {
	    goto transposed_to_bed_ret_READ_FAIL;
	  }
	  goto transposed_to_bed_ret_INVALID_FORMAT;
	}
	cptr2 = skip_initial_spaces(tbuf);
      }
      axptr = cptr2;
      axlen = strlen_se(cptr2);
      cptr2 = &(axptr[axlen]);
      if (!(*cptr2)) {
	if (!axlen) {
	  goto transposed_to_bed_ret_INVALID_FORMAT;
	}
	cptr3 = memcpya(allele_buf, axptr, axlen);
        axptr = allele_buf;
	do {
	  cptr2 = tbuf;
	  if (!fgets(tbuf, MAXLINELEN, infile)) {
	    if (ferror(infile)) {
	      goto transposed_to_bed_ret_READ_FAIL;
	    } else if (indiv_idx != indiv_ct - 1) {
              goto transposed_to_bed_ret_INVALID_FORMAT;
	    } else {
	      tbuf[0] = '\0';
	      break;
	    }
	  }
          if (!is_space_or_eoln(*cptr2)) {
	    cptr2 = item_endnn(cptr2);
	  }
	  if ((((uintptr_t)(cptr3 - allele_buf)) + ((uintptr_t)(cptr2 - tbuf))) >= NON_WKSPACE_MIN) {
	    goto transposed_to_bed_ret_NOMEM;
	  }
	  cptr3 = memcpya(cptr3, tbuf, cptr2 - tbuf);
	} while (!(*cptr2));
	axlen = (uintptr_t)(cptr3 - allele_buf);
      }
      if ((*axptr != missing_geno) || (axlen != 1)) {
	if (uii == 4) {
	  goto transposed_to_bed_ret_INVALID_FORMAT_4;
	}
	retval = update_tped_alleles_and_cts(&allele_tot, alleles, alens, allele_cts, axptr, axlen, &ujj);
	if (retval) {
	  if (retval == RET_INVALID_FORMAT) {
	    goto transposed_to_bed_ret_INVALID_FORMAT_6;
	  }
	  goto transposed_to_bed_ret_NOMEM;
	}
        prewritebuf[indiv_idx] = uii * 4 + ujj;
      } else {
	if (uii != 4) {
	  goto transposed_to_bed_ret_INVALID_FORMAT_4;
	}
	prewritebuf[indiv_idx] = 16;
      }
    }

    memcpy(salleles, alleles, 4 * sizeof(intptr_t));
    for (uii = 1; uii < 4; uii++) {
      ujj = allele_cts[uii];
      if (allele_cts[uii - 1] < ujj) {
	axptr = salleles[uii];
	ii = uii;
	do {
	  ii--;
	  salleles[((uint32_t)ii) + 1] = salleles[(uint32_t)ii];
	  allele_cts[((uint32_t)ii) + 1] = allele_cts[(uint32_t)ii];
	} while (ii && (allele_cts[((uint32_t)ii) - 1] < ujj));
	salleles[(uint32_t)ii] = axptr;
	allele_cts[(uint32_t)ii] = ujj;
      }
    }
    if (allele_cts[2]) {
      putchar('\r');
      sprintf(logbuf, "Note: Variant %" PRIuPTR " is %sallelic.  Setting rarest alleles to missing.\n", marker_ct - 1, allele_cts[3]? "quad" : "tri");
      logprintb();
      transposed_to_bed_print_pct(pct);
    }
    for (uii = 0; uii < 4; uii++) {
      axptr = alleles[uii];
      ucptr = &(writemap[4 * uii]);
      if (!axptr) {
	memset(ucptr, 1, 4);
      } else if (axptr == salleles[0]) {
        for (ujj = 0; ujj < 4; ujj++) {
	  axptr = alleles[ujj];
	  if (!axptr) {
	    *ucptr++ = 1;
	  } else if (axptr == salleles[0]) {
	    *ucptr++ = 3;
	  } else if (axptr == salleles[1]) {
	    *ucptr++ = 2;
	  } else {
	    *ucptr++ = 1;
	  }
	}
      } else if (axptr == salleles[1]) {
	for (ujj = 0; ujj < 4; ujj++) {
	  axptr = alleles[ujj];
	  if (!axptr) {
	    *ucptr++ = 1;
	  } else if (axptr == salleles[0]) {
	    *ucptr++ = 2;
	  } else if (axptr == salleles[1]) {
	    *ucptr++ = 0;
	  } else {
	    *ucptr++ = 1;
	  }
	}
      } else {
        memset(ucptr, 1, 4);
      }
    }
    uii = indiv_ct & (~3U);
    ucptr = writebuf;
    for (ujj = 0; ujj < uii; ujj += 4) {
      *ucptr++ = writemap[prewritebuf[ujj]] | (writemap[prewritebuf[ujj + 1]] << 2) | (writemap[prewritebuf[ujj + 2]] << 4) | (writemap[prewritebuf[ujj + 3]] << 6);
    }
    ucc = 0;
    ucptr2 = &(prewritebuf[uii]);
    uii = indiv_ct % 4;
    if (uii) {
      for (ujj = 0; ujj < uii; ujj++) {
        ucc |= (writemap[*ucptr2++]) << (ujj * 2);
      }
      *ucptr = ucc;
    }
    fwrite(writebuf, 1, indiv_ct4, outfile);
    if (!salleles[1]) {
      putc(missing_geno, bimfile);
    } else {
      uii = strlen(salleles[1]);
      if (uii >= max_marker_allele_len) {
	max_marker_allele_len = uii + 1;
      }
      fputs(salleles[1], bimfile);
    }
    putc('\t', bimfile);
    if (!salleles[0]) {
      putc(missing_geno, bimfile);
    } else {
      uii = strlen(salleles[0]);
      if (uii >= max_marker_allele_len) {
	max_marker_allele_len = uii + 1;
      }
      fputs(salleles[0], bimfile);
    }
    for (uii = 0; uii < allele_tot; uii++) {
      if (alleles[uii][1]) {
	free(alleles[uii]);
      }
    }
    allele_tot = 0;
    if (putc_checked('\n', bimfile)) {
      goto transposed_to_bed_ret_WRITE_FAIL;
    }
    if (no_extra_cols) {
      cptr2 = skip_initial_spaces(cptr2);
      while (cptr2 == &(tbuf[MAXLINELEN - 1])) {
	if (cptr2[-1] == '\n') {
	  break;
	}
	cptr2 = tbuf;
	if (!fgets(tbuf, MAXLINELEN, infile)) {
	  if (ferror(infile)) {
	    goto transposed_to_bed_ret_READ_FAIL;
	  }
	  tbuf[0] = '\0';
	  break;
	}
        cptr2 = skip_initial_spaces(cptr2);
      }
      if (!is_space_or_eoln(*cptr2)) {
	no_extra_cols = 0;
	putchar('\r');
	logprint("Note: Extra columns in .tped file.  Ignoring.\n");
	transposed_to_bed_print_pct(pct);
	goto transposed_to_bed_nextline;
      }
    } else {
    transposed_to_bed_nextline:
      cptr2 = (char*)memchr(cptr2, 0, MAXLINELEN - ((uintptr_t)(cptr2 - tbuf)));
      while (cptr2 == &(tbuf[MAXLINELEN - 1])) {
	if (cptr2[-1] == '\n') {
	  break;
	}
        if (!fgets(tbuf, MAXLINELEN, infile)) {
          if (ferror(infile)) {
	    goto transposed_to_bed_ret_READ_FAIL;
	  }
          break;
	}
	cptr2 = (char*)memchr(tbuf, 0, MAXLINELEN);
      }
    }
  }
  // topsize = 0;
  if (fclose_null(&infile)) {
    goto transposed_to_bed_ret_READ_FAIL;
  }
  if (fclose_null(&bimfile)) {
    goto transposed_to_bed_ret_WRITE_FAIL;
  }
  if (fclose_null(&outfile)) {
    goto transposed_to_bed_ret_WRITE_FAIL;
  }

  if (map_is_unsorted) {
    loadbuf_size = 2 * max_marker_allele_len + MAXLINELEN;

    // mapvals is already positioned here, this just updates wkspace_base and
    // wkspace_left
    mapvals = (int64_t*)wkspace_alloc(marker_ct * sizeof(int64_t));

    if (wkspace_alloc_ll_checked(&ll_buf, marker_ct * sizeof(int64_t)) ||
        wkspace_alloc_ui_checked(&pos_buf, marker_ct * sizeof(int32_t)) ||
        wkspace_alloc_c_checked(&marker_ids, marker_ct * max_marker_id_len) ||
	wkspace_alloc_d_checked(&marker_cms, marker_ct * sizeof(double)) ||
        wkspace_alloc_c_checked(&loadbuf, loadbuf_size)) {
      goto transposed_to_bed_ret_NOMEM;
    }
    marker_allele_ptrs = (char**)wkspace_alloc(marker_ct * 2 * sizeof(intptr_t));
    if (!marker_allele_ptrs) {
      goto transposed_to_bed_ret_NOMEM;
    }
    // prevent cleanup from failing
    memset(marker_allele_ptrs, 0, marker_ct * 2 * sizeof(intptr_t));

    for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
      pos_buf[marker_idx] = (uint32_t)((uint64_t)mapvals[marker_idx]);
      ll_buf[marker_idx] = (mapvals[marker_idx] & 0xffffffff00000000LLU) | marker_idx;
    }
    sort_marker_chrom_pos(ll_buf, marker_ct, pos_buf, chrom_start, chrom_id, &chrom_ct);

    memcpy(outname_end, ".bim.tmp", 9);
    if (fopen_checked(&infile, outname, "r")) {
      goto transposed_to_bed_ret_OPEN_FAIL;
    }
    outname_end[4] = '\0';
    if (fopen_checked(&outfile, outname, "w")) {
      goto transposed_to_bed_ret_OPEN_FAIL;
    }
    marker_idx = 0;
    while (fgets(loadbuf, loadbuf_size, infile)) {
      // .tmp file, guaranteed to be no spaces in front
      cptr = next_item(loadbuf);
      cptr2 = item_endl(cptr);
      if (!cptr2) {
	goto transposed_to_bed_ret_INVALID_FORMAT;
      }
      cptr3 = skip_initial_spaces(cptr2);
      cptr4 = next_item_mult(cptr3, 2);
      uii = cptr2 - cptr;
      memcpyx(&(marker_ids[marker_idx * max_marker_id_len]), cptr, uii, '\0');
      if (scan_double(cptr3, &(marker_cms[marker_idx]))) {
	goto transposed_to_bed_ret_INVALID_FORMAT;
      }
      uii = strlen_se(cptr4);
      if (allele_set(&(marker_allele_ptrs[2 * marker_idx]), cptr4, uii)) {
	goto transposed_to_bed_ret_NOMEM;
      }
      cptr4 = skip_initial_spaces(&(cptr4[uii + 1]));
      uii = strlen_se(cptr4);
      if (allele_set(&(marker_allele_ptrs[2 * marker_idx + 1]), cptr4, uii)) {
	goto transposed_to_bed_ret_NOMEM;
      }
      marker_idx++;
    }
    if (!feof(infile)) {
      goto transposed_to_bed_ret_READ_FAIL;
    }
    fclose_null(&infile);
    marker_idx = 0;
    map_reverse = (uint32_t*)mapvals;
    for (uii = 0; uii < chrom_ct; uii++) {
      cur_chrom = chrom_id[uii];
      ujj = chrom_start[uii + 1];
      cptr2 = chrom_name_write(&(tbuf[MAXLINELEN]), chrom_info_ptr, cur_chrom, 0);
      *cptr2++ = '\t';
      for (; marker_idx < ujj; marker_idx++) {
	marker_uidx = (uint32_t)ll_buf[marker_idx];
	fwrite(&(tbuf[MAXLINELEN]), 1, cptr2 - (&(tbuf[MAXLINELEN])), outfile);
	fputs(&(marker_ids[marker_uidx * max_marker_id_len]), outfile);
	tbuf[0] = '\t';
	cptr = uint32_writex(double_g_writex(&(tbuf[1]), marker_cms[marker_uidx], '\t'), (uint32_t)(ll_buf[marker_idx] >> 32), '\t');
	if (fwrite_checked(tbuf, (uintptr_t)(cptr - tbuf), outfile)) {
	  goto transposed_to_bed_ret_WRITE_FAIL;
	}
        fputs(marker_allele_ptrs[2 * marker_uidx], outfile);
        putc('\t', outfile);
        fputs(marker_allele_ptrs[2 * marker_uidx + 1], outfile);
        if (putc_checked('\n', outfile)) {
	  goto transposed_to_bed_ret_WRITE_FAIL;
	}
	map_reverse[marker_uidx] = marker_idx;
      }
    }
    if (fclose_null(&outfile)) {
      goto transposed_to_bed_ret_WRITE_FAIL;
    }

    wkspace_reset((unsigned char*)map_reverse);
    map_reverse = (uint32_t*)wkspace_alloc(marker_ct * sizeof(int32_t));
    outname_end[4] = '.';
    unlink(outname);

    outname_end[2] = 'e';
    outname_end[3] = 'd';
    if (fopen_checked(&infile, outname, "rb")) {
      goto transposed_to_bed_ret_OPEN_FAIL;
    }
    outname_end[4] = '\0';
    if (fopen_checked(&outfile, outname, "wb")) {
      goto transposed_to_bed_ret_OPEN_FAIL;
    }
    if (fwrite_checked("l\x1b\x01", 3, outfile)) {
      goto transposed_to_bed_ret_WRITE_FAIL;
    }
    uii = 0xfffffffeU; // last marker uidx
    for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
      marker_uidx = map_reverse[marker_idx];
      if (marker_uidx != uii + 1) {
        if (fseeko(infile, 3 + ((uint64_t)marker_uidx) * indiv_ct4, SEEK_SET)) {
	  goto transposed_to_bed_ret_READ_FAIL;
	}
      }
      if (fread(writebuf, 1, indiv_ct4, infile) < indiv_ct4) {
	goto transposed_to_bed_ret_READ_FAIL;
      }
      if (fwrite_checked(writebuf, indiv_ct4, outfile)) {
	goto transposed_to_bed_ret_WRITE_FAIL;
      }
      uii = marker_uidx;
    }
    fclose_null(&infile);
    outname_end[4] = '.';
    unlink(outname);
    outname_end[4] = '\0';
  } else {
    uii = (outname_end - outname);
    memcpy(outname_end, ".bim.tmp", 9);
    memcpy(tbuf, outname, 9 + uii);
    outname_end[4] = '\0';
    if (rename(tbuf, outname)) {
      goto transposed_to_bed_ret_WRITE_FAIL;
    }
    tbuf[uii + 2] = 'e';
    tbuf[uii + 3] = 'd';
    outname_end[2] = 'e';
    outname_end[3] = 'd';
    if (rename(tbuf, outname)) {
      goto transposed_to_bed_ret_WRITE_FAIL;
    }
  }
  fputs("\rProcessing .tped file... done.\n", stdout);
  sprintf(logbuf, "%s + .bim + .fam written.\n", outname);
  logprintb();

  while (0) {
  transposed_to_bed_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  transposed_to_bed_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  transposed_to_bed_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  transposed_to_bed_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  transposed_to_bed_ret_INVALID_FORMAT:
    putchar('\r');
    logprint("Error: Improperly formatted .tped file.\n");
  transposed_to_bed_ret_INVALID_FORMAT_2:
    retval = RET_INVALID_FORMAT;
    break;
  transposed_to_bed_ret_INVALID_FORMAT_4:
    sprintf(logbuf, "Error: Half-missing call at variant %" PRIuPTR ", indiv %" PRIuPTR " in .tped file.\n", marker_ct - 1, indiv_idx);
  transposed_to_bed_ret_INVALID_FORMAT_5:
    putchar('\r');
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  transposed_to_bed_ret_INVALID_FORMAT_6:
    putchar('\r');
    sprintf(logbuf, "Error: More than four alleles at variant %" PRIuPTR ".\n", marker_ct - 1);
    logprintb();
    break;
  }
 transposed_to_bed_ret_1:
  for (uii = 0; uii < allele_tot; uii++) {
    if (alleles[uii][1]) {
      free(alleles[uii]);
    }
  }
  if (marker_allele_ptrs && (max_marker_allele_len > 2)) {
    for (marker_idx = 0; marker_idx < marker_ct * 2; marker_idx++) {
      cptr = marker_allele_ptrs[marker_idx];
      if (cptr && ((cptr < g_one_char_strs) || (cptr >= (&(g_one_char_strs[512]))))) {
	free(cptr);
      }
    }
  }
  fclose_cond(infile);
  fclose_cond(bimfile);
  fclose_cond(outfile);
  wkspace_reset(wkspace_mark);
  return retval;
}

// oh, what the hell, may as well be liberal
#define MAX_VCF_ALT 65534

int32_t vcf_to_bed(char* vcfname, char* outname, char* outname_end, int32_t missing_pheno, uint64_t misc_flags, char* const_fid, char id_delim, double vcf_min_qual, char* vcf_filter_exceptions_flattened, Chrom_info* chrom_info_ptr) {
  unsigned char* wkspace_mark = wkspace_base;
  gzFile gz_infile = NULL;
  FILE* outfile = NULL;
  FILE* bimfile = NULL;
  FILE* skip3file = NULL;
  char* sorted_fexcepts = NULL;
  uint32_t double_id = (misc_flags / MISC_DOUBLE_ID) & 1;
  uint32_t check_qual = (vcf_min_qual == -INFINITY)? 0 : 1;
  uint32_t allow_extra_chroms = (misc_flags / MISC_ALLOW_EXTRA_CHROMS) & 1;
  uint32_t biallelic_only = (misc_flags / MISC_BIALLELIC_ONLY) & 1;
  uint32_t biallelic_strict = (misc_flags / MISC_BIALLELIC_ONLY_STRICT) & 1;
  uint32_t skip3_list = (misc_flags / MISC_BIALLELIC_ONLY_LIST) & 1;
  uint32_t marker_ct = 0;
  uint32_t marker_skip_ct = 0;
  uintptr_t fexcept_ct = 0;
  uintptr_t max_fexcept_len = 5;
  uintptr_t const_fid_len = 0;
  uintptr_t indiv_ct = 0;
  char missing_geno = *g_missing_geno_ptr;
  int32_t retval = 0;
  char fam_trailer[20];
  uint32_t* vcf_alt_cts;
  char* loadbuf;
  char* bufptr;
  char* bufptr2;
  char* bufptr3;
  char* chrom_ptr;
  char* marker_id;
  char* pos_str;
  char* alt_alleles;
  char* wptr;
  uintptr_t* base_bitfields;
  uintptr_t* alt_bitfield;
  uintptr_t* ref_ptr;
  uintptr_t* alt_ptr;
  uintptr_t final_mask;
  uintptr_t indiv_ctl2;
  uintptr_t indiv_ctv2;
  uintptr_t indiv_ct4;
  uintptr_t indiv_idx;
  uintptr_t loadbuf_size;
  uintptr_t fam_trailer_len;
  uintptr_t slen;
  uintptr_t alt_idx;
  uintptr_t alt_ct;
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t ulkk;
  uintptr_t alt_allele_idx;
  double dxx;
  uint32_t chrom_len;
  uint32_t marker_id_len;
  int32_t ii;
  char cc;
  if (gzopen_checked(&gz_infile, vcfname, "rb")) {
    goto vcf_to_bed_ret_OPEN_FAIL;
  }
  if (gzbuffer(gz_infile, 131072)) {
    goto vcf_to_bed_ret_NOMEM;
  }
  if (misc_flags & MISC_VCF_FILTER) {
    // automatically include "." and "PASS"
    fexcept_ct = 2;
    if (vcf_filter_exceptions_flattened) {
      bufptr = vcf_filter_exceptions_flattened;
      do {
        fexcept_ct++;
        slen = strlen(bufptr);
        if (slen >= max_fexcept_len) {
          max_fexcept_len = slen + 1;
	}
	bufptr = &(bufptr[slen + 1]);
      } while (*bufptr);
    }
    if (wkspace_alloc_c_checked(&sorted_fexcepts, fexcept_ct * max_fexcept_len)) {
      goto vcf_to_bed_ret_NOMEM;
    }
    memcpy(sorted_fexcepts, ".", 2);
    memcpy(&(sorted_fexcepts[max_fexcept_len]), "PASS", 5);
    if (vcf_filter_exceptions_flattened) {
      bufptr = vcf_filter_exceptions_flattened;
      for (ulii = 2; ulii < fexcept_ct; ulii++) {
        slen = strlen(bufptr) + 1;
        memcpy(&(sorted_fexcepts[ulii * max_fexcept_len]), bufptr, slen);
        bufptr = &(bufptr[slen]);
      }
      qsort(sorted_fexcepts, fexcept_ct, max_fexcept_len, strcmp_casted);
      fexcept_ct = collapse_duplicate_ids(sorted_fexcepts, fexcept_ct, max_fexcept_len, NULL);
      // there can't be many filter exceptions, so don't bother to free unused
      // memory in corner case
    }
  }
  bufptr = memcpya(fam_trailer, "\t0\t0\t0\t", 7);
  bufptr = int32_writex(bufptr, missing_pheno, '\n');
  fam_trailer_len = (uintptr_t)(bufptr - fam_trailer);

  loadbuf_size = wkspace_left;
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  } else if (loadbuf_size <= MAXLINELEN) {
    goto vcf_to_bed_ret_NOMEM;
  }
  
  loadbuf = (char*)wkspace_base;
  loadbuf[loadbuf_size - 1] = ' ';
  while (1) {
    if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
      goto vcf_to_bed_ret_READ_FAIL;
    }
    if (!loadbuf[loadbuf_size - 1]) {
      if (loadbuf_size == MAXLINEBUFLEN) {
        logprint("Error: Pathologically long line in .vcf file.\n");
        goto vcf_to_bed_ret_INVALID_FORMAT;
      }
      goto vcf_to_bed_ret_NOMEM;
    }
    bufptr = skip_initial_spaces(loadbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    if (*bufptr != '#') {
      logprint("Error: Missing header line in .vcf file.\n");
      goto vcf_to_bed_ret_INVALID_FORMAT;
    }
    if (bufptr[1] != '#') {
      break;
    }
  }
  if (memcmp(bufptr, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", 38)) {
    logprint("Error: Improperly formatted .vcf header line.\n");
    goto vcf_to_bed_ret_INVALID_FORMAT;
  }
  bufptr = &(bufptr[38]);
  if (memcmp(bufptr, "\tFORMAT\t", 8) || (((unsigned char)bufptr[8]) <= ' ')) {
    logprint("Error: No genotype data in .vcf file.\n");
    goto vcf_to_bed_ret_INVALID_FORMAT;
  }
  bufptr = &(bufptr[8]);
  memcpy(outname_end, ".fam", 5);
  if (fopen_checked(&outfile, outname, "w")) {
    goto vcf_to_bed_ret_OPEN_FAIL;
  }
  if (const_fid) {
    const_fid_len = strlen(const_fid);
  } else if ((!double_id) && (!id_delim)) {
    // default: --double-id + --id-delim
    double_id = 1;
    id_delim = '_';
  }
  do {
    indiv_ct++;
    bufptr2 = strchr(bufptr, '\t');
    if (bufptr2) {
      slen = (uintptr_t)(bufptr2 - bufptr);
    } else {
      slen = strlen_se(bufptr);
      bufptr2 = &(bufptr[slen]);
    }
    if (slen > 65535) {
      logprint("Error: --vcf does not support sample IDs longer than 65535 characters.\n");
      goto vcf_to_bed_ret_INVALID_FORMAT;
    }
    if ((*bufptr == '0') && (slen == 1)) {
      logprint("Error: Sample ID cannot be '0'.\n");
      goto vcf_to_bed_ret_INVALID_FORMAT;
    }
    if (id_delim) {
      if (*bufptr == id_delim) {
	sprintf(logbuf, "Error: '%c' at beginning of sample ID.\n", id_delim);
	goto vcf_to_bed_ret_INVALID_FORMAT_2;
      } else if (bufptr[slen - 1] == id_delim) {
	sprintf(logbuf, "Error: '%c' at end of sample ID.\n", id_delim);
	goto vcf_to_bed_ret_INVALID_FORMAT_2;
      }
      bufptr3 = (char*)memchr(bufptr, (unsigned char)id_delim, slen);
      if (!bufptr3) {
	if (double_id) {
	  goto vcf_to_bed_double_id;
	} else if (const_fid) {
	  goto vcf_to_bed_const_id;
	} else {
	  sprintf(logbuf, "Error: No '%c' in sample ID.\n", id_delim);
	  goto vcf_to_bed_ret_INVALID_FORMAT_2;
	}
      }
      if (memchr(&(bufptr3[1]), (unsigned char)id_delim, (uintptr_t)(bufptr2 - &(bufptr3[1])))) {
        sprintf(logbuf, "Error: Multiple instances of '%c' in sample ID.\n", id_delim);
        goto vcf_to_bed_ret_INVALID_FORMAT_2;
      }
      wptr = memcpyax(tbuf, bufptr, (uintptr_t)(bufptr3 - bufptr), '\t');
      bufptr3++;
      if ((*bufptr3 == '0') && (bufptr2 == &(bufptr3[1]))) {
        sprintf(logbuf, "Error: Sample ID ends with \"%c0\", which induces an invalid IID of '0'.\n", id_delim);
        goto vcf_to_bed_ret_INVALID_FORMAT_2;
      }
      wptr = memcpya(wptr, bufptr3, (uintptr_t)(bufptr2 - bufptr3));
    } else {
      if (double_id) {
      vcf_to_bed_double_id:
	wptr = memcpyax(tbuf, bufptr, (uintptr_t)(bufptr2 - bufptr), '\t');
      } else {
      vcf_to_bed_const_id:
        wptr = memcpyax(tbuf, const_fid, const_fid_len, '\t');
      }
      wptr = memcpya(wptr, bufptr, (uintptr_t)(bufptr2 - bufptr));
    }
    wptr = memcpya(wptr, fam_trailer, fam_trailer_len);
    if (fwrite_checked(tbuf, wptr - tbuf, outfile)) {
      goto vcf_to_bed_ret_WRITE_FAIL;
    }
    if (*bufptr2 != '\t') {
      break;
    }
    bufptr = &(bufptr2[1]);
  } while (((unsigned char)bufptr[0]) > ' ');
  if (fclose_null(&outfile)) {
    goto vcf_to_bed_ret_WRITE_FAIL;
  }

  indiv_ct4 = (indiv_ct + 3) / 4;
  indiv_ctl2 = (indiv_ct + BITCT2 - 1) / BITCT2;
  indiv_ctv2 = 2 * ((indiv_ct + BITCT - 1) / BITCT);
  final_mask = (~ZEROLU) >> (2 * ((0x7fffffe0 - indiv_ct) % BITCT2));
  if (wkspace_alloc_ul_checked(&base_bitfields, indiv_ctv2 * 10 * sizeof(intptr_t)) ||
      wkspace_alloc_ui_checked(&vcf_alt_cts, MAX_VCF_ALT * sizeof(int32_t))) {
    goto vcf_to_bed_ret_NOMEM;
  }
  memcpy(outname_end, ".bim", 5);
  if (fopen_checked(&bimfile, outname, "w")) {
    goto vcf_to_bed_ret_OPEN_FAIL;
  }
  memcpyl3(&(outname_end[2]), "ed");
  if (fopen_checked(&outfile, outname, "wb")) {
    goto vcf_to_bed_ret_OPEN_FAIL;
  }
  if (fwrite_checked("l\x1b\x01", 3, outfile)) {
    goto vcf_to_bed_ret_WRITE_FAIL;
  }
  loadbuf_size = wkspace_left;
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  } else if (loadbuf_size <= MAXLINELEN) {
    goto vcf_to_bed_ret_NOMEM;
  }
  
  loadbuf = (char*)wkspace_base;
  loadbuf[loadbuf_size - 1] = ' ';
  while (1) {
    if (!gzgets(gz_infile, loadbuf, loadbuf_size)) {
      if (!gzeof(gz_infile)) {
        goto vcf_to_bed_ret_READ_FAIL;
      }
      break;
    }
    if (!loadbuf[loadbuf_size - 1]) {
      if (loadbuf_size == MAXLINEBUFLEN) {
        logprint("\nError: Pathologically long line in .vcf file.\n");
        goto vcf_to_bed_ret_INVALID_FORMAT;
      }
      goto vcf_to_bed_ret_NOMEM;
    }
    bufptr = skip_initial_spaces(loadbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    bufptr2 = strchr(bufptr, '\t');
    if (!bufptr2) {
      goto vcf_to_bed_ret_INVALID_FORMAT_3;
    }
    ii = get_chrom_code(chrom_info_ptr, bufptr);
    if (ii == -1) {
      if (!allow_extra_chroms) {
	*bufptr2 = '\0';
        sprintf(logbuf, "\nError: Invalid chromosome code '%s' in .vcf file.\n(Use --allow-extra-chr to force it to be accepted.", bufptr);
	goto vcf_to_bed_ret_INVALID_FORMAT_2;
      }
      retval = resolve_or_add_chrom_name(chrom_info_ptr, bufptr, &ii);
      if (retval) {
	putchar('\n');
        goto vcf_to_bed_ret_1;
      }
    }
    if (!is_set(chrom_info_ptr->chrom_mask, ii)) {
      marker_skip_ct++;
      continue;
    }
    chrom_ptr = bufptr;
    chrom_len = (uintptr_t)(bufptr2 - bufptr);
    bufptr2++;
    pos_str = bufptr2;
    marker_id = strchr(bufptr2, '\t');
    if (!marker_id) {
      goto vcf_to_bed_ret_INVALID_FORMAT_3;
    }
    marker_id++;
    if ((((unsigned char)(*pos_str)) - '0') >= 10) {
      marker_id[-1] = '\0';
      sprintf(logbuf, "\nError: Invalid variant position '%s' in .vcf file.\n", bufptr2);
      goto vcf_to_bed_ret_INVALID_FORMAT_2;
    }
    bufptr3 = strchr(marker_id, '\t');
    if (!bufptr3) {
      goto vcf_to_bed_ret_INVALID_FORMAT_3;
    }
    marker_id_len = (uintptr_t)(bufptr3 - marker_id);
    bufptr3++;
    // now bufptr3 = ref allele
    bufptr = strchr(bufptr3, '\t');
    if (!bufptr) {
      goto vcf_to_bed_ret_INVALID_FORMAT_3;
    }
    bufptr++;
    alt_ct = 1;
    alt_alleles = bufptr;
    cc = *bufptr;
    // ',' < '.'
    while (1) {
      if ((unsigned char)cc <= ',') {
	logprint("\nError: Invalid alternate allele in .vcf file.\n");
	goto vcf_to_bed_ret_INVALID_FORMAT;
      }
      do {
	cc = *(++bufptr);
      } while ((unsigned char)cc > ',');
      if (cc != ',') {
	break;
      }
      cc = *(++bufptr);
      alt_ct++;      
    }
    if (cc != '\t') {
      goto vcf_to_bed_ret_INVALID_FORMAT_4;
    }
    if (biallelic_strict && (alt_ct > 1)) {
      goto vcf_to_bed_skip3;
    }
    bufptr++;
    bufptr2 = strchr(bufptr, '\t');
    if (!bufptr2) {
      goto vcf_to_bed_ret_INVALID_FORMAT_3;
    }
    if (check_qual) {
      if (*bufptr == '.') {
	marker_skip_ct++;
	continue;
      }
      if (scan_double(bufptr, &dxx)) {
        *bufptr2 = '\0';
        sprintf(logbuf, "\nError: Invalid QUAL value '%s' in .vcf file.\n", bufptr);
	goto vcf_to_bed_ret_INVALID_FORMAT_2;
      }
      if (dxx < vcf_min_qual) {
        marker_skip_ct++;
	continue;
      }
    }
    bufptr = &(bufptr2[1]);
    bufptr2 = strchr(bufptr, '\t');
    if (!bufptr2) {
      goto vcf_to_bed_ret_INVALID_FORMAT_3;
    }
    if (fexcept_ct) {
      *bufptr2 = '\0';
      ii = bsearch_str(bufptr, sorted_fexcepts, max_fexcept_len, 0, fexcept_ct - 1);
      if (ii == -1) {
	marker_skip_ct++;
	continue;
      }
    }
    bufptr = &(bufptr2[1]);
    bufptr2 = strchr(bufptr, '\t');
    if (!bufptr2) {
      goto vcf_to_bed_ret_INVALID_FORMAT_3;
    }
    bufptr = &(bufptr2[1]);
    bufptr2 = strchr(bufptr, '\t');
    if (!bufptr2) {
      goto vcf_to_bed_ret_INVALID_FORMAT_3;
    }
    if (memcmp(bufptr, "GT", 2)) {
      marker_skip_ct++;
      continue;
    }
    bufptr = &(bufptr2[1]);
    // okay, finally done with the line header
    geno_start = bufptr;
    if (alt_ct < 10) {
      // slightly faster parsing for the usual case
      // note that '.' is interpreted as an alternate allele since that doesn't
      // break anything
      fill_ulong_zero(base_bitfields, (alt_ct + 1) * indiv_ctv2);
      if ((!biallelic_only) || (alt_ct == 1)) {
	for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++, bufptr = &(bufptr2[1])) {
	  bufptr2 = strchr(bufptr, '\t');
	  if (!bufptr2) {
	    if (indiv_idx != indiv_ct - 1) {
	      goto vcf_to_bed_ret_INVALID_FORMAT_3;
	    }
	    bufptr2 = &(bufptr[strlen_se(bufptr)]);
	  }
	  if (*bufptr == '.') {
	    continue;
	  }
	  ulii = (unsigned char)(*bufptr) - '0';
	  cc = bufptr[1];
	  if ((cc != '/') && (cc != '|')) {
	    set_bit_ul(&(base_bitfields[ulii * indiv_ctv2]), indiv_idx * 2 + 1);
	  } else {
	    set_bit_ul(&(base_bitfields[ulii * indiv_ctv2]), indiv_idx * 2);
	    ulii = ((unsigned char)bufptr[2]) - '0';
	    base_bitfields[ulii * indiv_ctv2 + indiv_idx / BITCT2] += ONELU << (2 * (indiv_idx % BITCT2));
	  }
	}
	alt_allele_idx = 1;
	if (alt_ct > 1) {
	  ulii = popcount2_longs(&(base_bitfields[indiv_ctv2]), 0, indiv_ctl2);
	  for (alt_idx = 2; alt_idx <= alt_ct; alt_idx++) {
	    uljj = popcount2_longs(&(base_bitfields[indiv_ctv2 * alt_idx]), 0, indiv_ctl2);
	    if (uljj > ulii) {
	      ulii = uljj;
	      alt_allele_idx = alt_idx;
	    }
	  }
	}
      } else {
	// expect early termination in this case
	alt_allele_idx = 0;
	for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++, bufptr = &(bufptr2[1])) {
	  bufptr2 = strchr(bufptr, '\t');
	  if (!bufptr2) {
	    if (indiv_idx != indiv_ct - 1) {
	      goto vcf_to_bed_ret_INVALID_FORMAT_3;
	    }
	    bufptr2 = &(bufptr[strlen_se(bufptr)]);
	  }
	  if (*bufptr == '.') {
	    continue;
	  }
	  ulii = (unsigned char)(*bufptr) - '0';
	  if (ulii && (ulii != alt_allele_idx)) {
	    if (alt_allele_idx) {
	      goto vcf_to_bed_skip3;
	    }
	    alt_allele_idx = ulii;
	  }
	  cc = bufptr[1];
	  if ((cc != '/') && (cc != '|')) {
	    set_bit_ul(&(base_bitfields[ulii * indiv_ctv2]), indiv_idx * 2 + 1);
	  } else {
	    set_bit_ul(&(base_bitfields[ulii * indiv_ctv2]), indiv_idx * 2);
	    ulii = ((unsigned char)bufptr[2]) - '0';
	    if (ulii && (ulii != alt_allele_idx)) {
              if (alt_allele_idx) {
		goto vcf_to_bed_skip3;
	      }
	      alt_allele_idx = ulii;
	    }
	    base_bitfields[ulii * indiv_ctv2 + indiv_idx / BITCT2] += ONELU << (2 * (indiv_idx % BITCT2));
	  }
	}
	if (!alt_allele_idx) {
	  alt_allele_idx = 1;
	}
      }
      alt_bitfield = &(base_bitfields[alt_allele_idx * indiv_ctv2]);
    } else {
      // bleah, multi-digit genotype codes
      fill_uint_zero(vcf_alt_cts, alt_ct);
      logprint("\nError: --vcf does not yet support variants with 10+ alternate alleles.\n");
      goto vcf_to_bed_ret_1;
    }
    ref_ptr = base_bitfields;
    alt_ptr = alt_bitfield;
    for (indiv_idx = 0; indiv_idx < indiv_ctl2; indiv_idx++) {
      // take ref, then:
      // * if ref is nonzero, add 1
      // * if ref + alt is not two, force to 01
      ulii = *ref_ptr;
      uljj = *alt_ptr++;
      ulkk = (ulii + uljj) & AAAAMASK;
      uljj = ulii + ((ulii | (ulii >> 1)) & FIVEMASK);
      ulii = ulkk | (ulkk >> 1); // nonmissing?
      *ref_ptr++ = (uljj & ulii) | (((~ulkk) >> 1) & FIVEMASK);
    }
    ref_ptr[-1] &= final_mask;
    if (fwrite_checked(base_bitfields, indiv_ct4, outfile)) {
      goto vcf_to_bed_ret_WRITE_FAIL;
    }
    chrom_ptr[chrom_len] = '\0';
    fputs(chrom_ptr, bimfile);
    putc('\t', bimfile);
    fwrite(marker_id, 1, marker_id_len + 1, bimfile);
    putc('0', bimfile);
    putc('\t', bimfile);
    fwrite(pos_str, 1, marker_id - pos_str, bimfile);

    if (*alt_alleles == '.') {
      putc(missing_geno, bimfile);
    } else {
      bufptr = alt_alleles;
      for (alt_idx = 1; alt_idx < alt_allele_idx; alt_idx++) {
	bufptr = strchr(bufptr, ',');
	bufptr++;
      }
      bufptr2 = strchr(bufptr, (alt_allele_idx == alt_ct)? '\t' : ',');
      *bufptr2 = '\0';
      fputs(bufptr, bimfile);
    }
    putc('\t', bimfile);
    alt_alleles[-1] = '\n';
    *alt_alleles = '\0';
    if (fputs_checked(bufptr3, bimfile)) {
      goto vcf_to_bed_ret_WRITE_FAIL;
    }
    marker_ct++;
    if (!(marker_ct % 1000)) {
      printf("\r--vcf: %uk variants complete.", marker_ct / 1000);
      fflush(stdout);
    }
    continue;
  vcf_to_bed_skip3:
    if ((!marker_skip_ct) && skip3_list) {
      memcpy(outname_end, ".skip3", 7);
      if (fopen_checked(&skip3file, outname, "w")) {
	goto vcf_to_bed_ret_OPEN_FAIL;
      }
      memcpy(outname_end, ".bed", 5);
    }
    marker_id[marker_id_len] = '\0';
    if (fputs_checked(marker_id, skip3file)) {
      goto vcf_to_bed_ret_WRITE_FAIL;
    }
    putc('\n', skip3file);
    marker_skip_ct++;
  }
  putchar('\r');
  sprintf(logbuf, "--vcf: %s + .bim + .fam written.\n", outname);
  logprintb();
  if (marker_skip_ct) {
    sprintf(logbuf, "(%u variant%s skipped.)\n", marker_skip_ct, (marker_skip_ct == 1)? "" : "s");
    logprintb();
  }
  while (0) {
  vcf_to_bed_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  vcf_to_bed_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  vcf_to_bed_ret_READ_FAIL:
    putchar('\n');
    retval = RET_READ_FAIL;
    break;
  vcf_to_bed_ret_WRITE_FAIL:
    putchar('\n');
    retval = RET_WRITE_FAIL;
    break;
  vcf_to_bed_ret_INVALID_FORMAT_4:
    logprint("\nError: Improperly formatted .vcf file.\n");
    retval = RET_INVALID_FORMAT;
    break;
  vcf_to_bed_ret_INVALID_FORMAT_3:
    logprint("\nError: Fewer tokens than expected in .vcf line.\n");
    retval = RET_INVALID_FORMAT;
    break;
  vcf_to_bed_ret_INVALID_FORMAT_2:
    logprintb();
  vcf_to_bed_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 vcf_to_bed_ret_1:
  gzclose_cond(gz_infile);
  fclose_cond(outfile);
  fclose_cond(bimfile);
  fclose_cond(skip3file);
  wkspace_reset(wkspace_mark);
  return retval;
}

void announce_make_xylist_complete(uint32_t markers_left, char* outname, char* outname_end) {
  *outname_end = '\0';
  sprintf(logbuf, "--23file-make-xylist: %u variant ID%s written to %s.xylist.\n", markers_left, (markers_left == 1)? "" : "s", outname);
  logprintb();
}

uint32_t write_23_cached_chrom(char* write_cache, uint32_t markers_left, char chrom_second_char, FILE* outfile_bed, FILE* outfile_bim) {
  // N.B. must flush writebuf before calling this
  uint32_t uii;
  do {
    uii = strlen(write_cache);
    if (putc_checked('2', outfile_bim)) {
      return 1;
    }
    putc(chrom_second_char, outfile_bim);
    if (fwrite_checked(write_cache, uii, outfile_bim)) {
      return 1;
    }
    write_cache = &(write_cache[uii + 1]);
    if (putc_checked(*write_cache++, outfile_bed)) {
      return 1;
    }
  } while (--markers_left);
  return 0;
}

int32_t bed_from_23(char* infile_name, char* outname, char* outname_end, uint32_t modifier_23, char* fid_23, char* iid_23, double pheno_23, char* paternal_id_23, char* maternal_id_23, char* convert_xy_fname, Chrom_info* chrom_info_ptr) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* infile_23 = NULL;
  FILE* outfile_bed = NULL;
  FILE* outfile_txt = NULL;
  FILE* xylist_file = NULL;
  uint32_t is_male = modifier_23 & M23_MALE;
  uint32_t is_female = modifier_23 & M23_FEMALE;
  uint32_t do_convert_xy = modifier_23 & M23_CONVERT_XY;
  uint32_t make_xylist = modifier_23 & M23_MAKE_XYLIST;
  char* xylist = NULL;
  uint32_t y_present = 0;
  uint32_t nonmissing_y_present = 0;
  uintptr_t xylist_ct = 0;
  uintptr_t xylist_max_id_len = 0;
  uint32_t xychrom_written = 0;
  unsigned char* writebuf = (unsigned char*)(&(tbuf[MAXLINELEN]));
  char* writebuf2 = &(tbuf[MAXLINELEN * 2]);
  // xy_phase:
  // 0 = start of X chromosome, in initial pseudoautosomal region
  // 1 = in middle, out of region
  // 2 = at end, in final pseudoautosomal region
  uint32_t xy_phase = 0;
  int32_t retval = 0;
  uint32_t cur_chrom = 0;
  uint32_t chrom_mask_23 = (uint32_t)(chrom_info_ptr->chrom_mask[0]);
  uint32_t keep_x = ((uint32_t)(chrom_mask_23 >> 23)) & 1;
  uint32_t keep_xy = ((uint32_t)(chrom_mask_23 >> 25)) & 1;
  uint32_t indel_ct = 0;
  uint32_t xy_marker_ct = 0;
  uint32_t mid_x_marker_ct = 0;
  unsigned char* writebuf_cur;
  unsigned char* writebuf_end;
  char* new_xy_buf;
  char* new_xy_buf_cur;
  char* new_xy_buf_end;
  char* last_early_xy;
  char* id_start;
  char* chrom_start;
  char* pos_start;
  char* allele_start;
  char* writebuf2_cur;
  uintptr_t id_len;
  uint32_t allele_calls;
  uint32_t is_xy;
  uint32_t uii;
  int32_t ii;
  char cc;
  char cc2;
  unsigned char ucc;
  if (fopen_checked(&infile_23, infile_name, "r")) {
    goto bed_from_23_ret_OPEN_FAIL;
  }
  memcpy(outname_end, ".bim", 5);
  if (fopen_checked(&outfile_txt, outname, "w")) {
    goto bed_from_23_ret_OPEN_FAIL;
  }
  memcpy(&(outname_end[2]), "ed", 2);
  if (fopen_checked(&outfile_bed, outname, "wb")) {
    goto bed_from_23_ret_OPEN_FAIL;
  }
  if (convert_xy_fname) {
    retval = open_and_size_string_list(convert_xy_fname, &xylist_file, &xylist_ct, &xylist_max_id_len);
    if (retval) {
      goto bed_from_23_ret_1;
    }
#ifndef __LP64__
    if (((uint64_t)xylist_ct) * xylist_max_id_len > 0x7fffffffLLU) {
      goto bed_from_23_ret_NOMEM;
    }
#endif
    if (wkspace_alloc_c_checked(&xylist, xylist_ct * xylist_max_id_len)) {
      goto bed_from_23_ret_NOMEM;
    }
    retval = load_string_list(&xylist_file, xylist_max_id_len, xylist);
    if (retval) {
      goto bed_from_23_ret_1;
    }
    fclose_null(&xylist_file);
    qsort(xylist, xylist_ct, xylist_max_id_len, strcmp_casted);
  } else {
    if ((do_convert_xy || make_xylist) && (modifier_23 & M23_FORCE_MISSING_SEX)) {
      is_male = 1;
    }
    if (make_xylist) {
      memcpy(outname_end, ".xylist", 8);
      if (fopen_checked(&xylist_file, outname, "w")) {
        goto bed_from_23_ret_OPEN_FAIL;
      }
    }
  }
  if (wkspace_left < MAXLINELEN) {
    goto bed_from_23_ret_NOMEM;
  }
  new_xy_buf = (char*)wkspace_base;
  new_xy_buf_cur = new_xy_buf;
  last_early_xy = new_xy_buf;
  new_xy_buf_end = (char*)(&(wkspace_base[wkspace_left - MAXLINELEN]));
  writebuf_cur = (unsigned char*)memcpyl3a((char*)writebuf, "l\x1b\x01");
  writebuf_end = &(writebuf[MAXLINELEN]);
  tbuf[MAXLINELEN - 1] = ' ';
  while (fgets(tbuf, MAXLINELEN, infile_23)) {
    if (!tbuf[MAXLINELEN - 1]) {
      sprintf(logbuf, "Error: Pathologically long line in %s.\n", infile_name);
      logprintb();
      goto bed_from_23_ret_INVALID_FORMAT;
    }
    id_start = skip_initial_spaces(tbuf);
    cc = *id_start;
    if (is_eoln_kns(cc) || (cc == '#')) {
      continue;
    }
    chrom_start = item_endnn(id_start);
    id_len = (uintptr_t)(chrom_start - id_start);
    chrom_start = skip_initial_spaces(chrom_start);
    pos_start = next_item(chrom_start);
    allele_start = next_item(pos_start);
    if (no_more_items_kns(allele_start)) {
      logprint("Error: Fewer tokens than expected in --23file line.\n");
      goto bed_from_23_ret_INVALID_FORMAT;
    }
    allele_calls = strlen_se(allele_start);
    if (allele_calls > 2) {
      logprint("Error: Too many allele calls in --23file line.\n");
      goto bed_from_23_ret_INVALID_FORMAT;
    }
    ii = get_chrom_code(SPECIES_HUMAN, chrom_start);
    if (ii == -1) {
      logprint("Error: Invalid chromosome code in --23file line.\n");
      goto bed_from_23_ret_INVALID_FORMAT;
    }
    uii = (uint32_t)ii;
    if ((!(chrom_mask_23 & (1 << uii))) && ((uii != 23) || (!(keep_x || keep_xy)))) {
      continue;
    }
    if (uii == 25) {
      is_xy = 1;
      uii = 23;
      xy_phase = 2;
    } else {
      is_xy = 0;
    }
    if (uii < cur_chrom) {
      // allow 24 -> XY
      if ((!is_xy) || (cur_chrom == 26)) {
        sprintf(logbuf, "Error: Chromosomes in %s are out of order.\n", infile_name);
        goto bed_from_23_ret_INVALID_FORMAT;
      }
      cur_chrom = 25;
    } else if (uii > cur_chrom) {
      cur_chrom = uii;
      if ((uii == 26) && xy_marker_ct) {
	if (writebuf_cur != writebuf) {
	  if (fwrite_checked(writebuf, (uintptr_t)(writebuf_cur - writebuf), outfile_bed)) {
	    goto bed_from_23_ret_WRITE_FAIL;
	  }
	  writebuf_cur = writebuf;
	}
	if (xylist_file) {
	  announce_make_xylist_complete(xy_marker_ct, outname, outname_end);
	}
	if (do_convert_xy) {
          if (write_23_cached_chrom(new_xy_buf, xy_marker_ct, '5', outfile_bed, outfile_txt)) {
	    goto bed_from_23_ret_WRITE_FAIL;
	  }
	}
	xychrom_written = 1;
      }
    }
    cc2 = allele_start[0];
    if ((cur_chrom == 23) && (!is_xy) && (do_convert_xy || make_xylist)) {
      if (xylist) {
	// This should be robust to changes in the 23andMe marker set.
	// Strategy: We actually only care about the end of the beginning, and
	// the beginning of the end.  So we just look for the last hit in the
	// first half or so of the chromosome, and the first hit in the last
	// half.  While the precise bp length of the X chromosome varies
	// between reference builds, position 77000000 should be safely in the
	// middle.
	is_xy = 1;
	// if !last_early_xy, then we already saw "first hit in last half" and
	// force the marker to be on XY.
	if (last_early_xy) {
	  id_start[id_len] = '\0';
	  if (bsearch_str(id_start, xylist, xylist_max_id_len, 0, xylist_ct - 1) != -1) {
	    if (atoiz(pos_start, &ii)) {
	      pos_start[strlen_se(pos_start)] = '\0';
	      sprintf(logbuf, "Error: Invalid --23file variant position '%s'.\n", pos_start);
	      goto bed_from_23_ret_INVALID_FORMAT;
	    }
	    if (((uint32_t)ii) > 77000000) {
	      if (mid_x_marker_ct) {
		if (xy_marker_ct) {
		  last_early_xy = &(last_early_xy[strlen(last_early_xy) + 2]);
		}
		if (writebuf_cur != writebuf) {
		  if (fwrite_checked(writebuf, (uintptr_t)(writebuf_cur - writebuf), outfile_bed)) {
		    goto bed_from_23_ret_WRITE_FAIL;
		  }
		  writebuf_cur = writebuf;
		}
		if (write_23_cached_chrom(last_early_xy, mid_x_marker_ct, '3', outfile_bed, outfile_txt)) {
		  goto bed_from_23_ret_WRITE_FAIL;
		}
	      }
	      new_xy_buf_cur = last_early_xy;
	      last_early_xy = NULL;
	    } else {
	      last_early_xy = new_xy_buf_cur;
	      xy_marker_ct += mid_x_marker_ct;
	      mid_x_marker_ct = 0;
	      xy_phase = 0;
	    }
	  } else {
	    is_xy = 0;
	    mid_x_marker_ct++;
	  }
	}
      } else if (cc2 == '-') {
	// status quo
	if (xy_phase != 1) {
	  is_xy = 1;
	}
      } else if (allele_calls == 2) {
	is_xy = 1;
      }
      if (is_xy) {
        if (xy_phase == 1) {
	  xy_phase = 2;
	}
      } else if (xy_phase != 1) {
	if (xy_phase == 2) {
	  logprint("Error: Invalid X chromosome data (more than one haploid component).\n");
	  goto bed_from_23_ret_INVALID_FORMAT;
	}
	xy_phase = 1;
      }
    } else if ((cur_chrom == 24) && (!nonmissing_y_present)) {
      y_present = 1;
      if (cc2 != '-') {
	nonmissing_y_present = 1;
      }
    }
    if (cc2 == '-') {
      ucc = 1;
      cc = '0';
      cc2 = '0';
    } else if (allele_calls == 2) {
      cc = allele_start[1];
      if (cc == cc2) {
	cc = '0';
	ucc = 3;
      } else {
	ucc = 2;
      }
    } else {
      cc = '0';
      ucc = 3;
    }
    if ((cc2 == 'D') || (cc2 == 'I')) {
      indel_ct++;
      cc = (char)(((unsigned char)cc2) ^ 13); // swaps D and I
    }
    if (is_xy) {
      if (allele_calls != 2) {
	goto bed_from_23_ret_INVALID_FORMAT_2;
      }
      if (!keep_xy) {
	continue;
      }
      if (new_xy_buf_cur > new_xy_buf_end) {
	goto bed_from_23_ret_NOMEM;
      }
      if (xylist_file) {
        fwrite(id_start, 1, id_len, xylist_file);
        if (putc_checked('\n', xylist_file)) {
	  goto bed_from_23_ret_WRITE_FAIL;
	}
      }
      xy_marker_ct++;
    }
    if (do_convert_xy && (cur_chrom == 23) && (xylist || is_xy)) {
      *new_xy_buf_cur++ = '\t';
      new_xy_buf_cur = memcpya(new_xy_buf_cur, id_start, id_len);
      new_xy_buf_cur = memcpyl3a(new_xy_buf_cur, "\t0\t");
      new_xy_buf_cur = memcpyax(new_xy_buf_cur, pos_start, strlen_se(pos_start), '\t');
      *new_xy_buf_cur++ = cc;
      *new_xy_buf_cur++ = '\t';
      *new_xy_buf_cur++ = cc2;
      *new_xy_buf_cur++ = '\n';
      *new_xy_buf_cur++ = '\0';
      *new_xy_buf_cur++ = (char)ucc;
    } else {
      if ((cur_chrom == 23) && (!is_xy) && (!keep_x)) {
	continue;
      }
      if ((allele_calls == 1) && (cur_chrom <= 23) && (!is_xy)) {
	if ((cur_chrom == 23) && (!is_female)) {
	  is_male = 1;
	} else {
	  logprint("Error: Too few allele calls in --23file line.\n");
	  goto bed_from_23_ret_INVALID_FORMAT;
	}
      } else if ((cur_chrom == 24) && (cc2 != '0') && (!is_male)) {
	if (!is_female) {
	  is_male = 1;
	} else {
	  logprint("Error: Nonmissing female --23file allele call on Y chromosome.\n");
	  goto bed_from_23_ret_INVALID_FORMAT;
	}
      }
      writebuf2_cur = uint32_writex(writebuf2, cur_chrom, '\t');
      writebuf2_cur = memcpya(writebuf2_cur, id_start, id_len);
      writebuf2_cur = memcpyl3a(writebuf2_cur, "\t0\t");
      writebuf2_cur = memcpyax(writebuf2_cur, pos_start, strlen_se(pos_start), '\t');
      *writebuf2_cur++ = cc;
      *writebuf2_cur++ = '\t';
      *writebuf2_cur++ = cc2;
      *writebuf2_cur++ = '\n';
      if (fwrite_checked(writebuf2, (uintptr_t)(writebuf2_cur - writebuf2), outfile_txt)) {
	goto bed_from_23_ret_WRITE_FAIL;
      }
      if (writebuf_cur == writebuf_end) {
	if (fwrite_checked(writebuf, (uintptr_t)(writebuf_cur - writebuf), outfile_bed)) {
	  goto bed_from_23_ret_WRITE_FAIL;
	}
	writebuf_cur = writebuf;
      }
      *writebuf_cur++ = (char)ucc;
    }
  }
  if (!feof(infile_23)) {
    goto bed_from_23_ret_READ_FAIL;
  }
  if ((writebuf_cur == &(writebuf[3])) && (writebuf[0] == 'l')) {
    if (chrom_mask_23 == 0x7ffffff) {
      logprint("Error: No --23file variants.\n");
      goto bed_from_23_ret_INVALID_FORMAT;
    } else {
      logprint("Error: No --23file variants pass chromosome filter.\n");
      goto bed_from_23_ret_INVALID_CMDLINE;
    }
  }
  if (fwrite_checked(writebuf, (uintptr_t)(writebuf_cur - writebuf), outfile_bed)) {
    goto bed_from_23_ret_WRITE_FAIL;
  }
  if (xy_marker_ct && (!xychrom_written)) {
    if (do_convert_xy) {
      if (write_23_cached_chrom(new_xy_buf_cur, xy_marker_ct, '5', outfile_bed, outfile_txt)) {
	goto bed_from_23_ret_WRITE_FAIL;
      }
    }
    if (xylist_file) {
      announce_make_xylist_complete(xy_marker_ct, outname, outname_end);
    }
  }
  if (fclose_null(&outfile_txt)) {
    goto bed_from_23_ret_WRITE_FAIL;
  }
  memcpy(outname_end, ".fam", 5);
  if (fopen_checked(&outfile_txt, outname, "w")) {
    goto bed_from_23_ret_OPEN_FAIL;
  }
  if (fid_23) {
    fputs(fid_23, outfile_txt);
    putc(' ', outfile_txt);
  } else {
    fputs("FAM001 ", outfile_txt);
  }
  if (iid_23) {
    fputs(iid_23, outfile_txt);
    putc(' ', outfile_txt);
  } else {
    fputs("ID001 ", outfile_txt);
  }
  if (paternal_id_23) {
    fputs(paternal_id_23, outfile_txt);
    putc(' ', outfile_txt);
  } else {
    fputs("0 ", outfile_txt);
  }
  if (maternal_id_23) {
    fputs(maternal_id_23, outfile_txt);
  } else {
    putc('0', outfile_txt);
  }
  if (modifier_23 & M23_FORCE_MISSING_SEX) {
    cc = '0';
  } else if (is_male) {
    cc = '1';
  } else {
    cc = '2';
  }
  fprintf(outfile_txt, " %c %g\n", cc, pheno_23);
  if (fclose_null(&outfile_txt)) {
    goto bed_from_23_ret_WRITE_FAIL;
  }
  *outname_end = '\0';
  sprintf(logbuf, "--23file: binary fileset written to %s.bed + .bim + .fam.\n", outname);
  logprintb();
  if (indel_ct) {
    sprintf(logbuf, "%u variants with indel calls present.  --list-23-indels may be useful here.\n", indel_ct);
    logprintb();
  }
  if (!(modifier_23 & M23_SEX)) {
    sprintf(logbuf, "Inferred sex: %smale.\n", is_male? "" : "fe");
    logprintb();
  }
  if ((modifier_23 & M23_MALE) && y_present && (!nonmissing_y_present)) {
    if (keep_x && xy_marker_ct) {
      if (!xy_phase) {
	logprint("Warning: No explicit haploid calls on X chromosome, and no nonmissing calls on Y\nchromosome.  Double-check whether this is really a male sample.\n");
      }
    } else {
      logprint("Warning: No nonmissing calls on Y chromosome.  Double-check whether this is\nreally a male sample.\n");
    }
  }
  while (0) {
  bed_from_23_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  bed_from_23_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  bed_from_23_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  bed_from_23_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  bed_from_23_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  bed_from_23_ret_INVALID_FORMAT_2:
    logprint("Error: Too few allele calls in --23file line.\n");
  bed_from_23_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 bed_from_23_ret_1:
  fclose_cond(infile_23);
  fclose_cond(outfile_bed);
  fclose_cond(outfile_txt);
  fclose_cond(xylist_file);
  wkspace_reset(wkspace_mark);
  return retval;
}

int32_t generate_dummy(char* outname, char* outname_end, uint32_t flags, uintptr_t marker_ct, uintptr_t indiv_ct, double geno_mrate, double pheno_mrate) {
  FILE* outfile = NULL;
  unsigned char* wkspace_mark = wkspace_base;
  int32_t retval = 0;
  uintptr_t indiv_ct4 = (indiv_ct + 3) / 4;
  uint32_t dbl_indiv_mod4 = 2 * (indiv_ct % 4);
  uint32_t four_alleles = 0;
  uint32_t geno_m_check = (geno_mrate > 0.0);
  uint32_t geno_m32 = (uint32_t)(geno_mrate * 4294967296.0);
  uint32_t pheno_m_check = (pheno_mrate > 0.0);
  uint32_t pheno_m32 = (uint32_t)(geno_mrate * 4294967296.0);
  uintptr_t urand = 0;
  uint32_t saved_rnormal = 0;
  char wbuf[64];
  char* wptr = &(wbuf[5]);
  char* wptr2;
  double saved_rnormal_val;
  char alleles[13];
  unsigned char* writebuf;
  unsigned char* ucptr;
  unsigned char* ucptr2;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t pct;
  uint32_t loop_end;
  unsigned char ucc;
  unsigned char ucc2;
  unsigned char bmap[320];
  unsigned char* bmap2;
  uint64_t ullii;
  double dxx;
  fill_bmap_short(bmap, indiv_ct % 4);
  bmap2 = &(bmap[256]);
  if (flags & DUMMY_ACGT) {
    memcpy(alleles, "ACAGATCGCTGTA", 13);
    four_alleles = 1;
  } else if (flags & DUMMY_1234) {
    memcpy(alleles, "1213142324341", 13);
    four_alleles = 1;
  } else if (flags & DUMMY_12) {
    memcpyl3(alleles, "121");
  } else {
    memcpyl3(alleles, "ABA");
  }
  if (wkspace_alloc_uc_checked(&writebuf, indiv_ct4)) {
    goto generate_dummy_ret_NOMEM;
  }
  memcpy(outname_end, ".bim", 5);
  if (fopen_checked(&outfile, outname, "w")) {
    goto generate_dummy_ret_OPEN_FAIL;
  }
  memcpy(wbuf, "1\tsnp", 5);
  if (four_alleles) {
    for (uii = 0; uii < marker_ct; uii++) {
      if (!(uii % 8)) {
	do {
	  urand = sfmt_genrand_uint32(&sfmt);
	} while (urand < 425132032LU); // 2^32 - 12^8.  heck, why not
      }
      ukk = urand / 12U;
      ujj = urand - (ukk * 12U);
      urand = ukk;
      wptr2 = uint32_writex(memcpyl3a(uint32_write(wptr, uii), "\t0\t"), uii, '\t');
      wptr2[0] = alleles[ujj];
      wptr2[1] = '\t';
      wptr2[2] = alleles[ujj + 1];
      wptr2[3] = '\n';
      if (fwrite_checked(wbuf, 4 + (uintptr_t)(wptr2 - wbuf), outfile)) {
	goto generate_dummy_ret_WRITE_FAIL;
      }
    }
  } else {
    for (uii = 0; uii < marker_ct; uii++) {
      if (!(uii % 32)) {
	urand = sfmt_genrand_uint32(&sfmt);
      }
      ujj = urand & 1;
      urand >>= 1;
      wptr2 = uint32_writex(memcpyl3a(uint32_write(wptr, uii), "\t0\t"), uii, '\t');
      wptr2[0] = alleles[ujj];
      wptr2[1] = '\t';
      wptr2[2] = alleles[ujj + 1];
      wptr2[3] = '\n';
      if (fwrite_checked(wbuf, 4 + (uintptr_t)(wptr2 - wbuf), outfile)) {
	goto generate_dummy_ret_WRITE_FAIL;
      }
    }
  }
  if (fclose_null(&outfile)) {
    goto generate_dummy_ret_WRITE_FAIL;
  }
  memcpy(outname_end, ".fam", 5);
  if (fopen_checked(&outfile, outname, "w")) {
    goto generate_dummy_ret_OPEN_FAIL;
  }
  wptr = memcpyl3a(wbuf, "per");
  if (flags & DUMMY_SCALAR_PHENO) {
    for (uii = 0; uii < indiv_ct; uii++) {
      if (pheno_m_check && (sfmt_genrand_uint32(&sfmt) <= pheno_m32)) {
	dxx = -9;
      } else {
	if (saved_rnormal) {
	  dxx = saved_rnormal_val;
	  saved_rnormal = 0;
	} else {
	  dxx = rand_normal(&saved_rnormal_val);
	  saved_rnormal = 1;
	}
      }
      wptr2 = double_g_writex(memcpya(uint32_write(memcpya(uint32_write(wptr, uii), " per", 4), uii), " 0 0 2 ", 7), dxx, '\n');
      if (fwrite_checked(wbuf, wptr2 - wbuf, outfile)) {
	goto generate_dummy_ret_WRITE_FAIL;
      }
    }
  } else {
    for (uii = 0; uii < indiv_ct; uii++) {
      if (!(uii % 32)) {
	urand = sfmt_genrand_uint32(&sfmt);
      }
      wptr2 = uint32_write(memcpya(uint32_write(wptr, uii), " per", 4), uii);
      if (pheno_m_check && (sfmt_genrand_uint32(&sfmt) <= pheno_m32)) {
	wptr2 = memcpya(wptr2, " 0 0 2 -9\n", 10);
      } else {
	wptr2 = memcpya(wptr2, " 0 0 2 ", 7);
	*wptr2++ = (char)((urand & 1) + '1');
	*wptr2++ = '\n';
      }
      if (fwrite_checked(wbuf, wptr2 - wbuf, outfile)) {
	goto generate_dummy_ret_WRITE_FAIL;
      }
      urand >>= 1;
    }
  }
  if (fclose_null(&outfile)) {
    goto generate_dummy_ret_WRITE_FAIL;
  }
  memcpy(outname_end, ".bed", 5);
  if (fopen_checked(&outfile, outname, "wb")) {
    goto generate_dummy_ret_OPEN_FAIL;
  }
  if (fwrite_checked("l\x1b\x01", 3, outfile)) {
    goto generate_dummy_ret_WRITE_FAIL;
  }
  uii = 0;
  ullii = (3 * ONELU) + ((uint64_t)marker_ct) * indiv_ct4;
  if (ullii >= 10485760) {
    printf("Writing dummy .bed (%" PRIu64 " MB)... 0%%", ullii >> 20);
  } else {
    fputs("Writing dummy .bed... 0%", stdout);
  }
  fflush(stdout);
  for (pct = 1; pct <= 100; pct++) {
    loop_end = ((uint64_t)(pct * marker_ct)) / 100U;
    for (; uii < loop_end; uii++) {
      ucptr = writebuf;
      for (ujj = 0; ujj < indiv_ct4; ujj++) {
	if (!(ujj % 4)) {
	  urand = sfmt_genrand_uint32(&sfmt);
	}
	ucc = 0;
	for (ukk = 0; ukk < 8; ukk += 2) {
	  if (geno_m_check && (sfmt_genrand_uint32(&sfmt) < geno_m32)) {
	    ucc2 = 1;
	  } else {
	    ucc2 = urand & 3;
	    if (ucc2 == 1) {
	      ucc2 = 2;
	    }
	  }
	  ucc |= ucc2 << ukk;
	  urand >>= 2;
	}
	*ucptr++ = ucc;
      }
      if (dbl_indiv_mod4) {
	ucc = *(--ucptr);
	*ucptr = ucc >> dbl_indiv_mod4;
      }

      ujj = popcount_chars((uintptr_t*)writebuf, 0, indiv_ct4);
      if (ujj < indiv_ct) {
	ucptr = writebuf;
	ucptr2 = &(writebuf[indiv_ct / 4]);
	while (ucptr < ucptr2) {
	  *ucptr = bmap[*ucptr];
	  ucptr++;
	}
	if (dbl_indiv_mod4) {
	  *ucptr = bmap2[*ucptr];
	}
      }
      if (fwrite_checked(writebuf, indiv_ct4, outfile)) {
	putchar('\n');
	goto generate_dummy_ret_WRITE_FAIL;
      }
    }
    if (pct < 100) {
      if (pct > 10) {
	putchar('\b');
      }
      printf("\b\b%u%%", pct);
      fflush(stdout);
    }
  }
  putchar('\r');
  sprintf(logbuf, "Dummy data generated (%" PRIuPTR " %s, %" PRIuPTR " SNP%s).\n", indiv_ct, species_str(indiv_ct), marker_ct, (marker_ct == 1)? "" : "s");
  logprintb();
  while (0) {
  generate_dummy_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  generate_dummy_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  generate_dummy_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
  }
  fclose_cond(outfile);
  wkspace_reset(wkspace_mark);
  return retval;
}

void simulate_init_freqs_qt(uint32_t tags_or_haps, double dprime, double qt_var, double qt_dom, double missing_freq, double* freqs, uint64_t* thresholds, double* qt_adj) {
  // Initialize frequency table for current SNP.  Similar to instanceSNP_QT()
  // in PLINK 1.07 simul.cpp.  (The difference is that heterozygote frequencies
  // are not stored in a redundant manner.)
  // thresholds[0]: P(causal 11, marker 11) * 2^63
  // thresholds[1]: P(causal 11, marker 11 or missing) * 2^63
  // thresholds[2]: P(causal 11, marker 11 or missing or 12) * 2^63
  // ...
  // thresholds[14]: [1 - P(causal 22, marker 22)] * 2^63
  // We use 2^63 instead of 2^64 to avoid integer overflow headaches.
  double freq = freqs[0];
  double mfreq = freqs[1];
  double nonmissing_freq = 1 - missing_freq;
  double scaled_miss_freq = missing_freq * TWO_63;
  double scaled_nm_freq = nonmissing_freq * TWO_63;
  double ld = freq * (1 - mfreq);

  double qq = 1 - freq;
  double aa = sqrt(qt_var / (2 * freq * qq * ((1 + qt_dom * (qq - freq)) * (1 + qt_dom * (qq - freq)) + qt_dom * 2 * freq * qq * qt_dom)));
  double dd = qt_dom * aa;
  double dxx = aa * (1 - 2 * freq * (1 + qq * qt_dom));

  // Joint causal variant/marker probs, considering one allele at a time.
  // First digit indicates causal allele, second digit indicates marker.
  double h21 = qq * mfreq;
  double h11;
  double h12;
  double h22;

  // first 2 digits indicate causal variant genotype, last 2 indicate marker.
  // this is DIFFERENT from instanceSNP()'s usage.
  double h_11_11;
  double h_11_12;
  double h_11_22;
  double h_12_11;
  double h_12_12;
  double h_12_22;
  double h_22_11;
  double h_22_12;
  double h_22_22;

  // Constraints:
  // 1. mean = (qq^2) * qt_adj[3] + 2 * freq * qq * qt_adj[2] +
  //           (freq^2) * qt_adj[0]
  //         = 0
  // 2. variance = (qq^2) * (qt_adj[3]^2) + 2 * freq * qq * (qt_adj[2]^2) +
  //               (freq^2) * (qt_adj[0]^2)
  //             = qt_var
  // 3. qt_adj[2] = ((qt_adj[0] + qt_adj[3]) / 2) +
  //                qt_dom * ((qt_adj[0] - qt_adj[3]) / 2)
  // PLINK 1.07 computes the correct values, but incorrectly reverses their
  // positions, and also incorrectly uses sp[s].gAB in some places where .gBB
  // should be used and vice versa.
  qt_adj[0] = dxx + aa;
  qt_adj[2] = dxx + dd;
  qt_adj[3] = dxx - aa;

  if (h21 < ld) {
    ld = h21;
  }
  ld *= dprime;
  h11 = freq * mfreq + ld;
  h12 = freq * (1 - mfreq) - ld;
  h21 -= ld;
  h22 = qq * (1 - mfreq) + ld;
  h_11_11 = h11 * h11;
  h_11_12 = h11 * h12 * 2;
  h_11_22 = h12 * h12;
  h_12_11 = h21 * h11 * 2;
  h_12_12 = (h22 * h11 + h21 * h12) * 2;
  h_12_22 = h22 * h12 * 2;
  h_22_11 = h21 * h21;
  h_22_12 = h22 * h21 * 2;
  h_22_22 = h22 * h22;
  // determination of causal variant missing/nonmissing status is deferred, to
  // simplify phenotype updating
  if (!tags_or_haps) {
    thresholds[0] = (uint64_t)((h_11_11 + h_12_11 + h_22_11) * TWO_63);
    thresholds[1] = thresholds[0] + (uint64_t)((h_11_12 + h_12_12 + h_22_12) * TWO_63);
  } else {
    thresholds[0] = (uint64_t)(h_11_11 * scaled_nm_freq);
    thresholds[1] = thresholds[0] + (uint64_t)((h_11_11 + h_11_12 + h_11_22) * scaled_miss_freq);
    thresholds[2] = thresholds[1] + (uint64_t)(h_11_12 * scaled_nm_freq);
    thresholds[3] = thresholds[2] + (uint64_t)(h_11_22 * scaled_nm_freq);
    thresholds[4] = thresholds[3] + (uint64_t)(h_12_11 * scaled_nm_freq);
    thresholds[5] = thresholds[4] + (uint64_t)((h_12_11 + h_12_12 + h_12_22) * scaled_miss_freq);
    thresholds[6] = thresholds[5] + (uint64_t)(h_12_12 * scaled_nm_freq);
    thresholds[7] = thresholds[6] + (uint64_t)(h_12_22 * scaled_nm_freq);
    thresholds[8] = thresholds[7] + (uint64_t)(h_22_11 * scaled_nm_freq);
    thresholds[9] = thresholds[8] + (uint64_t)((h_22_11 + h_22_12 + h_22_22) * scaled_miss_freq);
    thresholds[10] = thresholds[9] + (uint64_t)(h_22_12 * scaled_nm_freq);
  }
}

void simulate_cc_get_conditional_probs(double prevalence, double g0, double g1, double g2, double het_odds, double hom0_odds, double* f0p, double* f1p, double* f2p) {
  // PLINK 1.07 interprets het_odds and hom0_odds as probability ratios instead
  // of odds ratios.  The two are nearly identical if prevalence is small, but
  // it's still better to avoid that approximation.  This requires solving a
  // cubic equation in X := odds(hom2).
  //
  // prevalence = g0 * (hom0_odds * X) / (1 + hom0_odds * X)
  //            + g1 * (het_odds * X) / (1 + het_odds * X)
  //            + g2 * X / (1 + X)
  //
  // prevalence * (1 + X) * (1 + het_odds * X) * (1 + hom0_odds * X) =
  //   X * (g0 * hom0_odds * (1 + het_odds * X) * (1 + X) +
  //        g1 * het_odds * (1 + hom0_odds * X) * (1 + X) +
  //        g2 * (1 + hom0_odds * X) * (1 + het_odds * X))
  //
  //   X^3 * prevalence * het_odds * hom0_odds
  // + X^2 * prevalence * (het_odds * hom0_odds + het_odds + hom0_odds)
  // + X * prevalence * (1 + het_odds + hom0_odds)
  // + prevalence =
  //   X^3 * het_odds * hom0_odds
  // + X^2 * (g0 * hom0_odds * (1 + het_odds) +
  //          g1 * het_odds * (1 + hom0_odds) +
  //          g2 * (hom0_odds + het_odds))
  // + X * (g0 * hom0_odds + g1 * het_odds + g2)
  //
  // 0 = X^3 * het_odds * hom0_odds * (1 - prevalence)
  //   + X^2 * (g0 * hom0_odds * (1 + het_odds) +
  //            g1 * het_odds * (1 + hom0_odds) +
  //            g2 * (hom0_odds + het_odds) -
  //            prevalence * (het_odds * hom0_odds + het_odds + hom0_odds))
  //   + X * (g0 * hom0_odds + g1 * het_odds + g2 -
  //          prevalence * (1 + het_odds + hom0_odds)
  //   - prevalence
  //
  // Cubic discriminant:
  //   27a^2d^2 + 4ac^3 + 4b^3d - b^2c^2 - 18abcd
  //   (Wikipedia article currently has the sign wrong!)
  //
  // ...eh, screw it, casework for applying the formula is a bit too annoying
  // for me to implement it for now.  Maybe if one or two other cubics come up,
  // I'll do it.  Meanwhile, here's a simple numeric binary search.
  uint32_t iters_left = 53;
  double min_f2 = 0.0;
  double max_f2 = 1.0;
  double cur_f2;
  double cur_f2_odds;
  double cur_f0;
  double cur_f1;
  double cur_prevalence;
  if ((prevalence == 0) || (prevalence == 1)) {
    *f0p = prevalence;
    *f1p = prevalence;
    *f2p = prevalence;
    return;
  }
  do {
    cur_f2 = 0.5 * (min_f2 + max_f2);
    cur_f2_odds = cur_f2 / (1 - cur_f2);
    // odds = p / (1 - p)
    // -> p = odds / (1 + odds)
    cur_f0 = cur_f2_odds * hom0_odds;
    cur_f1 = cur_f2_odds * het_odds;
    cur_f0 = cur_f0 / (1 + cur_f0);
    cur_f1 = cur_f1 / (1 + cur_f1);
    cur_prevalence = g0 * cur_f0 + g1 * cur_f1 + g2 * cur_f2;
    if (cur_prevalence > prevalence) {
      max_f2 = cur_f2;
    } else if (cur_prevalence < prevalence) {
      min_f2 = cur_f2;
    } else {
      break;
    }
  } while (--iters_left);
  *f0p = cur_f0;
  *f1p = cur_f1;
  *f2p = cur_f2;
}

void simulate_init_freqs_cc(uint32_t do_haps, double dprime, double* freqs, double prevalence, double het_odds, double hom0_odds, double missing_freq, uint64_t* ctrl_thresholds, uint64_t* case_thresholds) {
  // Similar to instanceSNP().
  double freq = freqs[0];
  double mfreq = freqs[1];
  double nonmissing_freq = 1 - missing_freq;
  double scaled_missing_freq = missing_freq * TWO_63;
  double scaled_nonmissing_freq = nonmissing_freq * TWO_63;

  // P(causal variant genotype)
  double g0 = freq * freq;
  double g1 = 2 * freq * (1 - freq);
  double g2 = 1 - g0 - g1;
  // P(marker genotype)
  double mg0 = mfreq * mfreq;
  double mg1 = 2 * mfreq * (1 - mfreq);
  double mg2 = 1 - mg0 - mg1;

  double ld = freq * (1 - mfreq);
  double h21 = (1 - freq) * mfreq;
  double h11;
  double h12;
  double h22;

  double f0;
  double f1;
  double f2;
  double mf0;
  double mf1;
  double mf2;
  // first 2 digits indicate causal variant genotype, last 2 indicate marker.
  // this is DIFFERENT from instanceSNP()'s usage.
  double h_11_11;
  double h_11_12;
  double h_11_22;
  double h_12_11;
  double h_12_12;
  double h_12_22;
  double h_22_11;
  double h_22_12;
  double h_22_22;

  double a0;
  double a1;
  double a2;
  double u0;
  double u1;
  double u2;

  double xh_11_11;
  double xh_11_12;
  double xh_11_22;
  double xh_12_11;
  double xh_12_12;
  double xh_12_22;
  double xh_22_11;
  double xh_22_12;
  double xh_22_22;

  double tot_recip;
  double tot_recip_nmsq;
  double tot_recip_1miss;

  if (h21 < ld) {
    ld = h21;
  }
  ld *= dprime;

  // Joint causal variant/marker probs, considering one allele at a time.
  h11 = freq * mfreq + ld;
  h12 = freq * (1 - mfreq) - ld;
  h21 -= ld;
  h22 = (1 - freq) * (1 - mfreq) + ld;


  // Now considering both alleles simultaneously.
  h_11_11 = h11 * h11;
  h_11_12 = h11 * h12 * 2;
  h_11_22 = h12 * h12;
  h_12_11 = h21 * h11 * 2;
  h_12_12 = (h22 * h11 + h21 * h12) * 2;
  h_12_22 = h22 * h12 * 2;
  h_22_11 = h21 * h21;
  h_22_12 = h22 * h21 * 2;
  h_22_22 = h22 * h22;

  // P(case | causal variant genotype)
  simulate_cc_get_conditional_probs(prevalence, g0, g1, g2, het_odds, hom0_odds, &f0, &f1, &f2);
  // P(case | marker genotype)
  mf0 = (f0 * h_11_11 + f1 * h_12_11 + f2 * h_22_11) / mg0;
  mf1 = (f0 * h_11_12 + f1 * h_12_12 + f2 * h_22_12) / mg1;
  mf2 = (f0 * h_11_22 + f1 * h_12_22 + f2 * h_22_22) / mg2;
  // P(marker genotype | affection status)
  a0 = mg0 * mf0;
  a1 = mg1 * mf1;
  a2 = mg2 * mf2;
  tot_recip = 1.0 / (a0 + a1 + a2);
  a0 *= tot_recip;
  a1 *= tot_recip;
  a2 *= tot_recip;
  u0 = mg0 * (1 - mf0);
  u1 = mg1 * (1 - mf1);
  u2 = mg2 * (1 - mf2);
  tot_recip = 1.0 / (u0 + u1 + u2);
  u0 *= tot_recip;
  u1 *= tot_recip;
  u2 *= tot_recip;

  if (!do_haps) {
    case_thresholds[0] = (uint64_t)(a0 * scaled_nonmissing_freq);
    case_thresholds[1] = case_thresholds[0] + (uint64_t)(scaled_missing_freq);
    case_thresholds[2] = case_thresholds[1] + (uint64_t)(a1 * scaled_nonmissing_freq);
    ctrl_thresholds[0] = (uint64_t)(u0 * scaled_nonmissing_freq);
    ctrl_thresholds[1] = ctrl_thresholds[0] + (uint64_t)(scaled_missing_freq);
    ctrl_thresholds[2] = ctrl_thresholds[1] + (uint64_t)(u1 * scaled_nonmissing_freq);
  } else {
    xh_11_11 = h_11_11 * f0;
    xh_11_12 = h_11_12 * f0;
    xh_11_22 = h_11_22 * f0;
    xh_12_11 = h_12_11 * f1;
    xh_12_12 = h_12_12 * f1;
    xh_12_22 = h_12_22 * f1;
    xh_22_11 = h_22_11 * f2;
    xh_22_12 = h_22_12 * f2;
    xh_22_22 = h_22_22 * f2;
    tot_recip = 1.0 / (xh_11_11 + xh_11_12 + xh_11_22 + xh_12_11 + xh_12_12 + xh_12_22 + xh_22_11 + xh_22_12 + xh_22_22);
    tot_recip_nmsq = tot_recip * nonmissing_freq * scaled_nonmissing_freq;
    tot_recip_1miss = tot_recip * nonmissing_freq * scaled_missing_freq;

    case_thresholds[0] = (uint64_t)(xh_11_11 * tot_recip_nmsq);
    case_thresholds[1] = case_thresholds[0] + (uint64_t)((xh_11_11 + xh_11_12 + xh_11_22) * tot_recip_1miss);
    case_thresholds[2] = case_thresholds[1] + (uint64_t)(xh_11_12 * tot_recip_nmsq);
    case_thresholds[3] = case_thresholds[2] + (uint64_t)(xh_11_22 * tot_recip_nmsq);
    case_thresholds[4] = case_thresholds[3] + (uint64_t)((xh_11_11 + xh_12_11 + xh_22_11) * tot_recip_1miss);
    case_thresholds[5] = case_thresholds[4] + (uint64_t)(missing_freq * scaled_missing_freq);
    case_thresholds[6] = case_thresholds[5] + (uint64_t)((xh_11_12 + xh_12_12 + xh_22_12) * tot_recip_1miss);
    case_thresholds[7] = case_thresholds[6] + (uint64_t)((xh_11_22 + xh_12_22 + xh_22_22) * tot_recip_1miss);
    case_thresholds[8] = case_thresholds[7] + (uint64_t)(xh_12_11 * tot_recip_nmsq);
    case_thresholds[9] = case_thresholds[8] + (uint64_t)((xh_12_11 + xh_12_12 + xh_12_22) * tot_recip_1miss);
    case_thresholds[10] = case_thresholds[9] + (uint64_t)(xh_12_12 * tot_recip_nmsq);
    case_thresholds[11] = case_thresholds[10] + (uint64_t)(xh_12_22 * tot_recip_nmsq);
    case_thresholds[12] = case_thresholds[11] + (uint64_t)(xh_22_11 * tot_recip_nmsq);
    case_thresholds[13] = case_thresholds[12] + (uint64_t)((xh_22_11 + xh_22_12 + xh_22_22) * tot_recip_1miss);
    case_thresholds[14] = case_thresholds[13] + (uint64_t)(xh_22_12 * tot_recip_nmsq);

    xh_11_11 = h_11_11 * (1 - f0);
    xh_11_12 = h_11_12 * (1 - f0);
    xh_11_22 = h_11_22 * (1 - f0);
    xh_12_11 = h_12_11 * (1 - f1);
    xh_12_12 = h_12_12 * (1 - f1);
    xh_12_22 = h_12_22 * (1 - f1);
    xh_22_11 = h_22_11 * (1 - f2);
    xh_22_12 = h_22_12 * (1 - f2);
    xh_22_22 = h_22_22 * (1 - f2);
    tot_recip = 1.0 / (xh_11_11 + xh_11_12 + xh_11_22 + xh_12_11 + xh_12_12 + xh_12_22 + xh_22_11 + xh_22_12 + xh_22_22);
    tot_recip_nmsq = tot_recip * nonmissing_freq * scaled_nonmissing_freq;
    tot_recip_1miss = tot_recip * nonmissing_freq * scaled_missing_freq;

    ctrl_thresholds[0] = (uint64_t)(xh_11_11 * tot_recip_nmsq);
    ctrl_thresholds[1] = ctrl_thresholds[0] + (uint64_t)((xh_11_11 + xh_11_12 + xh_11_22) * tot_recip_1miss);
    ctrl_thresholds[2] = ctrl_thresholds[1] + (uint64_t)(xh_11_12 * tot_recip_nmsq);
    ctrl_thresholds[3] = ctrl_thresholds[2] + (uint64_t)(xh_11_22 * tot_recip_nmsq);
    ctrl_thresholds[4] = ctrl_thresholds[3] + (uint64_t)((xh_11_11 + xh_12_11 + xh_22_11) * tot_recip_1miss);
    ctrl_thresholds[5] = ctrl_thresholds[4] + (uint64_t)(missing_freq * scaled_missing_freq);
    ctrl_thresholds[6] = ctrl_thresholds[5] + (uint64_t)((xh_11_12 + xh_12_12 + xh_22_12) * tot_recip_1miss);
    ctrl_thresholds[7] = ctrl_thresholds[6] + (uint64_t)((xh_11_22 + xh_12_22 + xh_22_22) * tot_recip_1miss);
    ctrl_thresholds[8] = ctrl_thresholds[7] + (uint64_t)(xh_12_11 * tot_recip_nmsq);
    ctrl_thresholds[9] = ctrl_thresholds[8] + (uint64_t)((xh_12_11 + xh_12_12 + xh_12_22) * tot_recip_1miss);
    ctrl_thresholds[10] = ctrl_thresholds[9] + (uint64_t)(xh_12_12 * tot_recip_nmsq);
    ctrl_thresholds[11] = ctrl_thresholds[10] + (uint64_t)(xh_12_22 * tot_recip_nmsq);
    ctrl_thresholds[12] = ctrl_thresholds[11] + (uint64_t)(xh_22_11 * tot_recip_nmsq);
    ctrl_thresholds[13] = ctrl_thresholds[12] + (uint64_t)((xh_22_11 + xh_22_12 + xh_22_22) * tot_recip_1miss);
    ctrl_thresholds[14] = ctrl_thresholds[13] + (uint64_t)(xh_22_12 * tot_recip_nmsq);
  }
}

int32_t simulate_dataset(char* outname, char* outname_end, uint32_t flags, char* simulate_fname, uint32_t case_ct, uint32_t ctrl_ct, double prevalence, uint32_t indiv_ct, double missing_freq, char* name_prefix) {
  FILE* infile = NULL;
  FILE* outfile_txt = NULL;
  FILE* outfile_simfreq = NULL;
  FILE* outfile_bed = NULL;
  unsigned char* wkspace_mark = wkspace_base;
  int32_t retval = 0;
  uint32_t do_tags = flags & SIMULATE_TAGS;
  uint32_t do_haps = flags & SIMULATE_HAPS;
  uint32_t tags_or_haps = do_tags | do_haps;
  uint32_t is_qt = flags & SIMULATE_QT;
  uint32_t randomize_alleles = flags & (SIMULATE_ACGT | SIMULATE_1234 | SIMULATE_12);
  uint32_t simulate_12 = flags & SIMULATE_12;
  uint32_t marker_pos = 1;
  double* qt_vals = NULL;
  double dxx = 0;
  double dyy = 0;
  double qt_totvar = 0;
  uint32_t cur_marker_ct = 0;
  char* cur_snp_label = &(tbuf[MAXLINELEN]);
  uint32_t snp_label_len = 0;
  uint32_t missing_thresh = (uint32_t)(missing_freq * 4294967296.0);
  double marker_freq_lb = 0;
  double marker_freq_ub = 0;
  double dprime = 0;
  double het_odds = 0;
  double hom0_odds = 0;
  double qt_var = 0;
  double qt_dom = 0;
  uint32_t marker_idx_offset = 0;
  uint32_t pct = 0;
  uint32_t name_prefix_len = 0;
  uint32_t zero_odds_ratio_warning_given = 0;
  char* marker_freq_lb_ptr = NULL;
  char* marker_ld_ptr = NULL;
  uintptr_t* writebuf2 = NULL;
  char alleles[13];
  char cur_alleles[4];
  double freqs[6];
  uint64_t thresholds[15];
  uint64_t case_thresholds[15];
  double qt_adj[4];
  uint32_t marker_ct;
  uintptr_t* writebuf;
  char* cptr;
  char* snp_label_ptr;
  char* freq_lb_ptr;
  char* penult_ptr; // het disease odds for C/C, additive genetic var for QT
  char* last_ptr; // homset disease odds for C/C, dominance/additive for QT
  char* wptr;
  uintptr_t* wbptr;
  uintptr_t* wbptr2;
  uintptr_t indiv_ct4;
  uintptr_t indiv_ctl2;
  sfmt_t* sfmt64p;
  uint32_t indiv_idx;
  uint32_t cur_marker_idx;
  uint32_t loop_end;
  double freq_lb;
  double freq_delta;
  int32_t ii;
  char cc;
  uint64_t ullii;
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t ulkk;
  uintptr_t ulmm;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  double dzz;
  if (!is_qt) {
    indiv_ct = case_ct + ctrl_ct;
  } else {
    if (wkspace_alloc_d_checked(&qt_vals, indiv_ct * sizeof(double))) {
      goto simulate_ret_NOMEM;
    }
    fill_double_zero(qt_vals, indiv_ct);
  }
  indiv_ct4 = (indiv_ct + 3) / 4;
  indiv_ctl2 = (indiv_ct + BITCT2 - 1) / BITCT2;
  if (randomize_alleles) {
    if (flags & SIMULATE_ACGT) {
      memcpy(alleles, "ACAGATCGCTGTA", 13);
    } else if (flags & SIMULATE_1234) {
      memcpy(alleles, "1213142324341", 13);
    } else {
      memcpyl3(alleles, "121");
    }
  } else {
    if (is_qt) {
      memcpy(alleles, "HLAB", 4);
    } else {
      memcpy(alleles, "DdAB", 4);
    }
  }
  if (wkspace_alloc_ul_checked(&writebuf, indiv_ctl2 * sizeof(intptr_t))) {
    goto simulate_ret_NOMEM;
  }
  if (do_haps) {
    if (wkspace_alloc_ul_checked(&writebuf2, indiv_ctl2 * sizeof(intptr_t))) {
      goto simulate_ret_NOMEM;
    }
  }
  if (fopen_checked(&infile, simulate_fname, "r")) {
    goto simulate_ret_OPEN_FAIL;
  }
  memcpy(outname_end, ".bim", 5);
  if (fopen_checked(&outfile_txt, outname, "w")) {
    goto simulate_ret_OPEN_FAIL;
  }
  memcpy(outname_end, ".simfreq", 9);
  if (fopen_checked(&outfile_simfreq, outname, "w")) {
    goto simulate_ret_OPEN_FAIL;
  }
  memcpy(outname_end, ".bed", 5);
  if (fopen_checked(&outfile_bed, outname, "wb")) {
    goto simulate_ret_OPEN_FAIL;
  }
  if (fwrite_checked("l\x1b\x01", 3, outfile_bed)) {
    goto simulate_ret_WRITE_FAIL;
  }
  sprintf(logbuf, "Writing --simulate%s dataset to %s + .bim + .fam...", is_qt? "-qt" : "", outname);
  logprintb();
  fputs(" 0%", stdout);
  fflush(stdout);
  sfmt64p = (sfmt_t*)wkspace_alloc(sizeof(sfmt_t));
  if (!sfmt64p) {
    goto simulate_ret_NOMEM;
  }
  init_sfmt64_from_sfmt32(&sfmt, sfmt64p);
  tbuf[MAXLINELEN - 1] = ' ';
  // just determine total marker ct in initial scan, for progress indicator
  ullii = 0;
  while (fgets(tbuf, MAXLINELEN, infile)) {
    if (!tbuf[MAXLINELEN - 1]) {
      sprintf(logbuf, "\nError: Excessively long line in --simulate%s input file.\n", is_qt? "-qt" : "");
      goto simulate_ret_INVALID_FORMAT_2;
    }
    cptr = skip_initial_spaces(tbuf);
    if (is_eoln_kns(*cptr)) {
      continue;
    }
    if (atoiz(cptr, &ii)) {
      sprintf(logbuf, "\nError: Invalid SNP count in --simulate%s input file.\n", is_qt? "-qt" : "");
      goto simulate_ret_INVALID_FORMAT_2;
    }
    ullii += ii;
  }
  if (!feof(infile)) {
    goto simulate_ret_READ_FAIL;
  }
  if (!ullii) {
    sprintf(logbuf, "\nError: --simulate%s input file specifies zero SNPs.\n", is_qt? "-qt" : "");
    goto simulate_ret_INVALID_FORMAT_2;
  } else if (ullii > (do_haps? 0x3fffffff : 0x7fffffff)) {
    sprintf(logbuf, "\nError: --simulate%s input file specifies too many SNPs.\n", is_qt? "-qt" : "");
    goto simulate_ret_INVALID_FORMAT_2;
  }
  marker_ct = ullii;
  loop_end = (marker_ct + 99) / 100;
  rewind(infile);
  while (fgets(tbuf, MAXLINELEN, infile)) {
    cptr = skip_initial_spaces(tbuf);
    if (is_eoln_kns(*cptr)) {
      continue;
    }
    snp_label_ptr = next_item(cptr);
    freq_lb_ptr = next_item(snp_label_ptr);
    if (tags_or_haps) {
      marker_freq_lb_ptr = next_item_mult(freq_lb_ptr, 2);
      marker_ld_ptr = next_item_mult(marker_freq_lb_ptr, 2);
      penult_ptr = next_item(marker_ld_ptr);
    } else {
      penult_ptr = next_item_mult(freq_lb_ptr, 2);
    }
    last_ptr = next_item(penult_ptr);
    if (no_more_items(last_ptr)) {
      sprintf(logbuf, "\nError: Missing token(s) in --simulate%s input file line.\n", is_qt? "-qt" : "");
      goto simulate_ret_INVALID_FORMAT_2;
    }
    if (!no_more_items(next_item(last_ptr))) {
      sprintf(logbuf, "\nError: Too many tokens in --simulate%s input file line.\n", is_qt? "-qt" : "");
      goto simulate_ret_INVALID_FORMAT_2;
    }
    atoiz(cptr, &ii);
    if (!ii) {
      continue;
    }
    cur_marker_ct = ii;
    snp_label_len = strlen_se(snp_label_ptr);
    memcpy(cur_snp_label, snp_label_ptr, snp_label_len);
    cur_snp_label[snp_label_len++] = '_';
    if (scan_two_doubles(freq_lb_ptr, &freq_lb, &freq_delta) || (freq_lb < 0) || (freq_delta < freq_lb) || (freq_delta > 1)) {
      sprintf(logbuf, "\nError: Invalid allele frequency bound in --simulate%s input file.\n", is_qt? "-qt" : "");
      goto simulate_ret_INVALID_FORMAT_2;
    }
    freq_delta -= freq_lb;
    if (tags_or_haps) {
      if (scan_two_doubles(marker_freq_lb_ptr, &marker_freq_lb, &marker_freq_ub) || (marker_freq_lb < 0) || (marker_freq_ub < marker_freq_lb) || (marker_freq_ub > 1)) {
	sprintf(logbuf, "\nError: Invalid marker allele frequency bound in --simulate%s input file.\n", is_qt? "-qt" : "");
	goto simulate_ret_INVALID_FORMAT_2;
      }
      if (scan_double(marker_ld_ptr, &dprime) || (dprime < 0) || (dprime > 1)) {
	sprintf(logbuf, "\nError: Invalid d-prime in --simulate%s input file.\n", is_qt? "-qt" : "");
	goto simulate_ret_INVALID_FORMAT_2;
      }
    } else {
      dprime = 1;
    }
    if (is_qt) {
      if (scan_double(penult_ptr, &qt_var) || (qt_var < 0) || (qt_var > 1)) {
	logprint("\nError: Invalid variance value in --simulate-qt input file.\n");
	goto simulate_ret_INVALID_FORMAT;
      }
      if ((qt_var > 0) && (((freq_delta == 0) && ((freq_lb == 0) || (freq_lb == 1))) || (tags_or_haps && (marker_freq_lb == marker_freq_ub) && ((marker_freq_lb == 0) || (marker_freq_lb == 1))))) {
	logprint("\nError: Nonzero variance with fixed 0/1 allele frequency in --simulate-qt input file.\n");
	goto simulate_ret_INVALID_FORMAT;
      }
      qt_totvar += ((intptr_t)cur_marker_ct) * qt_var;
      if (qt_totvar > 1 + EPSILON) {
	logprint("\nError: --simulate-qt input file specific QTL variance greater than 1.\n");
	goto simulate_ret_INVALID_FORMAT;
      }
      if (scan_double(last_ptr, &qt_dom)) {
	logprint("\nError: Invalid dominance deviation value in --simulate-qt input file.\n");
	goto simulate_ret_INVALID_FORMAT;
      }
    } else {
      if (scan_double(penult_ptr, &het_odds) || (het_odds < 0)) {
	logprint("\nError: Invalid heterozygote disease odds ratio in --simulate input file.\n");
	goto simulate_ret_INVALID_FORMAT;
      }
      if ((strlen_se(last_ptr) == 4) && match_upper_nt(last_ptr, "MULT", 4)) {
	hom0_odds = het_odds * het_odds;
      } else if (scan_double(last_ptr, &hom0_odds) || (hom0_odds < 0)) {
	logprint("\nError: Invalid homozygote disease odds ratio in --simulate input file.\n");
	goto simulate_ret_INVALID_FORMAT;
      }
      if ((!zero_odds_ratio_warning_given) && ((het_odds == 0) || (hom0_odds == 0))) {
        putchar('\r');
	logstr("\n");
	logprint("Warning: Zero odds ratio present in --simulate input file.  Did you mean\n--simulate-qt instead?\n");
	zero_odds_ratio_warning_given = 1;
        printf("%u%%", pct);
      }
    }
    tbuf[0] = '1';
    for (cur_marker_idx = 0; cur_marker_idx < cur_marker_ct; cur_marker_idx++) {
      freqs[0] = freq_lb + rand_unif() * freq_delta;
      if (tags_or_haps) {
	// force marker freq/c.v. freq and (1 - marker freq)/(1 - c.v. freq) to
	// not be smaller than dprime or larger than 1 / dprime, unless that
	// causes marker freq to be out of range.
	if (dprime > 0) {
	  dxx = 1 - (1 - freqs[0]) / dprime;
	  dzz = freqs[0] * dprime;
	  if (dxx < dzz) {
	    dxx = dzz;
	  }
	  if (dxx < marker_freq_lb) {
	    dxx = marker_freq_lb;
	  }
	  dyy = 1 - (1 - freqs[0]) * dprime;
	  dzz = freqs[0] / dprime;
	  if (dyy > dzz) {
	    dyy = dzz;
	  }
	  if (dyy > marker_freq_ub) {
	    dyy = marker_freq_ub;
	  }
	  if (dyy < dxx) {
	    if (dyy < marker_freq_lb) {
	      dyy = marker_freq_lb;
	    }
	    dxx = dyy;
	  }
	} else {
	  dxx = marker_freq_lb;
	  dyy = marker_freq_ub;
	}
	freqs[1] = dxx + rand_unif() * (dyy - dxx);
      } else {
	freqs[1] = freqs[0];
      }
      if (is_qt) {
	simulate_init_freqs_qt(do_haps, dprime, qt_var, qt_dom, missing_freq, freqs, thresholds, qt_adj);
      } else {
	simulate_init_freqs_cc(do_haps, dprime, freqs, prevalence, het_odds, hom0_odds, missing_freq, thresholds, case_thresholds);
      }
      wptr = &(tbuf[1]);
      *wptr++ = ' ';
      if (cur_marker_ct > 1) {
	wptr = memcpya(wptr, cur_snp_label, snp_label_len);
	wptr = uint32_write(wptr, cur_marker_idx);
      } else {
	wptr = memcpya(wptr, cur_snp_label, snp_label_len - 1);
      }
      *wptr++ = '\t';
      dxx = freqs[0];
      wptr = double_g_writex(wptr, dxx, ' ');
      wptr = double_g_writex(wptr, dxx, '\t');
      if (tags_or_haps) {
	dxx = freqs[1];
	wptr = double_g_writex(wptr, dxx, ' ');
	wptr = double_g_writex(wptr, dxx, '\t');
	wptr = double_g_writex(wptr, dprime, '\t');
      }
      if (is_qt) {
	wptr = double_g_writex(wptr, qt_var, '\t');
	wptr = double_g_writex(wptr, qt_dom, '\n');
      } else {
	wptr = double_g_writex(wptr, het_odds, '\t');
	wptr = double_g_writex(wptr, hom0_odds, '\n');
      }
      if (fwrite_checked(tbuf, (uintptr_t)(wptr - tbuf), outfile_simfreq)) {
	goto simulate_ret_WRITE_FAIL;
      }
      if (randomize_alleles) {
	if (!simulate_12) {
	  do {
	    uii = sfmt_genrand_uint32(&sfmt);
	  } while (uii >= 4294967184U); // largest multiple of 144 < 2^32
	  uii = uii % 144U;
	  ujj = uii / 12;
	  uii -= ujj * 12;
	} else {
	  uii = sfmt_genrand_uint32(&sfmt) & 3;
	  ujj = uii >> 1;
	  uii &= 1;
	}
	memcpy(cur_alleles, &(alleles[uii]), 2);
	memcpy(&(cur_alleles[2]), &(alleles[ujj]), 2);
      } else {
	memcpy(cur_alleles, alleles, 4);
      }
      ulii = 0;
      ukk = 0;
      wbptr = writebuf;
      if (!do_haps) {
	if (is_qt) {
	  for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
	    ullii = sfmt_genrand_uint64(sfmt64p) >> 1;
	    if (ullii > thresholds[1]) {
	      ulkk = 3;
	    } else if (ullii > thresholds[0]) {
	      ulkk = 2;
	    } else {
	      ulkk = 0;
	    }
	    qt_vals[indiv_idx] += qt_adj[ulkk];
	    if (sfmt_genrand_uint32(&sfmt) < missing_thresh) {
	      ulkk = 1;
	    }
	    ulii |= ulkk << ukk;
	    ukk += 2;
	    if (ukk == BITCT) {
	      *wbptr++ = ulii;
	      ulii = 0;
	      ukk = 0;
	    }
	  }
	} else {
	  for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
	    ullii = sfmt_genrand_uint64(sfmt64p) >> 1;
	    if (indiv_idx < case_ct) {
	      if (ullii > case_thresholds[1]) {
		if (ullii > case_thresholds[2]) {
		  ulkk = 3;
		} else {
		  ulkk = 2;
		}
	      } else if (ullii > case_thresholds[0]) {
		ulkk = 1;
	      } else {
		ulkk = 0;
	      }
	    } else {
	      if (ullii > thresholds[1]) {
		if (ullii > thresholds[2]) {
		  ulkk = 3;
		} else {
		  ulkk = 2;
		}
	      } else if (ullii > thresholds[0]) {
		ulkk = 1;
	      } else {
		ulkk = 0;
	      }
	    }
	    ulii |= ulkk << ukk;
	    ukk += 2;
	    if (ukk == BITCT) {
	      *wbptr++ = ulii;
	      ulii = 0;
	      ukk = 0;
	    }
	  }
	}
	if (ukk) {
	  *wbptr = ulii;
	}
      } else {
	uljj = 0;
	wbptr2 = writebuf2;
	if (is_qt) {
	  for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
	    ullii = sfmt_genrand_uint64(sfmt64p) >> 1;
	    ulkk = uint64arr_greater_than(thresholds, 11, ullii);
	    ulmm = ulkk & 3;
	    ulkk /= 4;
	    ulkk += (ulkk + 1) >> 1;
	    qt_vals[indiv_idx] += qt_adj[ulkk];
	    if (sfmt_genrand_uint32(&sfmt) < missing_thresh) {
	      ulkk = 1;
	    }
	    ulii |= ulkk << ukk;
	    uljj |= ulmm << ukk;
	    ukk += 2;
	    if (ukk == BITCT) {
	      *wbptr++ = ulii;
	      *wbptr2++ = uljj;
	      ulii = 0;
	      uljj = 0;
	      ukk = 0;
	    }
	  }
	} else {
	  for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
	    ullii = sfmt_genrand_uint64(sfmt64p) >> 1;
	    if (indiv_idx < case_ct) {
	      ulkk = uint64arr_greater_than(case_thresholds, 15, ullii);
	    } else {
	      ulkk = uint64arr_greater_than(thresholds, 15, ullii);
	    }
	    ulmm = ulkk & 3;
	    ulkk /= 4;
	    ulii |= ulkk << ukk;
	    uljj |= ulmm << ukk;
	    ukk += 2;
	    if (ukk == BITCT) {
	      *wbptr++ = ulii;
	      *wbptr2++ = uljj;
	      ulii = 0;
	      uljj = 0;
	      ukk = 0;
	    }
	  }
	}
	if (ukk) {
	  *wbptr = ulii;
	  *wbptr2 = uljj;
	}
      }
      if (popcount_longs(writebuf, 0, indiv_ctl2) < indiv_ct) {
	reverse_loadbuf((unsigned char*)writebuf, indiv_ct);
	cc = cur_alleles[0];
	cur_alleles[0] = cur_alleles[1];
	cur_alleles[1] = cc;
      }
      wptr = &(tbuf[1]);
      *wptr++ = '\t';
      if (cur_marker_ct > 1) {
	wptr = memcpya(wptr, cur_snp_label, snp_label_len);
	wptr = uint32_write(wptr, cur_marker_idx);
      } else {
	wptr = memcpya(wptr, cur_snp_label, snp_label_len - 1);
      }
      if (do_tags) {
	wptr = memcpya(wptr, "_M", 2);
      }
      wptr = memcpyl3a(wptr, "\t0\t");
      wptr = uint32_writex(wptr, marker_pos++, '\t');
      *wptr++ = cur_alleles[0];
      *wptr++ = '\t';
      *wptr++ = cur_alleles[1];
      *wptr++ = '\n';
      if (fwrite_checked(tbuf, (uintptr_t)(wptr - tbuf), outfile_txt)) {
	goto simulate_ret_WRITE_FAIL;
      }
      if (fwrite_checked((unsigned char*)writebuf, indiv_ct4, outfile_bed)) {
	goto simulate_ret_WRITE_FAIL;
      }
      if (do_haps) {
	if (popcount_longs(writebuf2, 0, indiv_ctl2) < indiv_ct) {
	  reverse_loadbuf((unsigned char*)writebuf2, indiv_ct);
	  cc = cur_alleles[2];
	  cur_alleles[2] = cur_alleles[3];
	  cur_alleles[3] = cc;
	}
	wptr = &(tbuf[2 + snp_label_len]);
	wptr = uint32_write(wptr, cur_marker_idx);
	wptr = memcpya(wptr, "_M\t0\t", 5);
	wptr = uint32_writex(wptr, marker_pos++, '\t');
	*wptr++ = cur_alleles[2];
	*wptr++ = '\t';
	*wptr++ = cur_alleles[3];
	*wptr++ = '\n';
	if (fwrite_checked(tbuf, (uintptr_t)(wptr - tbuf), outfile_txt)) {
	  goto simulate_ret_WRITE_FAIL;
	}
	if (fwrite_checked((unsigned char*)writebuf2, indiv_ct4, outfile_bed)) {
	  goto simulate_ret_WRITE_FAIL;
	}
      }
      if (cur_marker_idx >= loop_end) {
	if (pct > 9) {
	  putchar('\b');
	}
	pct = ((cur_marker_idx + marker_idx_offset) * 100LLU) / marker_ct;
	printf("\b\b%u%%", pct);
	fflush(stdout);
	loop_end = ((marker_ct * (pct + 1LLU) + 99LLU) / 100U) - marker_idx_offset;
      }
    }
    marker_idx_offset += cur_marker_ct;
    loop_end -= cur_marker_ct;
  }
  if (!feof(infile)) {
    goto simulate_ret_READ_FAIL;
  }
  if (fclose_null(&outfile_txt) || fclose_null(&outfile_bed)) {
    goto simulate_ret_WRITE_FAIL;
  }
  memcpy(outname_end, ".fam", 5);
  if (fopen_checked(&outfile_txt, outname, "w")) {
    goto simulate_ret_OPEN_FAIL;
  }
  wptr = tbuf;
  if (name_prefix) {
    name_prefix_len = strlen(name_prefix);
    wptr = memcpyax(wptr, name_prefix, name_prefix_len, '-');
    uii = 4 + name_prefix_len;
  } else {
    uii = 3;
  }
  memcpyl3(wptr, "per");
  if (is_qt) {
    if (qt_totvar < 1 - EPSILON) {
      dyy = sqrt(1 - qt_totvar);
    } else {
      dyy = 0;
    }
  }
  for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
    wptr = uint32_writex(&(tbuf[uii]), indiv_idx, ' ');
    if (name_prefix_len) {
      wptr = memcpyax(wptr, name_prefix, name_prefix_len, '-');
    }
    wptr = memcpyl3a(wptr, "per");
    wptr = uint32_write(wptr, indiv_idx);
    wptr = memcpya(wptr, " 0 0 2 ", 7);
    if (is_qt) {
      if (indiv_idx & 1) {
	dzz = qt_vals[indiv_idx] + dyy * dxx;
      } else {
	dzz = qt_vals[indiv_idx] + dyy * rand_normal(&dxx);
      }
      wptr = double_g_write(wptr, dzz);
    } else {
      if (indiv_idx < case_ct) {
	*wptr++ = '2';
      } else {
	*wptr++ = '1';
      }
    }
    *wptr++ = '\n';
    if (fwrite_checked(tbuf, (uintptr_t)(wptr - tbuf), outfile_txt)) {
      goto simulate_ret_WRITE_FAIL;
    }
  }
  *outname_end = '\0';
  if (pct > 9) {
    putchar('\b');
  }
  fputs("\b\b", stdout);
  if (!zero_odds_ratio_warning_given) {
    logstr(" ");
  }
  sprintf(logbuf, "done.\nRealized simulation parameters saved to %s.simfreq.\n", outname);
  logprintb();
  while (0) {
  simulate_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  simulate_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  simulate_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  simulate_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  simulate_ret_INVALID_FORMAT_2:
    logprintb();
  simulate_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
  fclose_cond(infile);
  fclose_cond(outfile_txt);
  fclose_cond(outfile_simfreq);
  fclose_cond(outfile_bed);
  wkspace_reset(wkspace_mark);
  return retval;
}

int32_t recode_allele_load(char* loadbuf, uintptr_t loadbuf_size, char* recode_allele_name, char*** allele_missing_ptr, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, uintptr_t* recode_allele_reverse, char* recode_allele_extra) {
  FILE* rafile = NULL;
  uint32_t missing_allele = 0;
  uintptr_t rae_size = 0;
  char* sorted_ids;
  uint32_t* id_map;
  char* bufptr;
  char* bufptr2;
  int32_t retval;
  uint32_t slen;
  uint32_t alen;
  int32_t ii;
  uintptr_t marker_uidx;
  if (fopen_checked(&rafile, recode_allele_name, "r")) {
    goto recode_allele_load_ret_OPEN_FAIL;
  }
  retval = sort_item_ids(&sorted_ids, &id_map, unfiltered_marker_ct, marker_exclude, unfiltered_marker_ct - marker_ct, marker_ids, max_marker_id_len, 0, 0, strcmp_deref);
  if (retval) {
    goto recode_allele_load_ret_1;
  }
  loadbuf[loadbuf_size - 1] = ' ';
  while (fgets(loadbuf, loadbuf_size, rafile)) {
    if (!loadbuf[loadbuf_size - 1]) {
      logprint("Error: Pathologically long line in --recode-allele file.\n");
      goto recode_allele_load_ret_INVALID_FORMAT;
    }
    bufptr = skip_initial_spaces(loadbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    slen = strlen_se(bufptr);
    bufptr2 = skip_initial_spaces(&(bufptr[slen]));
    if (is_eoln_kns(*bufptr2)) {
      logprint("Error: --recode-allele line has only one token.\n");
      goto recode_allele_load_ret_INVALID_FORMAT;
    }
    alen = strlen_se(bufptr2);
    bufptr[slen] = '\0';
    ii = bsearch_str(bufptr, sorted_ids, max_marker_id_len, 0, marker_ct - 1);
    if (ii != -1) {
      marker_uidx = id_map[(uint32_t)ii];
      bufptr2[alen++] = '\0';
      if (!strcmp(bufptr2, marker_allele_ptrs[2 * marker_uidx])) {
	CLEAR_BIT(recode_allele_reverse, marker_uidx);
      } else if (!strcmp(bufptr2, marker_allele_ptrs[2 * marker_uidx + 1])) {
	SET_BIT(recode_allele_reverse, marker_uidx);
      } else {
	if (rae_size + alen > wkspace_left) {
	  goto recode_allele_load_ret_NOMEM;
	}
	missing_allele = 1;
	(*allele_missing_ptr)[marker_uidx] = &(recode_allele_extra[rae_size]);
	memcpy(&(recode_allele_extra[rae_size]), bufptr2, alen);
	rae_size += alen;
      }
    }
  }
  if (!feof(rafile)) {
    goto recode_allele_load_ret_READ_FAIL;
  }
  while (0) {
  recode_allele_load_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  recode_allele_load_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  recode_allele_load_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  recode_allele_load_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
  }
 recode_allele_load_ret_1:
  fclose_cond(rafile);
  if (missing_allele) {
    recode_allele_extra = (char*)wkspace_alloc(rae_size);
  } else {
    wkspace_reset((unsigned char*)(*allele_missing_ptr));
    *allele_missing_ptr = NULL;
  }
  return retval;
}

uint32_t recode_load_to(unsigned char* loadbuf, FILE* bedfile, uintptr_t bed_offset, uintptr_t unfiltered_marker_ct, uintptr_t marker_idx, uintptr_t marker_idx_end, uintptr_t* marker_exclude, uintptr_t* marker_reverse, uintptr_t* marker_uidx_ptr, uintptr_t unfiltered_indiv_ct) {
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t marker_uidx = *marker_uidx_ptr;
  uintptr_t marker_uidx_start;
  uintptr_t marker_uidx_stop;
  uintptr_t ulii;
  if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
    return 1;
  }
  while (marker_idx < marker_idx_end) {
    if (IS_SET(marker_exclude, marker_uidx)) {
      marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
      if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	return 1;
      }
    }
    if (unfiltered_marker_ct - marker_uidx > marker_idx_end - marker_idx) {
      ulii = next_set_ul_unsafe(marker_exclude, marker_uidx) - marker_uidx;
    } else {
      ulii = unfiltered_marker_ct - marker_uidx;
    }
    marker_uidx_start = marker_uidx;
    marker_uidx_stop = marker_uidx + ulii;
    marker_idx += ulii;
    ulii *= unfiltered_indiv_ct4;
    if (fread(loadbuf, 1, ulii, bedfile) < ulii) {
      return 1;
    }
    while (1) {
      next_set_ul_ck(marker_reverse, &marker_uidx, marker_uidx_stop);
      if (marker_uidx == marker_uidx_stop) {
	break;
      }
      reverse_loadbuf(&(loadbuf[(marker_uidx - marker_uidx_start) * unfiltered_indiv_ct4]), unfiltered_indiv_ct);
      marker_uidx++;
    }
    loadbuf = &(loadbuf[ulii]);
  }
  *marker_uidx_ptr = marker_uidx;
  return 0;
}

static inline int32_t recode_write_first_cols(FILE* outfile, uintptr_t indiv_uidx, char delimiter, char* person_ids, uintptr_t max_person_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* sex_nm, uintptr_t* sex_male, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, char* output_missing_pheno) {
  char wbuf[16];
  char* cptr = &(person_ids[indiv_uidx * max_person_id_len]);
  uintptr_t ulii = strlen_se(cptr);
  fwrite(cptr, 1, ulii, outfile);
  putc(delimiter, outfile);
  fputs(&(cptr[ulii + 1]), outfile);
  putc(delimiter, outfile);
  fputs(paternal_ids? (&(paternal_ids[indiv_uidx * max_paternal_id_len])) : "0", outfile);
  putc(delimiter, outfile);
  fputs(maternal_ids? (&(maternal_ids[indiv_uidx * max_maternal_id_len])) : "0", outfile);
  putc(delimiter, outfile);
  putc(sexchar(sex_nm, sex_male, indiv_uidx), outfile);
  putc(delimiter, outfile);
  if (!IS_SET(pheno_nm, indiv_uidx)) {
    fputs(output_missing_pheno, outfile);
  } else if (pheno_c) {
    putc('1' + IS_SET(pheno_c, indiv_uidx), outfile);
  } else {
    cptr = double_g_write(wbuf, pheno_d[indiv_uidx]);
    fwrite(wbuf, 1, cptr - wbuf, outfile);
  }
  if (putc_checked(delimiter, outfile)) {
    return -1;
  }
  return 0;
}

void init_recode_cmax0(char* allele1, char* allele2, char** cur_mk_allelesx, uint32_t* cmalen, char delimiter, char delim2) {
  uintptr_t alen1 = strlen(allele1);
  uintptr_t alen2 = strlen(allele2);
  char* cmaptr = cur_mk_allelesx[0];
  *cmaptr++ = delimiter;
  cmaptr = memcpyax(cmaptr, allele1, alen1, delim2);
  memcpy(cmaptr, allele1, alen1);
  cmalen[0] = alen1 * 2 + 2;
  cmaptr = cur_mk_allelesx[2];
  *cmaptr++ = delimiter;
  cmaptr = memcpyax(cmaptr, allele1, alen1, delim2);
  memcpy(cmaptr, allele2, alen2);
  cmalen[2] = alen1 + alen2 + 2;
  cmaptr = cur_mk_allelesx[3];
  *cmaptr++ = delimiter;
  cmaptr = memcpyax(cmaptr, allele2, alen2, delim2);
  memcpy(cmaptr, allele2, alen2);
  cmalen[3] = alen2 * 2 + 2;
}

void init_recode_cmax(char* allele1, char* allele2, char** cur_mk_allelesx, uint32_t* cmalen, char delimiter, char delim2) {
  uintptr_t alen1 = strlen(allele1);
  uintptr_t alen2 = strlen(allele2);
  char* cmaptr = memcpyax(cur_mk_allelesx[0], allele1, alen1, delim2);
  memcpyx(cmaptr, allele1, alen1, delimiter);
  cmalen[0] = alen1 * 2 + 2;
  cmaptr = memcpyax(cur_mk_allelesx[2], allele1, alen1, delim2);
  memcpyx(cmaptr, allele2, alen2, delimiter);
  cmalen[2] = alen1 + alen2 + 2;
  cmaptr = memcpyax(cur_mk_allelesx[3], allele2, alen2, delim2);
  memcpyx(cmaptr, allele2, alen2, delimiter);
  cmalen[3] = alen2 * 2 + 2;
}

void init_cur_mk_allelesx(char* mk_alleles, uintptr_t max_marker_allele_len, uint32_t do_reverse, char** cur_mk_allelesx, uint32_t* cmalen, char delimiter, char delim2) {
  uint32_t alen = strlen_se(mk_alleles);
  memcpy(cur_mk_allelesx[0], mk_alleles, alen);
  cmalen[0] = alen;
  if (delimiter) {
    cmalen[0] += 1;
    cur_mk_allelesx[0][alen] = delim2;
    memcpy(cur_mk_allelesx[1], mk_alleles, alen);
    cur_mk_allelesx[1][alen] = delimiter;
  } else {
    cur_mk_allelesx[0][alen] = '\0';
  }
  mk_alleles = &(mk_alleles[max_marker_allele_len]);
  alen = strlen_se(mk_alleles);
  cmalen[1] = alen;
  if (delimiter) {
    cmalen[1] += 1;
    memcpy(cur_mk_allelesx[2], mk_alleles, alen);
    cur_mk_allelesx[2][alen] = delim2;
    memcpy(cur_mk_allelesx[3], mk_alleles, alen);
    cur_mk_allelesx[3][alen] = delimiter;
    if (do_reverse) {
      cur_mk_allelesx[4] = cur_mk_allelesx[2];
      cur_mk_allelesx[5] = cur_mk_allelesx[1];
    } else {
      cur_mk_allelesx[4] = cur_mk_allelesx[0];
      cur_mk_allelesx[5] = cur_mk_allelesx[3];
    }
  } else {
    memcpy(cur_mk_allelesx[1], mk_alleles, alen + 1);
    if (do_reverse) {
      cur_mk_allelesx[4] = cur_mk_allelesx[1];
      cur_mk_allelesx[5] = cur_mk_allelesx[0];
    } else {
      cur_mk_allelesx[4] = cur_mk_allelesx[0];
      cur_mk_allelesx[5] = cur_mk_allelesx[1];
    }
  }
  if (do_reverse) {
    cmalen[2] = cmalen[1];
    cmalen[3] = cmalen[0];
  } else {
    cmalen[2] = cmalen[0];
    cmalen[3] = cmalen[1];
  }
}

int32_t recode_beagle_new_chrom(char* outname, char* outname_end2, uintptr_t* marker_exclude, Chrom_info* chrom_info_ptr, uintptr_t* marker_uidx_ptr, uint32_t* chrom_fo_idx_ptr, uint32_t* chrom_idx_ptr, uint32_t* chrom_end_ptr, FILE* bedfile, uintptr_t bed_offset, uintptr_t unfiltered_indiv_ct4, FILE** datfile_ptr, FILE** mapfile_ptr, char* dat_header, uintptr_t dat_header_len) {
  uintptr_t marker_uidx = *marker_uidx_ptr;
  uint32_t chrom_fo_idx = *chrom_fo_idx_ptr;
  uint32_t chrom_end;
  uint32_t chrom_idx;
  uint32_t is_x;
  uint32_t is_y;
  uint32_t is_haploid;
  char* wbufptr;
  refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_haploid);
  chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
  if ((chrom_idx > 22) && (chrom_idx < 27)) {
    do {
      marker_uidx = next_unset_ul_unsafe(marker_exclude, chrom_end);
      chrom_fo_idx++;
      refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_haploid);
      chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
    } while ((chrom_idx > 22) && (chrom_idx < 27));
    if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
      return RET_READ_FAIL;
    }
  }

  wbufptr = chrom_name_write(outname_end2, chrom_info_ptr, chrom_idx, 0);
  memcpy(wbufptr, ".dat", 5);
  if (fopen_checked(datfile_ptr, outname, "w")) {
    return RET_OPEN_FAIL;
  }
  memcpy(wbufptr, ".map", 5);
  if (fopen_checked(mapfile_ptr, outname, "w")) {
    return RET_OPEN_FAIL;
  }
  if (fwrite_checked(dat_header, dat_header_len, *datfile_ptr)) {
    return RET_WRITE_FAIL;
  }
  *marker_uidx_ptr = marker_uidx;
  *chrom_fo_idx_ptr = chrom_fo_idx;
  *chrom_idx_ptr = chrom_idx;
  *chrom_end_ptr = chrom_end;
  return 0;
}

int32_t open_and_write_fastphase_header(FILE** outfile_ptr, char* outname, uintptr_t* marker_exclude, uint32_t* marker_pos, uint32_t marker_uidx, uint32_t chrom_size, uint32_t indiv_ct) {
  char wbuf[16];
  char* wptr;
  uint32_t marker_idx;
  if (fopen_checked(outfile_ptr, outname, "w")) {
    return RET_OPEN_FAIL;
  }
  wptr = uint32_writex(wbuf, indiv_ct, '\n');
  if (fwrite_checked(wbuf, wptr - wbuf, *outfile_ptr)) {
    return RET_WRITE_FAIL;
  }
  wptr = uint32_write(wbuf, chrom_size);
  fwrite(wbuf, 1, wptr - wbuf, *outfile_ptr);
  fputs("\nP ", *outfile_ptr);
  for (marker_idx = 0; marker_idx < chrom_size; marker_uidx++, marker_idx++) {
    next_unset_unsafe_ck(marker_exclude, &marker_uidx);
    wptr = uint32_writex(wbuf, marker_pos[marker_uidx], ' ');
    fwrite(wbuf, 1, wptr - wbuf, *outfile_ptr);
  }
  if (putc_checked('\n', *outfile_ptr)) {
    return RET_WRITE_FAIL;
  }
  return 0;
}

uint32_t write_ped_lines(FILE* outfile, unsigned char* loadbuf, uintptr_t* marker_exclude, uintptr_t marker_uidx_start, uintptr_t marker_ct, uintptr_t unfiltered_indiv_ct4, uintptr_t* indiv_exclude, uintptr_t* indiv_uidx_ptr, uintptr_t indiv_idx_start, uintptr_t indiv_idx_end, uint32_t recode_compound, char** mk_allele_ptrs, char delimiter, char delim2, char* person_ids, uintptr_t max_person_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* sex_nm, uintptr_t* sex_male, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, char output_missing_geno, char* output_missing_pheno, char* writebuf) {
  uintptr_t indiv_uidx = *indiv_uidx_ptr;
  char missing4[4];
  unsigned char* bufptr;
  char* wbufptr;
  char* aptr;
  uintptr_t marker_uidx;
  uintptr_t marker_idx;
  uintptr_t indiv_idx;
  uint32_t shiftval;
  uint32_t alen;
  unsigned char ucc;
  missing4[0] = output_missing_geno;
  missing4[1] = delim2;
  missing4[2] = output_missing_geno;
  missing4[3] = delimiter;
  for (indiv_idx = indiv_idx_start; indiv_idx < indiv_idx_end; indiv_uidx++, indiv_idx++) {
    next_unset_ul_unsafe_ck(indiv_exclude, &indiv_uidx);
    if (recode_write_first_cols(outfile, indiv_uidx, delimiter, person_ids, max_person_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, sex_nm, sex_male, pheno_nm, pheno_c, pheno_d, output_missing_pheno)) {
      return 1;
    }
    bufptr = &(loadbuf[indiv_uidx / 4]);
    wbufptr = writebuf;
    shiftval = (indiv_uidx % 4) * 2;
    marker_uidx = marker_uidx_start;
    if (recode_compound) {
      for (marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
	next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
	ucc = ((*bufptr) >> shiftval) & 3;
	if (ucc) {
	  if (ucc == 3) {
	    *wbufptr++ = mk_allele_ptrs[2 * marker_uidx + 1][0];
	    *wbufptr = mk_allele_ptrs[2 * marker_uidx + 1][0];
	  } else if (ucc == 2) {
	    *wbufptr++ = mk_allele_ptrs[2 * marker_uidx][0];
	    *wbufptr = mk_allele_ptrs[2 * marker_uidx + 1][0];
	  } else {
	    *wbufptr++ = output_missing_geno;
	    *wbufptr = output_missing_geno;
	  }
	} else {
	  *wbufptr++ = mk_allele_ptrs[2 * marker_uidx][0];
	  *wbufptr = mk_allele_ptrs[2 * marker_uidx][0];
	}
	bufptr = &(bufptr[unfiltered_indiv_ct4]);
	wbufptr = &(wbufptr[2]);
      }
      if (fwrite_checked(writebuf, marker_ct * 3, outfile)) {
	return 1;
      }
    } else {
      for (marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
	next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
	ucc = ((*bufptr) >> shiftval) & 3;
	if (ucc) {
	  if (ucc == 3) {
	    aptr = mk_allele_ptrs[2 * marker_uidx + 1];
	    alen = strlen(aptr);
	    memcpy(wbufptr, aptr, alen);
	    wbufptr[alen] = delim2;
	    memcpy(&(wbufptr[alen + 1]), aptr, alen);
	    wbufptr[2 * alen + 1] = delimiter;
	    wbufptr = &(wbufptr[2 * alen + 2]);
	  } else if (ucc == 2) {
	    wbufptr = strcpyax(wbufptr, mk_allele_ptrs[2 * marker_uidx], delim2);
	    wbufptr = strcpyax(wbufptr, mk_allele_ptrs[2 * marker_uidx + 1], delimiter);
	  } else {
	    wbufptr = memcpya(wbufptr, missing4, 4);
	  }
	} else {
	  aptr = mk_allele_ptrs[2 * marker_uidx];
	  alen = strlen(aptr);
	  memcpy(wbufptr, aptr, alen);
	  wbufptr[alen] = delim2;
	  memcpy(&(wbufptr[alen + 1]), aptr, alen);
	  wbufptr[2 * alen + 1] = delimiter;
	  wbufptr = &(wbufptr[2 * alen + 2]);
	}
	bufptr = &(bufptr[unfiltered_indiv_ct4]);
      }
      wbufptr[-1] = '\n';
      if (fwrite_checked(writebuf, wbufptr - writebuf, outfile)) {
	return 1;
      }
    }
  }
  *indiv_uidx_ptr = indiv_uidx;
  return 0;
}

uint32_t write_haploview_map(FILE* outfile, uintptr_t* marker_exclude, uintptr_t marker_uidx_start, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t* marker_pos) {
  char wbuf[16];
  char* wptr;
  uintptr_t marker_idx;
  for (marker_idx = 0; marker_idx < marker_ct; marker_uidx_start++, marker_idx++) {
    next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx_start);
    fputs(&(marker_ids[marker_uidx_start * max_marker_id_len]), outfile);
    putc('\t', outfile);
    wptr = uint32_writex(wbuf, marker_pos[marker_uidx_start], '\n');
    fwrite(wbuf, 1, wptr - wbuf, outfile);
  }
  if (ferror(outfile)) {
    return 1;
  }
  return 0;
}

int32_t recode(uint32_t recode_modifier, FILE* bedfile, uintptr_t bed_offset, FILE* famfile, char* outname, char* outname_end, char* recode_allele_name, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_ct, char* marker_ids, uintptr_t max_marker_id_len, double* marker_cms, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, uint32_t* marker_pos, uintptr_t* marker_reverse, char* person_ids, uintptr_t max_person_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* sex_nm, uintptr_t* sex_male, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, char* output_missing_pheno, uint64_t misc_flags, uint32_t hh_exists, Chrom_info* chrom_info_ptr) {
  FILE* outfile = NULL;
  FILE* outfile2 = NULL;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t indiv_ctv2 = 2 * ((indiv_ct + (BITCT - 1)) / BITCT);
  uintptr_t topsize = 0;
  unsigned char* wkspace_mark = wkspace_base;
  char delimiter = (recode_modifier & RECODE_TAB)? '\t' : ' ';
  uintptr_t* recode_allele_reverse = NULL;
  char** mk_allele_ptrs = marker_allele_ptrs;
  char** allele_missing = NULL;
  char* recode_allele_extra = NULL;
  char delim2 = delimiter;
  uintptr_t* indiv_include2 = NULL;
  uintptr_t* indiv_male_include2 = NULL;
  uint32_t lgen_ref = (recode_modifier & RECODE_LGEN_REF);
  uint32_t rlist = (recode_modifier & RECODE_RLIST);
  uint32_t vcf_not_fid = (recode_modifier & RECODE_VCF) && (!(recode_modifier & RECODE_FID));
  uint32_t vcf_not_iid = (recode_modifier & RECODE_VCF) && (!(recode_modifier & RECODE_IID));
  uint32_t vcf_two_ids = vcf_not_fid && vcf_not_iid;
  uint32_t recode_12 = recode_modifier & RECODE_12;
  uint32_t set_hh_missing = (misc_flags / MISC_SET_HH_MISSING) & 1;
  uint32_t zero_extra_chroms = (misc_flags / MISC_ZERO_EXTRA_CHROMS) & 1;
  uint32_t xmhh_exists_orig = hh_exists & XMHH_EXISTS;
  uintptr_t header_len = 0;
  uintptr_t max_chrom_size = 0;
  uintptr_t max_fid_len = 0;
  uintptr_t fid_ct = 0;
  uint32_t last_chrom_fo_idx = 0;
  uint32_t onechar_max = (chrom_info_ptr->max_code > 9)? 9 : chrom_info_ptr->max_code;
  uint32_t last_pos = 0;
  char missing_geno = *g_missing_geno_ptr;
  char output_missing_geno = *g_output_missing_geno_ptr;
  uintptr_t* loadbuf_collapsed = NULL;
  uintptr_t* loadbuf_collapsed_end = NULL;
  char* person_ids_collapsed = NULL;
  char* writebuf2 = NULL;
  char* writebuf3 = NULL;
  uint32_t* fid_map = NULL;
  char* cur_mk_allelesx_buf = NULL;
  int32_t retval = 0;
  char* writebufl[4];
  char* writebuflp[4];
  char* writebuflps[4];
  char* cur_mk_allelesx[6];
  uint32_t cmalen[4];
  uint32_t recode_compound;
  time_t rawtime;
  struct tm *loctime;
  uint32_t is_x;
  uint32_t is_y;
  uint32_t is_haploid;
  uint32_t chrom_fo_idx;
  uint32_t chrom_idx;
  uint32_t chrom_end;
  uintptr_t autosomal_marker_ct;
  uintptr_t marker_uidx_start;
  uintptr_t marker_uidx;
  uintptr_t marker_idx;
  uintptr_t indiv_uidx;
  uintptr_t indiv_idx;
  unsigned char* loadbuf;
  uintptr_t* ulptr;
  uintptr_t* ulptr2;
  uintptr_t* ulptr_end;
  char* writebuf;
  unsigned char* bufptr;
  char* wbufptr;
  char* cptr;
  char* aptr;
  char* aptr2;
  uint32_t alen;
  uint32_t alen2;
  unsigned char ucc;
  unsigned char ucc2;
  char cc;
  uint32_t pct;
  uint32_t uidx_cur;
  uint32_t uidx_stop;
  uintptr_t loop_end;
  uintptr_t cur_word;
  uintptr_t ref_word;
  uintptr_t ulii;
  uintptr_t uljj;
  uint32_t shiftval;
  uint32_t shiftmax;
  uint32_t cur_fid;
  uint32_t uii;
  int32_t ii;
  if (!hh_exists) {
    set_hh_missing = 0;
  }
  if (set_hh_missing) {
    if (recode_modifier & (RECODE_23 | RECODE_BEAGLE | RECODE_BIMBAM | RECODE_BIMBAM_1CHR | RECODE_LGEN | RECODE_LGEN_REF | RECODE_LIST | RECODE_RLIST | RECODE_TRANSPOSE | RECODE_VCF)) {
      // SNP-major and no need for indiv_uidx in inner loop, so we can use
      // collapsed representation
      if (alloc_collapsed_haploid_filters(unfiltered_indiv_ct, indiv_ct, hh_exists | ((recode_modifier & RECODE_VCF)? XMHH_EXISTS : 0), 0, indiv_exclude, sex_male, &indiv_include2, &indiv_male_include2)) {
	goto recode_ret_NOMEM;
      }
    } else {
      // individual-major output (in which case we load large blocks and use
      // haploid_fix_multiple())
      if (alloc_raw_haploid_filters(unfiltered_indiv_ct, hh_exists, 0, indiv_exclude, sex_male, &indiv_include2, &indiv_male_include2)) {
	goto recode_ret_NOMEM;
      }
    }
  }
  if (recode_12) {
    mk_allele_ptrs = (char**)wkspace_alloc(unfiltered_marker_ct * 2 * sizeof(intptr_t));
    if (!mk_allele_ptrs) {
      goto recode_ret_NOMEM;
    }
    for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
      mk_allele_ptrs[2 * marker_uidx] = (char*)(&(g_one_char_strs[98])); // "1"
      mk_allele_ptrs[2 * marker_uidx + 1] = (char*)(&(g_one_char_strs[100])); // "2"
    }
    max_marker_allele_len = 2;
  }
  if (recode_modifier & (RECODE_BEAGLE | RECODE_BIMBAM | RECODE_BIMBAM_1CHR | RECODE_LGEN | RECODE_LGEN_REF | RECODE_LIST | RECODE_RLIST | RECODE_TRANSPOSE | RECODE_VCF)) {
    if (wkspace_alloc_ul_checked(&loadbuf_collapsed, indiv_ctv2 * sizeof(intptr_t))) {
      goto recode_ret_NOMEM;
    }
    loadbuf_collapsed_end = &(loadbuf_collapsed[(indiv_ct + (BITCT2 - 1)) / BITCT2]);
    if (recode_modifier & (RECODE_LGEN | RECODE_LGEN_REF | RECODE_LIST | RECODE_RLIST)) {
      // need to collapse person_ids to remove need for indiv_uidx in inner
      // loop
      person_ids_collapsed = alloc_and_init_collapsed_arr(person_ids, max_person_id_len, unfiltered_indiv_ct, indiv_exclude, indiv_ct, 1);
      if (!person_ids_collapsed) {
        goto recode_ret_NOMEM;
      }
    }
  }
  if (recode_modifier & RECODE_VCF) {
    if (wkspace_alloc_c_checked(&writebuf, indiv_ct * 4)) {
      goto recode_ret_NOMEM;
    }
  } else if (recode_modifier & RECODE_BEAGLE) {
    // common header:
    // "P FID " + ... + "\n"
    // "I IID " + ... + "\n"
    // "? PHE " + ... + "\n"
    // ... space requirement is bounded above by
    //   (2 * indiv_ct * max_person_id_len) + 64 * indiv_ct
    //
    // per-marker:
    //   "M " + [marker name] + " " + ... + "\n"
    ulii = strlen(output_missing_pheno);
    if (wkspace_alloc_c_checked(&writebuf, 2 * ulii + 2) ||
        wkspace_alloc_c_checked(&writebuf2, 21 + indiv_ct * (2 * max_person_id_len + 64))) {
      goto recode_ret_NOMEM;
    }
    wbufptr = memcpyax(writebuf, output_missing_pheno, ulii, ' ');
    wbufptr = memcpyax(wbufptr, output_missing_pheno, ulii, ' ');
    ulii = 2 * ulii + 2;
    wbufptr = memcpya(writebuf2, "P FID ", 6);
    indiv_uidx = 0;
    for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_uidx++, indiv_idx++) {
      next_unset_ul_unsafe_ck(indiv_exclude, &indiv_uidx);
      cptr = &(person_ids[indiv_uidx * max_person_id_len]);
      aptr = (char*)memchr(cptr, '\t', max_person_id_len);
      ulii = (uintptr_t)(aptr - cptr);
      wbufptr = memcpyax(wbufptr, cptr, ulii, ' ');
      wbufptr = memcpyax(wbufptr, cptr, ulii, ' ');
    }
    wbufptr = memcpya(wbufptr, "\nI IID ", 7);
    indiv_uidx = 0;
    for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_uidx++, indiv_idx++) {
      next_unset_ul_unsafe_ck(indiv_exclude, &indiv_uidx);
      cptr = &(person_ids[indiv_uidx * max_person_id_len]);
      cptr = (char*)memchr(cptr, '\t', max_person_id_len);
      cptr++;
      ulii = strlen(cptr);
      wbufptr = memcpyax(wbufptr, cptr, ulii, ' ');
      wbufptr = memcpyax(wbufptr, cptr, ulii, ' ');
    }
    indiv_uidx = 0;
    if (pheno_c) {
      wbufptr = memcpya(wbufptr, "\nA PHE ", 7);
      for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_uidx++, indiv_idx++) {
	next_unset_ul_unsafe_ck(indiv_exclude, &indiv_uidx);
        if (IS_SET(pheno_nm, indiv_uidx)) {
	  cc = (unsigned char)('1' + IS_SET(pheno_c, indiv_uidx));
	  *wbufptr++ = cc;
	  *wbufptr++ = ' ';
	  *wbufptr++ = cc;
	  *wbufptr++ = ' ';
	} else {
	  wbufptr = memcpya(wbufptr, writebuf, ulii);
	}
      }
    } else {
      wbufptr = memcpya(wbufptr, "\nT PHE ", 7);
      for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_uidx++, indiv_idx++) {
	next_unset_ul_unsafe_ck(indiv_exclude, &indiv_uidx);
        if (IS_SET(pheno_nm, indiv_uidx)) {
	  cptr = double_g_writex(wbufptr, pheno_d[indiv_uidx], ' ');
          wbufptr = memcpya(cptr, wbufptr, (uintptr_t)(cptr - wbufptr));
	} else {
	  wbufptr = memcpya(wbufptr, writebuf, ulii);
	}
      }
    }
    *wbufptr++ = '\n';
    // free unused space, and save header length
    header_len = (uintptr_t)(wbufptr - writebuf2);
    wkspace_reset((unsigned char*)writebuf2);
    writebuf2 = (char*)wkspace_alloc(header_len);
    cmalen[1] = 4;
    ulii = 2 * max_marker_allele_len;
    if (wkspace_alloc_c_checked(&cur_mk_allelesx_buf, 4 * max_marker_allele_len)) {
      goto recode_ret_NOMEM;
    }
    cur_mk_allelesx[0] = cur_mk_allelesx_buf;
    cur_mk_allelesx[1] = &(cur_mk_allelesx_buf[ulii]);
    cur_mk_allelesx[2] = &(cur_mk_allelesx_buf[ulii * 2]);
    cur_mk_allelesx[3] = &(cur_mk_allelesx_buf[ulii * 3]);
    memcpy(cur_mk_allelesx[0], "    ", 4);
    memcpy(cur_mk_allelesx[1], " 0 0", 4);
    memcpy(cur_mk_allelesx[2], "    ", 4);
    memcpy(cur_mk_allelesx[3], "    ", 4);
  } else if (recode_modifier & (RECODE_BIMBAM | RECODE_BIMBAM_1CHR)) {
    if (max_marker_allele_len != 2) {
      logprint("Error: --recode bimbam cannot be used with multi-character allele names.\n");
      goto recode_ret_INVALID_FORMAT;
    }
    indiv_uidx = 0;
    for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_uidx++, indiv_idx++) {
      next_unset_ul_unsafe_ck(indiv_exclude, &indiv_uidx);
      cptr = (char*)memchr(&(person_ids[indiv_uidx * max_person_id_len]), '\t', max_person_id_len);
      if (strchr(&(cptr[1]), ',')) {
        logprint("Error: Comma present in individual ID during --recode bimbam run.\n");
        goto recode_ret_INVALID_FORMAT;
      }
    }
    marker_uidx = 0;
    for (marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
      next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
      if (strchr(&(marker_ids[marker_uidx * max_marker_id_len]), ',')) {
        logprint("Error: Comma present in SNP ID during --recode bimbam run.\n");
      }
    }
    // +1 because memcpyl3a() copies an extra character
    if (wkspace_alloc_c_checked(&writebuf, 3 * indiv_ct + 1) ||
        wkspace_alloc_c_checked(&writebuf2, 32)) {
      goto recode_ret_NOMEM;
    }
  } else if (recode_modifier & (RECODE_FASTPHASE | RECODE_FASTPHASE_1CHR)) {
    max_chrom_size = get_max_chrom_size(chrom_info_ptr, marker_exclude, &last_chrom_fo_idx);
    if ((recode_modifier & RECODE_FASTPHASE_1CHR) && (max_chrom_size != marker_ct)) {
      logprint("Error: --recode fastphase-1chr requires a single-chromosome dataset.  Did you\nmean '--recode fastphase'?  (Note the lack of a dash in the middle.)\n");
      goto recode_ret_INVALID_CMDLINE;
    }
    if (wkspace_left < ((uint64_t)unfiltered_indiv_ct4) * max_chrom_size + 2 * ((max_chrom_size + 63) & (~(63 * ONELU)))) {
      goto recode_ret_NO_MULTIPASS_YET;
    }
    if (wkspace_alloc_c_checked(&writebuf, max_chrom_size) ||
        wkspace_alloc_c_checked(&writebuf2, max_chrom_size)) {
      goto recode_ret_NOMEM;
    }
  } else if (recode_modifier & (RECODE_LGEN | RECODE_LGEN_REF)) {
    if (recode_modifier & RECODE_COMPOUND) {
      if (max_marker_allele_len != 2) {
        logprint("Error: --recode compound-genotypes cannot be used with multi-character allele\nnames.\n");
        goto recode_ret_INVALID_FORMAT;
      }
      ulii = 0;
    } else {
      ulii = 1;
    }
    ulii += 2 * max_marker_allele_len + max_marker_id_len + max_person_id_len;
    if (wkspace_alloc_c_checked(&writebuf, 4 * ulii)) {
      goto recode_ret_NOMEM;
    }
  } else if (recode_modifier & (RECODE_LIST | RECODE_RLIST)) {
    // --list:
    // 3 for chromosome and delim
    // + max_marker_id_len
    // + 3, or (2 * max_marker_allele_len - 1)
    // + indiv_ct * max_person_id_len + 1
    //
    // --rlist:
    // max_marker_id_len
    // + 4 for "HOM"/"HET"/"NIL" and delim
    // + 4, or (2 * max_marker_allele_len)
    // + indiv_ct * max_person_id_len + 1
    ulii = 3 + max_marker_id_len + 2 * max_marker_allele_len + indiv_ct * max_person_id_len;
    if ((recode_modifier & RECODE_LIST) && (max_marker_allele_len != 2)) {
      logprint("Error: --recode list cannot be used with multi-character allele names.\n");
      goto recode_ret_INVALID_FORMAT;
    }
    if (rlist) {
      ulii += 2;
    }
    if (wkspace_alloc_c_checked(&writebuf, ulii * 4)) {
      goto recode_ret_NOMEM;
    }
    writebufl[0] = writebuf;
    writebufl[1] = &(writebuf[ulii]);
    writebufl[2] = &(writebuf[ulii * 2]);
    writebufl[3] = &(writebuf[ulii * 3]);
  } else if (recode_modifier & RECODE_23) {
    if (indiv_ct != 1) {
      logprint("Error: --recode 23 can only be used on a file with exactly one individual.\n");
      goto recode_ret_INVALID_FORMAT;
    } else if (max_marker_allele_len != 2) {
      logprint("Error: --recode 23 cannot be used with multi-character allele names.\n");
      goto recode_ret_INVALID_FORMAT;
    }
    // chromosome code, marker position, single-char alleles
    if (wkspace_alloc_c_checked(&writebuf, 32)) {
      goto recode_ret_NOMEM;
    }
  } else if (recode_modifier & RECODE_STRUCTURE) {
    indiv_uidx = 0;
    for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_uidx++, indiv_idx++) {
      next_unset_ul_unsafe_ck(indiv_exclude, &indiv_uidx);
      cptr = &(person_ids[indiv_uidx * max_person_id_len]);
      aptr = (char*)memchr(cptr, '\t', max_person_id_len);
      ulii = (uintptr_t)(aptr - cptr);
      if (ulii >= max_fid_len) {
        max_fid_len = ulii + 1;
	if (max_fid_len == max_person_id_len - 2) {
	  break;
	}
      }
    }
    if (wkspace_alloc_c_checked(&writebuf3, max_fid_len * indiv_ct)) {
      goto recode_ret_NOMEM;
    }
    indiv_uidx = 0;
    for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_uidx++, indiv_idx++) {
      next_unset_ul_unsafe_ck(indiv_exclude, &indiv_uidx);
      cptr = &(person_ids[indiv_uidx * max_person_id_len]);
      aptr = (char*)memchr(cptr, '\t', max_fid_len);
      ulii = (uintptr_t)(aptr - cptr);
      memcpy(&(writebuf3[indiv_idx * max_fid_len]), cptr, ulii);
      writebuf3[indiv_idx * max_fid_len + ulii] = '\0';
    }
    qsort(writebuf3, indiv_ct, max_fid_len, strcmp_casted);
    // simpler collapse_duplicate_ids(), probably want to make this a function
    for (ulii = 1; ulii < indiv_ct; ulii++) {
      if (!strcmp(&(writebuf3[(ulii - 1) * max_fid_len]), &(writebuf3[ulii * max_fid_len]))) {
	break;
      }
    }
    fid_ct = ulii;
    while (++ulii < indiv_ct) {
      if (strcmp(&(writebuf3[(fid_ct - 1) * max_fid_len]), &(writebuf3[ulii * max_fid_len]))) {
        strcpy(&(writebuf3[fid_ct * max_fid_len]), &(writebuf3[ulii * max_fid_len]));
	fid_ct++;
      }
    }
    wkspace_reset((unsigned char*)writebuf3);
    writebuf3 = (char*)wkspace_alloc(fid_ct * max_fid_len);
    if (wkspace_alloc_ui_checked(&fid_map, fid_ct * sizeof(int32_t)) ||
        wkspace_alloc_c_checked(&writebuf, 4 * marker_ct) ||
        wkspace_alloc_c_checked(&writebuf2, 16)) {
      goto recode_ret_NOMEM;
    }
    fill_uint_zero(fid_map, fid_ct);
  } else {
    if (recode_modifier & (RECODE_A | RECODE_AD)) {
      recode_modifier &= ~RECODE_COMPOUND;
      if (wkspace_alloc_c_checked(&writebuf2, 32)) {
	goto recode_ret_NOMEM;
      }
      memcpy(writebuf2, "2 0     1 1 0 0 NA NA", 21);
      writebuf2[1] = delimiter;
      writebuf2[3] = delimiter;
      writebuf2[9] = delimiter;
      writebuf2[11] = delimiter;
      writebuf2[13] = delimiter;
      writebuf2[15] = delimiter;
      writebuf2[18] = delimiter;
      writebuf2[21] = delimiter;
    }
    max_chrom_size = marker_ct;
    if (recode_modifier & RECODE_AD) {
      ulii = 6;
    } else if (recode_modifier & (RECODE_A | RECODE_COMPOUND)) {
      if ((max_marker_allele_len != 2) && (recode_modifier & RECODE_COMPOUND)) {
	logprint("Error: --recode compound-genotypes cannot be used with multi-character allele\nnames.\n");
	goto recode_ret_INVALID_FORMAT;
      }
      ulii = 3;
    } else if (recode_modifier & (RECODE_HV | RECODE_HV_1CHR)) {
      max_chrom_size = get_max_chrom_size(chrom_info_ptr, marker_exclude, &last_chrom_fo_idx);
      if ((recode_modifier & RECODE_HV_1CHR) && (max_chrom_size != marker_ct)) {
	logprint("Error: --recode HV-1chr requires a single-chromosome dataset.  Did you mean\n'--recode HV'?  (Note the lack of a dash in the middle.)\n");
	goto recode_ret_INVALID_CMDLINE;
      }
      if (max_marker_allele_len == 2) {
        ulii = max_chrom_size * 4;
      } else {
        ulii = 0;
	for (chrom_fo_idx = 0; chrom_fo_idx < chrom_info_ptr->chrom_ct; chrom_fo_idx++) {
          uljj = 0;
	  chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx];
	  for (marker_uidx = next_unset_ul(marker_exclude, chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx], chrom_end); marker_uidx < chrom_end;) {
	    alen = strlen(mk_allele_ptrs[marker_uidx * 2]);
	    alen2 = strlen(mk_allele_ptrs[marker_uidx * 2 + 1]);
	    uljj += MAXV(alen, alen2) + 1;
	    marker_uidx++;
            next_unset_ul_ck(marker_exclude, &marker_uidx, chrom_end);
	  }
          if (uljj > ulii) {
	    ulii = uljj;
	  }
	}
	ulii *= 2;
      }
    } else {
      // all chromosomes at once
      // calculate maximum length of .ped line
      if (max_marker_allele_len == 2) {
	ulii = marker_ct * 4;
      } else {
	ulii = marker_ct;
	for (marker_uidx = 0, marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
	  next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
	  alen = strlen(mk_allele_ptrs[marker_uidx * 2]);
	  alen2 = strlen(mk_allele_ptrs[marker_uidx * 2 + 1]);
	  ulii += MAXV(alen, alen2);
	}
	ulii *= 2;
      }
    }
    if ((recode_modifier & (RECODE_TAB | RECODE_DELIMX)) == RECODE_TAB) {
      delim2 = ' ';
    }
    if (recode_modifier & (RECODE_A | RECODE_AD | RECODE_COMPOUND)) {
      if (wkspace_alloc_c_checked(&writebuf, max_chrom_size * ulii)) {
	goto recode_ret_NOMEM;
      }
      if (recode_modifier & RECODE_COMPOUND) {
        memset(writebuf, delimiter, max_chrom_size * 3 - 1);
	writebuf[max_chrom_size * 3 - 1] = '\n';
      }
      if (recode_allele_name) {
	ulii = ((unfiltered_marker_ct + (BITCT - 1)) / BITCT) * sizeof(intptr_t);
	if (wkspace_alloc_ul_checked(&recode_allele_reverse, ulii * sizeof(intptr_t))) {
	  goto recode_ret_NOMEM;
	}
	// this indicates when we want to report the A2 allele instead of the
	// A1.  (potential double negatives, bleah)
	fill_ulong_zero(recode_allele_reverse, ulii);
	allele_missing = (char**)wkspace_alloc(unfiltered_marker_ct * sizeof(char**));
	if (!allele_missing) {
	  goto recode_ret_NOMEM;
	}
	recode_allele_extra = (char*)wkspace_base;
	fill_ulong_zero((uintptr_t*)allele_missing, unfiltered_marker_ct);
	ulii = (max_marker_allele_len + MAXLINELEN + 15) & (~(15 * ONELU));
	loadbuf = (unsigned char*)top_alloc(&topsize, ulii);
	if (!loadbuf) {
	  goto recode_ret_NOMEM;
	}
	wkspace_left -= topsize;
	// When '12' and 'A'/'AD' are simultaneously present, most sensible
	// behavior is to match against real allele IDs and just apply '12' to
	// the output header line.  If that's not what the user wants, they can
	// do a two-step recode.
	// (--recode12 simply overrode --recodeA/--recodeAD in PLINK 1.07; no
	// need to replicate that.) 
	retval = recode_allele_load((char*)loadbuf, ulii, recode_allele_name, &allele_missing, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, marker_allele_ptrs, max_marker_allele_len, recode_allele_reverse, recode_allele_extra);
	wkspace_left += topsize;
	topsize = 0;
	if (retval) {
	  goto recode_ret_1;
	}
      }
    } else {
      // --recode, --recode HV
      if (wkspace_alloc_c_checked(&writebuf, ulii)) {
	goto recode_ret_NOMEM;
      }
    }
  }
  recode_compound = recode_modifier & RECODE_COMPOUND;
  if (!(recode_modifier & (RECODE_A | RECODE_AD | RECODE_BEAGLE | RECODE_FASTPHASE | RECODE_FASTPHASE_1CHR | RECODE_LGEN | RECODE_LGEN_REF | RECODE_VCF))) {
    if (wkspace_alloc_c_checked(&cur_mk_allelesx_buf, 8 * max_marker_allele_len)) {
      goto recode_ret_NOMEM;
    }
    cur_mk_allelesx[0] = cur_mk_allelesx_buf;
    cur_mk_allelesx[1] = &(cur_mk_allelesx_buf[max_marker_allele_len * 2]);
    cur_mk_allelesx[2] = &(cur_mk_allelesx_buf[max_marker_allele_len * 4]);
    cur_mk_allelesx[3] = &(cur_mk_allelesx_buf[max_marker_allele_len * 6]);
  } else if (recode_modifier & RECODE_VCF) {
    if (wkspace_alloc_c_checked(&cur_mk_allelesx_buf, 16)) {
      goto recode_ret_NOMEM;
    }
    cur_mk_allelesx[0] = cur_mk_allelesx_buf;
    cur_mk_allelesx[1] = &(cur_mk_allelesx_buf[4]);
    cur_mk_allelesx[2] = &(cur_mk_allelesx_buf[8]);
    cur_mk_allelesx[3] = &(cur_mk_allelesx_buf[12]);
  }
  if (fseeko(bedfile, bed_offset, SEEK_SET)) {
    goto recode_ret_READ_FAIL;
  }
  marker_uidx = 0;
  marker_idx = 0;
  if (wkspace_left < unfiltered_indiv_ct4) {
    goto recode_ret_NOMEM;
  }
  loadbuf = wkspace_base;
  chrom_fo_idx = 0;
  refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_haploid);
  chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
  if (recode_modifier & RECODE_TRANSPOSE) {
    strcpy(outname_end, ".tped");
    if (fopen_checked(&outfile, outname, "w")) {
      goto recode_ret_OPEN_FAIL;
    }
    sprintf(logbuf, "--recode transpose to %s + .tfam... ", outname);
    logprintb();
    fputs("0%", stdout);
    fflush(stdout);
    if ((recode_modifier & (RECODE_TAB | RECODE_DELIMX)) == RECODE_TAB) {
      delim2 = ' '; // within genotype
    }
    cur_mk_allelesx[1][0] = delimiter;
    cur_mk_allelesx[1][1] = output_missing_geno;
    cur_mk_allelesx[1][2] = delim2;
    cur_mk_allelesx[1][3] = output_missing_geno;
    cmalen[1] = 4;

    cptr = chrom_name_write(tbuf, chrom_info_ptr, chrom_idx, zero_extra_chroms);
    *cptr++ = delimiter;
    for (pct = 1; pct <= 100; pct++) {
      loop_end = (((uint64_t)pct) * marker_ct) / 100;
      for (; marker_idx < loop_end; marker_uidx++, marker_idx++) {
	if (IS_SET(marker_exclude, marker_uidx)) {
	  marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
	  if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	    goto recode_ret_READ_FAIL;
	  }
	}
	if (marker_uidx >= chrom_end) {
	  chrom_fo_idx++;
	  refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_haploid);
	  chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	  cptr = chrom_name_write(tbuf, chrom_info_ptr, chrom_idx, zero_extra_chroms);
	  *cptr++ = delimiter;
	}
	wbufptr = strcpyax(cptr, &(marker_ids[marker_uidx * max_marker_id_len]), delimiter);
	if (!marker_cms) {
	  *wbufptr++ = '0';
	} else {
	  wbufptr = double_g_write(wbufptr, marker_cms[marker_uidx]);
	}
	*wbufptr++ = delimiter;
	wbufptr = uint32_write(wbufptr, marker_pos[marker_uidx]);
        if (fwrite_checked(tbuf, wbufptr - tbuf, outfile)) {
	  goto recode_ret_WRITE_FAIL;
	}

	if (load_and_collapse(bedfile, (uintptr_t*)loadbuf, unfiltered_indiv_ct, loadbuf_collapsed, indiv_ct, indiv_exclude, IS_SET(marker_reverse, marker_uidx))) {
	  goto recode_ret_READ_FAIL;
	}
	if (is_haploid && set_hh_missing) {
          haploid_fix(hh_exists, indiv_include2, indiv_male_include2, indiv_ct, is_x, is_y, (unsigned char*)loadbuf_collapsed);
	}
	init_recode_cmax0(mk_allele_ptrs[2 * marker_uidx], mk_allele_ptrs[2 * marker_uidx + 1], cur_mk_allelesx, cmalen, delimiter, delim2);
	ulptr = loadbuf_collapsed;
	ulptr_end = &(loadbuf_collapsed[indiv_ct / BITCT2]);
	shiftmax = BITCT2;
	while (1) {
	  while (ulptr < ulptr_end) {
	    cur_word = *ulptr++;
	    for (shiftval = 0; shiftval < shiftmax; shiftval++) {
	      ulii = cur_word & 3;
	      fwrite(cur_mk_allelesx[ulii], 1, cmalen[ulii], outfile);
	      cur_word >>= 2;
	    }
	  }
	  if (ulptr == loadbuf_collapsed_end) {
	    break;
	  }
	  ulptr_end++;
	  shiftmax = indiv_ct % BITCT2;
	}
	if (putc_checked('\n', outfile)) {
	  goto recode_ret_WRITE_FAIL;
	}
      }
      if (pct < 100) {
	if (pct > 10) {
	  putchar('\b');
	}
	printf("\b\b%u%%", pct);
	fflush(stdout);
      }
    }
  } else if (recode_modifier & RECODE_VCF) {
    strcpy(outname_end, ".vcf");
    if (fopen_checked(&outfile, outname, "w")) {
      goto recode_ret_OPEN_FAIL;
    }
    if (fputs_checked(
"##fileformat=VCFv4.1\n"
"##filedate=", outfile)) {
      goto recode_ret_WRITE_FAIL;
    }
    time(&rawtime);
    loctime = localtime(&rawtime);
    strftime(tbuf, MAXLINELEN, "%Y%m%d", loctime);
    fputs(tbuf, outfile);
    fputs(
"\n##source=" PROG_NAME_CAPS "\n"
"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
, outfile);
    indiv_uidx = 0;
    shiftval = 0; // repurposed: underscore seen in ID?
    for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_uidx++, indiv_idx++) {
      next_unset_ul_unsafe_ck(indiv_exclude, &indiv_uidx);
      cptr = &(person_ids[indiv_uidx * max_person_id_len]);
      ulii = strlen_se(cptr);
      putc('\t', outfile);
      if (vcf_not_iid) {
	fwrite(cptr, 1, ulii, outfile);
	if (vcf_two_ids) {
	  if (!shiftval) {
	    if (strchr(cptr, '_')) {
	      shiftval = 1;
	      logprint("Warning: Underscore(s) present in individual IDs.\n");
	    }
	  }
	  putc('_', outfile);
	}
      }
      if (vcf_not_fid) {
	fputs(&(cptr[ulii + 1]), outfile);
      }
    }
    sprintf(logbuf, "--recode vcf%s to %s... ", vcf_not_iid? (vcf_not_fid? "" : "-fid") : "-iid", outname);
    logprintb();
    fputs("0%", stdout);
    tbuf[0] = '\n';
    if (set_hh_missing && is_haploid && (!is_x)) {
      uii = 2;
    } else {
      uii = 4;
    }
    cmalen[0] = uii;
    cmalen[1] = uii;
    cmalen[2] = uii;
    cmalen[3] = uii;
    memcpy(cur_mk_allelesx[0], "\t1/1", 4);
    memcpy(cur_mk_allelesx[1], "\t./.", 4);
    memcpy(cur_mk_allelesx[2], "\t0/1", 4);
    memcpy(cur_mk_allelesx[3], "\t0/0", 4);
    for (pct = 1; pct <= 100; pct++) {
      loop_end = (((uint64_t)pct) * marker_ct) / 100;
      for (; marker_idx < loop_end; marker_uidx++, marker_idx++) {
	if (IS_SET(marker_exclude, marker_uidx)) {
          marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
	  if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	    goto recode_ret_READ_FAIL;
	  }
	}
	if (marker_uidx >= chrom_end) {
	  chrom_fo_idx++;
	  refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_haploid);
	  chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	  // todo: only emit one for male X?
	  if (set_hh_missing && is_haploid && (!is_x)) {
	    uii = 2;
	  } else {
	    uii = 4;
	  }
	  cmalen[0] = uii;
	  cmalen[1] = uii;
	  cmalen[2] = uii;
	  cmalen[3] = uii;
	}
	wbufptr = chrom_name_write(&(tbuf[1]), chrom_info_ptr, chrom_idx, zero_extra_chroms);
	*wbufptr++ = '\t';
	wbufptr = uint32_writex(wbufptr, marker_pos[marker_uidx], '\t');
	wbufptr = strcpyax(wbufptr, &(marker_ids[marker_uidx * max_marker_id_len]), '\t');
	if (fwrite_checked(tbuf, wbufptr - tbuf, outfile)) {
	  goto recode_ret_WRITE_FAIL;
	}
	fputs(mk_allele_ptrs[2 * marker_uidx + 1], outfile);
	putc('\t', outfile);
	fputs(mk_allele_ptrs[2 * marker_uidx], outfile);
	fputs("\t.\t.\t.\tGT", outfile);

	if (load_and_collapse(bedfile, (uintptr_t*)loadbuf, unfiltered_indiv_ct, loadbuf_collapsed, indiv_ct, indiv_exclude, IS_SET(marker_reverse, marker_uidx))) {
	  goto recode_ret_READ_FAIL;
	}
	if (is_haploid && set_hh_missing) {
	  haploid_fix(hh_exists, indiv_include2, indiv_male_include2, indiv_ct, is_x, is_y, (unsigned char*)loadbuf_collapsed);
	}
	wbufptr = writebuf;
	ulptr = loadbuf_collapsed;
	ulptr_end = &(loadbuf_collapsed[indiv_ct / BITCT2]);
	shiftmax = BITCT2;
	if ((!set_hh_missing) || (!is_x)) {
	  while (1) {
	    while (ulptr < ulptr_end) {
	      cur_word = *ulptr++;
	      for (shiftval = 0; shiftval < shiftmax; shiftval++) {
		ulii = cur_word & 3;
		wbufptr = memcpya(wbufptr, cur_mk_allelesx[ulii], cmalen[ulii]);
		cur_word >>= 2;
	      }
	    }
	    if (ulptr == loadbuf_collapsed_end) {
	      break;
	    }
	    ulptr_end++;
	    shiftmax = indiv_ct % BITCT2;
	  }
	} else {
	  ulptr2 = indiv_male_include2;
	  while (1) {
	    while (ulptr < ulptr_end) {
	      cur_word = *ulptr++;
	      ref_word = (*ulptr2++) << 1;
	      for (shiftval = 0; shiftval < shiftmax; shiftval++) {
		ulii = cur_word & 3;
		uljj = ref_word & 3;
		wbufptr = memcpya(wbufptr, cur_mk_allelesx[ulii], 4 - uljj);
		cur_word >>= 2;
		ref_word >>= 2;
	      }
	    }
	    if (ulptr == loadbuf_collapsed_end) {
	      break;
	    }
	    ulptr_end++;
	    shiftmax = indiv_ct % BITCT2;
	  }
	}
	if (fwrite_checked(writebuf, wbufptr - writebuf, outfile)) {
	  goto recode_ret_WRITE_FAIL;
	}
      }
      if (pct < 100) {
	if (pct > 10) {
	  putchar('\b');
	}
	printf("\b\b%u%%", pct);
	fflush(stdout);
      }
    }
    if (putc_checked('\n', outfile)) {
      goto recode_ret_WRITE_FAIL;
    }
  } else if (recode_modifier & RECODE_23) {
    strcpy(outname_end, ".txt");
    if (fopen_checked(&outfile, outname, "w")) {
      goto recode_ret_OPEN_FAIL;
    }
    sprintf(logbuf, "--recode 23 to %s... ", outname);
    logprintb();
    time(&rawtime);
    fprintf(outfile, "# This data file generated by " PROG_NAME_CAPS " at: %s", ctime(&rawtime));
    if (fputs_checked(
"#\n"
"# Below is a text version of your data.  Fields are TAB-separated.\n"
"# Each line corresponds to a single SNP.  For each SNP, we provide its\n"
"# identifier, its location on a reference human genome, and the genotype call.\n"
"# For further information (e.g. which reference build was used), consult the\n"
"# original source of your data.\n"
"#\n"
"# rsid\tchromosome\tposition\tgenotype\n"
, outfile)) {
      goto recode_ret_WRITE_FAIL;
    }
    writebuf[0] = '\t';
    cptr = chrom_print_human(&(writebuf[1]), chrom_idx);
    *cptr++ = '\t';
    indiv_uidx = next_unset_unsafe(indiv_exclude, 0);
    ucc = IS_SET(sex_male, indiv_uidx);
    ucc2 = ((chrom_idx == 24) || (chrom_idx == 26) || (ucc && (chrom_idx == 23) && (!xmhh_exists_orig)))? 1 : 0;
    // todo: add a way to use output of --23file-make-xylist
    for (; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
      marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
      if (marker_uidx >= chrom_end) {
	chrom_fo_idx++;
	refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_haploid);
	chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
        cptr = chrom_print_human(&(writebuf[1]), chrom_idx);
	*cptr++ = '\t';
	ucc2 = ((chrom_idx == 24) || (chrom_idx == 26) || (ucc && (chrom_idx == 23) && (!xmhh_exists_orig)))? 1 : 0;
      }

      if (fseeko(bedfile, bed_offset + (indiv_uidx / 4) + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	goto recode_ret_READ_FAIL;
      }
      ii = fgetc(bedfile);
      if (ii == EOF) {
	goto recode_ret_READ_FAIL;
      }
      cur_word = (((uint32_t)ii) >> ((indiv_uidx % 4) * 2)) & 3;
      if (is_haploid && set_hh_missing) {
	haploid_fix(hh_exists, indiv_include2, indiv_male_include2, 1, is_x, is_y, (unsigned char*)(&cur_word));
      }
      if (IS_SET(marker_reverse, marker_uidx)) {
	cur_word = cur_word ^ (((~(cur_word ^ (cur_word >> 1))) & 1) * 3);
      }
      aptr = uint32_writex(cptr, marker_pos[marker_uidx], '\t');
      if (cur_word) {
	if (cur_word == 3) {
	  *aptr++ = mk_allele_ptrs[2 * marker_uidx + 1][0];
	  *aptr++ = mk_allele_ptrs[2 * marker_uidx + 1][0];
	} else if (cur_word == 2) {
	  *aptr++ = mk_allele_ptrs[2 * marker_uidx + 1][0];
	  *aptr++ = mk_allele_ptrs[2 * marker_uidx][0];
	} else {
	  aptr = memseta(aptr, '-', 2);
	}
      } else {
	*aptr++ = mk_allele_ptrs[2 * marker_uidx][0];
	*aptr++ = mk_allele_ptrs[2 * marker_uidx][0];
      }
      fputs(&(marker_ids[marker_uidx * max_marker_id_len]), outfile);
      if (ucc2 && (cur_word != 1)) {
	aptr--;
      }
      *aptr++ = '\n';
      if (fwrite_checked(writebuf, aptr - writebuf, outfile)) {
	goto recode_ret_WRITE_FAIL;
      }
    }
  } else if (recode_modifier & RECODE_BEAGLE) {
    // for backward compatibility, also exclude XY.  don't exclude custom name
    // chromosomes, though, since chromosome 0 was actually processed
    autosomal_marker_ct = marker_ct - count_non_autosomal_markers(chrom_info_ptr, marker_exclude, 1) - count_chrom_markers(chrom_info_ptr, 25, marker_exclude);
    if (!autosomal_marker_ct) {
      logprint("Error: No autosomal variants for --recode beagle.\n");
      goto recode_ret_ALL_MARKERS_EXCLUDED;
    }
    memcpy(outname_end, ".chr-", 6);
    sprintf(logbuf, "--recode beagle to %s*.dat + %s*.map... ", outname, outname);
    logprintb();
    fputs("0%", stdout);
    retval = recode_beagle_new_chrom(outname, &(outname_end[5]), marker_exclude, chrom_info_ptr, &marker_uidx, &chrom_fo_idx, &chrom_idx, &chrom_end, bedfile, bed_offset, unfiltered_indiv_ct4, &outfile, &outfile2, writebuf2, header_len);
    if (retval) {
      goto recode_ret_1;
    }
    for (pct = 1; pct <= 100; pct++) {
      loop_end = (((uint64_t)pct) * autosomal_marker_ct) / 100;
      for (; marker_idx < loop_end; marker_uidx++, marker_idx++) {
	if (IS_SET(marker_exclude, marker_uidx)) {
	  marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
          if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	    goto recode_ret_READ_FAIL;
	  }
	}
        if (marker_uidx >= chrom_end) {
	  if (fclose_null(&outfile) || fclose_null(&outfile2)) {
	    goto recode_ret_WRITE_FAIL;
	  }
	  chrom_fo_idx++;
          retval = recode_beagle_new_chrom(outname, &(outname_end[5]), marker_exclude, chrom_info_ptr, &marker_uidx, &chrom_fo_idx, &chrom_idx, &chrom_end, bedfile, bed_offset, unfiltered_indiv_ct4, &outfile, &outfile2, writebuf2, header_len);
	  if (retval) {
	    goto recode_ret_1;
	  }
	}
	if (load_and_collapse(bedfile, (uintptr_t*)loadbuf, unfiltered_indiv_ct, loadbuf_collapsed, indiv_ct, indiv_exclude, IS_SET(marker_reverse, marker_uidx))) {
	  goto recode_ret_READ_FAIL;
	}
	cptr = &(marker_ids[marker_uidx * max_marker_id_len]);
        if (fputs_checked("M ", outfile)) {
	  goto recode_ret_WRITE_FAIL;
	}
        fputs(cptr, outfile);
        putc(' ', outfile);
	if (fputs_checked(cptr, outfile2)) {
	  goto recode_ret_WRITE_FAIL;
	}
	putc('\t', outfile2);
        wbufptr = uint32_writex(tbuf, marker_pos[marker_uidx], '\t');
	fwrite(tbuf, 1, wbufptr - tbuf, outfile2);
	aptr = mk_allele_ptrs[2 * marker_uidx];
	aptr2 = mk_allele_ptrs[2 * marker_uidx + 1];
	alen = strlen(aptr);
	alen2 = strlen(aptr2);
	wbufptr = memcpyax(&(cur_mk_allelesx[0][1]), aptr, alen, ' ');
	memcpy(wbufptr, aptr, alen);
	cmalen[0] = 2 * alen + 2;
	wbufptr = memcpyax(&(cur_mk_allelesx[2][1]), aptr, alen, ' ');
	memcpy(wbufptr, aptr2, alen2);
	fputs(replace_if_zstr(aptr, "X"), outfile2);
	putc('\t', outfile2);
	fputs(replace_if_zstr(aptr2, "X"), outfile2);
	cmalen[2] = alen + alen2 + 2;
	wbufptr = memcpyax(&(cur_mk_allelesx[3][1]), aptr2, alen2, ' ');
	memcpy(wbufptr, aptr2, alen2);
	cmalen[3] = 2 * alen2 + 2;
	if (putc_checked('\n', outfile2)) {
	  goto recode_ret_WRITE_FAIL;
	}
	ulptr = loadbuf_collapsed;
	ulptr_end = &(loadbuf_collapsed[indiv_ct / BITCT2]);
        shiftmax = BITCT2;
        while (1) {
	  while (ulptr < ulptr_end) {
            cur_word = *ulptr++;
            for (shiftval = 0; shiftval < shiftmax; shiftval++) {
	      ulii = cur_word & 3;
	      fwrite(cur_mk_allelesx[ulii], 1, cmalen[ulii], outfile);
	      cur_word >>= 2;
	    }
	  }
	  if (ulptr == loadbuf_collapsed_end) {
	    break;
	  }
	  ulptr_end++;
	  shiftmax = indiv_ct % BITCT2;
	}
	if (putc_checked('\n', outfile)) {
	  goto recode_ret_WRITE_FAIL;
	}
      }
      if (pct < 100) {
	if (pct > 10) {
	  putchar('\b');
	}
	printf("\b\b%u%%", pct);
	fflush(stdout);
      }
    }
    if (fclose_null(&outfile2)) {
      goto recode_ret_WRITE_FAIL;
    }
  } else if (recode_modifier & (RECODE_BIMBAM | RECODE_BIMBAM_1CHR)) {
    if (recode_modifier & RECODE_BIMBAM_1CHR) {
      ii = single_chrom_start(chrom_info_ptr, unfiltered_marker_ct, marker_exclude);
      if (ii == -1) {
        logprint("Error: --recode bimbam-1chr requires a single-chromosome dataset.  Did you mean\n'--recode bimbam'?  (Note the lack of a dash in the middle.)\n");
        goto recode_ret_INVALID_CMDLINE;
      }
      marker_uidx = ((uint32_t)ii);
      if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	goto recode_ret_READ_FAIL;
      }
    }
    memcpy(outname_end, ".recode.", 9);
    sprintf(logbuf, "--recode bimbam%s to %s{geno,pheno,pos}.txt... ", (recode_modifier & RECODE_BIMBAM_1CHR)? "-1chr" : "", outname);
    logprintb();
    fputs("0%", stdout);
    memcpy(&(outname_end[8]), "pos.txt", 8);
    if (fopen_checked(&outfile, outname, "w")) {
      goto recode_ret_OPEN_FAIL;
    }
    writebuf2[0] = ' ';
    ulii = recode_modifier & RECODE_BIMBAM;
    for (marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
      next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
      fputs(&(marker_ids[marker_uidx * max_marker_id_len]), outfile);
      wbufptr = uint32_write(&(writebuf2[1]), marker_pos[marker_uidx]);
      if (ulii) {
	if (marker_uidx >= chrom_end) {
          chrom_idx = get_marker_chrom(chrom_info_ptr, marker_uidx);
          chrom_end = chrom_info_ptr->chrom_end[chrom_idx];
	}
        *wbufptr++ = ' ';
        wbufptr = chrom_name_write(wbufptr, chrom_info_ptr, chrom_idx, 1);
      }
      *wbufptr++ = '\n';
      if (fwrite_checked(writebuf2, wbufptr - writebuf2, outfile)) {
        goto recode_ret_WRITE_FAIL;
      }
    }
    if (fclose_null(&outfile)) {
      goto recode_ret_WRITE_FAIL;
    }
    memcpy(&(outname_end[8]), "pheno.txt", 10);
    if (fopen_checked(&outfile, outname, "w")) {
      goto recode_ret_OPEN_FAIL;
    }
    indiv_uidx = 0;
    if (pheno_c) {
      for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_uidx++, indiv_idx++) {
	next_unset_ul_unsafe_ck(indiv_exclude, &indiv_uidx);
	if (IS_SET(pheno_nm, indiv_uidx)) {
          putc('1' + IS_SET(pheno_c, indiv_uidx), outfile);
	} else {
	  fputs(output_missing_pheno, outfile);
	}
        if (putc_checked('\n', outfile)) {
          goto recode_ret_WRITE_FAIL;
	}
      }
    } else {
      for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_uidx++, indiv_idx++) {
	next_unset_ul_unsafe_ck(indiv_exclude, &indiv_uidx);
	if (IS_SET(pheno_nm, indiv_uidx)) {
          wbufptr = double_g_write(writebuf2, pheno_d[indiv_uidx]);
	  fwrite(writebuf2, 1, (uintptr_t)(wbufptr - writebuf2), outfile);
	} else {
          fputs(output_missing_pheno, outfile);
	}
	if (putc_checked('\n', outfile)) {
	  goto recode_ret_WRITE_FAIL;
	}
      }
    }
    if (fclose_null(&outfile)) {
      goto recode_ret_WRITE_FAIL;
    }
    memcpy(&(outname_end[8]), "geno.txt", 9);
    if (fopen_checked(&outfile, outname, "w")) {
      goto recode_ret_OPEN_FAIL;
    }
    wbufptr = uint32_writex(writebuf2, indiv_ct, '\n');
    wbufptr = uint32_write(wbufptr, marker_ct);
    wbufptr = memcpya(wbufptr, "\nIND", 4);
    if (fwrite_checked(writebuf2, wbufptr - writebuf2, outfile)) {
      goto recode_ret_WRITE_FAIL;
    }
    indiv_uidx = 0;
    for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_uidx++, indiv_idx++) {
      next_unset_ul_unsafe_ck(indiv_exclude, &indiv_uidx);
      cptr = (char*)memchr(&(person_ids[indiv_uidx * max_person_id_len]), '\t', max_person_id_len);
      putc(',', outfile);
      fputs(&(cptr[1]), outfile);
    }
    if (putc_checked('\n', outfile)) {
      goto recode_ret_WRITE_FAIL;
    }
    marker_uidx = 0;
    marker_idx = 0;
    chrom_fo_idx = 0;
    refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_haploid);
    chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
    writebuf2[0] = ',';
    memcpy(&(writebuf2[4]), ",??", 4);
    writebuf2[8] = ',';
    writebuf2[12] = ',';
    for (pct = 1; pct <= 100; pct++) {
      loop_end = (((uint64_t)pct) * marker_ct) / 100;
      for (; marker_idx < loop_end; marker_uidx++, marker_idx++) {
	if (IS_SET(marker_exclude, marker_uidx)) {
	  marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
          if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
            goto recode_ret_READ_FAIL;
	  }
        }
	if (marker_uidx >= chrom_end) {
          chrom_fo_idx++;
          refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_haploid);
          chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	}
	if (load_and_collapse(bedfile, (uintptr_t*)loadbuf, unfiltered_indiv_ct, loadbuf_collapsed, indiv_ct, indiv_exclude, IS_SET(marker_reverse, marker_uidx))) {
	  goto recode_ret_READ_FAIL;
	}
	if (is_haploid && set_hh_missing) {
	  haploid_fix(hh_exists, indiv_include2, indiv_male_include2, indiv_ct, is_x, is_y, (unsigned char*)loadbuf_collapsed);
	}
        ucc = mk_allele_ptrs[2 * marker_uidx][0];
        ucc2 = mk_allele_ptrs[2 * marker_uidx + 1][0];
        writebuf2[1] = ucc;
        writebuf2[2] = ucc;
        writebuf2[9] = ucc;
        writebuf2[10] = ucc2;
        writebuf2[13] = ucc2;
        writebuf2[14] = ucc2;
	if (fputs_checked(&(marker_ids[marker_uidx * max_marker_id_len]), outfile)) {
	  goto recode_ret_WRITE_FAIL;
	}
	wbufptr = writebuf;
	ulptr = loadbuf_collapsed;
	ulptr_end = &(loadbuf_collapsed[indiv_ct / BITCT2]);
        shiftmax = BITCT2;
        while (1) {
	  while (ulptr < ulptr_end) {
            cur_word = *ulptr++;
            for (shiftval = 0; shiftval < shiftmax; shiftval++) {
	      ulii = cur_word & 3;
	      wbufptr = memcpyl3a(wbufptr, &(writebuf2[4 * ulii]));
	      cur_word >>= 2;
	    }
	  }
	  if (ulptr == loadbuf_collapsed_end) {
	    break;
	  }
	  ulptr_end++;
	  shiftmax = indiv_ct % BITCT2;
	}
        fwrite(writebuf, 1, 3 * indiv_ct, outfile);
        if (putc_checked('\n', outfile)) {
	  goto recode_ret_WRITE_FAIL;
	}
      }
      if (pct < 100) {
	if (pct > 10) {
	  putchar('\b');
	}
	printf("\b\b%u%%", pct);
	fflush(stdout);
      }
    }
  } else if (recode_modifier & (RECODE_FASTPHASE | RECODE_FASTPHASE_1CHR)) {
    if (recode_modifier & RECODE_FASTPHASE) {
      memcpy(outname_end, ".chr-*", 7);
    } else {
      *outname_end = '\0';
    }
    sprintf(logbuf, "--recode fastphase%s to %s.recode.phase.inp... ", (recode_modifier & RECODE_FASTPHASE)? "" : "-1chr", outname);
    logprintb();
    fflush(stdout);
    chrom_fo_idx = 0xffffffffU; // exploit overflow for initialization
    if (recode_modifier & RECODE_FASTPHASE) {
      fputs("[chromosome   ", stdout);
    }
    do {
      chrom_fo_idx++;
      next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
      refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_haploid);
      chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
      ulii = count_chrom_markers(chrom_info_ptr, chrom_idx, marker_exclude);
      if (recode_modifier & RECODE_FASTPHASE) {
        wbufptr = chrom_name_write(&(outname_end[5]), chrom_info_ptr, chrom_idx, zero_extra_chroms);
        if (chrom_idx <= chrom_info_ptr->max_code) {
          printf("\b\b%u] \b", chrom_idx);
	} else {
	  // nonstandard chromosome name
          fputs("\b\b**] \b", stdout);
	}
        fflush(stdout);
      } else {
	wbufptr = outname_end;
      }
      memcpy(wbufptr, ".recode.phase.inp", 18);
      retval = open_and_write_fastphase_header(&outfile, outname, marker_exclude, marker_pos, marker_uidx, ulii, indiv_ct);
      if (retval) {
	goto recode_ret_1;
      }
      marker_uidx_start = marker_uidx;
      if (recode_load_to(loadbuf, bedfile, bed_offset, chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx + 1], 0, ulii, marker_exclude, marker_reverse, &marker_uidx, unfiltered_indiv_ct)) {
	goto recode_ret_READ_FAIL;
      }
      if (set_hh_missing) {
	haploid_fix_multiple(marker_exclude, marker_uidx_start, ulii, chrom_info_ptr, hh_exists, indiv_include2, indiv_male_include2, unfiltered_indiv_ct, unfiltered_indiv_ct4, loadbuf);
      }
      indiv_uidx = 0;
      for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_uidx++, indiv_idx++) {
        next_unset_ul_unsafe_ck(indiv_exclude, &indiv_uidx);
	if (fputs_checked("# ID ", outfile)) {
	  goto recode_ret_WRITE_FAIL;
	}
        cptr = (char*)memchr(&(person_ids[indiv_uidx * max_person_id_len]), '\t', max_person_id_len);
        fputs(&(cptr[1]), outfile);
	putc('\n', outfile);
	bufptr = &(loadbuf[indiv_uidx / 4]);
	shiftval = (indiv_uidx % 4) * 2;
	aptr = writebuf;
	aptr2 = writebuf2;
	uidx_cur = marker_uidx_start;
	marker_idx = 0;
        do {
          uidx_cur = next_unset_unsafe(marker_exclude, uidx_cur);
          uidx_stop = next_set(marker_exclude, marker_uidx, chrom_end);
	  marker_idx += uidx_stop - uidx_cur;
	  do {
            ucc = ((*bufptr) >> shiftval) & 3;
            *aptr++ = "0?01"[ucc];
	    *aptr2++ = "0?11"[ucc];
	    bufptr = &(bufptr[unfiltered_indiv_ct4]);
	  } while (++uidx_cur < uidx_stop);
	} while (marker_idx < ulii);
	fwrite(writebuf, 1, ulii, outfile);
	putc('\n', outfile);
        fwrite(writebuf2, 1, ulii, outfile);
	if (putc_checked('\n', outfile)) {
	  goto recode_ret_WRITE_FAIL;
	}
      }
      if (recode_modifier & RECODE_FASTPHASE) {
        if (chrom_idx > onechar_max) {
          putchar('\b');
	}
      }
      if (fclose_null(&outfile)) {
        goto recode_ret_WRITE_FAIL;
      }
    } while (chrom_fo_idx < last_chrom_fo_idx);
    if (recode_modifier & RECODE_FASTPHASE) {
      fputs("\b\b\b\b\b\b\b\b\b\b\b\b\b\b               \b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", stdout);
    }
  } else if (recode_modifier & (RECODE_LGEN | RECODE_LGEN_REF)) {
    if (lgen_ref) {
      strcpy(outname_end, ".ref");
      if (fopen_checked(&outfile2, outname, "w")) {
	goto recode_ret_OPEN_FAIL;
      }
    }
    strcpy(outname_end, ".lgen");
    if (fopen_checked(&outfile, outname, "w")) {
      goto recode_ret_OPEN_FAIL;
    }
    if (delimiter == ' ') {
      indiv_delim_convert(unfiltered_indiv_ct, indiv_exclude, indiv_ct, person_ids, max_person_id_len, '\t', ' ');
    } else {
      if (!(recode_modifier & RECODE_DELIMX)) {
	delim2 = ' ';
      }
    }
    sprintf(logbuf, "--recode lgen%s to %s + .map + .fam%s... ", lgen_ref? "-ref" : "", outname, lgen_ref? " + .ref" : "");
    logprintb();
    fputs("0%", stdout);
    fflush(stdout);
    ref_word = 4;
    writebuf[0] = delimiter;
    for (pct = 1; pct <= 100; pct++) {
      loop_end = (((uint64_t)pct) * marker_ct) / 100;
      for (; marker_idx < loop_end; marker_uidx++, marker_idx++) {
	if (IS_SET(marker_exclude, marker_uidx)) {
	  marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
	  if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	    goto recode_ret_READ_FAIL;
	  }
	}
	if (marker_uidx >= chrom_end) {
	  chrom_fo_idx++;
	  refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_haploid);
	}
        if (load_and_collapse(bedfile, (uintptr_t*)loadbuf, unfiltered_indiv_ct, loadbuf_collapsed, indiv_ct, indiv_exclude, IS_SET(marker_reverse, marker_uidx))) {
	  goto recode_ret_READ_FAIL;
	}
	if (is_haploid && set_hh_missing) {
	  haploid_fix(hh_exists, indiv_include2, indiv_male_include2, indiv_ct, is_x, is_y, (unsigned char*)loadbuf_collapsed);
	}
	wbufptr = &(marker_ids[marker_uidx * max_marker_id_len]);
	cptr = strcpya(&(writebuf[1]), wbufptr);
	cptr = memseta(cptr, delimiter, 2);
	ulii = (uintptr_t)(cptr - writebuf);
	alen = 2 * max_marker_allele_len + ulii;
	cur_mk_allelesx[0] = cptr;
	cur_mk_allelesx[1] = &(cptr[alen]);
	cur_mk_allelesx[2] = &(cptr[2 * alen]);
	cur_mk_allelesx[3] = &(cptr[3 * alen]);
	memcpy(&(writebuf[alen]), writebuf, ulii);
	memcpy(&(writebuf[2 * alen]), writebuf, ulii);
	memcpy(&(writebuf[3 * alen]), writebuf, ulii);
	cur_mk_allelesx[1][0] = output_missing_geno;
	cur_mk_allelesx[1][1] = delim2;
	cur_mk_allelesx[1][2] = output_missing_geno;
	cur_mk_allelesx[1][3] = '\n';
	init_recode_cmax(mk_allele_ptrs[2 * marker_uidx], mk_allele_ptrs[2 * marker_uidx + 1], cur_mk_allelesx, cmalen, '\n', delim2);
	cmalen[0] += ulii;
	cmalen[1] = 4 + ulii;
	cmalen[2] += ulii;
	cmalen[3] += ulii;
	if (lgen_ref) {
	  aptr = mk_allele_ptrs[2 * marker_uidx + 1];
	  aptr2 = mk_allele_ptrs[2 * marker_uidx];
	  if ((aptr[0] != missing_geno) || aptr[1] || (aptr2[0] != missing_geno) || aptr2[1]) {
	    fputs(wbufptr, outfile2);
	    if ((aptr[0] != missing_geno) || aptr[1]) {
	      putc(delimiter, outfile2);
	      fputs(aptr, outfile2);
	    }
	    if ((aptr2[0] != missing_geno) || aptr2[1]) {
	      putc(delimiter, outfile2);
	      fputs(aptr2, outfile2);
	    }
	    if (putc_checked('\n', outfile2)) {
	      goto recode_ret_WRITE_FAIL;
	    }
	  }
	  ref_word = 3;
	}

	ulptr = loadbuf_collapsed;
	ulptr_end = &(loadbuf_collapsed[indiv_ct / BITCT2]);
	shiftmax = BITCT2;
	indiv_idx = 0;
	while (1) {
	  while (ulptr < ulptr_end) {
	    cur_word = *ulptr++;
	    for (; indiv_idx < shiftmax; indiv_idx++) {
	      ulii = cur_word & 3;
	      if (ulii != ref_word) {
		fputs(&(person_ids_collapsed[indiv_idx * max_person_id_len]), outfile);
		fwrite(&(writebuf[ulii * alen]), 1, cmalen[ulii], outfile);
	      }
	      cur_word >>= 2;
	    }
	    shiftmax += BITCT2;
	  }
	  if (indiv_idx == indiv_ct) {
	    break;
	  }
	  ulptr_end++;
	  shiftmax = indiv_ct;
	}
      }
      if (pct < 100) {
	if (pct > 10) {
	  putchar('\b');
	}
	printf("\b\b%u%%", pct);
	fflush(stdout);
      }
    }
    indiv_delim_convert(unfiltered_indiv_ct, indiv_exclude, indiv_ct, person_ids, max_person_id_len, ' ', '\t');
  } else if (recode_modifier & (RECODE_A | RECODE_AD)) {
    memcpy(outname_end, ".raw", 5);
    if (wkspace_left < ((uint64_t)unfiltered_indiv_ct4) * marker_ct) {
      goto recode_ret_NO_MULTIPASS_YET;
    }
    if (fopen_checked(&outfile, outname, "w")) {
      goto recode_ret_OPEN_FAIL;
    }
    if (fputs_checked("FID IID PAT MAT SEX PHENOTYPE", outfile)) {
      goto recode_ret_WRITE_FAIL;
    }
    for (marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
      next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
      cptr = &(marker_ids[marker_uidx * max_marker_id_len]);
      if (allele_missing && (allele_missing[marker_uidx])) {
	aptr = allele_missing[marker_uidx];
      } else {
	aptr = mk_allele_ptrs[2 * marker_uidx + IS_NONNULL_AND_SET(recode_allele_reverse, marker_uidx)];
      }
      putc(' ', outfile);
      fputs(cptr, outfile);
      putc('_', outfile);
      fputs(aptr, outfile);
      if (recode_modifier & RECODE_AD) {
	putc(' ', outfile);
	fputs(cptr, outfile);
	fputs("_HET", outfile);
      }
      if (ferror(outfile)) {
	goto recode_ret_WRITE_FAIL;
      }
    }
    if (putc_checked('\n', outfile)) {
      goto recode_ret_WRITE_FAIL;
    }
    marker_uidx = 0;
    sprintf(logbuf, "--recode A%s to %s... ", (recode_modifier & RECODE_AD)? "D" : "", outname);
    logprintb();
    if (recode_load_to(loadbuf, bedfile, bed_offset, unfiltered_marker_ct, 0, marker_ct, marker_exclude, marker_reverse, &marker_uidx, unfiltered_indiv_ct)) {
      goto recode_ret_READ_FAIL;
    }
    if (set_hh_missing) {
      haploid_fix_multiple(marker_exclude, 0, marker_ct, chrom_info_ptr, hh_exists, indiv_include2, indiv_male_include2, unfiltered_indiv_ct, unfiltered_indiv_ct4, loadbuf);
    }
    fputs("0%", stdout);
    indiv_uidx = 0;
    indiv_idx = 0;
    for (pct = 1; pct <= 100; pct++) {
      loop_end = ((uint64_t)pct * indiv_ct) / 100;
      for (; indiv_idx < loop_end; indiv_uidx++, indiv_idx++) {
	indiv_uidx = next_unset_ul_unsafe(indiv_exclude, indiv_uidx);
	if (recode_write_first_cols(outfile, indiv_uidx, delimiter, person_ids, max_person_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, sex_nm, sex_male, pheno_nm, pheno_c, pheno_d, output_missing_pheno)) {
	  goto recode_ret_WRITE_FAIL;
	}
	bufptr = &(loadbuf[indiv_uidx / 4]);
	wbufptr = writebuf;
	shiftval = (indiv_uidx % 4) * 2;
	marker_uidx = 0;
	marker_idx = 0;
	if (recode_modifier & RECODE_A) {
	  do {
	    marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
	    ulii = next_set_ul(marker_exclude, marker_uidx, unfiltered_marker_ct);
	    marker_idx += ulii - marker_uidx;
	    do {
	      ucc = ((*bufptr) >> shiftval) & 3;
	      if (allele_missing && allele_missing[marker_uidx]) {
		*wbufptr++ = "0N00"[ucc];
	      } else {
		*wbufptr++ = "2N10"[ucc];
	      }
	      if (ucc == 1) {
		*wbufptr++ = 'A';
	      }
	      *wbufptr++ = delimiter;
	      bufptr = &(bufptr[unfiltered_indiv_ct4]);
	    } while (++marker_uidx < ulii);
	  } while (marker_idx < marker_ct);
	} else {
	  do {
	    marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
	    ulii = next_set_ul(marker_exclude, marker_uidx, unfiltered_marker_ct);
	    marker_idx += ulii - marker_uidx;
	    do {
	      ucc = ((*bufptr) >> shiftval) & 3;
	      if (ucc != 1) {
		wbufptr = memcpya(wbufptr, &(writebuf2[4 * ((allele_missing && allele_missing[marker_uidx])? 3 : ucc)]), 4);
	      } else {
		wbufptr = memcpya(wbufptr, &(writebuf2[16]), 6);
	      }
	      bufptr = &(bufptr[unfiltered_indiv_ct4]);
	    } while (++marker_uidx < ulii);
	  } while (marker_idx < marker_ct);
	}
	wbufptr[-1] = '\n';
	ulii = (uintptr_t)(wbufptr - writebuf);
	if (fwrite_checked(writebuf, ulii, outfile)) {
	  goto recode_ret_WRITE_FAIL;
	}
      }
      if (pct < 100) {
	if (pct > 10) {
	  putchar('\b');
	}
	printf("\b\b%u%%", pct);
	fflush(stdout);
      }
    }
  } else if (recode_modifier & (RECODE_LIST | RECODE_RLIST)) {
    // todo: fix
    strcpy(outname_end, rlist? ".rlist" : ".list");
    if (fopen_checked(&outfile, outname, "w")) {
      goto recode_ret_OPEN_FAIL;
    }
    if (delimiter != '\t') {
      indiv_delim_convert(unfiltered_indiv_ct, indiv_exclude, indiv_ct, person_ids, max_person_id_len, '\t', ' ');
    }
    if (rlist) {
      sprintf(logbuf, "--recode rlist to %s + .map + .fam... ", outname);
    } else {
      sprintf(logbuf, "--recode list to %s... ", outname);
    }
    logprintb();
    fputs("0%", stdout);
    fflush(stdout);
    cur_mk_allelesx[1][0] = missing_geno;
    cur_mk_allelesx[1][1] = delimiter;
    cur_mk_allelesx[1][2] = missing_geno;
    cmalen[1] = 3;
    for (pct = 1; pct <= 100; pct++) {
      loop_end = (((uint64_t)pct) * marker_ct) / 100;
      for (; marker_idx < loop_end; marker_uidx++, marker_idx++) {
	if (IS_SET(marker_exclude, marker_uidx)) {
	  marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
	  if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	    goto recode_ret_READ_FAIL;
	  }
	}
	if (marker_uidx >= chrom_end) {
	  chrom_fo_idx++;
	  refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_haploid);
	  chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	}
	if (load_and_collapse(bedfile, (uintptr_t*)loadbuf, unfiltered_indiv_ct, loadbuf_collapsed, indiv_ct, indiv_exclude, IS_SET(marker_reverse, marker_uidx))) {
	  goto recode_ret_READ_FAIL;
	}
	if (is_haploid && set_hh_missing) {
          haploid_fix(hh_exists, indiv_include2, indiv_male_include2, indiv_ct, is_x, is_y, (unsigned char*)loadbuf_collapsed);
	}
	init_recode_cmax(mk_allele_ptrs[2 * marker_uidx], mk_allele_ptrs[2 * marker_uidx + 1], cur_mk_allelesx, cmalen, '\0', delimiter);
	cmalen[0] -= 1;
	cmalen[2] -= 1;
	cmalen[3] -= 1;
	aptr = &(marker_ids[marker_uidx * max_marker_id_len]);
	alen = strlen(aptr);
	for (ulii = 0; ulii < 4; ulii++) {
	  wbufptr = writebufl[ulii];
	  if (!rlist) {
	    wbufptr = chrom_name_write(wbufptr, chrom_info_ptr, chrom_idx, zero_extra_chroms);
	    *wbufptr++ = delimiter;
	  }
	  wbufptr = memcpyax(wbufptr, aptr, alen, delimiter);
	  if (rlist) {
	    switch (ulii) {
	    case 0:
	    case 3:
	      wbufptr = memcpyl3a(wbufptr, "HOM");
	      break;
	    case 1:
	      wbufptr = memcpyl3a(wbufptr, "NIL");
	      break;
	    case 2:
	      wbufptr = memcpyl3a(wbufptr, "HET");
	      break;
	    }
	    *wbufptr++ = delimiter;
	    wbufptr = memcpya(wbufptr, cur_mk_allelesx[ulii], cmalen[ulii]);
	  } else {
	    *wbufptr++ = cur_mk_allelesx[ulii][0];
	    *wbufptr++ = cur_mk_allelesx[ulii][2];
	  }
	  writebuflp[ulii] = wbufptr;
	}
	ulptr = loadbuf_collapsed;
	ulptr_end = &(loadbuf_collapsed[indiv_ct / BITCT2]);
	indiv_idx = 0;
	indiv_uidx = BITCT2; // repurposed as stop value
	if (rlist) {
	  for (ulii = 0; ulii < 4; ulii++) {
	    writebuflps[ulii] = writebuflp[ulii];
	  }
	  while (1) {
            while (ulptr < ulptr_end) {
	      cur_word = *ulptr++;
	      for (; indiv_idx < indiv_uidx; indiv_idx++, cur_word >>= 2) {
                ulii = cur_word & 3;
		if (ulii != 3) {
		  *(writebuflp[ulii]++) = delimiter;
		  writebuflp[ulii] = strcpya(writebuflp[ulii], &(person_ids_collapsed[indiv_idx * max_person_id_len]));
		}
	      }
	      indiv_uidx += BITCT2;
	    }
	    if (ulptr == loadbuf_collapsed_end) {
	      break;
	    }
            ulptr_end++;
	    indiv_uidx = indiv_ct;
	  }
	  if (writebuflp[2] != writebuflps[2]) {
	    *(writebuflp[2]++) = '\n';
	    fwrite(writebufl[2], 1, writebuflp[2] - writebufl[2], outfile);
	  }
	  if (writebuflp[0] != writebuflps[0]) {
	    *(writebuflp[0]++) = '\n';
	    fwrite(writebufl[0], 1, writebuflp[0] - writebufl[0], outfile);
	  }
	  if (writebuflp[1] != writebuflps[1]) {
	    *(writebuflp[1]++) = '\n';
	    fwrite(writebufl[1], 1, writebuflp[1] - writebufl[1], outfile);
	  }
	} else {
	  while (1) {
            while (ulptr < ulptr_end) {
	      cur_word = *ulptr++;
	      for (; indiv_idx < indiv_uidx; indiv_idx++, cur_word >>= 2) {
                ulii = cur_word & 3;
		*(writebuflp[ulii]++) = delimiter;
		writebuflp[ulii] = strcpya(writebuflp[ulii], &(person_ids_collapsed[indiv_idx * max_person_id_len]));
	      }
	      indiv_uidx += BITCT2;
	    }
	    if (ulptr == loadbuf_collapsed_end) {
	      break;
	    }
            ulptr_end++;
	    indiv_uidx = indiv_ct;
	  }
	  *(writebuflp[0]++) = '\n';
	  *(writebuflp[1]++) = '\n';
	  *(writebuflp[2]++) = '\n';
	  *(writebuflp[3]++) = '\n';
	  fwrite(writebufl[0], 1, writebuflp[0] - writebufl[0], outfile);
	  fwrite(writebufl[2], 1, writebuflp[2] - writebufl[2], outfile);
	  fwrite(writebufl[3], 1, writebuflp[3] - writebufl[3], outfile);
	  fwrite(writebufl[1], 1, writebuflp[1] - writebufl[1], outfile);
	}
	if (ferror(outfile)) {
	  goto recode_ret_WRITE_FAIL;
	}
      }
      if (pct < 100) {
	if (pct > 10) {
	  putchar('\b');
	}
	printf("\b\b%u%%", pct);
	fflush(stdout);
      }
    }
    if (delimiter != '\t') {
      indiv_delim_convert(unfiltered_indiv_ct, indiv_exclude, indiv_ct, person_ids, max_person_id_len, ' ', '\t');
    }
  } else if (recode_modifier & (RECODE_HV | RECODE_HV_1CHR)) {
    if (wkspace_left < ((uint64_t)unfiltered_indiv_ct4) * max_chrom_size) {
      goto recode_ret_NO_MULTIPASS_YET;
    }
    if (recode_modifier & RECODE_HV) {
      memcpy(outname_end, ".chr-*", 7);
    } else {
      *outname_end = '\0';
    }
    sprintf(logbuf, "--recode HV%s to %s.ped + .info... ", (recode_modifier & RECODE_HV)? "" : "-1chr", outname);
    logprintb();
    fflush(stdout);
    if (recode_modifier & RECODE_HV) {
      fputs("[chromosome   ", stdout);
    }
    chrom_fo_idx = 0xffffffffU;
    do {
      chrom_fo_idx++;
      marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx);
      refresh_chrom_info(chrom_info_ptr, marker_uidx, &chrom_end, &chrom_fo_idx, &is_x, &is_y, &is_haploid);
      chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
      ulii = count_chrom_markers(chrom_info_ptr, chrom_idx, marker_exclude);
      if (recode_modifier & RECODE_HV) {
        wbufptr = chrom_name_write(&(outname_end[5]), chrom_info_ptr, chrom_idx, zero_extra_chroms);
        if (chrom_idx <= chrom_info_ptr->max_code) {
          printf("\b\b%u] \b", chrom_idx);
        } else {
          // nonstandard chromosome name
          fputs("\b\b**] \b", stdout);
        }
        fflush(stdout);
      } else {
        wbufptr = outname_end;
      }
      memcpy(wbufptr, ".ped", 5);
      if (fopen_checked(&outfile, outname, "w")) {
	goto recode_ret_OPEN_FAIL;
      }
      marker_uidx_start = marker_uidx;
      if (recode_load_to(loadbuf, bedfile, bed_offset, chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx + 1], 0, ulii, marker_exclude, marker_reverse, &marker_uidx, unfiltered_indiv_ct)) {
	goto recode_ret_READ_FAIL;
      }
      if (set_hh_missing) {
        haploid_fix_multiple(marker_exclude, marker_uidx_start, ulii, chrom_info_ptr, hh_exists, indiv_include2, indiv_male_include2, unfiltered_indiv_ct, unfiltered_indiv_ct4, loadbuf);
      }
      indiv_uidx = 0;
      if (recode_compound) {
        writebuf[ulii * 3 - 1] = '\n';
      }
      if (write_ped_lines(outfile, loadbuf, marker_exclude, marker_uidx_start, ulii, unfiltered_indiv_ct4, indiv_exclude, &indiv_uidx, 0, indiv_ct, recode_compound, mk_allele_ptrs, delimiter, delim2, person_ids, max_person_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, sex_nm, sex_male, pheno_nm, pheno_c, pheno_d, output_missing_geno, output_missing_pheno, writebuf)) {
	goto recode_ret_WRITE_FAIL;
      }
      if (recode_compound) {
	writebuf[ulii * 3 - 1] = delimiter;
      }
      if (fclose_null(&outfile)) {
	goto recode_ret_WRITE_FAIL;
      }
      memcpy(wbufptr, ".info", 6);
      if (fopen_checked(&outfile, outname, "w")) {
	goto recode_ret_OPEN_FAIL;
      }
      if (write_haploview_map(outfile, marker_exclude, marker_uidx_start, ulii, marker_ids, max_marker_id_len, marker_pos)) {
        goto recode_ret_WRITE_FAIL;
      }
      if (fclose_null(&outfile)) {
	goto recode_ret_WRITE_FAIL;
      }
      if (recode_modifier & RECODE_HV) {
        if (chrom_idx > onechar_max) {
          putchar('\b');
	}
      }
    } while (chrom_fo_idx < last_chrom_fo_idx);
    if (recode_modifier & RECODE_HV) {
      fputs("\b\b\b\b\b\b\b\b\b\b\b\b\b\b               \b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", stdout);
    }
  } else if (recode_modifier & RECODE_STRUCTURE) {
    memcpy(outname_end, ".recode.strct_in", 17);
    if (wkspace_left < ((uint64_t)unfiltered_indiv_ct4) * marker_ct) {
      goto recode_ret_NO_MULTIPASS_YET;
    }
    if (fopen_checked(&outfile, outname, "w")) {
      goto recode_ret_OPEN_FAIL;
    }
    sprintf(logbuf, "--recode structure to %s... ", outname);
    logprintb();
    for (marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
      next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
      fputs(&(marker_ids[marker_uidx * max_marker_id_len]), outfile);
      putc(' ', outfile);
    }
    if (putc_checked('\n', outfile)) {
      goto recode_ret_WRITE_FAIL;
    }
    marker_uidx = 0;
    chrom_fo_idx = 0xffffffffU; // exploit overflow
    chrom_end = 0;
    for (marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
      next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
      if (marker_uidx >= chrom_end) {
	do {
          chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[(++chrom_fo_idx) + 1];
	} while (marker_uidx >= chrom_end);
	fputs("-1 ", outfile);
      } else {
        wbufptr = uint32_writex(writebuf2, marker_pos[marker_uidx] - last_pos, ' ');
        fwrite(writebuf2, 1, wbufptr - writebuf2, outfile);
      }
      last_pos = marker_pos[marker_uidx];
    }
    if (putc_checked('\n', outfile)) {
      goto recode_ret_WRITE_FAIL;
    }
    marker_uidx = 0;
    if (recode_load_to(loadbuf, bedfile, bed_offset, unfiltered_marker_ct, 0, marker_ct, marker_exclude, marker_reverse, &marker_uidx, unfiltered_indiv_ct)) {
      goto recode_ret_READ_FAIL;
    }
    if (set_hh_missing) {
      haploid_fix_multiple(marker_exclude, 0, marker_ct, chrom_info_ptr, hh_exists, indiv_include2, indiv_male_include2, unfiltered_indiv_ct, unfiltered_indiv_ct4, loadbuf);
    }
    indiv_uidx = 0;
    indiv_idx = 0;
    fputs("0%", stdout);
    last_pos = 0; // now unique FIDs encountered so far
    memcpy(writebuf2, " 1 1 0 0 1 2 2 2", 16);
    for (pct = 1; pct <= 100; pct++) {
      loop_end = (((uint64_t)pct) * indiv_ct) / 100;
      for (; indiv_idx < loop_end; indiv_uidx++, indiv_idx++) {
	next_unset_ul_unsafe_ck(indiv_exclude, &indiv_uidx);
	bufptr = &(loadbuf[indiv_uidx / 4]);
	shiftval = (indiv_uidx % 4) * 2;
        cptr = &(person_ids[indiv_uidx * max_person_id_len]);
        aptr = (char*)memchr(cptr, '\t', max_person_id_len);
        fputs(&(aptr[1]), outfile);
	memcpy(tbuf, cptr, aptr - cptr);
        tbuf[(uintptr_t)(aptr - cptr)] = '\0';
        ii = bsearch_str(tbuf, writebuf3, max_fid_len, 0, fid_ct - 1);
	cur_fid = fid_map[(uint32_t)ii];
	if (!cur_fid) {
	  cur_fid = ++last_pos;
          fid_map[(uint32_t)ii] = last_pos;
	}
	tbuf[0] = ' ';
        wbufptr = uint32_write(&(tbuf[1]), cur_fid);
        fwrite(tbuf, 1, wbufptr - tbuf, outfile);
	marker_uidx = 0;
	wbufptr = writebuf;
        for (marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
	  next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
          ucc = ((*bufptr) >> shiftval) & 3;
          wbufptr = memcpya(wbufptr, &(writebuf2[4 * ucc]), 4);
	  bufptr = &(bufptr[unfiltered_indiv_ct4]);
	}
	fwrite(writebuf, 1, wbufptr - writebuf, outfile);
	if (putc_checked('\n', outfile)) {
          goto recode_ret_WRITE_FAIL;
	}
      }
      if (pct < 100) {
	if (pct > 10) {
	  putchar('\b');
	}
	printf("\b\b%u%%", pct);
	fflush(stdout);
      }
    }
  } else {
    memcpy(outname_end, ".ped", 5);
    if (wkspace_left < ((uint64_t)unfiltered_indiv_ct4) * marker_ct) {
      goto recode_ret_NO_MULTIPASS_YET;
    }
    if (fopen_checked(&outfile, outname, "w")) {
      goto recode_ret_OPEN_FAIL;
    }
    sprintf(logbuf, "--recode to %s + .map... ", outname);
    logprintb();
    if (recode_load_to(loadbuf, bedfile, bed_offset, unfiltered_marker_ct, 0, marker_ct, marker_exclude, marker_reverse, &marker_uidx, unfiltered_indiv_ct)) {
      goto recode_ret_READ_FAIL;
    }
    if (set_hh_missing) {
      haploid_fix_multiple(marker_exclude, 0, marker_ct, chrom_info_ptr, hh_exists, indiv_include2, indiv_male_include2, unfiltered_indiv_ct, unfiltered_indiv_ct4, loadbuf);
    }
    indiv_uidx = 0;
    indiv_idx = 0;
    fputs("0%", stdout);
    for (pct = 1; pct <= 100; pct++) {
      loop_end = (((uint64_t)pct) * indiv_ct) / 100;
      if (write_ped_lines(outfile, loadbuf, marker_exclude, 0, marker_ct, unfiltered_indiv_ct4, indiv_exclude, &indiv_uidx, indiv_idx, loop_end, recode_compound, mk_allele_ptrs, delimiter, delim2, person_ids, max_person_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, sex_nm, sex_male, pheno_nm, pheno_c, pheno_d, output_missing_geno, output_missing_pheno, writebuf)) {
	goto recode_ret_WRITE_FAIL;
      }
      indiv_idx = loop_end;
      if (pct < 100) {
	if (pct > 10) {
	  putchar('\b');
	}
	printf("\b\b%u%%", pct);
	fflush(stdout);
      }
    }
  }
  if (outfile) {
    if (fclose_null(&outfile)) {
      goto recode_ret_WRITE_FAIL;
    }
  }

  if (recode_modifier & (RECODE_TRANSPOSE | RECODE_LGEN | RECODE_LGEN_REF | RECODE_RLIST)) {
    if (recode_modifier & RECODE_TRANSPOSE) {
      strcpy(outname_end, ".tfam");
    } else {
      strcpy(outname_end, ".fam");
      if ((recode_modifier & (RECODE_TAB | RECODE_DELIMX)) == RECODE_TAB) {
	delimiter = ' ';
      }
    }
    retval = write_fam(outname, unfiltered_indiv_ct, indiv_exclude, indiv_ct, person_ids, max_person_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, sex_nm, sex_male, pheno_nm, pheno_c, pheno_d, output_missing_pheno, delimiter, NULL);
    if (retval) {
      goto recode_ret_1;
    }
  }

  if (!(recode_modifier & (RECODE_TRANSPOSE | RECODE_23 | RECODE_A | RECODE_AD | RECODE_BEAGLE | RECODE_BIMBAM | RECODE_BIMBAM_1CHR | RECODE_FASTPHASE | RECODE_FASTPHASE_1CHR | RECODE_HV | RECODE_HV_1CHR | RECODE_LIST | RECODE_STRUCTURE | RECODE_VCF))) {
    strcpy(outname_end, ".map");
    retval = write_map_or_bim(outname, marker_exclude, marker_ct, marker_ids, max_marker_id_len, marker_cms, marker_pos, NULL, ((recode_modifier & (RECODE_TAB | RECODE_DELIMX)) == RECODE_DELIMX)? ' ' : '\t', zero_extra_chroms, chrom_info_ptr);
    if (retval) {
      goto recode_ret_1;
    }
  }
  if (!(recode_modifier & (RECODE_23 | RECODE_FASTPHASE | RECODE_FASTPHASE_1CHR | RECODE_HV | RECODE_HV_1CHR))) {
    fputs("\b\b\b", stdout);
  }
  logprint("done.\n");
  while (0) {
  recode_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  recode_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  recode_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  recode_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  recode_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  recode_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  recode_ret_NO_MULTIPASS_YET:
    // probably want to implement this later
    logprint("Error: --recode does not yet support multipass recoding of very large files;\ncontact the " PROG_NAME_CAPS " developers if you need this.\nFor now, you can try using a machine with more memory, and/or split the file\ninto smaller pieces and recode them separately.\n");
    retval = RET_CALC_NOT_YET_SUPPORTED;
    break;
  recode_ret_ALL_MARKERS_EXCLUDED:
    retval = RET_ALL_MARKERS_EXCLUDED;
    break;
  }
 recode_ret_1:
  wkspace_reset(wkspace_mark);
  fclose_cond(outfile2);
  fclose_cond(outfile);
  return retval;
}

// .fam
typedef struct ll_entry_struct {
  struct ll_entry_struct* next;
  double pheno;
  uint32_t orig_order;
  char idstr[];
} Ll_entry;

// .bim
typedef struct ll_entry2_struct {
  struct ll_entry2_struct* next;
  int64_t pos;
  double cm;
  char* allele[2];
  char idstr[];
} Ll_entry2;

static inline int32_t idmatch(char* idtab, char* id0, uint32_t id0_len_p1, char* id1, uint32_t id1_len_p1) {
  return (!(memcmp(idtab, id0, id0_len_p1) || memcmp(&(idtab[id0_len_p1]), id1, id1_len_p1)));
}

static inline uint32_t hashval(char* id1, uint32_t id1_len, char* id2, uint32_t id2_len) {
  // just interpret as little-endian number and take modulo HASHSIZE_S for now.
  // may want to upgrade this later
  unsigned char* ucptr = (unsigned char*)id1;
  unsigned char* ucp_end = &(ucptr[id1_len]);
  uint32_t vv = *ucptr;
  while (++ucptr != ucp_end) {
    vv = ((vv << 8) + (*ucptr)) % HASHSIZE_S;
  }
  vv = ((vv << 8) + 9) % HASHSIZE_S;
  ucptr = (unsigned char*)id2;
  ucp_end = &(ucptr[id2_len]);
  do {
    vv = ((vv << 8) + (*ucptr++)) % HASHSIZE_S;
  } while (ucptr != ucp_end);
  return vv;
}

static inline uint32_t hashval2(char* idstr, uint32_t idlen) {
  unsigned char* ucptr = (unsigned char*)idstr;
  unsigned char* ucp_end = &(ucptr[idlen]);
  uint32_t vv = *ucptr;
  while (++ucptr != ucp_end) {
    vv = ((vv << 8) + (*ucptr)) % HASHSIZE;
  }
  return vv;
}

static inline Ll_entry* top_alloc_ll(uintptr_t* topsize_ptr, uint32_t size) {
  return (Ll_entry*)top_alloc(topsize_ptr, size + sizeof(Ll_entry));
}

static inline Ll_entry2* top_alloc_ll2(uintptr_t* topsize_ptr, uint32_t size) {
  return (Ll_entry2*)top_alloc(topsize_ptr, size + sizeof(Ll_entry2));
}

int32_t merge_fam_id_scan(char* bedname, char* famname, uintptr_t* max_person_id_len_ptr, uint32_t* max_person_full_len_ptr, uint32_t* is_dichot_pheno_ptr, Ll_entry** htable, uintptr_t* topsize_ptr, uint64_t* tot_indiv_ct_ptr, uint32_t* ped_buflen_ptr, uint32_t* cur_indiv_ct_ptr, uint32_t* orig_idx_ptr) {
  uintptr_t max_person_id_len = *max_person_id_len_ptr;
  uint32_t max_person_full_len = *max_person_full_len_ptr;
  uint32_t is_dichot_pheno = *is_dichot_pheno_ptr;
  uintptr_t topsize = *topsize_ptr;
  uint64_t tot_indiv_ct = *tot_indiv_ct_ptr;
  uint32_t orig_idx = *orig_idx_ptr;
  uint32_t cur_indiv_ct = 0;
  FILE* infile = NULL;
  uint32_t text_file = 0;
  int32_t retval = 0;
  uint32_t col1_len;
  uint32_t col2_len;
  uint32_t col3_len;
  uint32_t col4_len;
  uint32_t tot_len;
  uintptr_t ulii;
  uint32_t uii;
  Ll_entry** ll_pptr;
  Ll_entry* ll_ptr;
  char* col2_start_ptr;
  char* col3_start_ptr;
  char* col4_start_ptr;
  char* col5_start_ptr;
  char* col6_start_ptr;
  char* col1_end_ptr;
  char* col2_end_ptr;
  char* col1_start_ptr;
  char* wptr;
  double pheno;
  char cc;
  if (!famname) {
    famname = bedname;
    text_file = 1;
  }
  if (fopen_checked(&infile, famname, "r")) {
    goto merge_fam_id_scan_ret_OPEN_FAIL;
  }
  tbuf[MAXLINELEN - 1] = ' ';
  while (fgets(tbuf, MAXLINELEN, infile)) {
    col1_start_ptr = skip_initial_spaces(tbuf);
    cc = *col1_start_ptr;
    if (!is_eoln_or_comment(cc)) {
      col1_end_ptr = item_endnn(col1_start_ptr);
      col1_len = col1_end_ptr - col1_start_ptr;
      col2_start_ptr = skip_initial_spaces(col1_end_ptr);
      col2_end_ptr = item_endnn(col2_start_ptr);
      col2_len = col2_end_ptr - col2_start_ptr;
      col3_start_ptr = skip_initial_spaces(col2_end_ptr);
      col4_start_ptr = item_endnn(col3_start_ptr);
      col3_len = col4_start_ptr - col3_start_ptr;
      col4_start_ptr = skip_initial_spaces(col4_start_ptr);
      col5_start_ptr = item_endnn(col4_start_ptr);
      col4_len = col5_start_ptr - col4_start_ptr;
      col5_start_ptr = skip_initial_spaces(col5_start_ptr);
      col6_start_ptr = item_endnn(col5_start_ptr);
      uii = col6_start_ptr - col5_start_ptr;
      if (uii != 1) {
	*col5_start_ptr = '0';
      }
      col6_start_ptr = skip_initial_spaces(col6_start_ptr);
      if (no_more_items_kns(col6_start_ptr)) {
	goto merge_fam_id_scan_ret_INVALID_FORMAT;
      }
      *col1_end_ptr = '\t';
      *col2_end_ptr = '\t';
      uii = col1_len + col2_len + 2;
      if (uii > max_person_id_len) {
	max_person_id_len = uii;
      }
      tot_len = uii + col3_len + col4_len + 4;
      uii = hashval(col1_start_ptr, col1_len, col2_start_ptr, col2_len);
      ll_pptr = &(htable[uii]);
      ll_ptr = *ll_pptr;
      uii = 1;
      if (is_dichot_pheno) {
	is_dichot_pheno = eval_affection(col6_start_ptr, -9, 2, 0);
      }
      if (scan_double(col6_start_ptr, &pheno)) {
	pheno = -9;
      }
      while (ll_ptr) {
	if (idmatch(ll_ptr->idstr, col1_start_ptr, col1_len + 1, col2_start_ptr, col2_len + 1)) {
	  uii = 0;
	  /*
	  // possibly for future: add parental ID/sex merge (not in PLINK 1.07)
	  if (merge_mode == 1) {
	    if (fabs(pheno - ll_ptr->pheno) > PHENO_EPSILON) {
	      ll_ptr->pheno = -9;
	    }
	  } else if (merge_mode == 2) {
	    if (ll_ptr->pheno == -9) {
	      ll_ptr->pheno = pheno;
	    }
	  } else if ((merge_mode == 5) || ((merge_mode == 3) && (pheno != -9))) {
	    ll_ptr->pheno = pheno;
	  }
	  */
	  break;
	}
        ll_pptr = &(ll_ptr->next);
	ll_ptr = *ll_pptr;
      }
      if (uii) {
	if (tot_len > max_person_full_len) {
	  max_person_full_len = tot_len;
	}
	ll_ptr = top_alloc_ll(&topsize, tot_len);
	ll_ptr->next = NULL;
	ll_ptr->pheno = pheno;
	ll_ptr->orig_order = orig_idx++;
	wptr = memcpyax(memcpyax(memcpyax(memcpyax(ll_ptr->idstr, col1_start_ptr, col1_len, '\t'), col2_start_ptr, col2_len, '\t'), col3_start_ptr, col3_len, '\t'), col4_start_ptr, col4_len, '\t');
	*wptr = *col5_start_ptr;
	wptr[1] = '\0';
	*ll_pptr = ll_ptr;
	tot_indiv_ct++;
      }
      cur_indiv_ct++;
    }
    if (!tbuf[MAXLINELEN - 1]) {
      if (!text_file) {
	sprintf(logbuf, "Error: .fam line too long in %s.\n", famname);
	goto merge_fam_id_scan_ret_INVALID_FORMAT;
      }
      ulii = 0;
      do {
	tbuf[MAXLINELEN - 1] = ' ';
	if (tbuf[MAXLINELEN - 2] == '\n') {
	  break;
	}
	ulii += MAXLINELEN - 1;
#ifndef __LP64__
	if (ulii > 0x7fffffff) {
	  sprintf(logbuf, "Error: .ped line too long (>= 2GB) in %s for 32-bit " PROG_NAME_CAPS ".\n", famname);
	  goto merge_fam_id_scan_ret_INVALID_FORMAT;
	}
#endif
        if (!fgets(tbuf, MAXLINELEN, infile)) {
	  goto merge_fam_id_scan_ret_READ_FAIL;
	}
      } while (!tbuf[MAXLINELEN - 1]);
      ulii += strlen(tbuf) + 1;
      if (ulii > (*ped_buflen_ptr)) {
	*ped_buflen_ptr = ulii;
      }
    }
  }
  if (!feof(infile)) {
    goto merge_fam_id_scan_ret_READ_FAIL;
  }
  if (!cur_indiv_ct) {
    sprintf(logbuf, "Error: No %s in %s.\n", g_species_plural, famname);
    goto merge_fam_id_scan_ret_INVALID_FORMAT;
  }
  *max_person_id_len_ptr = max_person_id_len;
  *max_person_full_len_ptr = max_person_full_len;
  *is_dichot_pheno_ptr = is_dichot_pheno;
  *topsize_ptr = topsize;
  *tot_indiv_ct_ptr = tot_indiv_ct;
  *cur_indiv_ct_ptr = cur_indiv_ct;
  *orig_idx_ptr = orig_idx;
  while (0) {
  merge_fam_id_scan_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  merge_fam_id_scan_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  merge_fam_id_scan_ret_INVALID_FORMAT:
    logprintb();
    retval = RET_INVALID_FORMAT;
  }
  fclose_cond(infile);
  return retval;
}

int32_t merge_bim_scan(char* bimname, uint32_t is_binary, uintptr_t* max_marker_id_len_ptr, Ll_entry2** htable2, uintptr_t* topsize_ptr, uint32_t* max_bim_linelen_ptr, uint64_t* tot_marker_ct_ptr, uint32_t* cur_marker_ct_ptr, uint32_t* position_warning_ct_ptr, Ll_str** non_biallelics_ptr, uint32_t allow_extra_chroms, Chrom_info* chrom_info_ptr) {
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t max_marker_id_len = *max_marker_id_len_ptr;
  uintptr_t topsize = *topsize_ptr;
  uint32_t max_bim_linelen = *max_bim_linelen_ptr;
  uint64_t tot_marker_ct = *tot_marker_ct_ptr;
  uint32_t cur_marker_ct = *cur_marker_ct_ptr;
  uint32_t position_warning_ct = *position_warning_ct_ptr;
  uint32_t loadbuf_size = MAXLINELEN;
  double cm = 0.0;
  FILE* infile = NULL;
  int32_t retval = 0;
  uint32_t alen1 = 1;
  uint32_t alen2 = 1;
  uint32_t cm_col;
  int64_t llxx;
  char* loadbuf;
  char* bufptr;
  char* bufptr2;
  char* bufptr3;
  uint32_t name_match;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  int32_t ii;
  int32_t jj;
  Ll_entry2** ll_pptr;
  Ll_entry2* ll_ptr;
  Ll_str* ll_string_new;
  char* aptr1;
  char* aptr2;
  char* new_aptr;
  uint32_t allele_ct;
  char* cur_alleles[2];
  if (fopen_checked(&infile, bimname, "r")) {
    goto merge_bim_scan_ret_OPEN_FAIL;
  }
  if (is_binary) {
    if (wkspace_left - topsize > 0x7fffff7f) {
      loadbuf_size = 0x3fffffc0;
    } else if (wkspace_left - topsize >= MAXLINELEN * 2) {
      loadbuf_size = ((wkspace_left - topsize) / 2) & (~(CACHELINE - 1));
    } else {
      goto merge_bim_scan_ret_NOMEM;
    }
  }
  loadbuf = (char*)wkspace_alloc(loadbuf_size);
  loadbuf[loadbuf_size - 1] = ' ';
  if (check_cm_col(infile, loadbuf, is_binary, loadbuf_size, &cm_col)) {
    goto merge_bim_scan_ret_INVALID_FORMAT_2;
  }
  do {
    if (!loadbuf[loadbuf_size - 1]) {
      if ((loadbuf_size == 0x3fffffc0) || ((!is_binary) && (loadbuf_size == MAXLINELEN))) {
	sprintf(logbuf, "Error: Pathologically long line in %s.\n", bimname);
	goto merge_bim_scan_ret_INVALID_FORMAT;
      } else {
	goto merge_bim_scan_ret_NOMEM;
      }
    }
    uii = strlen(loadbuf);
    if (uii >= max_bim_linelen) {
      max_bim_linelen = uii + 1;
    }
    bufptr = skip_initial_spaces(loadbuf);
    if (is_eoln_or_comment(*bufptr)) {
      continue;
    }
    ii = get_chrom_code(chrom_info_ptr, bufptr);
    if (ii == -1) {
      if (!allow_extra_chroms) {
	sprintf(logbuf, "Error: Invalid chromosome index in %s.\n", bimname);
	goto merge_bim_scan_ret_INVALID_FORMAT;
      }
      retval = resolve_or_add_chrom_name(chrom_info_ptr, bufptr, &ii);
      if (retval) {
	goto merge_bim_scan_ret_1;
      }
    }
    bufptr = next_item(bufptr);
    bufptr2 = item_endl(bufptr);
    uii = bufptr2 - bufptr;
    bufptr2 = skip_initial_spaces(bufptr2);
    if (no_more_items_kns(bufptr2)) {
      goto merge_bim_scan_ret_INVALID_FORMAT_2;
    }
    if (cm_col) {
      if (scan_double(bufptr2, &cm)) {
	cm = 0;
      }
      bufptr2 = next_item(bufptr2);
      if (no_more_items_kns(bufptr2)) {
	goto merge_bim_scan_ret_INVALID_FORMAT_2;
      }
    }
    if (!atoiz2(bufptr2, &jj)) {
      if (is_binary) {
	aptr1 = next_item(bufptr2);
	aptr2 = next_item(aptr1);
	if (no_more_items_kns(aptr2)) {
	  goto merge_bim_scan_ret_INVALID_FORMAT_2;
	}
	alen1 = strlen_se(aptr1);
	alen2 = strlen_se(aptr2);
	if ((alen1 == 1) && (*aptr1 == '0')) {
	  aptr1 = NULL;
	}
	if (aptr1 && (alen1 == alen2) && (!memcmp(aptr1, aptr2, alen1))) {
	  sprintf(logbuf, "Error: A1 and A2 alleles identical for a variant in %s.\n", bimname);
	  goto merge_bim_scan_ret_INVALID_FORMAT;
	}
	if ((alen2 == 1) && (*aptr2 == '0')) {
	  aptr2 = NULL;
	}
      } else {
	aptr1 = NULL;
	aptr2 = NULL;
      }
      llxx = (((int64_t)ii) << 32) + jj;
      ujj = hashval2(bufptr, uii);
      ll_pptr = &(htable2[ujj]);
      ll_ptr = *ll_pptr;
      name_match = 0;
      bufptr[uii++] = '\0';
      while (ll_ptr) {
	if (!memcmp(ll_ptr->idstr, bufptr, uii)) {
	  if (is_binary) {
	    bufptr2 = ll_ptr->allele[0];
	    allele_ct = 0;
	    if (bufptr2) {
	      cur_alleles[0] = bufptr2;
	      allele_ct = 1;
	    }
	    bufptr3 = ll_ptr->allele[1];
	    if (bufptr3) {
	      cur_alleles[allele_ct++] = bufptr3;
	    }
	    if (aptr2) {
	      for (ukk = 0; ukk < allele_ct; ukk++) {
		if ((!memcmp(aptr2, cur_alleles[ukk], alen2)) && (!cur_alleles[ukk][alen2])) {
		  break;
		}
	      }
	      if (ukk == allele_ct) {
		if (allele_ct == 2) {
		  ll_string_new = top_alloc_llstr(&topsize, uii);
		  if (!ll_string_new) {
		    goto merge_bim_scan_ret_NOMEM;
		  }
		  ll_string_new->next = *non_biallelics_ptr;
		  memcpy(ll_string_new->ss, bufptr, uii);
		  *non_biallelics_ptr = ll_string_new;
		} else {
		  if (allele_set(&new_aptr, aptr2, alen2)) {
		    goto merge_bim_scan_ret_NOMEM;
		  }
		  if (!ll_ptr->allele[1]) {
		    ll_ptr->allele[1] = new_aptr;
		  } else {
		    ll_ptr->allele[0] = new_aptr;
		  }
		  cur_alleles[allele_ct++] = new_aptr;
		}
	      }
	    }
	    if (aptr1) {
	      for (ukk = 0; ukk < allele_ct; ukk++) {
		if ((!memcmp(aptr1, cur_alleles[ukk], alen1)) && (!cur_alleles[ukk][alen1])) {
		  break;
		}
	      }
	      if (ukk == allele_ct) {
		if (allele_ct == 2) {
		  ll_string_new = top_alloc_llstr(&topsize, uii);
		  if (!ll_string_new) {
		    goto merge_bim_scan_ret_NOMEM;
		  }
		  ll_string_new->next = *non_biallelics_ptr;
		  memcpy(ll_string_new->ss, bufptr, uii);
		  *non_biallelics_ptr = ll_string_new;
		} else {
		  if (allele_set(&new_aptr, aptr1, alen1)) {
		    goto merge_bim_scan_ret_NOMEM;
		  }
		  if (!ll_ptr->allele[1]) {
		    ll_ptr->allele[1] = new_aptr;
		  } else {
		    ll_ptr->allele[0] = new_aptr;
		  }
		  cur_alleles[allele_ct++] = new_aptr;
		}
	      }
	    }
	  }
	  if (ll_ptr->pos != llxx) {
	    if ((ll_ptr->pos >> 32) == (llxx >> 32)) {
	      sprintf(logbuf, "Warning: Multiple positions seen for variant %s.\n", bufptr);
	      if (position_warning_ct < 3) {
		position_warning_ct++;
		if (position_warning_ct < 3) {
		  logprintb();
		} else {
		  logstr(logbuf);
		  printf("More position warnings: see log file.\n");
		}
	      } else {
		logstr(logbuf);
	      }
	    } else {
	      sprintf(logbuf, "Warning: Multiple chromosomes seen for variant %s.\n", bufptr);
	      logprintb();
	    }
	  }
	  name_match = 1;
	  break;
	}
        ll_pptr = &(ll_ptr->next);
	ll_ptr = *ll_pptr;
      }
      if (!name_match) {
        if (uii > max_marker_id_len) {
	  max_marker_id_len = uii;
	}
	ll_ptr = top_alloc_ll2(&topsize, uii);
	if (!ll_ptr) {
	  goto merge_bim_scan_ret_NOMEM;
	}
	ll_ptr->next = NULL;
	ll_ptr->pos = llxx;
	ll_ptr->cm = cm;
	if (aptr1) {
	  if (allele_set(&(ll_ptr->allele[0]), aptr1, alen1)) {
	    goto merge_bim_scan_ret_NOMEM;
	  }
	} else {
	  ll_ptr->allele[0] = NULL;
	}
	if (aptr2) {
	  if (allele_set(&(ll_ptr->allele[1]), aptr2, alen2)) {
	    goto merge_bim_scan_ret_NOMEM;
	  }
	} else {
	  ll_ptr->allele[1] = NULL;
	}
	memcpy(ll_ptr->idstr, bufptr, uii);
	*ll_pptr = ll_ptr;
	tot_marker_ct++;
      }
      cur_marker_ct++;
    } else if (!atoi(bufptr2)) {
      sprintf(logbuf, "Error: Invalid base-pair position in %s.\n", bimname);
      goto merge_bim_scan_ret_INVALID_FORMAT;
    }
  } while (fgets(loadbuf, loadbuf_size, infile));
  if (!feof(infile)) {
    goto merge_bim_scan_ret_READ_FAIL;
  }
  if (!cur_marker_ct) {
    sprintf(logbuf, "Error: No variants in %s.\n", bimname);
    goto merge_bim_scan_ret_INVALID_FORMAT;
  }
  *max_marker_id_len_ptr = max_marker_id_len;
  *topsize_ptr = topsize;
  *max_bim_linelen_ptr = max_bim_linelen;
  *tot_marker_ct_ptr = tot_marker_ct;
  *cur_marker_ct_ptr = cur_marker_ct;
  *position_warning_ct_ptr = position_warning_ct;
  
  while (0) {
  merge_bim_scan_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  merge_bim_scan_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  merge_bim_scan_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  merge_bim_scan_ret_INVALID_FORMAT_2:
    sprintf(logbuf, "Error: %s is not a properly formatted .map/.bim file.\n", bimname);
  merge_bim_scan_ret_INVALID_FORMAT:
    logprintb();
    retval = RET_INVALID_FORMAT;
  }
 merge_bim_scan_ret_1:
  fclose_cond(infile);
  wkspace_reset(wkspace_mark);
  return retval;
}

int32_t report_non_biallelics(char* outname, char* outname_end, Ll_str* non_biallelics) {
  FILE* outfile = NULL;
  int32_t retval = 0;
  uintptr_t nbmarker_ct_dup = 0;
  uintptr_t nbmarker_ct = 1;
  uintptr_t max_nbmarker_id_len = 0;
  Ll_str* cur_ptr = non_biallelics;
  uintptr_t cur_slen;
  char* id_arr;
  char* id_arr_ptr;
  char* id_arr_ptr_old;
  char* id_arr_end;
  do {
    cur_slen = strlen(cur_ptr->ss);
    if (cur_slen >= max_nbmarker_id_len) {
      max_nbmarker_id_len = cur_slen + 1;
    }
    nbmarker_ct_dup++;
    cur_ptr = cur_ptr->next;
  } while (cur_ptr);
  if (wkspace_alloc_c_checked(&id_arr, nbmarker_ct_dup * max_nbmarker_id_len)) {
    goto report_non_biallelics_ret_NOMEM;
  }
  cur_ptr = non_biallelics;
  id_arr_ptr = id_arr;
  do {
    strcpy(id_arr_ptr, cur_ptr->ss);
    cur_ptr = cur_ptr->next;
    id_arr_ptr = &(id_arr_ptr[max_nbmarker_id_len]);
  } while (cur_ptr);
  qsort(id_arr, nbmarker_ct_dup, max_nbmarker_id_len, strcmp_casted);
  memcpy(outname_end, ".missnp", 8);
  if (fopen_checked(&outfile, outname, "w")) {
    goto report_non_biallelics_ret_OPEN_FAIL;
  }
  id_arr_ptr = id_arr;
  if (fputs_checked(id_arr_ptr, outfile)) {
    goto report_non_biallelics_ret_WRITE_FAIL;
  }
  putc('\n', outfile);
  id_arr_end = &(id_arr[(nbmarker_ct_dup - 1) * max_nbmarker_id_len]);
  while (id_arr_ptr != id_arr_end) {
    id_arr_ptr_old = id_arr_ptr;
    id_arr_ptr = &(id_arr_ptr[max_nbmarker_id_len]);
    if (strcmp(id_arr_ptr, id_arr_ptr_old)) {
      nbmarker_ct++;
      fputs(id_arr_ptr, outfile);
      if (putc_checked('\n', outfile)) {
        goto report_non_biallelics_ret_WRITE_FAIL;
      }
    }
  }
  if (fclose_null(&outfile)) {
    goto report_non_biallelics_ret_WRITE_FAIL;
  }
  sprintf(logbuf, "Error: %" PRIuPTR " variant%s with 3+ alleles present.  Use --flip with\n%s for now.\n", nbmarker_ct, (nbmarker_ct == 1)? "" : "s", outname);
  logprintb();
  while (0) {
  report_non_biallelics_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  report_non_biallelics_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  report_non_biallelics_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  fclose_cond(outfile);
  return retval;
}

void merge_alleles_update_str(char* marker_allele_ptr, char** allele_ptrs, uint32_t* distinct_allele_ct_ptr) {
  uint32_t distinct_allele_ct = *distinct_allele_ct_ptr;
  uint32_t allele_idx = 0;
  if ((marker_allele_ptr == g_missing_geno_ptr) || (distinct_allele_ct == 3)) {
    return;
  }
  while (allele_idx < distinct_allele_ct) {
    if (!strcmp(marker_allele_ptr, allele_ptrs[allele_idx])) {
      return;
    }
    allele_idx++;
  }
  *distinct_allele_ct_ptr = distinct_allele_ct + 1;
  if (distinct_allele_ct == 2) {
    return;
  }
  allele_ptrs[distinct_allele_ct] = marker_allele_ptr;
}

uint32_t merge_alleles(char** marker_allele_ptrs, uint32_t marker_uidx, uint32_t marker_uidx2) {
  uint32_t distinct_allele_ct = 0;
  char* allele_ptrs[2];
  allele_ptrs[0] = NULL;
  allele_ptrs[1] = NULL;
  // reverse order so --keep-allele-order works
  merge_alleles_update_str(marker_allele_ptrs[2 * marker_uidx + 1], allele_ptrs, &distinct_allele_ct);
  merge_alleles_update_str(marker_allele_ptrs[2 * marker_uidx], allele_ptrs, &distinct_allele_ct);
  merge_alleles_update_str(marker_allele_ptrs[2 * marker_uidx2 + 1], allele_ptrs, &distinct_allele_ct);
  merge_alleles_update_str(marker_allele_ptrs[2 * marker_uidx2], allele_ptrs, &distinct_allele_ct);
  if (distinct_allele_ct > 2) {
    return 1;
  }
  if (allele_ptrs[0]) {
    // todo: test if this write is inefficient enough that we want to guard it
    // with an if statement
    marker_allele_ptrs[2 * marker_uidx + 1] = allele_ptrs[0];
    if (allele_ptrs[1]) {
      marker_allele_ptrs[2 * marker_uidx] = allele_ptrs[1];
    }
  }
  return 0;
}

static inline uint32_t merge_post_msort_update_maps(char* marker_ids, uintptr_t max_marker_id_len, uint32_t* marker_map, uint32_t* pos_buf, int64_t* ll_buf, uint32_t* chrom_start, uint32_t* chrom_id, uint32_t chrom_ct, uint32_t* dedup_marker_ct_ptr, uint32_t merge_equal_pos, char** marker_allele_ptrs, Chrom_info* chrom_info_ptr) {
  // Input: ll_buf is a sequence of sorted arrays (one per chromosome) with
  // base-pair positions in high 32 bits, and pre-sort indices in low 32 bits.
  // Chromosome boundaries are stored in chrom_start[].
  // Pre-sort indices refer to marker ID ASCII ordering, which is how lookup
  // will be performed mid-merge.  There may be duplicate positions, and
  // markers that don't pass the chromosome filter.

  // Result: Duplicates have been collapsed, with chrom_start[] updated.
  // pos_buf contains sorted base-pair positions,
  // post-chromosome-filtering-and-duplicate-removal.
  // marker_map[n] is the post-filtering position in all other arrays of marker
  // ID n by ASCII ordering.
  uintptr_t* chrom_mask = chrom_info_ptr->chrom_mask;
  uint32_t read_pos = 0;
  uint32_t write_pos = 0; // may be lower than read_pos due to dups
  uint32_t chrom_idx;
  uint32_t chrom_read_end_idx;
  int64_t llxx;
  uint32_t prev_bp;
  uint32_t cur_bp;
  uint32_t presort_idx;
  for (chrom_idx = 0; chrom_idx < chrom_ct; chrom_idx++) {
    if (!IS_SET(chrom_mask, chrom_id[chrom_idx])) {
      read_pos = chrom_start[chrom_idx + 1];
      chrom_start[chrom_idx + 1] = write_pos;
      continue;
    }
    chrom_read_end_idx = chrom_start[chrom_idx + 1];
    // ll_buf has base-pair positions in high 32 bits, and pre-sort indices in
    // low 32 bits.
    llxx = ll_buf[read_pos++];
    prev_bp = (uint32_t)(llxx >> 32);
    pos_buf[write_pos] = prev_bp;
    marker_map[(uint32_t)llxx] = write_pos++;
    for (; read_pos < chrom_read_end_idx; read_pos++) {
      llxx = ll_buf[read_pos];
      presort_idx = (uint32_t)llxx;
      cur_bp = (uint32_t)(llxx >> 32);
      if (prev_bp == cur_bp) {
	if (merge_equal_pos && merge_alleles(marker_allele_ptrs, ((uint32_t)ll_buf[read_pos - 1]), presort_idx)) {
	  sprintf(logbuf, "Error: --merge-equal-pos failure.  Variants %s and %s\nhave the same position, but do not share the same alleles.\n", &(marker_ids[max_marker_id_len * presort_idx]), &(marker_ids[max_marker_id_len * ((uint32_t)ll_buf[read_pos - 1])]));
	  logprintb();
	  return 1;
	}
	sprintf(logbuf, "Warning: Variants %s and %s have the same position.\n", &(marker_ids[max_marker_id_len * presort_idx]), &(marker_ids[max_marker_id_len * ((uint32_t)ll_buf[read_pos - 1])]));
	logprintb();
	if (merge_equal_pos) {
	  marker_map[presort_idx] = write_pos - 1;
	  continue;
	}
      } else {
	prev_bp = cur_bp;
      }
      marker_map[presort_idx] = write_pos;
      pos_buf[write_pos++] = cur_bp;
    }
    read_pos = chrom_start[chrom_idx + 1];
    chrom_start[chrom_idx + 1] = write_pos;
  }
  *dedup_marker_ct_ptr = write_pos;
  return 0;
}

static inline int32_t merge_must_track_write(int32_t mm) {
  // modes 6 and 7 can be sped up with early pruning of nonoverlapping
  // markers, but not worth complicating code for this.
  return (mm == 1) || (mm > 5) || (mm == 4);
}

static inline int32_t merge_first_mode(int32_t mm, uint32_t merge_equal_pos) {
  if (merge_equal_pos) {
    return (mm > 5)? 4 : mm;
  } else {
    return merge_must_track_write(mm)? 4 : 5;
  }
}

int32_t merge_diff_print(FILE* outfile, char* idbuf, char* marker_id, char* person_id, unsigned char newval, unsigned char oldval, char** marker_allele_ptrs) {
  char* bufptr = item_endnn(person_id);
  uint32_t slen = strlen_se(marker_id);
  char* ma1p[4];
  char* ma2p[4];
  char wbuf[8];
  uint32_t slen2;
  memcpyx(idbuf, marker_id, slen, 0);
  fprintf(outfile, "%20s", idbuf);
  slen = (bufptr++) - person_id;
  memcpyx(idbuf, person_id, slen, 0);
  fprintf(outfile, " %20s %20s", idbuf, bufptr);
  ma1p[0] = marker_allele_ptrs[0];
  ma1p[1] = (char*)g_missing_geno_ptr;
  ma1p[2] = marker_allele_ptrs[0];
  ma1p[3] = marker_allele_ptrs[1];
  ma2p[0] = marker_allele_ptrs[0];
  ma2p[1] = ma1p[1];
  ma2p[2] = marker_allele_ptrs[1];
  ma2p[3] = marker_allele_ptrs[1];
  slen = strlen(ma1p[newval]);
  slen2 = strlen(ma2p[newval]);
  if (slen + slen2 > 6) {
    fprintf(outfile, " %s/%s", ma1p[newval], ma2p[newval]);
  } else {
    memcpy(wbuf, ma1p[newval], slen);
    wbuf[slen] = '/';
    memcpy(&(wbuf[slen + 1]), ma2p[newval], slen2 + 1);
    fprintf(outfile, " %8s", wbuf);
  }
  slen = strlen(ma1p[oldval]);
  slen2 = strlen(ma2p[oldval]);
  if (slen + slen2 > 6) {
    fprintf(outfile, " %s/%s \n", ma1p[oldval], ma2p[oldval]);
  } else {
    memcpy(wbuf, ma1p[oldval], slen);
    wbuf[slen] = '/';
    memcpy(&(wbuf[slen + 1]), ma2p[oldval], slen2 + 1);
    fprintf(outfile, " %8s \n", wbuf);
  }
  return ferror(outfile);
}

int32_t merge_main(char* bedname, char* bimname, char* famname, char* bim_loadbuf, uint32_t max_bim_linelen, uint32_t tot_indiv_ct, uint32_t tot_marker_ct, uint32_t dedup_marker_ct, uint32_t start_marker_idx, uint32_t marker_window_size, char** marker_allele_ptrs, char* marker_ids, uintptr_t max_marker_id_len, char* person_ids, uintptr_t max_person_id_len, uint32_t merge_nsort, uint32_t* indiv_nsmap, uint32_t* flex_map, uint32_t* marker_map, char* idbuf, unsigned char* readbuf, unsigned char* writebuf, uint32_t merge_mode, uintptr_t* markbuf, FILE* outfile, uint64_t* diff_total_overlap_ptr, uint64_t* diff_not_both_genotyped_ptr, uint64_t* diff_discordant_ptr, uint32_t ped_buflen, unsigned char* bmap_raw) {
  // flex_map maps individuals for binary filesets, and markers for text
  // filesets.
  uint32_t is_binary = famname? 1 : 0;
  FILE* bedfile = NULL;
  FILE* infile2 = NULL;
  int32_t retval = 0;
  // bugfix: there was a potential integer overflow back when these were
  // uint32_t
  uintptr_t tot_indiv_ct4 = (tot_indiv_ct + 3) / 4;
  uintptr_t tot_indiv_ctl = (tot_indiv_ct + (BITCT - 1)) / BITCT;
  uint32_t end_marker_idx = start_marker_idx + marker_window_size;
  uint32_t marker_in_idx = 0xffffffffU; // overflow to zero on first add
  uint32_t last_marker_in_idx = 0xfffffffeU;
  uint32_t cur_indiv_ct = 0;
  uintptr_t* mbufptr = NULL; // merge mode 1, 4, 6, 7
  const char* missing_geno_ptr = g_missing_geno_ptr;
  uint64_t diff_total_overlap = 0;
  uint64_t diff_not_both_genotyped = 0;
  uint64_t diff_discordant = 0;
  uint32_t cur_indiv_ct4 = 0;
  uint32_t cur_indiv_ctd4 = 0;
  uint32_t is_ped_compound = 1; // 0 = no, 1 = unresolved, 2 = yes
  uint32_t alen1 = 1;
  uint32_t alen2 = 1;
  uintptr_t uljj = 0;
  uintptr_t* mbufptr2;
  unsigned char* bmap;
  uint32_t cm_col;
  char* aptr1;
  char* aptr2;
  char* bufptr;
  char* bufptr2;
  char* bufptr3;
  char* bufptr4;
  char* bufptr5;
  unsigned char* rbufptr;
  unsigned char* wbufptr;
  unsigned char* wbufptr2;
  uintptr_t ulii;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  uint32_t unn;
  uintptr_t indiv_idx;
  uint32_t marker_out_idx;
  unsigned char ucc;
  unsigned char ucc2;
  unsigned char ucc3;
  unsigned char ucc4;
  char cc;
  int32_t ii;
  if (is_binary) {
    if (fopen_checked(&infile2, famname, "r")) {
      goto merge_main_ret_OPEN_FAIL;
    }
    while (fgets(tbuf, MAXLINELEN, infile2)) {
      bufptr = skip_initial_spaces(tbuf);
      if (is_eoln_kns(*bufptr)) {
	continue;
      }
      bufptr2 = item_endnn(bufptr);
      bufptr3 = skip_initial_spaces(bufptr2);
      bufptr4 = item_endnn(bufptr3);
      uii = (bufptr2 - bufptr);
      ujj = (bufptr4 - bufptr3);
      memcpyx(memcpyax(idbuf, bufptr, uii, '\t'), bufptr3, ujj, 0);
      if (merge_nsort) {
        ii = bsearch_str_natural(idbuf, person_ids, max_person_id_len, 0, tot_indiv_ct - 1);
      } else {
	ii = bsearch_str(idbuf, person_ids, max_person_id_len, 0, tot_indiv_ct - 1);
	if (indiv_nsmap && (ii != -1)) {
	  ii = indiv_nsmap[(uint32_t)ii];
	}
      }
      if (ii == -1) {
	// previously validated, so give read failure error code instead of
	// invalid format
	goto merge_main_ret_READ_FAIL;
      }
      flex_map[cur_indiv_ct++] = ii;
    }
    if (!feof(infile2)) {
      goto merge_main_ret_READ_FAIL;
    }
    fclose_null(&infile2);
    cur_indiv_ct4 = (cur_indiv_ct + 3) / 4;
    cur_indiv_ctd4 = cur_indiv_ct / 4;
  } else {
    bim_loadbuf = tbuf;
    max_bim_linelen = MAXLINELEN;
  }
  if (fopen_checked(&infile2, bimname, "r")) {
    goto merge_main_ret_OPEN_FAIL;
  }
  if (check_cm_col(infile2, bim_loadbuf, is_binary, max_bim_linelen, &cm_col)) {
    goto merge_main_ret_READ_FAIL;
  }
  if (fopen_checked(&bedfile, bedname, is_binary? "rb" : "r")) {
    goto merge_main_ret_OPEN_FAIL;
  }
  do {
    bufptr = skip_initial_spaces(bim_loadbuf);
    if (is_eoln_or_comment(*bufptr)) {
      continue;
    }
    ++marker_in_idx;
    bufptr = next_item(bufptr);
    bufptr2 = next_item_mult(bufptr, 1 + cm_col);
    if (!bufptr2) {
      goto merge_main_ret_READ_FAIL;
    }
    if (*bufptr2 == '-') {
      if (!is_binary) {
	flex_map[marker_in_idx] = 0xffffffffU;
      }
      continue;
    }
    bufptr3 = item_endnn(bufptr);
    uii = (bufptr3 - bufptr);
    memcpyx(idbuf, bufptr, uii, 0);
    ii = bsearch_str(idbuf, marker_ids, max_marker_id_len, 0, tot_marker_ct - 1);
    if (ii == -1) {
      goto merge_main_ret_READ_FAIL;
    }
    marker_out_idx = marker_map[(uint32_t)ii];
    if ((marker_out_idx < start_marker_idx) || (marker_out_idx >= end_marker_idx)) {
      if (!is_binary) {
	flex_map[marker_in_idx] = 0xffffffffU;
      }
      continue;
    }
    if (is_binary) {
      if (marker_in_idx != last_marker_in_idx + 1) {
	if (fseeko(bedfile, 3 + ((uint64_t)marker_in_idx) * cur_indiv_ct4, SEEK_SET)) {
	  goto merge_main_ret_READ_FAIL;
	}
      }
      bufptr2 = next_item(bufptr2);
      bufptr3 = next_item(bufptr2);
      if (no_more_items_kns(bufptr3)) {
	goto merge_main_ret_READ_FAIL;
      }
      alen1 = strlen_se(bufptr2);
      bufptr2[alen1] = '\0';
      alen2 = strlen_se(bufptr3);
      bufptr3[alen2] = '\0';
      bufptr4 = marker_allele_ptrs[((uint32_t)ii) * 2];
      bufptr5 = marker_allele_ptrs[((uint32_t)ii) * 2 + 1];
      if ((((*bufptr2 != '0') || (alen1 != 1)) && (!strcmp(bufptr2, bufptr5)))  || (((*bufptr3 != '0') || (alen2 != 1)) && (!strcmp(bufptr3, bufptr4)))) {
	bmap = &(bmap_raw[256]);
      } else {
	bmap = bmap_raw;
      }

      last_marker_in_idx = marker_in_idx;
      if (fread(readbuf, 1, cur_indiv_ct4, bedfile) < cur_indiv_ct4) {
	goto merge_main_ret_READ_FAIL;
      }
      rbufptr = readbuf;
      wbufptr = &(writebuf[(marker_out_idx - start_marker_idx) * tot_indiv_ct4]);
      if (merge_must_track_write(merge_mode)) {
	mbufptr = &(markbuf[(marker_out_idx - start_marker_idx) * tot_indiv_ctl]);
      }
      switch (merge_mode) {
      case 1: // difference -> missing
	indiv_idx = 0;
	for (uii = 0; uii < cur_indiv_ctd4; uii++) {
	  ucc = bmap[*rbufptr++];
	  do {
	    ujj = flex_map[indiv_idx++];
	    ukk = ujj / BITCT;
	    ulii = ONELU << (ujj % BITCT);
	    wbufptr2 = &(wbufptr[ujj / 4]);
	    umm = (ujj % 4) * 2;
	    unn = 3U << umm;
	    if (mbufptr[ukk] & ulii) {
	      ucc2 = *wbufptr2;
	      if ((ucc2 ^ ((ucc & 3) << umm)) & unn) {
		*wbufptr2 = (ucc2 & (~unn)) | (1U << umm);
	      }
	    } else {
	      mbufptr[ukk] |= ulii;
	      *wbufptr2 = ((*wbufptr2) & (~unn)) | ((ucc & 3) << umm);
	    }
	    ucc >>= 2;
	  } while (indiv_idx & 3);
	}
	ucc = bmap[*rbufptr];
	while (indiv_idx < cur_indiv_ct) {
	  ujj = flex_map[indiv_idx++];
	  ukk = ujj / BITCT;
	  ulii = ONELU << (ujj % BITCT);
	  wbufptr2 = &(wbufptr[ujj / 4]);
	  umm = (ujj % 4) * 2;
	  unn = 3U << umm;
	  if (mbufptr[ukk] & ulii) {
	    ucc2 = *wbufptr2;
	    if ((ucc2 ^ ((ucc & 3) << umm)) & unn) {
	      *wbufptr2 = (ucc2 & (~unn)) | (1U << umm);
	    }
	  } else {
	    mbufptr[ukk] |= ulii;
	    *wbufptr2 = ((*wbufptr2) & (~unn)) | ((ucc & 3) << umm);
	  }
	  ucc >>= 2;
	}
	break;
      case 2: // only overwrite originally missing
	indiv_idx = 0;
	for (uii = 0; uii < cur_indiv_ctd4; uii++) {
	  ucc = bmap[*rbufptr++];
	  do {
	    ujj = flex_map[indiv_idx++];
	    ukk = (ujj % 4) * 2;
	    wbufptr2 = &(wbufptr[ujj / 4]);
	    ucc2 = *wbufptr2;
	    if (((ucc2 >> ukk) & 3) == 1) {
	      *wbufptr2 = (ucc2 & (~(3U << ukk))) | ((ucc & 3) << ukk);
	    }
	    ucc >>= 2;
	  } while (indiv_idx & 3);
	}
	ucc = bmap[*rbufptr++];
	while (indiv_idx < cur_indiv_ct) {
	  ujj = flex_map[indiv_idx++];
	  ukk = (ujj % 4) * 2;
	  wbufptr2 = &(wbufptr[ujj / 4]);
	  ucc2 = *wbufptr2;
	  if (((ucc2 >> ukk) & 3) == 1) {
	    *wbufptr2 = (ucc2 & (~(3U << ukk))) | ((ucc & 3) << ukk);
	  }
	  ucc >>= 2;
	}
	break;
      case 3: // only overwrite if nonmissing in new file
	indiv_idx = 0;
	for (uii = 0; uii < cur_indiv_ctd4; uii++) {
	  ucc = bmap[*rbufptr++];
	  do {
	    ucc2 = ucc & 3;
	    if (ucc2 != 1) {
	      ujj = flex_map[indiv_idx];
	      ukk = (ujj % 4) * 2;
	      wbufptr2 = &(wbufptr[ujj / 4]);
	      *wbufptr2 = ((*wbufptr2) & (~(3U << ukk))) | (ucc2 << ukk);
	    }
	    ucc >>= 2;
	  } while ((++indiv_idx) & 3);
	}
	ucc = bmap[*rbufptr++];
	while (indiv_idx < cur_indiv_ct) {
	  ucc2 = ucc & 3;
	  if (ucc2 != 1) {
	    ujj = flex_map[indiv_idx];
	    ukk = (ujj % 4) * 2;
	    wbufptr2 = &(wbufptr[ujj / 4]);
	    *wbufptr2 = ((*wbufptr2) & (~(3U << ukk))) | (ucc2 << ukk);
	  }
	  ucc >>= 2;
	  indiv_idx++;
	}
	break;
      case 4: // never overwrite
	indiv_idx = 0;
	for (uii = 0; uii < cur_indiv_ctd4; uii++) {
	  ucc = bmap[*rbufptr++];
	  do {
	    ujj = flex_map[indiv_idx++];
	    ukk = ujj / BITCT;
	    ulii = ONELU << (ujj % BITCT);
	    if (!(mbufptr[ukk] & ulii)) {
	      mbufptr[ukk] |= ulii;
	      wbufptr2 = &(wbufptr[ujj / 4]);
	      ukk = (ujj % 4) * 2;
	      *wbufptr2 = ((*wbufptr2) & (~(3U << ukk))) | ((ucc & 3) << ukk);
	    }
	    ucc >>= 2;
	  } while (indiv_idx & 3);
	}
	ucc = bmap[*rbufptr];
	while (indiv_idx < cur_indiv_ct) {
	  ujj = flex_map[indiv_idx++];
	  ukk = ujj / BITCT;
	  ulii = ONELU << (ujj % BITCT);
	  if (!(mbufptr[ukk] & ulii)) {
	    mbufptr[ukk] |= ulii;
	    wbufptr2 = &(wbufptr[ujj / 4]);
	    ukk = (ujj % 4) * 2;
	    *wbufptr2 = ((*wbufptr2) & (~(3U << ukk))) | ((ucc & 3) << ukk);
	  }
	  ucc >>= 2;
	}
	break;
      case 5: // always overwrite
	indiv_idx = 0;
	for (uii = 0; uii < cur_indiv_ctd4; uii++) {
	  ucc = bmap[*rbufptr++];
	  do {
	    ujj = flex_map[indiv_idx++];
	    ukk = (ujj % 4) * 2;
	    wbufptr2 = &(wbufptr[ujj / 4]);
	    *wbufptr2 = ((*wbufptr2) & (~(3U << ukk))) | ((ucc & 3) << ukk);
	    ucc >>= 2;
	  } while (indiv_idx & 3);
	}
	ucc = bmap[*rbufptr];
	while (indiv_idx < cur_indiv_ct) {
	  ujj = flex_map[indiv_idx++];
	  ukk = (ujj % 4) * 2;
	  wbufptr2 = &(wbufptr[ujj / 4]);
	  *wbufptr2 = ((*wbufptr2) & (~(3U << ukk))) | ((ucc & 3) << ukk);
	  ucc >>= 2;
	}
	break;
      case 6: // report all mismatches
	indiv_idx = 0;
	for (uii = 0; uii < cur_indiv_ctd4; uii++) {
	  ucc = bmap[*rbufptr++];
	  do {
	    ujj = flex_map[indiv_idx++];
	    if (mbufptr[ujj / BITCT] & (ONELU << (ujj % BITCT))) {
	      // would prefer to do this by multiplying indiv overlap with
	      // marker overlap, but the same-position automerge screws
	      // with that
	      diff_total_overlap++;
	      ukk = (ujj % 4) * 2;
	      ucc2 = ucc & 3;
	      ucc3 = (wbufptr[ujj / 4] >> ukk) & 3;
	      umm = ((ucc2 == 1) || (ucc3 == 1));
	      if (umm) {
		diff_not_both_genotyped++;
	      }
	      if (ucc2 != ucc3) {
		if (!umm) {
		  diff_discordant++;
		}
		if (merge_diff_print(outfile, idbuf, bufptr, &(person_ids[ujj * max_person_id_len]), ucc2, ucc3, &(marker_allele_ptrs[((uint32_t)ii) * 2]))) {
		  goto merge_main_ret_WRITE_FAIL;
		}
	      }
	    }
	    ucc >>= 2;
	  } while (indiv_idx & 3);
	}
	ucc = bmap[*rbufptr];
	while (indiv_idx < cur_indiv_ct) {
	  ujj = flex_map[indiv_idx++];
	  if (mbufptr[ujj / BITCT] & (ONELU << (ujj % BITCT))) {
	    diff_total_overlap++;
	    ukk = (ujj % 4) * 2;
	    ucc2 = ucc & 3;
	    ucc3 = (wbufptr[ujj / 4] >> ukk) & 3;
	    umm = ((ucc2 == 1) || (ucc3 == 1));
	    if (umm) {
	      diff_not_both_genotyped++;
	    }
	    if (ucc2 != ucc3) {
	      if (!umm) {
		diff_discordant++;
	      }
	      if (merge_diff_print(outfile, idbuf, bufptr, &(person_ids[ujj * max_person_id_len]), ucc2, ucc3, &(marker_allele_ptrs[((uint32_t)ii) * 2]))) {
		goto merge_main_ret_WRITE_FAIL;
	      }
	    }
	  }
	  ucc >>= 2;
	}
	break;
      case 7: // report nonmissing mismatches
	indiv_idx = 0;
	for (uii = 0; uii < cur_indiv_ctd4; uii++) {
	  ucc = bmap[*rbufptr++];
	  do {
	    ujj = flex_map[indiv_idx++];
	    if (mbufptr[ujj / BITCT] & (ONELU << (ujj % BITCT))) {
	      diff_total_overlap++;
	      ukk = (ujj % 4) * 2;
	      ucc2 = ucc & 3;
	      ucc3 = (wbufptr[ujj / 4] >> ukk) & 3;
	      if ((ucc2 == 1) || (ucc3 == 1)) {
		diff_not_both_genotyped++;
	      } else if (ucc2 != ucc3) {
		diff_discordant++;
		if (merge_diff_print(outfile, idbuf, bufptr, &(person_ids[ujj * max_person_id_len]), ucc2, ucc3, &(marker_allele_ptrs[((uint32_t)ii) * 2]))) {
		  goto merge_main_ret_WRITE_FAIL;
		}
	      }
	    }
	    ucc >>= 2;
	  } while (indiv_idx & 3);
	}
	ucc = bmap[*rbufptr];
	while (indiv_idx < cur_indiv_ct) {
	  ujj = flex_map[indiv_idx++];
	  if (mbufptr[ujj / BITCT] & (ONELU << (ujj % BITCT))) {
	    diff_total_overlap++;
	    ukk = (ujj % 4) * 2;
	    ucc2 = ucc & 3;
	    ucc3 = (wbufptr[ujj / 4] >> ukk) & 3;
	    if ((ucc2 == 1) || (ucc3 == 1)) {
	      diff_not_both_genotyped++;
	    } else if (ucc2 != ucc3) {
	      diff_discordant++;
	      if (merge_diff_print(outfile, idbuf, bufptr, &(person_ids[ujj * max_person_id_len]), ucc2, ucc3, &(marker_allele_ptrs[((uint32_t)ii) * 2]))) {
		goto merge_main_ret_WRITE_FAIL;
	      }
	    }
	  }
	  ucc >>= 2;
	}
	break;
      }
    } else {
      flex_map[marker_in_idx] = ii;
    }
  } while (fgets(bim_loadbuf, max_bim_linelen, infile2));
  if (!feof(infile2)) {
    goto merge_main_ret_READ_FAIL;
  }
  if (!is_binary) {
    last_marker_in_idx = marker_in_idx + 1; // total count
    while (fgets((char*)readbuf, ped_buflen, bedfile)) {
      bufptr = skip_initial_spaces((char*)readbuf);
      cc = *bufptr;
      if (is_eoln_or_comment(cc)) {
	continue;
      }
      bufptr2 = item_endnn(bufptr);
      uii = (bufptr2 - bufptr);
      bufptr3 = skip_initial_spaces(bufptr2);
      bufptr2 = item_endnn(bufptr3);
      ujj = (bufptr2 - bufptr3);
      memcpyx(memcpyax(idbuf, bufptr, uii, '\t'), bufptr3, ujj, 0);
      if (merge_nsort) {
        ii = bsearch_str_natural(idbuf, person_ids, max_person_id_len, 0, tot_indiv_ct - 1);
      } else {
	ii = bsearch_str(idbuf, person_ids, max_person_id_len, 0, tot_indiv_ct - 1);
	if (indiv_nsmap && (ii != -1)) {
	  ii = indiv_nsmap[(uint32_t)ii];
	}
      }
      if (ii == -1) {
	goto merge_main_ret_READ_FAIL;
      }
      bufptr3 = next_item_mult(skip_initial_spaces(bufptr2), 4);
      if (!bufptr3) {
	goto merge_main_ret_READ_FAIL;
      }
      if (is_ped_compound == 1) {
	bufptr4 = bufptr3;
	for (marker_in_idx = 0; marker_in_idx < 2 * last_marker_in_idx; marker_in_idx++) {
          if (is_eoln_kns(*bufptr4)) {
	    is_ped_compound = 2;
	    break;
	  }
	  bufptr4 = skip_initial_spaces(item_endnn(bufptr4));
	}
	if (is_ped_compound == 1) {
	  is_ped_compound = 0;
	}
      }
      wbufptr = &(writebuf[((uint32_t)ii) / 4]);
      ujj = (((uint32_t)ii) % 4) * 2;
      ucc = ~(3U << ujj);
      ucc4 = (1U << ujj);
      if (merge_must_track_write(merge_mode)) {
	mbufptr = &(markbuf[((uint32_t)ii) / BITCT]);
	uljj = ONELU << (((uint32_t)ii) % BITCT);
      }
      for (marker_in_idx = 0; marker_in_idx < last_marker_in_idx; marker_in_idx++) {
	if (!is_ped_compound) {
	  if (is_eoln_kns(*bufptr3)) {
	    goto merge_main_ret_INVALID_FORMAT_2;
	  }
	  aptr1 = bufptr3;
	  bufptr3 = item_endnn(bufptr3);
	  alen1 = (uintptr_t)(bufptr3 - aptr1);
	  bufptr3 = skip_initial_spaces(bufptr3);
	  aptr1[alen1] = '\0';
	  if (is_eoln_kns(*bufptr3)) {
	    goto merge_main_ret_INVALID_FORMAT_2;
	  }
	  aptr2 = bufptr3;
	  bufptr3 = item_endnn(bufptr3);
	  alen2 = (uintptr_t)(bufptr3 - aptr2);
	  bufptr3 = skip_initial_spaces(bufptr3);
	  aptr2[alen2] = '\0';
	} else {
	  cc = *bufptr3;
	  if (is_eoln_kns(cc)) {
	    goto merge_main_ret_INVALID_FORMAT_2;
	  }
          aptr1 = (char*)(&(g_one_char_strs[((unsigned char)cc) * 2]));
          bufptr3 = skip_initial_spaces(&(bufptr3[1]));
          cc = *bufptr3;
          if (is_eoln_kns(cc)) {
	    goto merge_main_ret_INVALID_FORMAT_2;
	  }
	  aptr2 = (char*)(&(g_one_char_strs[((unsigned char)cc) * 2]));
	  bufptr3 = skip_initial_spaces(&(bufptr3[1]));
	}

	// lexicographic position (or 0xffffffffU skip indicator)
	uii = flex_map[marker_in_idx];
	if (uii == 0xffffffffU) {
	  continue;
	}
	// actual relative position
	marker_out_idx = marker_map[uii] - start_marker_idx;

	if ((*aptr1 == '0') && (alen1 == 1)) {
          if ((*aptr2 != '0') || (alen2 != 1)) {
	    goto merge_main_ret_INVALID_FORMAT_3;
	  }
	  ucc2 = 1; // final PLINK encoding
	} else if ((*aptr2 == '0') && (alen2 == 1)) {
	  goto merge_main_ret_INVALID_FORMAT_3;
	} else {
	  ucc2 = 0; // A2 count
	  if (!strcmp(aptr1, marker_allele_ptrs[uii * 2 + 1])) {
	    ucc2++;
	  } else if (strcmp(aptr1, marker_allele_ptrs[uii * 2])) {
	    // new allele code
	    if (marker_allele_ptrs[uii * 2 + 1] == missing_geno_ptr) {
	      // fill A2 first
	      ucc2++;
	      ukk = uii * 2 + 1;
	    } else if (marker_allele_ptrs[uii * 2] == missing_geno_ptr) {
	      ukk = uii * 2;
	    } else {
	      printf("%s %s %s\n", marker_allele_ptrs[uii * 2 + 1], marker_allele_ptrs[uii * 2], aptr1);
	      goto merge_main_ret_INVALID_FORMAT;
	    }
	    if (allele_set(&(marker_allele_ptrs[ukk]), aptr1, alen1)) {
	      goto merge_main_ret_NOMEM;
	    }
	  }
	  if (!strcmp(aptr2, marker_allele_ptrs[uii * 2 + 1])) {
	    ucc2++;
	  } else if (strcmp(aptr2, marker_allele_ptrs[uii * 2])) {
	    // put A2 check second since the only way A2 will be unset when A1
	    // is set at this point is if it was specified that way in an
	    // earlier binary file.
	    if (marker_allele_ptrs[uii * 2] == missing_geno_ptr) {
	      ukk = uii * 2;
	    } else if (marker_allele_ptrs[uii * 2 + 1] == missing_geno_ptr) {
              ukk = uii * 2 + 1;
	    } else {
	      goto merge_main_ret_INVALID_FORMAT;
	    }
	    if (allele_set(&(marker_allele_ptrs[ukk]), aptr2, alen2)) {
	      goto merge_main_ret_NOMEM;
	    }
	  }
	  // transform to final PLINK encoding
	  if (ucc2) {
	    ucc2++;
	  }
	}
	switch (merge_mode) {
	case 1:
	  mbufptr2 = &(mbufptr[marker_out_idx * tot_indiv_ctl]);
	  wbufptr2 = &(wbufptr[marker_out_idx * tot_indiv_ct4]);
	  if ((*mbufptr2) & uljj) {
	    ucc3 = *wbufptr2;
	    if (((ucc3 >> ujj) & 3) != ucc2) {
	      *wbufptr2 = (ucc3 & ucc) | ucc4;
	    }
	  } else {
	    *mbufptr2 |= uljj;
	    *wbufptr2 = ((*wbufptr2) & ucc) | (ucc2 << ujj);
	  }
	  break;
	case 2:
	  wbufptr2 = &(wbufptr[marker_out_idx * tot_indiv_ct4]);
	  ucc3 = *wbufptr2;
	  if (((ucc3 >> ujj) & 3) == 1) {
	    *wbufptr2 = (ucc3 & ucc) | (ucc2 << ujj);
	  }
	  break;
	case 3:
	  if (ucc2 == 1) {
	    break;
	  }
	  // fall through
	case 5:
	  wbufptr2 = &(wbufptr[marker_out_idx * tot_indiv_ct4]);
	  *wbufptr2 = ((*wbufptr2) & ucc) | (ucc2 << ujj);
	  break;
	case 4:
	  mbufptr2 = &(mbufptr[marker_out_idx * tot_indiv_ctl]);
	  if (!((*mbufptr2) & uljj)) {
	    *mbufptr2 |= uljj;
	    wbufptr2 = &(wbufptr[marker_out_idx * tot_indiv_ct4]);
	    *wbufptr2 = ((*wbufptr2) & ucc) | (ucc2 << ujj);
	  }
	  break;
	case 6:
	  if (mbufptr[marker_out_idx * tot_indiv_ctl] & uljj) {
	    diff_total_overlap++;
	    ucc3 = ((wbufptr[marker_out_idx * tot_indiv_ct4] >> ujj) & 3);
	    umm = ((ucc2 == 1) || (ucc3 == 1));
	    if (umm) {
	      diff_not_both_genotyped++;
	    }
	    if (ucc2 != ucc3) {
	      if (!umm) {
		diff_discordant++;
	      }
	      if (merge_diff_print(outfile, idbuf, &(marker_ids[uii * max_marker_id_len]), &(person_ids[((uint32_t)ii) * max_person_id_len]), ucc2, ucc3, &(marker_allele_ptrs[uii * 2]))) {
		goto merge_main_ret_WRITE_FAIL;
	      }
	    }
	  }
	  break;
	case 7:
	  if (mbufptr[marker_out_idx * tot_indiv_ctl] & uljj) {
	    diff_total_overlap++;
	    ucc3 = ((wbufptr[marker_out_idx * tot_indiv_ct4] >> ujj) & 3);
	    if ((ucc2 == 1) || (ucc3 == 1)) {
	      diff_not_both_genotyped++;
	    } else if (ucc2 != ucc3) {
	      diff_discordant++;
	      if (merge_diff_print(outfile, idbuf, &(marker_ids[uii * max_marker_id_len]), &(person_ids[((uint32_t)ii) * max_person_id_len]), ucc2, ucc3, &(marker_allele_ptrs[uii * 2]))) {
		goto merge_main_ret_WRITE_FAIL;
	      }
	    }
	  }
	}
      }
      if (!is_eoln_kns(*bufptr3)) {
	fill_idbuf_fam_indiv(idbuf, bufptr, ' ');
	sprintf(logbuf, "Error: %s line has too many tokens (indiv id %s).\n", bedname, idbuf);
	goto merge_main_ret_INVALID_FORMAT_4;
      }
    }
    if (!feof(bedfile)) {
      goto merge_main_ret_READ_FAIL;
    }
  }
  if (merge_mode > 5) {
    *diff_total_overlap_ptr += diff_total_overlap;
    *diff_not_both_genotyped_ptr += diff_not_both_genotyped;
    *diff_discordant_ptr += diff_discordant;
  }
  while (0) {
  merge_main_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  merge_main_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  merge_main_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  merge_main_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  merge_main_ret_INVALID_FORMAT:
    sprintf(logbuf, "Error: Variant %s is not biallelic.\n", &(marker_ids[uii * max_marker_id_len]));
    putchar('\n');
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  merge_main_ret_INVALID_FORMAT_3:
    fill_idbuf_fam_indiv(idbuf, bufptr, ' ');
    sprintf(logbuf, "Error: Half-missing call in %s (indiv id %s, variant %s).\n", bedname, idbuf, &(marker_ids[uii * max_marker_id_len]));
    putchar('\n');
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  merge_main_ret_INVALID_FORMAT_2:
    fill_idbuf_fam_indiv(idbuf, bufptr, ' ');
    sprintf(logbuf, "Error: Line too short in %s (indiv id %s).\n", bedname, idbuf);
  merge_main_ret_INVALID_FORMAT_4:
    putchar('\n');
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
  fclose_cond(bedfile);
  fclose_cond(infile2);
  return retval;
}

int32_t merge_datasets(char* bedname, char* bimname, char* famname, char* outname, char* outname_end, char* mergename1, char* mergename2, char* mergename3, uint64_t calculation_type, uint32_t merge_type, uint32_t indiv_sort, uint64_t misc_flags, Chrom_info* chrom_info_ptr) {
  FILE* mergelistfile = NULL;
  FILE* outfile = NULL;
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t max_person_id_len = 0;
  uint32_t max_person_full_len = 0;
  uintptr_t max_marker_id_len = 0;
  uint32_t keep_allele_order = (misc_flags / MISC_KEEP_ALLELE_ORDER) & 1;
  uint32_t allow_extra_chroms = (misc_flags / MISC_ALLOW_EXTRA_CHROMS) & 1;
  uint32_t zero_extra_chroms = (misc_flags / MISC_ZERO_EXTRA_CHROMS) & 1;
  uint32_t is_dichot_pheno = 1;
  uint32_t merge_mode = (merge_type & MERGE_MODE_MASK);
  uint32_t merge_nsort = ((!indiv_sort) || (indiv_sort == INDIV_SORT_NATURAL))? 1 : 0;
  uint32_t merge_equal_pos = (merge_type & MERGE_EQUAL_POS)? 1 : 0;
  Ll_entry** htable = (Ll_entry**)(&(wkspace_base[wkspace_left - HASHMEM_S]));
  Ll_entry2** htable2 = (Ll_entry2**)(&(wkspace_base[wkspace_left - HASHMEM]));
  Ll_str* non_biallelics = NULL;
  uint32_t ped_buflen = MAXLINELEN;
  uint32_t max_bim_linelen = 0;
  char* missing_geno_ptr = (char*)g_missing_geno_ptr;
  char* pheno_c_char = NULL;
  double* pheno_d = NULL;
  uint32_t* indiv_nsmap = NULL;
  uint32_t max_cur_indiv_ct = 0;
  uint32_t max_cur_marker_text_ct = 0;
  uintptr_t* markbuf = NULL; // needed for merge modes 1, 4, 6, 7
  uint64_t diff_total_overlap = 0;
  uint64_t diff_not_both_genotyped = 0;
  uint64_t diff_discordant = 0;
  uint32_t orig_idx = 0;
  uint32_t position_warning_ct = 0;
  uint32_t cur_marker_ct = 0;
  uint32_t tot_marker_ct = 0;
  uint32_t* map_reverse = NULL;
  uintptr_t* reversed = NULL;
  char* bim_loadbuf = NULL;
  // N.B. marker_allele_ptrs are ordered by marker_id instead of position
  char** marker_allele_ptrs = NULL;
  unsigned char bmap_raw[576];
  unsigned char* bmap;
  unsigned char* bmap2;
  unsigned char* ucptr;
  unsigned char* ucptr_end;
  uintptr_t* pcptr;
  uintptr_t markers_per_pass;
  uint32_t pass_ct;
  uintptr_t topsize;
  char* person_ids;
  char* person_fids;
  char* marker_ids;
  uint32_t* marker_map;
  uint32_t* flex_map;
  double* marker_cms;
  uint32_t* pos_buf;
  int64_t* ll_buf;
  uintptr_t mlpos;
  uintptr_t merge_ct;
  char* idbuf;
  char* mergelist_buf;
  char** mergelist_bed;
  char** mergelist_bim;
  char** mergelist_fam;
  uint32_t cur_indiv_ct;
  uint32_t tot_indiv_ct;
  uint32_t tot_indiv_ct4;
  uint32_t dedup_marker_ct;
  uint64_t ullxx;
  int64_t llxx;
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t ulkk;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  char* bufptr;
  char* bufptr2;
  char* bufptr3;
  char* bufptr4;
  Ll_entry* ll_ptr;
  Ll_entry2* ll_ptr2;
  uint32_t chrom_start[MAX_POSSIBLE_CHROM + 1];
  uint32_t chrom_id[MAX_POSSIBLE_CHROM];
  uint32_t chrom_ct;
  unsigned char* readbuf;
  unsigned char* writebuf;
  unsigned char* ubufptr;
  char cc;
  unsigned char ucc;
  int32_t retval;

  if (!merge_mode) {
    merge_mode = 1;
  }
  if (merge_type & MERGE_LIST) {
    if (fopen_checked(&mergelistfile, mergename1, "r")) {
      goto merge_datasets_ret_READ_FAIL;
    }
    merge_ct = 1;
    ullxx = 0;
    // first pass: determine merge_ct, mergelist_buf size, verify no lines have
    // > 3 entries
    tbuf[MAXLINELEN - 1] = ' ';
    while (fgets(tbuf, MAXLINELEN, mergelistfile)) {
      if (!tbuf[MAXLINELEN - 1]) {
	logprint("Error: Pathologically long line in --merge-list file.\n");
	goto merge_datasets_ret_INVALID_FORMAT;
      }
      bufptr = skip_initial_spaces(tbuf);
      if (no_more_items_kns(bufptr)) {
	continue;
      }
      bufptr2 = next_item_mult(bufptr, 3);
      if (!no_more_items_kns(bufptr2)) {
	logprint("Error: More than three tokens on --merge-list line.\n");
        goto merge_datasets_ret_INVALID_FORMAT;
      }
      if (no_more_items_kns(next_item(bufptr))) {
	bufptr2 = item_endnn(bufptr);
	ulii = bufptr2 - bufptr;
	if (ulii > FNAMESIZE - 5) {
	  logprint("Error: Excessively long fileset prefix in --merge-list file.\n");
	  goto merge_datasets_ret_INVALID_FORMAT;
	}
	ullxx += 3 * ulii + 15;
      } else {
	do {
	  bufptr2 = item_endnn(bufptr);
	  ulii = bufptr2 - bufptr;
	  if (ulii > FNAMESIZE - 1) {
	    logprint("Error: Excessively long filename in --merge-list file.\n");
	    goto merge_datasets_ret_INVALID_FORMAT;
	  }
	  ullxx += ulii + 1;
	  bufptr = skip_initial_spaces(bufptr2);
	} while (!no_more_items_kns(bufptr));
      }
      merge_ct++;
    }
    if (!feof(mergelistfile)) {
      goto merge_datasets_ret_READ_FAIL;
    }
    if (!merge_ct) {
      logprint("Error: Empty --merge-list file.\n");
      goto merge_datasets_ret_INVALID_FORMAT;
    }
#ifndef __LP64__
    if (ullxx > 0x7fffffff) {
      goto merge_datasets_ret_NOMEM;
    }
#endif
    mergelist_bed = (char**)wkspace_alloc(merge_ct * sizeof(intptr_t));
    mergelist_bim = (char**)wkspace_alloc(merge_ct * sizeof(intptr_t));
    mergelist_fam = (char**)wkspace_alloc(merge_ct * sizeof(intptr_t));
    if (wkspace_alloc_c_checked(&mergelist_buf, (uintptr_t)ullxx)) {
      goto merge_datasets_ret_NOMEM;
    }
    rewind(mergelistfile);
    bufptr4 = mergelist_buf;
    mlpos = 1;
    while (fgets(tbuf, MAXLINELEN, mergelistfile)) {
      bufptr = skip_initial_spaces(tbuf);
      if (no_more_items_kns(bufptr)) {
	continue;
      }
      bufptr2 = item_endnn(bufptr);
      ulii = (bufptr2 - bufptr);
      bufptr3 = skip_initial_spaces(bufptr2);
      if (no_more_items_kns(bufptr3)) {
	mergelist_bed[mlpos] = bufptr4;
	bufptr4 = memcpya(memcpya(bufptr4, bufptr, ulii), ".bed", 5);
	mergelist_bim[mlpos] = bufptr4;
	bufptr4 = memcpya(memcpya(bufptr4, bufptr, ulii), ".bim", 5);
	mergelist_fam[mlpos] = bufptr4;
	bufptr4 = memcpya(memcpya(bufptr4, bufptr, ulii), ".fam", 5);
      } else {
	mergelist_bed[mlpos] = bufptr4;
	bufptr4 = memcpyax(bufptr4, bufptr, ulii, 0);
	bufptr2 = item_endnn(bufptr3);
	ulii = bufptr2 - bufptr3;
	bufptr = skip_initial_spaces(bufptr2);
	mergelist_bim[mlpos] = bufptr4;
	bufptr4 = memcpyax(bufptr4, bufptr3, ulii, 0);
	if (no_more_items_kns(bufptr)) {
	  mergelist_fam[mlpos] = NULL;
	} else {
	  bufptr2 = item_endnn(bufptr);
	  ulii = bufptr2 - bufptr;
	  mergelist_fam[mlpos] = bufptr4;
	  bufptr4 = memcpyax(bufptr4, bufptr, ulii, 0);
	}
      }
      if (++mlpos == merge_ct) {
	break;
      }
    }
    if ((mlpos < merge_ct) && (!feof(mergelistfile))) {
      goto merge_datasets_ret_READ_FAIL;
    }
    fclose_null(&mergelistfile);
  } else {
    merge_ct = 2;
    mergelist_bed = (char**)wkspace_alloc(2 * sizeof(intptr_t));
    mergelist_bim = (char**)wkspace_alloc(2 * sizeof(intptr_t));
    mergelist_fam = (char**)wkspace_alloc(2 * sizeof(intptr_t));
    mergelist_bed[1] = mergename1;
    mergelist_bim[1] = mergename2;
    mergelist_fam[1] = (merge_type & MERGE_BINARY)? mergename3 : NULL;
  }
  mergelist_bed[0] = bedname;
  mergelist_bim[0] = bimname;
  mergelist_fam[0] = famname;

  // ID counting/duplicate detection strategy:
  // - We do NOT want to scan through .ped files any more times than absolutely
  // necessary.  So we actually use *gasp* a hash table here.
  // - The hash table is positioned at the FAR end of wkspace, automatically
  // sized to ~4MB (or ~2MB on 32-bit systems).  IDs are then stored
  // backwards from there.  This simplifies copying into a sorted list.
  if (wkspace_left < HASHSIZE_S * sizeof(intptr_t)) {
    goto merge_datasets_ret_NOMEM;
  }
  for (uii = 0; uii < HASHSIZE_S; uii++) {
    htable[uii] = NULL;
  }
  topsize = HASHMEM_S;

  ullxx = 0;
  mlpos = 0;
  do {
    retval = merge_fam_id_scan(mergelist_bed[mlpos], mergelist_fam[mlpos], &max_person_id_len, &max_person_full_len, &is_dichot_pheno, htable, &topsize, &ullxx, &ped_buflen, &cur_indiv_ct, &orig_idx);
    if (retval) {
      goto merge_datasets_ret_1;
    }
    if (cur_indiv_ct > max_cur_indiv_ct) {
      max_cur_indiv_ct = cur_indiv_ct;
    }
  } while (++mlpos < merge_ct);
#ifdef __LP64__
  if (ullxx > 0x7fffffff) {
    sprintf(logbuf, "Error: Too many %s (max 2147483647).\n", g_species_plural);
    goto merge_datasets_ret_INVALID_FORMAT_2;
  }
#else
  // avoid integer overflow in wkspace_alloc calls
  if (ullxx * max_person_full_len > 0x7fffffff) {
    sprintf(logbuf, "Error: Too many %s for 32-bit " PROG_NAME_CAPS ".\n", g_species_plural);
    goto merge_datasets_ret_INVALID_FORMAT_2;
  }
#endif
  tot_indiv_ct = ullxx;
  fill_bmap_raw(bmap_raw, tot_indiv_ct % 4);
  bmap = &(bmap_raw[256]);
  bmap2 = &(bmap_raw[512]);
  // "allocate" first hash table off far side of stack before making regular
  // stack allocations
  wkspace_left -= topsize;
  if (indiv_sort == INDIV_SORT_NONE) {
    if (wkspace_alloc_ui_checked(&indiv_nsmap, tot_indiv_ct * sizeof(int32_t))) {
      goto merge_datasets_ret_NOMEM2;
    }
  }
  if (wkspace_alloc_c_checked(&person_ids, max_person_id_len * tot_indiv_ct) ||
      wkspace_alloc_c_checked(&person_fids, max_person_full_len * tot_indiv_ct)) {
    goto merge_datasets_ret_NOMEM2;
  }
  if (is_dichot_pheno) {
    if (wkspace_alloc_c_checked(&pheno_c_char, tot_indiv_ct)) {
      goto merge_datasets_ret_NOMEM2;
    }
  } else {
    if (wkspace_alloc_d_checked(&pheno_d, tot_indiv_ct * sizeof(double))) {
      goto merge_datasets_ret_NOMEM2;
    }
  }
  if (indiv_sort == INDIV_SORT_NONE) {
    if (wkspace_alloc_ui_checked(&map_reverse, tot_indiv_ct * sizeof(int32_t))) {
      goto merge_datasets_ret_NOMEM2;
    }
    for (uii = 0; uii < HASHSIZE_S; uii++) {
      if (htable[uii]) {
	ll_ptr = htable[uii];
	do {
	  ujj = ll_ptr->orig_order;
	  strcpy(&(person_fids[ujj * max_person_full_len]), ll_ptr->idstr);
	  if (is_dichot_pheno) {
	    if (ll_ptr->pheno == -9) {
	      pheno_c_char[ujj] = -1;
	    } else {
	      pheno_c_char[ujj] = ll_ptr->pheno - 1;
	    }
	  } else {
	    pheno_d[ujj] = ll_ptr->pheno;
	  }
	  ll_ptr = ll_ptr->next;
	} while (ll_ptr);
      }
    }
    for (uii = 0; uii < tot_indiv_ct; uii++) {
      indiv_nsmap[uii] = uii;
    }
    if (qsort_ext(person_fids, tot_indiv_ct, max_person_full_len, strcmp_deref, (char*)indiv_nsmap, sizeof(int32_t))) {
      goto merge_datasets_ret_NOMEM2;
    }
  } else {
    ulii = 0;
    bufptr = person_fids;
    for (uii = 0; uii < HASHSIZE_S; uii++) {
      if (htable[uii]) {
	ll_ptr = htable[uii];
	do {
	  strcpy(bufptr, ll_ptr->idstr);
	  bufptr = &(bufptr[max_person_full_len]);
	  if (is_dichot_pheno) {
	    if (ll_ptr->pheno == -9) {
	      pheno_c_char[ulii] = -1;
	    } else {
	      pheno_c_char[ulii] = ll_ptr->pheno - 1;
	    }
	  } else {
	    pheno_d[ulii] = ll_ptr->pheno;
	  }
	  ulii++;
	  ll_ptr = ll_ptr->next;
	} while (ll_ptr);
      }
    }
    if (is_dichot_pheno) {
      if (qsort_ext(person_fids, tot_indiv_ct, max_person_full_len, merge_nsort? strcmp_natural_deref : strcmp_deref, pheno_c_char, 1)) {
	goto merge_datasets_ret_NOMEM2;
      }
    } else {
      if (qsort_ext(person_fids, tot_indiv_ct, max_person_full_len, merge_nsort? strcmp_natural_deref : strcmp_deref, (char*)pheno_d, sizeof(double))) {
	goto merge_datasets_ret_NOMEM2;
      }
    }
  }
  wkspace_left += topsize; // deallocate first hash table
  if (merge_mode < 6) {
    memcpy(outname_end, ".fam", 5);
    if (fopen_checked(&outfile, outname, "w")) {
      goto merge_datasets_ret_OPEN_FAIL;
    }
  }
  if (indiv_sort == INDIV_SORT_NONE) {
    for (ulii = 0; ulii < tot_indiv_ct; ulii++) {
      map_reverse[indiv_nsmap[ulii]] = ulii;
    }
    for (ulii = 0; ulii < tot_indiv_ct; ulii++) {
      ujj = map_reverse[ulii];
      bufptr = &(person_fids[ujj * max_person_full_len]);
      bufptr3 = &(person_ids[ujj * max_person_id_len]);
      bufptr2 = next_item_mult(bufptr, 2);
      uii = (bufptr2 - bufptr) - 1;
      memcpyx(bufptr3, bufptr, uii, 0);
      if (merge_mode < 6) {
	uii += strlen(bufptr2) + 1;
	if (fwrite_checked(bufptr, uii, outfile)) {
	  goto merge_datasets_ret_WRITE_FAIL;
	}
	if (is_dichot_pheno) {
	  cc = pheno_c_char[ulii];
	  fprintf(outfile, "\t%s\n", cc? ((cc == 1)? "2" : "-9") : "1");
	} else {
	  fprintf(outfile, "\t%g\n", pheno_d[ulii]);
	}
      }
    }
  } else {
    bufptr = person_fids;
    bufptr3 = person_ids;
    for (ulii = 0; ulii < tot_indiv_ct; ulii++) {
      bufptr2 = next_item_mult(bufptr, 2);
      uii = (bufptr2 - bufptr) - 1;
      memcpyx(bufptr3, bufptr, uii, 0);
      bufptr3 = &(bufptr3[max_person_id_len]);
      if (merge_mode < 6) {
	uii += strlen(bufptr2) + 1;
	fwrite(bufptr, 1, uii, outfile);
	if (is_dichot_pheno) {
	  cc = pheno_c_char[ulii];
	  fprintf(outfile, "\t%s\n", cc? ((cc == 1)? "2" : "-9") : "1");
	} else {
	  fprintf(outfile, "\t%g\n", pheno_d[ulii]);
	}
      }
      bufptr = &(bufptr[max_person_full_len]);
    }
  }
  if (merge_mode < 6) {
    if (ferror(outfile)) {
      goto merge_datasets_ret_WRITE_FAIL;
    }
    if (fclose_null(&outfile)) {
      goto merge_datasets_ret_WRITE_FAIL;
    }
  }
  wkspace_reset((unsigned char*)person_fids);
  for (uii = 0; uii < HASHSIZE; uii++) {
    htable2[uii] = NULL;
  }
  topsize = HASHMEM;

  ullxx = 0;
  mlpos = 0;
  do {
    retval = merge_bim_scan(mergelist_bim[mlpos], (mergelist_fam[mlpos])? 1 : 0, &max_marker_id_len, htable2, &topsize, &max_bim_linelen, &ullxx, &cur_marker_ct, &position_warning_ct, &non_biallelics, allow_extra_chroms, chrom_info_ptr);
    if (retval) {
      goto merge_datasets_ret_1;
    }
    if (!mergelist_fam[mlpos]) {
      if (cur_marker_ct > max_cur_marker_text_ct) {
        max_cur_marker_text_ct = cur_marker_ct;
      }
    }
  } while (++mlpos < merge_ct);
#ifdef __LP64__
  if (ullxx > 0x7fffffff) {
    logprint("Error: Too many variants (max 2147483647).\n");
    goto merge_datasets_ret_INVALID_FORMAT;
  }
#else
  if (ullxx * MAXV(max_marker_id_len, 8) > 0x7fffffff) {
    logprint("Error: Too many variants for 32-bit " PROG_NAME_CAPS ".\n");
    goto merge_datasets_ret_INVALID_FORMAT;
  }
#endif
  if (non_biallelics) {
    wkspace_reset(wkspace_mark);
    retval = report_non_biallelics(outname, outname_end, non_biallelics);
    if (retval) {
      goto merge_datasets_ret_1;
    }
    goto merge_datasets_ret_INVALID_FORMAT;
  }
  tot_marker_ct = ullxx;
  // "allocate" second hash table off far side of stack before making regular
  // stack allocations
  wkspace_left -= topsize;
  marker_allele_ptrs = (char**)wkspace_alloc(tot_marker_ct * 2 * sizeof(intptr_t));
  if (!marker_allele_ptrs) {
    goto merge_datasets_ret_NOMEM2;
  }
  for (uii = 0; uii < tot_marker_ct * 2; uii++) {
    marker_allele_ptrs[uii] = NULL;
  }
  if (max_bim_linelen) {
    max_bim_linelen++;
    if (wkspace_alloc_c_checked(&bim_loadbuf, max_bim_linelen)) {
      goto merge_datasets_ret_NOMEM2;
    }
  }
  if (wkspace_alloc_c_checked(&marker_ids, max_marker_id_len * tot_marker_ct) ||
      wkspace_alloc_ui_checked(&marker_map, tot_marker_ct * sizeof(int32_t)) ||
      wkspace_alloc_d_checked(&marker_cms, tot_marker_ct * sizeof(double)) ||
      wkspace_alloc_ui_checked(&pos_buf, tot_marker_ct * sizeof(int32_t)) ||
      wkspace_alloc_ll_checked(&ll_buf, tot_marker_ct * sizeof(int64_t))) {
    goto merge_datasets_ret_NOMEM2;
  }
  for (uii = 0; uii < tot_marker_ct; uii++) {
    pos_buf[uii] = uii;
  }
  ulii = 0;
  for (uii = 0; uii < HASHSIZE; uii++) {
    if (htable2[uii]) {
      ll_ptr2 = htable2[uii];
      do {
	strcpy(&(marker_ids[ulii * max_marker_id_len]), ll_ptr2->idstr);
        ulii++;
	ll_ptr2 = ll_ptr2->next;
      } while (ll_ptr2);
    }
  }
  if (qsort_ext(marker_ids, tot_marker_ct, max_marker_id_len, strcmp_deref, (char*)pos_buf, sizeof(int32_t))) {
    goto merge_datasets_ret_NOMEM2;
  }
  // pos_buf[n] contains the position of lexicographic marker #n in the hash
  // table.  invert this map, then traverse the hash table.
  for (uii = 0; uii < tot_marker_ct; uii++) {
    marker_map[pos_buf[uii]] = uii;
  }
  wkspace_left += topsize; // deallocate second hash table
  ulii = 0;
  for (uii = 0; uii < HASHSIZE; uii++) {
    if (htable2[uii]) {
      ll_ptr2 = htable2[uii];
      do {
	ujj = marker_map[ulii++];
	llxx = ll_ptr2->pos;
	pos_buf[ujj] = (uint32_t)llxx;
	bufptr = ll_ptr2->allele[0];
	if (bufptr) {
          marker_allele_ptrs[ujj * 2] = bufptr;
	} else {
	  marker_allele_ptrs[ujj * 2] = missing_geno_ptr;
	}
	bufptr = ll_ptr2->allele[1];
	if (bufptr) {
	  marker_allele_ptrs[ujj * 2 + 1] = bufptr;
	} else {
	  marker_allele_ptrs[ujj * 2 + 1] = missing_geno_ptr;
	}
	marker_cms[ujj] = ll_ptr2->cm;
	ll_buf[ujj] = (llxx & 0xffffffff00000000LL) | ujj;
	ll_ptr2 = ll_ptr2->next;
      } while (ll_ptr2);
    }
  }
  sort_marker_chrom_pos(ll_buf, tot_marker_ct, pos_buf, chrom_start, chrom_id, &chrom_ct);
  if (merge_post_msort_update_maps(marker_ids, max_marker_id_len, marker_map, pos_buf, ll_buf, chrom_start, chrom_id, chrom_ct, &dedup_marker_ct, merge_equal_pos, marker_allele_ptrs, chrom_info_ptr)) {
    goto merge_datasets_ret_INVALID_FORMAT;
  }
  if (!dedup_marker_ct) {
    logprint("Error: No variants in merged file.\n");
    goto merge_datasets_ret_INVALID_FORMAT;
  }
  wkspace_reset((unsigned char*)ll_buf);

  tot_indiv_ct4 = (tot_indiv_ct + 3) / 4;

  if (!keep_allele_order) {
    ulii = (tot_marker_ct + (BITCT - 1)) / BITCT;
    if (wkspace_alloc_ul_checked(&reversed, ulii * sizeof(intptr_t))) {
      goto merge_datasets_ret_NOMEM;
    }
    fill_ulong_zero(reversed, ulii);
  }
  if (wkspace_alloc_ui_checked(&flex_map, MAXV(max_cur_indiv_ct, max_cur_marker_text_ct) * sizeof(int32_t)) ||
      wkspace_alloc_c_checked(&idbuf, MAXV(max_marker_id_len, max_person_id_len))) {
    goto merge_datasets_ret_NOMEM;
  }

  if (tot_indiv_ct4 > ped_buflen) {
    ulii = tot_indiv_ct4;
  } else {
    ulii = ped_buflen;
  }
  if (wkspace_alloc_uc_checked(&readbuf, ulii)) {
    goto merge_datasets_ret_NOMEM;
  }
  if (merge_must_track_write(merge_mode)) {
    ulii = (tot_indiv_ct + (BITCT - 1)) / BITCT;
    markers_per_pass = wkspace_left / (3 * sizeof(intptr_t) * ulii);
    if (markers_per_pass > dedup_marker_ct) {
      markers_per_pass = dedup_marker_ct;
    }
    markbuf = (uintptr_t*)wkspace_alloc(markers_per_pass * ulii * sizeof(intptr_t));
  } else {
    markers_per_pass = wkspace_left / tot_indiv_ct4;
    if (markers_per_pass > dedup_marker_ct) {
      markers_per_pass = dedup_marker_ct;
    }
  }
  if (!markers_per_pass) {
    goto merge_datasets_ret_NOMEM;
  }
  pass_ct = 1 + ((dedup_marker_ct - 1) / markers_per_pass);

  writebuf = wkspace_base;
  pcptr = (uintptr_t*)wkspace_base;
  if (merge_mode < 6) {
    memcpy(outname_end, ".bed", 5);
    if (fopen_checked(&outfile, outname, "wb")) {
      goto merge_datasets_ret_OPEN_FAIL;
    }
    if (fwrite_checked("l\x1b\x01", 3, outfile)) {
      goto merge_datasets_ret_WRITE_FAIL;
    }
    if (pass_ct == 1) {
      sprintf(logbuf, "Performing single-pass merge (%u %s, %u variant%s).\n", tot_indiv_ct, species_str(tot_indiv_ct), dedup_marker_ct, (dedup_marker_ct == 1)? "" : "s");
    } else {
      sprintf(logbuf, "Performing %u-pass merge (%u %s, %" PRIuPTR "/%u variant%s per pass).\n", pass_ct, tot_indiv_ct, species_str(tot_indiv_ct), markers_per_pass, dedup_marker_ct, (dedup_marker_ct == 1)? "" : "s");
    }
  } else {
    memcpy(outname_end, ".diff", 6);
    if (fopen_checked(&outfile, outname, "w")) {
      goto merge_datasets_ret_OPEN_FAIL;
    }
    if (fputs_checked("                 SNP                  FID                  IID      NEW      OLD \n", outfile)) {
      goto merge_datasets_ret_WRITE_FAIL;
    }
    sprintf(logbuf, "Performing %u-pass diff (mode %u), writing results to %s.\n", pass_ct, merge_mode, outname);
  }
  logprintb();
  for (uii = 0; uii < pass_ct; uii++) {
    if (uii + 1 == pass_ct) {
      ujj = dedup_marker_ct - markers_per_pass * uii;
    } else {
      ujj = markers_per_pass;
    }
    if (tot_indiv_ct % 4) {
      umm = tot_indiv_ct / 4;
      ubufptr = writebuf;
      ucc = 0x55 >> ((tot_indiv_ct % 4) * 2);
      for (ukk = 0; ukk < ujj; ukk++) {
        memset(ubufptr, 0x55, umm);
        ubufptr[umm] = ucc;
	ubufptr = &(ubufptr[tot_indiv_ct4]);
      }
    } else {
      memset(writebuf, 0x55, ((uintptr_t)ujj) * tot_indiv_ct4);
    }
    if (merge_must_track_write(merge_mode)) {
      fill_ulong_zero(markbuf, ujj * ulii);
    }
    for (mlpos = 0; mlpos < merge_ct; mlpos++) {
      retval = merge_main(mergelist_bed[mlpos], mergelist_bim[mlpos], mergelist_fam[mlpos], bim_loadbuf, max_bim_linelen, tot_indiv_ct, tot_marker_ct, dedup_marker_ct, uii * markers_per_pass, ujj, marker_allele_ptrs, marker_ids, max_marker_id_len, person_ids, max_person_id_len, merge_nsort, indiv_nsmap, flex_map, marker_map, idbuf, readbuf, writebuf, mlpos? merge_mode : merge_first_mode(merge_mode, merge_equal_pos), markbuf, outfile, &diff_total_overlap, &diff_not_both_genotyped, &diff_discordant, ped_buflen, bmap_raw);
      if (retval) {
	goto merge_datasets_ret_1;
      }
      if (mlpos != merge_ct - 1) {
        printf("\rPass %u: fileset #%" PRIuPTR " complete.", uii + 1, mlpos + 1);
	fflush(stdout);
      }
    }
    if (merge_mode < 6) {
      if (!keep_allele_order) {
	for (ukk = 0; ukk < ujj; ukk++) {
	  uljj = ((uintptr_t)ukk) * tot_indiv_ct4;
	  umm = popcount_chars(pcptr, uljj, uljj + tot_indiv_ct4);
	  if (umm < tot_indiv_ct) {
	    ulkk = (uii * markers_per_pass) + ukk;
	    reversed[ulkk / BITCT] |= (ONELU << (ulkk % BITCT));
	    ucptr = &(writebuf[uljj]);
	    ucptr_end = &(writebuf[uljj + tot_indiv_ct / 4]);
	    while (ucptr < ucptr_end) {
	      *ucptr = bmap[*ucptr];
	      ucptr++;
	    }
	    if (tot_indiv_ct % 4) {
	      *ucptr = bmap2[*ucptr];
	    }
	  }
	}
      }
      if (fwrite_checked(writebuf, ((uintptr_t)ujj) * tot_indiv_ct4, outfile)) {
        goto merge_datasets_ret_WRITE_FAIL;
      }
    }
    fputs("\r                                              \r", stdout);
    if (uii + 1 != pass_ct) {
      sprintf(logbuf, "Pass %u complete.\n", uii + 1);
      logprintb();
    }
  }
  if (fclose_null(&outfile)) {
    goto merge_datasets_ret_WRITE_FAIL;
  }
  wkspace_reset((unsigned char*)flex_map);
  if (wkspace_alloc_ui_checked(&map_reverse, dedup_marker_ct * sizeof(int32_t))) {
    goto merge_datasets_ret_NOMEM;
  }
  if (merge_mode < 6) {
    memcpy(outname_end, ".bim", 5);
    if (fopen_checked(&outfile, outname, "w")) {
      goto merge_datasets_ret_OPEN_FAIL;
    }
    uii = tot_marker_ct;
    while (uii--) {
      map_reverse[marker_map[uii]] = uii;
    }
    for (ulii = 0; ulii < chrom_ct; ulii++) {
      uii = chrom_start[ulii + 1];
      ujj = chrom_start[ulii];
      ukk = chrom_id[ulii];
      for (; ujj < uii; ujj++) {
	umm = map_reverse[ujj];
	bufptr = chrom_name_write(tbuf, chrom_info_ptr, ukk, zero_extra_chroms);
	fwrite(tbuf, 1, bufptr - tbuf, outfile);
	if (keep_allele_order || (!IS_SET(reversed, ujj))) {
	  bufptr = marker_allele_ptrs[2 * umm];
	  bufptr2 = marker_allele_ptrs[2 * umm + 1];
	} else {
	  bufptr = marker_allele_ptrs[2 * umm + 1];
	  bufptr2 = marker_allele_ptrs[2 * umm];
	}
	fprintf(outfile, "\t%s\t%g\t%u\t%s\t%s\n", &(marker_ids[map_reverse[ujj] * max_marker_id_len]), marker_cms[ujj], pos_buf[ujj], bufptr, bufptr2);
      }
      if (ferror(outfile)) {
	goto merge_datasets_ret_WRITE_FAIL;
      }
    }
    if (fclose_null(&outfile)) {
      goto merge_datasets_ret_WRITE_FAIL;
    }
    memcpy(outname_end, ".bed", 5);
    sprintf(logbuf, "Merged fileset written to %s + .bim + .fam.\n", outname);
    logprintb();
  } else {
    // undo the "not"
    diff_not_both_genotyped = diff_total_overlap - diff_not_both_genotyped;
    sprintf(logbuf, "%" PRIu64 " overlapping variants, %" PRIu64 " genotyped in both filesets.\n%" PRIu64 " concordant, for a concordance rate of %g.\n", diff_total_overlap, diff_not_both_genotyped, diff_not_both_genotyped - diff_discordant, 1.0 - (((double)diff_discordant) / ((double)diff_not_both_genotyped)));
    logprintb();
  }

  forget_extra_chrom_names(chrom_info_ptr);
  while (0) {
  merge_datasets_ret_NOMEM2:
    wkspace_left += topsize;
  merge_datasets_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  merge_datasets_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  merge_datasets_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  merge_datasets_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  merge_datasets_ret_INVALID_FORMAT_2:
    logprintb();
  merge_datasets_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 merge_datasets_ret_1:
  if (marker_allele_ptrs) {
    for (uii = 0; uii < tot_marker_ct * 2; uii++) {
      bufptr = marker_allele_ptrs[uii];
      if (bufptr && ((bufptr < g_one_char_strs) || (bufptr >= (&(g_one_char_strs[512]))))) {
	free(bufptr);
      }
    }
  }
  fclose_cond(mergelistfile);
  fclose_cond(outfile);
  wkspace_reset(wkspace_mark);
  return retval;
}
