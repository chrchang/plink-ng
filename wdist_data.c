#include <fcntl.h>
#include <time.h>

#ifndef _WIN32
#include <sys/mman.h>
#endif

#include <sys/stat.h>
#include <sys/types.h>
#include "wdist_common.h"

#define PHENO_EPSILON 0.000030517578125

// last prime before 2^19
#define HASHSIZE 524287
#define HASHSIZE_S 524287

#ifdef __LP64__
#define HASHMEM 4194304
#define HASHMEM_S 4194304
#else
#define HASHMEM 2097152
#define HASHMEM_S 2097152
#endif

int32_t sort_item_ids_nx(char** sorted_ids_ptr, int32_t** id_map_ptr, uintptr_t item_ct, char* item_ids, uintptr_t max_id_len) {
  // Version of sort_item_ids() with no exclusion.
  uintptr_t ulii;
  uintptr_t uljj;
  char* sorted_ids;
  sorted_ids = (char*)wkspace_alloc(item_ct * max_id_len);
  if (!sorted_ids) {
    return RET_NOMEM;
  }
  *sorted_ids_ptr = sorted_ids;
  *id_map_ptr = (int32_t*)wkspace_alloc(item_ct * sizeof(int32_t));
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
  uljj = item_ct - 1;
  for (ulii = 0; ulii < uljj; ulii++) {
    if (!strcmp(&(sorted_ids[ulii * max_id_len]), &(sorted_ids[(ulii + 1) * max_id_len]))) {
      logprint("Error: Duplicate IDs.\n");
      return RET_INVALID_FORMAT;
    }
  }
  return 0;
}

int32_t sort_item_ids(char** sorted_ids_ptr, int32_t** id_map_ptr, uintptr_t unfiltered_ct, uintptr_t* exclude_arr, uintptr_t exclude_ct, char* item_ids, uintptr_t max_id_len, int(* comparator_deref)(const void*, const void*), int32_t* duplicate_fail_ptr) {
  // Allocates space for *id_map_ptr and *sorted_ids_ptr from wkspace; then
  // stores a lexicographically sorted list of IDs in *sorted_ids_ptr and
  // the raw positions of the corresponding markers/indivs in *id_map_ptr.
  // Does not include excluded markers/indivs in the list.
  // If *duplicate_fail_ptr is nonzero, this fails out if duplicate IDs are
  // found (RET_INVALID_FORMAT is returned).  If it is zero, it is set to one
  // if there are duplicates.
  uintptr_t item_ct = unfiltered_ct - exclude_ct;
  uint32_t uii;
  uint32_t ujj;
  char* sorted_ids;
  *id_map_ptr = (int32_t*)wkspace_alloc(item_ct * sizeof(int32_t));
  if (!(*id_map_ptr)) {
    return RET_NOMEM;
  }
  sorted_ids = (char*)wkspace_alloc(item_ct * max_id_len);
  if (!sorted_ids) {
    return RET_NOMEM;
  }
  if (!item_ct) {
    return 0;
  }
  *sorted_ids_ptr = sorted_ids;
  uii = 0;
  for (ujj = 0; ujj < item_ct; ujj++) {
    uii = next_non_set_unsafe(exclude_arr, uii);
    memcpy(&(sorted_ids[ujj * max_id_len]), &(item_ids[uii * max_id_len]), max_id_len);
    (*id_map_ptr)[ujj] = uii++;
  }
  if (qsort_ext(sorted_ids, item_ct, max_id_len, comparator_deref, (char*)(*id_map_ptr), sizeof(int32_t))) {
    return RET_NOMEM;
  }
  ujj = item_ct - 1;
  for (uii = 0; uii < ujj; uii++) {
    if (!strcmp(&(sorted_ids[uii * max_id_len]), &(sorted_ids[(uii + 1) * max_id_len]))) {
      if (*duplicate_fail_ptr) {
        logprint("Error: Duplicate IDs.\n");
        return RET_INVALID_FORMAT;
      }
      *duplicate_fail_ptr = 1;
      break;
    }
  }
  return 0;
}

#ifdef _WIN32
int32_t indiv_major_to_snp_major(char* indiv_major_fname, char* outname, FILE** outfile_ptr, uintptr_t unfiltered_marker_ct) {
  logprint("Error: Win32 WDIST does not yet support transposition of individual-major .bed\nfiles.  Contact the developers if you need this.\n");
  return RET_CALC_NOT_YET_SUPPORTED;
}
#else
int32_t indiv_major_to_snp_major(char* indiv_major_fname, char* outname, FILE** outfile_ptr, uintptr_t unfiltered_marker_ct) {
  // This implementation only handles large files on 64-bit Unix systems; a
  // more portable version needs to be written for Windows and 32-bit Unix.
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
  *outfile_ptr = NULL;
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
  if (fopen_checked(outfile_ptr, outname, "wb")) {
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
  if (fwrite_checked("l\x1b\x01", 3, *outfile_ptr)) {
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
    if (fwrite_checked(wkspace_base, ((int64_t)block_last_marker) * indiv_ct4, *outfile_ptr)) {
      goto indiv_major_to_snp_major_ret_WRITE_FAIL;
    }
    superblock_offset += max_4blocks_in_mem;
  }
  retval = 0;
  while (0) {
  indiv_major_to_snp_major_ret_INVALID_FORMAT:
    sprintf(logbuf, "Error: %s's file size is inconsistent with the marker count.\n", indiv_major_fname);
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
  if (*outfile_ptr) {
    fclose(*outfile_ptr);
    *outfile_ptr = NULL;
  }
  if (in_contents != MAP_FAILED) {
    munmap(in_contents, sb.st_size);
  }
  close(in_fd);
  return retval;
}
#endif

const char errstr_map_format[] = "Error: Improperly formatted .map file.\n";

int32_t load_map(FILE** mapfile_ptr, char* mapname, int32_t* map_cols_ptr, uintptr_t* unfiltered_marker_ct_ptr, uintptr_t* marker_exclude_ct_ptr, uintptr_t* max_marker_id_len_ptr, uintptr_t** marker_exclude_ptr, char** marker_ids_ptr, Chrom_info* chrom_info_ptr, uint32_t** marker_pos_ptr, int32_t* map_is_unsorted_ptr) {
  // todo: some cleanup
  uintptr_t unfiltered_marker_ct = 0;
  uintptr_t max_marker_id_len = 0;
  uint64_t loaded_chrom_mask = 0;
  int32_t last_chrom = -1;
  int32_t last_pos = 0;
  int32_t marker_pos_needed = 0;
  uint32_t species = chrom_info_ptr->species;
  uint32_t unfiltered_marker_ctl;
  char* bufptr;
  uintptr_t ulii;
  int32_t jj;
  int32_t cur_pos;
  int32_t chroms_encountered_m1 = -1;
  uintptr_t marker_uidx;
  if (fopen_checked(mapfile_ptr, mapname, "r")) {
    return RET_OPEN_FAIL;
  }
  // first pass: count columns, determine raw marker count, determine maximum
  // marker ID length if necessary.
  tbuf[MAXLINELEN - 6] = ' ';
  while (fgets(tbuf, MAXLINELEN - 5, *mapfile_ptr) != NULL) {
    if (!tbuf[MAXLINELEN - 6]) {
      sprintf(logbuf, "Error: Excessively long line in .map/.bim file (max %u chars).\n", MAXLINELEN - 8);
      logprintb();
      return RET_INVALID_FORMAT;
    }
    bufptr = skip_initial_spaces(tbuf);
    if (is_eoln_or_comment(*bufptr)) {
      continue;
    }
    bufptr = next_item(bufptr);
    if (no_more_items_kns(bufptr)) {
      logprint(errstr_map_format);
      return RET_INVALID_FORMAT;
    }
    ulii = strlen_se(bufptr) + 1;
    if (ulii > max_marker_id_len) {
      max_marker_id_len = ulii;
    }
    if (!unfiltered_marker_ct) {
      bufptr = next_item_mult(bufptr, 2);
      if (!bufptr) {
	logprint(errstr_map_format);
	return RET_INVALID_FORMAT;
      }
      if (*bufptr > ' ') {
	*map_cols_ptr = 4;
      }
    }
    unfiltered_marker_ct++;
  }
  if (!feof(*mapfile_ptr)) {
    return RET_READ_FAIL;
  }
  if (!unfiltered_marker_ct) {
    logprint("Error: No markers in .map/.bim file.");
    return RET_INVALID_FORMAT;
  }
  *unfiltered_marker_ct_ptr = unfiltered_marker_ct;
  *max_marker_id_len_ptr = max_marker_id_len;
  rewind(*mapfile_ptr);
  unfiltered_marker_ctl = (unfiltered_marker_ct + (BITCT - 1)) / BITCT;

  // unfiltered_marker_ct can be very large, so use wkspace for all allocations
  // that are a multiple of it

  // permanent stack allocation #1: marker_exclude
  if (wkspace_alloc_ul_checked(marker_exclude_ptr, unfiltered_marker_ctl * sizeof(intptr_t))) {
    return RET_NOMEM;
  }
  fill_ulong_zero(*marker_exclude_ptr, unfiltered_marker_ctl);
  fill_uint_zero(chrom_info_ptr->chrom_file_order, MAX_POSSIBLE_CHROM);
  fill_uint_zero(chrom_info_ptr->chrom_file_order_marker_idx, MAX_POSSIBLE_CHROM + 1);
  fill_uint_zero(chrom_info_ptr->chrom_start, MAX_POSSIBLE_CHROM);
  fill_uint_zero(chrom_info_ptr->chrom_end, MAX_POSSIBLE_CHROM);
  // permanent stack allocation #3, if needed: marker_pos
  if (marker_pos_needed) {
    *marker_pos_ptr = (uint32_t*)wkspace_alloc(unfiltered_marker_ct * sizeof(int32_t));
    if (!(*marker_pos_ptr)) {
      return RET_NOMEM;
    }
  }
  *marker_ids_ptr = (char*)wkspace_alloc(unfiltered_marker_ct * max_marker_id_len);
  if (!(*marker_ids_ptr)) {
    return RET_NOMEM;
  }

  // second pass: actually load stuff
  for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
    if (get_next_noncomment(*mapfile_ptr, &bufptr)) {
      return RET_READ_FAIL;
    }
    jj = marker_code(species, bufptr);
    if (jj == -1) {
      logprint("Error: Invalid chromosome index in .map/.bim file.\n");
      return RET_INVALID_FORMAT;
    }
    if (jj != last_chrom) {
      if (last_chrom != -1) {
	chrom_info_ptr->chrom_end[last_chrom] = marker_uidx;
      }
      if (jj < last_chrom) {
	if (!(*map_is_unsorted_ptr)) {
	  *map_is_unsorted_ptr = 1;
	}
      }
      last_chrom = jj;
      if (loaded_chrom_mask & (1LLU << jj)) {
	*map_is_unsorted_ptr = 2;
      } else {
	loaded_chrom_mask |= 1LLU << jj;
	chrom_info_ptr->chrom_start[jj] = marker_uidx;
	chrom_info_ptr->chrom_file_order[++chroms_encountered_m1] = jj;
	chrom_info_ptr->chrom_file_order_marker_idx[chroms_encountered_m1] = marker_uidx;
      }
      last_pos = 0;
    }

    if (!(chrom_info_ptr->chrom_mask & (1LLU << jj))) {
      set_bit(*marker_exclude_ptr, marker_uidx, marker_exclude_ct_ptr);
    } else {
      bufptr = next_item(bufptr);
      if (no_more_items_kns(bufptr)) {
	logprint(errstr_map_format);
	return RET_INVALID_FORMAT;
      }
      read_next_terminate(&((*marker_ids_ptr)[marker_uidx * max_marker_id_len]), bufptr);
      bufptr = next_item_mult(bufptr, *map_cols_ptr - 2);
      if (no_more_items_kns(bufptr)) {
        logprint(errstr_map_format);
        return RET_INVALID_FORMAT;
      }
      if (*bufptr == '-') {
	// negative marker positions now have same effect in .bim as .map
	set_bit(*marker_exclude_ptr, marker_uidx, marker_exclude_ct_ptr);
      } else {
	cur_pos = atoi(bufptr);
	if (cur_pos < last_pos) {
	  *map_is_unsorted_ptr = 2;
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
  return 0;
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
  do {
    if ((cur_pos >= sf_pos[cur_idx]) && (cur_pos <= sf_pos[cur_idx + 1])) {
      return 0;
    }
    cur_idx += 2;
  } while (cur_idx < end_idx);
  return 1;
}

int32_t load_bim(FILE** bimfile_ptr, char* mapname, int32_t* map_cols_ptr, uintptr_t* unfiltered_marker_ct_ptr, uintptr_t* marker_exclude_ct_ptr, uintptr_t* max_marker_id_len_ptr, uint32_t* plink_maxsnp_ptr, uintptr_t** marker_exclude_ptr, double** set_allele_freqs_ptr, char** marker_alleles_ptr, uintptr_t* max_marker_allele_len_ptr, char** marker_ids_ptr, Chrom_info* chrom_info_ptr, uint32_t** marker_pos_ptr, char* extractname, char* excludename, char* freqname, char* refalleles, uint64_t calculation_type, uint32_t recode_modifier, int32_t allelexxxx, int32_t marker_pos_start, int32_t marker_pos_end, uint32_t snp_window_size, char* markername_from, char* markername_to, char* markername_snp, char* sf_markers, unsigned char* sf_starts_range, uint32_t sf_ct, uint32_t sf_max_len, int32_t* map_is_unsorted_ptr) {
  unsigned char* wkspace_mark = wkspace_base;
  // may want to just always load
  uint32_t marker_alleles_needed = (freqname || (calculation_type & (CALC_FREQ | CALC_HARDY | CALC_MAKE_BED | CALC_RECODE | CALC_REGRESS_PCS | CALC_MODEL)) || allelexxxx);
  uintptr_t unfiltered_marker_ct = 0;
  uintptr_t max_marker_id_len = 0;
  uintptr_t max_marker_allele_len = 1;
  uint32_t plink_maxsnp = 4;
  uint64_t loaded_chrom_mask = 0;
  int32_t last_chrom = -1;
  int32_t last_pos = 0;
  int32_t marker_pos_needed = (((calculation_type & (CALC_GENOME | CALC_LD_PRUNE | CALC_REGRESS_PCS | CALC_MODEL)) != 0) || ((calculation_type & CALC_RECODE) && (recode_modifier & RECODE_23)));
  uint32_t species = chrom_info_ptr->species;
  uint32_t from_slen = markername_from? strlen(markername_from) : 0;
  uint32_t to_slen = markername_to? strlen(markername_to) : 0;
  uint32_t snp_slen = markername_snp? strlen(markername_snp) : 0;
  uint32_t slen_check = from_slen || to_slen || snp_slen || sf_ct;
  uint32_t from_chrom = MAX_POSSIBLE_CHROM;
  uint32_t to_chrom = MAX_POSSIBLE_CHROM;
  uint32_t snp_chrom = MAX_POSSIBLE_CHROM;
  int32_t chroms_encountered_m1 = -1;
  uint32_t snp_pos = 0;
  uint32_t mcm2 = 1;
  uint32_t* sf_pos = NULL;
  uint32_t sf_entry_ct;
  uint32_t* sf_str_chroms = NULL;
  uint32_t* sf_str_pos = NULL;
  uint32_t* sf_str_lens = NULL;
  uint32_t* sf_llbuf = NULL;
  char* marker_alleles = NULL;
  uint64_t sf_mask;
  uint32_t sf_lltop;
  uint32_t sf_start_idxs[MAX_POSSIBLE_CHROM + 1];
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
  int32_t ii;
  int32_t jj;
  int32_t cur_pos;
  uintptr_t marker_uidx;
  int32_t retval;
  if (sf_ct) {
    if (wkspace_alloc_ui_checked(&sf_str_chroms, sf_ct * sizeof(int32_t)) ||
	wkspace_alloc_ui_checked(&sf_str_pos, sf_ct * sizeof(int32_t)) ||
        wkspace_alloc_ui_checked(&sf_str_lens, sf_ct * sizeof(int32_t)) ||
        wkspace_alloc_ui_checked(&sf_llbuf, 3 * (species_max_code[species] + sf_ct) * sizeof(int32_t))) {
      goto load_bim_ret_NOMEM;
    }
    for (uii = 0; uii < sf_ct; uii++) {
      sf_str_chroms[uii] = MAX_POSSIBLE_CHROM;
      sf_str_lens[uii] = strlen(&(sf_markers[uii * sf_max_len]));
    }
  }
  if (fopen_checked(bimfile_ptr, mapname, "r")) {
    goto load_bim_ret_OPEN_FAIL;
  }
  // first pass: count columns, determine raw marker count,  determine maximum
  // marker ID length and/or marker allele length if necessary.
  tbuf[MAXLINELEN - 6] = ' ';
  while (fgets(tbuf, MAXLINELEN - 5, *bimfile_ptr) != NULL) {
    if (!tbuf[MAXLINELEN - 6]) {
      sprintf(logbuf, "Error: Excessively long line in .bim file (max %u chars).\n", MAXLINELEN - 8);
      goto load_bim_ret_INVALID_FORMAT_2;
    }
    // bufptr = col 1 start
    bufptr = skip_initial_spaces(tbuf);
    if (is_eoln_or_comment(*bufptr)) {
      continue;
    }
    ii = marker_code(species, bufptr);
    if (ii == -1) {
      jj = strlen_se(bufptr);
      bufptr[jj] = '\0';
      sprintf(logbuf, "Error: Invalid chromosome index '%s' in .bim file.\n", bufptr);
      goto load_bim_ret_INVALID_FORMAT_2;
    }

    // bufpter = col 2 start
    bufptr = next_item(bufptr);
    if (no_more_items_kns(bufptr)) {
      goto load_bim_ret_INVALID_FORMAT_5;
    }
    ulii = strlen_se(bufptr) + 1;
    if (ulii > max_marker_id_len) {
      max_marker_id_len = ulii;
    }
    if (ulii > (plink_maxsnp + 1)) {
      plink_maxsnp = ulii + 1;
    }
    if (!unfiltered_marker_ct) {
      // bufptr2 = col 5 start
      bufptr2 = next_item_mult(bufptr, 3);
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
      if (uii > max_marker_allele_len) {
	max_marker_allele_len = uii;
      }
      uii = strlen_se(bufptr3);
      if (uii > max_marker_allele_len) {
	max_marker_allele_len = uii;
      }
    }
    if (slen_check) {
      ulii--;
      if (sf_ct) {
	uii = 0;
	do {
	  if ((ulii == sf_str_lens[uii]) && (!memcmp(bufptr, &(sf_markers[uii * sf_max_len]), ulii))) {
	    if (sf_str_chroms[uii] != MAX_POSSIBLE_CHROM) {
	      goto load_bim_ret_INVALID_FORMAT_3;
	    }
	    sf_str_chroms[uii] = ii;
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
	  from_chrom = ii;
	  if (to_chrom != MAX_POSSIBLE_CHROM) {
	    if (from_chrom != to_chrom) {
	      goto load_bim_ret_INVALID_FORMAT_4;
	    }
	  }
	  chrom_info_ptr->chrom_mask = 1LLU << from_chrom;
	}
	if ((ulii == to_slen) && (!memcmp(bufptr, markername_to, ulii))) {
	  if (to_chrom != MAX_POSSIBLE_CHROM) {
	    goto load_bim_ret_INVALID_FORMAT_3;
	  }
	  if (load_bim_scan_position(&uii, bufptr, mcm2)) {
	    goto load_bim_ret_INVALID_FORMAT;
	  }
	  marker_pos_end = uii;
	  to_chrom = ii;
	  if (from_chrom != MAX_POSSIBLE_CHROM) {
	    if (to_chrom != from_chrom) {
	      goto load_bim_ret_INVALID_FORMAT_4;
	    }
	  }
	  chrom_info_ptr->chrom_mask = 1LLU << to_chrom;
	}
	if ((ulii == snp_slen) && (!memcmp(bufptr, markername_snp, ulii))) {
	  if (snp_chrom != MAX_POSSIBLE_CHROM) {
	    goto load_bim_ret_INVALID_FORMAT_3;
	  }
	  if (load_bim_scan_position(&snp_pos, bufptr, mcm2)) {
	    goto load_bim_ret_INVALID_FORMAT;
	  }
	  snp_chrom = ii;
	  chrom_info_ptr->chrom_mask = 1LLU << snp_chrom;
	}
      }
    }
    unfiltered_marker_ct++;
  }
  if (sf_ct) {
    for (uii = 0; uii < sf_ct; uii++) {
      if (sf_str_chroms[uii] == MAX_POSSIBLE_CHROM) {
	sprintf(logbuf, "Error: Marker '%s' not found in .bim file.\n", &(sf_markers[uii * sf_max_len]));
	goto load_bim_ret_INVALID_FORMAT_2;
      }
    }
    // effectively build out one linked list per chromosome
    sf_mask = chrom_info_ptr->chrom_mask;
    sf_entry_ct = 0;
    sf_lltop = 0;
    ujj = species_max_code[species];
    for (uii = 0; uii <= ujj; uii++) {
      sf_start_idxs[uii] = 1; // impossible (multiples of 3)
    }
    uii = 0;
    do {
      ujj = sf_str_chroms[uii];
      ukk = sf_str_pos[uii];
      if (sf_starts_range[uii]) {
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
	  if (sf_mask & (1LLU << ujj)) {
	    load_bim_sf_insert(ujj, ukk, 2147483647, sf_start_idxs, sf_llbuf, &sf_lltop, &sf_entry_ct);
	  }
	  for (uoo = ujj + 1; uoo < umm; uoo++) {
	    if (sf_mask & (1LLU << uoo)) {
	      load_bim_sf_insert(uoo, 0, 2147483647, sf_start_idxs, sf_llbuf, &sf_lltop, &sf_entry_ct);
	    }
	  }
	  if (sf_mask & (1LLU << umm)) {
	    load_bim_sf_insert(umm, 0, unn, sf_start_idxs, sf_llbuf, &sf_lltop, &sf_entry_ct);
	  }
	} else {
	  if (ukk > unn) {
            umm = ukk;
	    ukk = unn;
	    unn = umm;
	  }
	  if (sf_mask & (1LLU << ujj)) {
	    load_bim_sf_insert(ujj, ukk, unn, sf_start_idxs, sf_llbuf, &sf_lltop, &sf_entry_ct);
	  }
	}
	uii += 2;
      } else {
	if (sf_mask & (1LLU << ujj)) {
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
    ujj = species_max_code[species];
    ukk = 0;
    for (uii = 0; uii <= ujj; uii++) {
      if (sf_start_idxs[uii] == 1) {
	sf_mask &= ~(1LLU << uii);
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
    chrom_info_ptr->chrom_mask = sf_mask;
    wkspace_reset(wkspace_mark);
  }
  if (!feof(*bimfile_ptr)) {
    goto load_bim_ret_READ_FAIL;
  }
  if (!unfiltered_marker_ct) {
    logprint("Error: No markers in .bim file.");
    goto load_bim_ret_INVALID_FORMAT;
  }
  if (from_slen || to_slen) {
    if (from_slen && (from_chrom == MAX_POSSIBLE_CHROM)) {
      sprintf(logbuf, "Error: --from marker '%s' not found.\n", markername_from);
      goto load_bim_ret_INVALID_FORMAT_2;
    }
    if (to_slen && (to_chrom == MAX_POSSIBLE_CHROM)) {
      sprintf(logbuf, "Error: --to marker '%s' not found.\n", markername_to);
      goto load_bim_ret_INVALID_FORMAT_2;
    }
    if (marker_pos_start > marker_pos_end) {
      ii = marker_pos_start;
      marker_pos_start = marker_pos_end;
      marker_pos_end = ii;
    }
  }
  if (snp_slen) {
    if (snp_chrom == MAX_POSSIBLE_CHROM) {
      sprintf(logbuf, "Error: --snp marker '%s' not found.\n", markername_snp);
      goto load_bim_ret_INVALID_FORMAT_2;
    }
    if (snp_window_size > snp_pos) {
      marker_pos_start = 0;
    } else {
      marker_pos_start = snp_pos - snp_window_size;
    }
    if (snp_window_size > (2147483647 - snp_pos)) {
      marker_pos_end = 2147483647;
    } else {
      marker_pos_end = snp_pos + snp_window_size;
    }
  }

  *unfiltered_marker_ct_ptr = unfiltered_marker_ct;
  *max_marker_id_len_ptr = max_marker_id_len;
  *plink_maxsnp_ptr = plink_maxsnp;
  rewind(*bimfile_ptr);
  ii = (unfiltered_marker_ct + (BITCT - 1)) / BITCT;

  // unfiltered_marker_ct can be very large, so use wkspace for all allocations
  // that are a multiple of it

  // permanent stack allocation #1: marker_exclude
  // permanent stack allocation #2: set_allele_freqs
  if (wkspace_alloc_ul_checked(marker_exclude_ptr, ii * sizeof(intptr_t)) ||
      wkspace_alloc_d_checked(set_allele_freqs_ptr, unfiltered_marker_ct * sizeof(double))) {
    goto load_bim_ret_NOMEM;
  }
  fill_ulong_zero(*marker_exclude_ptr, ii);
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
  if (max_marker_allele_len > 1) {
    *max_marker_allele_len_ptr = ++max_marker_allele_len;
  }
  if (marker_alleles_needed) {
    if (wkspace_alloc_c_checked(marker_alleles_ptr, unfiltered_marker_ct * sizeof(char) * 2 * max_marker_allele_len)) {
      goto load_bim_ret_NOMEM;
    }
    marker_alleles = *marker_alleles_ptr;
    memset(marker_alleles, 0, unfiltered_marker_ct * sizeof(char) * 2 * max_marker_allele_len);
  }
  if (wkspace_alloc_c_checked(marker_ids_ptr, unfiltered_marker_ct * max_marker_id_len)) {
    goto load_bim_ret_NOMEM;
  }

  // second pass: actually load stuff
  for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
    if (get_next_noncomment(*bimfile_ptr, &bufptr)) {
      goto load_bim_ret_READ_FAIL;
    }
    jj = marker_code(species, bufptr);
    if (jj != last_chrom) {
      if (last_chrom != -1) {
	chrom_info_ptr->chrom_end[last_chrom] = marker_uidx;
      }
      if (jj < last_chrom) {
	if (!(*map_is_unsorted_ptr)) {
	  *map_is_unsorted_ptr = 1;
	}
      }
      last_chrom = jj;
      if (loaded_chrom_mask & (1LLU << jj)) {
	if (calculation_type != CALC_MAKE_BED) {
	  logprint("Error: .bim file is unsorted.  Use --make-bed by itself to remedy this.\n");
	  goto load_bim_ret_INVALID_FORMAT;
	}
	*map_is_unsorted_ptr = 2;
      } else {
	loaded_chrom_mask |= 1LLU << jj;
	chrom_info_ptr->chrom_start[jj] = marker_uidx;
	chrom_info_ptr->chrom_file_order[++chroms_encountered_m1] = jj;
	chrom_info_ptr->chrom_file_order_marker_idx[chroms_encountered_m1] = marker_uidx;
      }
      last_pos = 0;
    }

    if (!(chrom_info_ptr->chrom_mask & (1LLU << jj))) {
      set_bit(*marker_exclude_ptr, marker_uidx, marker_exclude_ct_ptr);
    } else {
      bufptr = next_item(bufptr);
      if (no_more_items_kns(bufptr)) {
	goto load_bim_ret_INVALID_FORMAT_5;
	return RET_INVALID_FORMAT;
      }
      read_next_terminate(&((*marker_ids_ptr)[marker_uidx * max_marker_id_len]), bufptr);
      bufptr = next_item_mult(bufptr, mcm2);
      if (no_more_items_kns(bufptr)) {
	goto load_bim_ret_INVALID_FORMAT_5;
      }
      if (*bufptr == '-') {
	// negative marker positions now have same effect in .bim as .map
	set_bit(*marker_exclude_ptr, marker_uidx, marker_exclude_ct_ptr);
      } else {
	cur_pos = atoi(bufptr);
	if (cur_pos < last_pos) {
	  *map_is_unsorted_ptr = 2;
	} else {
	  last_pos = cur_pos;
	}
	if ((sf_ct && sf_out_of_range((uint32_t)cur_pos, (uint32_t)jj, sf_start_idxs, sf_pos)) || ((marker_pos_start != -1) && ((cur_pos < marker_pos_start) || (cur_pos > marker_pos_end)))) {
	  set_bit(*marker_exclude_ptr, marker_uidx, marker_exclude_ct_ptr);
	} else {
	  if (marker_pos_needed && jj) {
	    (*marker_pos_ptr)[marker_uidx] = cur_pos;
	  }
	  if (marker_alleles_needed) {
	    bufptr = next_item(bufptr);
	    bufptr2 = next_item(bufptr);
	    if (max_marker_allele_len == 1) {
	      marker_alleles[marker_uidx * 2] = *bufptr;
	      marker_alleles[marker_uidx * 2 + 1] = *bufptr2;
	    } else {
	      uii = strlen_se(bufptr);
	      ulii = marker_uidx * 2 * max_marker_allele_len;
	      memcpy(&(marker_alleles[ulii]), bufptr, uii);
	      // no need to append \0 since was zeroed out earlier
	      uii = strlen_se(bufptr2);
	      memcpy(&(marker_alleles[ulii + max_marker_allele_len]), bufptr2, uii);
	    }
	  }
	}
      }
    }
  }
  if (unfiltered_marker_ct == *marker_exclude_ct_ptr) {
    logprint("Error: All markers excluded.\n");
    goto load_bim_ret_INVALID_FORMAT;
  }
  chrom_info_ptr->chrom_mask &= loaded_chrom_mask;
  chrom_info_ptr->chrom_end[last_chrom] = marker_uidx;
  chrom_info_ptr->chrom_ct = ++chroms_encountered_m1;
  chrom_info_ptr->chrom_file_order_marker_idx[chroms_encountered_m1] = marker_uidx;
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
    logprint("Error: Improperly formatted .bim file.\n");
    retval = RET_INVALID_FORMAT;
    break;
  load_bim_ret_INVALID_FORMAT_4:
    logprint("Error: --from and --to markers are not on the same chromosome.\n");
    retval = RET_INVALID_FORMAT;
    break;
  load_bim_ret_INVALID_FORMAT_3:
    jj = strlen_se(bufptr);
    bufptr[jj] = '\0';
    sprintf(logbuf, "Error: Duplicate marker ID '%s' in .bim file.\n", bufptr);
  load_bim_ret_INVALID_FORMAT_2:
    logprintb();
  load_bim_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
  }
  free_cond(sf_pos);
  return retval;
}

void sort_marker_chrom_pos(int64_t* ll_buf, uintptr_t marker_ct, int32_t* pos_buf, uint32_t* chrom_start, uint32_t* chrom_id, uint32_t* chrom_ct_ptr) {
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
  ll_buf[0] = ((int64_t)uii) | (((int64_t)pos_buf[uii]) << 32);
  for (marker_idx = 1; marker_idx < marker_ct; marker_idx++) {
    if ((ll_buf[marker_idx] >> 32) != cur_chrom) {
      cur_chrom = ll_buf[marker_idx] >> 32;
      chrom_start[++chrom_ct] = marker_idx;
      chrom_id[chrom_ct] = cur_chrom;
    }
    uii = (uint32_t)ll_buf[marker_idx];
    ll_buf[marker_idx] = ((int64_t)uii) | (((int64_t)pos_buf[uii]) << 32);
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

int32_t sort_and_write_bim(int32_t** map_reverse_ptr, FILE* mapfile, int32_t map_cols, FILE** bimfile_ptr, char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, uintptr_t max_marker_id_len, int32_t species, char* marker_alleles, uintptr_t max_marker_allele_len, uintptr_t* marker_reverse, int32_t compact_map_reverse) {
  int32_t tmp_map = (bimfile_ptr == NULL);
  int64_t* ll_buf = NULL;
  int32_t retval = 0;
  char zstr[2];
  FILE* map_outfile;
  char* marker_ids;
  double* gd_vals;
  int32_t* pos_buf;
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
  char cc;
  char cc2;
  char* aptr1;
  char* aptr2;
  // There can be a LOT of markers (one existing collaborator has a file with
  // ~35 million), so speeding up the sorting step over just calling
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
  if (wkspace_alloc_i_checked(map_reverse_ptr, (compact_map_reverse? marker_ct : unfiltered_marker_ct) * sizeof(int32_t)) ||
      wkspace_alloc_ll_checked(&ll_buf, marker_ct * sizeof(int64_t)) ||
      wkspace_alloc_c_checked(&marker_ids, marker_ct * max_marker_id_len) ||
      wkspace_alloc_d_checked(&gd_vals, marker_ct * sizeof(double)) ||
      wkspace_alloc_i_checked(&pos_buf, marker_ct * sizeof(int32_t)) ||
      wkspace_alloc_ui_checked(&unpack_map, marker_ct * sizeof(int32_t))) {
    goto sort_and_write_bim_ret_NOMEM;
  }
  rewind(mapfile);
  marker_idx = 0;
  for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
    if (get_next_noncomment(mapfile, &bufptr)) {
      goto sort_and_write_bim_ret_READ_FAIL;
    }
    if (is_set(marker_exclude, marker_uidx)) {
      continue;
    }
    ll_buf[marker_idx] = (((int64_t)marker_code(species, bufptr)) << 32) + marker_idx;
    bufptr = next_item(bufptr);
    uii = strlen_se(bufptr);
    memcpy(&(marker_ids[marker_idx * max_marker_id_len]), bufptr, uii);
    marker_ids[marker_idx * max_marker_id_len + uii] = '\0';
    bufptr = next_item(bufptr);
    if (map_cols == 4) {
      gd_vals[marker_idx] = atof(bufptr);
      bufptr = next_item(bufptr);
    } else {
      gd_vals[marker_idx] = 0.0;
    }
    unpack_map[marker_idx] = marker_uidx;
    pos_buf[marker_idx++] = atoi(bufptr);
  }
  sort_marker_chrom_pos(ll_buf, marker_ct, pos_buf, chrom_start, chrom_id, &chrom_ct);

  if (tmp_map) {
    strcpy(outname_end, ".map.tmp");
    bimfile_ptr = &map_outfile;
  } else {
    strcpy(outname_end, ".bim");
  }
  if (fopen_checked(bimfile_ptr, outname, "w")) {
    goto sort_and_write_bim_ret_OPEN_FAIL;
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
      if (tmp_map) {
        if (fprintf(*bimfile_ptr, "%u\t%s\t%g\t%u\n", cur_chrom, &(marker_ids[marker_idx2 * max_marker_id_len]), gd_vals[marker_idx2], (uint32_t)(ll_buf[marker_idx] >> 32)) < 0) {
	  goto sort_and_write_bim_ret_WRITE_FAIL;
        }
      } else {
	if (max_marker_allele_len == 1) {
	  if ((!marker_reverse) || (!is_set(marker_reverse, marker_uidx))) {
	    cc = marker_alleles[2 * marker_uidx];
	    cc2 = marker_alleles[2 * marker_uidx + 1];
	  } else {
	    cc = marker_alleles[2 * marker_uidx + 1];
	    cc2 = marker_alleles[2 * marker_uidx];
	  }
	  if (!cc) {
	    cc = '0';
	  }
	  if (!cc2) {
	    cc2 = '0';
	  }
	  if (fprintf(*bimfile_ptr, "%u\t%s\t%g\t%u\t%c\t%c\n", cur_chrom, &(marker_ids[marker_idx2 * max_marker_id_len]), gd_vals[marker_idx2], (uint32_t)(ll_buf[marker_idx] >> 32), cc, cc2) < 0) {
	    goto sort_and_write_bim_ret_WRITE_FAIL;
	  }
	} else {
	  if ((!marker_reverse) || (!is_set(marker_reverse, marker_uidx))) {
	    aptr1 = &(marker_alleles[2 * marker_uidx * max_marker_allele_len]);
	    aptr2 = &(marker_alleles[(2 * marker_uidx + 1) * max_marker_allele_len]);
	  } else {
	    aptr1 = &(marker_alleles[(2 * marker_uidx + 1) * max_marker_allele_len]);
	    aptr2 = &(marker_alleles[2 * marker_uidx * max_marker_allele_len]);
	  }
	  if (!(*aptr1)) {
	    aptr1 = zstr;
	  }
	  if (!(*aptr2)) {
	    aptr2 = zstr;
	  }
	  if (fprintf(*bimfile_ptr, "%u\t%s\t%g\t%u\t%s\t%s\n", cur_chrom, &(marker_ids[marker_idx2 * max_marker_id_len]), gd_vals[marker_idx2], (uint32_t)(ll_buf[marker_idx] >> 32), aptr1, aptr2) < 0) {
	    goto sort_and_write_bim_ret_WRITE_FAIL;
	  }
	}
      }
      (*map_reverse_ptr)[compact_map_reverse? marker_idx2 : marker_uidx] = marker_idx;
    }
  }
  if (tmp_map) {
    if (fclose_null(bimfile_ptr)) {
      goto sort_and_write_bim_ret_WRITE_FAIL;
    }
  }
  while (0) {
  sort_and_write_bim_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  sort_and_write_bim_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  sort_and_write_bim_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  sort_and_write_bim_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
  }
  if (ll_buf) {
    wkspace_reset((unsigned char*)ll_buf);
  }
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

void fill_bmap_hap(unsigned char* bmap_hap, uint32_t ct_mod4) {
  unsigned char imap[4];
  unsigned char* bmap = bmap_hap;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  uint32_t unn;
  uint32_t uoo;
  imap[0] = 0;
  imap[1] = 1;
  imap[2] = 1;
  imap[3] = 3;
  for (uii = 0; uii < 4; uii++) {
    umm = imap[uii] * 64;
    for (ujj = 0; ujj < 4; ujj++) {
      unn = umm + imap[ujj] * 16;
      for (ukk = 0; ukk < 4; ukk++) {
	uoo = unn + imap[ukk] * 4;
	*bmap++ = uoo;
	*bmap++ = uoo + 1;
	*bmap++ = uoo + 1;
	*bmap++ = uoo + 3;
      }
    }
  }
  imap[0] = 3;
  imap[3] = 0;
  for (uii = 0; uii < 4; uii++) {
    umm = imap[uii] * 64;
    for (ujj = 0; ujj < 4; ujj++) {
      unn = umm + imap[ujj] * 16;
      for (ukk = 0; ukk < 4; ukk++) {
	uoo = unn + imap[ukk] * 4;
	*bmap++ = uoo + 3;
	*bmap++ = uoo + 1;
	*bmap++ = uoo + 1;
	*bmap++ = uoo;
      }
    }
  }
  if (ct_mod4) {
    switch (ct_mod4) {
    case 1:
      bmap[0] = 3;
      bmap[1] = 1;
      bmap[2] = 1;
      bmap[3] = 0;
      return;
    case 2:
      for (uii = 0; uii < 4; uii++) {
	ukk = imap[uii] * 4;
	*bmap++ = ukk + 3;
	*bmap++ = ukk + 1;
	*bmap++ = ukk + 1;
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
	  *bmap++ = umm + 1;
	  *bmap++ = umm;
	}
      }
    }
  }
}

int32_t make_bed(FILE* bedfile, int32_t bed_offset, FILE* bimfile, int32_t map_cols, FILE** bedoutfile_ptr, FILE** famoutfile_ptr, FILE** bimoutfile_ptr, char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_alleles, uintptr_t max_marker_allele_len, uintptr_t* marker_reverse, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_ct, char* person_ids, uintptr_t max_person_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* sex_nm, uintptr_t* sex_male, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, double missing_phenod, char* output_missing_pheno, uintptr_t max_marker_id_len, int32_t map_is_unsorted, uint32_t* indiv_sort_map, uint32_t set_hh_missing, Chrom_info* chrom_info_ptr) {
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  uintptr_t indiv_ct4 = (indiv_ct + 3) / 4;
  unsigned char* wkspace_mark = wkspace_base;
  int32_t affection = (pheno_c != NULL);
  int32_t species = chrom_info_ptr->species;
  uintptr_t* is_male_arr;
  uintptr_t marker_uidx;
  uintptr_t marker_idx;
  uintptr_t indiv_uidx;
  uintptr_t indiv_uidx2;
  uintptr_t indiv_idx;
  uint32_t ii_mod4;
  uint32_t slen;
  uint32_t slen2;
  uint32_t chrom_fo_idx;
  uint32_t chrom_end;
  int32_t chrom_idx;
  unsigned char* loadbuf;
  unsigned char* writebuf;
  char* bufptr;
  char* sptr;
  char* sptr2;
  char* cptr;
  uint32_t clen;
  unsigned char cc;
  uint32_t pct;
  uint32_t loop_end;
  int32_t* map_reverse;
  unsigned char* writeptr;
  unsigned char bmap_raw[576];
  unsigned char bmap_hap[576];
  unsigned char imap[4];
  unsigned char imap_m[4];
  unsigned char* bmap_rawh;
  unsigned char* bmap;
  unsigned char* bmap2;
  int32_t retval;
  if (wkspace_alloc_uc_checked(&loadbuf, unfiltered_indiv_ct4)) {
    return RET_NOMEM;
  }
  fill_bmap_raw(bmap_raw, indiv_ct & 3);

  strcpy(outname_end, ".bed");
  if (fopen_checked(bedoutfile_ptr, outname, "wb")) {
    return RET_OPEN_FAIL;
  }

  marker_uidx = 0;
  marker_idx = 0;
  sprintf(logbuf, "--make-bed to %s + .bim + .fam... ", outname);
  logprintb();
  fputs("0%", stdout);
  if (fwrite_checked("l\x1b\x01", 3, *bedoutfile_ptr)) {
    return RET_WRITE_FAIL;
  }
  fflush(stdout);
  if (map_is_unsorted) {
    if (set_hh_missing) {
      logprint("\nError: --set-hh-missing cannot currently be used on an unsorted .bim file.  Use\n--make-bed without --set-hh-missing to sort by position first; then run\n--make-bed + --set-hh-missing on the new fileset.\n");
      return RET_CALC_NOT_YET_SUPPORTED;
    }
    retval = sort_and_write_bim(&map_reverse, bimfile, map_cols, bimoutfile_ptr, outname, outname_end, unfiltered_marker_ct, marker_exclude, marker_ct, max_marker_id_len, species, marker_alleles, max_marker_allele_len, marker_reverse, 0);
    if (retval) {
      return retval;
    }

    if (wkspace_alloc_uc_checked(&writebuf, marker_ct * indiv_ct4)) {
      logprint("\nError: Insufficient memory for current --make-bed implementation.  Try raising\nthe --memory value for now.\n");
      return RET_CALC_NOT_YET_SUPPORTED;
    } else {
      if (fseeko(bedfile, bed_offset, SEEK_SET)) {
	return RET_READ_FAIL;
      }
      for (pct = 1; pct <= 100; pct++) {
	loop_end = ((uint64_t)pct * marker_ct) / 100LU;
	for (; marker_idx < loop_end; marker_idx++) {
	  if (is_set(marker_exclude, marker_uidx)) {
	    marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx + 1);
	    if (fseeko(bedfile, bed_offset + (uint64_t)marker_uidx * unfiltered_indiv_ct4, SEEK_SET)) {
	      return RET_READ_FAIL;
	    }
	  }
	  if (is_set(marker_reverse, marker_uidx)) {
	    bmap = &(bmap_raw[256]);
	    bmap2 = &(bmap_raw[512]);
	  } else {
	    bmap = bmap_raw;
	    bmap2 = bmap_raw;
	  }
	  if (fread(loadbuf, 1, unfiltered_indiv_ct4, bedfile) < unfiltered_indiv_ct4) {
	    return RET_READ_FAIL;
	  }
	  writeptr = &(writebuf[indiv_ct4 * map_reverse[marker_uidx]]);
	  indiv_uidx = 0;
	  cc = 0;
	  if (indiv_sort_map) {
	    for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
	      do {
		indiv_uidx2 = indiv_sort_map[indiv_uidx++];
	      } while (is_set(indiv_exclude, indiv_uidx2));
	      ii_mod4 = indiv_idx & 3;
	      cc |= (((loadbuf[indiv_uidx2 / 4] >> ((indiv_uidx2 & 3) * 2)) & 3) << (ii_mod4 * 2));
	      if (ii_mod4 == 3) {
		writeptr[indiv_idx / 4] = bmap[cc];
		cc = 0;
	      }
	    }
	  } else {
	    for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
	      indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx);
	      ii_mod4 = indiv_idx & 3;
	      cc |= (((loadbuf[indiv_uidx / 4] >> ((indiv_uidx & 3) * 2)) & 3) << (ii_mod4 * 2));
	      if (ii_mod4 == 3) {
		writeptr[indiv_idx / 4] = bmap[cc];
		cc = 0;
	      }
	      indiv_uidx++;
	    }
	  }
	  if (indiv_ct & 3) {
	    writeptr[indiv_ct4 - 1] = bmap2[cc];
	  }
	  marker_uidx++;
	}
	if (pct < 100) {
	  if (pct > 10) {
	    putchar('\b');
	  }
	  printf("\b\b%u%%", pct);
	  fflush(stdout);
	}
      }
      if (fwrite_checked(writebuf, marker_ct * indiv_ct4, *bedoutfile_ptr)) {
	return RET_WRITE_FAIL;
      }
    }
  } else {
    if (wkspace_alloc_uc_checked(&writebuf, indiv_ct4)) {
      return RET_NOMEM;
    }
    fill_bmap_hap(bmap_hap, indiv_ct & 3);
    if (set_hh_missing) {
      if (wkspace_alloc_ul_checked(&is_male_arr, unfiltered_indiv_ctl * sizeof(intptr_t))) {
	return RET_NOMEM;
      }
      for (indiv_idx = 0; indiv_idx < unfiltered_indiv_ctl; indiv_idx++) {
	is_male_arr[indiv_idx] = sex_nm[indiv_idx] & sex_male[indiv_idx];
      }
      pct = 0;
      loop_end = marker_ct / 100;
      imap[1] = 1;
      imap[2] = 2;
      imap_m[1] = 1;
      imap_m[2] = 1;
      for (chrom_fo_idx = 0; chrom_fo_idx < chrom_info_ptr->chrom_ct; chrom_fo_idx++) {
	marker_uidx = chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx];
	chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[chrom_fo_idx + 1];
	if (fseeko(bedfile, bed_offset + (uint64_t)marker_uidx * unfiltered_indiv_ct4, SEEK_SET)) {
	  return RET_READ_FAIL;
	}
	chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	if (chrom_idx == species_x_code[species]) {
	  for (; marker_uidx < chrom_end; marker_uidx++) {
	    if (is_set(marker_exclude, marker_uidx)) {
	      marker_uidx = next_non_set(marker_exclude, marker_uidx + 1, chrom_end);
	      if (marker_uidx == chrom_end) {
		break;
	      }
	      if (fseeko(bedfile, bed_offset + (uint64_t)marker_uidx * unfiltered_indiv_ct4, SEEK_SET)) {
		return RET_READ_FAIL;
	      }
	    }
	    if (is_set(marker_reverse, marker_uidx)) {
	      imap[0] = 3;
	      imap[3] = 0;
	      imap_m[0] = 3;
	      imap_m[3] = 0;
	    } else {
	      imap[0] = 0;
	      imap[3] = 3;
	      imap_m[0] = 0;
	      imap_m[3] = 3;
	    }
	    if (fread(loadbuf, 1, unfiltered_indiv_ct4, bedfile) < unfiltered_indiv_ct4) {
	      return RET_READ_FAIL;
	    }
	    indiv_uidx = 0;
	    cc = 0;
	    if (indiv_sort_map) {
	      for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
		do {
		  indiv_uidx2 = indiv_sort_map[indiv_uidx++];
		} while (is_set(indiv_exclude, indiv_uidx2));
		ii_mod4 = indiv_idx & 3;
		if (is_set(is_male_arr, indiv_uidx2)) {
		  cc |= (imap_m[(loadbuf[indiv_uidx2 / 4] >> ((indiv_uidx2 & 3) * 2)) & 3] << (ii_mod4 * 2));
		} else {
		  cc |= (imap[(loadbuf[indiv_uidx2 / 4] >> ((indiv_uidx2 & 3) * 2)) & 3] << (ii_mod4 * 2));
		}
		if (ii_mod4 == 3) {
		  writebuf[indiv_idx / 4] = cc;
		  cc = 0;
		}
	      }
	    } else {
	      for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
		indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx);
		ii_mod4 = indiv_idx & 3;
		if (is_set(is_male_arr, indiv_uidx)) {
		  cc |= (imap_m[(loadbuf[indiv_uidx / 4] >> ((indiv_uidx & 3) * 2)) & 3] << (ii_mod4 * 2));
		} else {
		  cc |= (imap[(loadbuf[indiv_uidx / 4] >> ((indiv_uidx & 3) * 2)) & 3] << (ii_mod4 * 2));
		}
		if (ii_mod4 == 3) {
		  writebuf[indiv_idx / 4] = cc;
		  cc = 0;
		}
		indiv_uidx++;
	      }
	    }
	    if (indiv_ct & 3) {
	      writebuf[indiv_ct4 - 1] = cc;
	    }
	    if (fwrite_checked(writebuf, indiv_ct4, *bedoutfile_ptr)) {
	      return RET_WRITE_FAIL;
	    }
	    if ((++marker_idx) >= loop_end) {
	      if (marker_idx < marker_ct) {
		if (pct >= 10) {
		  putchar('\b');
		}
		pct = (marker_idx * 100LLU) / marker_ct;
		printf("\b\b%u%%", pct);
		fflush(stdout);
		loop_end = ((uint64_t)(pct + 1LLU) * marker_ct) / 100;
	      } else {
		chrom_fo_idx = chrom_info_ptr->chrom_ct - 1;
		break;
	      }
	    }
	  }
	} else {
	  if ((species_haploid_mask[species] >> chrom_idx) & 1LLU) {
	    bmap_rawh = bmap_hap;
	  } else {
	    bmap_rawh = bmap_raw;
	  }
	  for (; marker_uidx < chrom_end; marker_uidx++) {
	    if (is_set(marker_exclude, marker_uidx)) {
	      marker_uidx = next_non_set(marker_exclude, marker_uidx + 1, chrom_end);
	      if (marker_uidx == chrom_end) {
		break;
	      }
	      if (fseeko(bedfile, bed_offset + (uint64_t)marker_uidx * unfiltered_indiv_ct4, SEEK_SET)) {
		return RET_READ_FAIL;
	      }
	    }
	    if (is_set(marker_reverse, marker_uidx)) {
	      bmap = &(bmap_rawh[256]);
	      bmap2 = &(bmap_rawh[512]);
	    } else {
	      bmap = bmap_rawh;
	      bmap2 = bmap_rawh;
	    }
	    if (fread(loadbuf, 1, unfiltered_indiv_ct4, bedfile) < unfiltered_indiv_ct4) {
	      return RET_READ_FAIL;
	    }
	    indiv_uidx = 0;
	    cc = 0;
	    if (indiv_sort_map) {
	      for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
		do {
		  indiv_uidx2 = indiv_sort_map[indiv_uidx++];
		} while (is_set(indiv_exclude, indiv_uidx2));
		ii_mod4 = indiv_idx & 3;
		cc |= (((loadbuf[indiv_uidx2 / 4] >> ((indiv_uidx2 & 3) * 2)) & 3) << (ii_mod4 * 2));
		if (ii_mod4 == 3) {
		  writebuf[indiv_idx / 4] = bmap[cc];
		  cc = 0;
		}
	      }
	    } else {
	      for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
		indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx);
		ii_mod4 = indiv_idx & 3;
		cc |= (((loadbuf[indiv_uidx / 4] >> ((indiv_uidx & 3) * 2)) & 3) << (ii_mod4 * 2));
		if (ii_mod4 == 3) {
		  writebuf[indiv_idx / 4] = bmap[cc];
		  cc = 0;
		}
		indiv_uidx++;
	      }
	    }
	    if (indiv_ct & 3) {
	      writebuf[indiv_ct4 - 1] = bmap2[cc];
	    }
	    if (fwrite_checked(writebuf, indiv_ct4, *bedoutfile_ptr)) {
	      return RET_WRITE_FAIL;
	    }
	    if ((++marker_idx) >= loop_end) {
	      if (marker_idx < marker_ct) {
		if (pct >= 10) {
		  putchar('\b');
		}
		pct = (marker_idx * 100LLU) / marker_ct;
		printf("\b\b%u%%", pct);
		fflush(stdout);
		loop_end = ((uint64_t)(pct + 1LLU) * marker_ct) / 100;
	      } else {
		chrom_fo_idx = chrom_info_ptr->chrom_ct - 1;
		break;
	      }
	    }
	  }
	}
      }
    } else {
      if (fseeko(bedfile, bed_offset, SEEK_SET)) {
	return RET_READ_FAIL;
      }
      for (pct = 1; pct <= 100; pct++) {
	loop_end = ((uint64_t)pct * marker_ct) / 100;
	for (; marker_idx < loop_end; marker_idx++) {
	  if (is_set(marker_exclude, marker_uidx)) {
	    marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx + 1);
	    if (fseeko(bedfile, bed_offset + (uint64_t)marker_uidx * unfiltered_indiv_ct4, SEEK_SET)) {
	      return RET_READ_FAIL;
	    }
	  }
	  if (is_set(marker_reverse, marker_uidx)) {
	    bmap = &(bmap_raw[256]);
	    bmap2 = &(bmap_raw[512]);
	  } else {
	    bmap = bmap_raw;
	    bmap2 = bmap_raw;
	  }
	  if (fread(loadbuf, 1, unfiltered_indiv_ct4, bedfile) < unfiltered_indiv_ct4) {
	    return RET_READ_FAIL;
	  }
	  indiv_uidx = 0;
	  cc = 0;
	  if (indiv_sort_map) {
	    for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
	      do {
		indiv_uidx2 = indiv_sort_map[indiv_uidx++];
	      } while (is_set(indiv_exclude, indiv_uidx2));
	      ii_mod4 = indiv_idx & 3;
	      cc |= (((loadbuf[indiv_uidx2 / 4] >> ((indiv_uidx2 & 3) * 2)) & 3) << (ii_mod4 * 2));
	      if (ii_mod4 == 3) {
		writebuf[indiv_idx / 4] = bmap[cc];
		cc = 0;
	      }
	    }
	  } else {
	    for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
	      indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx);
	      ii_mod4 = indiv_idx & 3;
	      cc |= (((loadbuf[indiv_uidx / 4] >> ((indiv_uidx & 3) * 2)) & 3) << (ii_mod4 * 2));
	      if (ii_mod4 == 3) {
		writebuf[indiv_idx / 4] = bmap[cc];
		cc = 0;
	      }
	      indiv_uidx++;
	    }
	  }
	  if (indiv_ct & 3) {
	    writebuf[indiv_ct4 - 1] = bmap2[cc];
	  }
	  if (fwrite_checked(writebuf, indiv_ct4, *bedoutfile_ptr)) {
	    return RET_WRITE_FAIL;
	  }
	  marker_uidx++;
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
  }
  fclose_null(bedoutfile_ptr);

  strcpy(outname_end, ".fam");
  if (fopen_checked(famoutfile_ptr, outname, "w")) {
    return RET_OPEN_FAIL;
  }
  indiv_uidx = 0;
  indiv_uidx2 = 0;
  for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
    if (indiv_sort_map) {
      do {
	indiv_uidx = indiv_sort_map[indiv_uidx2++];
      } while (is_set(indiv_exclude, indiv_uidx));
    } else {
      indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx);
    }
    bufptr = tbuf;
    cptr = &(person_ids[indiv_uidx * max_person_id_len]);
    clen = strlen(cptr);
    memcpy(bufptr, cptr, clen);
    bufptr[strlen_se(cptr)] = ' ';
    bufptr += clen;
    *bufptr++ = ' ';
    if (paternal_ids) {
      cptr = &(paternal_ids[indiv_uidx * max_paternal_id_len]);
      clen = strlen(cptr);
      memcpy(bufptr, cptr, clen);
      bufptr += clen;
    } else {
      *bufptr++ = '0';
    }
    *bufptr++ = ' ';
    if (maternal_ids) {
      cptr = &(maternal_ids[indiv_uidx * max_maternal_id_len]);
      clen = strlen(cptr);
      memcpy(bufptr, cptr, clen);
      bufptr += clen;
    } else {
      *bufptr++ = '0';
    }
    *bufptr++ = ' ';
    *bufptr++ = sexchar(sex_nm, sex_male, indiv_uidx);
    *bufptr++ = ' ';
    if (!is_set(pheno_nm, indiv_uidx)) {
      sprintf(bufptr, "%s\n", output_missing_pheno);
    } else if (affection) {
      bufptr[0] = is_set(pheno_c, indiv_uidx)? '2' : '1';
      bufptr[1] = '\n';
      bufptr[2] = '\0';
    } else {
      sprintf(bufptr, "%g\n", pheno_d[indiv_uidx]);
    }
    if (fputs(tbuf, *famoutfile_ptr) == EOF) {
      return RET_WRITE_FAIL;
    }
    if (!indiv_sort_map) {
      indiv_uidx++;
    }
  }
  fclose_null(famoutfile_ptr);

  if (!map_is_unsorted) {
    strcpy(outname_end, ".bim");
    if (fopen_checked(bimoutfile_ptr, outname, "w")) {
      return RET_OPEN_FAIL;
    }
    rewind(bimfile);
    for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
      if (get_next_noncomment(bimfile, &bufptr)) {
	return RET_READ_FAIL;
      }
      if (is_set(marker_exclude, marker_uidx)) {
	continue;
      }
      cptr = next_item_mult(bufptr, map_cols);
      if (max_marker_allele_len == 1) {
	if (is_set(marker_reverse, marker_uidx)) {
	  *cptr = marker_alleles[2 * marker_uidx + 1];
	  cptr = next_item(cptr);
	  *cptr = marker_alleles[2 * marker_uidx];
	} else {
	  *cptr = marker_alleles[2 * marker_uidx];
	  cptr = next_item(cptr);
	  *cptr = marker_alleles[2 * marker_uidx + 1];
	}
      } else {
	sptr = &(marker_alleles[2 * marker_uidx * max_marker_allele_len]);
	sptr2 = &(sptr[max_marker_allele_len]);
	slen = strlen_se(sptr);
	slen2 = strlen_se(sptr2);
	if (is_set(marker_reverse, marker_uidx)) {
	  memcpy(cptr, sptr2, slen2);
	  cptr[slen2] = '\t';
	  memcpy(&(cptr[slen2 + 1]), sptr, slen);
	} else {
	  memcpy(cptr, sptr, slen);
	  cptr[slen] = '\t';
	  memcpy(&(cptr[slen + 1]), sptr2, slen2);
	}
	cptr[slen + slen2 + 1] = '\n';
	cptr[slen + slen2 + 2] = '\0';
      }
      if (fputs(bufptr, *bimoutfile_ptr) == EOF) {
	return RET_WRITE_FAIL;
      }
    }
    fclose_null(bimoutfile_ptr);
  }

  fputs("\b\b\b", stdout);
  logprint("done.\n");
  wkspace_reset(wkspace_mark);
  return 0;
}

const char errstr_fam_format[] = "Error: Improperly formatted .fam file.\n";

int32_t load_fam(FILE* famfile, uint32_t buflen, int32_t fam_col_1, int32_t fam_col_34, int32_t fam_col_5, int32_t fam_col_6, int32_t true_fam_col_6, int32_t missing_pheno, int32_t missing_pheno_len, int32_t affection_01, uintptr_t* unfiltered_indiv_ct_ptr, char** person_ids_ptr, uintptr_t* max_person_id_len_ptr, char** paternal_ids_ptr, uintptr_t* max_paternal_id_len_ptr, char** maternal_ids_ptr, uintptr_t* max_maternal_id_len_ptr, uintptr_t** sex_nm_ptr, uintptr_t** sex_male_ptr, int32_t* affection_ptr, uintptr_t** pheno_nm_ptr, uintptr_t** pheno_c_ptr, double** pheno_d_ptr, uintptr_t** founder_info_ptr, uintptr_t** indiv_exclude_ptr) {
  char* bufptr0;
  char* bufptr;
  uintptr_t unfiltered_indiv_ct = 0;
  uintptr_t max_person_id_len = 4;
  uintptr_t max_paternal_id_len = 2;
  uintptr_t max_maternal_id_len = 2;
  int32_t affection = 1;
  uint64_t last_tell = 0;
  uint32_t new_buflen = 0;
  unsigned char* wkspace_mark = wkspace_base;
  double missing_phenod = (double)missing_pheno;
  char* linebuf;
  char* person_ids;
  char* paternal_ids = NULL;
  char* maternal_ids = NULL;
  uintptr_t* pheno_c = NULL;
  double* pheno_d = NULL;
  char cc;
  uint64_t* line_locs;
  uint64_t* tmp_ullp;
  uint32_t max_people;
  uintptr_t tmp_len;
  uintptr_t tmp_len2;
  uintptr_t indiv_uidx;
  uintptr_t unfiltered_indiv_ctl;
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
      if (fam_col_1) {
	bufptr = next_item(bufptr0);
      } else {
	bufptr = bufptr0;
      }
      tmp_len = strlen_se(bufptr0) + strlen_se(bufptr) + 2;
      if (tmp_len > max_person_id_len) {
	max_person_id_len = tmp_len;
      }
      if (fam_col_34) {
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
      if (fam_col_5) {
	bufptr = next_item(bufptr);
      }
      if (fam_col_6) {
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
  if (wkspace_alloc_c_checked(person_ids_ptr, unfiltered_indiv_ct * max_person_id_len)) {
    return RET_NOMEM;
  }
  person_ids = *person_ids_ptr;
  if (fam_col_34) {
    // could make this conditional, but memory footprint of this is typically
    // negligible
    if (wkspace_alloc_c_checked(paternal_ids_ptr, unfiltered_indiv_ct * max_paternal_id_len)) {
      return RET_NOMEM;
    }
    paternal_ids = *paternal_ids_ptr;
    if (wkspace_alloc_c_checked(maternal_ids_ptr, unfiltered_indiv_ct * max_maternal_id_len)) {
      return RET_NOMEM;
    }
    maternal_ids = *maternal_ids_ptr;
  }
  unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  if (wkspace_alloc_ul_checked(sex_nm_ptr, unfiltered_indiv_ctl * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(sex_male_ptr, unfiltered_indiv_ctl * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(founder_info_ptr, unfiltered_indiv_ctl * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(indiv_exclude_ptr, unfiltered_indiv_ctl * sizeof(intptr_t)) ||
      wkspace_alloc_ul_checked(pheno_nm_ptr, unfiltered_indiv_ctl * sizeof(intptr_t))) {
    return RET_NOMEM;
  }

  if (fam_col_6) {
    if (affection) {
      pheno_c = (uintptr_t*)malloc(unfiltered_indiv_ctl * sizeof(intptr_t));
      if (!pheno_c) {
	return RET_NOMEM;
      }
      *pheno_c_ptr = pheno_c;
    } else {
      pheno_d = (double*)malloc(unfiltered_indiv_ct * sizeof(double));
      if (!pheno_d) {
	return RET_NOMEM;
      }
      *pheno_d_ptr = pheno_d;
    }
  }
  wkspace_mark = wkspace_base;
  if (new_buflen) {
    buflen = new_buflen;
  }
  if (wkspace_alloc_c_checked(&linebuf, buflen)) {
    return RET_NOMEM;
  }
  if (wkspace_alloc_ull_checked(&tmp_ullp, unfiltered_indiv_ct * sizeof(int64_t))) {
    return RET_NOMEM;
  }
  for (ii = unfiltered_indiv_ct - 1; ii >= 0; ii--) {
    tmp_ullp[ii] = line_locs[ii];
  }
  line_locs = tmp_ullp;
  if (fam_col_34) {
    fill_ulong_zero(*founder_info_ptr, unfiltered_indiv_ctl);
  } else {
    fill_ulong_one(*founder_info_ptr, unfiltered_indiv_ctl);
  }
  fill_ulong_zero(*sex_nm_ptr, unfiltered_indiv_ctl);
  fill_ulong_zero(*sex_male_ptr, unfiltered_indiv_ctl);
  fill_ulong_zero(*indiv_exclude_ptr, unfiltered_indiv_ctl);
  fill_ulong_zero(*pheno_nm_ptr, unfiltered_indiv_ctl);

  // ----- .fam read, second pass -----
  for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ct; indiv_uidx++) {
    if (fseeko(famfile, line_locs[indiv_uidx], SEEK_SET)) {
      return RET_READ_FAIL;
    }
    if (fgets(linebuf, buflen, famfile) == NULL) {
      return RET_READ_FAIL;
    }
    if (fam_col_1) {
      bufptr = next_item(linebuf);
    } else {
      bufptr = linebuf;
    }
    tmp_len = strlen_se(linebuf);
    tmp_len2 = strlen_se(bufptr);
    memcpy(&(person_ids[indiv_uidx * max_person_id_len]), linebuf, tmp_len);
    person_ids[indiv_uidx * max_person_id_len + tmp_len] = '\t';
    memcpy(&(person_ids[indiv_uidx * max_person_id_len + tmp_len + 1]), bufptr, tmp_len2);
    person_ids[indiv_uidx * max_person_id_len + tmp_len + tmp_len2 + 1] = '\0';
    if (fam_col_34) {
      bufptr = next_item(bufptr);
      cc = *bufptr;
      tmp_len = strlen_se(bufptr);
      memcpy(&(paternal_ids[indiv_uidx * max_paternal_id_len]), bufptr, tmp_len);
      paternal_ids[indiv_uidx * max_paternal_id_len + tmp_len] = '\0';
      bufptr = next_item(bufptr);
      tmp_len2 = strlen_se(bufptr);
      memcpy(&(maternal_ids[indiv_uidx * max_maternal_id_len]), bufptr, tmp_len2);
      maternal_ids[indiv_uidx * max_maternal_id_len + tmp_len] = '\0';
      if ((tmp_len == 1) && (tmp_len2 == 1) && (cc == '0') && (*bufptr == '0')) {
	set_bit_noct(*founder_info_ptr, indiv_uidx);
      }
    }
    if (fam_col_5) {
      bufptr = next_item(bufptr);
      if (strlen_se(bufptr) == 1) {
	if (*bufptr == '1') {
	  set_bit_noct(*sex_nm_ptr, indiv_uidx);
	  set_bit_noct(*sex_male_ptr, indiv_uidx);
	} else if (*bufptr == '2') {
	  set_bit_noct(*sex_nm_ptr, indiv_uidx);
	}
      }
    }
    if (fam_col_6) {
      bufptr = next_item(bufptr);
      if (affection) {
	if (!is_missing(bufptr, missing_pheno, missing_pheno_len, affection_01)) {
	  set_bit_noct(*pheno_nm_ptr, indiv_uidx);
	  if (*bufptr == (affection_01? '1' : '2')) {
	    set_bit_noct(pheno_c, indiv_uidx);
	  } else {
	    clear_bit_noct(pheno_c, indiv_uidx);
	  }
	}
      } else {
	if (sscanf(bufptr, "%lg", &(pheno_d[indiv_uidx])) == 1) {
	  if (pheno_d[indiv_uidx] != missing_phenod) {
	    set_bit_noct(*pheno_nm_ptr, indiv_uidx);
	  }
	}
      }
    }
    if (true_fam_col_6 && (!fam_col_6)) {
      bufptr = next_item(bufptr);
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
int32_t check_gd_col(FILE* bimfile, char* tbuf, uint32_t is_binary, uint32_t* gd_col_ptr) {
  char* bufptr;
  while (fgets(tbuf, MAXLINELEN - 5, bimfile)) {
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
  int32_t ii;
  for (ii = 0; ii < 4; ii++) {
    if (marker_alleles[ii] == '\0') {
      marker_alleles[ii] = cc;
      marker_allele_cts[ii] = 1;
      return 0;
    } else if (marker_alleles[ii] == cc) {
      marker_allele_cts[ii] += 1;
      return 0;
    }
  }
  return -1;
}

typedef struct ll_str_struct {
  struct ll_str_struct* next;
  char ss[];
} Ll_str;

typedef struct ll_str_fixed_struct {
  struct ll_str_struct* next;
#ifdef __LP64__
  char ss[8];
#else
  char ss[12];
#endif
} Ll_str_fixed;

static inline Ll_str* top_alloc_llstr(uintptr_t* topsize_ptr, uint32_t size) {
  return (Ll_str*)top_alloc(topsize_ptr, size + sizeof(Ll_str));
}

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
	cur_allele_name_start = llptr->ss;
      }
      memcpy(cur_allele_name_start, allele_name, an_len);
      cur_allele_name_start[an_len] = '\0';
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
	chars_left = ((2147483632 - sizeof(intptr_t)) - ((uintptr_t)(&(cur_allele_name_start[slen + 1]) - allele_list_start->ss))) & 15;
	if (chars_left > an_len) {
	  cur_allele_name_start[slen] = '\t';
	  memcpy(&(cur_allele_name_start[slen + 1]), allele_name, an_len);
	  cur_allele_name_start[slen + an_len + 1] = '\0';
	} else {
	  llptr = top_alloc_llstr(topsize_ptr, an_len + 1);
	  if (!llptr) {
	    return RET_NOMEM;
	  }
	  allele_list_start->next = llptr;
	  cur_allele_name_start = llptr->ss;
	  memcpy(cur_allele_name_start, allele_name, an_len);
	  cur_allele_name_start[an_len + 1] = '\0';
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

int32_t ped_to_bed_multichar_allele(uintptr_t max_marker_allele_len, FILE** pedfile_ptr, FILE** outfile_ptr, char* outname, char* outname_end, FILE** mapfile_ptr, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_alleles_f, uint32_t map_is_unsorted, int32_t fam_col_1, int32_t fam_col_34, int32_t fam_col_5, int32_t fam_col_6, uint32_t ped_col_skip, uint32_t gd_col, uint32_t* map_reverse, int64_t ped_size) {
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
  uint32_t pass_ct;
  uintptr_t indiv_ct4;
  char* loadbuf;
  char* col1_ptr;
  char* col2_ptr;
  char* bufptr;
  char* bufptr2;
  char* bufptr3;
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
  uint64_t ullii;
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
  rewind(*pedfile_ptr);
  fputs("Rescanning .ped file... 0%", stdout);
  fflush(stdout);
  while (1) {
    loadbuf_size = wkspace_left - topsize;
    if (loadbuf_size > 2147483584) {
      loadbuf_size = 2147483584;
    }
    loadbuf[loadbuf_size - 1] = ' ';
  ped_to_bed_multichar_allele_loop_1_start:
    if (!fgets(loadbuf, loadbuf_size, *pedfile_ptr)) {
      break;
    }
    if (!loadbuf[loadbuf_size - 1]) {
      putchar('\n');
      goto ped_to_bed_multichar_allele_ret_NOMEM;
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
    cur_slen_rdup = (cur_slen + CACHELINE) & (CACHELINE - 1);
    if (fam_col_1) {
      col2_ptr = next_item(col1_ptr);
    } else {
      col2_ptr = col1_ptr;
    }
    bufptr = next_item_mult(col2_ptr, ped_col_skip - 1);
    if (no_more_items_kns(bufptr)) {
      putchar('\n');
      sprintf(logbuf, "Error: Missing columns in .ped line: %s\n", col1_ptr);
      goto ped_to_bed_multichar_allele_ret_INVALID_FORMAT_2;
    }
    if ((bufptr - col1_ptr) > (MAXLINELEN / 2) - 4) {
      putchar('\n');
      logprint("Error: Pathologically long header item(s) in .ped file.\n");
      goto ped_to_bed_multichar_allele_ret_INVALID_FORMAT;
    }
    uii = strlen_se(col1_ptr);
    memcpy(tbuf, col1_ptr, uii);
    tbuf[uii++] = '\t';
    ujj = strlen_se(col2_ptr);
    memcpy(&(tbuf[uii]), col2_ptr, ujj);
    uii += ujj;
    tbuf[uii++] = '\t';
    bufptr2 = item_endnn(col2_ptr);
    if (fam_col_34) {
      copy_item(tbuf, &uii, &bufptr2);
      copy_item(tbuf, &uii, &bufptr2);
    } else {
      memcpy(&(tbuf[uii]), "0\t0\t", 4);
      uii += 4;
    }
    if (fam_col_5) {
      copy_item(tbuf, &uii, &bufptr2);
    } else {
      memcpy(&(tbuf[uii]), "0\t", 2);
      uii += 2;
    }
    if (fam_col_6) {
      copy_item(tbuf, &uii, &bufptr2);
      tbuf[uii - 1] = '\n';
    } else {
      memcpy(&(tbuf[uii]), "-9\n", 3);
      uii += 3;
    }
    if (fwrite_checked(tbuf, uii, *outfile_ptr)) {
      goto ped_to_bed_multichar_allele_ret_WRITE_FAIL;
    }
    wkspace_left -= cur_slen_rdup;
    marker_idx = 0;
    for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
      alen1 = strlen_se(bufptr);
      aptr1 = bufptr;
      bufptr = skip_initial_spaces(&(bufptr[alen1]));
      alen2 = strlen_se(bufptr);
      if (!alen2) {
	wkspace_left += cur_slen_rdup;
	goto ped_to_bed_multichar_allele_ret_INVALID_FORMAT_3;
      }
      aptr2 = bufptr;
      bufptr = skip_initial_spaces(&(bufptr[alen2]));
      if (is_set(marker_exclude, marker_uidx)) {
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
      if (alen1 > max_marker_allele_len) {
	max_marker_allele_len = alen1;
      }
      if (alen2 > max_marker_allele_len) {
	max_marker_allele_len = alen2;
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
    wkspace_left += cur_slen_rdup;
    if (!is_eoln_kns(*bufptr)) {
      putchar('\n');
      sprintf(logbuf, "Error: Too many entries in .ped line for indiv %" PRIuPTR ".\n", indiv_ct);
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
  fputs("\r.ped scan complete (for binary autoconversion).\n", stdout);
  if (!indiv_ct) {
    sprintf(logbuf, "Error: No %s in .ped file.\n", species_plural);
    goto ped_to_bed_multichar_allele_ret_INVALID_FORMAT_2;
  }
  if (fclose_null(outfile_ptr)) {
    goto ped_to_bed_multichar_allele_ret_WRITE_FAIL;
  }
  max_marker_allele_len++;
  ullii = marker_ct * ((uint64_t)max_marker_allele_len) * 2 * sizeof(char);
  if (topsize + ullii > ((uint64_t)wkspace_left)) {
    goto ped_to_bed_multichar_allele_ret_NOMEM;
  }
  marker_alleles_f = (char*)wkspace_alloc(ullii);
  memcpy(outname_end, ".bim", 5);
  if (fopen_checked(outfile_ptr, outname, "w")) {
    goto ped_to_bed_multichar_allele_ret_OPEN_FAIL;
  }
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
      ukk = marker_allele_cts[4 * marker_idx + 3];
      if (map_is_unsorted) {
        sprintf(logbuf, "Warning: Marker %u (post-sort/filter) %sallelic; setting rarest missing.\n", map_reverse[marker_idx] + 1, (ukk? "quad" : "tri"));
      } else {
        sprintf(logbuf, "Warning: Marker %" PRIuPTR " %sallelic; setting rarest alleles missing.\n", marker_idx + 1, (ukk? "quad" : "tri"));
      }
      get_top_two(&(marker_allele_cts[4 * marker_idx]), ukk? 4 : 3, &uii, &ujj);
      ukk = map_reverse[marker_idx];
    } else {
      uii = (marker_allele_cts[4 * marker_idx] < marker_allele_cts[4 * marker_idx + 1])? 1 : 0;
      ujj = uii ^ 1;
      ukk = marker_idx;
    }
    aptr1 = &(marker_alleles_f[2 * marker_idx * max_marker_allele_len]);
    aptr2 = &(aptr1[max_marker_allele_len]);
    copy_nse(aptr1, get_llstr((Ll_str*)(&(marker_alleles_tmp[ukk])), uii));
    copy_nse(aptr2, get_llstr((Ll_str*)(&(marker_alleles_tmp[ukk])), ujj));
    bufptr3 = &(tbuf[MAXLINELEN]);
    if (map_is_unsorted) {
      bufptr = (char*)memchr(tbuf, '\n', MAXLINELEN);
      ulii = (uintptr_t)(bufptr - tbuf);
      memcpy(bufptr3, tbuf, ulii);
      bufptr2 = &(bufptr3[ulii]);
    } else {
      uii = 0;
      copy_item(bufptr3, &uii, &bufptr);
      copy_item(bufptr3, &uii, &bufptr);
      if (gd_col) {
	copy_item(bufptr3, &uii, &bufptr);
      } else {
	memcpy(&(bufptr3[uii]), "0\t", 2);
	uii += 2;
      }
      bufptr = skip_initial_spaces(bufptr);
      ujj = strlen_se(bufptr);
      memcpy(&(bufptr3[uii]), bufptr, ujj);
      bufptr2 = &(bufptr3[uii + ujj]);
    }
    *bufptr2++ = '\t';
    alen1 = strlen_se(aptr1);
    if (alen1) {
      memcpy(bufptr2, aptr1, alen1);
      bufptr2 = &(bufptr2[alen1]);
    } else {
      *bufptr2++ = '0';
    }
    *bufptr2++ = '\t';
    alen2 = strlen_se(aptr2);
    if (alen2) {
      memcpy(bufptr2, aptr2, alen2);
      bufptr2 = &(bufptr2[alen2]);
    } else {
      *bufptr2++ = '0';
    }
    *bufptr2++ = '\n';
    if (fwrite_checked(bufptr3, bufptr2 - bufptr3, *outfile_ptr)) {
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
    sprintf(logbuf, "Performing single-pass .bed write (%" PRIuPTR " marker%s, %" PRIuPTR " %s).\n", marker_ct, (marker_ct == 1)? "" : "s", indiv_ct, species_str(indiv_ct));
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
    sprintf(logbuf, "Performing %u-pass .bed write (%u/%" PRIuPTR " marker%s/pass, %" PRIuPTR " %s).\n", pass_ct, markers_per_pass, marker_ct, (markers_per_pass == 1)? "" : "s", indiv_ct, species_str(indiv_ct));
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
	    aptr1[alen1++] = '\0';
	    aptr2 = bufptr;
	    bufptr = item_endnn(bufptr);
	    alen2 = (uintptr_t)(bufptr - aptr2);
	    bufptr = skip_initial_spaces(bufptr);
	    aptr2[alen2++] = '\0';
	    if (is_set(marker_exclude, marker_uidx)) {
	      continue;
	    }
	    ukk = map_reverse[umm++];
	    if ((ukk >= marker_start) && (ukk < marker_end)) {
	      ucc = 1;
	      if (!memcmp(aptr1, &(marker_alleles_f[(2 * ukk + 1) * max_marker_allele_len]), alen1)) {
		if (!memcmp(aptr1, aptr2, alen1)) {
		  ucc = 3;
		} else if (!memcmp(aptr2, &(marker_alleles_f[2 * ukk * max_marker_allele_len]), alen2)) {
		  ucc = 2;
		}
	      } else if (!memcmp(aptr1, &(marker_alleles_f[2 * ukk * max_marker_allele_len]), alen1)) {
		if (!memcmp(aptr1, aptr2, alen1)) {
		  ucc = 0;
		} else if (!memcmp(aptr2, &(marker_alleles_f[(2 * ukk + 1) * max_marker_allele_len]), alen2)) {
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
	    aptr1[alen1++] = '\0';
	    aptr2 = bufptr;
	    bufptr = item_endnn(bufptr);
	    alen2 = (uintptr_t)(bufptr - aptr2);
	    bufptr = skip_initial_spaces(bufptr);
	    aptr2[alen2++] = '\0';
	    if (is_set(marker_exclude, marker_uidx)) {
	      continue;
	    }
	    ucc = 1;
	    if (!memcmp(aptr1, &(marker_alleles_f[(2 * marker_idx + 1) * max_marker_allele_len]), alen1)) {
	      if (!memcmp(aptr1, aptr2, alen1)) {
		ucc = 3;
	      } else if (!memcmp(aptr2, &(marker_alleles_f[2 * marker_idx * max_marker_allele_len]), alen2)) {
		ucc = 2;
	      }
	    } else if (!memcmp(aptr1, &(marker_alleles_f[2 * marker_idx * max_marker_allele_len]), alen1)) {
	      if (!memcmp(aptr1, aptr2, alen1)) {
		ucc = 0;
	      } else if (!memcmp(aptr2, &(marker_alleles_f[(2 * marker_idx + 1) * max_marker_allele_len]), alen2)) {
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
    wkspace_left += cur_slen_rdup;
    putchar('\n');
    if (retval != RET_NOMEM) {
      sprintf(logbuf, "Error: More than 4 different alleles at marker %u%s.\n", uii + 1, map_is_unsorted? " (post-sort/filter)" : "");
      logprintb();
    }
    break;
  ped_to_bed_multichar_allele_ret_INVALID_FORMAT_4:
    wkspace_left += cur_slen_rdup;
    putchar('\n');
    sprintf(logbuf, "Error: Half-missing call in .ped file at marker %" PRIuPTR ", indiv %" PRIuPTR ".\n", marker_uidx + 1, indiv_ct + 1);
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  ped_to_bed_multichar_allele_ret_INVALID_FORMAT_3:
    putchar('\n');
    sprintf(logbuf, "Error: Not enough markers in .ped line %" PRIuPTR ".\n", indiv_ct + 1);
  ped_to_bed_multichar_allele_ret_INVALID_FORMAT_2:
    logprintb();
  ped_to_bed_multichar_allele_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
  return retval;
}

int32_t ped_to_bed(char* pedname, char* mapname, char* outname, char* outname_end, int32_t fam_col_1, int32_t fam_col_34, int32_t fam_col_5, int32_t fam_col_6, int32_t affection_01, int32_t missing_pheno, Chrom_info* chrom_info_ptr) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* mapfile = NULL;
  FILE* pedfile = NULL;
  FILE* outfile = NULL;
  uintptr_t unfiltered_marker_ct = 0;
  uintptr_t marker_exclude_ct = 0;
  uintptr_t* marker_exclude;
  uintptr_t max_marker_id_len = 0;
  uintptr_t marker_ct = 0;
  uint32_t map_is_unsorted = 0;
  int32_t last_chrom = 0;
  uint32_t last_mpos = 0;
  uintptr_t indiv_ct = 0;
  uint32_t ped_buflen = 1;
  int32_t retval = 0;
  uint32_t ped_col_skip = 1 + fam_col_1 + 2 * fam_col_34 + fam_col_5 + fam_col_6;
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
  uint32_t gd_col;
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
  char* bufptr3;
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
  if (check_gd_col(mapfile, tbuf, 0, &gd_col)) {
    sprintf(logbuf, "Error: Missing columns in .map line: %s", skip_initial_spaces(tbuf));
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
    bufptr = next_item_mult(col2_ptr, 1 + gd_col);
    if (no_more_items_kns(bufptr)) {
      sprintf(logbuf, "Error: Missing columns in .map line: %s", col1_ptr);
      goto ped_to_bed_ret_INVALID_FORMAT_2;
    }
    ii = marker_code(chrom_info_ptr->species, col1_ptr);
    if (ii == -1) {
      logprint("Error: Invalid chromosome index in .map file.\n");
      goto ped_to_bed_ret_INVALID_FORMAT;
    }
    if ((*bufptr == '-') || (!(chrom_info_ptr->chrom_mask & (1LLU << ii)))) {
      set_bit(marker_exclude, unfiltered_marker_ct, &marker_exclude_ct);
    } else {
      if ((*bufptr < '0') || (*bufptr > '9')) {
	logprint("Error: Non-numeric marker position in .map file.\n");
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
    if (unfiltered_marker_ct > 2147483647) {
      logprint("Error: Too many markers in .map file (max 2147483647).\n");
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
    logprint("Error: No markers in current analysis.\n");
    goto ped_to_bed_ret_INVALID_FORMAT;
  }
  marker_exclude = (uintptr_t*)wkspace_alloc(((unfiltered_marker_ct + (BITCT - 1)) / BITCT) * sizeof(intptr_t));

  // provisionally assume max_marker_allele_len == 1
  if (wkspace_alloc_c_checked(&marker_alleles_f, marker_ct * 2)) {
    goto ped_to_bed_ret_NOMEM;
  }
  if (map_is_unsorted) {
    retval = sort_and_write_bim((int32_t**)(&map_reverse), mapfile, 3 + gd_col, NULL, outname, outname_end, unfiltered_marker_ct, marker_exclude, marker_ct, max_marker_id_len, chrom_info_ptr->species, NULL, 1, NULL, 1);
    if (retval) {
      goto ped_to_bed_ret_1;
    }
    gd_col = 1;
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
  if (loadbuf_size > 2147483584) {
    loadbuf_size = 2147483584;
  }
  if (loadbuf_size < MAXLINELEN) {
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
      goto ped_to_bed_ret_NOMEM;
    }
    col1_ptr = skip_initial_spaces(loadbuf);
    if (is_eoln_or_comment(*col1_ptr)) {
      ulii = strlen(loadbuf) + 1;
      if (ulii > ped_buflen) {
	ped_buflen = ulii;
      }
      continue;
    }
    if (fam_col_1) {
      col2_ptr = next_item(col1_ptr);
    } else {
      col2_ptr = col1_ptr;
    }
    bufptr = next_item_mult(col2_ptr, ped_col_skip - 1);
    if (no_more_items_kns(bufptr)) {
      sprintf(logbuf, "\nError: Missing columns in .ped line: %s\n", col1_ptr);
      goto ped_to_bed_ret_INVALID_FORMAT_2;
    }
    if ((bufptr - col1_ptr) > (MAXLINELEN / 2) - 4) {
      logprint("\nError: Pathologically long header item(s) in .ped file.\n");
      goto ped_to_bed_ret_INVALID_FORMAT;
    }
    uii = strlen_se(col1_ptr);
    memcpy(tbuf, col1_ptr, uii);
    tbuf[uii++] = '\t';
    ujj = strlen_se(col2_ptr);
    memcpy(&(tbuf[uii]), col2_ptr, ujj);
    uii += ujj;
    tbuf[uii++] = '\t';
    bufptr2 = item_endnn(col2_ptr);
    if (fam_col_34) {
      copy_item(tbuf, &uii, &bufptr2);
      copy_item(tbuf, &uii, &bufptr2);
    } else {
      memcpy(&(tbuf[uii]), "0\t0\t", 4);
      uii += 4;
    }
    if (fam_col_5) {
      copy_item(tbuf, &uii, &bufptr2);
    } else {
      memcpy(&(tbuf[uii]), "0\t", 2);
      uii += 2;
    }
    if (fam_col_6) {
      copy_item(tbuf, &uii, &bufptr2);
      tbuf[uii - 1] = '\n';
    } else {
      memcpy(&(tbuf[uii]), "-9\n", 3);
      uii += 3;
    }
    if (fwrite_checked(tbuf, uii, outfile)) {
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
      if (is_set(marker_exclude, marker_uidx)) {
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
      sprintf(logbuf, "Possibly irregular .ped line.  Restarting scan, assuming multichar alleles.\n");
      logprintb();
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
      sprintf(logbuf, "\nError: No %s in .ped file.\n", species_plural);
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
	  sprintf(logbuf, "Warning: Marker %u (post-sort/filter) %sallelic; setting rarest missing.\n", map_reverse[marker_idx] + 1, (cc? "quad" : "tri"));
	} else {
	  sprintf(logbuf, "Warning: Marker %" PRIuPTR " %sallelic; setting rarest alleles missing.\n", marker_idx + 1, (cc? "quad" : "tri"));
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
      bufptr3 = &(tbuf[MAXLINELEN]);
      if (map_is_unsorted) {
	bufptr = (char*)memchr(tbuf, '\n', MAXLINELEN);
	ulii = (uintptr_t)(bufptr - tbuf);
	memcpy(bufptr3, tbuf, ulii);
	bufptr2 = &(bufptr3[ulii]);
      } else {
	uii = 0;
	copy_item(bufptr3, &uii, &bufptr);
	copy_item(bufptr3, &uii, &bufptr);
	if (gd_col) {
	  copy_item(bufptr3, &uii, &bufptr);
	} else {
	  memcpy(&(bufptr3[uii]), "0\t", 2);
	  uii += 2;
	}
	bufptr = skip_initial_spaces(bufptr);
	ujj = strlen_se(bufptr);
	memcpy(&(bufptr3[uii]), bufptr, ujj);
	bufptr2 = &(bufptr3[uii + ujj]);
      }
      *bufptr2++ = '\t';
      *bufptr2++ = cc;
      *bufptr2++ = '\t';
      *bufptr2++ = cc2;
      *bufptr2++ = '\n';
      if (fwrite_checked(bufptr3, bufptr2 - bufptr3, outfile)) {
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
      sprintf(logbuf, "Performing single-pass .bed write (%" PRIuPTR " marker%s, %" PRIuPTR " %s).\n", marker_ct, (marker_ct == 1)? "" : "s", indiv_ct, species_str(indiv_ct));
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
      sprintf(logbuf, "Performing %u-pass .bed write (%u/%" PRIuPTR " marker%s/pass, %" PRIuPTR " %s).\n", pass_ct, markers_per_pass, marker_ct, (markers_per_pass == 1)? "" : "s", indiv_ct, species_str(indiv_ct));
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
	      if (is_set(marker_exclude, marker_uidx)) {
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
	      if (is_set(marker_exclude, marker_uidx)) {
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
    retval = ped_to_bed_multichar_allele(max_marker_allele_len, &pedfile, &outfile, outname, outname_end, &mapfile, unfiltered_marker_ct, marker_exclude, marker_ct, marker_alleles_f, map_is_unsorted, fam_col_1, fam_col_34, fam_col_5, fam_col_6, ped_col_skip, gd_col, map_reverse, ped_size);
    if (retval) {
      goto ped_to_bed_ret_1;
    }
  }

  if (fclose_null(&outfile)) {
    goto ped_to_bed_ret_WRITE_FAIL_2;
  }
  putchar('\r');
  sprintf(logbuf, "Binary fileset written to %s + .bim + .fam.\n", outname);
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
    sprintf(logbuf, "Error: Not enough markers in .ped line %" PRIuPTR ".\n", indiv_ct + 1);
  ped_to_bed_ret_INVALID_FORMAT_2:
    logprintb();
  ped_to_bed_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 ped_to_bed_ret_1:
  fclose_cond(pedfile);
  fclose_cond(mapfile);
  fclose_cond(outfile);
  wkspace_reset(wkspace_mark);
  return retval;
}

// 256 - 32 = 224
char one_char_strs[448];

int32_t heapstr_alloc(char** new_sptr, char* new_str, uint32_t slen) {
  char* sptr;
  if (slen == 1) {
    *new_sptr = &(one_char_strs[2 * (((unsigned char)(*new_str)) - ' ')]);
    return 0;
  } else {
    sptr = (char*)malloc(slen + 1);
    if (!sptr) {
      return -1;
    }
    memcpy(sptr, new_str, slen + 1);
    *new_sptr = sptr;
    return 0;
  }
}

int32_t lgen_to_bed(char* lgen_namebuf, char* outname, char* outname_end, int32_t missing_pheno, int32_t affection_01, uint32_t lgen_modifier, char* lgen_reference_fname, Chrom_info* chrom_info_ptr) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* infile = NULL;
  FILE* outfile = NULL;
  char* name_end = (char*)memchr(lgen_namebuf, 0, FNAMESIZE);
  uint32_t lgen_allele_count = lgen_modifier & LGEN_ALLELE_COUNT;
  int32_t map_cols = 3;
  uintptr_t* marker_exclude = NULL;
  uintptr_t marker_exclude_ct = 0;
  uintptr_t max_marker_id_len = 0;
  uintptr_t unfiltered_marker_ct = 0;
  uintptr_t marker_ct = 0;
  char** marker_alleles = NULL;
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
  int32_t affection = 0;
  uintptr_t* founder_info = NULL;
  uintptr_t* indiv_exclude = NULL;
  int32_t map_is_unsorted = 0;
  int32_t duplicate_fail = 0;
  uintptr_t* pheno_nm = NULL;
  uintptr_t* pheno_c = NULL;
  double* pheno_d = NULL;
  char* sorted_marker_ids;
  int32_t* marker_id_map;
  int32_t* map_reverse;
  char* sorted_indiv_ids;
  int32_t* indiv_id_map;
  unsigned char* writebuf;
  uintptr_t indiv_ct4;
  uintptr_t max_person_id_len;
  uintptr_t marker_idx;
  unsigned char ucc;
  unsigned char* ucptr;
  unsigned char bmap_short[320];
  unsigned char* bmap2;
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
  cptr = one_char_strs;
  for (uii = 0; uii < 224; uii++) {
    *cptr++ = (char)(uii + ' ');
    *cptr++ = '\0';
  }

  memcpy(name_end, ".map", 5);
  retval = load_map(&infile, lgen_namebuf, &map_cols, &unfiltered_marker_ct, &marker_exclude_ct, &max_marker_id_len, &marker_exclude, &marker_ids, chrom_info_ptr, &marker_pos, &map_is_unsorted);
  if (retval) {
    goto lgen_to_bed_ret_1;
  }
  marker_ct = unfiltered_marker_ct - marker_exclude_ct;
  duplicate_fail = 1;
  retval = sort_item_ids(&sorted_marker_ids, &marker_id_map, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_ids, max_marker_id_len, strcmp_deref, &duplicate_fail);
  if (retval) {
    goto lgen_to_bed_ret_1;
  }
  if (map_is_unsorted) {
    // Writes a temporary .map which is read later, and then deleted.
    retval = sort_and_write_bim(&map_reverse, infile, map_cols, NULL, outname, outname_end, unfiltered_marker_ct, marker_exclude, marker_ct, max_marker_id_len, chrom_info_ptr->species, NULL, 1, NULL, 0);
    if (retval) {
      goto lgen_to_bed_ret_1;
    }
    for (uii = 0; uii < marker_ct; uii++) {
      marker_id_map[uii] = map_reverse[marker_id_map[uii]];
    }
  }
  // collapse
  if (wkspace_alloc_i_checked(&indiv_id_map, unfiltered_marker_ct * sizeof(int32_t))) {
    goto lgen_to_bed_ret_NOMEM;
  }
  marker_idx = 0;
  for (uii = 0; uii < marker_ct; uii++) {
    marker_idx = next_non_set_unsafe(marker_exclude, marker_idx);
    indiv_id_map[marker_idx++] = uii;
  }
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
  retval = load_fam(infile, MAXLINELEN, 1, 1, 1, 1, 1, missing_pheno, intlen(missing_pheno), affection_01, &indiv_ct, &person_ids, &max_person_id_len, &paternal_ids, &max_paternal_id_len, &maternal_ids, &max_maternal_id_len, &sex_nm, &sex_male, &affection, &pheno_nm, &pheno_c, &pheno_d, &founder_info, &indiv_exclude);
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
  marker_alleles = (char**)wkspace_alloc(2 * marker_ct * sizeof(char*));
  if (!marker_alleles) {
    goto lgen_to_bed_ret_NOMEM;
  }
  memset(marker_alleles, 0, 2 * marker_ct * sizeof(char*));
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
  if (lgen_modifier & LGEN_REFERENCE) {
    if (fopen_checked(&infile, lgen_reference_fname, "r")) {
      goto lgen_to_bed_ret_OPEN_FAIL;
    }
    while (fgets(tbuf, MAXLINELEN, infile)) {
      cptr = skip_initial_spaces(tbuf);
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
	memcpy(id_buf, cptr, a1len);
	id_buf[a1len] = '\0';
	ii = bsearch_str(id_buf, marker_ids, max_marker_id_len, 0, marker_ct - 1);
	if (ii != -1) {
	  marker_idx = marker_id_map[ii];
	  if (marker_alleles[2 * marker_idx + 1]) {
	    goto lgen_to_bed_ret_INVALID_FORMAT_4;
	  }
	  sptr = item_end(a1ptr);
	  a2ptr = skip_initial_spaces(sptr);
	  a1len = (uintptr_t)(sptr - a1ptr);
	  a1ptr[a1len] = '\0';
	  if (heapstr_alloc(&(marker_alleles[2 * marker_idx + 1]), a1ptr, a1len)) {
	    goto lgen_to_bed_ret_NOMEM;
	  }
	  if (no_more_items_kns(a2ptr)) {
	    if (lgen_allele_count) {
	      a1ptr[a1len++] = 'v';
	      a1ptr[a1len] = '\0';
	      if (heapstr_alloc(&(marker_alleles[2 * marker_idx]), a1ptr, a1len)) {
		goto lgen_to_bed_ret_NOMEM;
	      }
	    }
	  } else {
	    a2len = strlen_se(a2ptr);
	    a2ptr[a2len] = '\0';
	    if (heapstr_alloc(&(marker_alleles[2 * marker_idx]), a2ptr, a2len)) {
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
  // PLINK reports an error whenever there are 3+ alleles at one locus, so
  // backwards compatibility does not mandate that we worry about that case.
  // Thus we just use the obvious one-pass load, and save proper handling of
  // triallelic sites, etc. for the future VCF engine.
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
    while (fgets(tbuf, MAXLINELEN, infile)) {
      cptr = skip_initial_spaces(tbuf);
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
      if (no_more_items_kns(a2ptr)) {
	goto lgen_to_bed_ret_INVALID_FORMAT;
      }
      a2len = strlen_se(a2ptr);
      ii = bsearch_fam_indiv(id_buf, sorted_indiv_ids, max_person_id_len, indiv_ct, cptr, cptr2);
      if (ii == -1) {
	goto lgen_to_bed_ret_INVALID_FORMAT_2;
      }
      indiv_idx = indiv_id_map[ii];
      ulii = (uintptr_t)(cptr4 - cptr3);
      memcpy(id_buf, cptr3, ulii);
      id_buf[ulii] = '\0';
      ii = bsearch_str(id_buf, marker_ids, max_marker_id_len, 0, marker_ct - 1);
      if (ii != -1) {
	marker_idx = marker_id_map[ii];
	sptr = marker_alleles[2 * marker_idx + 1]; // existing A2
	a1ptr[a1len] = '\0';
	a2ptr[a2len] = '\0';
	if ((*a1ptr == '0') && (a1len == 1)) {
	  if ((*a2ptr == '0') && (a2len == 1)) {
	    uii = 1;
	  } else {
	    goto lgen_to_bed_ret_INVALID_FORMAT_5;
	  }
	} else if ((*a2ptr == '0') && (a2len == 1)) {
	  goto lgen_to_bed_ret_INVALID_FORMAT_5;
        } else {
          if (!sptr) {
	    if (heapstr_alloc(&(marker_alleles[2 * marker_idx + 1]), a1ptr, a1len)) {
	      goto lgen_to_bed_ret_NOMEM;
	    }
	    if (!strcmp(a1ptr, a2ptr)) {
	      uii = 2;
	    } else {
	      uii = 1;
	      if (heapstr_alloc(&(marker_alleles[2 * marker_idx]), a2ptr, a2len)) {
		goto lgen_to_bed_ret_NOMEM;
	      }
	    }
	  } else {
	    sptr2 = marker_alleles[2 * marker_idx];
	    if (!sptr2) {
	      if (!strcmp(a1ptr, sptr)) {
		if (!strcmp(a2ptr, sptr)) {
		  uii = 2;
		} else {
		  uii = 1;
		  if (heapstr_alloc(&(marker_alleles[2 * marker_idx]), a2ptr, a2len)) {
		    goto lgen_to_bed_ret_NOMEM;
		  }
		}
	      } else {
		if (heapstr_alloc(&(marker_alleles[2 * marker_idx]), a1ptr, a1len)) {
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
    while (fgets(tbuf, MAXLINELEN, infile)) {
      cptr = skip_initial_spaces(tbuf);
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
      indiv_idx = indiv_id_map[ii];
      ulii = (uintptr_t)(cptr4 - cptr3);
      memcpy(id_buf, cptr3, ulii);
      id_buf[ulii] = '\0';
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
      cptr = marker_alleles[uii * 2];
      marker_alleles[uii * 2] = marker_alleles[uii * 2 + 1];
      marker_alleles[uii * 2 + 1] = cptr;
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
    if (!marker_alleles[ujj]) {
      marker_alleles[ujj] = &(one_char_strs[32]);
    }
  }
  uii = 0;
  marker_idx = 0;
  while (fgets(tbuf, MAXLINELEN, infile)) {
    if (is_eoln_or_comment(*(skip_initial_spaces(tbuf)))) {
      continue;
    }
    if (is_set(marker_exclude, uii)) {
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
    sprintf(cptr, "\t%s\t%s\n", marker_alleles[marker_idx * 2], marker_alleles[marker_idx * 2 + 1]);
    ulii = strlen(cptr) + (uintptr_t)(cptr - tbuf);
    if (fwrite_checked(tbuf, ulii, outfile)) {
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
  if (strcmp(lgen_namebuf, outname)) {
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
  sprintf(logbuf, "%s + .bim + .fam written.\n", outname);
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
    retval = RET_INVALID_FORMAT;
    break;
  lgen_to_bed_ret_INVALID_FORMAT_3:
    sprintf(logbuf, "Error: Marker %s in .lgen file has 3+ different alleles.\n", id_buf);
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  lgen_to_bed_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  }
 lgen_to_bed_ret_1:
  if (marker_alleles) {
    ma_end = &(marker_alleles[2 * marker_ct]);
    cptr = &(one_char_strs[448]);
    while (marker_alleles < ma_end) {
      sptr = *marker_alleles++;
      if (sptr && (((uintptr_t)sptr < (uintptr_t)one_char_strs) || ((uintptr_t)sptr >= (uintptr_t)cptr))) {
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

inline uint32_t update_alleles_and_cts(char* alleles, uint32_t* allele_cts, char cc) {
  if (cc == alleles[0]) {
    allele_cts[0] += 1;
    return 0;
  } else if (!alleles[0]) {
    alleles[0] = cc;
    allele_cts[0] = 1;
    return 0;
  } else if (cc == alleles[1]) {
    allele_cts[1] += 1;
    return 1;
  } else if (!alleles[1]) {
    alleles[1] = cc;
    allele_cts[1] = 1;
    return 1;
  } else if (cc == alleles[2]) {
    allele_cts[2] += 1;
    return 2;
  } else if (!alleles[2]) {
    alleles[2] = cc;
    allele_cts[2] = 1;
    return 2;
  } else if (cc == alleles[3]) {
    allele_cts[3] += 1;
    return 3;
  } else if (!alleles[3]) {
    alleles[3] = cc;
    allele_cts[3] += 1;
    return 3;
  } else {
    return 4;
  }
}

void transposed_to_bed_print_pct(uint32_t pct) {
  printf("Processing .tped file... %u%%", pct);
  fflush(stdout);
}

int32_t transposed_to_bed(char* tpedname, char* tfamname, char* outname, char* outname_end, char missing_geno, Chrom_info* chrom_info_ptr) {
  FILE* infile = NULL;
  FILE* bimfile = NULL;
  FILE* outfile = NULL;
  uintptr_t indiv_ct = 0;
  uint32_t no_extra_cols = 1;
  int32_t retval = 0;
  uint32_t pct = 0;
  uintptr_t marker_idx = 0;
  int64_t last_mapval = 0;
  int32_t map_is_unsorted = 0;
  uintptr_t max_marker_id_len = 0;
  uintptr_t indiv_ct4;
  uintptr_t indiv_idx;
  uintptr_t ulii;
  uint32_t uii;
  uint32_t ujj;
  int32_t ii;
  // It's unnecessary, but we choose to include nice handling of triallelic and
  // quadallelic sites in this routine.
  char alleles[4];
  char orig_alleles[4];
  uint32_t allele_cts[4];
  unsigned char* writebuf;
  unsigned char* prewritebuf;
  unsigned char writemap[17];
  unsigned char ucc;
  uint32_t max_load;
  char* loadbuf;
  char* cptr;
  char* cptr2;
  char* cptr3;
  char* cptr4;
  unsigned char* ucptr;
  unsigned char* ucptr2;
  char cc;
  char cc2;
  int64_t tped_size;
  int64_t tped_next_thresh;
  int64_t cur_mapval;
  int64_t* mapvals;
  uintptr_t marker_ct;
  uintptr_t marker_uidx;
  uint32_t* map_reverse;
  int64_t* ll_buf;
  int32_t* pos_buf;
  char* marker_ids;
  uint32_t chrom_start[MAX_POSSIBLE_CHROM + 1];
  uint32_t chrom_id[MAX_POSSIBLE_CHROM];
  uint32_t cur_chrom;
  uint32_t chrom_ct;
  double* gd_vals;
  char* marker_alleles;

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
      goto transposed_to_bed_ret_WRITE_FAIL;
    }
    indiv_ct++;
  }
  if (!feof(infile)) {
    goto transposed_to_bed_ret_READ_FAIL;
  }
  if (!indiv_ct) {
    sprintf(logbuf, "Error: No %s in .tfam file.\n", species_plural);
    goto transposed_to_bed_ret_INVALID_FORMAT_5;
  }
  indiv_ct4 = (indiv_ct + 3) / 4;
  fclose_null(&infile);
  fclose_null(&outfile);
  if (fopen_checked(&infile, tpedname, "r")) {
    goto transposed_to_bed_ret_OPEN_FAIL;
  }
  if (fseeko(infile, 0, SEEK_END)) {
    goto transposed_to_bed_ret_READ_FAIL;
  }
  tped_size = ftello(infile);
  rewind(infile);
  tped_next_thresh = tped_size / 100;

  memcpy(outname_end, ".bim.tmp", 9);
  if (fopen_checked(&bimfile, outname, "w")) {
    goto transposed_to_bed_ret_OPEN_FAIL;
  }
  memcpy(outname_end, ".bed.tmp", 9);
  if (fopen_checked(&outfile, outname, "wb")) {
    goto transposed_to_bed_ret_OPEN_FAIL;
  }
  if (wkspace_alloc_uc_checked(&writebuf, indiv_ct4)) {
    goto transposed_to_bed_ret_NOMEM;
  }
  if (wkspace_alloc_uc_checked(&prewritebuf, indiv_ct)) {
    goto transposed_to_bed_ret_NOMEM;
  }
  mapvals = (int64_t*)wkspace_base;
  loadbuf = (char*)(&wkspace_base[sizeof(int64_t)]);
  if (wkspace_left > 2147483592) {
    max_load = 2147483584;
  } else {
    max_load = wkspace_left - 8;
  }
  writemap[16] = 1;
  if (fwrite_checked("l\x1b\x01", 3, outfile)) {
    goto transposed_to_bed_ret_WRITE_FAIL;
  }
  while (fgets(loadbuf, max_load, infile)) {
    cptr = skip_initial_spaces(loadbuf);
    if (is_eoln_kns(*cptr)) {
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
    ii = marker_code(chrom_info_ptr->species, cptr);
    if ((ii == -1) || (!(chrom_info_ptr->chrom_mask & (1LLU << ii)))) {
      continue;
    }
    uii = strlen_se(cptr2) + 1;
    if (uii > max_marker_id_len) {
      max_marker_id_len = uii;
    }
    if (*cptr3 == '-') {
      continue;
    }
    memset(alleles, 0, 4);
    allele_cts[0] = 0;
    allele_cts[1] = 0;
    allele_cts[2] = 0;
    allele_cts[3] = 0;
    cur_mapval = (((int64_t)ii) << 32) | atoi(cptr3);
    mapvals[marker_idx++] = cur_mapval;
    if (last_mapval > cur_mapval) {
      map_is_unsorted = 1;
    } else {
      last_mapval = cur_mapval;
    }
    for (uii = 0; uii < 3; uii++) {
      cptr2 = item_endnn(cptr);
      *cptr2++ = '\t';
      if (fwrite_checked(cptr, (cptr2 - cptr), bimfile)) {
	goto transposed_to_bed_ret_WRITE_FAIL;
      }
      cptr = skip_initial_spaces(cptr2);
    }
    cptr2 = item_endnn(cptr);
    *cptr2++ = '\t';
    if (fwrite_checked(cptr, (cptr2 - cptr), bimfile)) {
      goto transposed_to_bed_ret_WRITE_FAIL;
    }
    cptr2 = cptr4;
    for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
      cc = *cptr2++;
      if (!cc) {
	goto transposed_to_bed_ret_INVALID_FORMAT_3;
      }
      cptr2 = skip_initial_spaces(cptr2);
      cc2 = *cptr2++;
      if (!cc2) {
	goto transposed_to_bed_ret_INVALID_FORMAT_3;
      }
      cptr2 = skip_initial_spaces(cptr2);
      if (cc == missing_geno) {
	if (cc2 != missing_geno) {
	  goto transposed_to_bed_ret_INVALID_FORMAT_4;
	}
	prewritebuf[indiv_idx] = 16;
      } else if (cc2 == missing_geno) {
	goto transposed_to_bed_ret_INVALID_FORMAT_4;
      } else {
	uii = update_alleles_and_cts(alleles, allele_cts, cc);
	ujj = update_alleles_and_cts(alleles, allele_cts, cc2);
	if ((uii == 4) || (ujj == 4)) {
	  cptr[strlen_se(cptr)] = '\0';
	  putchar('\r');
	  sprintf(logbuf, "Error: More than four alleles at marker %s.\n", cptr);
	  logprintb();
	  goto transposed_to_bed_ret_INVALID_FORMAT_2;
	}
        prewritebuf[indiv_idx] = uii * 4 + ujj;
      }
    }

    if (no_extra_cols && (!is_eoln_kns(*cptr2))) {
      no_extra_cols = 0;
      putchar('\r');
      logprint("Note: Extra columns in .tped file.  Ignoring.\n");
      transposed_to_bed_print_pct(pct);
    }
    if (max_load == 1 + strlen(cptr2) + (cptr2 - loadbuf)) {
      goto transposed_to_bed_ret_NOMEM;
    }

    memcpy(orig_alleles, alleles, 4);
    for (uii = 1; uii < 4; uii++) {
      ujj = allele_cts[uii];
      if (allele_cts[uii - 1] < ujj) {
	cc = alleles[uii];
	ii = uii;
	do {
	  ii--;
	  alleles[ii + 1] = alleles[ii];
	  allele_cts[ii + 1] = allele_cts[ii];
	} while ((ii > 0) && (allele_cts[ii - 1] < ujj));
	alleles[ii] = cc;
	allele_cts[ii] = ujj;
      }
    }
    if (allele_cts[2]) {
      ulii = strlen_se(cptr);
      cc = cptr[ulii];
      cptr[ulii] = '\0';
      putchar('\r');
      sprintf(logbuf, "Note: Marker %s is %sallelic.  Setting rarest alleles to missing.\n", cptr, allele_cts[3]? "quad" : "tri");
      logprintb();
      transposed_to_bed_print_pct(pct);
      cptr[ulii] = cc;
    }
    for (uii = 0; uii < 4; uii++) {
      cc = orig_alleles[uii];
      ucptr = &(writemap[4 * uii]);
      if (cc == '\0') {
	memset(ucptr, 1, 4);
      } else if (cc == alleles[0]) {
        for (ujj = 0; ujj < 4; ujj++) {
	  cc = orig_alleles[ujj];
	  if (cc == '\0') {
	    *ucptr++ = 1;
	  } else if (cc == alleles[0]) {
	    *ucptr++ = 3;
	  } else if (cc == alleles[1]) {
	    *ucptr++ = 2;
	  } else {
	    *ucptr++ = 1;
	  }
	}
      } else if (cc == alleles[1]) {
	for (ujj = 0; ujj < 4; ujj++) {
	  cc = orig_alleles[ujj];
	  if (cc == '\0') {
	    *ucptr++ = 1;
	  } else if (cc == alleles[0]) {
	    *ucptr++ = 2;
	  } else if (cc == alleles[1]) {
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
    if (fwrite_checked(writebuf, indiv_ct4, outfile)) {
      goto transposed_to_bed_ret_WRITE_FAIL;
    }
    if (fprintf(bimfile, "%c\t%c\n", alleles[1]? alleles[1] : '0', alleles[0]? alleles[0] : '0') < 0) {
      goto transposed_to_bed_ret_WRITE_FAIL;
    }
    loadbuf = &(loadbuf[8]);
    max_load -= 8;
  }
  if (!feof(infile)) {
    goto transposed_to_bed_ret_READ_FAIL;
  }
  fclose_null(&infile);
  if (fclose_null(&bimfile)) {
    goto transposed_to_bed_ret_WRITE_FAIL;
  }
  if (fclose_null(&outfile)) {
    goto transposed_to_bed_ret_WRITE_FAIL;
  }

  if (map_is_unsorted) {
    marker_ct = marker_idx;

    // mapvals is already positioned here, this just updates wkspace_base and
    // wkspace_left
    mapvals = (int64_t*)wkspace_alloc(marker_ct * sizeof(int64_t));

    if (wkspace_alloc_ll_checked(&ll_buf, marker_ct * sizeof(int64_t)) ||
        wkspace_alloc_i_checked(&pos_buf, marker_ct * sizeof(int32_t)) ||
        wkspace_alloc_c_checked(&marker_ids, marker_ct * max_marker_id_len) ||
	wkspace_alloc_d_checked(&gd_vals, marker_ct * sizeof(double)) ||
        wkspace_alloc_c_checked(&marker_alleles, marker_ct * 2)) {
      goto transposed_to_bed_ret_NOMEM;
    }

    for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
      pos_buf[marker_idx] = (int)mapvals[marker_idx];
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
    while (fgets(tbuf, MAXLINELEN, infile)) {
      // .tmp file, guaranteed to be no spaces in front
      cptr = next_item(tbuf);
      cptr2 = item_endl(cptr);
      if (!cptr2) {
	goto transposed_to_bed_ret_INVALID_FORMAT;
      }
      cptr3 = skip_initial_spaces(cptr2);
      cptr4 = next_item_mult(cptr3, 2);
      uii = cptr2 - cptr;
      memcpy(&(marker_ids[marker_idx * max_marker_id_len]), cptr, uii);
      marker_ids[marker_idx * max_marker_id_len + uii] = '\0';
      if (sscanf(cptr3, "%lg", &(gd_vals[marker_idx])) != 1) {
	goto transposed_to_bed_ret_INVALID_FORMAT;
      }
      marker_alleles[2 * marker_idx] = *cptr4;
      marker_alleles[2 * marker_idx + 1] = *(next_item(cptr4));
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
      for (; marker_idx < ujj; marker_idx++) {
	marker_uidx = (uint32_t)ll_buf[marker_idx];
	if (fprintf(outfile, "%u\t%s\t%g\t%u\t%c\t%c\n", cur_chrom, &(marker_ids[marker_uidx * max_marker_id_len]), gd_vals[marker_uidx], (uint32_t)(ll_buf[marker_idx] >> 32), marker_alleles[2 * marker_uidx], marker_alleles[2 * marker_uidx + 1]) < 0) {
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
    uii = 4294967294U; // last marker uidx
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
  transposed_to_bed_ret_INVALID_FORMAT_3:
    cptr[strlen_se(cptr)] = '\0';
    sprintf(logbuf, "Error: Missing entries at marker %s in .tped file.\n", cptr);
  transposed_to_bed_ret_INVALID_FORMAT_5:
    putchar('\r');
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  transposed_to_bed_ret_INVALID_FORMAT_4:
    cptr[strlen_se(cptr)] = '\0';
    putchar('\r');
    sprintf(logbuf, "Error: half-missing call at marker %s, indiv %" PRIuPTR " in .tped file.\n", cptr, indiv_idx);
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
  if (infile) {
    fclose(infile);
  }
  if (bimfile) {
    fclose(bimfile);
  }
  if (outfile) {
    fclose(outfile);
  }
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
    memcpy(alleles, "121", 3);
  } else {
    memcpy(alleles, "ABA", 3);
  }
  if (wkspace_alloc_uc_checked(&writebuf, indiv_ct4)) {
    goto generate_dummy_ret_NOMEM;
  }
  memcpy(outname_end, ".bim", 5);
  if (fopen_checked(&outfile, outname, "w")) {
    goto generate_dummy_ret_OPEN_FAIL;
  }
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
      if (fprintf(outfile, "1\tsnp%u\t0\t%u\t%c\t%c\n", uii, uii, alleles[ujj], alleles[ujj + 1]) < 0) {
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
      if (fprintf(outfile, "1\tsnp%u\t0\t%u\t%c\t%c\n", uii, uii, alleles[ujj], alleles[ujj + 1]) < 0) {
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
      if (fprintf(outfile, "per%u per%u 0 0 2 %g\n", uii, uii, dxx) < 0) {
	goto generate_dummy_ret_OPEN_FAIL;
      }
    }
  } else {
    for (uii = 0; uii < indiv_ct; uii++) {
      if (!(uii % 32)) {
	urand = sfmt_genrand_uint32(&sfmt);
      }
      if (pheno_m_check && (sfmt_genrand_uint32(&sfmt) <= pheno_m32)) {
	if (fprintf(outfile, "per%u per%u 0 0 2 -9\n", uii, uii) < 0) {
	  goto generate_dummy_ret_OPEN_FAIL;
	}
      } else {
	if (fprintf(outfile, "per%u per%u 0 0 2 %c\n", uii, uii, (char)((urand & 1) + '1')) < 0) {
	  goto generate_dummy_ret_OPEN_FAIL;
	}
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
    loop_end = ((uint64_t)(pct * marker_ct)) / 100LLU;
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
  sprintf(logbuf, "Dummy data generated (%" PRIuPTR " %s, %" PRIuPTR " markers).\n", indiv_ct, species_str(indiv_ct), marker_ct);
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

int32_t recode_allele_load(char* recode_allele_name, char*** allele_missing_ptr, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, char* marker_alleles, uintptr_t max_marker_allele_len, uintptr_t* recode_allele_reverse, char* recode_allele_extra) {
  FILE* rafile = NULL;
  uint32_t missing_allele = 0;
  int32_t duplicate_fail = 1;
  uintptr_t rae_size = 0;
  char* sorted_ids;
  int32_t* id_map;
  char* bufptr;
  char* bufptr2;
  int32_t retval;
  uint32_t slen;
  uint32_t alen;
  int32_t ii;
  uintptr_t marker_uidx;
  char cc;
  if (fopen_checked(&rafile, recode_allele_name, "r")) {
    goto recode_allele_load_ret_OPEN_FAIL;
  }
  retval = sort_item_ids(&sorted_ids, &id_map, unfiltered_marker_ct, marker_exclude, unfiltered_marker_ct - marker_ct, marker_ids, max_marker_id_len, strcmp_deref, &duplicate_fail);
  if (retval) {
    goto recode_allele_load_ret_1;
  }
  tbuf[MAXLINELEN - 1] = ' ';
  while (fgets(tbuf, MAXLINELEN, rafile)) {
    if (!tbuf[MAXLINELEN - 1]) {
      logprint("Error: Pathologically long line in --recode-allele file.\n");
      goto recode_allele_load_ret_INVALID_FORMAT;
    }
    bufptr = skip_initial_spaces(tbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    slen = strlen_se(bufptr);
    bufptr2 = skip_initial_spaces(&(bufptr[slen]));
    if (is_eoln_kns(*bufptr2)) {
      logprint("Error: --recode-allele line has only one entry.\n");
      goto recode_allele_load_ret_INVALID_FORMAT;
    }
    alen = strlen_se(bufptr2);
    bufptr[slen] = '\0';
    ii = bsearch_str(bufptr, sorted_ids, max_marker_id_len, 0, marker_ct - 1);
    if (ii != -1) {
      marker_uidx = id_map[ii];
      if (max_marker_allele_len == 1) {
	cc = *bufptr2;
	if ((cc == marker_alleles[2 * marker_uidx]) && (alen == 1)) {
	  clear_bit_noct(recode_allele_reverse, marker_uidx);
	} else if ((cc == marker_alleles[2 * marker_uidx + 1]) && (alen == 1)) {
	  set_bit_noct(recode_allele_reverse, marker_uidx);
	} else {
	  if (rae_size + alen + 1 > wkspace_left) {
	    goto recode_allele_load_ret_NOMEM;
	  }
	  missing_allele = 1;
	  (*allele_missing_ptr)[marker_uidx] = &(recode_allele_extra[rae_size]);
	  memcpy(&(recode_allele_extra[rae_size]), bufptr2, alen);
	  rae_size += alen + 1;
	  recode_allele_extra[rae_size - 1] = '\0';
	}
      } else {
	bufptr2[alen++] = '\0';
	if ((alen <= max_marker_allele_len) && (!memcmp(bufptr2, &(marker_alleles[2 * marker_uidx * max_marker_allele_len]), alen))) {
	  clear_bit_noct(recode_allele_reverse, marker_uidx);
	} else if ((alen <= max_marker_allele_len) && (!memcmp(bufptr2, &(marker_alleles[(2 * marker_uidx + 1) * max_marker_allele_len]), alen))) {
	  set_bit_noct(recode_allele_reverse, marker_uidx);
	} else {
	  missing_allele = 1;
	  (*allele_missing_ptr)[marker_uidx] = &(recode_allele_extra[rae_size]);
	  memcpy(&(recode_allele_extra[rae_size]), bufptr2, alen);
	  rae_size += alen;
	}
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

int32_t recode_load_to(unsigned char* loadbuf, FILE* bedfile, int32_t bed_offset, uintptr_t unfiltered_marker_ct, uintptr_t marker_idx, uintptr_t marker_idx_end, uintptr_t* marker_exclude, uintptr_t* marker_uidx_ptr, uintptr_t unfiltered_indiv_ct4) {
  uintptr_t marker_uidx = *marker_uidx_ptr;
  uintptr_t ulii;
  while (marker_idx < marker_idx_end) {
    if (is_set(marker_exclude, marker_uidx)) {
      marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx + 1);
      if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	return RET_READ_FAIL;
      }
    }
    if (unfiltered_marker_ct - marker_uidx > marker_idx_end - marker_idx) {
      ulii = next_set_unsafe(marker_exclude, marker_uidx) - marker_uidx;
    } else {
      ulii = unfiltered_marker_ct - marker_uidx;
    }
    marker_uidx += ulii;
    marker_idx += ulii;
    ulii *= unfiltered_indiv_ct4;
    if (fread(loadbuf, 1, ulii, bedfile) < ulii) {
      return RET_READ_FAIL;
    }
    loadbuf = &(loadbuf[ulii]);
  }
  *marker_uidx_ptr = marker_uidx;
  return 0;
}

static inline int32_t recode_write_first_cols(FILE* outfile, uintptr_t indiv_uidx, char delimiter, char* person_ids, uintptr_t max_person_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* sex_nm, uintptr_t* sex_male, uintptr_t* pheno_nm, uintptr_t* pheno_c, char* output_missing_pheno, int32_t phenos_present, double* pheno_d, double missing_phenod) {
  char* cptr = &(person_ids[indiv_uidx * max_person_id_len]);
  uintptr_t ulii = strlen_se(cptr);
  if (fwrite_checked(cptr, ulii, outfile)) {
    return -1;
  }
  if (fprintf(outfile, "%c%s%c%s%c%s%c%c%c", delimiter, &(cptr[ulii + 1]), delimiter, paternal_ids? (&(paternal_ids[indiv_uidx * max_paternal_id_len])) : "0", delimiter, maternal_ids? (&(maternal_ids[indiv_uidx * max_maternal_id_len])) : "0", delimiter, sexchar(sex_nm, sex_male, indiv_uidx), delimiter) < 0) {
    return -1;
  }
  if (!is_set(pheno_nm, indiv_uidx)) {
    if (fprintf(outfile, "%s%c", output_missing_pheno, delimiter) < 0) {
      return -1;
    }
  } else if (pheno_c) {
    if (fprintf(outfile, "%c%c", is_set(pheno_c, indiv_uidx)? '2' : '1', delimiter) < 0) {
      return -1;
    }
  } else {
    if (fprintf(outfile, "%g%c", pheno_d[indiv_uidx], delimiter) < 0) {
      return -1;
    }
  }
  return 0;
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

int32_t recode(uint32_t recode_modifier, FILE* bedfile, int32_t bed_offset, FILE* famfile, FILE* bimfile, FILE** outfile_ptr, char* outname, char* outname_end, char* recode_allele_name, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_ct, char* marker_ids, uintptr_t max_marker_id_len, char* marker_alleles, uintptr_t max_marker_allele_len, uint32_t* marker_pos, uintptr_t* marker_reverse, char* person_ids, uintptr_t max_person_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* sex_nm, uintptr_t* sex_male, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, double missing_phenod, char output_missing_geno, char* output_missing_pheno, uint32_t set_hh_missing, uint32_t xmhh_exists, uint32_t nxmhh_exists, Chrom_info* chrom_info_ptr) {
  FILE* ref_file = NULL;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  unsigned char* wkspace_mark = wkspace_base;
  int32_t affection = (pheno_c != NULL);
  int32_t phenos_present = (affection || (pheno_d != NULL));
  char delimiter = (recode_modifier & RECODE_TAB)? '\t' : ' ';
  uintptr_t* recode_allele_reverse = NULL;
  char** allele_missing = NULL;
  char* recode_allele_extra = NULL;
  char delim2 = delimiter;
  uintptr_t* indiv_include2 = NULL;
  uintptr_t* indiv_male_include2 = NULL;
  uint32_t lgen_ref = (recode_modifier & RECODE_LGEN_REF);
  uint32_t rlist = (recode_modifier & RECODE_RLIST);
  int32_t retval = 0;
  time_t rawtime;
  FILE* outfile;
  uint32_t is_x;
  uint32_t is_haploid;
  uint32_t chrom_fo_idx;
  uint32_t chrom_idx;
  uint32_t chrom_end;
  uintptr_t marker_uidx;
  uintptr_t marker_idx;
  uintptr_t indiv_uidx;
  uintptr_t indiv_idx;
  unsigned char* loadbuf;
  char* writebuf;
  char* writebufl[4];
  char* writebuflp[4];
  char* writebuflps[4];
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
  char cur_mk_alleles[4];
  char* cur_mk_allelesx_buf;
  char* cur_mk_allelesx[6];
  uint32_t cmalen[4];
  uint32_t pct;
  uint32_t loop_end;
  uintptr_t ulii;
  uint32_t shiftval;
  char* mk_alleles;
  cur_mk_alleles[0] = '\0';
  cur_mk_alleles[1] = '\0';
  cur_mk_alleles[2] = '\0';
  cur_mk_alleles[3] = '\0';
  if ((!xmhh_exists) && (!nxmhh_exists)) {
    set_hh_missing = 0;
  } else if (!set_hh_missing) {
    xmhh_exists = 0;
    nxmhh_exists = 0;
  }
  if (set_hh_missing) {
    if (nxmhh_exists) {
      if (wkspace_alloc_ul_checked(&indiv_include2, unfiltered_indiv_ctl * 2 * sizeof(intptr_t))) {
	goto recode_ret_NOMEM;
      }
      exclude_to_vec_include(unfiltered_indiv_ct, indiv_include2, indiv_exclude);
    }
    if (xmhh_exists) {
      if (wkspace_alloc_ul_checked(&indiv_male_include2, unfiltered_indiv_ctl * 2 * sizeof(intptr_t))) {
        goto recode_ret_NOMEM;
      }
      if (nxmhh_exists) {
	memcpy(indiv_male_include2, indiv_include2, unfiltered_indiv_ctl * 2 * sizeof(intptr_t));
      } else {
        exclude_to_vec_include(unfiltered_indiv_ct, indiv_male_include2, indiv_exclude);
      }
      vec_include_mask_in(unfiltered_indiv_ct, indiv_male_include2, sex_nm);
      vec_include_mask_in(unfiltered_indiv_ct, indiv_male_include2, sex_male);
    }
  }
  if (recode_modifier & RECODE_TRANSPOSE) {
    if (max_marker_allele_len == 1) {
      ulii = 4;
    } else {
      ulii = 2 * max_marker_allele_len;
    }
    if (wkspace_alloc_c_checked(&writebuf, indiv_ct * ulii)) {
      goto recode_ret_NOMEM;
    }
  } else if (recode_modifier & (RECODE_LGEN | RECODE_LGEN_REF)) {
    if (wkspace_alloc_c_checked(&writebuf, 3)) {
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
    if (max_marker_allele_len == 1) {
      ulii += 2;
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
    } else if (max_marker_allele_len != 1) {
      logprint("Error: --recode 23 cannot be used with multi-character allele names.\n");
      goto recode_ret_INVALID_FORMAT;
    }
    // just the chromosome code, and space for single-char alleles
    if (wkspace_alloc_c_checked(&writebuf, 5)) {
      goto recode_ret_NOMEM;
    }
  } else {
    if (recode_modifier & RECODE_AD) {
      ulii = 6;
    } else if (recode_modifier & (RECODE_A | RECODE_COMPOUND)) {
      if ((max_marker_allele_len != 1) && (recode_modifier & RECODE_COMPOUND)) {
	logprint("Error: --recode compound-genotypes cannot be used with multi-character allele\nnames.\n");
	goto recode_ret_INVALID_FORMAT;
      }
      ulii = 3;
    } else {
      if (max_marker_allele_len == 1) {
        ulii = 4;
      } else {
	ulii = 2 * max_marker_allele_len;
      }
    }
    if (wkspace_alloc_c_checked(&writebuf, marker_ct * ulii)) {
      goto recode_ret_NOMEM;
    }
    if (recode_allele_name) {
      ulii = ((unfiltered_marker_ct + (BITCT - 1)) / BITCT) * sizeof(intptr_t);
      if (wkspace_alloc_ul_checked(&recode_allele_reverse, ulii * sizeof(intptr_t))) {
	goto recode_ret_NOMEM;
      }
      memcpy(recode_allele_reverse, marker_reverse, ulii);
      marker_reverse = recode_allele_reverse;
      allele_missing = (char**)wkspace_alloc(unfiltered_marker_ct * sizeof(char**));
      if (!allele_missing) {
	goto recode_ret_NOMEM;
      }
      recode_allele_extra = (char*)wkspace_base;
      fill_ulong_zero((uintptr_t*)allele_missing, unfiltered_marker_ct);
      retval = recode_allele_load(recode_allele_name, &allele_missing, unfiltered_marker_ct, marker_exclude, marker_ct, marker_ids, max_marker_id_len, marker_alleles, max_marker_allele_len, recode_allele_reverse, recode_allele_extra);
      if (retval) {
	goto recode_ret_1;
      }
    }
  }
  if (recode_modifier & RECODE_12) {
    if (wkspace_alloc_c_checked(&mk_alleles, unfiltered_marker_ct * max_marker_allele_len * 2)) {
      goto recode_ret_NOMEM;
    }
    marker_uidx = 0;
    for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
      marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx);
      if (max_marker_allele_len == 1) {
	mk_alleles[2 * marker_uidx] = '1';
	mk_alleles[2 * marker_uidx + 1] = '2';
      } else {
	mk_alleles[2 * marker_uidx * max_marker_allele_len] = '1';
	mk_alleles[2 * marker_uidx * max_marker_allele_len + 1] = '\0';
	mk_alleles[(2 * marker_uidx + 1) * max_marker_allele_len] = '2';
	mk_alleles[(2 * marker_uidx + 1) * max_marker_allele_len + 1] = '\0';
      }
      marker_uidx++;
    }
  } else {
    mk_alleles = marker_alleles;
  }
  if (max_marker_allele_len > 1) {
    if (wkspace_alloc_c_checked(&cur_mk_allelesx_buf, 4 * max_marker_allele_len)) {
      goto recode_ret_NOMEM;
    }
    cur_mk_allelesx[0] = cur_mk_allelesx_buf;
    cur_mk_allelesx[1] = &(cur_mk_allelesx_buf[max_marker_allele_len]);
    cur_mk_allelesx[2] = &(cur_mk_allelesx_buf[max_marker_allele_len * 2]);
    cur_mk_allelesx[3] = &(cur_mk_allelesx_buf[max_marker_allele_len * 3]);
  }
  if (fseeko(bedfile, bed_offset, SEEK_SET)) {
    goto recode_ret_READ_FAIL;
  }
  marker_uidx = 0;
  marker_idx = 0;
  loadbuf = wkspace_base;
  chrom_fo_idx = 0;
  refresh_chrom_info(chrom_info_ptr, marker_uidx, set_hh_missing, 0, &chrom_end, &chrom_fo_idx, &is_x, &is_haploid);
  chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
  if (recode_modifier & RECODE_TRANSPOSE) {
    strcpy(outname_end, ".tped");
    if (fopen_checked(outfile_ptr, outname, "w")) {
      goto recode_ret_OPEN_FAIL;
    }
    sprintf(logbuf, "--recode to %s + .tfam... ", outname);
    logprintb();
    fputs("0%", stdout);
    if (max_marker_allele_len == 1) {
      if ((recode_modifier & (RECODE_TAB | RECODE_DELIMX)) == RECODE_TAB) {
	loop_end = indiv_ct * 4;
	for (ulii = 1; ulii < loop_end; ulii += 4) {
	  writebuf[ulii] = ' ';
	  writebuf[ulii + 2] = '\t';
	}
      } else {
	memset(writebuf, delimiter, indiv_ct * 4 - 1);
      }
      writebuf[indiv_ct * 4 - 1] = '\n';
    } else {
      if ((recode_modifier & (RECODE_TAB | RECODE_DELIMX)) == RECODE_TAB) {
	delim2 = ' '; // within genotype
      }
    }
    rewind(bimfile);

    for (pct = 1; pct <= 100; pct++) {
      loop_end = (((uint64_t)pct) * marker_ct) / 100;
      for (; marker_idx < loop_end; marker_idx++) {
	if (is_set(marker_exclude, marker_uidx)) {
	  if (get_next_noncomment_excl(bimfile, &wbufptr, marker_exclude, &marker_uidx)) {
	    goto recode_ret_READ_FAIL;
	  }
	  if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	    goto recode_ret_READ_FAIL;
	  }
	} else {
	  if (get_next_noncomment(bimfile, &wbufptr)) {
	    goto recode_ret_READ_FAIL;
	  }
	}
	if (marker_uidx >= chrom_end) {
	  chrom_fo_idx++;
	  refresh_chrom_info(chrom_info_ptr, marker_uidx, set_hh_missing, 0, &chrom_end, &chrom_fo_idx, &is_x, &is_haploid);
	  chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	}

	for (indiv_idx = 0; indiv_idx < 3; indiv_idx++) {
	  cptr = item_endnn(wbufptr);
	  ulii = 1 + (uintptr_t)(cptr - wbufptr);
	  *cptr = delimiter;
	  if (fwrite_checked(wbufptr, ulii, *outfile_ptr)) {
	    goto recode_ret_WRITE_FAIL;
	  }
	  wbufptr = skip_initial_spaces(cptr);
	}
	ulii = strlen_se(wbufptr);
	wbufptr[ulii] = delimiter;
	if (fwrite_checked(wbufptr, ulii + 1, *outfile_ptr)) {
	  goto recode_ret_WRITE_FAIL;
	}

	if (fread(loadbuf, 1, unfiltered_indiv_ct4, bedfile) < unfiltered_indiv_ct4) {
	  goto recode_ret_READ_FAIL;
	}
	if (is_haploid) {
	  if (is_x) {
	    if (xmhh_exists) {
	      hh_reset(loadbuf, indiv_male_include2, unfiltered_indiv_ct);
	    }
	  } else if (nxmhh_exists) {
	    hh_reset(loadbuf, indiv_include2, unfiltered_indiv_ct);
	  }
	}
	indiv_uidx = 0;
	wbufptr = writebuf;
	if (max_marker_allele_len == 1) {
	  cur_mk_alleles[0] = mk_alleles[2 * marker_uidx];
	  cur_mk_alleles[1] = mk_alleles[2 * marker_uidx + 1];
	  if (is_set(marker_reverse, marker_uidx)) {
	    cur_mk_alleles[2] = cur_mk_alleles[1];
	    cur_mk_alleles[3] = *cur_mk_alleles;
	  } else {
	    cur_mk_alleles[2] = *cur_mk_alleles;
	    cur_mk_alleles[3] = cur_mk_alleles[1];
	  }
	  for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
	    indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx);
	    ucc = (loadbuf[indiv_uidx / 4] >> ((indiv_uidx % 4) * 2)) & 3;
	    if (ucc) {
	      if (ucc == 2) {
		*wbufptr = cur_mk_alleles[2];
		wbufptr[2] = cur_mk_alleles[3];
	      } else if (ucc == 3) {
		*wbufptr = cur_mk_alleles[1];
		wbufptr[2] = cur_mk_alleles[1];
	      } else {
		*wbufptr = output_missing_geno;
		wbufptr[2] = output_missing_geno;
	      }
	    } else {
	      *wbufptr = *cur_mk_alleles;
	      wbufptr[2] = *cur_mk_alleles;
	    }
	    wbufptr = &(wbufptr[4]);
	    indiv_uidx++;
	  }
	  if (fwrite_checked(writebuf, indiv_ct * 4, *outfile_ptr)) {
	    goto recode_ret_WRITE_FAIL;
	  }
	} else {
	  init_cur_mk_allelesx(&(mk_alleles[2 * marker_uidx * max_marker_allele_len]), max_marker_allele_len, is_set(marker_reverse, marker_uidx), cur_mk_allelesx, cmalen, delimiter, delim2);
	  for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
	    indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx);
	    ucc = (loadbuf[indiv_uidx / 4] >> ((indiv_uidx % 4) * 2)) & 3;
	    if (ucc) {
	      if (ucc == 2) {
		memcpy(wbufptr, cur_mk_allelesx[4], cmalen[2]);
		wbufptr = &(wbufptr[cmalen[2]]);
		memcpy(wbufptr, cur_mk_allelesx[5], cmalen[3]);
		wbufptr = &(wbufptr[cmalen[3]]);
	      } else if (ucc == 3) {
		memcpy(wbufptr, cur_mk_allelesx[2], cmalen[1]);
		wbufptr = &(wbufptr[cmalen[1]]);
		memcpy(wbufptr, cur_mk_allelesx[3], cmalen[1]);
		wbufptr = &(wbufptr[cmalen[1]]);
	      } else {
		*wbufptr++ = output_missing_geno;
		*wbufptr++ = delim2;
		*wbufptr++ = output_missing_geno;
		*wbufptr++ = delimiter;
	      }
	    } else {
	      memcpy(wbufptr, cur_mk_allelesx[0], cmalen[0]);
	      wbufptr = &(wbufptr[cmalen[0]]);
	      memcpy(wbufptr, cur_mk_allelesx[1], cmalen[0]);
	      wbufptr = &(wbufptr[cmalen[0]]);
	    }
	    indiv_uidx++;
	  }
	  wbufptr[-1] = '\n';
	  if (fwrite_checked(writebuf, (uintptr_t)(wbufptr - writebuf), *outfile_ptr)) {
	    goto recode_ret_WRITE_FAIL;
	  }
	}
	marker_uidx++;
      }
      if (pct < 100) {
	if (pct > 10) {
	  putchar('\b');
	}
	printf("\b\b%u%%", pct);
	fflush(stdout);
      }
    }
  } else if (recode_modifier & RECODE_23) {
    strcpy(outname_end, ".txt");
    if (fopen_checked(outfile_ptr, outname, "w")) {
      goto recode_ret_OPEN_FAIL;
    }
    sprintf(logbuf, "--recode 23 to %s... ", outname);
    logprintb();
    time(&rawtime);
    if (fprintf(*outfile_ptr, "# This data file generated by WDIST at: %s", ctime(&rawtime)) < 0) {
      goto recode_ret_WRITE_FAIL;
    }
    if (fputs(
"#\n"
"# Below is a text version of your data.  Fields are TAB-separated.\n"
"# Each line corresponds to a single SNP.  For each SNP, we provide its\n"
"# identifier, its location on a reference human genome, and the genotype call.\n"
"# For further information (e.g. which reference build was used), consult the\n"
"# original source of your data.\n"
"#\n"
"# rsid\tchromosome\tposition\tgenotype\n"
, *outfile_ptr) == EOF) {
      goto recode_ret_WRITE_FAIL;
    }
    chrom_print_human_terminate(writebuf, chrom_idx);
    indiv_uidx = next_non_set_unsafe(indiv_exclude, 0);
    // May want to add special handling of X/XY, depending on what our pipeline
    // generates.
    for (; marker_idx < marker_ct; marker_idx++) {
      if (is_set(marker_exclude, marker_uidx)) {
	if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	  goto recode_ret_READ_FAIL;
	}
      }
      if (marker_uidx >= chrom_end) {
	chrom_fo_idx++;
	refresh_chrom_info(chrom_info_ptr, marker_uidx, set_hh_missing, 0, &chrom_end, &chrom_fo_idx, &is_x, &is_haploid);
	chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	chrom_print_human_terminate(writebuf, chrom_idx);
      }

      if (fread(loadbuf, 1, unfiltered_indiv_ct4, bedfile) < unfiltered_indiv_ct4) {
	goto recode_ret_READ_FAIL;
      }
      if (is_haploid) {
	if (is_x) {
	  if (xmhh_exists) {
	    hh_reset(loadbuf, indiv_male_include2, unfiltered_indiv_ct);
	  }
	} else if (nxmhh_exists) {
	  hh_reset(loadbuf, indiv_include2, unfiltered_indiv_ct);
	}
      }

      cur_mk_alleles[0] = mk_alleles[2 * marker_uidx];
      cur_mk_alleles[1] = mk_alleles[2 * marker_uidx + 1];
      if (is_set(marker_reverse, marker_uidx)) {
	cur_mk_alleles[2] = cur_mk_alleles[1];
	cur_mk_alleles[3] = *cur_mk_alleles;
      } else {
	cur_mk_alleles[2] = *cur_mk_alleles;
	cur_mk_alleles[3] = cur_mk_alleles[1];
      }
      ucc = (loadbuf[indiv_uidx / 4] >> ((indiv_uidx % 4) * 2)) & 3;
      if (ucc) {
	if (ucc == 2) {
	  writebuf[3] = cur_mk_alleles[2];
	  writebuf[4] = cur_mk_alleles[3];
	} else if (ucc == 3) {
	  writebuf[3] = cur_mk_alleles[1];
	  writebuf[4] = cur_mk_alleles[1];
	} else {
	  writebuf[3] = output_missing_geno;
	  writebuf[4] = output_missing_geno;
	}
      } else {
	writebuf[3] = *cur_mk_alleles;
	writebuf[4] = *cur_mk_alleles;
      }
      if (fprintf(*outfile_ptr, "%s\t%s\t%u\t%c%c\n", &(marker_ids[marker_uidx * max_marker_id_len]), writebuf, marker_pos[marker_uidx], writebuf[3], writebuf[4]) < 0) {
	goto recode_ret_WRITE_FAIL;
      }
      marker_uidx++;
    }
  } else if (recode_modifier & (RECODE_LGEN | RECODE_LGEN_REF)) {
    strcpy(outname_end, ".ref");
    if (fopen_checked(&ref_file, outname, "w")) {
      goto recode_ret_OPEN_FAIL;
    }
    strcpy(outname_end, ".lgen");
    if (fopen_checked(outfile_ptr, outname, "w")) {
      goto recode_ret_OPEN_FAIL;
    }
    outfile = *outfile_ptr;
    if (delimiter == ' ') {
      indiv_delim_convert(unfiltered_indiv_ct, indiv_exclude, indiv_ct, person_ids, max_person_id_len, '\t', ' ');
    } else {
      if (!(recode_modifier & RECODE_DELIMX)) {
	delim2 = ' ';
      }
    }
    sprintf(logbuf, "--recode to %s + .map + .fam%s... ", outname, lgen_ref? " + .ref" : "");
    logprintb();
    fputs("0%", stdout);
    ucc2 = 1;
    tbuf[0] = delimiter;
    for (pct = 1; pct <= 100; pct++) {
      loop_end = (((uint64_t)pct) * marker_ct) / 100;
      for (; marker_idx < loop_end; marker_idx++) {
	if (is_set(marker_exclude, marker_uidx)) {
	  marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx + 1);
	  if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	    goto recode_ret_READ_FAIL;
	  }
	}
	if (marker_uidx >= chrom_end) {
	  chrom_fo_idx++;
	  refresh_chrom_info(chrom_info_ptr, marker_uidx, set_hh_missing, 0, &chrom_end, &chrom_fo_idx, &is_x, &is_haploid);
	}
	if (fread(loadbuf, 1, unfiltered_indiv_ct4, bedfile) < unfiltered_indiv_ct4) {
	  goto recode_ret_READ_FAIL;
	}
	if (is_haploid) {
	  if (is_x) {
	    if (xmhh_exists) {
	      hh_reset(loadbuf, indiv_male_include2, unfiltered_indiv_ct);
	    }
	  } else if (nxmhh_exists) {
	    hh_reset(loadbuf, indiv_include2, unfiltered_indiv_ct);
	  }
	}
	wbufptr = &(marker_ids[marker_uidx * max_marker_id_len]);
	alen = strlen(wbufptr);
	memcpy(&(tbuf[1]), wbufptr, alen);
	cptr = &(tbuf[alen + 1]);
	if (delimiter == ' ') {
	  *cptr++ = ' ';
	  *cptr++ = ' ';
	} else {
	  *cptr++ = '\t';
	}
	indiv_uidx = 0;
	if (lgen_ref) {
	  fputs(wbufptr, ref_file);
	  putc(delimiter, ref_file);
	}
	if (max_marker_allele_len == 1) {
	  cur_mk_alleles[0] = mk_alleles[2 * marker_uidx];
	  cur_mk_alleles[1] = mk_alleles[2 * marker_uidx + 1];
	  if (is_set(marker_reverse, marker_uidx)) {
	    cur_mk_alleles[2] = cur_mk_alleles[1];
	    cur_mk_alleles[3] = *cur_mk_alleles;
	  } else {
	    cur_mk_alleles[2] = *cur_mk_alleles;
	    cur_mk_alleles[3] = cur_mk_alleles[1];
	  }
	  if (lgen_ref) {
	    putc(cur_mk_alleles[3], ref_file);
	    putc(delimiter, ref_file);
	    putc(cur_mk_alleles[2], ref_file);
	    if (putc('\n', ref_file) == EOF) {
	      goto recode_ret_WRITE_FAIL;
	    }
	    ucc2 = is_set(marker_reverse, marker_uidx)? 0 : 3;
	  }
	  alen = 4 + (uintptr_t)(cptr - tbuf);
	  *cptr++ = *cur_mk_alleles;
	  *cptr++ = delim2;
	  *cptr++ = *cur_mk_alleles;
	  *cptr++ = '\n';
	  memcpy(cptr, tbuf, alen);
	  memcpy(&(cptr[alen]), tbuf, alen);
	  memcpy(&(cptr[alen * 2]), tbuf, alen);
	  cptr[alen - 4] = '0';
	  cptr[alen - 2] = '0';
	  cptr[2 * alen - 4] = cur_mk_alleles[2];
	  cptr[2 * alen - 2] = cur_mk_alleles[3];
	  cptr[3 * alen - 4] = cur_mk_alleles[1];
	  cptr[3 * alen - 2] = cur_mk_alleles[1];
	  for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
	    indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx);
	    ucc = (loadbuf[indiv_uidx / 4] >> ((indiv_uidx % 4) * 2)) & 3;
	    if (ucc != ucc2) {
	      fputs(&(person_ids[indiv_uidx * max_person_id_len]), outfile);
	      if (fwrite_checked(&(tbuf[ucc * alen]), alen, outfile)) {
		goto recode_ret_WRITE_FAIL;
	      }
	    }
	    indiv_uidx++;
	  }
	} else {
	  *cptr = '\0';
	  init_cur_mk_allelesx(&(mk_alleles[2 * marker_uidx * max_marker_allele_len]), max_marker_allele_len, is_set(marker_reverse, marker_uidx), cur_mk_allelesx, cmalen, '\0', '\0');
	  if (lgen_ref) {
	    fputs(cur_mk_allelesx[5], ref_file);
	    putc(delimiter, ref_file);
	    fputs(cur_mk_allelesx[4], ref_file);
	    if (putc('\n', ref_file) == EOF) {
	      goto recode_ret_WRITE_FAIL;
	    }
	    ucc2 = is_set(marker_reverse, marker_uidx)? 0 : 3;
	  }
	  for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
	    indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx);
	    ucc = (loadbuf[indiv_uidx / 4] >> ((indiv_uidx % 4) * 2)) & 3;
	    if (ucc != ucc2) {
	      fputs(&(person_ids[indiv_uidx * max_person_id_len]), outfile);
	      fputs(tbuf, outfile);
	      if (ucc == 3) {
		fputs(cur_mk_allelesx[1], outfile);
		putc(delim2, outfile);
		fputs(cur_mk_allelesx[1], outfile);
	      } else if (ucc == 2) {
		fputs(cur_mk_allelesx[4], outfile);
		putc(delim2, outfile);
		fputs(cur_mk_allelesx[5], outfile);
	      } else if (!ucc) {
		fputs(cur_mk_allelesx[0], outfile);
		putc(delim2, outfile);
		fputs(cur_mk_allelesx[0], outfile);
	      } else {
		putc('0', outfile);
		putc(delim2, outfile);
		putc('0', outfile);
	      }
	      if (putc('\n', outfile) == EOF) {
		goto recode_ret_WRITE_FAIL;
	      }
	    }
	    indiv_uidx++;
	  }
	}
	marker_uidx++;
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
    strcpy(outname_end, ".raw");
    if (wkspace_left >= (uint64_t)unfiltered_indiv_ct4 * marker_ct) {
      if (fopen_checked(outfile_ptr, outname, "w")) {
	goto recode_ret_OPEN_FAIL;
      }
      outfile = *outfile_ptr;
      if (fputs("FID IID PAT MAT SEX PHENOTYPE", outfile) == EOF) {
	goto recode_ret_WRITE_FAIL;
      }
      for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
	marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx);
	cptr = &(marker_ids[marker_uidx * max_marker_id_len]);
	if ((max_marker_allele_len == 1) && ((!allele_missing) || (!allele_missing[marker_uidx]))) {
	  cc = mk_alleles[2 * marker_uidx + is_set(marker_reverse, marker_uidx)];
	  putc(' ', outfile);
	  fputs(cptr, outfile);
	  putc('_', outfile);
	  if (putc(cc, outfile) == EOF) {
	    goto recode_ret_WRITE_FAIL;
	  }
	  if (recode_modifier & RECODE_AD) {
	    putc(' ', outfile);
	    fputs(cptr, outfile);
	    if (fputs("_HET", outfile) == EOF) {
	      goto recode_ret_WRITE_FAIL;
	    }
	  }
	} else {
	  if (allele_missing && (allele_missing[marker_uidx])) {
	    aptr = allele_missing[marker_uidx];
	  } else {
	    aptr = &(mk_alleles[(2 * marker_uidx + is_set(marker_reverse, marker_uidx)) * max_marker_allele_len]);
	  }
	  if (recode_modifier & RECODE_A) {
	    if (fprintf(outfile, " %s_%s", cptr, aptr) < 0) {
	      goto recode_ret_WRITE_FAIL;
	    }
	  } else {
	    if (fprintf(outfile, " %s_%s %s_HET", cptr, aptr, cptr) < 0) {
	      goto recode_ret_WRITE_FAIL;
	    }
	  }
	}
	marker_uidx++;
      }
      if (putc('\n', outfile) == EOF) {
	goto recode_ret_WRITE_FAIL;
      }
      marker_uidx = 0;
      marker_idx = 0;
      sprintf(logbuf, "--recode to %s... ", outname);
      logprintb();
      retval = recode_load_to(loadbuf, bedfile, bed_offset, unfiltered_marker_ct, 0, marker_ct, marker_exclude, &marker_uidx, unfiltered_indiv_ct4);
      if (retval) {
	goto recode_ret_1;
      }
      if (set_hh_missing) {
	marker_uidx = 0;
	for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
	  marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx);
	  if (marker_uidx >= chrom_end) {
	    chrom_fo_idx++;
	    refresh_chrom_info(chrom_info_ptr, marker_uidx, set_hh_missing, 0, &chrom_end, &chrom_fo_idx, &is_x, &is_haploid);
	  }
	  if (is_haploid) {
	    if (is_x) {
	      if (xmhh_exists) {
		hh_reset(&(loadbuf[marker_uidx * unfiltered_indiv_ct4]), indiv_male_include2, unfiltered_indiv_ct);
	      }
	    } else if (nxmhh_exists) {
	      hh_reset(&(loadbuf[marker_uidx * unfiltered_indiv_ct4]), indiv_include2, unfiltered_indiv_ct);
	    }
	  }
	  marker_uidx++;
	}
      }
      fputs("0%", stdout);
      indiv_uidx = 0;
      indiv_idx = 0;
      for (pct = 1; pct <= 100; pct++) {
	loop_end = ((uint64_t)pct * indiv_ct) / 100;
	for (; indiv_idx < loop_end; indiv_idx++) {
	  indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx);
	  if (recode_write_first_cols(outfile, indiv_uidx, delimiter, person_ids, max_person_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, sex_nm, sex_male, pheno_nm, pheno_c, output_missing_pheno, phenos_present, pheno_d, missing_phenod)) {
	    goto recode_ret_WRITE_FAIL;
	  }
	  bufptr = &(loadbuf[indiv_idx / 4]);
	  wbufptr = writebuf;
	  shiftval = (indiv_idx % 4) * 2;
	  marker_uidx = 0;
	  if (recode_modifier & RECODE_A) {
	    for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
	      marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx);
	      ucc = ((*bufptr) >> shiftval) & 3;
	      if (allele_missing && allele_missing[marker_uidx]) {
		if (ucc != 1) {
		  *wbufptr++ = '0';
		} else {
		  *wbufptr++ = 'N';
		  *wbufptr++ = 'A';
		}
	      } else {
		if (ucc) {
		  if (ucc == 2) {
		    *wbufptr++ = '1';
		  } else if (ucc == 3) {
		    *wbufptr++ = is_set(marker_reverse, marker_uidx)? '2' : '0';
		  } else {
		    *wbufptr++ = 'N';
		    *wbufptr++ = 'A';
		  }
		} else {
		  *wbufptr++ = is_set(marker_reverse, marker_uidx)? '0' : '2';
		}
	      }
	      *wbufptr++ = delimiter;
	      bufptr = &(bufptr[unfiltered_indiv_ct4]);
	      marker_uidx++;
	    }
	  } else {
	    for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
	      marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx);
	      ucc = ((*bufptr) >> shiftval) & 3;
	      if (allele_missing && allele_missing[marker_idx]) {
		if (ucc != 1) {
		  *wbufptr++ = '0';
		  *wbufptr++ = delimiter;
		  *wbufptr++ = '0';
		} else {
		  *wbufptr++ = 'N';
		  *wbufptr++ = 'A';
		  *wbufptr++ = delimiter;
		  *wbufptr++ = 'N';
		  *wbufptr++ = 'A';
		}
	      } else {
		if (ucc) {
		  if (ucc == 2) {
		    *wbufptr++ = '1';
		    *wbufptr++ = delimiter;
		    *wbufptr++ = '1';
		  } else if (ucc == 3) {
		    *wbufptr++ = is_set(marker_reverse, marker_uidx)? '2' : '0';
		    *wbufptr++ = delimiter;
		    *wbufptr++ = '0';
		  } else {
		    *wbufptr++ = 'N';
		    *wbufptr++ = 'A';
		    *wbufptr++ = delimiter;
		    *wbufptr++ = 'N';
		    *wbufptr++ = 'A';
		  }
		} else {
		  *wbufptr++ = is_set(marker_reverse, marker_uidx)? '0' : '2';
		  *wbufptr++ = delimiter;
		  *wbufptr++ = '0';
		}
	      }
	      *wbufptr++ = delimiter;
	      bufptr = &(bufptr[unfiltered_indiv_ct4]);
	      marker_uidx++;
	    }
	  }
	  wbufptr[-1] = '\n';
	  ulii = (uintptr_t)(wbufptr - writebuf);
	  if (fwrite_checked(writebuf, ulii, outfile)) {
	    goto recode_ret_WRITE_FAIL;
	  }
	  indiv_uidx++;
	}
	if (pct < 100) {
	  if (pct > 10) {
	    putchar('\b');
	  }
	  printf("\b\b%u%%", pct);
	  fflush(stdout);
	}
      }
      *outfile_ptr = outfile;
    } else {
      sprintf(logbuf, "Error: --recode A%s does not yet support multipass recoding of very large .bed\nfiles; contact the WDIST developers if you need this.\n", (recode_modifier & RECODE_AD)? "D" : "");
      logprintb();
      retval = RET_CALC_NOT_YET_SUPPORTED;
      goto recode_ret_1;
    }
  } else if (recode_modifier & (RECODE_LIST | RECODE_RLIST)) {
    strcpy(outname_end, rlist? ".rlist" : ".list");
    if (fopen_checked(outfile_ptr, outname, "w")) {
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
    for (pct = 1; pct <= 100; pct++) {
      loop_end = (((uint64_t)pct) * marker_ct) / 100;
      for (; marker_idx < loop_end; marker_idx++) {
	if (is_set(marker_exclude, marker_uidx)) {
	  marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx + 1);
	  if (fseeko(bedfile, bed_offset + ((uint64_t)marker_uidx) * unfiltered_indiv_ct4, SEEK_SET)) {
	    goto recode_ret_READ_FAIL;
	  }
	}
	if (marker_uidx >= chrom_end) {
	  chrom_fo_idx++;
	  refresh_chrom_info(chrom_info_ptr, marker_uidx, set_hh_missing, 0, &chrom_end, &chrom_fo_idx, &is_x, &is_haploid);
	  chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	}
	if (fread(loadbuf, 1, unfiltered_indiv_ct4, bedfile) < unfiltered_indiv_ct4) {
	  goto recode_ret_READ_FAIL;
	}
	if (is_haploid) {
	  if (is_x) {
	    if (xmhh_exists) {
	      hh_reset(loadbuf, indiv_male_include2, unfiltered_indiv_ct);
	    }
	  } else if (nxmhh_exists) {
	    hh_reset(loadbuf, indiv_include2, unfiltered_indiv_ct);
	  }
	}
	if (max_marker_allele_len == 1) {
	  cur_mk_alleles[0] = mk_alleles[2 * marker_uidx];
	  cur_mk_alleles[1] = mk_alleles[2 * marker_uidx + 1];
	  if (is_set(marker_reverse, marker_uidx)) {
	    cur_mk_alleles[2] = cur_mk_alleles[1];
	    cur_mk_alleles[3] = *cur_mk_alleles;
	  } else {
	    cur_mk_alleles[2] = *cur_mk_alleles;
	    cur_mk_alleles[3] = cur_mk_alleles[1];
	  }
	} else {
	  init_cur_mk_allelesx(&(mk_alleles[2 * marker_uidx * max_marker_allele_len]), max_marker_allele_len, is_set(marker_reverse, marker_uidx), cur_mk_allelesx, cmalen, '\0', '\0');
	}
	aptr = &(marker_ids[marker_uidx * max_marker_id_len]);
	alen = strlen(aptr);
	for (ulii = 0; ulii < 4; ulii++) {
	  wbufptr = writebufl[ulii];
	  if (!rlist) {
	    wbufptr = uint32_write(chrom_idx, wbufptr);
	    *wbufptr++ = delimiter;
	  }
	  memcpy(wbufptr, aptr, alen);
	  wbufptr += alen;
	  *wbufptr++ = delimiter;
	  if (rlist) {
	    switch (ulii) {
	    case 0:
	    case 3:
	      *wbufptr++ = 'H';
	      *wbufptr++ = 'O';
	      *wbufptr++ = 'M';
	      break;
	    case 1:
	      *wbufptr++ = 'N';
	      *wbufptr++ = 'I';
	      *wbufptr++ = 'L';
	      break;
	    case 2:
	      *wbufptr++ = 'H';
	      *wbufptr++ = 'E';
	      *wbufptr++ = 'T';
	      break;
	    }
	    *wbufptr++ = delimiter;
	  }
	  if (max_marker_allele_len == 1) {
	    switch (ulii) {
	    case 0:
	      *wbufptr++ = *cur_mk_alleles;
	      if (rlist) {
		*wbufptr++ = delimiter;
	      }
	      *wbufptr++ = *cur_mk_alleles;
	      break;
	    case 1:
	      *wbufptr++ = '0';
	      if (rlist) {
		*wbufptr++ = delimiter;
	      }
	      *wbufptr++ = '0';
	      break;
	    case 2:
	      *wbufptr++ = cur_mk_alleles[2];
	      if (rlist) {
		*wbufptr++ = delimiter;
	      }
	      *wbufptr++ = cur_mk_alleles[3];
	      break;
	    case 3:
	      *wbufptr++ = cur_mk_alleles[1];
	      if (rlist) {
		*wbufptr++ = delimiter;
	      }
	      *wbufptr++ = cur_mk_alleles[1];
	      break;
	    }
	  } else {
	    switch (ulii) {
	    case 0:
	      memcpy(wbufptr, cur_mk_allelesx[0], cmalen[0]);
	      wbufptr += cmalen[0];
	      if (rlist) {
		*wbufptr++ = delimiter;
	      }
	      memcpy(wbufptr, cur_mk_allelesx[0], cmalen[0]);
	      wbufptr += cmalen[0];
	      break;
	    case 1:
	      *wbufptr++ = '0';
	      if (rlist) {
		*wbufptr++ = delimiter;
	      }
	      *wbufptr++ = '0';
	      break;
	    case 2:
	      memcpy(wbufptr, cur_mk_allelesx[4], cmalen[2]);
	      wbufptr += cmalen[2];
	      if (rlist) {
	        *wbufptr++ = delimiter;
	      }
	      memcpy(wbufptr, cur_mk_allelesx[5], cmalen[3]);
	      wbufptr += cmalen[3];
	      break;
	    case 3:
	      memcpy(wbufptr, cur_mk_allelesx[1], cmalen[1]);
	      wbufptr += cmalen[1];
	      if (rlist) {
		*wbufptr++ = delimiter;
	      }
	      memcpy(wbufptr, cur_mk_allelesx[1], cmalen[1]);
	      wbufptr += cmalen[1];
	      break;
	    }
	  }
	  writebuflp[ulii] = wbufptr;
	}
	indiv_uidx = 0;
	if (rlist) {
	  for (ulii = 0; ulii < 4; ulii++) {
	    writebuflps[ulii] = writebuflp[ulii];
	  }
	  ucc2 = is_set(marker_reverse, marker_uidx)? 0 : 3;
	  for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
	    indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx);
	    ucc = (loadbuf[indiv_uidx / 4] >> ((indiv_uidx % 4) * 2)) & 3;
	    if (ucc != ucc2) {
	      *(writebuflp[ucc]++) = delimiter;
	      aptr2 = &(person_ids[indiv_uidx * max_person_id_len]);
	      alen2 = strlen(aptr2);
	      memcpy(writebuflp[ucc], &(person_ids[indiv_uidx * max_person_id_len]), alen2);
	      writebuflp[ucc] += alen2;
	    }
	    indiv_uidx++;
	  }
	  if (writebuflp[2] != writebuflps[2]) {
	    *(writebuflp[2]++) = '\n';
	    if (fwrite_checked(writebufl[2], (uintptr_t)(writebuflp[2] - writebufl[2]), *outfile_ptr)) {
	      goto recode_ret_WRITE_FAIL;
	    }
	  }
	  ucc2 = 3 - ucc2;
	  if (writebuflp[ucc2] != writebuflps[ucc2]) {
	    *(writebuflp[ucc2]++) = '\n';
	    if (fwrite_checked(writebufl[ucc2], (uintptr_t)(writebuflp[ucc2] - writebufl[ucc2]), *outfile_ptr)) {
	      goto recode_ret_WRITE_FAIL;
	    }
	  }
	  if (writebuflp[1] != writebuflps[1]) {
	    *(writebuflp[1]++) = '\n';
	    if (fwrite_checked(writebufl[1], (uintptr_t)(writebuflp[1] - writebufl[1]), *outfile_ptr)) {
	      goto recode_ret_WRITE_FAIL;
	    }
	  }
	} else {
	  for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
	    indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx);
	    ucc = (loadbuf[indiv_uidx / 4] >> ((indiv_uidx % 4) * 2)) & 3;
	    *(writebuflp[ucc]++) = delimiter;
	    aptr2 = &(person_ids[indiv_uidx * max_person_id_len]);
	    alen2 = strlen(aptr2);
	    memcpy(writebuflp[ucc], &(person_ids[indiv_uidx * max_person_id_len]), alen2);
	    writebuflp[ucc] += alen2;
	    indiv_uidx++;
	  }
	  *(writebuflp[0]++) = '\n';
	  *(writebuflp[1]++) = '\n';
	  *(writebuflp[2]++) = '\n';
	  *(writebuflp[3]++) = '\n';
	  ulii = is_set(marker_reverse, marker_uidx)? 3 : 0;
	  if (fwrite_checked(writebufl[ulii], (uintptr_t)(writebuflp[ulii] - writebufl[ulii]), *outfile_ptr)) {
	    goto recode_ret_WRITE_FAIL;
	  }
	  if (fwrite_checked(writebufl[2], (uintptr_t)(writebuflp[2] - writebufl[2]), *outfile_ptr)) {
	    goto recode_ret_WRITE_FAIL;
	  }
	  ulii = 3 - ulii;
	  if (fwrite_checked(writebufl[ulii], (uintptr_t)(writebuflp[ulii] - writebufl[ulii]), *outfile_ptr)) {
	    goto recode_ret_WRITE_FAIL;
	  }
	  if (fwrite_checked(writebufl[1], (uintptr_t)(writebuflp[1] - writebufl[1]), *outfile_ptr)) {
	    goto recode_ret_WRITE_FAIL;
	  }
	}
	marker_uidx++;
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
  } else {
    strcpy(outname_end, ".ped");
    if (wkspace_left >= (uint64_t)unfiltered_indiv_ct4 * marker_ct) {
      if (fopen_checked(outfile_ptr, outname, "w")) {
	goto recode_ret_OPEN_FAIL;
      }
      sprintf(logbuf, "--recode to %s + .map... ", outname);
      logprintb();
      retval = recode_load_to(loadbuf, bedfile, bed_offset, unfiltered_marker_ct, 0, marker_ct, marker_exclude, &marker_uidx, unfiltered_indiv_ct4);
      if (retval) {
	goto recode_ret_1;
      }
      if (set_hh_missing) {
	marker_uidx = 0;
	for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
	  marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx);
	  if (marker_uidx >= chrom_end) {
	    chrom_fo_idx++;
	    refresh_chrom_info(chrom_info_ptr, marker_uidx, set_hh_missing, 0, &chrom_end, &chrom_fo_idx, &is_x, &is_haploid);
	  }
	  if (is_haploid) {
	    if (is_x) {
	      if (xmhh_exists) {
		hh_reset(&(loadbuf[marker_uidx * unfiltered_indiv_ct4]), indiv_male_include2, unfiltered_indiv_ct);
	      }
	    } else if (nxmhh_exists) {
	      hh_reset(&(loadbuf[marker_uidx * unfiltered_indiv_ct4]), indiv_include2, unfiltered_indiv_ct);
	    }
	  }
	  marker_uidx++;
	}
      }
      indiv_uidx = 0;
      indiv_idx = 0;
      if (max_marker_allele_len == 1) {
	if (recode_modifier & RECODE_COMPOUND) {
	  memset(writebuf, delimiter, marker_ct * 3 - 1);
	  writebuf[marker_ct * 3 - 1] = '\n';
	} else {
	  if ((recode_modifier & (RECODE_TAB | RECODE_DELIMX)) == RECODE_TAB) {
	    // PLINK actually writes spaces between same-locus alleles even
            // under --tab.  We replicate this for compatibility.
	    loop_end = marker_ct * 4;
	    for (ulii = 1; ulii < loop_end; ulii += 4) {
	      writebuf[ulii] = ' ';
	      writebuf[ulii + 2] = '\t';
	    }
	  } else {
	    memset(writebuf, delimiter, marker_ct * 4 - 1);
	  }
	}
	writebuf[marker_ct * 4 - 1] = '\n';
      } else {
	if ((recode_modifier & (RECODE_TAB | RECODE_DELIMX)) == RECODE_TAB) {
	  delim2 = ' ';
	}
      }
      fputs("0%", stdout);
      for (pct = 1; pct <= 100; pct++) {
	loop_end = ((uint64_t)pct * indiv_ct) / 100;
	for (; indiv_idx < loop_end; indiv_idx++) {
	  indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx);
	  if (recode_write_first_cols(*outfile_ptr, indiv_uidx, delimiter, person_ids, max_person_id_len, paternal_ids, max_paternal_id_len, maternal_ids, max_maternal_id_len, sex_nm, sex_male, pheno_nm, pheno_c, output_missing_pheno, phenos_present, pheno_d, missing_phenod)) {
	    goto recode_ret_WRITE_FAIL;
	  }
	  bufptr = &(loadbuf[indiv_idx / 4]);
	  wbufptr = writebuf;
	  shiftval = (indiv_idx % 4) * 2;
	  marker_uidx = 0;
	  if (recode_modifier & RECODE_COMPOUND) {
	    for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
	      marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx);
	      ucc = ((*bufptr) >> shiftval) & 3;
	      if (ucc) {
		if (ucc == 2) {
		  ucc2 = is_set(marker_exclude, marker_uidx);
		  *wbufptr++ = mk_alleles[2 * marker_uidx + ucc2];
		  *wbufptr = mk_alleles[2 * marker_uidx + (1 ^ ucc2)];
		} else if (ucc == 3) {
		  *wbufptr++ = mk_alleles[2 * marker_uidx + 1];
		  *wbufptr = mk_alleles[2 * marker_uidx + 1];
		} else {
		  *wbufptr++ = output_missing_geno;
		  *wbufptr = output_missing_geno;
		}
	      } else {
		*wbufptr++ = mk_alleles[2 * marker_uidx];
		*wbufptr = mk_alleles[2 * marker_uidx];
	      }
	      bufptr = &(bufptr[unfiltered_indiv_ct4]);
	      wbufptr = &(wbufptr[2]);
	      marker_uidx++;
	    }
	    if (fwrite_checked(writebuf, marker_ct * 3, *outfile_ptr)) {
	      goto recode_ret_WRITE_FAIL;
	    }
	  } else if (max_marker_allele_len == 1) {
	    for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
	      marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx);
	      ucc = ((*bufptr) >> shiftval) & 3;
	      if (ucc) {
		if (ucc == 2) {
		  ucc2 = is_set(marker_exclude, marker_uidx);
		  *wbufptr = mk_alleles[2 * marker_uidx + ucc2];
		  wbufptr[2] = mk_alleles[2 * marker_uidx + (1 ^ ucc2)];
		} else if (ucc == 3) {
		  *wbufptr = mk_alleles[2 * marker_uidx + 1];
		  wbufptr[2] = mk_alleles[2 * marker_uidx + 1];
		} else {
		  *wbufptr = output_missing_geno;
		  wbufptr[2] = output_missing_geno;
		}
	      } else {
		*wbufptr = mk_alleles[2 * marker_uidx];
		wbufptr[2] = mk_alleles[2 * marker_uidx];
	      }
	      bufptr = &(bufptr[unfiltered_indiv_ct4]);
	      wbufptr = &(wbufptr[4]);
	      marker_uidx++;
	    }
	    if (fwrite_checked(writebuf, marker_ct * 4, *outfile_ptr)) {
	      goto recode_ret_WRITE_FAIL;
	    }
	  } else {
	    for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
	      marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx);
	      ucc = ((*bufptr) >> shiftval) & 3;
	      if (ucc) {
		if (ucc == 2) {
		  ucc2 = is_set(marker_exclude, marker_uidx);
		  aptr = &(mk_alleles[(2 * marker_uidx + ucc2) * max_marker_allele_len]);
		  alen = strlen(aptr);
		  aptr2 = &(mk_alleles[(2 * marker_uidx + (1 ^ ucc2)) * max_marker_allele_len]);
		  alen2 = strlen(aptr2);
		  memcpy(wbufptr, aptr, alen);
		  wbufptr[alen] = delim2;
		  memcpy(&(wbufptr[alen + 1]), aptr2, alen2);
		  wbufptr[alen + alen2 + 1] = delimiter;
		  wbufptr = &(wbufptr[alen + alen2 + 2]);
		} else if (ucc == 3) {
		  aptr = &(mk_alleles[(2 * marker_uidx + 1) * max_marker_allele_len]);
		  alen = strlen(aptr);
		  memcpy(wbufptr, aptr, alen);
		  wbufptr[alen] = delim2;
		  memcpy(&(wbufptr[alen + 1]), aptr, alen);
		  wbufptr[2 * alen + 1] = delimiter;
		  wbufptr = &(wbufptr[2 * alen + 2]);
		} else {
		  *wbufptr++ = output_missing_geno;
		  *wbufptr++ = delim2;
		  *wbufptr++ = output_missing_geno;
		  *wbufptr++ = delimiter;
		}
	      } else {
		aptr = &(mk_alleles[2 * marker_uidx * max_marker_allele_len]);
		alen = strlen(aptr);
		memcpy(wbufptr, aptr, alen);
		wbufptr[alen] = delim2;
		memcpy(&(wbufptr[alen + 1]), aptr, alen);
		wbufptr[2 * alen + 1] = delimiter;
		wbufptr = &(wbufptr[2 * alen + 2]);
	      }
	      bufptr = &(bufptr[unfiltered_indiv_ct4]);
	      marker_uidx++;
	    }
	    wbufptr[-1] = '\n';
	    if (fwrite_checked(writebuf, (uintptr_t)(wbufptr - writebuf), *outfile_ptr)) {
	      goto recode_ret_WRITE_FAIL;
	    }
	  }
	  indiv_uidx++;
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
      logprint("Error: --recode does not yet support multipass recoding of very large .bed\nfiles; contact the WDIST developers if you need this.\n");
      retval = RET_CALC_NOT_YET_SUPPORTED;
      goto recode_ret_1;
    }
  }
  fclose_null(outfile_ptr);

  if (recode_modifier & (RECODE_TRANSPOSE | RECODE_LGEN | RECODE_LGEN_REF | RECODE_RLIST)) {
    if (recode_modifier & RECODE_TRANSPOSE) {
      strcpy(outname_end, ".tfam");
    } else {
      strcpy(outname_end, ".fam");
      if ((recode_modifier & (RECODE_TAB | RECODE_DELIMX)) == RECODE_TAB) {
	delimiter = ' ';
      }
    }
    if (fopen_checked(outfile_ptr, outname, "w")) {
      goto recode_ret_OPEN_FAIL;
    }
    outfile = *outfile_ptr;
    rewind(famfile);
    for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ct; indiv_uidx++) {
      if (fgets(tbuf, MAXLINELEN, famfile) == NULL) {
	goto recode_ret_READ_FAIL;
      }
      if (is_set(indiv_exclude, indiv_uidx)) {
	continue;
      }
      wbufptr = tbuf;
      for (loop_end = 0; loop_end < 5; loop_end++) {
	cptr = item_endnn(wbufptr);
	ulii = 1 + (uintptr_t)(cptr - wbufptr);
	*cptr = delimiter;
	if (fwrite_checked(wbufptr, ulii, outfile)) {
	  goto recode_ret_WRITE_FAIL;
	}
	wbufptr = skip_initial_spaces(cptr);
      }
      if (!is_set(pheno_nm, indiv_uidx)) {
	fputs(output_missing_pheno, outfile);
      } else if (affection) {
	putc(is_set(pheno_c, indiv_uidx) + '1', outfile);
      } else {
	if (fprintf(outfile, "%g", pheno_d[indiv_uidx]) < 0) {
	  goto recode_ret_WRITE_FAIL;
	}
      }
      if (putc('\n', outfile) == EOF) {
	goto recode_ret_WRITE_FAIL;
      }
    }
    *outfile_ptr = outfile;
    fclose_null(outfile_ptr);
  }

  if (!(recode_modifier & (RECODE_TRANSPOSE | RECODE_23 | RECODE_A | RECODE_AD | RECODE_LIST))) {
    strcpy(outname_end, ".map");
    if (fopen_checked(outfile_ptr, outname, "w")) {
      goto recode_ret_OPEN_FAIL;
    }
    rewind(bimfile);
    if ((recode_modifier & (RECODE_TAB | RECODE_DELIMX)) == RECODE_DELIMX) {
      cc = ' ';
    } else {
      // PLINK does not use space delimiter here
      cc = '\t';
    }
    for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
      if (get_next_noncomment(bimfile, &wbufptr)) {
	goto recode_ret_READ_FAIL;
      }
      if (is_set(marker_exclude, marker_uidx)) {
	continue;
      }
      wbufptr = tbuf;
      for (loop_end = 0; loop_end < 3; loop_end++) {
	cptr = item_endnn(wbufptr);
	ulii = 1 + (uintptr_t)(cptr - wbufptr);
	*cptr = cc;
	if (fwrite_checked(wbufptr, ulii, *outfile_ptr)) {
	  goto recode_ret_WRITE_FAIL;
	}
	wbufptr = skip_initial_spaces(cptr);
      }
      ulii = strlen_se(wbufptr);
      wbufptr[ulii] = '\n';
      if (fwrite_checked(wbufptr, ulii + 1, *outfile_ptr)) {
	goto recode_ret_WRITE_FAIL;
      }
    }
    fclose_null(outfile_ptr);
  }
  if (!(recode_modifier & RECODE_23)) {
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
  recode_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
 recode_ret_1:
  wkspace_reset(wkspace_mark);
  fclose_cond(ref_file);
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
  double gd_val;
  char* allele[2];
  char idstr[];
} Ll_entry2;

static inline int32_t idmatch(char* idtab, char* id0, uint32_t id0_len_p1, char* id1, uint32_t id1_len_p1) {
  return (!(memcmp(idtab, id0, id0_len_p1) || memcmp(&(idtab[id0_len_p1]), id1, id1_len_p1)));
}

static inline uint32_t hashval(char* id1, uint32_t id1_len, char* id2, uint32_t id2_len) {
  // just interpret as little-endian number and take modulo HASHSIZE_S
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

int top_alloc_str(uintptr_t* topsize_ptr, char* ss, uint32_t slen, char** new_alloc_ptr) {
  if ((slen == 1) && (*((unsigned char*)ss) > ' ')) {
    *new_alloc_ptr = &(one_char_strs[2 * ((unsigned char)(*ss) - 32)]);
    return 0;
  } else {
    *new_alloc_ptr = (char*)top_alloc(topsize_ptr, slen + 1);
    if (!(*new_alloc_ptr)) {
      return -1;
    }
    memcpy(*new_alloc_ptr, ss, slen);
    (*new_alloc_ptr)[slen] = '\0';
    return 0;
  }
}

static inline Ll_entry* top_alloc_ll(uintptr_t* topsize_ptr, uint32_t size) {
  return (Ll_entry*)top_alloc(topsize_ptr, size + sizeof(Ll_entry));
}

static inline Ll_entry2* top_alloc_ll2(uintptr_t* topsize_ptr, uint32_t size) {
  return (Ll_entry2*)top_alloc(topsize_ptr, size + sizeof(Ll_entry2));
}

int32_t merge_fam_id_scan(char* bedname, char* famname, uintptr_t* max_person_id_len_ptr, uint32_t* max_person_full_len_ptr, int32_t* is_dichot_pheno_ptr, Ll_entry** htable, uintptr_t* topsize_ptr, uint64_t* tot_indiv_ct_ptr, uint32_t* ped_buflen_ptr, uint32_t* cur_indiv_ct_ptr, uint32_t* orig_idx_ptr) {
  uintptr_t max_person_id_len = *max_person_id_len_ptr;
  uint32_t max_person_full_len = *max_person_full_len_ptr;
  int32_t is_dichot_pheno = *is_dichot_pheno_ptr;
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
      if (sscanf(col6_start_ptr, "%lg", &pheno) != 1) {
	pheno = -9;
      }
      while (ll_ptr) {
	if (idmatch(ll_ptr->idstr, col1_start_ptr, col1_len + 1, col2_start_ptr, col2_len + 1)) {
	  uii = 0;
	  /*
	  // possibly for future: add parental ID/sex merge (not in PLINK)
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
	memcpy(ll_ptr->idstr, col1_start_ptr, col1_len);
	ll_ptr->idstr[col1_len] = '\t';
	memcpy(&(ll_ptr->idstr[col1_len + 1]), col2_start_ptr, col2_len);
	tot_len = col1_len + col2_len + 1;
	ll_ptr->idstr[tot_len++] = '\t';
	memcpy(&(ll_ptr->idstr[tot_len]), col3_start_ptr, col3_len);
        ll_ptr->idstr[tot_len + col3_len] = '\t';
	tot_len += col3_len + 1;
	memcpy(&(ll_ptr->idstr[tot_len]), col4_start_ptr, col4_len);
	ll_ptr->idstr[tot_len + col4_len] = '\t';
	tot_len += col4_len + 1;
	memcpy(&(ll_ptr->idstr[tot_len]), col5_start_ptr, 1);
	ll_ptr->idstr[tot_len + 1] = '\0';
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
	if (ulii > 2147483647) {
	  sprintf(logbuf, "Error: .ped line too long (>= 2GB) in %s for 32-bit WDIST.\n", famname);
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
    sprintf(logbuf, "Error: No %s in %s.\n", species_plural, famname);
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

int32_t merge_bim_scan(char* bimname, uint32_t is_binary, uintptr_t* max_marker_id_len_ptr, uintptr_t* max_marker_allele_len_ptr, Ll_entry2** htable2, uintptr_t* topsize_ptr, uint64_t* tot_marker_ct_ptr, uint32_t* cur_marker_ct_ptr, int32_t species) {
  uintptr_t max_marker_id_len = *max_marker_id_len_ptr;
  uintptr_t max_marker_allele_len = *max_marker_allele_len_ptr;
  uintptr_t topsize = *topsize_ptr;
  uint64_t tot_marker_ct = *tot_marker_ct_ptr;
  uint32_t cur_marker_ct = *cur_marker_ct_ptr;
  FILE* infile = NULL;
  int32_t retval = 0;
  uint32_t alen1 = 1;
  uint32_t alen2 = 1;
  uint32_t gd_col;
  int64_t llxx;
  char* bufptr;
  char* bufptr2;
  char* bufptr3;
  int32_t ii;
  int32_t jj;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  Ll_entry2** ll_pptr;
  Ll_entry2* ll_ptr;
  double gd_val;
  char* aptr1;
  char* aptr2;
  char* new_aptr;
  uint32_t allele_ct;
  char* cur_alleles[2];
  if (fopen_checked(&infile, bimname, "r")) {
    goto merge_bim_scan_ret_OPEN_FAIL;
  }
  if (check_gd_col(infile, tbuf, is_binary, &gd_col)) {
    goto merge_bim_scan_ret_INVALID_FORMAT_2;
  }
  do {
    bufptr = skip_initial_spaces(tbuf);
    if (is_eoln_or_comment(*bufptr)) {
      continue;
    }
    ii = marker_code(species, bufptr);
    if (ii == -1) {
      sprintf(logbuf, "Error: Invalid chromosome index in %s.\n", bimname);
      goto merge_bim_scan_ret_INVALID_FORMAT;
    }
    bufptr = next_item(bufptr);
    bufptr2 = item_endl(bufptr);
    uii = bufptr2 - bufptr;
    bufptr2 = skip_initial_spaces(bufptr2);
    if (no_more_items_kns(bufptr2)) {
      goto merge_bim_scan_ret_INVALID_FORMAT_2;
    }
    if (gd_col) {
      if (sscanf(bufptr2, "%lg", &gd_val) != 1) {
	gd_val = 0;
      }
      bufptr2 = next_item(bufptr2);
      if (no_more_items_kns(bufptr2)) {
	goto merge_bim_scan_ret_INVALID_FORMAT_2;
      }
    }
    jj = atoi(bufptr2);
    if (jj > 0) {
      if (is_binary) {
	aptr1 = next_item(bufptr2);
	aptr2 = next_item(aptr1);
	if (no_more_items_kns(aptr2)) {
	  goto merge_bim_scan_ret_INVALID_FORMAT_2;
	}
	alen1 = strlen_se(aptr1);
	alen2 = strlen_se(aptr2);
	if (alen1 > max_marker_allele_len) {
	  max_marker_allele_len = alen1;
	}
	if (alen2 > max_marker_allele_len) {
	  max_marker_allele_len = alen2;
	}
	if ((alen1 == 1) && (*aptr1 == '0')) {
	  aptr1 = NULL;
	}
	if (aptr1 && (alen1 == alen2) && (!memcmp(aptr1, aptr2, alen1))) {
	  sprintf(logbuf, "Error: A1 and A2 alleles identical for a marker in %s.\n", bimname);
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
      ukk = 1; // no match?
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
		  goto merge_bim_scan_ret_INVALID_FORMAT_3;
		}
		if (top_alloc_str(&topsize, aptr2, alen2, &new_aptr)) {
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
	    if (aptr1) {
	      for (ukk = 0; ukk < allele_ct; ukk++) {
		if ((!memcmp(aptr1, cur_alleles[ukk], alen1)) && (!cur_alleles[ukk][alen1])) {
		  break;
		}
	      }
	      if (ukk == allele_ct) {
		if (allele_ct == 2) {
		  goto merge_bim_scan_ret_INVALID_FORMAT_3;
		}
		if (top_alloc_str(&topsize, aptr1, alen1, &new_aptr)) {
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
	  if (ll_ptr->pos != llxx) {
	    if ((ll_ptr->pos >> 32) == (llxx >> 32)) {
	      sprintf(logbuf, "Warning: Multiple positions seen for marker %s.\n", bufptr);
	    } else {
	      sprintf(logbuf, "Warning: Multiple chromosomes seen for marker %s.\n", bufptr);
	    }
	    logprintb();
	  }
	  ukk = 0;
	  break;
	}
        ll_pptr = &(ll_ptr->next);
	ll_ptr = *ll_pptr;
      }
      if (ukk) {
        if (uii > max_marker_id_len) {
	  max_marker_id_len = uii;
	}
	ll_ptr = top_alloc_ll2(&topsize, uii);
	if (!ll_ptr) {
	  goto merge_bim_scan_ret_NOMEM;
	}
	ll_ptr->next = NULL;
	ll_ptr->pos = llxx;
	ll_ptr->gd_val = gd_val;
	if (aptr1) {
	  if (top_alloc_str(&topsize, aptr1, alen1, &(ll_ptr->allele[0]))) {
	    goto merge_bim_scan_ret_NOMEM;
	  }
	} else {
	  ll_ptr->allele[0] = NULL;
	}
	if (aptr2) {
	  if (top_alloc_str(&topsize, aptr2, alen2, &(ll_ptr->allele[1]))) {
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
    } else if (!jj) {
      sprintf(logbuf, "Error: Invalid base-pair position in %s.\n", bimname);
      goto merge_bim_scan_ret_INVALID_FORMAT;
    }
  } while (fgets(tbuf, MAXLINELEN, infile));
  if (!feof(infile)) {
    goto merge_bim_scan_ret_READ_FAIL;
  }
  if (!cur_marker_ct) {
    sprintf(logbuf, "Error: No markers in %s.\n", bimname);
    goto merge_bim_scan_ret_INVALID_FORMAT;
  }
  *max_marker_allele_len_ptr = max_marker_allele_len;
  *max_marker_id_len_ptr = max_marker_id_len;
  *topsize_ptr = topsize;
  *tot_marker_ct_ptr = tot_marker_ct;
  *cur_marker_ct_ptr = cur_marker_ct;
  
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
  merge_bim_scan_ret_INVALID_FORMAT_3:
    sprintf(logbuf, "Error: Marker %s is not biallelic.\n", bufptr);
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  merge_bim_scan_ret_INVALID_FORMAT_2:
    sprintf(logbuf, "Error: %s is not a properly formatted .map/.bim file.\n", bimname);
  merge_bim_scan_ret_INVALID_FORMAT:
    logprintb();
    retval = RET_INVALID_FORMAT;
  }
  fclose_cond(infile);
  return retval;
}

static inline void merge_post_msort_update_maps(char* marker_ids, uintptr_t max_marker_id_len, uint32_t* marker_map, uint32_t* pos_buf, int64_t* ll_buf, uint32_t* chrom_start, uint32_t* chrom_id, uint32_t chrom_ct, uint32_t* dedup_marker_ct_ptr, uint32_t merge_disallow_equal_pos, uint64_t chrom_mask) {
  // Input: ll_buf is a effectively sequence of sorted arrays (one per
  // chromosome) with base-pair positions in high 32 bits, and pre-sort indices
  // in low 32 bits.  Chromosome boundaries are stored in chrom_start[].
  // Pre-sort indices refer to marker ID ASCII ordering, which is how lookup
  // will be performed mid-merge.  There may be duplicate positions, and
  // markers that don't pass the chromosome filter.

  // Result: Duplicates have been collapsed, with chrom_start[] updated.
  // pos_buf contains sorted base-pair positions,
  // post-chromosome-filtering-and-duplicate-removal.
  // marker_map[n] is the post-filtering position in all other arrays of marker
  // ID n by ASCII ordering.
  uint32_t read_pos = 0;
  uint32_t write_pos = 0; // may be lower than read_pos due to dups
  uint32_t chrom_idx;
  uint32_t chrom_read_end_idx;
  int64_t llxx;
  uint32_t prev_bp;
  uint32_t cur_bp;
  uint32_t presort_idx;
  for (chrom_idx = 0; chrom_idx < chrom_ct; chrom_idx++) {
    if (!((chrom_mask >> chrom_id[chrom_idx]) & 1LLU)) {
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
	sprintf(logbuf, "Warning: Markers %s and %s have the same position.\n", &(marker_ids[max_marker_id_len * presort_idx]), &(marker_ids[max_marker_id_len * ((uint32_t)ll_buf[read_pos - 1])]));
	logprintb();
	if (merge_disallow_equal_pos) {
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
}

static inline int32_t merge_must_track_write(int32_t mm) {
  // modes 6 and 7 can be sped up with early pruning of nonoverlapping
  // markers, but not worth complicating code for this.
  return (mm == 1) || (mm > 5) || (mm == 4);
}

static inline int32_t merge_first_mode(int32_t mm, uint32_t merge_disallow_equal_pos) {
  if (merge_disallow_equal_pos) {
    return (mm > 5)? 4 : mm;
  } else {
    return merge_must_track_write(mm)? 4 : 5;
  }
}

int32_t merge_diff_print(FILE* outfile, char* idbuf, char* marker_id, char* person_id, unsigned char newval, unsigned char oldval, char* marker_alleles, uintptr_t max_marker_allele_len) {
  char* bufptr = item_endnn(person_id);
  uint32_t slen = strlen_se(marker_id);
  char ma1[4];
  char ma2[4];
  char* ma1p[4];
  char* ma2p[4];
  char wbuf[8];
  uint32_t slen2;
  memcpy(idbuf, marker_id, slen);
  idbuf[slen] = '\0';
  if (fprintf(outfile, "%20s", idbuf) < 0) {
    return -1;
  }
  slen = (bufptr++) - person_id;
  memcpy(idbuf, person_id, slen);
  idbuf[slen] = '\0';
  if (fprintf(outfile, " %20s %20s", idbuf, bufptr) < 0) {
    return -1;
  }
  if (max_marker_allele_len == 1) {
    ma1[0] = marker_alleles[0];
    ma1[1] = '0';
    ma1[2] = marker_alleles[0];
    ma1[3] = marker_alleles[1];
    ma2[0] = marker_alleles[0];
    ma2[1] = '0';
    ma2[2] = marker_alleles[1];
    ma2[3] = marker_alleles[1];
    return (fprintf(outfile, "      %c/%c      %c/%c \n", ma1[newval], ma2[newval], ma1[oldval], ma2[oldval]) < 0);
  }
  ma1p[0] = marker_alleles;
  ma1p[1] = &(one_char_strs[32]);
  ma1p[2] = marker_alleles;
  ma1p[3] = &(marker_alleles[max_marker_allele_len]);
  ma2p[0] = marker_alleles;
  ma2p[1] = &(one_char_strs[32]);
  ma2p[2] = &(marker_alleles[max_marker_allele_len]);
  ma2p[3] = &(marker_alleles[max_marker_allele_len]);
  slen = strlen(ma1p[newval]);
  slen2 = strlen(ma2p[newval]);
  if (slen + slen2 > 6) {
    if (fprintf(outfile, " %s/%s", ma1p[newval], ma2p[newval]) < 0) {
      return -1;
    }
  } else {
    memcpy(wbuf, ma1p[newval], slen);
    wbuf[slen] = '/';
    memcpy(&(wbuf[slen + 1]), ma2p[newval], slen2 + 1);
    if (fprintf(outfile, " %8s", wbuf) < 0) {
      return -1;
    }
  }
  slen = strlen(ma1p[oldval]);
  slen2 = strlen(ma2p[oldval]);
  if (slen + slen2 > 6) {
    return (fprintf(outfile, " %s/%s \n", ma1p[oldval], ma2p[oldval]) < 0);
  } else {
    memcpy(wbuf, ma1p[oldval], slen);
    wbuf[slen] = '/';
    memcpy(&(wbuf[slen + 1]), ma2p[oldval], slen2 + 1);
    return (fprintf(outfile, " %8s \n", wbuf) < 0);
  }
}

int32_t merge_main(char* bedname, char* bimname, char* famname, uint32_t tot_indiv_ct, uint32_t tot_marker_ct, uint32_t dedup_marker_ct, uint32_t start_marker_idx, uint32_t marker_window_size, char* marker_alleles, uintptr_t max_marker_allele_len, char* marker_ids, uintptr_t max_marker_id_len, char* person_ids, uintptr_t max_person_id_len, uint32_t merge_nsort, uint32_t* indiv_nsmap, uint32_t* flex_map, uint32_t* marker_map, uint32_t* chrom_start, uint32_t* chrom_id, uint32_t chrom_ct, char* idbuf, unsigned char* readbuf, unsigned char* writebuf, uint32_t merge_mode, uintptr_t* markbuf, FILE* outfile, uint64_t* diff_total_overlap_ptr, uint64_t* diff_not_both_genotyped_ptr, uint64_t* diff_discordant_ptr, uint32_t ped_buflen, unsigned char* bmap_raw) {
  // flex_map maps individuals for binary filesets, and markers for text
  // filesets.
  uint32_t is_binary = famname? 1 : 0;
  FILE* bedfile = NULL;
  FILE* infile2 = NULL;
  int32_t retval = 0;
  uint32_t tot_indiv_ct4 = (tot_indiv_ct + 3) / 4;
  uint32_t tot_indiv_ctl = (tot_indiv_ct + (BITCT - 1)) / BITCT;
  uint32_t end_marker_idx = start_marker_idx + marker_window_size;
  uint32_t marker_in_idx = 4294967295U; // overflow to zero on first add
  uint32_t last_marker_in_idx = 4294967294U;
  uint32_t cur_indiv_ct = 0;
  uintptr_t* mbufptr = NULL; // merge mode 1, 4, 6, 7
  uint64_t diff_total_overlap = 0;
  uint64_t diff_not_both_genotyped = 0;
  uint64_t diff_discordant = 0;
  uint32_t cur_indiv_ct4 = 0;
  uint32_t cur_indiv_ctd4 = 0;
  uintptr_t uljj = 0;
  uintptr_t* mbufptr2;
  unsigned char* bmap;
  uint32_t gd_col;
  char* bufptr;
  char* bufptr2;
  char* bufptr3;
  char* bufptr4;
  unsigned char* rbufptr;
  unsigned char* wbufptr;
  unsigned char* wbufptr2;
  uintptr_t ulii;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  uint32_t unn;
  uint32_t alen1;
  uint32_t alen2;
  uintptr_t indiv_idx;
  uint32_t marker_out_idx;
  unsigned char ucc;
  unsigned char ucc2;
  unsigned char ucc3;
  unsigned char ucc4;
  char cc;
  char cc2;
  char cc3;
  char cc4;
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
      memcpy(idbuf, bufptr, uii);
      idbuf[uii] = '\t';
      memcpy(&(idbuf[uii + 1]), bufptr3, ujj);
      idbuf[uii + ujj + 1] = '\0';
      if (merge_nsort) {
        ii = bsearch_str_natural(idbuf, person_ids, max_person_id_len, 0, tot_indiv_ct - 1);
      } else {
	ii = bsearch_str(idbuf, person_ids, max_person_id_len, 0, tot_indiv_ct - 1);
	if (indiv_nsmap && (ii != -1)) {
	  ii = indiv_nsmap[ii];
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
  }
  if (fopen_checked(&infile2, bimname, "r")) {
    goto merge_main_ret_OPEN_FAIL;
  }
  if (check_gd_col(infile2, tbuf, is_binary, &gd_col)) {
    goto merge_main_ret_READ_FAIL;
  }
  if (fopen_checked(&bedfile, bedname, is_binary? "rb" : "r")) {
    goto merge_main_ret_OPEN_FAIL;
  }
  do {
    bufptr = skip_initial_spaces(tbuf);
    if (is_eoln_or_comment(*bufptr)) {
      continue;
    }
    ++marker_in_idx;
    bufptr = next_item(bufptr);
    bufptr2 = next_item_mult(bufptr, 1 + gd_col);
    if (!bufptr2) {
      goto merge_main_ret_READ_FAIL;
    }
    if (*bufptr2 == '-') {
      if (!is_binary) {
	flex_map[marker_in_idx] = 4294967295U;
      }
      continue;
    }
    bufptr3 = item_endnn(bufptr);
    uii = (bufptr3 - bufptr);
    memcpy(idbuf, bufptr, uii);
    idbuf[uii] = '\0';
    ii = bsearch_str(idbuf, marker_ids, max_marker_id_len, 0, tot_marker_ct - 1);
    if (ii == -1) {
      goto merge_main_ret_READ_FAIL;
    }
    marker_out_idx = marker_map[ii];
    if ((marker_out_idx < start_marker_idx) || (marker_out_idx >= end_marker_idx)) {
      if (!is_binary) {
	flex_map[marker_in_idx] = 4294967295U;
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
      if (max_marker_allele_len == 1) {
	cc = *bufptr2;
	cc2 = *bufptr3;
	if (((cc != '0') && (cc == marker_alleles[ii * 2 + 1])) || ((cc2 == marker_alleles[ii * 2]) && (cc2 != '0'))) {
	  bmap = &(bmap_raw[256]);
	} else {
	  bmap = bmap_raw;
	}
      } else {
	alen1 = strlen_se(bufptr2);
	alen2 = strlen_se(bufptr3);
	bufptr4 = &(marker_alleles[((uintptr_t)(ii * 2)) * max_marker_allele_len]);
	if ((((*bufptr2 != '0') || (alen1 != 1)) && ((!memcmp(bufptr2, &(bufptr4[max_marker_allele_len]), alen1)) && (!bufptr4[max_marker_allele_len + alen1]))) || (((*bufptr3 != '0') || (alen2 != 1)) && ((!memcmp(bufptr3, bufptr4, alen2)) && (!bufptr4[alen2])))) {
	  bmap = &(bmap_raw[256]);
	} else {
	  bmap = bmap_raw;
	}
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
		if (merge_diff_print(outfile, idbuf, bufptr, &(person_ids[ujj * max_person_id_len]), ucc2, ucc3, &(marker_alleles[((uintptr_t)(ii * 2)) * max_marker_allele_len]), max_marker_allele_len)) {
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
	      if (merge_diff_print(outfile, idbuf, bufptr, &(person_ids[ujj * max_person_id_len]), ucc2, ucc3, &(marker_alleles[((uintptr_t)(ii * 2)) * max_marker_allele_len]), max_marker_allele_len)) {
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
		if (merge_diff_print(outfile, idbuf, bufptr, &(person_ids[ujj * max_person_id_len]), ucc2, ucc3, &(marker_alleles[((uintptr_t)(ii * 2)) * max_marker_allele_len]), max_marker_allele_len)) {
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
	      if (merge_diff_print(outfile, idbuf, bufptr, &(person_ids[ujj * max_person_id_len]), ucc2, ucc3, &(marker_alleles[((uintptr_t)(ii * 2)) * max_marker_allele_len]), max_marker_allele_len)) {
		goto merge_main_ret_WRITE_FAIL;
	      }
	    }
	  }
	  ucc >>= 2;
	}
	break;
      }
    } else {
      flex_map[marker_in_idx] = marker_out_idx - start_marker_idx;
    }
  } while (fgets(tbuf, MAXLINELEN, infile2));
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
      memcpy(idbuf, bufptr, uii);
      idbuf[uii] = '\t';
      bufptr3 = skip_initial_spaces(bufptr2);
      bufptr2 = item_endnn(bufptr3);
      ujj = (bufptr2 - bufptr3);
      memcpy(&(idbuf[uii + 1]), bufptr3, ujj);
      idbuf[uii + ujj + 1] = '\0';
      if (merge_nsort) {
        ii = bsearch_str_natural(idbuf, person_ids, max_person_id_len, 0, tot_indiv_ct - 1);
      } else {
	ii = bsearch_str(idbuf, person_ids, max_person_id_len, 0, tot_indiv_ct - 1);
	if (indiv_nsmap && (ii != -1)) {
	  ii = indiv_nsmap[ii];
	}
      }
      if (ii == -1) {
	goto merge_main_ret_READ_FAIL;
      }
      bufptr3 = next_item_mult(skip_initial_spaces(bufptr2), 4);
      if (!bufptr3) {
	goto merge_main_ret_READ_FAIL;
      }
      wbufptr = &(writebuf[ii / 4]);
      ujj = (ii % 4) * 2;
      ucc = ~(3U << ujj);
      ucc4 = (1U << ujj);
      if (merge_must_track_write(merge_mode)) {
	mbufptr = &(markbuf[ii / BITCT]);
	uljj = ONELU << (ii % BITCT);
      }
      for (marker_in_idx = 0; marker_in_idx < last_marker_in_idx; marker_in_idx++) {
	cc = *bufptr3;
	if (is_eoln_kns(cc)) {
	  goto merge_main_ret_INVALID_FORMAT_2;
	}
	bufptr3 = skip_initial_spaces(&(bufptr3[1]));
	cc2 = *bufptr3;
	if (is_eoln_kns(cc2)) {
	  goto merge_main_ret_INVALID_FORMAT_2;
	}
	bufptr3 = skip_initial_spaces(&(bufptr3[1]));
        uii = flex_map[marker_in_idx];
	if (uii == 4294967295U) {
	  continue;
	}
	if (cc == '0') {
	  if (cc2 != '0') {
	    goto merge_main_ret_INVALID_FORMAT_3;
	  }
          ucc2 = 1;
	} else if (cc2 == '0') {
	  goto merge_main_ret_INVALID_FORMAT_3;
	} else {
	  ucc2 = 0; // A2 count
	  if (max_marker_allele_len == 1) {
	    cc3 = marker_alleles[uii * 2];
	    cc4 = marker_alleles[uii * 2 + 1];
	  } else {
	    if (marker_alleles[uii * 2 * max_marker_allele_len + 1]) {
	      cc3 = ' ';
	    } else {
	      cc3 = marker_alleles[uii * 2 * max_marker_allele_len];
	    }
	    if (marker_alleles[(uii * 2 + 1) * max_marker_allele_len + 1]) {
	      cc4 = ' ';
	    } else {
	      cc4 = marker_alleles[(uii * 2 + 1) * max_marker_allele_len];
	    }
	  }
	  if (cc == cc4) {
	    ucc2++;
	  } else if (cc != cc3) {
	    if (cc4 == '0') {
	      cc4 = cc;
	      ucc2++;
	      marker_alleles[(uii * 2 + 1) * max_marker_allele_len] = cc4;
	    } else if (cc3 == '0') {
	      cc3 = cc;
	      marker_alleles[(uii * 2) * max_marker_allele_len] = cc3;
	    } else {
	      goto merge_main_ret_INVALID_FORMAT;
	    }
	  }
	  if (cc2 == cc4) {
	    ucc2++;
	  } else if (cc2 != cc3) {
	    if (cc3 == '0') {
	      marker_alleles[(uii * 2) * max_marker_allele_len] = cc2;
	    } else if (cc4 == '0') {
	      // put this second since only way cc3 != '0' and cc4 == '0' is if
	      // it was specified that way in earlier binary file, and cc
	      // matches cc3.
	      marker_alleles[(uii * 2 + 1) * max_marker_allele_len] = cc2;
	      ucc2++;
	    } else {
	      goto merge_main_ret_INVALID_FORMAT;
	    }
	  }
	  // transform to PLINK encoding
	  if (ucc2) {
	    ucc2++;
	  }
	}
	switch (merge_mode) {
	case 1:
	  mbufptr2 = &(mbufptr[uii * tot_indiv_ctl]);
	  wbufptr2 = &(wbufptr[uii * tot_indiv_ct4]);
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
	  wbufptr2 = &(wbufptr[uii * tot_indiv_ct4]);
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
	  wbufptr2 = &(wbufptr[uii * tot_indiv_ct4]);
	  *wbufptr2 = ((*wbufptr2) & ucc) | (ucc2 << ujj);
	  break;
	case 4:
	  mbufptr2 = &(mbufptr[uii * tot_indiv_ctl]);
	  if (!((*mbufptr2) & uljj)) {
	    *mbufptr2 |= uljj;
	    wbufptr2 = &(wbufptr[uii * tot_indiv_ct4]);
	    *wbufptr2 = ((*wbufptr2) & ucc) | (ucc2 << ujj);
	  }
	  break;
	case 6:
	  if (mbufptr[uii * tot_indiv_ctl] & uljj) {
	    diff_total_overlap++;
	    ucc3 = ((wbufptr[uii * tot_indiv_ct4] >> ujj) & 3);
	    umm = ((ucc2 == 1) || (ucc3 == 1));
	    if (umm) {
	      diff_not_both_genotyped++;
	    }
	    if (ucc2 != ucc3) {
	      if (!umm) {
		diff_discordant++;
	      }
	      if (merge_diff_print(outfile, idbuf, &(marker_ids[uii * max_marker_id_len]), &(person_ids[((uint32_t)ii) * max_person_id_len]), ucc2, ucc3, &(marker_alleles[uii * 2 * max_marker_allele_len]), max_marker_allele_len)) {
		goto merge_main_ret_WRITE_FAIL;
	      }
	    }
	  }
	  break;
	case 7:
	  if (mbufptr[uii * tot_indiv_ctl] & uljj) {
	    diff_total_overlap++;
	    ucc3 = ((wbufptr[uii * tot_indiv_ct4] >> ujj) & 3);
	    if ((ucc2 == 1) || (ucc3 == 1)) {
	      diff_not_both_genotyped++;
	    } else if (ucc2 != ucc3) {
	      diff_discordant++;
	      if (merge_diff_print(outfile, idbuf, &(marker_ids[uii * max_marker_id_len]), &(person_ids[((uint32_t)ii) * max_person_id_len]), ucc2, ucc3, &(marker_alleles[uii * 2 * max_marker_allele_len]), max_marker_allele_len)) {
		goto merge_main_ret_WRITE_FAIL;
	      }
	    }
	  }
	}
      }
      if (!is_eoln_kns(*bufptr3)) {
	fill_idbuf_fam_indiv(idbuf, bufptr, ' ');
	sprintf(logbuf, "Error: %s line has too many entries (indiv id %s).\nIf this file has multi-character alleles, convert to binary before merging.\n", bedname, idbuf);
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
    sprintf(logbuf, "Error: Marker %s is not biallelic.\n", &(marker_ids[uii * max_marker_id_len]));
    putchar('\n');
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  merge_main_ret_INVALID_FORMAT_3:
    fill_idbuf_fam_indiv(idbuf, bufptr, ' ');
    sprintf(logbuf, "Error: Half-missing call in %s (indiv id %s, marker %u)", bedname, idbuf, marker_in_idx);
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

int32_t merge_datasets(char* bedname, char* bimname, char* famname, char* outname, char* outname_end, char* mergename1, char* mergename2, char* mergename3, uint64_t calculation_type, int32_t merge_type, int32_t indiv_sort, int32_t keep_allele_order, Chrom_info* chrom_info_ptr) {
  FILE* mergelistfile = NULL;
  FILE* outfile = NULL;
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t max_person_id_len = 0;
  uint32_t max_person_full_len = 0;
  uintptr_t max_marker_id_len = 0;
  uintptr_t max_marker_allele_len = 1;
  int32_t is_dichot_pheno = 1;
  int32_t merge_mode = (merge_type & MERGE_MODE_MASK);
  uint32_t merge_nsort = ((!indiv_sort) || (indiv_sort == INDIV_SORT_NATURAL))? 1 : 0;
  uint32_t merge_disallow_equal_pos = (merge_type & MERGE_ALLOW_EQUAL_POS)? 0 : 1;
  Ll_entry** htable = (Ll_entry**)(&(wkspace_base[wkspace_left - HASHMEM_S]));
  Ll_entry2** htable2 = (Ll_entry2**)(&(wkspace_base[wkspace_left - HASHMEM]));
  uint32_t ped_buflen = MAXLINELEN;
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
  uint32_t* map_reverse = NULL;
  uintptr_t* reversed = NULL;
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
  // N.B. marker_alleles are ordered by marker_id instead of position
  char* marker_alleles;
  uint32_t* marker_map;
  uint32_t* flex_map;
  double* gd_vals;
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
  uint32_t cur_marker_ct;
  uint32_t tot_indiv_ct;
  uint32_t tot_indiv_ct4;
  uint32_t tot_marker_ct;
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
  char cc2;
  unsigned char ucc;
  int32_t retval;

  bufptr = one_char_strs;
  for (uii = 0; uii < 224; uii++) {
    *bufptr++ = (char)(uii + ' ');
    *bufptr++ = '\0';
  }

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
    while (fgets(tbuf, MAXLINELEN, mergelistfile)) {
      bufptr = skip_initial_spaces(tbuf);
      if (no_more_items_kns(bufptr)) {
	continue;
      }
      bufptr2 = next_item_mult(bufptr, 3);
      if (!no_more_items_kns(bufptr2)) {
	logprint("Error: More than three items on --merge-list line.\n");
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
    if (ullxx > 2147483647) {
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
        memcpy(bufptr4, bufptr, ulii);
	memcpy(&(bufptr4[ulii]), ".bed", 5);
	mergelist_bed[mlpos] = bufptr4;
	bufptr4 = &(bufptr4[ulii + 5]);
        memcpy(bufptr4, bufptr, ulii);
	memcpy(&(bufptr4[ulii]), ".bim", 5);
	mergelist_bim[mlpos] = bufptr4;
	bufptr4 = &(bufptr4[ulii + 5]);
        memcpy(bufptr4, bufptr, ulii);
	memcpy(&(bufptr4[ulii]), ".fam", 5);
	mergelist_fam[mlpos] = bufptr4;
	bufptr4 = &(bufptr4[ulii + 5]);
      } else {
	memcpy(bufptr4, bufptr, ulii);
	bufptr4[ulii] = '\0';
	mergelist_bed[mlpos] = bufptr4;
	bufptr4 = &(bufptr4[ulii + 1]);
	bufptr2 = item_endnn(bufptr3);
	ulii = bufptr2 - bufptr3;
	bufptr = skip_initial_spaces(bufptr2);
	memcpy(bufptr4, bufptr3, ulii);
	bufptr4[ulii] = '\0';
	mergelist_bim[mlpos] = bufptr4;
	bufptr4 = &(bufptr4[ulii + 1]);
	if (no_more_items_kns(bufptr)) {
	  mergelist_fam[mlpos] = NULL;
	} else {
	  bufptr2 = item_endnn(bufptr);
	  ulii = bufptr2 - bufptr;
	  memcpy(bufptr4, bufptr, ulii);
	  bufptr4[ulii] = '\0';
	  mergelist_fam[mlpos] = bufptr4;
	  bufptr4 = &(bufptr4[ulii + 1]);
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
  if (ullxx > 2147483647) {
    sprintf(logbuf, "Error: Too many %s (max 2147483647).\n", species_plural);
    goto merge_datasets_ret_INVALID_FORMAT_2;
  }
#else
  // avoid integer overflow in wkspace_alloc calls
  if (ullxx * max_person_full_len > 2147483647) {
    sprintf(logbuf, "Error: Too many %s for 32-bit WDIST.\n", species_plural);
    goto merge_datasets_ret_INVALID_FORMAT_2;
  }
#endif
  tot_indiv_ct = ullxx;
  fill_bmap_raw(bmap_raw, tot_indiv_ct % 4);
  bmap = &(bmap_raw[256]);
  bmap2 = &(bmap_raw[512]);
  if (indiv_sort == INDIV_SORT_NONE) {
    if (wkspace_alloc_ui_checked(&indiv_nsmap, tot_indiv_ct * sizeof(int32_t))) {
      goto merge_datasets_ret_NOMEM;
    }
  }
  if (wkspace_alloc_c_checked(&person_ids, max_person_id_len * tot_indiv_ct) ||
      wkspace_alloc_c_checked(&person_fids, max_person_full_len * tot_indiv_ct)) {
    goto merge_datasets_ret_NOMEM;
  }
  if (is_dichot_pheno) {
    if (wkspace_left < topsize + tot_indiv_ct) {
      goto merge_datasets_ret_NOMEM;
    }
    wkspace_alloc_c_checked(&pheno_c_char, tot_indiv_ct * sizeof(char));
  } else {
    if (wkspace_left < topsize + tot_indiv_ct * sizeof(double)) {
      goto merge_datasets_ret_NOMEM;
    }
    wkspace_alloc_d_checked(&pheno_d, tot_indiv_ct * sizeof(double));
  }
  if (indiv_sort == INDIV_SORT_NONE) {
    if (wkspace_alloc_ui_checked(&map_reverse, tot_indiv_ct * sizeof(int32_t))) {
      goto merge_datasets_ret_NOMEM;
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
    wkspace_left -= topsize;
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
    wkspace_left -= topsize;
    if (is_dichot_pheno) {
      if (qsort_ext(person_fids, tot_indiv_ct, max_person_full_len, merge_nsort? strcmp_natural_deref : strcmp_deref, pheno_c_char, 1)) {
	goto merge_datasets_ret_NOMEM2;
      }
    } else {
      if (qsort_ext(person_fids, tot_indiv_ct, max_person_full_len, merge_nsort? strcmp_natural_deref : strcmp_deref, (char*)pheno_d, sizeof(double))) {
	goto merge_datasets_ret_NOMEM2;
      }
    }
    wkspace_left += topsize;
  }
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
      memcpy(bufptr3, bufptr, uii);
      bufptr3[uii] = '\0';
      if (merge_mode < 6) {
	uii += strlen(bufptr2) + 1;
	if (fwrite_checked(bufptr, uii, outfile)) {
	  goto merge_datasets_ret_WRITE_FAIL;
	}
	if (is_dichot_pheno) {
	  cc = pheno_c_char[ulii];
	  if (fprintf(outfile, "\t%s\n", cc? ((cc == 1)? "2" : "-9") : "1") < 0) {
	    goto merge_datasets_ret_WRITE_FAIL;
	  }
	} else {
	  if (fprintf(outfile, "\t%g\n", pheno_d[ulii]) < 0) {
	    goto merge_datasets_ret_WRITE_FAIL;
	  }
	}
      }
    }
  } else {
    bufptr = person_fids;
    bufptr3 = person_ids;
    for (ulii = 0; ulii < tot_indiv_ct; ulii++) {
      bufptr2 = next_item_mult(bufptr, 2);
      uii = (bufptr2 - bufptr) - 1;
      memcpy(bufptr3, bufptr, uii);
      bufptr3[uii] = '\0';
      bufptr3 = &(bufptr3[max_person_id_len]);
      if (merge_mode < 6) {
	uii += strlen(bufptr2) + 1;
	if (fwrite_checked(bufptr, uii, outfile)) {
	  goto merge_datasets_ret_WRITE_FAIL;
	}
	if (is_dichot_pheno) {
	  cc = pheno_c_char[ulii];
	  if (fprintf(outfile, "\t%s\n", cc? ((cc == 1)? "2" : "-9") : "1") < 0) {
	    goto merge_datasets_ret_WRITE_FAIL;
	  }
	} else {
	  if (fprintf(outfile, "\t%g\n", pheno_d[ulii]) < 0) {
	    goto merge_datasets_ret_WRITE_FAIL;
	  }
	}
      }
      bufptr = &(bufptr[max_person_full_len]);
    }
  }
  if (merge_mode < 6) {
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
    retval = merge_bim_scan(mergelist_bim[mlpos], (mergelist_fam[mlpos])? 1 : 0, &max_marker_id_len, &max_marker_allele_len, htable2, &topsize, &ullxx, &cur_marker_ct, chrom_info_ptr->species);
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
  if (ullxx > 2147483647) {
    logprint("Error: Too many markers (max 2147483647).\n");
    goto merge_datasets_ret_INVALID_FORMAT;
  }
#else
  if (ullxx * MAXV(max_marker_id_len, 8) > 2147483647) {
    logprint("Error: Too many markers for 32-bit WDIST.\n");
    goto merge_datasets_ret_INVALID_FORMAT;
  }
#endif
  tot_marker_ct = ullxx;
  wkspace_left -= topsize;
  if (max_marker_allele_len > 1) {
    max_marker_allele_len++;
  }
  if (wkspace_alloc_c_checked(&marker_ids, max_marker_id_len * tot_marker_ct) ||
      wkspace_alloc_c_checked(&marker_alleles, (tot_marker_ct * max_marker_allele_len) * 2) ||
      wkspace_alloc_ui_checked(&marker_map, tot_marker_ct * sizeof(int32_t)) ||
      wkspace_alloc_d_checked(&gd_vals, tot_marker_ct * sizeof(double)) ||
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
  wkspace_left += topsize;
  ulii = 0;
  if (max_marker_allele_len == 1) {
    for (uii = 0; uii < HASHSIZE; uii++) {
      if (htable2[uii]) {
	ll_ptr2 = htable2[uii];
	do {
	  ujj = marker_map[ulii++];
	  llxx = ll_ptr2->pos;
	  pos_buf[ujj] = (uint32_t)llxx;
	  bufptr = ll_ptr2->allele[0];
	  marker_alleles[ujj * 2] = bufptr? (*bufptr) : '0';
	  bufptr = ll_ptr2->allele[1];
	  marker_alleles[ujj * 2 + 1] = bufptr? (*bufptr) : '0';
	  gd_vals[ujj] = ll_ptr2->gd_val;
	  ll_buf[ujj] = (llxx & 0xffffffff00000000LL) | ujj;
	  ll_ptr2 = ll_ptr2->next;
	} while (ll_ptr2);
      }
    }
  } else {
    for (uii = 0; uii < HASHSIZE; uii++) {
      if (htable2[uii]) {
	ll_ptr2 = htable2[uii];
	do {
	  ujj = marker_map[ulii++];
	  llxx = ll_ptr2->pos;
	  pos_buf[ujj] = (uint32_t)llxx;
	  bufptr = ll_ptr2->allele[0];
	  if (bufptr) {
	    strcpy(&(marker_alleles[(ujj * 2) * max_marker_allele_len]), bufptr);
	  } else {
	    marker_alleles[(ujj * 2) * max_marker_allele_len] = '0';
	    marker_alleles[(ujj * 2) * max_marker_allele_len + 1] = '\0';
	  }
	  bufptr = ll_ptr2->allele[1];
	  if (bufptr) {
	    strcpy(&(marker_alleles[(ujj * 2 + 1) * max_marker_allele_len]), bufptr);
	  } else {
	    marker_alleles[(ujj * 2 + 1) * max_marker_allele_len] = '0';
	    marker_alleles[(ujj * 2 + 1) * max_marker_allele_len + 1] = '\0';
	  }
	  gd_vals[ujj] = ll_ptr2->gd_val;
	  ll_buf[ujj] = (llxx & 0xffffffff00000000LL) | ujj;
	  ll_ptr2 = ll_ptr2->next;
	} while (ll_ptr2);
      }
    }
  }
  sort_marker_chrom_pos(ll_buf, tot_marker_ct, (int32_t*)pos_buf, chrom_start, chrom_id, &chrom_ct);
  merge_post_msort_update_maps(marker_ids, max_marker_id_len, marker_map, pos_buf, ll_buf, chrom_start, chrom_id, chrom_ct, &dedup_marker_ct, merge_disallow_equal_pos, chrom_info_ptr->chrom_mask);
  if (!dedup_marker_ct) {
    logprint("Error: No markers in merged file.\n");
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
      uii = dedup_marker_ct;
    } else {
      uii = markers_per_pass;
    }
    markbuf = (uintptr_t*)wkspace_alloc(uii * ulii * sizeof(intptr_t));
  } else {
    markers_per_pass = wkspace_left / tot_indiv_ct4;
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
      sprintf(logbuf, "Performing single-pass merge (%u %s, %u marker%s).\n", tot_indiv_ct, species_str(tot_indiv_ct), dedup_marker_ct, (dedup_marker_ct == 1)? "" : "s");
    } else {
      sprintf(logbuf, "Performing %u-pass merge (%u %s, %" PRIuPTR "/%u marker%s per pass).\n", pass_ct, tot_indiv_ct, species_str(tot_indiv_ct), markers_per_pass, dedup_marker_ct, (dedup_marker_ct == 1)? "" : "s");
    }
  } else {
    memcpy(outname_end, ".diff", 6);
    if (fopen_checked(&outfile, outname, "w")) {
      goto merge_datasets_ret_OPEN_FAIL;
    }
    if (fputs("                 SNP                  FID                  IID      NEW      OLD \n", outfile) == EOF) {
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
      retval = merge_main(mergelist_bed[mlpos], mergelist_bim[mlpos], mergelist_fam[mlpos], tot_indiv_ct, tot_marker_ct, dedup_marker_ct, uii * markers_per_pass, ujj, marker_alleles, max_marker_allele_len, marker_ids, max_marker_id_len, person_ids, max_person_id_len, merge_nsort, indiv_nsmap, flex_map, marker_map, chrom_start, chrom_id, chrom_ct, idbuf, readbuf, writebuf, mlpos? merge_mode : merge_first_mode(merge_mode, merge_disallow_equal_pos), markbuf, outfile, &diff_total_overlap, &diff_not_both_genotyped, &diff_discordant, ped_buflen, bmap_raw);
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
    if (max_marker_allele_len == 1) {
      for (ulii = 0; ulii < chrom_ct; ulii++) {
	uii = chrom_start[ulii + 1];
	ujj = chrom_start[ulii];
	ukk = chrom_id[ulii];
	for (; ujj < uii; ujj++) {
	  umm = map_reverse[ujj];
	  if (keep_allele_order || (!is_set(reversed, ujj))) {
	    cc = marker_alleles[2 * umm];
	    cc2 = marker_alleles[2 * umm + 1];
	  } else {
	    cc = marker_alleles[2 * umm + 1];
	    cc2 = marker_alleles[2 * umm];
	  }
	  if (fprintf(outfile, "%u\t%s\t%g\t%u\t%c\t%c\n", ukk, &(marker_ids[map_reverse[ujj] * max_marker_id_len]), gd_vals[ujj], pos_buf[ujj], cc, cc2) < 0) {
	    goto merge_datasets_ret_WRITE_FAIL;
	  }
	}
      }
    } else {
      for (ulii = 0; ulii < chrom_ct; ulii++) {
	uii = chrom_start[ulii + 1];
	ujj = chrom_start[ulii];
	ukk = chrom_id[ulii];
	for (; ujj < uii; ujj++) {
	  umm = map_reverse[ujj];
	  if (keep_allele_order || (!is_set(reversed, ujj))) {
	    bufptr = &(marker_alleles[(2 * umm) * max_marker_allele_len]);
	    bufptr2 = &(marker_alleles[(2 * umm + 1) * max_marker_allele_len]);
	  } else {
	    bufptr = &(marker_alleles[(2 * umm + 1) * max_marker_allele_len]);
	    bufptr2 = &(marker_alleles[(2 * umm) * max_marker_allele_len]);
	  }
	  if (fprintf(outfile, "%u\t%s\t%g\t%u\t%s\t%s\n", ukk, &(marker_ids[map_reverse[ujj] * max_marker_id_len]), gd_vals[ujj], pos_buf[ujj], bufptr, bufptr2) < 0) {
	    goto merge_datasets_ret_WRITE_FAIL;
	  }
	}
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
    sprintf(logbuf, "%" PRIu64 " overlapping markers, %" PRIu64 " genotyped in both filesets.\n%" PRIu64 " concordant, for a concordance rate of %g.\n", diff_total_overlap, diff_not_both_genotyped, diff_not_both_genotyped - diff_discordant, 1.0 - (((double)diff_discordant) / ((double)diff_not_both_genotyped)));
    logprintb();
  }

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
  fclose_cond(mergelistfile);
  fclose_cond(outfile);
  wkspace_reset(wkspace_mark);
  return retval;
}
