#include <fcntl.h>
#include <sys/mman.h>
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

int sort_item_ids_nx(char** sorted_ids_ptr, int** id_map_ptr, int item_ct, char* item_ids, unsigned long max_id_len) {
  // Version of sort_item_ids() with no exclusion.
  int ii;
  int jj;
  char* sorted_ids;
  sorted_ids = (char*)wkspace_alloc(item_ct * max_id_len);
  if (!sorted_ids) {
    return RET_NOMEM;
  }
  *sorted_ids_ptr = sorted_ids;
  *id_map_ptr = (int*)wkspace_alloc(item_ct * sizeof(int));
  if (!(*id_map_ptr)) {
    return RET_NOMEM;
  }
  for (ii = 0; ii < item_ct; ii++) {
    memcpy(&(sorted_ids[ii * max_id_len]), &(item_ids[ii * max_id_len]), max_id_len);
    (*id_map_ptr)[ii] = ii;
  }
  if (qsort_ext(sorted_ids, item_ct, max_id_len, strcmp_deref, (char*)(*id_map_ptr), sizeof(int))) {
    return RET_NOMEM;
  }
  jj = item_ct - 1;
  for (ii = 0; ii < jj; ii++) {
    if (!strcmp(&(sorted_ids[ii * max_id_len]), &(sorted_ids[(ii + 1) * max_id_len]))) {
      logprint("Error: Duplicate IDs.\n");
      return RET_INVALID_FORMAT;
    }
  }
  return 0;
}

int indiv_major_to_snp_major(char* indiv_major_fname, char* outname, FILE** outfile_ptr, int unfiltered_marker_ct) {
  // This implementation only handles large files on 64-bit Unix systems; a
  // more portable version needs to be written for Windows and 32-bit Unix.
  int in_fd = open(indiv_major_fname, O_RDONLY);
  unsigned char* in_contents = (unsigned char*)MAP_FAILED;
  int unfiltered_marker_ct4 = (unfiltered_marker_ct + 3) / 4;
  struct stat sb;
  unsigned char* icoff;
  int retval;
  int indiv_ct;
  int indiv_ct4;
  int indiv_ct4l;
  long max_4blocks_in_mem;
  int superblock_offset;
  int block_last_marker;
  int ii;
  int add_val;
  int rshift_val;
  int jj;
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
    ii = 3;
  } else if ((!in_contents[0]) && (!((sb.st_size - 1) % unfiltered_marker_ct4))) {
    ii = 1;
  } else {
    ii = 0;
  }
  icoff = &(in_contents[ii]);
  if ((sb.st_size - ii) % unfiltered_marker_ct4) {
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
    for (ii = 0; ii < block_last_marker; ii++) {
      rshift_val = (ii % 4) * 2;
      add_val = ii / 4;
      for (jj = 0; jj < indiv_ct4l; jj++) {
        *write_ptr++ = ((icoff[4 * jj * unfiltered_marker_ct4 + add_val] >> rshift_val) & 3) + (((icoff[(4 * jj + 1) * unfiltered_marker_ct4 + add_val] >> rshift_val) & 3) << 2) + (((icoff[(4 * jj + 2) * unfiltered_marker_ct4 + add_val] >> rshift_val) & 3) << 4) + (((icoff[(4 * jj + 3) * unfiltered_marker_ct4 + add_val] >> rshift_val) & 3) << 6);
      }
      if (indiv_ct % 4) {
	*write_ptr = 0;
	for (jj = 0; jj < (indiv_ct % 4); jj++) {
	  *write_ptr |= ((icoff[(jj + indiv_ct) * unfiltered_marker_ct4 + add_val] >> rshift_val) & 3) << (jj * 2);
	}
	write_ptr++;
      }
    }
    if (fwrite_checked(wkspace_base, ((long long)block_last_marker) * indiv_ct4, *outfile_ptr)) {
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

const char errstr_map_format[] = "Error: Improperly formatted .map file.\n";

int load_map_or_bim(FILE** mapfile_ptr, char* mapname, int binary_files, int* map_cols_ptr, unsigned int* unfiltered_marker_ct_ptr, unsigned int* marker_exclude_ct_ptr, unsigned long* max_marker_id_len_ptr, unsigned int* plink_maxsnp_ptr, unsigned long** marker_exclude_ptr, double** set_allele_freqs_ptr, char** marker_alleles_ptr, char** marker_ids_ptr, Chrom_info* chrom_info_ptr, unsigned int** marker_pos_ptr, char* extractname, char* excludename, char* freqname, int calculation_type, int recode_modifier, int* map_is_unsorted_ptr) {
  int marker_ids_needed = (extractname || excludename || freqname || (calculation_type & (CALC_FREQ | CALC_WRITE_SNPLIST | CALC_LD_PRUNE | CALC_REGRESS_PCS)) || ((calculation_type & CALC_RECODE) && (recode_modifier & RECODE_LGEN)));
  unsigned int unfiltered_marker_ct = 0;
  unsigned long max_marker_id_len = 0;
  unsigned int plink_maxsnp = 4;
  unsigned long long loaded_chrom_mask = 0;
  int last_chrom = -1;
  int last_pos = 0;
  int marker_pos_needed = calculation_type & (CALC_GENOME | CALC_LD_PRUNE | CALC_REGRESS_PCS);
  unsigned int species = chrom_info_ptr->species;
  char* bufptr;
  char* bufptr2;
  unsigned long ulii;
  int ii;
  int jj;
  int cur_pos;
  int chroms_encountered_m1 = -1;
  unsigned int marker_uidx;
  if (fopen_checked(mapfile_ptr, mapname, "r")) {
    return RET_OPEN_FAIL;
  }
  // first pass: count columns, determine raw marker count, determine maximum
  // marker ID length if necessary.
  tbuf[MAXLINELEN - 6] = ' ';
  while (fgets(tbuf, MAXLINELEN - 5, *mapfile_ptr) != NULL) {
    if (!tbuf[MAXLINELEN - 6]) {
      sprintf(logbuf, "Error: Excessively long line in .map/.bim file (max %d chars).\n", MAXLINELEN - 8);
      logprintb();
      return RET_INVALID_FORMAT;
    }
    if (tbuf[0] > ' ') {
      bufptr = next_item(tbuf);
      if (no_more_items(bufptr)) {
	logprint(errstr_map_format);
	return RET_INVALID_FORMAT;
      }
      ulii = strlen_se(bufptr) + 1;
      if (ulii > max_marker_id_len) {
	max_marker_id_len = ulii;
      }
      if (marker_ids_needed || (!unfiltered_marker_ct)) {
	if (ulii > (plink_maxsnp + 1)) {
	  plink_maxsnp = ulii + 1;
	}
	if (!unfiltered_marker_ct) {
	  bufptr = next_item_mult(bufptr, 2);
	  if (binary_files) {
	    bufptr = next_item_mult(bufptr, 2);
	  }
	  if (!bufptr) {
	    logprint(errstr_map_format);
	    return RET_INVALID_FORMAT;
	  }
	  if (*bufptr > ' ') {
	    *map_cols_ptr = 4;
	  }
	}
      }
      unfiltered_marker_ct++;
    }
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
  *plink_maxsnp_ptr = plink_maxsnp;
  rewind(*mapfile_ptr);
  ii = (unfiltered_marker_ct + (BITCT - 1)) / BITCT;

  // unfiltered_marker_ct can be very large, so use wkspace for all allocations
  // that are a multiple of it

  // permanent stack allocation #1: marker_exclude
  if (wkspace_alloc_ul_checked(marker_exclude_ptr, ii * sizeof(long))) {
    return RET_NOMEM;
  }
  fill_ulong_zero(*marker_exclude_ptr, ii);
  // permanent stack allocation #2: set_allele_freqs
  *set_allele_freqs_ptr = (double*)wkspace_alloc(unfiltered_marker_ct * sizeof(double));
  if (!(*set_allele_freqs_ptr)) {
    return RET_NOMEM;
  }
  for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
    (*set_allele_freqs_ptr)[marker_uidx] = -1.0;
  }
  fill_uint_zero(chrom_info_ptr->chrom_file_order, MAX_POSSIBLE_CHROM);
  fill_uint_zero(chrom_info_ptr->chrom_file_order_marker_idx, MAX_POSSIBLE_CHROM + 1);
  fill_uint_zero(chrom_info_ptr->chrom_start, MAX_POSSIBLE_CHROM);
  fill_uint_zero(chrom_info_ptr->chrom_end, MAX_POSSIBLE_CHROM);
  // permanent stack allocation #3, if needed: marker_pos
  if (marker_pos_needed) {
    *marker_pos_ptr = (unsigned int*)wkspace_alloc(unfiltered_marker_ct * sizeof(int));
    if (!(*marker_pos_ptr)) {
      return RET_NOMEM;
    }
  }
  if (binary_files) {
    if (freqname || (calculation_type & (CALC_FREQ | CALC_RECODE | CALC_REGRESS_PCS))) {
      *marker_alleles_ptr = (char*)wkspace_alloc(unfiltered_marker_ct * 2 * sizeof(char));
      if (!(*marker_alleles_ptr)) {
	return RET_NOMEM;
      }
      memset(*marker_alleles_ptr, 0, unfiltered_marker_ct * 2);
    }
  }
  if (marker_ids_needed) {
    *marker_ids_ptr = (char*)wkspace_alloc(unfiltered_marker_ct * max_marker_id_len);
    if (!(*marker_ids_ptr)) {
      return RET_NOMEM;
    }
  }

  // second pass: actually load stuff
  for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
    do {
      if (fgets(tbuf, MAXLINELEN, *mapfile_ptr) == NULL) {
        return RET_READ_FAIL;
      }
    } while (tbuf[0] <= ' ');
    jj = marker_code(species, tbuf);
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
	if (calculation_type != CALC_MAKE_BED) {
	  logprint("Error: .map/.bim file is unsorted.  Use --make-bed by itself to remedy this.\n");
	  return RET_INVALID_FORMAT;
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
      bufptr = next_item(tbuf);
      if (*marker_ids_ptr) {
        if (no_more_items(bufptr)) {
	  logprint(errstr_map_format);
          return RET_INVALID_FORMAT;
        }
        read_next_terminate(&((*marker_ids_ptr)[marker_uidx * max_marker_id_len]), bufptr);
      }
      bufptr = next_item(bufptr);
      if (*map_cols_ptr == 4) {
	bufptr = next_item(bufptr);
      }
      if (no_more_items(bufptr)) {
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
	if ((calculation_type & marker_pos_needed) && jj) {
	  (*marker_pos_ptr)[marker_uidx] = cur_pos;
	}
        if (binary_files && (*marker_alleles_ptr)) {
	  bufptr = next_item(bufptr);
	  bufptr2 = next_item(bufptr);
	  if (no_more_items(bufptr2)) {
	    logprint(errstr_map_format);
	    return RET_INVALID_FORMAT;
	  }
	  (*marker_alleles_ptr)[marker_uidx * 2] = *bufptr;
	  (*marker_alleles_ptr)[marker_uidx * 2 + 1] = *bufptr2;
	}
      }
    }
  }
  chrom_info_ptr->chrom_end[last_chrom] = marker_uidx;
  chrom_info_ptr->chrom_ct = ++chroms_encountered_m1;
  chrom_info_ptr->chrom_file_order_marker_idx[chroms_encountered_m1] = marker_uidx;
  return 0;
}

void sort_marker_chrom_pos(long long* ll_buf, unsigned int marker_ct, int* pos_buf, unsigned int* chrom_start, unsigned int* chrom_id, unsigned int* chrom_ct_ptr) {
  // assumes ll_buf is filled with chromosome idxs in high 32 bits, and
  // internal marker indices in low 32 bits.  pos_buf is expected to have
  // base-pair positions.
  unsigned int marker_idx;
  unsigned int uii;
  unsigned int cur_chrom;
  unsigned int chrom_ct;
#ifdef __cplusplus
  std::sort(ll_buf, &(ll_buf[marker_ct]));
#else
  qsort(ll_buf, marker_ct, sizeof(long long), llcmp);
#endif
  cur_chrom = ll_buf[0] >> 32;
  chrom_ct = 0;
  chrom_start[0] = 0;
  chrom_id[0] = cur_chrom;
  uii = (unsigned int)ll_buf[0];
  ll_buf[0] = ((long long)uii) | (((long long)pos_buf[uii]) << 32);
  for (marker_idx = 1; marker_idx < marker_ct; marker_idx++) {
    if ((ll_buf[marker_idx] >> 32) != cur_chrom) {
      cur_chrom = ll_buf[marker_idx] >> 32;
      chrom_start[++chrom_ct] = marker_idx;
      chrom_id[chrom_ct] = cur_chrom;
    }
    uii = (unsigned int)ll_buf[marker_idx];
    ll_buf[marker_idx] = ((long long)uii) | (((long long)pos_buf[uii]) << 32);
  }
  chrom_start[++chrom_ct] = marker_ct;
  for (uii = 0; uii < chrom_ct; uii++) {
#ifdef __cplusplus
    std::sort(&(ll_buf[chrom_start[uii]]), &(ll_buf[chrom_start[uii + 1]]));
#else
    qsort(&(ll_buf[chrom_start[uii]]), chrom_start[uii + 1] - chrom_start[uii], sizeof(long long), llcmp);
#endif
  }
  *chrom_ct_ptr = chrom_ct;
}

int sort_and_write_bim(int** map_reverse_ptr, FILE* mapfile, int map_cols, FILE** bimfile_ptr, char* outname, char* outname_end, unsigned int unfiltered_marker_ct, unsigned long* marker_exclude, unsigned int marker_ct, unsigned long max_marker_id_len, int species, char* marker_alleles) {
  int tmp_map = (bimfile_ptr == NULL);
  int load_markers = (marker_alleles == NULL) && (!tmp_map);
  FILE* map_outfile;
  long long* ll_buf;
  char* marker_ids;
  double* gd_vals;
  int* pos_buf;
  unsigned int marker_uidx;
  unsigned int marker_idx;
  char* bufptr;
  unsigned int uii;
  unsigned int ujj;
  unsigned int chrom_start[MAX_POSSIBLE_CHROM + 1];
  unsigned int chrom_id[MAX_POSSIBLE_CHROM];
  unsigned int cur_chrom;
  unsigned int chrom_ct;
  char cc;
  if (wkspace_alloc_i_checked(map_reverse_ptr, unfiltered_marker_ct * sizeof(int))) {
    return RET_NOMEM;
  }
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
  if (wkspace_alloc_ll_checked(&ll_buf, marker_ct * sizeof(long long)) ||
      wkspace_alloc_c_checked(&marker_ids, marker_ct * max_marker_id_len) ||
      wkspace_alloc_d_checked(&gd_vals, marker_ct * sizeof(double)) ||
      wkspace_alloc_i_checked(&pos_buf, marker_ct * sizeof(int))) {
    return RET_NOMEM;
  }
  if (load_markers) {
    if (wkspace_alloc_c_checked(&marker_alleles, unfiltered_marker_ct * 2)) {
      return RET_NOMEM;
    }
  }
  rewind(mapfile);
  marker_idx = 0;
  for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
    do {
      if (fgets(tbuf, MAXLINELEN, mapfile) == NULL) {
	return RET_READ_FAIL;
      }
    } while (tbuf[0] <= ' ');
    if (is_set(marker_exclude, marker_uidx)) {
      continue;
    }
    ll_buf[marker_idx] = (((long long)marker_code(species, tbuf)) << 32) + marker_uidx;
    bufptr = next_item(tbuf);
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
    pos_buf[marker_idx++] = atoi(bufptr);
    if (load_markers) {
      bufptr = next_item(bufptr);
      marker_alleles[2 * marker_uidx] = *bufptr;
      marker_alleles[2 * marker_uidx + 1] = *(next_item(bufptr));
    }
  }
  sort_marker_chrom_pos(ll_buf, marker_ct, pos_buf, chrom_start, chrom_id, &chrom_ct);

  if (tmp_map) {
    strcpy(outname_end, ".map.tmp");
    bimfile_ptr = &map_outfile;
  } else {
    strcpy(outname_end, ".bim");
  }
  if (fopen_checked(bimfile_ptr, outname, "w")) {
    return RET_OPEN_FAIL;
  }

  marker_idx = 0;
  for (uii = 0; uii < chrom_ct; uii++) {
    cur_chrom = chrom_id[uii];
    ujj = chrom_start[uii + 1];
    for (; marker_idx < ujj; marker_idx++) {
      marker_uidx = (unsigned int)ll_buf[marker_idx];
      if (tmp_map) {
        if (fprintf(*bimfile_ptr, "%u\t%s\t%g\t%u\n", cur_chrom, &(marker_ids[marker_uidx * max_marker_id_len]), gd_vals[marker_uidx], (unsigned int)(ll_buf[marker_idx] >> 32)) < 0) {
	  return RET_WRITE_FAIL;
        }
      } else {
        cc = marker_alleles[2 * marker_uidx];
        if (fprintf(*bimfile_ptr, "%u\t%s\t%g\t%u\t%c\t%c\n", cur_chrom, &(marker_ids[marker_uidx * max_marker_id_len]), gd_vals[marker_uidx], (unsigned int)(ll_buf[marker_idx] >> 32), cc? cc : '0', marker_alleles[2 * marker_uidx + 1]) < 0) {
	  return RET_WRITE_FAIL;
        }
      }
      (*map_reverse_ptr)[marker_uidx] = marker_idx;
    }
  }
  if (tmp_map) {
    if (fclose_null(bimfile_ptr)) {
      return RET_WRITE_FAIL;
    }
  }

  wkspace_reset((unsigned char*)ll_buf);
  return 0;
}

int text_to_bed(FILE** bedtmpfile_ptr, FILE** famtmpfile_ptr, FILE** bimtmpfile_ptr, FILE** outfile_ptr, char* outname, char* outname_end, FILE** pedfile_ptr, unsigned long long* line_locs, unsigned long long* line_mids, int pedbuflen, FILE** famfile_ptr, FILE* mapfile, int map_cols, unsigned int unfiltered_marker_ct, unsigned long* marker_exclude, unsigned int marker_ct, unsigned int unfiltered_indiv_ct, unsigned long* indiv_exclude, unsigned int indiv_ct, char* marker_alleles, char* person_ids, unsigned int max_person_id_len, char* paternal_ids, unsigned int max_paternal_id_len, char* maternal_ids, unsigned int max_maternal_id_len, unsigned char* sex_info, char* pheno_c, double* pheno_d, double missing_phenod, char* output_missing_pheno, unsigned long max_marker_id_len, int map_is_unsorted, unsigned int species) {
  unsigned int marker_ct4 = (marker_ct + 3) / 4;
  unsigned int indiv_ct4 = (indiv_ct + 3) / 4;
  unsigned char* wkspace_mark = wkspace_base;
  int affection = (pheno_c != NULL);
  int phenos_present = (affection || (pheno_d != NULL));
  int* map_reverse = NULL;
  unsigned int marker_uidx;
  unsigned int marker_idx;
  unsigned int indiv_uidx;
  unsigned int indiv_idx;
  char* loadbuf;
  char* bufptr;
  char* bufptr2;
  unsigned char* writebuf;
  char cc;
  char cc2;
  unsigned char newval;
  unsigned int uii;
  int retval;
  if (wkspace_alloc_c_checked(&loadbuf, pedbuflen)) {
    return RET_NOMEM;
  }
  if (map_is_unsorted) {
    retval = sort_and_write_bim(&map_reverse, mapfile, map_cols, bimtmpfile_ptr, outname, outname_end, unfiltered_marker_ct, marker_exclude, marker_ct, max_marker_id_len, species, marker_alleles);
    if (retval) {
      return retval;
    }
  }
  writebuf = wkspace_base;
  if (((long long)wkspace_left) >= (marker_ct * ((long long)indiv_ct4))) {
    strcpy(outname_end, ".bed");
    sprintf(logbuf, "Converting .ped file to SNP-major %s.\n", outname);
    logprintb();
    if (fopen_checked(bedtmpfile_ptr, outname, "wb")) {
      return RET_OPEN_FAIL;
    }
    strcpy(outname_end, ".fam");
    if (fopen_checked(famtmpfile_ptr, outname, "w")) {
      return RET_OPEN_FAIL;
    }
    if (fwrite_checked("l\x1b\x01", 3, *bedtmpfile_ptr)) {
      return RET_WRITE_FAIL;
    }
    memset(writebuf, 0, marker_ct * ((long long)indiv_ct4));
    indiv_idx = 0;
    rewind(*pedfile_ptr);
    for (indiv_uidx = 0; indiv_idx < indiv_ct; indiv_uidx++) {
      if (is_set(indiv_exclude, indiv_uidx)) {
	indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx + 1);
	if (fseeko(*pedfile_ptr, line_locs[indiv_uidx], SEEK_SET)) {
	  return RET_READ_FAIL;
	}
      }
      if (fgets(loadbuf, pedbuflen, *pedfile_ptr) == NULL) {
	return RET_READ_FAIL;
      }
      bufptr = tbuf;
      bufptr2 = &(person_ids[indiv_uidx * max_person_id_len]);
      bufptr += sprintf(tbuf, "%s %s %s %d ", bufptr2, paternal_ids? (&(paternal_ids[indiv_uidx * max_paternal_id_len])) : "0", maternal_ids? (&(maternal_ids[indiv_uidx * max_maternal_id_len])) : "0", sex_info? sex_info[indiv_uidx] : 0);
      tbuf[strlen_se(bufptr2)] = ' ';
      if (affection) {
	cc = pheno_c[indiv_uidx];
	if (cc == -1) {
	  sprintf(bufptr, "%s\n", output_missing_pheno);
	} else {
	  sprintf(bufptr, "%d\n", cc + 1);
	}
      } else if ((!phenos_present) || (pheno_d[indiv_uidx] == missing_phenod)) {
	sprintf(bufptr, "%s\n", output_missing_pheno);
      } else {
	sprintf(bufptr, "%g\n", pheno_d[indiv_uidx]);
      }
      if (fputs(tbuf, *famtmpfile_ptr) == EOF) {
	return RET_WRITE_FAIL;
      }
      bufptr = &(loadbuf[line_mids[indiv_uidx] - line_locs[indiv_uidx]]);
      marker_idx = 0;
      for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
	if (is_set(marker_exclude, marker_uidx)) {
	  bufptr++;
	  while ((*bufptr == ' ') || (*bufptr == '\t')) {
	    bufptr++;
	  }
	  bufptr++;
	  while ((*bufptr == ' ') || (*bufptr == '\t')) {
	    bufptr++;
	  }
	  continue;
	}
	cc = *bufptr++;
	while ((*bufptr == ' ') || (*bufptr == '\t')) {
	  bufptr++;
	}
	cc2 = *bufptr++;
	while ((*bufptr == ' ') || (*bufptr == '\t')) {
	  bufptr++;
	}
	newval = 1;
	if (cc == marker_alleles[2 * marker_uidx]) {
	  if (cc2 == cc) {
	    newval = 0;
	  } else if (cc2 == marker_alleles[2 * marker_uidx + 1]) {
	    newval = 2;
	  }
	} else if (cc == marker_alleles[2 * marker_uidx + 1]) {
	  if (cc2 == marker_alleles[2 * marker_uidx]) {
	    newval = 2;
	  } else if (cc2 == cc) {
	    newval = 3;
	  }
	}
	if (map_is_unsorted) {
	  writebuf[map_reverse[marker_uidx] * indiv_ct4 + (indiv_idx / 4)] |= newval << ((indiv_idx % 4) * 2);
	} else {
	  writebuf[marker_idx * indiv_ct4 + (indiv_idx / 4)] |= newval << ((indiv_idx % 4) * 2);
	}
	marker_idx++;
      }
      indiv_idx++;
    }
    if (fwrite_checked(writebuf, marker_ct * ((long long)indiv_ct4), *bedtmpfile_ptr)) {
      return RET_WRITE_FAIL;
    }
    if (fclose_null(bedtmpfile_ptr)) {
      return RET_WRITE_FAIL;
    }
  } else if (wkspace_left < marker_ct4) {
    return RET_NOMEM;
  } else {
    // this should be rewritten
    logprint("Converting very large .ped file to temporary individual-major .bed.\n");
    strcpy(outname_end, ".fam");
    if (fopen_checked(famtmpfile_ptr, outname, "w")) {
      return RET_OPEN_FAIL;
    }
    // .bed.tmp must be last since filename is copied before
    // indiv_major_to_snp_major() call
    strcpy(outname_end, ".bed.tmp");
    if (fopen_checked(bedtmpfile_ptr, outname, "wb")) {
      return RET_OPEN_FAIL;
    }
    if (fwrite_checked("l\x1b", 3, *bedtmpfile_ptr)) {
      return RET_WRITE_FAIL;
    }
    rewind(*pedfile_ptr);

    indiv_idx = 0;
    for (indiv_uidx = 0; indiv_idx < indiv_ct; indiv_uidx++) {
      if (is_set(indiv_exclude, indiv_uidx)) {
	indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx + 1);
	if (fseeko(*pedfile_ptr, line_locs[indiv_uidx], SEEK_SET)) {
	  return RET_READ_FAIL;
	}
      }
      if (fgets(loadbuf, pedbuflen, *pedfile_ptr) == NULL) {
	return RET_READ_FAIL;
      }
      bufptr = &(loadbuf[line_mids[indiv_uidx] - line_locs[indiv_uidx]]);
      bufptr2 = bufptr - 2;
      while ((*bufptr2 == ' ') || (*bufptr2 == '\t')) {
	bufptr2--;
      }
      cc = 0;
      bufptr2[1] = '\n';
      // note that if the .ped file is missing a few columns, the generated
      // .fam will also be missing the same columns
      if (fwrite_checked(loadbuf, 2 + (long)(bufptr2 - loadbuf), *famtmpfile_ptr)) {
	return RET_WRITE_FAIL;
      }
      marker_idx = 0;
      memset(writebuf, 0, marker_ct4);
      for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
        if (is_set(marker_exclude, marker_uidx)) {
	  bufptr++;
	  while ((*bufptr == ' ') || (*bufptr == '\t')) {
	    bufptr++;
	  }
	  bufptr++;
	  while ((*bufptr == ' ') || (*bufptr == '\t')) {
	    bufptr++;
	  }
	  continue;
	}
	cc = *bufptr++;
	while ((*bufptr == ' ') || (*bufptr == '\t')) {
	  bufptr++;
	}
	cc2 = *bufptr++;
	while ((*bufptr == ' ') || (*bufptr == '\t')) {
	  bufptr++;
	}
	newval = 1;
        if (cc == marker_alleles[2 * marker_uidx]) {
	  if (cc2 == cc) {
	    newval = 0;
	  } else if (cc2 == marker_alleles[2 * marker_uidx + 1]) {
	    newval = 2;
	  }
	} else if (cc == marker_alleles[2 * marker_uidx + 1]) {
	  if (cc2 == marker_alleles[2 * marker_uidx]) {
	    newval = 2;
	  } else if (cc2 == cc) {
	    newval = 3;
	  }
	}
	if (map_is_unsorted) {
	  uii = map_reverse[marker_uidx];
	  writebuf[uii / 4] |= newval << ((uii % 4) * 2);
	} else {
          writebuf[marker_idx / 4] |= newval << ((marker_idx % 4) * 2);
	}
        marker_idx++;
      }
      if (fwrite_checked(writebuf, marker_ct4, *bedtmpfile_ptr)) {
	return RET_WRITE_FAIL;
      }
      indiv_idx++;
    }
    if (fclose_null(bedtmpfile_ptr)) {
      return RET_WRITE_FAIL;
    }

    sprintf(logbuf, "Now transposing to SNP-major %s.\n", outname);
    logprintb();
    strcpy(tbuf, outname);
    strcpy(outname_end, ".bed");
    retval = indiv_major_to_snp_major(tbuf, outname, outfile_ptr, marker_ct);
    unlink(tbuf);
    if (retval) {
      return retval;
    }
  }

  // ----- .map -> .bim -----
  if (!map_is_unsorted) {
    strcpy(outname_end, ".bim");
    if (fopen_checked(bimtmpfile_ptr, outname, "w")) {
      return RET_OPEN_FAIL;
    }

    rewind(mapfile);
    for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
      do {
	if (fgets(tbuf, MAXLINELEN, mapfile) == NULL) {
	  return RET_READ_FAIL;
	}
      } while (tbuf[0] <= ' ');
      if (is_set(marker_exclude, marker_uidx)) {
	continue;
      }
      bufptr = next_item(tbuf);
      bufptr = next_item(bufptr);
      if (map_cols == 4) {
	bufptr = next_item(bufptr);
      }
      while (!is_space_or_eoln(*bufptr)) {
	bufptr++;
      }
      *bufptr++ = '\t';
      if (marker_alleles[2 * marker_uidx]) {
	*bufptr++ = marker_alleles[2 * marker_uidx];
      } else {
	*bufptr++ = '0';
      }
      *bufptr++ = '\t';
      *bufptr++ = marker_alleles[2 * marker_uidx + 1];
      *bufptr++ = '\n';
      if (fwrite_checked(tbuf, (int)(bufptr - tbuf), *bimtmpfile_ptr)) {
	return RET_WRITE_FAIL;
      }
    }
  }
  fclose(*pedfile_ptr);
  fclose_null(bimtmpfile_ptr);
  fclose_null(famtmpfile_ptr);
  strcpy(outname_end, ".bed");
  if (fopen_checked(pedfile_ptr, outname, "rb")) {
    return RET_OPEN_FAIL;
  }
  strcpy(outname_end, ".fam");
  if (fopen_checked(famfile_ptr, outname, "r")) {
    return RET_OPEN_FAIL;
  }
  logprint("Automatic --make-bed complete.\n");
  wkspace_reset(wkspace_mark);
  return 0;
}

void fill_bmap_raw(unsigned char* bmap_raw, unsigned int ct_mod4) {
  // possibilities:
  // 0. A1 -> A1, A2 -> A2 [0..255]
  // 1. A1 -> A2, A2 -> A1 [256..511]
  // 1 last char. [512..575]
  unsigned char imap[4] = {3, 1, 2, 0};
  unsigned char* bmap = bmap_raw;
  unsigned int uii = 0;
  unsigned int ujj;
  unsigned int ukk;
  unsigned int umm;
  unsigned int unn;
  unsigned int extra_ct;
  do {
    *bmap++ = uii++;
  } while (uii < 256);
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

int make_bed(FILE* bedfile, int bed_offset, FILE* bimfile, int map_cols, FILE** bedoutfile_ptr, FILE** famoutfile_ptr, FILE** bimoutfile_ptr, char* outname, char* outname_end, unsigned int unfiltered_marker_ct, unsigned long* marker_exclude, unsigned int marker_ct, unsigned long* marker_reverse, unsigned int unfiltered_indiv_ct, unsigned long* indiv_exclude, unsigned int indiv_ct, char* person_ids, unsigned int max_person_id_len, char* paternal_ids, unsigned int max_paternal_id_len, char* maternal_ids, unsigned int max_maternal_id_len, unsigned char* sex_info, char* pheno_c, double* pheno_d, double missing_phenod, char* output_missing_pheno, unsigned long max_marker_id_len, int map_is_unsorted, unsigned int* indiv_sort_map, int species) {
  unsigned int unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  unsigned long indiv_ct4 = (indiv_ct + 3) / 4;
  unsigned char* wkspace_mark = wkspace_base;
  int affection = (pheno_c != NULL);
  int phenos_present = (affection || (pheno_d != NULL));
  unsigned int marker_uidx;
  unsigned int marker_idx;
  unsigned int indiv_uidx;
  unsigned int indiv_uidx2;
  unsigned int indiv_idx;
  unsigned int ii_mod4;
  unsigned char* loadbuf;
  unsigned char* writebuf;
  char* bufptr;
  char* cptr;
  unsigned char cc;
  char cc2;
  unsigned int pct;
  unsigned int loop_end;
  int* map_reverse;
  unsigned char* writeptr;
  unsigned char bmap_raw[576];
  unsigned char* bmap;
  unsigned char* bmap2;
  int retval;
  if (wkspace_alloc_uc_checked(&loadbuf, unfiltered_indiv_ct4)) {
    return RET_NOMEM;
  }
  fill_bmap_raw(bmap_raw, indiv_ct & 3);

  if (fseeko(bedfile, bed_offset, SEEK_SET)) {
    return RET_READ_FAIL;
  }
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
    retval = sort_and_write_bim(&map_reverse, bimfile, map_cols, bimoutfile_ptr, outname, outname_end, unfiltered_marker_ct, marker_exclude, marker_ct, max_marker_id_len, species, NULL);
    if (retval) {
      return retval;
    }

    if (wkspace_alloc_uc_checked(&writebuf, marker_ct * indiv_ct4)) {
      logprint("\nError: Insufficient memory for current --make-bed implementation.  Try raising\nthe --memory value for now.\n");
      return RET_CALC_NOT_YET_SUPPORTED;
    } else {
      for (pct = 1; pct <= 100; pct++) {
	loop_end = ((unsigned long long)pct * marker_ct) / 100;
	for (; marker_idx < loop_end; marker_idx++) {
	  if (is_set(marker_exclude, marker_uidx)) {
	    marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx + 1);
	    if (fseeko(bedfile, bed_offset + (unsigned long long)marker_uidx * unfiltered_indiv_ct4, SEEK_SET)) {
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
	  writeptr = &(writebuf[indiv_ct4 * map_reverse[marker_idx]]);
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
	    printf("\b\b\b%u%%", pct);
	  } else {
	    printf("\b\b%u%%", pct);
	  }
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
    for (pct = 1; pct <= 100; pct++) {
      loop_end = ((unsigned long long)pct * marker_ct) / 100;
      for (; marker_idx < loop_end; marker_idx++) {
	if (is_set(marker_exclude, marker_uidx)) {
	  marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx + 1);
	  if (fseeko(bedfile, bed_offset + (unsigned long long)marker_uidx * unfiltered_indiv_ct4, SEEK_SET)) {
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
	  printf("\b\b\b%u%%", pct);
	} else {
	  printf("\b\b%u%%", pct);
	}
	fflush(stdout);
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
    bufptr += sprintf(tbuf, "%s %s %s %d ", cptr, paternal_ids? (&(paternal_ids[indiv_uidx * max_paternal_id_len])) : "0", maternal_ids? (&(maternal_ids[indiv_uidx * max_maternal_id_len])) : "0", sex_info? sex_info[indiv_uidx] : 0);
    tbuf[strlen_se(cptr)] = ' ';
    if (affection) {
      cc2 = pheno_c[indiv_uidx];
      if (cc2 == -1) {
        sprintf(bufptr, "%s\n", output_missing_pheno);
      } else {
        sprintf(bufptr, "%d\n", cc2 + 1);
      }
    } else if ((!phenos_present) || (pheno_d[indiv_uidx] == missing_phenod)) {
      sprintf(bufptr, "%s\n", output_missing_pheno);
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
      do {
	if (fgets(tbuf, MAXLINELEN, bimfile) == NULL) {
	  return RET_READ_FAIL;
	}
      } while (tbuf[0] <= ' ');
      if (is_set(marker_exclude, marker_uidx)) {
	continue;
      }
      if (fputs(tbuf, *bimoutfile_ptr) == EOF) {
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

const char errstr_fam_format[] = "Error: Improperly formatted .fam/.ped file.\n";

int load_fam(FILE* famfile, unsigned long buflen, int fam_col_1, int fam_col_34, int fam_col_5, int fam_col_6, int true_fam_col_6, int missing_pheno, int missing_pheno_len, int affection_01, unsigned int* unfiltered_indiv_ct_ptr, char** person_ids_ptr, unsigned int* max_person_id_len_ptr, char** paternal_ids_ptr, unsigned int* max_paternal_id_len_ptr, char** maternal_ids_ptr, unsigned int* max_maternal_id_len_ptr, unsigned char** sex_info_ptr, int* affection_ptr, char** pheno_c_ptr, double** pheno_d_ptr, unsigned long** founder_info_ptr, unsigned long** indiv_exclude_ptr, int binary_files, unsigned long long** line_locs_ptr, unsigned long long** line_mids_ptr, int* pedbuflen_ptr) {
  char* bufptr;
  unsigned int unfiltered_indiv_ct = 0;
  unsigned int max_person_id_len = 4;
  unsigned int max_paternal_id_len = 2;
  unsigned int max_maternal_id_len = 2;
  int affection = 1;
  unsigned long long last_tell = 0;
  unsigned long new_buflen = 0;
  unsigned char* wkspace_mark = wkspace_base;
  char* linebuf;
  char* person_ids;
  char* paternal_ids = NULL;
  char* maternal_ids = NULL;
  unsigned char* sex_info = NULL;
  char* pheno_c = NULL;
  double* pheno_d = NULL;
  char cc;
  unsigned long long* line_locs;
  unsigned long long* tmp_ullp;
  unsigned int max_people;
  unsigned long tmp_len;
  unsigned long tmp_len2;
  unsigned int indiv_uidx;
  unsigned int unfiltered_indiv_ctl;
  int ii;
  char* fgets_return;
  if (wkspace_alloc_c_checked(&linebuf, buflen)) {
    return RET_NOMEM;
  }
  linebuf[buflen - 1] = ' ';
  line_locs = (unsigned long long*)wkspace_base;
  max_people = wkspace_left / sizeof(long long);
  // ----- .fam/[.ped first columns] read, first pass -----
  // count number of people, determine maximum person/father/mother ID lengths,
  // affection status, verify all floating point phenotype values are valid
  while (fgets(linebuf, buflen, famfile) != NULL) {
    if (linebuf[0] > ' ') {
      if (linebuf[0] != '#') {
	if (fam_col_1) {
	  bufptr = next_item(linebuf);
	} else {
	  bufptr = linebuf;
	}
	tmp_len = strlen_se(linebuf) + strlen_se(bufptr) + 2;
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
	  if (no_more_items(bufptr)) {
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
    logprint("Error: Nobody in .fam/.ped file.\n");
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
  if (fam_col_5) {
    if (wkspace_alloc_uc_checked(sex_info_ptr, unfiltered_indiv_ct)) {
      return RET_NOMEM;
    }
    sex_info = *sex_info_ptr;
  }
  unfiltered_indiv_ctl = (unfiltered_indiv_ct + (BITCT - 1)) / BITCT;
  if (wkspace_alloc_ul_checked(founder_info_ptr, unfiltered_indiv_ctl * sizeof(long))) {
    return RET_NOMEM;
  }
  if (wkspace_alloc_ul_checked(indiv_exclude_ptr, unfiltered_indiv_ctl * sizeof(long))) {
    return RET_NOMEM;
  }

  if (fam_col_6) {
    if (affection) {
      pheno_c = (char*)malloc(unfiltered_indiv_ct * sizeof(char));
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
  if (wkspace_alloc_ull_checked(&tmp_ullp, unfiltered_indiv_ct * sizeof(long long))) {
    return RET_NOMEM;
  }
  for (ii = unfiltered_indiv_ct - 1; ii >= 0; ii--) {
    tmp_ullp[ii] = line_locs[ii];
  }
  line_locs = tmp_ullp;
  if (!binary_files) {
    *line_locs_ptr = (unsigned long long*)malloc(unfiltered_indiv_ct * sizeof(long long));
    if (!(*line_locs_ptr)) {
      return RET_NOMEM;
    }
    *line_mids_ptr = (unsigned long long*)malloc(unfiltered_indiv_ct * sizeof(long long));
    if (!(*line_mids_ptr)) {
      return RET_NOMEM;
    }
    *pedbuflen_ptr = buflen;
  }
  if (fam_col_34) {
    fill_ulong_zero(*founder_info_ptr, unfiltered_indiv_ctl);
  } else {
    fill_ulong_one(*founder_info_ptr, unfiltered_indiv_ctl);
  }
  if (fam_col_5) {
    memset(sex_info, 0, unfiltered_indiv_ct);
  }
  fill_ulong_zero(*indiv_exclude_ptr, unfiltered_indiv_ctl);

  // ----- .fam/[.ped first columns] read, second pass -----
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
	  sex_info[indiv_uidx] = 1;
	} else if (*bufptr == '2') {
	  sex_info[indiv_uidx] = 2;
	}
      }
    }
    if (fam_col_6) {
      bufptr = next_item(bufptr);
      if (affection) {
	if (is_missing(bufptr, missing_pheno, missing_pheno_len, affection_01)) {
	  pheno_c[indiv_uidx] = -1;
	} else if (affection_01) {
	  pheno_c[indiv_uidx] = *bufptr - '0';
	} else {
	  pheno_c[indiv_uidx] = *bufptr - '1';
	}
      } else {
	if (sscanf(bufptr, "%lg", &(pheno_d[indiv_uidx])) != 1) {
	  pheno_d[indiv_uidx] = (double)missing_pheno;
	}
      }
    }
    if (true_fam_col_6 && (!fam_col_6)) {
      bufptr = next_item(bufptr);
    }
    if (!binary_files) {
      bufptr = next_item(bufptr);
      if (no_more_items(bufptr)) {
	logprint(errstr_fam_format);
	return RET_INVALID_FORMAT;
      }
      (*line_locs_ptr)[indiv_uidx] = line_locs[indiv_uidx];
      (*line_mids_ptr)[indiv_uidx] = line_locs[indiv_uidx] + (unsigned long)(bufptr - linebuf);
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

int lgen_to_bed(char* lgen_namebuf, char* outname, int missing_pheno, int affection_01, Chrom_info* chrom_info_ptr) {
  FILE* infile = NULL;
  FILE* outfile = NULL;
  char* name_end = (char*)memchr(lgen_namebuf, 0, FNAMESIZE);
  char* outname_end = (char*)memchr(outname, 0, FNAMESIZE);
  int map_cols = 3;
  unsigned long* marker_exclude = NULL;
  unsigned int marker_exclude_ct = 0;
  unsigned long max_marker_id_len = 0;
  double* set_allele_freqs = NULL;
  unsigned int unfiltered_marker_ct = 0;
  unsigned int marker_ct = 0;
  char* marker_alleles = NULL;
  unsigned int* marker_allele_cts;
  char* marker_ids = NULL;
  unsigned int* marker_pos = NULL;
  unsigned int indiv_ct = 0;
  char* person_ids = NULL;
  char* paternal_ids = NULL;
  unsigned int max_paternal_id_len = 2;
  char* maternal_ids = NULL;
  unsigned int max_maternal_id_len = 2;
  unsigned char* sex_info = NULL;
  int affection = 0;
  unsigned long* founder_info = NULL;
  unsigned long* indiv_exclude = NULL;
  int map_is_unsorted = 0;
  char* pheno_c = NULL;
  double* pheno_d = NULL;
  char* sorted_marker_ids;
  int* marker_id_map;
  int* map_reverse;
  char* sorted_indiv_ids;
  int* indiv_id_map;
  unsigned char* writebuf;
  unsigned int indiv_ct4;
  unsigned int plink_maxsnp;
  unsigned int max_person_id_len;
  unsigned long marker_idx;
  unsigned char ucc;
  unsigned char* ucptr;
  unsigned char itable_short[4];
  unsigned char itable[256];
  unsigned char itable2[64];
  char* id_buf;
  char* cptr;
  char* cptr2;
  char* cptr3;
  char* cptr4;
  char* cptr5;
  char* cptr6;
  unsigned int indiv_idx;
  unsigned long ulii;
  unsigned int uii;
  unsigned int ujj;
  unsigned int ukk;
  unsigned int umm;
  long long lgen_size;
  unsigned int pct;
  long long lgen_next_thresh;
  char cc;
  char cc2;
  int ii;
  int jj;
  int retval;
  memcpy(name_end, ".map", 5);
  // use CALC_WRITE_SNPLIST to force loading of marker IDs, without other
  // unwanted side-effects
  retval = load_map_or_bim(&infile, lgen_namebuf, 0, &map_cols, &unfiltered_marker_ct, &marker_exclude_ct, &max_marker_id_len, &plink_maxsnp, &marker_exclude, &set_allele_freqs, NULL, &marker_ids, chrom_info_ptr, &marker_pos, '\0', '\0', '\0', CALC_WRITE_SNPLIST, 0, &map_is_unsorted);
  if (retval) {
    goto lgen_to_bed_ret_1;
  }
  marker_ct = unfiltered_marker_ct - marker_exclude_ct;
  retval = sort_item_ids_nx(&sorted_marker_ids, &marker_id_map, unfiltered_marker_ct, marker_ids, max_marker_id_len);
  if (retval) {
    goto lgen_to_bed_ret_1;
  }
  if (map_is_unsorted) {
    // Writes a temporary .map which is read later, and then deleted.
    retval = sort_and_write_bim(&map_reverse, infile, map_cols, NULL, outname, outname_end, unfiltered_marker_ct, marker_exclude, marker_ct, max_marker_id_len, chrom_info_ptr->species, NULL);
    if (retval) {
      goto lgen_to_bed_ret_1;
    }
    for (uii = 0; uii < unfiltered_marker_ct; uii++) {
      if (marker_id_map[uii] != -1) {
        marker_id_map[uii] = map_reverse[marker_id_map[uii]];
      }
    }
    wkspace_reset((unsigned char*)map_reverse);
  } else if (marker_exclude_ct) {
    // temporary storage of missing collapse info
    if (wkspace_alloc_i_checked(&indiv_id_map, unfiltered_marker_ct * sizeof(int))) {
      goto lgen_to_bed_ret_NOMEM;
    }
    marker_idx = 0;
    for (uii = 0; uii < marker_ct; uii++) {
      marker_idx = next_non_set_unsafe(marker_exclude, marker_idx);
      indiv_id_map[marker_idx] = uii;
      marker_idx++;
    }
    for (uii = 0; uii < unfiltered_marker_ct; uii++) {
      ii = marker_id_map[uii];
      if (is_set(marker_exclude, ii)) {
	marker_id_map[uii] = -1;
      } else {
	marker_id_map[uii] = indiv_id_map[ii];
      }
    }
    wkspace_reset((unsigned char*)indiv_id_map);
  }
  fclose_null(&infile);

  memcpy(name_end, ".fam", 5);
  if (fopen_checked(&infile, lgen_namebuf, "r")) {
    goto lgen_to_bed_ret_OPEN_FAIL;
  }
  retval = load_fam(infile, MAXLINELEN, 1, 1, 1, 1, 1, missing_pheno, intlen(missing_pheno), affection_01, &indiv_ct, &person_ids, &max_person_id_len, &paternal_ids, &max_paternal_id_len, &maternal_ids, &max_maternal_id_len, &sex_info, &affection, &pheno_c, &pheno_d, &founder_info, &indiv_exclude, 1, NULL, NULL, NULL);
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
  if (wkspace_alloc_c_checked(&marker_alleles, 2 * unfiltered_marker_ct)) {
    goto lgen_to_bed_ret_NOMEM;
  }
  if (wkspace_alloc_ui_checked(&marker_allele_cts, 2 * unfiltered_marker_ct * sizeof(int))) {
    goto lgen_to_bed_ret_NOMEM;
  }
  memset(marker_alleles, 0, 2 * marker_ct);
  fill_uint_zero(marker_allele_cts, 2 * marker_ct);
  indiv_ct4 = (indiv_ct + 3) / 4;
  if (wkspace_alloc_uc_checked(&writebuf, ((unsigned long)marker_ct) * indiv_ct4)) {
    logprint("Error: Very large .lgen -> .bed conversions are not yet supported.  Try this\nwith more memory (use --memory and/or a better machine).\n");
    goto lgen_to_bed_ret_CALC_NOT_YET_SUPPORTED;
  }
  if (indiv_ct % 4) {
    ucc = 0x15 >> (6 - 2 * indiv_ct);
    for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
      memset(&(writebuf[marker_idx * indiv_ct4]), 0x55, indiv_ct4 - 1);
      writebuf[(marker_idx + 1) * indiv_ct4 - 1] = ucc;
    }
  } else {
    memset(writebuf, 0x55, marker_ct * indiv_ct4);
  }
  // PLINK reports an error whenever there are 3+ alleles at one locus, so
  // backwards compatibility does not mandate that we worry about that case.
  // Thus we just use the obvious one-pass load.
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
  while (fgets(tbuf, MAXLINELEN, infile)) {
    if (is_eoln(*tbuf)) {
      continue;
    }
    cptr = skip_initial_spaces(tbuf);
    cptr2 = next_item(cptr);
    cptr3 = next_item(cptr2);
    cptr4 = item_end(cptr3);
    if (!cptr4) {
      goto lgen_to_bed_ret_INVALID_FORMAT;
    }
    cptr5 = skip_initial_spaces(cptr4);
    cptr6 = skip_initial_spaces(&(cptr5[1]));
    if (no_more_items(cptr6)) {
      goto lgen_to_bed_ret_INVALID_FORMAT;
    }
    ii = bsearch_fam_indiv(id_buf, sorted_indiv_ids, max_person_id_len, indiv_ct, cptr, cptr2);
    if (ii == -1) {
      cptr[strlen_se(cptr)] = '\0';
      cptr2[strlen_se(cptr2)] = '\0';
      sprintf(logbuf, "Error: Person %s %s in .lgen file but missing from .fam.\n", cptr, cptr2);
      logprintb();
      goto lgen_to_bed_ret_INVALID_FORMAT_2;
    }
    indiv_idx = indiv_id_map[ii];
    ulii = (unsigned long)(cptr4 - cptr3);
    memcpy(id_buf, cptr3, ulii);
    id_buf[ulii] = '\0';
    ii = bsearch_str(id_buf, sorted_marker_ids, max_marker_id_len, 0, unfiltered_marker_ct - 1);
    if (ii == -1) {
      sprintf(logbuf, "Error: Marker %s in .lgen file but missing from .map.\n", id_buf);
      logprintb();
      goto lgen_to_bed_ret_INVALID_FORMAT_2;
    }
    jj = marker_id_map[ii];
    if (jj != -1) {
      marker_idx = jj;
      cc = marker_alleles[2 * marker_idx + 1]; // A2
      if (cc == '\0') {
	cc = *cptr5;
	marker_alleles[2 * marker_idx + 1] = cc;
	if (*cptr6 == cc) {
	  uii = 2;
	} else {
	  uii = 1;
	  marker_alleles[2 * marker_idx] = *cptr6;
	}
      } else {
	cc2 = marker_alleles[2 * marker_idx];
	if (cc2 == '\0') {
	  if (*cptr5 == cc) {
	    if (*cptr6 == cc) {
	      uii = 2;
	    } else {
	      uii = 1;
	      marker_alleles[2 * marker_idx] = *cptr6;
	    }
	  } else {
	    cc2 = *cptr5;
	    marker_alleles[2 * marker_idx] = cc2;
	    if (*cptr6 == cc) {
	      uii = 1;
	    } else if (*cptr6 == cc2) {
	      uii = 0;
	    } else {
	      goto lgen_to_bed_ret_INVALID_FORMAT_3;
	    }
	  }
	} else {
	  if (*cptr5 == cc) {
	    uii = 1;
	  } else if (*cptr5 == cc2) {
	    uii = 0;
	  } else {
	    goto lgen_to_bed_ret_INVALID_FORMAT_3;
	  }
	  if (*cptr6 == cc) {
	    uii++;
	  } else if (*cptr6 != cc2) {
	    goto lgen_to_bed_ret_INVALID_FORMAT_3;
	  }
	}
      }
      marker_allele_cts[2 * marker_idx] += 2 - uii;
      marker_allele_cts[2 * marker_idx + 1] += uii;
      if (uii) {
	uii++;
      }
      ulii = marker_idx * indiv_ct4 + (indiv_idx / 4);
      ujj = (indiv_idx % 4) * 2;
      writebuf[ulii] = (writebuf[ulii] & (~(3 << ujj))) | (uii << ujj);
    }
    if (ftello(infile) >= lgen_next_thresh) {
      uii = (ftello(infile) * 100) / lgen_size;
      if (pct < 10) {
	printf("\b\b%u%%", uii);
      } else {
	printf("\b\b\b%u%%", uii);
      }
      fflush(stdout);
      pct = uii;
      lgen_next_thresh = ((pct + 1) * lgen_size) / 100;
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
  itable_short[0] = 3;
  itable_short[1] = 1;
  itable_short[2] = 2;
  itable_short[3] = 0;
  ucptr = itable;
  for (uii = 0; uii < 4; uii++) {
    for (ujj = 0; ujj < 4; ujj++) {
      for (ukk = 0; ukk < 4; ukk++) {
	umm = itable_short[uii] * 64 + itable_short[ujj] * 16 + itable_short[ukk] * 4;
	*ucptr++ = umm + 3;
	*ucptr++ = umm + 1;
	*ucptr++ = umm + 2;
	*ucptr++ = umm;
      }
    }
  }
  ukk = indiv_ct / 4;
  umm = indiv_ct % 4;
  ucptr = itable2;
  if (umm == 3) {
    for (uii = 0; uii < 4; uii++) {
      for (ujj = 0; ujj < 4; ujj++) {
	ukk = itable_short[uii] * 16 + itable_short[ujj] * 4;
	*ucptr++ = ukk + 3;
	*ucptr++ = ukk + 1;
	*ucptr++ = ukk + 2;
	*ucptr++ = ukk;
      }
    }
  } else if (umm == 2) {
    for (uii = 0; uii < 4; uii++) {
      ujj = itable_short[uii] * 4;
      *ucptr++ = ujj + 3;
      *ucptr++ = ujj + 1;
      *ucptr++ = ujj + 2;
      *ucptr++ = ujj;
    }
  } else if (umm == 1) {
    *ucptr++ = 3;
    *ucptr++ = 1;
    *ucptr++ = 2;
    *ucptr++ = 0;
  }
  for (uii = 0; uii < marker_ct; uii++) {
    if (marker_allele_cts[uii * 2] > marker_allele_cts[uii * 2 + 1]) {
      ucptr = &(writebuf[uii * indiv_ct4]);
      for (ujj = 0; ujj < ukk; ujj++) {
	*ucptr = itable[*ucptr];
	ucptr++;
      }
      if (umm) {
        *ucptr = itable2[*ucptr];
	ucptr++;
      }
      cc = marker_alleles[uii * 2];
      marker_alleles[uii * 2] = marker_alleles[uii * 2 + 1];
      marker_alleles[uii * 2 + 1] = cc;
    }
  }
  if (fwrite_checked(writebuf, ((unsigned long)marker_ct) * indiv_ct4, outfile)) {
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
      marker_alleles[ujj] = '0';
    }
  }
  uii = 0;
  marker_idx = 0;
  while (fgets(tbuf, MAXLINELEN, infile)) {
    if (is_eoln(*tbuf)) {
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
    sprintf(cptr, "\t%c\t%c\n", marker_alleles[marker_idx * 2], marker_alleles[marker_idx * 2 + 1]);
    ulii = strlen(cptr) + (unsigned long)(cptr - tbuf);
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
      if (is_eoln(*tbuf)) {
	continue;
      }
      ulii = strlen(tbuf);
      if (tbuf[ulii - 1] != '\n') {
	tbuf[ulii++] = '\n';
	tbuf[ulii] = '\0';
      }
      if (fwrite_checked(tbuf, ulii, outfile)) {
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
  lgen_to_bed_ret_INVALID_FORMAT:
    logprint("Error: Improperly formatted .lgen file.\n");
  lgen_to_bed_ret_INVALID_FORMAT_2:
    retval = RET_INVALID_FORMAT;
    break;
  lgen_to_bed_ret_INVALID_FORMAT_3:
    sprintf(logbuf, "Error: Marker %s in .lgen file has 3+ different alleles.\n", id_buf);
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
 lgen_to_bed_ret_1:
  if (infile) {
    fclose(infile);
  }
  if (outfile) {
    fclose(outfile);
  }
  return retval;
}

inline unsigned int update_alleles_and_cts(char* alleles, unsigned int* allele_cts, char cc) {
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

void transposed_to_bed_print_pct(unsigned int pct) {
  printf("Processing .tped file... %u%%", pct);
  fflush(stdout);
}

int transposed_to_bed(char* tpedname, char* tfamname, char* outname, char missing_geno, int species) {
  FILE* infile = NULL;
  FILE* bimfile = NULL;
  FILE* outfile = NULL;
  char* outname_end = (char*)memchr(outname, 0, FNAMESIZE);
  unsigned int indiv_ct = 0;
  unsigned int no_extra_cols = 1;
  int retval = 0;
  unsigned int pct = 0;
  unsigned int marker_idx = 0;
  long long last_mapval = 0;
  int map_is_unsorted = 0;
  unsigned long max_marker_id_len = 0;
  unsigned int indiv_ct4;
  unsigned int indiv_idx;
  unsigned long ulii;
  unsigned int uii;
  unsigned int ujj;
  int ii;
  // may as well not error out on triallelic (or quadallelic) input.  Instead,
  // print notes when it happens, and truncate.
  char alleles[4];
  char orig_alleles[4];
  unsigned int allele_cts[4];
  unsigned char* writebuf;
  unsigned char* prewritebuf;
  unsigned char writemap[17];
  unsigned char ucc;
  unsigned int max_load;
  char* loadbuf;
  char* cptr;
  char* cptr2;
  char* cptr3;
  char* cptr4;
  unsigned char* ucptr;
  unsigned char* ucptr2;
  char cc;
  char cc2;
  long long tped_size;
  long long tped_next_thresh;
  long long cur_mapval;
  long long* mapvals;
  unsigned int marker_ct;
  unsigned int marker_uidx;
  unsigned int* map_reverse;
  long long* ll_buf;
  int* pos_buf;
  char* marker_ids;
  unsigned int chrom_start[MAX_POSSIBLE_CHROM + 1];
  unsigned int chrom_id[MAX_POSSIBLE_CHROM];
  unsigned int cur_chrom;
  unsigned int chrom_ct;
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
    if (is_eoln(*tbuf)) {
      continue;
    }
    ulii = strlen(tbuf);
    if (tbuf[ulii - 1] != '\n') {
      tbuf[ulii++] = '\n';
    }
    if (fwrite_checked(tbuf, ulii, outfile)) {
      goto transposed_to_bed_ret_WRITE_FAIL;
    }
    indiv_ct++;
  }
  if (!feof(infile)) {
    goto transposed_to_bed_ret_READ_FAIL;
  }
  if (!indiv_ct) {
    logprint("\rError: No individuals in .tfam file.\n");
    goto transposed_to_bed_ret_INVALID_FORMAT_2;
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
  mapvals = (long long*)wkspace_base;
  loadbuf = (char*)(&wkspace_base[sizeof(long long)]);
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
    if (is_eoln(*loadbuf)) {
      continue;
    }
    memset(alleles, 0, 4);
    allele_cts[0] = 0;
    allele_cts[1] = 0;
    allele_cts[2] = 0;
    allele_cts[3] = 0;
    cptr = skip_initial_spaces(loadbuf);
    cptr2 = next_item(cptr);
    cptr3 = next_item_mult(cptr2, 2);
    cptr4 = next_item(cptr3);
    if (no_more_items(cptr4)) {
      goto transposed_to_bed_ret_INVALID_FORMAT;
    }
    ii = atoi(cptr);
    uii = strlen_se(cptr2) + 1;
    if (uii > max_marker_id_len) {
      max_marker_id_len = uii;
    }
    if (*cptr3 == '-') {
      continue;
    }
    cur_mapval = (((long long)ii) << 32) | atoi(cptr3);
    mapvals[marker_idx++] = cur_mapval;
    if (last_mapval > cur_mapval) {
      map_is_unsorted = 1;
    } else {
      last_mapval = cur_mapval;
    }
    for (uii = 0; uii < 3; uii++) {
      cptr2 = item_endnn2(cptr);
      *cptr2++ = '\t';
      if (fwrite_checked(cptr, (cptr2 - cptr), bimfile)) {
	goto transposed_to_bed_ret_WRITE_FAIL;
      }
      cptr = skip_initial_spaces(cptr2);
    }
    cptr2 = item_endnn2(cptr);
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

    if (no_extra_cols && (!is_eoln(*cptr2)) && (*cptr2 != '\0')) {
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
      if (allele_cts[3]) {
	sprintf(logbuf, "Note: Marker %s is quadallelic.  Truncating.\n", cptr);
      } else {
	sprintf(logbuf, "Note: Marker %s is triallelic.  Truncating.\n", cptr);
      }
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
    if (ftello(infile) >= tped_next_thresh) {
      uii = (ftello(infile) * 100) / tped_size;
      if (pct < 10) {
	printf("\b\b%u%%", uii);
      } else {
	printf("\b\b\b%u%%", uii);
      }
      fflush(stdout);
      pct = uii;
      tped_next_thresh = ((pct + 1) * tped_size) / 100;
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
    mapvals = (long long*)wkspace_alloc(marker_ct * sizeof(long long));

    if (wkspace_alloc_ll_checked(&ll_buf, marker_ct * sizeof(long long)) ||
        wkspace_alloc_i_checked(&pos_buf, marker_ct * sizeof(int)) ||
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
    map_reverse = (unsigned int*)mapvals;
    for (uii = 0; uii < chrom_ct; uii++) {
      cur_chrom = chrom_id[uii];
      ujj = chrom_start[uii + 1];
      for (; marker_idx < ujj; marker_idx++) {
	marker_uidx = (unsigned int)ll_buf[marker_idx];
	if (fprintf(outfile, "%u\t%s\t%g\t%u\t%c\t%c\n", cur_chrom, &(marker_ids[marker_uidx * max_marker_id_len]), gd_vals[marker_uidx], (unsigned int)(ll_buf[marker_idx] >> 32), marker_alleles[2 * marker_uidx], marker_alleles[2 * marker_uidx + 1]) < 0) {
	  goto transposed_to_bed_ret_WRITE_FAIL;
	}
	map_reverse[marker_uidx] = marker_idx;
      }
    }
    if (fclose_null(&outfile)) {
      goto transposed_to_bed_ret_WRITE_FAIL;
    }

    wkspace_reset((unsigned char*)map_reverse);
    map_reverse = (unsigned int*)wkspace_alloc(marker_ct * sizeof(int));
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
        if (fseeko(infile, 3 + ((unsigned long long)marker_uidx) * indiv_ct4, SEEK_SET)) {
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
    putchar('\r');
    sprintf(logbuf, "Error: Missing entries at marker %s in .tped file.\n", cptr);
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  transposed_to_bed_ret_INVALID_FORMAT_4:
    cptr[strlen_se(cptr)] = '\0';
    putchar('\r');
    sprintf(logbuf, "Error: half-missing call at marker %s, indiv %u in .tped file.\n", cptr, indiv_idx);
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

int recode(int recode_modifier, FILE* bedfile, int bed_offset, FILE* famfile, FILE* bimfile, FILE** outfile_ptr, char* outname, char* outname_end, unsigned int unfiltered_marker_ct, unsigned long* marker_exclude, unsigned int marker_ct, unsigned int unfiltered_indiv_ct, unsigned long* indiv_exclude, unsigned int indiv_ct, char* marker_ids, unsigned long max_marker_id_len, char* marker_alleles, char* person_ids, unsigned int max_person_id_len, char* paternal_ids, unsigned int max_paternal_id_len, char* maternal_ids, unsigned int max_maternal_id_len, unsigned char* sex_info, char* pheno_c, double* pheno_d, double missing_phenod, char output_missing_geno, char* output_missing_pheno, double* set_allele_freqs, int binary_files) {
  unsigned int unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  unsigned char* wkspace_mark = wkspace_base;
  int affection = (pheno_c != NULL);
  int phenos_present = (affection || (pheno_d != NULL));
  char delimiter = (recode_modifier & RECODE_TAB)? '\t' : ' ';
  char delim2;
  unsigned int marker_uidx;
  unsigned int marker_idx;
  unsigned int indiv_uidx;
  unsigned int indiv_idx;
  unsigned char* loadbuf;
  char* writebuf;
  unsigned char* bufptr;
  char* wbufptr;
  char* cptr;
  unsigned char ucc;
  unsigned char ucc2;
  char cc;
  char cc2;
  char cur_mk_alleles[2];
  unsigned int pct;
  unsigned int loop_end;
  unsigned long ulii;
  unsigned int shiftval;
  char* mk_alleles;
  char* het_reverse;
  if (recode_modifier & RECODE_TRANSPOSE) {
    if (wkspace_alloc_c_checked(&writebuf, indiv_ct * 4)) {
      return RET_NOMEM;
    }
  } else if (recode_modifier & RECODE_LGEN) {
    if (wkspace_alloc_c_checked(&writebuf, 3)) {
      return RET_NOMEM;
    }
    if ((!binary_files) && ((recode_modifier & (RECODE_TAB | RECODE_DELIMX)) == (RECODE_TAB | RECODE_DELIMX))) {
      // corner case: would generate two different .fam files with colliding
      // names
      strcpy(outname_end, "-lf");
      outname_end = &(outname_end[3]);
    }
  } else {
    if (wkspace_alloc_c_checked(&writebuf, marker_ct * 4)) {
      return RET_NOMEM;
    }
  }
  if (wkspace_alloc_c_checked(&het_reverse, marker_ct)) {
    return RET_NOMEM;
  }
  if (recode_modifier & RECODE_12) {
    if (wkspace_alloc_c_checked(&mk_alleles, unfiltered_marker_ct * 2)) {
      return RET_NOMEM;
    }
  } else {
    mk_alleles = marker_alleles;
  }
  marker_uidx = 0;
  for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
    marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx);
    if (set_allele_freqs[marker_uidx] < 0.5) {
      het_reverse[marker_idx] = 1;
    } else {
      het_reverse[marker_idx] = 0;
    }
    if (recode_modifier & RECODE_12) {
      mk_alleles[2 * marker_uidx] = '1';
      mk_alleles[2 * marker_uidx + 1] = '2';
    }
    marker_uidx++;
  }
  if (fseeko(bedfile, bed_offset, SEEK_SET)) {
    return RET_READ_FAIL;
  }
  fflush(stdout);
  marker_uidx = 0;
  marker_idx = 0;
  loadbuf = wkspace_base;
  if (recode_modifier & RECODE_TRANSPOSE) {
    strcpy(outname_end, ".tped");
    if (fopen_checked(outfile_ptr, outname, "w")) {
      return RET_OPEN_FAIL;
    }
    sprintf(logbuf, "--recode to %s + .tfam... ", outname);
    logprintb();
    fputs("0%", stdout);
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
    rewind(bimfile);
    for (pct = 1; pct <= 100; pct++) {
      loop_end = (((unsigned long long)pct) * marker_ct) / 100;
      for (; marker_idx < loop_end; marker_idx++) {
	do {
	  if (fgets(tbuf, MAXLINELEN, bimfile) == NULL) {
	    return RET_READ_FAIL;
	  }
	} while (tbuf[0] <= ' ');
	if (is_set(marker_exclude, marker_uidx)) {
	  do {
	    do {
	      if (fgets(tbuf, MAXLINELEN, bimfile) == NULL) {
		return RET_READ_FAIL;
	      }
	    } while (tbuf[0] <= ' ');
	  } while (is_set(marker_exclude, ++marker_uidx));
	  if (fseeko(bedfile, bed_offset + (marker_uidx * unfiltered_indiv_ct4), SEEK_SET)) {
	    return RET_READ_FAIL;
	  }
	}

	wbufptr = tbuf;
	for (indiv_idx = 0; indiv_idx < 3; indiv_idx++) {
	  cptr = item_endnn2(wbufptr);
	  ulii = 1 + (unsigned long)(cptr - wbufptr);
	  *cptr = delimiter;
	  if (fwrite_checked(wbufptr, ulii, *outfile_ptr)) {
	    return RET_WRITE_FAIL;
	  }
	  wbufptr = skip_initial_spaces(cptr);
	}
	ulii = strlen_se(wbufptr);
	wbufptr[ulii] = delimiter;
	if (fwrite_checked(wbufptr, ulii + 1, *outfile_ptr)) {
	  return RET_WRITE_FAIL;
	}

	if (fread(loadbuf, 1, unfiltered_indiv_ct4, bedfile) < unfiltered_indiv_ct4) {
	  return RET_READ_FAIL;
	}
	indiv_uidx = 0;
	wbufptr = writebuf;
	cur_mk_alleles[0] = mk_alleles[2 * marker_uidx];
	cur_mk_alleles[1] = mk_alleles[2 * marker_uidx + 1];
	for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
	  indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx);
	  ucc = (loadbuf[indiv_uidx / 4] >> ((indiv_uidx % 4) * 2)) & 3;
	  if (ucc) {
	    if (ucc == 2) {
	      ucc2 = het_reverse[marker_idx];
	      *wbufptr = cur_mk_alleles[ucc2];
	      wbufptr[2] = cur_mk_alleles[1 ^ ucc2];
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
	  return RET_WRITE_FAIL;
	}
	marker_uidx++;
      }
      if (pct < 100) {
	if (pct > 10) {
	  printf("\b\b\b%u%%", pct);
	} else {
	  printf("\b\b%u%%", pct);
	}
	fflush(stdout);
      }
    }
  } else if (recode_modifier & RECODE_LGEN) {
    strcpy(outname_end, ".lgen");
    if (fopen_checked(outfile_ptr, outname, "w")) {
      return RET_OPEN_FAIL;
    }
    sprintf(logbuf, "--recode to %s + .map + .fam... ", outname);
    logprintb();
    fputs("0%", stdout);
    delim2 = delimiter;
    if (delimiter == ' ') {
      memcpy(writebuf, "  ", 3);
    } else {
      memcpy(writebuf, "\t", 2);
      if (!(recode_modifier & RECODE_DELIMX)) {
	delim2 = ' ';
      }
    }
    for (pct = 1; pct <= 100; pct++) {
      loop_end = (((unsigned long long)pct) * marker_ct) / 100;
      for (; marker_idx < loop_end; marker_idx++) {
	if (is_set(marker_exclude, marker_uidx)) {
	  marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx + 1);
	  if (fseeko(bedfile, bed_offset + (marker_uidx * unfiltered_indiv_ct4), SEEK_SET)) {
	    return RET_READ_FAIL;
	  }
	}
	if (fread(loadbuf, 1, unfiltered_indiv_ct4, bedfile) < unfiltered_indiv_ct4) {
	  return RET_READ_FAIL;
	}
	wbufptr = &(marker_ids[marker_uidx * max_marker_id_len]);
	indiv_uidx = 0;
	cur_mk_alleles[0] = mk_alleles[2 * marker_uidx];
	cur_mk_alleles[1] = mk_alleles[2 * marker_uidx + 1];
	for (indiv_idx = 0; indiv_idx < indiv_ct; indiv_idx++) {
	  indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx);
	  ucc = (loadbuf[indiv_uidx / 4] >> ((indiv_uidx % 4) * 2)) & 3;
	  if (ucc != 1) {
	    if (ucc == 2) {
	      ucc2 = het_reverse[marker_idx];
	      cc = cur_mk_alleles[ucc2];
	      cc2 = cur_mk_alleles[1 ^ ucc2];
	    } else if (ucc == 3) {
	      cc = cur_mk_alleles[1];
	      cc2 = cc;
	    } else {
	      cc = *cur_mk_alleles;
	      cc2 = *cur_mk_alleles;
	    }
	    cptr = &(person_ids[indiv_uidx * max_person_id_len]);
	    ulii = strlen_se(cptr);
	    if (fwrite_checked(cptr, ulii, *outfile_ptr)) {
	      return RET_WRITE_FAIL;
	    }
	    if (fprintf(*outfile_ptr, "%c%s%c%s%s%c%c%c\n", delimiter, &(cptr[ulii + 1]), delimiter, wbufptr, writebuf, cc, delim2, cc2) < 0) {
	      return RET_WRITE_FAIL;
	    }
	  }
	  indiv_uidx++;
	}
	marker_uidx++;
      }
      if (pct < 100) {
	if (pct > 10) {
	  printf("\b\b\b%u%%", pct);
	} else {
	  printf("\b\b%u%%", pct);
	}
	fflush(stdout);
      }
    }
  } else {
    strcpy(outname_end, ".ped");
    if ((long long)wkspace_left >= (long long)unfiltered_indiv_ct4 * marker_ct) {
      if (fopen_checked(outfile_ptr, outname, "w")) {
	return RET_OPEN_FAIL;
      }
      sprintf(logbuf, "--recode to %s + .map... ", outname);
      logprintb();
      fputs("0%", stdout);
      while (marker_idx < marker_ct) {
	if (is_set(marker_exclude, marker_uidx)) {
	  marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx + 1);
	  if (fseeko(bedfile, bed_offset + (marker_uidx * unfiltered_indiv_ct4), SEEK_SET)) {
	    return RET_READ_FAIL;
	  }
	}
	if (unfiltered_marker_ct - marker_uidx > marker_ct - marker_idx) {
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
      loadbuf = wkspace_base;
      indiv_uidx = 0;
      indiv_idx = 0;
      if ((recode_modifier & (RECODE_TAB | RECODE_DELIMX)) == RECODE_TAB) {
	// PLINK actually writes spaces between same-locus alleles even under
        // --tab.  We replicate this for compatibility.
	loop_end = marker_ct * 4;
	for (ulii = 1; ulii < loop_end; ulii += 4) {
	  writebuf[ulii] = ' ';
	  writebuf[ulii + 2] = '\t';
	}
      } else {
        memset(writebuf, delimiter, marker_ct * 4 - 1);
      }
      writebuf[marker_ct * 4 - 1] = '\n';
      for (pct = 1; pct <= 100; pct++) {
	loop_end = ((unsigned long long)pct * indiv_ct) / 100;
	for (; indiv_idx < loop_end; indiv_idx++) {
	  indiv_uidx = next_non_set_unsafe(indiv_exclude, indiv_uidx);
	  cptr = &(person_ids[indiv_uidx * max_person_id_len]);
	  ulii = strlen_se(cptr);
	  if (fwrite_checked(cptr, ulii, *outfile_ptr)) {
	    return RET_WRITE_FAIL;
	  }
	  if (fprintf(*outfile_ptr, "%c%s%c%s%c%s%c%d%c", delimiter, &(cptr[ulii + 1]), delimiter, paternal_ids? (&(paternal_ids[indiv_uidx * max_paternal_id_len])) : "0", delimiter,  maternal_ids? (&(maternal_ids[indiv_uidx * max_maternal_id_len])) : "0", delimiter, sex_info? sex_info[indiv_uidx] : 0, delimiter) < 0) {
	    return RET_WRITE_FAIL;
	  }
	  if (affection) {
	    cc = pheno_c[indiv_uidx];
	    if (cc == -1) {
	      if (fprintf(*outfile_ptr, "%s%c", output_missing_pheno, delimiter) < 0) {
		return RET_WRITE_FAIL;
	      }
	    } else {
	      if (fprintf(*outfile_ptr, "%d%c", cc + 1, delimiter) < 0) {
		return RET_WRITE_FAIL;
	      }
	    }
	  } else if ((!phenos_present) || (pheno_d[indiv_uidx] == missing_phenod)) {
	    if (fprintf(*outfile_ptr, "%s%c", output_missing_pheno, delimiter) < 0) {
	      return RET_WRITE_FAIL;
	    }
	  } else {
	    if (fprintf(*outfile_ptr, "%g%c", pheno_d[indiv_uidx], delimiter) < 0) {
	      return RET_WRITE_FAIL;
	    }
	  }
	  bufptr = &(loadbuf[indiv_idx / 4]);
	  wbufptr = writebuf;
	  shiftval = (indiv_idx % 4) * 2;
	  marker_uidx = 0;
	  for (marker_idx = 0; marker_idx < marker_ct; marker_idx++) {
	    marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx);
	    ucc = ((*bufptr) >> shiftval) & 3;
	    if (ucc) {
	      if (ucc == 2) {
		ucc2 = het_reverse[marker_idx];
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
	    return RET_WRITE_FAIL;
	  }
	  indiv_uidx++;
	}
	if (pct < 100) {
	  if (pct > 10) {
	    printf("\b\b\b%u%%", pct);
	  } else {
	    printf("\b\b%u%%", pct);
	  }
	  fflush(stdout);
	}
      }
    } else {
      logprint("Error: --recode does not yet support very large .bed files.\n");
      return RET_CALC_NOT_YET_SUPPORTED;
    }
  }
  fclose_null(outfile_ptr);

  if ((recode_modifier & RECODE_TRANSPOSE) || ((recode_modifier & RECODE_LGEN) && (binary_files || ((recode_modifier & (RECODE_TAB | RECODE_DELIMX)) == (RECODE_TAB | RECODE_DELIMX))))) {
    if (recode_modifier & RECODE_TRANSPOSE) {
      strcpy(outname_end, ".tfam");
    } else {
      strcpy(outname_end, ".fam");
      if ((recode_modifier & (RECODE_TAB | RECODE_DELIMX)) == RECODE_TAB) {
	delimiter = ' ';
      }
    }
    if (fopen_checked(outfile_ptr, outname, "w")) {
      return RET_OPEN_FAIL;
    }
    rewind(famfile);
    for (indiv_uidx = 0; indiv_uidx < unfiltered_indiv_ct; indiv_uidx++) {
      if (fgets(tbuf, MAXLINELEN, famfile) == NULL) {
	return RET_READ_FAIL;
      }
      if (is_set(indiv_exclude, indiv_uidx)) {
	continue;
      }
      wbufptr = tbuf;
      for (loop_end = 0; loop_end < 5; loop_end++) {
	cptr = item_endnn2(wbufptr);
	ulii = 1 + (unsigned long)(cptr - wbufptr);
	*cptr = delimiter;
	if (fwrite_checked(wbufptr, ulii, *outfile_ptr)) {
	  return RET_WRITE_FAIL;
	}
	wbufptr = skip_initial_spaces(cptr);
      }
      if (affection) {
	cc = pheno_c[indiv_uidx];
	if (cc == -1) {
	  if (fprintf(*outfile_ptr, "%s\n", output_missing_pheno) < 0) {
	    return RET_WRITE_FAIL;
	  }
	} else {
	  if (fprintf(*outfile_ptr, "%d\n", cc + 1) < 0) {
	    return RET_WRITE_FAIL;
	  }
	}
      } else if ((!phenos_present) || (pheno_d[indiv_uidx] == missing_phenod)) {
	if (fprintf(*outfile_ptr, "%s\n", output_missing_pheno) < 0) {
	  return RET_WRITE_FAIL;
	}
      } else {
	if (fprintf(*outfile_ptr, "%g\n", pheno_d[indiv_uidx]) < 0) {
	  return RET_WRITE_FAIL;
	}
      }
    }
    fclose_null(outfile_ptr);
  }

  if (!(recode_modifier & RECODE_TRANSPOSE)) {
    strcpy(outname_end, ".map");
    if (fopen_checked(outfile_ptr, outname, "w")) {
      return RET_OPEN_FAIL;
    }
    rewind(bimfile);
    if ((recode_modifier & (RECODE_TAB | RECODE_DELIMX)) == RECODE_DELIMX) {
      cc = ' ';
    } else {
      // PLINK does not use space delimiter here
      cc = '\t';
    }
    for (marker_uidx = 0; marker_uidx < unfiltered_marker_ct; marker_uidx++) {
      do {
	fgets(tbuf, MAXLINELEN, bimfile);
      } while (tbuf[0] <= ' ');
      if (is_set(marker_exclude, marker_uidx)) {
	continue;
      }
      wbufptr = tbuf;
      for (loop_end = 0; loop_end < 3; loop_end++) {
	cptr = item_endnn2(wbufptr);
	ulii = 1 + (unsigned long)(cptr - wbufptr);
	*cptr = cc;
	if (fwrite_checked(wbufptr, ulii, *outfile_ptr)) {
	  return RET_WRITE_FAIL;
	}
	wbufptr = skip_initial_spaces(cptr);
      }
      ulii = strlen_se(wbufptr);
      wbufptr[ulii] = '\n';
      if (fwrite_checked(wbufptr, ulii + 1, *outfile_ptr)) {
	return RET_WRITE_FAIL;
      }
    }
    fclose_null(outfile_ptr);
  }
  fputs("\b\b\b", stdout);
  logprint("done.\n");

  wkspace_reset(wkspace_mark);
  return 0;
}

typedef struct ll_entry_struct {
  struct ll_entry_struct* next;
  double pheno;
  unsigned int orig_order;
  char idstr[];
} Ll_entry;

typedef struct ll_entry2_struct {
  struct ll_entry2_struct* next;
  long long pos;
  double gd_val;
  char idstr[];
} Ll_entry2;

static inline int idmatch(char* idtab, char* id0, unsigned int id0_len_p1, char* id1, unsigned int id1_len_p1) {
  return (!(memcmp(idtab, id0, id0_len_p1) || memcmp(&(idtab[id0_len_p1]), id1, id1_len_p1)));
}

static inline unsigned int hashval(char* id1, unsigned int id1_len, char* id2, unsigned int id2_len) {
  // just interpret as little-endian number and take modulo HASHSIZE_S
  unsigned char* ucptr = (unsigned char*)id1;
  unsigned char* ucp_end = &(ucptr[id1_len]);
  unsigned int vv = *ucptr;
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

static inline unsigned int hashval2(char* idstr, unsigned int idlen) {
  unsigned char* ucptr = (unsigned char*)idstr;
  unsigned char* ucp_end = &(ucptr[idlen]);
  unsigned int vv = *ucptr;
  while (++ucptr != ucp_end) {
    vv = ((vv << 8) + (*ucptr)) % HASHSIZE;
  }
  return vv;
}

static inline Ll_entry* alloc_ll(unsigned long* topsize_ptr, unsigned int size) {
  unsigned long ts = *topsize_ptr + ((size + sizeof(Ll_entry) + 15) & (~15LU));
  if (ts > wkspace_left) {
    return NULL;
  } else {
    *topsize_ptr = ts;
    return (Ll_entry*)(&(wkspace_base[wkspace_left - ts]));
  }
}

static inline Ll_entry2* alloc_ll2(unsigned long* topsize_ptr, unsigned int size) {
  unsigned long ts = *topsize_ptr + ((size + sizeof(Ll_entry2) + 15) & (~15LU));
  if (ts > wkspace_left) {
    return NULL;
  } else {
    *topsize_ptr = ts;
    return (Ll_entry2*)(&(wkspace_base[wkspace_left - ts]));
  }
}

int merge_fam_id_scan(char* bedname, char* famname, unsigned int* max_person_id_len_ptr, unsigned int* max_person_full_len_ptr, int* is_dichot_pheno_ptr, Ll_entry** htable, unsigned long* topsize_ptr, unsigned long long* tot_indiv_ct_ptr, unsigned long* ped_buflen_ptr, unsigned int* cur_indiv_ct_ptr, unsigned int* orig_idx_ptr) {
  unsigned int max_person_id_len = *max_person_id_len_ptr;
  unsigned int max_person_full_len = *max_person_full_len_ptr;
  int is_dichot_pheno = *is_dichot_pheno_ptr;
  unsigned long topsize = *topsize_ptr;
  unsigned long long tot_indiv_ct = *tot_indiv_ct_ptr;
  unsigned int orig_idx = *orig_idx_ptr;
  unsigned int cur_indiv_ct = 0;
  FILE* infile = NULL;
  unsigned int text_file = 0;
  int retval = 0;
  unsigned int col1_len;
  unsigned int col2_len;
  unsigned int col3_len;
  unsigned int col4_len;
  unsigned int tot_len;
  unsigned long ulii;
  unsigned int uii;
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
    if ((!is_eoln(cc)) && (cc != '#')) {
      col1_end_ptr = item_endnn2(col1_start_ptr);
      col1_len = col1_end_ptr - col1_start_ptr;
      col2_start_ptr = skip_initial_spaces(col1_end_ptr);
      col2_end_ptr = item_endnn2(col2_start_ptr);
      col2_len = col2_end_ptr - col2_start_ptr;
      col3_start_ptr = skip_initial_spaces(col2_end_ptr);
      col4_start_ptr = item_endnn2(col3_start_ptr);
      col3_len = col4_start_ptr - col3_start_ptr;
      col4_start_ptr = skip_initial_spaces(col4_start_ptr);
      col5_start_ptr = item_endnn2(col4_start_ptr);
      col4_len = col5_start_ptr - col4_start_ptr;
      col5_start_ptr = skip_initial_spaces(col5_start_ptr);
      col6_start_ptr = item_endnn2(col5_start_ptr);
      uii = col6_start_ptr - col5_start_ptr;
      if (uii != 1) {
	*col5_start_ptr = '0';
      }
      col6_start_ptr = skip_initial_spaces(col6_start_ptr);
      if (no_more_items(col6_start_ptr)) {
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
	  // for future: add parental ID/sex merge (PLINK doesn't have it)
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
	ll_ptr = alloc_ll(&topsize, tot_len);
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
    sprintf(logbuf, "Error: No individuals in %s.\n", famname);
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

// side effect: initializes tbuf to first nonempty line of .map/.bim
int check_gd_col(FILE* bimfile, char* tbuf, unsigned int is_binary, unsigned int* gd_col_ptr) {
  char* bufptr;
  while (fgets(tbuf, MAXLINELEN, bimfile)) {
    if (is_eoln(*tbuf)) {
      continue;
    }
    bufptr = next_item_mult(skip_initial_spaces(tbuf), 2 + 2 * is_binary);
    if (no_more_items(bufptr)) {
      return -1;
    }
    if (no_more_items(next_item(bufptr))) {
      *gd_col_ptr = 0;
    }
    *gd_col_ptr = 1;
    return 0;
  }
  return -1;
}

int merge_bim_scan(char* bimname, unsigned int is_binary, unsigned long* max_marker_id_len_ptr, Ll_entry2** htable2, unsigned long* topsize_ptr, unsigned long long* tot_marker_ct_ptr, unsigned int* cur_marker_ct_ptr, int species) {
  unsigned long max_marker_id_len = *max_marker_id_len_ptr;
  unsigned long topsize = *topsize_ptr;
  unsigned long long tot_marker_ct = *tot_marker_ct_ptr;
  unsigned int cur_marker_ct = *cur_marker_ct_ptr;
  FILE* infile = NULL;
  int retval = 0;
  unsigned int gd_col;
  long long llxx;
  char* bufptr;
  char* bufptr2;
  char* bufptr3;
  int ii;
  int jj;
  unsigned int uii;
  unsigned int ujj;
  unsigned int ukk;
  Ll_entry2** ll_pptr;
  Ll_entry2* ll_ptr;
  double gd_val;
  char cc;
  char cc2;
  unsigned int allele_ct;
  char alleles[2];
  if (fopen_checked(&infile, bimname, "r")) {
    goto merge_bim_scan_ret_OPEN_FAIL;
  }
  if (check_gd_col(infile, tbuf, is_binary, &gd_col)) {
    goto merge_bim_scan_ret_INVALID_FORMAT_2;
  }
  do {
    if (is_eoln(*tbuf)) {
      continue;
    }
    bufptr = skip_initial_spaces(tbuf);
    ii = marker_code(species, bufptr);
    if (ii == -1) {
      sprintf(logbuf, "Error: Invalid chromosome index in %s.\n", bimname);
      goto merge_bim_scan_ret_INVALID_FORMAT;
    }
    bufptr = next_item(bufptr);
    bufptr2 = item_endl(bufptr);
    uii = bufptr2 - bufptr;
    bufptr2 = skip_initial_spaces(bufptr2);
    if (no_more_items(bufptr2)) {
      goto merge_bim_scan_ret_INVALID_FORMAT_2;
    }
    if (gd_col) {
      if (sscanf(bufptr2, "%lg", &gd_val) != 1) {
	gd_val = 0;
      }
      bufptr2 = next_item(bufptr2);
      if (no_more_items(bufptr2)) {
	goto merge_bim_scan_ret_INVALID_FORMAT_2;
      }
    }
    jj = atoi(bufptr2);
    if (jj > 0) {
      if (is_binary) {
	bufptr2 = next_item(bufptr2);
	bufptr3 = next_item(bufptr2);
	if (no_more_items(bufptr3)) {
	  goto merge_bim_scan_ret_INVALID_FORMAT_2;
	}
	cc = *bufptr2;
	cc2 = *bufptr3;
	if ((cc == cc2) && (cc == '0')) {
	  sprintf(logbuf, "Error: A1 and A2 alleles identical for a marker in %s.\n", bimname);
	  goto merge_bim_scan_ret_INVALID_FORMAT;
	}
      } else {
	cc = '0';
	cc2 = '0';
      }
      llxx = (((long long)ii) << 32) + jj;
      ujj = hashval2(bufptr, uii);
      ll_pptr = &(htable2[ujj]);
      ll_ptr = *ll_pptr;
      ukk = 1;
      bufptr[uii++] = '\0';
      while (ll_ptr) {
	if (!memcmp(&(ll_ptr->idstr[2]), bufptr, uii)) {
	  if (is_binary) {
	    bufptr2 = ll_ptr->idstr;
	    allele_ct = 0;
	    if (bufptr2[0] != '0') {
	      alleles[0] = bufptr2[0];
	      allele_ct++;
	    }
	    if (bufptr2[1] != '0') {
	      alleles[allele_ct++] = bufptr2[1];
	    }
	    if (cc2 != '0') {
	      for (ukk = 0; ukk < allele_ct; ukk++) {
		if (alleles[ukk] == cc2) {
		  break;
		}
	      }
	      if (ukk == allele_ct) {
		if (!allele_ct) {
		  bufptr2[1] = cc2;
		} else if (allele_ct == 1) {
		  if (bufptr2[0] == '0') {
		    bufptr2[0] = cc2;
		  } else {
		    bufptr2[1] = cc2;
		  }
		} else {
		  goto merge_bim_scan_ret_INVALID_FORMAT_3;
		}
		alleles[allele_ct++] = cc2;
	      }
	    }
	    if (cc != '0') {
	      for (ukk = 0; ukk < allele_ct; ukk++) {
		if (alleles[ukk] == cc) {
		  break;
		}
	      }
	      if (ukk == allele_ct) {
		if (!allele_ct) {
		  bufptr2[1] = cc;
		} else if (allele_ct == 1) {
		  if (bufptr2[0] == '0') {
		    bufptr2[0] = cc;
		  } else {
		    bufptr2[1] = cc;
		  }
		} else {
		  goto merge_bim_scan_ret_INVALID_FORMAT_3;
		}
		alleles[allele_ct++] = cc;
	      }
	    }
	    ukk = 1;
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
	ll_ptr = alloc_ll2(&topsize, uii + 2);
	ll_ptr->next = NULL;
	ll_ptr->pos = llxx;
	ll_ptr->gd_val = gd_val;
	ll_ptr->idstr[0] = cc;
	ll_ptr->idstr[1] = cc2;
	memcpy(&(ll_ptr->idstr[2]), bufptr, uii);
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
  *max_marker_id_len_ptr = max_marker_id_len;
  *topsize_ptr = topsize;
  *tot_marker_ct_ptr = tot_marker_ct;
  *cur_marker_ct_ptr = cur_marker_ct;
  
  while (0) {
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

static inline void merge_post_msort_update_maps(char* marker_ids, unsigned long max_marker_id_len, unsigned int* marker_map, unsigned int* pos_buf, long long* ll_buf, unsigned int* chrom_start, unsigned int chrom_ct, unsigned int* dedup_marker_ct_ptr, unsigned int merge_disallow_equal_pos) {
  // Input: ll_buf is a effectively sequence of sorted arrays (one per
  // chromosome) with base-pair positions in high 32 bits, and pre-sort indices
  // in low 32 bits.  Chromosome boundaries are stored in chrom_start[].
  // Pre-sort indices refer to marker ID ASCII ordering, which is how lookup
  // will be performed mid-merge.  There may be duplicate positions.

  // Result: Duplicates have been collapsed, with chrom_start[] updated.
  // pos_buf contains sorted base-pair positions, post-duplicate-removal.
  // marker_map[n] is the post-duplicate-removal position in all other arrays
  // of marker ID n by ASCII ordering.
  unsigned int read_pos = 0;
  unsigned int write_pos = 0; // may be lower than read_pos due to dups
  unsigned int chrom_idx;
  unsigned int chrom_read_end_idx;
  long long llxx;
  unsigned int prev_bp;
  unsigned int cur_bp;
  unsigned int presort_idx;
  for (chrom_idx = 0; chrom_idx < chrom_ct; chrom_idx++) {
    chrom_read_end_idx = chrom_start[chrom_idx + 1];
    // ll_buf has base-pair positions in high 32 bits, and pre-sort indices in
    // low 32 bits.
    llxx = ll_buf[read_pos++];
    prev_bp = (unsigned int)(llxx >> 32);
    pos_buf[write_pos] = prev_bp;
    marker_map[(unsigned int)llxx] = write_pos++;
    for (; read_pos < chrom_read_end_idx; read_pos++) {
      llxx = ll_buf[read_pos];
      presort_idx = (unsigned int)llxx;
      cur_bp = (unsigned int)(llxx >> 32);
      if (prev_bp == cur_bp) {
	sprintf(logbuf, "Warning: Markers %s and %s have the same position.\n", &(marker_ids[max_marker_id_len * presort_idx]), &(marker_ids[max_marker_id_len * ((unsigned int)ll_buf[read_pos - 1])]));
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

static inline int merge_must_track_write(int mm) {
  // modes 6 and 7 can be sped up with early pruning of nonoverlapping
  // markers/SNPs, but not worth complicating code for this.
  return (mm == 1) || (mm > 5) || (mm == 4);
}

static inline int merge_first_mode(int mm, unsigned int merge_disallow_equal_pos) {
  if (merge_disallow_equal_pos) {
    return (mm > 5)? 4 : mm;
  } else {
    return merge_must_track_write(mm)? 4 : 5;
  }
}

int merge_diff_print(FILE* outfile, char* idbuf, char* marker_id, char* person_id, unsigned char newval, unsigned char oldval, char* marker_alleles) {
  char* bufptr = item_endnn2(marker_id);
  unsigned int slen = bufptr - marker_id;
  char ma1[4];
  char ma2[4];
  memcpy(idbuf, marker_id, slen);
  idbuf[slen] = '\0';
  if (fprintf(outfile, "%20s", idbuf) < 0) {
    return -1;
  }
  ma1[0] = marker_alleles[0];
  ma1[1] = '0';
  ma1[2] = marker_alleles[0];
  ma1[3] = marker_alleles[1];
  ma2[0] = marker_alleles[0];
  ma2[1] = '0';
  ma2[2] = marker_alleles[1];
  ma2[3] = marker_alleles[1];
  bufptr = item_endnn2(person_id);
  slen = (bufptr++) - person_id;
  memcpy(idbuf, person_id, slen);
  idbuf[slen] = '\0';
  return (fprintf(outfile, " %20s %20s      %c/%c      %c/%c \n", idbuf, bufptr, ma1[newval], ma2[newval], ma1[oldval], ma2[oldval]) < 0);
}

int merge_main(char* bedname, char* bimname, char* famname, unsigned int tot_indiv_ct, unsigned int tot_marker_ct, unsigned int dedup_marker_ct, unsigned int start_marker_idx, unsigned int marker_window_size, char* marker_alleles, char* marker_ids, unsigned long max_marker_id_len, char* person_ids, unsigned int max_person_id_len, unsigned int merge_nsort, unsigned int* indiv_nsmap, unsigned int* flex_map, unsigned int* marker_map, unsigned int* chrom_start, unsigned int* chrom_id, unsigned int chrom_ct, char* idbuf, unsigned char* readbuf, unsigned char* writebuf, unsigned int merge_mode, unsigned long* markbuf, FILE* outfile, unsigned long long* diff_total_overlap_ptr, unsigned long long* diff_not_both_genotyped_ptr, unsigned long long* diff_discordant_ptr, unsigned long ped_buflen, unsigned char* bmap_raw) {
  // flex_map maps individuals for binary filesets, and markers for text
  // filesets.
  unsigned int is_binary = famname? 1 : 0;
  FILE* bedfile = NULL;
  FILE* infile2 = NULL;
  int retval = 0;
  unsigned int tot_indiv_ct4 = (tot_indiv_ct + 3) / 4;
  unsigned int tot_indiv_ctl = (tot_indiv_ct + (BITCT - 1)) / BITCT;
  unsigned int end_marker_idx = start_marker_idx + marker_window_size;
  unsigned int marker_in_idx = 4294967295U; // overflow to zero on first add
  unsigned int last_marker_in_idx = 4294967294U;
  unsigned int cur_indiv_ct = 0;
  unsigned long* mbufptr = NULL; // merge mode 1, 4, 6, 7
  unsigned long long diff_total_overlap = 0;
  unsigned long long diff_not_both_genotyped = 0;
  unsigned long long diff_discordant = 0;
  unsigned int cur_indiv_ct4 = 0;
  unsigned int cur_indiv_ctd4 = 0;
  unsigned long uljj = 0;
  unsigned long* mbufptr2;
  unsigned char* bmap;
  unsigned int gd_col;
  char* bufptr;
  char* bufptr2;
  char* bufptr3;
  char* bufptr4;
  unsigned char* rbufptr;
  unsigned char* wbufptr;
  unsigned char* wbufptr2;
  unsigned long ulii;
  unsigned int uii;
  unsigned int ujj;
  unsigned int ukk;
  unsigned int umm;
  unsigned int unn;
  unsigned int indiv_idx;
  unsigned int marker_out_idx;
  unsigned char ucc;
  unsigned char ucc2;
  unsigned char ucc3;
  unsigned char ucc4;
  char cc;
  char cc2;
  char cc3;
  char cc4;
  int ii;
  if (is_binary) {
    if (fopen_checked(&infile2, famname, "r")) {
      goto merge_main_ret_OPEN_FAIL;
    }
    while (fgets(tbuf, MAXLINELEN, infile2)) {
      bufptr = skip_initial_spaces(tbuf);
      if (is_eoln(*bufptr)) {
	continue;
      }
      bufptr2 = item_endnn2(bufptr);
      bufptr3 = skip_initial_spaces(bufptr2);
      bufptr4 = item_endnn2(bufptr3);
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
    if (is_eoln(*tbuf)) {
      continue;
    }
    ++marker_in_idx;
    bufptr = next_item(skip_initial_spaces(tbuf));
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
    bufptr3 = item_endnn2(bufptr);
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
	if (fseeko(bedfile, 3 + (marker_in_idx * cur_indiv_ct4), SEEK_SET)) {
	  goto merge_main_ret_READ_FAIL;
	}
      }
      bufptr2 = next_item(bufptr2);
      bufptr3 = next_item(bufptr2);
      if (no_more_items(bufptr3)) {
	goto merge_main_ret_READ_FAIL;
      }
      ucc = *bufptr2;
      ucc2 = *bufptr3;
      if (((ucc != '0') && (ucc == marker_alleles[ii * 2 + 1])) || ((ucc2 == marker_alleles[ii * 2]) && (ucc2 != '0'))) {
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
	    ulii = 1LU << (ujj % BITCT);
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
	  ulii = 1LU << (ujj % BITCT);
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
	    ulii = 1LU << (ujj % BITCT);
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
	  ulii = 1LU << (ujj % BITCT);
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
	    if (mbufptr[ujj / BITCT] & (1LU << (ujj % BITCT))) {
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
		if (merge_diff_print(outfile, idbuf, bufptr, &(person_ids[ujj * max_person_id_len]), ucc2, ucc3, &(marker_alleles[ii * 2]))) {
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
	  if (mbufptr[ujj / BITCT] & (1LU << (ujj % BITCT))) {
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
	      if (merge_diff_print(outfile, idbuf, bufptr, &(person_ids[ujj * max_person_id_len]), ucc2, ucc3, &(marker_alleles[ii * 2]))) {
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
	    if (mbufptr[ujj / BITCT] & (1LU << (ujj % BITCT))) {
	      diff_total_overlap++;
	      ukk = (ujj % 4) * 2;
	      ucc2 = ucc & 3;
	      ucc3 = (wbufptr[ujj / 4] >> ukk) & 3;
	      if ((ucc2 == 1) || (ucc3 == 1)) {
		diff_not_both_genotyped++;
	      } else if (ucc2 != ucc3) {
		diff_discordant++;
		if (merge_diff_print(outfile, idbuf, bufptr, &(person_ids[ujj * max_person_id_len]), ucc2, ucc3, &(marker_alleles[ii * 2]))) {
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
	  if (mbufptr[ujj / BITCT] & (1LU << (ujj % BITCT))) {
	    diff_total_overlap++;
	    ukk = (ujj % 4) * 2;
	    ucc2 = ucc & 3;
	    ucc3 = (wbufptr[ujj / 4] >> ukk) & 3;
	    if ((ucc2 == 1) || (ucc3 == 1)) {
	      diff_not_both_genotyped++;
	    } else if (ucc2 != ucc3) {
	      diff_discordant++;
	      if (merge_diff_print(outfile, idbuf, bufptr, &(person_ids[ujj * max_person_id_len]), ucc2, ucc3, &(marker_alleles[ii * 2]))) {
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
      if ((is_eoln(cc)) || (cc == '#')) {
	continue;
      }
      bufptr2 = item_endnn2(bufptr);
      uii = (bufptr2 - bufptr);
      memcpy(idbuf, bufptr, uii);
      idbuf[uii] = '\t';
      bufptr3 = skip_initial_spaces(bufptr2);
      bufptr2 = item_endnn2(bufptr3);
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
	uljj = 1LU << (ii % BITCT);
      }
      for (marker_in_idx = 0; marker_in_idx < last_marker_in_idx; marker_in_idx++) {
	cc = *bufptr3;
	if (is_eoln(cc)) {
	  goto merge_main_ret_INVALID_FORMAT_2;
	}
	bufptr3 = skip_initial_spaces(&(bufptr3[1]));
	cc2 = *bufptr3;
	if (is_eoln(cc2)) {
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
	  cc3 = marker_alleles[uii * 2];
	  cc4 = marker_alleles[uii * 2 + 1];
	  if (cc == cc4) {
	    ucc2++;
	  } else if (cc != cc3) {
	    if (cc4 == '0') {
	      cc4 = cc;
	      ucc2++;
	      marker_alleles[uii * 2 + 1] = cc4;
	    } else if (cc3 == '0') {
	      cc3 = cc;
	      marker_alleles[uii * 2] = cc3;
	    } else {
	      goto merge_main_ret_INVALID_FORMAT;
	    }
	  }
	  if (cc2 == cc4) {
	    ucc2++;
	  } else if (cc2 != cc3) {
	    if (cc3 == '0') {
	      marker_alleles[uii * 2] = cc2;
	    } else if (cc4 == '0') {
	      // put this second since only way cc3 != '0' and cc4 == '0' is if
	      // it was specified that way in earlier binary file, and cc
	      // matches cc3.
	      marker_alleles[uii * 2 + 1] = cc2;
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
	      if (merge_diff_print(outfile, idbuf, &(marker_ids[uii * max_marker_id_len]), &(person_ids[((unsigned int)ii) * max_person_id_len]), ucc2, ucc3, &(marker_alleles[uii * 2]))) {
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
	      if (merge_diff_print(outfile, idbuf, &(marker_ids[uii * max_marker_id_len]), &(person_ids[((unsigned int)ii) * max_person_id_len]), ucc2, ucc3, &(marker_alleles[uii * 2]))) {
		goto merge_main_ret_WRITE_FAIL;
	      }
	    }
	  }
	}
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
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  merge_main_ret_INVALID_FORMAT_3:
    fill_idbuf_fam_indiv(idbuf, bufptr, ' ');
    sprintf(logbuf, "Error: Half-missing call in %s (indiv id %s, marker %u)", bedname, idbuf, marker_in_idx);
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  merge_main_ret_INVALID_FORMAT_2:
    fill_idbuf_fam_indiv(idbuf, bufptr, ' ');
    sprintf(logbuf, "Error: Line too short in %s (indiv id %s).\n", bedname, idbuf);
    logprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
  fclose_cond(bedfile);
  fclose_cond(infile2);
  return retval;
}

int merge_datasets(char* bedname, char* bimname, char* famname, char* outname, char* outname_end, char* mergename1, char* mergename2, char* mergename3, int calculation_type, int merge_type, int indiv_sort, int keep_allele_order, int species) {
  FILE* mergelistfile = NULL;
  FILE* outfile = NULL;
  unsigned char* wkspace_mark = wkspace_base;
  unsigned int max_person_id_len = 0;
  unsigned int max_person_full_len = 0;
  unsigned long max_marker_id_len = 0;
  int is_dichot_pheno = 1;
  int merge_mode = (merge_type & MERGE_MODE_MASK);
  unsigned int merge_nsort = ((!indiv_sort) || (indiv_sort == INDIV_SORT_NATURAL))? 1 : 0;
  unsigned int merge_disallow_equal_pos = (merge_type & MERGE_ALLOW_EQUAL_POS)? 0 : 1;
  Ll_entry** htable = (Ll_entry**)(&(wkspace_base[wkspace_left - HASHMEM_S]));
  Ll_entry2** htable2 = (Ll_entry2**)(&(wkspace_base[wkspace_left - HASHMEM]));
  unsigned long ped_buflen = MAXLINELEN;
  char* pheno_c = NULL;
  double* pheno_d = NULL;
  unsigned int* indiv_nsmap = NULL;
  unsigned int max_cur_indiv_ct = 0;
  unsigned int max_cur_marker_text_ct = 0;
  unsigned long* markbuf = NULL; // needed for merge modes 1, 4, 6, 7
  unsigned long long diff_total_overlap = 0;
  unsigned long long diff_not_both_genotyped = 0;
  unsigned long long diff_discordant = 0;
  unsigned int orig_idx = 0;
  unsigned int* map_reverse = NULL;
  unsigned long* reversed = NULL;
  unsigned char bmap_raw[576];
  unsigned char* bmap;
  unsigned char* bmap2;
  unsigned char* ucptr;
  unsigned char* ucptr_end;
  unsigned long* pcptr;
  unsigned long markers_per_pass;
  unsigned int pass_ct;
  unsigned long topsize;
  char* person_ids;
  char* person_fids;
  char* marker_ids;
  // N.B. marker_alleles are ordered by marker_id instead of position
  char* marker_alleles;
  unsigned int* marker_map;
  unsigned int* flex_map;
  double* gd_vals;
  unsigned int* pos_buf;
  long long* ll_buf;
  unsigned long mlpos;
  unsigned long merge_ct;
  char* idbuf;
  char* mergelist_buf;
  char** mergelist_bed;
  char** mergelist_bim;
  char** mergelist_fam;
  unsigned int cur_indiv_ct;
  unsigned int cur_marker_ct;
  unsigned int tot_indiv_ct;
  unsigned int tot_indiv_ct4;
  unsigned int tot_marker_ct;
  unsigned int dedup_marker_ct;
  unsigned long long ullxx;
  long long llxx;
  unsigned long ulii;
  unsigned long uljj;
  unsigned long ulkk;
  unsigned int uii;
  unsigned int ujj;
  unsigned int ukk;
  unsigned int umm;
  char* bufptr;
  char* bufptr2;
  char* bufptr3;
  char* bufptr4;
  Ll_entry* ll_ptr;
  Ll_entry2* ll_ptr2;
  unsigned int chrom_start[MAX_POSSIBLE_CHROM + 1];
  unsigned int chrom_id[MAX_POSSIBLE_CHROM];
  unsigned int chrom_ct;
  unsigned char* readbuf;
  unsigned char* writebuf;
  unsigned char* ubufptr;
  char cc;
  char cc2;
  unsigned char ucc;
  int retval;

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
      if (no_more_items(bufptr)) {
	continue;
      }
      bufptr2 = next_item_mult(bufptr, 3);
      if (!no_more_items(bufptr2)) {
	logprint("Error: More than three items on --merge-list line.\n");
        goto merge_datasets_ret_INVALID_FORMAT;
      }
      if (no_more_items(next_item(bufptr))) {
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
	} while (!no_more_items(bufptr));
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
    mergelist_bed = (char**)wkspace_alloc(merge_ct * sizeof(long));
    mergelist_bim = (char**)wkspace_alloc(merge_ct * sizeof(long));
    mergelist_fam = (char**)wkspace_alloc(merge_ct * sizeof(long));
    if (wkspace_alloc_c_checked(&mergelist_buf, (unsigned long)ullxx)) {
      goto merge_datasets_ret_NOMEM;
    }
    rewind(mergelistfile);
    bufptr4 = mergelist_buf;
    mlpos = 1;
    while (fgets(tbuf, MAXLINELEN, mergelistfile)) {
      bufptr = skip_initial_spaces(tbuf);
      if (no_more_items(bufptr)) {
	continue;
      }
      bufptr2 = item_endnn(bufptr);
      ulii = (bufptr2 - bufptr);
      bufptr3 = skip_initial_spaces(bufptr2);
      if (no_more_items(bufptr3)) {
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
	if (no_more_items(bufptr)) {
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
    mergelist_bed = (char**)wkspace_alloc(2 * sizeof(long));
    mergelist_bim = (char**)wkspace_alloc(2 * sizeof(long));
    mergelist_fam = (char**)wkspace_alloc(2 * sizeof(long));
    mergelist_bed[1] = mergename1;
    mergelist_bim[1] = mergename2;
    mergelist_fam[1] = (merge_type & MERGE_BINARY)? mergename3 : NULL;
  }
  mergelist_bed[0] = bedname;
  mergelist_bim[0] = bimname;
  mergelist_fam[0] = (famname[0])? famname : NULL;

  // ID counting/duplicate detection strategy:
  // - We do NOT want to scan through .ped files any more times than absolutely
  // necessary.  So we actually use *gasp* a hash table here.
  // - The hash table is positioned at the FAR end of wkspace, automatically
  // sized to ~4MB (or ~2MB on 32-bit systems).  IDs are then stored
  // backwards from there.  This simplifies copying into a sorted list.
  if (wkspace_left < HASHSIZE_S * sizeof(long)) {
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
    logprint("Error: Too many individuals (max 2147483647).\n");
    goto merge_datasets_ret_INVALID_FORMAT;
  }
#else
  // avoid integer overflow in wkspace_alloc calls
  if (ullxx * max_person_full_len > 2147483647) {
    logprint("Error: Too many individuals for 32-bit WDIST.\n");
    goto merge_datasets_ret_INVALID_FORMAT;
  }
#endif
  tot_indiv_ct = ullxx;
  fill_bmap_raw(bmap_raw, tot_indiv_ct % 4);
  bmap = &(bmap_raw[256]);
  bmap2 = &(bmap_raw[512]);
  if (indiv_sort == INDIV_SORT_NONE) {
    if (wkspace_alloc_ui_checked(&indiv_nsmap, tot_indiv_ct * sizeof(int))) {
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
    wkspace_alloc_c_checked(&pheno_c, tot_indiv_ct);
  } else {
    if (wkspace_left < topsize + tot_indiv_ct * sizeof(double)) {
      goto merge_datasets_ret_NOMEM;
    }
    wkspace_alloc_d_checked(&pheno_d, tot_indiv_ct * sizeof(double));
  }
  if (indiv_sort == INDIV_SORT_NONE) {
    if (wkspace_alloc_ui_checked(&map_reverse, tot_indiv_ct * sizeof(int))) {
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
	      pheno_c[ujj] = -1;
	    } else {
	      pheno_c[ujj] = ((ll_ptr->pheno == 1)? 0 : 1);
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
    if (qsort_ext(person_fids, tot_indiv_ct, max_person_full_len, strcmp_deref, (char*)indiv_nsmap, sizeof(int))) {
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
	      pheno_c[ulii] = -1;
	    } else if (ll_ptr->pheno == 1) {
	      pheno_c[ulii] = 0;
	    } else {
	      pheno_c[ulii] = 1;
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
      if (qsort_ext(person_fids, tot_indiv_ct, max_person_full_len, merge_nsort? strcmp_natural_deref : strcmp_deref, pheno_c, 1)) {
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
	  cc = pheno_c[ulii];
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
	  cc = pheno_c[ulii];
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
    retval = merge_bim_scan(mergelist_bim[mlpos], (mergelist_fam[mlpos])? 1 : 0, &max_marker_id_len, htable2, &topsize, &ullxx, &cur_marker_ct, species);
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
  if (wkspace_alloc_c_checked(&marker_ids, max_marker_id_len * tot_marker_ct) ||
      wkspace_alloc_c_checked(&marker_alleles, tot_marker_ct * 2) ||
      wkspace_alloc_ui_checked(&marker_map, tot_marker_ct * sizeof(int)) ||
      wkspace_alloc_d_checked(&gd_vals, tot_marker_ct * sizeof(double)) ||
      wkspace_alloc_ui_checked(&pos_buf, tot_marker_ct * sizeof(int)) ||
      wkspace_alloc_ll_checked(&ll_buf, tot_marker_ct * sizeof(long long))) {
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
	strcpy(&(marker_ids[ulii * max_marker_id_len]), &(ll_ptr2->idstr[2]));
        ulii++;
	ll_ptr2 = ll_ptr2->next;
      } while (ll_ptr2);
    }
  }
  if (qsort_ext(marker_ids, tot_marker_ct, max_marker_id_len, strcmp_deref, (char*)pos_buf, sizeof(int))) {
    goto merge_datasets_ret_NOMEM2;
  }
  // pos_buf[n] contains the position of lexicographic marker #n in the hash
  // table.  invert this map, then traverse the hash table.
  for (uii = 0; uii < tot_marker_ct; uii++) {
    marker_map[pos_buf[uii]] = uii;
  }
  wkspace_left += topsize;
  ulii = 0;
  for (uii = 0; uii < HASHSIZE; uii++) {
    if (htable2[uii]) {
      ll_ptr2 = htable2[uii];
      do {
        ujj = marker_map[ulii++];
	llxx = ll_ptr2->pos;
	pos_buf[ujj] = (unsigned int)llxx;
	memcpy(&(marker_alleles[ujj * 2]), ll_ptr2->idstr, 2);
	gd_vals[ujj] = ll_ptr2->gd_val;
	ll_buf[ujj] = (llxx & 0xffffffff00000000LL) | ujj;
	ll_ptr2 = ll_ptr2->next;
      } while (ll_ptr2);
    }
  }
  sort_marker_chrom_pos(ll_buf, tot_marker_ct, (int*)pos_buf, chrom_start, chrom_id, &chrom_ct);
  merge_post_msort_update_maps(marker_ids, max_marker_id_len, marker_map, pos_buf, ll_buf, chrom_start, chrom_ct, &dedup_marker_ct, merge_disallow_equal_pos);
  wkspace_reset((unsigned char*)ll_buf);

  tot_indiv_ct4 = (tot_indiv_ct + 3) / 4;

  if (!keep_allele_order) {
    ulii = (tot_marker_ct + (BITCT - 1)) / BITCT;
    if (wkspace_alloc_ul_checked(&reversed, ulii * sizeof(long))) {
      goto merge_datasets_ret_NOMEM;
    }
    fill_ulong_zero(reversed, ulii);
  }
  if (wkspace_alloc_ui_checked(&flex_map, MAXV(max_cur_indiv_ct, max_cur_marker_text_ct) * sizeof(int)) ||
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
    markers_per_pass = wkspace_left / (3 * sizeof(long) * ulii);
    if (markers_per_pass > dedup_marker_ct) {
      uii = dedup_marker_ct;
    } else {
      uii = markers_per_pass;
    }
    markbuf = (unsigned long*)wkspace_alloc(uii * ulii * sizeof(long));
  } else {
    markers_per_pass = wkspace_left / tot_indiv_ct4;
  }
  pass_ct = 1 + ((dedup_marker_ct - 1) / markers_per_pass);

  writebuf = wkspace_base;
  pcptr = (unsigned long*)wkspace_base;
  if (merge_mode < 6) {
    memcpy(outname_end, ".bed", 5);
    if (fopen_checked(&outfile, outname, "wb")) {
      goto merge_datasets_ret_OPEN_FAIL;
    }
    if (fwrite_checked("l\x1b\x01", 3, outfile)) {
      goto merge_datasets_ret_WRITE_FAIL;
    }
    if (pass_ct == 1) {
      sprintf(logbuf, "Performing single-pass merge (%u individuals, %u markers).\n", tot_indiv_ct, dedup_marker_ct);
    } else {
      sprintf(logbuf, "Performing %u-pass merge (%u indivs, %lu/%u markers per pass).\n", pass_ct, tot_indiv_ct, markers_per_pass, dedup_marker_ct);
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
      memset(writebuf, 0x55, ((unsigned long)ujj) * tot_indiv_ct4);
    }
    if (merge_must_track_write(merge_mode)) {
      fill_ulong_zero(markbuf, ujj * ulii);
    }
    for (mlpos = 0; mlpos < merge_ct; mlpos++) {
      retval = merge_main(mergelist_bed[mlpos], mergelist_bim[mlpos], mergelist_fam[mlpos], tot_indiv_ct, tot_marker_ct, dedup_marker_ct, uii * markers_per_pass, ujj, marker_alleles, marker_ids, max_marker_id_len, person_ids, max_person_id_len, merge_nsort, indiv_nsmap, flex_map, marker_map, chrom_start, chrom_id, chrom_ct, idbuf, readbuf, writebuf, mlpos? merge_mode : merge_first_mode(merge_mode, merge_disallow_equal_pos), markbuf, outfile, &diff_total_overlap, &diff_not_both_genotyped, &diff_discordant, ped_buflen, bmap_raw);
      if (retval) {
	goto merge_datasets_ret_1;
      }
      if (mlpos != merge_ct - 1) {
        printf("\rPass %u: fileset #%lu complete.", uii + 1, mlpos + 1);
	fflush(stdout);
      }
    }
    if (merge_mode < 6) {
      if (!keep_allele_order) {
	for (ukk = 0; ukk < ujj; ukk++) {
	  uljj = ((unsigned long)ukk) * tot_indiv_ct4;
	  umm = popcount_chars(pcptr, uljj, uljj + tot_indiv_ct4);
	  if (umm < tot_indiv_ct) {
	    ulkk = (uii * markers_per_pass) + ukk;
	    reversed[ulkk / BITCT] |= (1LU << (ulkk % BITCT));
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
      if (fwrite_checked(writebuf, ((unsigned long)ujj) * tot_indiv_ct4, outfile)) {
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
  if (wkspace_alloc_ui_checked(&map_reverse, dedup_marker_ct * sizeof(int))) {
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
    if (fclose_null(&outfile)) {
      goto merge_datasets_ret_WRITE_FAIL;
    }
    memcpy(outname_end, ".bed", 5);
    sprintf(logbuf, "Merged fileset written to %s + .bim + .fam.\n", outname);
    logprintb();
  } else {
    // undo the "not"
    diff_not_both_genotyped = diff_total_overlap - diff_not_both_genotyped;
    sprintf(logbuf, "%llu overlapping SNPs, %llu genotyped in both filesets.\n%llu concordant, for a concordance rate of %g.\n", diff_total_overlap, diff_not_both_genotyped, diff_not_both_genotyped - diff_discordant, 1.0 - (((double)diff_discordant) / ((double)diff_not_both_genotyped)));
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
