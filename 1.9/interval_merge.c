#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>

#define MAXLINELEN 256

char tbuf[MAXLINELEN];
char gene_name_buf[MAXLINELEN];
char chrom_name_buf[MAXLINELEN];

uint32_t scan_uint_capped(char* ss, uint32_t* valp, uint32_t cap_div_10, uint32_t cap_mod_10) {
  // Reads an integer in [0, cap].  Assumes first character is nonspace. 
  uint32_t val = (uint32_t)((unsigned char)*ss) - 48;
  uint32_t cur_digit;
  if (val < 10) {
    while (1) {
    scan_uint_capped_main_loop:
      cur_digit = (uint32_t)((unsigned char)(*(++ss))) - 48;
      if (cur_digit >= 10) {
	*valp = val;
	return 0;
      }
      if ((val >= cap_div_10) && ((val > cap_div_10) || (cur_digit > cap_mod_10))) {
	return 1;
      }
      val = val * 10 + cur_digit;
    }
  }
  // '-' has ascii code 45, so unsigned 45 - 48 = 0xfffffffdU
  ss++;
  if (val != 0xfffffffdU) {
    if (val == 0xfffffffbU) {
      val = (uint32_t)((unsigned char)(*ss)) - 48;
      if (val < 10) {
	goto scan_uint_capped_main_loop;
      }
    }
    return 1;
  }
  // accept "-0", "-00", etc.
  if (*ss != '0') {
    return 1;
  }
  while (*(++ss) == '0');
  *valp = 0;
  return ((uint32_t)((unsigned char)(*ss)) - 48) < 10;
}

static inline uint32_t scan_uint_defcap(char* ss, uint32_t* valp) {
  return scan_uint_capped(ss, valp, 0x7ffffffe / 10, 0x7ffffffe % 10);
}

static inline int32_t is_eoln_kns(char cc) {
  return ((unsigned char)cc) < 32;
}

static inline char* skip_initial_spaces(char* sptr) {
  while ((*sptr == ' ') || (*sptr == '\t')) {
    sptr++;
  }
  return sptr;
}

static inline int32_t is_space_or_eoln(char cc) {
  return ((unsigned char)cc) <= 32;
}

char* next_token(char* sptr) {
  if (!sptr) {
    return NULL;
  }
  while ((*sptr != ' ') && (*sptr != '\t')) {
    if (!(*sptr)) {
      return NULL;
    }
    sptr++;
  }
  return skip_initial_spaces(sptr);
}

static inline char* token_endnn(char* sptr) {
  while (!is_space_or_eoln(*(++sptr)));
  return sptr;
}

static inline void memcpyx(char* target, const void* source, uint32_t ct, const char extra_char) {
  memcpy(target, source, ct);
  target[ct] = extra_char;
}

int32_t main(int32_t argc, char** argv) {
  uint32_t first_pos = 0;
  uint32_t last_pos = 0;
  uint32_t gene_len = 0;
  uint32_t chrom_len = 0;
  uint32_t par_match = 0;
  uint32_t par_first_pos = 0;
  uint32_t par_last_pos = 0;
  int32_t retval = 0;
  char* chrom_ptr;
  char* gene_ptr;
  char* sptr;
  uint32_t cur_first_pos;
  uint32_t cur_last_pos;
  uint32_t cur_gene_len;
  uint32_t cur_chrom_len;
  *chrom_name_buf = '\0';
  *gene_name_buf = '\0';
  tbuf[MAXLINELEN - 1] = ' ';
  while (fgets(tbuf, MAXLINELEN, stdin)) {
    if (!tbuf[MAXLINELEN - 1]) {
      goto main_ret_LONG_LINE;
    }
    gene_ptr = skip_initial_spaces(tbuf);
    if (is_eoln_kns(*gene_ptr)) {
      continue;
    }
    chrom_ptr = next_token(gene_ptr);
    sptr = next_token(chrom_ptr);
    if ((!sptr) || scan_uint_defcap(sptr, &cur_first_pos)) {
      goto main_ret_INVALID_TOKEN;
    }
    sptr = next_token(sptr);
    if ((!sptr) || scan_uint_defcap(sptr, &cur_last_pos)) {
      goto main_ret_INVALID_TOKEN;
    }
    sptr = token_endnn(gene_ptr);
    cur_gene_len = (uintptr_t)(sptr - gene_ptr);
    sptr = token_endnn(chrom_ptr);
    cur_chrom_len = (uintptr_t)(sptr - chrom_ptr);
    if ((gene_len == cur_gene_len) && (chrom_len == cur_chrom_len) && (!memcmp(gene_name_buf, gene_ptr, gene_len)) && (!memcmp(chrom_name_buf, chrom_ptr, chrom_len)) && (cur_first_pos <= last_pos + 1)) {
      if (cur_first_pos < first_pos) {
	first_pos = cur_first_pos;
      }
      if (cur_last_pos > last_pos) {
	last_pos = cur_last_pos;
      }
    } else {
      if (gene_len) {
        printf("%s %u %u %s\n", chrom_name_buf, first_pos, last_pos, gene_name_buf);
	if ((chrom_len == 1) && (*chrom_name_buf == 'X') && (*chrom_ptr == 'Y') && (gene_len == cur_gene_len) && (!memcmp(gene_name_buf, gene_ptr, gene_len))) {
	  par_first_pos = first_pos;
          par_last_pos = last_pos;
	  par_match = 1;
	} else if (par_match) {
	  printf("XY %u %u %s\n", par_first_pos, par_last_pos, gene_name_buf);
	  par_match = 0;
	}
      }
      memcpyx(gene_name_buf, gene_ptr, cur_gene_len, '\0');
      memcpyx(chrom_name_buf, chrom_ptr, cur_chrom_len, '\0');
      gene_len = cur_gene_len;
      chrom_len = cur_chrom_len;
      first_pos = cur_first_pos;
      last_pos = cur_last_pos;
    }
  }
  printf("%s %u %u %s\n", chrom_name_buf, first_pos, last_pos, gene_name_buf);
  if (par_match) {
    printf("XY %u %u %s\n", par_first_pos, par_last_pos, gene_name_buf);
  }
  while (0) {
  main_ret_LONG_LINE:
    fputs("Error: Excessively long line.\n", stderr);
    retval = 1;
    break;
  main_ret_INVALID_TOKEN:
    fputs("Error: Missing or invalid token.\n", stderr);
    retval = 2;
    break;
  }
  return 0;
}
