// This file is part of PLINK 1.90, copyright (C) 2005-2020 Shaun Purcell,
// Christopher Chang.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


#include "plink_common.h"

// #include "pigz.h"

// no leading \n since this is used in LOGPRINTFWW expressions
const char g_errstr_fopen[] = "Error: Failed to open %s.\n";

const char g_cmdline_format_str[] = "\n  " PROG_NAME_STR " <input flag(s)...> [command flag(s)...] [other flag(s)...]\n  " PROG_NAME_STR " --help [flag name(s)...]\n\n";

char g_textbuf[TEXTBUF_SIZE];

// note that \xxx character constants are interpreted in octal.
// technically no need to represent 0-31, but 64 extra bytes of data is
// probably cheaper than the code to subtract 32 everywhere.
const char g_one_char_strs[] = "\0\0\1\0\2\0\3\0\4\0\5\0\6\0\7\0\10\0\11\0\12\0\13\0\14\0\15\0\16\0\17\0\20\0\21\0\22\0\23\0\24\0\25\0\26\0\27\0\30\0\31\0\32\0\33\0\34\0\35\0\36\0\37\0\40\0\41\0\42\0\43\0\44\0\45\0\46\0\47\0\50\0\51\0\52\0\53\0\54\0\55\0\56\0\57\0\60\0\61\0\62\0\63\0\64\0\65\0\66\0\67\0\70\0\71\0\72\0\73\0\74\0\75\0\76\0\77\0\100\0\101\0\102\0\103\0\104\0\105\0\106\0\107\0\110\0\111\0\112\0\113\0\114\0\115\0\116\0\117\0\120\0\121\0\122\0\123\0\124\0\125\0\126\0\127\0\130\0\131\0\132\0\133\0\134\0\135\0\136\0\137\0\140\0\141\0\142\0\143\0\144\0\145\0\146\0\147\0\150\0\151\0\152\0\153\0\154\0\155\0\156\0\157\0\160\0\161\0\162\0\163\0\164\0\165\0\166\0\167\0\170\0\171\0\172\0\173\0\174\0\175\0\176\0\177\0\200\0\201\0\202\0\203\0\204\0\205\0\206\0\207\0\210\0\211\0\212\0\213\0\214\0\215\0\216\0\217\0\220\0\221\0\222\0\223\0\224\0\225\0\226\0\227\0\230\0\231\0\232\0\233\0\234\0\235\0\236\0\237\0\240\0\241\0\242\0\243\0\244\0\245\0\246\0\247\0\250\0\251\0\252\0\253\0\254\0\255\0\256\0\257\0\260\0\261\0\262\0\263\0\264\0\265\0\266\0\267\0\270\0\271\0\272\0\273\0\274\0\275\0\276\0\277\0\300\0\301\0\302\0\303\0\304\0\305\0\306\0\307\0\310\0\311\0\312\0\313\0\314\0\315\0\316\0\317\0\320\0\321\0\322\0\323\0\324\0\325\0\326\0\327\0\330\0\331\0\332\0\333\0\334\0\335\0\336\0\337\0\340\0\341\0\342\0\343\0\344\0\345\0\346\0\347\0\350\0\351\0\352\0\353\0\354\0\355\0\356\0\357\0\360\0\361\0\362\0\363\0\364\0\365\0\366\0\367\0\370\0\371\0\372\0\373\0\374\0\375\0\376\0\377";
const char* g_missing_geno_ptr = &(g_one_char_strs[96]);
const char* g_output_missing_geno_ptr = &(g_one_char_strs[96]);

uintptr_t g_failed_alloc_attempt_size = 0;

sfmt_t g_sfmt;

FILE* g_logfile = nullptr;

char g_logbuf[MAXLINELEN * 2];

uint32_t g_debug_on = 0;
uint32_t g_log_failed = 0;
uint32_t g_thread_ct;

uint32_t aligned_malloc(uintptr_t size, uintptr_t** aligned_pp) {
#if defined __LP64__ && !defined __APPLE__
  // Avoid random segfaults on 64-bit machines which have 8-byte- instead of
  // 16-byte-aligned malloc().  (Slightly different code is needed if malloc()
  // does not even guarantee 8-byte alignment.)
  uintptr_t* malloc_ptr = (uintptr_t*)malloc(size + VEC_BYTES);
  if (!malloc_ptr) {
    g_failed_alloc_attempt_size = size + VEC_BYTES;
    return 1;
  }
  *aligned_pp = (uintptr_t*)((((uintptr_t)malloc_ptr) + VEC_BYTES) & (~(VEC_BYTES_M1 * ONELU)));
  (*aligned_pp)[-1] = (uintptr_t)malloc_ptr;
#else
  // no SSE2 concerns here
  *aligned_pp = (uintptr_t*)malloc(size);
  if (!(*aligned_pp)) {
    g_failed_alloc_attempt_size = size;
    return 1;
  }
#endif
  return 0;
}

void aligned_free(uintptr_t* aligned_pp) {
#if defined __LP64__ && !defined __APPLE__
  free((uintptr_t*)(aligned_pp[-1]));
#else
  free(aligned_pp);
#endif
}

uint32_t push_ll_str(const char* ss, Ll_str** ll_stack_ptr) {
  uintptr_t str_bytes = strlen(ss) + 1;
  Ll_str* new_ll_str = (Ll_str*)malloc(sizeof(Ll_str) + str_bytes);
  if (!new_ll_str) {
    g_failed_alloc_attempt_size = sizeof(Ll_str) + str_bytes;
    return 1;
  }
  new_ll_str->next = *ll_stack_ptr;
  memcpy(new_ll_str->ss, ss, str_bytes);
  *ll_stack_ptr = new_ll_str;
  return 0;
}

void logstr(const char* ss) {
  if (!g_debug_on) {
    fputs(ss, g_logfile);
    if (ferror(g_logfile)) {
      putc_unlocked('\n', stdout);
      fflush(stdout);
      fprintf(stderr, "Warning: Logging failure on:\n%s\nFurther logging will not be attempted in this run.\n", ss);
      g_log_failed = 1;
    }
  } else {
    if (g_log_failed) {
      fflush(stdout);
      fputs(ss, stderr);
    } else {
      fputs(ss, g_logfile);
      if (ferror(g_logfile)) {
	putc_unlocked('\n', stdout);
	fflush(stdout);
        fprintf(stderr, "Error: Debug logging failure.  Dumping to stderr:\n%s", ss);
	g_log_failed = 1;
      } else {
	fflush(g_logfile);
      }
    }
  }
}

void logprint(const char* ss) {
  logstr(ss);
  fputs(ss, stdout);
}

void logerrprint(const char* ss) {
  logstr(ss);
  fflush(stdout);
  fputs(ss, stderr);
}

void logprintb() {
  logstr(g_logbuf);
  fputs(g_logbuf, stdout);
}

void logerrprintb() {
  logstr(g_logbuf);
  fflush(stdout);
  fputs(g_logbuf, stderr);
}

void wordwrap(uint32_t suffix_len, char* ss) {
  // Input: A null-terminated string with no intermediate newlines.  If
  //        suffix_len is zero, there should be a terminating \n; otherwise,
  //        the last character should be a space.
  // Effect: Spaces are replaced with newlines in a manner that plays well with
  //         80 column terminal windows.  (Multi-space blocks are never
  //         collapsed.)
  char* token_start = ss;
  char* line_end = &(ss[79]);
  char* token_end;
  while (1) {
    while (*token_start == ' ') {
      token_start++;
    }
    if (token_start > line_end) {
      do {
	*line_end = '\n';
	line_end = &(line_end[80]);
      } while (token_start > line_end);
    }
    token_end = strchr(token_start, ' ');
    if (!token_end) {
      if (&(token_start[79]) == line_end) {
	return;
      }
      token_end = strchr(token_start, '\0');
      if (!suffix_len) {
	if (token_end <= &(line_end[1])) {
	  // okay if end-of-string is one past the end, because function
	  // assumes last character is \n in suffix_len == 0 case
	  assert(token_end[-1] == '\n');
	  return;
	}
      } else {
        if (&(token_end[suffix_len]) <= line_end) {
	  return;
	}
	// because of terminal space assumption, token_start actually points
	// to the end of the string
	assert(token_start[-1] == ' ');
      }
      token_start[-1] = '\n';
      return;
    }
    if (token_end > line_end) {
      if (&(token_start[79]) != line_end) {
	token_start[-1] = '\n';
        line_end = &(token_start[79]);
	if (token_end > line_end) {
	  // single really long token, can't do anything beyond putting it on
	  // its own line
          *token_end = '\n';
	  line_end = &(token_end[80]);
	}
      } else {
	// single really long token, *and* previous token was either
	// nonexistent or long
	*token_end = '\n';
	line_end = &(token_end[80]);
      }
    }
    token_start = &(token_end[1]);
  }
}

void wordwrapb(uint32_t suffix_len) {
  wordwrap(suffix_len, g_logbuf);
}

int32_t fopen_checked(const char* fname, const char* mode, FILE** target_ptr) {
  *target_ptr = fopen(fname, mode);
  if (!(*target_ptr)) {
    LOGERRPRINTFWW(g_errstr_fopen, fname);
    return -1;
  }
  return 0;
}

int32_t fwrite_checked(const void* buf, size_t len, FILE* outfile) {
  while (len > 0x7ffff000) {
    // OS X can't perform 2GB+ writes
    // typical disk block size is 4kb, so 0x7ffff000 is the largest sensible
    // write size
    fwrite(buf, 1, 0x7ffff000, outfile);
    buf = &(((unsigned char*)buf)[0x7ffff000]);
    len -= 0x7ffff000;
  }
  fwrite(buf, 1, len, outfile);
  return ferror(outfile);
}

int32_t gzopen_read_checked(const char* fname, gzFile* gzf_ptr) {
  *gzf_ptr = gzopen(fname, FOPEN_RB);
  if (!(*gzf_ptr)) {
    LOGERRPRINTFWW(g_errstr_fopen, fname);
    return RET_OPEN_FAIL;
  }
  if (gzbuffer(*gzf_ptr, 131072)) {
    return RET_NOMEM;
  }
  return 0;
}

// manually managed, very large stack
unsigned char* g_bigstack_base;
unsigned char* g_bigstack_end;

unsigned char* bigstack_alloc(uintptr_t size) {
  unsigned char* alloc_ptr;
  size = round_up_pow2(size, CACHELINE);
  if (bigstack_left() < size) {
    g_failed_alloc_attempt_size = size;
    return nullptr;
  }
  alloc_ptr = g_bigstack_base;
  g_bigstack_base += size;
  return alloc_ptr;
}

void bigstack_shrink_top(const void* rebase, uintptr_t new_size) {
  uintptr_t freed_bytes = ((uintptr_t)(g_bigstack_base - ((unsigned char*)rebase))) - round_up_pow2(new_size, CACHELINE);
  g_bigstack_base -= freed_bytes;
}

unsigned char* bigstack_end_alloc_presized(uintptr_t size) {
  assert(!(size & END_ALLOC_CHUNK_M1));
  uintptr_t cur_bigstack_left = bigstack_left();
  if (size > cur_bigstack_left) {
    g_failed_alloc_attempt_size = size;
    return nullptr;
  } else {
    g_bigstack_end -= size;
    return g_bigstack_end;
  }
}

uint32_t match_upper(const char* ss, const char* fixed_str) {
  char cc = *fixed_str++;
  do {
    if ((((unsigned char)(*ss++)) & 0xdf) != ((unsigned char)cc)) {
      return 0;
    }
    cc = *fixed_str++;
  } while (cc);
  return !(*ss);
}

uint32_t match_upper_counted(const char* ss, const char* fixed_str, uint32_t ct) {
  do {
    if ((((unsigned char)(*ss++)) & 0xdf) != ((unsigned char)(*fixed_str++))) {
      return 0;
    }
  } while (--ct);
  return 1;
}

#ifdef __LP64__
static inline uint32_t scan_uint_capped_finish(const char* ss, uint64_t cap, uint32_t* valp) {
  uint64_t val = *valp;
  while (1) {
    // a little bit of unrolling seems to help
    const uint64_t cur_digit = (uint64_t)((unsigned char)(*ss++)) - 48;
    if (cur_digit >= 10) {
      break;
    }
    // val = val * 10 + cur_digit;
    const uint64_t cur_digit2 = (uint64_t)((unsigned char)(*ss++)) - 48;
    if (cur_digit2 >= 10) {
      val = val * 10 + cur_digit;
      if (val > cap) {
	return 1;
      }
      break;
    }
    val = val * 100 + cur_digit * 10 + cur_digit2;
    if (val > cap) {
      return 1;
    }
  }
  *valp = val;
  return 0;
}

uint32_t scan_posint_capped(const char* ss, uint64_t cap, uint32_t* valp) {
  // '0' has ascii code 48
  *valp = (uint32_t)((unsigned char)(*ss++)) - 48;
  if (*valp >= 10) {
    // permit leading '+' (ascii 43), but not '++' or '+-'
    if (*valp != 0xfffffffbU) {
      return 1;
    }
    *valp = (uint32_t)((unsigned char)(*ss++)) - 48;
    if (*valp >= 10) {
      return 1;
    }
  }
  while (!(*valp)) {
    *valp = (uint32_t)((unsigned char)(*ss++)) - 48;
    if ((*valp) >= 10) {
      return 1;
    }
  }
  return scan_uint_capped_finish(ss, cap, valp);
}

uint32_t scan_uint_capped(const char* ss, uint64_t cap, uint32_t* valp) {
  // Reads an integer in [0, cap].  Assumes first character is nonspace.
  uint32_t val = (uint32_t)((unsigned char)(*ss++)) - 48;
  if (val >= 10) {
    if (val != 0xfffffffbU) {
      // '-' has ascii code 45, so unsigned 45 - 48 = 0xfffffffdU
      if ((val != 0xfffffffdU) || (*ss != '0')) {
	return 1;
      }
      // accept "-0", "-00", etc.
      while (*(++ss) == '0');
      *valp = 0;
      return ((uint32_t)((unsigned char)(*ss)) - 48) < 10;
    }
    // accept leading '+'
    val = (uint32_t)((unsigned char)(*ss++)) - 48;
    if (val >= 10) {
      return 1;
    }
  }
  *valp = val;
  return scan_uint_capped_finish(ss, cap, valp);
}

uint32_t scan_int_abs_bounded(const char* ss, uint64_t bound, int32_t* valp) {
  // Reads an integer in [-bound, bound].  Assumes first character is nonspace.
  *valp = (uint32_t)((unsigned char)(*ss++)) - 48;
  int32_t sign = 1;
  if (((uint32_t)*valp) >= 10) {
    if (*valp == -3) {
      sign = -1;
    } else if (*valp != -5) {
      return 1;
    }
    *valp = (uint32_t)((unsigned char)(*ss++)) - 48;
    if (((uint32_t)*valp) >= 10) {
      return 1;
    }
  }
  if (scan_uint_capped_finish(ss, bound, (uint32_t*)valp)) {
    return 1;
  }
  *valp *= sign;
  return 0;
}
#else // not __LP64__
uint32_t scan_posint_capped32(const char* ss, uint32_t cap_div_10, uint32_t cap_mod_10, uint32_t* valp) {
  // '0' has ascii code 48
  uint32_t val = (uint32_t)((unsigned char)(*ss++)) - 48;
  if (val >= 10) {
    if (val != 0xfffffffbU) {
      return 1;
    }
    val = (uint32_t)((unsigned char)(*ss++)) - 48;
    if (val >= 10) {
      return 1;
    }
  }
  while (!val) {
    val = (uint32_t)((unsigned char)(*ss++)) - 48;
    if (val >= 10) {
      return 1;
    }
  }
  while (1) {
    const uint32_t cur_digit = (uint32_t)((unsigned char)(*ss++)) - 48;
    if (cur_digit >= 10) {
      *valp = val;
      return 0;
    }
    // avoid integer overflow in middle of computation
    if ((val >= cap_div_10) && ((val > cap_div_10) || (cur_digit > cap_mod_10))) {
      return 1;
    }
    val = val * 10 + cur_digit;
  }
}

uint32_t scan_uint_capped32(const char* ss, uint32_t cap_div_10, uint32_t cap_mod_10, uint32_t* valp) {
  // Reads an integer in [0, cap].  Assumes first character is nonspace.
  uint32_t val = (uint32_t)((unsigned char)(*ss++)) - 48;
  if (val >= 10) {
    if (val != 0xfffffffbU) {
      if ((val != 0xfffffffdU) || (*ss != '0')) {
	return 1;
      }
      while (*(++ss) == '0');
      *valp = 0;
      return ((uint32_t)((unsigned char)(*ss)) - 48) < 10;
    }
    val = (uint32_t)((unsigned char)(*ss++)) - 48;
    if (val >= 10) {
      return 1;
    }
  }
  while (1) {
    const uint32_t cur_digit = (uint32_t)((unsigned char)(*ss++)) - 48;
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

uint32_t scan_int_abs_bounded32(const char* ss, uint32_t bound_div_10, uint32_t bound_mod_10, int32_t* valp) {
  // Reads an integer in [-bound, bound].  Assumes first character is nonspace.
  uint32_t val = (uint32_t)((unsigned char)(*ss++)) - 48;
  int32_t sign = 1;
  if (val >= 10) {
    if (val == 0xfffffffdU) {
      sign = -1;
    } else if (val != 0xfffffffbU) {
      return 1;
    }
    val = (uint32_t)((unsigned char)(*ss++)) - 48;
    if (val >= 10) {
      return 1;
    }
  }
  while (1) {
    const uint32_t cur_digit = (uint32_t)((unsigned char)(*ss++)) - 48;
    if (cur_digit >= 10) {
      *valp = sign * ((int32_t)val);
      return 0;
    }
    if ((val >= bound_div_10) && ((val > bound_div_10) || (cur_digit > bound_mod_10))) {
      return 1;
    }
    val = val * 10 + cur_digit;
  }
}
#endif

uint32_t scan_posintptr(const char* ss, uintptr_t* valp) {
  // Reads an integer in [1, 2^BITCT - 1].  Assumes first character is
  // nonspace.
  uintptr_t val = (uintptr_t)((unsigned char)(*ss++)) - 48;
  if (val >= 10) {
#ifdef __LP64__
    if (val != 0xfffffffffffffffbLLU) {
      return 1;
    }
#else
    if (val != 0xfffffffbU) {
      return 1;
    }
#endif
    val = (uintptr_t)((unsigned char)(*ss++)) - 48;
    if (val >= 10) {
      return 1;
    }
  }
  while (!val) {
    val = (uintptr_t)((unsigned char)(*ss++)) - 48;
    if (val >= 10) {
      return 1;
    }
  }
  // limit is 20 digits, we've already read one
#ifdef __LP64__
  const char* ss_limit = &(ss[20]);
#else
  const char* ss_limit = &(ss[10]);
#endif
  while (1) {
    const uintptr_t cur_digit = (uintptr_t)((unsigned char)(*ss++)) - 48;
    if (cur_digit >= 10) {
      *valp = val;
      return 0;
    }
    const uintptr_t cur_digit2 = (uintptr_t)((unsigned char)(*ss++)) - 48;
    if (ss == ss_limit) {
      if ((cur_digit2 < 10) || ((val >= (~ZEROLU) / 10) && ((val > (~ZEROLU) / 10) || (cur_digit > (~ZEROLU) % 10)))) {
	return 1;
      }
      *valp = val * 10 + cur_digit;
      return 0;
    }
    if (cur_digit2 >= 10) {
      *valp = val * 10 + cur_digit;
      return 0;
    }
    val = val * 100 + cur_digit * 10 + cur_digit2;
  }
}

/*
uint32_t scan_uintptr(char* ss, uintptr_t* valp) {
  // [0, 2^BITCT - 1].
  uintptr_t val = (uint32_t)((unsigned char)*ss) - 48;
  uintptr_t cur_digit;
  if (val < 10) {
    while (1) {
    scan_uintptr_main_loop:
      cur_digit = (uint32_t)((unsigned char)(*(++ss))) - 48;
      if (cur_digit >= 10) {
	*valp = val;
	return 0;
      }
      if ((val >= (~ZEROLU) / 10) && ((val > (~ZEROLU) / 10) || (cur_digit > (~ZEROLU) % 10))) {
	return 1;
      }
      val = val * 10 + cur_digit;
    }
  }
  ss++;
  if (val != 0xfffffffdU) {
    if (val == 0xfffffffbU) {
      val = (uint32_t)((unsigned char)(*ss)) - 48;
      if (val < 10) {
	goto scan_uintptr_main_loop;
      }
    }
    return 1;
  }
  if (*ss != '0') {
    return 1;
  }
  while (*(++ss) == '0');
  *valp = 0;
  return ((uint32_t)((unsigned char)(*ss)) - 48) < 10;
}
*/

uint32_t scan_posint_cappedx(const char* ss, uint64_t cap, uint32_t* valp) {
  double val;
  if (scan_doublex(ss, &val) || (val < 1.0) || (val > ((double)cap))) {
    return 1;
  }
  *valp = (uint32_t)val;
  return (val != ((double)(*valp)));
}

uint32_t scan_uint_cappedx(const char* ss, uint64_t cap, uint32_t* valp) {
  double val;
  if (scan_doublex(ss, &val) || (val < 0.0) || (val > ((double)cap))) {
    return 1;
  }
  *valp = (uint32_t)val;
  return (val != ((double)(*valp)));
}

uint32_t scan_int_abs_boundedx(const char* ss, uint64_t bound, int32_t* valp) {
  const double bound_d = (double)bound;
  double val;
  if (scan_doublex(ss, &val) || (val < -bound_d) || (val > bound_d)) {
    return 1;
  }
  *valp = (int32_t)val;
  return (val != ((double)(*valp)));
}

uint32_t scan_posintptrx(const char* ss, uintptr_t* valp) {
  double val;
  if (scan_doublex(ss, &val) || (val < 1.0) || (val > ((double)(~ZEROLU)))) {
    return 1;
  }
  *valp = (uintptr_t)val;
  return (val != ((double)(*valp)));
}


uint32_t scan_two_doubles(char* ss, double* __restrict val1p, double* __restrict val2p) {
  char* ss2;
  *val1p = strtod(ss, &ss2);
  if (ss == ss2) {
    return 1;
  }
  ss = skip_initial_spaces(ss2);
  *val2p = strtod(ss, &ss2);
  return (ss == ss2)? 1 : 0;
}

int32_t scan_token_ct_len(uintptr_t half_bufsize, FILE* infile, char* buf, uintptr_t* __restrict token_ct_ptr, uintptr_t* __restrict max_token_len_ptr) {
  // buf must be of size >= (2 * half_bufsize + 2)
  // max_token_len includes trailing null
  uintptr_t full_bufsize = half_bufsize * 2;
  uintptr_t curtoklen = 0;
  uintptr_t token_ct = *token_ct_ptr;
  uintptr_t max_token_len = *max_token_len_ptr;
  char* midbuf = &(buf[half_bufsize]);
  char* bufptr;
  char* bufptr2;
  char* buf_end;
  uintptr_t bufsize;
  while (1) {
    if (fread_checked(midbuf, half_bufsize, infile, &bufsize)) {
      return RET_READ_FAIL;
    }
    if (!bufsize) {
      if (curtoklen) {
        // corner case
        if (curtoklen >= max_token_len) {
	  max_token_len = curtoklen + 1;
	}
	token_ct++;
      }
      break;
    }
    buf_end = &(midbuf[bufsize]);
    *buf_end = ' ';
    buf_end[1] = '0';
    bufptr = &(buf[half_bufsize - curtoklen]);
    bufptr2 = midbuf;
    if (curtoklen) {
      goto scan_token_ct_len_tok_start;
    }
    while (1) {
      while (*bufptr <= ' ') {
	bufptr++;
      }
      if (bufptr >= buf_end) {
	curtoklen = 0;
	break;
      }
      bufptr2 = &(bufptr[1]);
    scan_token_ct_len_tok_start:
      while (*bufptr2 > ' ') {
	bufptr2++;
      }
      curtoklen = (uintptr_t)(bufptr2 - bufptr);
      if ((bufptr2 == buf_end) && (buf_end == &(buf[full_bufsize]))) {
	if (curtoklen >= half_bufsize) {
	  return RET_INVALID_FORMAT;
	}
	break;
      }
      if (curtoklen >= max_token_len) {
	if (curtoklen >= half_bufsize) {
	  return RET_INVALID_FORMAT;
	}
	max_token_len = curtoklen + 1;
      }
      token_ct++;
      bufptr = &(bufptr2[1]);
    }
  }
  if (!feof(infile)) {
    return RET_READ_FAIL;
  }
  *max_token_len_ptr = max_token_len;
  *token_ct_ptr = token_ct;
  return 0;
}

int32_t read_tokens(uintptr_t half_bufsize, uintptr_t token_ct, uintptr_t max_token_len, FILE* infile, char* __restrict buf, char* __restrict token_name_buf) {
  // buf must be of size >= (2 * half_bufsize + 2).
  // max_token_len includes trailing null
  uintptr_t full_bufsize = half_bufsize * 2;
  uintptr_t curtoklen = 0;
  uintptr_t token_idx = 0;
  char* midbuf = &(buf[half_bufsize]);
  char* bufptr = midbuf;
  char* bufptr2;
  char* bufptr3;
  char* buf_end;
  uintptr_t bufsize;
  while (1) {
    if (fread_checked(midbuf, half_bufsize, infile, &bufsize)) {
      return RET_READ_FAIL;
    }
    if (!bufsize) {
      if (curtoklen) {
        if (token_idx + 1 == token_ct) {
          memcpyx(&(token_name_buf[token_idx * max_token_len]), bufptr, curtoklen, '\0');
	  return 0;
        }
      }
      // something very strange has to happen to get here
      return RET_READ_FAIL;
    }
    buf_end = &(midbuf[bufsize]);
    *buf_end = ' ';
    buf_end[1] = '0';
    bufptr2 = midbuf;
    if (curtoklen) {
      goto read_tokens_tok_start;
    }
    while (1) {
      while (*bufptr <= ' ') {
	bufptr++;
      }
      if (bufptr >= buf_end) {
        curtoklen = 0;
	bufptr = midbuf;
	break;
      }
      bufptr2 = &(bufptr[1]);
    read_tokens_tok_start:
      while (*bufptr2 > ' ') {
	bufptr2++;
      }
      curtoklen = (uintptr_t)(bufptr2 - bufptr);
      if ((bufptr2 == buf_end) && (buf_end == &(buf[full_bufsize]))) {
	bufptr3 = &(buf[half_bufsize - curtoklen]);
        memcpy(bufptr3, bufptr, curtoklen);
        bufptr = bufptr3;
	break;
      }
      memcpyx(&(token_name_buf[token_idx * max_token_len]), bufptr, curtoklen, '\0');
      if (++token_idx == token_ct) {
	return 0;
      }
      bufptr = &(bufptr2[1]);
    }
  }
}

int32_t gzputs_w4(gzFile gz_outfile, const char* ss) {
  if (!ss[1]) {
    if (gzputs(gz_outfile, "   ") == -1) {
      return -1;
    }
    return gzputc(gz_outfile, ss[0]);
  }
  if (!ss[2]) {
    if (gzputs(gz_outfile, "  ") == -1) {
      return -1;
    }
  } else if (!ss[3]) {
    if (gzputc(gz_outfile, ' ') == -1) {
      return -1;
    }
  }
  return gzputs(gz_outfile, ss);
}

int32_t get_next_noncomment(FILE* fptr, char** lptr_ptr, uintptr_t* line_idx_ptr) {
  char* lptr;
  do {
    if (!fgets(g_textbuf, MAXLINELEN, fptr)) {
      return -1;
    }
    *line_idx_ptr += 1;
    lptr = skip_initial_spaces(g_textbuf);
  } while (is_eoln_or_comment_kns(*lptr));
  *lptr_ptr = lptr;
  return 0;
}

int32_t get_next_noncomment_excl(const uintptr_t* __restrict marker_exclude, FILE* fptr, char** lptr_ptr, uintptr_t* __restrict line_idx_ptr, uintptr_t* __restrict marker_uidx_ptr) {
  while (!get_next_noncomment(fptr, lptr_ptr, line_idx_ptr)) {
    if (!is_set_ul(marker_exclude, *marker_uidx_ptr)) {
      return 0;
    }
    *marker_uidx_ptr += 1;
  }
  return -1;
}

void get_top_two_ui(const uint32_t* __restrict uint_arr, uintptr_t uia_size, uintptr_t* __restrict top_idx_ptr, uintptr_t* __restrict second_idx_ptr) {
  assert(uia_size > 1);
  uintptr_t top_idx = (uint_arr[1] > uint_arr[0])? 1 : 0;
  uintptr_t second_idx = 1 ^ top_idx;
  uint32_t top_val = uint_arr[top_idx];
  uint32_t second_val = uint_arr[second_idx];
  uintptr_t cur_idx;
  uintptr_t cur_val;
  for (cur_idx = 2; cur_idx < uia_size; ++cur_idx) {
    cur_val = uint_arr[cur_idx];
    if (cur_val > second_val) {
      if (cur_val > top_val) {
	second_val = top_val;
	second_idx = top_idx;
	top_val = cur_val;
	top_idx = cur_idx;
      } else {
	second_val = cur_val;
	second_idx = cur_idx;
      }
    }
  }
  *top_idx_ptr = top_idx;
  *second_idx_ptr = second_idx;
}

uint32_t intlen(int32_t num) {
  int32_t retval = 1;
  uint32_t absnum;
  if (num < 0) {
    absnum = -num;
    retval++;
  } else {
    absnum = num;
  }
  while (absnum > 99) {
    // division by a constant is faster for unsigned ints
    absnum /= 100;
    retval += 2;
  }
  if (absnum > 9) {
    retval++;
  }
  return retval;
}

int32_t strcmp_se(const char* s_read, const char* s_const, uint32_t s_const_len) {
  return memcmp(s_read, s_const, s_const_len) || (!is_space_or_eoln(s_read[s_const_len]));
}

char* next_token(char* sptr) {
  if (!sptr) {
    return nullptr;
  }
  unsigned char ucc = *sptr;
  while (ucc > 32) {
    ucc = *(++sptr);
  }
  while ((ucc == ' ') || (ucc == '\t')) {
    ucc = *(++sptr);
  }
  return (ucc > 32)? sptr : nullptr;
}

char* next_token_mult(char* sptr, uint32_t ct) {
  assert(ct);
  if (!sptr) {
    return nullptr;
  }
  unsigned char ucc = *sptr;
  do {
    while (ucc > 32) {
      ucc = *(++sptr);
    }
    while ((ucc == ' ') || (ucc == '\t')) {
      ucc = *(++sptr);
    }
    if (ucc <= 32) {
      return nullptr;
    }
  } while (--ct);
  return sptr;
}

uint32_t count_tokens(const char* bufptr) {
  uint32_t token_ct = 0;
  while ((*bufptr == ' ') || (*bufptr == '\t')) {
    bufptr++;
  }
  while (!is_eoln_kns(*bufptr)) {
    token_ct++;
    while (!is_space_or_eoln(*(++bufptr)));
    while ((*bufptr == ' ') || (*bufptr == '\t')) {
      bufptr++;
    }
  }
  return token_ct;
}

uint32_t count_and_measure_multistr(const char* multistr, uintptr_t* max_slen_ptr) {
  // max_slen includes null terminator
  // assumes multistr is nonempty
  uint32_t ct = 0;
  uintptr_t max_slen = *max_slen_ptr;
  uintptr_t slen;
  do {
    slen = strlen(multistr) + 1;
    if (slen > max_slen) {
      max_slen = slen;
    }
    multistr = &(multistr[slen]);
    ct++;
  } while (*multistr);
  *max_slen_ptr = max_slen;
  return ct;
}

// number-to-string encoders

static const char digit2_table[200] = {
  '0', '0', '0', '1', '0', '2', '0', '3', '0', '4',
  '0', '5', '0', '6', '0', '7', '0', '8', '0', '9',
  '1', '0', '1', '1', '1', '2', '1', '3', '1', '4',
  '1', '5', '1', '6', '1', '7', '1', '8', '1', '9',
  '2', '0', '2', '1', '2', '2', '2', '3', '2', '4',
  '2', '5', '2', '6', '2', '7', '2', '8', '2', '9',
  '3', '0', '3', '1', '3', '2', '3', '3', '3', '4',
  '3', '5', '3', '6', '3', '7', '3', '8', '3', '9',
  '4', '0', '4', '1', '4', '2', '4', '3', '4', '4',
  '4', '5', '4', '6', '4', '7', '4', '8', '4', '9',
  '5', '0', '5', '1', '5', '2', '5', '3', '5', '4',
  '5', '5', '5', '6', '5', '7', '5', '8', '5', '9',
  '6', '0', '6', '1', '6', '2', '6', '3', '6', '4',
  '6', '5', '6', '6', '6', '7', '6', '8', '6', '9',
  '7', '0', '7', '1', '7', '2', '7', '3', '7', '4',
  '7', '5', '7', '6', '7', '7', '7', '8', '7', '9',
  '8', '0', '8', '1', '8', '2', '8', '3', '8', '4',
  '8', '5', '8', '6', '8', '7', '8', '8', '8', '9',
  '9', '0', '9', '1', '9', '2', '9', '3', '9', '4',
  '9', '5', '9', '6', '9', '7', '9', '8', '9', '9'};

char* uint32toa(uint32_t uii, char* start) {
  // Memory-efficient fast integer writer.  (You can do a bit better sometimes
  // by using a larger lookup table, but on average I doubt that pays off.)
  // Returns a pointer to the end of the integer (not null-terminated).
  uint32_t quotient;
  if (uii < 1000) {
    if (uii < 10) {
      *start++ = '0' + uii;
      return start;
    }
    if (uii < 100) {
      goto uint32toa_2;
    }
    quotient = uii / 100;
    *start++ = '0' + quotient;
  } else {
    if (uii < 10000000) {
      if (uii >= 100000) {
	if (uii < 1000000) {
	  goto uint32toa_6;
	}
	quotient = uii / 1000000;
	*start++ = '0' + quotient;
	goto uint32toa_6b;
      }
      if (uii < 10000) {
	goto uint32toa_4;
      }
      quotient = uii / 10000;
      *start++ = '0' + quotient;
    } else {
      if (uii >= 100000000) {
	quotient = uii / 100000000;
	if (uii >= 1000000000) {
	  start = memcpya(start, &(digit2_table[quotient * 2]), 2);
	} else {
	  *start++ = '0' + quotient;
	}
	uii -= 100000000 * quotient;
      }
      quotient = uii / 1000000;
      start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    uint32toa_6b:
      uii -= 1000000 * quotient;
    uint32toa_6:
      quotient = uii / 10000;
      start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    }
    uii -= 10000 * quotient;
  uint32toa_4:
    // could make a uitoa_z4() call here, but that's slightly slower
    quotient = uii / 100;
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
  }
  uii -= 100 * quotient;
 uint32toa_2:
  return memcpya(start, &(digit2_table[uii * 2]), 2);
}

char* int32toa(int32_t ii, char* start) {
  uint32_t uii = ii;
  if (ii < 0) {
    // -INT_MIN is undefined, but negating the unsigned int equivalent works
    *start++ = '-';
    uii = -uii;
  }
  return uint32toa(uii, start);
}

char* uitoa_z4(uint32_t uii, char* start) {
  uint32_t quotient = uii / 100;
  assert(quotient < 100);
  uii -= 100 * quotient;
  start = memcpya(start, &(digit2_table[quotient * 2]), 2);
  return memcpya(start, &(digit2_table[uii * 2]), 2);
}

char* uitoa_z6(uint32_t uii, char* start) {
  uint32_t quotient = uii / 10000;
  start = memcpya(start, &(digit2_table[quotient * 2]), 2);
  return uitoa_z4(uii - 10000 * quotient, start);
}

char* uitoa_z8(uint32_t uii, char* start) {
  uint32_t quotient = uii / 1000000;
  start = memcpya(start, &(digit2_table[quotient * 2]), 2);
  return uitoa_z6(uii - 1000000 * quotient, start);
}

char* int64toa(int64_t llii, char* start) {
  uint64_t ullii = llii;
  uint64_t top_digits;
  uint32_t bottom_eight;
  uint32_t middle_eight;
  if (llii < 0) {
    *start++ = '-';
    ullii = -ullii;
  }
  if (ullii <= 0xffffffffLLU) {
    return uint32toa((uint32_t)ullii, start);
  }
  top_digits = ullii / 100000000;
  bottom_eight = (uint32_t)(ullii - (top_digits * 100000000));
  if (top_digits <= 0xffffffffLLU) {
    start = uint32toa((uint32_t)top_digits, start);
    return uitoa_z8(bottom_eight, start);
  }
  ullii = top_digits / 100000000;
  middle_eight = (uint32_t)(top_digits - (ullii * 100000000));
  start = uint32toa((uint32_t)ullii, start);
  start = uitoa_z8(middle_eight, start);
  return uitoa_z8(bottom_eight, start);
}

char* uint32toa_w4(uint32_t uii, char* start) {
  uint32_t quotient;
  if (uii < 1000) {
    if (uii < 10) {
      // assumes little-endian
      *((uint32_t*)start) = 0x30202020 + (uii << 24);
      return &(start[4]);
    }
    if (uii < 100) {
      memset(start, 32, 2);
    } else {
      quotient = uii / 100;
      *start++ = ' ';
      *start++ = '0' + quotient;
      uii -= quotient * 100;
    }
    return memcpya(start, &(digit2_table[uii * 2]), 2);
  } else {
    // presumably the field width is 4 for a reason; don't bother optimizing
    // this
    return uint32toa(uii, start);
  }
}

char* uint32toa_w6(uint32_t uii, char* start) {
  uint32_t quotient;
  if (uii < 1000) {
    if (uii < 10) {
      start = memseta(start, 32, 5);
      *start++ = '0' + uii;
      return start;
    }
    if (uii < 100) {
      start = memseta(start, 32, 4);
      goto uint32toa_w6_2;
    }
    quotient = uii / 100;
    // the little-endian trick doesn't seem to help here.  possibly relevant
    // differences from uint32toa_w4() and _w8(): sequential dependence on
    // quotient, need to interpret pointer as a char* again
    start = memseta(start, 32, 3);
    *start++ = '0' + quotient;
  } else {
    if (uii < 10000000) {
      if (uii >= 100000) {
	if (uii < 1000000) {
	  goto uint32toa_w6_6;
	}
	quotient = uii / 1000000;
	*start++ = '0' + quotient;
	goto uint32toa_w6_6b;
      } else if (uii >= 10000) {
	*start++ = ' ';
	quotient = uii / 10000;
	*start++ = '0' + quotient;
      } else {
	start = memseta(start, 32, 2);
	goto uint32toa_w6_4;
      }
    } else {
      if (uii >= 100000000) {
	quotient = uii / 100000000;
	if (uii >= 1000000000) {
	  start = memcpya(start, &(digit2_table[quotient * 2]), 2);
	} else {
	  *start++ = '0' + quotient;
	}
	uii -= 100000000 * quotient;
      }
      quotient = uii / 1000000;
      start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    uint32toa_w6_6b:
      uii -= 1000000 * quotient;
    uint32toa_w6_6:
      quotient = uii / 10000;
      start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    }
    uii -= 10000 * quotient;
  uint32toa_w6_4:
    quotient = uii / 100;
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
  }
  uii -= 100 * quotient;
 uint32toa_w6_2:
  return memcpya(start, &(digit2_table[uii * 2]), 2);
}

char* uint32toa_w7(uint32_t uii, char* start) {
  uint32_t quotient;
  if (uii < 1000) {
    if (uii < 10) {
      start = memseta(start, 32, 6);
      *start++ = '0' + uii;
      return start;
    }
    if (uii < 100) {
      start = memseta(start, 32, 5);
      goto uint32toa_w7_2;
    }
    quotient = uii / 100;
    start = memseta(start, 32, 4);
    *start++ = '0' + quotient;
  } else {
    if (uii < 10000000) {
      if (uii >= 100000) {
	if (uii >= 1000000) {
	  quotient = uii / 1000000;
	  *start++ = '0' + quotient;
	  goto uint32toa_w7_6b;
	}
	*start++ = ' ';
	goto uint32toa_w7_6;
      } else if (uii >= 10000) {
	start = memseta(start, 32, 2);
	quotient = uii / 10000;
	*start++ = '0' + quotient;
      } else {
	start = memseta(start, 32, 3);
	goto uint32toa_w7_4;
      }
    } else {
      if (uii >= 100000000) {
	quotient = uii / 100000000;
	if (uii >= 1000000000) {
	  start = memcpya(start, &(digit2_table[quotient * 2]), 2);
	} else {
	  *start++ = '0' + quotient;
	}
	uii -= 100000000 * quotient;
      }
      quotient = uii / 1000000;
      start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    uint32toa_w7_6b:
      uii -= 1000000 * quotient;
    uint32toa_w7_6:
      quotient = uii / 10000;
      start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    }
    uii -= 10000 * quotient;
  uint32toa_w7_4:
    quotient = uii / 100;
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
  }
  uii -= 100 * quotient;
 uint32toa_w7_2:
  return memcpya(start, &(digit2_table[uii * 2]), 2);
}

char* uint32toa_w8(uint32_t uii, char* start) {
  uint32_t quotient;
  if (uii < 1000) {
    if (uii < 10) {
#ifdef __LP64__
      *((uintptr_t*)start) = 0x3020202020202020LLU + (((uintptr_t)uii) << 56);
      return &(start[8]);
#else
      start = memseta(start, 32, 7);
      *start++ = '0' + uii;
      return start;
#endif
    }
    if (uii < 100) {
      start = memseta(start, 32, 6);
      goto uint32toa_w8_2;
    }
    quotient = uii / 100;
    start = memseta(start, 32, 5);
    *start++ = '0' + quotient;
  } else {
    if (uii < 10000000) {
      if (uii >= 100000) {
	if (uii < 1000000) {
	  start = memseta(start, 32, 2);
	  goto uint32toa_w8_6;
	}
	quotient = uii / 1000000;
	*start = ' ';
	start[1] = '0' + quotient;
	start += 2;
	goto uint32toa_w8_6b;
      } else if (uii < 10000) {
	start = memseta(start, 32, 4);
	goto uint32toa_w8_4;
      }
      memset(start, 32, 3);
      quotient = uii / 10000;
      start[3] = '0' + quotient;
      start += 4;
    } else {
      if (uii >= 100000000) {
	quotient = uii / 100000000;
	if (uii >= 1000000000) {
	  start = memcpya(start, &(digit2_table[quotient * 2]), 2);
	} else {
	  *start++ = '0' + quotient;
	}
	uii -= 100000000 * quotient;
      }
      quotient = uii / 1000000;
      start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    uint32toa_w8_6b:
      uii -= 1000000 * quotient;
    uint32toa_w8_6:
      quotient = uii / 10000;
      start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    }
    uii -= 10000 * quotient;
  uint32toa_w8_4:
    quotient = uii / 100;
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
  }
  uii -= 100 * quotient;
 uint32toa_w8_2:
  return memcpya(start, &(digit2_table[uii * 2]), 2);
}

char* uint32toa_w10(uint32_t uii, char* start) {
  // if we decide to reduce code size and optimize only one field width, this
  // should be it
  uint32_t quotient;
  if (uii < 1000) {
    if (uii < 10) {
      start = memseta(start, 32, 9);
      *start++ = '0' + uii;
      return start;
    }
    if (uii < 100) {
      start = memseta(start, 32, 8);
      goto uint32toa_w10_2;
    }
    quotient = uii / 100;
    start = memseta(start, 32, 7);
    *start++ = '0' + quotient;
  } else {
    if (uii < 10000000) {
      if (uii >= 100000) {
	if (uii < 1000000) {
	  start = memseta(start, 32, 4);
	  goto uint32toa_w10_6;
	}
	quotient = uii / 1000000;
	memset(start, 32, 3);
	start[3] = '0' + quotient;
	start += 4;
	goto uint32toa_w10_6b;
      } else if (uii < 10000) {
	start = memseta(start, 32, 6);
	goto uint32toa_w10_4;
      }
      memset(start, 32, 5);
      quotient = uii / 10000;
      start[5] = '0' + quotient;
      start += 6;
    } else {
      if (uii >= 100000000) {
	quotient = uii / 100000000;
	if (uii >= 1000000000) {
	  memcpy(start, &(digit2_table[quotient * 2]), 2);
	} else {
	  *start = ' ';
	  start[1] = '0' + quotient;
	}
	uii -= 100000000 * quotient;
      } else {
	memset(start, 32, 2);
      }
      quotient = uii / 1000000;
      memcpy(&(start[2]), &(digit2_table[quotient * 2]), 2);
      start += 4;
    uint32toa_w10_6b:
      uii -= 1000000 * quotient;
    uint32toa_w10_6:
      quotient = uii / 10000;
      start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    }
    uii -= 10000 * quotient;
  uint32toa_w10_4:
    quotient = uii / 100;
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
  }
  uii -= 100 * quotient;
 uint32toa_w10_2:
  return memcpya(start, &(digit2_table[uii * 2]), 2);
}

static inline char* uitoa_trunc2(uint32_t uii, char* start) {
  // Given 0 < uii < 100, writes uii without *trailing* zeroes.  (I.e. this is
  // for floating-point encoder use.)
  memcpy(start, &(digit2_table[uii * 2]), 2);
  if (start[1] != '0') {
    return &(start[2]);
  }
  return &(start[1]);
}

static inline char* uitoa_trunc3(uint32_t uii, char* start) {
  *start++ = '0' + (uii / 100);
  uii %= 100;
  if (!uii) {
    return start;
  }
  memcpy(start, &(digit2_table[uii * 2]), 2);
  if (start[1] != '0') {
    return &(start[2]);
  }
  return &(start[1]);
}

static inline char* uitoa_trunc4(uint32_t uii, char* start) {
  uint32_t quotient = uii / 100;
  memcpy(start, &(digit2_table[quotient * 2]), 2);
  uii -= 100 * quotient;
  if (uii) {
    start += 2;
    memcpy(start, &(digit2_table[uii * 2]), 2);
  }
  if (start[1] != '0') {
    return &(start[2]);
  }
  return &(start[1]);
}

static inline char* uitoa_trunc6(uint32_t uii, char* start) {
  uint32_t quotient = uii / 10000;
  memcpy(start, &(digit2_table[quotient * 2]), 2);
  uii -= 10000 * quotient;
  if (uii) {
    quotient = uii / 100;
    start += 2;
    memcpy(start, &(digit2_table[quotient * 2]), 2);
    uii -= 100 * quotient;
    if (uii) {
      start += 2;
      memcpy(start, &(digit2_table[uii * 2]), 2);
    }
  }
  if (start[1] != '0') {
    return &(start[2]);
  }
  return &(start[1]);
}

static inline char* uitoa_trunc8(uint32_t uii, char* start) {
  uint32_t quotient = uii / 1000000;
  memcpy(start, &(digit2_table[quotient * 2]), 2);
  uii -= 1000000 * quotient;
  if (uii) {
    quotient = uii / 10000;
    start += 2;
    memcpy(start, &(digit2_table[quotient * 2]), 2);
    uii -= 10000 * quotient;
    if (uii) {
      quotient = uii / 100;
      start += 2;
      memcpy(start, &(digit2_table[quotient * 2]), 2);
      uii -= 100 * quotient;
      if (uii) {
	start += 2;
	memcpy(start, &(digit2_table[uii * 2]), 2);
      }
    }
  }
  if (start[1] != '0') {
    return &(start[2]);
  }
  return &(start[1]);
}

static inline char* qrtoa_1p1(uint32_t quotient, uint32_t remainder, char* start) {
  start[0] = '0' + quotient;
  if (!remainder) {
    return &(start[1]);
  }
  start[1] = '.';
  start[2] = '0' + remainder;
  return &(start[3]);
}

static inline char* qrtoa_1p2(uint32_t quotient, uint32_t remainder, char* start) {
  *start++ = '0' + quotient;
  if (!remainder) {
    return start;
  }
  *start++ = '.';
  memcpy(start, &(digit2_table[remainder * 2]), 2);
  if (start[1] != '0') {
    return &(start[2]);
  }
  return &(start[1]);
}

static inline char* qrtoa_1p3(uint32_t quotient, uint32_t remainder, char* start) {
  // quotient = (int32_t)dxx;
  // remainder = ((int32_t)(dxx * 1000)) - (quotient * 1000);
  *start++ = '0' + quotient;
  if (!remainder) {
    return start;
  }
  *start++ = '.';
  quotient = remainder / 10;
  memcpy(start, &(digit2_table[quotient * 2]), 2);
  remainder -= 10 * quotient;
  if (remainder) {
    start[2] = '0' + remainder;
    return &(start[3]);
  }
  if (start[1] != '0') {
    return &(start[2]);
  }
  return &(start[1]);
}

static inline char* qrtoa_1p5(uint32_t quotient, uint32_t remainder, char* start) {
  *start++ = '0' + quotient;
  if (!remainder) {
    return start;
  }
  *start++ = '.';
  quotient = remainder / 1000;
  memcpy(start, &(digit2_table[quotient * 2]), 2);
  remainder -= 1000 * quotient;
  if (remainder) {
    quotient = remainder / 10;
    start += 2;
    memcpy(start, &(digit2_table[quotient * 2]), 2);
    remainder -= 10 * quotient;
    if (remainder) {
      start[2] = '0' + remainder;
      return &(start[3]);
    }
  }
  if (start[1] != '0') {
    return &(start[2]);
  }
  return &(start[1]);
}

static inline char* qrtoa_1p7(uint32_t quotient, uint32_t remainder, char* start) {
  *start++ = '0' + quotient;
  if (!remainder) {
    return start;
  }
  *start++ = '.';
  quotient = remainder / 100000;
  memcpy(start, &(digit2_table[quotient * 2]), 2);
  remainder -= 100000 * quotient;
  if (remainder) {
    quotient = remainder / 1000;
    start += 2;
    memcpy(start, &(digit2_table[quotient * 2]), 2);
    remainder -= 1000 * quotient;
    if (remainder) {
      quotient = remainder / 10;
      start += 2;
      memcpy(start, &(digit2_table[quotient * 2]), 2);
      remainder -= 10 * quotient;
      if (remainder) {
	start[2] = '0' + remainder;
	return &(start[3]);
      }
    }
  }
  if (start[1] != '0') {
    return &(start[2]);
  }
  return &(start[1]);
}

// Okay, time to do banker's rounding when printing doubles.  14 digits of
// precision are used in judging equality to 0.5 (actual precision of doubles
// is 15-17 digits); the intention is to capture all directly loaded or exactly
// computed edge cases (so enough tolerance is needed to survive the internal
// multiplications by powers of 10, etc.), while rounding a negligible number
// of honest-to-god 0.4999999s up and 0.5000001s down.
// To avoid inadvertent printing of an extra digit, there's a deliberate gap
// between the 99.9994999...-type bounds and the largest numbers that would
// actually round down.
static const double banker_round5[] = {0.499995, 0.500005};
static const double banker_round6[] = {0.4999995, 0.5000005};
static const double banker_round7[] = {0.49999995, 0.50000005};
static const double banker_round8[] = {0.499999995, 0.500000005};
static const double banker_round9[] = {0.4999999995, 0.5000000005};
static const double banker_round10[] = {0.49999999995, 0.50000000005};
static const double banker_round11[] = {0.499999999995, 0.500000000005};
static const double banker_round12[] = {0.4999999999995, 0.5000000000005};

static inline uint32_t double_bround(double dxx, const double* banker_round) {
  uint32_t result = (int32_t)dxx;
  return result + (int32_t)((dxx - ((int32_t)result)) + banker_round[result & 1]);
}

// These are separate functions so the compiler can optimize the integer
// divisions.
static inline void double_bround1(double dxx, const double* banker_round, uint32_t* quotientp, uint32_t* remainderp) {
  dxx *= 10;
  uint32_t remainder = (int32_t)dxx;
  remainder += (int32_t)((dxx - ((int32_t)remainder)) + banker_round[remainder & 1]);
  *quotientp = remainder / 10;
  *remainderp = remainder - (*quotientp) * 10;
}

static inline void double_bround2(double dxx, const double* banker_round, uint32_t* quotientp, uint32_t* remainderp) {
  dxx *= 100;
  uint32_t remainder = (int32_t)dxx;
  remainder += (int32_t)((dxx - ((int32_t)remainder)) + banker_round[remainder & 1]);
  *quotientp = remainder / 100;
  *remainderp = remainder - (*quotientp) * 100;
}

static inline void double_bround3(double dxx, const double* banker_round, uint32_t* quotientp, uint32_t* remainderp) {
  dxx *= 1000;
  uint32_t remainder = (int32_t)dxx;
  remainder += (int32_t)((dxx - ((int32_t)remainder)) + banker_round[remainder & 1]);
  *quotientp = remainder / 1000;
  *remainderp = remainder - (*quotientp) * 1000;
}

static inline void double_bround4(double dxx, const double* banker_round, uint32_t* quotientp, uint32_t* remainderp) {
  dxx *= 10000;
  uint32_t remainder = (int32_t)dxx;
  remainder += (int32_t)((dxx - ((int32_t)remainder)) + banker_round[remainder & 1]);
  *quotientp = remainder / 10000;
  *remainderp = remainder - (*quotientp) * 10000;
}

static inline void double_bround5(double dxx, const double* banker_round, uint32_t* quotientp, uint32_t* remainderp) {
  dxx *= 100000;
  uint32_t remainder = (int32_t)dxx;
  remainder += (int32_t)((dxx - ((int32_t)remainder)) + banker_round[remainder & 1]);
  *quotientp = remainder / 100000;
  *remainderp = remainder - (*quotientp) * 100000;
}

static inline void double_bround6(double dxx, const double* banker_round, uint32_t* quotientp, uint32_t* remainderp) {
  dxx *= 1000000;
  uint32_t remainder = (int32_t)dxx;
  remainder += (int32_t)((dxx - ((int32_t)remainder)) + banker_round[remainder & 1]);
  *quotientp = remainder / 1000000;
  *remainderp = remainder - (*quotientp) * 1000000;
}

static inline void double_bround7(double dxx, const double* banker_round, uint32_t* quotientp, uint32_t* remainderp) {
  dxx *= 10000000;
  uint32_t remainder = (int32_t)dxx;
  remainder += (int32_t)((dxx - ((int32_t)remainder)) + banker_round[remainder & 1]);
  *quotientp = remainder / 10000000;
  *remainderp = remainder - (*quotientp) * 10000000;
}

char* dtoa_so6(double dxx, char* start) {
  // 6 sig fig number, 0.999995 <= dxx < 999999.5
  // 'so' = "significand only"
  // Just hardcoding all six cases, in the absence of a better approach...
  uint32_t uii;
  uint32_t quotient;
  uint32_t remainder;
  if (dxx < 99.999949999999) {
    if (dxx < 9.9999949999999) {
      double_bround5(dxx, banker_round8, &quotient, &remainder);
      return qrtoa_1p5(quotient, remainder, start);
    }
    double_bround4(dxx, banker_round8, &quotient, &remainder);
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    if (!remainder) {
      return start;
    }
    *start++ = '.';
    quotient = remainder / 100;
    memcpy(start, &(digit2_table[quotient * 2]), 2);
    remainder -= 100 * quotient;
    if (remainder) {
      start += 2;
    dtoa_so6_pretail:
      memcpy(start, &(digit2_table[remainder * 2]), 2);
    }
  dtoa_so6_tail:
    if (start[1] != '0') {
      return &(start[2]);
    }
    return &(start[1]);
  } else if (dxx < 9999.9949999999) {
    if (dxx < 999.99949999999) {
      double_bround3(dxx, banker_round8, &uii, &remainder);
      quotient = uii / 100;
      *start++ = '0' + quotient;
      quotient = uii - 100 * quotient;
      start = memcpya(start, &(digit2_table[quotient * 2]), 2);
      if (!remainder) {
	return start;
      }
      *start++ = '.';
      quotient = remainder / 10;
      memcpy(start, &(digit2_table[quotient * 2]), 2);
      remainder -= quotient * 10;
      if (!remainder) {
        goto dtoa_so6_tail;
      }
      start[2] = '0' + remainder;
      return &(start[3]);
    }
    double_bround2(dxx, banker_round8, &uii, &remainder);
    quotient = uii / 100;
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    quotient = uii - (100 * quotient);
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    if (!remainder) {
      return start;
    }
    *start++ = '.';
    goto dtoa_so6_pretail;
  } else if (dxx < 99999.949999999) {
    double_bround1(dxx, banker_round8, &uii, &remainder);
    quotient = uii / 10000;
    *start = '0' + quotient;
    uii -= 10000 * quotient;
    quotient = uii / 100;
    start = memcpya(&(start[1]), &(digit2_table[quotient * 2]), 2);
    uii = uii - 100 * quotient;
    start = memcpya(start, &(digit2_table[uii * 2]), 2);
    if (!remainder) {
      return start;
    }
    *start++ = '.';
    *start = '0' + remainder;
    return &(start[1]);
  } else {
    return uitoa_z6(double_bround(dxx, banker_round8), start);
  }
}

// Briefly had banker's rounding for floats, but then I realized that the only
// float-printing function calls are --make-grm related, they all request 6-7
// digits of precision, and at that point it's impossible to distinguish exact
// 0.5-matches in the remainder.  So we just have generic rounding functions
// here, with similar interfaces to the double-rounding functions to minimize
// the need for separate reasoning about this code.
static inline uint32_t float_round(float fxx) {
  return (uint32_t)((int32_t)(fxx + 0.5));
}

static inline void float_round1(float fxx, uint32_t* quotientp, uint32_t* remainderp) {
  uint32_t remainder = float_round(fxx * 10);
  *quotientp = remainder / 10;
  *remainderp = remainder - (*quotientp) * 10;
}

static inline void float_round2(float fxx, uint32_t* quotientp, uint32_t* remainderp) {
  uint32_t remainder = float_round(fxx * 100);
  *quotientp = remainder / 100;
  *remainderp = remainder - (*quotientp) * 100;
}

static inline void float_round3(float fxx, uint32_t* quotientp, uint32_t* remainderp) {
  uint32_t remainder = float_round(fxx * 1000);
  *quotientp = remainder / 1000;
  *remainderp = remainder - (*quotientp) * 1000;
}

static inline void float_round4(float fxx, uint32_t* quotientp, uint32_t* remainderp) {
  uint32_t remainder = float_round(fxx * 10000);
  *quotientp = remainder / 10000;
  *remainderp = remainder - (*quotientp) * 10000;
}

static inline void float_round5(float fxx, uint32_t* quotientp, uint32_t* remainderp) {
  uint32_t remainder = float_round(fxx * 100000);
  *quotientp = remainder / 100000;
  *remainderp = remainder - (*quotientp) * 100000;
}

static inline void float_round6(float fxx, uint32_t* quotientp, uint32_t* remainderp) {
  uint32_t remainder = float_round(fxx * 1000000);
  *quotientp = remainder / 1000000;
  *remainderp = remainder - (*quotientp) * 1000000;
}

char* ftoa_so6(float fxx, char* start) {
  uint32_t uii;
  uint32_t quotient;
  uint32_t remainder;
  // difference between consecutive floats near 10 can be as large as
  // 10 * 2^{-23}, which is just under 1.2e-6.  So, to avoid printing an extra
  // digit, we have to set this bound to be robust to an addition error of size
  // 6e-7.
  // (possible todo: just brute-force test this on all <2^32 possible floats
  // and look for a better threshold)
  if (fxx < 99.999944) {
    if (fxx < 9.9999944) {
      float_round5(fxx, &quotient, &remainder);
      return qrtoa_1p5(quotient, remainder, start);
    }
    float_round4(fxx, &quotient, &remainder);
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    if (!remainder) {
      return start;
    }
    *start++ = '.';
    quotient = remainder / 100;
    memcpy(start, &(digit2_table[quotient * 2]), 2);
    remainder -= 100 * quotient;
    if (remainder) {
      start += 2;
    ftoa_so6_pretail:
      memcpy(start, &(digit2_table[remainder * 2]), 2);
    }
  ftoa_so6_tail:
    if (start[1] != '0') {
      return &(start[2]);
    }
    return &(start[1]);
  } else if (fxx < 9999.9944) {
    if (fxx < 999.99944) {
      float_round3(fxx, &uii, &remainder);
      quotient = uii / 100;
      *start = '0' + quotient;
      quotient = uii - 100 * quotient;
      start = memcpya(&(start[1]), &(digit2_table[quotient * 2]), 2);
      if (!remainder) {
	return start;
      }
      *start++ = '.';
      quotient = remainder / 10;
      memcpy(start, &(digit2_table[quotient * 2]), 2);
      remainder -= quotient * 10;
      if (!remainder) {
        goto ftoa_so6_tail;
      }
      start[2] = '0' + remainder;
      return &(start[3]);
    }
    float_round2(fxx, &uii, &remainder);
    quotient = uii / 100;
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    quotient = uii - (100 * quotient);
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    if (!remainder) {
      return start;
    }
    *start++ = '.';
    goto ftoa_so6_pretail;
  } else if (fxx < 99999.944) {
    float_round1(fxx, &uii, &remainder);
    quotient = uii / 10000;
    *start = '0' + quotient;
    uii -= 10000 * quotient;
    quotient = uii / 100;
    start = memcpya(&(start[1]), &(digit2_table[quotient * 2]), 2);
    uii = uii - 100 * quotient;
    start = memcpya(start, &(digit2_table[uii * 2]), 2);
    if (!remainder) {
      return start;
    }
    *start = '.';
    start[1] = '0' + remainder;
    return &(start[2]);
  } else {
    return uitoa_z6(float_round(fxx), start);
  }
}

char* dtoa_so2(double dxx, char* start) {
  // 2 sig fig number, 0.95 <= dxx < 99.5
  uint32_t quotient;
  uint32_t remainder;
  if (dxx < 9.9499999999999) {
    double_bround1(dxx, banker_round12, &quotient, &remainder);
    return qrtoa_1p1(quotient, remainder, start);
  }
  return memcpya(start, &(digit2_table[(double_bround(dxx, banker_round12)) * 2]), 2);
}

char* dtoa_so3(double dxx, char* start) {
  // 3 sig fig number, 0.995 <= dxx < 999.5
  uint32_t quotient;
  uint32_t remainder;
  if (dxx < 99.949999999999) {
    if (dxx < 9.9949999999999) {
      double_bround2(dxx, banker_round11, &quotient, &remainder);
      return qrtoa_1p2(quotient, remainder, start);
    }
    double_bround1(dxx, banker_round11, &quotient, &remainder);
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    if (!remainder) {
      return start;
    }
    *start++ = '.';
  } else {
    quotient = double_bround(dxx, banker_round11);
    start = memcpya(start, &(digit2_table[(quotient / 10) * 2]), 2);
    remainder = quotient % 10;
  }
  *start = '0' + remainder;
  return &(start[1]);
}

char* dtoa_so4(double dxx, char* start) {
  // 4 sig fig number, 0.9995 <= dxx < 9999.5
  uint32_t uii;
  uint32_t quotient;
  uint32_t remainder;
  if (dxx < 99.994999999999) {
    if (dxx < 9.9994999999999) {
      double_bround3(dxx, banker_round10, &quotient, &remainder);
      return qrtoa_1p3(quotient, remainder, start);
    }
    double_bround2(dxx, banker_round10, &quotient, &remainder);
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    if (!remainder) {
      return start;
    }
    *start++ = '.';
    memcpy(start, &(digit2_table[remainder * 2]), 2);
    if (start[1] != '0') {
      return &(start[2]);
    }
    return &(start[1]);
  } else if (dxx < 999.94999999999) {
    double_bround1(dxx, banker_round10, &uii, &remainder);
    quotient = uii / 100;
    *start = '0' + quotient;
    quotient = uii - 100 * quotient;
    start = memcpya(&(start[1]), &(digit2_table[quotient * 2]), 2);
    if (!remainder) {
      return start;
    }
    *start = '.';
    start[1] = '0' + remainder;
    return &(start[2]);
  } else {
    uitoa_z4(double_bround(dxx, banker_round10), start);
    return &(start[4]);
  }
}

char* dtoa_so8(double dxx, char* start) {
  // 8 sig fig number, 0.99999995 <= dxx < 99999999.5
  uint32_t uii;
  uint32_t quotient;
  uint32_t remainder;
  if (dxx < 99.999999499999) {
    if (dxx < 9.9999999499999) {
      double_bround7(dxx, banker_round6, &quotient, &remainder);
      return qrtoa_1p7(quotient, remainder, start);
    }
    double_bround6(dxx, banker_round6, &quotient, &remainder);
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    if (!remainder) {
      return start;
    }
    *start++ = '.';
    quotient = remainder / 10000;
    memcpy(start, &(digit2_table[quotient * 2]), 2);
    remainder -= 10000 * quotient;
    if (remainder) {
      start += 2;
    dtoa_so8_pretail4:
      quotient = remainder / 100;
      memcpy(start, &(digit2_table[quotient * 2]), 2);
      remainder -= 100 * quotient;
      if (remainder) {
	start += 2;
      dtoa_so8_pretail2:
        memcpy(start, &(digit2_table[remainder * 2]), 2);
      }
    }
  dtoa_so8_tail:
    if (start[1] != '0') {
      return &(start[2]);
    }
    return &(start[1]);
  } else if (dxx < 9999.9999499999) {
    if (dxx < 999.99999499999) {
      double_bround5(dxx, banker_round6, &uii, &remainder);
      quotient = uii / 100;
      *start++ = '0' + quotient;
      quotient = uii - 100 * quotient;
      start = memcpya(start, &(digit2_table[quotient * 2]), 2);
      if (!remainder) {
	return start;
      }
      *start++ = '.';
      quotient = remainder / 1000;
      memcpy(start, &(digit2_table[quotient * 2]), 2);
      remainder -= quotient * 1000;
      if (!remainder) {
        goto dtoa_so8_tail;
      }
      start += 2;
    dtoa_so8_pretail3:
      quotient = remainder / 10;
      memcpy(start, &(digit2_table[quotient * 2]), 2);
      remainder -= quotient * 10;
      if (!remainder) {
	goto dtoa_so8_tail;
      }
      start[2] = '0' + remainder;
      return &(start[3]);
    }
    double_bround4(dxx, banker_round6, &uii, &remainder);
    quotient = uii / 100;
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    quotient = uii - (100 * quotient);
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    if (!remainder) {
      return start;
    }
    *start++ = '.';
    goto dtoa_so8_pretail4;
  } else if (dxx < 999999.99499999) {
    if (dxx < 99999.999499999) {
      double_bround3(dxx, banker_round6, &uii, &remainder);
      quotient = uii / 10000;
      *start = '0' + quotient;
      uii -= 10000 * quotient;
      quotient = uii / 100;
      start = memcpya(&(start[1]), &(digit2_table[quotient * 2]), 2);
      uii -= 100 * quotient;
      start = memcpya(start, &(digit2_table[uii * 2]), 2);
      if (!remainder) {
	return start;
      }
      *start++ = '.';
      goto dtoa_so8_pretail3;
    }
    double_bround2(dxx, banker_round6, &uii, &remainder);
    quotient = uii / 10000;
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    uii -= 10000 * quotient;
    quotient = uii / 100;
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    uii -= 100 * quotient;
    start = memcpya(start, &(digit2_table[uii * 2]), 2);
    if (!remainder) {
      return start;
    }
    *start++ = '.';
    goto dtoa_so8_pretail2;
  } else if (dxx < 9999999.9499999) {
    double_bround1(dxx, banker_round6, &uii, &remainder);
    quotient = uii / 1000000;
    *start = '0' + quotient;
    uii -= 1000000 * quotient;
    quotient = uii / 10000;
    start = memcpya(&(start[1]), &(digit2_table[quotient * 2]), 2);
    uii -= 10000 * quotient;
    quotient = uii / 100;
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    uii -= 100 * quotient;
    start = memcpya(start, &(digit2_table[uii * 2]), 2);
    if (!remainder) {
      return start;
    }
    *start = '.';
    start[1] = '0' + remainder;
    return &(start[2]);
  } else {
    return uitoa_z8(double_bround(dxx, banker_round6), start);
  }
}

char* dtoa_e(double dxx, char* start) {
  uint32_t xp10 = 0;
  uint32_t quotient;
  uint32_t remainder;
  char sign;
  if (dxx != dxx) {
    // do this first to avoid generating exception
    return memcpyl3a(start, "nan");
  } else if (dxx < 0) {
    *start++ = '-';
    dxx = -dxx;
  }
  if (dxx >= 9.9999994999999e-1) {
    if (dxx >= 9.9999994999999e7) {
      if (dxx >= 9.9999994999999e127) {
	if (dxx == INFINITY) {
	  return memcpyl3a(start, "inf");
	} else if (dxx >= 9.9999994999999e255) {
	  dxx *= 1.0e-256;
	  xp10 |= 256;
	} else {
	  dxx *= 1.0e-128;
	  xp10 |= 128;
	}
      }
      if (dxx >= 9.9999994999999e63) {
	dxx *= 1.0e-64;
	xp10 |= 64;
      }
      if (dxx >= 9.9999994999999e31) {
	dxx *= 1.0e-32;
	xp10 |= 32;
      }
      if (dxx >= 9.9999994999999e15) {
	dxx *= 1.0e-16;
	xp10 |= 16;
      }
      if (dxx >= 9.9999994999999e7) {
	dxx *= 1.0e-8;
	xp10 |= 8;
      }
    }
    if (dxx >= 9.9999994999999e3) {
      dxx *= 1.0e-4;
      xp10 |= 4;
    }
    if (dxx >= 9.9999994999999e1) {
      dxx *= 1.0e-2;
      xp10 |= 2;
    }
    if (dxx >= 9.9999994999999) {
      dxx *= 1.0e-1;
      xp10++;
    }
    sign = '+';
  } else {
    if (dxx < 9.9999994999999e-8) {
      // general case
      if (dxx < 9.9999994999999e-128) {
	if (dxx == 0.0) {
	  return memcpya(start, "0.000000e+00", 12);
	}
	if (dxx < 9.9999994999999e-256) {
	  dxx *= 1.0e256;
	  xp10 |= 256;
	} else {
	  dxx *= 1.0e128;
	  xp10 |= 128;
	}
      }
      if (dxx < 9.9999994999999e-64) {
	dxx *= 1.0e64;
	xp10 |= 64;
      }
      if (dxx < 9.9999994999999e-32) {
	dxx *= 1.0e32;
	xp10 |= 32;
      }
      if (dxx < 9.9999994999999e-16) {
	dxx *= 1.0e16;
	xp10 |= 16;
      }
      if (dxx < 9.9999994999999e-8) {
	dxx *= 100000000;
	xp10 |= 8;
      }
    }
    if (dxx < 9.999994999999e-4) {
      dxx *= 10000;
      xp10 |= 4;
    }
    if (dxx < 9.9999994999999e-2) {
      dxx *= 100;
      xp10 |= 2;
    }
    if (dxx < 9.9999994999999e-1) {
      dxx *= 10;
      xp10++;
    }
    sign = '-';
  }
  double_bround6(dxx, banker_round7, &quotient, &remainder);
  *start++ = '0' + quotient;
  *start++ = '.';
  start = uitoa_z6(remainder, start);
  *start++ = 'e';
  *start++ = sign;
  if (xp10 >= 100) {
    quotient = xp10 / 100;
    *start++ = '0' + quotient;
    xp10 -= quotient * 100;
  }
  return memcpya(start, &(digit2_table[xp10 * 2]), 2);
}

char* ftoa_e(float fxx, char* start) {
  uint32_t xp10 = 0;
  uint32_t quotient;
  uint32_t remainder;
  char sign;
  if (fxx != fxx) {
    // do this first to avoid generating exception
    return memcpyl3a(start, "nan");
  } else if (fxx < 0) {
    *start++ = '-';
    fxx = -fxx;
  }
  if (fxx >= 9.9999995e-1) {
    if (fxx >= 9.9999995e15) {
      if (fxx == INFINITY) {
	return memcpyl3a(start, "inf");
      } else if (fxx >= 9.9999995e31) {
	fxx *= 1.0e-32;
	xp10 |= 32;
      } else {
	fxx *= 1.0e-16;
	xp10 |= 16;
      }
    }
    if (fxx >= 9.9999995e7) {
      fxx *= 1.0e-8;
      xp10 |= 8;
    }
    if (fxx >= 9.9999995e3) {
      fxx *= 1.0e-4;
      xp10 |= 4;
    }
    if (fxx >= 9.9999995e1) {
      fxx *= 1.0e-2;
      xp10 |= 2;
    }
    if (fxx >= 9.9999995) {
      fxx *= 1.0e-1;
      xp10++;
    }
    sign = '+';
  } else {
    if (fxx < 9.9999995e-16) {
      if (fxx == 0.0) {
	return memcpya(start, "0.000000e+00", 12);
      } else if (fxx < 9.9999995e-32) {
	fxx *= 1.0e32;
	xp10 |= 32;
      } else {
	fxx *= 1.0e16;
	xp10 |= 16;
      }
    }
    if (fxx < 9.9999995e-8) {
      fxx *= 100000000;
      xp10 |= 8;
    }
    if (fxx < 9.9999995e-4) {
      fxx *= 10000;
      xp10 |= 4;
    }
    if (fxx < 9.9999995e-2) {
      fxx *= 100;
      xp10 |= 2;
    }
    if (fxx < 9.9999995e-1) {
      fxx *= 10;
      xp10++;
    }
    sign = '-';
  }
  float_round6(fxx, &quotient, &remainder);
  *start++ = '0' + quotient;
  *start++ = '.';
  start = uitoa_z6(remainder, start);
  *start++ = 'e';
  *start++ = sign;
  return memcpya(start, &(digit2_table[xp10 * 2]), 2);
}

char* dtoa_f_p2(double dxx, char* start) {
  const double* br_ptr;
  uint32_t quotient;
  uint32_t remainder;
  if (dxx != dxx) {
    return memcpyl3a(start, "nan");
  } else if (dxx < 9.9949999999999) {
    if (dxx < 0) {
      *start++ = '-';
      dxx = -dxx;
      if (dxx >= 9.9949999999999) {
        goto dtoa_f_p2_10;
      }
    }
    double_bround2(dxx, banker_round11, &quotient, &remainder);
    *start++ = '0' + quotient;
  dtoa_f_p2_dec:
    *start++ = '.';
    return memcpya(start, &(digit2_table[remainder * 2]), 2);
  }
 dtoa_f_p2_10:
  if (dxx < 9999999.9949999) {
    if (dxx < 999.99499999999) {
      if (dxx < 99.994999999999) {
	br_ptr = banker_round10;
      } else {
        br_ptr = banker_round9;
      }
    } else if (dxx < 99999.994999999) {
      if (dxx < 9999.9949999999) {
	br_ptr = banker_round8;
      } else {
	br_ptr = banker_round7;
      }
    } else if (dxx < 999999.99499999) {
      br_ptr = banker_round6;
    } else {
      br_ptr = banker_round5;
    }
    double_bround2(dxx, br_ptr, &quotient, &remainder);
    start = uint32toa(quotient, start);
    goto dtoa_f_p2_dec;
  }
  if (dxx == INFINITY) {
    return memcpyl3a(start, "inf");
  }
  // just punt larger numbers to glibc for now, this isn't a bottleneck
  start += sprintf(start, "%.2f", dxx);
  return start;
}

char* dtoa_f_p3(double dxx, char* start) {
  const double* br_ptr;
  uint32_t quotient;
  uint32_t remainder;
  if (dxx != dxx) {
    return memcpyl3a(start, "nan");
  } else if (dxx < 9.9994999999999) {
    if (dxx < 0) {
      *start++ = '-';
      dxx = -dxx;
      if (dxx >= 9.9994999999999) {
        goto dtoa_f_p3_10;
      }
    }
    double_bround3(dxx, banker_round10, &quotient, &remainder);
    *start++ = '0' + quotient;
  dtoa_f_p3_dec:
    *start++ = '.';
    quotient = remainder / 100;
    remainder -= 100 * quotient;
    *start++ = '0' + quotient;
    return memcpya(start, &(digit2_table[remainder * 2]), 2);
  }
 dtoa_f_p3_10:
  if (dxx < 999999.99949999) {
    if (dxx < 999.99949999999) {
      if (dxx < 99.999499999999) {
	br_ptr = banker_round9;
      } else {
        br_ptr = banker_round8;
      }
    } else if (dxx < 99999.999499999) {
      if (dxx < 9999.9994999999) {
	br_ptr = banker_round7;
      } else {
	br_ptr = banker_round6;
      }
    } else {
      br_ptr = banker_round5;
    }
    double_bround3(dxx, br_ptr, &quotient, &remainder);
    start = uint32toa(quotient, start);
    goto dtoa_f_p3_dec;
  }
  if (dxx == INFINITY) {
    return memcpyl3a(start, "inf");
  }
  start += sprintf(start, "%.3f", dxx);
  return start;
}

char* dtoa_f_w9p6(double dxx, char* start) {
  uint32_t quotient;
  uint32_t remainder;
  if (dxx != dxx) {
    return memcpya(start, "      nan", 9);
  } else if (dxx < 9.9999994999999) {
    if (dxx < 0) {
      *start++ = '-';
      dxx = -dxx;
      if (dxx >= 9.9999994999999) {
	goto dtoa_f_w9p6_10;
      }
    } else {
      *start++ = ' ';
    }
    double_bround6(dxx, banker_round7, &quotient, &remainder);
    *start++ = '0' + quotient;
  dtoa_f_w9p6_dec:
    *start++ = '.';
    return uitoa_z6(remainder, start);
  }
 dtoa_f_w9p6_10:
  if (dxx < 999.99999949999) {
    double_bround6(dxx, (dxx < 99.999999499999)? banker_round6 : banker_round5, &quotient, &remainder);
    start = uint32toa(quotient, start);
    goto dtoa_f_w9p6_dec;
  }
  if (dxx == INFINITY) {
    return memcpya(start, "      inf", 9);
  }
  start += sprintf(start, "%.6f", dxx);
  return start;
}

char* dtoa_f_w7p4(double dxx, char* start) {
  const double* br_ptr;
  uint32_t quotient;
  uint32_t remainder;
  if (dxx != dxx) {
    return memcpya(start, "    nan", 7);
  } else if (dxx < 9.9999499999999) {
    if (dxx < 0) {
      *start++ = '-';
      dxx = -dxx;
      if (dxx >= 9.9999499999999) {
	goto dtoa_f_w7p4_10;
      }
    } else {
      *start++ = ' ';
    }
    double_bround4(dxx, banker_round9, &quotient, &remainder);
    *start++ = '0' + quotient;
  dtoa_f_w7p4_dec:
    *start++ = '.';
    quotient = remainder / 100;
    remainder -= 100 * quotient;
    return memcpya(memcpya(start, &(digit2_table[quotient * 2]), 2), &(digit2_table[remainder * 2]), 2);
  }
 dtoa_f_w7p4_10:
  if (dxx < 99999.999949999) {
    if (dxx < 999.99994999999) {
      if (dxx < 99.999949999999) {
	br_ptr = banker_round8;
      } else {
	br_ptr = banker_round7;
      }
    } else if (dxx < 9999.9999499999) {
      br_ptr = banker_round6;
    } else {
      br_ptr = banker_round5;
    }
    double_bround4(dxx, br_ptr, &quotient, &remainder);
    start = uint32toa(quotient, start);
    goto dtoa_f_w7p4_dec;
  }
  if (dxx == INFINITY) {
    return memcpya(start, "    inf", 7);
  }
  start += sprintf(start, "%.4f", dxx);
  return start;
}

char* dtoa_f_w9p6_spaced(double dxx, char* start) {
  // Prettier fixed-width decimal: removes trailing zero(es) if and only if the
  // match appears to be exact.
  // Does not detect exact matches when abs(dxx) > 2^31 / 10^5.
  double dyy = dxx * 100000 + 0.00000005;
  start = dtoa_f_w9p6(dxx, start);
  if (dyy - ((double)((int32_t)dyy)) >= 0.0000001) {
    return start;
  }
  trailing_zeroes_to_spaces(start);
  return start;
}

char* dtoa_f_w9p6_clipped(double dxx, char* start) {
  // same conditions as _spaced()
  double dyy = dxx * 100000 + 0.00000005;
  start = dtoa_f_w9p6(dxx, start);
  if (dyy - ((double)((int32_t)dyy)) >= 0.0000001) {
    return start;
  }
  return clip_trailing_zeroes(start);
}

char* dtoa_g(double dxx, char* start) {
  uint32_t xp10 = 0;
  uint32_t quotient;
  uint32_t remainder;
  if (dxx != dxx) {
    return memcpyl3a(start, "nan");
  } else if (dxx < 0) {
    *start++ = '-';
    dxx = -dxx;
  }
  if (dxx < 9.9999949999999e-5) {
    // 6 sig fig exponential notation, small
    if (dxx < 9.9999949999999e-16) {
      if (dxx < 9.9999949999999e-128) {
	if (dxx == 0.0) {
	  *start = '0';
	  return &(start[1]);
	} else if (dxx < 9.9999949999999e-256) {
	  dxx *= 1.0e256;
	  xp10 |= 256;
	} else {
	  dxx *= 1.0e128;
	  xp10 |= 128;
	}
      }
      if (dxx < 9.9999949999999e-64) {
	dxx *= 1.0e64;
	xp10 |= 64;
      }
      if (dxx < 9.9999949999999e-32) {
	dxx *= 1.0e32;
	xp10 |= 32;
      }
      if (dxx < 9.9999949999999e-16) {
	dxx *= 1.0e16;
	xp10 |= 16;
      }
    }
    if (dxx < 9.9999949999999e-8) {
      dxx *= 100000000;
      xp10 |= 8;
    }
    if (dxx < 9.9999949999999e-4) {
      dxx *= 10000;
      xp10 |= 4;
    }
    if (dxx < 9.9999949999999e-2) {
      dxx *= 100;
      xp10 |= 2;
    }
    if (dxx < 9.9999949999999e-1) {
      dxx *= 10;
      xp10++;
    }
    double_bround5(dxx, banker_round8, &quotient, &remainder);
    start = memcpya(qrtoa_1p5(quotient, remainder, start), "e-", 2);
    if (xp10 >= 100) {
      quotient = xp10 / 100;
      *start++ = '0' + quotient;
      xp10 -= 100 * quotient;
    }
    return memcpya(start, &(digit2_table[xp10 * 2]), 2);
  } else if (dxx >= 999999.49999999) {
    // 6 sig fig exponential notation, large
    if (dxx >= 9.9999949999999e15) {
      if (dxx >= 9.9999949999999e127) {
	if (dxx == INFINITY) {
	  return memcpyl3a(start, "inf");
	} else if (dxx >= 9.9999949999999e255) {
	  dxx *= 1.0e-256;
	  xp10 |= 256;
	} else {
	  dxx *= 1.0e-128;
	  xp10 |= 128;
	}
      }
      if (dxx >= 9.9999949999999e63) {
	dxx *= 1.0e-64;
	xp10 |= 64;
      }
      if (dxx >= 9.9999949999999e31) {
	dxx *= 1.0e-32;
	xp10 |= 32;
      }
      if (dxx >= 9.9999949999999e15) {
	dxx *= 1.0e-16;
	xp10 |= 16;
      }
    }
    if (dxx >= 9.9999949999999e7) {
      dxx *= 1.0e-8;
      xp10 |= 8;
    }
    if (dxx >= 9.9999949999999e3) {
      dxx *= 1.0e-4;
      xp10 |= 4;
    }
    if (dxx >= 9.9999949999999e1) {
      dxx *= 1.0e-2;
      xp10 |= 2;
    }
    if (dxx >= 9.9999949999999e0) {
      dxx *= 1.0e-1;
      xp10++;
    }
    double_bround5(dxx, banker_round8, &quotient, &remainder);
    start = memcpya(qrtoa_1p5(quotient, remainder, start), "e+", 2);
    if (xp10 >= 100) {
      quotient = xp10 / 100;
      *start++ = '0' + quotient;
      xp10 -= 100 * quotient;
    }
    return memcpya(start, &(digit2_table[xp10 * 2]), 2);
  } else if (dxx >= 0.99999949999999) {
    return dtoa_so6(dxx, start);
  } else {
    // 6 sig fig decimal, no less than ~0.0001
    start = memcpya(start, "0.", 2);
    if (dxx < 9.9999949999999e-3) {
      dxx *= 100;
      start = memcpya(start, "00", 2);
    }
    if (dxx < 9.9999949999999e-2) {
      dxx *= 10;
      *start++ = '0';
    }
    return uitoa_trunc6(double_bround(dxx * 1000000, banker_round8), start);
  }
}

char* ftoa_g(float fxx, char* start) {
  uint32_t xp10 = 0;
  uint32_t quotient;
  uint32_t remainder;
  if (fxx != fxx) {
    return memcpyl3a(start, "nan");
  } else if (fxx < 0) {
    *start++ = '-';
    fxx = -fxx;
  }
  if (fxx < 9.9999944e-5) {
    if (fxx < 9.9999944e-16) {
      if (fxx == 0.0) {
	*start = '0';
	return &(start[1]);
      } else if (fxx < 9.9999944e-32) {
	fxx *= 1.0e32;
	xp10 |= 32;
      } else {
	fxx *= 1.0e16;
	xp10 |= 16;
      }
    }
    if (fxx < 9.9999944e-8) {
      fxx *= 100000000;
      xp10 |= 8;
    }
    if (fxx < 9.9999944e-4) {
      fxx *= 10000;
      xp10 |= 4;
    }
    if (fxx < 9.9999944e-2) {
      fxx *= 100;
      xp10 |= 2;
    }
    if (fxx < 9.9999944e-1) {
      fxx *= 10;
      xp10++;
    }
    float_round5(fxx, &quotient, &remainder);
    return memcpya(memcpya(qrtoa_1p5(quotient, remainder, start), "e-", 2), &(digit2_table[xp10 * 2]), 2);
  } else if (fxx >= 999999.44) {
    if (fxx >= 9.9999944e15) {
      if (fxx == INFINITY) {
	return memcpyl3a(start, "inf");
      } else if (fxx >= 9.9999944e31) {
	fxx *= 1.0e-32;
	xp10 |= 32;
      } else {
	fxx *= 1.0e-16;
	xp10 |= 16;
      }
    }
    if (fxx >= 9.9999944e7) {
      fxx *= 1.0e-8;
      xp10 |= 8;
    }
    if (fxx >= 9.9999944e3) {
      fxx *= 1.0e-4;
      xp10 |= 4;
    }
    if (fxx >= 9.9999944e1) {
      fxx *= 1.0e-2;
      xp10 |= 2;
    }
    if (fxx >= 9.9999944e0) {
      fxx *= 1.0e-1;
      xp10++;
    }
    float_round5(fxx, &quotient, &remainder);
    return memcpya(memcpya(qrtoa_1p5(quotient, remainder, start), "e+", 2), &(digit2_table[xp10 * 2]), 2);
  } else if (fxx >= 0.99999944) {
    return ftoa_so6(fxx, start);
  } else {
    // 6 sig fig decimal, no less than ~0.0001
    start = memcpya(start, "0.", 2);
    if (fxx < 9.9999944e-3) {
      fxx *= 100;
      start = memcpya(start, "00", 2);
    }
    if (fxx < 9.9999944e-2) {
      fxx *= 10;
      *start++ = '0';
    }
    return uitoa_trunc6(float_round(fxx * 1000000), start);
  }
}

char* dtoa_g_wxp2(double dxx, uint32_t min_width, char* start) {
  assert(min_width >= 5);
  uint32_t xp10 = 0;
  char wbuf[16];
  char* wpos = wbuf;
  uint32_t quotient;
  uint32_t remainder;
  if (dxx != dxx) {
    memcpy(memseta(start, 32, min_width - 4), " nan", 4);
    return &(start[min_width]);
  } else if (dxx < 0) {
    *wpos++ = '-';
    dxx = -dxx;
  }
  if (dxx < 9.9499999999999e-5) {
    // 2 sig fig exponential notation, small
    if (dxx < 9.9499999999999e-16) {
      if (dxx < 9.9499999999999e-128) {
	if (dxx == 0.0) {
          memset(start, 32, min_width - 1);
	  start[min_width - 1] = '0';
	  return &(start[min_width]);
        } else if (dxx < 9.9499999999999e-256) {
	  dxx *= 1.0e256;
	  xp10 |= 256;
	} else {
	  dxx *= 1.0e128;
	  xp10 |= 128;
	}
      }
      if (dxx < 9.9499999999999e-64) {
	dxx *= 1.0e64;
	xp10 |= 64;
      }
      if (dxx < 9.9499999999999e-32) {
	dxx *= 1.0e32;
	xp10 |= 32;
      }
      if (dxx < 9.9499999999999e-16) {
	dxx *= 1.0e16;
	xp10 |= 16;
      }
    }
    if (dxx < 9.9499999999999e-8) {
      dxx *= 100000000;
      xp10 |= 8;
    }
    if (dxx < 9.9499999999999e-4) {
      dxx *= 10000;
      xp10 |= 4;
    }
    if (dxx < 9.9499999999999e-2) {
      dxx *= 100;
      xp10 |= 2;
    }
    if (dxx < 9.9499999999999e-1) {
      dxx *= 10;
      xp10++;
    }
    double_bround1(dxx, banker_round12, &quotient, &remainder);
    wpos = qrtoa_1p1(quotient, remainder, wpos);
    remainder = wpos - wbuf;
    if (xp10 >= 100) {
      if (remainder < min_width - 5) {
	memcpy(memseta(start, 32, min_width - 5 - remainder), wbuf, remainder);
	start = &(start[min_width - 5]);
      } else {
	start = memcpya(start, wbuf, remainder);
      }
      quotient = xp10 / 100;
      start = memcpyax(start, "e-", 2, '0' + quotient);
      xp10 -= 100 * quotient;
    } else {
      if (remainder < min_width - 4) {
	memcpy(memseta(start, 32, min_width - 4 - remainder), wbuf, remainder);
	start = &(start[min_width - 4]);
      } else {
	start = memcpya(start, wbuf, remainder);
      }
      start = memcpya(start, "e-", 2);
    }
    return memcpya(start, &(digit2_table[xp10 * 2]), 2);
  } else if (dxx >= 99.499999999999) {
    // 2 sig fig exponential notation, large
    if (dxx >= 9.9499999999999e15) {
      if (dxx >= 9.9499999999999e127) {
	if (dxx == INFINITY) {
	  start = memseta(start, 32, min_width - 4);
	  if (wpos == wbuf) {
	    return memcpya(start, " inf", 4);
	  } else {
	    return memcpya(start, "-inf", 4);
	  }
	} else if (dxx >= 9.9499999999999e255) {
	  dxx *= 1.0e-256;
	  xp10 |= 256;
	} else {
	  dxx *= 1.0e-128;
	  xp10 |= 128;
	}
      }
      if (dxx >= 9.9499999999999e63) {
	dxx *= 1.0e-64;
	xp10 |= 64;
      }
      if (dxx >= 9.9499999999999e31) {
	dxx *= 1.0e-32;
	xp10 |= 32;
      }
      if (dxx >= 9.9499999999999e15) {
	dxx *= 1.0e-16;
	xp10 |= 16;
      }
    }
    if (dxx >= 9.9499999999999e7) {
      dxx *= 1.0e-8;
      xp10 |= 8;
    }
    if (dxx >= 9.9499999999999e3) {
      dxx *= 1.0e-4;
      xp10 |= 4;
    }
    if (dxx >= 9.9499999999999e1) {
      dxx *= 1.0e-2;
      xp10 |= 2;
    }
    if (dxx >= 9.9499999999999e0) {
      dxx *= 1.0e-1;
      xp10++;
    }
    double_bround1(dxx, banker_round12, &quotient, &remainder);
    wpos = qrtoa_1p1(quotient, remainder, wpos);
    remainder = wpos - wbuf;
    if (xp10 >= 100) {
      if (remainder < min_width - 5) {
	memcpy(memseta(start, 32, min_width - 5 - remainder), wbuf, remainder);
	start = &(start[min_width - 5]);
      } else {
	start = memcpya(start, wbuf, remainder);
      }
      quotient = xp10 / 100;
      start = memcpyax(start, "e+", 2, '0' + quotient);
      xp10 -= 100 * quotient;
    } else {
      if (remainder < min_width - 4) {
	memcpy(memseta(start, 32, min_width - 4 - remainder), wbuf, remainder);
	start = &(start[min_width - 4]);
      } else {
	start = memcpya(start, wbuf, remainder);
      }
      start = memcpya(start, "e+", 2);
    }
    return memcpya(start, &(digit2_table[xp10 * 2]), 2);
  } else {
    if (dxx >= 0.99499999999999) {
      wpos = dtoa_so2(dxx, wpos);
    } else {
      // 2 sig fig decimal, no less than ~0.0001
      wpos = memcpya(wpos, "0.", 2);
      if (dxx < 9.9499999999999e-3) {
	dxx *= 100;
	wpos = memcpya(wpos, "00", 2);
      }
      if (dxx < 9.9499999999999e-2) {
	dxx *= 10;
	*wpos++ = '0';
      }
      wpos = uitoa_trunc2(double_bround(dxx * 100, banker_round12), wpos);
    }
    remainder = wpos - wbuf;
    if (remainder < min_width) {
      memcpy(memseta(start, 32, min_width - remainder), wbuf, remainder);
      return &(start[min_width]);
    } else {
      return memcpya(start, wbuf, remainder);
    }
  }
}

char* dtoa_g_wxp3(double dxx, uint32_t min_width, char* start) {
  assert(min_width >= 5);
  uint32_t xp10 = 0;
  char wbuf[16];
  char* wpos = wbuf;
  uint32_t quotient;
  uint32_t remainder;
  if (dxx != dxx) {
    memcpy(memseta(start, 32, min_width - 4), " nan", 4);
    return &(start[min_width]);
  } else if (dxx < 0) {
    *wpos++ = '-';
    dxx = -dxx;
  }
  if (dxx < 9.9949999999999e-5) {
    // 3 sig fig exponential notation, small
    if (dxx < 9.9949999999999e-16) {
      if (dxx < 9.9949999999999e-128) {
	if (dxx == 0.0) {
          memset(start, 32, min_width - 1);
	  start[min_width - 1] = '0';
	  return &(start[min_width]);
        } else if (dxx < 9.9949999999999e-256) {
	  dxx *= 1.0e256;
	  xp10 |= 256;
	} else {
	  dxx *= 1.0e128;
	  xp10 |= 128;
	}
      }
      if (dxx < 9.9949999999999e-64) {
	dxx *= 1.0e64;
	xp10 |= 64;
      }
      if (dxx < 9.9949999999999e-32) {
	dxx *= 1.0e32;
	xp10 |= 32;
      }
      if (dxx < 9.9949999999999e-16) {
	dxx *= 1.0e16;
	xp10 |= 16;
      }
    }
    if (dxx < 9.9949999999999e-8) {
      dxx *= 100000000;
      xp10 |= 8;
    }
    if (dxx < 9.9949999999999e-4) {
      dxx *= 10000;
      xp10 |= 4;
    }
    if (dxx < 9.9949999999999e-2) {
      dxx *= 100;
      xp10 |= 2;
    }
    if (dxx < 9.9949999999999e-1) {
      dxx *= 10;
      xp10++;
    }
    double_bround2(dxx, banker_round11, &quotient, &remainder);
    wpos = qrtoa_1p2(quotient, remainder, wpos);
    remainder = wpos - wbuf;
    if (xp10 >= 100) {
      if (remainder < min_width - 5) {
	memcpy(memseta(start, 32, min_width - 5 - remainder), wbuf, remainder);
	start = &(start[min_width - 5]);
      } else {
	start = memcpya(start, wbuf, remainder);
      }
      quotient = xp10 / 100;
      start = memcpyax(start, "e-", 2, '0' + quotient);
      xp10 -= 100 * quotient;
    } else {
      if (remainder < min_width - 4) {
	memcpy(memseta(start, 32, min_width - 4 - remainder), wbuf, remainder);
	start = &(start[min_width - 4]);
      } else {
	start = memcpya(start, wbuf, remainder);
      }
      start = memcpya(start, "e-", 2);
    }
    return memcpya(start, &(digit2_table[xp10 * 2]), 2);
  } else if (dxx >= 999.49999999999) {
    // 3 sig fig exponential notation, large
    if (dxx >= 9.9949999999999e15) {
      if (dxx >= 9.9949999999999e127) {
	if (dxx == INFINITY) {
	  start = memseta(start, 32, min_width - 4);
	  if (wpos == wbuf) {
	    return memcpya(start, " inf", 4);
	  } else {
	    return memcpya(start, "-inf", 4);
	  }
	} else if (dxx >= 9.9949999999999e255) {
	  dxx *= 1.0e-256;
	  xp10 |= 256;
	} else {
	  dxx *= 1.0e-128;
	  xp10 |= 128;
	}
      }
      if (dxx >= 9.9949999999999e63) {
	dxx *= 1.0e-64;
	xp10 |= 64;
      }
      if (dxx >= 9.9949999999999e31) {
	dxx *= 1.0e-32;
	xp10 |= 32;
      }
      if (dxx >= 9.9949999999999e15) {
	dxx *= 1.0e-16;
	xp10 |= 16;
      }
    }
    if (dxx >= 9.9949999999999e7) {
      dxx *= 1.0e-8;
      xp10 |= 8;
    }
    if (dxx >= 9.9949999999999e3) {
      dxx *= 1.0e-4;
      xp10 |= 4;
    }
    if (dxx >= 9.9949999999999e1) {
      dxx *= 1.0e-2;
      xp10 |= 2;
    }
    if (dxx >= 9.9949999999999e0) {
      dxx *= 1.0e-1;
      xp10++;
    }
    double_bround2(dxx, banker_round11, &quotient, &remainder);
    wpos = qrtoa_1p2(quotient, remainder, wpos);
    remainder = wpos - wbuf;
    if (xp10 >= 100) {
      if (remainder < min_width - 5) {
	memcpy(memseta(start, 32, min_width - 5 - remainder), wbuf, remainder);
	start = &(start[min_width - 5]);
      } else {
	start = memcpya(start, wbuf, remainder);
      }
      quotient = xp10 / 100;
      start = memcpyax(start, "e+", 2, '0' + quotient);
      xp10 -= 100 * quotient;
    } else {
      if (remainder < min_width - 4) {
	memcpy(memseta(start, 32, min_width - 4 - remainder), wbuf, remainder);
	start = &(start[min_width - 4]);
      } else {
	start = memcpya(start, wbuf, remainder);
      }
      start = memcpya(start, "e+", 2);
    }
    return memcpya(start, &(digit2_table[xp10 * 2]), 2);
  } else {
    if (dxx >= 0.99949999999999) {
      wpos = dtoa_so3(dxx, wpos);
    } else {
      // 3 sig fig decimal, no less than ~0.001
      wpos = memcpya(wpos, "0.", 2);
      if (dxx < 9.9949999999999e-3) {
	dxx *= 100;
	wpos = memcpya(wpos, "00", 2);
      }
      if (dxx < 9.9949999999999e-2) {
	dxx *= 10;
	*wpos++ = '0';
      }
      wpos = uitoa_trunc3(double_bround(dxx * 1000, banker_round11), wpos);
    }
    remainder = wpos - wbuf;
    if (remainder < min_width) {
      memcpy(memseta(start, 32, min_width - remainder), wbuf, remainder);
      return &(start[min_width]);
    } else {
      return memcpya(start, wbuf, remainder);
    }
  }
}

char* dtoa_g_wxp4(double dxx, uint32_t min_width, char* start) {
  uint32_t xp10 = 0;
  char wbuf[16];
  char* wpos = wbuf;
  uint32_t quotient;
  uint32_t remainder;
  if (dxx != dxx) {
    if (min_width > 3) {
      start = memseta(start, 32, min_width - 3);
    }
    return memcpyl3a(start, "nan");
  } else if (dxx < 0) {
    *wpos++ = '-';
    dxx = -dxx;
  }
  if (dxx < 9.9994999999999e-5) {
    // 4 sig fig exponential notation, small
    if (dxx < 9.9994999999999e-16) {
      if (dxx < 9.9994999999999e-128) {
	if (dxx == 0.0) {
          memset(start, 32, min_width - 1);
	  start[min_width - 1] = '0';
	  return &(start[min_width]);
        } else if (dxx < 9.9994999999999e-256) {
	  dxx *= 1.0e256;
	  xp10 |= 256;
	} else {
	  dxx *= 1.0e128;
	  xp10 |= 128;
	}
      }
      if (dxx < 9.9994999999999e-64) {
	dxx *= 1.0e64;
	xp10 |= 64;
      }
      if (dxx < 9.9994999999999e-32) {
	dxx *= 1.0e32;
	xp10 |= 32;
      }
      if (dxx < 9.9994999999999e-16) {
	dxx *= 1.0e16;
	xp10 |= 16;
      }
    }
    if (dxx < 9.9994999999999e-8) {
      dxx *= 100000000;
      xp10 |= 8;
    }
    if (dxx < 9.9994999999999e-4) {
      dxx *= 10000;
      xp10 |= 4;
    }
    if (dxx < 9.9994999999999e-2) {
      dxx *= 100;
      xp10 |= 2;
    }
    if (dxx < 9.9994999999999e-1) {
      dxx *= 10;
      xp10++;
    }
    double_bround3(dxx, banker_round10, &quotient, &remainder);
    wpos = qrtoa_1p3(quotient, remainder, wpos);
    remainder = wpos - wbuf;
    if (xp10 >= 100) {
      if (remainder + 5 < min_width) {
	memcpy(memseta(start, 32, min_width - (remainder + 5)), wbuf, remainder);
	start = &(start[min_width - 5]);
      } else {
	start = memcpya(start, wbuf, remainder);
      }
      quotient = xp10 / 100;
      start = memcpyax(start, "e-", 2, '0' + quotient);
      xp10 -= 100 * quotient;
    } else {
      if (remainder + 4 < min_width) {
	memcpy(memseta(start, 32, min_width - (remainder + 4)), wbuf, remainder);
	start = &(start[min_width - 4]);
      } else {
	start = memcpya(start, wbuf, remainder);
      }
      start = memcpya(start, "e-", 2);
    }
    return memcpya(start, &(digit2_table[xp10 * 2]), 2);
  } else if (dxx >= 9999.4999999999) {
    // 4 sig fig exponential notation, large
    if (dxx >= 9.9994999999999e15) {
      if (dxx >= 9.9994999999999e127) {
	if (dxx == INFINITY) {
	  if (min_width > 4) {
	    start = memseta(start, 32, min_width - 4);
	  }
	  if (wpos == wbuf) {
	    return memcpya(start, " inf", 4);
	  } else {
	    return memcpya(start, "-inf", 4);
	  }
	} else if (dxx >= 9.9994999999999e255) {
	  dxx *= 1.0e-256;
	  xp10 |= 256;
	} else {
	  dxx *= 1.0e-128;
	  xp10 |= 128;
	}
      }
      if (dxx >= 9.9994999999999e63) {
	dxx *= 1.0e-64;
	xp10 |= 64;
      }
      if (dxx >= 9.9994999999999e31) {
	dxx *= 1.0e-32;
	xp10 |= 32;
      }
      if (dxx >= 9.9994999999999e15) {
	dxx *= 1.0e-16;
	xp10 |= 16;
      }
    }
    if (dxx >= 9.9994999999999e7) {
      dxx *= 1.0e-8;
      xp10 |= 8;
    }
    if (dxx >= 9.9994999999999e3) {
      dxx *= 1.0e-4;
      xp10 |= 4;
    }
    if (dxx >= 9.9994999999999e1) {
      dxx *= 1.0e-2;
      xp10 |= 2;
    }
    if (dxx >= 9.9994999999999e0) {
      dxx *= 1.0e-1;
      xp10++;
    }
    double_bround3(dxx, banker_round10, &quotient, &remainder);
    wpos = qrtoa_1p3(quotient, remainder, wpos);
    remainder = wpos - wbuf;
    if (xp10 >= 100) {
      if (remainder + 5 < min_width) {
	memcpy(memseta(start, 32, min_width - (remainder + 5)), wbuf, remainder);
	start = &(start[min_width - 5]);
      } else {
	start = memcpya(start, wbuf, remainder);
      }
      quotient = xp10 / 100;
      start = memcpyax(start, "e+", 2, '0' + quotient);
      xp10 -= 100 * quotient;
    } else {
      if (remainder + 4 < min_width) {
	memcpy(memseta(start, 32, min_width - (remainder + 4)), wbuf, remainder);
	start = &(start[min_width - 4]);
      } else {
	start = memcpya(start, wbuf, remainder);
      }
      start = memcpya(start, "e+", 2);
    }
    return memcpya(start, &(digit2_table[xp10 * 2]), 2);
  } else {
    if (dxx >= 0.99994999999999) {
      wpos = dtoa_so4(dxx, wpos);
    } else {
      // 4 sig fig decimal, no less than ~0.0001
      wpos = memcpya(wpos, "0.", 2);
      if (dxx < 9.9994999999999e-3) {
	dxx *= 100;
	wpos = memcpya(wpos, "00", 2);
      }
      if (dxx < 9.9994999999999e-2) {
	dxx *= 10;
	*wpos++ = '0';
      }
      wpos = uitoa_trunc4(double_bround(dxx * 10000, banker_round10), wpos);
    }
    remainder = wpos - wbuf;
    if (remainder < min_width) {
      memcpy(memseta(start, 32, min_width - remainder), wbuf, remainder);
      return &(start[min_width]);
    } else {
      return memcpya(start, wbuf, remainder);
    }
  }
}

char* dtoa_g_wxp8(double dxx, uint32_t min_width, char* start) {
  uint32_t xp10 = 0;
  char wbuf[16];
  char* wpos = wbuf;
  uint32_t quotient;
  uint32_t remainder;
  if (dxx != dxx) {
    if (min_width > 3) {
      start = memseta(start, 32, min_width - 3);
    }
    return memcpyl3a(start, "nan");
  } else if (dxx < 0) {
    *wpos++ = '-';
    dxx = -dxx;
  }
  if (dxx < 9.9999999499999e-5) {
    // 8 sig fig exponential notation, small
    if (dxx < 9.9999999499999e-16) {
      if (dxx < 9.9999999499999e-128) {
	if (dxx == 0.0) {
          memset(start, 32, min_width - 1);
	  start[min_width - 1] = '0';
	  return &(start[min_width]);
        } else if (dxx < 9.9999999499999e-256) {
	  dxx *= 1.0e256;
	  xp10 |= 256;
	} else {
	  dxx *= 1.0e128;
	  xp10 |= 128;
	}
      }
      if (dxx < 9.9999999499999e-64) {
	dxx *= 1.0e64;
	xp10 |= 64;
      }
      if (dxx < 9.9999999499999e-32) {
	dxx *= 1.0e32;
	xp10 |= 32;
      }
      if (dxx < 9.9999999499999e-16) {
	dxx *= 1.0e16;
	xp10 |= 16;
      }
    }
    if (dxx < 9.9999999499999e-8) {
      dxx *= 100000000;
      xp10 |= 8;
    }
    if (dxx < 9.9999999499999e-4) {
      dxx *= 10000;
      xp10 |= 4;
    }
    if (dxx < 9.9999999499999e-2) {
      dxx *= 100;
      xp10 |= 2;
    }
    if (dxx < 9.9999999499999e-1) {
      dxx *= 10;
      xp10++;
    }
    double_bround7(dxx, banker_round6, &quotient, &remainder);
    wpos = qrtoa_1p7(quotient, remainder, wpos);
    remainder = wpos - wbuf;
    if (xp10 >= 100) {
      if (remainder + 5 < min_width) {
	memcpy(memseta(start, 32, min_width - (remainder + 5)), wbuf, remainder);
	start = &(start[min_width - 5]);
      } else {
	start = memcpya(start, wbuf, remainder);
      }
      quotient = xp10 / 100;
      start = memcpyax(start, "e-", 2, '0' + quotient);
      xp10 -= 100 * quotient;
    } else {
      if (remainder + 4 < min_width) {
	memcpy(memseta(start, 32, min_width - (remainder + 4)), wbuf, remainder);
	start = &(start[min_width - 4]);
      } else {
	start = memcpya(start, wbuf, remainder);
      }
      start = memcpya(start, "e-", 2);
    }
    return memcpya(start, &(digit2_table[xp10 * 2]), 2);
  } else if (dxx >= 99999999.499999) {
    // 8 sig fig exponential notation, large
    if (dxx >= 9.9999999499999e15) {
      if (dxx >= 9.9999999499999e127) {
	if (dxx == INFINITY) {
	  if (min_width > 4) {
	    start = memseta(start, 32, min_width - 4);
	  }
	  if (wpos == wbuf) {
	    return memcpya(start, " inf", 4);
	  } else {
	    return memcpya(start, "-inf", 4);
	  }
	} else if (dxx >= 9.9999999499999e255) {
	  dxx *= 1.0e-256;
	  xp10 |= 256;
	} else {
	  dxx *= 1.0e-128;
	  xp10 |= 128;
	}
      }
      if (dxx >= 9.9999999499999e63) {
	dxx *= 1.0e-64;
	xp10 |= 64;
      }
      if (dxx >= 9.9999999499999e31) {
	dxx *= 1.0e-32;
	xp10 |= 32;
      }
      if (dxx >= 9.9999999499999e15) {
	dxx *= 1.0e-16;
	xp10 |= 16;
      }
    }
    if (dxx >= 9.9999999499999e7) {
      dxx *= 1.0e-8;
      xp10 |= 8;
    }
    if (dxx >= 9.9999999499999e3) {
      dxx *= 1.0e-4;
      xp10 |= 4;
    }
    if (dxx >= 9.9999999499999e1) {
      dxx *= 1.0e-2;
      xp10 |= 2;
    }
    if (dxx >= 9.9999999499999e0) {
      dxx *= 1.0e-1;
      xp10++;
    }
    double_bround7(dxx, banker_round6, &quotient, &remainder);
    wpos = qrtoa_1p7(quotient, remainder, wpos);
    remainder = wpos - wbuf;
    if (xp10 >= 100) {
      if (remainder + 5 < min_width) {
	memcpy(memseta(start, 32, min_width - (remainder + 5)), wbuf, remainder);
	start = &(start[min_width - 5]);
      } else {
	start = memcpya(start, wbuf, remainder);
      }
      quotient = xp10 / 100;
      start = memcpyax(start, "e+", 2, '0' + quotient);
      xp10 -= 100 * quotient;
    } else {
      if (remainder + 4 < min_width) {
	memcpy(memseta(start, 32, min_width - (remainder + 4)), wbuf, remainder);
	start = &(start[min_width - 4]);
      } else {
	start = memcpya(start, wbuf, remainder);
      }
      start = memcpya(start, "e+", 2);
    }
    return memcpya(start, &(digit2_table[xp10 * 2]), 2);
  } else {
    if (dxx >= 0.99999999499999) {
      wpos = dtoa_so8(dxx, wpos);
    } else {
      // 8 sig fig decimal, no less than ~0.0001
      wpos = memcpya(wpos, "0.", 2);
      if (dxx < 9.9999999499999e-3) {
	dxx *= 100;
	wpos = memcpya(wpos, "00", 2);
      }
      if (dxx < 9.9999999499999e-2) {
	dxx *= 10;
	*wpos++ = '0';
      }
      wpos = uitoa_trunc8(double_bround(dxx * 100000000, banker_round6), wpos);
    }
    remainder = wpos - wbuf;
    if (remainder < min_width) {
      memcpy(memseta(start, 32, min_width - remainder), wbuf, remainder);
      return &(start[min_width]);
    } else {
      return memcpya(start, wbuf, remainder);
    }
  }
}

char* chrom_print_human(uint32_t num, char* buf) {
  uint32_t n10;
  if (num < 10) {
    *buf = '0' + num;
    return &(buf[1]);
  } else if (num < 23) {
    n10 = num / 10;
    *buf = '0' + n10;
    buf[1] = '0' + (num - 10 * n10);
    return &(buf[2]);
  } else if (num < 25) {
    // X is 24th letter of alphabet, and 23rd chromosome
    *buf = 'A' + num;
    return &(buf[1]);
  } else if (num > 26) {
    // --allow-extra-chr 0
    *buf = '0';
    return &(buf[1]);
  } else if (num == 25) {
    memcpy(buf, "XY", 2);
    return &(buf[2]);
  } else {
    memcpy(buf, "MT", 2);
    return &(buf[2]);
  }
}

void magic_num(uint32_t divisor, uint64_t* multp, uint32_t* __restrict pre_shiftp, uint32_t* __restrict post_shiftp, uint32_t* __restrict incrp) {
  // Enables fast integer division by a constant not known until runtime.  See
  // http://ridiculousfish.com/blog/posts/labor-of-division-episode-iii.html .
  // Assumes divisor is not zero, of course.
  // May want to populate a struct instead.
  uint32_t down_multiplier = 0;
  uint32_t down_exponent = 0;
  uint32_t has_magic_down = 0;
  uint32_t quotient;
  uint32_t remainder;
  uint32_t ceil_log_2_d;
  uint32_t exponent;
  uint32_t uii;
  if (divisor & (divisor - 1)) {
    quotient = 0x80000000U / divisor;
    remainder = 0x80000000U - (quotient * divisor);
    ceil_log_2_d = 32 - __builtin_clz(divisor);
    for (exponent = 0; ; exponent++) {
      if (remainder >= divisor - remainder) {
        quotient = quotient * 2 + 1;
	remainder = remainder * 2 - divisor;
      } else {
	quotient = quotient * 2;
	remainder = remainder * 2;
      }
      if ((exponent >= ceil_log_2_d) || (divisor - remainder) <= (1U << exponent)) {
	break;
      }
      if ((!has_magic_down) && (remainder <= (1U << exponent))) {
	has_magic_down = 1;
	down_multiplier = quotient;
	down_exponent = exponent;
      }
    }
    if (exponent < ceil_log_2_d) {
      *multp = quotient + 1;
      *pre_shiftp = 0;
      *post_shiftp = 32 + exponent;
      *incrp = 0;
      return;
    } else if (divisor & 1) {
      *multp = down_multiplier;
      *pre_shiftp = 0;
      *post_shiftp = 32 + down_exponent;
      *incrp = 1;
      return;
    } else {
      *pre_shiftp = __builtin_ctz(divisor);
      magic_num(divisor >> (*pre_shiftp), multp, &uii, post_shiftp, incrp);
      return;
    }
  } else {
    // power of 2
    *multp = 1;
    *pre_shiftp = 0;
    *post_shiftp = __builtin_ctz(divisor);
    *incrp = 0;
  }
}

void fill_bits(uintptr_t loc_start, uintptr_t len, uintptr_t* bitarr) {
  assert(len);
  uintptr_t maj_start = loc_start / BITCT;
  uintptr_t maj_end = (loc_start + len) / BITCT;
  uintptr_t minor;
  if (maj_start == maj_end) {
    bitarr[maj_start] |= (ONELU << ((loc_start + len) % BITCT)) - (ONELU << (loc_start % BITCT));
  } else {
    bitarr[maj_start] |= ~((ONELU << (loc_start % BITCT)) - ONELU);
    fill_ulong_one(maj_end - maj_start - 1, &(bitarr[maj_start + 1]));
    minor = (loc_start + len) % BITCT;
    if (minor) {
      bitarr[maj_end] |= (ONELU << minor) - ONELU;
    }
  }
}

void clear_bits(uintptr_t loc_start, uintptr_t len, uintptr_t* bitarr) {
  assert(len);
  uintptr_t maj_start = loc_start / BITCT;
  uintptr_t maj_end = (loc_start + len) / BITCT;
  uintptr_t minor;
  if (maj_start == maj_end) {
    bitarr[maj_start] &= ~((ONELU << ((loc_start + len) % BITCT)) - (ONELU << (loc_start % BITCT)));
  } else {
    bitarr[maj_start] &= ((ONELU << (loc_start % BITCT)) - ONELU);
    fill_ulong_zero(maj_end - maj_start - 1, &(bitarr[maj_start + 1]));
    minor = (loc_start + len) % BITCT;
    if (minor) {
      bitarr[maj_end] &= ~((ONELU << minor) - ONELU);
    }
  }
}

uint32_t next_unset_unsafe(const uintptr_t* bitarr, uint32_t loc) {
  const uintptr_t* bitarr_ptr = &(bitarr[loc / BITCT]);
  uintptr_t ulii = (~(*bitarr_ptr)) >> (loc % BITCT);
  if (ulii) {
    return loc + CTZLU(ulii);
  }
  do {
    ulii = *(++bitarr_ptr);
  } while (ulii == ~ZEROLU);
  return ((uintptr_t)(bitarr_ptr - bitarr)) * BITCT + CTZLU(~ulii);
}

#ifdef __LP64__
uintptr_t next_unset_ul_unsafe(const uintptr_t* bitarr, uintptr_t loc) {
  const uintptr_t* bitarr_ptr = &(bitarr[loc / BITCT]);
  uintptr_t ulii = (~(*bitarr_ptr)) >> (loc % BITCT);
  if (ulii) {
    return loc + CTZLU(ulii);
  }
  do {
    ulii = *(++bitarr_ptr);
  } while (ulii == ~ZEROLU);
  return (((uintptr_t)(bitarr_ptr - bitarr)) * BITCT + CTZLU(~ulii));
}
#endif

uint32_t next_unset(const uintptr_t* bitarr, uint32_t loc, uint32_t ceil) {
  // safe version.
  assert(ceil >= 1);
  const uintptr_t* bitarr_ptr = &(bitarr[loc / BITCT]);
  uintptr_t ulii = (~(*bitarr_ptr)) >> (loc % BITCT);
  const uintptr_t* bitarr_last;
  if (ulii) {
    loc += CTZLU(ulii);
    return MINV(loc, ceil);
  }
  bitarr_last = &(bitarr[(ceil - 1) / BITCT]);
  do {
    if (bitarr_ptr >= bitarr_last) {
      return ceil;
    }
    ulii = *(++bitarr_ptr);
  } while (ulii == ~ZEROLU);
  loc = ((uintptr_t)(bitarr_ptr - bitarr)) * BITCT + CTZLU(~ulii);
  return MINV(loc, ceil);
}

#ifdef __LP64__
uintptr_t next_unset_ul(const uintptr_t* bitarr, uintptr_t loc, uintptr_t ceil) {
  const uintptr_t* bitarr_ptr = &(bitarr[loc / BITCT]);
  uintptr_t ulii = (~(*bitarr_ptr)) >> (loc % BITCT);
  const uintptr_t* bitarr_last;
  if (ulii) {
    ulii = loc + CTZLU(ulii);
    return MINV(ulii, ceil);
  }
  bitarr_last = &(bitarr[(ceil - 1) / BITCT]);
  do {
    if (bitarr_ptr >= bitarr_last) {
      return ceil;
    }
    ulii = *(++bitarr_ptr);
  } while (ulii == ~ZEROLU);
  ulii = ((uintptr_t)(bitarr_ptr - bitarr)) * BITCT + CTZLU(~ulii);
  return MINV(ulii, ceil);
}
#endif

uint32_t next_set_unsafe(const uintptr_t* bitarr, uint32_t loc) {
  const uintptr_t* bitarr_ptr = &(bitarr[loc / BITCT]);
  uintptr_t ulii = (*bitarr_ptr) >> (loc % BITCT);
  if (ulii) {
    return loc + CTZLU(ulii);
  }
  do {
    ulii = *(++bitarr_ptr);
  } while (!ulii);
  return ((uintptr_t)(bitarr_ptr - bitarr)) * BITCT + CTZLU(ulii);
}

#ifdef __LP64__
uintptr_t next_set_ul_unsafe(const uintptr_t* bitarr, uintptr_t loc) {
  const uintptr_t* bitarr_ptr = &(bitarr[loc / BITCT]);
  uintptr_t ulii = (*bitarr_ptr) >> (loc % BITCT);
  if (ulii) {
    return loc + CTZLU(ulii);
  }
  do {
    ulii = *(++bitarr_ptr);
  } while (!ulii);
  return ((uintptr_t)(bitarr_ptr - bitarr)) * BITCT + CTZLU(ulii);
}
#endif

uint32_t next_set(const uintptr_t* bitarr, uint32_t loc, uint32_t ceil) {
  const uintptr_t* bitarr_ptr = &(bitarr[loc / BITCT]);
  uintptr_t ulii = (*bitarr_ptr) >> (loc % BITCT);
  const uintptr_t* bitarr_last;
  uint32_t rval;
  if (ulii) {
    rval = loc + CTZLU(ulii);
    return MINV(rval, ceil);
  }
  bitarr_last = &(bitarr[(ceil - 1) / BITCT]);
  do {
    if (bitarr_ptr >= bitarr_last) {
      return ceil;
    }
    ulii = *(++bitarr_ptr);
  } while (!ulii);
  rval = ((uintptr_t)(bitarr_ptr - bitarr)) * BITCT + CTZLU(ulii);
  return MINV(rval, ceil);
}

#ifdef __LP64__
uintptr_t next_set_ul(const uintptr_t* bitarr, uintptr_t loc, uintptr_t ceil) {
  const uintptr_t* bitarr_ptr = &(bitarr[loc / BITCT]);
  uintptr_t ulii = (*bitarr_ptr) >> (loc % BITCT);
  const uintptr_t* bitarr_last;
  if (ulii) {
    ulii = loc + CTZLU(ulii);
    return MINV(ulii, ceil);
  }
  bitarr_last = &(bitarr[(ceil - 1) / BITCT]);
  do {
    if (bitarr_ptr >= bitarr_last) {
      return ceil;
    }
    ulii = *(++bitarr_ptr);
  } while (!ulii);
  ulii = ((uintptr_t)(bitarr_ptr - bitarr)) * BITCT + CTZLU(ulii);
  return MINV(ulii, ceil);
}
#endif

int32_t last_set_bit(const uintptr_t* bitarr, uint32_t word_ct) {
  const uintptr_t* bitarr_ptr = &(bitarr[word_ct]);
  uintptr_t ulii;
  do {
    ulii = *(--bitarr_ptr);
    if (ulii) {
      return ((uintptr_t)(bitarr_ptr - bitarr)) * BITCT + BITCT - 1 - CLZLU(ulii);
    }
  } while (bitarr_ptr > bitarr);
  return -1;
}

int32_t last_clear_bit(const uintptr_t* bitarr, uint32_t ceil) {
  // can return ceil or any lower number
  const uintptr_t* bitarr_ptr = &(bitarr[ceil / BITCT]);
  uint32_t remainder = ceil % BITCT;
  uintptr_t ulii;
  if (remainder) {
    ulii = (~(*bitarr_ptr)) & ((ONELU << remainder) - ONELU);
    if (ulii) {
      return (ceil | (BITCT - 1)) - CLZLU(ulii);
    }
  }
  while (bitarr_ptr > bitarr) {
    ulii = ~(*(--bitarr_ptr));
    if (ulii) {
      return ((uintptr_t)(bitarr_ptr - bitarr)) * BITCT + BITCT - 1 - CLZLU(ulii);
    }
  }
  return -1;
}

uint32_t prev_unset_unsafe(const uintptr_t* bitarr, uint32_t loc) {
  // unlike the next_{un}set family, this always returns a STRICTLY earlier
  // position
  const uintptr_t* bitarr_ptr = &(bitarr[loc / BITCT]);
  uint32_t remainder = loc % BITCT;
  uintptr_t ulii;
  if (remainder) {
    ulii = (~(*bitarr_ptr)) & ((ONELU << remainder) - ONELU);
    if (ulii) {
      return (loc | (BITCT - 1)) - CLZLU(ulii);
    }
  }
  do {
    ulii = ~(*(--bitarr_ptr));
  } while (!ulii);
  return ((uintptr_t)(bitarr_ptr - bitarr)) * BITCT + BITCT - 1 - CLZLU(ulii);
}

/*
uint32_t prev_unset(uintptr_t* bitarr, uint32_t loc, uint32_t floor) {
  uintptr_t* bitarr_ptr = &(bitarr[loc / BITCT]);
  uint32_t remainder = loc % BITCT;
  uintptr_t* bitarr_first;
  uintptr_t ulii;
  if (remainder) {
    ulii = (~(*bitarr_ptr)) & ((ONELU << remainder) - ONELU);
    if (ulii) {
      loc = (loc | (BITCT - 1)) - CLZLU(ulii);
      return MAXV(loc, floor);
    }
  }
  bitarr_first = &(bitarr[floor / BITCT]);
  do {
    if (bitarr_ptr == bitarr_first) {
      return floor;
    }
    ulii = ~(*(--bitarr_ptr));
  } while (!ulii);
  loc = ((uintptr_t)(bitarr_ptr - bitarr)) * BITCT + BITCT - 1 - CLZLU(ulii);
  return MAXV(loc, floor);
}
*/


int32_t bigstack_calloc_uc(uintptr_t ct, unsigned char** ucp_ptr) {
  *ucp_ptr = (unsigned char*)bigstack_alloc(ct);
  if (!(*ucp_ptr)) {
    return 1;
  }
  memset(*ucp_ptr, 0, ct);
  return 0;
}

int32_t bigstack_calloc_d(uintptr_t ct, double** dp_ptr) {
  *dp_ptr = (double*)bigstack_alloc(ct * sizeof(double));
  if (!(*dp_ptr)) {
    return 1;
  }
  fill_double_zero(ct, *dp_ptr);
  return 0;
}

int32_t bigstack_calloc_f(uintptr_t ct, float** fp_ptr) {
  *fp_ptr = (float*)bigstack_alloc(ct * sizeof(float));
  if (!(*fp_ptr)) {
    return 1;
  }
  fill_float_zero(ct, *fp_ptr);
  return 0;
}

int32_t bigstack_calloc_ui(uintptr_t ct, uint32_t** uip_ptr) {
  *uip_ptr = (uint32_t*)bigstack_alloc(ct * sizeof(int32_t));
  if (!(*uip_ptr)) {
    return 1;
  }
  fill_uint_zero(ct, *uip_ptr);
  return 0;
}

int32_t bigstack_calloc_ul(uintptr_t ct, uintptr_t** ulp_ptr) {
  *ulp_ptr = (uintptr_t*)bigstack_alloc(ct * sizeof(intptr_t));
  if (!(*ulp_ptr)) {
    return 1;
  }
  fill_ulong_zero(ct, *ulp_ptr);
  return 0;
}

int32_t bigstack_calloc_ull(uintptr_t ct, uint64_t** ullp_ptr) {
  *ullp_ptr = (uint64_t*)bigstack_alloc(ct * sizeof(int64_t));
  if (!(*ullp_ptr)) {
    return 1;
  }
  fill_ull_zero(ct, *ullp_ptr);
  return 0;
}

int32_t bigstack_end_calloc_uc(uintptr_t ct, unsigned char** ucp_ptr) {
  *ucp_ptr = (unsigned char*)bigstack_end_alloc(ct);
  if (!(*ucp_ptr)) {
    return 1;
  }
  memset(*ucp_ptr, 0, ct);
  return 0;
}

int32_t bigstack_end_calloc_d(uintptr_t ct, double** dp_ptr) {
  *dp_ptr = (double*)bigstack_end_alloc(ct * sizeof(double));
  if (!(*dp_ptr)) {
    return 1;
  }
  fill_double_zero(ct, *dp_ptr);
  return 0;
}

int32_t bigstack_end_calloc_f(uintptr_t ct, float** fp_ptr) {
  *fp_ptr = (float*)bigstack_end_alloc(ct * sizeof(float));
  if (!(*fp_ptr)) {
    return 1;
  }
  fill_float_zero(ct, *fp_ptr);
  return 0;
}

int32_t bigstack_end_calloc_ui(uintptr_t ct, uint32_t** uip_ptr) {
  *uip_ptr = (uint32_t*)bigstack_end_alloc(ct * sizeof(int32_t));
  if (!(*uip_ptr)) {
    return 1;
  }
  fill_uint_zero(ct, *uip_ptr);
  return 0;
}

int32_t bigstack_end_calloc_ul(uintptr_t ct, uintptr_t** ulp_ptr) {
  *ulp_ptr = (uintptr_t*)bigstack_end_alloc(ct * sizeof(intptr_t));
  if (!(*ulp_ptr)) {
    return 1;
  }
  fill_ulong_zero(ct, *ulp_ptr);
  return 0;
}

int32_t bigstack_end_calloc_ull(uintptr_t ct, uint64_t** ullp_ptr) {
  *ullp_ptr = (uint64_t*)bigstack_end_alloc(ct * sizeof(int64_t));
  if (!(*ullp_ptr)) {
    return 1;
  }
  fill_ull_zero(ct, *ullp_ptr);
  return 0;
}


// MurmurHash3, from
// https://code.google.com/p/smhasher/source/browse/trunk/MurmurHash3.cpp

static inline uint32_t rotl32(uint32_t x, int8_t r) {
  return (x << r) | (x >> (32 - r));
}

static inline uint32_t getblock32(const uint32_t* p, int i) {
  return p[i];
}

//-----------------------------------------------------------------------------
// Finalization mix - force all bits of a hash block to avalanche

static inline uint32_t fmix32(uint32_t h) {
  h ^= h >> 16;
  h *= 0x85ebca6b;
  h ^= h >> 13;
  h *= 0xc2b2ae35;
  h ^= h >> 16;

  return h;
}

uint32_t murmurhash3_32(const void* key, uint32_t len) {
  const uint8_t* data = (const uint8_t*)key;
  const int32_t nblocks = len / 4;

  uint32_t h1 = 0;
  // uint32_t h1 = seed;

  const uint32_t c1 = 0xcc9e2d51;
  const uint32_t c2 = 0x1b873593;

  //----------
  // body

  const uint32_t* blocks = (const uint32_t*)(data + nblocks*4);

  int32_t i;
  uint32_t k1;
  for(i = -nblocks; i; i++) {
      k1 = getblock32(blocks,i);

      k1 *= c1;
      k1 = rotl32(k1,15);
      k1 *= c2;

      h1 ^= k1;
      h1 = rotl32(h1,13);
      h1 = h1*5+0xe6546b64;
  }

  //----------
  // tail

  const uint8_t* tail = (const uint8_t*)(data + nblocks*4);

  k1 = 0;

  switch(len & 3) {
    case 3:
      k1 ^= tail[2] << 16;
      // fall through
    case 2:
      k1 ^= tail[1] << 8;
      // fall through
    case 1:
      k1 ^= tail[0];
      k1 *= c1;
      k1 = rotl32(k1,15);
      k1 *= c2;
      h1 ^= k1;
  }

  //----------
  // finalization

  h1 ^= len;

  return fmix32(h1);
}

uint32_t is_composite6(uintptr_t num) {
  // assumes num is congruent to 1 or 5 mod 6.
  // can speed this up by ~50% by hardcoding avoidance of multiples of 5/7,
  // but this isn't currently a bottleneck so I'll keep this simple
  uintptr_t divisor = 5;
  while (divisor * divisor <= num) {
    if (!(num % divisor)) {
      return 1;
    }
    divisor += 2;
    if (!(num % divisor)) {
      return 1;
    }
    divisor += 4;
  }
  return 0;
}

uintptr_t geqprime(uintptr_t floor) {
  // assumes floor is odd and greater than 1.  Returns 5 if floor = 3,
  // otherwise returns the first prime >= floor.
  uintptr_t ulii = floor % 3;
  if (!ulii) {
    floor += 2;
  } else if (ulii == 1) {
    goto geqprime_1mod6;
  }
  while (is_composite6(floor)) {
    floor += 2;
  geqprime_1mod6:
    if (!is_composite6(floor)) {
      return floor;
    }
    floor += 4;
  }
  return floor;
}

int32_t populate_id_htable(uintptr_t unfiltered_ct, const uintptr_t* exclude_arr, uintptr_t item_ct, const char* item_ids, uintptr_t max_id_len, uint32_t store_dups, uint32_t id_htable_size, uint32_t* id_htable) {
  // While unique IDs are normally assumed (and enforced) here, --extract and
  // --exclude are an exception, since we want to be able to e.g. exclude all
  // variants named '.'.  Since there could be millions of them, ordinary
  // O(n^2) hash table duplicate resolution is unacceptably slow; instead, we
  // allocate additional linked lists past the end of id_htable to track all
  // unfiltered indexes of duplicate names.  (This requires the
  // alloc_and_populate_id_htable interface; bigstack_end_alloc doesn't work
  // there.)
  uintptr_t item_uidx = 0;
  uint32_t extra_alloc = 0;
  uint32_t prev_llidx = 0;
  // needs to be synced with extract_exclude_flag_norange()
  uint32_t* extra_alloc_base = (uint32_t*)g_bigstack_base;
  uint32_t item_idx = 0;
  const char* sptr;
  uintptr_t prev_uidx;
  uintptr_t cur_bigstack_left;
  uint32_t max_extra_alloc;
  uint32_t slen;
  uint32_t hashval;
  uint32_t next_incr;
  uint32_t top_diff;
  uint32_t hash_result;
  uint32_t cur_dup;
  fill_uint_one(id_htable_size, id_htable);
  if (!store_dups) {
    for (; item_idx < item_ct; item_uidx++, item_idx++) {
      next_unset_ul_unsafe_ck(exclude_arr, &item_uidx);
      sptr = &(item_ids[item_uidx * max_id_len]);
      slen = strlen(sptr);
      hashval = murmurhash3_32(sptr, slen) % id_htable_size;
      next_incr = 1;
      while (1) {
	hash_result = id_htable[hashval];
	if (hash_result == 0xffffffffU) {
	  id_htable[hashval] = item_uidx;
	  break;
	} else if (!memcmp(sptr, &(item_ids[hash_result * max_id_len]), slen + 1)) {
	  // could add an allow_dups parameter which controls whether this is
	  // an error
	  LOGERRPRINTFWW("Error: Duplicate ID '%s'.\n", sptr);
	  return RET_INVALID_FORMAT;
	}
	// defend against overflow
	top_diff = id_htable_size - hashval;
	if (top_diff > next_incr) {
	  hashval += next_incr;
	} else {
	  hashval = next_incr - top_diff;
	}
	next_incr += 2; // quadratic probing
      }
    }
  } else {
    cur_bigstack_left = bigstack_left();
#ifdef __LP64__
    if (cur_bigstack_left >= 0x400000000LLU) {
      max_extra_alloc = 0xfffffffeU;
    } else {
      max_extra_alloc = cur_bigstack_left / sizeof(int32_t);
    }
#else
    max_extra_alloc = cur_bigstack_left / sizeof(int32_t);
#endif
    for (; item_idx < item_ct; item_uidx++, item_idx++) {
      next_unset_ul_unsafe_ck(exclude_arr, &item_uidx);
      sptr = &(item_ids[item_uidx * max_id_len]);
      slen = strlen(sptr);
      hashval = murmurhash3_32(sptr, slen) % id_htable_size;
      next_incr = 1;
      while (1) {
	hash_result = id_htable[hashval];
	if (hash_result == 0xffffffffU) {
	  id_htable[hashval] = item_uidx;
	  break;
        } else {
	  cur_dup = hash_result >> 31;
          if (cur_dup) {
	    prev_llidx = hash_result << 1;
	    prev_uidx = extra_alloc_base[prev_llidx];
          } else {
	    prev_uidx = hash_result;
          }
          if (!memcmp(sptr, &(item_ids[prev_uidx * max_id_len]), slen + 1)) {
	    if (extra_alloc + 4 > max_extra_alloc) {
	      return RET_NOMEM;
	    }
	    // point to linked list entry instead
	    if (!cur_dup) {
	      extra_alloc_base[extra_alloc] = hash_result;
	      extra_alloc_base[extra_alloc + 1] = 0xffffffffU; // list end
	      prev_llidx = extra_alloc;
	      extra_alloc += 2;
	    }
	    extra_alloc_base[extra_alloc] = item_uidx;
	    extra_alloc_base[extra_alloc + 1] = prev_llidx;
	    id_htable[hashval] = 0x80000000U | (extra_alloc >> 1);
	    extra_alloc += 2;
	    break; // bugfix
          }
	}
	top_diff = id_htable_size - hashval;
	if (top_diff > next_incr) {
	  hashval += next_incr;
	} else {
	  hashval = next_incr - top_diff;
	}
	next_incr += 2;
      }
    }
    if (extra_alloc) {
      bigstack_alloc(extra_alloc * sizeof(int32_t));
    }
  }
  return 0;
}

uint32_t id_htable_find(const char* id_buf, uintptr_t cur_id_len, const uint32_t* id_htable, uint32_t id_htable_size, const char* item_ids, uintptr_t max_id_len) {
  // assumes no duplicate entries, and nonzero id_htable_size
  // returns 0xffffffffU on failure
  if (cur_id_len >= max_id_len) {
    return 0xffffffffU;
  }
  uint32_t hashval = murmurhash3_32(id_buf, cur_id_len) % id_htable_size;
  uint32_t next_incr = 1;
  const char* sptr;
  uint32_t hash_result;
  uint32_t top_diff;
  while (1) {
    hash_result = id_htable[hashval];
    if (hash_result == 0xffffffffU) {
      return 0xffffffffU;
    }
    sptr = &(item_ids[hash_result * max_id_len]);
    if ((!memcmp(id_buf, sptr, cur_id_len)) && (!sptr[cur_id_len])) {
      return hash_result;
    }
    top_diff = id_htable_size - hashval;
    if (top_diff > next_incr) {
      hashval += next_incr;
    } else {
      hashval = next_incr - top_diff;
    }
    next_incr += 2;
  }
}

void fill_idx_to_uidx(const uintptr_t* exclude_arr, uintptr_t unfiltered_item_ct, uintptr_t item_ct, uint32_t* idx_to_uidx) {
  uint32_t* idx_to_uidx_end = &(idx_to_uidx[item_ct]);
  uint32_t item_uidx = 0;
  uint32_t item_uidx_stop;
  while (idx_to_uidx < idx_to_uidx_end) {
    item_uidx = next_unset_unsafe(exclude_arr, item_uidx);
    item_uidx_stop = next_set(exclude_arr, item_uidx, unfiltered_item_ct);
    do {
      *idx_to_uidx++ = item_uidx++;
    } while (item_uidx < item_uidx_stop);
  }
}

void fill_idx_to_uidx_incl(const uintptr_t* include_arr, uintptr_t unfiltered_item_ct, uintptr_t item_ct, uint32_t* idx_to_uidx) {
  uint32_t* idx_to_uidx_end = &(idx_to_uidx[item_ct]);
  uint32_t item_uidx = 0;
  uint32_t item_uidx_stop;
  while (idx_to_uidx < idx_to_uidx_end) {
    item_uidx = next_set_unsafe(include_arr, item_uidx);
    item_uidx_stop = next_unset(include_arr, item_uidx, unfiltered_item_ct);
    do {
      *idx_to_uidx++ = item_uidx++;
    } while (item_uidx < item_uidx_stop);
  }
}

void fill_uidx_to_idx(const uintptr_t* exclude_arr, uint32_t unfiltered_item_ct, uint32_t item_ct, uint32_t* uidx_to_idx) {
  uint32_t item_uidx = 0;
  uint32_t item_idx = 0;
  uint32_t* uidx_to_idx_ptr;
  uint32_t* uidx_to_idx_stop;
  while (item_idx < item_ct) {
    item_uidx = next_unset_unsafe(exclude_arr, item_uidx);
    uidx_to_idx_ptr = &(uidx_to_idx[item_uidx]);
    item_uidx = next_set(exclude_arr, item_uidx, unfiltered_item_ct);
    uidx_to_idx_stop = &(uidx_to_idx[item_uidx]);
    do {
      *uidx_to_idx_ptr++ = item_idx++;
    } while (uidx_to_idx_ptr < uidx_to_idx_stop);
  }
}

void fill_uidx_to_idx_incl(const uintptr_t* include_arr, uint32_t unfiltered_item_ct, uint32_t item_ct, uint32_t* uidx_to_idx) {
  uint32_t item_uidx = 0;
  uint32_t item_idx = 0;
  uint32_t* uidx_to_idx_ptr;
  uint32_t* uidx_to_idx_stop;
  while (item_idx < item_ct) {
    item_uidx = next_set_unsafe(include_arr, item_uidx);
    uidx_to_idx_ptr = &(uidx_to_idx[item_uidx]);
    item_uidx = next_unset(include_arr, item_uidx, unfiltered_item_ct);
    uidx_to_idx_stop = &(uidx_to_idx[item_uidx]);
    do {
      *uidx_to_idx_ptr++ = item_idx++;
    } while (uidx_to_idx_ptr < uidx_to_idx_stop);
  }
}

void fill_midx_to_idx(const uintptr_t* exclude_arr_orig, const uintptr_t* exclude_arr, uint32_t item_ct, uint32_t* midx_to_idx) {
  // Assumes item_ct is nonzero.

  // May want to switch to alternate behavior: when current midx is excluded,
  // fill midx_to_idx[] with the next item_idx.
  uint32_t item_uidx = next_unset_unsafe(exclude_arr_orig, 0);
  uint32_t item_idx = 0;
  uint32_t item_midx;
  for (item_midx = 0; item_idx < item_ct; item_uidx++, item_midx++) {
    next_unset_unsafe_ck(exclude_arr_orig, &item_uidx);
    if (!IS_SET(exclude_arr, item_uidx)) {
      midx_to_idx[item_midx] = item_idx++;
    }
  }
}

void fill_quatervec_55(uint32_t ct, uintptr_t* quatervec) {
  uint32_t rem = ct & (BITCT - 1);
#ifdef __LP64__
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  __m128i* vecp = (__m128i*)quatervec;
  __m128i* vec_end = (__m128i*)(&(quatervec[2 * (ct / BITCT)]));
  uintptr_t* second_to_last;
  while (vecp < vec_end) {
    *vecp++ = m1;
  }
  if (rem) {
    second_to_last = (uintptr_t*)vecp;
    if (rem > BITCT2) {
      second_to_last[0] = FIVEMASK;
      second_to_last[1] = FIVEMASK >> ((BITCT - rem) * 2);
    } else {
      second_to_last[0] = FIVEMASK >> ((BITCT2 - rem) * 2);
      second_to_last[1] = 0;
    }
  }
#else
  uintptr_t* vec_end = &(quatervec[2 * (ct / BITCT)]);
  while (quatervec < vec_end) {
    *quatervec++ = FIVEMASK;
  }
  if (rem) {
    if (rem > BITCT2) {
      quatervec[0] = FIVEMASK;
      quatervec[1] = FIVEMASK >> ((BITCT - rem) * 2);
    } else {
      quatervec[0] = FIVEMASK >> ((BITCT2 - rem) * 2);
      quatervec[1] = 0;
    }
  }
#endif
}

void quaterarr_collapse_init(const uintptr_t* __restrict unfiltered_bitarr, uint32_t unfiltered_ct, const uintptr_t* __restrict filter_bitarr, uint32_t filtered_ct, uintptr_t* __restrict output_quaterarr) {
  // Used to unpack e.g. unfiltered sex_male to a filtered quaterarr usable as
  // a raw input bitmask.
  // Assumes output_quaterarr is sized to a multiple of 16 bytes.
  uintptr_t cur_write = 0;
  uint32_t item_uidx = 0;
  uint32_t write_bit = 0;
  uint32_t item_idx = 0;
  uint32_t item_uidx_stop;
  while (item_idx < filtered_ct) {
    item_uidx = next_set_unsafe(filter_bitarr, item_uidx);
    item_uidx_stop = next_unset(filter_bitarr, item_uidx, unfiltered_ct);
    item_idx += item_uidx_stop - item_uidx;
    do {
      cur_write |= ((unfiltered_bitarr[item_uidx / BITCT] >> (item_uidx % BITCT)) & 1) << (write_bit * 2);
      if (++write_bit == BITCT2) {
	*output_quaterarr++ = cur_write;
        cur_write = 0;
	write_bit = 0;
      }
    } while (++item_uidx < item_uidx_stop);
  }
  if (write_bit) {
    *output_quaterarr++ = cur_write;
  }
  if ((filtered_ct + (BITCT2 - 1)) & BITCT2) {
    *output_quaterarr = 0;
  }
}

void quaterarr_collapse_init_exclude(const uintptr_t* __restrict unfiltered_bitarr, uint32_t unfiltered_ct, const uintptr_t* __restrict filter_exclude_bitarr, uint32_t filtered_ct, uintptr_t* __restrict output_quaterarr) {
  uintptr_t cur_write = 0;
  uint32_t item_uidx = 0;
  uint32_t write_bit = 0;
  uint32_t item_idx = 0;
  uint32_t item_uidx_stop;
  while (item_idx < filtered_ct) {
    item_uidx = next_unset_unsafe(filter_exclude_bitarr, item_uidx);
    item_uidx_stop = next_set(filter_exclude_bitarr, item_uidx, unfiltered_ct);
    item_idx += item_uidx_stop - item_uidx;
    do {
      cur_write |= ((unfiltered_bitarr[item_uidx / BITCT] >> (item_uidx % BITCT)) & 1) << (write_bit * 2);
      if (++write_bit == BITCT2) {
	*output_quaterarr++ = cur_write;
        cur_write = 0;
	write_bit = 0;
      }
    } while (++item_uidx < item_uidx_stop);
  }
  if (write_bit) {
    *output_quaterarr++ = cur_write;
  }
  if ((filtered_ct + (BITCT2 - 1)) & BITCT2) {
    *output_quaterarr = 0;
  }
}

uint32_t alloc_collapsed_haploid_filters(const uintptr_t* __restrict sample_bitarr, const uintptr_t* __restrict sex_male, uint32_t unfiltered_sample_ct, uint32_t sample_ct, uint32_t hh_exists, uint32_t is_include, uintptr_t** sample_include_quatervec_ptr, uintptr_t** sample_male_include_quatervec_ptr) {
  uintptr_t sample_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(sample_ct);
  if (hh_exists & (Y_FIX_NEEDED | NXMHH_EXISTS)) {
    // if already allocated, we assume this is fully initialized
    if (!(*sample_include_quatervec_ptr)) {
      if (bigstack_alloc_ul(sample_ctv2, sample_include_quatervec_ptr)) {
	return 1;
      }
      fill_quatervec_55(sample_ct, *sample_include_quatervec_ptr);
    }
  }
  if (hh_exists & (XMHH_EXISTS | Y_FIX_NEEDED)) {
    // if already allocated, we assume it's been bigstack_end_alloc'd but not
    // initialized
    if (!(*sample_male_include_quatervec_ptr)) {
      if (bigstack_alloc_ul(sample_ctv2, sample_male_include_quatervec_ptr)) {
	return 1;
      }
    }
    if (is_include) {
      quaterarr_collapse_init(sex_male, unfiltered_sample_ct, sample_bitarr, sample_ct, *sample_male_include_quatervec_ptr);
    } else {
      quaterarr_collapse_init_exclude(sex_male, unfiltered_sample_ct, sample_bitarr, sample_ct, *sample_male_include_quatervec_ptr);
    }
  }
  return 0;
}

void sample_delim_convert(uintptr_t unfiltered_sample_ct, const uintptr_t* sample_exclude, uint32_t sample_ct, uintptr_t max_sample_id_len, char oldc, char newc, char* sample_ids) {
  // assumes there is exactly one delimiter to convert per name
  uintptr_t sample_uidx = 0;
  uint32_t sample_idx;
  char* nptr;
  for (sample_idx = 0; sample_idx < sample_ct; sample_uidx++, sample_idx++) {
    next_unset_ul_unsafe_ck(sample_exclude, &sample_uidx);
    nptr = (char*)memchr(&(sample_ids[sample_uidx * max_sample_id_len]), (unsigned char)oldc, max_sample_id_len);
    *nptr = newc;
  }
}

void get_set_wrange_align(const uintptr_t* __restrict bitarr, uintptr_t word_ct, uintptr_t* __restrict firstw_ptr, uintptr_t* __restrict wlen_ptr) {
  const uintptr_t* bitarr_ptr = bitarr;
  const uintptr_t* bitarr_end = &(bitarr[word_ct]);
#ifdef __LP64__
  const uintptr_t* bitarr_end2 = &(bitarr[word_ct & (~ONELU)]);
  while (bitarr_ptr < bitarr_end2) {
    if (bitarr_ptr[0] || bitarr_ptr[1]) {
      *firstw_ptr = (uintptr_t)(bitarr_ptr - bitarr);
      while (!(*(--bitarr_end)));
      *wlen_ptr = 1 + (uintptr_t)(bitarr_end - bitarr_ptr);
      return;
    }
    bitarr_ptr = &(bitarr_ptr[2]);
  }
  if ((bitarr_end2 != bitarr_end) && (*bitarr_end2)) {
    *firstw_ptr = word_ct - 1;
    *wlen_ptr = 1;
    return;
  }
#else
  while (bitarr_ptr < bitarr_end) {
    if (*bitarr_ptr) {
      *firstw_ptr = (uintptr_t)(bitarr_ptr - bitarr);
      while (!(*(--bitarr_end)));
      *wlen_ptr = 1 + (uintptr_t)(bitarr_end - bitarr_ptr);
      return;
    }
    bitarr_ptr++;
  }
#endif
  *firstw_ptr = 0;
  *wlen_ptr = 0;
}

// hashval computation left to caller since this is frequently used with
// chromosome IDs, where the compiler can optimize the integer modulus
// operation since the hash table size is preset
uint32_t unklen_id_htable_find(const char* cur_id, const char* const* item_ids, const uint32_t* id_htable, uint32_t hashval, uint32_t id_htable_size) {
  // returns 0xffffffffU on failure
  uint32_t next_incr = 1;
  while (1) {
    const uint32_t hash_result = id_htable[hashval];
    if (hash_result == 0xffffffffU) {
      return 0xffffffffU;
    }
    const char* htable_entry = item_ids[hash_result];
    if (!strcmp(cur_id, htable_entry)) {
      return hash_result;
    }
    const uint32_t top_diff = id_htable_size - hashval;
    if (top_diff > next_incr) {
      hashval += next_incr;
    } else {
      hashval = next_incr - top_diff;
    }
    next_incr += 2;
  }
}

static inline uint32_t nonstd_chrom_name_htable_find(const char* chrom_name, const char* const* nonstd_names, const uint32_t* nonstd_id_htable, uint32_t name_slen) {
  const uint32_t hashval = murmurhash3_32(chrom_name, name_slen) % CHROM_NAME_HTABLE_SIZE;
  return unklen_id_htable_find(chrom_name, nonstd_names, nonstd_id_htable, hashval, CHROM_NAME_HTABLE_SIZE);
}


// Global since species_str() may be called by functions which don't actually
// care about chrom_info.  (chrom_info is really a global variable too, but I
// find it easier to maintain this code when chrom_info dependencies are made
// explicit in the function signatures; in contrast, g_species_singular and
// g_species_plural are just for pretty printing and lend no insight into what
// the functions which reference them are doing.)
const char* g_species_singular = nullptr;
const char* g_species_plural = nullptr;

int32_t init_chrom_info(Chrom_info* chrom_info_ptr) {
  // "constructor".  initializes with maximum capacity.  doesn't use bigstack.
  // chrom_mask, haploid_mask: bits
  // chrom_file_order, chrom_idx_to_foidx: int32s
  // chrom_fo_vidx_start: int32s, with an extra trailing element
  // nonstd_names: intptr_ts
  // nonstd_id_htable: CHROM_NAME_HTABLE_SIZE int32s

  assert(!(MAX_POSSIBLE_CHROM % VEC_BYTES));
  const uintptr_t vecs_required = 2 * BITCT_TO_VECCT(MAX_POSSIBLE_CHROM) + 3 * (MAX_POSSIBLE_CHROM / VEC_INT32) + 1 + (MAX_POSSIBLE_CHROM / VEC_WORDS) + (CHROM_NAME_HTABLE_SIZE + (VEC_INT32 - 1)) / VEC_INT32;

  // needed for proper cleanup
  chrom_info_ptr->name_ct = 0;
  chrom_info_ptr->incl_excl_name_stack = nullptr;
  if (aligned_malloc(vecs_required * VEC_BYTES, &(chrom_info_ptr->chrom_mask))) {
    return RET_NOMEM;
  }
  uintptr_t* alloc_iter = &(chrom_info_ptr->chrom_mask[BITCT_TO_VECCT(MAX_POSSIBLE_CHROM) * VEC_WORDS]);
  chrom_info_ptr->haploid_mask = alloc_iter;
  alloc_iter = &(alloc_iter[BITCT_TO_VECCT(MAX_POSSIBLE_CHROM) * VEC_WORDS]);
  chrom_info_ptr->chrom_file_order = (uint32_t*)alloc_iter;
  alloc_iter = &(alloc_iter[(MAX_POSSIBLE_CHROM / VEC_INT32) * VEC_WORDS]);
  chrom_info_ptr->chrom_fo_vidx_start = (uint32_t*)alloc_iter;
  alloc_iter = &(alloc_iter[((MAX_POSSIBLE_CHROM / VEC_INT32) + 1) * VEC_WORDS]);
  chrom_info_ptr->chrom_idx_to_foidx = (uint32_t*)alloc_iter;
  alloc_iter = &(alloc_iter[(MAX_POSSIBLE_CHROM / VEC_INT32) * VEC_WORDS]);
  chrom_info_ptr->nonstd_names = (char**)alloc_iter;
  alloc_iter = &(alloc_iter[MAX_POSSIBLE_CHROM]);
  chrom_info_ptr->nonstd_id_htable = (uint32_t*)alloc_iter;
  // alloc_iter = &(alloc_iter[((CHROM_NAME_HTABLE_SIZE + (VEC_INT32 - 1)) / VEC_INT32) * VEC_WORDS]);
  // postpone nonstd_id_htable initialization until first nonstandard ID is
  // loaded
  // fill_uint_one(CHROM_NAME_HTABLE_SIZE, chrom_info_ptr->nonstd_id_htable);
  return 0;
}

// if these are defined within init_species(), they may not persist after
// function exit
static const char species_singular_constants[][7] = {"person", "cow", "dog", "horse", "mouse", "plant", "sheep", "sample"};
static const char species_plural_constants[][8] = {"people", "cattle", "dogs", "horses", "mice", "plants", "sheep", "samples"};

void init_species(uint32_t species_code, Chrom_info* chrom_info_ptr) {
  // human: 22, X, Y, XY, MT
  // cow: 29, X, Y, MT
  // dog: 38, X, Y, XY, MT
  // horse: 31, X, Y
  // mouse: 19, X, Y
  // rice: 12
  // sheep: 26, X, Y
  const int32_t species_xymt_codes[] = {
    23, 24, 25, 26,
    30, 31, -2, 33,
    39, 40, 41, 42,
    32, 33, -2, -2,
    20, 21, -2, -2,
    -2, -2, -2, -2,
    27, 28, -2, -2};
  const uint32_t species_autosome_ct[] = {22, 29, 38, 31, 19, 12, 26};
  const uint32_t species_max_code[] = {26, 33, 42, 33, 21, 12, 28};
  fill_ulong_zero(CHROM_MASK_WORDS, chrom_info_ptr->chrom_mask);
  chrom_info_ptr->output_encoding = 0;
  chrom_info_ptr->zero_extra_chroms = 0;
  chrom_info_ptr->species = species_code;
  chrom_info_ptr->is_include_stack = 0;
  g_species_singular = species_singular_constants[species_code];
  g_species_plural = species_plural_constants[species_code];
  if (species_code != SPECIES_UNKNOWN) {
    // these are assumed to be already initialized in the SPECIES_UNKNOWN case

    // bugfix: haploid_mask was being cleared in --chr-set case
    fill_ulong_zero(CHROM_MASK_WORDS, chrom_info_ptr->haploid_mask);
    memcpy(chrom_info_ptr->xymt_codes, &(species_xymt_codes[species_code * XYMT_OFFSET_CT]), XYMT_OFFSET_CT * sizeof(int32_t));
    chrom_info_ptr->autosome_ct = species_autosome_ct[species_code];
    chrom_info_ptr->max_code = species_max_code[species_code];
    switch (species_code) {
    case SPECIES_HUMAN:
      chrom_info_ptr->haploid_mask[0] = 0x1800000;
      break;
    case SPECIES_COW:
      chrom_info_ptr->haploid_mask[0] = 0xc0000000LU;
      break;
    case SPECIES_DOG:
#ifdef __LP64__
      chrom_info_ptr->haploid_mask[0] = 0x18000000000LLU;
#else
      chrom_info_ptr->haploid_mask[1] = 0x180;
#endif
      break;
    case SPECIES_HORSE:
#ifdef __LP64__
      chrom_info_ptr->haploid_mask[0] = 0x300000000LLU;
#else
      chrom_info_ptr->haploid_mask[1] = 3;
#endif
      break;
    case SPECIES_MOUSE:
      chrom_info_ptr->haploid_mask[0] = 0x300000;
      break;
    case SPECIES_RICE:
      chrom_info_ptr->haploid_mask[0] = 0x1fff;
      break;
    case SPECIES_SHEEP:
      chrom_info_ptr->haploid_mask[0] = 0x18000000;
      break;
    }
  }
  fill_uint_one(chrom_info_ptr->max_code + 1, chrom_info_ptr->chrom_idx_to_foidx);
}

void init_default_chrom_mask(Chrom_info* chrom_info_ptr) {
  if (chrom_info_ptr->species != SPECIES_UNKNOWN) {
    fill_all_bits(chrom_info_ptr->max_code + 1, chrom_info_ptr->chrom_mask);
  } else {
    fill_all_bits(chrom_info_ptr->autosome_ct + 1, chrom_info_ptr->chrom_mask);
    // --chr-set support
    for (uint32_t xymt_idx = 0; xymt_idx < XYMT_OFFSET_CT; ++xymt_idx) {
      int32_t cur_code = chrom_info_ptr->xymt_codes[xymt_idx];
      if (cur_code != -2) {
	set_bit(chrom_info_ptr->xymt_codes[xymt_idx], chrom_info_ptr->chrom_mask);
      }
    }
  }
}

void forget_extra_chrom_names(uint32_t reinitialize, Chrom_info* chrom_info_ptr) {
  const uint32_t name_ct = chrom_info_ptr->name_ct;
  // guard against init_species() not being called yet
  if (name_ct) {
    char** nonstd_names = chrom_info_ptr->nonstd_names;
    const uint32_t chrom_idx_last = chrom_info_ptr->max_code + name_ct;
    for (uint32_t chrom_idx = chrom_info_ptr->max_code + 1; chrom_idx <= chrom_idx_last; ++chrom_idx) {
      free(nonstd_names[chrom_idx]);
      nonstd_names[chrom_idx] = nullptr;
    }
    if (reinitialize) {
      fill_uint_one(CHROM_NAME_HTABLE_SIZE, chrom_info_ptr->nonstd_id_htable);
      chrom_info_ptr->name_ct = 0;
    }
  }
}

int32_t finalize_chrom_info(Chrom_info* chrom_info_ptr) {
  const uint32_t chrom_ct = chrom_info_ptr->chrom_ct;
  const uint32_t name_ct = chrom_info_ptr->name_ct;
  const uint32_t chrom_code_end = chrom_info_ptr->max_code + 1 + name_ct;
  const uint32_t chrom_code_bitvec_ct = BITCT_TO_VECCT(chrom_code_end);
  const uint32_t chrom_ct_int32vec_ct = (chrom_ct + (VEC_INT32 - 1)) / VEC_INT32;
  const uint32_t chrom_ct_p1_int32vec_ct = 1 + (chrom_ct / VEC_INT32);
  const uint32_t chrom_code_end_int32vec_ct = (chrom_code_end + (VEC_INT32 - 1)) / VEC_INT32;
  const uint32_t chrom_code_end_wordvec_ct = (chrom_code_end + (VEC_WORDS - 1)) / VEC_WORDS;
  uint32_t final_vecs_required = 2 * chrom_code_bitvec_ct + chrom_ct_int32vec_ct + chrom_ct_p1_int32vec_ct + chrom_code_end_int32vec_ct;
  if (name_ct) {
    final_vecs_required += chrom_code_end_wordvec_ct + (CHROM_NAME_HTABLE_SIZE + (VEC_INT32 - 1)) / VEC_INT32;
  }
  uintptr_t* new_alloc;
  if (aligned_malloc(final_vecs_required * VEC_BYTES, &new_alloc)) {
    return RET_NOMEM;
  }
  uintptr_t* old_alloc = chrom_info_ptr->chrom_mask;
  uintptr_t* new_alloc_iter = new_alloc;

  memcpy(new_alloc_iter, chrom_info_ptr->chrom_mask, chrom_code_bitvec_ct * VEC_BYTES);
  chrom_info_ptr->chrom_mask = new_alloc_iter;
  new_alloc_iter = &(new_alloc_iter[chrom_code_bitvec_ct * VEC_WORDS]);

  memcpy(new_alloc_iter, chrom_info_ptr->haploid_mask, chrom_code_bitvec_ct * VEC_BYTES);
  chrom_info_ptr->haploid_mask = new_alloc_iter;
  new_alloc_iter = &(new_alloc_iter[chrom_code_bitvec_ct * VEC_WORDS]);

  memcpy(new_alloc_iter, chrom_info_ptr->chrom_file_order, chrom_ct_int32vec_ct * VEC_BYTES);
  chrom_info_ptr->chrom_file_order = (uint32_t*)new_alloc_iter;
  new_alloc_iter = &(new_alloc_iter[chrom_ct_int32vec_ct * VEC_WORDS]);

  memcpy(new_alloc_iter, chrom_info_ptr->chrom_fo_vidx_start, chrom_ct_p1_int32vec_ct * VEC_BYTES);
  chrom_info_ptr->chrom_fo_vidx_start = (uint32_t*)new_alloc_iter;
  new_alloc_iter = &(new_alloc_iter[chrom_ct_p1_int32vec_ct * VEC_WORDS]);

  memcpy(new_alloc_iter, chrom_info_ptr->chrom_idx_to_foidx, chrom_code_end_int32vec_ct * VEC_BYTES);
  chrom_info_ptr->chrom_idx_to_foidx = (uint32_t*)new_alloc_iter;

  if (!name_ct) {
    chrom_info_ptr->nonstd_names = nullptr;
    chrom_info_ptr->nonstd_id_htable = nullptr;
  } else {
    new_alloc_iter = &(new_alloc_iter[chrom_code_end_int32vec_ct * VEC_WORDS]);

    memcpy(new_alloc_iter, chrom_info_ptr->nonstd_names, chrom_code_end_wordvec_ct * VEC_BYTES);
    chrom_info_ptr->nonstd_names = (char**)new_alloc_iter;
    new_alloc_iter = &(new_alloc_iter[chrom_code_end_wordvec_ct * VEC_WORDS]);

    memcpy(new_alloc_iter, chrom_info_ptr->nonstd_id_htable, CHROM_NAME_HTABLE_SIZE * sizeof(int32_t));
    chrom_info_ptr->nonstd_id_htable = (uint32_t*)new_alloc_iter;
  }
  aligned_free(old_alloc);
  return 0;
}

void cleanup_chrom_info(Chrom_info* chrom_info_ptr) {
  if (chrom_info_ptr->chrom_mask) {
    // bugfix: this must happened before aligned_free() call
    forget_extra_chrom_names(0, chrom_info_ptr);

    aligned_free(chrom_info_ptr->chrom_mask);
    chrom_info_ptr->chrom_mask = nullptr;
  }
  Ll_str* ll_str_ptr = chrom_info_ptr->incl_excl_name_stack;
  while (ll_str_ptr) {
    Ll_str* next_ptr = ll_str_ptr->next;
    free(ll_str_ptr);
    ll_str_ptr = next_ptr;
  }
  chrom_info_ptr->incl_excl_name_stack = nullptr;
}

char* chrom_name_std(const Chrom_info* chrom_info_ptr, uint32_t chrom_idx, char* buf) {
  const uint32_t output_encoding = chrom_info_ptr->output_encoding;
  if (output_encoding & (CHR_OUTPUT_PREFIX | CHR_OUTPUT_0M)) {
    if (output_encoding == CHR_OUTPUT_0M) {
      // force two chars
      if (chrom_idx <= chrom_info_ptr->autosome_ct) {
	buf = (char*)memcpya(buf, &(digit2_table[chrom_idx * 2]), 2);
      } else if ((int32_t)chrom_idx == chrom_info_ptr->xymt_codes[XY_OFFSET]) {
	buf = (char*)memcpya(buf, "XY", 2);
      } else {
	*buf++ = '0';
	if ((int32_t)chrom_idx == chrom_info_ptr->xymt_codes[X_OFFSET]) {
	  *buf++ = 'X';
	} else {
	  // assumes only X/Y/XY/MT defined
	  *buf++ = ((int32_t)chrom_idx == chrom_info_ptr->xymt_codes[Y_OFFSET])? 'Y' : 'M';
	}
      }
      return buf;
    }
    buf = memcpyl3a(buf, "chr");
  }
  if ((!(output_encoding & (CHR_OUTPUT_M | CHR_OUTPUT_MT))) || (chrom_idx <= chrom_info_ptr->autosome_ct)) {
    return uint32toa(chrom_idx, buf);
  } else if ((int32_t)chrom_idx == chrom_info_ptr->xymt_codes[X_OFFSET]) {
    *buf++ = 'X';
  } else if ((int32_t)chrom_idx == chrom_info_ptr->xymt_codes[Y_OFFSET]) {
    *buf++ = 'Y';
  } else if ((int32_t)chrom_idx == chrom_info_ptr->xymt_codes[XY_OFFSET]) {
    buf = (char*)memcpya(buf, "XY", 2);
  } else {
    *buf++ = 'M';
    if (output_encoding & CHR_OUTPUT_MT) {
      *buf++ = 'T';
    }
  }
  return buf;
}

char* chrom_name_write(const Chrom_info* chrom_info_ptr, uint32_t chrom_idx, char* buf) {
  // assumes chrom_idx is valid
  if (!chrom_idx) {
    *buf++ = '0';
    return buf;
  } else if (chrom_idx <= chrom_info_ptr->max_code) {
    return chrom_name_std(chrom_info_ptr, chrom_idx, buf);
  } else if (chrom_info_ptr->zero_extra_chroms) {
    *buf++ = '0';
    return buf;
  } else {
    return strcpya(buf, chrom_info_ptr->nonstd_names[chrom_idx]);
  }
}

char* chrom_name_buf5w4write(const Chrom_info* chrom_info_ptr, uint32_t chrom_idx, uint32_t* chrom_name_len_ptr, char* buf5) {
  uint32_t slen;
  *chrom_name_len_ptr = 4;
  if (!chrom_idx) {
    memcpy(buf5, "   0", 4);
  } else if (chrom_idx <= chrom_info_ptr->max_code) {
    if (chrom_info_ptr->output_encoding & CHR_OUTPUT_PREFIX) {
      *chrom_name_len_ptr = (uintptr_t)(chrom_name_std(chrom_info_ptr, chrom_idx, buf5) - buf5);
    } else {
      width_force(4, buf5, chrom_name_std(chrom_info_ptr, chrom_idx, buf5));
    }
  } else if (chrom_info_ptr->zero_extra_chroms) {
    memcpy(buf5, "   0", 4);
  } else {
    slen = strlen(chrom_info_ptr->nonstd_names[chrom_idx]);
    if (slen < 4) {
      fw_strcpyn(4, slen, chrom_info_ptr->nonstd_names[chrom_idx], buf5);
    } else {
      *chrom_name_len_ptr = slen;
      return chrom_info_ptr->nonstd_names[chrom_idx];
    }
  }
  return buf5;
}

uint32_t get_max_chrom_slen(const Chrom_info* chrom_info_ptr) {
  // does not include trailing null
  // can be overestimate
  // if more functions start calling this, it should just be built into
  // load_bim() instead
  if (chrom_info_ptr->zero_extra_chroms) {
    return 3 + MAX_CHROM_TEXTNUM_SLEN;
  }
  const uint32_t chrom_ct = chrom_info_ptr->chrom_ct;
  const uint32_t max_code = chrom_info_ptr->max_code;
  uint32_t max_chrom_slen = 3 + MAX_CHROM_TEXTNUM_SLEN;
  for (uint32_t chrom_fo_idx = 0; chrom_fo_idx < chrom_ct; chrom_fo_idx++) {
    const uint32_t chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
    if (!is_set(chrom_info_ptr->chrom_mask, chrom_idx)) {
      continue;
    }
    if (chrom_idx > max_code) {
      const uint32_t name_slen = strlen(chrom_info_ptr->nonstd_names[chrom_idx]);
      if (name_slen > max_chrom_slen) {
	max_chrom_slen = name_slen;
      }
    }
  }
  return max_chrom_slen;
}

uint32_t haploid_chrom_present(const Chrom_info* chrom_info_ptr) {
  const uintptr_t* chrom_mask = chrom_info_ptr->chrom_mask;
  const uintptr_t* haploid_mask = chrom_info_ptr->haploid_mask;
  for (uint32_t widx = 0; widx < CHROM_MASK_INITIAL_WORDS; widx++) {
    if (chrom_mask[widx] & haploid_mask[widx]) {
      return 1;
    }
  }
  return 0;
}

static inline int32_t single_letter_chrom(uint32_t letter) {
  letter &= 0xdf;
  if (letter == 'X') {
    return CHROM_X;
  } else if (letter == 'Y') {
    return CHROM_Y;
  } else if (letter == 'M') {
    return CHROM_MT;
  } else {
    return -1;
  }
}

int32_t get_chrom_code_raw(const char* sptr) {
  // any character <= ' ' is considered a terminator
  // note that char arithmetic tends to be compiled to int32 operations, so we
  // mostly work with ints here
  // assumes MAX_CHROM_TEXTNUM_SLEN == 2
  uint32_t first_char_code = (unsigned char)sptr[0];
  uint32_t second_char_code = (unsigned char)sptr[1];
  if ((first_char_code & 0xdf) == 'C') {
    if (((second_char_code & 0xdf) == 'H') && ((((unsigned char)sptr[2]) & 0xdf) == 'R')) {
      sptr = &(sptr[3]);
      first_char_code = (unsigned char)sptr[0];
      second_char_code = (unsigned char)sptr[1];
    } else {
      return -1;
    }
  }
  if (second_char_code > ' ') {
    if (sptr[2] > ' ') {
      return -1;
    }
    const uint32_t first_char_toi = first_char_code - '0';
    if (first_char_toi < 10) {
      const uint32_t second_char_toi = second_char_code - '0';
      if (second_char_toi < 10) {
	return first_char_toi * 10 + second_char_toi;
      } else if (!first_char_toi) {
	// accept '0X', '0Y', '0M' emitted by Oxford software
	return single_letter_chrom(second_char_code);
      }
    } else {
      first_char_code &= 0xdf;
      if (first_char_code == 'X') {
        if ((second_char_code == 'Y') || (second_char_code == 'y')) {
	  return CHROM_XY;
	}
      } else if (first_char_code == 'M') {
        if ((second_char_code == 'T') || (second_char_code == 't')) {
	  return CHROM_MT;
	}
      }
    }
  } else {
    const uint32_t first_char_toi = first_char_code - '0';
    if (first_char_toi < 10) {
      return first_char_toi;
    } else {
      return single_letter_chrom(first_char_code);
    }
  }
  return -1;
}

int32_t get_chrom_code(const char* chrom_name, const Chrom_info* chrom_info_ptr, uint32_t name_slen) {
  // requires chrom_name to be null-terminated
  // in practice, name_slen will usually already be known, may as well avoid
  // redundant strlen() calls even though this uglifies the interface
  // does not perform exhaustive error-checking
  // -1 = --allow-extra-chr ok, -2 = total fail
  const int32_t chrom_code_raw = get_chrom_code_raw(chrom_name);
  if (((const uint32_t)chrom_code_raw) <= chrom_info_ptr->max_code) {
    return chrom_code_raw;
  }
  if (chrom_code_raw != -1) {
    if (chrom_code_raw >= MAX_POSSIBLE_CHROM) {
      return chrom_info_ptr->xymt_codes[chrom_code_raw - MAX_POSSIBLE_CHROM];
    }
    return -2;
  }
  if (!chrom_info_ptr->name_ct) {
    return -1;
  }
  // 0xffffffffU gets casted to -1
  return (int32_t)nonstd_chrom_name_htable_find(chrom_name, (const char* const*)chrom_info_ptr->nonstd_names, chrom_info_ptr->nonstd_id_htable, name_slen);
}

int32_t get_chrom_code_counted(const Chrom_info* chrom_info_ptr, uint32_t name_slen, char* chrom_name) {
  // when the chromosome name isn't null-terminated
  char* s_end = &(chrom_name[name_slen]);
  const char tmpc = *s_end;
  *s_end = '\0';
  const int32_t retval = get_chrom_code(chrom_name, chrom_info_ptr, name_slen);
  *s_end = tmpc;
  return retval;
}

uint32_t get_variant_chrom_fo_idx(const Chrom_info* chrom_info_ptr, uintptr_t variant_uidx) {
  const uint32_t* variant_binsearch = chrom_info_ptr->chrom_fo_vidx_start;
  uint32_t chrom_fo_min = 0;
  uint32_t chrom_ct = chrom_info_ptr->chrom_ct;
  while (chrom_ct - chrom_fo_min > 1) {
    const uint32_t chrom_fo_cur = (chrom_ct + chrom_fo_min) / 2;
    if (variant_binsearch[chrom_fo_cur] > variant_uidx) {
      chrom_ct = chrom_fo_cur;
    } else {
      chrom_fo_min = chrom_fo_cur;
    }
  }
  return chrom_fo_min;
}

void chrom_error(const char* chrom_name, const char* file_descrip, const Chrom_info* chrom_info_ptr, uintptr_t line_idx, int32_t error_code) {
  // assumes chrom_name is null-terminated
  const int32_t raw_code = get_chrom_code_raw(chrom_name);
  logprint("\n");
  if (line_idx) {
    LOGERRPRINTFWW("Error: Invalid chromosome code '%s' on line %" PRIuPTR " of %s.\n", chrom_name, line_idx, file_descrip);
  } else {
    LOGERRPRINTFWW("Error: Invalid chromosome code '%s' in %s.\n", chrom_name, file_descrip);
  }
  if ((raw_code > ((int32_t)chrom_info_ptr->max_code)) && ((raw_code <= MAX_CHROM_TEXTNUM + XYMT_OFFSET_CT) || (raw_code >= MAX_POSSIBLE_CHROM))) {
    if (chrom_info_ptr->species != SPECIES_UNKNOWN) {
      if (chrom_info_ptr->species == SPECIES_HUMAN) {
	logerrprint("(This is disallowed for humans.  Check if the problem is with your data, or if\nyou forgot to define a different chromosome set with e.g. --chr-set.)\n");
      } else {
	logerrprint("(This is disallowed by the PLINK 1.07 species flag you used.  You can\ntemporarily work around this restriction with --chr-set; contact the developers\nif you want the flag to be permanently redefined.)\n");
      }
    } else {
      logerrprint("(This is disallowed by your --chr-set/--autosome-num parameters.  Check if the\nproblem is with your data, or your command line.)\n");
    }
  } else if (error_code == -1) {
    logerrprint("(Use --allow-extra-chr to force it to be accepted.)\n");
  }
}

int32_t try_to_add_chrom_name(const char* chrom_name, const char* file_descrip, uintptr_t line_idx, uint32_t name_slen, uint32_t allow_extra_chroms, int32_t* chrom_idx_ptr, Chrom_info* chrom_info_ptr) {
  // assumes chrom_name is nonstandard (i.e. not "2", "chr2", "chrX", etc.)
  // requires chrom_name to be null-terminated
  // assumes chrom_idx currently has the return value of get_chrom_code()
  if ((!allow_extra_chroms) || ((*chrom_idx_ptr) == -2)) {
    chrom_error(chrom_name, file_descrip, chrom_info_ptr, line_idx, *chrom_idx_ptr);
    return RET_MALFORMED_INPUT;
  }

  // quasi-bugfix: remove redundant hash table check

  if (chrom_name[0] == '#') {
    // redundant with some of the comment-skipping loaders, but this isn't
    // performance-critical
    logprint("\n");
    logerrprint("Error: Chromosome/contig names may not begin with '#'.\n");
    return RET_MALFORMED_INPUT;
  }
  if (name_slen > MAX_ID_SLEN) {
    logprint("\n");
    if (line_idx) {
      LOGERRPRINTFWW("Error: Line %" PRIuPTR " of %s has an excessively long chromosome/contig name. (The " PROG_NAME_CAPS " limit is " MAX_ID_SLEN_STR " characters.)\n", line_idx, file_descrip);
    } else {
      LOGERRPRINTFWW("Error: Excessively long chromosome/contig name in %s. (The " PROG_NAME_CAPS " limit is " MAX_ID_SLEN_STR " characters.)\n", file_descrip);
    }
    return RET_MALFORMED_INPUT;
  }
  const uint32_t max_code_p1 = chrom_info_ptr->max_code + 1;
  const uint32_t name_ct = chrom_info_ptr->name_ct;
  const uint32_t chrom_code_end = max_code_p1 + name_ct;
  if (chrom_code_end == MAX_POSSIBLE_CHROM) {
    logprint("\n");
    logerrprint("Error: Too many distinct nonstandard chromosome/contig names.\n");
    return RET_MALFORMED_INPUT;
  }
  if (!name_ct) {
    // lazy initialization
    fill_uint_one(CHROM_NAME_HTABLE_SIZE, chrom_info_ptr->nonstd_id_htable);
  }
  char** nonstd_names = chrom_info_ptr->nonstd_names;
  nonstd_names[chrom_code_end] = (char*)malloc(name_slen + 1);
  if (!nonstd_names[chrom_code_end]) {
    return RET_NOMEM;
  }
  Ll_str* name_stack_ptr = chrom_info_ptr->incl_excl_name_stack;
  uint32_t in_name_stack = 0;
  while (name_stack_ptr) {
    // there shouldn't be many of these, so sorting is unimportant
    if (!strcmp(chrom_name, name_stack_ptr->ss)) {
      in_name_stack = 1;
      break;
    }
    name_stack_ptr = name_stack_ptr->next;
  }
  if ((in_name_stack && chrom_info_ptr->is_include_stack) || ((!in_name_stack) && (!chrom_info_ptr->is_include_stack))) {
    SET_BIT(chrom_code_end, chrom_info_ptr->chrom_mask);
    if (chrom_info_ptr->haploid_mask[0] & 1) {
      SET_BIT(chrom_code_end, chrom_info_ptr->haploid_mask);
    }
  }
  memcpy(nonstd_names[chrom_code_end], chrom_name, name_slen + 1);
  *chrom_idx_ptr = (int32_t)chrom_code_end;
  chrom_info_ptr->name_ct = name_ct + 1;
  uint32_t* id_htable = chrom_info_ptr->nonstd_id_htable;
  uint32_t hashval = murmurhash3_32(chrom_name, name_slen) % CHROM_NAME_HTABLE_SIZE;
  uint32_t next_incr = 1;
  while (1) {
    if (id_htable[hashval] == 0xffffffffU) {
      id_htable[hashval] = chrom_code_end;
      return 0;
    }
    // no overflow danger here
    hashval += next_incr;
    if (hashval >= CHROM_NAME_HTABLE_SIZE) {
      hashval -= CHROM_NAME_HTABLE_SIZE;
    }
    next_incr += 2; // quadratic probing
  }
}

uint32_t allele_set(const char* newval, uint32_t slen, char** allele_ptr) {
  char* newptr;
  if (slen == 1) {
    newptr = (char*)(&(g_one_char_strs[((unsigned char)*newval) * 2]));
  } else {
    newptr = (char*)malloc(slen + 1);
    if (!newptr) {
      return 1;
    }
    memcpyx(newptr, newval, slen, '\0');
  }
  *allele_ptr = newptr;
  return 0;
}

uint32_t allele_reset(const char* newval, uint32_t slen, char** allele_ptr) {
  char* newptr;
  if (slen == 1) {
    newptr = (char*)(&(g_one_char_strs[((uint8_t)*newval) * 2]));
  } else {
    newptr = (char*)malloc(slen + 1);
    if (!newptr) {
      return 1;
    }
    memcpyx(newptr, newval, slen, '\0');
  }
  if (allele_ptr[0][1]) {
    free(*allele_ptr);
  }
  *allele_ptr = newptr;
  return 0;
}

void cleanup_allele_storage(uint32_t max_allele_slen, uintptr_t allele_storage_entry_ct, char** allele_storage) {
  if (allele_storage && (max_allele_slen > 1)) {
    const uintptr_t one_char_strs_addr = (uintptr_t)g_one_char_strs;
    for (uintptr_t idx = 0; idx < allele_storage_entry_ct; ++idx) {
      char* cur_entry = allele_storage[idx];
      assert(cur_entry);
      // take advantage of unsigned wraparound
      if ((((uintptr_t)cur_entry) - one_char_strs_addr) >= 512) {
	free(cur_entry);
      }
    }
  }
}

void cleanup_allele_storage2(uintptr_t allele_storage_entry_ct, char** allele_storage) {
  if (allele_storage) {
    const uintptr_t one_char_strs_addr = (uintptr_t)g_one_char_strs;
    for (uintptr_t idx = 0; idx < allele_storage_entry_ct;) {
      char* cur_entry = allele_storage[idx];
      if (!cur_entry) {
        // --merge-equal-pos hacked entry
        idx += 2;
        continue;
      }
      // take advantage of unsigned wraparound
      if ((((uintptr_t)cur_entry) - one_char_strs_addr) >= 512) {
	free(cur_entry);
      }
      ++idx;
    }
  }
}

void refresh_chrom_info(const Chrom_info* chrom_info_ptr, uintptr_t marker_uidx, uint32_t* __restrict chrom_end_ptr, uint32_t* __restrict chrom_fo_idx_ptr, uint32_t* __restrict is_x_ptr, uint32_t* __restrict is_y_ptr, uint32_t* __restrict is_mt_ptr, uint32_t* __restrict is_haploid_ptr) {
  // assumes we are at the end of the chromosome denoted by chrom_fo_idx.  Ok
  // for chrom_fo_idx == 0xffffffffU.
  // assumes marker_uidx < unfiltered_marker_ct
  *chrom_end_ptr = chrom_info_ptr->chrom_fo_vidx_start[(*chrom_fo_idx_ptr) + 1];
  while (marker_uidx >= (*chrom_end_ptr)) {
    *chrom_end_ptr = chrom_info_ptr->chrom_fo_vidx_start[(++(*chrom_fo_idx_ptr)) + 1];
  }
  const int32_t chrom_idx = chrom_info_ptr->chrom_file_order[*chrom_fo_idx_ptr];
  *is_x_ptr = (chrom_idx == chrom_info_ptr->xymt_codes[X_OFFSET]);
  *is_y_ptr = (chrom_idx == chrom_info_ptr->xymt_codes[Y_OFFSET]);
  *is_mt_ptr = (chrom_idx == chrom_info_ptr->xymt_codes[MT_OFFSET]);
  *is_haploid_ptr = is_set(chrom_info_ptr->haploid_mask, chrom_idx);
}

int32_t single_chrom_start(const Chrom_info* chrom_info_ptr, const uintptr_t* marker_exclude, uint32_t unfiltered_marker_ct) {
  // Assumes there is at least one marker, and there are no split chromosomes.
  // Returns first marker_uidx in chromosome if there is only one, or -1 if
  // there's more than one chromosome.
  uint32_t first_marker_uidx = next_unset_unsafe(marker_exclude, 0);
  uint32_t last_marker_chrom = get_variant_chrom(chrom_info_ptr, last_clear_bit(marker_exclude, unfiltered_marker_ct));
  return (get_variant_chrom(chrom_info_ptr, first_marker_uidx) == last_marker_chrom)? first_marker_uidx : -1;
}

double get_dmedian(const double* sorted_arr, uintptr_t len) {
  if (len) {
    if (len % 2) {
      return sorted_arr[len / 2];
    } else {
      return (sorted_arr[len / 2] + sorted_arr[(len / 2) - 1]) * 0.5;
    }
  } else {
    return 0.0;
  }
}

#ifdef __cplusplus
double destructive_get_dmedian(uintptr_t len, double* unsorted_arr) {
  if (!len) {
    return 0.0;
  }
  uintptr_t len_d2 = len / 2;
  std::nth_element(unsorted_arr, &(unsorted_arr[len_d2]), &(unsorted_arr[len]));
  if (!(len % 2)) {
    std::nth_element(unsorted_arr, &(unsorted_arr[len_d2 - 1]), &(unsorted_arr[len_d2]));
    return (unsorted_arr[len_d2 - 1] + unsorted_arr[len_d2]) * 0.5;
  } else {
    return unsorted_arr[len_d2];
  }
}
#else
double destructive_get_dmedian(uintptr_t len, double* unsorted_arr) {
  // no, I'm not gonna bother reimplementing introselect just for folks who
  // insist on using gcc over g++
  qsort(unsorted_arr, len, sizeof(double), double_cmp);
  return get_dmedian(unsorted_arr, len);
}
#endif

int32_t strcmp_casted(const void* s1, const void* s2) {
  return strcmp((char*)s1, (char*)s2);
}

// PLINK 2's natural sort uses the following logic:
// - All alphabetic characters act as if they are capitalized, except for
// tiebreaking purposes (where ASCII is used).
// - Numbers are compared by magnitude, with the exception of...
// - Numbers with leading zero(es).  If you're putting extraneous zeroes in
// front of IDs, we assume they're there to force particular items to be sorted
// earlier, rather than just appearing at random.  So, unlike many natural sort
// implementations, we sort 00200 < 021 < 20: all numbers with n leading zeroes
// are sorted before all numbers with (n-1) leading zeroes; magnitude only
// applies if the leading zero counts match.  This handles e.g. subbasement
// room numbering properly.
//
// This won't always do what you want if your IDs have variable-length decimals
// in them (e.g. it yields 0.99 < 0.101); if you don't want to fall back on
// ASCII sort, enforce a fixed number of digits after the decimal point.  Also
// note that ASCII sort is outright better for e.g. numbers represented in
// hexadecimal or base 36.  In principle, it's possible to reliably autodetect
// some of these cases (especially hexadecimal numbers beginning with "0x"),
// but that'll never be perfect so we just let the user toggle the sort method.
int32_t strcmp_natural_scan_forward(const unsigned char* s1, const unsigned char* s2) {
  // assumes s1 and s2 currently point to the middle of a mismatching number,
  // where s1 < s2.
  unsigned char c1;
  unsigned char c2;
  do {
    c1 = *(++s1);
    c2 = *(++s2);
    if (is_not_digit(c1)) {
      return -1;
    }
  } while (is_digit(c2));
  return 1;
}

// We have the following major states:
//   0 (initial): strings perfectly match so far, last char (if any) is
//                nonnumeric.
//   1: strings perfectly match so far, last char is numeric.
//   2: strings match except for capitalization, last char is nonnumeric.
//   3: strings match except for capitalization, last char is numeric.
// strcmp_natural_tiebroken() expresses the logic for states 2 and 3, while
// strcmp_natural_uncasted() handles states 0 and 1.
int32_t strcmp_natural_tiebroken(const unsigned char* s1, const unsigned char* s2) {
  // assumes ties should be broken in favor of s2.
  unsigned char c1 = *(++s1);
  unsigned char c2 = *(++s2);
  while (is_not_nzdigit(c1) && is_not_nzdigit(c2)) {
    // state 2
  strcmp_natural_tiebroken_state_2:
    if (c1 != c2) {
      if ((c1 >= 'a') && (c1 <= 'z')) {
	c1 -= 32;
      }
      if ((c2 >= 'a') && (c2 <= 'z')) {
	c2 -= 32;
      }
      if (c1 < c2) {
	return -1;
      } else if (c1 > c2) {
	return 1;
      }
    } else if (!c1) {
      return -1;
    }
    c1 = *(++s1);
    c2 = *(++s2);
  }
  if (is_not_nzdigit(c1) || is_not_nzdigit(c2)) {
    return (c1 < c2)? -1 : 1;
  }
  do {
    // state 3
    if (c1 != c2) {
      if (is_digit(c2)) {
	if (c1 < c2) {
	  return strcmp_natural_scan_forward(s1, s2);
	} else {
	  return -strcmp_natural_scan_forward(s2, s1);
	}
      }
      return 1;
    }
    c1 = *(++s1);
    c2 = *(++s2);
  } while (is_digit(c1));
  if (is_digit(c2)) {
    return -1;
  }
  // skip the while (is_not_digit...) check
  goto strcmp_natural_tiebroken_state_2;
}

static inline int32_t strcmp_natural_uncasted(const unsigned char* s1, const unsigned char* s2) {
  unsigned char c1 = *s1;
  unsigned char c2 = *s2;
  while (is_not_nzdigit(c1) && is_not_nzdigit(c2)) {
    // state 0
  strcmp_natural_uncasted_state_0:
    if (c1 != c2) {
      if ((c1 >= 'a') && (c1 <= 'z')) {
	if (c2 + 32 == c1) {
	  return -strcmp_natural_tiebroken(s2, s1);
	} else if ((c2 < 'a') || (c2 > 'z')) {
	  c1 -= 32;
	}
      } else if ((c2 >= 'a') && (c2 <= 'z')) {
	c2 -= 32;
	if (c1 == c2) {
	  return strcmp_natural_tiebroken(s1, s2);
	}
      }
      return (c1 < c2)? -1 : 1;
    } else if (!c1) {
      return 0;
    }
    c1 = *(++s1);
    c2 = *(++s2);
  }
  if (is_not_nzdigit(c1) || is_not_nzdigit(c2)) {
    return (c1 < c2)? -1 : 1;
  }
  do {
    // state 1
    if (c1 != c2) {
      if (is_digit(c2)) {
	if (c1 < c2) {
	  return strcmp_natural_scan_forward(s1, s2);
	} else {
	  return -strcmp_natural_scan_forward(s2, s1);
	}
      }
      return 1;
    }
    c1 = *(++s1);
    c2 = *(++s2);
  } while (is_digit(c1));
  if (is_digit(c2)) {
    return -1;
  }
  goto strcmp_natural_uncasted_state_0;
}

int32_t strcmp_natural(const void* s1, const void* s2) {
  return strcmp_natural_uncasted((unsigned char*)s1, (unsigned char*)s2);
}

int32_t strcmp_deref(const void* s1, const void* s2) {
  return strcmp(*(char**)s1, *(char**)s2);
}

int32_t strcmp_natural_deref(const void* s1, const void* s2) {
  return strcmp_natural_uncasted(*(unsigned char**)s1, *(unsigned char**)s2);
}

int32_t get_uidx_from_unsorted(const char* idstr, const uintptr_t* exclude_arr, uint32_t id_ct, const char* unsorted_ids, uintptr_t max_id_len) {
  uintptr_t id_uidx = 0;
  uintptr_t slen_p1 = strlen(idstr) + 1;
  uint32_t id_idx;
  if (slen_p1 > max_id_len) {
    return -1;
  }
  for (id_idx = 0; id_idx < id_ct; id_uidx++, id_idx++) {
    id_uidx = next_unset_ul_unsafe(exclude_arr, id_uidx);
    if (!memcmp(idstr, &(unsorted_ids[id_uidx * max_id_len]), slen_p1)) {
      return (int32_t)((uint32_t)id_uidx);
    }
  }
  return -1;
}

char* scan_for_duplicate_ids(char* sorted_ids, uintptr_t id_ct, uintptr_t max_id_len) {
  uintptr_t id_idx;
  id_ct--;
  for (id_idx = 0; id_idx < id_ct; id_idx++) {
    if (!strcmp(&(sorted_ids[id_idx * max_id_len]), &(sorted_ids[(id_idx + 1) * max_id_len]))) {
      return &(sorted_ids[id_idx * max_id_len]);
    }
  }
  return nullptr;
}

char* scan_for_duplicate_or_overlap_ids(char* sorted_ids, uintptr_t id_ct, uintptr_t max_id_len, const char* sorted_nonoverlap_ids, uintptr_t nonoverlap_id_ct, uintptr_t max_nonoverlap_id_len) {
  // extended scan_for_duplicate_ids() which also verifies that no entry in
  // sorted_ids matches any entry in sorted_nonoverlap_ids.
  // nonoverlap_id_ct == 0 and sorted_nonoverlap_ids == nullptr ok.  id_ct
  // cannot be zero, though.
  uintptr_t nonoverlap_id_idx = 0;
  uintptr_t id_idx = 0;
  char* cur_id_ptr = sorted_ids;
  const char* nonoverlap_id_ptr;
  char* other_id_ptr;
  int32_t ii;
  while (1) {
    if (nonoverlap_id_idx == nonoverlap_id_ct) {
      return scan_for_duplicate_ids(cur_id_ptr, id_ct - id_idx, max_id_len);
    }
    nonoverlap_id_ptr = &(sorted_nonoverlap_ids[nonoverlap_id_idx * max_nonoverlap_id_len]);
    ii = strcmp(cur_id_ptr, nonoverlap_id_ptr);
    if (ii < 0) {
      if (++id_idx == id_ct) {
	return nullptr;
      }
      other_id_ptr = &(cur_id_ptr[max_id_len]);
      if (!strcmp(cur_id_ptr, other_id_ptr)) {
	return cur_id_ptr;
      }
      cur_id_ptr = other_id_ptr;
      continue;
    } else if (!ii) {
      return cur_id_ptr;
    }
    nonoverlap_id_idx++;
  }
}

int32_t eval_affection(const char* bufptr, double missing_phenod) {
  // turns out --1 had the side-effect of *forcing* case/control
  // interpretation in 1.07.  We replicate that for backward compatibility, and
  // no longer call this function in that context.
  char* ss;
  double dxx;
  // this used to be an integer read, but that could do the wrong thing if e.g.
  // all phenotypes were -9.xxx...
  dxx = strtod(bufptr, &ss);
  if ((ss == bufptr) || (dxx == missing_phenod)) {
    return 1;
  }
  return ((bufptr[0] == '0') || (bufptr[0] == '1') || (bufptr[0] == '2')) && is_space_or_eoln(bufptr[1]);
}

uint32_t triangle_divide(int64_t cur_prod, int32_t modif) {
  // return smallest integer vv for which (vv * (vv + modif)) is no smaller
  // than cur_prod, and neither term in the product is negative.  (Note the
  // lack of a divide by two; cur_prod should also be double its "true" value
  // as a result.)
  int64_t vv;
  if (cur_prod == 0) {
    if (modif < 0) {
      return -modif;
    } else {
      return 0;
    }
  }
  vv = (int64_t)sqrt((double)cur_prod);
  while ((vv - 1) * (vv + modif - 1) >= cur_prod) {
    vv--;
  }
  while (vv * (vv + modif) < cur_prod) {
    vv++;
  }
  return vv;
}

void parallel_bounds(uint32_t ct, int32_t start, uint32_t parallel_idx, uint32_t parallel_tot, int32_t* __restrict bound_start_ptr, int32_t* __restrict bound_end_ptr) {
  int32_t modif = 1 - start * 2;
  int64_t ct_tot = ((int64_t)ct) * (ct + modif);
  *bound_start_ptr = triangle_divide((ct_tot * parallel_idx) / parallel_tot, modif);
  *bound_end_ptr = triangle_divide((ct_tot * (parallel_idx + 1)) / parallel_tot, modif);
}

// this might belong in plink_calc instead, not being used anywhere else
// set align to 1 for no alignment
void triangle_fill(uint32_t ct, uint32_t pieces, uint32_t parallel_idx, uint32_t parallel_tot, uint32_t start, uint32_t align, uint32_t* target_arr) {
  int32_t modif = 1 - start * 2;
  uint32_t cur_piece = 1;
  int64_t ct_tr;
  int64_t cur_prod;
  int32_t lbound;
  int32_t ubound;
  uint32_t uii;
  uint32_t align_m1;
  parallel_bounds(ct, start, parallel_idx, parallel_tot, &lbound, &ubound);
  // x(x+1)/2 is divisible by y iff (x % (2y)) is 0 or (2y - 1).
  align *= 2;
  align_m1 = align - 1;
  target_arr[0] = lbound;
  target_arr[pieces] = ubound;
  cur_prod = ((int64_t)lbound) * (lbound + modif);
  ct_tr = (((int64_t)ubound) * (ubound + modif) - cur_prod) / pieces;
  while (cur_piece < pieces) {
    cur_prod += ct_tr;
    lbound = triangle_divide(cur_prod, modif);
    uii = (lbound - ((int32_t)start)) & align_m1;
    if ((uii) && (uii != align_m1)) {
      lbound = start + ((lbound - ((int32_t)start)) | align_m1);
    }
    // lack of this check caused a nasty bug earlier
    if (((uint32_t)lbound) > ct) {
      lbound = ct;
    }
    target_arr[cur_piece++] = lbound;
  }
}

int32_t relationship_req(uint64_t calculation_type) {
  return (calculation_type & (CALC_RELATIONSHIP | CALC_UNRELATED_HERITABILITY | CALC_REL_CUTOFF | CALC_REGRESS_REL | CALC_PCA))? 1 : 0;
}

int32_t distance_req(const char* read_dists_fname, uint64_t calculation_type) {
  return ((calculation_type & CALC_DISTANCE) || ((calculation_type & (CALC_PLINK1_DISTANCE_MATRIX | CALC_PLINK1_IBS_MATRIX)) && (!(calculation_type & CALC_GENOME))) || ((!read_dists_fname) && (calculation_type & (CALC_IBS_TEST | CALC_GROUPDIST | CALC_REGRESS_DISTANCE))));
}

int32_t double_cmp(const void* aa, const void* bb) {
  double cc = *((const double*)aa) - *((const double*)bb);
  if (cc > 0.0) {
    return 1;
  } else if (cc < 0.0) {
    return -1;
  } else {
    return 0;
  }
}

int32_t double_cmp_decr(const void* aa, const void* bb) {
  double cc = *((const double*)aa) - *((const double*)bb);
  if (cc < 0.0) {
    return 1;
  } else if (cc > 0.0) {
    return -1;
  } else {
    return 0;
  }
}

int32_t double_cmp_deref(const void* aa, const void* bb) {
  double cc = **((const double**)aa) - **((const double**)bb);
  if (cc > 0.0) {
    return 1;
  } else if (cc < 0.0) {
    return -1;
  } else {
    return 0;
  }
}

int32_t char_cmp_deref(const void* aa, const void* bb) {
  return (int32_t)(**((const char**)aa) - **((const char**)bb));
}

int32_t intcmp(const void* aa, const void* bb) {
  return *((const int32_t*)aa) - *((const int32_t*)bb);
}

int32_t uintcmp(const void* aa, const void* bb) {
  if (*((const uint32_t*)aa) < *((const uint32_t*)bb)) {
    return -1;
  } else {
    return (*((const uint32_t*)aa) > *((const uint32_t*)bb));
  }
}

int32_t intcmp2(const void* aa, const void* bb) {
  if (*((const int32_t*)aa) < *((const int32_t*)bb)) {
    return -1;
  } else {
    return (*((const int32_t*)aa) > *((const int32_t*)bb));
  }
}

int32_t intcmp3_decr(const void* aa, const void* bb) {
  int32_t ii = *((const int32_t*)bb) - *((const int32_t*)aa);
  if (ii) {
    return ii;
  }
  ii = ((const int32_t*)bb)[1] - ((const int32_t*)aa)[1];
  if (ii) {
    return ii;
  }
  return ((const int32_t*)bb)[2] - ((const int32_t*)aa)[2];
}

#ifndef __cplusplus
int32_t llcmp(const void* aa, const void* bb) {
  int64_t diff = *((const int64_t*)aa) - *((const int64_t*)bb);
  if (diff > 0) {
    return 1;
  } else if (diff < 0) {
    return -1;
  } else {
    return 0;
  }
}
#endif

// alas, qsort_r not available on some Linux distributions

// Normally use qsort_ext(), but this version is necessary before g_bigstack
// has been allocated.
void qsort_ext2(char* main_arr, uintptr_t arr_length, uintptr_t item_length, int(* comparator_deref)(const void*, const void*), char* secondary_arr, uintptr_t secondary_item_len, char* proxy_arr, uintptr_t proxy_len) {
  uintptr_t ulii;
  for (ulii = 0; ulii < arr_length; ulii++) {
    *(char**)(&(proxy_arr[ulii * proxy_len])) = &(main_arr[ulii * item_length]);
    memcpy(&(proxy_arr[ulii * proxy_len + sizeof(void*)]), &(secondary_arr[ulii * secondary_item_len]), secondary_item_len);
  }
  qsort(proxy_arr, arr_length, proxy_len, comparator_deref);
  for (ulii = 0; ulii < arr_length; ulii++) {
    memcpy(&(secondary_arr[ulii * secondary_item_len]), &(proxy_arr[ulii * proxy_len + sizeof(void*)]), secondary_item_len);
    memcpy(&(proxy_arr[ulii * proxy_len]), *(char**)(&(proxy_arr[ulii * proxy_len])), item_length);
  }
  for (ulii = 0; ulii < arr_length; ulii++) {
    memcpy(&(main_arr[ulii * item_length]), &(proxy_arr[ulii * proxy_len]), item_length);
  }
}

// This actually tends to be faster than just sorting an array of indices,
// because of memory locality issues.
int32_t qsort_ext(char* main_arr, uintptr_t arr_length, uintptr_t item_length, int(* comparator_deref)(const void*, const void*), char* secondary_arr, intptr_t secondary_item_len) {
  // main_arr = packed array of equal-length items to sort
  // arr_length = number of items
  // item_length = byte count of each main_arr item
  // comparator_deref = returns positive if *first > *second, 0 if equal,
  //                    negative if *first < *second.  Note the extra
  //                    dereference.
  // secondary_arr = packed array of fixed-length records associated with the
  //                 main_arr items, to be resorted in the same way.  (e.g.
  //                 if one is building an index, this could start as a sorted
  //                 0..(n-1) sequence of integers; then, post-sort, this would
  //                 be a lookup table for the original position of each
  //                 main_arr item.)
  // secondary_item_len = byte count of each secondary_arr item
  uintptr_t proxy_len = secondary_item_len + sizeof(void*);
  unsigned char* bigstack_mark = g_bigstack_base;
  char* proxy_arr;
  if (!arr_length) {
    return 0;
  }
  if (proxy_len < item_length) {
    proxy_len = item_length;
  }
  if (bigstack_alloc_c(arr_length * proxy_len, &proxy_arr)) {
    return -1;
  }
  qsort_ext2(main_arr, arr_length, item_length, comparator_deref, secondary_arr, secondary_item_len, proxy_arr, proxy_len);
  bigstack_reset(bigstack_mark);
  return 0;
}

int32_t sort_item_ids_noalloc(uintptr_t unfiltered_ct, const uintptr_t* exclude_arr, uintptr_t item_ct, const char* __restrict item_ids, uintptr_t max_id_len, uint32_t allow_dups, uint32_t collapse_idxs, int(* comparator_deref)(const void*, const void*), char* __restrict sorted_ids, uint32_t* id_map) {
  // Stores a lexicographically sorted list of IDs in sorted_ids and the raw
  // positions of the corresponding markers/samples in *id_map_ptr.  Does not
  // include excluded markers/samples in the list.
  // Assumes sorted_ids and id_map have been allocated; use the sort_item_ids()
  // wrapper if they haven't been.
  // Note that this DOES still perform a "stack" allocation (in the qsort_ext()
  // call).
  uint32_t uii = 0;
  char* dup_id;
  char* tptr;
  uint32_t ujj;
  if (!item_ct) {
    return 0;
  }
  if (!collapse_idxs) {
    for (ujj = 0; ujj < item_ct; uii++, ujj++) {
      next_unset_unsafe_ck(exclude_arr, &uii);
      memcpy(&(sorted_ids[ujj * max_id_len]), &(item_ids[uii * max_id_len]), max_id_len);
      id_map[ujj] = uii;
    }
  } else {
    for (ujj = 0; ujj < item_ct; uii++, ujj++) {
      next_unset_unsafe_ck(exclude_arr, &uii);
      memcpy(&(sorted_ids[ujj * max_id_len]), &(item_ids[uii * max_id_len]), max_id_len);
      id_map[ujj] = ujj;
    }
  }
  if (qsort_ext(sorted_ids, item_ct, max_id_len, comparator_deref, (char*)id_map, sizeof(int32_t))) {
    return RET_NOMEM;
  }
  if (!allow_dups) {
    dup_id = scan_for_duplicate_ids(sorted_ids, item_ct, max_id_len);
    if (dup_id) {
      tptr = strchr(dup_id, '\t');
      if (tptr) {
        *tptr = ' ';
      }
      LOGERRPRINTFWW("Error: Duplicate ID '%s'.\n", dup_id);
      return RET_INVALID_FORMAT;
    }
  }
  return 0;
}

int32_t sort_item_ids(uintptr_t unfiltered_ct, const uintptr_t* exclude_arr, uintptr_t exclude_ct, const char* __restrict item_ids, uintptr_t max_id_len, uint32_t allow_dups, uint32_t collapse_idxs, int(* comparator_deref)(const void*, const void*), char** sorted_ids_ptr, uint32_t** id_map_ptr) {
  uintptr_t item_ct = unfiltered_ct - exclude_ct;
  // id_map on bottom because --indiv-sort frees *sorted_ids_ptr
  if (bigstack_alloc_ui(item_ct, id_map_ptr) ||
      bigstack_alloc_c(item_ct * max_id_len, sorted_ids_ptr)) {
    return RET_NOMEM;
  }
  return sort_item_ids_noalloc(unfiltered_ct, exclude_arr, item_ct, item_ids, max_id_len, allow_dups, collapse_idxs, comparator_deref, *sorted_ids_ptr, *id_map_ptr);
}

uint32_t uint32arr_greater_than(const uint32_t* sorted_uint32_arr, uint32_t arr_length, uint32_t uii) {
  // assumes arr_length is nonzero, and sorted_uint32_arr is in nondecreasing
  // order.  (useful for searching marker_pos.)
  // uii guaranteed to be larger than sorted_uint32_arr[min_idx - 1] if it
  // exists, but NOT necessarily sorted_uint32_arr[min_idx].
  int32_t min_idx = 0;
  // similarly, uii guaranteed to be no greater than
  // sorted_uint32_arr[max_idx + 1] if it exists, but not necessarily
  // sorted_uint32_arr[max_idx].  Signed integer since it could become -1.
  int32_t max_idx = arr_length - 1;
  uint32_t mid_idx;
  while (min_idx < max_idx) {
    mid_idx = (((uint32_t)min_idx) + ((uint32_t)max_idx)) / 2;
    if (uii > sorted_uint32_arr[mid_idx]) {
      min_idx = mid_idx + 1;
    } else {
      max_idx = mid_idx - 1;
    }
  }
  if (uii > sorted_uint32_arr[((uint32_t)min_idx)]) {
    return (min_idx + 1);
  } else {
    return min_idx;
  }
}

uint32_t int32arr_greater_than(const int32_t* sorted_int32_arr, uint32_t arr_length, int32_t ii) {
  int32_t min_idx = 0;
  int32_t max_idx = arr_length - 1;
  uint32_t mid_idx;
  while (min_idx < max_idx) {
    mid_idx = (((uint32_t)min_idx) + ((uint32_t)max_idx)) / 2;
    if (ii > sorted_int32_arr[mid_idx]) {
      min_idx = mid_idx + 1;
    } else {
      max_idx = mid_idx - 1;
    }
  }
  if (ii > sorted_int32_arr[((uint32_t)min_idx)]) {
    return (min_idx + 1);
  } else {
    return min_idx;
  }
}

uintptr_t uint64arr_greater_than(const uint64_t* sorted_uint64_arr, uintptr_t arr_length, uint64_t ullii) {
  intptr_t min_idx = 0;
  intptr_t max_idx = arr_length - 1;
  uintptr_t mid_idx;
  while (min_idx < max_idx) {
    mid_idx = (((uintptr_t)min_idx) + ((uintptr_t)max_idx)) / 2;
    if (ullii > sorted_uint64_arr[mid_idx]) {
      min_idx = mid_idx + 1;
    } else {
      max_idx = mid_idx - 1;
    }
  }
  if (ullii > sorted_uint64_arr[((uintptr_t)min_idx)]) {
    return (min_idx + 1);
  } else {
    return min_idx;
  }
}

uintptr_t doublearr_greater_than(const double* sorted_dbl_arr, uintptr_t arr_length, double dxx) {
  // returns number of items in sorted_dbl_arr which dxx is greater than.
  // assumes array is nonempty and sorted in nondecreasing order
  intptr_t min_idx = 0;
  intptr_t max_idx = arr_length - 1;
  uintptr_t mid_idx;
  while (min_idx < max_idx) {
    mid_idx = (((uintptr_t)min_idx) + ((uintptr_t)max_idx)) / 2;
    if (dxx > sorted_dbl_arr[mid_idx]) {
      min_idx = mid_idx + 1;
    } else {
      max_idx = mid_idx - 1;
    }
  }
  if (dxx > sorted_dbl_arr[((uintptr_t)min_idx)]) {
    return (min_idx + 1);
  } else {
    return min_idx;
  }
}

uintptr_t nonincr_doublearr_leq_stride(const double* nonincr_dbl_arr, uintptr_t arr_length, uintptr_t stride, double dxx) {
  // assumes relevant elements of array are sorted in nonincreasing order
  // instead, and they are spaced stride units apart
  intptr_t min_idx = 0;
  intptr_t max_idx = arr_length - 1;
  uintptr_t mid_idx;
  while (min_idx < max_idx) {
    mid_idx = (((uintptr_t)min_idx) + ((uintptr_t)max_idx)) / 2;
    if (dxx <= nonincr_dbl_arr[mid_idx * stride]) {
      min_idx = mid_idx + 1;
    } else {
      max_idx = mid_idx - 1;
    }
  }
  if (dxx <= nonincr_dbl_arr[((uintptr_t)min_idx) * stride]) {
    return (min_idx + 1);
  } else {
    return min_idx;
  }
}

int32_t bsearch_str(const char* id_buf, uintptr_t cur_id_len, const char* lptr, uintptr_t max_id_len, uintptr_t end_idx) {
  // does not assume null-terminated id_buf, or nonempty array.
  // N.B. max_id_len includes null terminator as usual, while cur_id_len does
  // NOT.
  uintptr_t start_idx = 0;
  uintptr_t mid_idx;
  int32_t ii;
  if (cur_id_len >= max_id_len) {
    return -1;
  }
  while (start_idx < end_idx) {
    mid_idx = (start_idx + end_idx) / 2;
    ii = memcmp(id_buf, &(lptr[mid_idx * max_id_len]), cur_id_len);
    if (ii > 0) {
      start_idx = mid_idx + 1;
    } else if ((ii < 0) || lptr[mid_idx * max_id_len + cur_id_len]) {
      end_idx = mid_idx;
    } else {
      return ((uint32_t)mid_idx);
    }
  }
  return -1;
}

int32_t bsearch_str_natural(const char* id_buf, const char* lptr, uintptr_t max_id_len, uintptr_t end_idx) {
  // unlike bsearch_str(), caller is responsible for slen > max_id_len check
  // if appropriate here
  uintptr_t start_idx = 0;
  uintptr_t mid_idx;
  int32_t ii;
  while (start_idx < end_idx) {
    mid_idx = (start_idx + end_idx) / 2;
    ii = strcmp_natural(id_buf, &(lptr[mid_idx * max_id_len]));
    if (ii > 0) {
      start_idx = mid_idx + 1;
    } else if (ii < 0) {
      end_idx = mid_idx;
    } else {
      return ((uint32_t)mid_idx);
    }
  }
  return -1;
}

uintptr_t bsearch_str_lb(const char* id_buf, uintptr_t cur_id_len, const char* lptr, uintptr_t max_id_len, uintptr_t end_idx) {
  // returns number of elements in lptr[] less than id_buf.
  uintptr_t start_idx = 0;
  uintptr_t mid_idx;
  if (cur_id_len > max_id_len) {
    cur_id_len = max_id_len;
  }
  while (start_idx < end_idx) {
    mid_idx = (start_idx + end_idx) / 2;
    if (memcmp(id_buf, &(lptr[mid_idx * max_id_len]), cur_id_len) > 0) {
      start_idx = mid_idx + 1;
    } else {
      end_idx = mid_idx;
    }
  }
  return start_idx;
}

uint32_t bsearch_read_fam_indiv(char* __restrict read_ptr, const char* __restrict lptr, uintptr_t max_id_len, uintptr_t filter_line_ct, char** read_pp_new, int32_t* retval_ptr, char* __restrict id_buf) {
  // id_buf = workspace
  // lptr = packed, sorted list of ID strings to search over
  // read_ptr is assumed to point to beginning of FID.  FID is terminated by
  // any space/eoln character, then IID is assumed to follow it (and is also
  // terminated by any space/eoln).  Nonzero error value is returned if IID
  // does not exist.
  char* iid_ptr;
  uintptr_t slen_fid;
  uintptr_t slen_iid;
  uintptr_t slen_final;
  slen_fid = strlen_se(read_ptr);
  iid_ptr = skip_initial_spaces(&(read_ptr[slen_fid]));
  if (is_eoln_kns(*iid_ptr)) {
    return 1;
  }
  slen_iid = strlen_se(iid_ptr);
  if (read_pp_new) {
    *read_pp_new = skip_initial_spaces(&(iid_ptr[slen_iid]));
  }
  slen_final = slen_fid + slen_iid + 1;
  if (slen_final >= max_id_len) {
    // avoid buffer overflow
    *retval_ptr = -1;
    return 0;
  }
  // error message bugfix: null-terminate this string
  memcpyx(memcpyax(id_buf, read_ptr, slen_fid, '\t'), iid_ptr, slen_iid, '\0');
  *retval_ptr = bsearch_str(id_buf, slen_final, lptr, max_id_len, filter_line_ct);
  return 0;
}

void bsearch_fam(const char* __restrict fam_id, const char* __restrict lptr, uintptr_t max_id_len, uint32_t filter_line_ct, uint32_t* __restrict first_idx_ptr, uint32_t* __restrict last_idx_ptr, char* __restrict id_buf) {
  uint32_t slen;
  uint32_t fidx;
  uint32_t loff;
  if (!filter_line_ct) {
    goto bsearch_fam_ret_null;
  }
  slen = strlen_se(fam_id);
  if (slen + 3 > max_id_len) {
    goto bsearch_fam_ret_null;
  }
  memcpy(id_buf, fam_id, slen);
  id_buf[slen] = '\t';
  fidx = bsearch_str_lb(id_buf, slen + 1, lptr, max_id_len, filter_line_ct);
  if (fidx == filter_line_ct) {
    goto bsearch_fam_ret_null;
  }
  id_buf[slen] = ' ';
  loff = bsearch_str_lb(id_buf, slen + 1, &(lptr[fidx * max_id_len]), max_id_len, filter_line_ct - fidx);
  if (!loff) {
    goto bsearch_fam_ret_null;
  }
  *first_idx_ptr = fidx;
  *last_idx_ptr = fidx + loff;
  return;
 bsearch_fam_ret_null:
  *first_idx_ptr = 0;
  *last_idx_ptr = 0;
}

void bitarr_invert(uintptr_t bit_ct, uintptr_t* bitarr) {
  uintptr_t* bitarr_stop = &(bitarr[bit_ct / BITCT]);
  while (bitarr < bitarr_stop) {
    *bitarr = ~(*bitarr);
    bitarr++;
  }
  if (bit_ct % BITCT) {
    *bitarr = (~(*bitarr)) & ((ONELU << (bit_ct % BITCT)) - ONELU);
  }
}

void bitarr_invert_copy(const uintptr_t* input_bitarr, uintptr_t bit_ct, uintptr_t* output_bitarr) {
  const uintptr_t* input_stop = &(input_bitarr[bit_ct / BITCT]);
  while (input_bitarr < input_stop) {
    *output_bitarr++ = ~(*input_bitarr++);
  }
  if (bit_ct % BITCT) {
    *output_bitarr = (~(*input_bitarr)) & ((ONELU << (bit_ct % BITCT)) - ONELU);
  }
}

void bitvec_and(const uintptr_t* __restrict arg_bitvec, uintptr_t word_ct, uintptr_t* __restrict main_bitvec) {
  // main_bitvec := main_bitvec AND arg_bitvec
#ifdef __LP64__
  __m128i* vv128 = (__m128i*)main_bitvec;
  const __m128i* iv128 = (const __m128i*)arg_bitvec;
  __m128i* vv128_end = &(vv128[word_ct / 2]);
  while (vv128 < vv128_end) {
    *vv128 = _mm_and_si128(*iv128++, *vv128);
    vv128++;
  }
  if (word_ct & 1) {
    word_ct--;
    main_bitvec[word_ct] &= arg_bitvec[word_ct];
  }
#else
  uintptr_t* bitvec_end = &(main_bitvec[word_ct]);
  while (main_bitvec < bitvec_end) {
    *main_bitvec++ &= *arg_bitvec++;
  }
#endif
}

void bitvec_andnot(const uintptr_t* __restrict exclude_bitvec, uintptr_t word_ct, uintptr_t* __restrict main_bitvec) {
  // main_bitvec := main_bitvec ANDNOT exclude_bitvec
  // note that this is the reverse of the _mm_andnot() operand order
#ifdef __LP64__
  __m128i* vv128 = (__m128i*)main_bitvec;
  const __m128i* ev128 = (const __m128i*)exclude_bitvec;
  __m128i* vv128_end = &(vv128[word_ct / 2]);
  while (vv128 < vv128_end) {
    *vv128 = _mm_andnot_si128(*ev128++, *vv128);
    vv128++;
  }
  if (word_ct & 1) {
    word_ct--;
    main_bitvec[word_ct] &= ~(exclude_bitvec[word_ct]);
  }
#else
  uintptr_t* bitvec_end = &(main_bitvec[word_ct]);
  while (main_bitvec < bitvec_end) {
    *main_bitvec++ &= ~(*exclude_bitvec++);
  }
#endif
}

void bitvec_andnot_reversed_args(const uintptr_t* __restrict include_bitvec, uintptr_t word_ct, uintptr_t* __restrict main_bitvec) {
  // main_bitvec := (~main_bitvec) AND include_bitvec
#ifdef __LP64__
  __m128i* vv128 = (__m128i*)main_bitvec;
  const __m128i* iv128 = (const __m128i*)include_bitvec;
  __m128i* vv128_end = &(vv128[word_ct / 2]);
  while (vv128 < vv128_end) {
    *vv128 = _mm_andnot_si128(*vv128, *iv128++);
    vv128++;
  }
  if (word_ct & 1) {
    word_ct--;
    main_bitvec[word_ct] = (~main_bitvec[word_ct]) & include_bitvec[word_ct];
  }
#else
  uintptr_t* bitvec_end = &(main_bitvec[word_ct]);
  while (main_bitvec < bitvec_end) {
    *main_bitvec = (~(*main_bitvec)) & (*include_bitvec++);
    main_bitvec++;
  }
#endif
}

void bitvec_or(const uintptr_t* __restrict arg_bitvec, uintptr_t word_ct, uintptr_t* main_bitvec) {
  // main_bitvec := main_bitvec OR arg_bitvec
#ifdef __LP64__
  __m128i* vv128 = (__m128i*)main_bitvec;
  const __m128i* ov128 = (const __m128i*)arg_bitvec;
  __m128i* vv128_end = &(vv128[word_ct / 2]);
  while (vv128 < vv128_end) {
    *vv128 = _mm_or_si128(*ov128++, *vv128);
    vv128++;
  }
  if (word_ct & 1) {
    word_ct--;
    main_bitvec[word_ct] |= arg_bitvec[word_ct];
  }
#else
  uintptr_t* vec_end = &(main_bitvec[word_ct]);
  while (main_bitvec < vec_end) {
    *main_bitvec++ |= *arg_bitvec++;
  }
#endif
}

void bitvec_ornot(const uintptr_t* __restrict inverted_or_bitvec, uintptr_t word_ct, uintptr_t* __restrict main_bitvec) {
  // main_bitvec := main_bitvec OR (~inverted_or_bitvec)
#ifdef __LP64__
#ifdef __APPLE__
  const __m128i all1 = {0xffffffffffffffffLLU, 0xffffffffffffffffLLU};
#else
  const __m128i all1 = {-1LL, -1LL};
#endif
  __m128i* vv128 = (__m128i*)main_bitvec;
  const __m128i* ev128 = (const __m128i*)inverted_or_bitvec;
  __m128i* vv128_end = &(vv128[word_ct / 2]);
  while (vv128 < vv128_end) {
    *vv128 = _mm_or_si128(_mm_xor_si128(*ev128++, all1), *vv128);
    vv128++;
  }
  if (word_ct & 1) {
    word_ct--;
    main_bitvec[word_ct] |= ~(inverted_or_bitvec[word_ct]);
  }
#else
  uintptr_t* vec_end = &(main_bitvec[word_ct]);
  while (main_bitvec < vec_end) {
    *main_bitvec++ |= ~(*inverted_or_bitvec++);
  }
#endif
}

void bitvec_xor(const uintptr_t* __restrict arg_bitvec, uintptr_t word_ct, uintptr_t* __restrict main_bitvec) {
  // main_bitvec := main_bitvec XOR xor_bitvec
#ifdef __LP64__
  __m128i* bitv128 = (__m128i*)main_bitvec;
  __m128i* xorv128 = (__m128i*)arg_bitvec;
  __m128i* bitv128_end = &(bitv128[word_ct / 2]);
  while (bitv128 < bitv128_end) {
    *bitv128 = _mm_xor_si128(*xorv128++, *bitv128);
    bitv128++;
  }
  if (word_ct & 1) {
    word_ct--;
    main_bitvec[word_ct] ^= arg_bitvec[word_ct];
  }
#else
  uintptr_t* main_bitvec_end = &(main_bitvec[word_ct]);
  while (main_bitvec < main_bitvec_end) {
    *main_bitvec++ ^= *arg_bitvec++;
  }
#endif
}

uint32_t is_monomorphic_a2(const uintptr_t* geno_arr, uint32_t sample_ct) {
  const uintptr_t* loop_end = &(geno_arr[sample_ct / BITCT2]);
  uint32_t sample_rem = sample_ct % BITCT2;
  for (; geno_arr < loop_end; geno_arr++) {
    if ((~(*geno_arr)) & FIVEMASK) {
      return 0;
    }
  }
  return (sample_rem && ((~(*geno_arr)) & (FIVEMASK >> (BITCT - sample_rem * 2))))? 0 : 1;
}

uint32_t is_monomorphic(const uintptr_t* geno_arr, uint32_t sample_ct) {
  uint32_t sample_ctd2 = sample_ct / BITCT2;
  uint32_t sample_rem = sample_ct % BITCT2;
  uintptr_t ulii;
  uintptr_t uljj;
  while (sample_ctd2) {
    ulii = *geno_arr++;
    uljj = (ulii >> 1) & FIVEMASK;
    ulii = ~ulii;
    // now ulii & FIVEMASK = low bit zero, uljj = high bit one
    if (uljj) {
      if (uljj & ulii) {
        // heterozygote observed
        return 0;
      }
      // homozyg A2 observed
      while (1) {
	// 00 and 10 now both demonstrate marker is polymorphic
	if (ulii & FIVEMASK) {
	  return 0;
	}
	if (!(--sample_ctd2)) {
	  return (sample_rem && ((~(*geno_arr)) & (FIVEMASK >> (BITCT - sample_rem * 2))))? 0 : 1;
	}
	ulii = ~(*geno_arr++);
      }
    } else if (ulii & FIVEMASK) {
      do {
        if (!(--sample_ctd2)) {
          return (sample_rem && ((*geno_arr) & (AAAAMASK >> (BITCT - sample_rem * 2))))? 0 : 1;
	}
	ulii = *geno_arr++;
      } while (!(ulii & AAAAMASK));
      return 0;
    }
    sample_ctd2--;
  }
  if (sample_rem) {
    ulii = *geno_arr;
    uljj = (ulii >> 1) & FIVEMASK;
    ulii = ~ulii;
    if ((uljj & ulii) || (uljj && (ulii & (~uljj) & (FIVEMASK >> (BITCT - sample_rem * 2))))) {
      return 0;
    }
  }
  return 1;
}

uint32_t less_than_two_genotypes(const uintptr_t* geno_arr, uint32_t sample_ct) {
  uint32_t sample_ctd2 = sample_ct / BITCT2;
  uint32_t sample_rem = sample_ct % BITCT2;
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t ulkk;
  uint32_t distinct_genotype_ct;
  while (sample_ctd2) {
    ulii = *geno_arr++;
    uljj = (ulii >> 1) & FIVEMASK;
    ulkk = ~ulii;
    if (uljj) {
      if (uljj & ulii) {
	// homozygote major observed; either 00 or 10 now demonstrate marker
	// is polymorphic
	while (1) {
	  if (ulkk & FIVEMASK) {
	    return 0;
	  }
	  if (!(--sample_ctd2)) {
	    return (sample_rem && ((~(*geno_arr)) & (FIVEMASK >> (BITCT - sample_rem * 2))))? 0 : 1;
	  }
	  ulkk = ~(*geno_arr++);
	}
      } else {
        // heterozygote observed; either 00 or 11 now means we have 2+
	// genotypes
	while (1) {
	  ulii = ~(*geno_arr++);
	  if (!(--sample_ctd2)) {
	    return (sample_rem && (((~ulii) ^ (ulii >> 1)) & (FIVEMASK >> (BITCT - sample_rem * 2))))? 0 : 1;
	  }
	  if (((~ulii) ^ (ulii >> 1)) & FIVEMASK) {
	    return 0;
	  }
	}
      }
    } else if (ulkk & FIVEMASK) {
      // homozygous minor observed; either 10 or 11 now demonstrate marker is
      // polymorphic
      do {
        if (!(--sample_ctd2)) {
          return (sample_rem && ((*geno_arr) & (AAAAMASK >> (BITCT - sample_rem * 2))))? 0 : 1;
	}
	ulii = *geno_arr++;
      } while (!(ulii & AAAAMASK));
      return 0;
    }
    sample_ctd2--;
  }
  if (sample_rem) {
    ulii = *geno_arr;
    uljj = (ulii >> 1) & FIVEMASK;
    ulkk = ~ulii;
    // homozygous minor present?
    distinct_genotype_ct = (ulkk & (~uljj) & (FIVEMASK >> (BITCT - sample_rem * 2)))? 1 : 0;
    // heterozygous present?
    distinct_genotype_ct += (uljj & ulkk)? 1 : 0;
    // homozygous major present?
    distinct_genotype_ct += (uljj & ulii)? 1 : 0;
    if (distinct_genotype_ct > 1) {
      return 0;
    }
  }
  return 1;
}

/*
uint32_t has_three_genotypes(uintptr_t* lptr, uint32_t sample_ct) {
  uintptr_t* lptr_end = &(lptr[sample_ct / BITCT2]);
  uint32_t sample_rem = sample_ct % BITCT2;
  uintptr_t* cur_lptr;
  uintptr_t ulii;
  uintptr_t uljj;
  cur_lptr = lptr;
  while (1) {
    ulii = ~(*cur_lptr);
    uljj = ulii & (ulii >> 1) & FIVEMASK;
    if (cur_lptr == lptr_end) {
      if ((!sample_rem) || (!(uljj << (BITCT - sample_rem * 2)))) {
	return 0;
      }
      break;
    }
    if (uljj) {
      // found hom A1
      break;
    }
    cur_lptr++;
  }
  cur_lptr = lptr;
  // zero-padding is benign for het and hom A2 checks
  lptr_end = &(lptr[QUATERCT_TO_WORDCT(sample_ct)]);
  while (1) {
    ulii = *cur_lptr;
    uljj = (ulii >> 1) & FIVEMASK;
    if ((~ulii) & uljj) {
      break;
    }
    if (++cur_lptr == lptr_end) {
      return 0;
    }
  }
  cur_lptr = lptr;
  do {
    ulii = *cur_lptr;
    uljj = (ulii >> 1) & FIVEMASK;
    if (ulii & uljj) {
      return 1;
    }
  } while (++cur_lptr < lptr_end);
  return 0;
}
*/

#ifdef __LP64__
// Basic SSE2 implementation of Lauradoux/Walisch popcount.
static inline uintptr_t popcount_vecs(const __m128i* vptr, uintptr_t ct) {
  // popcounts vptr[0..(ct-1)].  Assumes ct is a multiple of 3 (0 ok).
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  uintptr_t tot = 0;
  const __m128i* vend;
  __m128i count1;
  __m128i count2;
  __m128i half1;
  __m128i half2;
  __univec acc;

  while (ct >= 30) {
    ct -= 30;
    vend = &(vptr[30]);
  popcount_vecs_main_loop:
    acc.vi = _mm_setzero_si128();
    do {
      count1 = *vptr++;
      count2 = *vptr++;
      half1 = *vptr++;
      half2 = _mm_and_si128(_mm_srli_epi64(half1, 1), m1);
      half1 = _mm_and_si128(half1, m1);
      // Two bits can represent values from 0-3, so make each pair in count1
      // count2 store a partial bitcount covering themselves AND another bit
      // from elsewhere.
      count1 = _mm_sub_epi64(count1, _mm_and_si128(_mm_srli_epi64(count1, 1), m1));
      count2 = _mm_sub_epi64(count2, _mm_and_si128(_mm_srli_epi64(count2, 1), m1));
      count1 = _mm_add_epi64(count1, half1);
      count2 = _mm_add_epi64(count2, half2);
      // Four bits represent 0-15, so we can safely add four 0-3 partial
      // bitcounts together.
      count1 = _mm_add_epi64(_mm_and_si128(count1, m2), _mm_and_si128(_mm_srli_epi64(count1, 2), m2));
      count1 = _mm_add_epi64(count1, _mm_add_epi64(_mm_and_si128(count2, m2), _mm_and_si128(_mm_srli_epi64(count2, 2), m2)));
      // Accumulator stores sixteen 0-255 counts in parallel.
      acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(count1, m4), _mm_and_si128(_mm_srli_epi64(count1, 4), m4)));
    } while (vptr < vend);
    acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
    tot += ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
  }
  if (ct) {
    vend = &(vptr[ct]);
    ct = 0;
    goto popcount_vecs_main_loop;
  }
  return tot;
}

static inline uintptr_t popcount2_vecs(const __m128i* vptr, uintptr_t ct) {
  // assumes ct is a multiple of 6.
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  uintptr_t tot = 0;
  const __m128i* vend;
  __m128i loader1;
  __m128i loader2;
  __m128i count1;
  __m128i count2;
  __univec acc;

  while (ct >= 30) {
    ct -= 30;
    vend = &(vptr[30]);
  popcount2_vecs_main_loop:
    acc.vi = _mm_setzero_si128();
    do {
      loader1 = *vptr++;
      loader2 = *vptr++;
      count1 = _mm_add_epi64(_mm_and_si128(loader1, m2), _mm_and_si128(_mm_srli_epi64(loader1, 2), m2));
      count2 = _mm_add_epi64(_mm_and_si128(loader2, m2), _mm_and_si128(_mm_srli_epi64(loader2, 2), m2));

      loader1 = *vptr++;
      loader2 = *vptr++;
      count1 = _mm_add_epi64(count1, _mm_add_epi64(_mm_and_si128(loader1, m2), _mm_and_si128(_mm_srli_epi64(loader1, 2), m2)));
      count2 = _mm_add_epi64(count2, _mm_add_epi64(_mm_and_si128(loader2, m2), _mm_and_si128(_mm_srli_epi64(loader2, 2), m2)));

      loader1 = *vptr++;
      loader2 = *vptr++;
      count1 = _mm_add_epi64(count1, _mm_add_epi64(_mm_and_si128(loader1, m2), _mm_and_si128(_mm_srli_epi64(loader1, 2), m2)));
      count2 = _mm_add_epi64(count2, _mm_add_epi64(_mm_and_si128(loader2, m2), _mm_and_si128(_mm_srli_epi64(loader2, 2), m2)));

      acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(count1, m4), _mm_and_si128(_mm_srli_epi64(count1, 4), m4)));
      acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(count2, m4), _mm_and_si128(_mm_srli_epi64(count2, 4), m4)));
    } while (vptr < vend);
    acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
    tot += ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
  }
  if (ct) {
    vend = &(vptr[ct]);
    ct = 0;
    goto popcount2_vecs_main_loop;
  }
  return tot;
}

static inline uintptr_t popcount_vecs_exclude(const __m128i* __restrict vptr, const __m128i* __restrict exclude_ptr, uintptr_t ct) {
  // popcounts vptr ANDNOT exclude_ptr[0..(ct-1)].  ct is a multiple of 3.
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  uintptr_t tot = 0;
  const __m128i* vend;
  __m128i count1, count2, half1, half2;
  __univec acc;

  while (ct >= 30) {
    ct -= 30;
    vend = &(vptr[30]);
  popcount_vecs_exclude_main_loop:
    acc.vi = _mm_setzero_si128();
    do {
      // nots the FIRST value
      count1 = _mm_andnot_si128(*exclude_ptr++, *vptr++);
      count2 = _mm_andnot_si128(*exclude_ptr++, *vptr++);
      half1 = _mm_andnot_si128(*exclude_ptr++, *vptr++);
      half2 = _mm_and_si128(_mm_srli_epi64(half1, 1), m1);
      half1 = _mm_and_si128(half1, m1);
      count1 = _mm_sub_epi64(count1, _mm_and_si128(_mm_srli_epi64(count1, 1), m1));
      count2 = _mm_sub_epi64(count2, _mm_and_si128(_mm_srli_epi64(count2, 1), m1));
      count1 = _mm_add_epi64(count1, half1);
      count2 = _mm_add_epi64(count2, half2);
      count1 = _mm_add_epi64(_mm_and_si128(count1, m2), _mm_and_si128(_mm_srli_epi64(count1, 2), m2));
      count1 = _mm_add_epi64(count1, _mm_add_epi64(_mm_and_si128(count2, m2), _mm_and_si128(_mm_srli_epi64(count2, 2), m2)));
      acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(count1, m4), _mm_and_si128(_mm_srli_epi64(count1, 4), m4)));
    } while (vptr < vend);
    acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
    tot += ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
  }
  if (ct) {
    vend = &(vptr[ct]);
    ct = 0;
    goto popcount_vecs_exclude_main_loop;
  }
  return tot;
}

static inline uintptr_t popcount_vecs_intersect(const __m128i* __restrict vptr1, const __m128i* __restrict vptr2, uintptr_t ct) {
  // popcounts vptr1 AND vptr2[0..(ct-1)].  ct is a multiple of 3.
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  uintptr_t tot = 0;
  const __m128i* vend1;
  __m128i count1, count2, half1, half2;
  __univec acc;

  while (ct >= 30) {
    ct -= 30;
    vend1 = &(vptr1[30]);
  popcount_vecs_intersect_main_loop:
    acc.vi = _mm_setzero_si128();
    do {
      count1 = _mm_and_si128(*vptr2++, *vptr1++);
      count2 = _mm_and_si128(*vptr2++, *vptr1++);
      half1 = _mm_and_si128(*vptr2++, *vptr1++);
      half2 = _mm_and_si128(_mm_srli_epi64(half1, 1), m1);
      half1 = _mm_and_si128(half1, m1);
      count1 = _mm_sub_epi64(count1, _mm_and_si128(_mm_srli_epi64(count1, 1), m1));
      count2 = _mm_sub_epi64(count2, _mm_and_si128(_mm_srli_epi64(count2, 1), m1));
      count1 = _mm_add_epi64(count1, half1);
      count2 = _mm_add_epi64(count2, half2);
      count1 = _mm_add_epi64(_mm_and_si128(count1, m2), _mm_and_si128(_mm_srli_epi64(count1, 2), m2));
      count1 = _mm_add_epi64(count1, _mm_add_epi64(_mm_and_si128(count2, m2), _mm_and_si128(_mm_srli_epi64(count2, 2), m2)));
      acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(count1, m4), _mm_and_si128(_mm_srli_epi64(count1, 4), m4)));
    } while (vptr1 < vend1);
    acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
    tot += ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
  }
  if (ct) {
    vend1 = &(vptr1[ct]);
    ct = 0;
    goto popcount_vecs_intersect_main_loop;
  }
  return tot;
}
#endif

uintptr_t popcount_longs(const uintptr_t* lptr, uintptr_t word_ct) {
  // Efficiently popcounts lptr[0..(word_ct - 1)].  In the 64-bit case, lptr[]
  // must be 16-byte aligned.
  // The popcount_longs_nzbase() wrapper takes care of starting from a later
  // index.
  uintptr_t tot = 0;
  const uintptr_t* lptr_end = &(lptr[word_ct]);
#ifdef __LP64__
  uintptr_t six_ct;
  const __m128i* vptr;
  vptr = (const __m128i*)lptr;
  six_ct = word_ct / 6;
  tot += popcount_vecs(vptr, six_ct * 3);
  lptr = &(lptr[six_ct * 6]);
#else
  // The humble 16-bit lookup table actually beats
  // http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel
  // on my development machine by a hair.
  // However, if we take the hint from Lauradoux/Walisch and postpone the
  // multiply and right shift, this is no longer true.  Ah well.
  const uintptr_t* lptr_six_end;
  uintptr_t tmp_stor;
  uintptr_t loader;
  uintptr_t ulii;
  uintptr_t uljj;
  lptr_six_end = &(lptr[word_ct - (word_ct % 6)]);
  while (lptr < lptr_six_end) {
    loader = *lptr++;
    ulii = loader - ((loader >> 1) & FIVEMASK);
    loader = *lptr++;
    uljj = loader - ((loader >> 1) & FIVEMASK);
    loader = *lptr++;
    ulii += (loader >> 1) & FIVEMASK;
    uljj += loader & FIVEMASK;
    ulii = (ulii & 0x33333333) + ((ulii >> 2) & 0x33333333);
    ulii += (uljj & 0x33333333) + ((uljj >> 2) & 0x33333333);
    tmp_stor = (ulii & 0x0f0f0f0f) + ((ulii >> 4) & 0x0f0f0f0f);

    loader = *lptr++;
    ulii = loader - ((loader >> 1) & FIVEMASK);
    loader = *lptr++;
    uljj = loader - ((loader >> 1) & FIVEMASK);
    loader = *lptr++;
    ulii += (loader >> 1) & FIVEMASK;
    uljj += loader & FIVEMASK;
    ulii = (ulii & 0x33333333) + ((ulii >> 2) & 0x33333333);
    ulii += (uljj & 0x33333333) + ((uljj >> 2) & 0x33333333);
    tmp_stor += (ulii & 0x0f0f0f0f) + ((ulii >> 4) & 0x0f0f0f0f);

    // Each 8-bit slot stores a number in 0..48.  Multiplying by 0x01010101 is
    // equivalent to the left-shifts and adds we need to sum those four 8-bit
    // numbers in the high-order slot.
    tot += (tmp_stor * 0x01010101) >> 24;
  }
#endif
  while (lptr < lptr_end) {
    tot += popcount_long(*lptr++);
  }
  return tot;
}

uintptr_t popcount2_longs(const uintptr_t* lptr, uintptr_t word_ct) {
  // treats lptr[] as an array of two-bit instead of one-bit numbers
  uintptr_t tot = 0;
  const uintptr_t* lptr_end = &(lptr[word_ct]);
#ifdef __LP64__
  uintptr_t twelve_ct;
  const __m128i* vptr;
  vptr = (const __m128i*)lptr;
  twelve_ct = word_ct / 12;
  tot += popcount2_vecs(vptr, twelve_ct * 6);
  lptr = &(lptr[twelve_ct * 12]);
#else
  const uintptr_t* lptr_six_end;
  uintptr_t loader1;
  uintptr_t loader2;
  uintptr_t ulii;
  uintptr_t uljj;
  lptr_six_end = &(lptr[word_ct - (word_ct % 6)]);
  while (lptr < lptr_six_end) {
    loader1 = *lptr++;
    loader2 = *lptr++;
    ulii = (loader1 & 0x33333333) + ((loader1 >> 2) & 0x33333333);
    uljj = (loader2 & 0x33333333) + ((loader2 >> 2) & 0x33333333);
    loader1 = *lptr++;
    loader2 = *lptr++;
    ulii += (loader1 & 0x33333333) + ((loader1 >> 2) & 0x33333333);
    uljj += (loader2 & 0x33333333) + ((loader2 >> 2) & 0x33333333);
    loader1 = *lptr++;
    loader2 = *lptr++;
    ulii += (loader1 & 0x33333333) + ((loader1 >> 2) & 0x33333333);
    uljj += (loader2 & 0x33333333) + ((loader2 >> 2) & 0x33333333);
    ulii = (ulii & 0x0f0f0f0f) + ((ulii >> 4) & 0x0f0f0f0f);
    ulii += (uljj & 0x0f0f0f0f) + ((uljj >> 4) & 0x0f0f0f0f);

    // Each 8-bit slot stores a number in 0..48.  Multiplying by 0x01010101 is
    // equivalent to the left-shifts and adds we need to sum those four 8-bit
    // numbers in the high-order slot.
    tot += (ulii * 0x01010101) >> 24;
  }
#endif
  while (lptr < lptr_end) {
    tot += popcount2_long(*lptr++);
  }
  return tot;
}

uintptr_t popcount_bit_idx(const uintptr_t* lptr, uintptr_t start_idx, uintptr_t end_idx) {
  uintptr_t start_idxl = start_idx / BITCT;
  uintptr_t start_idxlr = start_idx & (BITCT - 1);
  uintptr_t end_idxl = end_idx / BITCT;
  uintptr_t end_idxlr = end_idx & (BITCT - 1);
  uintptr_t ct = 0;
  if (start_idxl == end_idxl) {
    return popcount_long(lptr[start_idxl] & ((ONELU << end_idxlr) - (ONELU << start_idxlr)));
  }
  if (start_idxlr) {
    ct = popcount_long(lptr[start_idxl++] >> start_idxlr);
  }
  if (end_idxl > start_idxl) {
    ct += popcount_longs_nzbase(lptr, start_idxl, end_idxl);
  }
  if (end_idxlr) {
    ct += popcount_long(lptr[end_idxl] & ((ONELU << end_idxlr) - ONELU));
  }
  return ct;
}

uint32_t chrom_window_max(const uint32_t* marker_pos, const uintptr_t* marker_exclude, const Chrom_info* chrom_info_ptr, uint32_t chrom_idx, uint32_t ct_max, uint32_t bp_max, uint32_t cur_window_max) {
  // okay, it's absurd to keep rewriting this from scratch, especially given
  // that makes it likely that some reimplementations suck (--indep{-pairwise}
  // version was O(n^2) instead of O(n); sure, it didn't really matter because
  // the main calculation was more expensive, but still, ugh).

  if (cur_window_max >= ct_max) {
    return ct_max;
  }
  // assumes chrom_idx exists
  uint32_t chrom_fo_idx = chrom_info_ptr->chrom_idx_to_foidx[chrom_idx];
  uint32_t chrom_end = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx + 1];
  uint32_t marker_uidx = next_unset(marker_exclude, chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx], chrom_end);
  uint32_t marker_ct = chrom_end - marker_uidx - popcount_bit_idx(marker_exclude, marker_uidx, chrom_end);
  if (marker_ct <= cur_window_max) {
    return cur_window_max;
  }
  uint32_t window_idx_first = 0;
  uint32_t window_uidx_first = marker_uidx;
  uint32_t window_pos_first = marker_pos[marker_uidx];
  uint32_t marker_idx;
  uint32_t marker_pos_thresh;
  for (marker_idx = 0; marker_idx < marker_ct; marker_uidx++, marker_idx++) {
    next_unset_unsafe_ck(marker_exclude, &marker_uidx);
    marker_pos_thresh = marker_pos[marker_uidx];
    if (marker_pos_thresh < bp_max) {
      marker_pos_thresh = 0;
    } else {
      marker_pos_thresh -= bp_max;
    }
    if (marker_pos_thresh > window_pos_first) {
      do {
        window_uidx_first++;
        next_unset_unsafe_ck(marker_exclude, &window_uidx_first);
        window_pos_first = marker_pos[window_uidx_first];
        window_idx_first++;
      } while (marker_pos_thresh > window_pos_first);
    } else if (marker_idx - window_idx_first == cur_window_max) {
      if (++cur_window_max == ct_max) {
	return cur_window_max;
      }
    }
  }
  return cur_window_max;
}

uint32_t window_back(const uint32_t* __restrict marker_pos, const double* __restrict marker_cms, const uintptr_t* marker_exclude, uint32_t marker_uidx_min, uint32_t marker_uidx_start, uint32_t count_max, uint32_t bp_max, double cm_max, uint32_t* __restrict window_trail_ct_ptr) {
  // Finds the earliest location which is within count_max sites, bp_max bps,
  // and (if marker_cms != nullptr) cm_max centimorgans.
  // count_max must be positive.
  if (marker_uidx_min == marker_uidx_start) {
    // special-case this since it happens frequently
    *window_trail_ct_ptr = 0;
    return marker_uidx_min;
  }
  double min_cm = marker_cms? (marker_cms[marker_uidx_start] - cm_max) : 0.0;
  uint32_t min_pos = 0;
  uint32_t marker_uwidx_cur = marker_uidx_start / BITCT;
  uint32_t uii = marker_uidx_start % BITCT;
  uint32_t marker_uidx_last = marker_uidx_start;
  uint32_t remaining_count = count_max;
  const uintptr_t* marker_exclude_cur = &(marker_exclude[marker_uwidx_cur]);
  uintptr_t cur_word;
  uint32_t ujj;
  uint32_t ukk;
  marker_uwidx_cur *= BITCT;
  if (bp_max <= marker_pos[marker_uidx_start]) {
    min_pos = marker_pos[marker_uidx_start] - bp_max;
  }
  if (!uii) {
    goto window_back_zstart;
  }
  cur_word = (~(*marker_exclude_cur)) & ((ONELU << uii) - ONELU);
  while (1) {
    if (marker_uwidx_cur <= marker_uidx_min) {
      cur_word &= ~((ONELU << (marker_uidx_min % BITCT)) - ONELU);
      marker_uwidx_cur = marker_uidx_min;
      uii = popcount_long(cur_word);
      if (uii >= remaining_count) {
	goto window_back_count;
      } else if ((marker_pos[marker_uwidx_cur] < min_pos) || (marker_cms && (marker_cms[marker_uwidx_cur] < min_cm))) {
	goto window_back_find_first_pos;
      } else {
	goto window_back_min;
      }
    }
    uii = popcount_long(cur_word);
    if (uii >= remaining_count) {
    window_back_count:
      uii -= remaining_count; // now a count of number of bits to advance
      while (uii) {
	cur_word &= cur_word - 1;
        uii--;
      }
      // bugfix (7 May 2017): forgot to round marker_uwidx_cur down to word
      //   boundary, before adding CTZLU(cur_word) offset
      marker_uwidx_cur = (marker_uwidx_cur & (~(BITCT - ONELU))) + CTZLU(cur_word);
      if ((marker_pos[marker_uwidx_cur] < min_pos) || (marker_cms && (marker_cms[marker_uwidx_cur] < min_cm))) {
	goto window_back_find_first_pos;
      }
      *window_trail_ct_ptr = count_max;
      return marker_uwidx_cur;
    }
    if ((marker_pos[marker_uwidx_cur] < min_pos) || (marker_cms && (marker_cms[marker_uwidx_cur] < min_cm))) {
    window_back_find_first_pos:
      ujj = uint32arr_greater_than(&(marker_pos[marker_uwidx_cur]), marker_uidx_last - marker_uwidx_cur, min_pos);
      if (marker_cms) {
	ukk = doublearr_greater_than(&(marker_cms[marker_uwidx_cur]), marker_uidx_last - marker_uwidx_cur, min_cm);
	if (ujj < ukk) {
	  ujj = ukk;
	}
      }
      marker_uwidx_cur += ujj;
      if (marker_uwidx_cur > marker_uidx_min) {
	next_unset_unsafe_ck(marker_exclude, &marker_uwidx_cur);
      }
    window_back_min:
      *window_trail_ct_ptr = marker_uidx_start - marker_uwidx_cur - popcount_bit_idx(marker_exclude, marker_uwidx_cur, marker_uidx_start);
      return marker_uwidx_cur;
    }
    remaining_count -= uii;
    marker_uidx_last = marker_uwidx_cur;
  window_back_zstart:
    cur_word = ~(*(--marker_exclude_cur));
    marker_uwidx_cur -= BITCT;
  }
}

uint32_t window_forward(const uint32_t* __restrict marker_pos, const double* __restrict marker_cms, const uintptr_t* marker_exclude, uint32_t marker_uidx_start, uint32_t marker_uidx_last, uint32_t count_max, uint32_t bp_max, double cm_max, uint32_t* __restrict window_lead_ct_ptr) {
  // window_lead_ct_ptr currently cannot be nullptr
  if (marker_uidx_start == marker_uidx_last) {
    *window_lead_ct_ptr = 0;
    return marker_uidx_start;
  }
  double max_cm = marker_cms? (cm_max + marker_cms[marker_uidx_start]) : 0.0;
  uint32_t marker_uwidx_prev = marker_uidx_start;
  uint32_t max_pos = bp_max + marker_pos[marker_uidx_start];
  uint32_t marker_uwidx_cur = (marker_uidx_start + 1) / BITCT;
  uint32_t uii = (marker_uidx_start + 1) % BITCT;
  uint32_t remaining_count = count_max;
  const uintptr_t* marker_exclude_cur = &(marker_exclude[marker_uwidx_cur]);
  uintptr_t cur_word;
  uint32_t ujj;
  uint32_t ukk;
  marker_uwidx_cur *= BITCT;
  cur_word = ~((*marker_exclude_cur) | ((ONELU << uii) - ONELU));
  while (1) {
    uii = popcount_long(cur_word);
    if (uii >= remaining_count) {
      while (--remaining_count) {
	cur_word &= cur_word - 1;
      }
      marker_uwidx_cur += CTZLU(cur_word);
      if (marker_uwidx_cur <= marker_uidx_last) {
	if ((marker_pos[marker_uwidx_cur] > max_pos) || (marker_cms && (marker_cms[marker_uwidx_cur] > max_cm))) {
	  break;
	}
	*window_lead_ct_ptr = count_max;
	return marker_uwidx_cur;
      }
      if ((marker_pos[marker_uidx_last] <= max_pos) && ((!marker_cms) || (marker_cms[marker_uidx_last] <= max_cm))) {
	marker_uwidx_prev = marker_uidx_last;
	goto window_forward_return;
      }
      marker_uwidx_cur = marker_uidx_last;
      break;
    }
    marker_uwidx_cur += BITCT;
    if (marker_uwidx_cur >= marker_uidx_last) {
      if ((marker_pos[marker_uidx_last] <= max_pos) && ((!marker_cms) || (marker_cms[marker_uidx_last] <= max_cm))) {
	marker_uwidx_prev = marker_uidx_last;
	goto window_forward_return;
      }
      marker_uwidx_cur = marker_uidx_last;
      break;
    } else if ((marker_pos[marker_uwidx_cur] > max_pos) || (marker_cms && (marker_cms[marker_uwidx_cur] > max_cm))) {
      break;
    }
    marker_uwidx_prev = marker_uwidx_cur;
    remaining_count -= uii;
    cur_word = ~(*(++marker_exclude_cur));
  }
  ujj = uint32arr_greater_than(&(marker_pos[marker_uwidx_prev]), marker_uwidx_cur - marker_uwidx_prev, max_pos + 1);
  if (marker_cms) {
    ukk = doublearr_greater_than(&(marker_cms[marker_uwidx_prev]), marker_uwidx_cur - marker_uwidx_prev, max_cm * (1 + SMALL_EPSILON));
    if (ujj > ukk) {
      ujj = ukk;
    }
  }
  marker_uwidx_prev += ujj;
  prev_unset_unsafe_ck(marker_exclude, &marker_uwidx_prev);
 window_forward_return:
  *window_lead_ct_ptr = marker_uwidx_prev - marker_uidx_start - popcount_bit_idx(marker_exclude, marker_uidx_start, marker_uwidx_prev);
  return marker_uwidx_prev;
}

uintptr_t jump_forward_unset_unsafe(const uintptr_t* bitvec, uintptr_t cur_pos, uintptr_t forward_ct) {
  // advances forward_ct unset bits; forward_ct must be positive.  (stays put
  // if forward_ct == 1 and current bit is unset.  may want to tweak this
  // interface, easy to introduce off-by-one bugs...)
  // In usual 64-bit case, also assumes bitvec is 16-byte aligned and the end
  // of the trailing 16-byte block can be safely read from.
  uintptr_t widx = cur_pos / BITCT;
  uintptr_t ulii = cur_pos % BITCT;
  const uintptr_t* bptr = &(bitvec[widx]);
  uintptr_t uljj;
  uintptr_t ulkk;
#ifdef __LP64__
  const __m128i* vptr;
#endif
  if (ulii) {
    uljj = (~(*bptr)) >> ulii;
    ulkk = popcount_long(uljj);
    if (ulkk >= forward_ct) {
    jump_forward_unset_unsafe_finish:
      ulkk = CTZLU(uljj);
      while (--forward_ct) {
        uljj &= uljj - 1;
        ulkk = CTZLU(uljj);
      }
      return widx * BITCT + ulii + ulkk;
    }
    forward_ct -= ulkk;
    widx++;
    bptr++;
  }
  ulii = 0;
#ifdef __LP64__
  if (widx & 1) {
    uljj = ~(*bptr);
    ulkk = popcount_long(uljj);
    if (ulkk >= forward_ct) {
      goto jump_forward_unset_unsafe_finish;
    }
    forward_ct -= ulkk;
    bptr++;
  }
  vptr = (const __m128i*)bptr;
  while (forward_ct > BITCT * 6) {
    uljj = ((forward_ct - 1) / (BITCT * 6)) * 3;
    ulkk = popcount_vecs(vptr, uljj);
    vptr = &(vptr[uljj]);
    forward_ct -= uljj * BITCT * 2 - ulkk;
  }
  bptr = (const uintptr_t*)vptr;
  while (forward_ct > BITCT) {
    forward_ct -= popcount_long(~(*bptr++));
  }
#else
  while (forward_ct > BITCT) {
    uljj = (forward_ct - 1) / BITCT;
    ulkk = popcount_longs(bptr, uljj);
    bptr = &(bptr[uljj]);
    forward_ct -= uljj * BITCT - ulkk;
  }
#endif
  while (1) {
    uljj = ~(*bptr);
    ulkk = popcount_long(uljj);
    if (ulkk >= forward_ct) {
      widx = (uintptr_t)(bptr - bitvec);
      goto jump_forward_unset_unsafe_finish;
    }
    forward_ct -= ulkk;
    bptr++;
  }
}

uintptr_t popcount_longs_exclude(const uintptr_t* __restrict lptr, const uintptr_t* __restrict exclude_arr, uintptr_t end_idx) {
  // popcounts lptr ANDNOT exclude_arr[0..(end_idx-1)].
  // N.B. on 64-bit systems, assumes lptr and exclude_arr are 16-byte aligned.
  uintptr_t tot = 0;
  const uintptr_t* lptr_end = &(lptr[end_idx]);
#ifdef __LP64__
  uintptr_t six_ct = end_idx / 6;
  tot += popcount_vecs_exclude((const __m128i*)lptr, (const __m128i*)exclude_arr, six_ct * 3);
  lptr = &(lptr[six_ct * 6]);
  exclude_arr = &(exclude_arr[six_ct * 6]);
#else
  const uintptr_t* lptr_six_end;
  uintptr_t tmp_stor;
  uintptr_t loader;
  uintptr_t ulii;
  uintptr_t uljj;
  lptr_six_end = &(lptr[end_idx - (end_idx % 6)]);
  while (lptr < lptr_six_end) {
    loader = (*lptr++) & (~(*exclude_arr++));
    ulii = loader - ((loader >> 1) & FIVEMASK);
    loader = (*lptr++) & (~(*exclude_arr++));
    uljj = loader - ((loader >> 1) & FIVEMASK);
    loader = (*lptr++) & (~(*exclude_arr++));
    ulii += (loader >> 1) & FIVEMASK;
    uljj += loader & FIVEMASK;
    ulii = (ulii & 0x33333333) + ((ulii >> 2) & 0x33333333);
    ulii += (uljj & 0x33333333) + ((uljj >> 2) & 0x33333333);
    tmp_stor = (ulii & 0x0f0f0f0f) + ((ulii >> 4) & 0x0f0f0f0f);

    loader = (*lptr++) & (~(*exclude_arr++));
    ulii = loader - ((loader >> 1) & FIVEMASK);
    loader = (*lptr++) & (~(*exclude_arr++));
    uljj = loader - ((loader >> 1) & FIVEMASK);
    loader = (*lptr++) & (~(*exclude_arr++));
    ulii += (loader >> 1) & FIVEMASK;
    uljj += loader & FIVEMASK;
    ulii = (ulii & 0x33333333) + ((ulii >> 2) & 0x33333333);
    ulii += (uljj & 0x33333333) + ((uljj >> 2) & 0x33333333);
    tmp_stor += (ulii & 0x0f0f0f0f) + ((ulii >> 4) & 0x0f0f0f0f);

    // Each 8-bit slot stores a number in 0..48.  Multiplying by 0x01010101 is
    // equivalent to the left-shifts and adds we need to sum those four 8-bit
    // numbers in the high-order slot.
    tot += (tmp_stor * 0x01010101) >> 24;
  }
#endif
  while (lptr < lptr_end) {
    tot += popcount_long((*lptr++) & (~(*exclude_arr++)));
  }
  return tot;
}

uintptr_t popcount_longs_intersect(const uintptr_t* __restrict lptr1, const uintptr_t* __restrict lptr2, uintptr_t word_ct) {
  uintptr_t tot = 0;
  const uintptr_t* lptr1_end = &(lptr1[word_ct]);
#ifdef __LP64__
  uintptr_t six_ct = word_ct / 6;
  tot += popcount_vecs_intersect((const __m128i*)lptr1, (const __m128i*)lptr2, six_ct * 3);
  lptr1 = &(lptr1[six_ct * 6]);
  lptr2 = &(lptr2[six_ct * 6]);
#else
  const uintptr_t* lptr1_six_end;
  uintptr_t tmp_stor;
  uintptr_t loader;
  uintptr_t ulii;
  uintptr_t uljj;
  lptr1_six_end = &(lptr1[word_ct - (word_ct % 6)]);
  while (lptr1 < lptr1_six_end) {
    loader = (*lptr1++) & (*lptr2++);
    ulii = loader - ((loader >> 1) & FIVEMASK);
    loader = (*lptr1++) & (*lptr2++);
    uljj = loader - ((loader >> 1) & FIVEMASK);
    loader = (*lptr1++) & (*lptr2++);
    ulii += (loader >> 1) & FIVEMASK;
    uljj += loader & FIVEMASK;
    ulii = (ulii & 0x33333333) + ((ulii >> 2) & 0x33333333);
    ulii += (uljj & 0x33333333) + ((uljj >> 2) & 0x33333333);
    tmp_stor = (ulii & 0x0f0f0f0f) + ((ulii >> 4) & 0x0f0f0f0f);

    loader = (*lptr1++) & (*lptr2++);
    ulii = loader - ((loader >> 1) & FIVEMASK);
    loader = (*lptr1++) & (*lptr2++);
    uljj = loader - ((loader >> 1) & FIVEMASK);
    loader = (*lptr1++) & (*lptr2++);
    ulii += (loader >> 1) & FIVEMASK;
    uljj += loader & FIVEMASK;
    ulii = (ulii & 0x33333333) + ((ulii >> 2) & 0x33333333);
    ulii += (uljj & 0x33333333) + ((uljj >> 2) & 0x33333333);
    tmp_stor += (ulii & 0x0f0f0f0f) + ((ulii >> 4) & 0x0f0f0f0f);

    // Each 8-bit slot stores a number in 0..48.  Multiplying by 0x01010101 is
    // equivalent to the left-shifts and adds we need to sum those four 8-bit
    // numbers in the high-order slot.
    tot += (tmp_stor * 0x01010101) >> 24;
  }
#endif
  while (lptr1 < lptr1_end) {
    tot += popcount_long((*lptr1++) & (*lptr2++));
  }
  return tot;
}

void vertical_bitct_subtract(const uintptr_t* bitarr, uint32_t item_ct, uint32_t* sum_arr) {
  // assumes trailing bits are zeroed out
  uintptr_t cur_word;
  uint32_t idx_offset;
  uint32_t last_set_bit;
  for (idx_offset = 0; idx_offset < item_ct; idx_offset += BITCT) {
    cur_word = *bitarr++;
    while (cur_word) {
      last_set_bit = CTZLU(cur_word);
      sum_arr[idx_offset + last_set_bit] -= 1;
      cur_word &= cur_word - ONELU;
    }
  }
}

#ifdef __LP64__
void count_2freq_dbl_960b(const VECITYPE* geno_vvec, const VECITYPE* geno_vvec_end, const VECITYPE* __restrict mask1vp, const VECITYPE* __restrict mask2vp, uint32_t* __restrict ct1abp, uint32_t* __restrict ct1cp, uint32_t* __restrict ct2abp, uint32_t* __restrict ct2cp) {
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  __m128i loader;
  __m128i loader2;
  __m128i loader3;
  __m128i to_ct1_ab;
  __m128i to_ct_abtmp;
  __m128i to_ct1_c;
  __m128i to_ct2_ab;
  __m128i to_ct2_c;
  __univec acc1_ab;
  __univec acc1_c;
  __univec acc2_ab;
  __univec acc2_c;

  acc1_ab.vi = _mm_setzero_si128();
  acc1_c.vi = _mm_setzero_si128();
  acc2_ab.vi = _mm_setzero_si128();
  acc2_c.vi = _mm_setzero_si128();
  do {
    loader = *geno_vvec++;
    loader2 = *mask1vp++;
    loader3 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
    loader2 = _mm_and_si128(loader2, loader);
    to_ct1_ab = _mm_add_epi64(loader3, loader2);
    to_ct1_c = _mm_andnot_si128(loader3, loader2);
    loader2 = *mask2vp++;
    loader3 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
    loader2 = _mm_and_si128(loader2, loader);
    to_ct2_ab = _mm_add_epi64(loader3, loader2);
    to_ct2_c = _mm_andnot_si128(loader3, loader2);
    to_ct1_ab = _mm_add_epi64(_mm_and_si128(to_ct1_ab, m2), _mm_and_si128(_mm_srli_epi64(to_ct1_ab, 2), m2));
    to_ct2_ab = _mm_add_epi64(_mm_and_si128(to_ct2_ab, m2), _mm_and_si128(_mm_srli_epi64(to_ct2_ab, 2), m2));

    loader = *geno_vvec++;
    loader2 = *mask1vp++;
    loader3 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
    loader2 = _mm_and_si128(loader2, loader);
    to_ct_abtmp = _mm_add_epi64(loader3, loader2);
    to_ct1_c = _mm_add_epi64(to_ct1_c, _mm_andnot_si128(loader3, loader2));
    to_ct1_ab = _mm_add_epi64(to_ct1_ab, _mm_add_epi64(_mm_and_si128(to_ct_abtmp, m2), _mm_and_si128(_mm_srli_epi64(to_ct_abtmp, 2), m2)));
    loader2 = *mask2vp++;
    loader3 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
    loader2 = _mm_and_si128(loader2, loader);
    to_ct_abtmp = _mm_add_epi64(loader3, loader2);
    to_ct2_c = _mm_add_epi64(to_ct2_c, _mm_andnot_si128(loader3, loader2));
    to_ct2_ab = _mm_add_epi64(to_ct2_ab, _mm_add_epi64(_mm_and_si128(to_ct_abtmp, m2), _mm_and_si128(_mm_srli_epi64(to_ct_abtmp, 2), m2)));

    loader = *geno_vvec++;
    loader2 = *mask1vp++;
    loader3 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
    loader2 = _mm_and_si128(loader2, loader);
    to_ct_abtmp = _mm_add_epi64(loader3, loader2);
    to_ct1_c = _mm_add_epi64(to_ct1_c, _mm_andnot_si128(loader3, loader2));
    to_ct1_ab = _mm_add_epi64(to_ct1_ab, _mm_add_epi64(_mm_and_si128(to_ct_abtmp, m2), _mm_and_si128(_mm_srli_epi64(to_ct_abtmp, 2), m2)));
    loader2 = *mask2vp++;
    loader3 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
    loader2 = _mm_and_si128(loader2, loader);
    to_ct_abtmp = _mm_add_epi64(loader3, loader2);
    to_ct2_c = _mm_add_epi64(to_ct2_c, _mm_andnot_si128(loader3, loader2));
    to_ct2_ab = _mm_add_epi64(to_ct2_ab, _mm_add_epi64(_mm_and_si128(to_ct_abtmp, m2), _mm_and_si128(_mm_srli_epi64(to_ct_abtmp, 2), m2)));

    to_ct1_c = _mm_add_epi64(_mm_and_si128(to_ct1_c, m2), _mm_and_si128(_mm_srli_epi64(to_ct1_c, 2), m2));
    to_ct2_c = _mm_add_epi64(_mm_and_si128(to_ct2_c, m2), _mm_and_si128(_mm_srli_epi64(to_ct2_c, 2), m2));

    acc1_ab.vi = _mm_add_epi64(acc1_ab.vi, _mm_add_epi64(_mm_and_si128(to_ct1_ab, m4), _mm_and_si128(_mm_srli_epi64(to_ct1_ab, 4), m4)));
    acc1_c.vi = _mm_add_epi64(acc1_c.vi, _mm_add_epi64(_mm_and_si128(to_ct1_c, m4), _mm_and_si128(_mm_srli_epi64(to_ct1_c, 4), m4)));
    acc2_ab.vi = _mm_add_epi64(acc2_ab.vi, _mm_add_epi64(_mm_and_si128(to_ct2_ab, m4), _mm_and_si128(_mm_srli_epi64(to_ct2_ab, 4), m4)));
    acc2_c.vi = _mm_add_epi64(acc2_c.vi, _mm_add_epi64(_mm_and_si128(to_ct2_c, m4), _mm_and_si128(_mm_srli_epi64(to_ct2_c, 4), m4)));
  } while (geno_vvec < geno_vvec_end);
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  acc1_ab.vi = _mm_add_epi64(_mm_and_si128(acc1_ab.vi, m8), _mm_and_si128(_mm_srli_epi64(acc1_ab.vi, 8), m8));
  acc1_c.vi = _mm_and_si128(_mm_add_epi64(acc1_c.vi, _mm_srli_epi64(acc1_c.vi, 8)), m8);
  acc2_ab.vi = _mm_add_epi64(_mm_and_si128(acc2_ab.vi, m8), _mm_and_si128(_mm_srli_epi64(acc2_ab.vi, 8), m8));
  acc2_c.vi = _mm_and_si128(_mm_add_epi64(acc2_c.vi, _mm_srli_epi64(acc2_c.vi, 8)), m8);
  *ct1abp += ((acc1_ab.u8[0] + acc1_ab.u8[1]) * 0x1000100010001LLU) >> 48;
  *ct1cp += ((acc1_c.u8[0] + acc1_c.u8[1]) * 0x1000100010001LLU) >> 48;
  *ct2abp += ((acc2_ab.u8[0] + acc2_ab.u8[1]) * 0x1000100010001LLU) >> 48;
  *ct2cp += ((acc2_c.u8[0] + acc2_c.u8[1]) * 0x1000100010001LLU) >> 48;
}

void count_3freq_1920b(const VECITYPE* geno_vvec, const VECITYPE* geno_vvec_end, const VECITYPE* __restrict maskvp, uint32_t* __restrict even_ctp, uint32_t* __restrict odd_ctp, uint32_t* __restrict homset_ctp) {
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  __m128i loader;
  __m128i loader2;
  __m128i loader3;
  __m128i even1;
  __m128i odd1;
  __m128i homset1;
  __m128i even2;
  __m128i odd2;
  __m128i homset2;
  __univec acc_even;
  __univec acc_odd;
  __univec acc_homset;

  acc_even.vi = _mm_setzero_si128();
  acc_odd.vi = _mm_setzero_si128();
  acc_homset.vi = _mm_setzero_si128();
  do {
    loader = *geno_vvec++;
    loader2 = *maskvp++;
    odd1 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
    even1 = _mm_and_si128(loader2, loader);
    homset1 = _mm_and_si128(odd1, loader);
    loader = *geno_vvec++;
    loader2 = *maskvp++;
    loader3 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
    even1 = _mm_add_epi64(even1, _mm_and_si128(loader2, loader));
    odd1 = _mm_add_epi64(odd1, loader3);
    homset1 = _mm_add_epi64(homset1, _mm_and_si128(loader3, loader));
    loader = *geno_vvec++;
    loader2 = *maskvp++;
    loader3 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
    even1 = _mm_add_epi64(even1, _mm_and_si128(loader2, loader));
    odd1 = _mm_add_epi64(odd1, loader3);
    homset1 = _mm_add_epi64(homset1, _mm_and_si128(loader3, loader));

    even1 = _mm_add_epi64(_mm_and_si128(even1, m2), _mm_and_si128(_mm_srli_epi64(even1, 2), m2));
    odd1 = _mm_add_epi64(_mm_and_si128(odd1, m2), _mm_and_si128(_mm_srli_epi64(odd1, 2), m2));
    homset1 = _mm_add_epi64(_mm_and_si128(homset1, m2), _mm_and_si128(_mm_srli_epi64(homset1, 2), m2));

    loader = *geno_vvec++;
    loader2 = *maskvp++;
    odd2 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
    even2 = _mm_and_si128(loader2, loader);
    homset2 = _mm_and_si128(odd2, loader);
    loader = *geno_vvec++;
    loader2 = *maskvp++;
    loader3 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
    even2 = _mm_add_epi64(even2, _mm_and_si128(loader2, loader));
    odd2 = _mm_add_epi64(odd2, loader3);
    homset2 = _mm_add_epi64(homset2, _mm_and_si128(loader3, loader));
    loader = *geno_vvec++;
    loader2 = *maskvp++;
    loader3 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
    even2 = _mm_add_epi64(even2, _mm_and_si128(loader2, loader));
    odd2 = _mm_add_epi64(odd2, loader3);
    homset2 = _mm_add_epi64(homset2, _mm_and_si128(loader3, loader));

    even1 = _mm_add_epi64(even1, _mm_add_epi64(_mm_and_si128(even2, m2), _mm_and_si128(_mm_srli_epi64(even2, 2), m2)));
    odd1 = _mm_add_epi64(odd1, _mm_add_epi64(_mm_and_si128(odd2, m2), _mm_and_si128(_mm_srli_epi64(odd2, 2), m2)));
    homset1 = _mm_add_epi64(homset1, _mm_add_epi64(_mm_and_si128(homset2, m2), _mm_and_si128(_mm_srli_epi64(homset2, 2), m2)));

    acc_even.vi = _mm_add_epi64(acc_even.vi, _mm_add_epi64(_mm_and_si128(even1, m4), _mm_and_si128(_mm_srli_epi64(even1, 4), m4)));
    acc_odd.vi = _mm_add_epi64(acc_odd.vi, _mm_add_epi64(_mm_and_si128(odd1, m4), _mm_and_si128(_mm_srli_epi64(odd1, 4), m4)));
    acc_homset.vi = _mm_add_epi64(acc_homset.vi, _mm_add_epi64(_mm_and_si128(homset1, m4), _mm_and_si128(_mm_srli_epi64(homset1, 4), m4)));
  } while (geno_vvec < geno_vvec_end);
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  acc_even.vi = _mm_add_epi64(_mm_and_si128(acc_even.vi, m8), _mm_and_si128(_mm_srli_epi64(acc_even.vi, 8), m8));
  acc_odd.vi = _mm_add_epi64(_mm_and_si128(acc_odd.vi, m8), _mm_and_si128(_mm_srli_epi64(acc_odd.vi, 8), m8));
  acc_homset.vi = _mm_add_epi64(_mm_and_si128(acc_homset.vi, m8), _mm_and_si128(_mm_srli_epi64(acc_homset.vi, 8), m8));
  *even_ctp += ((acc_even.u8[0] + acc_even.u8[1]) * 0x1000100010001LLU) >> 48;
  *odd_ctp += ((acc_odd.u8[0] + acc_odd.u8[1]) * 0x1000100010001LLU) >> 48;
  *homset_ctp += ((acc_homset.u8[0] + acc_homset.u8[1]) * 0x1000100010001LLU) >> 48;
}
#else
void count_2freq_dbl_24b(const uintptr_t* __restrict geno_vec, const uintptr_t* __restrict mask1p, const uintptr_t* __restrict mask2p, uint32_t* __restrict ct1abp, uint32_t* __restrict ct1cp, uint32_t* __restrict ct2abp, uint32_t* __restrict ct2cp) {
  uintptr_t loader = *geno_vec++;
  uintptr_t loader2 = *mask1p++;
  uintptr_t loader3 = (loader >> 1) & loader2;
  uintptr_t to_ct1_ab;
  uintptr_t to_ct1_c;
  uintptr_t to_ct2_ab;
  uintptr_t to_ct2_c;
  uintptr_t to_ct_abtmp;
  uintptr_t partial1_ab;
  uintptr_t partial1_c;
  uintptr_t partial2_ab;
  uintptr_t partial2_c;
  loader2 &= loader;
  to_ct1_ab = loader2 + loader3;
  to_ct1_c = loader2 & (~loader3);
  loader2 = *mask2p++;
  loader3 = (loader >> 1) & loader2;
  loader2 &= loader;
  to_ct2_ab = loader2 + loader3;
  to_ct2_c = loader2 & (~loader3);

  to_ct1_ab = (to_ct1_ab & 0x33333333) + ((to_ct1_ab >> 2) & 0x33333333);
  to_ct2_ab = (to_ct2_ab & 0x33333333) + ((to_ct2_ab >> 2) & 0x33333333);

  loader = *geno_vec++;
  loader2 = *mask1p++;
  loader3 = (loader >> 1) & loader2;
  loader2 &= loader;
  to_ct_abtmp = loader2 + loader3;
  to_ct1_c += loader2 & (~loader3);
  to_ct1_ab += (to_ct_abtmp & 0x33333333) + ((to_ct_abtmp >> 2) & 0x33333333);
  loader2 = *mask2p++;
  loader3 = (loader >> 1) & loader2;
  loader2 &= loader;
  to_ct_abtmp = loader2 + loader3;
  to_ct2_c += loader2 & (~loader3);
  to_ct2_ab += (to_ct_abtmp & 0x33333333) + ((to_ct_abtmp >> 2) & 0x33333333);

  loader = *geno_vec++;
  loader2 = *mask1p++;
  loader3 = (loader >> 1) & loader2;
  loader2 &= loader;
  to_ct_abtmp = loader2 + loader3;
  to_ct1_c += loader2 & (~loader3);
  to_ct1_ab += (to_ct_abtmp & 0x33333333) + ((to_ct_abtmp >> 2) & 0x33333333);
  loader2 = *mask2p++;
  loader3 = (loader >> 1) & loader2;
  loader2 &= loader;
  to_ct_abtmp = loader2 + loader3;
  to_ct2_c += loader2 & (~loader3);
  to_ct2_ab += (to_ct_abtmp & 0x33333333) + ((to_ct_abtmp >> 2) & 0x33333333);

  partial1_ab = (to_ct1_ab & 0x0f0f0f0f) + ((to_ct1_ab >> 4) & 0x0f0f0f0f);
  partial1_c = (to_ct1_c & 0x33333333) + ((to_ct1_c >> 2) & 0x33333333);
  partial2_ab = (to_ct2_ab & 0x0f0f0f0f) + ((to_ct2_ab >> 4) & 0x0f0f0f0f);
  partial2_c = (to_ct2_c & 0x33333333) + ((to_ct2_c >> 2) & 0x33333333);

  loader = *geno_vec++;
  loader2 = *mask1p++;
  loader3 = (loader >> 1) & loader2;
  loader2 &= loader;
  to_ct1_ab = loader2 + loader3;
  to_ct1_c = loader2 & (~loader3);
  loader2 = *mask2p++;
  loader3 = (loader >> 1) & loader2;
  loader2 &= loader;
  to_ct2_ab = loader2 + loader3;
  to_ct2_c = loader2 & (~loader3);

  to_ct1_ab = (to_ct1_ab & 0x33333333) + ((to_ct1_ab >> 2) & 0x33333333);
  to_ct2_ab = (to_ct2_ab & 0x33333333) + ((to_ct2_ab >> 2) & 0x33333333);

  loader = *geno_vec++;
  loader2 = *mask1p++;
  loader3 = (loader >> 1) & loader2;
  loader2 &= loader;
  to_ct_abtmp = loader2 + loader3;
  to_ct1_c += loader2 & (~loader3);
  to_ct1_ab += (to_ct_abtmp & 0x33333333) + ((to_ct_abtmp >> 2) & 0x33333333);
  loader2 = *mask2p++;
  loader3 = (loader >> 1) & loader2;
  loader2 &= loader;
  to_ct_abtmp = loader2 + loader3;
  to_ct2_c += loader2 & (~loader3);
  to_ct2_ab += (to_ct_abtmp & 0x33333333) + ((to_ct_abtmp >> 2) & 0x33333333);

  loader = *geno_vec++;
  loader2 = *mask1p++;
  loader3 = (loader >> 1) & loader2;
  loader2 &= loader;
  to_ct_abtmp = loader2 + loader3;
  to_ct1_c += loader2 & (~loader3);
  to_ct1_ab += (to_ct_abtmp & 0x33333333) + ((to_ct_abtmp >> 2) & 0x33333333);
  loader2 = *mask2p++;
  loader3 = (loader >> 1) & loader2;
  loader2 &= loader;
  to_ct_abtmp = loader2 + loader3;
  to_ct2_c += loader2 & (~loader3);
  to_ct2_ab += (to_ct_abtmp & 0x33333333) + ((to_ct_abtmp >> 2) & 0x33333333);

  partial1_ab += (to_ct1_ab & 0x0f0f0f0f) + ((to_ct1_ab >> 4) & 0x0f0f0f0f);
  partial1_c += (to_ct1_c & 0x33333333) + ((to_ct1_c >> 2) & 0x33333333);
  partial2_ab += (to_ct2_ab & 0x0f0f0f0f) + ((to_ct2_ab >> 4) & 0x0f0f0f0f);
  partial2_c += (to_ct2_c & 0x33333333) + ((to_ct2_c >> 2) & 0x33333333);

  partial1_c = (partial1_c & 0x0f0f0f0f) + ((partial1_c >> 4) & 0x0f0f0f0f);
  partial2_c = (partial2_c & 0x0f0f0f0f) + ((partial2_c >> 4) & 0x0f0f0f0f);

  *ct1abp += (partial1_ab * 0x01010101) >> 24;
  *ct1cp += (partial1_c * 0x01010101) >> 24;
  *ct2abp += (partial2_ab * 0x01010101) >> 24;
  *ct2cp += (partial2_c * 0x01010101) >> 24;
}

void count_3freq_48b(const uintptr_t* __restrict geno_vec, const uintptr_t* __restrict maskp, uint32_t* __restrict ctap, uint32_t* __restrict ctbp, uint32_t* __restrict ctcp) {
  uintptr_t loader = *geno_vec++;
  uintptr_t loader2 = *maskp++;
  uint32_t to_ct_a1 = loader & loader2;
  uint32_t to_ct_b1 = (loader >> 1) & loader2;
  uint32_t to_ct_c1 = loader & to_ct_b1;
  uintptr_t loader3;
  uint32_t to_ct_a2;
  uint32_t to_ct_b2;
  uint32_t to_ct_c2;
  uintptr_t partial_a;
  uintptr_t partial_b;
  uintptr_t partial_c;
  loader = *geno_vec++;
  loader2 = *maskp++;
  loader3 = (loader >> 1) & loader2;
  to_ct_a1 += loader & loader2;
  to_ct_b1 += loader3;
  to_ct_c1 += loader & loader3;
  loader = *geno_vec++;
  loader2 = *maskp++;
  loader3 = (loader >> 1) & loader2;
  to_ct_a1 += loader & loader2;
  to_ct_b1 += loader3;
  to_ct_c1 += loader & loader3;

  loader = *geno_vec++;
  loader2 = *maskp++;
  to_ct_a2 = loader & loader2;
  to_ct_b2 = (loader >> 1) & loader2;
  to_ct_c2 = loader & to_ct_b2;
  loader = *geno_vec++;
  loader2 = *maskp++;
  loader3 = (loader >> 1) & loader2;
  to_ct_a2 += loader & loader2;
  to_ct_b2 += loader3;
  to_ct_c2 += loader & loader3;
  loader = *geno_vec++;
  loader2 = *maskp++;
  loader3 = (loader >> 1) & loader2;
  to_ct_a2 += loader & loader2;
  to_ct_b2 += loader3;
  to_ct_c2 += loader & loader3;

  to_ct_a1 = (to_ct_a1 & 0x33333333) + ((to_ct_a1 >> 2) & 0x33333333);
  to_ct_a1 += (to_ct_a2 & 0x33333333) + ((to_ct_a2 >> 2) & 0x33333333);
  partial_a = (to_ct_a1 & 0x0f0f0f0f) + ((to_ct_a1 >> 4) & 0x0f0f0f0f);
  to_ct_b1 = (to_ct_b1 & 0x33333333) + ((to_ct_b1 >> 2) & 0x33333333);
  to_ct_b1 += (to_ct_b2 & 0x33333333) + ((to_ct_b2 >> 2) & 0x33333333);
  partial_b = (to_ct_b1 & 0x0f0f0f0f) + ((to_ct_b1 >> 4) & 0x0f0f0f0f);
  to_ct_c1 = (to_ct_c1 & 0x33333333) + ((to_ct_c1 >> 2) & 0x33333333);
  to_ct_c1 += (to_ct_c2 & 0x33333333) + ((to_ct_c2 >> 2) & 0x33333333);
  partial_c = (to_ct_c1 & 0x0f0f0f0f) + ((to_ct_c1 >> 4) & 0x0f0f0f0f);

  loader = *geno_vec++;
  loader2 = *maskp++;
  to_ct_a1 = loader & loader2;
  to_ct_b1 = (loader >> 1) & loader2;
  to_ct_c1 = loader & to_ct_b1;
  loader = *geno_vec++;
  loader2 = *maskp++;
  loader3 = (loader >> 1) & loader2;
  to_ct_a1 += loader & loader2;
  to_ct_b1 += loader3;
  to_ct_c1 += loader & loader3;
  loader = *geno_vec++;
  loader2 = *maskp++;
  loader3 = (loader >> 1) & loader2;
  to_ct_a1 += loader & loader2;
  to_ct_b1 += loader3;
  to_ct_c1 += loader & loader3;

  loader = *geno_vec++;
  loader2 = *maskp++;
  to_ct_a2 = loader & loader2;
  to_ct_b2 = (loader >> 1) & loader2;
  to_ct_c2 = loader & to_ct_b2;
  loader = *geno_vec++;
  loader2 = *maskp++;
  loader3 = (loader >> 1) & loader2;
  to_ct_a2 += loader & loader2;
  to_ct_b2 += loader3;
  to_ct_c2 += loader & loader3;
  loader = *geno_vec;
  loader2 = *maskp;
  loader3 = (loader >> 1) & loader2;
  to_ct_a2 += loader & loader2;
  to_ct_b2 += loader3;
  to_ct_c2 += loader & loader3;

  to_ct_a1 = (to_ct_a1 & 0x33333333) + ((to_ct_a1 >> 2) & 0x33333333);
  to_ct_a1 += (to_ct_a2 & 0x33333333) + ((to_ct_a2 >> 2) & 0x33333333);
  partial_a += (to_ct_a1 & 0x0f0f0f0f) + ((to_ct_a1 >> 4) & 0x0f0f0f0f);
  to_ct_b1 = (to_ct_b1 & 0x33333333) + ((to_ct_b1 >> 2) & 0x33333333);
  to_ct_b1 += (to_ct_b2 & 0x33333333) + ((to_ct_b2 >> 2) & 0x33333333);
  partial_b += (to_ct_b1 & 0x0f0f0f0f) + ((to_ct_b1 >> 4) & 0x0f0f0f0f);
  to_ct_c1 = (to_ct_c1 & 0x33333333) + ((to_ct_c1 >> 2) & 0x33333333);
  to_ct_c1 += (to_ct_c2 & 0x33333333) + ((to_ct_c2 >> 2) & 0x33333333);
  partial_c += (to_ct_c1 & 0x0f0f0f0f) + ((to_ct_c1 >> 4) & 0x0f0f0f0f);

  *ctap += (partial_a * 0x01010101) >> 24;
  *ctbp += (partial_b * 0x01010101) >> 24;
  *ctcp += (partial_c * 0x01010101) >> 24;
}
#endif

#ifdef __LP64__
void count_set_freq_60v(const __m128i* vptr, const __m128i* vend, const __m128i* __restrict include_vec, uint32_t* __restrict set_ctp, uint32_t* __restrict missing_ctp) {
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  __m128i loader;
  __m128i loader2;
  __m128i loader3;
  __m128i odds;
  __m128i evens;
  __m128i missings;
  __univec acc;
  __univec accm;
  acc.vi = _mm_setzero_si128();
  accm.vi = _mm_setzero_si128();
  do {
    loader = *vptr++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader3 = *include_vec++;
    odds = _mm_and_si128(loader2, loader3);
    evens = _mm_and_si128(odds, loader);
    missings = _mm_and_si128(loader, _mm_andnot_si128(loader2, loader3));

    loader = *vptr++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader3 = *include_vec++;
    odds = _mm_add_epi64(odds, _mm_and_si128(loader2, loader3));
    loader3 = _mm_and_si128(loader, loader3);
    evens = _mm_add_epi64(evens, _mm_and_si128(loader2, loader3));
    missings = _mm_add_epi64(missings, _mm_andnot_si128(loader2, loader3));

    loader = *vptr++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader3 = *include_vec++;
    odds = _mm_add_epi64(odds, _mm_and_si128(loader2, loader3));
    loader3 = _mm_and_si128(loader, loader3);
    evens = _mm_add_epi64(evens, _mm_and_si128(loader2, loader3));
    missings = _mm_add_epi64(missings, _mm_andnot_si128(loader2, loader3));

    odds = _mm_add_epi64(_mm_and_si128(odds, m2), _mm_and_si128(_mm_srli_epi64(odds, 2), m2));
    missings = _mm_add_epi64(_mm_and_si128(missings, m2), _mm_and_si128(_mm_srli_epi64(missings, 2), m2));
    odds = _mm_add_epi64(odds, _mm_add_epi64(_mm_and_si128(evens, m2), _mm_and_si128(_mm_srli_epi64(evens, 2), m2)));

    // each 4-bit value here <= 6, so safe to add before m4 mask
    accm.vi = _mm_add_epi64(accm.vi, _mm_and_si128(_mm_add_epi64(missings, _mm_srli_epi64(missings, 4)), m4));

    acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(odds, m4), _mm_and_si128(_mm_srli_epi64(odds, 4), m4)));
  } while (vptr < vend);
  // and each 8-bit value here <= 120
  accm.vi = _mm_and_si128(_mm_add_epi64(accm.vi, _mm_srli_epi64(accm.vi, 8)), m8);

  acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
  *set_ctp += ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
  *missing_ctp += ((accm.u8[0] + accm.u8[1]) * 0x1000100010001LLU) >> 48;
}

void count_set_freq_hap_120v(const __m128i* vptr, const __m128i* vend, const __m128i* __restrict include_vec, uint32_t* __restrict set_ctp, uint32_t* __restrict missing_ctp) {
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  __univec acc;
  __univec accm;
  __m128i loader;
  __m128i loader2;
  __m128i loader3;
  __m128i partial;
  __m128i partialm;
  __m128i partial2;
  __m128i partial2m;
  acc.vi = _mm_setzero_si128();
  accm.vi = _mm_setzero_si128();
  do {
    loader = *vptr++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader3 = *include_vec++;
    partial = _mm_and_si128(loader3, _mm_and_si128(loader, loader2));
    partialm = _mm_and_si128(loader3, _mm_xor_si128(loader, loader2));
    loader = *vptr++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader3 = *include_vec++;
    partial = _mm_add_epi64(partial, _mm_and_si128(loader3, _mm_and_si128(loader, loader2)));
    partialm = _mm_add_epi64(partialm, _mm_and_si128(loader3, _mm_and_si128(loader, loader2)));
    loader = *vptr++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader3 = *include_vec++;
    partial = _mm_add_epi64(partial, _mm_and_si128(loader3, _mm_and_si128(loader, loader2)));
    partialm = _mm_add_epi64(partialm, _mm_and_si128(loader3, _mm_and_si128(loader, loader2)));
    partial2 = _mm_add_epi64(_mm_and_si128(partial, m2), _mm_and_si128(_mm_srli_epi64(partial, 2), m2));
    partial2m = _mm_add_epi64(_mm_and_si128(partialm, m2), _mm_and_si128(_mm_srli_epi64(partialm, 2), m2));

    loader = *vptr++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader3 = *include_vec++;
    partial = _mm_and_si128(loader3, _mm_and_si128(loader, loader2));
    partialm = _mm_and_si128(loader3, _mm_xor_si128(loader, loader2));
    loader = *vptr++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader3 = *include_vec++;
    partial = _mm_add_epi64(partial, _mm_and_si128(loader3, _mm_and_si128(loader, loader2)));
    partialm = _mm_add_epi64(partialm, _mm_and_si128(loader3, _mm_and_si128(loader, loader2)));
    loader = *vptr++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader3 = *include_vec++;
    partial = _mm_add_epi64(partial, _mm_and_si128(loader3, _mm_and_si128(loader, loader2)));
    partialm = _mm_add_epi64(partialm, _mm_and_si128(loader3, _mm_and_si128(loader, loader2)));
    partial2 = _mm_add_epi64(partial2, _mm_add_epi64(_mm_and_si128(partial, m2), _mm_and_si128(_mm_srli_epi64(partial, 2), m2)));
    partial2m = _mm_add_epi64(partial2m, _mm_add_epi64(_mm_and_si128(partialm, m2), _mm_and_si128(_mm_srli_epi64(partialm, 2), m2)));
    acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(partial2, m4), _mm_and_si128(_mm_srli_epi64(partial2, 4), m4)));
    accm.vi = _mm_add_epi64(accm.vi, _mm_add_epi64(_mm_and_si128(partial2m, m4), _mm_and_si128(_mm_srli_epi64(partial2m, 4), m4)));
  } while (vptr < vend);
  acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
  accm.vi = _mm_add_epi64(_mm_and_si128(accm.vi, m8), _mm_and_si128(_mm_srli_epi64(accm.vi, 8), m8));
  *set_ctp += ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
  *missing_ctp += ((accm.u8[0] + accm.u8[1]) * 0x1000100010001LLU) >> 48;
}

void count_set_freq_x_60v(const __m128i* vptr, const __m128i* vend, const __m128i* __restrict include_vec, const __m128i* __restrict male_vec, uint32_t* __restrict set_ctp, uint32_t* __restrict missing_ctp) {
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  __m128i loader;
  __m128i loader2;
  __m128i loader3;
  __m128i loader4;
  __m128i set_odds;
  __m128i set_evens;
  __m128i missings_nm;
  __m128i missings_m;
  __m128i males;
  __univec acc;
  __univec accm;
  acc.vi = _mm_setzero_si128();
  accm.vi = _mm_setzero_si128();
  do {
    loader = *vptr++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader3 = *include_vec++;
    loader4 = _mm_andnot_si128(*male_vec, loader3);
    set_evens = _mm_and_si128(loader, loader4); // subtract missings_nm later
    set_odds = _mm_and_si128(loader2, loader4);
    missings_nm = _mm_andnot_si128(loader2, set_evens);
    males = _mm_and_si128(loader3, *male_vec++);
    set_evens = _mm_or_si128(set_evens, _mm_and_si128(_mm_and_si128(loader, loader2), males));
    missings_m = _mm_and_si128(_mm_xor_si128(loader, loader2), males);

    loader = *vptr++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader3 = *include_vec++;
    loader4 = _mm_andnot_si128(*male_vec, loader3);
    set_odds = _mm_add_epi64(set_odds, _mm_and_si128(loader2, loader4));
    loader4 = _mm_and_si128(loader, loader4);
    set_evens = _mm_add_epi64(set_evens, loader4);
    missings_nm = _mm_add_epi64(missings_nm, _mm_andnot_si128(loader2, loader4));
    loader4 = _mm_and_si128(loader3, *male_vec++);
    set_evens = _mm_add_epi64(set_evens, _mm_and_si128(_mm_and_si128(loader, loader2), loader4));
    missings_m = _mm_add_epi64(missings_m, _mm_and_si128(_mm_xor_si128(loader, loader2), loader4));
    males = _mm_add_epi64(males, loader4);

    loader = *vptr++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader3 = *include_vec++;
    loader4 = _mm_andnot_si128(*male_vec, loader3);
    set_odds = _mm_add_epi64(set_odds, _mm_and_si128(loader2, loader4));
    loader4 = _mm_and_si128(loader, loader4);
    set_evens = _mm_add_epi64(set_evens, loader4);
    missings_nm = _mm_add_epi64(missings_nm, _mm_andnot_si128(loader2, loader4));
    loader4 = _mm_and_si128(loader3, *male_vec++);
    set_evens = _mm_add_epi64(set_evens, _mm_and_si128(_mm_and_si128(loader, loader2), loader4));
    missings_m = _mm_add_epi64(missings_m, _mm_and_si128(_mm_xor_si128(loader, loader2), loader4));
    males = _mm_add_epi64(males, loader4);

    set_evens = _mm_sub_epi64(set_evens, missings_nm);
    missings_nm = _mm_slli_epi64(_mm_add_epi64(_mm_and_si128(missings_nm, m2), _mm_and_si128(_mm_srli_epi64(missings_nm, 2), m2)), 1);
    set_odds = _mm_add_epi64(_mm_and_si128(set_odds, m2), _mm_and_si128(_mm_srli_epi64(set_odds, 2), m2));
    missings_nm = _mm_add_epi64(missings_nm, _mm_add_epi64(_mm_and_si128(missings_m, m2), _mm_and_si128(_mm_srli_epi64(missings_m, 2), m2)));
    set_odds = _mm_add_epi64(set_odds, _mm_add_epi64(_mm_and_si128(set_evens, m2), _mm_and_si128(_mm_srli_epi64(set_evens, 2), m2)));
    missings_nm = _mm_add_epi64(missings_nm, _mm_add_epi64(_mm_and_si128(males, m2), _mm_and_si128(_mm_srli_epi64(males, 2), m2)));
    acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(set_odds, m4), _mm_and_si128(_mm_srli_epi64(set_odds, 4), m4)));
    accm.vi = _mm_add_epi64(accm.vi, _mm_add_epi64(_mm_and_si128(missings_nm, m4), _mm_and_si128(_mm_srli_epi64(missings_nm, 4), m4)));
  } while (vptr < vend);
  acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
  accm.vi = _mm_add_epi64(_mm_and_si128(accm.vi, m8), _mm_and_si128(_mm_srli_epi64(accm.vi, 8), m8));
  *set_ctp += ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
  *missing_ctp += ((accm.u8[0] + accm.u8[1]) * 0x1000100010001LLU) >> 48;
}

void count_set_freq_y_120v(const __m128i* vptr, const __m128i* vend, const __m128i* __restrict include_vec, const __m128i* __restrict nonmale_vec, uint32_t* __restrict set_ctp, uint32_t* __restrict missing_ctp) {
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  __m128i loader;
  __m128i loader2;
  __m128i loader3;
  __m128i loader4;
  __m128i sets1;
  __m128i missings1;
  __m128i sets2;
  __m128i missings2;
  __univec acc;
  __univec accm;
  acc.vi = _mm_setzero_si128();
  accm.vi = _mm_setzero_si128();
  do {
    loader = *vptr++;
    loader3 = *include_vec++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader4 = *nonmale_vec++;
    sets1 = _mm_and_si128(_mm_andnot_si128(loader4, loader3), _mm_and_si128(loader, loader2));
    missings1 = _mm_and_si128(loader3, _mm_or_si128(loader4, _mm_xor_si128(loader, loader2)));

    loader = *vptr++;
    loader3 = *include_vec++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader4 = *nonmale_vec++;
    sets1 = _mm_add_epi64(sets1, _mm_and_si128(_mm_andnot_si128(loader4, loader3), _mm_and_si128(loader, loader2)));
    missings1 = _mm_add_epi64(missings1, _mm_and_si128(loader3, _mm_or_si128(loader4, _mm_xor_si128(loader, loader2))));

    loader = *vptr++;
    loader3 = *include_vec++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader4 = *nonmale_vec++;
    sets1 = _mm_add_epi64(sets1, _mm_and_si128(_mm_andnot_si128(loader4, loader3), _mm_and_si128(loader, loader2)));
    missings1 = _mm_add_epi64(missings1, _mm_and_si128(loader3, _mm_or_si128(loader4, _mm_xor_si128(loader, loader2))));
    sets1 = _mm_add_epi64(_mm_and_si128(sets1, m2), _mm_and_si128(_mm_srli_epi64(sets1, 2), m2));
    missings1 = _mm_add_epi64(_mm_and_si128(missings1, m2), _mm_and_si128(_mm_srli_epi64(missings1, 2), m2));

    loader = *vptr++;
    loader3 = *include_vec++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader4 = *nonmale_vec++;
    sets2 = _mm_and_si128(_mm_andnot_si128(loader4, loader3), _mm_and_si128(loader, loader2));
    missings2 = _mm_and_si128(loader3, _mm_or_si128(loader4, _mm_xor_si128(loader, loader2)));

    loader = *vptr++;
    loader3 = *include_vec++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader4 = *nonmale_vec++;
    sets2 = _mm_add_epi64(sets2, _mm_and_si128(_mm_andnot_si128(loader4, loader3), _mm_and_si128(loader, loader2)));
    missings2 = _mm_add_epi64(missings2, _mm_and_si128(loader3, _mm_or_si128(loader4, _mm_xor_si128(loader, loader2))));

    loader = *vptr++;
    loader3 = *include_vec++;
    loader2 = _mm_srli_epi64(loader, 1);
    loader4 = *nonmale_vec++;
    sets2 = _mm_add_epi64(sets2, _mm_and_si128(_mm_andnot_si128(loader4, loader3), _mm_and_si128(loader, loader2)));
    missings2 = _mm_add_epi64(missings2, _mm_and_si128(loader3, _mm_or_si128(loader4, _mm_xor_si128(loader, loader2))));
    sets1 = _mm_add_epi64(sets1, _mm_add_epi64(_mm_and_si128(sets2, m2), _mm_and_si128(_mm_srli_epi64(sets2, 2), m2)));
    missings1 = _mm_add_epi64(missings1, _mm_add_epi64(_mm_and_si128(missings2, m2), _mm_and_si128(_mm_srli_epi64(missings2, 2), m2)));
    acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(sets1, m4), _mm_and_si128(_mm_srli_epi64(sets1, 4), m4)));
    accm.vi = _mm_add_epi64(accm.vi, _mm_add_epi64(_mm_and_si128(missings1, m4), _mm_and_si128(_mm_srli_epi64(missings1, 4), m4)));
  } while (vptr < vend);
  acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
  accm.vi = _mm_add_epi64(_mm_and_si128(accm.vi, m8), _mm_and_si128(_mm_srli_epi64(accm.vi, 8), m8));
  *set_ctp += ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
  *missing_ctp += ((accm.u8[0] + accm.u8[1]) * 0x1000100010001LLU) >> 48;
}

uintptr_t count_01_vecs(const __m128i* vptr, uintptr_t vct) {
  // counts number of aligned 01s (i.e. PLINK missing genotypes) in
  // [vptr, vend).  Assumes number of words in interval is a multiple of 12.
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  uintptr_t tot = 0;
  const __m128i* vend;
  __m128i loader1;
  __m128i loader2;
  __m128i count1;
  __m128i count2;
  __univec acc;

  while (vct >= 60) {
    vct -= 60;
    vend = &(vptr[60]);
  count_01_vecs_main_loop:
    acc.vi = _mm_setzero_si128();
    do {
      loader1 = *vptr++;
      loader2 = *vptr++;
      count1 = _mm_and_si128(_mm_andnot_si128(_mm_srli_epi64(loader1, 1), loader1), m1);
      count2 = _mm_and_si128(_mm_andnot_si128(_mm_srli_epi64(loader2, 1), loader2), m1);
      loader1 = *vptr++;
      loader2 = *vptr++;
      count1 = _mm_add_epi64(count1, _mm_and_si128(_mm_andnot_si128(_mm_srli_epi64(loader1, 1), loader1), m1));
      count2 = _mm_add_epi64(count2, _mm_and_si128(_mm_andnot_si128(_mm_srli_epi64(loader2, 1), loader2), m1));
      loader1 = *vptr++;
      loader2 = *vptr++;
      count1 = _mm_add_epi64(count1, _mm_and_si128(_mm_andnot_si128(_mm_srli_epi64(loader1, 1), loader1), m1));
      count2 = _mm_add_epi64(count2, _mm_and_si128(_mm_andnot_si128(_mm_srli_epi64(loader2, 1), loader2), m1));
      count1 = _mm_add_epi64(_mm_and_si128(count1, m2), _mm_and_si128(_mm_srli_epi64(count1, 2), m2));
      count1 = _mm_add_epi64(count1, _mm_add_epi64(_mm_and_si128(count2, m2), _mm_and_si128(_mm_srli_epi64(count2, 2), m2)));
      acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(count1, m4), _mm_and_si128(_mm_srli_epi64(count1, 4), m4)));
    } while (vptr < vend);
    acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
    tot += ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
  }
  if (vct) {
    vend = &(vptr[vct]);
    vct = 0;
    goto count_01_vecs_main_loop;
  }
  return tot;
}

#else
void count_set_freq_6(const uintptr_t* __restrict lptr, const uintptr_t* __restrict include_vec, uint32_t* __restrict set_ctp, uint32_t* __restrict missing_ctp) {
  uintptr_t loader = *lptr++;
  uintptr_t loader2 = loader >> 1;
  uintptr_t loader3 = *include_vec++;
  uintptr_t odds = loader2 & loader3;
  uintptr_t evens = odds & loader;
  uintptr_t missings = (~loader2) & loader3 & loader;
  uintptr_t acc;
  uintptr_t accm;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  odds += loader2 & loader3;
  loader3 &= loader;
  evens += loader2 & loader3;
  missings += (~loader2) & loader3;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  odds += loader2 & loader3;
  loader3 &= loader;
  evens += loader2 & loader3;
  missings += (~loader2) & loader3;

  odds = (odds & 0x33333333) + ((odds >> 2) & 0x33333333);
  odds += (evens & 0x33333333) + ((evens >> 2) & 0x33333333);
  accm = (missings & 0x33333333) + ((missings >> 2) & 0x33333333);
  acc = (odds & 0x0f0f0f0f) + ((odds >> 4) & 0x0f0f0f0f);

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  odds = loader2 & loader3;
  evens = odds & loader;
  missings = (~loader2) & loader3 & loader;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  odds += loader2 & loader3;
  loader3 &= loader;
  evens += loader2 & loader3;
  missings += (~loader2) & loader3;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  odds += loader2 & loader3;
  loader3 &= loader;
  evens += loader2 & loader3;
  missings += (~loader2) & loader3;

  odds = (odds & 0x33333333) + ((odds >> 2) & 0x33333333);
  accm += (missings & 0x33333333) + ((missings >> 2) & 0x33333333);
  odds += (evens & 0x33333333) + ((evens >> 2) & 0x33333333);
  accm = (accm & 0x0f0f0f0f) + ((accm >> 4) & 0x0f0f0f0f);
  acc += (odds & 0x0f0f0f0f) + ((odds >> 4) & 0x0f0f0f0f);
  *set_ctp += (acc * 0x01010101) >> 24;
  *missing_ctp += (accm * 0x01010101) >> 24;
}

void count_set_freq_hap_12(const uintptr_t* __restrict lptr, const uintptr_t* __restrict include_vec, uint32_t* __restrict set_ctp, uint32_t* __restrict missing_ctp) {
  uintptr_t loader = *lptr++;
  uintptr_t loader2 = loader >> 1;
  uintptr_t loader3 = *include_vec++;
  uintptr_t partial = loader & loader2 & loader3;
  uintptr_t partialm = (loader ^ loader2) & loader3;
  uintptr_t partial2;
  uintptr_t partial2m;
  uintptr_t acc;
  uintptr_t accm;
  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  partial += loader & loader2 & loader3;
  partialm += (loader ^ loader2) & loader3;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  partial += loader & loader2 & loader3;
  partialm += (loader ^ loader2) & loader3;
  partial2 = (partial & 0x33333333) + ((partial >> 2) & 0x33333333);
  partial2m = (partialm & 0x33333333) + ((partialm >> 2) & 0x33333333);

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  partial = loader & loader2 & loader3;
  partialm = (loader ^ loader2) & loader3;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  partial += loader & loader2 & loader3;
  partialm += (loader ^ loader2) & loader3;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  partial += loader & loader2 & loader3;
  partialm += (loader ^ loader2) & loader3;
  partial2 += (partial & 0x33333333) + ((partial >> 2) & 0x33333333);
  partial2m += (partialm & 0x33333333) + ((partialm >> 2) & 0x33333333);
  acc = (partial2 & 0x0f0f0f0f) + ((partial2 >> 4) & 0x0f0f0f0f);
  accm = (partial2m & 0x0f0f0f0f) + ((partial2m >> 4) & 0x0f0f0f0f);

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  partial = loader & loader2 & loader3;
  partialm = (loader ^ loader2) & loader3;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  partial += loader & loader2 & loader3;
  partialm += (loader ^ loader2) & loader3;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  partial += loader & loader2 & loader3;
  partialm += (loader ^ loader2) & loader3;
  partial2 = (partial & 0x33333333) + ((partial >> 2) & 0x33333333);
  partial2m = (partialm & 0x33333333) + ((partialm >> 2) & 0x33333333);

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  partial = loader & loader2 & loader3;
  partialm = (loader ^ loader2) & loader3;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  partial += loader & loader2 & loader3;
  partialm += (loader ^ loader2) & loader3;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  partial += loader & loader2 & loader3;
  partialm += (loader ^ loader2) & loader3;
  partial2 += (partial & 0x33333333) + ((partial >> 2) & 0x33333333);
  partial2m += (partialm & 0x33333333) + ((partialm >> 2) & 0x33333333);
  acc += (partial2 & 0x0f0f0f0f) + ((partial2 >> 4) & 0x0f0f0f0f);
  accm += (partial2m & 0x0f0f0f0f) + ((partial2m >> 4) & 0x0f0f0f0f);
  *set_ctp += (acc * 0x01010101) >> 24;
  *missing_ctp += (accm * 0x01010101) >> 24;
}

void count_set_freq_x_6(const uintptr_t* __restrict lptr, const uintptr_t* __restrict include_vec, const uintptr_t* __restrict male_vec, uint32_t* __restrict set_ctp, uint32_t* __restrict missing_ctp) {
  uintptr_t loader = *lptr++;
  uintptr_t loader2 = loader >> 1;
  uintptr_t loader3 = *include_vec++;
  uintptr_t loader4 = loader3 & (~(*male_vec));
  uintptr_t set_odds = loader2 & loader4;
  uintptr_t set_evens = loader & loader4;
  uintptr_t missings_nm = set_evens & (~loader2);
  uintptr_t missings_m;
  uintptr_t males;
  uintptr_t acc;
  uintptr_t accm;
  males = loader3 & (*male_vec++);
  set_evens |= loader & loader2 & males;
  missings_m = (loader ^ loader2) & males;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = loader3 & (~(*male_vec));
  set_odds += loader2 & loader4;
  loader4 &= loader;
  set_evens += loader4;
  missings_nm += loader4 & (~loader2);
  loader4 = loader3 & (*male_vec++);
  set_evens += loader & loader2 & loader4;
  missings_m += (loader ^ loader2) & loader4;
  males += loader4;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = loader3 & (~(*male_vec));
  set_odds += loader2 & loader4;
  loader4 &= loader;
  set_evens += loader4;
  missings_nm += loader4 & (~loader2);
  loader4 = loader3 & (*male_vec++);
  set_evens += loader & loader2 & loader4;
  missings_m += (loader ^ loader2) & loader4;
  males += loader4;

  set_evens -= missings_nm;
  set_odds = (set_odds & 0x33333333) + ((set_odds >> 2) & 0x33333333);
  set_odds += (set_evens & 0x33333333) + ((set_evens >> 2) & 0x33333333);
  missings_nm = ((missings_nm & 0x33333333) + ((missings_nm >> 2) & 0x33333333)) * 2;
  missings_nm += (missings_m & 0x33333333) + ((missings_m >> 2) & 0x33333333);
  missings_nm += (males & 0x33333333) + ((males >> 2) & 0x33333333);
  acc = (set_odds & 0x0f0f0f0f) + ((set_odds >> 4) & 0x0f0f0f0f);
  accm = (missings_nm & 0x0f0f0f0f) + ((missings_nm >> 4) & 0x0f0f0f0f);

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = loader3 & (~(*male_vec));
  set_odds = loader2 & loader4;
  set_evens = loader & loader4;
  missings_nm = set_evens & (~loader2);
  males = loader3 & (*male_vec++);
  set_evens |= loader & loader2 & males;
  missings_m = (loader ^ loader2) & males;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = loader3 & (~(*male_vec));
  set_odds += loader2 & loader4;
  loader4 &= loader;
  set_evens += loader4;
  missings_nm += loader4 & (~loader2);
  loader4 = loader3 & (*male_vec++);
  set_evens += loader & loader2 & loader4;
  missings_m += (loader ^ loader2) & loader4;
  males += loader4;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = loader3 & (~(*male_vec));
  set_odds += loader2 & loader4;
  loader4 &= loader;
  set_evens += loader4;
  missings_nm += loader4 & (~loader2);
  loader4 = loader3 & (*male_vec++);
  set_evens += loader & loader2 & loader4;
  missings_m += (loader ^ loader2) & loader4;
  males += loader4;

  set_evens -= missings_nm;
  set_odds = (set_odds & 0x33333333) + ((set_odds >> 2) & 0x33333333);
  set_odds += (set_evens & 0x33333333) + ((set_evens >> 2) & 0x33333333);
  missings_nm = ((missings_nm & 0x33333333) + ((missings_nm >> 2) & 0x33333333)) * 2;
  missings_nm += (missings_m & 0x33333333) + ((missings_m >> 2) & 0x33333333);
  missings_nm += (males & 0x33333333) + ((males >> 2) & 0x33333333);
  acc += (set_odds & 0x0f0f0f0f) + ((set_odds >> 4) & 0x0f0f0f0f);
  accm += (missings_nm & 0x0f0f0f0f) + ((missings_nm >> 4) & 0x0f0f0f0f);
  *set_ctp += (acc * 0x01010101) >> 24;
  *missing_ctp += (accm * 0x01010101) >> 24;
}

void count_set_freq_y_12(const uintptr_t* __restrict lptr, const uintptr_t* __restrict include_vec, const uintptr_t* __restrict nonmale_vec, uint32_t* __restrict set_ctp, uint32_t* __restrict missing_ctp) {
  uintptr_t loader = *lptr++;
  uintptr_t loader2 = loader >> 1;
  uintptr_t loader3 = *include_vec++;
  uintptr_t loader4 = *nonmale_vec++;
  uintptr_t sets1 = loader3 & loader & loader2 & (~loader4);
  uintptr_t missings1 = loader3 & (loader4 | (loader ^ loader2));
  uintptr_t sets2;
  uintptr_t missings2;
  uintptr_t acc;
  uintptr_t accm;

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = *nonmale_vec++;
  sets1 += loader3 & loader & loader2 & (~loader4);
  missings1 += loader3 & (loader4 | (loader ^ loader2));

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = *nonmale_vec++;
  sets1 += loader3 & loader & loader2 & (~loader4);
  missings1 += loader3 & (loader4 | (loader ^ loader2));
  sets1 = (sets1 & 0x33333333) + ((sets1 >> 2) & 0x33333333);
  missings1 = (missings1 & 0x33333333) + ((missings1 >> 2) & 0x33333333);

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = *nonmale_vec++;
  sets2 = loader3 & loader & loader2 & (~loader4);
  missings2 = loader3 & (loader4 | (loader ^ loader2));

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = *nonmale_vec++;
  sets2 += loader3 & loader & loader2 & (~loader4);
  missings2 += loader3 & (loader4 | (loader ^ loader2));

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = *nonmale_vec++;
  sets2 += loader3 & loader & loader2 & (~loader4);
  missings2 += loader3 & (loader4 | (loader ^ loader2));
  sets1 += (sets2 & 0x33333333) + ((sets2 >> 2) & 0x33333333);
  missings1 += (missings2 & 0x33333333) + ((missings2 >> 2) & 0x33333333);
  acc = (sets1 & 0x0f0f0f0f) + ((sets1 >> 4) & 0x0f0f0f0f);
  accm = (missings1 & 0x0f0f0f0f) + ((missings1 >> 4) & 0x0f0f0f0f);

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = *nonmale_vec++;
  sets1 = loader3 & loader & loader2 & (~loader4);
  missings1 = loader3 & (loader4 | (loader ^ loader2));

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = *nonmale_vec++;
  sets1 += loader3 & loader & loader2 & (~loader4);
  missings1 += loader3 & (loader4 | (loader ^ loader2));

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = *nonmale_vec++;
  sets1 += loader3 & loader & loader2 & (~loader4);
  missings1 += loader3 & (loader4 | (loader ^ loader2));
  sets1 = (sets1 & 0x33333333) + ((sets1 >> 2) & 0x33333333);
  missings1 = (missings1 & 0x33333333) + ((missings1 >> 2) & 0x33333333);

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = *nonmale_vec++;
  sets2 = loader3 & loader & loader2 & (~loader4);
  missings2 = loader3 & (loader4 | (loader ^ loader2));

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = *nonmale_vec++;
  sets2 += loader3 & loader & loader2 & (~loader4);
  missings2 += loader3 & (loader4 | (loader ^ loader2));

  loader = *lptr++;
  loader2 = loader >> 1;
  loader3 = *include_vec++;
  loader4 = *nonmale_vec++;
  sets2 += loader3 & loader & loader2 & (~loader4);
  missings2 += loader3 & (loader4 | (loader ^ loader2));
  sets1 += (sets2 & 0x33333333) + ((sets2 >> 2) & 0x33333333);
  missings1 += (missings2 & 0x33333333) + ((missings2 >> 2) & 0x33333333);
  acc += (sets1 & 0x0f0f0f0f) + ((sets1 >> 4) & 0x0f0f0f0f);
  accm += (missings1 & 0x0f0f0f0f) + ((missings1 >> 4) & 0x0f0f0f0f);
  *set_ctp += (acc * 0x01010101) >> 24;
  *missing_ctp += (accm * 0x01010101) >> 24;
}

uintptr_t count_01_12(const uintptr_t* lptr) {
  uintptr_t loader1 = *lptr++;
  uintptr_t loader2 = *lptr++;
  uintptr_t count1 = loader1 & (~(loader1 >> 1)) & FIVEMASK;
  uintptr_t count2 = loader2 & (~(loader2 >> 1)) & FIVEMASK;
  uintptr_t partial1;
  uintptr_t partial2;
  loader1 = *lptr++;
  loader2 = *lptr++;
  count1 += loader1 & (~(loader1 >> 1)) & FIVEMASK;
  count2 += loader2 & (~(loader2 >> 1)) & FIVEMASK;
  loader1 = *lptr++;
  loader2 = *lptr++;
  count1 += loader1 & (~(loader1 >> 1)) & FIVEMASK;
  count2 += loader2 & (~(loader2 >> 1)) & FIVEMASK;
  partial1 = (count1 & 0x33333333) + ((count1 >> 2) & 0x33333333);
  partial2 = (count2 & 0x33333333) + ((count2 >> 2) & 0x33333333);

  loader1 = *lptr++;
  loader2 = *lptr++;
  count1 = loader1 & (~(loader1 >> 1)) & FIVEMASK;
  count2 = loader2 & (~(loader2 >> 1)) & FIVEMASK;
  loader1 = *lptr++;
  loader2 = *lptr++;
  count1 += loader1 & (~(loader1 >> 1)) & FIVEMASK;
  count2 += loader2 & (~(loader2 >> 1)) & FIVEMASK;
  loader1 = *lptr++;
  loader2 = *lptr++;
  count1 += loader1 & (~(loader1 >> 1)) & FIVEMASK;
  count2 += loader2 & (~(loader2 >> 1)) & FIVEMASK;
  partial1 += (count1 & 0x33333333) + ((count1 >> 2) & 0x33333333);
  partial2 += (count2 & 0x33333333) + ((count2 >> 2) & 0x33333333);

  partial1 = (partial1 & 0x0f0f0f0f) + ((partial1 >> 4) & 0x0f0f0f0f);
  partial1 += (partial2 & 0x0f0f0f0f) + ((partial2 >> 4) & 0x0f0f0f0f);
  return (partial1 * 0x01010101) >> 24;
}
#endif

void genovec_set_freq(const uintptr_t* __restrict geno_vec, const uintptr_t* __restrict include_quatervec, uintptr_t sample_ctl2, uint32_t* __restrict set_ctp, uint32_t* __restrict missing_ctp) {
  // Assuming include_quatervec describes e.g. cases, and an autosomal marker,
  // this counts the number of case set alleles loaded in geno_vec[], as well
  // as the number of cases with missing genotype info.
  // See single_marker_freqs_and_hwe() for discussion.
  // missing count: popcount2(genotype & (~(genotype >> 1)) & 0x5555...)
  // set allele count: popcount(genotype) - missing count
  const uintptr_t* geno_vec_end = &(geno_vec[sample_ctl2]);
  uintptr_t loader;
  uintptr_t loader2;
  uintptr_t missing_incr;
  uint32_t acc = 0;
  uint32_t accm = 0;
#ifdef __LP64__
  uintptr_t cur_decr = 60;
  const uintptr_t* geno_vec_6x_end;
  sample_ctl2 -= sample_ctl2 % 6;
  while (sample_ctl2 >= 60) {
  genovec_set_freq_loop:
    geno_vec_6x_end = &(geno_vec[cur_decr]);
    count_set_freq_60v((const __m128i*)geno_vec, (const __m128i*)geno_vec_6x_end, (const __m128i*)include_quatervec, &acc, &accm);
    geno_vec = geno_vec_6x_end;
    include_quatervec = &(include_quatervec[cur_decr]);
    sample_ctl2 -= cur_decr;
  }
  if (sample_ctl2) {
    cur_decr = sample_ctl2;
    goto genovec_set_freq_loop;
  }
#else
  const uintptr_t* geno_vec_six_end = &(geno_vec[sample_ctl2 - (sample_ctl2 % 6)]);
  while (geno_vec < geno_vec_six_end) {
    count_set_freq_6(geno_vec, include_quatervec, &acc, &accm);
    geno_vec = &(geno_vec[6]);
    include_quatervec = &(include_quatervec[6]);
  }
#endif
  while (geno_vec < geno_vec_end) {
    loader = *geno_vec++;
    loader2 = *include_quatervec++;
    missing_incr = popcount2_long(loader & (~(loader >> 1)) & loader2);
    accm += missing_incr;
    acc += popcount_long(loader & (loader2 * 3)) - missing_incr;
  }
  *set_ctp = acc;
  *missing_ctp = accm;
}

void genovec_set_freq_x(const uintptr_t* __restrict geno_vec, const uintptr_t* __restrict include_quatervec, const uintptr_t* __restrict male_quatervec, uintptr_t sample_ctl2, uint32_t* __restrict set_ctp, uint32_t* __restrict missing_ctp) {
  // diploid counting for nonmales, haploid counting for males
  // missing_ct := male_obs + male_missing + 2 * female_missing
  const uintptr_t* geno_vec_end = &(geno_vec[sample_ctl2]);
  uintptr_t loader;
  uintptr_t loader2;
  uintptr_t loader3;
  uintptr_t loader4;
  uintptr_t missing_incr;
  uint32_t acc = 0;
  uint32_t accm = 0;
#ifdef __LP64__
  uintptr_t cur_decr = 60;
  const uintptr_t* geno_vec_6x_end;
  sample_ctl2 -= sample_ctl2 % 6;
  while (sample_ctl2 >= 60) {
  genovec_set_freq_x_loop:
    geno_vec_6x_end = &(geno_vec[cur_decr]);
    count_set_freq_x_60v((const __m128i*)geno_vec, (const __m128i*)geno_vec_6x_end, (const __m128i*)include_quatervec, (const __m128i*)male_quatervec, &acc, &accm);
    geno_vec = geno_vec_6x_end;
    include_quatervec = &(include_quatervec[cur_decr]);
    male_quatervec = &(male_quatervec[cur_decr]);
    sample_ctl2 -= cur_decr;
  }
  if (sample_ctl2) {
    cur_decr = sample_ctl2;
    goto genovec_set_freq_x_loop;
  }
#else
  const uintptr_t* geno_vec_six_end = &(geno_vec[sample_ctl2 - (sample_ctl2 % 6)]);
  while (geno_vec < geno_vec_six_end) {
    count_set_freq_x_6(geno_vec, include_quatervec, male_quatervec, &acc, &accm);
    geno_vec = &(geno_vec[6]);
    include_quatervec = &(include_quatervec[6]);
    male_quatervec = &(male_quatervec[6]);
  }
#endif
  while (geno_vec < geno_vec_end) {
    loader = *geno_vec++;
    loader2 = loader >> 1;
    loader3 = *include_quatervec++;
    loader4 = loader3 & (~(*male_quatervec));
    missing_incr = popcount2_long(loader & (~loader2) & loader4);
    accm += 2 * missing_incr;
    acc += popcount_long(loader & (loader4 * 3)) - missing_incr;

    loader4 = loader3 & (*male_quatervec++);
    acc += popcount2_long(loader & loader2 & loader4);
    accm += popcount_long(((loader ^ loader2) & loader4) | (loader4 << 1));
  }
  *set_ctp = acc;
  *missing_ctp = accm;
}

void genovec_set_freq_y(const uintptr_t* __restrict geno_vec, const uintptr_t* __restrict include_quatervec, const uintptr_t* __restrict nonmale_quatervec, uintptr_t sample_ctl2, uint32_t* __restrict set_ctp, uint32_t* __restrict missing_ctp) {
  // all nonmales contribute to missing_ct here
  const uintptr_t* geno_vec_end = &(geno_vec[sample_ctl2]);
  uintptr_t loader;
  uintptr_t loader2;
  uintptr_t loader3;
  uintptr_t loader4;
  uint32_t acc = 0;
  uint32_t accm = 0;
#ifdef __LP64__
  uintptr_t cur_decr = 120;
  const uintptr_t* geno_vec_12x_end;
  sample_ctl2 -= sample_ctl2 % 12;
  while (sample_ctl2 >= 120) {
  genovec_set_freq_y_loop:
    geno_vec_12x_end = &(geno_vec[cur_decr]);
    count_set_freq_y_120v((__m128i*)geno_vec, (__m128i*)geno_vec_12x_end, (__m128i*)include_quatervec, (__m128i*)nonmale_quatervec, &acc, &accm);
    geno_vec = geno_vec_12x_end;
    include_quatervec = &(include_quatervec[cur_decr]);
    nonmale_quatervec = &(nonmale_quatervec[cur_decr]);
    sample_ctl2 -= cur_decr;
  }
  if (sample_ctl2) {
    cur_decr = sample_ctl2;
    goto genovec_set_freq_y_loop;
  }
#else
  const uintptr_t* geno_vec_twelve_end = &(geno_vec[sample_ctl2 - (sample_ctl2 % 12)]);
  while (geno_vec < geno_vec_twelve_end) {
    count_set_freq_y_12(geno_vec, include_quatervec, nonmale_quatervec, &acc, &accm);
    geno_vec = &(geno_vec[12]);
    include_quatervec = &(include_quatervec[12]);
    nonmale_quatervec = &(nonmale_quatervec[12]);
  }
#endif
  while (geno_vec < geno_vec_end) {
    loader = *geno_vec++;
    loader2 = loader >> 1;
    loader3 = *include_quatervec++;
    loader4 = *nonmale_quatervec++;
    acc += popcount2_long(loader & loader2 & loader3 & (~loader4));
    accm += popcount2_long(loader3 & ((loader ^ loader2) | loader4));
  }
  *set_ctp = acc;
  *missing_ctp = accm;
}

void genovec_3freq(const uintptr_t* __restrict geno_vec, const uintptr_t* __restrict include_quatervec, uintptr_t sample_ctl2, uint32_t* __restrict missing_ctp, uint32_t* __restrict het_ctp, uint32_t* __restrict homset_ctp) {
  // generic routine for getting all counts.
  const uintptr_t* geno_vec_end = &(geno_vec[sample_ctl2]);
  uintptr_t loader;
  uintptr_t loader2;
  uintptr_t loader3;
  uint32_t acc_even = 0;
  uint32_t acc_odd = 0;
  uint32_t acc_and = 0;
#ifdef __LP64__
  uintptr_t cur_decr = 120;
  const uintptr_t* geno_vec_12x_end;
  sample_ctl2 -= sample_ctl2 % 12;
  while (sample_ctl2 >= 120) {
  genovec_3freq_loop:
    geno_vec_12x_end = &(geno_vec[cur_decr]);
    count_3freq_1920b((const __m128i*)geno_vec, (const __m128i*)geno_vec_12x_end, (const __m128i*)include_quatervec, &acc_even, &acc_odd, &acc_and);
    geno_vec = geno_vec_12x_end;
    include_quatervec = &(include_quatervec[cur_decr]);
    sample_ctl2 -= cur_decr;
  }
  if (sample_ctl2) {
    cur_decr = sample_ctl2;
    goto genovec_3freq_loop;
  }
#else
  const uintptr_t* geno_vec_twelve_end = &(geno_vec[sample_ctl2 - (sample_ctl2 % 12)]);
  while (geno_vec < geno_vec_twelve_end) {
    count_3freq_48b(geno_vec, include_quatervec, &acc_even, &acc_odd, &acc_and);
    geno_vec = &(geno_vec[12]);
    include_quatervec = &(include_quatervec[12]);
  }
#endif
  while (geno_vec < geno_vec_end) {
    loader = *geno_vec++;
    loader2 = *include_quatervec++;
    loader3 = loader2 & (loader >> 1);
    acc_even += popcount2_long(loader & loader2);
    acc_odd += popcount2_long(loader3);
    acc_and += popcount2_long(loader & loader3);
  }
  *missing_ctp = acc_even - acc_and;
  *het_ctp = acc_odd - acc_and;
  *homset_ctp = acc_and;
}

uintptr_t count_01(const uintptr_t* quatervec, uintptr_t word_ct) {
  // really just for getting a missing count
  // unlike popcount01_longs, this does not assume quatervec[] has no 11s
  const uintptr_t* quatervec_end = &(quatervec[word_ct]);
  uintptr_t loader;
#ifdef __LP64__
  uintptr_t acc;
  word_ct -= word_ct % 12;
  acc = count_01_vecs((__m128i*)quatervec, word_ct / 2);
  quatervec = &(quatervec[word_ct]);
#else
  const uintptr_t* quatervec_twelve_end = &(quatervec[word_ct - (word_ct % 12)]);
  uintptr_t acc = 0;
  while (quatervec < quatervec_twelve_end) {
    acc += count_01_12(quatervec);
    quatervec = &(quatervec[12]);
  }
#endif
  while (quatervec < quatervec_end) {
    loader = *quatervec++;
    acc += popcount2_long(loader & (~(loader >> 1)) & FIVEMASK);
  }
  return acc;
}

void fill_all_bits(uintptr_t ct, uintptr_t* bitarr) {
  // leaves bits beyond the end unset
  // ok for ct == 0
  uintptr_t quotient = ct / BITCT;
  uintptr_t remainder = ct % BITCT;
  fill_ulong_one(quotient, bitarr);
  if (remainder) {
    bitarr[quotient] = (ONELU << remainder) - ONELU;
  }
}

uint32_t numeric_range_list_to_bitarr(const Range_list* range_list_ptr, uint32_t item_ct, uint32_t offset, uint32_t ignore_overflow, uintptr_t* bitarr) {
  // bitarr assumed to be initialized
  const char* names = range_list_ptr->names;
  const unsigned char* starts_range = range_list_ptr->starts_range;
  uint32_t name_ct = range_list_ptr->name_ct;
  uint32_t name_max_len = range_list_ptr->name_max_len;
  uint32_t idx_max = item_ct + offset;
  uint32_t name_idx;
  uint32_t idx1;
  uint32_t idx2;
  for (name_idx = 0; name_idx < name_ct; name_idx++) {
    if (scan_uint_capped(&(names[name_idx * name_max_len]), idx_max, &idx1)) {
      if (ignore_overflow) {
	continue;
      }
      return 1;
    }
    if (starts_range[name_idx]) {
      name_idx++;
      if (scan_uint_capped(&(names[name_idx * name_max_len]), idx_max, &idx2)) {
	if (!ignore_overflow) {
	  return 1;
	}
        idx2 = idx_max - 1;
      }
      fill_bits(idx1 - offset, (idx2 - idx1) + 1, bitarr);
    } else {
      set_bit(idx1 - offset, bitarr);
    }
  }
  return 0;
}

int32_t string_range_list_to_bitarr(char* header_line, uint32_t item_ct, uint32_t fixed_len, const Range_list* range_list_ptr, const char* __restrict sorted_ids, const uint32_t* __restrict id_map, const char* __restrict range_list_flag, const char* __restrict file_descrip, uintptr_t* bitarr, int32_t* __restrict seen_idxs) {
  // bitarr assumed to be initialized
  // if fixed_len is zero, header_line is assumed to be a list of
  // space-delimited unequal-length names
  uintptr_t max_id_len = range_list_ptr->name_max_len;
  uintptr_t name_ct = range_list_ptr->name_ct;
  uint32_t item_idx = 0;
  int32_t retval = 0;
  char* bufptr;
  uint32_t cmdline_pos;
  int32_t ii;
  while (1) {
    bufptr = token_endnn(header_line);
    ii = bsearch_str(header_line, (uintptr_t)(bufptr - header_line), sorted_ids, max_id_len, name_ct);
    if (ii != -1) {
      cmdline_pos = id_map[(uint32_t)ii];
      if (seen_idxs[cmdline_pos] != -1) {
	sprintf(g_logbuf, "Error: Duplicate --%s token in %s.\n", range_list_flag, file_descrip);
        goto string_range_list_to_bitarr_ret_INVALID_FORMAT_2;
      }
      seen_idxs[cmdline_pos] = item_idx;
      if (cmdline_pos && range_list_ptr->starts_range[cmdline_pos - 1]) {
        if (seen_idxs[cmdline_pos - 1] == -1) {
          LOGPREPRINTFWW("Error: Second element of --%s range appears before first element in %s.\n", range_list_flag, file_descrip);
          goto string_range_list_to_bitarr_ret_INVALID_CMDLINE_2;
	}
	fill_bits(seen_idxs[cmdline_pos - 1], (item_idx - seen_idxs[cmdline_pos - 1]) + 1, bitarr);
      } else if (!(range_list_ptr->starts_range[cmdline_pos])) {
	SET_BIT(item_idx, bitarr);
      }
    }
    if (++item_idx == item_ct) {
      break;
    }
    if (fixed_len) {
      header_line = &(header_line[fixed_len]);
    } else {
      header_line = skip_initial_spaces(&(bufptr[1]));
    }
  }
  for (cmdline_pos = 0; cmdline_pos < name_ct; cmdline_pos++) {
    if (seen_idxs[cmdline_pos] == -1) {
      goto string_range_list_to_bitarr_ret_INVALID_CMDLINE_3;
    }
  }
  while (0) {
  string_range_list_to_bitarr_ret_INVALID_CMDLINE_3:
    sprintf(g_logbuf, "Error: Missing --%s token in %s.\n", range_list_flag, file_descrip);
  string_range_list_to_bitarr_ret_INVALID_CMDLINE_2:
    logerrprintb();
    retval = RET_INVALID_CMDLINE;
    break;
  string_range_list_to_bitarr_ret_INVALID_FORMAT_2:
    logerrprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
  return retval;
}

int32_t string_range_list_to_bitarr_alloc(char* header_line, uint32_t item_ct, uint32_t fixed_len, const Range_list* range_list_ptr, const char* __restrict range_list_flag, const char* __restrict file_descrip, uintptr_t** bitarr_ptr) {
  // wrapper for string_range_list_to_bitarr which allocates the bitfield and
  // temporary buffers on the heap
  uintptr_t item_ctl = BITCT_TO_WORDCT(item_ct);
  uintptr_t name_ct = range_list_ptr->name_ct;
  int32_t retval = 0;
  int32_t* seen_idxs;
  char* sorted_ids;
  uint32_t* id_map;
  if (bigstack_calloc_ul(item_ctl, bitarr_ptr) ||
      bigstack_alloc_i(name_ct, &seen_idxs)) {
    return RET_NOMEM;
  }
  // kludge to use sort_item_ids()
  fill_ulong_zero(BITCT_TO_WORDCT(name_ct), (uintptr_t*)seen_idxs);
  if (sort_item_ids(name_ct, (uintptr_t*)seen_idxs, 0, range_list_ptr->names, range_list_ptr->name_max_len, 0, 0, strcmp_deref, &sorted_ids, &id_map)) {
    return RET_NOMEM;
  }
  fill_int_one(name_ct, seen_idxs);
  retval = string_range_list_to_bitarr(header_line, item_ct, fixed_len, range_list_ptr, sorted_ids, id_map, range_list_flag, file_descrip, *bitarr_ptr, seen_idxs);
  bigstack_reset(seen_idxs);
  return retval;
}

int32_t string_range_list_to_bitarr2(const char* __restrict sorted_ids, const uint32_t* id_map, uintptr_t item_ct, uintptr_t max_id_len, const Range_list* __restrict range_list_ptr, const char* __restrict range_list_flag, uintptr_t* bitfield_excl) {
  // sorted_ids/id_map is for e.g. marker IDs instead of command line
  // parameters.  bitfield_excl is assumed to be initialized (since its length
  // is not known by this function).
  char* names = range_list_ptr->names;
  const unsigned char* starts_range = range_list_ptr->starts_range;
  uintptr_t name_max_len = range_list_ptr->name_max_len;
  uint32_t name_ct = range_list_ptr->name_ct;
  int32_t retval = 0;
  uint32_t param_idx;
  char* bufptr;
  uint32_t item_uidx;
  uint32_t item_uidx2;
  int32_t ii;
  for (param_idx = 0; param_idx < name_ct; param_idx++) {
    bufptr = &(names[param_idx * name_max_len]);
    ii = bsearch_str_nl(bufptr, sorted_ids, max_id_len, item_ct);
    if (ii == -1) {
      goto string_range_list_to_bitarr2_ret_INVALID_CMDLINE_3;
    }
    item_uidx = id_map[(uint32_t)ii];
    if (starts_range[param_idx]) {
      param_idx++;
      bufptr = &(names[param_idx * name_max_len]);
      ii = bsearch_str_nl(bufptr, sorted_ids, max_id_len, item_ct);
      if (ii == -1) {
        goto string_range_list_to_bitarr2_ret_INVALID_CMDLINE_3;
      }
      item_uidx2 = id_map[(uint32_t)ii];
      if (item_uidx2 < item_uidx) {
	sprintf(g_logbuf, "Error: Second element of --%s range appears before first.\n", range_list_flag);
	goto string_range_list_to_bitarr2_ret_INVALID_CMDLINE_2;
      }
      clear_bits(item_uidx, item_uidx2 - item_uidx + 1, bitfield_excl);
    } else {
      clear_bit(item_uidx, bitfield_excl);
    }
  }
  while (0) {
  string_range_list_to_bitarr2_ret_INVALID_CMDLINE_3:
    sprintf(g_logbuf, "Error: --%s ID not found.\n", range_list_flag);
  string_range_list_to_bitarr2_ret_INVALID_CMDLINE_2:
    logerrprintb();
    retval = RET_INVALID_CMDLINE;
    break;
  }
  return retval;
}

uint32_t count_non_autosomal_markers(const Chrom_info* chrom_info_ptr, const uintptr_t* marker_exclude, uint32_t count_x, uint32_t count_mt) {
  // for backward compatibility, unplaced markers are considered to be
  // autosomal here
  const int32_t x_code = chrom_info_ptr->xymt_codes[X_OFFSET];
  const int32_t y_code = chrom_info_ptr->xymt_codes[Y_OFFSET];
  const int32_t mt_code = chrom_info_ptr->xymt_codes[MT_OFFSET];
  uint32_t ct = 0;
  if (count_x && (x_code != -2)) {
    ct += count_chrom_markers(chrom_info_ptr, marker_exclude, x_code);
  }
  if (y_code != -2) {
    ct += count_chrom_markers(chrom_info_ptr, marker_exclude, y_code);
  }
  if (count_mt && (mt_code != -2)) {
    ct += count_chrom_markers(chrom_info_ptr, marker_exclude, mt_code);
  }
  return ct;
}

int32_t conditional_allocate_non_autosomal_markers(const Chrom_info* chrom_info_ptr, uintptr_t unfiltered_marker_ct, const uintptr_t* marker_exclude_orig, uint32_t marker_ct, uint32_t count_x, uint32_t count_mt, const char* calc_descrip, uintptr_t** marker_exclude_ptr, uint32_t* newly_excluded_ct_ptr) {
  // if all markers are autosomal (or pseudoautosomal) diploid, nothing
  // happens.  otherwise, this creates a marker_exclude copy with
  // non-{autosomal diploid} markers excluded for the caller to use.
  const uintptr_t unfiltered_marker_ctl = BITCT_TO_WORDCT(unfiltered_marker_ct);
  const int32_t* xymt_codes = chrom_info_ptr->xymt_codes;
  uint32_t xymt_cts[XYMT_OFFSET_CT];
  fill_uint_zero(XYMT_OFFSET_CT, xymt_cts);
  if (is_set(chrom_info_ptr->haploid_mask, 0)) {
    *newly_excluded_ct_ptr = marker_ct;
  } else {
    if (count_x && (xymt_codes[X_OFFSET] != -2)) {
      xymt_cts[X_OFFSET] = count_chrom_markers(chrom_info_ptr, marker_exclude_orig, xymt_codes[X_OFFSET]);
    }
    if (xymt_codes[Y_OFFSET] != -2) {
      xymt_cts[Y_OFFSET] = count_chrom_markers(chrom_info_ptr, marker_exclude_orig, xymt_codes[Y_OFFSET]);
    }
    if (count_mt && (xymt_codes[MT_OFFSET] != -2)) {
      xymt_cts[MT_OFFSET] = count_chrom_markers(chrom_info_ptr, marker_exclude_orig, xymt_codes[MT_OFFSET]);
    }
    *newly_excluded_ct_ptr = xymt_cts[X_OFFSET] + xymt_cts[Y_OFFSET] + xymt_cts[MT_OFFSET];
  }
  if (*newly_excluded_ct_ptr) {
    LOGPRINTF("Excluding %u variant%s on non-autosomes from %s.\n", *newly_excluded_ct_ptr, (*newly_excluded_ct_ptr == 1)? "" : "s", calc_descrip);
  }
  if (*newly_excluded_ct_ptr == marker_ct) {
    logerrprint("Error: No variants remaining.\n");
    return RET_INVALID_CMDLINE;
  }
  if (!(*newly_excluded_ct_ptr)) {
    return 0;
  }
  if (bigstack_alloc_ul(unfiltered_marker_ctl, marker_exclude_ptr)) {
    return RET_NOMEM;
  }
  memcpy(*marker_exclude_ptr, marker_exclude_orig, unfiltered_marker_ctl * sizeof(intptr_t));
  for (uint32_t xymt_idx = 0; xymt_idx < XYMT_OFFSET_CT; ++xymt_idx) {
    if (xymt_cts[xymt_idx]) {
      const uint32_t chrom_fo_idx = chrom_info_ptr->chrom_idx_to_foidx[xymt_codes[xymt_idx]];
      fill_bits(chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx], chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx + 1] - chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx], *marker_exclude_ptr);
    }
  }
  return 0;
}

uint32_t get_max_chrom_size(const Chrom_info* chrom_info_ptr, const uintptr_t* marker_exclude, uint32_t* last_chrom_fo_idx_ptr) {
  const uint32_t chrom_ct = chrom_info_ptr->chrom_ct;
  uint32_t max_chrom_size = 0;
  uint32_t last_chrom_fo_idx = 0;
  for (uint32_t chrom_fo_idx = 0; chrom_fo_idx < chrom_ct; chrom_fo_idx++) {
    const uint32_t cur_chrom_size = count_chrom_markers(chrom_info_ptr, marker_exclude, chrom_info_ptr->chrom_file_order[chrom_fo_idx]);
    if (cur_chrom_size) {
      last_chrom_fo_idx = chrom_fo_idx;
      if (cur_chrom_size > max_chrom_size) {
        max_chrom_size = cur_chrom_size;
      }
    }
  }
  if (last_chrom_fo_idx_ptr) {
    *last_chrom_fo_idx_ptr = last_chrom_fo_idx;
  }
  return max_chrom_size;
}

void count_genders(const uintptr_t* __restrict sex_nm, const uintptr_t* __restrict sex_male, const uintptr_t* __restrict sample_exclude, uintptr_t unfiltered_sample_ct, uint32_t* __restrict male_ct_ptr, uint32_t* __restrict female_ct_ptr, uint32_t* __restrict unk_ct_ptr) {
  // unfiltered_sample_ct can be zero
  uint32_t male_ct = 0;
  uint32_t female_ct = 0;
  uint32_t unk_ct = 0;
  uint32_t unfiltered_sample_ctld = unfiltered_sample_ct / BITCT;
  uint32_t unfiltered_sample_ct_rem = unfiltered_sample_ct & (BITCT - 1);
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t sample_bidx;
  for (sample_bidx = 0; sample_bidx < unfiltered_sample_ctld; sample_bidx++) {
    ulii = ~(*sample_exclude++);
  count_genders_last_loop:
    uljj = *sex_nm++;
    unk_ct += popcount_long(ulii & (~uljj));
    ulii &= uljj;
    uljj = *sex_male++;
    male_ct += popcount_long(ulii & uljj);
    female_ct += popcount_long(ulii & (~uljj));
  }
  if (unfiltered_sample_ct_rem) {
    ulii = (~(*sample_exclude)) & ((ONELU << unfiltered_sample_ct_rem) - ONELU);
    unfiltered_sample_ct_rem = 0;
    goto count_genders_last_loop;
  }
  *male_ct_ptr = male_ct;
  *female_ct_ptr = female_ct;
  *unk_ct_ptr = unk_ct;
}

void reverse_loadbuf(uintptr_t unfiltered_sample_ct, unsigned char* loadbuf) {
  // unfiltered_sample_ct can be zero
  uintptr_t sample_bidx = 0;
  unsigned char* loadbuf_end = &(loadbuf[(unfiltered_sample_ct + 3) / 4]);
  unsigned char ucc;
  unsigned char ucc2;
  uintptr_t unfiltered_sample_ctd;
  uint32_t* loadbuf_alias32;
  uint32_t uii;
  uint32_t ujj;
#ifdef __LP64__
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  __m128i* loadbuf_alias;
  __m128i vii;
  __m128i vjj;
  // todo: use this vector loop even when loadbuf is unaligned, so stuff like
  // recode_load_to() is faster
  if (!(((uintptr_t)loadbuf) & 15)) {
    loadbuf_alias = (__m128i*)loadbuf;
    unfiltered_sample_ctd = unfiltered_sample_ct / 64;
    for (; sample_bidx < unfiltered_sample_ctd; sample_bidx++) {
      vii = *loadbuf_alias;
      // we want to exchange 00 and 11, and leave 01/10 untouched.  So make
      // vjj := 11 iff vii is 00/11, and vjj := 00 otherwise; then xor.
      vjj = _mm_andnot_si128(_mm_xor_si128(vii, _mm_srli_epi64(vii, 1)), m1);
      vjj = _mm_or_si128(vjj, _mm_slli_epi64(vjj, 1));
      *loadbuf_alias++ = _mm_xor_si128(vii, vjj);
    }
    loadbuf = (unsigned char*)loadbuf_alias;
  } else if (!(((uintptr_t)loadbuf) & 3)) {
    loadbuf_alias32 = (uint32_t*)loadbuf;
    unfiltered_sample_ctd = unfiltered_sample_ct / BITCT2;
    for (; sample_bidx < unfiltered_sample_ctd; sample_bidx++) {
      uii = *loadbuf_alias32;
      ujj = 0x55555555 & (~(uii ^ (uii >> 1)));
      ujj *= 3;
      *loadbuf_alias32++ = uii ^ ujj;
    }
    loadbuf = (unsigned char*)loadbuf_alias32;
  }
#else
  if (!(((uintptr_t)loadbuf) & 3)) {
    loadbuf_alias32 = (uint32_t*)loadbuf;
    unfiltered_sample_ctd = unfiltered_sample_ct / BITCT2;
    for (; sample_bidx < unfiltered_sample_ctd; sample_bidx++) {
      uii = *loadbuf_alias32;
      ujj = 0x55555555 & (~(uii ^ (uii >> 1)));
      ujj *= 3;
      *loadbuf_alias32++ = uii ^ ujj;
    }
    loadbuf = (unsigned char*)loadbuf_alias32;
  }
#endif
  for (; loadbuf < loadbuf_end;) {
    ucc = *loadbuf;
    ucc2 = 0x55 & (~(ucc ^ (ucc >> 1)));
    ucc2 *= 3;
    *loadbuf++ = ucc ^ ucc2;
  }
  uii = unfiltered_sample_ct & 3;
  if (uii) {
    loadbuf[-1] &= (0xff >> (8 - 2 * uii));
  }
}

// deprecated, try to just use copy_quaterarr_nonempty_subset()
void copy_quaterarr_nonempty_subset_excl(const uintptr_t* __restrict raw_quaterarr, const uintptr_t* __restrict subset_excl, uint32_t raw_quaterarr_size, uint32_t subset_size, uintptr_t* __restrict output_quaterarr) {
  assert(subset_size);
  assert(raw_quaterarr_size >= subset_size);
  uintptr_t cur_output_word = 0;
  uintptr_t* output_quaterarr_last = &(output_quaterarr[subset_size / BITCT2]);
  const uint32_t word_write_halfshift_end = subset_size % BITCT2;
  uint32_t word_write_halfshift = 0;
  // if < 2/3-filled, use sparse copy algorithm
  if (subset_size * (3 * ONELU) < raw_quaterarr_size * (2 * ONELU)) {
    const uint32_t subset_excl_widx_last = raw_quaterarr_size / BITCT;
    uint32_t subset_excl_widx = 0;
    while (1) {
      uintptr_t cur_include_word = ~subset_excl[subset_excl_widx];

      // this, kiddies, is why exclude masks were a mistake.
      if (subset_excl_widx == subset_excl_widx_last) {
	cur_include_word &= (ONELU << (raw_quaterarr_size % BITCT)) - ONELU;
      }

      if (cur_include_word) {
	uint32_t wordhalf_idx = 0;
#ifdef __LP64__
	uint32_t cur_include_halfword = (uint32_t)cur_include_word;
#else
	uint32_t cur_include_halfword = (uint16_t)cur_include_word;
#endif
	while (1) {
	  if (cur_include_halfword) {
	    uintptr_t raw_quaterarr_word = raw_quaterarr[subset_excl_widx * 2 + wordhalf_idx];
	    do {
	      uint32_t rqa_idx_lowbits = __builtin_ctz(cur_include_halfword);
	      cur_output_word |= ((raw_quaterarr_word >> (rqa_idx_lowbits * 2)) & 3) << (word_write_halfshift * 2);
	      if (++word_write_halfshift == BITCT2) {
		*output_quaterarr++ = cur_output_word;
		word_write_halfshift = 0;
		cur_output_word = 0;
	      }
	      cur_include_halfword &= cur_include_halfword - 1;
	    } while (cur_include_halfword);
	  }
	  if (wordhalf_idx) {
	    break;
	  }
	  wordhalf_idx++;
#ifdef __LP64__
	  cur_include_halfword = cur_include_word >> 32;
#else
	  cur_include_halfword = cur_include_word >> 16;
#endif
	}
	if (output_quaterarr == output_quaterarr_last) {
	  if (word_write_halfshift == word_write_halfshift_end) {
            if (word_write_halfshift_end) {
	      *output_quaterarr_last = cur_output_word;
	    }
	    return;
	  }
	}
      }
      subset_excl_widx++;
    }
  }
  // blocked copy
  const uintptr_t* subset_excl_last = &(subset_excl[raw_quaterarr_size / BITCT]);
  while (1) {
    uintptr_t cur_include_word = ~(*subset_excl);
    if (subset_excl == subset_excl_last) {
      cur_include_word &= (ONELU << (raw_quaterarr_size % BITCT)) - ONELU;
    }
    subset_excl++;
    uint32_t wordhalf_idx = 0;
#ifdef __LP64__
    uintptr_t cur_include_halfword = (uint32_t)cur_include_word;
#else
    uint32_t cur_include_halfword = (uint16_t)cur_include_word;
#endif
    while (1) {
      uintptr_t raw_quaterarr_word = *raw_quaterarr++;
      while (cur_include_halfword) {
	uint32_t rqa_idx_lowbits = CTZLU(cur_include_halfword);
	uintptr_t halfword_invshifted = (~cur_include_halfword) >> rqa_idx_lowbits;
	uintptr_t raw_quaterarr_curblock_unmasked = raw_quaterarr_word >> (rqa_idx_lowbits * 2);
	uint32_t rqa_block_len = CTZLU(halfword_invshifted);
	uint32_t block_len_limit = BITCT2 - word_write_halfshift;
	cur_output_word |= raw_quaterarr_curblock_unmasked << (2 * word_write_halfshift);
	if (rqa_block_len < block_len_limit) {
	  word_write_halfshift += rqa_block_len;
	  cur_output_word &= (ONELU << (word_write_halfshift * 2)) - ONELU;
	} else {
	  // no need to mask, extra bits vanish off the high end
	  *output_quaterarr++ = cur_output_word;
	  word_write_halfshift = rqa_block_len - block_len_limit;
	  cur_output_word = (raw_quaterarr_curblock_unmasked >> (2 * block_len_limit)) & ((ONELU << (2 * word_write_halfshift)) - ONELU);
	}
	cur_include_halfword &= (~(ONELU << (rqa_block_len + rqa_idx_lowbits))) + ONELU;
      }
      if (wordhalf_idx) {
	break;
      }
      wordhalf_idx++;
#ifdef __LP64__
      cur_include_halfword = cur_include_word >> 32;
#else
      cur_include_halfword = cur_include_word >> 16;
#endif
    }
    if (output_quaterarr == output_quaterarr_last) {
      if (word_write_halfshift == word_write_halfshift_end) {
	if (word_write_halfshift_end) {
	  *output_quaterarr_last = cur_output_word;
	}
	return;
      }
    }
  }
}

uint32_t load_and_collapse(uint32_t unfiltered_sample_ct, uint32_t sample_ct, const uintptr_t* __restrict sample_exclude, uintptr_t final_mask, uint32_t do_reverse, FILE* bedfile, uintptr_t* __restrict rawbuf, uintptr_t* __restrict mainbuf) {
  assert(unfiltered_sample_ct);
  uint32_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  if (unfiltered_sample_ct == sample_ct) {
    rawbuf = mainbuf;
  }
  if (load_raw(unfiltered_sample_ct4, bedfile, rawbuf)) {
    return RET_READ_FAIL;
  }
  if (unfiltered_sample_ct != sample_ct) {
    copy_quaterarr_nonempty_subset_excl(rawbuf, sample_exclude, unfiltered_sample_ct, sample_ct, mainbuf);
  } else {
    rawbuf[(unfiltered_sample_ct - 1) / BITCT2] &= final_mask;
  }
  if (do_reverse) {
    reverse_loadbuf(sample_ct, (unsigned char*)mainbuf);
  }
  return 0;
}

void copy_quaterarr_nonempty_subset(const uintptr_t* __restrict raw_quaterarr, const uintptr_t* __restrict subset_mask, uint32_t raw_quaterarr_size, uint32_t subset_size, uintptr_t* __restrict output_quaterarr) {
  // in plink 2.0, we probably want (0-based) bit raw_quaterarr_size of
  // subset_mask to be always allocated and unset.  This removes a few special
  // cases re: iterating past the end of arrays.
  assert(subset_size);
  assert(raw_quaterarr_size >= subset_size);
  uintptr_t cur_output_word = 0;
  uintptr_t* output_quaterarr_last = &(output_quaterarr[subset_size / BITCT2]);
  const uint32_t word_write_halfshift_end = subset_size % BITCT2;
  uint32_t word_write_halfshift = 0;
  // if < 2/3-filled, use sparse copy algorithm
  if (subset_size * (3 * ONELU) < raw_quaterarr_size * (2 * ONELU)) {
    uint32_t subset_mask_widx = 0;
    while (1) {
      const uintptr_t cur_include_word = subset_mask[subset_mask_widx];
      if (cur_include_word) {
	uint32_t wordhalf_idx = 0;
#ifdef __LP64__
	uint32_t cur_include_halfword = (uint32_t)cur_include_word;
#else
	uint32_t cur_include_halfword = (uint16_t)cur_include_word;
#endif
	while (1) {
	  if (cur_include_halfword) {
	    uintptr_t raw_quaterarr_word = raw_quaterarr[subset_mask_widx * 2 + wordhalf_idx];
	    do {
	      uint32_t rqa_idx_lowbits = __builtin_ctz(cur_include_halfword);
	      cur_output_word |= ((raw_quaterarr_word >> (rqa_idx_lowbits * 2)) & 3) << (word_write_halfshift * 2);
	      if (++word_write_halfshift == BITCT2) {
		*output_quaterarr++ = cur_output_word;
		word_write_halfshift = 0;
		cur_output_word = 0;
	      }
	      cur_include_halfword &= cur_include_halfword - 1;
	    } while (cur_include_halfword);
	  }
	  if (wordhalf_idx) {
	    break;
	  }
	  wordhalf_idx++;
#ifdef __LP64__
	  cur_include_halfword = cur_include_word >> 32;
#else
	  cur_include_halfword = cur_include_word >> 16;
#endif
	}
	if (output_quaterarr == output_quaterarr_last) {
	  if (word_write_halfshift == word_write_halfshift_end) {
            if (word_write_halfshift_end) {
	      *output_quaterarr_last = cur_output_word;
	    }
	    return;
	  }
	}
      }
      subset_mask_widx++;
    }
  }
  // blocked copy
  while (1) {
    const uintptr_t cur_include_word = *subset_mask++;
    uint32_t wordhalf_idx = 0;
#ifdef __LP64__
    uintptr_t cur_include_halfword = (uint32_t)cur_include_word;
#else
    uint32_t cur_include_halfword = (uint16_t)cur_include_word;
#endif
    while (1) {
      uintptr_t raw_quaterarr_word = *raw_quaterarr++;
      while (cur_include_halfword) {
	uint32_t rqa_idx_lowbits = CTZLU(cur_include_halfword);
	uintptr_t halfword_invshifted = (~cur_include_halfword) >> rqa_idx_lowbits;
	uintptr_t raw_quaterarr_curblock_unmasked = raw_quaterarr_word >> (rqa_idx_lowbits * 2);
	uint32_t rqa_block_len = CTZLU(halfword_invshifted);
	uint32_t block_len_limit = BITCT2 - word_write_halfshift;
	cur_output_word |= raw_quaterarr_curblock_unmasked << (2 * word_write_halfshift);
	if (rqa_block_len < block_len_limit) {
	  word_write_halfshift += rqa_block_len;
	  cur_output_word &= (ONELU << (word_write_halfshift * 2)) - ONELU;
	} else {
	  // no need to mask, extra bits vanish off the high end
	  *output_quaterarr++ = cur_output_word;
	  word_write_halfshift = rqa_block_len - block_len_limit;
	  if (word_write_halfshift) {
	    cur_output_word = (raw_quaterarr_curblock_unmasked >> (2 * block_len_limit)) & ((ONELU << (2 * word_write_halfshift)) - ONELU);
	  } else {
	    // avoid potential right-shift-64
	    cur_output_word = 0;
	  }
	}
	cur_include_halfword &= (~(ONELU << (rqa_block_len + rqa_idx_lowbits))) + ONELU;
      }
      if (wordhalf_idx) {
	break;
      }
      wordhalf_idx++;
#ifdef __LP64__
      cur_include_halfword = cur_include_word >> 32;
#else
      cur_include_halfword = cur_include_word >> 16;
#endif
    }
    if (output_quaterarr == output_quaterarr_last) {
      if (word_write_halfshift == word_write_halfshift_end) {
	if (word_write_halfshift_end) {
	  *output_quaterarr_last = cur_output_word;
	}
	return;
      }
    }
  }
}

/*
void inplace_quaterarr_proper_subset(const uintptr_t* __restrict subset_mask, uint32_t orig_quaterarr_size, uint32_t subset_size, uintptr_t* __restrict main_quaterarr) {
  assert(orig_quaterarr_size > subset_size);
  // worthwhile to special-case this since we get to entirely skip
  // reading/writing these words
  if (!(~subset_mask[0])) {
    const uintptr_t* subset_mask_initial = subset_mask;
    // guaranteed to terminate since orig_quaterarr_size > subset_size.
    do {
      subset_mask++;
    } while (!(~subset_mask[0]));
    const uint32_t quaterarr_word_skip_ct = 2 * ((uintptr_t)(subset_mask - subset_mask_initial));
    main_quaterarr = &(main_quaterarr[quaterarr_word_skip_ct]);
    const uint32_t item_skip_ct = quaterarr_word_skip_ct * BITCT2;
    orig_quaterarr_size -= item_skip_ct;
    subset_size -= item_skip_ct;
  }
  uintptr_t cur_output_word = 0;
  uintptr_t* main_quaterarr_writer = main_quaterarr;
  uintptr_t* main_quaterarr_write_last = &(main_quaterarr[subset_size / BITCT2]);
  const uint32_t word_write_halfshift_end = subset_size % BITCT2;
  uint32_t word_write_halfshift = 0;
  // if <= 2/3-filled, use sparse copy algorithm
  if (subset_size * (3 * ONELU) <= orig_quaterarr_size * (2 * ONELU)) {
    uint32_t subset_mask_widx = 0;
    while (1) {
      const uintptr_t cur_include_word = subset_mask[subset_mask_widx];
      if (cur_include_word) {
	uint32_t wordhalf_idx = 0;
#ifdef __LP64__
	uint32_t cur_include_halfword = (uint32_t)cur_include_word;
#else
	uint32_t cur_include_halfword = (uint16_t)cur_include_word;
#endif
	while (1) {
	  if (cur_include_halfword) {
	    uintptr_t orig_quaterarr_word = main_quaterarr[subset_mask_widx * 2 + wordhalf_idx];
	    do {
	      uint32_t rqa_idx_lowbits = __builtin_ctz(cur_include_halfword);
	      cur_output_word |= ((orig_quaterarr_word >> (rqa_idx_lowbits * 2)) & 3) << (word_write_halfshift * 2);
	      if (++word_write_halfshift == BITCT2) {
		*main_quaterarr_writer++ = cur_output_word;
		word_write_halfshift = 0;
		cur_output_word = 0;
	      }
	      cur_include_halfword &= cur_include_halfword - 1;
	    } while (cur_include_halfword);
	  }
	  if (wordhalf_idx) {
	    break;
	  }
	  wordhalf_idx++;
#ifdef __LP64__
	  cur_include_halfword = cur_include_word >> 32;
#else
	  cur_include_halfword = cur_include_word >> 16;
#endif
	}
	if (main_quaterarr_writer == main_quaterarr_write_last) {
	  if (word_write_halfshift == word_write_halfshift_end) {
            if (word_write_halfshift_end) {
	      *main_quaterarr_writer = cur_output_word;
	    }
	    return;
	  }
	}
      }
      subset_mask_widx++;
    }
  }
  // blocked copy
  while (1) {
    const uintptr_t cur_include_word = *subset_mask++;
    uint32_t wordhalf_idx = 0;
#ifdef __LP64__
    uintptr_t cur_include_halfword = (uint32_t)cur_include_word;
#else
    uint32_t cur_include_halfword = (uint16_t)cur_include_word;
#endif
    while (1) {
      uintptr_t orig_quaterarr_word = *main_quaterarr++;
      while (cur_include_halfword) {
	uint32_t rqa_idx_lowbits = CTZLU(cur_include_halfword);
	uintptr_t halfword_invshifted = (~cur_include_halfword) >> rqa_idx_lowbits;
	uintptr_t orig_quaterarr_curblock_unmasked = orig_quaterarr_word >> (rqa_idx_lowbits * 2);
	uint32_t rqa_block_len = CTZLU(halfword_invshifted);
	uint32_t block_len_limit = BITCT2 - word_write_halfshift;
	cur_output_word |= orig_quaterarr_curblock_unmasked << (2 * word_write_halfshift);
	if (rqa_block_len < block_len_limit) {
	  word_write_halfshift += rqa_block_len;
	  cur_output_word &= (ONELU << (word_write_halfshift * 2)) - ONELU;
	} else {
	  // no need to mask, extra bits vanish off the high end

	  *main_quaterarr_writer++ = cur_output_word;
	  word_write_halfshift = rqa_block_len - block_len_limit;
	  if (word_write_halfshift) {
	    cur_output_word = (orig_quaterarr_curblock_unmasked >> (2 * block_len_limit)) & ((ONELU << (2 * word_write_halfshift)) - ONELU);
	  } else {
	    cur_output_word = 0;
	  }
	}
	cur_include_halfword &= (~(ONELU << (rqa_block_len + rqa_idx_lowbits))) + ONELU;
      }
      if (wordhalf_idx) {
	break;
      }
      wordhalf_idx++;
#ifdef __LP64__
      cur_include_halfword = cur_include_word >> 32;
#else
      cur_include_halfword = cur_include_word >> 16;
#endif
    }
    if (main_quaterarr_writer == main_quaterarr_write_last) {
      if (word_write_halfshift == word_write_halfshift_end) {
	if (word_write_halfshift_end) {
	  *main_quaterarr_writer = cur_output_word;
	}
	return;
      }
    }
  }
}
*/

uint32_t load_and_collapse_incl(uint32_t unfiltered_sample_ct, uint32_t sample_ct, const uintptr_t* __restrict sample_include, uintptr_t final_mask, uint32_t do_reverse, FILE* bedfile, uintptr_t* __restrict rawbuf, uintptr_t* __restrict mainbuf) {
  assert(unfiltered_sample_ct);
  uint32_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  if (unfiltered_sample_ct == sample_ct) {
    rawbuf = mainbuf;
  }
  if (load_raw(unfiltered_sample_ct4, bedfile, rawbuf)) {
    return RET_READ_FAIL;
  }
  if (unfiltered_sample_ct != sample_ct) {
    copy_quaterarr_nonempty_subset(rawbuf, sample_include, unfiltered_sample_ct, sample_ct, mainbuf);
  } else {
    mainbuf[(unfiltered_sample_ct - 1) / BITCT2] &= final_mask;
  }
  if (do_reverse) {
    reverse_loadbuf(sample_ct, (unsigned char*)mainbuf);
  }
  return 0;
}

/*
uint32_t load_and_collapse_incl_inplace(const uintptr_t* __restrict sample_include, uint32_t unfiltered_sample_ct, uint32_t sample_ct, uintptr_t final_mask, uint32_t do_reverse, FILE* bedfile, uintptr_t* __restrict mainbuf) {
  // mainbuf must be large enough to store unfiltered data
  uint32_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  if (load_raw(unfiltered_sample_ct4, bedfile, mainbuf)) {
    return RET_READ_FAIL;
  }
  if (unfiltered_sample_ct == sample_ct) {
    mainbuf[(unfiltered_sample_ct - 1) / BITCT2] &= final_mask;
  } else {
    inplace_quaterarr_proper_subset(sample_include, unfiltered_sample_ct, sample_ct, mainbuf);
  }
  if (do_reverse) {
    reverse_loadbuf(sample_ct, (unsigned char*)mainbuf);
  }
  return 0;
}
*/

uint32_t load_and_split(uint32_t unfiltered_sample_ct, const uintptr_t* __restrict pheno_nm, const uintptr_t* __restrict pheno_c, FILE* bedfile, uintptr_t* __restrict rawbuf, uintptr_t* __restrict casebuf, uintptr_t* __restrict ctrlbuf) {
  // add do_reverse later if needed
  uintptr_t* rawbuf_end = &(rawbuf[unfiltered_sample_ct / BITCT2]);
  uintptr_t case_word = 0;
  uintptr_t ctrl_word = 0;
  uint32_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
  uint32_t case_shift2 = 0;
  uint32_t ctrl_shift2 = 0;
  uint32_t read_shift_max = BITCT2;
  uint32_t sample_uidx = 0;
  uint32_t read_shift;
  uintptr_t read_word;
  uintptr_t ulii;
  if (load_raw(unfiltered_sample_ct4, bedfile, rawbuf)) {
    return RET_READ_FAIL;
  }
  while (1) {
    while (rawbuf < rawbuf_end) {
      read_word = *rawbuf++;
      for (read_shift = 0; read_shift < read_shift_max; sample_uidx++, read_shift++) {
	if (is_set(pheno_nm, sample_uidx)) {
	  ulii = read_word & 3;
	  if (is_set(pheno_c, sample_uidx)) {
	    case_word |= ulii << case_shift2;
	    case_shift2 += 2;
	    if (case_shift2 == BITCT) {
	      *casebuf++ = case_word;
	      case_word = 0;
	      case_shift2 = 0;
	    }
	  } else {
	    ctrl_word |= ulii << ctrl_shift2;
	    ctrl_shift2 += 2;
	    if (ctrl_shift2 == BITCT) {
	      *ctrlbuf++ = ctrl_word;
	      ctrl_word = 0;
	      ctrl_shift2 = 0;
	    }
	  }
	}
	read_word >>= 2;
      }
    }
    if (sample_uidx == unfiltered_sample_ct) {
      if (case_shift2) {
	*casebuf = case_word;
      }
      if (ctrl_shift2) {
	*ctrlbuf = ctrl_word;
      }
      return 0;
    }
    rawbuf_end++;
    read_shift_max = unfiltered_sample_ct % BITCT2;
  }
}

void init_quaterarr_from_bitarr(const uintptr_t* __restrict bitarr, uintptr_t unfiltered_sample_ct, uintptr_t* __restrict new_quaterarr) {
  // allows unfiltered_sample_ct == 0
  uint32_t unfiltered_sample_ctl = BITCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t ulkk;
  uintptr_t ulmm;
  uint32_t bit_idx;
  while (unfiltered_sample_ctl) {
    ulii = ~(*bitarr++);
    ulkk = FIVEMASK;
    ulmm = FIVEMASK;
    if (ulii) {
      uljj = ulii >> BITCT2;
#ifdef __LP64__
      ulii &= 0xffffffffLLU;
#else
      ulii &= 0xffffLU;
#endif
      if (ulii) {
	do {
	  bit_idx = CTZLU(ulii);
	  ulkk &= ~(ONELU << (bit_idx * 2));
	  ulii &= ulii - 1;
	} while (ulii);
      }
      if (uljj) {
	do {
	  bit_idx = CTZLU(uljj);
	  ulmm &= ~(ONELU << (bit_idx * 2));
	  uljj &= uljj - 1;
	} while (uljj);
      }
    }
    *new_quaterarr++ = ulkk;
    *new_quaterarr++ = ulmm;
    --unfiltered_sample_ctl;
  }
  ulii = unfiltered_sample_ct & (BITCT - 1);
  if (ulii) {
    new_quaterarr--;
    if (ulii < BITCT2) {
      *new_quaterarr-- = 0;
    } else {
      ulii -= BITCT2;
    }
    *new_quaterarr &= (ONELU << (ulii * 2)) - ONELU;
  }
}

void init_quaterarr_from_inverted_bitarr(const uintptr_t* __restrict inverted_bitarr, uintptr_t unfiltered_sample_ct, uintptr_t* __restrict new_quaterarr) {
  // allows unfiltered_sample_ct == 0
  uint32_t unfiltered_sample_ctl = BITCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t ulkk;
  uintptr_t ulmm;
  uint32_t bit_idx;
  while (unfiltered_sample_ctl) {
    ulii = *inverted_bitarr++;
    ulkk = FIVEMASK;
    ulmm = FIVEMASK;
    if (ulii) {
      uljj = ulii >> BITCT2;
#ifdef __LP64__
      ulii &= 0xffffffffLLU;
#else
      ulii &= 0xffffLU;
#endif
      if (ulii) {
	do {
	  bit_idx = CTZLU(ulii);
	  ulkk &= ~(ONELU << (bit_idx * 2));
	  ulii &= ulii - 1;
	} while (ulii);
      }
      if (uljj) {
	do {
	  bit_idx = CTZLU(uljj);
	  ulmm &= ~(ONELU << (bit_idx * 2));
	  uljj &= uljj - 1;
	} while (uljj);
      }
    }
    *new_quaterarr++ = ulkk;
    *new_quaterarr++ = ulmm;
    --unfiltered_sample_ctl;
  }
  ulii = unfiltered_sample_ct & (BITCT - 1);
  if (ulii) {
    new_quaterarr--;
    if (ulii < BITCT2) {
      *new_quaterarr-- = 0;
    } else {
      ulii -= BITCT2;
    }
    *new_quaterarr &= (ONELU << (ulii * 2)) - ONELU;
  }
}

void quatervec_01_init_invert(const uintptr_t* __restrict source_quatervec, uintptr_t entry_ct, uintptr_t* __restrict target_quatervec) {
  // Initializes a quatervec as the inverse of another.
  // Some modifications needed for AVX2.
  uint32_t vec_wsize = QUATERCT_TO_ALIGNED_WORDCT(entry_ct);
  uint32_t rem = entry_ct & (BITCT - 1);
#ifdef __LP64__
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  __m128i* tptr = (__m128i*)target_quatervec;
  __m128i* sptr = (__m128i*)source_quatervec;
  __m128i* tptr_end = (__m128i*)(&(target_quatervec[vec_wsize]));
  uintptr_t* second_to_last;
  while (tptr < tptr_end) {
    *tptr++ = _mm_andnot_si128(*sptr++, m1);
  }
  if (rem) {
    second_to_last = &(((uintptr_t*)tptr_end)[-2]);
    if (rem > BITCT2) {
      second_to_last[1] &= (~ZEROLU) >> ((BITCT - rem) * 2);
    } else {
      *second_to_last &= (~ZEROLU) >> ((BITCT2 - rem) * 2);
      second_to_last[1] = 0;
    }
  }
#else
  uintptr_t* tptr_end = &(target_quatervec[vec_wsize]);
  while (target_quatervec < tptr_end) {
    *target_quatervec++ = FIVEMASK & (~(*source_quatervec++));
  }
  if (rem) {
    if (rem > BITCT2) {
      target_quatervec[-1] &= (~ZEROLU) >> ((BITCT - rem) * 2);
    } else {
      target_quatervec[-2] &= (~ZEROLU) >> ((BITCT2 - rem) * 2);
      target_quatervec[-1] = 0;
    }
  }

#endif
}

void bitvec_andnot_copy(const uintptr_t* __restrict source_vec, const uintptr_t* __restrict exclude_vec, uintptr_t word_ct, uintptr_t* __restrict target_vec) {
  // target_vec := source_vec ANDNOT exclude_vec
  // may write an extra word
  assert(word_ct);
#ifdef __LP64__
  __m128i* tptr = (__m128i*)target_vec;
  __m128i* sptr = (__m128i*)source_vec;
  __m128i* xptr = (__m128i*)exclude_vec;
  __m128i* tptr_end = (__m128i*)(&(target_vec[round_up_pow2(word_ct, VEC_WORDS)]));
  do {
    *tptr++ = _mm_andnot_si128(*xptr++, *sptr++);
  } while (tptr < tptr_end);
#else
  uintptr_t* tptr_end = &(target_vec[word_ct]);
  do {
    *target_vec++ = (*source_vec++) & (~(*exclude_vec++));
  } while (target_vec < tptr_end);
#endif
}

void apply_bitarr_mask_to_quaterarr_01(const uintptr_t* __restrict mask_bitarr, uintptr_t unfiltered_sample_ct, uintptr_t* main_quaterarr) {
  // allows unfiltered_sample_ct == 0
  uint32_t unfiltered_sample_ctl = BITCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t ulkk;
  uintptr_t ulmm;
  uint32_t bit_idx;
  while (unfiltered_sample_ctl) {
    ulii = ~(*mask_bitarr++);
    ulkk = *main_quaterarr;
    ulmm = main_quaterarr[1];
    if (ulii) {
      uljj = ulii >> BITCT2;
#ifdef __LP64__
      ulii &= 0xffffffffLLU;
#else
      ulii &= 0xffffLU;
#endif
      if (ulii) {
	do {
	  bit_idx = CTZLU(ulii);
	  ulkk &= ~(ONELU << (bit_idx * 2));
	  ulii &= ulii - 1;
	} while (ulii);
      }
      if (uljj) {
	do {
	  bit_idx = CTZLU(uljj);
	  ulmm &= ~(ONELU << (bit_idx * 2));
	  uljj &= uljj - 1;
	} while (uljj);
      }
    }
    *main_quaterarr++ = ulkk;
    *main_quaterarr++ = ulmm;
    --unfiltered_sample_ctl;
  }
}

void apply_bitarr_excl_to_quaterarr_01(const uintptr_t* __restrict excl_bitarr, uintptr_t unfiltered_sample_ct, uintptr_t* __restrict main_quaterarr) {
  assert(unfiltered_sample_ct);
  uint32_t unfiltered_sample_ctl = BITCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t ulkk;
  uintptr_t ulmm;
  uint32_t bit_idx;
  do {
    ulii = *excl_bitarr++;
    ulkk = *main_quaterarr;
    ulmm = main_quaterarr[1];
    if (ulii) {
      uljj = ulii >> BITCT2;
#ifdef __LP64__
      ulii &= 0xffffffffLLU;
#else
      ulii &= 0xffffLU;
#endif
      if (ulii) {
	do {
	  bit_idx = CTZLU(ulii);
	  ulkk &= ~(ONELU << (bit_idx * 2));
	  ulii &= ulii - 1;
	} while (ulii);
      }
      if (uljj) {
	do {
	  bit_idx = CTZLU(uljj);
	  ulmm &= ~(ONELU << (bit_idx * 2));
	  uljj &= uljj - 1;
	} while (uljj);
      }
    }
    *main_quaterarr++ = ulkk;
    *main_quaterarr++ = ulmm;
  } while (--unfiltered_sample_ctl);
}

void apply_excl_intersect_to_quaterarr_01(const uintptr_t* __restrict excl_bitarr_1, const uintptr_t* __restrict excl_bitarr_2, uintptr_t unfiltered_sample_ct, uintptr_t* __restrict main_quaterarr) {
  assert(unfiltered_sample_ct);
  uint32_t unfiltered_sample_ctl = BITCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t ulii;
  uintptr_t uljj;
  uintptr_t ulkk;
  uintptr_t ulmm;
  uint32_t bit_idx;
  do {
    ulii = (*excl_bitarr_1++) & (*excl_bitarr_2++);
    ulkk = *main_quaterarr;
    ulmm = main_quaterarr[1];
    if (ulii) {
      uljj = ulii >> BITCT2;
#ifdef __LP64__
      ulii &= 0xffffffffLLU;
#else
      ulii &= 0xffffLU;
#endif
      if (ulii) {
	do {
	  bit_idx = CTZLU(ulii);
	  ulkk &= ~(ONELU << (bit_idx * 2));
	  ulii &= ulii - 1;
	} while (ulii);
      }
      if (uljj) {
	do {
	  bit_idx = CTZLU(uljj);
	  ulmm &= ~(ONELU << (bit_idx * 2));
	  uljj &= uljj - 1;
	} while (uljj);
      }
    }
    *main_quaterarr++ = ulkk;
    *main_quaterarr++ = ulmm;
  } while (--unfiltered_sample_ctl);
}

void quatervec_copy_only_01(const uintptr_t* __restrict input_quatervec, uintptr_t unfiltered_sample_ct, uintptr_t* __restrict output_quatervec) {
  // initializes result_ptr bits 01 iff input_quatervec bits are 01
  assert(unfiltered_sample_ct);
#ifdef __LP64__
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  __m128i* vec2_read = (__m128i*)input_quatervec;
  __m128i* read_end = &(vec2_read[QUATERCT_TO_VECCT(unfiltered_sample_ct)]);
  __m128i* vec2_write = (__m128i*)output_quatervec;
  __m128i loader;
  do {
    loader = *vec2_read++;
    *vec2_write++ = _mm_and_si128(_mm_andnot_si128(_mm_srli_epi64(loader, 1), loader), m1);
  } while (vec2_read < read_end);
#else
  const uintptr_t* read_end = &(input_quatervec[QUATERCT_TO_ALIGNED_WORDCT(unfiltered_sample_ct)]);
  uintptr_t loader;
  do {
    loader = *input_quatervec++;
    *output_quatervec++ = loader & (~(loader >> 1)) & FIVEMASK;
  } while (input_quatervec < read_end);
#endif
}

void quatervec_01_invert(uintptr_t unfiltered_sample_ct, uintptr_t* main_quatervec) {
  uintptr_t* vec2_last = &(main_quatervec[unfiltered_sample_ct / BITCT2]);
  uint32_t remainder = unfiltered_sample_ct & (BITCT2 - 1);
#ifdef __LP64__
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  __m128i* vec2_128 = (__m128i*)main_quatervec;
  __m128i* vec2_last128 = &(vec2_128[unfiltered_sample_ct / BITCT]);
  while (vec2_128 < vec2_last128) {
    *vec2_128 = _mm_xor_si128(*vec2_128, m1);
    vec2_128++;
  }
  main_quatervec = (uintptr_t*)vec2_128;
  if (main_quatervec != vec2_last) {
    *main_quatervec = (*main_quatervec) ^ FIVEMASK;
    main_quatervec++;
  }
#else
  while (main_quatervec != vec2_last) {
    *main_quatervec = (*main_quatervec) ^ FIVEMASK;
    main_quatervec++;
  }
#endif
  if (remainder) {
    *vec2_last = *vec2_last ^ (FIVEMASK >> (2 * (BITCT2 - remainder)));
  }
}

void vec_datamask(uintptr_t unfiltered_sample_ct, uint32_t matchval, uintptr_t* data_ptr, uintptr_t* mask_ptr, uintptr_t* result_ptr) {
  // vec_ptr assumed to be standard 00/01 bit vector
  // sets result_vec bits to 01 iff data_ptr bits are equal to matchval and
  // vec_ptr bit is set, 00 otherwise.
  // currently assumes matchval is not 1.
  assert(unfiltered_sample_ct);
#ifdef __LP64__
  __m128i* data_read = (__m128i*)data_ptr;
  __m128i* mask_read = (__m128i*)mask_ptr;
  __m128i* data_read_end = &(data_read[QUATERCT_TO_VECCT(unfiltered_sample_ct)]);
  __m128i* writer = (__m128i*)result_ptr;
  __m128i loader;
#else
  uintptr_t* data_read_end = &(data_ptr[QUATERCT_TO_ALIGNED_WORDCT(unfiltered_sample_ct)]);
  uintptr_t loader;
#endif
  if (matchval) {
    if (matchval == 2) {
#ifdef __LP64__
      do {
        loader = *data_read++;
        *writer++ = _mm_and_si128(_mm_andnot_si128(loader, _mm_srli_epi64(loader, 1)), *mask_read++);
      } while (data_read < data_read_end);
#else
      do {
	loader = *data_ptr++;
        *result_ptr++ = (~loader) & (loader >> 1) & (*mask_ptr++);
      } while (data_ptr < data_read_end);
#endif
    } else {
#ifdef __LP64__
      do {
        loader = *data_read++;
        *writer++ = _mm_and_si128(_mm_and_si128(loader, _mm_srli_epi64(loader, 1)), *mask_read++);
      } while (data_read < data_read_end);
#else
      do {
        loader = *data_ptr++;
        *result_ptr++ = loader & (loader >> 1) & (*mask_ptr++);
      } while (data_ptr < data_read_end);
#endif
    }
  } else {
#ifdef __LP64__
    do {
      loader = *data_read++;
      *writer++ = _mm_andnot_si128(_mm_or_si128(loader, _mm_srli_epi64(loader, 1)), *mask_read++);
    } while (data_read < data_read_end);
#else
    do {
      loader = *data_ptr++;
      *result_ptr++ = (~(loader | (loader >> 1))) & (*mask_ptr++);
    } while (data_ptr < data_read_end);
#endif
  }
}

/*
void vec_rotate_plink1_to_plink2(uintptr_t* lptr, uint32_t word_ct) {
#ifdef __LP64__
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  __m128i* vptr = (__m128i*)lptr;
  __m128i* vend = (__m128i*)(&(lptr[word_ct]));
  __m128i vii;
  __m128i vjj;
  do {
    // new high bit set iff old low bit was set
    // new low bit set iff old bits differed
    vii = *vptr;
    vjj = _mm_and_si128(vii, m1); // old low bit
    vii = _mm_and_si128(_mm_srli_epi64(vii, 1), m1); // old high bit, shifted
    *vptr = _mm_or_si128(_mm_slli_epi64(vjj, 1), _mm_xor_si128(vii, vjj));
  } while (++vptr != vend);
#else
  uintptr_t* lend = &(lptr[word_ct]);
  uintptr_t ulii;
  uintptr_t uljj;
  do {
    ulii = *lptr;
    uljj = ulii & FIVEMASK;
    ulii = (ulii >> 1) & FIVEMASK;
    *lptr = ulii ^ (uljj * 3);
  } while (++lptr != lend);
#endif
}
*/

// this was "rotate_plink1_to_plink2_...", until I noticed that the plink2
// format should store alt allele counts instead of ref allele counts.
void rotate_plink1_to_a2ct_and_copy(uintptr_t* loadbuf, uintptr_t* writebuf, uintptr_t word_ct) {
  // assumes positive word_ct
  uintptr_t* loadbuf_end = &(loadbuf[word_ct]);
  uintptr_t ulii;
  uintptr_t uljj;
  do {
    ulii = *loadbuf++;
    uljj = ulii & FIVEMASK;
    ulii = (ulii >> 1) & FIVEMASK;
    *writebuf++ = ulii ^ (uljj * 3);
  } while (loadbuf < loadbuf_end);
}

void extract_collapsed_missing_bitfield(uintptr_t* lptr, uintptr_t unfiltered_sample_ct, uintptr_t* sample_include_quaterarr, uintptr_t sample_ct, uintptr_t* missing_bitfield) {
  uint32_t word_ct = QUATERCT_TO_WORDCT(unfiltered_sample_ct);
  uintptr_t sample_idx;
  uintptr_t cur_word;
  uintptr_t cur_mask;
  uintptr_t cur_write;
  uint32_t woffset;
  uint32_t widx;
  uint32_t uii;
  if (unfiltered_sample_ct == sample_ct) {
    cur_write = 0;
    woffset = 0;
    for (widx = 0; widx < word_ct; widx++) {
      cur_word = *lptr++;
      cur_word = cur_word & ((~cur_word) >> 1) & (*sample_include_quaterarr++);
      while (cur_word) {
        uii = CTZLU(cur_word) / 2;
        cur_write |= ONELU << (woffset + uii);
	cur_word &= cur_word - 1;
      }
      if (woffset) {
        *missing_bitfield++ = cur_write;
        cur_write = 0;
        woffset = 0;
      } else {
	woffset = BITCT2;
      }
    }
    if (woffset) {
      *missing_bitfield++ = cur_write;
    }
  } else {
    fill_ulong_zero(BITCT_TO_WORDCT(sample_ct), missing_bitfield);
    sample_idx = 0;
    for (widx = 0; sample_idx < sample_ct; widx++, lptr++) {
      cur_mask = *sample_include_quaterarr++;
      if (cur_mask) {
        cur_word = *lptr;
        cur_word = cur_word & ((~cur_word) >> 1) & cur_mask;
	if (cur_mask == FIVEMASK) {
          if (cur_word) {
	    uii = sample_idx;
            do {
              set_bit((CTZLU(cur_word) / 2) + uii, missing_bitfield);
              cur_word &= cur_word - 1;
	    } while (cur_word);
	  }
	  sample_idx += BITCT2;
	} else {
	  if (cur_word) {
	    do {
	      uii = CTZLU(cur_mask);
	      if ((cur_word >> uii) & 1) {
                set_bit_ul(sample_idx, missing_bitfield);
	      }
	      sample_idx++;
	      cur_mask &= cur_mask - 1;
	    } while (cur_mask);
	  } else {
            sample_idx += popcount2_long(cur_mask);
	  }
        }
      }
    }
  }
}

void hh_reset(unsigned char* loadbuf, uintptr_t* sample_include_quaterarr, uintptr_t unfiltered_sample_ct) {
  uintptr_t sample_bidx = 0;
  unsigned char* loadbuf_end = &(loadbuf[(unfiltered_sample_ct + 3) / 4]);
  unsigned char* iicp;
  unsigned char ucc;
  unsigned char ucc2;
  uintptr_t unfiltered_sample_ctd;
  uint32_t* loadbuf_alias32;
  uint32_t uii;
  uint32_t ujj;
#ifdef __LP64__
  uint32_t* sample_include_quaterarr_alias32;
  __m128i* loadbuf_alias;
  __m128i* iivp;
  __m128i vii;
  __m128i vjj;
  if (!(((uintptr_t)loadbuf) & 15)) {
    loadbuf_alias = (__m128i*)loadbuf;
    iivp = (__m128i*)sample_include_quaterarr;
    unfiltered_sample_ctd = unfiltered_sample_ct / 64;
    for (; sample_bidx < unfiltered_sample_ctd; sample_bidx++) {
      vii = *loadbuf_alias;
      vjj = _mm_and_si128(_mm_andnot_si128(vii, _mm_srli_epi64(vii, 1)), *iivp++);
      *loadbuf_alias++ = _mm_sub_epi64(vii, vjj);
    }
    loadbuf = (unsigned char*)loadbuf_alias;
    iicp = (unsigned char*)iivp;
  } else if (!(((uintptr_t)loadbuf) & 3)) {
    loadbuf_alias32 = (uint32_t*)loadbuf;
    sample_include_quaterarr_alias32 = (uint32_t*)sample_include_quaterarr;
    unfiltered_sample_ctd = unfiltered_sample_ct / BITCT2;
    for (; sample_bidx < unfiltered_sample_ctd; sample_bidx++) {
      uii = *loadbuf_alias32;
      ujj = ((uii >> 1) & (~uii)) & (*sample_include_quaterarr_alias32++);
      *loadbuf_alias32++ = uii - ujj;
    }
    loadbuf = (unsigned char*)loadbuf_alias32;
    iicp = (unsigned char*)sample_include_quaterarr_alias32;
  } else {
    iicp = (unsigned char*)sample_include_quaterarr;
  }
#else
  if (!(((uintptr_t)loadbuf) & 3)) {
    loadbuf_alias32 = (uint32_t*)loadbuf;
    unfiltered_sample_ctd = unfiltered_sample_ct / BITCT2;
    for (; sample_bidx < unfiltered_sample_ctd; sample_bidx++) {
      uii = *loadbuf_alias32;
      ujj = ((uii >> 1) & (~uii)) & (*sample_include_quaterarr++);
      *loadbuf_alias32++ = uii - ujj;
    }
    loadbuf = (unsigned char*)loadbuf_alias32;
  }
  iicp = (unsigned char*)sample_include_quaterarr;
#endif
  for (; loadbuf < loadbuf_end;) {
    ucc = *loadbuf;
    ucc2 = ((ucc >> 1) & (~ucc)) & (*iicp++);
    *loadbuf++ = ucc - ucc2;
  }
}

void hh_reset_y(unsigned char* loadbuf, uintptr_t* sample_include_quaterarr, uintptr_t* sample_male_include_quaterarr, uintptr_t unfiltered_sample_ct) {
  uintptr_t sample_bidx = 0;
  unsigned char* loadbuf_end = &(loadbuf[(unfiltered_sample_ct + 3) / 4]);
  unsigned char* iicp;
  unsigned char* imicp;
  unsigned char ucc;
  unsigned char ucc2;
  unsigned char ucc3;
  uintptr_t unfiltered_sample_ctd;
  uint32_t* loadbuf_alias32;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
#ifdef __LP64__
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  uint32_t* sample_include_quaterarr_alias32;
  uint32_t* sample_male_include_quaterarr_alias32;
  __m128i* loadbuf_alias;
  __m128i* iivp;
  __m128i* imivp;
  __m128i vii;
  __m128i vjj;
  __m128i vkk;
  if (!(((uintptr_t)loadbuf) & 15)) {
    loadbuf_alias = (__m128i*)loadbuf;
    iivp = (__m128i*)sample_include_quaterarr;
    imivp = (__m128i*)sample_male_include_quaterarr;
    unfiltered_sample_ctd = unfiltered_sample_ct / 64;
    for (; sample_bidx < unfiltered_sample_ctd; sample_bidx++) {
      // sample_include_quaterarr & ~sample_male_include_quaterarr: force to 01
      // sample_male_include_quaterarr: convert 10 to 01, keep everything else
      vii = *imivp++;
      vjj = *iivp++;
      vkk = _mm_and_si128(*loadbuf_alias, _mm_or_si128(vii, _mm_slli_epi64(vii, 1)));
      *loadbuf_alias++ = _mm_or_si128(_mm_andnot_si128(vii, vjj), _mm_sub_epi64(vkk, _mm_and_si128(_mm_andnot_si128(vkk, _mm_srli_epi64(vkk, 1)), m1)));
    }
    loadbuf = (unsigned char*)loadbuf_alias;
    iicp = (unsigned char*)iivp;
    imicp = (unsigned char*)imivp;
  } else if (!(((uintptr_t)loadbuf) & 3)) {
    loadbuf_alias32 = (uint32_t*)loadbuf;
    sample_include_quaterarr_alias32 = (uint32_t*)sample_include_quaterarr;
    sample_male_include_quaterarr_alias32 = (uint32_t*)sample_male_include_quaterarr;
    unfiltered_sample_ctd = unfiltered_sample_ct / 16;
    for (; sample_bidx < unfiltered_sample_ctd; sample_bidx++) {
      uii = *sample_male_include_quaterarr_alias32++;
      ujj = *sample_include_quaterarr_alias32++;
      ukk = (*loadbuf_alias32) & (uii * 3);
      *loadbuf_alias32++ = ((~uii) & ujj) | (ukk - ((~ukk) & (ukk >> 1) & 0x55555555));
    }
    loadbuf = (unsigned char*)loadbuf_alias32;
    iicp = (unsigned char*)sample_include_quaterarr_alias32;
    imicp = (unsigned char*)sample_male_include_quaterarr_alias32;
  } else {
    iicp = (unsigned char*)sample_include_quaterarr;
    imicp = (unsigned char*)sample_male_include_quaterarr;
  }
#else
  if (!(((uintptr_t)loadbuf) & 3)) {
    loadbuf_alias32 = (uint32_t*)loadbuf;
    unfiltered_sample_ctd = unfiltered_sample_ct / 16;
    for (; sample_bidx < unfiltered_sample_ctd; sample_bidx++) {
      uii = *sample_male_include_quaterarr++;
      ujj = *sample_include_quaterarr++;
      ukk = (*loadbuf_alias32) & (uii * 3);
      *loadbuf_alias32++ = ((~uii) & ujj) | (ukk - ((~ukk) & (ukk >> 1) & 0x55555555));
    }
    loadbuf = (unsigned char*)loadbuf_alias32;
  }
  iicp = (unsigned char*)sample_include_quaterarr;
  imicp = (unsigned char*)sample_male_include_quaterarr;
#endif
  for (; loadbuf < loadbuf_end;) {
    ucc = *imicp++;
    ucc2 = *iicp++;
    ucc3 = (*loadbuf) & (ucc * 3);
    *loadbuf++ = ((~ucc) & ucc2) | (ucc3 - ((~ucc3) & (ucc3 >> 1) & 0x55));
  }
}

uint32_t alloc_raw_haploid_filters(uint32_t unfiltered_sample_ct, uint32_t hh_exists, uint32_t is_include, uintptr_t* sample_bitarr, uintptr_t* sex_male, uintptr_t** sample_raw_include_quatervec_ptr, uintptr_t** sample_raw_male_include_quatervec_ptr) {
  uintptr_t unfiltered_sample_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(unfiltered_sample_ct);
  uintptr_t* sample_raw_male_include_quatervec;
  if (hh_exists & (Y_FIX_NEEDED | NXMHH_EXISTS)) {
    if (bigstack_alloc_ul(unfiltered_sample_ctv2, sample_raw_include_quatervec_ptr)) {
      return 1;
    }
    if (is_include) {
      init_quaterarr_from_bitarr(sample_bitarr, unfiltered_sample_ct, *sample_raw_include_quatervec_ptr);
    } else {
      init_quaterarr_from_inverted_bitarr(sample_bitarr, unfiltered_sample_ct, *sample_raw_include_quatervec_ptr);
    }
  }
  if (hh_exists & (XMHH_EXISTS | Y_FIX_NEEDED)) {
    if (bigstack_alloc_ul(unfiltered_sample_ctv2, sample_raw_male_include_quatervec_ptr)) {
      return 1;
    }
    sample_raw_male_include_quatervec = *sample_raw_male_include_quatervec_ptr;
    if (hh_exists & (Y_FIX_NEEDED | NXMHH_EXISTS)) {
      memcpy(sample_raw_male_include_quatervec, *sample_raw_include_quatervec_ptr, unfiltered_sample_ctv2 * sizeof(intptr_t));
    } else {
      if (is_include) {
	init_quaterarr_from_bitarr(sample_bitarr, unfiltered_sample_ct, sample_raw_male_include_quatervec);
      } else {
	init_quaterarr_from_inverted_bitarr(sample_bitarr, unfiltered_sample_ct, sample_raw_male_include_quatervec);
      }
    }
    apply_bitarr_mask_to_quaterarr_01(sex_male, unfiltered_sample_ct, sample_raw_male_include_quatervec);
  }
  return 0;
}

void haploid_fix_multiple(uintptr_t* marker_exclude, uintptr_t marker_uidx_start, uintptr_t marker_ct, Chrom_info* chrom_info_ptr, uint32_t hh_exists, uint32_t set_hh_missing, uint32_t set_mixed_mt_missing, uintptr_t* sample_raw_include2, uintptr_t* sample_raw_male_include2, uintptr_t unfiltered_sample_ct, uintptr_t byte_ct_per_marker, unsigned char* loadbuf) {
  uintptr_t marker_idx = 0;
  uintptr_t marker_uidx = next_unset_ul_unsafe(marker_exclude, marker_uidx_start);
  uint32_t chrom_fo_idx = get_variant_chrom_fo_idx(chrom_info_ptr, marker_uidx);
  uint32_t chrom_idx;
  uint32_t is_x;
  uint32_t is_y;
  uint32_t is_mt;
  uint32_t is_haploid;
  uintptr_t chrom_end;
  uintptr_t marker_idx_chrom_end;

  while (marker_idx < marker_ct) {
    chrom_idx = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
    chrom_end = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx + 1];
    is_x = (chrom_info_ptr->xymt_codes[X_OFFSET] == (int32_t)chrom_idx);
    is_y = (chrom_info_ptr->xymt_codes[Y_OFFSET] == (int32_t)chrom_idx);
    is_mt = (chrom_info_ptr->xymt_codes[MT_OFFSET] == (int32_t)chrom_idx);
    is_haploid = IS_SET(chrom_info_ptr->haploid_mask, chrom_idx);
    marker_idx_chrom_end = marker_idx + chrom_end - marker_uidx - popcount_bit_idx(marker_exclude, marker_uidx, chrom_end);
    if (marker_idx_chrom_end > marker_ct) {
      marker_idx_chrom_end = marker_ct;
    }
    if (is_haploid && set_hh_missing) {
      if (is_x) {
	if (hh_exists & XMHH_EXISTS) {
	  for (; marker_idx < marker_idx_chrom_end; marker_idx++) {
	    hh_reset(&(loadbuf[marker_idx * byte_ct_per_marker]), sample_raw_male_include2, unfiltered_sample_ct);
	  }
	}
      } else if (is_y) {
	if (hh_exists & Y_FIX_NEEDED) {
	  for (; marker_idx < marker_idx_chrom_end; marker_idx++) {
	    hh_reset_y(&(loadbuf[marker_idx * byte_ct_per_marker]), sample_raw_include2, sample_raw_male_include2, unfiltered_sample_ct);
	  }
	}
      } else if (hh_exists & NXMHH_EXISTS) {
	for (; marker_idx < marker_idx_chrom_end; marker_idx++) {
	  hh_reset(&(loadbuf[marker_idx * byte_ct_per_marker]), sample_raw_include2, unfiltered_sample_ct);
	}
      }
    } else if (is_mt && set_mixed_mt_missing) {
      for (; marker_idx < marker_idx_chrom_end; marker_idx++) {
	hh_reset(&(loadbuf[marker_idx * byte_ct_per_marker]), sample_raw_include2, unfiltered_sample_ct);
      }
    }
    marker_idx = marker_idx_chrom_end;
    chrom_fo_idx++;
  }
}

void force_missing(unsigned char* loadbuf, uintptr_t* force_missing_include2, uintptr_t unfiltered_sample_ct) {
  uintptr_t sample_bidx = 0;
  unsigned char* loadbuf_end = &(loadbuf[(unfiltered_sample_ct + 3) / 4]);
  unsigned char* fmicp;
  unsigned char ucc;
  unsigned char ucc2;
  uintptr_t unfiltered_sample_ctd;
  uint32_t* loadbuf_alias32;
  uint32_t uii;
  uint32_t ujj;
#ifdef __LP64__
  uint32_t* force_missing_include2_alias32;
  __m128i* loadbuf_alias;
  __m128i* fmivp;
  __m128i vii;
  __m128i vjj;
  if (!(((uintptr_t)loadbuf) & 15)) {
    loadbuf_alias = (__m128i*)loadbuf;
    fmivp = (__m128i*)force_missing_include2;
    unfiltered_sample_ctd = unfiltered_sample_ct / 64;
    for (; sample_bidx < unfiltered_sample_ctd; sample_bidx++) {
      vii = *loadbuf_alias;
      vjj = *fmivp++;
      vii = _mm_or_si128(vii, vjj);
      vjj = _mm_slli_epi64(vjj, 1);
      *loadbuf_alias++ = _mm_andnot_si128(vjj, vii);
    }
    loadbuf = (unsigned char*)loadbuf_alias;
    fmicp = (unsigned char*)fmivp;
  } else if (!(((uintptr_t)loadbuf) & 3)) {
    loadbuf_alias32 = (uint32_t*)loadbuf;
    force_missing_include2_alias32 = (uint32_t*)force_missing_include2;
    unfiltered_sample_ctd = unfiltered_sample_ct / BITCT2;
    for (; sample_bidx < unfiltered_sample_ctd; sample_bidx++) {
      uii = *loadbuf_alias32;
      ujj = *force_missing_include2_alias32++;
      uii |= ujj;
      ujj <<= 1;
      *loadbuf_alias32++ = uii & (~ujj);
    }
    loadbuf = (unsigned char*)loadbuf_alias32;
    fmicp = (unsigned char*)force_missing_include2_alias32;
  } else {
    fmicp = (unsigned char*)force_missing_include2;
  }
#else
  if (!(((uintptr_t)loadbuf) & 3)) {
    loadbuf_alias32 = (uint32_t*)loadbuf;
    unfiltered_sample_ctd = unfiltered_sample_ct / BITCT2;
    for (; sample_bidx < unfiltered_sample_ctd; sample_bidx++) {
      uii = *loadbuf_alias32;
      ujj = *force_missing_include2++;
      uii |= ujj;
      ujj <<= 1;
      *loadbuf_alias32++ = uii & (~ujj);
    }
    loadbuf = (unsigned char*)loadbuf_alias32;
  }
  fmicp = (unsigned char*)force_missing_include2;
#endif
  for (; loadbuf < loadbuf_end;) {
    ucc = *loadbuf;
    ucc2 = *fmicp++;
    ucc |= ucc2;
    ucc2 <<= 1;
    *loadbuf++ = ucc & (~ucc2);
  }
}

int32_t open_and_size_string_list(char* fname, FILE** infile_ptr, uintptr_t* list_len_ptr, uintptr_t* max_str_len_ptr) {
  // assumes file is not open yet, and g_textbuf is safe to clobber
  uint32_t max_len = 0;
  uintptr_t line_idx = 0;
  uintptr_t list_len = 0;
  int32_t retval = 0;
  char* bufptr;
  uint32_t cur_len;
  if (fopen_checked(fname, "r", infile_ptr)) {
    goto open_and_size_string_list_ret_OPEN_FAIL;
  }
  g_textbuf[MAXLINELEN - 1] = ' ';
  while (fgets(g_textbuf, MAXLINELEN, *infile_ptr)) {
    line_idx++;
    if (!g_textbuf[MAXLINELEN - 1]) {
      LOGERRPRINTFWW("Error: Line %" PRIuPTR " of %s is pathologically long.\n", line_idx, fname);
      goto open_and_size_string_list_ret_INVALID_FORMAT;
    }
    bufptr = skip_initial_spaces(g_textbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    // don't complain about more than one entry on a line for now
    list_len++;
    cur_len = strlen_se(bufptr);
    if (cur_len >= max_len) {
      max_len = cur_len + 1;
    }
  }
  if (!feof(*infile_ptr)) {
    goto open_and_size_string_list_ret_READ_FAIL;
  }
  *list_len_ptr = list_len;
  *max_str_len_ptr = max_len;
  while (0) {
  open_and_size_string_list_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  open_and_size_string_list_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  open_and_size_string_list_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  }
  return retval;
}

int32_t load_string_list(FILE** infile_ptr, uintptr_t max_str_len, char* str_list) {
  // assumes file is open (probably by open_and_size_string_list), and
  // g_textbuf is safe to clobber
  int32_t retval = 0;
  char* bufptr;
  uint32_t cur_len;
  rewind(*infile_ptr);
  while (fgets(g_textbuf, MAXLINELEN, *infile_ptr)) {
    bufptr = skip_initial_spaces(g_textbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    cur_len = strlen_se(bufptr);
    memcpy(str_list, bufptr, cur_len);
    str_list[cur_len] = '\0';
    str_list = &(str_list[max_str_len]);
  }
  if (!feof(*infile_ptr)) {
    goto load_string_list_ret_READ_FAIL;
  }
  while (0) {
  load_string_list_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  }
  return retval;
}

int32_t open_and_skip_first_lines(FILE** infile_ptr, char* fname, char* loadbuf, uintptr_t loadbuf_size, uint32_t lines_to_skip) {
  uint32_t line_idx;
  loadbuf[loadbuf_size - 1] = ' ';
  if (fopen_checked(fname, "r", infile_ptr)) {
    return RET_OPEN_FAIL;
  }
  for (line_idx = 1; line_idx <= lines_to_skip; line_idx++) {
    if (!fgets(loadbuf, loadbuf_size, *infile_ptr)) {
      if (feof(*infile_ptr)) {
	LOGERRPRINTFWW("Error: Fewer lines than expected in %s.\n", fname);
	return RET_INVALID_FORMAT;
      } else {
	return RET_READ_FAIL;
      }
    }
    if (!(loadbuf[loadbuf_size - 1])) {
      if ((loadbuf_size == MAXLINELEN) || (loadbuf_size == MAXLINEBUFLEN)) {
	LOGERRPRINTFWW("Error: Line %u of %s is pathologically long.\n", line_idx, fname);
	return RET_INVALID_FORMAT;
      } else {
        return RET_NOMEM;
      }
    }
  }
  return 0;
}

int32_t load_to_first_token(FILE* infile, uintptr_t loadbuf_size, char comment_char, const char* file_descrip, char* loadbuf, char** bufptr_ptr, uintptr_t* line_idx_ptr) {
  uintptr_t line_idx = 0;
  while (fgets(loadbuf, loadbuf_size, infile)) {
    line_idx++;
    if (!(loadbuf[loadbuf_size - 1])) {
      // PLINK 1.9 has two text line loading modes: "regular" and "long".
      // * "Regular" mode limits lines to about MAXLINELEN (about 128k as of
      //   this writing) characters.
      // * "Long" mode theoretically accepts lines up to about MAXLINEBUFLEN
      //   (~2 GB) characters but degrades gracefully if less memory is
      //   available (in that case, an out-of-memory instead of an
      //   invalid-format error is reported on fgets overflow).  Any long
      //   buffer size larger than MAXLINELEN should work properly with
      //   plink_common.
      if ((loadbuf_size == MAXLINELEN) || (loadbuf_size == MAXLINEBUFLEN)) {
	LOGERRPRINTF("Error: Line %" PRIuPTR " of %s is pathologically long.\n", line_idx, file_descrip);
	return RET_INVALID_FORMAT;
      } else {
	return RET_NOMEM;
      }
    }
    *bufptr_ptr = skip_initial_spaces(loadbuf);
    if (!is_eoln_kns(**bufptr_ptr)) {
      if ((**bufptr_ptr) != comment_char) {
	*line_idx_ptr = line_idx;
        return 0;
      }
    }
  }
  if (!feof(infile)) {
    return RET_READ_FAIL;
  }
  LOGERRPRINTF("Error: Empty %s.\n", file_descrip);
  return RET_INVALID_FORMAT;
}

int32_t open_and_load_to_first_token(FILE** infile_ptr, char* fname, uintptr_t loadbuf_size, char comment_char, const char* file_descrip, char* loadbuf, char** bufptr_ptr, uintptr_t* line_idx_ptr) {
  loadbuf[loadbuf_size - 1] = ' ';
  if (fopen_checked(fname, "r", infile_ptr)) {
    return RET_OPEN_FAIL;
  }
  return load_to_first_token(*infile_ptr, loadbuf_size, comment_char, file_descrip, loadbuf, bufptr_ptr, line_idx_ptr);
}

int32_t scan_max_strlen(char* fname, uint32_t colnum, uint32_t colnum2, uint32_t headerskip, char skipchar, uintptr_t* max_str_len_ptr, uintptr_t* max_str2_len_ptr) {
  // colnum and colnum2 are 1-based indices.  If colnum2 is zero, only colnum
  // is scanned.
  // Includes terminating null in lengths.
  FILE* infile = nullptr;
  uintptr_t loadbuf_size = bigstack_left();
  uintptr_t max_str_len = *max_str_len_ptr;
  uintptr_t max_str2_len = 0;
  char* loadbuf = (char*)g_bigstack_base;
  uint32_t colmin;
  uint32_t coldiff;
  char* str1_ptr;
  char* str2_ptr;
  char cc;
  uintptr_t cur_str_len;
  uintptr_t line_idx;
  int32_t retval;
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  } else if (loadbuf_size <= MAXLINELEN) {
    goto scan_max_strlen_ret_NOMEM;
  }
  retval = open_and_skip_first_lines(&infile, fname, loadbuf, loadbuf_size, headerskip);
  if (retval) {
    goto scan_max_strlen_ret_1;
  }
  if (colnum < colnum2) {
    max_str2_len = *max_str2_len_ptr;
    colmin = colnum - 1;
    coldiff = colnum2 - colnum;
  } else if (colnum2) {
    max_str2_len = max_str_len;
    max_str_len = *max_str2_len_ptr;
    colmin = colnum2 - 1;
    coldiff = colnum - colnum2;
  } else {
    colmin = colnum - 1;
    coldiff = 0;
    colnum2 = 0xffffffffU;
  }
  line_idx = headerskip;
  while (fgets(loadbuf, loadbuf_size, infile)) {
    line_idx++;
    if (!(loadbuf[loadbuf_size - 1])) {
      if (loadbuf_size == MAXLINEBUFLEN) {
        LOGPREPRINTFWW("Error: Line %" PRIuPTR " of %s is pathologically long.\n", line_idx, fname);
	goto scan_max_strlen_ret_INVALID_FORMAT_2;
      } else {
        goto scan_max_strlen_ret_NOMEM;
      }
    }
    str1_ptr = skip_initial_spaces(loadbuf);
    cc = *str1_ptr;
    if (is_eoln_kns(cc) || (cc == skipchar)) {
      continue;
    }
    str1_ptr = next_token_multz(str1_ptr, colmin);
    str2_ptr = next_token_multz(str1_ptr, coldiff);
    if (no_more_tokens_kns(str2_ptr)) {
      // probably want option for letting this slide in the future
      LOGPREPRINTFWW("Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, fname);
      goto scan_max_strlen_ret_INVALID_FORMAT_2;
    }
    cur_str_len = strlen_se(str1_ptr);
    if (cur_str_len >= max_str_len) {
      max_str_len = cur_str_len + 1;
    }
    if (coldiff) {
      cur_str_len = strlen_se(str2_ptr);
      if (cur_str_len >= max_str2_len) {
	max_str2_len = cur_str_len + 1;
      }
    }
  }
  if (!feof(infile)) {
    goto scan_max_strlen_ret_READ_FAIL;
  }
  if (colnum < colnum2) {
    *max_str_len_ptr = max_str_len;
    if (coldiff) {
      *max_str2_len_ptr = max_str2_len;
    }
  } else {
    *max_str_len_ptr = max_str2_len;
    *max_str2_len_ptr = max_str_len;
  }
  while (0) {
  scan_max_strlen_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  scan_max_strlen_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  scan_max_strlen_ret_INVALID_FORMAT_2:
    logerrprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
 scan_max_strlen_ret_1:
  fclose_cond(infile);
  return retval;
}

int32_t scan_max_fam_indiv_strlen(char* fname, uint32_t colnum, uintptr_t* max_sample_id_len_ptr) {
  // colnum is a 1-based index with the FID column number; IID column is
  // assumed to follow.
  // Includes terminating null in lengths.
  FILE* infile = nullptr;
  uintptr_t loadbuf_size = bigstack_left();
  uintptr_t max_sample_id_len = *max_sample_id_len_ptr;
  uintptr_t line_idx = 0;
  char* loadbuf = (char*)g_bigstack_base;
  char* bufptr;
  char* bufptr2;
  uintptr_t cur_sample_id_len;
  int32_t retval;
  colnum--;
  if (loadbuf_size > MAXLINEBUFLEN) {
    loadbuf_size = MAXLINEBUFLEN;
  } else if (loadbuf_size <= MAXLINELEN) {
    goto scan_max_fam_indiv_strlen_ret_NOMEM;
  }
  retval = open_and_skip_first_lines(&infile, fname, loadbuf, loadbuf_size, 0);
  if (retval) {
    goto scan_max_fam_indiv_strlen_ret_1;
  }
  while (fgets(loadbuf, loadbuf_size, infile)) {
    line_idx++;
    if (!(loadbuf[loadbuf_size - 1])) {
      if (loadbuf_size == MAXLINEBUFLEN) {
	LOGPREPRINTFWW("Error: Line %" PRIuPTR " of %s is pathologically long.\n", line_idx, fname);
	goto scan_max_fam_indiv_strlen_ret_INVALID_FORMAT_2;
      } else {
        goto scan_max_fam_indiv_strlen_ret_NOMEM;
      }
    }
    bufptr = skip_initial_spaces(loadbuf);
    if (is_eoln_kns(*bufptr)) {
      continue;
    }
    bufptr = next_token_multz(bufptr, colnum);
    bufptr2 = next_token(bufptr);
    if (no_more_tokens_kns(bufptr2)) {
      LOGPREPRINTFWW("Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, fname);
      goto scan_max_fam_indiv_strlen_ret_INVALID_FORMAT_2;
    }
    cur_sample_id_len = strlen_se(bufptr) + strlen_se(bufptr2) + 2;
    if (cur_sample_id_len > max_sample_id_len) {
      max_sample_id_len = cur_sample_id_len;
    }
  }
  if (!feof(infile)) {
    goto scan_max_fam_indiv_strlen_ret_READ_FAIL;
  }
  *max_sample_id_len_ptr = max_sample_id_len;
  while (0) {
  scan_max_fam_indiv_strlen_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  scan_max_fam_indiv_strlen_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  scan_max_fam_indiv_strlen_ret_INVALID_FORMAT_2:
    logerrprintb();
    retval = RET_INVALID_FORMAT;
    break;
  }
 scan_max_fam_indiv_strlen_ret_1:
  fclose_cond(infile);
  return retval;
}

/*
void inplace_collapse_uint32(uint32_t* item_arr, uint32_t unfiltered_ct, uintptr_t* exclude_arr, uint32_t filtered_ct) {
  if (unfiltered_ct == filtered_ct) {
    return;
  }
  uint32_t item_uidx = next_set_unsafe(exclude_arr, 0);
  uint32_t item_idx = item_uidx;
  for (; item_idx < filtered_ct; item_idx++, item_uidx++) {
    next_unset_unsafe_ck(exclude_arr, &item_uidx);
    item_arr[item_idx] = item_arr[item_uidx];
  }
}
*/

void inplace_collapse_uint32_incl(uint32_t* item_arr, uint32_t unfiltered_ct, uintptr_t* incl_arr, uint32_t filtered_ct) {
  if (unfiltered_ct == filtered_ct) {
    return;
  }
  uint32_t item_uidx = next_unset_unsafe(incl_arr, 0);
  uint32_t item_idx = item_uidx;
  for (; item_idx < filtered_ct; item_idx++, item_uidx++) {
    next_set_unsafe_ck(incl_arr, &item_uidx);
    item_arr[item_idx] = item_arr[item_uidx];
  }
}

char* alloc_and_init_collapsed_arr(char* item_arr, uintptr_t item_len, uintptr_t unfiltered_ct, uintptr_t* exclude_arr, uintptr_t filtered_ct, uint32_t read_only) {
  uint32_t item_uidx = 0;
  char* new_arr;
  char* wptr;
  char* wptr_end;
  uintptr_t item_uidx_stop;
  uintptr_t delta;
  if (read_only && (unfiltered_ct == filtered_ct)) {
    return item_arr;
  }
  if (bigstack_alloc_c(filtered_ct * item_len, &new_arr)) {
    return nullptr;
  }
  wptr = new_arr;
  wptr_end = &(new_arr[filtered_ct * item_len]);
  while (wptr < wptr_end) {
    item_uidx = next_unset_ul_unsafe(exclude_arr, item_uidx);
    item_uidx_stop = next_set_ul(exclude_arr, item_uidx, unfiltered_ct);
    delta = item_uidx_stop - item_uidx;
    memcpy(wptr, &(item_arr[item_uidx * item_len]), delta * item_len);
    wptr = &(wptr[delta * item_len]);
    item_uidx = item_uidx_stop;
  }
  return new_arr;
}

char* alloc_and_init_collapsed_arr_incl(char* item_arr, uintptr_t item_len, uintptr_t unfiltered_ct, uintptr_t* include_arr, uintptr_t filtered_ct, uint32_t read_only) {
  uint32_t item_uidx = 0;
  char* new_arr;
  char* wptr;
  char* wptr_end;
  uintptr_t item_uidx_stop;
  uintptr_t delta;
  if (read_only && (unfiltered_ct == filtered_ct)) {
    return item_arr;
  }
  if (bigstack_alloc_c(filtered_ct * item_len, &new_arr)) {
    return nullptr;
  }
  wptr = new_arr;
  wptr_end = &(new_arr[filtered_ct * item_len]);
  do {
    item_uidx = next_set_ul_unsafe(include_arr, item_uidx);
    item_uidx_stop = next_unset_ul(include_arr, item_uidx, unfiltered_ct);
    delta = item_uidx_stop - item_uidx;
    memcpy(wptr, &(item_arr[item_uidx * item_len]), delta * item_len);
    wptr = &(wptr[delta * item_len]);
    item_uidx = item_uidx_stop;
  } while (wptr < wptr_end);
  return new_arr;
}

void inplace_delta_collapse_arr(char* item_arr, uintptr_t item_len, uintptr_t filtered_ct_orig, uintptr_t filtered_ct_new, uintptr_t* exclude_orig, uintptr_t* exclude_new) {
  // if this sort of collapse function is ever in an important loop, check
  // whether specialized 4-byte and 8-byte versions are much faster
  uintptr_t* exclude_orig_start = exclude_orig;
  char* write_end = &(item_arr[filtered_ct_new * item_len]);
  uintptr_t read_idx = 1;
  uint32_t uii = 0;
  char* write_ptr;
  uintptr_t ulii;
  uintptr_t uljj;
  uint32_t read_uidx;
  uint32_t ujj;
  if (filtered_ct_new == filtered_ct_orig) {
    return;
  }
  // find location of first newly excluded item
  while (1) {
    ulii = *exclude_orig;
    uljj = *exclude_new;
    if (ulii != uljj) {
      break;
    }
    uii += popcount_long(ulii);
    exclude_orig++;
    exclude_new++;
  }
  exclude_new -= ((uintptr_t)(exclude_orig - exclude_orig_start));
  read_uidx = BITCT * ((uintptr_t)(exclude_orig - exclude_orig_start));
  ujj = CTZLU(ulii ^ uljj);
  read_uidx += ujj;
  uii += popcount_long(ulii & ((ONELU << ujj) - ONELU));
  uii = read_uidx - uii; // now equal to # initial filtered indices skipped
  filtered_ct_new -= uii;
  item_arr = &(item_arr[uii * item_len]);
  write_ptr = item_arr;
  read_uidx++;
  for (; write_ptr < write_end; read_uidx++, read_idx++) {
    next_unset_unsafe_ck(exclude_orig_start, &read_uidx);
    if (IS_SET(exclude_new, read_uidx)) {
      continue;
    }
    memcpy(write_ptr, &(item_arr[read_idx * item_len]), item_len);
    write_ptr = &(write_ptr[item_len]);
  }
}

void inplace_delta_collapse_bitfield(uintptr_t* read_ptr, uint32_t filtered_ct_new, uintptr_t* exclude_orig, uintptr_t* exclude_new) {
  // only guaranteed to zero out trailing bits up to the nearest 16-byte
  // boundary on 64-bit systems
  uintptr_t* write_ptr = read_ptr;
  uintptr_t readw = *read_ptr++;
  uintptr_t writew = 0;
  uint32_t item_uidx = 0;
  uint32_t item_mwidx = 0;
  uint32_t item_idx = 0;
  for (; item_idx < filtered_ct_new; item_uidx++) {
    next_unset_unsafe_ck(exclude_orig, &item_uidx);
    if (!is_set(exclude_new, item_uidx)) {
      if ((readw >> item_mwidx) & 1) {
	writew |= ONELU << (item_idx % BITCT);
      }
      if (!((++item_idx) % BITCT)) {
	*write_ptr++ = writew;
	writew = 0;
      }
    }
    if (++item_mwidx == BITCT) {
      item_mwidx = 0;
      readw = *read_ptr++;
    }
  }
  if (write_ptr < read_ptr) {
    *write_ptr++ = writew;
    if (write_ptr < read_ptr) {
      *write_ptr = 0;
    }
  }
}

void copy_bitarr_subset_excl(const uintptr_t* __restrict raw_bitarr, const uintptr_t* __restrict subset_excl, uint32_t raw_bitarr_size, uint32_t subset_size, uintptr_t* __restrict output_bitarr) {
  uintptr_t cur_write = 0;
  uint32_t item_uidx = 0;
  uint32_t write_bit = 0;
  uint32_t item_idx = 0;
  uint32_t item_uidx_stop;
  if (!subset_excl[0]) {
    item_uidx = next_set(subset_excl, 0, raw_bitarr_size & (~(BITCT - 1))) & (~(BITCT - 1));
    memcpy(output_bitarr, raw_bitarr, item_uidx / 8);
    item_idx = item_uidx;
    output_bitarr = &(output_bitarr[item_uidx / BITCT]);
  }
  while (item_idx < subset_size) {
    item_uidx = next_unset_unsafe(subset_excl, item_uidx);
    item_uidx_stop = next_set(subset_excl, item_uidx, raw_bitarr_size);
    item_idx += item_uidx_stop - item_uidx;
    do {
      cur_write |= ((raw_bitarr[item_uidx / BITCT] >> (item_uidx % BITCT)) & 1) << write_bit;
      if (++write_bit == BITCT) {
	*output_bitarr++ = cur_write;
        cur_write = 0;
	write_bit = 0;
      }
    } while (++item_uidx < item_uidx_stop);
  }
  if (write_bit) {
    *output_bitarr = cur_write;
  }
}

void copy_bitarr_subset(const uintptr_t* __restrict raw_bitarr, const uintptr_t* __restrict subset_mask, uint32_t raw_bitarr_size, uint32_t subset_size, uintptr_t* __restrict output_bitarr) {
  // full-blown blocked copy not worth it due to undefined CTZLU(0), >> 64,
  // << 64
  uintptr_t cur_output_word = 0;
  uint32_t item_uidx = 0;
  uint32_t word_write_shift = 0;
  uint32_t item_idx = 0;
  uint32_t item_uidx_stop;
  if (!(~subset_mask[0])) {
    item_uidx = next_unset(subset_mask, 0, raw_bitarr_size & (~(BITCT - 1))) & (~(BITCT - 1));
    memcpy(output_bitarr, raw_bitarr, item_uidx / 8);
    item_idx = item_uidx;
    output_bitarr = &(output_bitarr[item_uidx / BITCT]);
  }
  while (item_idx < subset_size) {
    item_uidx = next_set_unsafe(subset_mask, item_uidx);

    // can speed this up a bit once we have a guaranteed unset bit at the end
    item_uidx_stop = next_unset(subset_mask, item_uidx, raw_bitarr_size);

    item_idx += item_uidx_stop - item_uidx;
    do {
      cur_output_word |= ((raw_bitarr[item_uidx / BITCT] >> (item_uidx % BITCT)) & 1) << word_write_shift;
      if (++word_write_shift == BITCT) {
	*output_bitarr++ = cur_output_word;
        cur_output_word = 0;
	word_write_shift = 0;
      }
    } while (++item_uidx < item_uidx_stop);
  }
  if (word_write_shift) {
    *output_bitarr = cur_output_word;
  }
}

void uncollapse_copy_flip_include_arr(uintptr_t* collapsed_include_arr, uintptr_t unfiltered_ct, uintptr_t* exclude_arr, uintptr_t* output_exclude_arr) {
  uintptr_t unfiltered_ctl = BITCT_TO_WORDCT(unfiltered_ct);
  uintptr_t* output_exclude_true_end = &(output_exclude_arr[unfiltered_ctl]);
  uintptr_t* output_exclude_end = &(output_exclude_arr[unfiltered_ct / BITCT]);
  uintptr_t cea_read = 0;
  uint32_t read_bit = BITCT;
  uint32_t write_bit;
  uintptr_t cur_write;
  uintptr_t cur_read = 0;
  if (!exclude_arr[0]) {
    // copy-with-possible-offset is substantially slower, so treat initial lack
    // of offset as a special case
    for (cur_read = 0; cur_read < unfiltered_ctl; cur_read++) {
      *output_exclude_arr++ = ~(*collapsed_include_arr++);
      if (*(++exclude_arr)) {
	break;
      }
    }
  }
  while (output_exclude_arr < output_exclude_end) {
    cur_write = *exclude_arr++;
    // want efficient handling of all-zeroes and all-ones here
    if (cur_write) {
      cur_read = ~cur_write;
    uncollapse_copy_flip_include_arr_loop:
      while (cur_read) {
        write_bit = CTZLU(cur_read);
        if (read_bit == BITCT) {
          cea_read = ~(*collapsed_include_arr++);
	  read_bit = 0;
        }
        cur_write |= (cea_read & ONELU) << write_bit;
        cea_read >>= 1;
        read_bit++;
        cur_read &= cur_read - ONELU;
      }
      *output_exclude_arr = cur_write;
    } else {
      if (read_bit == BITCT) {
        *output_exclude_arr = ~(*collapsed_include_arr++);
      } else {
        cur_write = cea_read;
        cea_read = ~(*collapsed_include_arr++);
        *output_exclude_arr = cur_write | (cea_read << (BITCT - read_bit));
	cea_read >>= read_bit;
      }
    }
    output_exclude_arr++;
  }
  if (output_exclude_arr < output_exclude_true_end) {
    cur_write = *exclude_arr++;
    cur_read = (~cur_write) & ((ONELU << (unfiltered_ct % BITCT)) - ONELU);
    goto uncollapse_copy_flip_include_arr_loop;
  }
}

void copy_when_nonmissing(uintptr_t* loadbuf, char* source, uintptr_t elem_size, uintptr_t unfiltered_sample_ct, uintptr_t missing_ct, char* dest) {
  uintptr_t* loadbuf_end = &(loadbuf[QUATERCT_TO_WORDCT(unfiltered_sample_ct)]);
  uintptr_t last_missing_p1 = 0;
  uintptr_t sample_idx_offset = 0;
  uintptr_t cur_word;
  uintptr_t new_missing_idx;
  uintptr_t diff;
  if (!missing_ct) {
    memcpy(dest, source, unfiltered_sample_ct * elem_size);
    return;
  }
  do {
    cur_word = *loadbuf++;
    cur_word = cur_word & (~(cur_word >> 1)) & FIVEMASK;
    while (cur_word) {
      new_missing_idx = sample_idx_offset + (CTZLU(cur_word) / 2);
      diff = new_missing_idx - last_missing_p1;
      if (diff) {
	dest = memcpya(dest, &(source[last_missing_p1 * elem_size]), diff * elem_size);
      }
      last_missing_p1 = new_missing_idx + 1;
      cur_word &= cur_word - 1;
    }
    sample_idx_offset += BITCT2;
  } while (loadbuf < loadbuf_end);
  diff = unfiltered_sample_ct - last_missing_p1;
  if (diff) {
    memcpy(dest, &(source[last_missing_p1 * elem_size]), diff * elem_size);
  }
}

uint32_t collapse_duplicate_ids(char* sorted_ids, uintptr_t id_ct, uintptr_t max_id_len, uint32_t* id_starts) {
  // Collapses array of sorted IDs to remove duplicates, and writes
  // pre-collapse positions to id_starts (so e.g. duplication count of any
  // sample ID can be determined via subtraction) if it isn't nullptr.
  // Returns id_ct of collapsed array.
  uintptr_t read_idx;
  uintptr_t write_idx;
  if (!id_ct) {
    return 0;
  }
  if (id_starts) {
    id_starts[0] = 0;
    for (read_idx = 1; read_idx < id_ct; read_idx++) {
      if (!strcmp(&(sorted_ids[(read_idx - 1) * max_id_len]), &(sorted_ids[read_idx * max_id_len]))) {
	break;
      }
      id_starts[read_idx] = read_idx;
    }
    write_idx = read_idx;
    while (++read_idx < id_ct) {
      if (strcmp(&(sorted_ids[(write_idx - 1) * max_id_len]), &(sorted_ids[read_idx * max_id_len]))) {
	strcpy(&(sorted_ids[write_idx * max_id_len]), &(sorted_ids[read_idx * max_id_len]));
	id_starts[write_idx++] = read_idx;
      }
    }
  } else {
    for (read_idx = 1; read_idx < id_ct; read_idx++) {
      if (!strcmp(&(sorted_ids[(read_idx - 1) * max_id_len]), &(sorted_ids[read_idx * max_id_len]))) {
	break;
      }
    }
    write_idx = read_idx;
    while (++read_idx < id_ct) {
      if (strcmp(&(sorted_ids[(write_idx - 1) * max_id_len]), &(sorted_ids[read_idx * max_id_len]))) {
	strcpy(&(sorted_ids[write_idx * max_id_len]), &(sorted_ids[read_idx * max_id_len]));
	write_idx++;
      }
    }
  }
  return write_idx;
}

void range_list_init(Range_list* range_list_ptr) {
  range_list_ptr->names = nullptr;
  range_list_ptr->starts_range = nullptr;
  range_list_ptr->name_ct = 0;
  range_list_ptr->name_max_len = 0;
}

void free_range_list(Range_list* range_list_ptr) {
  free_cond(range_list_ptr->names);
  free_cond(range_list_ptr->starts_range);
}

// implementation used in PLINK 1.07 stats.cpp
// probably want to remove this function and use erf() calls in the future
double normdist(double zz) {
  double sqrt2pi = 2.50662827463;
  double t0;
  double z1;
  double p0;
  t0 = 1 / (1 + 0.2316419 * fabs(zz));
  z1 = exp(-0.5 * zz * zz) / sqrt2pi;
  p0 = z1 * t0 * (0.31938153 + t0 * (-0.356563782 + t0 * (1.781477937 + t0 * (-1.821255978 + 1.330274429 * t0))));
  return zz >= 0 ? 1 - p0 : p0;
}

double rand_normal(double* secondval_ptr) {
  // N(0, 1)
  double dxx = sqrt(-2 * log(rand_unif()));
  double dyy = 2 * PI * rand_unif();
  *secondval_ptr = dxx * cos(dyy);
  return dxx * sin(dyy);
}

void init_sfmt64_from_sfmt32(sfmt_t* sfmt32, sfmt_t* sfmt64) {
  // sfmt_genrand_uint64() is not supposed to be called after
  // sfmt_genrand_uint32() is called on the same generator.  To work around
  // this, we initialize a new sfmt64 generator with this function when
  // necessary, and stick to genrand_uint32() calls with the main generator.
  uint32_t init_arr[4];
  uint32_t uii;
  for (uii = 0; uii < 4; uii++) {
    init_arr[uii] = sfmt_genrand_uint32(sfmt32);
  }
  sfmt_init_by_array(sfmt64, init_arr, 4);
}

void generate_perm1_interleaved(uint32_t tot_ct, uint32_t set_ct, uintptr_t perm_idx, uintptr_t perm_ct, uintptr_t* perm_buf) {
  uintptr_t tot_ctl = BITCT_TO_WORDCT(tot_ct);
  uintptr_t tot_rem = tot_ct & (BITCT - 1);
  uint32_t tot_quotient = (uint32_t)(0x100000000LLU / tot_ct);
  uint32_t upper_bound = tot_ct * tot_quotient - 1;
  uintptr_t uljj = perm_ct - perm_idx;
  uint32_t totq_preshift;
  uint64_t totq_magic;
  uint32_t totq_postshift;
  uint32_t totq_incr;
  uintptr_t* pbptr;
  uint32_t num_set;
  uint32_t urand;
  uintptr_t ulii;
  // seeing as how we're gonna divide by the same number a billion times or so,
  // it just might be worth optimizing that division...
  magic_num(tot_quotient, &totq_magic, &totq_preshift, &totq_postshift, &totq_incr);
  if (set_ct * 2 < tot_ct) {
    for (ulii = 0; ulii < tot_ctl; ulii++) {
      fill_ulong_zero(uljj, &(perm_buf[perm_idx + (ulii * perm_ct)]));
    }
    for (; perm_idx < perm_ct; perm_idx++) {
      pbptr = &(perm_buf[perm_idx]);
      for (num_set = 0; num_set < set_ct; num_set++) {
	do {
	  do {
	    urand = sfmt_genrand_uint32(&g_sfmt);
	  } while (urand > upper_bound);
	  // this is identical to ulii = urand / tot_quotient
	  ulii = (totq_magic * ((urand >> totq_preshift) + totq_incr)) >> totq_postshift;
	  uljj = ulii / BITCT;
	  ulii &= (BITCT - 1);
	} while ((pbptr[uljj * perm_ct] >> ulii) & 1);
	pbptr[uljj * perm_ct] |= (ONELU << ulii);
      }
    }
  } else {
    for (ulii = 0; ulii < tot_ctl; ulii++) {
      fill_ulong_one(uljj, &(perm_buf[perm_idx + (ulii * perm_ct)]));
    }
    // "set" has reversed meaning here
    set_ct = tot_ct - set_ct;
    for (; perm_idx < perm_ct; perm_idx++) {
      pbptr = &(perm_buf[perm_idx]);
      for (num_set = 0; num_set < set_ct; num_set++) {
	do {
	  do {
	    urand = sfmt_genrand_uint32(&g_sfmt);
	  } while (urand > upper_bound);
	  ulii = (totq_magic * ((urand >> totq_preshift) + totq_incr)) >> totq_postshift;
	  uljj = ulii / BITCT;
	  ulii &= (BITCT - 1);
	} while (!((pbptr[uljj * perm_ct] >> ulii) & 1));
	pbptr[uljj * perm_ct] &= ~(ONELU << ulii);
      }
    }
    if (tot_rem) {
      uljj = (~ZEROLU) >> (BITCT - tot_rem);
      pbptr = &(perm_buf[(tot_ctl - 1) * perm_ct + perm_idx]);
      for (ulii = perm_idx; ulii < perm_ct; ulii++) {
	*pbptr &= uljj;
	pbptr++;
      }
    }
  }
}

uint32_t cubic_real_roots(double coef_a, double coef_b, double coef_c, double* solutions) {
  // Analytically finds all real roots of x^3 + ax^2 + bx + c, saving them in
  // solutions[] (sorted from smallest to largest), and returning the count.
  // Multiple roots are only returned/counted once.
  // Additional research into numerical stability may be in order here.
  double a2 = coef_a * coef_a;
  double qq = (a2 - 3 * coef_b) * (1.0 / 9.0);
  double rr = (2 * a2 * coef_a - 9 * coef_a * coef_b + 27 * coef_c) * (1.0 / 54.0);
  double r2 = rr * rr;
  double q3 = qq * qq * qq;
  double adiv3 = coef_a * (1.0 / 3.0);
  double sq;
  double dxx;
  if (r2 < q3) {
    // three real roots
    sq = sqrt(qq);
    dxx = acos(rr / (qq * sq)) * (1.0 / 3.0);
    sq *= -2;
    solutions[0] = sq * cos(dxx) - adiv3;
    solutions[1] = sq * cos(dxx + (2.0 * PI / 3.0)) - adiv3;
    solutions[2] = sq * cos(dxx - (2.0 * PI / 3.0)) - adiv3;
    // now sort and check for within-epsilon equality
    if (solutions[0] > solutions[1]) {
      dxx = solutions[0];
      solutions[0] = solutions[1];
      if (dxx > solutions[2]) {
        solutions[1] = solutions[2];
	solutions[2] = dxx;
      } else {
	solutions[1] = dxx;
      }
      if (solutions[0] > solutions[1]) {
	dxx = solutions[0];
	solutions[0] = solutions[1];
	solutions[1] = dxx;
      }
    } else if (solutions[1] > solutions[2]) {
      dxx = solutions[1];
      solutions[1] = solutions[2];
      solutions[2] = dxx;
    }
    if (solutions[1] - solutions[0] < EPSILON) {
      solutions[1] = solutions[2];
      return (solutions[1] - solutions[0] < EPSILON)? 1 : 2;
    }
    return (solutions[2] - solutions[1] < EPSILON)? 2 : 3;
  }
  dxx = -pow(fabs(rr) + sqrt(r2 - q3), 1.0 / 3.0);
  if (dxx == 0.0) {
    solutions[0] = -adiv3;
    return 1;
  }
  if (rr < 0.0) {
    dxx = -dxx;
  }
  sq = qq / dxx;
  solutions[0] = dxx + sq - adiv3;
  // use of regular epsilon here has actually burned us
  if (fabs(dxx - sq) >= (EPSILON * 8)) {
    return 1;
  }
  if (dxx >= 0.0) {
    solutions[1] = solutions[0];
    solutions[0] = -dxx - adiv3;
  } else {
    solutions[1] = -dxx - adiv3;
  }
  return 2;
}

void join_threads(pthread_t* threads, uint32_t ctp1) {
  if (!(--ctp1)) {
    return;
  }
#ifdef _WIN32
  WaitForMultipleObjects(ctp1, threads, 1, INFINITE);
  for (uint32_t uii = 0; uii < ctp1; ++uii) {
    CloseHandle(threads[uii]);
  }
#else
  for (uint32_t uii = 0; uii < ctp1; uii++) {
    pthread_join(threads[uii], nullptr);
  }
#endif
}

#ifdef _WIN32
int32_t spawn_threads(pthread_t* threads, unsigned (__stdcall *start_routine)(void*), uintptr_t ct)
#else
int32_t spawn_threads(pthread_t* threads, void* (*start_routine)(void*), uintptr_t ct)
#endif
{
  uintptr_t ulii;
  if (ct == 1) {
    return 0;
  }
  for (ulii = 1; ulii < ct; ulii++) {
#ifdef _WIN32
    threads[ulii - 1] = (HANDLE)_beginthreadex(nullptr, 4096, start_routine, (void*)ulii, 0, nullptr);
    if (!threads[ulii - 1]) {
      join_threads(threads, ulii);
      return -1;
    }
#else
    if (pthread_create(&(threads[ulii - 1]), nullptr, start_routine, (void*)ulii)) {
      join_threads(threads, ulii);
      return -1;
    }
#endif
  }
  return 0;
}

// Okay, it's time to bite the bullet and stop creating and destroying threads
// like crazy, at least in the small-block-size GRM calculation; Intel
// MKL-powered GCTA 1.24 blew away our code on the NIH 512-core test machine
// when the maximum number of threads was used.  Mostly because threads were
// actually costing much more in creation/destruction time than they saved;
// much better wall-clock times would have resulted from manually setting
// --threads to a low number.  That's not cool.
//
// New framework:
// * On all operating systems, g_is_last_thread_block indicates whether all
//   threads should terminate upon completion of the current block.  (Initially
//   had this volatile, then realized that the presence of the sync-wait should
//   be enough to force the global variable to be reread.)
// * On Linux and OS X, if we aren't dealing with the final block,
//   spawn_threads2() also reinitializes g_thread_active_ct.
// * On Linux and OS X, spawn_threads2() checks if g_thread_mutex_initialized
//   is set.  If not, it, it is set, g_thread_sync_mutex,
//   g_thread_cur_block_done_condvar and g_thread_start_next_condvar are
//   initialized, then threads are launched.
//   If it has, pthread_cond_broadcast() acts on g_thread_start_next_condvar.
// * On Windows, spawn_threads2() checks if g_thread_mutex_initialized is set.
//   If it has not, it, along with g_thread_start_next_event[] and
//   g_thread_cur_block_done_events[], are initialized, then the threads are
//   launched.  If it has, SetEvent() acts on g_thread_start_next_event[].
//   (It used to act on only one event; then I realized that safely dealing
//   with a manual-reset event could be a pain if the first thread finishes
//   before the last one wakes up...)
// * Thread functions are expected to be of the form
//     THREAD_RET_TYPE function_name(void* arg) {
//       uintptr_t tidx = (uintptr_t)arg;
//       ...
//       while (1) {
//         ... // process current block
//         if ((!tidx) || g_is_last_thread_block) {
//           THREAD_RETURN;
//         }
//         THREAD_BLOCK_FINISH(tidx);
//       }
//     }
// * On Linux and OS X, THREAD_BLOCK_FINISH() acquires a mutex, decrements
//   g_thread_active_ct, calls pthread_cond_signal() on
//   g_thread_cur_block_done_condvar iff g_thread_active_ct is now zero, then
//   unconditionally calls pthread_cond_wait on g_thread_start_next_condvar and
//   the mutex.
// * On Windows, THREAD_BLOCK_FINISH() calls SetEvent() on
//   g_thread_cur_block_done_events[tidx - 1], then waits on
//   g_thread_start_next_event[tidx - 1].
// * If the termination variable is set, join_threads2() waits for all threads
//   to complete, then cleans up all multithreading objects.  Otherwise, on
//   Linux and OS X, it acquires the mutex and calls pthread_cond_wait() on
//   g_thread_cur_block_done_condvar and the mutex; and on Windows, it calls
//   WaitForMultipleObjects() on g_thread_cur_block_done_events[].
//   WaitForMultipleObjects has a 64 object limit, and for now it doesn't seem
//   too important to use a for loop to handle more objects?... well, we can
//   add that if anyone wants it, but for now the Windows thread limit is 65
//   (the main thread isn't part of the wait).
//
// This is only very slightly better than the original approach on my old
// MacBook Pro (since threading overhead was never high to begin with, there
// being only 2 cores...), but the impact should be more noticeable on heavily
// multicore machines.
//
// The next performance improvement to make is double-buffering; tricky to
// estimate how much (if any) "consumption" the main I/O thread should be
// doing, though, so it may want a job queue to go with it.

uintptr_t g_thread_spawn_ct;
uint32_t g_is_last_thread_block = 0;
#ifdef _WIN32
HANDLE g_thread_start_next_event[MAX_THREADS];
HANDLE g_thread_cur_block_done_events[MAX_THREADS];
#else
static pthread_mutex_t g_thread_sync_mutex;
static pthread_cond_t g_thread_cur_block_done_condvar;
static pthread_cond_t g_thread_start_next_condvar;
uint32_t g_thread_active_ct;

void THREAD_BLOCK_FINISH(uintptr_t tidx) {
  uintptr_t initial_spawn_ct = g_thread_spawn_ct;
  pthread_mutex_lock(&g_thread_sync_mutex);
  if (!(--g_thread_active_ct)) {
    pthread_cond_signal(&g_thread_cur_block_done_condvar);
  }
  while (g_thread_spawn_ct == initial_spawn_ct) {
    // spurious wakeup guard
    pthread_cond_wait(&g_thread_start_next_condvar, &g_thread_sync_mutex);
  }
  pthread_mutex_unlock(&g_thread_sync_mutex);
}
#endif
static uint32_t g_thread_mutex_initialized = 0;

void join_threads2(pthread_t* threads, uint32_t ctp1, uint32_t is_last_block) {
  uint32_t uii;
  if (!(--ctp1)) {
    if (is_last_block) {
      // allow another multithreaded function to be called later
      g_thread_mutex_initialized = 0;
    }
    return;
  }
#ifdef _WIN32
  if (!is_last_block) {
    WaitForMultipleObjects(ctp1, g_thread_cur_block_done_events, 1, INFINITE);
  } else {
    WaitForMultipleObjects(ctp1, threads, 1, INFINITE);
    for (uii = 0; uii < ctp1; uii++) {
      CloseHandle(threads[uii]);
      CloseHandle(g_thread_start_next_event[uii]);
      CloseHandle(g_thread_cur_block_done_events[uii]);
    }
    g_thread_mutex_initialized = 0;
  }
#else
  if (!is_last_block) {
    pthread_mutex_lock(&g_thread_sync_mutex);
    while (g_thread_active_ct) {
      pthread_cond_wait(&g_thread_cur_block_done_condvar, &g_thread_sync_mutex);
    }
    // keep mutex until next block loaded
  } else {
    for (uii = 0; uii < ctp1; uii++) {
      pthread_join(threads[uii], nullptr);
    }
    // slightly inefficient if there are multiple multithreaded commands being
    // run, but if different commands require different numbers of threads,
    // optimizing this sort of thing away could introduce bugs...
    pthread_mutex_destroy(&g_thread_sync_mutex);
    pthread_cond_destroy(&g_thread_cur_block_done_condvar);
    pthread_cond_destroy(&g_thread_start_next_condvar);
    g_thread_mutex_initialized = 0;
  }
#endif
}

#ifdef _WIN32
int32_t spawn_threads2(pthread_t* threads, unsigned (__stdcall *start_routine)(void*), uintptr_t ct, uint32_t is_last_block)
#else
int32_t spawn_threads2(pthread_t* threads, void* (*start_routine)(void*), uintptr_t ct, uint32_t is_last_block)
#endif
{
  uintptr_t ulii;
  // this needs to go before the ct == 1 check since start_routine() might need
  // it
  if (g_is_last_thread_block != is_last_block) {
    // might save us an unnecessary memory write that confuses the cache
    // coherency logic?
    g_is_last_thread_block = is_last_block;
  }
#ifdef _WIN32
  if (!g_thread_mutex_initialized) {
    g_thread_spawn_ct = 0;
    g_thread_mutex_initialized = 1;
    if (ct == 1) {
      return 0;
    }
    for (ulii = 1; ulii < ct; ulii++) {
      g_thread_start_next_event[ulii - 1] = CreateEvent(nullptr, FALSE, FALSE, nullptr);
      g_thread_cur_block_done_events[ulii - 1] = CreateEvent(nullptr, FALSE, FALSE, nullptr);
    }
    for (ulii = 1; ulii < ct; ulii++) {
      threads[ulii - 1] = (HANDLE)_beginthreadex(nullptr, 4096, start_routine, (void*)ulii, 0, nullptr);
      if (!threads[ulii - 1]) {
	if (ulii > 1) {
	  join_threads2(threads, ulii, is_last_block);
	  if (!is_last_block) {
	    for (uintptr_t uljj = 0; uljj < ulii - 1; ++uljj) {
	      CloseHandle(threads[uljj]);
	    }
	  }
	}
	if ((!is_last_block) || (ulii == 1)) {
	  for (uint32_t uii = 0; uii < ct - 1; ++uii) {
	    CloseHandle(g_thread_start_next_event[uii]);
	    CloseHandle(g_thread_cur_block_done_events[uii]);
	  }
	  g_thread_mutex_initialized = 0;
	}
	return -1;
      }
    }
  } else {
    g_thread_spawn_ct++;
    for (ulii = 1; ulii < ct; ulii++) {
      SetEvent(g_thread_start_next_event[ulii - 1]);
    }
  }
#else
  if (!is_last_block) {
    g_thread_active_ct = ct - 1;
  }
  if (!g_thread_mutex_initialized) {
    g_thread_spawn_ct = 0; // tidx 0 may need to know modulus
    g_thread_mutex_initialized = 1;
    if (ct == 1) {
      return 0;
    }
    if (pthread_mutex_init(&g_thread_sync_mutex, nullptr) ||
        pthread_cond_init(&g_thread_cur_block_done_condvar, nullptr) ||
        pthread_cond_init(&g_thread_start_next_condvar, nullptr)) {
      return -1;
    }
    for (ulii = 1; ulii < ct; ulii++) {
      if (pthread_create(&(threads[ulii - 1]), nullptr, start_routine, (void*)ulii)) {
	if (ulii > 1) {
	  join_threads2(threads, ulii, is_last_block);
	  if (!is_last_block) {
	    for (uintptr_t uljj = 0; uljj < ulii - 1; ++uljj) {
	      pthread_cancel(threads[uljj]);
	    }
	  }
	}
	if ((!is_last_block) || (ulii == 1)) {
	  pthread_mutex_destroy(&g_thread_sync_mutex);
	  pthread_cond_destroy(&g_thread_cur_block_done_condvar);
	  pthread_cond_destroy(&g_thread_start_next_condvar);
	  g_thread_mutex_initialized = 0;
	}
	return -1;
      }
    }
  } else {
    g_thread_spawn_ct++;
    if (ct == 1) {
      return 0;
    }
    // still holding mutex
    pthread_mutex_unlock(&g_thread_sync_mutex);
    pthread_cond_broadcast(&g_thread_start_next_condvar);
  }
#endif
  return 0;
}

sfmt_t** g_sfmtp_arr;

uint32_t bigstack_init_sfmtp(uint32_t thread_ct) {
  uint32_t uibuf[4];
  uint32_t tidx;
  uint32_t uii;
  g_sfmtp_arr = (sfmt_t**)bigstack_alloc(thread_ct * sizeof(intptr_t));
  if (!g_sfmtp_arr) {
    return 1;
  }
  g_sfmtp_arr[0] = &g_sfmt;
  if (thread_ct > 1) {
    for (tidx = 1; tidx < thread_ct; tidx++) {
      g_sfmtp_arr[tidx] = (sfmt_t*)bigstack_alloc(sizeof(sfmt_t));
      if (!g_sfmtp_arr[tidx]) {
	return 1;
      }
      for (uii = 0; uii < 4; uii++) {
	uibuf[uii] = sfmt_genrand_uint32(&g_sfmt);
      }
      sfmt_init_by_array(g_sfmtp_arr[tidx], uibuf, 4);
    }
  }
  return 0;
}
