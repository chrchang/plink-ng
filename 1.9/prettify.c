#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <unistd.h>

#ifdef _WIN32
#include <windows.h>
#endif

// Simple pretty printer.  (Christopher Chang, chrchang@alumni.caltech.edu)
//
// Unix "expand", "unexpand", and OS X "tab2space" have similar functions.
// Sadly, this cannot use Unix pipes since the entire input file must be
// scanned to determine column widths before the first output line can be
// written (presumably this is the main reason there is no preexisting Unix
// utility for this exact function).

#define RET_HELP 1
#define RET_NOMEM 2
#define RET_OPEN_FAIL 3
#define RET_INVALID_CMDLINE 4
#define RET_READ_FAIL 5
#define RET_INVALID_FORMAT 6
#define RET_WRITE_FAIL 7

const char errstr_fopen[] = "Error: Failed to open %s.\n";

#define FLAG_INPLACE 1
#define FLAG_RJUSTIFY 2
#define FLAG_SPACES_BEFORE_FIRST 4
#define FLAG_PAD 8
#define FLAG_SPACES_AFTER_LAST 0x10
#define FLAG_FINAL_EOLN 0x20
#define FLAG_STRIP_BLANK 0x40

#define FNAMESIZE 4096

char pathbuf[FNAMESIZE * 2 + 128];

void disp_usage(FILE* stream) {
  fputs(
"  prettify [flag(s)...] <input filename> [output filename]\n\n"
"  -i, --inplace      : Replace the input instead of writing to a new file.\n"
"  -s, --spacing <ct> : Set number of spaces between columns (default 2).\n"
"  -r, --ralign       : Make right sides of columns line up, instead of left.\n"
"  -l, --leading      : Add space(s) before the first column.\n"
"  -e, --extend-short : Use spaces to extend lines with fewer columns.\n"
"  -t, --trailing     : Add space(s) after the last column.\n"
"  -f, --force-eoln   : Force last line to be terminated by a newline.\n"
"  -n, --noblank      : Remove blank lines.\n\n"
"If no output filename is provided (and --inplace isn't in effect), results are\n"
"dumped to standard output.\n"
, stream);
}

void dispmsg(int32_t retval) {
  switch (retval) {
  case RET_NOMEM:
    fputs("\nError: Out of memory.\n", stderr);
    break;
  case RET_READ_FAIL:
    fputs("\nError: File read failure.\n", stderr);
    break;
  case RET_WRITE_FAIL:
    fputs("\nError: File write failure.\n", stderr);
    break;
  }
}

void free_cond(void* memptr) {
  if (memptr) {
    free(memptr);
  }
}

int32_t fopen_checked(FILE** target_ptr, const char* fname, const char* mode) {
  *target_ptr = fopen(fname, mode);
  if (!(*target_ptr)) {
    fprintf(stderr, errstr_fopen, fname);
    return -1;
  }
  return 0;
}

void fclose_cond(FILE* fptr) {
  if (fptr) {
    fclose(fptr);
  }
}

int32_t fclose_null(FILE** fptr_ptr) {
  int32_t ii = ferror(*fptr_ptr);
  int32_t jj = fclose(*fptr_ptr);
  *fptr_ptr = NULL;
  return ii || jj;
}

static inline void fill_ulong_zero(uintptr_t* ularr, size_t size) {
  uintptr_t* ulptr = &(ularr[size]);
  while (ularr < ulptr) {
    *ularr++ = 0;
  }
}

static inline uint32_t skip_spaces_ck(unsigned char** str_ptr, unsigned char* buf_end) {
  unsigned char* ss;
  for (ss = *str_ptr; ss < buf_end; ss++) {
    // conveniently, this treats the \r in \r\n as a space
    if ((*ss) > ' ') {
      *str_ptr = ss;
      return 0;
    }
  }
  return 1;
}

static inline unsigned char* get_token_end_ck(unsigned char* token_start, unsigned char* buf_end) {
  // assumes previous character was in a token, and does not support nonspace
  // characters below ASCII 33
  while ((token_start != buf_end) && ((*token_start) > ' ')) {
    token_start++;
  }
  return token_start;
}

#define BUFSIZE 131072
#define INITIAL_COLS 65536

unsigned char g_readbuf[BUFSIZE];

void handle_last_column(uintptr_t* col_widths, uintptr_t cur_col_idx, uintptr_t cur_col_width, uintptr_t* col_ct_ptr) {
  if (cur_col_width > col_widths[cur_col_idx]) {
    col_widths[cur_col_idx] = cur_col_width;
  }
  if (cur_col_idx > (*col_ct_ptr)) {
    *col_ct_ptr = cur_col_idx;
  }
}

int32_t scan_column_widths(FILE* infile, uintptr_t column_sep, uintptr_t** col_widths_ptr, uintptr_t* col_ct_ptr, unsigned char** spacebuf_ptr, unsigned char** rjustify_buf_ptr) {
  uintptr_t malloc_size = INITIAL_COLS * sizeof(intptr_t);
  uintptr_t* col_widths = (uintptr_t*)malloc(malloc_size);
  uintptr_t col_ct = 0;
  uintptr_t max_col_ct = INITIAL_COLS; // not a hard limit

  // actually a one-based index, to simplify distinguishing between
  // beginning-of-line-and-not-in-column from beginning-of-line-and-in-column
  // first element of col_widths[] is essentially unused as a result
  uintptr_t cur_col_idx = 0;

  uintptr_t cur_col_width = 0;
  uintptr_t line_idx = 1;
  int32_t retval = 0;
  unsigned char* readptr;
  unsigned char* line_end;
  unsigned char* readbuf_end;
  unsigned char* token_end;
  uintptr_t* new_col_widths;
  uintptr_t cur_read;
  if (!col_widths) {
    goto scan_column_widths_ret_NOMEM;
  }
  fill_ulong_zero(col_widths, max_col_ct);
  cur_read = fread(g_readbuf, 1, BUFSIZE, infile);
  if (ferror(infile)) {
    goto scan_column_widths_ret_READ_FAIL;
  }
  if (!cur_read) {
    goto scan_column_widths_ret_INVALID_FORMAT;
  }
  readptr = g_readbuf;
  readbuf_end = &(g_readbuf[cur_read]);
  while (1) {
    line_end = (unsigned char*)memchr(readptr, '\n', (uintptr_t)(readbuf_end - readptr));
    if (!line_end) {
      line_end = readbuf_end;
    }
    while (readptr < line_end) {
      if (!cur_col_width) {
	if (skip_spaces_ck(&readptr, line_end)) {
	  break;
	}
	if (++cur_col_idx == max_col_ct) {
	  malloc_size *= 2;
	  new_col_widths = (uintptr_t*)realloc(col_widths, malloc_size);
	  if (!new_col_widths) {
	    goto scan_column_widths_ret_READ_FAIL;
	  }
          col_widths = new_col_widths;
	  fill_ulong_zero(&(col_widths[max_col_ct]), max_col_ct);
	  max_col_ct *= 2;
	}
      }
      token_end = get_token_end_ck(readptr, line_end);
      cur_col_width += (uintptr_t)(token_end - readptr);
      if (token_end == line_end) {
	break;
      }
      if (cur_col_width > col_widths[cur_col_idx]) {
	col_widths[cur_col_idx] = cur_col_width;
      }
      cur_col_width = 0;
      readptr = token_end;
    }
    if ((line_end < readbuf_end) || (!cur_read)) {
      handle_last_column(col_widths, cur_col_idx, cur_col_width, &col_ct);
      if (line_end == readbuf_end) {
	// EOF
	break;
      }
      readptr = &(line_end[1]);
      line_idx++;
      cur_col_idx = 0;
      cur_col_width = 0;
      continue;
    }
    // in middle of line
    cur_read = fread(g_readbuf, 1, BUFSIZE, infile);
    if (ferror(infile)) {
      goto scan_column_widths_ret_READ_FAIL;
    }
    readptr = g_readbuf;
    readbuf_end = &(g_readbuf[cur_read]);
  }
  cur_col_width = 0;
  for (cur_col_idx = 1; cur_col_idx <= col_ct; cur_col_idx++) {
    if (col_widths[cur_col_idx] > cur_col_width) {
      cur_col_width = col_widths[cur_col_idx];
    }
  }
  *spacebuf_ptr = (unsigned char*)malloc(cur_col_width + column_sep);
  if (!(*spacebuf_ptr)) {
    goto scan_column_widths_ret_NOMEM;
  }
  memset(*spacebuf_ptr, 32, cur_col_width + column_sep);
  if (rjustify_buf_ptr) {
    *rjustify_buf_ptr = (unsigned char*)malloc(cur_col_width);
    if (!(*rjustify_buf_ptr)) {
      goto scan_column_widths_ret_NOMEM;
    }
  }
  while (0) {
  scan_column_widths_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  scan_column_widths_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  scan_column_widths_ret_INVALID_FORMAT:
    fputs("Error: Empty input file.\n", stderr);
    retval = RET_INVALID_FORMAT;
    break;
  }
  *col_widths_ptr = col_widths;
  *col_ct_ptr = col_ct;
  return retval;
}

int32_t pretty_write(FILE* infile, char* outname, uint32_t flags, uintptr_t column_sep, uintptr_t* col_widths, uintptr_t col_ct, unsigned char* spacebuf, unsigned char* rjustify_buf) {
  FILE* outfile = NULL;
  uintptr_t cur_col_idx = 0;
  uintptr_t cur_col_width = 0;
  uintptr_t prev_col_width = 0;
  uintptr_t rjbuf_len = 0;
  unsigned char* token_end = NULL;
  uint32_t spaces_before_first = flags & FLAG_SPACES_BEFORE_FIRST;
  uint32_t pad = flags & FLAG_PAD;
  uint32_t spaces_after_last = flags & FLAG_SPACES_AFTER_LAST;
  uint32_t final_eoln = flags & FLAG_FINAL_EOLN;
  uint32_t strip_blank = flags & FLAG_STRIP_BLANK;
  uint32_t no_final_newline = 0;
  int32_t retval = 0;
  unsigned char* readptr;
  unsigned char* line_end;
  unsigned char* readbuf_end;
  uintptr_t cur_read;
  if (!outname) {
    outfile = stdout;
  } else {
    if (fopen_checked(&outfile, outname, "w")) {
      goto pretty_write_ret_OPEN_FAIL;
    }
  }
  rewind(infile);
  cur_read = fread(g_readbuf, 1, BUFSIZE, infile);
  if (ferror(infile)) {
    goto pretty_write_ret_READ_FAIL;
  }
  readptr = g_readbuf;
  readbuf_end = &(g_readbuf[cur_read]);
  while (1) {
    line_end = (unsigned char*)memchr(readptr, '\n', (uintptr_t)(readbuf_end - readptr));
    if (!line_end) {
      if (readptr != readbuf_end) {
        no_final_newline = 1;
      }
      line_end = readbuf_end;
    }
    while (readptr < line_end) {
      if (!cur_col_width) {
	if (skip_spaces_ck(&readptr, line_end)) {
	  break;
	}
	if (cur_col_idx || spaces_before_first) {
	  if (!rjustify_buf) {
	    fwrite(spacebuf, 1, col_widths[cur_col_idx] + column_sep - prev_col_width, outfile);
	  } else {
	    fwrite(spacebuf, 1, column_sep, outfile);
	  }
	}
	cur_col_idx++;
      }
      token_end = get_token_end_ck(readptr, line_end);
      cur_col_width += (uintptr_t)(token_end - readptr);
      if (token_end == line_end) {
	break;
      }
      if (!rjustify_buf) {
        fwrite(readptr, 1, token_end - readptr, outfile);
        prev_col_width = cur_col_width;
      } else {
        fwrite(spacebuf, 1, col_widths[cur_col_idx] - cur_col_width, outfile);
	if (rjbuf_len) {
          fwrite(rjustify_buf, 1, rjbuf_len, outfile);
	  rjbuf_len = 0;
	}
        fwrite(readptr, 1, token_end - readptr, outfile);
      }
      cur_col_width = 0;
      readptr = token_end;
    }
    if ((line_end < readbuf_end) || (!cur_read)) {
      if (cur_col_idx) {
	if (!rjustify_buf) {
	  if (cur_col_width) {
	    // last column not dumped yet
	    fwrite(readptr, 1, token_end - readptr, outfile);
	    prev_col_width = cur_col_width;
	  }
	  if (pad || spaces_after_last) {
	    fwrite(spacebuf, 1, col_widths[cur_col_idx] - prev_col_width, outfile);
	  }
	} else {
	  fwrite(spacebuf, 1, col_widths[cur_col_idx] - cur_col_width, outfile);
	  if (rjbuf_len) {
	    fwrite(rjustify_buf, 1, rjbuf_len, outfile);
	    rjbuf_len = 0;
	  }
	  fwrite(readptr, 1, token_end - readptr, outfile);
	}
	if (pad) {
	  while (cur_col_idx < col_ct) {
	    fwrite(spacebuf, 1, col_widths[++cur_col_idx] + column_sep, outfile);
	  }
	}
	if (spaces_after_last) {
	  fwrite(spacebuf, 1, column_sep, outfile);
	}
      }
      if (!cur_read) {
	// EOF
	if (final_eoln && no_final_newline) {
	  putc('\n', outfile);
	}
	break;
      }
      if (cur_col_idx || (!strip_blank)) {
	putc('\n', outfile);
	if (ferror(outfile)) {
	  goto pretty_write_ret_WRITE_FAIL;
	}
      }
      readptr = &(line_end[1]);
      cur_col_idx = 0;
      cur_col_width = 0;
      prev_col_width = 0;
      continue;
    }
    // in middle of line
    if (cur_col_width) {
      if (!rjustify_buf) {
	fwrite(readptr, 1, token_end - readptr, outfile);
      } else {
	memcpy(&(rjustify_buf[rjbuf_len]), readptr, token_end - readptr);
	rjbuf_len += (uintptr_t)(token_end - readptr);
      }
    }

    cur_read = fread(g_readbuf, 1, BUFSIZE, infile);
    if (ferror(infile)) {
      goto pretty_write_ret_READ_FAIL;
    }
    if (cur_read) {
      no_final_newline = 0;
    }
    readptr = g_readbuf;
    readbuf_end = &(g_readbuf[cur_read]);
  }

  if (outname) {
    if (fclose_null(&outfile)) {
      goto pretty_write_ret_WRITE_FAIL;
    }
  }
  while (0) {
  pretty_write_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  pretty_write_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  pretty_write_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  if (outname) {
    fclose_cond(outfile);
  }
  return retval;
}

int32_t main(int32_t argc, char** argv) {
  FILE* infile = NULL;
  char* outname = NULL;
  uintptr_t column_sep = 2;
  uint32_t flags = 0;
  int32_t retval = 0;
  uintptr_t* col_widths = NULL;
  unsigned char* spacebuf = NULL;
  unsigned char* rjustify_buf = NULL;
  uintptr_t col_ct = 0;
  uint32_t infile_param_idx = 0;
  char* param_ptr;
#ifndef _WIN32
  char* cptr;
#endif
  uint32_t param_idx;
  uint32_t uii;
  int32_t ii;
  char cc;
  if (argc == 1) {
    goto main_ret_HELP;
  }
  for (param_idx = 1; param_idx < (uint32_t)argc; param_idx++) {
    if ((!strcmp(argv[param_idx], "--help")) || (!strcmp(argv[param_idx], "-help")) || (!strcmp(argv[param_idx], "-?")) || (!strcmp(argv[param_idx], "-h"))) {
      goto main_ret_HELP;
    }
  }

  if (argc > 10) {
    fputs("Error: Too many parameters.\n\n", stderr);
    goto main_ret_INVALID_CMDLINE_2;
  }
  for (param_idx = 1; param_idx < (uint32_t)argc; param_idx++) {
    if (argv[param_idx][0] != '-') {
      if (!infile_param_idx) {
	infile_param_idx = param_idx;
      } else if (!outname) {
	if (flags & FLAG_INPLACE) {
	  goto main_ret_INVALID_CMDLINE_3;
	}
        outname = argv[param_idx];
      } else {
	fputs("Error: Invalid parameter sequence.\n\n", stderr);
	goto main_ret_INVALID_CMDLINE_2;
      }
      continue;
    }
    param_ptr = &(argv[param_idx][1]);
    if (*param_ptr == '-') {
      // allow both single- and double-dash
      param_ptr++;
    }
    if (!strcmp(param_ptr, "inplace")) {
      if (outname) {
	goto main_ret_INVALID_CMDLINE_3;
      }
      flags |= FLAG_INPLACE;
    } else if ((!strcmp(param_ptr, "spacing")) || (!strcmp(param_ptr, "s"))) {
      if (++param_idx == (uint32_t)argc) {
	fputs("Error: Missing --spacing parameter.\n", stderr);
	goto main_ret_INVALID_CMDLINE;
      }
      ii = atoi(argv[param_idx]);
      if (ii < 1) {
	fprintf(stderr, "Error: Invalid --spacing parameter '%s'.\n", argv[param_idx]);
	goto main_ret_INVALID_CMDLINE;
      }
      column_sep = (uint32_t)ii;
    } else if (!strcmp(param_ptr, "ralign")) {
      flags |= FLAG_RJUSTIFY;
    } else if (!strcmp(param_ptr, "leading")) {
      flags |= FLAG_SPACES_BEFORE_FIRST;
    } else if (!strcmp(param_ptr, "extend-short")) {
      flags |= FLAG_PAD;
    } else if (!strcmp(param_ptr, "trailing")) {
      flags |= FLAG_SPACES_AFTER_LAST;
    } else if (!strcmp(param_ptr, "force-eoln")) {
      flags |= FLAG_FINAL_EOLN;
    } else if (!strcmp(param_ptr, "noblank")) {
      flags |= FLAG_STRIP_BLANK;
    } else {
      if ((argv[param_idx][1] != '-') && argv[param_idx][1]) {
	// permit abbreviated style
	while (1) {
	  cc = *param_ptr++;
	  if (!cc) {
	    break;
	  }
	  switch (cc) {
	  case 'i':
	    if (outname) {
	      goto main_ret_INVALID_CMDLINE_3;
	    }
	    flags |= FLAG_INPLACE;
	    break;
	  case 'r':
	    flags |= FLAG_RJUSTIFY;
	    break;
	  case 'l':
	    flags |= FLAG_SPACES_BEFORE_FIRST;
	    break;
	  case 'e':
	    flags |= FLAG_PAD;
	    break;
	  case 't':
	    flags |= FLAG_SPACES_AFTER_LAST;
	    break;
	  case 'f':
	    flags |= FLAG_FINAL_EOLN;
	    break;
	  case 'n':
	    flags |= FLAG_STRIP_BLANK;
	    break;
	  default:
            fprintf(stderr, "Error: Invalid flag '%s'.\n\n", argv[param_idx]);
	    goto main_ret_INVALID_CMDLINE_2;
	  }
	}
      } else {
	fprintf(stderr, "Error: Invalid flag '%s'.\n\n", argv[param_idx]);
	goto main_ret_INVALID_CMDLINE_2;
      }
    }
  }
  if (!infile_param_idx) {
    fputs("Error: No input filename.\n\n", stderr);
    goto main_ret_INVALID_CMDLINE_2;
  }
  if (flags & FLAG_INPLACE) {
    uii = strlen(argv[infile_param_idx]);
    outname = (char*)malloc(uii + 11);
    if (!outname) {
      goto main_ret_NOMEM;
    }
    memcpy(outname, argv[infile_param_idx], uii);
    memcpy(&(outname[uii]), "-temporary", 11);
  } else if (outname) {
#ifdef _WIN32
    uii = GetFullPathName(argv[infile_param_idx], FNAMESIZE, pathbuf, NULL);
    if ((!uii) || (uii > FNAMESIZE))
#else
    if (!realpath(argv[infile_param_idx], pathbuf))
#endif
    {
      fprintf(stderr, "Error: Failed to open %s.\n", argv[infile_param_idx]);
      goto main_ret_OPEN_FAIL;
    }
#ifdef _WIN32
    uii = GetFullPathName(outname, FNAMESIZE, &(pathbuf[FNAMESIZE + 64]), NULL);
    if (uii && (uii <= FNAMESIZE) && (!strcmp(pathbuf, &(pathbuf[FNAMESIZE + 64]))))
#else
    cptr = realpath(outname, &(pathbuf[FNAMESIZE + 64]));
    if (cptr && (!strcmp(pathbuf, &(pathbuf[FNAMESIZE + 64]))))
#endif
    {
      fputs("Error: Input and output files match.  Use --inplace instead.\n", stderr);
      goto main_ret_INVALID_CMDLINE;
    }
  }
  if (fopen_checked(&infile, argv[infile_param_idx], "rb")) {
    goto main_ret_OPEN_FAIL;
  }
  retval = scan_column_widths(infile, column_sep, &col_widths, &col_ct, &spacebuf, (flags & FLAG_RJUSTIFY)? (&rjustify_buf) : NULL);
  if (retval) {
    goto main_ret_1;
  }
  retval = pretty_write(infile, outname, flags, column_sep, col_widths, col_ct, spacebuf, rjustify_buf);
  if (retval) {
    goto main_ret_1;
  }
  fclose_null(&infile);
  if (flags & FLAG_INPLACE) {
    unlink(argv[infile_param_idx]);
    if (rename(outname, argv[infile_param_idx])) {
      fprintf(stderr, "Error: File rename failed.  Output is in %s instead of %s.\n", outname, argv[infile_param_idx]);
      goto main_ret_OPEN_FAIL;
    }
  }
  while (0) {
  main_ret_HELP:
    fputs(
"prettify v1.05 (5 Mar 2019)    Christopher Chang (chrchang@alumni.caltech.edu)\n\n"
"Takes a tab-and/or-space-delimited text table, and generates a space-delimited\n"
"pretty-printed version.  Multibyte character encodings are not currently\n"
"supported.\n\n"
, stdout);
    disp_usage(stdout);
    fputs(
"\nTo perform the simplest reverse conversion (multiple spaces to one tab), you\n"
"can use\n"
"  cat <input filename> | tr -s ' ' '\\t' > <output filename>\n"
"For one-to-one conversion between spaces and tabs instead, omit the \"-s\".  And\n"
"to strip leading and trailing tabs and spaces, try\n"
"  cat <in> | sed 's/^[[:space:]]*//g' | sed 's/[[:space:]]*$//g' > <out>\n"
, stdout);
    retval = RET_HELP;
    break;
  main_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  main_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  main_ret_INVALID_CMDLINE_3:
    fputs("Error: --inplace cannot be used with an output filename.\n", stderr);
    retval = RET_INVALID_CMDLINE;
    break;
  main_ret_INVALID_CMDLINE_2:
    disp_usage(stderr);
  main_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  }
 main_ret_1:
  free_cond(col_widths);
  fclose_cond(infile);
  dispmsg(retval);
  return retval;
}
