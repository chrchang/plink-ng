#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>

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

void disp_usage(FILE* stream) {
  fputs(
"  prettify {flag(s)...} [input filename] {output filename}\n\n"
"  --spacing [ct] : Change number of spaces between columns (default 2).\n"
"  --first        : Add spaces before the first column.\n"
"  --last         : Add spaces after the last column.\n"
"  --final-eoln   : Force last line to be terminated by a newline.\n"
"  --strip-blank  : Remove blank lines.\n\n"
"If no output filename is provided, standard output is used.\n"
, stream);
}

void free_cond(void* memptr) {
  if (memptr) {
    free(memptr);
  }
}

int32_t fopen_checked(FILE** target_ptr, const char* fname, const char* mode) {
  *target_ptr = fopen(fname, mode);
  if (!(*target_ptr)) {
    printf(errstr_fopen, fname);
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
  unsigned char ucc;
  for (ss = *str_ptr; ss < buf_end; ss++) {
    ucc = *ss;
    if ((ucc != ' ') && (ucc != '\t')) {
      *str_ptr = ss;
      return 0;
    }
  }
  return 1;
}

static inline unsigned char* get_token_end_ck(unsigned char* token_start, unsigned char* buf_end) {
  // assumes we are currently in a token, and does not support nonspace
  // characters below ASCII 33
  do {
    token_start++;
  } while ((token_start != buf_end) && ((*token_start) > ' '));
  return token_start;
}

#define BUFSIZE 131072
#define INITIAL_COLS 65536

unsigned char g_readbuf[BUFSIZE];

uint32_t malloc_double(uintptr_t* malloc_size_ptr, unsigned char** old_ptr_ptr) {
  uintptr_t malloc_size = *malloc_size_ptr;
  unsigned char* new_ptr;
  new_ptr = (unsigned char*)malloc(malloc_size * 2);
  if (!new_ptr) {
    return 1;
  }
  memcpy(new_ptr, *old_ptr_ptr, malloc_size);
  free(*old_ptr_ptr);
  *old_ptr_ptr = new_ptr;
  return 0;
}

void handle_last_column(uintptr_t* col_widths, uintptr_t cur_col_idx, uintptr_t cur_col_width, uintptr_t* col_ct_ptr) {
  if (cur_col_width > col_widths[cur_col_idx]) {
    col_widths[cur_col_idx] = cur_col_width;
  }
  if (cur_col_idx > (*col_ct_ptr)) {
    *col_ct_ptr = cur_col_idx;
  }
}

int32_t scan_column_widths(FILE* infile, uintptr_t** col_widths_ptr, uintptr_t* col_ct_ptr) {
  int32_t retval = 0;
  uintptr_t* col_widths = NULL;
  uintptr_t col_ct = 0;
  uintptr_t max_col_ct = INITIAL_COLS; // not a hard limit
  uintptr_t malloc_size = INITIAL_COLS * sizeof(intptr_t);

  // actually a one-based index, to simplify distinguishing between
  // beginning-of-line-and-not-in-column from beginning-of-line-and-in-column
  // first element of col_widths[] is essentially unused as a result
  uintptr_t cur_col_idx = 0;

  uintptr_t cur_col_width = 0;
  unsigned char* readptr;
  unsigned char* line_end;
  unsigned char* readbuf_end;
  unsigned char* token_end;
  uintptr_t cur_read;
  col_widths = (uintptr_t*)malloc(malloc_size);
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
      }
      if (++cur_col_idx == max_col_ct) {
        if (malloc_double(&malloc_size, (unsigned char**)(&col_widths))) {
          goto scan_column_widths_ret_READ_FAIL;
	}
        fill_ulong_zero(&(col_widths[max_col_ct]), max_col_ct);
	max_col_ct *= 2;
      }
      token_end = get_token_end_ck(readptr, line_end);
      cur_col_width += (uintptr_t)(token_end - readptr);
      if (readptr == line_end) {
	break;
      }
      if (cur_col_width > col_widths[cur_col_idx]) {
	col_widths[cur_col_idx] = cur_col_width;
      }
      cur_col_width = 0;
      readptr = token_end;
    }
    if (line_end < readbuf_end) {
      handle_last_column(col_widths, cur_col_idx, cur_col_width, &col_ct);
      readptr = &(line_end[1]);
      cur_col_idx = 0;
      cur_col_width = 0;
      continue;
    }
    if (feof(infile)) {
      handle_last_column(col_widths, cur_col_idx, cur_col_width, &col_ct);
      break;
    }
    // in middle of line
    cur_read = fread(g_readbuf, 1, BUFSIZE, infile);
    if (ferror(infile)) {
      goto scan_column_widths_ret_READ_FAIL;
    }
    readptr = g_readbuf;
    readbuf_end = &(g_readbuf[cur_read]);
  }
  while (0) {
  scan_column_widths_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  scan_column_widths_ret_READ_FAIL:
    fputs("Error: File read failure.\n", stderr);
    retval = RET_READ_FAIL;
    break;
  scan_column_widths_ret_INVALID_FORMAT:
    fputs("Error: Empty input file.\n", stderr);
    retval = RET_INVALID_FORMAT;
    break;
  }
  *col_widths_ptr = col_widths;
  return retval;
}

#define FLAG_SPACES_BEFORE_FIRST 1
#define FLAG_SPACES_AFTER_LAST 2
#define FLAG_FINAL_EOLN 4
#define FLAG_STRIP_BLANK 8

int32_t pretty_write(FILE* infile, char* outname, uint32_t flags, uint32_t column_sep, uintptr_t* col_widths, uintptr_t col_ct) {
  FILE* outfile = NULL;
  int32_t retval = 0;
  if (!outname) {
    outfile = stdout;
  } else {
    if (fopen_checked(&outfile, outname, "w")) {
      goto pretty_write_ret_OPEN_FAIL;
    }
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
  uint32_t flags = 0;
  uint32_t column_sep = 2;
  uint32_t param_idx = 1;
  int32_t retval = 0;
  uintptr_t* col_widths = NULL;
  uintptr_t col_ct = 0;
  uint32_t infile_param_idx;
  char* param_ptr;
  int32_t ii;
  if ((argc == 1) || ((argc == 2) && ((!strcmp(argv[1], "--help")) || (!strcmp(argv[1], "-help"))))) {
    goto main_ret_HELP;
  } else if (argc > 8) {
    fputs("Error: Too many parameters.\n\n", stderr);
    goto main_ret_INVALID_CMDLINE_2;
  }
  infile_param_idx = argc - 2;
  if ((!infile_param_idx) || (argv[infile_param_idx][0] == '-')) {
    infile_param_idx++;
  } else {
    outname = argv[infile_param_idx + 1];
  }
  for (param_idx = 1; param_idx < infile_param_idx; param_idx++) {
    if (argv[param_idx][0] != '-') {
      if (outname) {
        fputs("Error: All parameters before the last two must be flags.\n\n", stderr);
      } else {
	fputs("Error: Input filename currently must come after all flags.\n\n", stderr);
      }
      goto main_ret_INVALID_CMDLINE_2;
    }
    param_ptr = &(argv[param_idx][1]);
    if (*param_ptr == '-') {
      // allow both single- and double-dash
      param_ptr++;
    }
    if (!strcmp(param_ptr, "spacing")) {
      ii = atoi(argv[++param_idx]);
      if (ii < 1) {
	if (param_idx == infile_param_idx) {
	  fputs("Error: Missing --spacing parameter.\n\n", stderr);
	} else {
	  printf("Error: Invalid --spacing parameter '%s'.\n\n", argv[param_idx]);
	}
	goto main_ret_INVALID_CMDLINE;
      } else if (param_idx == infile_param_idx) {
	fputs("Error: Missing input filename.\n\n", stderr);
	goto main_ret_INVALID_CMDLINE_2;
      }
      column_sep = ii;
    } else if (!strcmp(param_ptr, "first")) {
      flags |= FLAG_SPACES_BEFORE_FIRST;
    } else if (!strcmp(param_ptr, "last")) {
      flags |= FLAG_SPACES_AFTER_LAST;
    } else if (!strcmp(param_ptr, "final-eoln")) {
      flags |= FLAG_FINAL_EOLN;
    } else if (!strcmp(param_ptr, "strip-blank")) {
      flags |= FLAG_STRIP_BLANK;
    } else if (!strcmp(param_ptr, "help")) {
      goto main_ret_HELP;
    } else {
      printf("Error: Invalid flag '%s'.\n\n", argv[param_idx]);
      goto main_ret_INVALID_CMDLINE_2;
    }
  }
  if (argv[infile_param_idx][0] == '-') {
    fputs("Error: Input filename cannot begin with '-'.\n\n", stderr);
    goto main_ret_INVALID_CMDLINE_2;
  }
  if (fopen_checked(&infile, argv[infile_param_idx], "rb")) {
    goto main_ret_OPEN_FAIL;
  }
  retval = scan_column_widths(infile, &col_widths, &col_ct);
  if (retval) {
    goto main_ret_1;
  }
  rewind(infile);
  retval = pretty_write(infile, argv[param_idx], flags, column_sep, col_widths, col_ct);
  if (retval) {
    goto main_ret_1;
  }
  fclose_null(&infile);
  while (0) {
  main_ret_HELP:
    fputs(
"prettify v1.0 (26 Sep 2013)   Christopher Chang (chrchang@alumni.caltech.edu)\n\n"
"Takes a tab-and/or-space-delimited text file, and generates a pretty-printed\n"
"space-delimited text file from it.  Multibyte character encodings are not\n"
"currently supported.\n\n"
, stdout);
    disp_usage(stdout);
    fputs(
"\nTo perform the simplest reverse conversion (spaces to one tab), you can use\n"
"  cat [input filename] | tr -s ' ' '\\t' > [output filename]\n"
"(and for one-to-one conversion between spaces and tabs, omit the \"-s\").\n"
, stdout);
    retval = RET_HELP;
    break;
  main_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
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
  return retval;
}
