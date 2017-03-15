#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Simple natural sort utility.  NOT industrial-strength: lines are limited to
// 256 chars, and the total number of lines is limited to 1048576.

static inline int is_digit(char cc) {
  return (cc <= '9') && (cc >= '0');
}

static inline int is_not_digit(char cc) {
  return (cc > '9') || (cc < '0');
}

static inline int is_not_nzdigit(char cc) {
  return (cc > '9') || (cc <= '0');
}

// PLINK's individual ID natural sort uses the following logic:
// - All alphabetic characters act as if they are capitalized, except for
// tiebreaking purposes (where ASCII is used).
// - Numbers in family and individual IDs are compared by magnitude, with the
// exception of...
// - Numbers with leading zero(es).  If you're putting extraneous zeroes in
// front of IDs, we assume they're there to force particular IDs to be sorted
// earlier, rather than just appearing at random.  So, unlike many natural sort
// implementations, we sort 00200 < 021 < 20: all numbers with n leading zeroes
// are sorted before all numbers with (n-1) leading zeroes; magnitude only
// applies if the leading zero counts match.  This handles e.g. subbasement
// room numbering properly.
//
// This won't always do what you want if your IDs have variable-length decimals
// in them (e.g. it yields 0.99 < 0.101); if you don't want to fall back on
// --merge-ascii-sort, enforce a fixed number of digits after the decimal
// point.
int strcmp_natural_scan_forward(const char* s1, const char* s2) {
  // assumes s1 and s2 currently point to the middle of a mismatching number,
  // where s1 < s2.
  char c1;
  char c2;
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
int strcmp_natural_tiebroken(const char* s1, const char* s2) {
  // assumes ties should be broken in favor of s2.
  char c1 = *(++s1);
  char c2 = *(++s2);
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

static inline int strcmp_natural_uncasted(const char* s1, const char* s2) {
  char c1 = *s1;
  char c2 = *s2;
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

int strcmp_natural(const void* s1, const void* s2) {
  return strcmp_natural_uncasted((char*)s1, (char*)s2);
}

char main_buf[268435456];

int main() {
  char tbuf[256];
  char* cptr = main_buf;
  unsigned int line_ct = 0;
  unsigned int uii;
  while (fgets(tbuf, 256, stdin)) {
    strcpy(cptr, tbuf);
    cptr = &(cptr[256]);
    line_ct++;
    if (line_ct == 1048576) {
      break;
    }
  }
  qsort(main_buf, line_ct, 256, strcmp_natural);
  for (uii = 0; uii < line_ct; uii++) {
    fputs(&(main_buf[uii * 256]), stdout);
  }
  return 0;
}
