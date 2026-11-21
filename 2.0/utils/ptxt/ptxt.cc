// Prototype plink2-library text tools.  Input may be uncompressed, gzip,
// BGZF (decompressed in parallel), or zstd.  Lines must end in LF.
//
// The user-facing documentation of every subcommand is the help text below
// (kSummaryHelp, kModeHelp; printed by 'ptxt' and 'ptxt SUB -h'), mirrored in
// README.md.  Keep the two in sync.

// Correctly-rounded decimal parsing.  plink2's ScanadvDouble may be off by an
// ulp or two, which flips comparisons against thresholds like "<= 0.99" when
// a field is exactly on the boundary (e.g. "0.9900").
#define FFC_IMPL
// ffc.h is C; it uses C-style casts throughout.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#include "../../include/ffc.h"
#pragma GCC diagnostic pop
#include "../../include/plink2_htable.h"
#include "../../include/plink2_text.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <string>
#include <algorithm>
#include <chrono>
#include <thread>
#include <vector>
using namespace plink2;

static const uint32_t kOutBufSize = 1 << 22;
static char* g_outbuf;
static char* g_outw;
static inline void Flush() { fwrite(g_outbuf, 1, g_outw - g_outbuf, stdout); g_outw = g_outbuf; }
static inline void Emit(const char* s, uintptr_t len) {
  if (g_outw + len > g_outbuf + kOutBufSize) {
    Flush();
    if (len > kOutBufSize) { fwrite(s, 1, len, stdout); return; }
  }
  memcpy(g_outw, s, len); g_outw += len;
}
static inline void EmitCh(char c) {
  if (g_outw == g_outbuf + kOutBufSize) Flush();
  *g_outw++ = c;
}

static void Die(const char* msg) { fprintf(stderr, "ptxt: %s\n", msg); exit(1); }

// ---- help ----
//
// 'ptxt' alone prints the subcommand summary; 'ptxt SUB -h' (or 'ptxt help
// SUB') prints everything about one subcommand; 'ptxt -v' prints kVersion.
// Keep lines <= 80 columns.

static const char kVersion[] = "ptxt v0.1.0 (21 Nov 2026)";

static const char kSummaryHelp[] = R"HELP(
  ptxt <subcommand> [options...] <arguments...>
  ptxt <subcommand> -h
  ptxt -v

Fast replacements for common text-processing one-liners.  Input may be
uncompressed, gzip, BGZF (decompressed on several threads) or zstd; use
/dev/stdin as FILE to read from a pipe.  Output goes to stdout.  Columns are
1-based numbers or header names.

Subcommands:
  -- Columns
     select    columns in any order, repeats allowed    awk '{print $3,$1}'
     cut       columns in ascending order (tab only)    cut -f
     paste     two files side by side                   paste
     prefix    add a leading string to a column         awk '{$1="chr"$1}1'
     unprefix  remove a leading string from a column    sed 's/^chr//'
     tsv       whitespace-delimited to tab-delimited    awk '{$1=$1}1' OFS=\t

  -- Rows
     filter    lines satisfying an expression           awk '$P < 5e-8'
     keep      lines whose column is in an ID list      grep -F -w -f
     remove    lines whose column isn't in an ID list   grep -v -F -w -f
     uniq      first line per distinct key              awk '!seen[$3]++'
     sort      stable multi-key sort                    sort -s -k1,1V -k2,2n
     join      hash join on key column(s)               join

  -- Summaries
     count     number of lines per distinct value       cut | sort | uniq -c
     sum       sum of a numeric column                  awk '{s+=$5} END...'

"ptxt <subcommand> -h" describes one subcommand in full; "ptxt -v" (or -V)
prints just the version.
)HELP";

// Options shared by most subcommands; printed after a subcommand's own help.
static const char kColumnsHelp[] = R"HELP(
Columns are 1-based numbers (e.g. 3) or header names (e.g. P); using a name
implies -H.
)HELP";

static const char kInputHelp[] = R"HELP(
Input may be uncompressed, gzip, BGZF or zstd; use /dev/stdin to read from a
pipe.  Lines must end in LF.
)HELP";

struct ModeHelp {
  const char* name;
  const char* body;        // usage line(s) and description
  const char* examples;
  const char* header_doc;  // what -H does; nullptr if [opts] isn't accepted
  uint32_t odelim;         // -O accepted
  uint32_t ranges;         // column ranges accepted
};

static const ModeHelp kModeHelp[] = {
  {"cut", R"HELP(
  ptxt cut LIST FILE

Same as 'cut -f LIST' on tab-delimited input.  LIST is a comma-separated list
of 1-based column numbers and ranges, in ascending order (e.g. 1-8 or 2,11).
Lines with no tab are printed whole, as cut does.  Contiguous ranges are
copied as one block, so this is fastest when you keep runs of columns.

Use 'ptxt select' to reorder or repeat columns, pick them by header name, or
read whitespace-delimited input.

)HELP", R"HELP(
Example:
  ptxt cut 1-5,8 calls.vcf.gz
)HELP", nullptr, 0, 0},

  {"paste", R"HELP(
  ptxt paste FILE1 FILE2

Same as 'paste FILE1 FILE2': prints each FILE1 line, a tab, then the
corresponding FILE2 line.  It is an error for the files to have different
numbers of lines.

)HELP", R"HELP(
Example:
  ptxt paste ids.txt scores.txt
)HELP", nullptr, 0, 0},

  {"select", R"HELP(
  ptxt select [options] COLS FILE

Prints the COLS columns of each line, like awk '{print $3, $1, ...}'.  COLS
is a comma-separated list; columns can be in any order and can be repeated.
Missing fields are printed as empty strings.

)HELP", R"HELP(
Examples:
  ptxt select SNP,P,BETA ss.tsv
  ptxt select -w -O , 2,1,5-7 plink.assoc
)HELP", "It is printed like any other line.", 1, 1},

  {"filter", R"HELP(
  ptxt filter [options] EXPR FILE

Prints the lines that satisfy EXPR, an awk-like expression.  All arguments
before FILE are joined with spaces, so EXPR doesn't have to be one shell
argument (but quote it anyway, to protect $, <, >, ( and &).

Syntax:
  Fields      $N, $NAME, or $"NAME" (names imply -H)
  Literals    numbers (any decimal or scientific notation), "strings"
  Operators   in increasing precedence:
                ||  (also |, or)
                &&  (also &, and)
                <  <=  >  >=  ==  (also =)  !=
                +  -
                *  /  %
                unary !  (also not),  unary -  +
  Functions   abs  log  log10  exp  sqrt
  Grouping    ( )

How values compare (as in awk, except for missing values):
  - A field is a number if the whole field parses as one.
  - Field vs "string" compares as strings.
  - Field vs field compares numerically if both are numbers, as strings
    otherwise.
  - Field vs number compares numerically if the field is a number.
    Otherwise the comparison is false, and != is true.  (This is plink2's
    --extract-if-info convention; awk would compare as strings.)
  - Arithmetic on a non-numeric field gives a missing value: every
    comparison against it is false (!= is true), and it is false as a
    condition.
  - A bare value used as a condition is true if it is a nonzero number or a
    nonempty string.
  - A field past the end of the line is an empty string.
  - Numbers are parsed with correct rounding, so boundary comparisons such as
    0.9900 <= 0.99 agree with awk.

)HELP", R"HELP(
Examples:
  ptxt filter '$P < 5e-8' ss.tsv
  ptxt filter '$P < 5e-8 && ($1 == "22" || $CHR == "X")' ss.tsv.gz
  ptxt filter 'abs($BETA / $SE) > 5 && $A1_FREQ >= 0.01' ss.tsv
  ptxt filter -w '-log10($9) > 8' plink.assoc
)HELP", "It is copied through.", 0, 0},

  {"keep", R"HELP(
  ptxt keep [options] COL IDFILE FILE
  ptxt remove [options] COL IDFILE FILE

'keep' prints the lines whose COL value is listed in IDFILE; 'remove' prints
the lines whose COL value isn't.  Each ID is the first whitespace-delimited
token on an IDFILE line; blank IDFILE lines are ignored.  Matching is exact
string equality, via a hash table, so IDFILE can hold millions of IDs.

This replaces awk 'NR==FNR {a[$1]; next} $3 in a' and 'grep -F -w -f', but
only looks at one column (grep -w matches anywhere on the line).

)HELP", R"HELP(
Examples:
  ptxt keep SNP hits.txt ss.tsv.gz
  ptxt remove -w 2 bad_ids.txt plink.assoc
)HELP", "It is copied through.", 0, 0},

  {"count", R"HELP(
  ptxt count [options] COL FILE

Prints "VALUE<TAB>COUNT" for each distinct value in column COL, in the order
the values first appear.  This replaces 'cut -f COL | sort | uniq -c' and
awk '{c[$1]++} END {for (k in c) print k, c[k]}'; the input doesn't need to
be sorted.  Lines that are missing COL aren't counted.

Pipe the output to 'ptxt sort 2:nr /dev/stdin' to order it by count.

)HELP", R"HELP(
Example:
  ptxt count CHR ss.tsv
)HELP", "It is skipped.", 0, 0},

  {"prefix", R"HELP(
  ptxt prefix [options] COL STR FILE
  ptxt unprefix [options] COL STR FILE

'prefix' adds STR to the start of column COL.  'unprefix' removes STR from
the start of COL wherever it is present.  The rest of each line, including
delimiters, is copied unchanged.

)HELP", R"HELP(
Examples:
  ptxt prefix 1 chr calls_b37.tsv
  ptxt unprefix CHR chr ss.tsv.gz
)HELP", "It is copied through.", 0, 0},

  {"tsv", R"HELP(
  ptxt tsv FILE

Converts whitespace-delimited lines, such as PLINK 1 reports, to
tab-delimited.  Leading and trailing spaces and tabs are dropped and each run
of spaces and tabs becomes one tab, as with awk '{$1=$1}1' OFS='\t'.  This is
what 'tr -s " " "\t"' is usually meant to do, except that tr leaves a leading
tab on padded lines.

)HELP", R"HELP(
Example:
  ptxt tsv plink.assoc > plink.assoc.tsv
)HELP", nullptr, 0, 0},

  {"sum", R"HELP(
  ptxt sum [options] COL FILE

Prints the sum of the numeric values in column COL.  Non-numeric and missing
values are skipped.

)HELP", R"HELP(
Example:
  ptxt sum -w 5 plink.frq.counts
)HELP", "It is skipped.", 0, 0},

  {"sort", R"HELP(
  ptxt sort [options] KEY[:TYPE][,KEY[:TYPE]...] FILE

Stable sort on one or more keys: lines with equal keys stay in input order,
as with 'sort -s'.  Each KEY is a single column; later keys break ties in
earlier ones (at most 32 keys).  TYPE is one of:
  s   string, in byte order (default)
  n   number, in any decimal or scientific notation ('g' is an alias).
      Non-numeric and missing values always sort last, in either direction.
  c   chromosome: 1, 2, ..., 22 (and higher numbers), X, Y, XY, MT (or M),
      with or without a chr prefix (any case), then other contigs in byte
      order.
Append r to a TYPE (or use just r, for a reversed string key) to sort that
key in descending order.

The whole input is loaded into memory.  The sort uses all cores, up to 16;
set PTXT_TIMING=1 to print per-phase timings to stderr.

)HELP", R"HELP(
Examples:
  ptxt sort CHR:c,POS:n ss.tsv > ss.sorted.tsv
  ptxt sort P:n ss.tsv.gz | head
  ptxt sort -H 7:nr,3 ss.tsv
)HELP", "It stays first.", 0, 0},

  {"join", R"HELP(
  ptxt join [options] [-a] [-e STR] KEYS1 FILE1 KEYS2 FILE2

Hash join.  For each FILE1 line, prints that line once for every FILE2 row
whose KEYS2 columns exactly match the line's KEYS1 columns, followed by the
FILE2 row's non-key columns.
  - Output keeps FILE1's order; several matching FILE2 rows come out in
    FILE2 order.  Neither file needs to be sorted.
  - KEYS1 and KEYS2 are comma-separated lists of the same length, so
    composite keys (e.g. CHR,POS) work.  Ranges aren't allowed.
  - FILE2 is held in memory and FILE1 is streamed, so put the smaller file
    second.
  - Short FILE2 rows are padded with empty fields.  FILE2 rows missing a key
    column, and blank FILE2 lines, are skipped.

Join options:
  -a        Also print FILE1 lines that have no match, padded with -e values.
  -e STR    Fill value for -a (default: NA).

With -H (or any column name), both files have header lines.  The output
header is FILE1's header followed by FILE2's non-key column names.

)HELP", R"HELP(
Examples:
  ptxt join SNP sumstats.tsv ID annot.tsv
  ptxt join -a -e . CHR,POS a.tsv chrom,pos b.tsv.gz
)HELP", "Both files must have one.", 1, 0},

  {"uniq", R"HELP(
  ptxt uniq [options] [-c] [-r|-u] [COLS] FILE

Prints the first line for each distinct key, in the order keys first appear.
The key is the whole line, or the COLS columns if COLS is given (a
comma-separated list; a missing field counts as empty, as in awk).  The input
doesn't need to be sorted: with no flags this is awk '!seen[$3]++', and on
sorted input the output matches uniq's.

Uniq options:
  -c        Put each key's count and a tab in front of its line.
  -r        Only print keys that occur more than once ('uniq -d'; -d is
            ptxt's delimiter option).
  -u        Only print keys that occur exactly once.
-c, -r and -u print their output at the end, keeping each key's first line in
memory.  Without them, lines are printed as soon as they're read.

Sorted input is especially fast: while keys arrive in ascending order no
hash table is needed.

)HELP", R"HELP(
Examples:
  ptxt uniq SNP ss.tsv                  # drop repeated SNP IDs
  ptxt uniq -c -r CHR,POS ss.tsv.gz     # list duplicated positions
  sort ids.txt | ptxt uniq /dev/stdin
)HELP", "It is copied through, with\n"
        "            \"count<TAB>\" in front under -c.", 0, 0},
};

static const ModeHelp* FindModeHelp(const char* mode) {
  if (!strcmp(mode, "remove")) mode = "keep";
  else if (!strcmp(mode, "unprefix")) mode = "prefix";
  for (const ModeHelp& h : kModeHelp) {
    if (!strcmp(h.name, mode)) return &h;
  }
  return nullptr;
}

static void PrintBanner() {
  // Right-align the tagline at column 80, as plink2 does.
  static const char kTagline[] = "text tools built on the plink2 library";
  printf("%s%*s\n", kVersion, static_cast<int>(80 - strlen(kVersion)), kTagline);
}

[[noreturn]] static void PrintSummaryHelp(int exit_code) {
  PrintBanner();
  fputs(kSummaryHelp, stdout);
  exit(exit_code);
}

[[noreturn]] static void PrintModeHelp(const ModeHelp* h, int exit_code) {
  fputs(h->body + 1, stdout);  // skip the raw string's leading newline
  if (h->header_doc) {
    fputs("Common options:\n"
          "  -w        Whitespace-delimited input (awk's default field splitting).\n"
          "  -d C      Input delimiter character (default: tab).\n", stdout);
    if (h->odelim) fputs("  -O STR    Output delimiter (default: tab).\n", stdout);
    printf("  -H        The first line is a header.  %s\n", h->header_doc);
    fputs("  --        End of options; needed only if the next argument is\n"
          "            literally one of the options above.\n", stdout);
    fputs(kColumnsHelp, stdout);
    if (h->ranges) fputs("Ranges such as 3-5 (or 5-3, descending) are allowed.\n", stdout);
  }
  // The body ends with a blank line; the sections after it don't.
  if (h->header_doc) putchar('\n');
  fputs(kInputHelp + 1, stdout);
  fputs(h->examples, stdout);
  exit(exit_code);
}

static const char* g_mode = "";

static bool IsHelpFlag(const char* a) { return !strcmp(a, "-h") || !strcmp(a, "--help"); }

// For argument-count errors: the one-line synopsis plus a pointer to -h.
[[noreturn]] static void DieUsage(const char* synopsis) {
  fprintf(stderr, "ptxt: usage: %s\n(\"ptxt %s -h\" for details.)\n", synopsis, g_mode);
  exit(1);
}

static void OpenTs(const char* fname, TextStream* txsp) {
  PreinitTextStream(txsp);
  if (TextStreamOpen(fname, txsp)) Die("failed to open input file");
}

static void FinishTs(TextStream* txsp) {
  if (TextStreamErrcode(txsp)) Die("read failure");
}

static inline uint32_t IsTokChar(char c) { return static_cast<unsigned char>(c) > 32; }

// Returns 1 iff [s, e) is entirely a finite decimal number (optional sign,
// digits/decimal point, optional exponent), parsed with correct rounding.
// inf/nan/hex are rejected, as awk does.
static inline uint32_t ParseNumExact(const char* s, const char* e, double* dp) {
  const char* p = s;
  const char* q;
  if (p != e && *p == '+') {
    q = ++p;
  } else {
    q = (p != e && *p == '-') ? p + 1 : p;
  }
  if (q == e) return 0;
  if (!((*q >= '0' && *q <= '9') || *q == '.')) return 0;
  if (*q == '0' && q + 1 != e && (q[1] == 'x' || q[1] == 'X')) return 0;
  const ffc_result r = ffc_from_chars_double(p, e, dp);
  return (r.outcome == FFC_OUTCOME_OK) && (r.ptr == e);
}

// ---- field location ----
// delim == '\0' means whitespace mode.  lend points at the line's '\n'.
// Returns the start of 0-based field col, or nullptr if the line has fewer
// fields.
static inline char* FieldStart(char* line, char* lend, uint32_t col, char delim) {
  if (delim) {
    if (!col) return line;
    char* p = AdvToNthDelimChecked(line, lend, col, delim);
    return p ? p + 1 : nullptr;
  }
  char* p = FirstNonTspace(line);
  if (!IsTokChar(*p)) return nullptr;
  return col ? NextTokenMult(p, col) : p;
}

static inline char* FieldEnd(char* fstart, char* lend, char delim) {
  if (delim) {
    char* e = static_cast<char*>(memchr(fstart, delim, lend - fstart));
    return e ? e : lend;
  }
  return CurTokenEnd(fstart);
}

// ---- options and column specs ----
struct Opts {
  char delim = '\t';
  const char* odelim = "\t";
  uint32_t header = 0;
  uint32_t left_outer = 0;  // join -a
  const char* fill = "NA";  // join -e
};

// Parses leading options; returns index of first positional arg.
// join_opts enables -a and -e (join only).
static int ParseOpts(int argc, char** argv, int start, Opts* o, uint32_t join_opts = 0) {
  int i = start;
  // Only exact option names are consumed, so a filter expression may start
  // with '-' (e.g. '-log10($P) > 8').  '--' ends option parsing.
  for (; i < argc; ++i) {
    const char* a = argv[i];
    if (!strcmp(a, "--")) { ++i; break; }
    if (IsHelpFlag(a)) PrintModeHelp(FindModeHelp(g_mode), 0);
    const uint32_t is_join_opt = join_opts && (!strcmp(a, "-a") || !strcmp(a, "-e"));
    if (strcmp(a, "-w") && strcmp(a, "-H") && strcmp(a, "-d") && strcmp(a, "-O") && !is_join_opt) break;
    if (!strcmp(a, "-w")) o->delim = '\0';
    else if (!strcmp(a, "-a")) o->left_outer = 1;
    else if (!strcmp(a, "-H")) o->header = 1;
    else if (i + 1 == argc) Die("option is missing its argument");
    else if (!strcmp(a, "-d")) {
      o->delim = argv[++i][0];
      if (o->delim == '\n' || !o->delim) Die("invalid -d delimiter");
    } else if (!strcmp(a, "-e")) {
      o->fill = argv[++i];
    } else {
      o->odelim = argv[++i];
    }
  }
  return i;
}

struct ColTok { uint32_t a, b; const char* name; };

static uint32_t IsRangeTok(const char* s, uint32_t len) {
  uint32_t i = 0, digits = 0;
  while (i < len && s[i] >= '0' && s[i] <= '9') { ++i; ++digits; }
  if (!digits) return 0;
  if (i == len) return 1;
  if (s[i] != '-') return 0;
  ++i; digits = 0;
  while (i < len && s[i] >= '0' && s[i] <= '9') { ++i; ++digits; }
  return digits && i == len;
}

// Returns 1 iff any column is given by name.
static uint32_t ParseCols(char* spec, std::vector<ColTok>* toks) {
  uint32_t has_names = 0;
  for (char* p = spec; ; ) {
    char* e = strchr(p, ',');
    uint32_t len = e ? e - p : strlen(p);
    if (!len) Die("empty column in list");
    ColTok t = {0, 0, nullptr};
    if (IsRangeTok(p, len)) {
      char* q;
      t.a = strtoul(p, &q, 10); t.b = t.a;
      if (*q == '-') t.b = strtoul(q + 1, &q, 10);
      if (!t.a || !t.b) Die("columns are 1-based");
    } else {
      if (e) *e = '\0';
      t.name = p; has_names = 1;
    }
    toks->push_back(t);
    if (!e) break;
    p = e + 1;
  }
  return has_names;
}

// Expands tokens into 0-based column indices, resolving names against the
// header line (may be nullptr if no names are present).
static std::vector<uint32_t> ResolveCols(const std::vector<ColTok>& toks, char* hline, char* hlend, char delim) {
  std::vector<uint32_t> cols;
  for (const ColTok& t : toks) {
    if (!t.name) {
      if (t.a <= t.b) for (uint32_t c = t.a; c <= t.b; ++c) cols.push_back(c - 1);
      else for (uint32_t c = t.a; c >= t.b; --c) cols.push_back(c - 1);
      continue;
    }
    const uint32_t nlen = strlen(t.name);
    uint32_t idx = 0;
    char* f = FieldStart(hline, hlend, 0, delim);
    for (; f; ++idx) {
      char* fe = FieldEnd(f, hlend, delim);
      if (static_cast<uint32_t>(fe - f) == nlen && !memcmp(f, t.name, nlen)) break;
      if (delim) f = (fe == hlend) ? nullptr : fe + 1;
      else { f = FirstNonTspace(fe); if (!IsTokChar(*f)) f = nullptr; }
    }
    if (!f) { fprintf(stderr, "ptxt: column '%s' not found in header\n", t.name); exit(1); }
    cols.push_back(idx);
  }
  return cols;
}

// Opens FILE, reads the header line if needed, and resolves COLS.
struct Input {
  TextStream txs;
  char* hline = nullptr;
  char* hlend = nullptr;
  std::vector<uint32_t> cols;
};

static void OpenInput(const char* fname, char* colspec, Opts* o, Input* in) {
  std::vector<ColTok> toks;
  if (colspec && ParseCols(colspec, &toks)) o->header = 1;
  OpenTs(fname, &in->txs);
  if (o->header) {
    if (TextNextLine(&in->txs, &in->hline)) Die("missing header line");
    in->hlend = TextLineEnd(&in->txs) - 1;
  }
  if (colspec) in->cols = ResolveCols(toks, in->hline, in->hlend, o->delim);
}

static uint32_t SingleCol(const Input& in) {
  if (in.cols.size() != 1) Die("exactly one column expected");
  return in.cols[0];
}

static inline void EmitHeader(const Input& in) {
  if (in.hline) Emit(in.hline, in.hlend + 1 - in.hline);
}

// ---- string set (open addressing, stores lengths so no overread) ----
struct StrSet {
  // Each slot holds (hash << 32 | index), so mismatched probes are rejected
  // without touching the entry; entries hold (pointer, length, hash) together
  // so a hash match costs one more cache line, not three.
  struct Entry {
    const char* p;
    uint32_t len;
    uint32_t h;
  };
  static constexpr uint64_t kEmpty = ~0ULL;
  std::vector<uint64_t> table;
  std::vector<Entry> entries;
  std::vector<char*> chunks;
  char* arena_cur = nullptr;
  uintptr_t arena_left = 0;

  StrSet() { table.assign(1 << 16, kEmpty); }
  ~StrSet() { for (char* c : chunks) free(c); }

  inline uint32_t Size() const { return entries.size(); }
  inline const char* Ptr(uint32_t idx) const { return entries[idx].p; }
  inline uint32_t Len(uint32_t idx) const { return entries[idx].len; }

  static inline uint32_t Slot(uint32_t h, uint32_t size) {
    return (static_cast<uint64_t>(h) * size) >> 32;
  }
  // Returns index, or UINT32_MAX.
  inline uint32_t Find(const char* s, uint32_t len, uint32_t h) const {
    const uint32_t size = table.size();
    for (uint32_t slot = Slot(h, size); ; ) {
      const uint64_t v = table[slot];
      if (v == kEmpty) return UINT32_MAX;
      if (static_cast<uint32_t>(v >> 32) == h) {
        const uint32_t idx = static_cast<uint32_t>(v);
        const Entry& e = entries[idx];
        if (e.len == len && !memcmp(e.p, s, len)) return idx;
      }
      if (++slot == size) slot = 0;
    }
  }
  // Returns index of s, inserting it if absent.  *added set to 1 iff new.
  uint32_t Insert(const char* s, uint32_t len, uint32_t h, uint32_t* added) {
    uint32_t idx = Find(s, len, h);
    *added = (idx == UINT32_MAX);
    if (!*added) return idx;
    if (len > arena_left) {
      const uintptr_t csize = len > (1 << 20) ? len : (1 << 20);
      char* c = static_cast<char*>(malloc(csize));
      if (!c) Die("out of memory");
      chunks.push_back(c); arena_cur = c; arena_left = csize;
    }
    memcpy(arena_cur, s, len);
    idx = entries.size();
    if (idx == UINT32_MAX) Die("too many distinct strings");
    entries.push_back({arena_cur, len, h});
    arena_cur += len; arena_left -= len;
    if (entries.size() * 4 > table.size()) {
      table.assign(table.size() * 2, kEmpty);
      for (uint32_t i = 0; i < entries.size(); ++i) Place(i);
    } else {
      Place(idx);
    }
    return idx;
  }
  // Appends s without indexing it; the caller guarantees it's new.  Used
  // while input is sorted.  BuildIndex() must be called before the next
  // Find/Insert.
  uint32_t AppendUnindexed(const char* s, uint32_t len) {
    if (len > arena_left) {
      const uintptr_t csize = len > (1 << 20) ? len : (1 << 20);
      char* c = static_cast<char*>(malloc(csize));
      if (!c) Die("out of memory");
      chunks.push_back(c); arena_cur = c; arena_left = csize;
    }
    memcpy(arena_cur, s, len);
    const uint32_t idx = entries.size();
    if (idx == UINT32_MAX) Die("too many distinct strings");
    entries.push_back({arena_cur, len, 0});
    arena_cur += len; arena_left -= len;
    return idx;
  }
  void BuildIndex() {
    uintptr_t size = table.size();
    while (entries.size() * 4 > size) size *= 2;
    table.assign(size, kEmpty);
    for (uint32_t i = 0; i < entries.size(); ++i) {
      entries[i].h = Hash32(entries[i].p, entries[i].len);
      Place(i);
    }
  }
  void Place(uint32_t idx) {
    const uint32_t size = table.size();
    const uint32_t h = entries[idx].h;
    for (uint32_t slot = Slot(h, size); ; ) {
      if (table[slot] == kEmpty) { table[slot] = (static_cast<uint64_t>(h) << 32) | idx; return; }
      if (++slot == size) slot = 0;
    }
  }
};

static inline uint32_t HashStr(const char* s, uint32_t len) { return Hash32(s, len); }

// ---- modes ----
static void ModeCut(int argc, char** argv) {
  if (argc != 4) DieUsage("ptxt cut LIST FILE");
  std::vector<uint32_t> cols;
  for (char* p = argv[2]; *p; ) {
    uint32_t a = strtoul(p, &p, 10), b = a;
    if (*p == '-') b = strtoul(p + 1, &p, 10);
    if (!a || b < a) Die("invalid LIST");
    for (uint32_t c = a; c <= b; ++c) cols.push_back(c - 1);
    if (*p == ',') ++p;
  }
  // Merge into contiguous ranges [rstart[r], rend[r]] (0-based, inclusive).
  std::vector<uint32_t> rstart, rend;
  for (uint32_t c : cols) {
    if (!rstart.empty() && c <= rend.back() + 1) {
      if (c > rend.back()) rend.back() = c;
    } else {
      rstart.push_back(c); rend.push_back(c);
    }
  }
  const uint32_t nr = rstart.size();
  TextStream txs; OpenTs(argv[3], &txs);
  char* line;
  while (!TextNextLine(&txs, &line)) {
    char* lend = TextLineEnd(&txs) - 1;
    // 'it' points at the start of field 'cur'.
    char* it = line; uint32_t cur = 0;
    for (uint32_t r = 0; r < nr; ++r) {
      uint32_t skip = rstart[r] - cur;
      if (skip) {
        it = AdvToNthDelimChecked(it, lend, skip, '\t');
        if (!it) {
          if (!r && !memchr(line, '\t', lend - line)) Emit(line, lend - line);
          break;
        }
        ++it;
      }
      // Find the delimiter terminating field rend[r]; copy the whole span.
      char* fend = AdvToNthDelimChecked(it, lend, rend[r] - rstart[r] + 1, '\t');
      if (!fend) fend = lend;
      if (r) EmitCh('\t');
      Emit(it, fend - it);
      if (fend == lend) break;
      it = fend + 1; cur = rend[r] + 1;
    }
    EmitCh('\n');
  }
  FinishTs(&txs);
}

static void ModePaste(int argc, char** argv) {
  if (argc != 4) DieUsage("ptxt paste FILE1 FILE2");
  TextStream t1, t2; OpenTs(argv[2], &t1); OpenTs(argv[3], &t2);
  char *l1, *l2;
  while (!TextNextLine(&t1, &l1)) {
    if (TextNextLine(&t2, &l2)) Die("FILE2 has fewer lines than FILE1");
    Emit(l1, TextLineEnd(&t1) - 1 - l1);
    EmitCh('\t');
    Emit(l2, TextLineEnd(&t2) - l2);
  }
  FinishTs(&t1);
  if (!TextNextLine(&t2, &l2)) Die("FILE2 has more lines than FILE1");
}

static void ModeSelect(int argc, char** argv) {
  Opts o;
  int i = ParseOpts(argc, argv, 2, &o);
  if (argc - i != 2) DieUsage("ptxt select [options] COLS FILE");
  Input in; OpenInput(argv[i + 1], argv[i], &o, &in);
  const std::vector<uint32_t>& req = in.cols;
  // Sorted unique needed columns, and map from request slot to needed slot.
  std::vector<uint32_t> need(req);
  std::sort(need.begin(), need.end());
  need.erase(std::unique(need.begin(), need.end()), need.end());
  std::vector<uint32_t> slot_of(req.size());
  for (uint32_t k = 0; k < req.size(); ++k) {
    slot_of[k] = std::lower_bound(need.begin(), need.end(), req[k]) - need.begin();
  }
  const uint32_t nn = need.size();
  std::vector<char*> fs(nn), fe(nn);
  const char delim = o.delim;
  const char* odelim = o.odelim;
  const uint32_t odlen = strlen(odelim);
  char* line = in.hline;
  char* lend = in.hlend;
  if (!line) {
    if (TextNextLine(&in.txs, &line)) { FinishTs(&in.txs); return; }
    lend = TextLineEnd(&in.txs) - 1;
  }
  do {
    // Locate needed fields in one left-to-right pass.
    uint32_t found = 0;
    char* it = FieldStart(line, lend, 0, delim);
    uint32_t cur = 0;
    if (it) {
      for (; found < nn; ++found) {
        const uint32_t skip = need[found] - cur;
        if (skip) {
          if (delim) {
            it = AdvToNthDelimChecked(it, lend, skip, delim);
            if (!it) break;
            ++it;
          } else {
            it = NextTokenMult(it, skip);
            if (!it) break;
          }
          cur = need[found];
        }
        fs[found] = it;
        fe[found] = FieldEnd(it, lend, delim);
      }
    }
    for (uint32_t k = 0; k < req.size(); ++k) {
      if (k) Emit(odelim, odlen);
      const uint32_t s = slot_of[k];
      if (s < found) Emit(fs[s], fe[s] - fs[s]);
    }
    EmitCh('\n');
    if (TextNextLine(&in.txs, &line)) break;
    lend = TextLineEnd(&in.txs) - 1;
  } while (1);
  FinishTs(&in.txs);
}

// ---- filter: awk-like boolean expressions ----
//
// Grammar (awk precedence, lowest first):
//   expr    := and ( ('||' | '|' | 'or') and )*
//   and     := cmp ( ('&&' | '&' | 'and') cmp )*
//   cmp     := add [ ('<' | '<=' | '>' | '>=' | '==' | '=' | '!=') add ]
//   add     := mul ( ('+' | '-') mul )*
//   mul     := unary ( ('*' | '/' | '%') unary )*
//   unary   := ('!' | 'not' | '-' | '+') unary | primary
//   primary := NUMBER | "STRING" | $N | $NAME | $"NAME" | FUNC '(' expr ')'
//            | '(' expr ')'
//   FUNC    := abs | log | log10 | exp | sqrt
//
// Value semantics follow awk, except that missing values never satisfy a
// comparison (plink2 --extract-if-info convention):
// - A field is "strnum": numeric if the whole field parses as a number.
// - field vs "string": string comparison.  field vs field: numeric if both
//   are numeric, string otherwise.
// - field vs number (or arithmetic result): numeric if the field is numeric;
//   otherwise the comparison is false ('!=' is true).
// - Arithmetic on a non-numeric field yields a missing value, which makes any
//   comparison false ('!=' true) and is false as a condition.
// - A bare value used as a condition is true iff it is nonzero (numbers) or
//   nonempty (strings), as in awk.
// - A field past the end of the line is the empty string.

enum ExprNodeType {
  kEnNum, kEnStr, kEnField, kEnNeg, kEnNot, kEnAdd, kEnSub, kEnMul, kEnDiv,
  kEnMod, kEnLt, kEnLe, kEnGt, kEnGe, kEnEq, kEnNe, kEnAnd, kEnOr, kEnFunc,
  // Specialized comparisons, rewritten after parsing: field OP number and
  // field OP "string" (op stored in ExprNode::func).
  kEnCmpFN, kEnCmpFS
};

enum ExprFunc { kEfAbs, kEfLog, kEfLog10, kEfExp, kEfSqrt };

struct ExprNode {
  ExprNodeType t;
  uint32_t a = 0, b = 0;  // child node indices
  double num = 0;
  std::string str;        // string literal, or column name for kEnField
  uint32_t col = 0;       // 0-based column (kEnField), resolved after parse
  uint32_t slot = 0;      // index into the sorted needed-column list
  uint32_t func = 0;
};

class ExprParser {
 public:
  ExprParser(const std::string& src, std::vector<ExprNode>* nodes) : src_(src), p_(0), nodes_(nodes) {}

  uint32_t Parse() {
    const uint32_t root = ParseOr();
    SkipWs();
    if (p_ != src_.size()) Fail("unexpected input");
    return root;
  }

 private:
  const std::string& src_;
  size_t p_;
  std::vector<ExprNode>* nodes_;

  [[noreturn]] void Fail(const char* msg) {
    fprintf(stderr, "ptxt: invalid filter expression (%s):\n  %s\n  %*s^\n", msg, src_.c_str(), static_cast<int>(p_), "");
    exit(1);
  }
  void SkipWs() {
    while (p_ < src_.size() && static_cast<unsigned char>(src_[p_]) <= 32) ++p_;
  }
  char Peek(size_t off = 0) const {
    return (p_ + off < src_.size()) ? src_[p_ + off] : '\0';
  }
  // Matches a case-insensitive keyword not followed by an identifier char.
  bool Keyword(const char* kw) {
    const size_t len = strlen(kw);
    if (p_ + len > src_.size() || strncasecmp(&src_[p_], kw, len)) return false;
    const char c = Peek(len);
    if (isalnum(static_cast<unsigned char>(c)) || c == '_') return false;
    p_ += len;
    return true;
  }
  uint32_t Add(ExprNode n) {
    nodes_->push_back(std::move(n));
    return nodes_->size() - 1;
  }
  uint32_t Bin(ExprNodeType t, uint32_t a, uint32_t b) {
    ExprNode n; n.t = t; n.a = a; n.b = b;
    return Add(std::move(n));
  }

  uint32_t ParseOr() {
    uint32_t lhs = ParseAnd();
    while (1) {
      SkipWs();
      if (Peek() == '|') {
        p_ += (Peek(1) == '|') ? 2 : 1;
      } else if (!Keyword("or")) {
        return lhs;
      }
      lhs = Bin(kEnOr, lhs, ParseAnd());
    }
  }
  uint32_t ParseAnd() {
    uint32_t lhs = ParseCmp();
    while (1) {
      SkipWs();
      if (Peek() == '&') {
        p_ += (Peek(1) == '&') ? 2 : 1;
      } else if (!Keyword("and")) {
        return lhs;
      }
      lhs = Bin(kEnAnd, lhs, ParseCmp());
    }
  }
  uint32_t ParseCmp() {
    const uint32_t lhs = ParseAdd();
    SkipWs();
    const char c0 = Peek(), c1 = Peek(1);
    ExprNodeType t;
    if (c0 == '<') { t = (c1 == '=') ? kEnLe : kEnLt; p_ += 1 + (c1 == '='); }
    else if (c0 == '>') { t = (c1 == '=') ? kEnGe : kEnGt; p_ += 1 + (c1 == '='); }
    else if (c0 == '=') { t = kEnEq; p_ += 1 + (c1 == '='); }
    else if (c0 == '!' && c1 == '=') { t = kEnNe; p_ += 2; }
    else if (c0 == '~' || (c0 == '!' && c1 == '~')) Fail("regex matching is not supported");
    else return lhs;
    const uint32_t rhs = ParseAdd();
    SkipWs();
    const char d0 = Peek(), d1 = Peek(1);
    if (d0 == '<' || d0 == '>' || (d0 == '=' && d1 != '=') || (d0 == '=' && d1 == '=') || (d0 == '!' && d1 == '=')) {
      Fail("comparisons can't be chained; use && (e.g. 0.01 <= $6 && $6 <= 0.99)");
    }
    return Bin(t, lhs, rhs);
  }
  uint32_t ParseAdd() {
    uint32_t lhs = ParseMul();
    while (1) {
      SkipWs();
      const char c = Peek();
      if (c != '+' && c != '-') return lhs;
      ++p_;
      lhs = Bin((c == '+') ? kEnAdd : kEnSub, lhs, ParseMul());
    }
  }
  uint32_t ParseMul() {
    uint32_t lhs = ParseUnary();
    while (1) {
      SkipWs();
      const char c = Peek();
      if (c != '*' && c != '/' && c != '%') return lhs;
      ++p_;
      lhs = Bin((c == '*') ? kEnMul : ((c == '/') ? kEnDiv : kEnMod), lhs, ParseUnary());
    }
  }
  uint32_t ParseUnary() {
    SkipWs();
    const char c = Peek();
    if ((c == '!' && Peek(1) != '=') || Keyword("not")) {
      if (c == '!') ++p_;
      return Bin(kEnNot, ParseUnary(), 0);
    }
    if (c == '-' || c == '+') {
      ++p_;
      const uint32_t child = ParseUnary();
      if (c == '+') return child;
      if ((*nodes_)[child].t == kEnNum) {
        (*nodes_)[child].num = -(*nodes_)[child].num;
        return child;
      }
      return Bin(kEnNeg, child, 0);
    }
    return ParsePrimary();
  }
  std::string ParseQuoted() {
    // assumes Peek() == '"'
    ++p_;
    std::string s;
    while (1) {
      if (p_ >= src_.size()) Fail("unterminated string");
      char c = src_[p_++];
      if (c == '"') return s;
      if (c == '\\') {
        if (p_ >= src_.size()) Fail("unterminated string");
        c = src_[p_++];
        if (c == 't') c = '\t';
        else if (c == 'n') c = '\n';
      }
      s.push_back(c);
    }
  }
  uint32_t ParsePrimary() {
    SkipWs();
    const char c = Peek();
    if (c == '(') {
      ++p_;
      const uint32_t inner = ParseOr();
      SkipWs();
      if (Peek() != ')') Fail("expected ')'");
      ++p_;
      return inner;
    }
    if (c == '"') {
      ExprNode n; n.t = kEnStr; n.str = ParseQuoted();
      return Add(std::move(n));
    }
    if (c == '$') {
      ++p_;
      ExprNode n; n.t = kEnField;
      if (Peek() == '"') {
        n.str = ParseQuoted();
      } else if (isdigit(static_cast<unsigned char>(Peek()))) {
        uint64_t v = 0;
        while (isdigit(static_cast<unsigned char>(Peek()))) {
          v = v * 10 + (src_[p_++] - '0');
          if (v > 0x7fffffff) Fail("column number too large");
        }
        if (!v) Fail("$0 is not supported; columns are 1-based");
        n.col = v - 1;
      } else {
        // Name runs until an operator, paren, quote, comma, or whitespace.
        const size_t start = p_;
        while (p_ < src_.size() && !strchr(" \t\n\r!\"&()<=>|~+-*/%,", src_[p_])) ++p_;
        if (p_ == start) Fail("expected column number or name after '$'");
        n.str = src_.substr(start, p_ - start);
      }
      return Add(std::move(n));
    }
    if (isdigit(static_cast<unsigned char>(c)) || (c == '.' && isdigit(static_cast<unsigned char>(Peek(1))))) {
      double v;
      const char* start = &src_[p_];
      const ffc_result r = ffc_from_chars_double(start, src_.data() + src_.size(), &v);
      if (r.outcome != FFC_OUTCOME_OK) Fail("invalid number");
      const char* e = r.ptr;
      p_ += e - start;
      if (isalpha(static_cast<unsigned char>(Peek())) || Peek() == '_') Fail("invalid number");
      ExprNode n; n.t = kEnNum; n.num = v;
      return Add(std::move(n));
    }
    static const char* kFuncNames[] = {"abs", "log", "log10", "exp", "sqrt"};
    // Try longer names first so 'log10' isn't parsed as 'log'.
    static const uint32_t kFuncOrder[] = {kEfLog10, kEfAbs, kEfLog, kEfExp, kEfSqrt};
    for (uint32_t fi : kFuncOrder) {
      const size_t save = p_;
      if (Keyword(kFuncNames[fi])) {
        SkipWs();
        if (Peek() != '(') { p_ = save; break; }
        ++p_;
        const uint32_t arg = ParseOr();
        SkipWs();
        if (Peek() != ')') Fail("expected ')'");
        ++p_;
        ExprNode n; n.t = kEnFunc; n.a = arg; n.func = fi;
        return Add(std::move(n));
      }
    }
    if (isalpha(static_cast<unsigned char>(c)) || c == '_') {
      Fail("unknown identifier (refer to columns as $3 or $NAME, strings as \"...\")");
    }
    Fail("expected a value");
  }
};

// Evaluation state for one line.
struct ExprLine {
  // Configuration.
  const std::vector<uint32_t>* need;  // sorted unique 0-based columns
  char delim;
  // Per-line state.
  char* line;
  char* lend;
  uint32_t located;     // need[0..located) have been located
  uint32_t exhausted;   // line ended before need[located]
  char* it;             // start of field it_col, if located > 0
  uint32_t it_col;
  std::vector<char*> fs, fe;
  std::vector<uint8_t> numstate;  // 0 = unparsed, 1 = numeric, 2 = not
  std::vector<double> nums;

  void Init(const std::vector<uint32_t>* need_, char delim_) {
    need = need_; delim = delim_;
    fs.resize(need->size()); fe.resize(need->size());
    numstate.resize(need->size()); nums.resize(need->size());
  }
  void Reset(char* line_, char* lend_) {
    line = line_; lend = lend_; located = 0; exhausted = 0;
  }
  // Locates need[slot] (and everything before it), scanning only as far as
  // needed.  Missing fields become empty strings.
  inline void Locate(uint32_t slot) {
    if (slot < located) return;
    LocateSlow(slot);
  }
  // Returns 1 and sets *dp iff field slot is numeric (parsed at most once).
  inline uint32_t FieldNum(uint32_t slot, double* dp) {
    Locate(slot);
    uint8_t st = numstate[slot];
    if (!st) {
      st = ParseNumExact(fs[slot], fe[slot], &nums[slot]) ? 1 : 2;
      numstate[slot] = st;
    }
    *dp = nums[slot];
    return st == 1;
  }
  void LocateSlow(uint32_t slot) {
    static char empty[1] = {'\0'};
    while (located <= slot) {
      const uint32_t col = (*need)[located];
      char* f = nullptr;
      if (!exhausted) {
        if (!located) {
          f = FieldStart(line, lend, col, delim);
        } else {
          const uint32_t skip = col - it_col;
          if (delim) {
            f = AdvToNthDelimChecked(it, lend, skip, delim);
            if (f) ++f;
          } else {
            f = NextTokenMult(it, skip);
          }
        }
      }
      if (f) {
        it = f; it_col = col;
        fs[located] = f;
        fe[located] = FieldEnd(f, lend, delim);
      } else {
        exhausted = 1;
        fs[located] = empty; fe[located] = empty;
      }
      numstate[located] = 0;
      ++located;
    }
  }
};

enum ValKind : uint8_t { kVkNum, kVkStr, kVkMissing };

struct Val {
  ValKind k;
  uint8_t is_field;  // string form available from the line
  double num;
  const char* s;
  uint32_t len;
};

static inline uint32_t Truthy(const Val& v) {
  if (v.k == kVkNum) return v.num != 0;
  if (v.k == kVkStr) return v.len != 0;
  return 0;
}

static inline Val NumVal(double d) {
  Val v; v.k = kVkNum; v.is_field = 0; v.num = d; v.s = nullptr; v.len = 0;
  return v;
}

static inline Val MissingVal() {
  Val v; v.k = kVkMissing; v.is_field = 0; v.num = 0; v.s = nullptr; v.len = 0;
  return v;
}

// Converts to a number for arithmetic; returns 0 if missing/non-numeric.
static inline uint32_t ToNum(const Val& v, double* dp) {
  if (v.k == kVkNum) { *dp = v.num; return 1; }
  if (v.k == kVkStr && !v.is_field && v.len) {
    // string literal: numeric iff it fully parses
    return ParseNumExact(v.s, v.s + v.len, dp);
  }
  return 0;
}

static inline int StrCmp3(const char* a, uint32_t alen, const char* b, uint32_t blen) {
  const int r = memcmp(a, b, alen < blen ? alen : blen);
  if (r) return r;
  return (alen > blen) - (alen < blen);
}

static inline uint32_t CmpResult(uint32_t op, int c) {
  switch (op) {
    case kEnLt: return c < 0;
    case kEnLe: return c <= 0;
    case kEnGt: return c > 0;
    case kEnGe: return c >= 0;
    case kEnEq: return c == 0;
    default: return c != 0;
  }
}

static Val EvalExpr(const std::vector<ExprNode>& nodes, uint32_t idx, ExprLine* lc) {
  const ExprNode& n = nodes[idx];
  switch (n.t) {
    case kEnNum:
      return NumVal(n.num);
    case kEnStr: {
      Val v; v.k = kVkStr; v.is_field = 0; v.num = 0;
      v.s = n.str.data(); v.len = n.str.size();
      return v;
    }
    case kEnField: {
      const uint32_t slot = n.slot;
      double d;
      const uint32_t is_num = lc->FieldNum(slot, &d);
      Val v; v.is_field = 1; v.s = lc->fs[slot];
      v.len = lc->fe[slot] - lc->fs[slot];
      if (is_num) { v.k = kVkNum; v.num = d; }
      else { v.k = kVkStr; v.num = 0; }
      return v;
    }
    case kEnCmpFN: {
      double d;
      if (!lc->FieldNum(n.slot, &d) || d != d) return NumVal(n.func == kEnNe);
      return NumVal(CmpResult(n.func, (d > n.num) - (d < n.num)));
    }
    case kEnCmpFS: {
      const uint32_t slot = n.slot;
      lc->Locate(slot);
      return NumVal(CmpResult(n.func, StrCmp3(lc->fs[slot], lc->fe[slot] - lc->fs[slot], n.str.data(), n.str.size())));
    }
    case kEnAnd:
      return NumVal(Truthy(EvalExpr(nodes, n.a, lc)) && Truthy(EvalExpr(nodes, n.b, lc)));
    case kEnOr:
      return NumVal(Truthy(EvalExpr(nodes, n.a, lc)) || Truthy(EvalExpr(nodes, n.b, lc)));
    case kEnNot:
      return NumVal(!Truthy(EvalExpr(nodes, n.a, lc)));
    case kEnNeg: {
      double d;
      if (!ToNum(EvalExpr(nodes, n.a, lc), &d)) return MissingVal();
      return NumVal(-d);
    }
    case kEnFunc: {
      double d;
      if (!ToNum(EvalExpr(nodes, n.a, lc), &d)) return MissingVal();
      switch (n.func) {
        case kEfAbs: d = fabs(d); break;
        case kEfLog: d = log(d); break;
        case kEfLog10: d = log10(d); break;
        case kEfExp: d = exp(d); break;
        default: d = sqrt(d); break;
      }
      // log(0), sqrt(-1), etc. propagate as NaN/inf, which compares like
      // missing (NaN) or as a real infinity.
      return NumVal(d);
    }
    case kEnAdd: case kEnSub: case kEnMul: case kEnDiv: case kEnMod: {
      double x, y;
      if (!ToNum(EvalExpr(nodes, n.a, lc), &x)) return MissingVal();
      if (!ToNum(EvalExpr(nodes, n.b, lc), &y)) return MissingVal();
      switch (n.t) {
        case kEnAdd: return NumVal(x + y);
        case kEnSub: return NumVal(x - y);
        case kEnMul: return NumVal(x * y);
        case kEnDiv: return y != 0 ? NumVal(x / y) : MissingVal();
        default: return y != 0 ? NumVal(fmod(x, y)) : MissingVal();
      }
    }
    default: break;
  }
  // Comparison.
  const Val x = EvalExpr(nodes, n.a, lc);
  const Val y = EvalExpr(nodes, n.b, lc);
  int c;
  if (x.k == kVkMissing || y.k == kVkMissing) return NumVal(n.t == kEnNe);
  if (x.k == kVkNum && y.k == kVkNum) {
    if (x.num != x.num || y.num != y.num) return NumVal(n.t == kEnNe);  // NaN
    c = (x.num > y.num) - (x.num < y.num);
  } else if ((x.k == kVkNum && !x.is_field && y.is_field) || (y.k == kVkNum && !y.is_field && x.is_field)) {
    // Number vs non-numeric field: missing.
    return NumVal(n.t == kEnNe);
  } else {
    // String comparison; a non-field number (vs a string literal) is
    // formatted like awk's CONVFMT.
    char bx[32], by[32];
    const char* sx = x.s; uint32_t lx = x.len;
    const char* sy = y.s; uint32_t ly = y.len;
    if (x.k == kVkNum && !x.is_field) { lx = snprintf(bx, sizeof bx, "%.6g", x.num); sx = bx; }
    if (y.k == kVkNum && !y.is_field) { ly = snprintf(by, sizeof by, "%.6g", y.num); sy = by; }
    c = StrCmp3(sx, lx, sy, ly);
  }
  switch (n.t) {
    case kEnLt: return NumVal(c < 0);
    case kEnLe: return NumVal(c <= 0);
    case kEnGt: return NumVal(c > 0);
    case kEnGe: return NumVal(c >= 0);
    case kEnEq: return NumVal(c == 0);
    default: return NumVal(c != 0);
  }
}

static void ModeFilter(int argc, char** argv) {
  Opts o;
  int i = ParseOpts(argc, argv, 2, &o);
  if (argc - i < 2) DieUsage("ptxt filter [options] EXPR FILE");
  // All args before FILE are joined with spaces, so quoting is optional
  // where the shell permits.
  std::string src;
  for (int k = i; k < argc - 1; ++k) {
    if (k > i) src.push_back(' ');
    src += argv[k];
  }
  std::vector<ExprNode> nodes;
  ExprParser parser(src, &nodes);
  const uint32_t root = parser.Parse();
  uint32_t any_names = 0;
  for (const ExprNode& n : nodes) any_names |= (n.t == kEnField && !n.str.empty());
  if (any_names) o.header = 1;
  Input in; OpenInput(argv[argc - 1], nullptr, &o, &in);
  const char delim = o.delim;
  // Resolve names; collect needed columns.
  std::vector<uint32_t> need;
  for (ExprNode& n : nodes) {
    if (n.t != kEnField) continue;
    if (!n.str.empty()) {
      ColTok t = {0, 0, n.str.c_str()};
      std::vector<ColTok> toks(1, t);
      n.col = ResolveCols(toks, in.hline, in.hlend, delim)[0];
    }
    need.push_back(n.col);
  }
  std::sort(need.begin(), need.end());
  need.erase(std::unique(need.begin(), need.end()), need.end());
  for (ExprNode& n : nodes) {
    if (n.t == kEnField) n.slot = std::lower_bound(need.begin(), need.end(), n.col) - need.begin();
  }
  // Rewrite field-vs-literal comparisons into specialized nodes.  A literal on
  // the left is moved to the right by mirroring the operator.
  for (ExprNode& n : nodes) {
    if (n.t < kEnLt || n.t > kEnNe) continue;
    uint32_t fa = n.a, la = n.b;
    ExprNodeType op = n.t;
    if (nodes[fa].t != kEnField) {
      std::swap(fa, la);
      if (op == kEnLt) op = kEnGt;
      else if (op == kEnGt) op = kEnLt;
      else if (op == kEnLe) op = kEnGe;
      else if (op == kEnGe) op = kEnLe;
    }
    if (nodes[fa].t != kEnField) continue;
    const ExprNodeType lt = nodes[la].t;
    if (lt != kEnNum && lt != kEnStr) continue;
    n.slot = nodes[fa].slot;
    n.func = op;
    if (lt == kEnNum) {
      n.t = kEnCmpFN;
      n.num = nodes[la].num;
    } else {
      n.t = kEnCmpFS;
      n.str = nodes[la].str;
    }
  }
  ExprLine lc;
  lc.Init(&need, delim);
  EmitHeader(in);
  char* line;
  while (!TextNextLine(&in.txs, &line)) {
    char* lnext = TextLineEnd(&in.txs);
    lc.Reset(line, lnext - 1);
    if (Truthy(EvalExpr(nodes, root, &lc))) Emit(line, lnext - line);
  }
  FinishTs(&in.txs);
}

static void ModeKeep(int argc, char** argv, uint32_t invert) {
  Opts o;
  int i = ParseOpts(argc, argv, 2, &o);
  if (argc - i != 3) DieUsage("ptxt keep|remove [options] COL IDFILE FILE");
  StrSet ids;
  {
    TextStream its; OpenTs(argv[i + 1], &its);
    char* line;
    while (!TextNextLine(&its, &line)) {
      char* t = FirstNonTspace(line);
      if (!IsTokChar(*t)) continue;
      char* te = CurTokenEnd(t);
      uint32_t added;
      ids.Insert(t, te - t, HashStr(t, te - t), &added);
    }
    FinishTs(&its);
    CleanupTextStream(&its, nullptr);
  }
  Input in; OpenInput(argv[i + 2], argv[i], &o, &in);
  const uint32_t col = SingleCol(in);
  const char delim = o.delim;
  EmitHeader(in);
  char* line;
  while (!TextNextLine(&in.txs, &line)) {
    char* lnext = TextLineEnd(&in.txs);
    char* lend = lnext - 1;
    char* f = FieldStart(line, lend, col, delim);
    uint32_t hit = 0;
    if (f) {
      char* fe = FieldEnd(f, lend, delim);
      const uint32_t len = fe - f;
      hit = (ids.Find(f, len, HashStr(f, len)) != UINT32_MAX);
    }
    if (hit != invert) Emit(line, lnext - line);
  }
  FinishTs(&in.txs);
}

static void ModeCount(int argc, char** argv) {
  Opts o;
  int i = ParseOpts(argc, argv, 2, &o);
  if (argc - i != 2) DieUsage("ptxt count [options] COL FILE");
  Input in; OpenInput(argv[i + 1], argv[i], &o, &in);
  const uint32_t col = SingleCol(in);
  const char delim = o.delim;
  StrSet vals;
  std::vector<uint64_t> counts;
  // Most columns worth counting have few distinct values, often repeated in
  // long runs (e.g. CHR); cache the previous value to skip hashing.
  const char* prev = nullptr; uint32_t prev_len = 0; uint32_t prev_idx = 0;
  char* line;
  while (!TextNextLine(&in.txs, &line)) {
    char* lend = TextLineEnd(&in.txs) - 1;
    char* f = FieldStart(line, lend, col, delim);
    if (!f) continue;
    char* fe = FieldEnd(f, lend, delim);
    const uint32_t len = fe - f;
    if (prev && len == prev_len && !memcmp(f, prev, len)) {
      ++counts[prev_idx];
      continue;
    }
    uint32_t added;
    const uint32_t idx = vals.Insert(f, len, HashStr(f, len), &added);
    if (added) counts.push_back(0);
    ++counts[idx];
    prev = vals.Ptr(idx); prev_len = len; prev_idx = idx;
  }
  FinishTs(&in.txs);
  for (uint32_t k = 0; k < vals.Size(); ++k) {
    char numbuf[24];
    const int nlen = snprintf(numbuf, sizeof numbuf, "\t%llu\n", static_cast<unsigned long long>(counts[k]));
    Emit(vals.Ptr(k), vals.Len(k));
    Emit(numbuf, nlen);
  }
}

static void ModePrefix(int argc, char** argv, uint32_t remove) {
  Opts o;
  int i = ParseOpts(argc, argv, 2, &o);
  if (argc - i != 3) DieUsage("ptxt prefix|unprefix [options] COL STR FILE");
  const char* str = argv[i + 1];
  const uint32_t slen = strlen(str);
  Input in; OpenInput(argv[i + 2], argv[i], &o, &in);
  const uint32_t col = SingleCol(in);
  const char delim = o.delim;
  EmitHeader(in);
  char* line;
  while (!TextNextLine(&in.txs, &line)) {
    char* lnext = TextLineEnd(&in.txs);
    char* f = FieldStart(line, lnext - 1, col, delim);
    if (f) {
      if (!remove) {
        Emit(line, f - line);
        Emit(str, slen);
        Emit(f, lnext - f);
        continue;
      }
      // FieldEnd isn't needed: if STR matches, it can't extend past the
      // field unless STR itself contains a delimiter.
      if (static_cast<uintptr_t>(lnext - 1 - f) >= slen && !memcmp(f, str, slen)) {
        Emit(line, f - line);
        Emit(f + slen, lnext - f - slen);
        continue;
      }
    }
    Emit(line, lnext - line);
  }
  FinishTs(&in.txs);
}

static void ModeTsv(int argc, char** argv) {
  if (argc != 3) DieUsage("ptxt tsv FILE");
  TextStream txs; OpenTs(argv[2], &txs);
  char* line;
  while (!TextNextLine(&txs, &line)) {
    char* p = FirstNonTspace(line);
    if (IsTokChar(*p)) {
      while (1) {
        char* e = CurTokenEnd(p);
        Emit(p, e - p);
        p = FirstNonTspace(e);
        if (!IsTokChar(*p)) break;
        EmitCh('\t');
      }
    }
    EmitCh('\n');
  }
  FinishTs(&txs);
}

static void ModeSum(int argc, char** argv) {
  Opts o;
  int i = ParseOpts(argc, argv, 2, &o);
  if (argc - i != 2) DieUsage("ptxt sum [options] COL FILE");
  Input in; OpenInput(argv[i + 1], argv[i], &o, &in);
  const uint32_t col = SingleCol(in);
  const char delim = o.delim;
  double tot = 0;
  char* line;
  while (!TextNextLine(&in.txs, &line)) {
    char* lend = TextLineEnd(&in.txs) - 1;
    char* f = FieldStart(line, lend, col, delim);
    if (!f) continue;
    double v;
    if (ParseNumExact(f, FieldEnd(f, lend, delim), &v)) tot += v;
  }
  FinishTs(&in.txs);
  char buf[32];
  Emit(buf, snprintf(buf, sizeof buf, "%.10g\n", tot));
}

// ---- join ----
//
// FILE2 (the lookup table) is loaded into a hash table keyed by its key
// columns; each row's non-key columns are stored as a ready-to-emit payload.
// FILE1 is then streamed: each line is printed once per matching FILE2 row
// (in FILE2 order), followed by that row's payload.

// Simple bump allocator.
struct Arena {
  std::vector<char*> chunks;
  char* cur = nullptr;
  uintptr_t left = 0;
  ~Arena() { for (char* c : chunks) free(c); }
  char* Alloc(uintptr_t len) {
    if (len > left) {
      const uintptr_t csize = len > (16 << 20) ? len : (16 << 20);
      cur = static_cast<char*>(malloc(csize));
      if (!cur) Die("out of memory");
      chunks.push_back(cur);
      left = csize;
    }
    char* r = cur;
    cur += len; left -= len;
    return r;
  }
};

// Splits a line into fields (all of them).  Absent trailing fields are just
// not present.
static void SplitFields(char* line, char* lend, char delim, std::vector<char*>* fs, std::vector<char*>* fe) {
  fs->clear(); fe->clear();
  char* f = FieldStart(line, lend, 0, delim);
  while (f) {
    char* e = FieldEnd(f, lend, delim);
    fs->push_back(f); fe->push_back(e);
    if (delim) {
      f = (e == lend) ? nullptr : e + 1;
    } else {
      f = FirstNonTspace(e);
      if (!IsTokChar(*f)) f = nullptr;
    }
  }
}

static void ModeJoin(int argc, char** argv) {
  Opts o;
  int i = ParseOpts(argc, argv, 2, &o, 1);
  if (argc - i != 4) DieUsage("ptxt join [options] [-a] [-e STR] KEYS1 FILE1 KEYS2 FILE2");
  std::vector<ColTok> toks1, toks2;
  if (ParseCols(argv[i], &toks1)) o.header = 1;
  if (ParseCols(argv[i + 2], &toks2)) o.header = 1;
  for (const std::vector<ColTok>* t : {&toks1, &toks2}) {
    for (const ColTok& c : *t) {
      if (!c.name && c.a != c.b) Die("join keys must be single columns, not ranges");
    }
  }
  if (toks1.size() != toks2.size()) Die("KEYS1 and KEYS2 must have the same number of columns");
  const uint32_t nk = toks1.size();
  const char delim = o.delim;
  const char* odelim = o.odelim;
  const uint32_t odlen = strlen(odelim);

  // ---- load FILE2 ----
  TextStream t2; OpenTs(argv[i + 3], &t2);
  char* line;
  std::vector<char*> fs, fe;
  std::string header2_payload;
  std::vector<uint32_t> kcols2;
  std::vector<uint8_t> is_key2;
  uint32_t fill_ct = UINT32_MAX;  // payload column count (from first line)
  auto resolve2 = [&](char* hline, char* hlend) {
    std::vector<uint32_t> cols = ResolveCols(toks2, hline, hlend, delim);
    kcols2 = cols;
    for (uint32_t c : cols) {
      if (c >= is_key2.size()) is_key2.resize(c + 1, 0);
      if (is_key2[c]) Die("duplicate column in KEYS2");
      is_key2[c] = 1;
    }
  };
  // Appends the non-key fields of the current split line, each preceded by
  // odelim, to *dst.
  // Returns the number of fields appended.
  auto payload = [&](std::string* dst) {
    uint32_t ct = 0;
    for (uint32_t c = 0; c < fs.size(); ++c) {
      if (c < is_key2.size() && is_key2[c]) continue;
      dst->append(odelim, odlen);
      dst->append(fs[c], fe[c] - fs[c]);
      ++ct;
    }
    return ct;
  };
  if (o.header) {
    if (TextNextLine(&t2, &line)) Die("FILE2 is missing its header line");
    char* lend = TextLineEnd(&t2) - 1;
    resolve2(line, lend);
    SplitFields(line, lend, delim, &fs, &fe);
    payload(&header2_payload);
    fill_ct = fs.size() - nk;
  } else {
    resolve2(nullptr, nullptr);
  }
  StrSet keys;
  std::vector<uint32_t> head, tail;     // per distinct key: first/last row
  std::vector<uint32_t> next;           // per row: next row with same key
  std::vector<const char*> pay;
  std::vector<uint32_t> paylen;
  Arena arena;
  std::string keybuf, paybuf;
  while (!TextNextLine(&t2, &line)) {
    char* lend = TextLineEnd(&t2) - 1;
    if (lend == line) continue;  // blank line (not a record)
    SplitFields(line, lend, delim, &fs, &fe);
    if (fs.empty()) continue;  // whitespace-only line under -w
    if (fill_ct == UINT32_MAX) fill_ct = (fs.size() >= nk) ? fs.size() - nk : 0;
    // Rows lacking a key column can't match anything.
    uint32_t ok = 1;
    for (uint32_t c : kcols2) ok &= (c < fs.size());
    if (!ok) continue;
    const char* kp;
    uint32_t klen;
    if (nk == 1) {
      kp = fs[kcols2[0]]; klen = fe[kcols2[0]] - fs[kcols2[0]];
    } else {
      keybuf.clear();
      for (uint32_t j = 0; j < nk; ++j) {
        if (j) keybuf.push_back('\0');
        keybuf.append(fs[kcols2[j]], fe[kcols2[j]] - fs[kcols2[j]]);
      }
      kp = keybuf.data(); klen = keybuf.size();
    }
    paybuf.clear();
    // Rows with fewer columns than the first are padded with empty fields,
    // so every output line has the same column count.
    for (uint32_t pc = payload(&paybuf); pc < fill_ct; ++pc) paybuf.append(odelim, odlen);
    char* pdst = arena.Alloc(paybuf.size());
    memcpy(pdst, paybuf.data(), paybuf.size());
    const uint32_t row = pay.size();
    if (row == UINT32_MAX) Die("FILE2 has too many rows");
    pay.push_back(pdst); paylen.push_back(paybuf.size()); next.push_back(UINT32_MAX);
    uint32_t added;
    const uint32_t kidx = keys.Insert(kp, klen, HashStr(kp, klen), &added);
    if (added) {
      head.push_back(row); tail.push_back(row);
    } else {
      next[tail[kidx]] = row;
      tail[kidx] = row;
    }
  }
  FinishTs(&t2);
  CleanupTextStream(&t2, nullptr);
  if (fill_ct == UINT32_MAX) fill_ct = 0;
  std::string fill;
  for (uint32_t k = 0; k < fill_ct; ++k) { fill.append(odelim, odlen); fill += o.fill; }

  // ---- stream FILE1 ----
  Input in; OpenInput(argv[i + 1], nullptr, &o, &in);
  std::vector<uint32_t> kcols1 = ResolveCols(toks1, in.hline, in.hlend, delim);
  if (in.hline) {
    Emit(in.hline, in.hlend - in.hline);
    Emit(header2_payload.data(), header2_payload.size());
    EmitCh('\n');
  }
  // Locate key fields in one ascending pass (they may be listed in any
  // order).
  std::vector<uint32_t> need(kcols1);
  std::sort(need.begin(), need.end());
  need.erase(std::unique(need.begin(), need.end()), need.end());
  std::vector<uint32_t> kslot(nk);
  for (uint32_t j = 0; j < nk; ++j) kslot[j] = std::lower_bound(need.begin(), need.end(), kcols1[j]) - need.begin();
  const uint32_t nn = need.size();
  std::vector<char*> ks(nn), ke(nn);
  while (!TextNextLine(&in.txs, &line)) {
    char* lnext = TextLineEnd(&in.txs);
    char* lend = lnext - 1;
    uint32_t kidx = UINT32_MAX;
    char* it = FieldStart(line, lend, 0, delim);
    uint32_t found = 0;
    if (it) {
      uint32_t cur = 0;
      for (; found < nn; ++found) {
        const uint32_t skip = need[found] - cur;
        if (skip) {
          if (delim) {
            it = AdvToNthDelimChecked(it, lend, skip, delim);
            if (!it) break;
            ++it;
          } else {
            it = NextTokenMult(it, skip);
            if (!it) break;
          }
          cur = need[found];
        }
        ks[found] = it;
        ke[found] = FieldEnd(it, lend, delim);
      }
    }
    if (found == nn) {
      if (nk == 1) {
        const uint32_t klen = ke[0] - ks[0];
        kidx = keys.Find(ks[0], klen, HashStr(ks[0], klen));
      } else {
        keybuf.clear();
        for (uint32_t j = 0; j < nk; ++j) {
          if (j) keybuf.push_back('\0');
          keybuf.append(ks[kslot[j]], ke[kslot[j]] - ks[kslot[j]]);
        }
        kidx = keys.Find(keybuf.data(), keybuf.size(), HashStr(keybuf.data(), keybuf.size()));
      }
    }
    if (kidx == UINT32_MAX) {
      if (o.left_outer) {
        Emit(line, lend - line);
        Emit(fill.data(), fill.size());
        EmitCh('\n');
      }
      continue;
    }
    for (uint32_t row = head[kidx]; row != UINT32_MAX; row = next[row]) {
      Emit(line, lend - line);
      Emit(pay[row], paylen[row]);
      EmitCh('\n');
    }
  }
  FinishTs(&in.txs);
}

// ---- sort ----
//
// Every key is parsed once into a uint64 whose unsigned order is the key
// order (reversal and missing-last already applied).  Where 64 bits can't
// decide equality (strings longer than 8 bytes, nonstandard contig names), a
// per-key flag bit sends ties to a full comparison.  Remaining ties are
// broken by input order, so the sort is stable.

enum SortKeyType : uint8_t { kSkStr, kSkNum, kSkChr };

struct SortKey {
  ColTok tok;
  uint32_t col;
  SortKeyType type;
  uint8_t rev;
};

// Keys are encoded as a sequence of uint64 words (1 per numeric/chromosome
// key, 2 per string key: a 16-byte big-endian prefix).  The first
// kSortInlineWords live in the record itself; any others go in a side array.
static const uint32_t kSortInlineWords = 3;

struct SortRec {
  uint64_t w[kSortInlineWords];
  uint32_t idx;     // input line index
  uint32_t flags;   // bit j: key j needs a full comparison on equality
};

static const uint64_t kSortMissing = ~0ULL;

static inline uint64_t EncodeDouble(double d) {
  if (d == 0) d = 0;  // -0 == 0
  uint64_t u;
  memcpy(&u, &d, sizeof u);
  return (u >> 63) ? ~u : (u | (1ULL << 63));
}

// Big-endian 16-byte prefix, zero-padded: (hi, lo) compares like the
// strings' first 16 bytes.
static inline void StrPrefix16(const char* s, uint32_t len, uint64_t* hi, uint64_t* lo) {
  unsigned char b[16] = {0};
  memcpy(b, s, len < 16 ? len : 16);
  uint64_t h = 0, l = 0;
  for (uint32_t i = 0; i < 8; ++i) h = (h << 8) | b[i];
  for (uint32_t i = 8; i < 16; ++i) l = (l << 8) | b[i];
  *hi = h; *lo = l;
}

// Natural chromosome order: 1, 2, ..., 22, ..., X, Y, XY, MT (optional "chr"
// prefix, case-insensitive); anything else sorts after these, by name.
// Returns UINT32_MAX for "anything else".
static uint32_t ChrRank(const char* s, uint32_t len) {
  if (len > 3 && !strncasecmp(s, "chr", 3)) { s += 3; len -= 3; }
  if (!len) return UINT32_MAX;
  if (s[0] >= '1' && s[0] <= '9' && len <= 6) {
    uint32_t v = 0;
    for (uint32_t i = 0; i < len; ++i) {
      if (s[i] < '0' || s[i] > '9') return UINT32_MAX;
      v = v * 10 + (s[i] - '0');
    }
    return v;
  }
  if (len == 1) {
    const char c = s[0] & 0xdf;
    if (c == 'X') return 1000001;
    if (c == 'Y') return 1000002;
    if (c == 'M') return 1000004;
    return UINT32_MAX;
  }
  if (len == 2) {
    const char c0 = s[0] & 0xdf, c1 = s[1] & 0xdf;
    if (c0 == 'X' && c1 == 'Y') return 1000003;
    if (c0 == 'M' && c1 == 'T') return 1000004;
  }
  return UINT32_MAX;
}

struct SortCtx {
  std::vector<SortKey> keys;
  std::vector<uint32_t> word_start, word_ct;    // per key
  uint32_t total_words = 0;
  uint32_t extra_stride = 0;                     // words beyond the inline ones
  std::vector<uint64_t> extra;                   // n * extra_stride
  std::vector<std::vector<const char*>> fptr;    // string/chr keys
  std::vector<std::vector<uint32_t>> flen;

  inline uint64_t Word(const SortRec& r, uint32_t w) const {
    return (w < kSortInlineWords) ? r.w[w] : extra[static_cast<uintptr_t>(r.idx) * extra_stride + (w - kSortInlineWords)];
  }
};

// Full comparison of key j for two lines whose encoded keys are equal.
static inline int SortFullCmp(const SortCtx& c, uint32_t j, uint32_t a, uint32_t b) {
  const SortKey& k = c.keys[j];
  if (k.type == kSkNum) return 0;
  const char* sa = c.fptr[j][a]; const uint32_t la = c.flen[j][a];
  const char* sb = c.fptr[j][b]; const uint32_t lb = c.flen[j][b];
  const int r = StrCmp3(sa, la, sb, lb);
  return k.rev ? -r : r;
}

static inline bool SortLess(const SortCtx& c, const SortRec& x, const SortRec& y) {
  // Common case: decided by the first inline word.
  if (x.w[0] != y.w[0]) return x.w[0] < y.w[0];
  const uint32_t fl = x.flags | y.flags;
  const uint32_t nk = c.keys.size();
  uint32_t w = 1;
  for (uint32_t j = 0; j < nk; ++j) {
    for (const uint32_t wend = c.word_start[j] + c.word_ct[j]; w < wend; ++w) {
      const uint64_t a = c.Word(x, w), b = c.Word(y, w);
      if (a != b) return a < b;
    }
    if (fl & (1U << j)) {
      const int r = SortFullCmp(c, j, x.idx, y.idx);
      if (r) return r < 0;
    }
  }
  return x.idx < y.idx;
}

static void ParseSortKeys(char* spec, std::vector<SortKey>* keys) {
  for (char* p = spec; ; ) {
    char* e = strchr(p, ',');
    if (e) *e = '\0';
    SortKey k;
    k.type = kSkStr; k.rev = 0; k.col = 0;
    char* colon = strchr(p, ':');
    if (colon) {
      *colon = '\0';
      for (const char* t = colon + 1; *t; ++t) {
        if (*t == 's') k.type = kSkStr;
        else if (*t == 'n' || *t == 'g') k.type = kSkNum;
        else if (*t == 'c') k.type = kSkChr;
        else if (*t == 'r') k.rev = 1;
        else Die("sort key types are s (string), n (numeric), c (chromosome), plus optional r (reverse)");
      }
    }
    std::vector<ColTok> toks;
    ParseCols(p, &toks);
    if (toks.size() != 1 || (!toks[0].name && toks[0].a != toks[0].b)) {
      Die("each sort key must be a single column");
    }
    k.tok = toks[0];
    keys->push_back(k);
    if (!e) break;
    p = e + 1;
  }
  if (keys->size() > 32) Die("at most 32 sort keys are supported");
}

template <typename F>
static void ParallelFor(uint32_t thread_ct, uintptr_t n, F fn) {
  if (thread_ct <= 1 || n < 65536) { fn(0, n); return; }
  std::vector<std::thread> threads;
  for (uint32_t t = 0; t < thread_ct; ++t) {
    const uintptr_t lo = n * t / thread_ct, hi = n * (t + 1) / thread_ct;
    threads.emplace_back([=]() { fn(lo, hi); });
  }
  for (std::thread& th : threads) th.join();
}

static void ModeSort(int argc, char** argv) {
  Opts o;
  int i = ParseOpts(argc, argv, 2, &o);
  if (argc - i != 2) DieUsage("ptxt sort [options] KEY[:TYPE][,KEY[:TYPE]...] FILE");
  SortCtx ctx;
  ParseSortKeys(argv[i], &ctx.keys);
  const uint32_t nk = ctx.keys.size();
  for (const SortKey& k : ctx.keys) if (k.tok.name) o.header = 1;
  Input in; OpenInput(argv[i + 1], nullptr, &o, &in);
  const char delim = o.delim;
  for (SortKey& k : ctx.keys) {
    std::vector<ColTok> toks(1, k.tok);
    k.col = ResolveCols(toks, in.hline, in.hlend, delim)[0];
  }
  // The header points into the stream's buffer, which loading overwrites.
  const std::string header = in.hline ? std::string(in.hline, in.hlend + 1 - in.hline) : std::string();

  // Phase timing on stderr when PTXT_TIMING is set.
  const bool timing = getenv("PTXT_TIMING") != nullptr;
  auto t_last = std::chrono::steady_clock::now();
  auto phase = [&](const char* name) {
    if (!timing) return;
    const auto now = std::chrono::steady_clock::now();
    fprintf(stderr, "%-8s %.3fs\n", name, std::chrono::duration<double>(now - t_last).count());
    t_last = now;
  };
  // Load all lines (including their '\n') into chunked storage.
  std::vector<char*> chunks;
  char* chunk_cur = nullptr;
  uintptr_t chunk_left = 0;
  std::vector<char*> lines;
  std::vector<uint32_t> lens;
  char* line;
  while (!TextNextLine(&in.txs, &line)) {
    const uintptr_t len = TextLineEnd(&in.txs) - line;
    if (len > 0xffffffffU) Die("line too long");
    if (len > chunk_left) {
      const uintptr_t csize = len > (64 << 20) ? len : (64 << 20);
      chunk_cur = static_cast<char*>(malloc(csize));
      if (!chunk_cur) Die("out of memory");
      chunks.push_back(chunk_cur);
      chunk_left = csize;
    }
    memcpy(chunk_cur, line, len);
    lines.push_back(chunk_cur);
    lens.push_back(len);
    chunk_cur += len; chunk_left -= len;
  }
  FinishTs(&in.txs);
  const uintptr_t n = lines.size();
  if (n > 0xffffffffU) Die("too many lines");

  phase("load");
  // Extract keys in parallel.
  uint32_t thread_ct = std::thread::hardware_concurrency();
  if (!thread_ct) thread_ct = 1;
  if (thread_ct > 16) thread_ct = 16;
  std::vector<uint32_t> need;
  for (const SortKey& k : ctx.keys) need.push_back(k.col);
  std::sort(need.begin(), need.end());
  need.erase(std::unique(need.begin(), need.end()), need.end());
  std::vector<uint32_t> key_slot(nk);
  for (uint32_t j = 0; j < nk; ++j) {
    key_slot[j] = std::lower_bound(need.begin(), need.end(), ctx.keys[j].col) - need.begin();
  }
  ctx.fptr.resize(nk); ctx.flen.resize(nk);
  for (uint32_t j = 0; j < nk; ++j) {
    ctx.word_start.push_back(ctx.total_words);
    const uint32_t wc = (ctx.keys[j].type == kSkStr) ? 2 : 1;
    ctx.word_ct.push_back(wc);
    ctx.total_words += wc;
    if (ctx.keys[j].type != kSkNum) { ctx.fptr[j].resize(n); ctx.flen[j].resize(n); }
  }
  if (ctx.total_words > kSortInlineWords) {
    ctx.extra_stride = ctx.total_words - kSortInlineWords;
    ctx.extra.resize(static_cast<uintptr_t>(n) * ctx.extra_stride);
  }
  std::vector<SortRec> recs(n);
  ParallelFor(thread_ct, n, [&](uintptr_t lo, uintptr_t hi) {
    const uint32_t nn = need.size();
    std::vector<char*> fs(nn), fe(nn);
    static char empty[1] = {'\0'};
    for (uintptr_t li = lo; li < hi; ++li) {
      char* ln = lines[li];
      char* lend = ln + lens[li] - 1;  // the '\n'
      // Locate needed fields left to right; absent ones are empty.
      char* it = FieldStart(ln, lend, 0, delim);
      uint32_t cur = 0, found = 0;
      if (it) {
        for (; found < nn; ++found) {
          const uint32_t skip = need[found] - cur;
          if (skip) {
            if (delim) {
              it = AdvToNthDelimChecked(it, lend, skip, delim);
              if (!it) break;
              ++it;
            } else {
              it = NextTokenMult(it, skip);
              if (!it) break;
            }
            cur = need[found];
          }
          fs[found] = it;
          fe[found] = FieldEnd(it, lend, delim);
        }
      }
      for (uint32_t s = found; s < nn; ++s) { fs[s] = empty; fe[s] = empty; }
      SortRec& r = recs[li];
      r.idx = li; r.flags = 0;
      for (uint32_t w = 0; w < kSortInlineWords; ++w) r.w[w] = 0;
      uint64_t* extra = ctx.extra_stride ? &ctx.extra[li * ctx.extra_stride] : nullptr;
      auto put = [&](uint32_t w, uint64_t v) {
        if (w < kSortInlineWords) r.w[w] = v;
        else extra[w - kSortInlineWords] = v;
      };
      for (uint32_t j = 0; j < nk; ++j) {
        const SortKey& k = ctx.keys[j];
        const char* f = fs[key_slot[j]];
        const uint32_t flen = fe[key_slot[j]] - fs[key_slot[j]];
        const uint32_t w0 = ctx.word_start[j];
        if (k.type == kSkNum) {
          double d;
          uint64_t v = kSortMissing;
          if (ParseNumExact(f, f + flen, &d)) {
            v = EncodeDouble(d);
            if (k.rev) v = ~v;
          }
          put(w0, v);
          continue;
        }
        ctx.fptr[j][li] = f; ctx.flen[j][li] = flen;
        if (k.type == kSkStr) {
          uint64_t pre_hi, pre_lo;
          StrPrefix16(f, flen, &pre_hi, &pre_lo);
          if (flen > 16) r.flags |= 1U << j;
          put(w0, k.rev ? ~pre_hi : pre_hi);
          put(w0 + 1, k.rev ? ~pre_lo : pre_lo);
        } else {
          const uint32_t rank = ChrRank(f, flen);
          if (rank == UINT32_MAX) r.flags |= 1U << j;
          put(w0, k.rev ? ~static_cast<uint64_t>(rank) : rank);
        }
      }
    }
  });

  phase("keys");
  // Parallel sort: sort thread_ct chunks, then merge pairwise.
  auto less = [&ctx](const SortRec& x, const SortRec& y) { return SortLess(ctx, x, y); };
  uint32_t part_ct = (n >= 65536) ? thread_ct : 1;
  std::vector<uintptr_t> bounds(part_ct + 1);
  for (uint32_t t = 0; t <= part_ct; ++t) bounds[t] = n * t / part_ct;
  {
    std::vector<std::thread> threads;
    for (uint32_t t = 0; t < part_ct; ++t) {
      threads.emplace_back([&, t]() { std::sort(recs.begin() + bounds[t], recs.begin() + bounds[t + 1], less); });
    }
    for (std::thread& th : threads) th.join();
  }
  phase("sort");
  if (part_ct > 1) {
    std::vector<SortRec> tmp(n);
    std::vector<SortRec>* src = &recs;
    std::vector<SortRec>* dst = &tmp;
    // Each round merges adjacent pairs of runs.  Every merge is split across
    // all threads: output position k of a merge of runs A and B takes
    // exactly i elements from A, where i is found by binary search ("merge
    // path"); the comparator is a strict total order, so i is unique.
    while (bounds.size() > 2) {
      std::vector<uintptr_t> nb;
      std::vector<std::thread> threads;
      for (size_t b = 0; b + 1 < bounds.size(); b += 2) {
        nb.push_back(bounds[b]);
        if (b + 2 == bounds.size()) {
          std::copy(src->begin() + bounds[b], src->begin() + bounds[b + 1], dst->begin() + bounds[b]);
          continue;
        }
        const SortRec* a = src->data() + bounds[b];
        const uintptr_t na = bounds[b + 1] - bounds[b];
        const SortRec* bb = src->data() + bounds[b + 1];
        const uintptr_t nbb = bounds[b + 2] - bounds[b + 1];
        SortRec* out = dst->data() + bounds[b];
        auto corank = [=, &less](uintptr_t k) {
          uintptr_t lo = (k > nbb) ? k - nbb : 0;
          uintptr_t hi = (k < na) ? k : na;
          while (lo < hi) {
            const uintptr_t mid = (lo + hi) / 2;
            if (less(a[mid], bb[k - mid - 1])) lo = mid + 1;
            else hi = mid;
          }
          return lo;
        };
        const uintptr_t total = na + nbb;
        for (uint32_t t = 0; t < thread_ct; ++t) {
          threads.emplace_back([=, &less]() {
            const uintptr_t k0 = total * t / thread_ct, k1 = total * (t + 1) / thread_ct;
            const uintptr_t i0 = corank(k0), i1 = corank(k1);
            std::merge(a + i0, a + i1, bb + (k0 - i0), bb + (k1 - i1), out + k0, less);
          });
        }
      }
      nb.push_back(n);
      for (std::thread& th : threads) th.join();
      bounds.swap(nb);
      std::swap(src, dst);
    }
    if (src != &recs) recs.swap(*src);
  }

  phase("merge");
  Emit(header.data(), header.size());
  Flush();
  // Gather lines in sorted order, in batches: each batch's per-thread
  // segments are copied in parallel into one buffer, then written at once.
  {
    const uintptr_t kBatchRecs = 1 << 20;
    std::vector<char> batch;
    std::vector<uintptr_t> seg_off(thread_ct + 1);
    for (uintptr_t lo = 0; lo < n; lo += kBatchRecs) {
      const uintptr_t hi = std::min<uintptr_t>(n, lo + kBatchRecs);
      // Byte offsets of each thread's segment.
      std::vector<uintptr_t> seg_bytes(thread_ct, 0);
      {
        std::vector<std::thread> threads;
        for (uint32_t t = 0; t < thread_ct; ++t) {
          threads.emplace_back([&, t]() {
            const uintptr_t s0 = lo + (hi - lo) * t / thread_ct, s1 = lo + (hi - lo) * (t + 1) / thread_ct;
            uintptr_t bytes = 0;
            for (uintptr_t q = s0; q < s1; ++q) bytes += lens[recs[q].idx];
            seg_bytes[t] = bytes;
          });
        }
        for (std::thread& th : threads) th.join();
      }
      seg_off[0] = 0;
      for (uint32_t t = 0; t < thread_ct; ++t) seg_off[t + 1] = seg_off[t] + seg_bytes[t];
      if (batch.size() < seg_off[thread_ct]) batch.resize(seg_off[thread_ct]);
      {
        std::vector<std::thread> threads;
        for (uint32_t t = 0; t < thread_ct; ++t) {
          threads.emplace_back([&, t]() {
            const uintptr_t s0 = lo + (hi - lo) * t / thread_ct, s1 = lo + (hi - lo) * (t + 1) / thread_ct;
            char* w = batch.data() + seg_off[t];
            for (uintptr_t q = s0; q < s1; ++q) {
              const uint32_t idx = recs[q].idx;
              memcpy(w, lines[idx], lens[idx]);
              w += lens[idx];
            }
          });
        }
        for (std::thread& th : threads) th.join();
      }
      if (fwrite(batch.data(), 1, seg_off[thread_ct], stdout) != seg_off[thread_ct]) Die("write failure");
    }
  }
  phase("output");
  for (char* c : chunks) free(c);
}

// ---- uniq ----
//
// Hash-based, so input need not be sorted; output is in first-seen order.
// On sorted input this is identical to uniq (with -r for -d), and
// without flags it's awk '!seen[KEY]++'.

static void ModeUniq(int argc, char** argv) {
  Opts o;
  uint32_t count = 0, dups = 0, singles = 0;
  int i = 2;
  while (1) {
    i = ParseOpts(argc, argv, i, &o);
    if (i > 2 && !strcmp(argv[i - 1], "--")) break;
    if (i == argc) break;
    const char* a = argv[i];
    if (!strcmp(a, "-c")) count = 1;
    else if (!strcmp(a, "-r")) dups = 1;
    else if (!strcmp(a, "-u")) singles = 1;
    else break;
    ++i;
  }
  if (argc - i != 1 && argc - i != 2) DieUsage("ptxt uniq [options] [-c] [-r|-u] [COLS] FILE");
  if (dups && singles) Die("-r and -u are mutually exclusive");
  const uint32_t whole = (argc - i == 1);
  Input in; OpenInput(argv[argc - 1], whole ? nullptr : argv[i], &o, &in);
  const std::vector<uint32_t>& cols = in.cols;
  const uint32_t nk = cols.size();
  const char delim = o.delim;
  const uint32_t deferred = count || dups || singles;
  if (in.hline) {
    if (count) Emit("count\t", 6);
    EmitHeader(in);
  }
  StrSet keys;
  std::vector<uint64_t> counts;  // only if deferred
  // First line of each key, if deferred and the key isn't the whole line.
  Arena lines;
  std::vector<const char*> lptr;
  std::vector<uint32_t> llen;
  std::string kbuf;
  // Sorted input has long runs of equal keys; skip hashing those.  And while
  // keys arrive in ascending order, each new key can't have been seen
  // before, so the hash table is only built if the order ever breaks.
  const char* prev = nullptr; uint32_t prev_len = 0; uint32_t prev_idx = 0;
  uint32_t ascending = 1;
  char* line;
  while (!TextNextLine(&in.txs, &line)) {
    char* lnext = TextLineEnd(&in.txs);
    char* lend = lnext - 1;
    const char* k;
    uint32_t klen;
    if (whole) {
      k = line; klen = lend - line;
    } else if (nk == 1) {
      // A missing field is an empty key, as in awk.
      char* f = FieldStart(line, lend, cols[0], delim);
      if (f) { k = f; klen = FieldEnd(f, lend, delim) - f; }
      else { k = lend; klen = 0; }
    } else {
      kbuf.clear();
      for (uint32_t j = 0; j < nk; ++j) {
        if (j) kbuf.push_back('\0');
        char* f = FieldStart(line, lend, cols[j], delim);
        if (f) kbuf.append(f, FieldEnd(f, lend, delim) - f);
      }
      k = kbuf.data(); klen = kbuf.size();
    }
    uint32_t idx;
    if (prev && klen == prev_len && !memcmp(k, prev, klen)) {
      idx = prev_idx;
    } else {
      uint32_t added = 1;
      if (ascending && prev && StrCmp3(k, klen, prev, prev_len) < 0) {
        ascending = 0;
        keys.BuildIndex();
      }
      if (ascending) idx = keys.AppendUnindexed(k, klen);
      else idx = keys.Insert(k, klen, HashStr(k, klen), &added);
      prev = keys.Ptr(idx); prev_len = klen; prev_idx = idx;
      if (added) {
        if (!deferred) {
          Emit(line, lnext - line);
          continue;
        }
        counts.push_back(0);
        if (!whole) {
          const uint32_t len = lend - line;
          char* c = lines.Alloc(len);
          memcpy(c, line, len);
          lptr.push_back(c); llen.push_back(len);
        }
      }
    }
    if (deferred) ++counts[idx];
  }
  FinishTs(&in.txs);
  if (!deferred) return;
  for (uint32_t j = 0; j < keys.Size(); ++j) {
    const uint64_t ct = counts[j];
    if ((dups && ct < 2) || (singles && ct > 1)) continue;
    if (count) {
      char numbuf[24];
      const int nlen = snprintf(numbuf, sizeof numbuf, "%llu\t", static_cast<unsigned long long>(ct));
      Emit(numbuf, nlen);
    }
    if (whole) Emit(keys.Ptr(j), keys.Len(j));
    else Emit(lptr[j], llen[j]);
    EmitCh('\n');
  }
}

int main(int argc, char** argv) {
  // 'ptxt' alone and 'ptxt SUB' alone print help and fail, as plink2 does
  // with no arguments; explicit help requests succeed.
  if (argc == 1) PrintSummaryHelp(1);
  if (!strcmp(argv[1], "-v") || !strcmp(argv[1], "-V") || !strcmp(argv[1], "--version")) {
    if (argc > 2) {
      fprintf(stderr, "ptxt: %s takes no arguments\n", argv[1]);
      return 1;
    }
    printf("%s\n", kVersion);
    return 0;
  }
  if (IsHelpFlag(argv[1]) || !strcmp(argv[1], "help")) {
    if (argc == 2) PrintSummaryHelp(0);
    const ModeHelp* h = FindModeHelp(argv[2]);
    if (!h) {
      fprintf(stderr, "ptxt: unrecognized subcommand '%s'; run 'ptxt' with no arguments for a list\n", argv[2]);
      return 1;
    }
    PrintModeHelp(h, 0);
  }
  g_mode = argv[1];
  const ModeHelp* mode_help = FindModeHelp(g_mode);
  if (!mode_help) {
    fprintf(stderr, "ptxt: unrecognized subcommand '%s'; run 'ptxt' with no arguments for a list\n", g_mode);
    return 1;
  }
  if (argc == 2) PrintModeHelp(mode_help, 1);
  if (IsHelpFlag(argv[2])) PrintModeHelp(mode_help, 0);
  g_outbuf = static_cast<char*>(malloc(kOutBufSize)); g_outw = g_outbuf;
  const char* mode = argv[1];
  if (!strcmp(mode, "cut")) ModeCut(argc, argv);
  else if (!strcmp(mode, "paste")) ModePaste(argc, argv);
  else if (!strcmp(mode, "select")) ModeSelect(argc, argv);
  else if (!strcmp(mode, "filter")) ModeFilter(argc, argv);
  else if (!strcmp(mode, "keep")) ModeKeep(argc, argv, 0);
  else if (!strcmp(mode, "remove")) ModeKeep(argc, argv, 1);
  else if (!strcmp(mode, "count")) ModeCount(argc, argv);
  else if (!strcmp(mode, "prefix")) ModePrefix(argc, argv, 0);
  else if (!strcmp(mode, "unprefix")) ModePrefix(argc, argv, 1);
  else if (!strcmp(mode, "tsv")) ModeTsv(argc, argv);
  else if (!strcmp(mode, "sum")) ModeSum(argc, argv);
  else if (!strcmp(mode, "sort")) ModeSort(argc, argv);
  else if (!strcmp(mode, "join")) ModeJoin(argc, argv);
  else if (!strcmp(mode, "uniq")) ModeUniq(argc, argv);
  else Die("unrecognized subcommand");
  Flush();
  return 0;
}
