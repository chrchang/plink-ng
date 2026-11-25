// Prototype plink2-library text tools.  Input may be uncompressed, gzip,
// BGZF (decompressed in parallel), or zstd.  Lines must end in LF.
//
// The user-facing documentation of every subcommand is the help text below
// (kSummaryHelp, kModeHelp; printed by 'ptxt' and 'ptxt SUB -h'), mirrored in
// README.md.  Keep the two in sync.
//
// Like the rest of plink2, this compiles as either C++ (the default) or C99
// ('make C99=1').  The only C++-specific code is std::sort in 'ptxt sort';
// the C build falls back on qsort() there, so 'ptxt sort' is slower in C.

// Correctly-rounded decimal parsing.  plink2's ScanadvDouble may be off by an
// ulp or two, which flips comparisons against thresholds like "<= 0.99" when
// a field is exactly on the boundary (e.g. "0.9900").
#define FFC_IMPL
#ifdef __cplusplus
// ffc.h is C; it uses C-style casts throughout.
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wold-style-cast"
#  pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#endif
#include "../../include/ffc.h"
#ifdef __cplusplus
#  pragma GCC diagnostic pop
#endif
#include "../../include/plink2_htable.h"
#include "../../include/plink2_text.h"
#include "../../include/plink2_thread.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <pthread.h>
#include <time.h>
#ifdef __cplusplus
#  include <algorithm>
using namespace plink2;
#endif

#define PTXT_NORETURN __attribute__((noreturn))

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

PTXT_NORETURN static void Die(const char* msg) { fprintf(stderr, "ptxt: %s\n", msg); exit(1); }

// ---- memory ----

static void* XMalloc(uintptr_t size) {
  void* p = malloc(size ? size : 1);
  if (!p) Die("out of memory");
  return p;
}

static void* XCalloc(uintptr_t ct, uintptr_t elem_size) {
  void* p = calloc(ct ? ct : 1, elem_size);
  if (!p) Die("out of memory");
  return p;
}

// Growable arrays: a struct with members p (element pointer), n (element
// count) and cap (capacity).  VEC_RESERVE and VEC_PUSH may move p.
#define VEC_DEF(tname, etype) typedef struct tname ## Struct { etype* p; uintptr_t n; uintptr_t cap; } tname

static void* VecGrowRaw(void* p, uintptr_t* capp, uintptr_t min_cap, uintptr_t elem_size) {
  uintptr_t cap = *capp * 2;
  if (cap < 16) cap = 16;
  if (cap < min_cap) cap = min_cap;
  void* r = realloc(p, cap * elem_size);
  if (!r) Die("out of memory");
  *capp = cap;
  return r;
}

#ifdef __cplusplus
template <typename T> static inline T* VecGrow(T* p, uintptr_t* capp, uintptr_t min_cap) {
  return static_cast<T*>(VecGrowRaw(p, capp, min_cap, sizeof(T)));
}
#  define VEC_RESERVE(v, ct) do { if ((ct) > (v).cap) (v).p = VecGrow((v).p, &(v).cap, (ct)); } while (0)
#else
#  define VEC_RESERVE(v, ct) do { if ((ct) > (v).cap) (v).p = VecGrowRaw((v).p, &(v).cap, (ct), sizeof(*(v).p)); } while (0)
#endif

#define VEC_PUSH(v, x) do { VEC_RESERVE((v), (v).n + 1); (v).p[(v).n++] = (x); } while (0)

VEC_DEF(U32Vec, uint32_t);
VEC_DEF(U64Vec, uint64_t);
VEC_DEF(CharVec, char);
VEC_DEF(CcpVec, const char*);
VEC_DEF(CharpVec, char*);

static inline void CvAppend(CharVec* v, const char* s, uintptr_t len) {
  VEC_RESERVE(*v, v->n + len);
  memcpy(&(v->p[v->n]), s, len);
  v->n += len;
}

static inline void CvPushCh(CharVec* v, char c) { VEC_PUSH(*v, c); }

// Simple bump allocator.
typedef struct ArenaStruct {
  CharpVec chunks;
  char* cur;
  uintptr_t left;
  uintptr_t chunk_size;
} Arena;

static void InitArena(uintptr_t chunk_size, Arena* a) {
  a->chunks.p = nullptr; a->chunks.n = 0; a->chunks.cap = 0;
  a->cur = nullptr;
  a->left = 0;
  a->chunk_size = chunk_size;
}

static void FreeArena(Arena* a) {
  for (uintptr_t k = 0; k < a->chunks.n; ++k) free(a->chunks.p[k]);
  free(a->chunks.p);
}

static inline char* ArenaAlloc(uintptr_t len, Arena* a) {
  if (len > a->left) {
    const uintptr_t csize = len > a->chunk_size ? len : a->chunk_size;
    a->cur = S_CAST(char*, XMalloc(csize));
    VEC_PUSH(a->chunks, a->cur);
    a->left = csize;
  }
  char* r = a->cur;
  a->cur += len; a->left -= len;
  return r;
}

// Sorts arr[0..ct) and removes duplicates; returns the new count.
static uint32_t SortUniqU32(uint32_t* arr, uint32_t ct) {
  if (ct < 2) return ct;
  STD_SORT(ct, u32cmp, arr);
  uint32_t w = 1;
  for (uint32_t r = 1; r < ct; ++r) {
    if (arr[r] != arr[w - 1]) arr[w++] = arr[r];
  }
  return w;
}

// Index of the first element of sorted arr[0..ct) that is >= val.
static uint32_t LowerBoundU32(const uint32_t* arr, uint32_t ct, uint32_t val) {
  uint32_t lo = 0, hi = ct;
  while (lo < hi) {
    const uint32_t mid = (lo + hi) / 2;
    if (arr[mid] < val) lo = mid + 1;
    else hi = mid;
  }
  return lo;
}

// ---- threads ----
//
// Each task runs fn(shared, lo, hi, tidx) on its own thread.
typedef void (*RangeFn)(void* shared, uintptr_t lo, uintptr_t hi, uint32_t tidx);

typedef struct RangeTaskStruct {
  RangeFn fn;
  void* shared;
  uintptr_t lo;
  uintptr_t hi;
  uint32_t tidx;
} RangeTask;

static void* RangeTaskMain(void* arg) {
  RangeTask* t = S_CAST(RangeTask*, arg);
  t->fn(t->shared, t->lo, t->hi, t->tidx);
  return nullptr;
}

static void RunTasks(RangeTask* tasks, uint32_t ct) {
  pthread_t* tids = S_CAST(pthread_t*, XMalloc(ct * sizeof(pthread_t)));
  for (uint32_t t = 0; t < ct; ++t) {
    if (pthread_create(&(tids[t]), nullptr, RangeTaskMain, &(tasks[t]))) Die("failed to create thread");
  }
  for (uint32_t t = 0; t < ct; ++t) pthread_join(tids[t], nullptr);
  free(tids);
}

// Splits [lo, hi) evenly across thread_ct threads (always spawned).
static void SpawnRanges(uint32_t thread_ct, uintptr_t lo, uintptr_t hi, RangeFn fn, void* shared) {
  RangeTask tasks[16];
  for (uint32_t t = 0; t < thread_ct; ++t) {
    tasks[t].fn = fn; tasks[t].shared = shared; tasks[t].tidx = t;
    tasks[t].lo = lo + (hi - lo) * t / thread_ct;
    tasks[t].hi = lo + (hi - lo) * (t + 1) / thread_ct;
  }
  RunTasks(tasks, thread_ct);
}

// Like SpawnRanges over [0, n), but runs on the calling thread when n is small.
static void ParallelFor(uint32_t thread_ct, uintptr_t n, RangeFn fn, void* shared) {
  if (thread_ct <= 1 || n < 65536) { fn(shared, 0, n, 0); return; }
  SpawnRanges(thread_ct, 0, n, fn, shared);
}

// ---- help ----
//
// 'ptxt' alone prints the subcommand summary; 'ptxt SUB -h' (or 'ptxt help
// SUB') prints everything about one subcommand; 'ptxt -v' prints kVersion.
// Keep lines <= 80 columns.

static const char kVersion[] = "ptxt v0.1.0 (21 Nov 2026)";

static const char kSummaryHelp[] = "\n"
"  ptxt <subcommand> [options...] <arguments...>\n"
"  ptxt <subcommand> -h\n"
"  ptxt -v\n"
"\n"
"Fast replacements for common text-processing one-liners.  Input may be\n"
"uncompressed, gzip, BGZF (decompressed on several threads) or zstd; use\n"
"/dev/stdin as FILE to read from a pipe.  Output goes to stdout.  Columns are\n"
"1-based numbers or header names.\n"
"\n"
"Subcommands:\n"
"  -- Columns\n"
"     select    columns in any order, repeats allowed    awk '{print $3,$1}'\n"
"     cut       columns in ascending order (tab only)    cut -f\n"
"     paste     two files side by side                   paste\n"
"     prefix    add a leading string to a column         awk '{$1=\"chr\"$1}1'\n"
"     unprefix  remove a leading string from a column    sed 's/^chr//'\n"
"     tsv       whitespace-delimited to tab-delimited    awk '{$1=$1}1' OFS=\\t\n"
"\n"
"  -- Rows\n"
"     filter    lines satisfying an expression           awk '$P < 5e-8'\n"
"     keep      lines whose column is in an ID list      grep -F -w -f\n"
"     remove    lines whose column isn't in an ID list   grep -v -F -w -f\n"
"     uniq      first line per distinct key              awk '!seen[$3]++'\n"
"     sort      stable multi-key sort                    sort -s -k1,1V -k2,2n\n"
"     join      hash join on key column(s)               join\n"
"\n"
"  -- Summaries\n"
"     count     number of lines per distinct value       cut | sort | uniq -c\n"
"     sum       sum of a numeric column                  awk '{s+=$5} END...'\n"
"\n"
"\"ptxt <subcommand> -h\" describes one subcommand in full; \"ptxt -v\" (or -V)\n"
"prints just the version.\n";

// Options shared by most subcommands; printed after a subcommand's own help.
static const char kColumnsHelp[] = "\n"
"Columns are 1-based numbers (e.g. 3) or header names (e.g. P); using a name\n"
"implies -H.\n";

static const char kInputHelp[] = "\n"
"Input may be uncompressed, gzip, BGZF or zstd; use /dev/stdin to read from a\n"
"pipe.  Lines must end in LF.\n";

typedef struct ModeHelpStruct {
  const char* name;
  const char* body;        // usage line(s) and description
  const char* examples;
  const char* header_doc;  // what -H does; nullptr if [opts] isn't accepted
  uint32_t odelim;         // -O accepted
  uint32_t ranges;         // column ranges accepted
} ModeHelp;

static const ModeHelp kModeHelp[] = {
  {"cut", "\n"
"  ptxt cut LIST FILE\n"
"\n"
"Same as 'cut -f LIST' on tab-delimited input.  LIST is a comma-separated list\n"
"of 1-based column numbers and ranges, in ascending order (e.g. 1-8 or 2,11).\n"
"Lines with no tab are printed whole, as cut does.  Contiguous ranges are\n"
"copied as one block, so this is fastest when you keep runs of columns.\n"
"\n"
"Use 'ptxt select' to reorder or repeat columns, pick them by header name, or\n"
"read whitespace-delimited input.\n"
"\n", "\n"
"Example:\n"
"  ptxt cut 1-5,8 calls.vcf.gz\n", nullptr, 0, 0},

  {"paste", "\n"
"  ptxt paste FILE1 FILE2\n"
"\n"
"Same as 'paste FILE1 FILE2': prints each FILE1 line, a tab, then the\n"
"corresponding FILE2 line.  It is an error for the files to have different\n"
"numbers of lines.\n"
"\n", "\n"
"Example:\n"
"  ptxt paste ids.txt scores.txt\n", nullptr, 0, 0},

  {"select", "\n"
"  ptxt select [options] COLS FILE\n"
"\n"
"Prints the COLS columns of each line, like awk '{print $3, $1, ...}'.  COLS\n"
"is a comma-separated list; columns can be in any order and can be repeated.\n"
"Missing fields are printed as empty strings.\n"
"\n", "\n"
"Examples:\n"
"  ptxt select SNP,P,BETA ss.tsv\n"
"  ptxt select -w -O , 2,1,5-7 plink.assoc\n", "It is printed like any other line.", 1, 1},

  {"filter", "\n"
"  ptxt filter [options] EXPR FILE\n"
"\n"
"Prints the lines that satisfy EXPR, an awk-like expression.  All arguments\n"
"before FILE are joined with spaces, so EXPR doesn't have to be one shell\n"
"argument (but quote it anyway, to protect $, <, >, ( and &).\n"
"\n"
"Syntax:\n"
"  Fields      $N, $NAME, or $\"NAME\" (names imply -H)\n"
"  Literals    numbers (any decimal or scientific notation), \"strings\"\n"
"  Operators   in increasing precedence:\n"
"                ||  (also |, or)\n"
"                &&  (also &, and)\n"
"                <  <=  >  >=  ==  (also =)  !=\n"
"                +  -\n"
"                *  /  %\n"
"                unary !  (also not),  unary -  +\n"
"  Functions   abs  log  log10  exp  sqrt\n"
"  Grouping    ( )\n"
"\n"
"How values compare (as in awk, except for missing values):\n"
"  - A field is a number if the whole field parses as one.\n"
"  - Field vs \"string\" compares as strings.\n"
"  - Field vs field compares numerically if both are numbers, as strings\n"
"    otherwise.\n"
"  - Field vs number compares numerically if the field is a number.\n"
"    Otherwise the comparison is false, and != is true.  (This is plink2's\n"
"    --extract-if-info convention; awk would compare as strings.)\n"
"  - Arithmetic on a non-numeric field gives a missing value: every\n"
"    comparison against it is false (!= is true), and it is false as a\n"
"    condition.\n"
"  - A bare value used as a condition is true if it is a nonzero number or a\n"
"    nonempty string.\n"
"  - A field past the end of the line is an empty string.\n"
"  - Numbers are parsed with correct rounding, so boundary comparisons such as\n"
"    0.9900 <= 0.99 agree with awk.\n"
"\n", "\n"
"Examples:\n"
"  ptxt filter '$P < 5e-8' ss.tsv\n"
"  ptxt filter '$P < 5e-8 && ($1 == \"22\" || $CHR == \"X\")' ss.tsv.gz\n"
"  ptxt filter 'abs($BETA / $SE) > 5 && $A1_FREQ >= 0.01' ss.tsv\n"
"  ptxt filter -w '-log10($9) > 8' plink.assoc\n", "It is copied through.", 0, 0},

  {"keep", "\n"
"  ptxt keep [options] COL IDFILE FILE\n"
"  ptxt remove [options] COL IDFILE FILE\n"
"\n"
"'keep' prints the lines whose COL value is listed in IDFILE; 'remove' prints\n"
"the lines whose COL value isn't.  Each ID is the first whitespace-delimited\n"
"token on an IDFILE line; blank IDFILE lines are ignored.  Matching is exact\n"
"string equality, via a hash table, so IDFILE can hold millions of IDs.\n"
"\n"
"This replaces awk 'NR==FNR {a[$1]; next} $3 in a' and 'grep -F -w -f', but\n"
"only looks at one column (grep -w matches anywhere on the line).\n"
"\n", "\n"
"Examples:\n"
"  ptxt keep SNP hits.txt ss.tsv.gz\n"
"  ptxt remove -w 2 bad_ids.txt plink.assoc\n", "It is copied through.", 0, 0},

  {"count", "\n"
"  ptxt count [options] COL FILE\n"
"\n"
"Prints \"VALUE<TAB>COUNT\" for each distinct value in column COL, in the order\n"
"the values first appear.  This replaces 'cut -f COL | sort | uniq -c' and\n"
"awk '{c[$1]++} END {for (k in c) print k, c[k]}'; the input doesn't need to\n"
"be sorted.  Lines that are missing COL aren't counted.\n"
"\n"
"Pipe the output to 'ptxt sort 2:nr /dev/stdin' to order it by count.\n"
"\n", "\n"
"Example:\n"
"  ptxt count CHR ss.tsv\n", "It is skipped.", 0, 0},

  {"prefix", "\n"
"  ptxt prefix [options] COL STR FILE\n"
"  ptxt unprefix [options] COL STR FILE\n"
"\n"
"'prefix' adds STR to the start of column COL.  'unprefix' removes STR from\n"
"the start of COL wherever it is present.  The rest of each line, including\n"
"delimiters, is copied unchanged.\n"
"\n", "\n"
"Examples:\n"
"  ptxt prefix 1 chr calls_b37.tsv\n"
"  ptxt unprefix CHR chr ss.tsv.gz\n", "It is copied through.", 0, 0},

  {"tsv", "\n"
"  ptxt tsv FILE\n"
"\n"
"Converts whitespace-delimited lines, such as PLINK 1 reports, to\n"
"tab-delimited.  Leading and trailing spaces and tabs are dropped and each run\n"
"of spaces and tabs becomes one tab, as with awk '{$1=$1}1' OFS='\\t'.  This is\n"
"what 'tr -s \" \" \"\\t\"' is usually meant to do, except that tr leaves a leading\n"
"tab on padded lines.\n"
"\n", "\n"
"Example:\n"
"  ptxt tsv plink.assoc > plink.assoc.tsv\n", nullptr, 0, 0},

  {"sum", "\n"
"  ptxt sum [options] COL FILE\n"
"\n"
"Prints the sum of the numeric values in column COL.  Non-numeric and missing\n"
"values are skipped.\n"
"\n", "\n"
"Example:\n"
"  ptxt sum -w 5 plink.frq.counts\n", "It is skipped.", 0, 0},

  {"sort", "\n"
"  ptxt sort [options] KEY[:TYPE][,KEY[:TYPE]...] FILE\n"
"\n"
"Stable sort on one or more keys: lines with equal keys stay in input order,\n"
"as with 'sort -s'.  Each KEY is a single column; later keys break ties in\n"
"earlier ones (at most 32 keys).  TYPE is one of:\n"
"  s   string, in byte order (default)\n"
"  n   number, in any decimal or scientific notation ('g' is an alias).\n"
"      Non-numeric and missing values always sort last, in either direction.\n"
"  c   chromosome: 1, 2, ..., 22 (and higher numbers), X, Y, XY, MT (or M),\n"
"      with or without a chr prefix (any case), then other contigs in byte\n"
"      order.\n"
"Append r to a TYPE (or use just r, for a reversed string key) to sort that\n"
"key in descending order.\n"
"\n"
"The whole input is loaded into memory.  The sort uses all cores, up to 16;\n"
"set PTXT_TIMING=1 to print per-phase timings to stderr.\n"
"\n", "\n"
"Examples:\n"
"  ptxt sort CHR:c,POS:n ss.tsv > ss.sorted.tsv\n"
"  ptxt sort P:n ss.tsv.gz | head\n"
"  ptxt sort -H 7:nr,3 ss.tsv\n", "It stays first.", 0, 0},

  {"join", "\n"
"  ptxt join [options] [-a] [-e STR] KEYS1 FILE1 KEYS2 FILE2\n"
"\n"
"Hash join.  For each FILE1 line, prints that line once for every FILE2 row\n"
"whose KEYS2 columns exactly match the line's KEYS1 columns, followed by the\n"
"FILE2 row's non-key columns.\n"
"  - Output keeps FILE1's order; several matching FILE2 rows come out in\n"
"    FILE2 order.  Neither file needs to be sorted.\n"
"  - KEYS1 and KEYS2 are comma-separated lists of the same length, so\n"
"    composite keys (e.g. CHR,POS) work.  Ranges aren't allowed.\n"
"  - FILE2 is held in memory and FILE1 is streamed, so put the smaller file\n"
"    second.\n"
"  - Short FILE2 rows are padded with empty fields.  FILE2 rows missing a key\n"
"    column, and blank FILE2 lines, are skipped.\n"
"\n"
"Join options:\n"
"  -a        Also print FILE1 lines that have no match, padded with -e values.\n"
"  -e STR    Fill value for -a (default: NA).\n"
"\n"
"With -H (or any column name), both files have header lines.  The output\n"
"header is FILE1's header followed by FILE2's non-key column names.\n"
"\n", "\n"
"Examples:\n"
"  ptxt join SNP sumstats.tsv ID annot.tsv\n"
"  ptxt join -a -e . CHR,POS a.tsv chrom,pos b.tsv.gz\n", "Both files must have one.", 1, 0},

  {"uniq", "\n"
"  ptxt uniq [options] [-c] [-r|-u] [COLS] FILE\n"
"\n"
"Prints the first line for each distinct key, in the order keys first appear.\n"
"The key is the whole line, or the COLS columns if COLS is given (a\n"
"comma-separated list; a missing field counts as empty, as in awk).  The input\n"
"doesn't need to be sorted: with no flags this is awk '!seen[$3]++', and on\n"
"sorted input the output matches uniq's.\n"
"\n"
"Uniq options:\n"
"  -c        Put each key's count and a tab in front of its line.\n"
"  -r        Only print keys that occur more than once ('uniq -d'; -d is\n"
"            ptxt's delimiter option).\n"
"  -u        Only print keys that occur exactly once.\n"
"-c, -r and -u print their output at the end, keeping each key's first line in\n"
"memory.  Without them, lines are printed as soon as they're read.\n"
"\n"
"Sorted input is especially fast: while keys arrive in ascending order no\n"
"hash table is needed.\n"
"\n", "\n"
"Examples:\n"
"  ptxt uniq SNP ss.tsv                  # drop repeated SNP IDs\n"
"  ptxt uniq -c -r CHR,POS ss.tsv.gz     # list duplicated positions\n"
"  sort ids.txt | ptxt uniq /dev/stdin\n", "It is copied through, with\n"
        "            \"count<TAB>\" in front under -c.", 0, 0},
};

static const ModeHelp* FindModeHelp(const char* mode) {
  if (!strcmp(mode, "remove")) mode = "keep";
  else if (!strcmp(mode, "unprefix")) mode = "prefix";
  for (uint32_t k = 0; k < sizeof(kModeHelp) / sizeof(kModeHelp[0]); ++k) {
    if (!strcmp(kModeHelp[k].name, mode)) return &(kModeHelp[k]);
  }
  return nullptr;
}

static void PrintBanner() {
  // Right-align the tagline at column 80, as plink2 does.
  static const char kTagline[] = "text tools built on the plink2 library";
  printf("%s%*s\n", kVersion, S_CAST(int, 80 - strlen(kVersion)), kTagline);
}

PTXT_NORETURN static void PrintSummaryHelp(int exit_code) {
  PrintBanner();
  fputs(kSummaryHelp, stdout);
  exit(exit_code);
}

PTXT_NORETURN static void PrintModeHelp(const ModeHelp* h, int exit_code) {
  fputs(h->body + 1, stdout);  // skip the leading newline
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

static uint32_t IsHelpFlag(const char* a) { return !strcmp(a, "-h") || !strcmp(a, "--help"); }

// For argument-count errors: the one-line synopsis plus a pointer to -h.
PTXT_NORETURN static void DieUsage(const char* synopsis) {
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

static inline uint32_t IsTokChar(char c) { return S_CAST(unsigned char, c) > 32; }

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
    char* e = S_CAST(char*, memchr(fstart, delim, lend - fstart));
    return e ? e : lend;
  }
  return CurTokenEnd(fstart);
}

// ---- options and column specs ----
typedef struct OptsStruct {
  char delim;
  const char* odelim;
  uint32_t header;
  uint32_t left_outer;  // join -a
  const char* fill;     // join -e
} Opts;

static void InitOpts(Opts* o) {
  o->delim = '\t';
  o->odelim = "\t";
  o->header = 0;
  o->left_outer = 0;
  o->fill = "NA";
}

// Parses leading options; returns index of first positional arg.
// join_opts enables -a and -e (join only).
static int ParseOpts(int argc, char** argv, int start, Opts* o, uint32_t join_opts) {
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

typedef struct ColTokStruct { uint32_t a, b; const char* name; } ColTok;
VEC_DEF(ColTokVec, ColTok);

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
static uint32_t ParseCols(char* spec, ColTokVec* toks) {
  uint32_t has_names = 0;
  for (char* p = spec; ; ) {
    char* e = strchr(p, ',');
    uint32_t len = e ? e - p : strlen(p);
    if (!len) Die("empty column in list");
    ColTok t;
    t.a = 0; t.b = 0; t.name = nullptr;
    if (IsRangeTok(p, len)) {
      char* q;
      t.a = strtoul(p, &q, 10); t.b = t.a;
      if (*q == '-') t.b = strtoul(q + 1, &q, 10);
      if (!t.a || !t.b) Die("columns are 1-based");
    } else {
      if (e) *e = '\0';
      t.name = p; has_names = 1;
    }
    VEC_PUSH(*toks, t);
    if (!e) break;
    p = e + 1;
  }
  return has_names;
}

// Returns the 0-based index of the header column with the given name.
static uint32_t ResolveNamedCol(const char* name, char* hline, char* hlend, char delim) {
  const uint32_t nlen = strlen(name);
  uint32_t idx = 0;
  char* f = FieldStart(hline, hlend, 0, delim);
  for (; f; ++idx) {
    char* fe = FieldEnd(f, hlend, delim);
    if (S_CAST(uint32_t, fe - f) == nlen && !memcmp(f, name, nlen)) break;
    if (delim) f = (fe == hlend) ? nullptr : fe + 1;
    else { f = FirstNonTspace(fe); if (!IsTokChar(*f)) f = nullptr; }
  }
  if (!f) { fprintf(stderr, "ptxt: column '%s' not found in header\n", name); exit(1); }
  return idx;
}

// Expands tokens into 0-based column indices (appended to *cols), resolving
// names against the header line (may be nullptr if no names are present).
static void ResolveCols(const ColTok* toks, uintptr_t tok_ct, char* hline, char* hlend, char delim, U32Vec* cols) {
  for (uintptr_t k = 0; k < tok_ct; ++k) {
    const ColTok* t = &(toks[k]);
    if (!t->name) {
      if (t->a <= t->b) for (uint32_t c = t->a; c <= t->b; ++c) VEC_PUSH(*cols, c - 1);
      else for (uint32_t c = t->a; c >= t->b; --c) VEC_PUSH(*cols, c - 1);
      continue;
    }
    VEC_PUSH(*cols, ResolveNamedCol(t->name, hline, hlend, delim));
  }
}

// Opens FILE, reads the header line if needed, and resolves COLS.
typedef struct InputStruct {
  TextStream txs;
  char* hline;
  char* hlend;
  U32Vec cols;
} Input;

static void OpenInput(const char* fname, char* colspec, Opts* o, Input* in) {
  ColTokVec toks = {nullptr, 0, 0};
  in->hline = nullptr;
  in->hlend = nullptr;
  in->cols.p = nullptr; in->cols.n = 0; in->cols.cap = 0;
  if (colspec && ParseCols(colspec, &toks)) o->header = 1;
  OpenTs(fname, &in->txs);
  if (o->header) {
    if (TextNextLine(&in->txs, &in->hline)) Die("missing header line");
    in->hlend = TextLineEnd(&in->txs) - 1;
  }
  if (colspec) ResolveCols(toks.p, toks.n, in->hline, in->hlend, o->delim, &in->cols);
  free(toks.p);
}

static uint32_t SingleCol(const Input* in) {
  if (in->cols.n != 1) Die("exactly one column expected");
  return in->cols.p[0];
}

static inline void EmitHeader(const Input* in) {
  if (in->hline) Emit(in->hline, in->hlend + 1 - in->hline);
}

// ---- string set (open addressing, stores lengths so no overread) ----
//
// Each slot holds (hash << 32 | index), so mismatched probes are rejected
// without touching the entry; entries hold (pointer, length, hash) together
// so a hash match costs one more cache line, not three.
typedef struct StrEntryStruct {
  const char* p;
  uint32_t len;
  uint32_t h;
} StrEntry;
VEC_DEF(StrEntryVec, StrEntry);

static const uint64_t kStrSetEmpty = ~0ULL;

typedef struct StrSetStruct {
  uint64_t* table;
  uint32_t table_size;
  StrEntryVec entries;
  Arena arena;
} StrSet;

// Allocates a table of the given size, all slots empty.
static void StrSetResetTable(uint32_t size, StrSet* ss) {
  free(ss->table);
  ss->table = S_CAST(uint64_t*, XMalloc(S_CAST(uintptr_t, size) * sizeof(uint64_t)));
  memset(ss->table, 0xff, S_CAST(uintptr_t, size) * sizeof(uint64_t));
  ss->table_size = size;
}

static void InitStrSet(StrSet* ss) {
  ss->table = nullptr;
  StrSetResetTable(1 << 16, ss);
  ss->entries.p = nullptr; ss->entries.n = 0; ss->entries.cap = 0;
  InitArena(1 << 20, &ss->arena);
}

static void FreeStrSet(StrSet* ss) {
  free(ss->table);
  free(ss->entries.p);
  FreeArena(&ss->arena);
}

static inline uint32_t StrSetSize(const StrSet* ss) { return ss->entries.n; }
static inline const char* StrSetPtr(const StrSet* ss, uint32_t idx) { return ss->entries.p[idx].p; }
static inline uint32_t StrSetLen(const StrSet* ss, uint32_t idx) { return ss->entries.p[idx].len; }

static inline uint32_t StrSetSlot(uint32_t h, uint32_t size) {
  return (S_CAST(uint64_t, h) * size) >> 32;
}

// Returns index, or UINT32_MAX.
static inline uint32_t StrSetFind(const StrSet* ss, const char* s, uint32_t len, uint32_t h) {
  const uint32_t size = ss->table_size;
  const uint64_t* table = ss->table;
  for (uint32_t slot = StrSetSlot(h, size); ; ) {
    const uint64_t v = table[slot];
    if (v == kStrSetEmpty) return UINT32_MAX;
    if (S_CAST(uint32_t, v >> 32) == h) {
      const uint32_t idx = S_CAST(uint32_t, v);
      const StrEntry* e = &(ss->entries.p[idx]);
      if (e->len == len && !memcmp(e->p, s, len)) return idx;
    }
    if (++slot == size) slot = 0;
  }
}

static void StrSetPlace(uint32_t idx, StrSet* ss) {
  const uint32_t size = ss->table_size;
  const uint32_t h = ss->entries.p[idx].h;
  uint64_t* table = ss->table;
  for (uint32_t slot = StrSetSlot(h, size); ; ) {
    if (table[slot] == kStrSetEmpty) { table[slot] = (S_CAST(uint64_t, h) << 32) | idx; return; }
    if (++slot == size) slot = 0;
  }
}

// Copies s into the arena and appends an entry for it; returns its index.
static uint32_t StrSetAppend(const char* s, uint32_t len, uint32_t h, StrSet* ss) {
  char* dst = ArenaAlloc(len, &ss->arena);
  memcpy(dst, s, len);
  const uintptr_t idx = ss->entries.n;
  if (idx == UINT32_MAX) Die("too many distinct strings");
  StrEntry e;
  e.p = dst; e.len = len; e.h = h;
  VEC_PUSH(ss->entries, e);
  return idx;
}

// Returns index of s, inserting it if absent.  *added set to 1 iff new.
static uint32_t StrSetInsert(const char* s, uint32_t len, uint32_t h, StrSet* ss, uint32_t* added) {
  uint32_t idx = StrSetFind(ss, s, len, h);
  *added = (idx == UINT32_MAX);
  if (!*added) return idx;
  idx = StrSetAppend(s, len, h, ss);
  if (ss->entries.n * 4 > ss->table_size) {
    StrSetResetTable(ss->table_size * 2, ss);
    for (uint32_t i = 0; i < ss->entries.n; ++i) StrSetPlace(i, ss);
  } else {
    StrSetPlace(idx, ss);
  }
  return idx;
}

// Appends s without indexing it; the caller guarantees it's new.  Used while
// input is sorted.  StrSetBuildIndex() must be called before the next
// Find/Insert.
static uint32_t StrSetAppendUnindexed(const char* s, uint32_t len, StrSet* ss) {
  return StrSetAppend(s, len, 0, ss);
}

static void StrSetBuildIndex(StrSet* ss) {
  uintptr_t size = ss->table_size;
  while (ss->entries.n * 4 > size) size *= 2;
  StrSetResetTable(size, ss);
  for (uint32_t i = 0; i < ss->entries.n; ++i) {
    ss->entries.p[i].h = Hash32(ss->entries.p[i].p, ss->entries.p[i].len);
    StrSetPlace(i, ss);
  }
}

static inline uint32_t HashStr(const char* s, uint32_t len) { return Hash32(s, len); }

// ---- modes ----
static void ModeCut(int argc, char** argv) {
  if (argc != 4) DieUsage("ptxt cut LIST FILE");
  U32Vec cols = {nullptr, 0, 0};
  for (char* p = argv[2]; *p; ) {
    uint32_t a = strtoul(p, &p, 10), b = a;
    if (*p == '-') b = strtoul(p + 1, &p, 10);
    if (!a || b < a) Die("invalid LIST");
    for (uint32_t c = a; c <= b; ++c) VEC_PUSH(cols, c - 1);
    if (*p == ',') ++p;
  }
  // Merge into contiguous ranges [rstart[r], rend[r]] (0-based, inclusive).
  U32Vec rstart = {nullptr, 0, 0};
  U32Vec rend = {nullptr, 0, 0};
  for (uintptr_t k = 0; k < cols.n; ++k) {
    const uint32_t c = cols.p[k];
    if (rstart.n && c <= rend.p[rend.n - 1] + 1) {
      if (c > rend.p[rend.n - 1]) rend.p[rend.n - 1] = c;
    } else {
      VEC_PUSH(rstart, c); VEC_PUSH(rend, c);
    }
  }
  const uint32_t nr = rstart.n;
  TextStream txs; OpenTs(argv[3], &txs);
  char* line;
  while (!TextNextLine(&txs, &line)) {
    char* lend = TextLineEnd(&txs) - 1;
    // 'it' points at the start of field 'cur'.
    char* it = line; uint32_t cur = 0;
    for (uint32_t r = 0; r < nr; ++r) {
      uint32_t skip = rstart.p[r] - cur;
      if (skip) {
        it = AdvToNthDelimChecked(it, lend, skip, '\t');
        if (!it) {
          if (!r && !memchr(line, '\t', lend - line)) Emit(line, lend - line);
          break;
        }
        ++it;
      }
      // Find the delimiter terminating field rend[r]; copy the whole span.
      char* fend = AdvToNthDelimChecked(it, lend, rend.p[r] - rstart.p[r] + 1, '\t');
      if (!fend) fend = lend;
      if (r) EmitCh('\t');
      Emit(it, fend - it);
      if (fend == lend) break;
      it = fend + 1; cur = rend.p[r] + 1;
    }
    EmitCh('\n');
  }
  FinishTs(&txs);
  free(cols.p);
  free(rstart.p);
  free(rend.p);
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
  Opts o; InitOpts(&o);
  int i = ParseOpts(argc, argv, 2, &o, 0);
  if (argc - i != 2) DieUsage("ptxt select [options] COLS FILE");
  Input in; OpenInput(argv[i + 1], argv[i], &o, &in);
  const uint32_t* req = in.cols.p;
  const uint32_t req_ct = in.cols.n;
  // Sorted unique needed columns, and map from request slot to needed slot.
  uint32_t* need = S_CAST(uint32_t*, XMalloc(req_ct * sizeof(int32_t)));
  memcpy(need, req, req_ct * sizeof(int32_t));
  const uint32_t nn = SortUniqU32(need, req_ct);
  uint32_t* slot_of = S_CAST(uint32_t*, XMalloc(req_ct * sizeof(int32_t)));
  for (uint32_t k = 0; k < req_ct; ++k) slot_of[k] = LowerBoundU32(need, nn, req[k]);
  char** fs = S_CAST(char**, XMalloc(nn * sizeof(intptr_t)));
  char** fe = S_CAST(char**, XMalloc(nn * sizeof(intptr_t)));
  const char delim = o.delim;
  const char* odelim = o.odelim;
  const uint32_t odlen = strlen(odelim);
  char* line = in.hline;
  char* lend = in.hlend;
  if (!line) {
    if (TextNextLine(&in.txs, &line)) line = nullptr;
    else lend = TextLineEnd(&in.txs) - 1;
  }
  while (line) {
    // Locate needed fields in one left-to-right pass.
    uint32_t found = 0;
    char* it = FieldStart(line, lend, 0, delim);
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
        fs[found] = it;
        fe[found] = FieldEnd(it, lend, delim);
      }
    }
    for (uint32_t k = 0; k < req_ct; ++k) {
      if (k) Emit(odelim, odlen);
      const uint32_t s = slot_of[k];
      if (s < found) Emit(fs[s], fe[s] - fs[s]);
    }
    EmitCh('\n');
    if (TextNextLine(&in.txs, &line)) break;
    lend = TextLineEnd(&in.txs) - 1;
  }
  FinishTs(&in.txs);
  free(in.cols.p);
  free(need);
  free(slot_of);
  free(fs);
  free(fe);
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

typedef enum {
  kEnNum, kEnStr, kEnField, kEnNeg, kEnNot, kEnAdd, kEnSub, kEnMul, kEnDiv,
  kEnMod, kEnLt, kEnLe, kEnGt, kEnGe, kEnEq, kEnNe, kEnAnd, kEnOr, kEnFunc,
  // Specialized comparisons, rewritten after parsing: field OP number and
  // field OP "string" (op stored in ExprNode::func).
  kEnCmpFN, kEnCmpFS
} ExprNodeType;

enum { kEfAbs, kEfLog, kEfLog10, kEfExp, kEfSqrt };

typedef struct ExprNodeStruct {
  ExprNodeType t;
  uint32_t a, b;     // child node indices
  double num;
  const char* str;   // string literal, or column name for kEnField
                     // (null-terminated)
  uint32_t slen;     // strlen(str)
  uint32_t col;      // 0-based column (kEnField), resolved after parse
  uint32_t slot;     // index into the sorted needed-column list
  uint32_t func;
} ExprNode;
VEC_DEF(ExprNodeVec, ExprNode);

typedef struct ExprParserStruct {
  const char* src;  // null-terminated
  uintptr_t len;
  uintptr_t p;
  ExprNodeVec* nodes;
} ExprParser;

PTXT_NORETURN static void ExprFail(const ExprParser* ps, const char* msg) {
  fprintf(stderr, "ptxt: invalid filter expression (%s):\n  %s\n  %*s^\n", msg, ps->src, S_CAST(int, ps->p), "");
  exit(1);
}

static void ExprSkipWs(ExprParser* ps) {
  while (ps->p < ps->len && S_CAST(unsigned char, ps->src[ps->p]) <= 32) ++ps->p;
}

static char ExprPeek(const ExprParser* ps, uintptr_t off) {
  return (ps->p + off < ps->len) ? ps->src[ps->p + off] : '\0';
}

// Matches a case-insensitive keyword not followed by an identifier char.
static uint32_t ExprKeyword(ExprParser* ps, const char* kw) {
  const uintptr_t len = strlen(kw);
  if (ps->p + len > ps->len || strncasecmp(&(ps->src[ps->p]), kw, len)) return 0;
  const char c = ExprPeek(ps, len);
  if (isalnum(S_CAST(unsigned char, c)) || c == '_') return 0;
  ps->p += len;
  return 1;
}

// Appends a zero-initialized node of type t; returns its index.
static uint32_t ExprNew(ExprParser* ps, ExprNodeType t) {
  ExprNode n;
  memset(&n, 0, sizeof(n));
  n.t = t;
  VEC_PUSH(*ps->nodes, n);
  return ps->nodes->n - 1;
}

static uint32_t ExprBin(ExprParser* ps, ExprNodeType t, uint32_t a, uint32_t b) {
  const uint32_t idx = ExprNew(ps, t);
  ps->nodes->p[idx].a = a;
  ps->nodes->p[idx].b = b;
  return idx;
}

static uint32_t ExprParseOr(ExprParser* ps);

// Assumes ExprPeek(ps, 0) == '"'.  Returns a malloc'd null-terminated copy
// (with \t and \n unescaped); its length is stored in *len_ptr.
static char* ExprParseQuoted(ExprParser* ps, uint32_t* len_ptr) {
  ++ps->p;
  CharVec s = {nullptr, 0, 0};
  while (1) {
    if (ps->p >= ps->len) ExprFail(ps, "unterminated string");
    char c = ps->src[ps->p++];
    if (c == '"') break;
    if (c == '\\') {
      if (ps->p >= ps->len) ExprFail(ps, "unterminated string");
      c = ps->src[ps->p++];
      if (c == 't') c = '\t';
      else if (c == 'n') c = '\n';
    }
    CvPushCh(&s, c);
  }
  *len_ptr = s.n;
  CvPushCh(&s, '\0');
  return s.p;
}

static uint32_t ExprParsePrimary(ExprParser* ps) {
  ExprSkipWs(ps);
  const char c = ExprPeek(ps, 0);
  if (c == '(') {
    ++ps->p;
    const uint32_t inner = ExprParseOr(ps);
    ExprSkipWs(ps);
    if (ExprPeek(ps, 0) != ')') ExprFail(ps, "expected ')'");
    ++ps->p;
    return inner;
  }
  if (c == '"') {
    uint32_t slen;
    char* str = ExprParseQuoted(ps, &slen);
    const uint32_t idx = ExprNew(ps, kEnStr);
    ps->nodes->p[idx].str = str;
    ps->nodes->p[idx].slen = slen;
    return idx;
  }
  if (c == '$') {
    ++ps->p;
    const char* str = nullptr;
    uint32_t slen = 0;
    uint32_t col = 0;
    if (ExprPeek(ps, 0) == '"') {
      str = ExprParseQuoted(ps, &slen);
    } else if (isdigit(S_CAST(unsigned char, ExprPeek(ps, 0)))) {
      uint64_t v = 0;
      while (isdigit(S_CAST(unsigned char, ExprPeek(ps, 0)))) {
        v = v * 10 + (ps->src[ps->p++] - '0');
        if (v > 0x7fffffff) ExprFail(ps, "column number too large");
      }
      if (!v) ExprFail(ps, "$0 is not supported; columns are 1-based");
      col = v - 1;
    } else {
      // Name runs until an operator, paren, quote, comma, or whitespace.
      const uintptr_t start = ps->p;
      while (ps->p < ps->len && !strchr(" \t\n\r!\"&()<=>|~+-*/%,", ps->src[ps->p])) ++ps->p;
      if (ps->p == start) ExprFail(ps, "expected column number or name after '$'");
      slen = ps->p - start;
      char* name = S_CAST(char*, XMalloc(slen + 1));
      memcpy(name, &(ps->src[start]), slen);
      name[slen] = '\0';
      str = name;
    }
    const uint32_t idx = ExprNew(ps, kEnField);
    ps->nodes->p[idx].str = str;
    ps->nodes->p[idx].slen = slen;
    ps->nodes->p[idx].col = col;
    return idx;
  }
  if (isdigit(S_CAST(unsigned char, c)) || (c == '.' && isdigit(S_CAST(unsigned char, ExprPeek(ps, 1))))) {
    double v;
    const char* start = &(ps->src[ps->p]);
    const ffc_result r = ffc_from_chars_double(start, ps->src + ps->len, &v);
    if (r.outcome != FFC_OUTCOME_OK) ExprFail(ps, "invalid number");
    const char* e = r.ptr;
    ps->p += e - start;
    if (isalpha(S_CAST(unsigned char, ExprPeek(ps, 0))) || ExprPeek(ps, 0) == '_') ExprFail(ps, "invalid number");
    const uint32_t idx = ExprNew(ps, kEnNum);
    ps->nodes->p[idx].num = v;
    return idx;
  }
  static const char* kFuncNames[] = {"abs", "log", "log10", "exp", "sqrt"};
  // Try longer names first so 'log10' isn't parsed as 'log'.
  static const uint32_t kFuncOrder[] = {kEfLog10, kEfAbs, kEfLog, kEfExp, kEfSqrt};
  for (uint32_t k = 0; k < sizeof(kFuncOrder) / sizeof(kFuncOrder[0]); ++k) {
    const uint32_t fi = kFuncOrder[k];
    const uintptr_t save = ps->p;
    if (ExprKeyword(ps, kFuncNames[fi])) {
      ExprSkipWs(ps);
      if (ExprPeek(ps, 0) != '(') { ps->p = save; break; }
      ++ps->p;
      const uint32_t arg = ExprParseOr(ps);
      ExprSkipWs(ps);
      if (ExprPeek(ps, 0) != ')') ExprFail(ps, "expected ')'");
      ++ps->p;
      const uint32_t idx = ExprBin(ps, kEnFunc, arg, 0);
      ps->nodes->p[idx].func = fi;
      return idx;
    }
  }
  if (isalpha(S_CAST(unsigned char, c)) || c == '_') {
    ExprFail(ps, "unknown identifier (refer to columns as $3 or $NAME, strings as \"...\")");
  }
  ExprFail(ps, "expected a value");
}

static uint32_t ExprParseUnary(ExprParser* ps) {
  ExprSkipWs(ps);
  const char c = ExprPeek(ps, 0);
  if ((c == '!' && ExprPeek(ps, 1) != '=') || ExprKeyword(ps, "not")) {
    if (c == '!') ++ps->p;
    return ExprBin(ps, kEnNot, ExprParseUnary(ps), 0);
  }
  if (c == '-' || c == '+') {
    ++ps->p;
    const uint32_t child = ExprParseUnary(ps);
    if (c == '+') return child;
    if (ps->nodes->p[child].t == kEnNum) {
      ps->nodes->p[child].num = -ps->nodes->p[child].num;
      return child;
    }
    return ExprBin(ps, kEnNeg, child, 0);
  }
  return ExprParsePrimary(ps);
}

static uint32_t ExprParseMul(ExprParser* ps) {
  uint32_t lhs = ExprParseUnary(ps);
  while (1) {
    ExprSkipWs(ps);
    const char c = ExprPeek(ps, 0);
    if (c != '*' && c != '/' && c != '%') return lhs;
    ++ps->p;
    const uint32_t rhs = ExprParseUnary(ps);
    lhs = ExprBin(ps, (c == '*') ? kEnMul : ((c == '/') ? kEnDiv : kEnMod), lhs, rhs);
  }
}

static uint32_t ExprParseAdd(ExprParser* ps) {
  uint32_t lhs = ExprParseMul(ps);
  while (1) {
    ExprSkipWs(ps);
    const char c = ExprPeek(ps, 0);
    if (c != '+' && c != '-') return lhs;
    ++ps->p;
    const uint32_t rhs = ExprParseMul(ps);
    lhs = ExprBin(ps, (c == '+') ? kEnAdd : kEnSub, lhs, rhs);
  }
}

static uint32_t ExprParseCmp(ExprParser* ps) {
  const uint32_t lhs = ExprParseAdd(ps);
  ExprSkipWs(ps);
  const char c0 = ExprPeek(ps, 0), c1 = ExprPeek(ps, 1);
  ExprNodeType t;
  if (c0 == '<') { t = (c1 == '=') ? kEnLe : kEnLt; ps->p += 1 + (c1 == '='); }
  else if (c0 == '>') { t = (c1 == '=') ? kEnGe : kEnGt; ps->p += 1 + (c1 == '='); }
  else if (c0 == '=') { t = kEnEq; ps->p += 1 + (c1 == '='); }
  else if (c0 == '!' && c1 == '=') { t = kEnNe; ps->p += 2; }
  else if (c0 == '~' || (c0 == '!' && c1 == '~')) ExprFail(ps, "regex matching is not supported");
  else return lhs;
  const uint32_t rhs = ExprParseAdd(ps);
  ExprSkipWs(ps);
  const char d0 = ExprPeek(ps, 0), d1 = ExprPeek(ps, 1);
  if (d0 == '<' || d0 == '>' || (d0 == '=' && d1 != '=') || (d0 == '=' && d1 == '=') || (d0 == '!' && d1 == '=')) {
    ExprFail(ps, "comparisons can't be chained; use && (e.g. 0.01 <= $6 && $6 <= 0.99)");
  }
  return ExprBin(ps, t, lhs, rhs);
}

static uint32_t ExprParseAnd(ExprParser* ps) {
  uint32_t lhs = ExprParseCmp(ps);
  while (1) {
    ExprSkipWs(ps);
    if (ExprPeek(ps, 0) == '&') {
      ps->p += (ExprPeek(ps, 1) == '&') ? 2 : 1;
    } else if (!ExprKeyword(ps, "and")) {
      return lhs;
    }
    const uint32_t rhs = ExprParseCmp(ps);
    lhs = ExprBin(ps, kEnAnd, lhs, rhs);
  }
}

static uint32_t ExprParseOr(ExprParser* ps) {
  uint32_t lhs = ExprParseAnd(ps);
  while (1) {
    ExprSkipWs(ps);
    if (ExprPeek(ps, 0) == '|') {
      ps->p += (ExprPeek(ps, 1) == '|') ? 2 : 1;
    } else if (!ExprKeyword(ps, "or")) {
      return lhs;
    }
    const uint32_t rhs = ExprParseAnd(ps);
    lhs = ExprBin(ps, kEnOr, lhs, rhs);
  }
}

// Parses src (null-terminated) into *nodes; returns the root index.
static uint32_t ExprParse(const char* src, uintptr_t len, ExprNodeVec* nodes) {
  ExprParser ps;
  ps.src = src; ps.len = len; ps.p = 0; ps.nodes = nodes;
  const uint32_t root = ExprParseOr(&ps);
  ExprSkipWs(&ps);
  if (ps.p != len) ExprFail(&ps, "unexpected input");
  return root;
}

// Evaluation state for one line.
typedef struct ExprLineStruct {
  // Configuration.
  const uint32_t* need;  // sorted unique 0-based columns
  char delim;
  // Per-line state.
  char* line;
  char* lend;
  uint32_t located;     // need[0..located) have been located
  uint32_t exhausted;   // line ended before need[located]
  char* it;             // start of field it_col, if located > 0
  uint32_t it_col;
  char** fs;
  char** fe;
  uint8_t* numstate;    // 0 = unparsed, 1 = numeric, 2 = not
  double* nums;
} ExprLine;

static void InitExprLine(const uint32_t* need, uint32_t nn, char delim, ExprLine* lc) {
  lc->need = need; lc->delim = delim;
  lc->fs = S_CAST(char**, XCalloc(nn, sizeof(intptr_t)));
  lc->fe = S_CAST(char**, XCalloc(nn, sizeof(intptr_t)));
  lc->numstate = S_CAST(uint8_t*, XCalloc(nn, 1));
  lc->nums = S_CAST(double*, XCalloc(nn, sizeof(double)));
}

static void FreeExprLine(ExprLine* lc) {
  free(lc->fs);
  free(lc->fe);
  free(lc->numstate);
  free(lc->nums);
}

static inline void ExprLineReset(char* line, char* lend, ExprLine* lc) {
  lc->line = line; lc->lend = lend; lc->located = 0; lc->exhausted = 0;
}

static void ExprLineLocateSlow(uint32_t slot, ExprLine* lc) {
  static char empty[1] = {'\0'};
  while (lc->located <= slot) {
    const uint32_t located = lc->located;
    const uint32_t col = lc->need[located];
    char* f = nullptr;
    if (!lc->exhausted) {
      if (!located) {
        f = FieldStart(lc->line, lc->lend, col, lc->delim);
      } else {
        const uint32_t skip = col - lc->it_col;
        if (lc->delim) {
          f = AdvToNthDelimChecked(lc->it, lc->lend, skip, lc->delim);
          if (f) ++f;
        } else {
          f = NextTokenMult(lc->it, skip);
        }
      }
    }
    if (f) {
      lc->it = f; lc->it_col = col;
      lc->fs[located] = f;
      lc->fe[located] = FieldEnd(f, lc->lend, lc->delim);
    } else {
      lc->exhausted = 1;
      lc->fs[located] = empty; lc->fe[located] = empty;
    }
    lc->numstate[located] = 0;
    lc->located = located + 1;
  }
}

// Locates need[slot] (and everything before it), scanning only as far as
// needed.  Missing fields become empty strings.
static inline void ExprLineLocate(uint32_t slot, ExprLine* lc) {
  if (slot < lc->located) return;
  ExprLineLocateSlow(slot, lc);
}

// Returns 1 and sets *dp iff field slot is numeric (parsed at most once).
static inline uint32_t ExprLineFieldNum(uint32_t slot, ExprLine* lc, double* dp) {
  ExprLineLocate(slot, lc);
  uint8_t st = lc->numstate[slot];
  if (!st) {
    st = ParseNumExact(lc->fs[slot], lc->fe[slot], &(lc->nums[slot])) ? 1 : 2;
    lc->numstate[slot] = st;
  }
  *dp = lc->nums[slot];
  return st == 1;
}

enum { kVkNum, kVkStr, kVkMissing };

typedef struct ValStruct {
  uint8_t k;         // kVk...
  uint8_t is_field;  // string form available from the line
  double num;
  const char* s;
  uint32_t len;
} Val;

static inline uint32_t Truthy(Val v) {
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
static inline uint32_t ToNum(Val v, double* dp) {
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

static Val EvalExpr(const ExprNode* nodes, uint32_t idx, ExprLine* lc) {
  const ExprNode* n = &(nodes[idx]);
  switch (n->t) {
    case kEnNum:
      return NumVal(n->num);
    case kEnStr: {
      Val v; v.k = kVkStr; v.is_field = 0; v.num = 0;
      v.s = n->str; v.len = n->slen;
      return v;
    }
    case kEnField: {
      const uint32_t slot = n->slot;
      double d;
      const uint32_t is_num = ExprLineFieldNum(slot, lc, &d);
      Val v; v.is_field = 1; v.s = lc->fs[slot];
      v.len = lc->fe[slot] - lc->fs[slot];
      if (is_num) { v.k = kVkNum; v.num = d; }
      else { v.k = kVkStr; v.num = 0; }
      return v;
    }
    case kEnCmpFN: {
      double d;
      if (!ExprLineFieldNum(n->slot, lc, &d) || d != d) return NumVal(n->func == kEnNe);
      return NumVal(CmpResult(n->func, (d > n->num) - (d < n->num)));
    }
    case kEnCmpFS: {
      const uint32_t slot = n->slot;
      ExprLineLocate(slot, lc);
      return NumVal(CmpResult(n->func, StrCmp3(lc->fs[slot], lc->fe[slot] - lc->fs[slot], n->str, n->slen)));
    }
    case kEnAnd:
      return NumVal(Truthy(EvalExpr(nodes, n->a, lc)) && Truthy(EvalExpr(nodes, n->b, lc)));
    case kEnOr:
      return NumVal(Truthy(EvalExpr(nodes, n->a, lc)) || Truthy(EvalExpr(nodes, n->b, lc)));
    case kEnNot:
      return NumVal(!Truthy(EvalExpr(nodes, n->a, lc)));
    case kEnNeg: {
      double d;
      if (!ToNum(EvalExpr(nodes, n->a, lc), &d)) return MissingVal();
      return NumVal(-d);
    }
    case kEnFunc: {
      double d;
      if (!ToNum(EvalExpr(nodes, n->a, lc), &d)) return MissingVal();
      switch (n->func) {
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
      if (!ToNum(EvalExpr(nodes, n->a, lc), &x)) return MissingVal();
      if (!ToNum(EvalExpr(nodes, n->b, lc), &y)) return MissingVal();
      switch (n->t) {
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
  const Val x = EvalExpr(nodes, n->a, lc);
  const Val y = EvalExpr(nodes, n->b, lc);
  int c;
  if (x.k == kVkMissing || y.k == kVkMissing) return NumVal(n->t == kEnNe);
  if (x.k == kVkNum && y.k == kVkNum) {
    if (x.num != x.num || y.num != y.num) return NumVal(n->t == kEnNe);  // NaN
    c = (x.num > y.num) - (x.num < y.num);
  } else if ((x.k == kVkNum && !x.is_field && y.is_field) || (y.k == kVkNum && !y.is_field && x.is_field)) {
    // Number vs non-numeric field: missing.
    return NumVal(n->t == kEnNe);
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
  switch (n->t) {
    case kEnLt: return NumVal(c < 0);
    case kEnLe: return NumVal(c <= 0);
    case kEnGt: return NumVal(c > 0);
    case kEnGe: return NumVal(c >= 0);
    case kEnEq: return NumVal(c == 0);
    default: return NumVal(c != 0);
  }
}

static void ModeFilter(int argc, char** argv) {
  Opts o; InitOpts(&o);
  int i = ParseOpts(argc, argv, 2, &o, 0);
  if (argc - i < 2) DieUsage("ptxt filter [options] EXPR FILE");
  // All args before FILE are joined with spaces, so quoting is optional
  // where the shell permits.
  CharVec src = {nullptr, 0, 0};
  for (int k = i; k < argc - 1; ++k) {
    if (k > i) CvPushCh(&src, ' ');
    CvAppend(&src, argv[k], strlen(argv[k]));
  }
  const uintptr_t src_len = src.n;
  CvPushCh(&src, '\0');
  ExprNodeVec nv = {nullptr, 0, 0};
  const uint32_t root = ExprParse(src.p, src_len, &nv);
  ExprNode* nodes = nv.p;
  const uint32_t node_ct = nv.n;
  uint32_t any_names = 0;
  for (uint32_t k = 0; k < node_ct; ++k) any_names |= (nodes[k].t == kEnField && nodes[k].slen);
  if (any_names) o.header = 1;
  Input in; OpenInput(argv[argc - 1], nullptr, &o, &in);
  const char delim = o.delim;
  // Resolve names; collect needed columns.
  uint32_t* need = S_CAST(uint32_t*, XMalloc(node_ct * sizeof(int32_t)));
  uint32_t need_ct = 0;
  for (uint32_t k = 0; k < node_ct; ++k) {
    ExprNode* n = &(nodes[k]);
    if (n->t != kEnField) continue;
    if (n->slen) n->col = ResolveNamedCol(n->str, in.hline, in.hlend, delim);
    need[need_ct++] = n->col;
  }
  need_ct = SortUniqU32(need, need_ct);
  for (uint32_t k = 0; k < node_ct; ++k) {
    if (nodes[k].t == kEnField) nodes[k].slot = LowerBoundU32(need, need_ct, nodes[k].col);
  }
  // Rewrite field-vs-literal comparisons into specialized nodes.  A literal on
  // the left is moved to the right by mirroring the operator.
  for (uint32_t k = 0; k < node_ct; ++k) {
    ExprNode* n = &(nodes[k]);
    if (n->t < kEnLt || n->t > kEnNe) continue;
    uint32_t fa = n->a, la = n->b;
    ExprNodeType op = n->t;
    if (nodes[fa].t != kEnField) {
      const uint32_t tmp = fa; fa = la; la = tmp;
      if (op == kEnLt) op = kEnGt;
      else if (op == kEnGt) op = kEnLt;
      else if (op == kEnLe) op = kEnGe;
      else if (op == kEnGe) op = kEnLe;
    }
    if (nodes[fa].t != kEnField) continue;
    const ExprNodeType lt = nodes[la].t;
    if (lt != kEnNum && lt != kEnStr) continue;
    n->slot = nodes[fa].slot;
    n->func = op;
    if (lt == kEnNum) {
      n->t = kEnCmpFN;
      n->num = nodes[la].num;
    } else {
      n->t = kEnCmpFS;
      n->str = nodes[la].str;
      n->slen = nodes[la].slen;
    }
  }
  ExprLine lc;
  InitExprLine(need, need_ct, delim, &lc);
  EmitHeader(&in);
  char* line;
  while (!TextNextLine(&in.txs, &line)) {
    char* lnext = TextLineEnd(&in.txs);
    ExprLineReset(line, lnext - 1, &lc);
    if (Truthy(EvalExpr(nodes, root, &lc))) Emit(line, lnext - line);
  }
  FinishTs(&in.txs);
  FreeExprLine(&lc);
  // Each literal/name string is owned by the kEnStr/kEnField node that parsed
  // it; rewritten kEnCmpFS nodes only borrow theirs.
  for (uint32_t k = 0; k < node_ct; ++k) {
    if (nodes[k].t == kEnStr || nodes[k].t == kEnField) free(K_CAST(char*, nodes[k].str));
  }
  free(nodes);
  free(need);
  free(src.p);
  free(in.cols.p);
}

static void ModeKeep(int argc, char** argv, uint32_t invert) {
  Opts o; InitOpts(&o);
  int i = ParseOpts(argc, argv, 2, &o, 0);
  if (argc - i != 3) DieUsage("ptxt keep|remove [options] COL IDFILE FILE");
  StrSet ids; InitStrSet(&ids);
  {
    TextStream its; OpenTs(argv[i + 1], &its);
    char* line;
    while (!TextNextLine(&its, &line)) {
      char* t = FirstNonTspace(line);
      if (!IsTokChar(*t)) continue;
      char* te = CurTokenEnd(t);
      uint32_t added;
      StrSetInsert(t, te - t, HashStr(t, te - t), &ids, &added);
    }
    FinishTs(&its);
    CleanupTextStream(&its, nullptr);
  }
  Input in; OpenInput(argv[i + 2], argv[i], &o, &in);
  const uint32_t col = SingleCol(&in);
  const char delim = o.delim;
  EmitHeader(&in);
  char* line;
  while (!TextNextLine(&in.txs, &line)) {
    char* lnext = TextLineEnd(&in.txs);
    char* lend = lnext - 1;
    char* f = FieldStart(line, lend, col, delim);
    uint32_t hit = 0;
    if (f) {
      char* fe = FieldEnd(f, lend, delim);
      const uint32_t len = fe - f;
      hit = (StrSetFind(&ids, f, len, HashStr(f, len)) != UINT32_MAX);
    }
    if (hit != invert) Emit(line, lnext - line);
  }
  FinishTs(&in.txs);
  FreeStrSet(&ids);
  free(in.cols.p);
}

static void ModeCount(int argc, char** argv) {
  Opts o; InitOpts(&o);
  int i = ParseOpts(argc, argv, 2, &o, 0);
  if (argc - i != 2) DieUsage("ptxt count [options] COL FILE");
  Input in; OpenInput(argv[i + 1], argv[i], &o, &in);
  const uint32_t col = SingleCol(&in);
  const char delim = o.delim;
  StrSet vals; InitStrSet(&vals);
  U64Vec counts = {nullptr, 0, 0};
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
      ++counts.p[prev_idx];
      continue;
    }
    uint32_t added;
    const uint32_t idx = StrSetInsert(f, len, HashStr(f, len), &vals, &added);
    if (added) VEC_PUSH(counts, 0);
    ++counts.p[idx];
    prev = StrSetPtr(&vals, idx); prev_len = len; prev_idx = idx;
  }
  FinishTs(&in.txs);
  for (uint32_t k = 0; k < StrSetSize(&vals); ++k) {
    char numbuf[24];
    const int nlen = snprintf(numbuf, sizeof numbuf, "\t%llu\n", S_CAST(unsigned long long, counts.p[k]));
    Emit(StrSetPtr(&vals, k), StrSetLen(&vals, k));
    Emit(numbuf, nlen);
  }
  FreeStrSet(&vals);
  free(counts.p);
  free(in.cols.p);
}

static void ModePrefix(int argc, char** argv, uint32_t remove) {
  Opts o; InitOpts(&o);
  int i = ParseOpts(argc, argv, 2, &o, 0);
  if (argc - i != 3) DieUsage("ptxt prefix|unprefix [options] COL STR FILE");
  const char* str = argv[i + 1];
  const uint32_t slen = strlen(str);
  Input in; OpenInput(argv[i + 2], argv[i], &o, &in);
  const uint32_t col = SingleCol(&in);
  const char delim = o.delim;
  EmitHeader(&in);
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
      if (S_CAST(uintptr_t, lnext - 1 - f) >= slen && !memcmp(f, str, slen)) {
        Emit(line, f - line);
        Emit(f + slen, lnext - f - slen);
        continue;
      }
    }
    Emit(line, lnext - line);
  }
  FinishTs(&in.txs);
  free(in.cols.p);
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
  Opts o; InitOpts(&o);
  int i = ParseOpts(argc, argv, 2, &o, 0);
  if (argc - i != 2) DieUsage("ptxt sum [options] COL FILE");
  Input in; OpenInput(argv[i + 1], argv[i], &o, &in);
  const uint32_t col = SingleCol(&in);
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
  free(in.cols.p);
  char buf[32];
  Emit(buf, snprintf(buf, sizeof buf, "%.10g\n", tot));
}

// ---- join ----
//
// FILE2 (the lookup table) is loaded into a hash table keyed by its key
// columns; each row's non-key columns are stored as a ready-to-emit payload.
// FILE1 is then streamed: each line is printed once per matching FILE2 row
// (in FILE2 order), followed by that row's payload.

// Splits a line into fields (all of them).  Absent trailing fields are just
// not present.
static void SplitFields(char* line, char* lend, char delim, CharpVec* fs, CharpVec* fe) {
  fs->n = 0; fe->n = 0;
  char* f = FieldStart(line, lend, 0, delim);
  while (f) {
    char* e = FieldEnd(f, lend, delim);
    VEC_PUSH(*fs, f); VEC_PUSH(*fe, e);
    if (delim) {
      f = (e == lend) ? nullptr : e + 1;
    } else {
      f = FirstNonTspace(e);
      if (!IsTokChar(*f)) f = nullptr;
    }
  }
}

// Appends the non-key fields of the split line, each preceded by odelim, to
// *dst.  Returns the number of fields appended.
static uint32_t JoinPayload(const CharpVec* fs, const CharpVec* fe, const uint8_t* is_key, uintptr_t is_key_n, const char* odelim, uint32_t odlen, CharVec* dst) {
  uint32_t ct = 0;
  for (uint32_t c = 0; c < fs->n; ++c) {
    if (c < is_key_n && is_key[c]) continue;
    CvAppend(dst, odelim, odlen);
    CvAppend(dst, fs->p[c], fe->p[c] - fs->p[c]);
    ++ct;
  }
  return ct;
}

static void ModeJoin(int argc, char** argv) {
  Opts o; InitOpts(&o);
  int i = ParseOpts(argc, argv, 2, &o, 1);
  if (argc - i != 4) DieUsage("ptxt join [options] [-a] [-e STR] KEYS1 FILE1 KEYS2 FILE2");
  ColTokVec toks1 = {nullptr, 0, 0};
  ColTokVec toks2 = {nullptr, 0, 0};
  if (ParseCols(argv[i], &toks1)) o.header = 1;
  if (ParseCols(argv[i + 2], &toks2)) o.header = 1;
  for (uint32_t tt = 0; tt < 2; ++tt) {
    const ColTokVec* t = tt ? &toks2 : &toks1;
    for (uintptr_t k = 0; k < t->n; ++k) {
      if (!t->p[k].name && t->p[k].a != t->p[k].b) Die("join keys must be single columns, not ranges");
    }
  }
  if (toks1.n != toks2.n) Die("KEYS1 and KEYS2 must have the same number of columns");
  const uint32_t nk = toks1.n;
  const char delim = o.delim;
  const char* odelim = o.odelim;
  const uint32_t odlen = strlen(odelim);

  // ---- load FILE2 ----
  TextStream t2; OpenTs(argv[i + 3], &t2);
  char* line;
  CharpVec fs = {nullptr, 0, 0};
  CharpVec fe = {nullptr, 0, 0};
  CharVec header2_payload = {nullptr, 0, 0};
  U32Vec kcols2 = {nullptr, 0, 0};
  uint8_t* is_key2 = nullptr;
  uint32_t is_key2_n = 0;
  uint32_t fill_ct = UINT32_MAX;  // payload column count (from first line)
  char* hline2 = nullptr;
  char* hlend2 = nullptr;
  if (o.header) {
    if (TextNextLine(&t2, &hline2)) Die("FILE2 is missing its header line");
    hlend2 = TextLineEnd(&t2) - 1;
  }
  ResolveCols(toks2.p, toks2.n, hline2, hlend2, delim, &kcols2);
  for (uint32_t k = 0; k < nk; ++k) {
    const uint32_t c = kcols2.p[k];
    if (c >= is_key2_n) is_key2_n = c + 1;
  }
  is_key2 = S_CAST(uint8_t*, XCalloc(is_key2_n, 1));
  for (uint32_t k = 0; k < nk; ++k) {
    const uint32_t c = kcols2.p[k];
    if (is_key2[c]) Die("duplicate column in KEYS2");
    is_key2[c] = 1;
  }
  if (o.header) {
    SplitFields(hline2, hlend2, delim, &fs, &fe);
    JoinPayload(&fs, &fe, is_key2, is_key2_n, odelim, odlen, &header2_payload);
    fill_ct = fs.n - nk;
  }
  StrSet keys; InitStrSet(&keys);
  U32Vec head = {nullptr, 0, 0};  // per distinct key: first/last row
  U32Vec tail = {nullptr, 0, 0};
  U32Vec next = {nullptr, 0, 0};  // per row: next row with same key
  CcpVec pay = {nullptr, 0, 0};
  U32Vec paylen = {nullptr, 0, 0};
  Arena arena; InitArena(16 << 20, &arena);
  CharVec keybuf = {nullptr, 0, 0};
  CharVec paybuf = {nullptr, 0, 0};
  while (!TextNextLine(&t2, &line)) {
    char* lend = TextLineEnd(&t2) - 1;
    if (lend == line) continue;  // blank line (not a record)
    SplitFields(line, lend, delim, &fs, &fe);
    if (!fs.n) continue;  // whitespace-only line under -w
    if (fill_ct == UINT32_MAX) fill_ct = (fs.n >= nk) ? fs.n - nk : 0;
    // Rows lacking a key column can't match anything.
    uint32_t ok = 1;
    for (uint32_t k = 0; k < nk; ++k) ok &= (kcols2.p[k] < fs.n);
    if (!ok) continue;
    const char* kp;
    uint32_t klen;
    if (nk == 1) {
      kp = fs.p[kcols2.p[0]]; klen = fe.p[kcols2.p[0]] - fs.p[kcols2.p[0]];
    } else {
      keybuf.n = 0;
      for (uint32_t j = 0; j < nk; ++j) {
        if (j) CvPushCh(&keybuf, '\0');
        CvAppend(&keybuf, fs.p[kcols2.p[j]], fe.p[kcols2.p[j]] - fs.p[kcols2.p[j]]);
      }
      kp = keybuf.p; klen = keybuf.n;
    }
    paybuf.n = 0;
    // Rows with fewer columns than the first are padded with empty fields,
    // so every output line has the same column count.
    for (uint32_t pc = JoinPayload(&fs, &fe, is_key2, is_key2_n, odelim, odlen, &paybuf); pc < fill_ct; ++pc) CvAppend(&paybuf, odelim, odlen);
    char* pdst = ArenaAlloc(paybuf.n, &arena);
    memcpy(pdst, paybuf.p, paybuf.n);
    const uintptr_t row = pay.n;
    if (row == UINT32_MAX) Die("FILE2 has too many rows");
    VEC_PUSH(pay, pdst); VEC_PUSH(paylen, paybuf.n); VEC_PUSH(next, UINT32_MAX);
    uint32_t added;
    const uint32_t kidx = StrSetInsert(kp, klen, HashStr(kp, klen), &keys, &added);
    if (added) {
      VEC_PUSH(head, row); VEC_PUSH(tail, row);
    } else {
      next.p[tail.p[kidx]] = row;
      tail.p[kidx] = row;
    }
  }
  FinishTs(&t2);
  CleanupTextStream(&t2, nullptr);
  if (fill_ct == UINT32_MAX) fill_ct = 0;
  CharVec fill = {nullptr, 0, 0};
  for (uint32_t k = 0; k < fill_ct; ++k) { CvAppend(&fill, odelim, odlen); CvAppend(&fill, o.fill, strlen(o.fill)); }

  // ---- stream FILE1 ----
  Input in; OpenInput(argv[i + 1], nullptr, &o, &in);
  U32Vec kcols1 = {nullptr, 0, 0};
  ResolveCols(toks1.p, toks1.n, in.hline, in.hlend, delim, &kcols1);
  if (in.hline) {
    Emit(in.hline, in.hlend - in.hline);
    Emit(header2_payload.p, header2_payload.n);
    EmitCh('\n');
  }
  // Locate key fields in one ascending pass (they may be listed in any
  // order).
  uint32_t* need = S_CAST(uint32_t*, XMalloc(nk * sizeof(int32_t)));
  memcpy(need, kcols1.p, nk * sizeof(int32_t));
  const uint32_t nn = SortUniqU32(need, nk);
  uint32_t* kslot = S_CAST(uint32_t*, XMalloc(nk * sizeof(int32_t)));
  for (uint32_t j = 0; j < nk; ++j) kslot[j] = LowerBoundU32(need, nn, kcols1.p[j]);
  char** ks = S_CAST(char**, XMalloc(nn * sizeof(intptr_t)));
  char** ke = S_CAST(char**, XMalloc(nn * sizeof(intptr_t)));
  while (!TextNextLine(&in.txs, &line)) {
    char* lnext = TextLineEnd(&in.txs);
    char* lend = lnext - 1;
    uint32_t kidx = UINT32_MAX;
    uint32_t found = 0;
    char* it = FieldStart(line, lend, 0, delim);
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
        kidx = StrSetFind(&keys, ks[0], klen, HashStr(ks[0], klen));
      } else {
        keybuf.n = 0;
        for (uint32_t j = 0; j < nk; ++j) {
          if (j) CvPushCh(&keybuf, '\0');
          CvAppend(&keybuf, ks[kslot[j]], ke[kslot[j]] - ks[kslot[j]]);
        }
        kidx = StrSetFind(&keys, keybuf.p, keybuf.n, HashStr(keybuf.p, keybuf.n));
      }
    }
    if (kidx == UINT32_MAX) {
      if (o.left_outer) {
        Emit(line, lend - line);
        Emit(fill.p, fill.n);
        EmitCh('\n');
      }
      continue;
    }
    for (uint32_t row = head.p[kidx]; row != UINT32_MAX; row = next.p[row]) {
      Emit(line, lend - line);
      Emit(pay.p[row], paylen.p[row]);
      EmitCh('\n');
    }
  }
  FinishTs(&in.txs);
  free(toks1.p); free(toks2.p);
  free(fs.p); free(fe.p);
  free(header2_payload.p);
  free(kcols2.p); free(is_key2);
  FreeStrSet(&keys);
  free(head.p); free(tail.p); free(next.p);
  free(pay.p); free(paylen.p);
  FreeArena(&arena);
  free(keybuf.p); free(paybuf.p); free(fill.p);
  free(kcols1.p); free(need); free(kslot); free(ks); free(ke);
  free(in.cols.p);
}

// ---- sort ----
//
// Every key is parsed once into a uint64 whose unsigned order is the key
// order (reversal and missing-last already applied).  Where 64 bits can't
// decide equality (strings longer than 8 bytes, nonstandard contig names), a
// per-key flag bit sends ties to a full comparison.  Remaining ties are
// broken by input order, so the sort is stable.

enum { kSkStr, kSkNum, kSkChr };

typedef struct SortKeyStruct {
  ColTok tok;
  uint32_t col;
  uint8_t type;  // kSk...
  uint8_t rev;
} SortKey;

CONSTI32(kMaxSortKeys, 32);

// Keys are encoded as a sequence of uint64 words (1 per numeric/chromosome
// key, 2 per string key: a 16-byte big-endian prefix).  The first
// kSortInlineWords live in the record itself; any others go in a side array.
CONSTI32(kSortInlineWords, 3);

typedef struct SortRecStruct {
  uint64_t w[kSortInlineWords];
  uint32_t idx;     // input line index
  uint32_t flags;   // bit j: key j needs a full comparison on equality
} SortRec;

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

typedef struct SortCtxStruct {
  SortKey keys[kMaxSortKeys];
  uint32_t key_ct;
  uint32_t word_start[kMaxSortKeys];
  uint32_t word_ct[kMaxSortKeys];
  uint32_t total_words;
  uint32_t extra_stride;               // words beyond the inline ones
  uint64_t* extra;                     // n * extra_stride
  const char** fptr[kMaxSortKeys];     // string/chr keys
  uint32_t* flen[kMaxSortKeys];
} SortCtx;

static inline uint64_t SortWord(const SortCtx* c, const SortRec* r, uint32_t w) {
  return (w < kSortInlineWords) ? r->w[w] : c->extra[S_CAST(uintptr_t, r->idx) * c->extra_stride + (w - kSortInlineWords)];
}

// Full comparison of key j for two lines whose encoded keys are equal.
static inline int SortFullCmp(const SortCtx* c, uint32_t j, uint32_t a, uint32_t b) {
  const SortKey* k = &(c->keys[j]);
  if (k->type == kSkNum) return 0;
  const char* sa = c->fptr[j][a]; const uint32_t la = c->flen[j][a];
  const char* sb = c->fptr[j][b]; const uint32_t lb = c->flen[j][b];
  const int r = StrCmp3(sa, la, sb, lb);
  return k->rev ? -r : r;
}

static inline uint32_t SortLess(const SortCtx* c, const SortRec* x, const SortRec* y) {
  // Common case: decided by the first inline word.
  if (x->w[0] != y->w[0]) return x->w[0] < y->w[0];
  const uint32_t fl = x->flags | y->flags;
  const uint32_t nk = c->key_ct;
  uint32_t w = 1;
  for (uint32_t j = 0; j < nk; ++j) {
    for (const uint32_t wend = c->word_start[j] + c->word_ct[j]; w < wend; ++w) {
      const uint64_t a = SortWord(c, x, w), b = SortWord(c, y, w);
      if (a != b) return a < b;
    }
    if (fl & (1U << j)) {
      const int r = SortFullCmp(c, j, x->idx, y->idx);
      if (r) return r < 0;
    }
  }
  return x->idx < y->idx;
}

#ifndef __cplusplus
// qsort() has no context argument.
static const SortCtx* g_sort_ctx;

static int SortRecCmp(const void* aa, const void* bb) {
  const SortRec* x = S_CAST(const SortRec*, aa);
  const SortRec* y = S_CAST(const SortRec*, bb);
  if (SortLess(g_sort_ctx, x, y)) return -1;
  return SortLess(g_sort_ctx, y, x);
}
#endif

static void ParseSortKeys(char* spec, SortCtx* ctx) {
  uint32_t key_ct = 0;
  for (char* p = spec; ; ) {
    char* e = strchr(p, ',');
    if (e) *e = '\0';
    // Keys past the limit are still parsed (into a scratch slot), so that
    // syntax errors anywhere in the list take precedence.
    SortKey scratch;
    SortKey* k = (key_ct < kMaxSortKeys) ? &(ctx->keys[key_ct]) : &scratch;
    ++key_ct;
    k->type = kSkStr; k->rev = 0; k->col = 0;
    char* colon = strchr(p, ':');
    if (colon) {
      *colon = '\0';
      for (const char* t = colon + 1; *t; ++t) {
        if (*t == 's') k->type = kSkStr;
        else if (*t == 'n' || *t == 'g') k->type = kSkNum;
        else if (*t == 'c') k->type = kSkChr;
        else if (*t == 'r') k->rev = 1;
        else Die("sort key types are s (string), n (numeric), c (chromosome), plus optional r (reverse)");
      }
    }
    ColTokVec toks = {nullptr, 0, 0};
    ParseCols(p, &toks);
    if (toks.n != 1 || (!toks.p[0].name && toks.p[0].a != toks.p[0].b)) {
      Die("each sort key must be a single column");
    }
    k->tok = toks.p[0];
    free(toks.p);
    if (!e) break;
    p = e + 1;
  }
  if (key_ct > kMaxSortKeys) Die("at most 32 sort keys are supported");
  ctx->key_ct = key_ct;
}

// Shared state for the parallel phases of ModeSort.
typedef struct SortJobStruct {
  SortCtx* ctx;
  char** lines;
  uint32_t* lens;
  const uint32_t* need;
  const uint32_t* key_slot;
  uint32_t nn;
  char delim;
  SortRec* recs;
  // Output gathering.
  uint32_t thread_ct;
  char* batch;
  uintptr_t* seg_bytes;
  uintptr_t* seg_off;
} SortJob;

// Stores encoded key word w of a record.
static inline void SortPut(uint32_t w, uint64_t v, SortRec* r, uint64_t* extra) {
  if (w < kSortInlineWords) r->w[w] = v;
  else extra[w - kSortInlineWords] = v;
}

static void SortExtractKeys(void* shared, uintptr_t lo, uintptr_t hi, __maybe_unused uint32_t tidx) {
  SortJob* job = S_CAST(SortJob*, shared);
  SortCtx* ctx = job->ctx;
  const uint32_t nn = job->nn;
  const uint32_t nk = ctx->key_ct;
  const uint32_t* need = job->need;
  const uint32_t* key_slot = job->key_slot;
  const char delim = job->delim;
  char** fs = S_CAST(char**, XMalloc(nn * sizeof(intptr_t)));
  char** fe = S_CAST(char**, XMalloc(nn * sizeof(intptr_t)));
  static char empty[1] = {'\0'};
  for (uintptr_t li = lo; li < hi; ++li) {
    char* ln = job->lines[li];
    char* lend = ln + job->lens[li] - 1;  // the '\n'
    // Locate needed fields left to right; absent ones are empty.
    uint32_t found = 0;
    char* it = FieldStart(ln, lend, 0, delim);
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
        fs[found] = it;
        fe[found] = FieldEnd(it, lend, delim);
      }
    }
    for (uint32_t s = found; s < nn; ++s) { fs[s] = empty; fe[s] = empty; }
    SortRec* r = &(job->recs[li]);
    r->idx = li; r->flags = 0;
    for (uint32_t w = 0; w < kSortInlineWords; ++w) r->w[w] = 0;
    uint64_t* extra = ctx->extra_stride ? &(ctx->extra[li * ctx->extra_stride]) : nullptr;
    for (uint32_t j = 0; j < nk; ++j) {
      const SortKey* k = &(ctx->keys[j]);
      const char* f = fs[key_slot[j]];
      const uint32_t flen = fe[key_slot[j]] - fs[key_slot[j]];
      const uint32_t w0 = ctx->word_start[j];
      if (k->type == kSkNum) {
        double d;
        uint64_t v = kSortMissing;
        if (ParseNumExact(f, f + flen, &d)) {
          v = EncodeDouble(d);
          if (k->rev) v = ~v;
        }
        SortPut(w0, v, r, extra);
        continue;
      }
      ctx->fptr[j][li] = f; ctx->flen[j][li] = flen;
      if (k->type == kSkStr) {
        uint64_t pre_hi, pre_lo;
        StrPrefix16(f, flen, &pre_hi, &pre_lo);
        if (flen > 16) r->flags |= 1U << j;
        SortPut(w0, k->rev ? ~pre_hi : pre_hi, r, extra);
        SortPut(w0 + 1, k->rev ? ~pre_lo : pre_lo, r, extra);
      } else {
        const uint32_t rank = ChrRank(f, flen);
        if (rank == UINT32_MAX) r->flags |= 1U << j;
        SortPut(w0, k->rev ? ~S_CAST(uint64_t, rank) : rank, r, extra);
      }
    }
  }
  free(fs);
  free(fe);
}

static void SortRange(void* shared, uintptr_t lo, uintptr_t hi, __maybe_unused uint32_t tidx) {
  SortJob* job = S_CAST(SortJob*, shared);
  SortRec* recs = job->recs;
#ifdef __cplusplus
  const SortCtx* ctxp = job->ctx;
  std::sort(&(recs[lo]), &(recs[hi]), [ctxp](const SortRec& x, const SortRec& y) { return SortLess(ctxp, &x, &y) != 0; });
#else
  qsort(&(recs[lo]), hi - lo, sizeof(SortRec), SortRecCmp);
#endif
}

// One slice of a two-run merge: output positions [lo, hi).
typedef struct MergeJobStruct {
  const SortCtx* ctx;
  const SortRec* a;
  uintptr_t na;
  const SortRec* b;
  uintptr_t nb;
  SortRec* out;
} MergeJob;

// Number of elements the first k merged outputs take from run a ("merge
// path"); unique since the comparator is a strict total order.
static uintptr_t MergeCorank(const MergeJob* m, uintptr_t k) {
  uintptr_t lo = (k > m->nb) ? k - m->nb : 0;
  uintptr_t hi = (k < m->na) ? k : m->na;
  while (lo < hi) {
    const uintptr_t mid = (lo + hi) / 2;
    if (SortLess(m->ctx, &(m->a[mid]), &(m->b[k - mid - 1]))) lo = mid + 1;
    else hi = mid;
  }
  return lo;
}

static void MergeSlice(void* shared, uintptr_t k0, uintptr_t k1, __maybe_unused uint32_t tidx) {
  const MergeJob* m = S_CAST(const MergeJob*, shared);
  const SortCtx* ctx = m->ctx;
  const uintptr_t i0 = MergeCorank(m, k0), i1 = MergeCorank(m, k1);
  const SortRec* a = &(m->a[i0]);
  const SortRec* a_end = &(m->a[i1]);
  const SortRec* b = &(m->b[k0 - i0]);
  const SortRec* b_end = &(m->b[k1 - i1]);
  SortRec* out = &(m->out[k0]);
  // Ties go to a, as with std::merge.
  while ((a != a_end) && (b != b_end)) {
    if (SortLess(ctx, b, a)) *out++ = *b++;
    else *out++ = *a++;
  }
  memcpy(out, a, (a_end - a) * sizeof(SortRec));
  out += a_end - a;
  memcpy(out, b, (b_end - b) * sizeof(SortRec));
}

static void SortCountBytes(void* shared, uintptr_t s0, uintptr_t s1, uint32_t tidx) {
  SortJob* job = S_CAST(SortJob*, shared);
  uintptr_t bytes = 0;
  for (uintptr_t q = s0; q < s1; ++q) bytes += job->lens[job->recs[q].idx];
  job->seg_bytes[tidx] = bytes;
}

static void SortGather(void* shared, uintptr_t s0, uintptr_t s1, uint32_t tidx) {
  SortJob* job = S_CAST(SortJob*, shared);
  char* w = &(job->batch[job->seg_off[tidx]]);
  for (uintptr_t q = s0; q < s1; ++q) {
    const uint32_t idx = job->recs[q].idx;
    memcpy(w, job->lines[idx], job->lens[idx]);
    w += job->lens[idx];
  }
}

static double MonotonicSec() {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return S_CAST(double, ts.tv_sec) + S_CAST(double, ts.tv_nsec) * 1e-9;
}

// Phase timing on stderr when PTXT_TIMING is set.
static void SortPhase(uint32_t timing, const char* name, double* t_last) {
  if (!timing) return;
  const double now = MonotonicSec();
  fprintf(stderr, "%-8s %.3fs\n", name, now - *t_last);
  *t_last = now;
}

static void ModeSort(int argc, char** argv) {
  Opts o; InitOpts(&o);
  int i = ParseOpts(argc, argv, 2, &o, 0);
  if (argc - i != 2) DieUsage("ptxt sort [options] KEY[:TYPE][,KEY[:TYPE]...] FILE");
  SortCtx ctx;
  memset(&ctx, 0, sizeof(ctx));
  ParseSortKeys(argv[i], &ctx);
  const uint32_t nk = ctx.key_ct;
  for (uint32_t j = 0; j < nk; ++j) if (ctx.keys[j].tok.name) o.header = 1;
  Input in; OpenInput(argv[i + 1], nullptr, &o, &in);
  const char delim = o.delim;
  for (uint32_t j = 0; j < nk; ++j) {
    const ColTok* t = &(ctx.keys[j].tok);
    ctx.keys[j].col = t->name ? ResolveNamedCol(t->name, in.hline, in.hlend, delim) : t->a - 1;
  }
  // The header points into the stream's buffer, which loading overwrites.
  char* header = nullptr;
  uintptr_t header_len = 0;
  if (in.hline) {
    header_len = in.hlend + 1 - in.hline;
    header = S_CAST(char*, XMalloc(header_len));
    memcpy(header, in.hline, header_len);
  }

  const uint32_t timing = getenv("PTXT_TIMING") != nullptr;
  double t_last = timing ? MonotonicSec() : 0;
  // Load all lines (including their '\n') into chunked storage.
  Arena chunks; InitArena(64 << 20, &chunks);
  CharpVec lines = {nullptr, 0, 0};
  U32Vec lens = {nullptr, 0, 0};
  char* line;
  while (!TextNextLine(&in.txs, &line)) {
    const uintptr_t len = TextLineEnd(&in.txs) - line;
    if (len > 0xffffffffU) Die("line too long");
    char* dst = ArenaAlloc(len, &chunks);
    memcpy(dst, line, len);
    VEC_PUSH(lines, dst);
    VEC_PUSH(lens, len);
  }
  FinishTs(&in.txs);
  const uintptr_t n = lines.n;
  if (n > 0xffffffffU) Die("too many lines");

  SortPhase(timing, "load", &t_last);
  // Extract keys in parallel.
  uint32_t thread_ct = NumCpu(nullptr);
  if (!thread_ct) thread_ct = 1;
  if (thread_ct > 16) thread_ct = 16;
  uint32_t need[kMaxSortKeys];
  for (uint32_t j = 0; j < nk; ++j) need[j] = ctx.keys[j].col;
  const uint32_t nn = SortUniqU32(need, nk);
  uint32_t key_slot[kMaxSortKeys];
  for (uint32_t j = 0; j < nk; ++j) key_slot[j] = LowerBoundU32(need, nn, ctx.keys[j].col);
  for (uint32_t j = 0; j < nk; ++j) {
    ctx.word_start[j] = ctx.total_words;
    const uint32_t wc = (ctx.keys[j].type == kSkStr) ? 2 : 1;
    ctx.word_ct[j] = wc;
    ctx.total_words += wc;
    if (ctx.keys[j].type != kSkNum) {
      ctx.fptr[j] = S_CAST(const char**, XMalloc(n * sizeof(intptr_t)));
      ctx.flen[j] = S_CAST(uint32_t*, XMalloc(n * sizeof(int32_t)));
    }
  }
  if (ctx.total_words > kSortInlineWords) {
    ctx.extra_stride = ctx.total_words - kSortInlineWords;
    ctx.extra = S_CAST(uint64_t*, XMalloc(n * ctx.extra_stride * sizeof(int64_t)));
  }
  SortRec* recs = S_CAST(SortRec*, XMalloc(n * sizeof(SortRec)));
  SortJob job;
  job.ctx = &ctx;
  job.lines = lines.p;
  job.lens = lens.p;
  job.need = need;
  job.key_slot = key_slot;
  job.nn = nn;
  job.delim = delim;
  job.recs = recs;
  job.thread_ct = thread_ct;
  ParallelFor(thread_ct, n, SortExtractKeys, &job);

  SortPhase(timing, "keys", &t_last);
  // Parallel sort: sort thread_ct chunks, then merge pairwise.
#ifndef __cplusplus
  g_sort_ctx = &ctx;
#endif
  const uint32_t part_ct = (n >= 65536) ? thread_ct : 1;
  ParallelFor(thread_ct, n, SortRange, &job);
  SortPhase(timing, "sort", &t_last);
  if (part_ct > 1) {
    SortRec* src = recs;
    SortRec* dst = S_CAST(SortRec*, XMalloc(n * sizeof(SortRec)));
    uintptr_t bounds[17];
    uintptr_t bound_ct = part_ct + 1;
    for (uint32_t t = 0; t <= part_ct; ++t) bounds[t] = n * t / part_ct;
    MergeJob mjobs[8];
    RangeTask* tasks = S_CAST(RangeTask*, XMalloc(8 * thread_ct * sizeof(RangeTask)));
    // Each round merges adjacent pairs of runs.  Every merge is split across
    // all threads: output position k of a merge of runs A and B takes
    // exactly i elements from A, where i is found by binary search ("merge
    // path").
    while (bound_ct > 2) {
      uintptr_t nb[17];
      uintptr_t nb_ct = 0;
      uint32_t task_ct = 0;
      for (uintptr_t b = 0; b + 1 < bound_ct; b += 2) {
        nb[nb_ct++] = bounds[b];
        if (b + 2 == bound_ct) {
          memcpy(&(dst[bounds[b]]), &(src[bounds[b]]), (bounds[b + 1] - bounds[b]) * sizeof(SortRec));
          continue;
        }
        MergeJob* m = &(mjobs[b / 2]);
        m->ctx = &ctx;
        m->a = &(src[bounds[b]]);
        m->na = bounds[b + 1] - bounds[b];
        m->b = &(src[bounds[b + 1]]);
        m->nb = bounds[b + 2] - bounds[b + 1];
        m->out = &(dst[bounds[b]]);
        const uintptr_t total = m->na + m->nb;
        for (uint32_t t = 0; t < thread_ct; ++t) {
          RangeTask* task = &(tasks[task_ct++]);
          task->fn = MergeSlice; task->shared = m; task->tidx = t;
          task->lo = total * t / thread_ct;
          task->hi = total * (t + 1) / thread_ct;
        }
      }
      nb[nb_ct++] = n;
      RunTasks(tasks, task_ct);
      memcpy(bounds, nb, nb_ct * sizeof(intptr_t));
      bound_ct = nb_ct;
      SortRec* tmp = src; src = dst; dst = tmp;
    }
    free(tasks);
    free(dst);
    recs = src;
    job.recs = recs;
  }

  SortPhase(timing, "merge", &t_last);
  Emit(header, header_len);
  Flush();
  // Gather lines in sorted order, in batches: each batch's per-thread
  // segments are copied in parallel into one buffer, then written at once.
  {
    const uintptr_t kBatchRecs = 1 << 20;
    uintptr_t seg_bytes[16];
    uintptr_t seg_off[17];
    job.seg_bytes = seg_bytes;
    job.seg_off = seg_off;
    job.batch = nullptr;
    uintptr_t batch_cap = 0;
    for (uintptr_t lo = 0; lo < n; lo += kBatchRecs) {
      const uintptr_t hi = (n < lo + kBatchRecs) ? n : lo + kBatchRecs;
      // Byte offsets of each thread's segment.
      SpawnRanges(thread_ct, lo, hi, SortCountBytes, &job);
      seg_off[0] = 0;
      for (uint32_t t = 0; t < thread_ct; ++t) seg_off[t + 1] = seg_off[t] + seg_bytes[t];
      if (batch_cap < seg_off[thread_ct]) {
        free(job.batch);
        batch_cap = seg_off[thread_ct];
        job.batch = S_CAST(char*, XMalloc(batch_cap));
      }
      SpawnRanges(thread_ct, lo, hi, SortGather, &job);
      if (fwrite(job.batch, 1, seg_off[thread_ct], stdout) != seg_off[thread_ct]) Die("write failure");
    }
    free(job.batch);
  }
  SortPhase(timing, "output", &t_last);
  free(recs);
  free(ctx.extra);
  for (uint32_t j = 0; j < nk; ++j) {
    free(ctx.fptr[j]);
    free(ctx.flen[j]);
  }
  free(lines.p);
  free(lens.p);
  FreeArena(&chunks);
  free(header);
  free(in.cols.p);
}

// ---- uniq ----
//
// Hash-based, so input need not be sorted; output is in first-seen order.
// On sorted input this is identical to uniq (with -r for -d), and
// without flags it's awk '!seen[KEY]++'.

static void ModeUniq(int argc, char** argv) {
  Opts o; InitOpts(&o);
  uint32_t count = 0, dups = 0, singles = 0;
  int i = 2;
  while (1) {
    i = ParseOpts(argc, argv, i, &o, 0);
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
  const uint32_t* cols = in.cols.p;
  const uint32_t nk = in.cols.n;
  const char delim = o.delim;
  const uint32_t deferred = count || dups || singles;
  if (in.hline) {
    if (count) Emit("count\t", 6);
    EmitHeader(&in);
  }
  StrSet keys; InitStrSet(&keys);
  U64Vec counts = {nullptr, 0, 0};  // only if deferred
  // First line of each key, if deferred and the key isn't the whole line.
  Arena lines; InitArena(16 << 20, &lines);
  CcpVec lptr = {nullptr, 0, 0};
  U32Vec llen = {nullptr, 0, 0};
  CharVec kbuf = {nullptr, 0, 0};
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
      kbuf.n = 0;
      for (uint32_t j = 0; j < nk; ++j) {
        if (j) CvPushCh(&kbuf, '\0');
        char* f = FieldStart(line, lend, cols[j], delim);
        if (f) CvAppend(&kbuf, f, FieldEnd(f, lend, delim) - f);
      }
      k = kbuf.p; klen = kbuf.n;
    }
    uint32_t idx;
    if (prev && klen == prev_len && !memcmp(k, prev, klen)) {
      idx = prev_idx;
    } else {
      uint32_t added = 1;
      if (ascending && prev && StrCmp3(k, klen, prev, prev_len) < 0) {
        ascending = 0;
        StrSetBuildIndex(&keys);
      }
      if (ascending) idx = StrSetAppendUnindexed(k, klen, &keys);
      else idx = StrSetInsert(k, klen, HashStr(k, klen), &keys, &added);
      prev = StrSetPtr(&keys, idx); prev_len = klen; prev_idx = idx;
      if (added) {
        if (!deferred) {
          Emit(line, lnext - line);
          continue;
        }
        VEC_PUSH(counts, 0);
        if (!whole) {
          const uint32_t len = lend - line;
          char* c = ArenaAlloc(len, &lines);
          memcpy(c, line, len);
          VEC_PUSH(lptr, c); VEC_PUSH(llen, len);
        }
      }
    }
    if (deferred) ++counts.p[idx];
  }
  FinishTs(&in.txs);
  for (uint32_t j = 0; deferred && (j < StrSetSize(&keys)); ++j) {
    const uint64_t ct = counts.p[j];
    if ((dups && ct < 2) || (singles && ct > 1)) continue;
    if (count) {
      char numbuf[24];
      const int nlen = snprintf(numbuf, sizeof numbuf, "%llu\t", S_CAST(unsigned long long, ct));
      Emit(numbuf, nlen);
    }
    if (whole) Emit(StrSetPtr(&keys, j), StrSetLen(&keys, j));
    else Emit(lptr.p[j], llen.p[j]);
    EmitCh('\n');
  }
  FreeStrSet(&keys);
  free(counts.p);
  FreeArena(&lines);
  free(lptr.p); free(llen.p); free(kbuf.p);
  free(in.cols.p);
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
  g_outbuf = S_CAST(char*, XMalloc(kOutBufSize)); g_outw = g_outbuf;
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
