# ptxt

`ptxt` is a prototype of fast replacements for the Unix text tools most often
used in bioinformatics scripts (`cut`, `paste`, `awk` one-liners, `sort`,
`join`, `uniq`, `tr -s`). It is built on the plink2 text-reading library.

- Input can be uncompressed, gzip, BGZF or zstd. BGZF is decompressed on
  several threads.
- Lines must end in LF.
- Use `/dev/stdin` as FILE to read from a pipe.
- Each subcommand aims to be at least 2× faster than every common macOS
  implementation (BSD tools, GNU coreutils, gawk, mawk) on typical inputs:
  summary statistics, VCFs and PLINK reports. See [Performance](#performance).

## Building and testing

```
make          # builds ptxt and gen (the benchmark/test data generator)
make test     # runs tests/run_tests.sh
```

The tests check output byte for byte against `cut`, `paste`, POSIX `awk`,
`sed`, `sort -s` and `uniq`. Where ptxt deliberately differs from those tools,
the tests use literal expected output instead. The tests pass with BSD or GNU
tools and with any of BSD awk, gawk or mawk. BGZF and zstd tests are skipped
when `bgzip` or `zstd` is missing.

`gen ss NLINE` writes a summary-statistics TSV with a header line.
`gen vcf NSAMPLE NVAR` writes a VCF-like file with GT-only genotype columns.
The output is deterministic.

## Common options and column syntax

These options work with select, filter, keep, remove, count, prefix,
unprefix, sum, sort, join and uniq:

| Option | Meaning |
|---|---|
| `-w` | Whitespace-delimited input (awk's default field splitting). |
| `-d C` | Input delimiter character (default: tab). |
| `-O STR` | Output delimiter for select and join (default: tab). |
| `-H` | The first line is a header (see below). |
| `--` | End of options; only needed if FILE or EXPR is literally `-w` etc. |

What each subcommand does with a `-H` header line:
- filter, keep, remove, prefix and unprefix copy it through.
- sort and uniq keep it first.
- count and sum skip it.
- select treats it like any other line.

Columns can be:
- 1-based numbers, e.g. `3`
- ranges, e.g. `3-5` (select only)
- comma-separated lists, e.g. `CHR,POS`
- header names, e.g. `P`. Using a name implies `-H`.

Numbers in fields are parsed with correct rounding (`ffc.h`), so comparisons
against boundary values such as `0.9900 <= 0.99` give the same answers as
awk.

## Subcommands

### cut

```
ptxt cut LIST FILE
```

Same as `cut -f LIST` on tab-delimited input. LIST must be ascending, e.g.
`1-8` or `2,11`. Contiguous ranges are copied as one block.

### paste

```
ptxt paste FILE1 FILE2
```

Same as `paste FILE1 FILE2`.

### select

```
ptxt select [opts] COLS FILE
```

Same as `awk '{print $3, $1, ...}'`. Columns can be in any order and can be
repeated. Output fields are separated by `-O` (default: tab).

### filter

```
ptxt filter [opts] EXPR FILE
```

Prints the lines that satisfy an awk-like expression. For example:

```
ptxt filter '$P < 5e-8 && ($1 == "22" || $CHR == "X") && abs($BETA/$SE) > 5' ss.tsv
ptxt filter '-log10($P) > 8' ss.tsv.gz
```

- Fields are written `$N`, `$NAME` or `$"NAME"`.
- Operators, in awk precedence:
  - `||` (also `|`, `or`)
  - `&&` (also `&`, `and`)
  - comparisons: `<`, `<=`, `>`, `>=`, `==` (also `=`), `!=`
  - `+`, `-`
  - `*`, `/`, `%`
  - unary `!` (also `not`), `-`, `+`
- Parentheses and the functions `abs`, `log`, `log10`, `exp` and `sqrt` are
  supported.
- All arguments before FILE are joined with spaces, so the expression doesn't
  have to be quoted as a single argument.

How values compare follows awk, except for missing values:
- A field counts as a number if the whole field parses as one.
- Field vs `"string"` compares as strings.
- Field vs field compares numerically if both are numbers, as strings
  otherwise.
- Field vs number:
  - If the field is numeric, they compare numerically.
  - Otherwise the comparison is false, and `!=` is true. This is plink2's
    `--extract-if-info` convention; awk would fall back to a string
    comparison.
- Arithmetic on a non-numeric field gives a missing value. Any comparison
  against it is false (`!=` is true), and it is false as a condition.
- A bare value used as a condition is true when it is nonzero or a nonempty
  string.
- A field past the end of the line is an empty string.

### keep / remove

```
ptxt keep [opts] COL IDFILE FILE
ptxt remove [opts] COL IDFILE FILE
```

`keep` prints the lines whose COL value is listed in IDFILE. `remove` prints
the lines whose COL value isn't. IDs are the first whitespace-delimited token
on each IDFILE line. This is the hash-lookup form of `grep -F -w -f` or
`awk 'NR==FNR{a[$1];next} $3 in a'`, but limited to a single column.

### count

```
ptxt count [opts] COL FILE
```

Prints `value<TAB>count` for each distinct value in COL, in the order the
values first appear. This replaces `cut -f COL | sort | uniq -c` and
`awk '{c[$1]++} END {...}'`, and the input doesn't need to be sorted.

### prefix / unprefix

```
ptxt prefix [opts] COL STR FILE
ptxt unprefix [opts] COL STR FILE
```

`prefix` adds STR to the start of COL. `unprefix` removes STR from the start
of COL where it is present. The rest of each line is copied unchanged.
Example: `ptxt prefix 1 chr ss.tsv`.

### tsv

```
ptxt tsv FILE
```

Converts whitespace-delimited lines, such as PLINK 1 reports, to TSV. This is
what `awk '{$1=$1}1' OFS='\t'` does: leading and trailing whitespace is
dropped and runs of spaces and tabs become single tabs. It is also what
`tr -s ' ' '\t'` is usually used for, but `tr` leaves a leading tab when
lines are padded.

### sum

```
ptxt sum [opts] COL FILE
```

Prints the sum of the numeric values in COL.

### sort

```
ptxt sort [opts] KEY[:TYPE][,KEY[:TYPE]...] FILE
```

A stable sort on one or more keys: lines with equal keys stay in input order,
as with `sort -s`. TYPE is one of:

- `s`: string, in byte order (the default)
- `n`: numeric, in any decimal or scientific notation. `g` is an alias.
  Values that aren't numbers or are missing always sort last, in either
  direction.
- `c`: chromosome. The order is 1–22 and up, X, Y, XY, MT (with or without a
  `chr` prefix), then other contigs by name.

Add `r` after a TYPE to reverse that key. Examples:

```
ptxt sort CHR:c,POS:n ss.tsv
ptxt sort P:n ss.tsv.gz
ptxt sort -H BETA:nr ss.tsv
```

The whole input is loaded into memory, and the sort uses all cores, up to 16.

### join

```
ptxt join [opts] KEYS1 FILE1 KEYS2 FILE2
```

A hash join:
- For each FILE1 line, ptxt prints that line once for every FILE2 row whose
  KEYS2 columns exactly match its KEYS1 columns, followed by the row's
  non-key columns.
- Output keeps FILE1's order. When several FILE2 rows match, they come out in
  FILE2 order.
- KEYS1 and KEYS2 are comma-separated lists of the same length, so composite
  keys work. They may name the columns in a different order in each file.
- FILE2 is held in memory and FILE1 is streamed, so neither file needs to be
  sorted.

Options:

| Option | Meaning |
|---|---|
| `-a` | Also print FILE1 lines with no match, padded with fill values. |
| `-e STR` | Fill value for `-a` (default: `NA`). |

- With `-H` or column names, both files have headers. The output header is
  FILE1's header followed by FILE2's non-key column names.
- FILE2 rows that are short are padded with empty fields.
- FILE2 rows that lack a key column are skipped, and so are blank FILE2
  lines.

```
ptxt join -a SNP sumstats.tsv ID annot.tsv
ptxt join CHR,POS a.tsv chrom,pos b.tsv
```

### uniq

```
ptxt uniq [opts] [-c] [-r|-u] [COLS] FILE
```

Prints the first line for each distinct key, in the order the keys first
appear. The key is the whole line, or the COLS fields if COLS is given. A
missing field counts as an empty key, as in awk.

The input doesn't need to be sorted. With no flags this is
`awk '!seen[$3]++'`. On sorted input the output is the same as `uniq`'s.

| Option | Meaning |
|---|---|
| `-c` | Put the key's count and a tab before each line. |
| `-r` | Print only keys that occur more than once (`uniq -d`; `-d` is ptxt's delimiter option). |
| `-u` | Print only keys that occur exactly once. |

- `-c`, `-r` and `-u` print their output at the end and keep each key's first
  line in memory.
- With `-H` the header is copied through, with `count<TAB>` added in front
  when `-c` is used.
- While keys arrive in ascending order, ptxt doesn't use a hash table, which
  makes it fast on sorted input. If a key arrives out of order, ptxt builds
  the table at that point.

## Performance

Measurements were taken on an M2 Pro with `LC_ALL=C` and a warm cache, with
output written to a file. Inputs were about 10M lines (0.8–1 GB). Each
speedup is against the fastest alternative tried.

| Task | Fastest alternative | Speedup |
|---|---|---|
| filter, compound expression | mawk | 4–7.5× |
| keep 10k / 1M IDs | ripgrep -F -w / awk | >5× |
| sort chr,pos / P-value / rsID | GNU sort (best of default and --parallel) | 4.7× / 52× / 3.0× |
| join, 1M-row lookup table | gawk; GNU join on pre-sorted input | 7.9×; 4.4× |
| join -a | GNU join -a1 on pre-sorted input | 9.3× |
| uniq by column, unsorted (`!s[$3]++`) | gawk | 9× |
| uniq whole line, unsorted | mawk | 12× |
| uniq on sorted full lines (plain/-c/-d) | GNU uniq | 2.2–3.5× |
| tsv (vs `tr -s ' ' '\t'`) | GNU tr | 2.75× |

Two cases come in under 2×: `cut -f1-8` on an uncompressed wide VCF (1.3×
against GNU cut), and string-equality filters (1.1–1.7× against ripgrep).
