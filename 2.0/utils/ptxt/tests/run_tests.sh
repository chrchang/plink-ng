#!/bin/bash

# ptxt, against the standard tools it replaces (cut, paste, awk, sed) and, where
# ptxt deliberately differs from them, against literal expected output.
#
# Usage: ./run_tests.sh <directory containing the ptxt and gen binaries>
#
# All oracle awk programs stick to POSIX awk, so any of mawk/gawk/BSD awk
# works.  bgzip and zstd tests are skipped when those tools are absent.

set -eo pipefail
export LC_ALL=C

BINDIR=$(cd "$1" && pwd)
P=$BINDIR/ptxt
GEN=$BINDIR/gen

WORKDIR=$(mktemp -d)
trap 'rm -rf "${WORKDIR:?}"' EXIT
cd "$WORKDIR"

PASS=0

# check NAME ORACLE_CMD PTXT_CMD: outputs must be byte-identical.
check() {
  if ! cmp -s <(bash -c "$2") <(bash -c "$3"); then
    echo "FAIL: $1"
    echo "  oracle: $2"
    echo "  ptxt:   $3"
    diff <(bash -c "$2") <(bash -c "$3") | head -10
    exit 1
  fi
  PASS=$((PASS + 1))
}

# expect NAME PTXT_CMD EXPECTED: output must equal EXPECTED (printf-escaped).
expect() {
  if ! cmp -s <(printf "$3") <(bash -c "$2"); then
    echo "FAIL: $1"
    echo "  cmd: $2"
    echo "  expected:"; printf "$3" | od -c | head -5
    echo "  got:"; bash -c "$2" | od -c | head -5
    exit 1
  fi
  PASS=$((PASS + 1))
}

# fails NAME PTXT_CMD MESSAGE_SUBSTRING: must exit nonzero, mentioning MESSAGE.
fails() {
  local out
  if out=$(bash -c "$2" 2>&1); then
    echo "FAIL: $1 (expected an error)"; echo "  cmd: $2"; exit 1
  fi
  if [[ "$out" != *"$3"* ]]; then
    echo "FAIL: $1 (wrong error)"; echo "  cmd: $2"
    echo "  expected substring: $3"; echo "  got: $out"; exit 1
  fi
  PASS=$((PASS + 1))
}

# ---- fixtures ----
"$GEN" ss 20000 > ss.tsv
"$GEN" vcf 300 400 > wide.vcf
# PLINK 1-style whitespace-padded report.
awk -F'\t' '{printf "%3s %12s %11s %3s %3s %8s %7s %9s %9s %8s %11s %6s\n", $1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12}' ss.tsv > ws.assoc
awk -F'\t' 'NR > 1 && NR % 7 == 0 {print $3}' ss.tsv > ids.txt
awk -F'\t' -v OFS='\t' 'NR > 1 {$1 = "chr" $1} 1' ss.tsv > sschr.tsv
printf 'a\tb\t\tdd\n\n  x  y   z \nsolo\nq\t\t\t\n' > edge.tsv
printf 'a\tb\nc\td' > nonl.tsv   # no trailing newline

# ---- cut ----
for L in 1 2 1-2 1,3 3-5 1-8 2,11 12 1,3-12 9-30 20; do
  check "cut -f$L ss" "cut -f$L ss.tsv" "$P cut $L ss.tsv"
done
for L in 1-8 1,3-309 10-160 300-309 1,5,9-12,309 308-320; do
  check "cut -f$L wide" "cut -f$L wide.vcf" "$P cut $L wide.vcf"
done
for L in 1 2 1-2 2-3 1,3 1-4 3-9 1,2,4 5; do
  check "cut -f$L edge" "cut -f$L edge.tsv" "$P cut $L edge.tsv"
done
expect "cut no trailing newline" "$P cut 2 nonl.tsv" 'b\nd\n'
expect "filter no trailing newline" "$P filter '\$1 == \"c\"' nonl.tsv" 'c\td\n'

# ---- paste ----
check "paste" "paste ss.tsv ws.assoc" "$P paste ss.tsv ws.assoc"
fails "paste short FILE2" "$P paste ss.tsv ids.txt" "fewer lines"
fails "paste long FILE2" "$P paste ids.txt ss.tsv" "more lines"

# ---- stdin as - ----
check "cut -" "cut -f2-3 ss.tsv" "cat ss.tsv | $P cut 2-3 -"
check "paste - FILE2" "paste ss.tsv ws.assoc" "cat ss.tsv | $P paste - ws.assoc"
check "paste FILE1 -" "paste ss.tsv ws.assoc" "cat ws.assoc | $P paste ss.tsv -"
check "filter -" "awk -F'\t' 'NR == 1 || \$11 < 0.5' ss.tsv" "cat ss.tsv | $P filter '\$P < 0.5' -"
check "filter - after --" "awk -F'\t' 'NR == 1 || \$11 < 0.5' ss.tsv" "cat ss.tsv | $P filter -H -- '\$P < 0.5' -"
check "count -" "$P count CHR ss.tsv" "cat ss.tsv | $P count CHR -"
check "tsv -" "awk -v OFS='\t' '{\$1 = \$1} 1' ws.assoc" "cat ws.assoc | $P tsv -"
check "uniq whole-line -" "awk '!s[\$0]++' ss.tsv" "cat ss.tsv | $P uniq -"
check "uniq -c -" "$P uniq -c 1 ss.tsv" "cat ss.tsv | $P uniq -c 1 -"
check "remove -" "awk -F'\t' 'NR == FNR {k[\$1]; next} !(\$3 in k)' ids.txt ss.tsv" "cat ss.tsv | $P remove 3 ids.txt -"
check "./- is a file" "cut -f1 ss.tsv" "cp ss.tsv ./- && $P cut 1 ./- </dev/null; rm -f ./-"
fails "paste - -" "$P paste - - </dev/null" "only one input file can be stdin"
fails "keep - /dev/stdin" "$P keep 1 - /dev/stdin </dev/null" "only one input file can be stdin"
fails "join - -" "$P join 1 - 1 - </dev/null" "only one input file can be stdin"

# ---- select ----
check "select reorder" "awk -F'\t' -v OFS='\t' '{print \$3, \$1, \$2, \$11}' ss.tsv" "$P select 3,1,2,11 ss.tsv"
check "select names" "awk -F'\t' -v OFS='\t' '{print \$3, \$1, \$11}' ss.tsv" "$P select SNP,CHR,P ss.tsv"
check "select repeat" "awk -F'\t' -v OFS='\t' '{print \$2, \$2, \$1}' ss.tsv" "$P select 2,2,1 ss.tsv"
check "select reversed range" "awk -F'\t' -v OFS='\t' '{print \$12, \$11, \$10}' ss.tsv" "$P select 12-10 ss.tsv"
check "select -w -O" "awk '{print \$2, \$11}' ws.assoc" "$P select -w -O ' ' 2,11 ws.assoc"
check "select -w names" "awk -v OFS='\t' '{print \$2, \$11}' ws.assoc" "$P select -w SNP,P ws.assoc"
check "select missing fields" "awk -F'\t' -v OFS='\t' '{print \$4, \$1, \$3}' edge.tsv" "$P select 4,1,3 edge.tsv"
check "select -w missing fields" "awk -v OFS='\t' '{print \$3, \$1}' edge.tsv" "$P select -w 3,1 edge.tsv"
check "select -d ," "awk -F'\t' -v OFS=',' '{print \$1, \$2, \$3}' ss.tsv | awk -F, -v OFS='\t' '{print \$3, \$1}'" "awk -F'\t' -v OFS=',' '{print \$1, \$2, \$3}' ss.tsv | $P select -d , 3,1 -"
fails "select unknown name" "$P select NOPE ss.tsv" "column 'NOPE' not found"

# ---- filter: agreement with awk ----
f() {  # f NAME AWK_COND PTXT_ARGS...: header passes through on both sides.
  local name=$1 cond=$2; shift 2
  check "filter $name" "awk -F'\t' 'NR == 1 || ($cond)' ss.tsv" "$P filter -H $* ss.tsv"
}
f "numeric lt" '$11 < 5e-8' "'\$11 < 5e-8'"
f "name" '$11 < 0.001' "'\$P < 0.001'"
f "literal on left" '$11 < 0.001' "'0.001 > \$P'"
f "le ge" '$6 >= 0.25 && $6 <= 0.75' "'\$A1_FREQ >= 0.25 && \$A1_FREQ <= 0.75'"
f "and chain" '$11 < 0.01 && $6 >= 0.01 && $6 <= 0.99' "'\$P < 0.01 && \$A1_FREQ >= 0.01 && \$A1_FREQ <= 0.99'"
f "or and parens" '($1 == "22" || $1 == "3") && $11 < 0.05' "'(\$CHR == \"22\" || \$CHR == \"3\") && \$P < 0.05'"
f "keywords" '!($4 == "A" || $4 == "T")' "'not (\$A1 == \"A\" or \$A1 == \"T\")'"
f "single-char ops" '$4 == "A" && $5 != "C"' "'\$A1 == \"A\" & \$A2 != \"C\"'"
f "field vs field str" '$4 != $5' "'\$A1 != \$A2'"
f "field vs field num" '$8 > $9' "'\$BETA > \$SE'"
f "string order" '$3 < "rs2"' "'\$SNP < \"rs2\"'"
f "numeric eq" '$7 == 400500' "'\$N == 400500'"
f "numeric ne" '$7 != 400500' "'\$N != 400500'"
f "numeric string vs number" '$1 == 22' "'\$CHR == 22'"
f "arith" '$7 - 400000 < 3 * 2 + 1' "'\$N - 400000 < 3 * 2 + 1'"
f "mod" '$7 % 100 == 0' "'\$N % 100 == 0'"
f "unary minus" '-$8 > 0.05' "'-\$BETA > 0.05'"
f "abs" '($8 / $9 > 20) || ($8 / $9 < -20)' "'abs(\$BETA / \$SE) > 20'"
f "log10" '-log($11) / log(10) > 3' "'-log10(\$P) > 3'"
f "exp sqrt" 'sqrt($6) > exp(-1)' "'sqrt(\$A1_FREQ) > exp(-1)'"
f "joined args" '$11 < 0.01' "'\$P' '<' '0.01'"
check "filter -w" "awk 'NR == 1 || (\$11 < 1e-3 && \$1 == 22)' ws.assoc" "$P filter -w '\$P < 1e-3 && \$CHR == 22' ws.assoc"
check "filter no -H" "awk -F'\t' '\$1 > \$12' ss.tsv" "$P filter '\$1 > \$12' ss.tsv"
check "filter --" "awk -F'\t' 'NR == 1 || \$11 < 0.01' ss.tsv" "$P filter -- '\$P < 0.01' ss.tsv"


# ---- filter: correctly-rounded parsing ----
# Each value sits exactly on its threshold in decimal; an inexact parser
# (plink2's ScanadvDouble is off by an ulp on e.g. "0.9900" and "0.00140")
# flips some of these.
printf 'v\tbeta\tse\n0.9900\t-0.07000\t0.00140\n0.99\t0.06866\t0.00245\n0.990\t0.0700\t0.0014\n0.98999\t0.1\t0.002\n' > bound.tsv
check "boundary le" "awk -F'\t' 'NR == 1 || \$1 <= 0.99' bound.tsv" "$P filter -H '\$v <= 0.99' bound.tsv"
expect "boundary le literal" "$P filter -H '\$v <= 0.99' bound.tsv | cut -f1" 'v\n0.9900\n0.99\n0.990\n0.98999\n'
expect "boundary eq" "$P filter -H '\$v == 0.99' bound.tsv | cut -f1" 'v\n0.9900\n0.99\n0.990\n'
expect "boundary ratio" "$P filter -H 'abs(\$beta / \$se) >= 50' bound.tsv | cut -f2" 'beta\n-0.07000\n0.0700\n0.1\n'
check "boundary ratio awk" "awk -F'\t' 'NR == 1 || \$2 / \$3 >= 50 || \$2 / \$3 <= -50' bound.tsv" "$P filter -H 'abs(\$beta / \$se) >= 50' bound.tsv"

# ---- filter: semantics that intentionally differ from awk ----
printf 'x\ty\tz\n1\tNA\t3\n2\t\t5\n0\t0.5\n10\t1e3\tfoo\n-3\t.5\t-0\n' > ex.tsv
expect "NA lt false" "$P filter -H '\$y < 1' ex.tsv" 'x\ty\tz\n0\t0.5\n-3\t.5\t-0\n'
expect "NA ne true" "$P filter -H '\$y != 1' ex.tsv | cut -f1" 'x\n1\n2\n0\n10\n-3\n'
expect "NA arith missing" "$P filter -H '\$y + 1 > 1' ex.tsv | cut -f1" 'x\n0\n10\n-3\n'
expect "empty field" "$P filter -H '\$y == \"\"' ex.tsv | cut -f1" 'x\n2\n'
expect "past end is empty" "$P filter -H '\$z == \"\"' ex.tsv | cut -f1" 'x\n0\n'
expect "truthy field" "$P filter -H '\$z' ex.tsv | cut -f1" 'x\n1\n2\n10\n'
expect "falsy zero" "$P filter -H '!\$x' ex.tsv | cut -f1" 'x\n0\n'
expect "number vs string literal" "$P filter -H '\$x == \"10.0\"' ex.tsv | cut -f1" 'x\n'
expect "numeric equality" "$P filter -H '\$x == 10.0' ex.tsv | cut -f1" 'x\n10\n'
expect "negative zero" "$P filter -H '\$z == 0' ex.tsv | cut -f1" 'x\n-3\n'
printf 'x\n+5\n-5\n+-5\n0x10\ninf\nnan\n1e400\n.5\n5.\n' > num.tsv
expect "numeric forms" "$P filter -H '\$x == \$x + 0' num.tsv" 'x\n+5\n-5\n.5\n5.\n'
expect "-d ," "printf '1,a\n3,b\n' | $P filter -d , '\$1 > 2' /dev/stdin" '3,b\n'

# ---- filter: errors ----
fails "incomplete" "$P filter '\$1 <' ss.tsv" "expected a value"
fails "regex" "$P filter '\$1 ~ \"a\"' ss.tsv" "regex matching is not supported"
fails "chained" "$P filter '0 < \$1 < 5' ss.tsv" "can't be chained"
fails "bare identifier" "$P filter 'P < 3' ss.tsv" "unknown identifier"
fails "unbalanced paren" "$P filter '(\$1 > 2' ss.tsv" "expected ')'"
fails "unterminated string" "$P filter '\$1 > \"abc' ss.tsv" "unterminated string"
fails "unknown column" "$P filter '\$NOPE > 1' ss.tsv" "column 'NOPE' not found"
fails "option missing arg" "$P filter -d" "missing its argument"

# ---- filter: len() and substr() ----
# Every 3rd/5th line gets a multi-base A1/A2, so SNP and indel filters both
# select something.
awk -F'\t' -v OFS='\t' 'NR > 1 && NR % 3 == 0 {$4 = $4 "TG"} NR > 1 && NR % 5 == 0 {$5 = $5 "A"} 1' ss.tsv > ssindel.tsv
fi2() {  # fi2 NAME AWK_COND PTXT_ARGS...: like f, on ssindel.tsv.
  local name=$1 cond=$2; shift 2
  check "filter $name" "awk -F'\t' 'NR == 1 || ($cond)' ssindel.tsv" "$P filter -H $* ssindel.tsv"
}
fi2 "len snp only" 'length($4) == 1 && length($5) == 1' "'len(\$A1) == 1 && len(\$A2) == 1'"
fi2 "len indel" 'length($4) != length($5)' "'len(\$A1) != len(\$A2)'"
fi2 "length synonym" 'length($4) > 1 && $11 < 0.5' "'length(\$A1) > 1 && \$P < 0.5'"
fi2 "len literal on left" '1 < length($5)' "'1 < len(\$A2)'"
fi2 "len ge le" 'length($4) >= 2 && length($5) <= 1' "'len(\$4) >= 2 && len(\$5) <= 1'"
fi2 "len arith" 'length($4) + length($5) == 4' "'len(\$A1) + len(\$A2) == 4'"
fi2 "len of number field" 'length($11) < 6' "'len(\$P) < 6'"
fi2 "len of computed number" 'length($7 + 0.5) == 8' "'len(\$N + 0.5) == 8'"
fi2 "len no -H" 'length($4) == 3' "'len(\$4) == 3'"
check "len -w" "awk 'length(\$4) == 1 && length(\$5) == 1' ws.assoc" "$P filter -w 'len(\$4) == 1 && len(\$5) == 1' ws.assoc"
f "substr prefix" 'substr($3, 1, 3) == "rs1"' "'substr(\$SNP, 1, 3) == \"rs1\"'"
f "substr literal on left" '"rs1" == substr($3, 1, 3)' "'\"rs1\" == substr(\$SNP, 1, 3)'"
f "substr no n" 'substr($3, 3) < "5"' "'substr(\$SNP, 3) < \"5\"'"
f "substr of number field" 'substr($11, 1, 4) == "0.12"' "'substr(\$P, 1, 4) == \"0.12\"'"
f "substr field vs field" 'substr($4, 1, 1) == substr($5, 1, 1)' "'substr(\$A1, 1, 1) == substr(\$A2, 1, 1)'"
f "substr computed args" 'substr($3, $1 % 3 + 1, 2) == "s1"' "'substr(\$SNP, \$CHR % 3 + 1, 2) == \"s1\"'"
f "substr nested" 'length(substr($3, 3, 2)) == 2' "'len(substr(\$SNP, 3, 2)) == 2'"
f "substr then arith" 'substr($3, 3) + 0 < 500' "'substr(\$SNP, 3) + 0 < 500'"
check "substr chr strip" "awk -F'\t' 'NR == 1 || substr(\$1, 4) + 0 > 20' sschr.tsv" "$P filter -H 'substr(\$CHR, 4) + 0 > 20' sschr.tsv"
printf 'w\thello\nx\t\ny\n' > sub.tsv
# awk writes the expected substring as column 3; every line must match it.
# (-1, 3) is left out: mawk disagrees with gawk and BSD awk there.
for args in '2, 3' '1' '0, 2' '0' '5, 100' '6' '7, 1' '2, 0' '2, -1' '1.7, 1' '2.7, 1' '1, 1.7' '-0.5, 2' '"1", "2"'; do
  awk -F'\t' -v OFS='\t' "{print \$1, \$2, substr(\$2, $args)}" sub.tsv > subx.tsv
  check "substr($args)" "cat subx.tsv" "$P filter 'substr(\$2, $args) == \$3' subx.tsv"
  # Literal comparisons take a specialized path.
  exp=$(awk "BEGIN {print substr(\"hello\", $args)}")
  for op in '==' '!=' '<' '>='; do
    check "substr($args) $op literal" "awk -F'\t' 'substr(\$2, $args) $op \"$exp\"' sub.tsv" "$P filter 'substr(\$2, $args) $op \"$exp\"' sub.tsv"
  done
done
expect "substr start below 1" "$P filter 'substr(\$2, -1, 3) == \"hel\"' sub.tsv" 'w\thello\n'
expect "substr is a string" "$P filter -H 'substr(\$x, 1) == 10' ex.tsv | cut -f1" 'x\n10\n'
expect "substr string vs number" "printf '10.0\n10\n' | $P filter 'substr(\$1, 1) == 10' -" '10\n'
expect "substr of number literal" "printf 'a\n' | $P filter 'substr(1e7, 2, 3) == \"000\" && substr(1/4, 1) == \"0.25\" && len(1/3) == 8' -" 'a\n'
expect "number vs string converts integer" "printf 'a\n' | $P filter '1e7 == \"10000000\"' -" 'a\n'
expect "len past end" "$P filter -H 'len(\$z) == 0' ex.tsv | cut -f1" 'x\n0\n'
expect "len empty field" "$P filter -H 'len(\$y) == 0' ex.tsv | cut -f1" 'x\n2\n'
expect "len NA field" "$P filter -H 'len(\$y) == 2' ex.tsv | cut -f1" 'x\n1\n-3\n'
expect "len of missing" "$P filter -H 'len(\$y + 1) >= 0' ex.tsv | cut -f1" 'x\n0\n10\n-3\n'
expect "substr missing arg" "$P filter -H 'substr(\$z, \$y) == \"3\"' ex.tsv | cut -f1" 'x\n'
expect "substr truthy" "$P filter -H 'substr(\$z, 2)' ex.tsv | cut -f1" 'x\n10\n-3\n'
# ---- filter: $0 ----
check "len \$0" "awk 'length(\$0) > 70' ss.tsv" "$P filter 'len(\$0) > 70' ss.tsv"
check "len \$0 eq" "awk 'length(\$0) == 70' ss.tsv" "$P filter 'length(\$0) == 70' ss.tsv"
check "len \$0 with fields" "awk -F'\t' 'NR == 1 || (length(\$0) > 70 && \$11 < 0.5 && \$4 == \"A\")' ss.tsv" "$P filter -H 'len(\$0) > 70 && \$P < 0.5 && \$A1 == \"A\"' ss.tsv"
check "len \$0 -w" "awk 'length(\$0) > 90' ws.assoc" "$P filter -w 'len(\$0) > 90' ws.assoc"
check "substr \$0" "awk 'substr(\$0, 1, 3) == \"22\t\"' ss.tsv" "$P filter 'substr(\$0, 1, 3) == \"22\t\"' ss.tsv"
check "\$0 vs field" "awk -F'\t' 'substr(\$0, 1, length(\$1)) == \$1' edge.tsv" "$P filter 'substr(\$0, 1, len(\$1)) == \$1' edge.tsv"
expect "\$0 numeric" "printf '5\n10\n1e3\nx\n' | $P filter '\$0 > 6' -" '10\n1e3\n'
expect "\$0 string" "printf 'ab\nabc\n' | $P filter '\$0 == \"abc\"' -" 'abc\n'
expect "\$0 no trailing newline" "printf 'abc' | $P filter 'len(\$0) == 3' -" 'abc\n'
expect "\$0 keeps CR" "printf 'abc\\r\\nde\\n' | $P filter 'len(\$0) == 4' -" 'abc\r\n'
expect "\$0 -w keeps leading space" "printf '  a b\n' | $P filter -w 'len(\$0) == 5 && \$1 == \"a\"' -" '  a b\n'
fails "len no parens" "$P filter 'len > 1' ss.tsv" "expected '('"
fails "len too many" "$P filter 'len(\$1, 2) > 1' ss.tsv" "too many arguments to len()"
fails "substr too few" "$P filter 'substr(\$1) == \"a\"' ss.tsv" "too few arguments to substr()"
fails "substr too many" "$P filter 'substr(\$1, 1, 2, 3) == \"a\"' ss.tsv" "too many arguments to substr()"
fails "top-level comma" "$P filter '\$1 > 1, \$2' ss.tsv" "unexpected input"

# ---- keep / remove ----
check "keep" "awk -F'\t' 'NR == FNR {k[\$1]; next} FNR == 1 || (\$3 in k)' ids.txt ss.tsv" "$P keep SNP ids.txt ss.tsv"
check "remove" "awk -F'\t' 'NR == FNR {k[\$1]; next} FNR == 1 || !(\$3 in k)' ids.txt ss.tsv" "$P remove -H 3 ids.txt ss.tsv"
check "keep -w" "awk 'NR == FNR {k[\$1]; next} FNR == 1 || (\$2 in k)' ids.txt ws.assoc" "$P keep -w SNP ids.txt ws.assoc"
check "keep ID file with extra columns" "awk -F'\t' 'NR == FNR {k[\$1]; next} (\$3 in k)' ids.txt ss.tsv" "awk '{print \$1, \"extra\", \"cols\"}' ids.txt | $P keep 3 - ss.tsv"
# Blank ID lines are ignored (awk would add the empty string as a key).
printf 'b\nzz\n\ny\n' > eids.txt
expect "keep blank IDs" "$P keep -w 2 eids.txt edge.tsv" 'a\tb\t\tdd\n  x  y   z \n'
expect "keep empty field" "$P keep 2 eids.txt edge.tsv" 'a\tb\t\tdd\n'

# ---- count ----
check "count" "awk -F'\t' 'NR > 1 {if (!(\$1 in c)) o[++n] = \$1; c[\$1]++} END {for (i = 1; i <= n; i++) print o[i] \"\t\" c[o[i]]}' ss.tsv" "$P count CHR ss.tsv"
check "count high cardinality" "awk -F'\t' '{if (!(\$7 in c)) o[++n] = \$7; c[\$7]++} END {for (i = 1; i <= n; i++) print o[i] \"\t\" c[o[i]]}' ss.tsv" "$P count 7 ss.tsv"
check "count -w" "awk 'NR > 1 {if (!(\$4 in c)) o[++n] = \$4; c[\$4]++} END {for (i = 1; i <= n; i++) print o[i] \"\t\" c[o[i]]}' ws.assoc" "$P count -w A1 ws.assoc"

# ---- prefix / unprefix ----
check "prefix -H" "cat sschr.tsv" "$P prefix CHR chr ss.tsv"
check "prefix all lines" "awk '{print \"chr\" \$0}' ss.tsv" "$P prefix 1 chr ss.tsv"
check "prefix middle column" "awk -F'\t' -v OFS='\t' '{\$3 = \"X\" \$3} 1' ss.tsv" "$P prefix 3 X ss.tsv"
check "unprefix" "cat ss.tsv" "$P unprefix 1 chr sschr.tsv"
check "unprefix roundtrip" "cat ws.assoc" "$P prefix -w 2 id_ ws.assoc | $P unprefix -w 2 id_ /dev/stdin"
check "prefix past end" "cat edge.tsv" "$P prefix 9 chr edge.tsv"

# ---- tsv ----
check "tsv" "awk -v OFS='\t' '{\$1 = \$1} 1' ws.assoc" "$P tsv ws.assoc"
check "tsv edge" "awk -v OFS='\t' '{\$1 = \$1} 1' edge.tsv" "$P tsv edge.tsv"

# ---- sum ----
check "sum" "awk -F'\t' 'NR > 1 {s += \$8} END {printf \"%.10g\n\", s}' ss.tsv" "$P sum -H 8 ss.tsv"
check "sum skips non-numeric" "echo 15" "printf '5\nNA\n10\n\n' | $P sum 1 -"

# ---- sort ----
# Oracle: the system sort with -s (stable), which BSD and GNU sort both
# support.  100k lines is above ptxt's threshold (65536) for its parallel
# sort/merge paths.
"$GEN" ss 100000 | tail -n +2 | awk 'BEGIN {srand(7)} {print rand() "\t" $0}' | sort -k1,1 | cut -f2- > shuf.tsv
T=$'\t'
s() {  # s NAME SORT_KEYS PTXT_KEYS [PTXT_FILE]
  check "sort $1" "sort -s -t'$T' $2 shuf.tsv" "$P sort $3 ${4:-shuf.tsv}"
}
s "int keys" "-k1,1n -k2,2n" "1:n,2:n"
s "float" "-k11,11g" "11:n"
s "float reverse" "-k11,11gr" "11:nr"
s "string" "-k3,3" "3"
s "string reverse" "-k3,3r" "3:sr"
s "string then float" "-k4,4 -k8,8g" "4,8:n"
s "mixed directions" "-k4,4 -k5,5r -k7,7nr" "4,5:sr,7:nr"
s "keys beyond inline words" "-k4,4 -k5,5 -k6,6g -k12,12gr" "4,5,6:n,12:nr"
s "g alias" "-k6,6g" "6:g"
(head -1 ss.tsv; cat shuf.tsv) > shufh.tsv
check "sort header names" "head -1 ss.tsv; sort -s -t'$T' -k1,1n -k2,2n shuf.tsv" "$P sort CHR:c,POS:n shufh.tsv"
check "sort -H" "head -1 ss.tsv; sort -s -t'$T' -k11,11g shuf.tsv" "$P sort -H 11:n shufh.tsv"
check "sort stdin" "head -1 ss.tsv; sort -s -t'$T' -k11,11g shuf.tsv" "cat shufh.tsv | $P sort P:n -"
gzip -c shufh.tsv > shufh.tsv.gz
check "sort gzip" "head -1 ss.tsv; sort -s -t'$T' -k3,3 shuf.tsv" "$P sort SNP shufh.tsv.gz"
check "sort -w" "head -1 ws.assoc; tail -n +2 ws.assoc | sort -s -k11,11g" "$P sort -w P:n ws.assoc"
check "sort small (serial path)" "sort -s -t'$T' -k11,11g ss.tsv | grep -v CHR" "grep -v CHR ss.tsv | $P sort 11:n /dev/stdin"
check "sort already sorted" "tail -n +2 ss.tsv" "tail -n +2 ss.tsv | $P sort 1:n,2:n /dev/stdin"
# Strings sharing a 16-byte prefix need the full-comparison fallback.
printf 'abcdefghijklmnopqr\nabcdefghijklmnop\nabcdefghijklmnopq\nabcdefgh\nabc\nabcdefghijklmnopqr2\nabcdefghijklmnopqr\n\n' > str.txt
check "sort long-prefix strings" "sort -s -k1,1 str.txt" "$P sort 1 str.txt"
check "sort long-prefix strings reverse" "sort -s -k1,1r str.txt" "$P sort 1:r str.txt"
check "sort missing fields" "sort -s -t'$T' -k2,2 edge.tsv" "$P sort 2 edge.tsv"
# Chromosome order: numeric, X, Y, XY, MT (any chr prefix/case), then other
# names in byte order.
printf 'chrX\t5\nchr2\t1\nchrUn_gl1\t3\nchr10\t1\nMT\t2\nchrY\t1\nchr1\t9\nchrXY\t1\nchrM\t4\nchr2\t0\nchr_random\t1\nx\t3\n' > chr.tsv
expect "sort chromosome order" "$P sort 1:c,2:n chr.tsv | cut -f1 | tr '\n' ' '" 'chr1 chr2 chr2 chr10 x chrX chrY chrXY MT chrM chrUn_gl1 chr_random '
# One key only, so chrX/x and MT/chrM ties keep input order.
expect "sort chromosome reverse" "$P sort 1:cr chr.tsv | cut -f1 | tr '\n' ' '" 'chr_random chrUn_gl1 MT chrM chrXY chrY chrX x chr10 chr2 chr2 chr1 '
# Missing/non-numeric values go last in both directions; -0 == 0 (stable).
printf 'a\t0.5\nb\tNA\nc\t-1\nd\t\ne\t1e-300\nf\t-0\ng\t0\nh\t3\n' > na.tsv
expect "sort numeric missing last" "$P sort 2:n na.tsv | cut -f1 | tr -d '\n'" 'cfgeahbd'
expect "sort numeric reverse missing last" "$P sort 2:nr na.tsv | cut -f1 | tr -d '\n'" 'haefgcbd'
expect "sort empty input" "$P sort 1 /dev/null" ''
fails "sort bad type" "$P sort 1:x chr.tsv" "sort key types"
fails "sort range key" "$P sort 1-3 chr.tsv" "single column"
fails "sort unknown name" "$P sort NOPE chr.tsv" "column 'NOPE' not found"

# ---- join ----
# Lookup tables built from the fixtures: every 3rd SNP gets an annotation.
awk -F'\t' -v OFS='\t' 'NR > 1 && NR % 3 == 0 {print $3, "GENE" (NR % 97), NR / 1000}' ss.tsv > annot.tsv
(printf 'ID\tGENE\tSCORE\n'; cat annot.tsv) > annoth.tsv
J_INNER='NR == FNR {a[$1] = $2 OFS $3; next}'
check "join inner" "awk -F'\t' -v OFS='\t' '$J_INNER FNR > 1 && (\$3 in a) {print \$0, a[\$3]}' annot.tsv ss.tsv" "tail -n +2 ss.tsv | $P join 3 - 1 annot.tsv"
check "join header names" "awk -F'\t' -v OFS='\t' '$J_INNER FNR == 1 {print \$0, \"GENE\", \"SCORE\"; next} (\$3 in a) {print \$0, a[\$3]}' annot.tsv ss.tsv" "$P join SNP ss.tsv ID annoth.tsv"
check "join -a" "awk -F'\t' -v OFS='\t' '$J_INNER FNR == 1 {print \$0, \"GENE\", \"SCORE\"; next} {print \$0, ((\$3 in a) ? a[\$3] : \"NA\" OFS \"NA\")}' annot.tsv ss.tsv" "$P join -a SNP ss.tsv ID annoth.tsv"
check "join -a -e" "awk -F'\t' -v OFS='\t' '$J_INNER FNR == 1 {print \$0, \"GENE\", \"SCORE\"; next} {print \$0, ((\$3 in a) ? a[\$3] : \".\" OFS \".\")}' annot.tsv ss.tsv" "$P join -a -e . SNP ss.tsv ID annoth.tsv"
check "join -O" "awk -F'\t' -v OFS='\t' '$J_INNER FNR > 1 && (\$3 in a) {split(a[\$3], v, \"\t\"); print \$0 \",\" v[1] \",\" v[2]}' annot.tsv ss.tsv" "tail -n +2 ss.tsv | $P join -O , 3 /dev/stdin 1 annot.tsv"
# Key column in the middle of FILE2; only its non-key columns are appended.
awk -F'\t' -v OFS='\t' '{print $2, $1, $3}' annot.tsv > annot_mid.tsv
check "join key mid FILE2" "awk -F'\t' -v OFS='\t' 'NR == FNR {a[\$2] = \$1 OFS \$3; next} FNR > 1 && (\$3 in a) {print \$0, a[\$3]}' annot_mid.tsv ss.tsv" "tail -n +2 ss.tsv | $P join 3 /dev/stdin 2 annot_mid.tsv"
# Composite key, listed in different column orders in the two files.
awk -F'\t' -v OFS='\t' 'NR > 1 && NR % 5 == 0 {print $2, $1, "v" NR}' ss.tsv > cp.tsv
check "join composite" "awk -F'\t' -v OFS='\t' 'NR == FNR {a[\$2 SUBSEP \$1] = \$3; next} FNR > 1 && ((\$1 SUBSEP \$2) in a) {print \$0, a[\$1 SUBSEP \$2]}' cp.tsv ss.tsv" "tail -n +2 ss.tsv | $P join 1,2 /dev/stdin 2,1 cp.tsv"
check "join -w" "awk -v OFS='\t' '$J_INNER FNR > 1 && (\$2 in a) {print \$0, a[\$2]}' annot.tsv ws.assoc" "tail -n +2 ws.assoc | $P join -w 2 /dev/stdin 1 annot.tsv"
gzip -c annoth.tsv > annoth.tsv.gz
check "join gzip FILE2" "$P join -a SNP ss.tsv ID annoth.tsv" "$P join -a SNP ss.tsv ID annoth.tsv.gz"
# One-to-many: FILE1 order kept, matches in FILE2 order; unmatched/short
# FILE1 lines are filled under -a; FILE2 rows lacking the key column are
# skipped; blank FILE2 lines are ignored.
printf 'A\tx\nB\ty\n\nA\tz\nA\tw\nC\n' > dup.tsv
printf 'A\t1\nB\t2\nD\t3\n\nA\t4\nC\t5\n' > main.tsv
expect "join one-to-many" "$P join 1 main.tsv 1 dup.tsv" 'A\t1\tx\nA\t1\tz\nA\t1\tw\nB\t2\ty\nA\t4\tx\nA\t4\tz\nA\t4\tw\nC\t5\t\n'
expect "join one-to-many -a" "$P join -a 2 main.tsv 1 dup.tsv" 'A\t1\tNA\nB\t2\tNA\nD\t3\tNA\n\tNA\nA\t4\tNA\nC\t5\tNA\n'
printf 'x\tA\nB\n' > shortkey.tsv
expect "join FILE2 short rows skipped" "$P join 1 main.tsv 2 shortkey.tsv" 'A\t1\tx\nA\t4\tx\n'
fails "join key count mismatch" "$P join 1,2 main.tsv 1 dup.tsv" "same number of columns"
fails "join range key" "$P join 1-2 main.tsv 1-2 dup.tsv" "single columns"
fails "join unknown name" "$P join SNP ss.tsv NOPE annoth.tsv" "column 'NOPE' not found"
fails "join usage" "$P join -a" "usage: ptxt join"
fails "join -a outside join" "$P filter -a '\$1' main.tsv" "invalid filter expression"

# ---- uniq ----
# Unsorted input with many duplicate lines/keys, and its sorted version.
awk -F'\t' -v OFS='\t' 'NR > 1 {print $1, $4, $5}' ss.tsv > dupl.tsv
sort dupl.tsv > dupl_sorted.tsv
# Sorted, then out of order: exercises the switch from the ascending fast
# path to hashing.
(cat dupl_sorted.tsv; cat dupl.tsv) > dupl_mixed.tsv
U_FIRST='{if (!(k in c)) {o[++n] = k; l[k] = $0} c[k]++}'
U_END='END {for (i = 1; i <= n; i++) {k = o[i]; if (c[k] >= MIN && c[k] <= MAX) print (CT ? c[k] "\t" : "") l[k]}}'
check "uniq whole line" "awk '!s[\$0]++' dupl.tsv" "$P uniq dupl.tsv"
check "uniq whole line mixed order" "awk '!s[\$0]++' dupl_mixed.tsv" "$P uniq dupl_mixed.tsv"
check "uniq sorted" "uniq dupl_sorted.tsv" "$P uniq dupl_sorted.tsv"
check "uniq sorted -r" "uniq -d dupl_sorted.tsv" "$P uniq -r dupl_sorted.tsv"
check "uniq sorted -u" "uniq -u dupl_sorted.tsv" "$P uniq -u dupl_sorted.tsv"
check "uniq sorted -c" "uniq -c dupl_sorted.tsv | awk '{c = \$1; sub(/^ *[0-9]+ /, \"\"); print c \"\t\" \$0}'" "$P uniq -c dupl_sorted.tsv"
check "uniq -c mixed order" "awk -v MIN=1 -v MAX=1e18 -v CT=1 '{k = \$0} $U_FIRST $U_END' dupl_mixed.tsv" "$P uniq -c dupl_mixed.tsv"
check "uniq column" "awk -F'\t' '!s[\$4]++' ss.tsv" "$P uniq 4 ss.tsv"
check "uniq column name" "awk -F'\t' 'NR == 1 || !s[\$4]++' ss.tsv" "$P uniq A1 ss.tsv"
check "uniq composite" "awk -F'\t' '!s[\$5, \$4]++' ss.tsv" "$P uniq 5,4 ss.tsv"
check "uniq composite mixed order" "awk -F'\t' '!s[\$1, \$3]++' dupl_mixed.tsv" "$P uniq 1,3 dupl_mixed.tsv"
check "uniq -w" "awk '!s[\$4]++' ws.assoc" "$P uniq -w 4 ws.assoc"
check "uniq -d" "tr '\t' , < ss.tsv | awk -F, '!s[\$5]++'" "tr '\t' , < ss.tsv | $P uniq -d , 5 -"
check "uniq column -c" "awk -F'\t' -v MIN=1 -v MAX=1e18 -v CT=1 'NR == 1 {print \"count\t\" \$0; next} {k = \$1 SUBSEP \$4} $U_FIRST $U_END' ss.tsv" "$P uniq -c CHR,A1 ss.tsv"
check "uniq column -r" "awk -F'\t' -v MIN=2 -v MAX=1e18 -v CT=0 '{k = \$1 SUBSEP \$2} $U_FIRST $U_END' dupl_mixed.tsv" "$P uniq -r 1,2 dupl_mixed.tsv"
check "uniq column -c -u" "awk -F'\t' -v MIN=1 -v MAX=1 -v CT=1 '{k = \$7} $U_FIRST $U_END' ss.tsv" "$P uniq -c -u 7 ss.tsv"
# A missing field is an empty key, as in awk; blank lines are keys too.
printf 'a\tb\nc\nd\t\ne\n\nf\tb\n\n' > ukeys.tsv
expect "uniq missing field" "$P uniq 2 ukeys.tsv" 'a\tb\nc\n'
expect "uniq blank lines" "$P uniq -c ukeys.tsv" '1\ta\tb\n1\tc\n1\td\t\n1\te\n2\t\n1\tf\tb\n'
fails "uniq -r -u" "$P uniq -r -u dupl.tsv" "mutually exclusive"
fails "uniq usage" "$P uniq -c" "usage: ptxt uniq"
fails "uniq unknown name" "$P uniq NOPE ss.tsv" "column 'NOPE' not found"

# ---- compressed input ----
gzip -c ss.tsv > ss.tsv.gz
check "gzip" "$P filter '\$P < 0.01' ss.tsv" "$P filter '\$P < 0.01' ss.tsv.gz"
check "gzip cut" "cut -f2,11 ss.tsv" "$P cut 2,11 ss.tsv.gz"
if command -v bgzip > /dev/null; then
  bgzip -c ss.tsv > ss.bgz.gz
  bgzip -c wide.vcf > wide.vcf.gz
  check "bgzf filter" "$P filter '\$P < 0.01' ss.tsv" "$P filter '\$P < 0.01' ss.bgz.gz"
  check "bgzf cut" "cut -f1-8 wide.vcf" "$P cut 1-8 wide.vcf.gz"
  check "bgzf keep" "$P keep 3 ids.txt ss.tsv" "$P keep 3 ids.txt ss.bgz.gz"
else
  echo "(bgzip not found; skipping BGZF tests)"
fi
if command -v zstd > /dev/null; then
  zstd -q -c ss.tsv > ss.tsv.zst
  check "zstd count" "$P count 1 ss.tsv" "$P count 1 ss.tsv.zst"
  check "zstd uniq" "$P uniq -c 1,4 ss.tsv" "$P uniq -c 1,4 ss.tsv.zst"
else
  echo "(zstd not found; skipping zstd tests)"
fi

# ---- help and usage errors ----
fails "unknown subcommand" "$P frobnicate x" "unrecognized subcommand"
fails "help for unknown subcommand" "$P help frobnicate" "unrecognized subcommand"
# 'ptxt' alone prints the summary but exits nonzero, like plink2.
fails "no arguments" "$P" "ptxt <subcommand> -h"
if ! out=$("$P" --help) || [[ "$out" != *"Subcommands:"* ]]; then
  echo "FAIL: ptxt --help"; exit 1
fi
PASS=$((PASS + 1))
# Every subcommand in the summary has full help: 'SUB -h', 'help SUB' and
# '-h' after other options all give the same text, starting with SUB's
# synopsis, and it fits in 80 columns.
SUBS=$("$P" --help | awk '/^     [a-z]+  / {print $1}')
[[ $(wc -w <<< "$SUBS") -eq 14 ]] || { echo "FAIL: summary lists $(wc -w <<< "$SUBS") subcommands"; exit 1; }
for sub in $SUBS; do
  check "$sub -h = help $sub" "$P help $sub" "$P $sub -h"
  if ! out=$("$P" $sub -h) || [[ "$out" != *"ptxt $sub "* ]]; then
    echo "FAIL: $sub -h doesn't show its synopsis"; exit 1
  fi
  if awk 'length > 80 {bad = 1} END {exit !bad}' <<< "$out"; then
    echo "FAIL: $sub -h has lines over 80 columns"; exit 1
  fi
  PASS=$((PASS + 1))
done
awk 'length > 80 {bad = 1} END {exit !bad}' <<< "$("$P" --help)" && { echo "FAIL: summary has lines over 80 columns"; exit 1; }
# The version string isn't hardcoded here, so bumping it doesn't break tests.
VERSION=$("$P" -v)
[[ "$VERSION" =~ ^ptxt\ v[0-9]+\.[0-9]+\.[0-9]+\ \([0-9]{1,2}\ [A-Z][a-z]{2}\ [0-9]{4}\)$ ]] || { echo "FAIL: -v output malformed: $VERSION"; exit 1; }
PASS=$((PASS + 1))
for v in -V --version; do
  expect "$v" "$P $v" "$VERSION\\n"
done
[[ $("$P" --help | head -1) == "$VERSION"* ]] || { echo "FAIL: banner lacks version"; exit 1; }
PASS=$((PASS + 1))
fails "-v with arguments" "$P -v sort" "takes no arguments"
check "-h after options" "$P sort -h" "$P sort -H -w -h"
check "-h after uniq flags" "$P uniq -h" "$P uniq -c -h"
check "--help" "$P join -h" "$P join -a --help"
fails "subcommand alone" "$P sort" "ptxt sort [options]"
fails "usage hint" "$P cut 1" "\"ptxt cut -h\" for details"
fails "missing file" "$P cut 1 does_not_exist.tsv" "failed to open"

echo "All $PASS ptxt tests passed."
