#!/bin/bash

set -exo pipefail

# Deterministic pseudorandom inputs, built with a plain LCG so that every awk
# implementation produces the same files.  The low bits of an LCG have short
# periods, so every value is taken from bits 9-30.  Only autosomes are used, so
# that plink 1.9's numeric chromosome codes can be compared against plink 2.0's
# directly.
awk 'BEGIN {
  s = 4242
  for (i = 0; i < 250; ++i) {
    s = (1103515245 * s + 12345) % 2147483648
    chrom = 1 + (int(s / 512) % 22)
    s = (1103515245 * s + 12345) % 2147483648
    start = 1 + (int(s / 512) % 900000)
    s = (1103515245 * s + 12345) % 2147483648
    end = start + 1 + (int(s / 512) % 200000)
    s = (1103515245 * s + 12345) % 2147483648
    printf "%d %d %d R%d\n", chrom, start, end, int(s / 512) % 40
  }
}' > ranges.txt

awk 'BEGIN {
  s = 777
  for (i = 0; i < 40; ++i) {
    s = (1103515245 * s + 12345) % 2147483648
    chrom = 1 + (int(s / 512) % 22)
    s = (1103515245 * s + 12345) % 2147483648
    start = 1 + (int(s / 512) % 900000)
    s = (1103515245 * s + 12345) % 2147483648
    printf "%d %d %d F%d\n", chrom, start, start + 1 + (int(s / 512) % 300000), i
  }
}' > filter.txt

awk 'BEGIN {
  split("coding utr intron promoter enhancer", a, " ")
  s = 31337
  for (i = 0; i < 3000; ++i) {
    s = (1103515245 * s + 12345) % 2147483648
    chrom = 1 + (int(s / 512) % 22)
    s = (1103515245 * s + 12345) % 2147483648
    bp = 1 + (int(s / 512) % 1100000)
    s = (1103515245 * s + 12345) % 2147483648
    p = (int(s / 512) % 1000000) / 1000000.0
    printf "%d %d rs%d %g\n", chrom, bp, i, p > "report.txt"
    s = (1103515245 * s + 12345) % 2147483648
    if ((int(s / 512) % 10) < 6) {
      s = (1103515245 * s + 12345) % 2147483648
      j = 1 + (int(s / 512) % 5)
      s = (1103515245 * s + 12345) % 2147483648
      k = 1 + (int(s / 512) % 5)
      if (j == k) {
        printf "rs%d %s\n", i, a[j] > "attrib.txt"
      } else {
        printf "rs%d %s %s\n", i, a[j], a[k] > "attrib.txt"
      }
    }
    s = (1103515245 * s + 12345) % 2147483648
    if ((int(s / 512) % 3) == 0) {
      printf "rs%d\n", i > "snps.txt"
    }
  }
  print "CHR BP SNP P" > "header.txt"
  for (i = 0; i < 40; ++i) {
    if ((i % 4) == 0) {
      printf "R%d\n", i > "subset.txt"
    }
  }
}'
cat header.txt report.txt > report_full.txt
mv report_full.txt report.txt

# plink 1.9 pads DIST/SGN to fixed widths and separates columns with spaces;
# plink 2.0 emits tabs and a '#'-prefixed header.  Reduce both to the same
# whitespace-normalized form.
norm19() {
    awk '{ $1 = $1; print }' OFS='\t' "$1" > "$2"
}
norm20() {
    sed '1s/^#//' "$1" | awk '{ $1 = $1; print }' OFS='\t' > "$2"
}

check() {
    norm19 plink_annot.annot n19.txt
    norm20 plink2_annot.annot n20.txt
    diff -q n19.txt n20.txt
    # Guard against a vacuous comparison.
    test $(wc -l < n20.txt) -ge 100
}

run_both() {
    plink --annotate report.txt "$@" --out plink_annot
    $D/plink2 $T1 $T2 --annotate report.txt "$@" --out plink2_annot
    check
}

D=$1
T1=$2
T2=$3

run_both ranges=ranges.txt
run_both attrib=attrib.txt
run_both ranges=ranges.txt attrib=attrib.txt
run_both ranges=ranges.txt attrib=attrib.txt block
run_both ranges=ranges.txt NA
run_both ranges=ranges.txt prune
run_both ranges=ranges.txt minimal
run_both ranges=ranges.txt filter=filter.txt
run_both ranges=ranges.txt attrib=attrib.txt snps=snps.txt
run_both ranges=ranges.txt subset=subset.txt

# --pfilter, and the border flags (named --border in plink 1.9).
plink --annotate report.txt ranges=ranges.txt attrib=attrib.txt --pfilter 0.1 --out plink_annot
$D/plink2 $T1 $T2 --annotate report.txt ranges=ranges.txt attrib=attrib.txt --pfilter 0.1 --out plink2_annot
check

# The parenthesized distances inside ANNOT are rounded to 4 significant digits
# by plink 1.9, so 'minimal' is used here; the DIST/SGN columns are written at
# full precision by both, and are compared.
plink --annotate report.txt ranges=ranges.txt distance minimal --border 40 --out plink_annot
$D/plink2 $T1 $T2 --annotate report.txt ranges=ranges.txt distance minimal --annotate-border 40 --out plink2_annot
check

plink --annotate report.txt ranges=ranges.txt attrib=attrib.txt block distance --border 25 --out plink_annot
$D/plink2 $T1 $T2 --annotate report.txt ranges=ranges.txt attrib=attrib.txt block distance --annotate-border 25 --out plink2_annot
check

# plink 2.0 extras: 0-based interval input, explicit field names, and Zstd
# output must all reproduce the default output.
$D/plink2 $T1 $T2 --annotate report.txt ranges=ranges.txt attrib=attrib.txt --out plink2_annot
awk '{print $1, $2 - 1, $3, $4}' ranges.txt > ranges0.txt
$D/plink2 $T1 $T2 --annotate report.txt ranges=ranges0.txt attrib=attrib.txt 0based --out plink2_zero
diff -q plink2_zero.annot plink2_annot.annot

awk 'NR > 1 {print $1, $2, $3, $4}' report.txt > report_alt.txt
sed -i.bak '1i\
C POSITION VARIANT PVAL
' report_alt.txt
$D/plink2 $T1 $T2 --annotate report_alt.txt ranges=ranges.txt attrib=attrib.txt --annotate-chr-field C --annotate-pos-field POSITION --annotate-id-field VARIANT --annotate-p-field PVAL --out plink2_alt
cut -f5 plink2_alt.annot > annot_alt.txt
cut -f5 plink2_annot.annot > annot_def.txt
diff -q annot_alt.txt annot_def.txt

$D/plink2 $T1 $T2 --annotate report.txt ranges=ranges.txt attrib=attrib.txt zs --out plink2_zs
$D/plink2 $T1 $T2 --zst-decompress plink2_zs.annot.zst > plink2_zs.annot
diff -q plink2_zs.annot plink2_annot.annot
