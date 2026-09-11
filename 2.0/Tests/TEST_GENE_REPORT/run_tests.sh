#!/bin/bash

set -exo pipefail

# Deterministic pseudorandom gene list and association report, built with a
# plain LCG so that every awk implementation produces the same files.  The low
# bits of an LCG have short periods, so every value is taken from bits 9-30.
# Only autosomes are used, so that plink 1.9's numeric chromosome codes can be
# compared against plink 2.0's directly.
awk 'BEGIN {
  s = 12345
  for (i = 0; i < 400; ++i) {
    s = (1103515245 * s + 12345) % 2147483648
    chrom = 1 + (int(s / 512) % 22)
    s = (1103515245 * s + 12345) % 2147483648
    start = 1 + (int(s / 512) % 900000)
    s = (1103515245 * s + 12345) % 2147483648
    end = start + 1 + (int(s / 512) % 200000)
    s = (1103515245 * s + 12345) % 2147483648
    printf "%d %d %d GENE%d\n", chrom, start, end, int(s / 512) % 80
  }
}' > genes.txt

awk 'BEGIN {
  s = 999
  print "CHR BP SNP P" > "report19.txt"
  print "#CHROM\tPOS\tID\tP" > "report20.txt"
  for (i = 0; i < 5000; ++i) {
    s = (1103515245 * s + 12345) % 2147483648
    chrom = 1 + (int(s / 512) % 22)
    s = (1103515245 * s + 12345) % 2147483648
    bp = 1 + (int(s / 512) % 1100000)
    s = (1103515245 * s + 12345) % 2147483648
    p = (int(s / 512) % 1000000) / 1000000.0
    printf "%d %d rs%d %g\n", chrom, bp, i, p > "report19.txt"
    printf "%d\t%d\trs%d\t%g\n", chrom, bp, i, p > "report20.txt"
  }
}'

printf 'GENE1\nGENE7\nGENE13\nGENE44\n' > subset.txt

check() {
    awk -f flatten19.awk plink_gr.range.report > flat19.txt
    awk 'NR > 1 {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $7}' plink2_gr.gene.report > flat20.txt
    diff -q flat19.txt flat20.txt
    # Guard against a vacuous comparison.
    test $(wc -l < flat20.txt) -ge 100
}

plink --gene-report report19.txt genes.txt --out plink_gr
$1/plink2 $2 $3 --gene-report report20.txt genes.txt --out plink2_gr
check

plink --gene-report report19.txt genes.txt --gene-list-border 25 --out plink_gr
$1/plink2 $2 $3 --gene-report report20.txt genes.txt --gene-list-border 25 --out plink2_gr
check

plink --gene-report report19.txt genes.txt --gene-subset subset.txt --out plink_gr
$1/plink2 $2 $3 --gene-report report20.txt genes.txt --gene-subset subset.txt --out plink2_gr
check

plink --gene-report report19.txt genes.txt --pfilter 0.05 --out plink_gr
$1/plink2 $2 $3 --gene-report report20.txt genes.txt --pfilter 0.05 --out plink2_gr
check

plink --gene-report report19.txt genes.txt --gene-list-border 10 --gene-subset subset.txt --pfilter 0.2 --out plink_gr
$1/plink2 $2 $3 --gene-report report20.txt genes.txt --gene-list-border 10 --gene-subset subset.txt --pfilter 0.2 --out plink2_gr
check

# plink 1.9 field names, and a report with no p-value column.
awk 'NR > 1 {print $1, $2, $3}' report19.txt > report_nop.txt
sed -i.bak '1i\
CHR BP SNP
' report_nop.txt
$1/plink2 $2 $3 --gene-report report_nop.txt genes.txt --out plink2_nop
# The P column must be dropped from the output when the input lacks one.
head -n 1 plink2_nop.gene.report | grep -qv 'P$'

# Explicit field names, and a 0-based (BED-style) gene range file.
awk '{print $1, $2, $3, $4}' report20.txt | tail -n +2 > report_alt.txt
sed -i.bak '1i\
C POSITION VARIANT PVAL
' report_alt.txt
$1/plink2 $2 $3 --gene-report report_alt.txt genes.txt --gene-report-chr-field C --gene-report-pos-field POSITION --gene-report-id-field VARIANT --gene-report-p-field PVAL --out plink2_alt
awk 'NR > 1 {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $7}' plink2_alt.gene.report > flat_alt.txt
$1/plink2 $2 $3 --gene-report report20.txt genes.txt --out plink2_gr
awk 'NR > 1 {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $7}' plink2_gr.gene.report > flat20.txt
diff -q flat_alt.txt flat20.txt

awk '{print $1, $2 - 1, $3, $4}' genes.txt > genes0.txt
$1/plink2 $2 $3 --gene-report report20.txt genes0.txt 0based --out plink2_zero
diff -q plink2_zero.gene.report plink2_gr.gene.report

# Column subset, and Zstd output.
$1/plink2 $2 $3 --gene-report report20.txt genes.txt cols=pos --out plink2_cols
head -n 1 plink2_cols.gene.report | grep -qx '#GENE	ID	POS'
awk 'NR > 1 {print $1 "\t" $2}' plink2_cols.gene.report > flat_cols.txt
awk 'NR > 1 {print $1 "\t" $7}' plink2_gr.gene.report > flat_gid.txt
diff -q flat_cols.txt flat_gid.txt

$1/plink2 $2 $3 --gene-report report20.txt genes.txt zs --out plink2_zs
$1/plink2 $2 $3 --zst-decompress plink2_zs.gene.report.zst > plink2_zs.gene.report
diff -q plink2_zs.gene.report plink2_gr.gene.report
