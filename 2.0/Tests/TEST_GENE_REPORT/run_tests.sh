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
    awk 'NR > 1 {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $8}' plink2_gr.gene.report > flat20.txt
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

# --extract, which restricts the report rows rather than a variant_include.
# The ID list has to be large enough that the report's matches would collide
# with it if it were allocated in the wrong part of the arena.
awk 'NR > 1 && (NR % 3) == 0 {print $3}' report19.txt > extract.txt
test $(wc -l < extract.txt) -ge 1000

plink --gene-report report19.txt genes.txt --extract extract.txt --out plink_gr
$1/plink2 $2 $3 --gene-report report20.txt genes.txt --extract extract.txt --out plink2_gr
check

plink --gene-report report19.txt genes.txt --extract extract.txt --pfilter 0.1 --gene-list-border 25 --out plink_gr
$1/plink2 $2 $3 --gene-report report20.txt genes.txt --extract extract.txt --pfilter 0.1 --gene-list-border 25 --out plink2_gr
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
awk 'NR > 1 {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $8}' plink2_alt.gene.report > flat_alt.txt
$1/plink2 $2 $3 --gene-report report20.txt genes.txt --out plink2_gr
awk 'NR > 1 {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $8}' plink2_gr.gene.report > flat20.txt
diff -q flat_alt.txt flat20.txt

awk '{print $1, $2 - 1, $3, $4}' genes.txt > genes0.txt
$1/plink2 $2 $3 --gene-report report20.txt genes0.txt 0based --out plink2_zero
diff -q plink2_zero.gene.report plink2_gr.gene.report

# Column subset, and Zstd output.
$1/plink2 $2 $3 --gene-report report20.txt genes.txt cols=pos --out plink2_cols
head -n 1 plink2_cols.gene.report | grep -qx '#GENE	POS	ID'
awk 'NR > 1 {print $1 "\t" $3}' plink2_cols.gene.report > flat_cols.txt
awk 'NR > 1 {print $1 "\t" $8}' plink2_gr.gene.report > flat_gid.txt
diff -q flat_cols.txt flat_gid.txt

$1/plink2 $2 $3 --gene-report report20.txt genes.txt zs --out plink2_zs
$1/plink2 $2 $3 --zst-decompress plink2_zs.gene.report.zst > plink2_zs.gene.report
diff -q plink2_zs.gene.report plink2_gr.gene.report

# Nonstandard contigs, including two 15000-character names that differ only in
# their last character.  No main dataset is loaded, so the contig names must be
# taken from the gene range file.  One contig appears only in the gene list and
# another only in the report, and some p-values are NA.
awk 'BEGIN {
  long = "q"
  for (i = 1; i < 15000; ++i) long = long substr("ACGTxyz_", 1 + (i * 7) % 8, 1)
  longb = substr(long, 1, 14999) "Q"
  n = split("ctgA " long " " longb " ctgGenesOnly", c, " ")
  s = 4242
  for (i = 0; i < 120; ++i) {
    s = (1103515245 * s + 12345) % 2147483648
    k = 1 + (int(s / 512) % n)
    s = (1103515245 * s + 12345) % 2147483648
    start = 1 + (int(s / 512) % 900000)
    s = (1103515245 * s + 12345) % 2147483648
    end = start + 1 + (int(s / 512) % 200000)
    s = (1103515245 * s + 12345) % 2147483648
    printf "%s %d %d XGENE%d\n", c[k], start, end, int(s / 512) % 30 > "genes_extra.txt"
  }
  c[4] = "ctgReportOnly"
  for (i = 0; i < 2000; ++i) {
    s = (1103515245 * s + 12345) % 2147483648
    k = 1 + (int(s / 512) % n)
    s = (1103515245 * s + 12345) % 2147483648
    bp = 1 + (int(s / 512) % 1100000)
    s = (1103515245 * s + 12345) % 2147483648
    p = (int(s / 512) % 1000000) / 1000000.0
    if (i % 97 == 0) p = "NA"
    printf "%s %d xs%d %s\n", c[k], bp, i, p > "report_rows_extra19.txt"
    printf "%s\t%d\txs%d\t%s\n", c[k], bp, i, p > "report_rows_extra20.txt"
  }
}'
cat genes.txt genes_extra.txt > genes_all.txt
cat report19.txt report_rows_extra19.txt > report_extra19.txt
cat report20.txt report_rows_extra20.txt > report_extra20.txt

plink --gene-report report_extra19.txt genes_all.txt --allow-extra-chr --out plink_gr
$1/plink2 $2 $3 --gene-report report_extra20.txt genes_all.txt --out plink2_gr
check
# Both long contigs must have matched.
test $(awk -F'\t' 'length($2) == 15000' flat20.txt | cut -f 2 | sort -u | wc -l) -eq 2
# NA p-values must be written as NA, not as exp(kLnPvalError).
test $(awk -F'\t' 'NR > 1 && $NF == "NA"' plink2_gr.gene.report | wc -l) -ge 10
test $(awk -F'\t' 'NR > 1 && $NF != "NA" && $NF > 1' plink2_gr.gene.report | wc -l) -eq 0

plink --gene-report report_extra19.txt genes_all.txt --allow-extra-chr --gene-list-border 25 --pfilter 0.5 --out plink_gr
$1/plink2 $2 $3 --gene-report report_extra20.txt genes_all.txt --gene-list-border 25 --pfilter 0.5 --out plink2_gr
check

plink --gene-report report_extra19.txt genes_all.txt --allow-extra-chr --chr ctgA --out plink_gr
$1/plink2 $2 $3 --gene-report report_extra20.txt genes_all.txt --chr ctgA --out plink2_gr
diff -q <(awk -f flatten19.awk plink_gr.range.report) <(awk 'NR > 1 {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $8}' plink2_gr.gene.report)
test $(awk 'NR > 1' plink2_gr.gene.report | wc -l) -ge 100
test $(awk 'NR > 1 {print $2}' plink2_gr.gene.report | sort -u) = ctgA

# --strict-extra-chr still rejects them.
if $1/plink2 $2 $3 --strict-extra-chr --gene-report report_extra20.txt genes_all.txt --out plink2_strict; then
    exit 1
fi
