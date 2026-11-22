#!/bin/bash

# --adjust-file 'log10' with --output-min-p.
#
# -log10(p) columns used to print the --output-min-p string (a p-value, not a
# -log10 value) for p-values below it.

set -exo pipefail

plink2="$1/plink2 $2 $3"

printf '#CHROM\tPOS\tID\tA1\tTEST\tP\n1\t100\tv1\tG\tADD\t1e-30\n1\t200\tv2\tG\tADD\t0.5\n1\t300\tv3\tG\tADD\t0.01\n' > tmp_assoc.txt
$plink2 --adjust-file tmp_assoc.txt log10 test=ADD --output-min-p 1e-10 --out tmp_adj > /dev/null
# v1's unadjusted -log10(p) must be 30, not "1e-10".
awk -F '\t' 'NR == 1 {for (i = 1; i <= NF; ++i) {if ($i == "NEG_LOG10_UNADJ") {c = i}}; next} $2 == "v1" {print $c}' tmp_adj.adjusted | grep -qx '30'
