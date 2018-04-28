#!/bin/bash

set -exo pipefail

plink --dummy 4321 333 0.05 --out tmp_data
plink --bfile tmp_data --mind 0.05 --make-bed --keep-allele-order --out tmp_filtered

cat tmp_filtered.fam | cut -d ' ' -f 1-2 > keep.txt
$1/plink2 $2 $3 --bfile tmp_data --keep keep.txt --make-bed --out bed_filtered

cat tmp_filtered.fam | tr ' ' '\t' > tmp_filtered_tabs.fam

diff -q tmp_filtered_tabs.fam bed_filtered.fam
diff -q tmp_filtered.bed bed_filtered.bed
$1/pgen_compress tmp_data.bed tmp_data.pgen 4321
$1/plink2 $2 $3 --bpfile tmp_data --keep keep.txt --make-bed --out pgen_filtered
diff -q bed_filtered.fam pgen_filtered.fam
diff -q bed_filtered.bed pgen_filtered.bed
