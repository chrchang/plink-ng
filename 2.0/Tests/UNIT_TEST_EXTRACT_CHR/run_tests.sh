#!/bin/bash

set -exo pipefail

plink --simulate simulate.txt --out tmp_data

cat tmp_data.bim | head -n 147 > tmp_data1.bim
cat tmp_data.bim | tail -n 150 | sed 's/^1/2/g' > tmp_data24.bim
cat tmp_data24.bim | head -n 100 > tmp_data2.bim
cat tmp_data24.bim | tail -n 50 | sed 's/^2/4/g' > tmp_data4.bim

cat tmp_data1.bim tmp_data2.bim tmp_data4.bim > tmp_data.bim
$1/pgen_compress tmp_data.bed tmp_data.pgen 2000

for c in 1 2 4
do
    plink --bfile tmp_data --chr $c --make-bed --out plink_c$c
    $1/plink2 $2 $3 --bfile tmp_data --chr $c --make-bed --out plink2_c$c
    diff -q plink_c$c.bim plink2_c$c.bim
    diff -q plink_c$c.bed plink2_c$c.bed
    $1/plink2 $2 $3 --bpfile tmp_data --chr $c --make-bed --out plink2x_c$c
    diff -q plink2x_c$c.bed plink2_c$c.bed
    diff -q plink2x_c$c.bim plink2_c$c.bim
done
