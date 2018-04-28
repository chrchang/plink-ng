#!/bin/bash

set -exo pipefail

plink --simulate simulate.txt --simulate-missing 0.01 --out tmp_data
plink --bfile tmp_data --a2-allele tmp_data.bim 5 2 --freqx
$1/plink2 $2 $3 --bfile tmp_data --geno-counts

cat plink.frqx | tail -n +2 > plink.frqx.headerless

cat plink2.gcount | tail -n +2 > plink2.gcount.headerless1
diff -q plink.frqx.headerless plink2.gcount.headerless1
$1/pgen_compress tmp_data.bed tmp_data.pgen 2000
$1/plink2 $2 $3 --bpfile tmp_data --geno-counts

cat plink2.gcount | tail -n +2 > plink2.gcount.headerless2
diff -q plink.frqx.headerless plink2.gcount.headerless2
