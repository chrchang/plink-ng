#!/bin/bash

set -exo pipefail

plink --simulate simulate.txt --simulate-missing 0.01 --out tmp_data

# Force last few variants to chrX/chrY.  (chrM is expected to be different.)
# (Unfortunately, "head -n -2" doesn't work on macOS.)
head -n 295 tmp_data.bim > tmp_data2.bim
tail -n 2 tmp_data.bim > last_two.txt
head -n 1 last_two.txt | sed 's/^1/X/g' >> tmp_data2.bim
tail -n 1 last_two.txt | sed 's/^1/Y/g' >> tmp_data2.bim
cp tmp_data2.bim tmp_data.bim
# Force some samples to male or no-sex.
head -n 1234 tmp_data.fam > tmp_data2.fam
head -n 1357 tmp_data.fam | tail -n +1235 | sed 's/2 1$/0 1/g' >> tmp_data2.fam
tail -n +1358 tmp_data.fam | sed 's/2 1$/1 1/g' >> tmp_data2.fam
cp tmp_data2.fam tmp_data.fam

plink --bfile tmp_data --a2-allele tmp_data.bim 5 2 --freqx
$1/plink2 $2 $3 --bfile tmp_data --output-chr 26 --geno-counts

cat plink.frqx | tail -n +2 > plink.frqx.headerless

cat plink2.gcount | tail -n +2 > plink2.gcount.headerless1
diff -q plink.frqx.headerless plink2.gcount.headerless1
$1/pgen_compress tmp_data.bed tmp_data.pgen 2000
$1/plink2 $2 $3 --bpfile tmp_data --output-chr 26 --geno-counts

cat plink2.gcount | tail -n +2 > plink2.gcount.headerless2
diff -q plink.frqx.headerless plink2.gcount.headerless2
