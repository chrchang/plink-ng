#!/bin/bash

set -exo pipefail

plink --simulate simulate.txt --simulate-missing 0.098 --out tmp_data
plink --bfile tmp_data --keep-allele-order --geno 0.1 --maf 0.01 --make-bed --out plink_filtered
$1/pgen_compress tmp_data.bed tmp_data.pgen 2000
$1/plink2 $2 $3 --bfile tmp_data --geno 0.1 --maf 0.01 --make-bpgen --out new_filtered
diff -q plink_filtered.bim new_filtered.bim
$1/pgen_compress -u new_filtered.pgen new_filtered.bed
diff -q plink_filtered.bed new_filtered.bed
$1/plink2 $2 $3 --bpfile tmp_data --geno 0.1 --maf 0.01 --make-bpgen --out new_filtered2
diff -q plink_filtered.bim new_filtered2.bim
diff -q new_filtered.pgen new_filtered2.pgen
