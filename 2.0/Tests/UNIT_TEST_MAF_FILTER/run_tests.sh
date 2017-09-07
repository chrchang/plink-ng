#!/bin/bash

# http://stackoverflow.com/questions/5195607/checking-bash-exit-status-of-several-commands-efficiently
function run {
    "$@"
    local status=$?
    if [ $status -ne 0 ]; then
	exit $status
    fi
    return $status
}


run plink --simulate simulate.txt --simulate-missing 0.098 --out tmp_data
run plink --bfile tmp_data --keep-allele-order --geno 0.1 --maf 0.01 --make-bed --out plink_filtered
run $1/pgen_compress tmp_data.bed tmp_data.pgen 2000
run $1/plink2 $2 --bfile tmp_data --geno 0.1 --maf 0.01 --make-bpgen --out new_filtered
run diff -q plink_filtered.bim new_filtered.bim
run $1/pgen_compress -u new_filtered.pgen new_filtered.bed
run diff -q plink_filtered.bed new_filtered.bed
run $1/plink2 $2 --bpfile tmp_data --geno 0.1 --maf 0.01 --make-bpgen --out new_filtered2
run diff -q plink_filtered.bim new_filtered2.bim
run diff -q new_filtered.pgen new_filtered2.pgen
