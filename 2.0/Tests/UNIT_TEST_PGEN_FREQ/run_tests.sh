#!/bin/bash

# shamelessly stolen from kickstart.sh
# http://stackoverflow.com/questions/5195607/checking-bash-exit-status-of-several-commands-efficiently
function run {
    "$@"
    local status=$?
    if [ $status -ne 0 ]; then
	exit $status
    fi
    return $status
}


run plink --simulate simulate.txt --simulate-missing 0.01 --out tmp_data
run plink --bfile tmp_data --a2-allele tmp_data.bim 5 2 --freqx
run $1/plink2 $2 --bfile tmp_data --geno-counts

run cat plink.frqx | tail -n +2 > plink.frqx.headerless

run cat plink2.gcount | tail -n +2 > plink2.gcount.headerless1
run diff -q plink.frqx.headerless plink2.gcount.headerless1
run $1/pgen_compress tmp_data.bed tmp_data.pgen 2000
run $1/plink2 $2 --bpfile tmp_data --geno-counts

run cat plink2.gcount | tail -n +2 > plink2.gcount.headerless2
run diff -q plink.frqx.headerless plink2.gcount.headerless2
