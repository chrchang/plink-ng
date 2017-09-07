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

run plink --dummy 4321 333 0.05 --out tmp_data
run plink --bfile tmp_data --mind 0.05 --make-bed --keep-allele-order --out tmp_filtered

run cat tmp_filtered.fam | cut -d ' ' -f 1-2 > keep.txt
run $1/plink2 $2 --bfile tmp_data --keep keep.txt --make-bed --out bed_filtered

cat tmp_filtered.fam | tr ' ' '\t' > tmp_filtered_tabs.fam

run diff -q tmp_filtered_tabs.fam bed_filtered.fam
run diff -q tmp_filtered.bed bed_filtered.bed
run $1/pgen_compress tmp_data.bed tmp_data.pgen 4321
run $1/plink2 $2 --bpfile tmp_data --keep keep.txt --make-bed --out pgen_filtered
run diff -q bed_filtered.fam pgen_filtered.fam
run diff -q bed_filtered.bed pgen_filtered.bed
