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


run plink --simulate simulate.txt --out tmp_data

run cat tmp_data.bim | head -n 147 > tmp_data1.bim
run cat tmp_data.bim | tail -n 150 | sed 's/^1/2/g' > tmp_data24.bim
run cat tmp_data24.bim | head -n 100 > tmp_data2.bim
run cat tmp_data24.bim | tail -n 50 | sed 's/^2/4/g' > tmp_data4.bim

run cat tmp_data1.bim tmp_data2.bim tmp_data4.bim > tmp_data.bim
run $1/pgen_compress tmp_data.bed tmp_data.pgen 2000

for c in 1 2 4
do
    run plink --bfile tmp_data --chr $c --make-bed --out plink_c$c
    run $1/plink2 --bfile tmp_data --chr $c --make-bed --out plink2_c$c
    run diff -q plink_c$c.bim plink2_c$c.bim
    run diff -q plink_c$c.bed plink2_c$c.bed
    run $1/plink2 --bpfile tmp_data --chr $c --make-bed --out plink2x_c$c
    run diff -q plink2x_c$c.bed plink2_c$c.bed
    run diff -q plink2x_c$c.bim plink2_c$c.bim
done
