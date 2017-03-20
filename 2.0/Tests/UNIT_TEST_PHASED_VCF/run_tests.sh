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


# 1kg_phase3_chr21_start.vcf.gz is the header (253 lines) + first 600
# biallelic lines from the official data release (grep -v "MULTI_ALLELIC")
# (was many more lines, but the test was too slow.  should manually test with
# full chr21 file before every merge with develop.)
# todo: support multiallelic variants
run plink --vcf 1kg_phase3_chr21_start.vcf.gz --out plink1_data
run $1/plink2 --vcf 1kg_phase3_chr21_start.vcf.gz --out plink2_data
run $1/plink2 --pfile plink2_data --export vcf --out plink2_data

rm -f 1kg_phase3_chr21_start.vcf

run gunzip -k 1kg_phase3_chr21_start.vcf.gz

run cat 1kg_phase3_chr21_start.vcf | sed '/^#/ d' > 1kg_noheader.txt

# INFO column is now preserved.
run cat plink2_data.vcf | sed '/^#/ d' > plink2_noheader.txt
run diff -q 1kg_noheader.txt plink2_noheader.txt

# 10th column is the first genotype column, so this removes samples 492-569 and
# sample 667
run cat 1kg_noheader.txt | cut -f 1-500,579-675,677- > 1kg_noheader_pruned.txt

run cat plink1_data.fam | head -n 569 | tail -n 78 | cut -f 1-2 > remove.txt
run cat plink1_data.fam | head -n 667 | tail -n 1 | cut -f 1-2 >> remove.txt
run $1/plink2 --pfile plink2_data --remove remove.txt --export vcf --out plink2_pruned

# there's a tiny chance that fileDate will change between this and the next VCF
# export command.  it's the second header line, so remove the top two header
# lines before comparing.
run cat plink2_pruned.vcf | tail -n +3 > plink2_pruned_nodate.txt

run cat plink2_pruned.vcf | sed '/^#/ d' > plink2_noheader_pruned.txt
run diff -q 1kg_noheader_pruned.txt plink2_noheader_pruned.txt
run $1/plink2 --pfile plink2_data --remove remove.txt --make-pgen --out plink2_pruned

# this exercises different subsetting code paths
run $1/plink2 --pfile plink2_pruned --export vcf --out plink2_pruned2

run cat plink2_pruned2.vcf | tail -n +3 > plink2_pruned2_nodate.txt
run diff -q plink2_pruned_nodate.txt plink2_pruned2_nodate.txt

run plink --bfile plink1_data --thin 0.1 --thin-indiv 0.1 --make-bed --out sparse

run cat sparse.bim | cut -f 2 > extract.txt
run $1/plink2 --pfile plink2_data --keep sparse.fam --extract extract.txt --export vcf --out plink2_sparse

# test phased genotype sorting code path
run $1/plink2 --pfile plink2_data --keep sparse.fam --indiv-sort f sparse.fam --extract extract.txt --make-pgen --out plink2_sparse2
run $1/plink2 --vcf plink2_sparse.vcf --make-pgen --out plink2_sparse
run diff -q plink2_sparse.pgen plink2_sparse2.pgen
run diff -q plink2_sparse.pvar plink2_sparse2.pvar
run diff -q plink2_sparse.psam plink2_sparse2.psam
