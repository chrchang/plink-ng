#!/bin/bash

set -exo pipefail

# 1kg_phase3_chr21_start.vcf.gz is the header (253 lines) + first 600
# biallelic lines from the official data release (grep -v "MULTI_ALLELIC").
# This previously included many more lines, but the test was too slow.
# Ideally, we want to test the full chr21 file, or even better, chrX/Y/MT,
# before significant releases.
# todo: support multiallelic variants
plink --vcf 1kg_phase3_chr21_start.vcf.gz --out plink1_data
$1/plink2 $2 $3 --vcf 1kg_phase3_chr21_start.vcf.gz --double-id --out plink2_data
$1/plink2 $2 $3 --pfile plink2_data --export vcf --out plink2_data

rm -f 1kg_phase3_chr21_start.vcf

gunzip -k 1kg_phase3_chr21_start.vcf.gz

cat 1kg_phase3_chr21_start.vcf | sed '/^#/ d' > 1kg_noheader.txt

# INFO column is now preserved.
cat plink2_data.vcf | sed '/^#/ d' > plink2_noheader.txt
diff -q 1kg_noheader.txt plink2_noheader.txt

# 10th column is the first genotype column, so this removes samples 492-569 and
# sample 667
cat 1kg_noheader.txt | cut -f 1-500,579-675,677- > 1kg_noheader_pruned.txt

cat plink1_data.fam | head -n 569 | tail -n 78 | cut -f 1-2 > remove.txt
cat plink1_data.fam | head -n 667 | tail -n 1 | cut -f 1-2 >> remove.txt
$1/plink2 $2 $3 --pfile plink2_data --remove remove.txt --export vcf --out plink2_pruned

# there's a tiny chance that fileDate will change between this and the next VCF
# export command.  it's the second header line, so remove the top two header
# lines before comparing.
cat plink2_pruned.vcf | tail -n +3 > plink2_pruned_nodate.txt

cat plink2_pruned.vcf | sed '/^#/ d' > plink2_noheader_pruned.txt
diff -q 1kg_noheader_pruned.txt plink2_noheader_pruned.txt
$1/plink2 $2 $3 --pfile plink2_data --remove remove.txt --make-pgen --out plink2_pruned

# this exercises different subsetting code paths
$1/plink2 $2 $3 --pfile plink2_pruned --export vcf --out plink2_pruned2

cat plink2_pruned2.vcf | tail -n +3 > plink2_pruned2_nodate.txt
diff -q plink2_pruned_nodate.txt plink2_pruned2_nodate.txt

# N.B. to reproduce a failure, it may be necessary to reuse the
# sparse.bim/sparse.fam generated here (or reuse the random seed in
# sparse.log).
plink --bfile plink1_data --thin 0.1 --thin-indiv 0.1 --make-bed --out sparse

cat sparse.bim | cut -f 2 > extract.txt
$1/plink2 $2 $3 --pfile plink2_data --keep sparse.fam --extract extract.txt --export vcf --out plink2_sparse

# test phased genotype sorting code path
$1/plink2 $2 $3 --pfile plink2_data --keep sparse.fam --indiv-sort f sparse.fam --extract extract.txt --make-pgen --out plink2_sparse2
$1/plink2 $2 $3 --vcf plink2_sparse.vcf --id-delim --make-pgen --out plink2_sparse
diff -q plink2_sparse.pgen plink2_sparse2.pgen
diff -q plink2_sparse.pvar plink2_sparse2.pvar
diff -q plink2_sparse.psam plink2_sparse2.psam

# test GRM/PCA, with tolerance for floating point errors and sign differences
# (work in progress, comparison routine not written yet)
plink --bfile plink1_data --maf 0.02 --pca 5 header tabs var-wts --out plink1_pca
$1/plink2 $2 $3 --bfile plink1_data --maf 0.02 --pca 5 var-wts --out plink2_pca
python pca_compare.py -1 plink1_pca -2 plink2_pca -t 0.000002

# note that this run depends on the random seed.
$1/plink2 $2 $3 --bfile plink1_data --maf 0.02 --pca 5 approx var-wts --out plink2_pca_approx
python pca_compare.py -1 plink1_pca -2 plink2_pca_approx -t 0.1
