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
# Verify BCF round-trip.
$1/plink2 $2 $3 --pfile plink2_data --export bcf --out plink2_roundtrip
$1/plink2 $2 $3 --bcf plink2_roundtrip.bcf --out plink2_roundtrip
diff -q plink2_data.pgen plink2_roundtrip.pgen

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
# reduce to 3 PCs since that's the most that --pca approx can handle with 29
# variants
plink --bfile plink1_data --maf 0.02 --pca 3 header tabs var-wts --out plink1_pca
$1/plink2 $2 $3 --bfile plink1_data --maf 0.02 --pca 3 biallelic-var-wts --out plink2_pca
python3 pca_compare.py -1 plink1_pca -2 plink2_pca -t 0.000002

# while we're at it, test all --read-freq input formats.
# (todo: add multiallelic --read-freq test when --score gets better
# multiallelic-variant support.)
plink --bfile plink1_data --freq
$1/plink2 $2 $3 --bfile plink1_data --read-freq plink.frq --maf 0.02 --pca 3 biallelic-var-wts --out plink2_pca_rfreq
# .frq only has 4-digit precision
python3 pca_compare.py -1 plink2_pca -2 plink2_pca_rfreq -t 0.0002
plink --bfile plink1_data --freqx
$1/plink2 $2 $3 --bfile plink1_data --read-freq plink.frqx --maf 0.02 --pca 3 biallelic-var-wts --out plink2_pca_rfreq
diff -q plink2_pca.eigenvec plink2_pca_rfreq.eigenvec
$1/plink2 --bfile plink1_data --freq
$1/plink2 $2 $3 --bfile plink1_data --read-freq plink2.afreq --maf 0.02 --pca 3 biallelic-var-wts --out plink2_pca_rfreq
python3 pca_compare.py -1 plink2_pca -2 plink2_pca_rfreq -t 0.000002
$1/plink2 --bfile plink1_data --geno-counts
$1/plink2 $2 $3 --bfile plink1_data --read-freq plink2.gcount --maf 0.02 --pca 3 biallelic-var-wts --out plink2_pca_rfreq
diff -q plink2_pca.eigenvec plink2_pca_rfreq.eigenvec

# note that this run depends on the random seed.
$1/plink2 $2 $3 --bfile plink1_data --maf 0.02 --pca 3 approx biallelic-var-wts --out plink2_pca_approx
python3 pca_compare.py -1 plink1_pca -2 plink2_pca_approx -t 0.0064

# Test --glm.
# Generate random binary and quantitative phenotypes for 1kg_phase3_chr21
# samples, and verify regression results are consistent with plink 1.9 within a
# tolerance (tbd: see if we can generate new phenotypes each time, or if it's
# better to just hardcode one pair)
#
# conditions:
# - no covariates
# - no covariates, 'genotypic'
# - top PCs as covariates
# - top PCs, 'genotypic'
# - top PCs, --tests
#
# TODO: Also test --condition/--condition-list, local-covar.
# TODO: TEST_MULTIALLELIC_VCF fixture, etc.

$1/plink2 --dummy 2504 1 --seed 1 --out pheno_cc
$1/plink2 --dummy 2504 1 0.01 0.01 scalar-pheno --seed 1 --out pheno_qt
cat pheno_cc.psam | tail -n +2 | cut -f 3 > pheno_cc_col.txt
cat pheno_qt.psam | tail -n +2 | cut -f 3 > pheno_qt_col.txt
cat plink1_data.fam | cut -f 1-2 > sample_ids.txt
paste sample_ids.txt pheno_cc_col.txt > pheno_cc.txt
paste sample_ids.txt pheno_qt_col.txt > pheno_qt.txt

plink --bfile plink1_data --maf 0.02 --pheno pheno_cc.txt --logistic --allow-no-sex --out plink1_glm
$1/plink2 $2 $3 --bfile plink1_data --maf 0.02 --pheno pheno_cc.txt --glm allow-no-covars no-firth --out plink2_glm
python3 glm_compare.py -1 plink1_glm.assoc.logistic -2 plink2_glm.PHENO1.glm.logistic -t 0.1
plink --bfile plink1_data --maf 0.02 --pheno pheno_qt.txt --linear --allow-no-sex --out plink1_glm
$1/plink2 $2 $3 --bfile plink1_data --maf 0.02 --pheno pheno_qt.txt --glm allow-no-covars --out plink2_glm
python3 glm_compare.py -1 plink1_glm.assoc.linear -2 plink2_glm.PHENO1.glm.linear -t 0.1

plink --bfile plink1_data --maf 0.02 --pheno pheno_cc.txt --logistic genotypic --allow-no-sex --out plink1_glm
$1/plink2 $2 $3 --bfile plink1_data --maf 0.02 --pheno pheno_cc.txt --glm allow-no-covars no-firth genotypic --out plink2_glm
python3 glm_compare.py -1 plink1_glm.assoc.logistic -2 plink2_glm.PHENO1.glm.logistic -t 0.1
plink --bfile plink1_data --maf 0.02 --pheno pheno_qt.txt --linear genotypic --allow-no-sex --out plink1_glm
$1/plink2 $2 $3 --bfile plink1_data --maf 0.02 --pheno pheno_qt.txt --glm allow-no-covars genotypic --out plink2_glm
python3 glm_compare.py -1 plink1_glm.assoc.linear -2 plink2_glm.PHENO1.glm.linear -t 0.1

plink --bfile plink1_data --maf 0.02 --pheno pheno_cc.txt --logistic --covar plink1_pca.eigenvec --allow-no-sex --out plink1_glm
$1/plink2 $2 $3 --bfile plink1_data --maf 0.02 --pheno pheno_cc.txt --glm no-firth --covar plink2_pca.eigenvec --out plink2_glm
python3 glm_compare.py -1 plink1_glm.assoc.logistic -2 plink2_glm.PHENO1.glm.logistic -t 0.1
plink --bfile plink1_data --maf 0.02 --pheno pheno_qt.txt --linear --covar plink1_pca.eigenvec --allow-no-sex --out plink1_glm
$1/plink2 $2 $3 --bfile plink1_data --maf 0.02 --pheno pheno_qt.txt --glm --covar plink2_pca.eigenvec --out plink2_glm
python3 glm_compare.py -1 plink1_glm.assoc.linear -2 plink2_glm.PHENO1.glm.linear -t 0.1

plink --bfile plink1_data --maf 0.02 --pheno pheno_cc.txt --logistic genotypic --covar plink1_pca.eigenvec --allow-no-sex --out plink1_glm
$1/plink2 $2 $3 --bfile plink1_data --maf 0.02 --pheno pheno_cc.txt --glm no-firth genotypic --covar plink2_pca.eigenvec --out plink2_glm
python3 glm_compare.py -1 plink1_glm.assoc.logistic -2 plink2_glm.PHENO1.glm.logistic -t 0.1
plink --bfile plink1_data --maf 0.02 --pheno pheno_qt.txt --linear genotypic --covar plink1_pca.eigenvec --allow-no-sex --out plink1_glm
$1/plink2 $2 $3 --bfile plink1_data --maf 0.02 --pheno pheno_qt.txt --glm genotypic --covar plink2_pca.eigenvec --out plink2_glm
python3 glm_compare.py -1 plink1_glm.assoc.linear -2 plink2_glm.PHENO1.glm.linear -t 0.1

plink --bfile plink1_data --maf 0.02 --pheno pheno_cc.txt --logistic --covar plink1_pca.eigenvec --tests 1,3,4 --allow-no-sex --out plink1_glm
$1/plink2 $2 $3 --bfile plink1_data --maf 0.02 --pheno pheno_cc.txt --glm no-firth --covar plink2_pca.eigenvec --tests 1,3,4 --out plink2_glm
python3 glm_compare.py -1 plink1_glm.assoc.logistic -2 plink2_glm.PHENO1.glm.logistic -t 0.1
# Note that 5 is "out of range" here.  This was handled incorrectly by ALL
# versions of plink before this test!!
plink --bfile plink1_data --maf 0.02 --pheno pheno_qt.txt --linear --covar plink1_pca.eigenvec --tests 1,3,5 --allow-no-sex --out plink1_glm
$1/plink2 $2 $3 --bfile plink1_data --maf 0.02 --pheno pheno_qt.txt --glm --covar plink2_pca.eigenvec --tests 1,3,5 --out plink2_glm
python3 glm_compare.py -1 plink1_glm.assoc.linear -2 plink2_glm.PHENO1.glm.linear -t 0.1
