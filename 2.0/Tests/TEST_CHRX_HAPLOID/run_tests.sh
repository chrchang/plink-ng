#!/bin/bash

set -exo pipefail

# chrX fileset where 20% of the male calls are het haploid.  Those calls are
# malformed input, and every per-variant count is expected to treat them as
# missing.
python3 make_chrx_vcf.py 1 200 300

$1/plink2 $2 $3 --vcf chrx.vcf --update-sex chrx_sex.txt --make-pgen --out plink2_chrx

$1/plink2 $2 $3 --pfile plink2_chrx --freq counts --out plink2_acount
$1/plink2 $2 $3 --pfile plink2_chrx --geno-counts --out plink2_gcount
python3 check_counts.py plink2_acount.acount plink2_gcount.gcount

# --set-invalid-haploid-missing removes the het haploid calls up front, so it
# must not change any of these counts.
$1/plink2 $2 $3 --pfile plink2_chrx --set-invalid-haploid-missing --make-pgen --out plink2_chrx_clean
$1/plink2 $2 $3 --pfile plink2_chrx_clean --freq counts --out plink2_acount_clean
$1/plink2 $2 $3 --pfile plink2_chrx_clean --geno-counts --out plink2_gcount_clean
diff -q plink2_acount.acount plink2_acount_clean.acount
diff -q plink2_gcount.gcount plink2_gcount_clean.gcount

# Same check on a sample subset, to cover the sample_ct != raw_sample_ct code
# path.
head -n 100 plink2_chrx.psam | tail -n +2 | cut -f 1 > subset.id
$1/plink2 $2 $3 --pfile plink2_chrx --keep subset.id --freq counts --out plink2_acount_sub
$1/plink2 $2 $3 --pfile plink2_chrx --keep subset.id --geno-counts --out plink2_gcount_sub
python3 check_counts.py plink2_acount_sub.acount plink2_gcount_sub.gcount

# And single-threaded, since the counts are accumulated per thread.
$1/plink2 $2 --threads 1 --pfile plink2_chrx --freq counts --out plink2_acount_st
diff -q plink2_acount.acount plink2_acount_st.acount
