#!/bin/bash

# This generates a few datasets for our tests to run on.  (It is not efficient
# to just store them in GitHub.)  It is assumed that the PLINK 1.9 build to be
# tested is named/symlinked as "plink2", and a PLINK 1.07 build is
# named/symlinked as "plink1".

# When the number of samples isn't divisible by 4, special-case logic is
# sometimes needed to handle the last data byte for each variant.  This has
# been a source of bugs, so it's important to have a test case for each sample
# count remainder mod 4.

# 385 samples x 1089 markers should be a large enough dataset to test generic
# advancement loops, and the lower the runtime of the test suite, the more
# frequently it can be run.

# This procedure has a small chance of generating a fileset that causes a test
# to fail due to floating point error.  It's the responsibility of the tester
# to figure out when that has happened (in which case they can just rerun this
# script), vs. when an actual bug is involved.

plink19 --silent --dummy 513 1423 0.02 --out dummy_cc1
plink19 --silent --dummy 512 1234 0.04 --out dummy_cc2
plink19 --silent --dummy 387 1112 0.03 scalar-pheno --out dummy1
plink19 --silent --dummy 478 1111 0.05 scalar-pheno --out dummy2
plink19 --silent --dummy 3 9999 0.01 --out trio_tmp
rm trio_tmp.fam
echo "fam1 dad 0 0 1 2" > trio_tmp.fam
echo "fam1 mom 0 0 2 1" >> trio_tmp.fam
echo "fam1 son dad mom 1 1" >> trio_tmp.fam
plink19 --silent --bfile trio_tmp --geno 0.6 --make-bed --out trio
rm trio_tmp*
rm set.txt
echo "1 10 20 set1" > set.txt
echo "1 15 40 set2" >> set.txt
rm score.txt
echo "snp1 A 0.1" > score.txt
echo "snp3 A 0.04" >> score.txt
