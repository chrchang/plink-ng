#!/bin/bash

set -exo pipefail

# For files that can't be re-imported, we compare output against
# manually-sanity-checked ground truth files.

$1/plink2 $2 $3 --pfile all_ploidy vzs --snps-only --export phylip --out plink2_test
diff -q plink2_test.phy plink2_test.phy.want
$1/plink2 $2 $3 --pfile all_ploidy vzs --chr PAR1,X,PAR2 --keep-females --snps-only --export phylip-phased --out plink2_test.phased
diff -q plink2_test.phased.phy plink2_test.phased.phy.want
