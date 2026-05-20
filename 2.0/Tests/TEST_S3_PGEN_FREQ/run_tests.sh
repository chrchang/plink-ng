#!/bin/bash

# S3 streaming test for plink2 --pfile and --bfile.
# Gated behind PLINK_TESTS_S3_SUPPORT=1; skips gracefully otherwise.

set -exo pipefail

if [ -z "$PLINK_TESTS_S3_SUPPORT" ]; then
  echo "SKIPPED: TEST_S3_PGEN_FREQ (set PLINK_TESTS_S3_SUPPORT=1 to enable)"
  exit 0
fi

export AWS_REGION="${AWS_REGION:-eu-west-1}"

S3_PREFIX="s3://ngi-igenomes/testdata/nf-core/modules/genomics/homo_sapiens/popgen/plink_simulated"

# Test 1: --pfile from S3 (pgen/psam/pvar)
$1/plink2 $2 $3 --pfile $S3_PREFIX --geno-counts --out s3_pfile
diff -q expected/pfile.gcount s3_pfile.gcount

# Test 2: --bfile from S3 (bed/bim/fam)
$1/plink2 $2 $3 --bfile $S3_PREFIX --geno-counts --out s3_bfile
diff -q expected/bfile.gcount s3_bfile.gcount
