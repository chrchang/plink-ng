#!/bin/bash

# Generate gold-standard expected output files for TEST_S3_PGEN_FREQ.
#
# Downloads the test data from the public nf-core/igenomes S3 bucket and runs
# plink2 locally to produce the reference .gcount files committed under
# expected/.
#
# Prerequisites:
#   - AWS CLI installed (aws s3 cp --no-sign-request)
#   - plink2 binary (pass build dir as $1, or defaults to ../../build_dynamic)
#
# Usage:
#   ./generate_expected.sh [plink2_build_dir]

set -exo pipefail

if [[ $# -eq 0 ]]; then
    PLINK2=../../build_dynamic/plink2
else
    PLINK2=$1/plink2
fi

export AWS_DEFAULT_REGION="${AWS_DEFAULT_REGION:-eu-west-1}"

S3_BUCKET="s3://ngi-igenomes/testdata/nf-core/modules/genomics/homo_sapiens/popgen"
LOCAL_DIR="tmp_local_data"

mkdir -p "$LOCAL_DIR"
mkdir -p expected

# Download all required files
for ext in pgen psam pvar bed bim fam; do
  aws s3 cp --no-sign-request "${S3_BUCKET}/plink_simulated.${ext}" "${LOCAL_DIR}/plink_simulated.${ext}"
done

# Generate expected output from pfile (pgen/psam/pvar)
$PLINK2 --pfile "${LOCAL_DIR}/plink_simulated" --geno-counts --out "${LOCAL_DIR}/pfile_result"
cp "${LOCAL_DIR}/pfile_result.gcount" expected/pfile.gcount

# Generate expected output from bfile (bed/bim/fam)
$PLINK2 --bfile "${LOCAL_DIR}/plink_simulated" --geno-counts --out "${LOCAL_DIR}/bfile_result"
cp "${LOCAL_DIR}/bfile_result.gcount" expected/bfile.gcount

# Clean up downloaded data
rm -rf "$LOCAL_DIR"

echo "Gold-standard files generated in expected/"
echo "  expected/pfile.gcount"
echo "  expected/bfile.gcount"
