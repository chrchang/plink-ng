#!/bin/bash

# Usage: ./run_tests.sh {plink2 + pgen_compress build dir}
#   {up to 2 args, e.g. --randmem, "--threads 1"}
# Requires plink to be in the system PATH.

set -exo pipefail

if [[ $# -eq 0 ]]; then
    d=../../build_dynamic
else
    # doesn't always work, but should be good enough
    DIR=$1
    if [ "${DIR:0:1}" = "/" ]; then
        d=$DIR
    else
        d=../$DIR
    fi
fi

cd TEST_EXTRACT_CHR
./run_tests.sh $d $2 $3 > TEST_EXTRACT_CHR.log
cd ..
echo "TEST_EXTRACT_CHR passed."

cd TEST_MAF_FILTER
./run_tests.sh $d $2 $3 > TEST_MAF_FILTER.log
cd ..
echo "TEST_MAF_FILTER passed."

cd TEST_PGEN_FREQ
./run_tests.sh $d $2 $3 > TEST_PGEN_FREQ.log
cd ..
echo "TEST_PGEN_FREQ passed."

cd TEST_PHASED_VCF
./run_tests.sh $d $2 $3 > TEST_PHASED_VCF.log
cd ..
echo "TEST_PHASED_VCF passed."

cd TEST_SAMPLE_SUBSET
./run_tests.sh $d $2 $3 > TEST_SAMPLE_SUBSET.log
cd ..
echo "TEST_SAMPLE_SUBSET passed."

cd TEST_DOSAGE_ROUND_TRIP
./run_tests.sh $d $2 $3 > TEST_DOSAGE_ROUND_TRIP.log
cd ..
echo "TEST_DOSAGE_ROUND_TRIP passed."

echo "All tests passed."
