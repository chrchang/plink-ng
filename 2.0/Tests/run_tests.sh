#!/bin/bash

# Usage: ./run_tests.sh {plink2 + pgen_compress build dir}
# Requires plink to be in the system PATH.

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

# http://stackoverflow.com/questions/5195607/checking-bash-exit-status-of-several-commands-efficiently
function run {
    "$@"
    local status=$?
    if [ $status -ne 0 ]; then
	exit $status
    fi
    return $status
}

run cd UNIT_TEST_EXTRACT_CHR
run ./run_tests.sh $d > UNIT_TEST_EXTRACT_CHR.log
run cd ..
echo "UNIT_TEST_EXTRACT_CHR passed."

run cd UNIT_TEST_MAF_FILTER
run ./run_tests.sh $d > UNIT_TEST_MAF_FILTER.log
run cd ..
echo "UNIT_TEST_MAF_FILTER passed."

run cd UNIT_TEST_PGEN_FREQ
run ./run_tests.sh $d > UNIT_TEST_PGEN_FREQ.log
run cd ..
echo "UNIT_TEST_PGEN_FREQ passed."

run cd UNIT_TEST_PHASED_VCF
run ./run_tests.sh $d > UNIT_TEST_PHASED_VCF.log
run cd ..
echo "UNIT_TEST_PHASED_VCF passed."

run cd UNIT_TEST_SAMPLE_SUBSET
run ./run_tests.sh $d > UNIT_TEST_SAMPLE_SUBSET.log
run cd ..
echo "UNIT_TEST_SAMPLE_SUBSET passed."

echo "All tests passed."
