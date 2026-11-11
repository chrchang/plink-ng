#!/bin/bash

# Regression tests for PLINK 1.9.
#
# Usage: ./run_tests.sh {plink2 build dir} {up to 2 plink2 args}
# Requires plink (1.9) to be in the system PATH.  Some tests use plink2 as a
# reference; the plink2 build dir defaults to ../../2.0/build_dynamic.

set -exo pipefail

if [[ $# -eq 0 ]]; then
    d=../../../2.0/build_dynamic
else
    # doesn't always work, but should be good enough
    DIR=$1
    if [ "${DIR:0:1}" = "/" ]; then
        d=$DIR
    else
        d=../$DIR
    fi
fi

cd TEST_19_REVIEW_FIXES
./run_tests.sh > TEST_19_REVIEW_FIXES.log
cd ..
echo "TEST_19_REVIEW_FIXES passed."

cd TEST_19_ASSOC_LD_FIXES
./run_tests.sh $d $2 $3 > TEST_19_ASSOC_LD_FIXES.log
cd ..
echo "TEST_19_ASSOC_LD_FIXES passed."

cd TEST_19_PLINK107_FIXES
./run_tests.sh > TEST_19_PLINK107_FIXES.log
cd ..
echo "TEST_19_PLINK107_FIXES passed."

cd TEST_19_CONVERSION_FIXES
./run_tests.sh $d $2 $3 > TEST_19_CONVERSION_FIXES.log
cd ..
echo "TEST_19_CONVERSION_FIXES passed."

cd TEST_19_SET_FILTER_FIXES
./run_tests.sh $d $2 $3 > TEST_19_SET_FILTER_FIXES.log
cd ..
echo "TEST_19_SET_FILTER_FIXES passed."

cd TEST_19_DFAM_MDS_FIXES
./run_tests.sh > TEST_19_DFAM_MDS_FIXES.log
cd ..
echo "TEST_19_DFAM_MDS_FIXES passed."

cd TEST_19_PERM_FIXES
./run_tests.sh > TEST_19_PERM_FIXES.log
cd ..
echo "TEST_19_PERM_FIXES passed."

cd TEST_19_BCF_CHECKS
./run_tests.sh > TEST_19_BCF_CHECKS.log
cd ..
echo "TEST_19_BCF_CHECKS passed."

cd TEST_19_LONG_CONTIG
./run_tests.sh > TEST_19_LONG_CONTIG.log
cd ..
echo "TEST_19_LONG_CONTIG passed."

echo "All tests passed."
