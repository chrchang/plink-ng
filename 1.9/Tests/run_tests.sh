#!/bin/bash

# Regression tests for PLINK 1.9.
#
# Usage: ./run_tests.sh {plink 1.9 dir} {plink2 build dir}
#   {up to 2 plink2 args}
# plink (1.9) is taken from {plink 1.9 dir}, or from the system PATH if that
# is empty or omitted.  Some tests use plink2 as a reference; the plink2 build
# dir defaults to ../../2.0/build_dynamic.  Tests that compare against PLINK
# 1.07 take it from $PLINK107, defaulting to ../../1.07/plink, and skip those
# checks if it isn't there.

set -exo pipefail

# Relative paths are given relative to this directory, while the tests run
# one level down.
subdir_path() {
    if [ "${1:0:1}" = "/" ]; then
        echo "$1"
    else
        echo "../$1"
    fi
}

p19=
if [[ -n "${1:-}" ]]; then
    p19=$(subdir_path "$1")
fi
if [[ -z "${2:-}" ]]; then
    d=../../../2.0/build_dynamic
else
    d=$(subdir_path "$2")
fi

cd TEST_19_REVIEW_FIXES
./run_tests.sh "$p19" > TEST_19_REVIEW_FIXES.log
cd ..
echo "TEST_19_REVIEW_FIXES passed."

cd TEST_19_ASSOC_LD_FIXES
./run_tests.sh "$p19" $d $3 $4 > TEST_19_ASSOC_LD_FIXES.log
cd ..
echo "TEST_19_ASSOC_LD_FIXES passed."

cd TEST_19_PLINK107_FIXES
./run_tests.sh "$p19" > TEST_19_PLINK107_FIXES.log
cd ..
echo "TEST_19_PLINK107_FIXES passed."

cd TEST_19_CONVERSION_FIXES
./run_tests.sh "$p19" $d $3 $4 > TEST_19_CONVERSION_FIXES.log
cd ..
echo "TEST_19_CONVERSION_FIXES passed."

cd TEST_19_SET_FILTER_FIXES
./run_tests.sh "$p19" $d $3 $4 > TEST_19_SET_FILTER_FIXES.log
cd ..
echo "TEST_19_SET_FILTER_FIXES passed."

cd TEST_19_DFAM_MDS_FIXES
./run_tests.sh "$p19" > TEST_19_DFAM_MDS_FIXES.log
cd ..
echo "TEST_19_DFAM_MDS_FIXES passed."

cd TEST_19_PERM_FIXES
./run_tests.sh "$p19" > TEST_19_PERM_FIXES.log
cd ..
echo "TEST_19_PERM_FIXES passed."

cd TEST_19_BCF_CHECKS
./run_tests.sh "$p19" > TEST_19_BCF_CHECKS.log
cd ..
echo "TEST_19_BCF_CHECKS passed."

cd TEST_19_LONG_CONTIG
./run_tests.sh "$p19" > TEST_19_LONG_CONTIG.log
cd ..
echo "TEST_19_LONG_CONTIG passed."

cd TEST_19_FORMATS
./run_tests.sh "$p19" $d $3 $4 > TEST_19_FORMATS.log
cd ..
echo "TEST_19_FORMATS passed."

cd TEST_19_MERGE
./run_tests.sh "$p19" $d $3 $4 > TEST_19_MERGE.log
cd ..
echo "TEST_19_MERGE passed."

cd TEST_19_QC_STATS
./run_tests.sh "$p19" $d $3 $4 > TEST_19_QC_STATS.log
cd ..
echo "TEST_19_QC_STATS passed."

cd TEST_19_ASSOC_BASIC
./run_tests.sh "$p19" $d $3 $4 > TEST_19_ASSOC_BASIC.log
cd ..
echo "TEST_19_ASSOC_BASIC passed."

cd TEST_19_TDT_EPI_CLUSTER
./run_tests.sh "$p19" > TEST_19_TDT_EPI_CLUSTER.log
cd ..
echo "TEST_19_TDT_EPI_CLUSTER passed."

cd TEST_19_PLINK107_COMPARE
./run_tests.sh "$p19" > TEST_19_PLINK107_COMPARE.log
cd ..
echo "TEST_19_PLINK107_COMPARE passed."

cd TEST_19_FILTERS_VS_PLINK2
./run_tests.sh "$p19" $d $3 $4 > TEST_19_FILTERS_VS_PLINK2.log
cd ..
echo "TEST_19_FILTERS_VS_PLINK2 passed."

cd TEST_19_EDIT_IO_VS_PLINK2
./run_tests.sh "$p19" $d $3 $4 > TEST_19_EDIT_IO_VS_PLINK2.log
cd ..
echo "TEST_19_EDIT_IO_VS_PLINK2 passed."

cd TEST_19_PHENO_CHRSET_VS_PLINK2
./run_tests.sh "$p19" $d $3 $4 > TEST_19_PHENO_CHRSET_VS_PLINK2.log
cd ..
echo "TEST_19_PHENO_CHRSET_VS_PLINK2 passed."

cd TEST_19_LD_FAMILY_VS_PLINK2
./run_tests.sh "$p19" $d $3 $4 > TEST_19_LD_FAMILY_VS_PLINK2.log
cd ..
echo "TEST_19_LD_FAMILY_VS_PLINK2 passed."

cd TEST_19_PLINK107_OPTIONS
./run_tests.sh "$p19" > TEST_19_PLINK107_OPTIONS.log
cd ..
echo "TEST_19_PLINK107_OPTIONS passed."

cd TEST_19_MISC_OPTIONS
./run_tests.sh "$p19" > TEST_19_MISC_OPTIONS.log
cd ..
echo "TEST_19_MISC_OPTIONS passed."

cd TEST_19_RECALC
./run_tests.sh "$p19" > TEST_19_RECALC.log
cd ..
echo "TEST_19_RECALC passed."

cd TEST_19_RUNTIME_FLAGS
./run_tests.sh "$p19" > TEST_19_RUNTIME_FLAGS.log
cd ..
echo "TEST_19_RUNTIME_FLAGS passed."

echo "All tests passed."
