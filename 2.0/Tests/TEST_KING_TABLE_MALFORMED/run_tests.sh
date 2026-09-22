#!/bin/bash

# A .kin0 data line that ends after the first sample ID must be rejected as
# malformed by both --king-table-subset and --king-cutoff-table.
#
# Both readers skipped the whitespace after the first ID and handed whatever
# came next to XidRead(), which asserts that it is not at end of line.  So the
# truncated line aborted plink2 instead of producing an error.

set -exo pipefail

$1/plink2 $2 $3 --dummy 6 200 --seed 1 --make-pgen --out tmp_data
$1/plink2 $2 $3 --pfile tmp_data --make-king-table --king-table-filter -9 --out tmp_data

# Sanity check: the unmodified table is accepted by both readers.
$1/plink2 $2 $3 --pfile tmp_data --king-table-subset tmp_data.kin0 --make-king-table --out plink2_ok
$1/plink2 $2 $3 --pfile tmp_data --king-cutoff-table tmp_data.kin0 0.2 --make-just-psam --out plink2_ok

# Line 3 keeps only its first sample ID.
awk 'NR == 3 { print $1; next } { print }' tmp_data.kin0 > tmp_short.kin0

if $1/plink2 $2 $3 --pfile tmp_data --king-table-subset tmp_short.kin0 --make-king-table --out plink2_subset 2> tmp_err_subset.txt; then
    echo "expected --king-table-subset to reject tmp_short.kin0"
    exit 1
fi
grep -q "Line 3 of --king-table-subset file has fewer tokens than expected" tmp_err_subset.txt

if $1/plink2 $2 $3 --pfile tmp_data --king-cutoff-table tmp_short.kin0 0.2 --make-just-psam --out plink2_cutoff 2> tmp_err_cutoff.txt; then
    echo "expected --king-cutoff-table to reject tmp_short.kin0"
    exit 1
fi
grep -q "Fewer tokens than expected on line 3 of tmp_short.kin0" tmp_err_cutoff.txt
