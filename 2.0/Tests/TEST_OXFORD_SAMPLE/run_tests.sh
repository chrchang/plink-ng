#!/bin/bash

# Numeric categorical (discrete, type 'D') columns in an Oxford .sample file:
# values up to 2^31 - 1 are accepted, and larger ones are rejected without
# undefined behavior.
#
# OxSampleToPsam() converted the parsed double to int32_t before checking
# its range, so a value such as 4294967296 was an out-of-range
# floating-point-to-integer conversion.  The error message was right, but the
# undefined-behavior sanitizer flags the conversion; the check below for
# "runtime error" catches it in a sanitizer build.

set -exo pipefail

$1/plink2 $2 $3 --dummy 4 5 --seed 1 --export bgen-1.2 --out tmp_data

# Append a categorical column, with value $1 for the last sample.
make_sample() {
    awk -v v="$1" 'NR == 1 { print $0, "CAT"; next } NR == 2 { print $0, "D"; next } { print $0, (NR == 6)? v : NR - 2 }' tmp_data.sample
}

make_sample 2147483647 > tmp_ok.sample
$1/plink2 $2 $3 --bgen tmp_data.bgen ref-first --sample tmp_ok.sample --make-just-psam --out plink2_ok 2> tmp_err_ok.txt
if grep -q 'runtime error' tmp_err_ok.txt; then
    exit 1
fi
test "$(awk '$1 == "per3" { print $NF }' plink2_ok.psam)" = C2147483647

for v in 2147483648 4294967296 1e300; do
    make_sample $v > tmp_bad.sample
    if $1/plink2 $2 $3 --bgen tmp_data.bgen ref-first --sample tmp_bad.sample --make-just-psam --out plink2_bad 2> tmp_err_bad.txt; then
        echo "expected --sample to reject categorical value $v"
        exit 1
    fi
    grep -q "Invalid categorical phenotype value '$v' on line 6, column 6" tmp_err_bad.txt
    if grep -q 'runtime error' tmp_err_bad.txt; then
        exit 1
    fi
done
