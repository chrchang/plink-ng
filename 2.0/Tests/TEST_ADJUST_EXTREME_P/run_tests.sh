#!/bin/bash

# --adjust-file on extreme p-values.
#
# 1. A p-value of exactly 0 (or a -log10(p) too large for ln(p) to be finite)
#    used to overflow LnPToChisq() and abort on an assertion; it is now
#    truncated to log(DBL_MIN), as the 'INF' string already was.
# 2. The Sidak corrections computed 1 - exp(c log(1-p)), which is exactly 0
#    once cp falls below ~2^{-53}: p = 1e-20 over 3 tests came out as 0
#    instead of 3e-20.

set -exo pipefail

plink2="$1/plink2 $2 $3"

# 1.
printf '#CHROM\tPOS\tID\tA1\tP\n1\t1\ta\tA\t0.01\n1\t2\tz\tA\t0\n1\t3\tc\tA\t0.3\n' > tmp_zero.txt
printf '#CHROM\tPOS\tID\tA1\tP\n1\t1\ta\tA\t0.01\n1\t2\tz\tA\tINF\n1\t3\tc\tA\t0.3\n' > tmp_inf.txt
for m in "" gc log10
do
$plink2 --adjust-file tmp_zero.txt $m --out tmp_zero > /dev/null
$plink2 --adjust-file tmp_inf.txt $m --out tmp_inf > /dev/null
diff -q tmp_zero.adjusted tmp_inf.adjusted
if grep -qiE 'nan|2147483647' tmp_zero.adjusted; then
    echo "unexpected nan/2147483647 in tmp_zero.adjusted"
    exit 1
fi
done
printf '#CHROM\tPOS\tID\tA1\tLOG10_P\n1\t1\ta\tA\t2\n1\t2\tz\tA\t1e308\n1\t3\tc\tA\t0.5\n' > tmp_huge.txt
$plink2 --adjust-file tmp_inf.txt --out tmp_inf > /dev/null
$plink2 --adjust-file tmp_huge.txt input-log10 --out tmp_huge > /dev/null
diff <(grep -w z tmp_huge.adjusted) <(grep -w z tmp_inf.adjusted)

# 2. Unadjusted p-values from 1e-5 to 1e-40 over 3 tests: SIDAK_SS and
#    SIDAK_SD must be 3p and 3p (smallest p first), within printed precision.
for e in 5 10 17 20 25 30 40
do
printf '#CHROM\tPOS\tID\tA1\tP\n1\t1\ta\tA\t1e-%s\n1\t2\tb\tA\t0.5\n1\t3\tc\tA\t0.9\n' $e > tmp_sidak.txt
$plink2 --adjust-file tmp_sidak.txt --out tmp_sidak > /dev/null
awk -F '\t' -v p=1e-$e '$2 == "a" { for (i = 7; i <= 8; ++i) { r = $i / (3 * p); if ((r < 0.9999) || (r > 1.0001)) { print "column " i " is " $i " for p = " p; exit 1 } } }' tmp_sidak.adjusted
done
