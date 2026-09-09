#!/bin/bash

# --write-var-ranges, checked against PLINK 1.9.

set -exo pipefail

plink --simulate simulate.txt --make-bed --out tmp_data > /dev/null

# Block counts from 1 up to the variant count, so both extremes are covered.
for n in 1 2 3 7 10 97
do
    plink --bfile tmp_data --write-var-ranges $n --allow-no-sex --out plink19
    $1/plink2 $2 $3 --bfile tmp_data --write-var-ranges $n --out plink2_vr
    # PLINK 1.9's header is "FIRST LAST"; plink2 prefixes it with '#'.
    diff -q <(tail -n +2 plink19.var.ranges) <(tail -n +2 plink2_vr.var.ranges)
    test "$(tail -n +2 plink2_vr.var.ranges | wc -l)" -eq $n
done

# One block per variant.
variant_ct=$(wc -l < tmp_data.bim)
plink --bfile tmp_data --write-var-ranges $variant_ct --allow-no-sex --out plink19_all
$1/plink2 $2 $3 --bfile tmp_data --write-var-ranges $variant_ct --out plink2_all
diff -q <(tail -n +2 plink19_all.var.ranges) <(tail -n +2 plink2_all.var.ranges)

# Zstd output round-trips.
$1/plink2 $2 $3 --bfile tmp_data --write-var-ranges 10 zs --out plink2_zs
$1/plink2 $2 $3 --zst-decompress plink2_zs.var.ranges.zst > plink2_zs.var.ranges
$1/plink2 $2 $3 --bfile tmp_data --write-var-ranges 10 --out plink2_plain
diff -q plink2_zs.var.ranges plink2_plain.var.ranges

# More blocks than variants is an error.
if $1/plink2 $2 $3 --bfile tmp_data --write-var-ranges $((variant_ct + 1)) --out plink2_bad 2> tmp_err.txt; then
    echo "expected --write-var-ranges to reject a block count above the variant count"
    exit 1
fi
grep -q "exceeds the number of variants" tmp_err.txt
