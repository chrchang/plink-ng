#!/bin/bash

# --twolocus, checked against PLINK 1.9.
#
# PLINK 1.9 writes a fixed-width count matrix with marginal totals per group;
# plink2 writes one row per (group, genotype, genotype) cell.  compare.awk
# reads both and matches the cells.

set -exo pipefail

BUILD=$1
EXTRA1=$2
EXTRA2=$3

awk -f make_vcf.awk > tmp_h.vcf
$BUILD/plink2 $EXTRA1 $EXTRA2 --vcf tmp_h.vcf --double-id --make-bed --out tmp_data
# PLINK 1.9 ignores the phenotypes of missing-sex samples, which would drop its
# CASE/CTRL sections, so the fixture assigns sex as well.
awk 'BEGIN{OFS=" "} {print $1, $2, 0, 0, (NR % 2) + 1, ((NR % 3)? 2 : 1)}' tmp_data.fam > tmp_ped.fam
mv tmp_ped.fam tmp_data.fam

compare_pair() {
    plink --bfile tmp_data --twolocus "$1" "$2" --out plink19
    $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --twolocus "$1" "$2" --out plink2
    awk -v id1="$1" -v id2="$2" -f compare.awk plink19.twolocus plink2.twolocus
}

# 1. A pair inside one haplotype block, a pair spanning two, and a pair
#    involving a variant with missing calls.
compare_pair b0v0 b0v3
compare_pair b0v1 b1v7
compare_pair b5v2 b5v9

# 2. Without a case/control phenotype there is only the ALL group.
awk 'BEGIN{OFS=" "} {print $1, $2, 0, 0, (NR % 2) + 1, -9}' tmp_data.fam > tmp_nopheno.fam
cp tmp_data.bed tmp_nopheno.bed
cp tmp_data.bim tmp_nopheno.bim
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_nopheno --twolocus b0v0 b0v3 --out plink2_nopheno
test "$(awk '!/^#/ {print $1}' plink2_nopheno.twolocus | sort -u | tr -d '\n')" = "ALL"

# 3. Every cell count is present, so the marginals PLINK 1.9 prints can be
#    recovered: the ALL group has to sum to the sample count.
test "$(awk '!/^#/ && $1 == "ALL" {s += $6} END {print s}' plink2.twolocus)" -eq "$(wc -l < tmp_data.fam)"

# 4. A duplicate or unknown variant ID is an error rather than a wrong report.
if $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --twolocus b0v0 nosuchvariant --out plink2_bad 2> tmp_bad_err.txt; then
    echo "--twolocus unexpectedly accepted an unknown variant ID"
    exit 1
fi

# 5. 'zs' is the same report, compressed.
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --twolocus zs b5v2 b5v9 --out plink2_zs
$BUILD/plink2 $EXTRA1 $EXTRA2 --zst-decompress plink2_zs.twolocus.zst > plink2_zs.twolocus
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --twolocus b5v2 b5v9 --out plink2_plain
diff -q plink2_plain.twolocus plink2_zs.twolocus
