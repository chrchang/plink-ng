#!/bin/bash

# --clump-index-first and --clump-replicate, checked against PLINK 1.9.
#
# Two association reports over the same fileset -- one on all samples, one on
# half of them -- give the two-file setup both flags are about.

set -exo pipefail

BUILD=$1
EXTRA1=$2
EXTRA2=$3

awk -f make_vcf.awk > tmp_h.vcf
$BUILD/plink2 $EXTRA1 $EXTRA2 --vcf tmp_h.vcf --double-id --make-bed --out tmp_data
awk 'BEGIN{OFS=" "} {print $1, $2, 0, 0, (NR % 2) + 1, (NR % 2) + 1}' tmp_data.fam > tmp_ped.fam
mv tmp_ped.fam tmp_data.fam

plink --bfile tmp_data --assoc --out tmp_a1 > /dev/null
head -n 200 tmp_data.fam | cut -d' ' -f1,2 > tmp_half.txt
plink --bfile tmp_data --keep tmp_half.txt --assoc --out tmp_a2 > /dev/null

# The phenotype is unrelated to the genotypes, so the thresholds have to be
# loose enough for clumps to form at all.
P1="--clump-p1 0.5 --clump-p2 0.8"

compare() {
    plink --bfile tmp_data --clump tmp_a1.assoc,tmp_a2.assoc $P1 $1 --out plink19
    $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --clump tmp_a1.assoc,tmp_a2.assoc $P1 $1 --out plink2
    awk -f compare.awk plink19.clumped plink2.clumps
}

# 1. Plain --clump over two files, as the baseline.
compare ""

# 2. --clump-index-first: every index variant comes from the first file, and
#    its reported p-value and F are that file's.
compare "--clump-index-first"
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --clump tmp_a1.assoc,tmp_a2.assoc $P1 --clump-index-first --out plink2_if
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --clump tmp_a1.assoc,tmp_a2.assoc $P1 --out plink2_plain
test "$(awk '!/^#/ && NF && $4 != 1' plink2_if.clumps | wc -l)" -eq 0
test "$(awk '!/^#/ && NF && $4 != 1' plink2_plain.clumps | wc -l)" -gt 0

# 3. --clump-replicate: clumps whose below-p2 results all come from one file
#    are dropped, so it can only shrink the report.
compare "--clump-replicate"
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --clump tmp_a1.assoc,tmp_a2.assoc $P1 --clump-replicate --out plink2_rep
test "$(grep -vc '^#' plink2_rep.clumps)" -le "$(grep -vc '^#' plink2_plain.clumps)"

# 4. Both together.
compare "--clump-index-first --clump-replicate"

# 5. A single file with --clump-replicate leaves nothing, since no clump can
#    have results from two files.
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --clump tmp_a1.assoc $P1 --clump-replicate --out plink2_one 2> tmp_one_err.txt || true
test "$(awk '!/^#/ && NF' plink2_one.clumps 2> /dev/null | wc -l)" -eq 0
