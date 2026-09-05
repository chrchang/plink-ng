#!/bin/bash

# --test-mishap.
#
# The fixture has 400 samples, which is what exposed an unaligned genotype
# buffer: the three ring slots were spaced by the packed word count, so
# PgrGetInv1()'s in-place inversion tripped its alignment assertion.  Keeping
# a sample count whose packed and aligned word counts differ is the regression
# test for that.
#
# The HETERO rows are compared against PLINK 1.9.  Its haplotype rows are not
# comparable in general: the two programs label haplotypes by different
# alleles, and 1.9's EM step emits numerically degenerate rows (counts like
# 9.99e-16/2.02e-14) where plink2 reports 0/0.

set -exo pipefail

BUILD=$1
EXTRA1=$2
EXTRA2=$3

awk -f make_vcf.awk > tmp_h.vcf
$BUILD/plink2 $EXTRA1 $EXTRA2 --vcf tmp_h.vcf --double-id --make-bed --out tmp_data
test "$(wc -l < tmp_data.fam)" -eq 400

# 1. HETERO rows against PLINK 1.9.
plink --bfile tmp_data --test-mishap --out plink19
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --test-mishap --out plink2
awk -f compare_hetero.awk plink19.missing.hap plink2.missing.hap

# 2. Same, with --maf, which changes which variants have flanking neighbors.
plink --bfile tmp_data --maf 0.1 --test-mishap --out plink19_maf
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --maf 0.1 --test-mishap --out plink2_maf
awk -f compare_hetero.awk plink19_maf.missing.hap plink2_maf.missing.hap

# 3. Structure: every checked locus gets one HETERO row plus one row per
#    observed flanking haplotype (four when all four occur, two when only two
#    do), and the two counts in each M_H field are consistent.
awk '
    !/^#/ {
        ++rows[$1];
        if ($2 == "HETERO") { ++hetero[$1] }
        split($5, a, "/"); split($6, b, "/");
        if (a[1] > a[2] || b[1] > b[2]) { print "count exceeds its total on " $1; failed = 1; exit 1 }
    }
    END {
        if (failed) { exit 1 }
        for (id in rows) {
            if (rows[id] != 5 && rows[id] != 3) { print id " has " rows[id] " rows, expected 3 or 5"; exit 1 }
            if (hetero[id] != 1) { print id " has " hetero[id] " HETERO rows"; exit 1 }
            ++n;
        }
        if (!n) { print "no loci checked"; exit 1 }
        print n " loci have the expected row structure"
    }
' plink2.missing.hap

# 4. The result must not depend on --threads.
for t in 1 3 8
do
    $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --test-mishap --threads $t --out plink2_t$t
    diff -q plink2.missing.hap plink2_t$t.missing.hap
done

# 5. 'zs' is the same report, compressed.
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --test-mishap zs --out plink2_zs
$BUILD/plink2 $EXTRA1 $EXTRA2 --zst-decompress plink2_zs.missing.hap.zst > plink2_zs.missing.hap
diff -q plink2.missing.hap plink2_zs.missing.hap
