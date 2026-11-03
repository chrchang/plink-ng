#!/bin/bash

# --test-mishap.
#
# The fixture has 400 samples, which is what exposed an unaligned genotype
# buffer: the three ring slots were spaced by the packed word count, so
# PgrGetInv1()'s in-place inversion tripped its alignment assertion.  Keeping
# a sample count whose packed and aligned word counts differ is the regression
# test for that.
#
# The HETERO rows are compared against PLINK 1.9 numerically, and every
# haplotype row by its label (compare_labels.awk): plink2 used to name each
# flanking allele after the major allele where PLINK 1.x names the nonmajor
# one, so all haplotype rows but the symmetric ones carried the wrong label.
# Numerically degenerate EM rows (counts like 9.99e-16/2.02e-14 on either side)
# are skipped.

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
awk -f compare_labels.awk plink19.missing.hap plink2.missing.hap

# 2. Same, with --maf, which changes which variants have flanking neighbors.
plink --bfile tmp_data --maf 0.1 --test-mishap --out plink19_maf
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --maf 0.1 --test-mishap --out plink2_maf
awk -f compare_hetero.awk plink19_maf.missing.hap plink2_maf.missing.hap
awk -f compare_labels.awk plink19_maf.missing.hap plink2_maf.missing.hap

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

# 6. Multiallelic variants are handled as major vs. rest: a fileset whose
#    nonmajor alleles have been merged into one gives the same numbers.  Only
#    HAPLOTYPE differs ('.' for the pooled side, and the row order when pooling
#    makes the rest outnumber the old major allele), so the rows are compared
#    as sorted sets without that column.
awk -f make_multi_vcf.awk
$BUILD/plink2 $EXTRA1 $EXTRA2 --vcf tmp_multi.vcf --test-mishap --out plink2_multi
$BUILD/plink2 $EXTRA1 $EXTRA2 --vcf tmp_pooled.vcf --test-mishap --out plink2_pooled
test "$(grep -c '^v' plink2_multi.missing.hap)" -gt 100
grep -q '\.' <(cut -f 2 plink2_multi.missing.hap)
diff <(cut -f 1,3- plink2_multi.missing.hap | sort) <(cut -f 1,3- plink2_pooled.missing.hap | sort)

# 6. Known phase.  A flanking double-het phased at both flanks is a known pair
#    of haplotypes, as in --r2-phased, so a dataset with every sample written
#    twice must give the same haplotype rows as one with every haplotype
#    written as a homozygous sample: both have the same haplotype counts and
#    no double-hets of unknown phase.  (The HETERO rows count heterozygous
#    samples, so they differ by construction and are left out.)  Erasing the
#    phase must change the haplotype rows on this fixture, or the check above
#    proves nothing.
for m in phased dup haphom
do
    awk -v mode=$m -f make_vcf_phased.awk > tmp_$m.vcf
    $BUILD/plink2 $EXTRA1 $EXTRA2 --vcf tmp_$m.vcf --make-pgen --out tmp_$m > /dev/null
    $BUILD/plink2 $EXTRA1 $EXTRA2 --pfile tmp_$m --test-mishap --out plink2_$m
    grep -v HETERO plink2_$m.missing.hap > tmp_$m.haprows
done
test "$(wc -l < tmp_dup.haprows)" -gt 100
diff -q tmp_dup.haprows tmp_haphom.haprows
$BUILD/plink2 $EXTRA1 $EXTRA2 --pfile tmp_phased --make-pgen erase-phase --out tmp_unphased > /dev/null
$BUILD/plink2 $EXTRA1 $EXTRA2 --pfile tmp_unphased --test-mishap --out plink2_unphased
if diff -q plink2_phased.missing.hap plink2_unphased.missing.hap > /dev/null; then
    echo "erasing phase did not change the --test-mishap report"
    exit 1
fi
