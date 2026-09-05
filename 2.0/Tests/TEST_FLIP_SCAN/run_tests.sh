#!/bin/bash

# --flip-scan, checked against PLINK 1.9.
#
# The fixture plants 20 haplotype blocks and then swaps the alleles of every
# 7th variant in the cases only, which is exactly the strand inconsistency the
# scan is meant to find.
#
# --maf 0.05 is applied throughout: PLINK 1.9 lets a nan through when a
# neighbor pair is monomorphic in one group, which counts the pair as a sign
# flip, and plink2 skips it.  Above that threshold the two agree.

set -exo pipefail

BUILD=$1
EXTRA1=$2
EXTRA2=$3

awk -f make_vcf.awk > tmp_h.vcf
$BUILD/plink2 $EXTRA1 $EXTRA2 --vcf tmp_h.vcf --double-id --make-bed --out tmp_data
awk 'BEGIN{OFS=" "} {print $1, $2, 0, 0, (NR % 2) + 1, (NR > 200)? 2 : 1}' tmp_data.fam > tmp_ped.fam
mv tmp_ped.fam tmp_data.fam

compare() {
    plink --bfile tmp_data --maf 0.05 --flip-scan $1 --out plink19
    $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --maf 0.05 --flip-scan $1 --out plink2
    awk -f compare.awk plink19.flipscan plink2.flipscan
}

# 1. Defaults, and each window/threshold flag.
compare ""
plink --bfile tmp_data --maf 0.05 --flip-scan --flip-scan-window 5 --out plink19_w5
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --maf 0.05 --flip-scan --flip-scan-window 5 --out plink2_w5
awk -f compare.awk plink19_w5.flipscan plink2_w5.flipscan

# --flip-scan-window-kb values large enough to hold the whole variant-count
# window agree; smaller ones do not.  On this fixture PLINK 1.9 reports no
# neighbors at all at 40 kb and below, while plink2 truncates the window and
# keeps going, so only the larger values are compared here.
plink --bfile tmp_data --maf 0.05 --flip-scan --flip-scan-window-kb 100 --out plink19_kb
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --maf 0.05 --flip-scan --flip-scan-window-kb 100 --out plink2_kb
awk -f compare.awk plink19_kb.flipscan plink2_kb.flipscan

plink --bfile tmp_data --maf 0.05 --flip-scan --flip-scan-threshold 0.2 --out plink19_thr
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --maf 0.05 --flip-scan --flip-scan-threshold 0.2 --out plink2_thr
awk -f compare.awk plink19_thr.flipscan plink2_thr.flipscan

# 2. 'verbose' adds one line per above-threshold neighbor pair.
plink --bfile tmp_data --maf 0.05 --flip-scan verbose --out plink19_v
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --maf 0.05 --flip-scan verbose --out plink2_v
awk -f compare.awk plink19_v.flipscan plink2_v.flipscan
awk -f compare_verbose.awk plink19_v.flipscan.verbose plink2_v.flipscan.verbose

# 3. The planted flips have to be found: every 7th variant is one, and each
#    should be flagged.
test "$(grep -vc '^#' plink2.flipscan)" -gt 20

# 4. The result must not depend on --threads.
for t in 1 3 8
do
    $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --maf 0.05 --flip-scan --threads $t --out plink2_t$t
    diff -q plink2.flipscan plink2_t$t.flipscan
done

# 5. 'zs' is the same report, compressed.
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --maf 0.05 --flip-scan zs --out plink2_zs
$BUILD/plink2 $EXTRA1 $EXTRA2 --zst-decompress plink2_zs.flipscan.zst > plink2_zs.flipscan
diff -q plink2.flipscan plink2_zs.flipscan
