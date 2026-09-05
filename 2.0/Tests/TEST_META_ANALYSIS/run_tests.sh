#!/bin/bash

# --meta-analysis, checked against PLINK 1.9.
#
# Two simulated case/control studies of the same 400 variants, with different
# sample sizes, run through PLINK 1.9's --logistic (with --ci, which is what
# adds the SE column the meta-analysis needs).

set -exo pipefail

BUILD=$1
EXTRA1=$2
EXTRA2=$3

# One simulated fileset split into two studies by sample, so both reports
# cover the same variants with the same alleles; two independent --simulate
# runs would disagree on A1/A2 for some variants and the comparison would then
# be about allele-mismatch handling rather than the meta-analysis itself.
plink --simulate simulate.txt --simulate-ncases 800 --simulate-ncontrols 800 --simulate-missing 0.02 --out tmp_all > /dev/null
head -n 800 tmp_all.fam | cut -f1,2 -d' ' > tmp_studyA.txt
tail -n 800 tmp_all.fam | cut -f1,2 -d' ' > tmp_studyB.txt
plink --bfile tmp_all --keep tmp_studyA.txt --make-bed --out tmp_s1 > /dev/null
plink --bfile tmp_all --keep tmp_studyB.txt --make-bed --out tmp_s2 > /dev/null
plink --bfile tmp_s1 --logistic --ci 0.95 --out tmp_r1 > /dev/null
plink --bfile tmp_s2 --logistic --ci 0.95 --out tmp_r2 > /dev/null
plink --bfile tmp_s1 --logistic beta --ci 0.95 --out tmp_b1 > /dev/null
plink --bfile tmp_s2 --logistic beta --ci 0.95 --out tmp_b2 > /dev/null

compare() {
    awk -f compare.awk "$1" "$2"
}

# 1. Defaults.
plink --meta-analysis tmp_r1.assoc.logistic tmp_r2.assoc.logistic --out plink19
$BUILD/plink2 $EXTRA1 $EXTRA2 --meta-analysis tmp_r1.assoc.logistic tmp_r2.assoc.logistic --out plink2
compare plink19.meta plink2.meta

# 2. Each mode that changes the row set or the effect-size column.
for mode in "report-all" "no-map" "no-allele" "study"
do
    plink --meta-analysis tmp_r1.assoc.logistic tmp_r2.assoc.logistic + $mode --out plink19_$mode
    $BUILD/plink2 $EXTRA1 $EXTRA2 --meta-analysis tmp_r1.assoc.logistic tmp_r2.assoc.logistic + $mode --out plink2_$mode
    compare plink19_$mode.meta plink2_$mode.meta
done

# 3. 'logscale', against reports that carry BETA instead of OR.
plink --meta-analysis tmp_b1.assoc.logistic tmp_b2.assoc.logistic + logscale --out plink19_log
$BUILD/plink2 $EXTRA1 $EXTRA2 --meta-analysis tmp_b1.assoc.logistic tmp_b2.assoc.logistic + logscale --out plink2_log
compare plink19_log.meta plink2_log.meta

# 4. 'weighted-z' adds two columns; PLINK 1.9 computes the same Z.
plink --meta-analysis tmp_r1.assoc.logistic tmp_r2.assoc.logistic + weighted-z --out plink19_wz
$BUILD/plink2 $EXTRA1 $EXTRA2 --meta-analysis tmp_r1.assoc.logistic tmp_r2.assoc.logistic + weighted-z --out plink2_wz
compare plink19_wz.meta plink2_wz.meta

# 5. The field-name flags, on reports whose columns have been renamed so the
#    defaults cannot match.
awk 'NR == 1 {sub(/ CHR /, " CHROMOSOME "); sub(/ SNP /, " MARKER "); sub(/ BP /, " POSITION "); sub(/ SE /, " STDERR ")} {print}' tmp_r1.assoc.logistic > tmp_renamed1.txt
awk 'NR == 1 {sub(/ CHR /, " CHROMOSOME "); sub(/ SNP /, " MARKER "); sub(/ BP /, " POSITION "); sub(/ SE /, " STDERR ")} {print}' tmp_r2.assoc.logistic > tmp_renamed2.txt
$BUILD/plink2 $EXTRA1 $EXTRA2 --meta-analysis tmp_renamed1.txt tmp_renamed2.txt \
    --meta-analysis-chr-field CHROMOSOME --meta-analysis-snp-field MARKER \
    --meta-analysis-bp-field POSITION --meta-analysis-se-field STDERR --out plink2_renamed
compare plink19.meta plink2_renamed.meta

# 6. Without those flags the renamed reports must be rejected, not silently
#    mis-parsed.
if $BUILD/plink2 $EXTRA1 $EXTRA2 --meta-analysis tmp_renamed1.txt tmp_renamed2.txt --out plink2_bad 2> tmp_bad_err.txt; then
    echo "--meta-analysis unexpectedly accepted a report with no recognized CHR field"
    exit 1
fi

# 7. 'zs' is the same report, compressed.
$BUILD/plink2 $EXTRA1 $EXTRA2 --meta-analysis tmp_r1.assoc.logistic tmp_r2.assoc.logistic + zs --out plink2_zs
$BUILD/plink2 $EXTRA1 $EXTRA2 --zst-decompress plink2_zs.meta.zst > plink2_zs.meta
diff -q plink2.meta plink2_zs.meta
