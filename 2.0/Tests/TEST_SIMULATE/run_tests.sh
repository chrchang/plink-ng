#!/bin/bash

# --simulate and --simulate-qt, checked against PLINK 1.9.
#
# Both are supposed to produce exactly the dataset PLINK 1.9 does given the
# same --seed, so the comparison is byte-for-byte rather than statistical.
# plink2 writes a .pgen fileset, so each run is converted to .bed first.

set -exo pipefail

BUILD=$1
EXTRA1=$2
EXTRA2=$3

compare() {
    # $1: PLINK 1.9 output prefix, $2: plink2 output prefix
    $BUILD/plink2 $EXTRA1 $EXTRA2 --pfile "$2" --make-bed --out "$2"_bed
    cmp "$1".bed "$2"_bed.bed
    diff -q "$1".bim "$2"_bed.bim
    diff -q <(tr -s ' \t' ' ' < "$1".fam) <(tr -s ' \t' ' ' < "$2"_bed.fam)
    diff -q "$1".simfreq "$2".simfreq
}

# 1. --simulate, with case/control counts and missingness.
plink --simulate simulate.txt --simulate-ncases 200 --simulate-ncontrols 300 --simulate-missing 0.03 --seed 12345 --out plink19_cc
$BUILD/plink2 $EXTRA1 $EXTRA2 --simulate simulate.txt --simulate-ncases 200 --simulate-ncontrols 300 --simulate-missing 0.03 --seed 12345 --out plink2_cc
compare plink19_cc plink2_cc

# 2. --simulate-qt.
plink --simulate-qt qt.txt --simulate-n 400 --simulate-missing 0.02 --seed 777 --out plink19_qt
$BUILD/plink2 $EXTRA1 $EXTRA2 --simulate-qt qt.txt --simulate-n 400 --simulate-missing 0.02 --seed 777 --out plink2_qt
compare plink19_qt plink2_qt

# 3. Each allele-coding modifier.
for coding in acgt 1234 12
do
    plink --simulate simulate.txt $coding --simulate-ncases 100 --simulate-ncontrols 100 --seed 5 --out plink19_$coding
    $BUILD/plink2 $EXTRA1 $EXTRA2 --simulate simulate.txt $coding --simulate-ncases 100 --simulate-ncontrols 100 --seed 5 --out plink2_$coding
    compare plink19_$coding plink2_$coding
done

# 4. 'tags' and 'haps', which take the nine-column parameter file.
for mode in tags haps
do
    plink --simulate simulate_tags.txt $mode --simulate-ncases 100 --simulate-ncontrols 100 --seed 9 --out plink19_$mode
    $BUILD/plink2 $EXTRA1 $EXTRA2 --simulate simulate_tags.txt $mode --simulate-ncases 100 --simulate-ncontrols 100 --seed 9 --out plink2_$mode
    compare plink19_$mode plink2_$mode
done

plink --simulate-qt qt_tags.txt tags --simulate-n 300 --seed 21 --out plink19_qttags
$BUILD/plink2 $EXTRA1 $EXTRA2 --simulate-qt qt_tags.txt tags --simulate-n 300 --seed 21 --out plink2_qttags
compare plink19_qttags plink2_qttags

# 5. --simulate-prevalence and --simulate-label.
plink --simulate simulate.txt --simulate-ncases 150 --simulate-ncontrols 150 --simulate-prevalence 0.05 --simulate-label FAM --seed 31 --out plink19_prev
$BUILD/plink2 $EXTRA1 $EXTRA2 --simulate simulate.txt --simulate-ncases 150 --simulate-ncontrols 150 --simulate-prevalence 0.05 --simulate-label FAM --seed 31 --out plink2_prev
compare plink19_prev plink2_prev

# 6. The seed is what makes this reproducible: same seed twice agrees, a
#    different seed does not.
$BUILD/plink2 $EXTRA1 $EXTRA2 --simulate simulate.txt --simulate-ncases 200 --simulate-ncontrols 300 --simulate-missing 0.03 --seed 12345 --out plink2_again
cmp plink2_cc.pgen plink2_again.pgen
$BUILD/plink2 $EXTRA1 $EXTRA2 --simulate simulate.txt --simulate-ncases 200 --simulate-ncontrols 300 --simulate-missing 0.03 --seed 54321 --out plink2_other
if cmp -s plink2_cc.pgen plink2_other.pgen; then
    echo "a different --seed produced the same dataset"
    exit 1
fi
