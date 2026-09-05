#!/bin/bash

# "--mgf", the BIMBAM mean genotype reader, checked as the inverse of
# "--export mgf".

set -exo pipefail

BUILD=$1
EXTRA1=$2
EXTRA2=$3

# 1. Hardcalls.  A round trip through the text format has to reproduce the
#    .mgf and .pos.txt byte for byte, and the .pvar apart from its header.
fails_with() {
    local expected=$1
    shift
    if "$@" > /dev/null 2> tmp_err.txt; then
        echo "expected failure: $*"
        exit 1
    fi
    # Error messages are word-wrapped, so the newlines have to come out before
    # a phrase can be matched across them.
    tr '\n' ' ' < tmp_err.txt | tr -s ' ' | grep -q "$expected"
}

plink --simulate simulate.txt --simulate-missing 0.03 --simulate-ncases 60 --simulate-ncontrols 60 --out tmp_hard > /dev/null
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_hard --export mgf --out tmp_hard_e
$BUILD/plink2 $EXTRA1 $EXTRA2 --mgf tmp_hard_e.mgf tmp_hard_e.pos.txt --out tmp_hard_r
$BUILD/plink2 $EXTRA1 $EXTRA2 --pfile tmp_hard_r --export mgf --out tmp_hard_e2
diff -q tmp_hard_e.mgf tmp_hard_e2.mgf
diff -q tmp_hard_e.pos.txt tmp_hard_e2.pos.txt

# The .pvar has to name the same variants at the same positions, with REF and
# ALT on the same sides: the .mgf's dosages count its second column, which is
# ALT, so a swap here would silently invert every call.
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_hard --make-just-pvar --out tmp_hard_p
diff -q <(grep -v '^#' tmp_hard_p.pvar) <(grep -v '^#' tmp_hard_r.pvar)

# Sample IDs are not in the format, so they are synthesized as with --dummy.
diff -q <(printf '#IID\n'; seq 0 119 | sed 's/^/per/') tmp_hard_r.psam

# 2. Dosages.  --export mgf writes five significant digits, so a second trip
#    through it is exact even though the first one rounds.
$BUILD/plink2 $EXTRA1 $EXTRA2 --dummy 120 400 0.03 dosage-freq=0.7 --seed 1 --out tmp_dose
$BUILD/plink2 $EXTRA1 $EXTRA2 --pfile tmp_dose --export mgf --out tmp_dose_e
$BUILD/plink2 $EXTRA1 $EXTRA2 --mgf tmp_dose_e.mgf tmp_dose_e.pos.txt --out tmp_dose_r
$BUILD/plink2 $EXTRA1 $EXTRA2 --pfile tmp_dose_r --export mgf --out tmp_dose_e2
diff -q tmp_dose_e.mgf tmp_dose_e2.mgf
# The fixture has to actually exercise the dosage path.
grep -q '\.' tmp_dose_e.mgf

# 3. Commas are separators in both files, and 'NA', '-9', '?', '??' and '.'
#    all mark a missing call.
tr '\t' ',' < tmp_dose_e.mgf > tmp_comma.mgf
tr '\t' ',' < tmp_dose_e.pos.txt > tmp_comma.pos.txt
$BUILD/plink2 $EXTRA1 $EXTRA2 --mgf tmp_comma.mgf tmp_comma.pos.txt --out tmp_comma_r
$BUILD/plink2 $EXTRA1 $EXTRA2 --pfile tmp_comma_r --export mgf --out tmp_comma_e
diff -q tmp_dose_e.mgf tmp_comma_e.mgf

for code in -9 '?' '??' '.'; do
    sed "s/\bNA\b/${code}/g" tmp_dose_e.mgf > tmp_miss.mgf
    $BUILD/plink2 $EXTRA1 $EXTRA2 --mgf tmp_miss.mgf tmp_dose_e.pos.txt --out tmp_miss_r
    $BUILD/plink2 $EXTRA1 $EXTRA2 --pfile tmp_miss_r --export mgf --out tmp_miss_e
    diff -q tmp_dose_e.mgf tmp_miss_e.mgf
done

# 4. 'pheno=' loads one row per sample; the columns become PHENO1, PHENO2, ...
awk 'BEGIN {print "#IID\tqt1\tqt2"; for (i = 0; i < 120; ++i) {print "per" i "\t" (i / 8) "\t" (120 - i)}}' > tmp_pheno.txt
$BUILD/plink2 $EXTRA1 $EXTRA2 --pfile tmp_dose --pheno tmp_pheno.txt --export mgf --out tmp_ph_e
test -s tmp_ph_e.pheno.txt
$BUILD/plink2 $EXTRA1 $EXTRA2 --mgf tmp_ph_e.mgf tmp_ph_e.pos.txt pheno=tmp_ph_e.pheno.txt --out tmp_ph_r
head -n 1 tmp_ph_r.psam | grep -q '^#IID	PHENO1	PHENO2$'
diff -q <(cut -f 2- tmp_ph_r.psam | tail -n +2) tmp_ph_e.pheno.txt

# Without 'pheno=' there are no phenotype columns at all.
head -n 1 tmp_dose_r.psam | grep -q '^#IID$'

# 5. --psam/--fam supplies real sample IDs in place of the synthesized ones,
#    without touching the genotypes.
awk 'BEGIN {print "#IID\tSEX\tPHENO1"; for (i = 0; i < 120; ++i) {print "s" i "\t" (i % 2 + 1) "\t" (i % 2 + 1)}}' > tmp_real.psam
$BUILD/plink2 $EXTRA1 $EXTRA2 --mgf tmp_dose_e.mgf tmp_dose_e.pos.txt --psam tmp_real.psam --out tmp_ids
diff -q <(cut -f 1 tmp_ids.psam | tail -n +2) <(seq 0 119 | sed 's/^/s/')
$BUILD/plink2 $EXTRA1 $EXTRA2 --pfile tmp_ids --export mgf --out tmp_ids_e
diff -q tmp_dose_e.mgf tmp_ids_e.mgf

# The PLINK 1.x .fam spelling works too, and brings FID along.
awk 'BEGIN {for (i = 0; i < 120; ++i) {print "f" i, "s" i, 0, 0, (i % 2 + 1), (i % 2 + 1)}}' > tmp_real.fam
$BUILD/plink2 $EXTRA1 $EXTRA2 --mgf tmp_dose_e.mgf tmp_dose_e.pos.txt --fam tmp_real.fam --out tmp_ids2
head -n 1 tmp_ids2.psam | grep -q '^#FID	IID'
diff -q <(cut -f 2 tmp_ids2.psam | tail -n +2) <(seq 0 119 | sed 's/^/s/')

# A row count that disagrees with the .mgf is an error, and so is combining
# the two phenotype sources.
head -n 6 tmp_real.psam > tmp_short.psam
fails_with "tmp_short.psam has 5 samples, while tmp_dose_e.mgf implies 120" \
    $BUILD/plink2 $EXTRA1 $EXTRA2 --mgf tmp_dose_e.mgf tmp_dose_e.pos.txt --psam tmp_short.psam --out tmp_bad7
fails_with "cannot be used with --psam/--fam" \
    $BUILD/plink2 $EXTRA1 $EXTRA2 --mgf tmp_ph_e.mgf tmp_ph_e.pos.txt pheno=tmp_ph_e.pheno.txt --psam tmp_real.psam --out tmp_bad8

# 6. A block-gzipped .mgf reads back the same as the plain one.
$BUILD/plink2 $EXTRA1 $EXTRA2 --pfile tmp_dose --export mgf bgz --out tmp_bgz
$BUILD/plink2 $EXTRA1 $EXTRA2 --mgf tmp_bgz.mgf.gz tmp_bgz.pos.txt --out tmp_bgz_r
$BUILD/plink2 $EXTRA1 $EXTRA2 --pfile tmp_bgz_r --export mgf --out tmp_bgz_e
diff -q tmp_dose_e.mgf tmp_bgz_e.mgf

# 7. Malformed input has to be rejected, naming the offending line.

# A variant with no position.
head -n 3 tmp_hard_e.pos.txt > tmp_short.pos.txt
fails_with "line 4 of tmp_hard_e.mgf is absent from tmp_short.pos.txt" \
    $BUILD/plink2 $EXTRA1 $EXTRA2 --mgf tmp_hard_e.mgf tmp_short.pos.txt --out tmp_bad1

# A dosage that is neither a number in [0, 2] nor a missing code.  Reading it
# as missing would quietly drop data.
awk 'NR == 2 {$5 = "xyz"} {print}' OFS='\t' tmp_hard_e.mgf > tmp_badtok.mgf
fails_with "Invalid dosage 'xyz' on line 2" \
    $BUILD/plink2 $EXTRA1 $EXTRA2 --mgf tmp_badtok.mgf tmp_hard_e.pos.txt --out tmp_bad2
awk 'NR == 3 {$6 = "2.5"} {print}' OFS='\t' tmp_hard_e.mgf > tmp_badrange.mgf
fails_with "Invalid dosage '2.5' on line 3" \
    $BUILD/plink2 $EXTRA1 $EXTRA2 --mgf tmp_badrange.mgf tmp_hard_e.pos.txt --out tmp_bad3

# A row with the wrong number of samples.
awk 'NR == 5 {NF = NF - 1} {print}' OFS='\t' tmp_hard_e.mgf > tmp_ragged.mgf
fails_with "Line 5 of tmp_ragged.mgf" \
    $BUILD/plink2 $EXTRA1 $EXTRA2 --mgf tmp_ragged.mgf tmp_hard_e.pos.txt --out tmp_bad4

# A phenotype file whose length disagrees with the .mgf.
head -n 5 tmp_ph_e.pheno.txt > tmp_short.pheno.txt
fails_with "5 rows, while tmp_ph_e.mgf implies 120 samples" \
    $BUILD/plink2 $EXTRA1 $EXTRA2 --mgf tmp_ph_e.mgf tmp_ph_e.pos.txt pheno=tmp_short.pheno.txt --out tmp_bad5

# A duplicate ID in the position file, which would make the lookup ambiguous.
{ cat tmp_hard_e.pos.txt; head -n 1 tmp_hard_e.pos.txt; } > tmp_dup.pos.txt
fails_with "Duplicate variant ID" \
    $BUILD/plink2 $EXTRA1 $EXTRA2 --mgf tmp_hard_e.mgf tmp_dup.pos.txt --out tmp_bad6
