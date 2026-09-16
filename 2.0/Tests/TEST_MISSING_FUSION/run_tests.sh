#!/bin/bash

# --missing's per-sample counts used to come from a dedicated pass over the
# genotypes.  They are now tallied inside the allele-frequency pass for the
# autosomal biallelic part of a dosage-free file, with chrX/chrY/haploid left
# to the old pass.  Every case below pins the two against each other, since a
# mistake there would silently change QC numbers rather than fail loudly.
#
# --mind forces the old path: its filter has to run before allele frequencies
# are computed, so the counts cannot ride along with them.  The threshold is
# set so that nothing is actually removed, which makes the two runs comparable.

set -exo pipefail

UNFUSED="--mind 0.99999"
SCOLS="scols=maybefid,maybesid,nmissdosage,nmiss,nmisshh,hethap,nobs,fmissdosage,fmiss,fmisshh"

awk -f make_vcf.awk > tmp_h.vcf
awk -v multi=1 -f make_vcf.awk > tmp_hm.vcf
for i in $(seq 0 119); do
    s=$(printf 's%03d' $i)
    if [ $(( i % 17 )) -eq 0 ]; then
        echo -e "$s\t$s\t0"
    elif [ $(( i % 2 )) -eq 0 ]; then
        echo -e "$s\t$s\t1"
    else
        echo -e "$s\t$s\t2"
    fi
done > tmp_sex.txt
awk 'NR % 2 == 0 {print $1, $2}' tmp_sex.txt > tmp_keep.txt

BUILD=$1
EXTRA1=$2
EXTRA2=$3

$1/plink2 $2 $3 --vcf tmp_h.vcf --double-id --update-sex tmp_sex.txt --split-par b37 --make-pgen --out tmp_all
# autosomes only: the separate pass disappears entirely
$1/plink2 $2 $3 --pfile tmp_all --chr 1,2 --make-pgen --out tmp_auto
# a multiallelic file, which is not eligible and must fall back
$BUILD/plink2 $EXTRA1 $EXTRA2 --vcf tmp_hm.vcf --double-id --update-sex tmp_sex.txt --split-par b37 --make-pgen --out tmp_multi

check() {
    local ds="$1"
    shift
    $BUILD/plink2 $EXTRA1 $EXTRA2 --pfile $ds "$@" --out tmp_fused
    $BUILD/plink2 $EXTRA1 $EXTRA2 --pfile $ds "$@" $UNFUSED --out tmp_unfused
    # the threshold must not actually drop anyone, or the reports would cover
    # different samples
    grep -q "^0 samples removed due to missing genotype data" tmp_unfused.log
    diff -q tmp_fused.smiss tmp_unfused.smiss
    diff -q tmp_fused.vmiss tmp_unfused.vmiss
}

for ds in tmp_all tmp_auto tmp_multi; do
    check $ds --missing
    check $ds --missing $SCOLS
    check $ds --missing --freq --hardy --geno-counts
    check $ds --missing --nonfounders
    check $ds --missing $SCOLS --y-nosex-missing-stats
    check $ds --missing --geno 0.15
    # a sample subset makes the fused path derive the genotype counts itself
    check $ds --missing --keep tmp_keep.txt
    check $ds --missing --keep tmp_keep.txt --freq --hardy
done

# The autosome-only file needs no separate pass at all; the mixed one still
# runs a small one for chrX/chrY/chrMT.
$BUILD/plink2 $EXTRA1 $EXTRA2 --pfile tmp_auto --missing --out tmp_a > tmp_a.stdout
test "$(grep -c 'Calculating sample missingness' tmp_a.log)" -eq 0
$BUILD/plink2 $EXTRA1 $EXTRA2 --pfile tmp_all --missing --out tmp_m > tmp_m.stdout
test "$(grep -c 'Calculating sample missingness' tmp_m.log)" -eq 1
# ...and a multiallelic file is not eligible, so it keeps the full pass.
$BUILD/plink2 $EXTRA1 $EXTRA2 --pfile tmp_multi --missing --out tmp_x > tmp_x.stdout
test "$(grep -c 'Calculating sample missingness' tmp_x.log)" -eq 1

# --mind still filters on the same counts it always did.
$BUILD/plink2 $EXTRA1 $EXTRA2 --pfile tmp_all --mind 0.08 --missing --out tmp_mind
$BUILD/plink2 $EXTRA1 $EXTRA2 --pfile tmp_all --mind 0.08 --write-samples --out tmp_mind2
diff -q <(cut -f 1,2 tmp_mind.smiss | tail -n +2) <(tail -n +2 tmp_mind2.id)

# Cross-check the autosomal counts against PLINK 1.9, which computes them its
# own way.
$BUILD/plink2 $EXTRA1 $EXTRA2 --pfile tmp_auto --export ped --out tmp_auto19
plink --file tmp_auto19 --missing --out tmp_19
$BUILD/plink2 $EXTRA1 $EXTRA2 --pfile tmp_auto --missing --out tmp_p2
awk 'NR > 1 {print $2, $3}' tmp_p2.smiss | sort > tmp_p2_cts.txt
awk 'NR > 1 {print $2, $4}' tmp_19.imiss | sort > tmp_19_cts.txt
diff -q tmp_p2_cts.txt tmp_19_cts.txt
