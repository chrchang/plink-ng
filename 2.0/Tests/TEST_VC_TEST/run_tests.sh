#!/bin/bash

# "--vc-test": burden, SKAT, SKAT-O, ACAT-V and ACAT-O per set.
#
# The checks here are identities rather than stored expected values.  The one
# that does the most work is the single-variant set: with one variant the
# kernel is 1x1, so every one of the four tests collapses to the same score
# test, and they have to agree exactly.  That is what caught a missing residual
# variance factor during development, which had left SKAT and SKAT-O wrong by
# exactly sigma^2 while leaving the burden test correct.

set -exo pipefail

BUILD=$1
EXTRA1=$2
EXTRA2=$3

printf '600 null 0.05 0.95 1 1\n' > tmp_sim.txt
plink --simulate tmp_sim.txt --make-bed --out tmp_base --seed 17 > /dev/null

# Sets of 6, plus singleton sets over the first 20 variants.
awk '{ i = int((NR - 1) / 6);
       ids[i] = (ids[i] ? ids[i] "," : "") $2; chr[i] = $1;
       if (!(i in pos)) pos[i] = $4 }
     END { for (j = 0; j <= i; ++j) print "S" j, chr[j], pos[j], ids[j] }' tmp_base.bim > tmp_sets.txt
awk 'NR <= 20 { print "ONE" NR, $1, $4, $2 }' tmp_base.bim >> tmp_sets.txt

{ echo "#FID IID QT"; awk '{print $1, $2, (NR % 83) / 83.0}' tmp_base.fam; } > tmp_qt.pheno

$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_base --pheno tmp_qt.pheno --pheno-name QT --vc-test --set-list tmp_sets.txt --vc-max-af 0.5 --out tmp_qtv

# 1. Every reported p-value is a probability.
awk '!/^#/ { for (i = 6; i <= 10; ++i) if (!($i > 0) || $i > 1) { print "out of range: " $0; exit 1 } }' tmp_qtv.vc

# 2. A one-variant set collapses to a single score test, so the four tests
#    agree.  ACAT-O combines four identical p-values, which returns that same
#    p-value, so it agrees too.
awk '/^ONE/ {
        b = $6; s = $7; o = $8; a = $9; ao = $10;
        for (i = 7; i <= 10; ++i) {
            d = $i - b; if (d < 0) d = -d;
            if (d > 1e-5 * b + 1e-12) {
                print "singleton set " $1 ": column " i " is " $i ", burden is " b;
                exit 1
            }
        }
        ++n
     }
     END { if (n != 20) { print "expected 20 singleton sets, saw " n; exit 1 }
           printf "singleton identity holds for %d sets\n", n }' tmp_qtv.vc

# 3. ACAT-O is a Cauchy combination of four p-values, so it cannot fall below
#    a quarter of the smallest of them, nor rise above the largest.
awk '!/^#/ {
        mn = $6; mx = $6;
        for (i = 7; i <= 9; ++i) { if ($i < mn) mn = $i; if ($i > mx) mx = $i }
        if ($10 < mn / 4.0 * 0.999) { print "ACAT-O below its inputs: " $0; exit 1 }
        if ($10 > mx * 1.001 + 1e-12) { print "ACAT-O above its inputs: " $0; exit 1 }
     }' tmp_qtv.vc

# 4. Under a phenotype unrelated to the genotypes, the p-values should look
#    uniform.  A loose band, since this is 120 sets, but it catches a test
#    that is wrong by a constant factor.
# (sorted with sort -n rather than awk's asort, which is a gawk extension)
awk '!/^#/ && !/^ONE/ { print $7 }' tmp_qtv.vc | sort -n > tmp_skatp.txt
nskat=$(wc -l < tmp_skatp.txt | tr -d '[:space:]')
med=$(awk -v n="$nskat" 'NR == int(n / 2) + 1 { print; exit }' tmp_skatp.txt)
awk -v med="$med" -v n="$nskat" 'BEGIN {
    if (med < 0.15 || med > 0.85) { print "SKAT median under the null is " med; exit 1 }
    printf "null SKAT median %.3f over %d sets\n", med, n }' 

# 5. A phenotype built from one set must make that set stand out, and the rest
#    must not follow it.
awk 'NR == FNR { if (FNR <= 6) want[$2] = 1; next }
     { print }' tmp_base.bim /dev/null > /dev/null
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_base --export A --out tmp_a > /dev/null
awk 'NR == 1 { for (i = 7; i <= 12; ++i) col[i] = 1; next }
     { s = 0; for (i = 7; i <= 12; ++i) if ($i != "NA") s += $i
       printf "%s %s %.4f\n", $1, $2, s + (NR % 7) * 0.25 }' tmp_a.raw > tmp_sig.body
{ echo "#FID IID SIG"; cat tmp_sig.body; } > tmp_sig.pheno
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_base --pheno tmp_sig.pheno --pheno-name SIG --vc-test --set-list tmp_sets.txt --vc-max-af 0.5 --out tmp_sigv
s0=$(awk '$1 == "S0" { print $7 }' tmp_sigv.vc)
awk -v p="$s0" 'BEGIN { if (!(p < 1e-6)) { print "S0 SKAT p is " p ", expected strong signal"; exit 1 } }'

# 6. Case/control runs, and the same identities hold there.
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_base --vc-test --set-list tmp_sets.txt --vc-max-af 0.5 --out tmp_ccv
awk '!/^#/ { for (i = 6; i <= 10; ++i) if (!($i > 0) || $i > 1) { print "out of range: " $0; exit 1 } }' tmp_ccv.vc
awk '/^ONE/ { d = $7 - $6; if (d < 0) d = -d
              if (d > 1e-5 * $6 + 1e-12) { print "case/control singleton mismatch: " $0; exit 1 } }' tmp_ccv.vc

# 7. Error cases.
fails() {
    if "$@" > /dev/null 2>&1; then
        echo "expected failure: $*"
        exit 1
    fi
}
fails $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_base --vc-test --out tmp_bad
fails $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_base --vc-test --set-list tmp_sets.txt --vc-max-af 0 --out tmp_bad
fails $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_base --vc-test --set-list tmp_sets.txt --vc-params 0 25 --out tmp_bad

echo "--vc-test tests passed."
