#!/bin/bash

# LD pruning and LD Score options with no other test coverage:
# --indep-pairphase, --indep-order, --indep-preferred, --ld-score-founders,
# and --ld-score-window-cm.
# Requires plink (1.9) to be in the system PATH.
#
# Regression check: --ld-score's default 1 cM window was stored as exactly
# 1.0, while an explicit --ld-score-window-cm value is stretched by a small
# epsilon, so variants exactly 1 cM away were dropped by the default window
# but kept by an explicit "--ld-score-window-cm 1".

set -exo pipefail

plink2="$1/plink2 $2 $3"

python3 make_data.py
$plink2 --vcf tmp_data.vcf --make-pgen --out tmp_data > /dev/null
$plink2 --pfile tmp_data --make-bed --out tmp_data > /dev/null

# --- --indep-order ---
# '1' replicates PLINK 1.x exactly; the default ('2') prunes differently
# with these parameters.
plink --bfile tmp_data --indep-pairwise 30 3 0.2 --out plink19_prune > /dev/null
$plink2 --pfile tmp_data --indep-pairwise 30 3 0.2 --indep-order 1 --out tmp_o1 > /dev/null
$plink2 --pfile tmp_data --indep-pairwise 30 3 0.2 --out tmp_o2 > /dev/null
$plink2 --pfile tmp_data --indep-pairwise 30 3 0.2 --indep-order 2 --out tmp_o2b > /dev/null
diff -q tmp_o1.prune.in plink19_prune.prune.in
diff -q tmp_o1.prune.out plink19_prune.prune.out
diff -q tmp_o2.prune.in tmp_o2b.prune.in
if diff -q tmp_o1.prune.in tmp_o2.prune.in > /dev/null; then
    exit 1
fi
if $plink2 --pfile tmp_data --indep-pairwise 30 3 0.2 --indep-order 3 --out tmp_bad > /dev/null 2>&1; then
    exit 1
fi

# --- --indep-preferred ---
# Four variants the default run prunes, far enough apart not to compete
# with each other, are all kept when preferred.  An unknown ID is ignored.
awk -F 'v' 'NR == 1 {print "v" $2; last = $2; next} $2 - last >= 60 {print "v" $2; last = $2}' tmp_o2.prune.out | head -n 4 > tmp_pref.txt
test "$(wc -l < tmp_pref.txt)" -eq 4
test "$(grep -cxFf tmp_pref.txt tmp_o2.prune.in || true)" -eq 0
for order in 1 2; do
    $plink2 --pfile tmp_data --indep-pairwise 30 3 0.2 --indep-order $order --indep-preferred tmp_pref.txt --out tmp_pref$order > /dev/null
    test "$(grep -cxFf tmp_pref.txt tmp_pref$order.prune.in)" -eq 4
    # prune.in and prune.out still partition the variants.
    test "$(cat tmp_pref$order.prune.in tmp_pref$order.prune.out | sort -u | wc -l)" -eq 300
done
{ cat tmp_pref.txt; echo nosuch; } > tmp_pref_extra.txt
$plink2 --pfile tmp_data --indep-pairwise 30 3 0.2 --indep-preferred tmp_pref_extra.txt --out tmp_pref_extra > /dev/null
diff -q tmp_pref2.prune.in tmp_pref_extra.prune.in

# --- --indep-pairphase ---
# v5swap has the same genotypes as v5 (genotype r^2 = 1), but opposite phase
# at every het, so haplotype r^2 is low: --indep-pairwise prunes one of the
# two, --indep-pairphase keeps both.
$plink2 --vcf tmp_swap.vcf --make-pgen --out tmp_swap > /dev/null
$plink2 --pfile tmp_swap --indep-pairwise 10 1 0.5 --out tmp_pw > /dev/null
$plink2 --pfile tmp_swap --indep-pairphase 10 1 0.5 --out tmp_pp > /dev/null
test "$(grep -cxE 'v5|v5swap' tmp_pw.prune.in)" -eq 1
test "$(grep -cxE 'v5|v5swap' tmp_pp.prune.in)" -eq 2
$plink2 --pfile tmp_swap --indep-pairphase 10kb 0.5 --out tmp_ppkb > /dev/null
test "$(grep -cxE 'v5|v5swap' tmp_ppkb.prune.in)" -eq 2
# It requires phased genotypes.
if $plink2 --bfile tmp_data --indep-pairphase 10 1 0.5 --out tmp_bad > /dev/null 2>&1; then
    exit 1
fi

# --- --ld-score-founders / --ld-score-window-cm ---
# cM = bp / 10000, so 1 cM = 10 kb; every fifth sample is a nonfounder.
awk 'BEGIN {OFS = "\t"} /^##/ {print; next} /^#CHROM/ {print $1, $2, $3, $4, $5, "CM"; next} {print $1, $2, $3, $4, $5, $2 / 10000}' tmp_data.pvar > tmp_cm.pvar
awk 'BEGIN {OFS = "\t"} NR == 1 {print "#IID", "PAT", "MAT", "SEX"; next} {i = NR - 2; if (i % 5 == 4) {print $1, "S" (i - 4), "S" (i - 3), "NA"} else {print $1, "0", "0", "NA"}}' tmp_data.psam > tmp_cm.psam
cp tmp_data.pgen tmp_cm.pgen
# $1: output prefix; remaining args: flags.  Prints the scores' hash.
ldscore_hash() {
    local out=$1
    shift
    $plink2 --pfile tmp_cm --ld-score cols=chrom,pos,nobsi,l2 "$@" --out $out > /dev/null
    tail -n +2 $out.ldscore | cksum
}
# Nonfounders must be handled explicitly.
if $plink2 --pfile tmp_cm --ld-score --ld-score-window-kb 20 --out tmp_bad > /dev/null 2>&1; then
    exit 1
fi
founders_kb20=$(ldscore_hash tmp_f20 --ld-score-founders --ld-score-window-kb 20)
test "$(ldscore_hash tmp_kf20 --keep-founders --ld-score-window-kb 20)" = "$founders_kb20"
test "$(ldscore_hash tmp_nf20 --nonfounders --ld-score-window-kb 20)" != "$founders_kb20"
# cM windows match the equivalent kb windows.
test "$(ldscore_hash tmp_cm2 --ld-score-founders --ld-score-window-cm 2)" = "$founders_kb20"
founders_kb10=$(ldscore_hash tmp_f10 --ld-score-founders --ld-score-window-kb 10)
test "$(ldscore_hash tmp_cm1 --ld-score-founders --ld-score-window-cm 1)" = "$founders_kb10"
# The default window is 1 cM, including variants exactly 1 cM away.
test "$(ldscore_hash tmp_def --ld-score-founders)" = "$founders_kb10"
test "$(awk 'NR == 12 {print $4}' tmp_def.ldscore)" -eq 21
# With several limits, the narrowest applies.
test "$(ldscore_hash tmp_kb20cm1 --ld-score-founders --ld-score-window-kb 20 --ld-score-window-cm 1)" = "$founders_kb10"
test "$(ldscore_hash tmp_kb5cm1 --ld-score-founders --ld-score-window-kb 5 --ld-score-window-cm 1)" = "$(ldscore_hash tmp_kb5 --ld-score-founders --ld-score-window-kb 5)"
test "$(ldscore_hash tmp_w7cm1 --ld-score-founders --ld-score-window 7 --ld-score-window-cm 1)" = "$(ldscore_hash tmp_w7 --ld-score-founders --ld-score-window 7)"
# Without cM coordinates, the default window is an error.
if $plink2 --pfile tmp_data --ld-score --out tmp_bad > /dev/null 2>&1; then
    exit 1
fi
if $plink2 --pfile tmp_cm --ld-score --ld-score-founders --ld-score-window-cm -1 --out tmp_bad > /dev/null 2>&1; then
    exit 1
fi
