#!/bin/bash

# --qfam / --qfam-parents / --qfam-between / --qfam-total, against PLINK 1.9.
#
# BETA, STAT, RAW_P and NIND are deterministic and are compared exactly (to
# PLINK 1.9's 4-significant-digit output precision).  EMP1 comes from a
# different random stream in each program, so it is only checked for agreement
# within Monte Carlo error.

set -exo pipefail

python3 make_qfam.py 5

plink --vcf qfam.vcf --keep-allele-order --make-bed --out tmp_p19
tail -n +2 qfam.psam > tmp_p19.fam
$1/plink2 $2 $3 --vcf qfam.vcf --psam qfam.psam --make-pgen --out tmp_data

# 1. The four tests against PLINK 1.9, with a fixed permutation count so that
#    the EMP1 comparison has a known Monte Carlo error (sd ~ 0.05 at p = 0.5
#    for 200 permutations of each program).
for spec in ":within" "-parents:parents" "-between:between" "-total:total"; do
    flag_suffix=${spec%%:*}
    out_suffix=${spec#*:}
    plink --bfile tmp_p19 --qfam$flag_suffix --mperm 200 --seed 7 \
        --out plink19_$out_suffix
    $1/plink2 $2 $3 --pfile tmp_data --qfam$flag_suffix mperm=200 --seed 7 \
        --out plink2_$out_suffix
    python3 compare_qfam.py \
        plink19_$out_suffix.qfam.$out_suffix \
        plink19_$out_suffix.qfam.$out_suffix.perm \
        plink2_$out_suffix.qfam.$out_suffix 1e-3 0.02
done

# All four kinds of group have to be present, or the comparison is weaker than
# it looks: 50 families, 20 sibships, and the rest singletons.
grep -q "Permuting 130 families/sibships/singletons" plink2_within.log

# 2. chrX and chrMT are excluded.
grep -q "Excluding 2 haploid/MT variants" plink2_within.log
test "$(grep -cv '^#' plink2_within.qfam.within)" -eq 400

# 3. Adaptive permutation, which is the default, stops early on most variants
#    and keeps going on a few.
$1/plink2 $2 $3 --pfile tmp_data --qfam --seed 7 --out plink2_adapt
awk 'NR > 1 { print $11 }' plink2_adapt.qfam.within | sort -n > tmp_np.txt
test "$(head -1 tmp_np.txt)" -lt 20
test "$(tail -1 tmp_np.txt)" -gt 100
# The deterministic columns must not depend on the permutation settings.
awk 'NR > 1 { print $3, $6, $7, $8, $9 }' plink2_within.qfam.within > tmp_fixed.txt
awk 'NR > 1 { print $3, $6, $7, $8, $9 }' plink2_adapt.qfam.within > tmp_fixed2.txt
diff -q tmp_fixed.txt tmp_fixed2.txt

# 4. 'perm-count' reports counts, and they line up with the p-values.
$1/plink2 $2 $3 --pfile tmp_data --qfam mperm=200 perm-count --seed 7 \
    --out plink2_count
paste plink2_within.qfam.within plink2_count.qfam.within | awk 'NR > 1 {
    if ($3 != $14) { print "variant order differs at " $3; exit 1 }
    if ($10 == "NA") { next }
    # EMP1 = (2 * count + 2) / (2 * (NP + 1))
    want = (2 * $21 + 2) / (2 * ($11 + 1))
    d = $10 - want; if (d < 0) d = -d
    if (d > 1e-6) { print "count/p-value mismatch at " $3 ": " $10 " vs " want; exit 1 }
}'

# 5. 'emp-se' adds two columns, and 'zs' round-trips.
$1/plink2 $2 $3 --pfile tmp_data --qfam mperm=200 emp-se --seed 7 \
    --out plink2_empse
head -1 plink2_empse.qfam.within | grep -q 'EMP_BETA	EMP_SE'
$1/plink2 $2 $3 --pfile tmp_data --qfam mperm=200 zs --seed 7 --out plink2_zs
$1/plink2 $2 $3 --zst-decompress plink2_zs.qfam.within.zst > plink2_zs.qfam.within
diff -q plink2_within.qfam.within plink2_zs.qfam.within

# 6. Column sets.
$1/plink2 $2 $3 --pfile tmp_data --qfam mperm=200 cols=chrom,pos,ref,alt,test,beta,stat,emp1 --seed 7 --out plink2_cols
head -1 plink2_cols.qfam.within | grep -qx '#CHROM	POS	ID	REF	ALT	TEST	BETA	STAT	EMP1'

# 7. Error paths.
awk 'BEGIN { OFS = "\t"; print "#FID", "IID", "CC" }
NR > 1 { print $1, $2, 1 + (NR % 2) }' qfam.psam > tmp_cc.txt
if $1/plink2 $2 $3 --pfile tmp_data --pheno tmp_cc.txt --no-psam-pheno \
       --qfam mperm=20 --out plink2_bad 2> tmp_err.txt; then
    echo "expected --qfam to require a quantitative phenotype"
    exit 1
fi
grep -q "quantitative phenotype" tmp_err.txt

$1/plink2 $2 $3 --pfile tmp_data --chr MT --make-pgen --out tmp_mt
if $1/plink2 $2 $3 --pfile tmp_mt --qfam mperm=20 --out plink2_bad2 2> tmp_err2.txt; then
    echo "expected --qfam to reject an MT-only dataset"
    exit 1
fi
grep -q "No variants remaining" tmp_err2.txt

if $1/plink2 $2 $3 --pfile tmp_data --qfam --qfam-total --out plink2_bad3 2> tmp_err3.txt; then
    echo "expected two --qfam flags to be rejected"
    exit 1
fi
grep -q "Only one --qfam" tmp_err3.txt
