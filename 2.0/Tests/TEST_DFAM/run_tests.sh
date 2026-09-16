#!/bin/bash

# --dfam, against an independent reference implementation and against
# PLINK 1.9.
#
# PLINK 2.0's statistic differs from PLINK 1.x's in one place: PLINK 1.07 and
# 1.9 add the magnitude of Cov(hom-ALT count, het count) to the sibship
# variance where it has to be subtracted, which inflates the variance and makes
# the statistic depend on which allele is counted.  ref_dfam.py implements both
# versions, so the tests below pin PLINK 2.0 to the corrected one and PLINK 1.9
# to the PLINK 1.x one, which together show that this is the only difference.

set -exo pipefail

python3 make_dfam.py 3

# 1. Against the reference implementation.
$1/plink2 $2 $3 --vcf dfam.vcf --psam dfam.psam --make-pgen --out tmp_data
$1/plink2 $2 $3 --pfile tmp_data --dfam --out plink2_d
python3 ref_dfam.py dfam.psam dfam.vcf > tmp_ref.txt
python3 compare_dfam.py plink2 plink2_d.dfam tmp_ref.txt 1e-5

# All four kinds of group have to be present, or the comparison above is
# weaker than it looks.
grep -q "50 informative families, 20 mixed sibships, 1 cluster of unrelateds" plink2_d.log

# 2. Against PLINK 1.9, through the same reference implementation in
#    PLINK 1.x-compatible mode.  PLINK 1.9 counts the minor allele, and its
#    statistic isn't invariant under that choice, so the orientation it
#    reports is fed back in.
plink --vcf dfam.vcf --keep-allele-order --make-bed --out tmp_p19
tail -n +2 dfam.psam > tmp_p19.fam
plink --bfile tmp_p19 --dfam --out plink19_d
awk 'NR > 1 { print $2, $3 }' plink19_d.dfam > tmp_a1.txt
python3 ref_dfam.py dfam.psam dfam.vcf plink19 tmp_a1.txt > tmp_ref19.txt
# PLINK 1.9 prints 4 significant digits.
python3 compare_dfam.py plink19 plink19_d.dfam tmp_ref19.txt 1e-3

# 3. Swapping REF and ALT must not change the statistic.  (It does in
#    PLINK 1.x; see the comment at the top.)
$1/plink2 $2 $3 --vcf dfam_flipped.vcf --psam dfam.psam --make-pgen --out tmp_flipped
$1/plink2 $2 $3 --pfile tmp_flipped --dfam --out plink2_flipped
awk 'NR > 1 { print $3, $6, $7 }' plink2_d.dfam > tmp_stats.txt
awk 'NR > 1 { print $3, $6, $7 }' plink2_flipped.dfam > tmp_stats_flipped.txt
# Cancellation can leave a statistic that is zero on one side and 1e-32 on the
# other, so this is a numeric comparison rather than a byte comparison.
paste tmp_stats.txt tmp_stats_flipped.txt | awk '{
    if ($1 != $4) { print "variant order differs at " $1; exit 1 }
    for (i = 2; i <= 3; i++) {
        a = $i; b = $(i + 3)
        if (a == "NA" || b == "NA") {
            if (a != b) { print "NA mismatch at " $1; exit 1 }
            continue
        }
        d = a - b; if (d < 0) d = -d
        m = (a < 0)? -a : a; if (m < 1e-12) m = 1e-12
        if (d > 1e-6 * m && d > 1e-12) {
            print "REF/ALT swap changed " $1 ": " a " vs " b
            exit 1
        }
    }
}'

# 4. chrX and chrMT are excluded.
grep -q "Excluding 2 haploid/MT variants" plink2_d.log
test "$(grep -cv '^#' plink2_d.dfam)" -eq 400

# 5. 'no-unrelateds' drops the unrelated component, so the statistic changes.
$1/plink2 $2 $3 --pfile tmp_data --dfam no-unrelateds --out plink2_nounrel
grep -q "0 clusters of unrelateds" plink2_nounrel.log
if diff -q plink2_d.dfam plink2_nounrel.dfam > /dev/null; then
    echo "no-unrelateds did not change the output"
    exit 1
fi

# 6. Column sets and zs.
$1/plink2 $2 $3 --pfile tmp_data --dfam cols=chrom,pos,ref,alt,obs,exp,chisq,p --out plink2_cols
head -1 plink2_cols.dfam | grep -qx '#CHROM	POS	ID	REF	ALT	OBS_CT	EXP_CT	CHISQ	P'
$1/plink2 $2 $3 --pfile tmp_data --dfam cols=chisq --out plink2_chisq
head -1 plink2_chisq.dfam | grep -qx '#ID	CHISQ'
$1/plink2 $2 $3 --pfile tmp_data --dfam zs --out plink2_zs
$1/plink2 $2 $3 --zst-decompress plink2_zs.dfam.zst > plink2_zs.dfam
diff -q plink2_d.dfam plink2_zs.dfam

# 7. A quantitative phenotype is rejected.
awk 'BEGIN { OFS = "\t"; print "#FID", "IID", "QT" }
NR > 1 { print $1, $2, NR * 0.25 }' dfam.psam > tmp_qt.txt
if $1/plink2 $2 $3 --pfile tmp_data --pheno tmp_qt.txt --no-psam-pheno --dfam --out plink2_bad 2> tmp_err.txt; then
    echo "expected --dfam to require a case/control phenotype"
    exit 1
fi
grep -q "case/control phenotype" tmp_err.txt

# 8. A dataset with no usable variants is an error.
$1/plink2 $2 $3 --pfile tmp_data --chr MT --make-pgen --out tmp_mt
if $1/plink2 $2 $3 --pfile tmp_mt --dfam --out plink2_bad2 2> tmp_err2.txt; then
    echo "expected --dfam to reject an MT-only dataset"
    exit 1
fi
grep -q "No variants remaining" tmp_err2.txt
