#!/bin/bash

# --tucc, against PLINK 1.9 and against an independently computed answer.

set -exo pipefail

python3 make_trios.py 1

# 1. The pseudo-case/pseudo-control genotypes.  The expected values are
#    derived from the genotypes the fixture generator drew, not from another
#    PLINK run.
$1/plink2 $2 $3 --vcf trios.vcf --psam trios.psam --make-pgen --out tmp_data
$1/plink2 $2 $3 --pfile tmp_data --tucc --out plink2_t
awk 'NR > 1 { print $3, $5 }' plink2_t.tucc.pvar > tmp_alleles.txt
$1/plink2 $2 $3 --pfile plink2_t.tucc --export A --export-allele tmp_alleles.txt --out plink2_raw
python3 check_tucc.py plink2_raw.raw expected.txt

# 2. The .psam: two samples per trio, the child's FID/sex, '_T'/'_U' suffixes,
#    case then control.  The two families whose child isn't part of a full
#    trio must be absent.
awk 'BEGIN { OFS = "\t"; print "#FID", "IID", "SEX", "PHENO1" }
NR > 1 && $3 != "0" && $4 != "0" {
    print $1, $2 "_T", $5, 2
    print $1, $2 "_U", $5, 1
}' trios.psam > expected_psam.txt
diff -q expected_psam.txt plink2_t.tucc.psam
test "$(grep -c '' plink2_t.tucc.psam)" -eq 77

# 3. Chromosome and multiallelic filtering: chrX, chrMT and the multiallelic
#    variant are dropped, everything else is kept, in input order.
awk 'NR > 1 { print $3 }' plink2_t.tucc.pvar > tmp_kept.txt
awk '!/^#/ && $1 == "1" && $5 == "G" { print $3 }' trios.vcf > expected_kept.txt
diff -q expected_kept.txt tmp_kept.txt

# 4. Against PLINK 1.9.  It cannot read the multiallelic variant, so that one
#    is dropped on import; --tucc drops it on the plink2 side too, so the two
#    datasets still line up.  A .ped comparison needs the two heterozygote
#    allele orders normalized, since A1/A2 assignment differs.
plink --vcf trios.vcf --biallelic-only strict --make-bed --out tmp_p19
tail -n +2 trios.psam > tmp_p19.fam
plink --bfile tmp_p19 --tucc write-bed --out plink19_t
plink --bfile plink19_t.tucc --recode --out plink19_ped
$1/plink2 $2 $3 --pfile plink2_t.tucc --export ped --out plink2_ped
normalize() {
    awk '{
        printf "%s %s %s %s %s %s", $1, $2, $3, $4, $5, $6
        for (i = 7; i < NF; i += 2) {
            if ($i <= $(i + 1)) { printf " %s %s", $i, $(i + 1) }
            else { printf " %s %s", $(i + 1), $i }
        }
        printf "\n"
    }' "$1" > "$2"
}
normalize plink19_ped.ped plink19_norm.ped
normalize plink2_ped.ped plink2_norm.ped
diff -q plink19_norm.ped plink2_norm.ped
# Non-vacuous: the fileset has to have real genotypes in it.
test "$(grep -c '' plink2_norm.ped)" -eq 76
grep -q ' G ' plink2_norm.ped

# 5. 'vzs' compresses the .pvar, and nothing else changes.
$1/plink2 $2 $3 --pfile tmp_data --tucc vzs --out plink2_vzs
$1/plink2 $2 $3 --zst-decompress plink2_vzs.tucc.pvar.zst > plink2_vzs.tucc.pvar
diff -q plink2_t.tucc.pvar plink2_vzs.tucc.pvar
cmp plink2_t.tucc.pgen plink2_vzs.tucc.pgen
diff -q plink2_t.tucc.psam plink2_vzs.tucc.psam

# 6. No trios: a warning, not an error, and no output fileset.
awk 'BEGIN { OFS = "\t" } NR == 1 { print; next } { $3 = "0"; $4 = "0"; print }' trios.psam > tmp_nofam.psam
$1/plink2 $2 $3 --pfile tmp_data --psam tmp_nofam.psam --tucc --out plink2_notrio 2> tmp_notrio.txt
grep -q "no trios" tmp_notrio.txt
test ! -e plink2_notrio.tucc.pgen

# 7. No usable variants is an error.
$1/plink2 $2 $3 --pfile tmp_data --chr MT --make-pgen --out tmp_mt
if $1/plink2 $2 $3 --pfile tmp_mt --tucc --out plink2_bad 2> tmp_err.txt; then
    echo "expected --tucc to reject an MT-only dataset"
    exit 1
fi
grep -q "No variants remaining" tmp_err.txt
