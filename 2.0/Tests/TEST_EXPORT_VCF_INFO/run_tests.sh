#!/bin/bash

# --export vcf/bcf 'vcf-info=AC,AN'.
# input.vcf covers biallelic and multiallelic variants, missing calls, phased
# and unphased calls, a variant with no ALT allele, and chrX (males,
# females, and an unknown-sex sample; including a male het), chrY, and MT.
# Its INFO column carries stale AC/AN entries (with a conflicting AC header
# line) which must be replaced rather than duplicated.
# check_ac_an.py recomputes AC and AN from the exported GT field, so the
# counts are checked against the ploidy each call was actually written with.

set -exo pipefail

$1/plink2 $2 $3 --vcf input.vcf --update-sex sex.txt --make-pgen --out tmp_data

$1/plink2 $2 $3 --pfile tmp_data --export vcf vcf-info=AC,AN --out tmp_vcf
python3 check_ac_an.py tmp_vcf.vcf AC,AN
grep -v '^#' tmp_vcf.vcf | cut -f 8 > tmp_info.txt
diff -q expected_info.txt tmp_info.txt
# The stale AC and AN header lines must be gone; the others must be kept.
test "$(grep -c '^##INFO=<ID=A[CN],' tmp_vcf.vcf)" -eq 2
grep -q '^##INFO=<ID=DP,' tmp_vcf.vcf
grep -q '^##INFO=<ID=DB,' tmp_vcf.vcf

# BCF: re-import, and verify that INFO decodes to the same text and that the
# genotypes are unchanged.
$1/plink2 $2 $3 --pfile tmp_data --export bcf vcf-info=AC,AN --out tmp_bcf
$1/plink2 $2 $3 --bcf tmp_bcf.bcf --update-sex sex.txt --make-pgen --out tmp_bcf_rt
grep -v '^#' tmp_bcf_rt.pvar | cut -f 6 > tmp_bcf_info.txt
diff -q expected_info.txt tmp_bcf_info.txt
diff -q tmp_data.pgen tmp_bcf_rt.pgen

# Sample subsets exercise the other ploidy code paths: all-male chrX (plain
# haploid), no males (plain diploid), and unphased data.
for subset in "--keep-males" "--keep-females" "--keep-nosex"
do
    $1/plink2 $2 $3 --pfile tmp_data $subset --export vcf vcf-info=AC,AN --out tmp_sub
    python3 check_ac_an.py tmp_sub.vcf AC,AN
    $1/plink2 $2 $3 --pfile tmp_data $subset --export bcf vcf-info=AC,AN --out tmp_sub
    $1/plink2 $2 $3 --bcf tmp_sub.bcf --update-sex sex.txt --make-pgen --out tmp_sub_rt
    diff -q <(grep -v '^#' tmp_sub.vcf | cut -f 8) <(grep -v '^#' tmp_sub_rt.pvar | cut -f 6)
done
sed 's/|/\//g' input.vcf > tmp_unphased.vcf
$1/plink2 $2 $3 --vcf tmp_unphased.vcf --update-sex sex.txt --make-pgen --out tmp_unphased
$1/plink2 $2 $3 --pfile tmp_unphased --export vcf vcf-info=AC,AN --out tmp_unphased
python3 check_ac_an.py tmp_unphased.vcf AC,AN
$1/plink2 $2 $3 --pfile tmp_unphased --export bcf vcf-info=AC,AN --out tmp_unphased
$1/plink2 $2 $3 --bcf tmp_unphased.bcf --update-sex sex.txt --make-pgen --out tmp_unphased_rt
diff -q <(grep -v '^#' tmp_unphased.vcf | cut -f 8) <(grep -v '^#' tmp_unphased_rt.pvar | cut -f 6)

# Swapped REF/ALT (--ref-allele), with the counts following the new ALT.
printf "v1\tG\nx1\tG\nm1\tG\n" > tmp_ref.txt
$1/plink2 $2 $3 --pfile tmp_data --ref-allele force tmp_ref.txt 2 1 --export vcf vcf-info=AC,AN --out tmp_swap
python3 check_ac_an.py tmp_swap.vcf AC,AN
test "$(awk -F '\t' '$3 == "v1" {print $8}' tmp_swap.vcf)" = "DP=10;AC=4;AN=10"

# Only one of the two keys: the other .pvar entry is left alone.
$1/plink2 $2 $3 --pfile tmp_data --export vcf vcf-info=AN --out tmp_an
python3 check_ac_an.py tmp_an.vcf AN
test "$(awk -F '\t' '$3 == "v1" {print $8}' tmp_an.vcf)" = "AC=7;DP=10;AN=10"
grep -q '^##INFO=<ID=AC,Number=1,Type=Float,' tmp_an.vcf
$1/plink2 $2 $3 --pfile tmp_data --export bcf vcf-info=AC --out tmp_ac
$1/plink2 $2 $3 --bcf tmp_ac.bcf --update-sex sex.txt --make-pgen --out tmp_ac_rt
test "$(awk -F '\t' '$3 == "v1" {print $6}' tmp_ac_rt.pvar)" = "DP=10;AN=99;AC=6"

# Larger counts need int16 and int32 BCF encodings.  (Only consistency is
# checked here, since the genotypes are random.)
$1/plink2 $2 $3 --dummy 17000 20 0.01 acgt --seed 1 --make-pgen --out tmp_big
$1/plink2 $2 $3 --pfile tmp_big --export vcf vcf-info=AC,AN --out tmp_big
python3 check_ac_an.py tmp_big.vcf AC,AN
$1/plink2 $2 $3 --pfile tmp_big --export bcf vcf-info=AC,AN --out tmp_big
$1/plink2 $2 $3 --bcf tmp_big.bcf --make-pgen --out tmp_big_rt
diff -q <(grep -v '^#' tmp_big.vcf | cut -f 8) <(grep -v '^#' tmp_big_rt.pvar | cut -f 6)

# Invalid combinations.
if $1/plink2 $2 $3 --pfile tmp_data --export ped vcf-info=AC --out tmp_bad 2> tmp_err.txt; then
    exit 1
fi
grep -q "only applies to --export's vcf and bcf" tmp_err.txt
if $1/plink2 $2 $3 --pfile tmp_data --export vcf vcf-info=AC,AF --out tmp_bad 2> tmp_err.txt; then
    exit 1
fi
grep -q "Invalid --export vcf-info= argument" tmp_err.txt
if $1/plink2 $2 $3 --pfile tmp_data --export vcf vcf-dosage=DS-only vcf-info=AN --out tmp_bad 2> tmp_err.txt; then
    exit 1
fi
grep -q "cannot be used with vcf-dosage=DS-only" tmp_err.txt
