#!/bin/bash

set -exo pipefail

# A variant whose REF equals an ALT, or with a repeated ALT, is rejected when
# the .pvar/.bim is loaded.  Missing codes ('.', or '0' in a .bim) may repeat.

dup_vcf() {
    printf '##fileformat=VCFv4.2\n##contig=<ID=1>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n1\t100\tv1\tA\tC\t.\t.\t.\tGT\t0/1\n1\t200\tv2\t%s\t%s\t.\t.\t.\tGT\t0/1\n' "$1" "$2" > tmp_dup.vcf
}

for alleles in "A A" "A C,C" "AT C,AT"; do
    dup_vcf $alleles
    if $1/plink2 $2 $3 --vcf tmp_dup.vcf --make-pgen --out tmp_dup; then
        echo "duplicate alleles $alleles were accepted"
        exit 1
    fi
    grep -q "Duplicate allele code" tmp_dup.log
done

for alleles in "A C,G" ". ." "A ."; do
    dup_vcf $alleles
    $1/plink2 $2 $3 --vcf tmp_dup.vcf --make-pgen --out tmp_ok
done

printf 'f\ti\t0\t0\t1\t1\n' > tmp_dup.fam
printf 'l\x1b\x01\x02' > tmp_dup.bed
printf '1\tv1\t0\t100\t0\t0\n' > tmp_dup.bim
$1/plink2 $2 $3 --bfile tmp_dup --freq --out tmp_ok
printf '1\tv1\t0\t100\tG\tG\n' > tmp_dup.bim
if $1/plink2 $2 $3 --bfile tmp_dup --freq --out tmp_dup; then
    echo "duplicate .bim alleles were accepted"
    exit 1
fi
grep -q "Duplicate allele code" tmp_dup.log
