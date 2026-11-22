#!/bin/bash

set -exo pipefail

# --validate rejects a variant whose REF equals an ALT, or that repeats an ALT.

for alleles in "A A" "AT C,AT,G"; do
    printf '##fileformat=VCFv4.2\n##contig=<ID=1>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n1\t100\tv1\tA\tC\t.\t.\t.\tGT\t0/1\n1\t200\tv2\t%s\t%s\t.\t.\t.\tGT\t0/1\n' $alleles > tmp_dup.vcf
    $1/plink2 $2 $3 --vcf tmp_dup.vcf --make-pgen --out tmp_dup
    if $1/plink2 $2 $3 --pfile tmp_dup --validate --out tmp_val; then
        echo "--validate accepted duplicate alleles $alleles"
        exit 1
    fi
    grep -q "Duplicate allele code in variant 'v2'" tmp_val.log
done

printf '##fileformat=VCFv4.2\n##contig=<ID=1>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n1\t100\tv1\tA\tC\t.\t.\t.\tGT\t0/1\n1\t200\tv2\tAT\tC,G\t.\t.\t.\tGT\t1/2\n' > tmp_ok.vcf
$1/plink2 $2 $3 --vcf tmp_ok.vcf --make-pgen --out tmp_ok
$1/plink2 $2 $3 --pfile tmp_ok --validate --out tmp_val
