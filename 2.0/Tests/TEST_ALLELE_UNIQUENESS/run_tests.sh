#!/bin/bash

set -exo pipefail

# Commands which depend on allele uniqueness (--ref-allele, --score, etc.)
# reject a variant with a repeated allele code.  The repeated codes need not
# be adjacent to the alphabetically first allele: "A C,C" sorts to A, C, C.

for alleles in "A C" "A A" "A C,C" "G T,C,T" "AT C,G,AT"; do
    printf '##fileformat=VCFv4.2\n##contig=<ID=1>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n1\t100\tv1\tA\tC\t.\t.\t.\tGT\t0/1\n1\t200\tv2\t%s\t%s\t.\t.\t.\tGT\t0/1\n' $alleles > tmp_data.vcf
    $1/plink2 $2 $3 --vcf tmp_data.vcf --make-pgen --out tmp_data
    printf 'v1\tC\n' > tmp_ref.txt
    if [ "$alleles" = "A C" ]; then
        $1/plink2 $2 $3 --pfile tmp_data --ref-allele force tmp_ref.txt 2 1 --make-pgen --out tmp_out
    else
        if $1/plink2 $2 $3 --pfile tmp_data --ref-allele force tmp_ref.txt 2 1 --make-pgen --out tmp_out; then
            echo "duplicate alleles $alleles were accepted"
            exit 1
        fi
        grep -q "Duplicate allele code in variant 'v2'" tmp_out.log
    fi
done
