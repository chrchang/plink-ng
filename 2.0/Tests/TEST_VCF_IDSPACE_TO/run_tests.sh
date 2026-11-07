#!/bin/bash

# --idspace-to on VCF import without --id-delim.
#
# Spaces in VCF sample IDs were only checked/converted when an ID delimiter
# was in effect, so without one they ended up in the .psam unconverted.

set -exo pipefail

plink2="$1/plink2 $2 $3"

printf '##fileformat=VCFv4.2\n##contig=<ID=1,length=1000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample one\tsample two\n1\t100\tv1\tA\tG\t.\t.\t.\tGT\t0/1\t0/0\n' > tmp_data.vcf
$plink2 --vcf tmp_data.vcf --idspace-to _ --make-pgen --out tmp_data > /dev/null
[[ "$(tail -n +2 tmp_data.psam | cut -f 1 | tr '\n' ' ')" == "sample_one sample_two " ]] || exit 1
# Without --idspace-to, this is an error.
if $plink2 --vcf tmp_data.vcf --make-pgen --out tmp_err > /dev/null 2>&1; then
    echo "space in sample ID accepted"
    exit 1
fi
