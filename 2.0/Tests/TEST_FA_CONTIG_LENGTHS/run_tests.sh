#!/bin/bash

# Contig lengths scraped from --fa for VCF export.
#
# The last contig in the .fa file never had its length recorded.

set -exo pipefail

plink2="$1/plink2 $2 $3"

printf '>1\nACGTACGTACGTACGTACGT\n>2\nACGTACGTAC\nACG\n' > tmp_ref.fa
cat > tmp_data.vcf <<'VCF'
##fileformat=VCFv4.2
##contig=<ID=1>
##contig=<ID=2>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
1	3	v1	G	T	.	.	.	GT	0/1	0/0
2	3	v2	G	T	.	.	.	GT	0/1	1/1
VCF
$plink2 --vcf tmp_data.vcf --make-pgen --out tmp_data > /dev/null
$plink2 --pfile tmp_data --fa tmp_ref.fa --export vcf --out tmp_out > /dev/null
[[ "$(grep '^##contig' tmp_out.vcf | tr '\n' ' ')" == "##contig=<ID=1,length=20> ##contig=<ID=2,length=13> " ]] || exit 1
