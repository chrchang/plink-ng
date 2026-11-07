#!/bin/bash

# --flip permissive + --flip-subset with a non-SNP in the --flip list.
#
# The non-SNP was skipped when flipping alleles, but stayed marked for the
# sample-subset genotype rewrite.

set -exo pipefail

plink2="$1/plink2 $2 $3"

cat > tmp_data.vcf <<'VCF'
##fileformat=VCFv4.2
##contig=<ID=1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
1	100	v1	A	T	.	.	.	GT	0/0	0/1	1/1	0/0
1	200	v2	A	AT	.	.	.	GT	0/0	0/1	1/1	0/1
1	300	v3	G	C	.	.	.	GT	0/1	1/1	0/0	0/0
VCF
$plink2 --vcf tmp_data.vcf --make-pgen --out tmp_data > /dev/null
printf 'v1\nv2\nv3\n' > tmp_flip.txt
printf 'S1\nS2\n' > tmp_subset.txt
$plink2 --pfile tmp_data --flip tmp_flip.txt permissive --flip-subset tmp_subset.txt --make-pgen --out tmp_flipped > /dev/null
$plink2 --pfile tmp_flipped --export vcf --out tmp_flipped > /dev/null
# The indel's calls must be unchanged, and the A/T and C/G SNPs' calls are
# swapped for S1 and S2 only.
[[ "$(awk '!/^#/ && $3 == "v2" {print $10, $11, $12, $13}' tmp_flipped.vcf)" == "0/0 0/1 1/1 0/1" ]] || exit 1
[[ "$(awk '!/^#/ && $3 == "v1" {print $10, $11, $12, $13}' tmp_flipped.vcf)" == "1/1 0/1 1/1 0/0" ]] || exit 1
[[ "$(awk '!/^#/ && $3 == "v3" {print $10, $11, $12, $13}' tmp_flipped.vcf)" == "0/1 0/0 0/0 0/0" ]] || exit 1
