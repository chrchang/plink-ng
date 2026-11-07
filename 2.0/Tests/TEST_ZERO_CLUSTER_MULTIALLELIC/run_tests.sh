#!/bin/bash

# --zero-cluster on multiallelic variants.
#
# Calls set to missing kept their ALT2+ patch entries, so they came back as
# nonmissing.

set -exo pipefail

plink2="$1/plink2 $2 $3"

cat > tmp_data.vcf <<'VCF'
##fileformat=VCFv4.2
##contig=<ID=1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5	S6
1	100	v1	A	C,G	.	.	.	GT	0/2	1/2	2/2	0/2	1/1	0/1
1	200	v2	A	T	.	.	.	GT	0/1	1/1	0/0	0/1	0/0	1/1
VCF
printf '#IID\tSEX\tCLUST\nS1\t1\tc1\nS2\t2\tc1\nS3\t1\tc1\nS4\t2\tc2\nS5\t1\tc2\nS6\t2\tc2\n' > tmp_data.psam
$plink2 --vcf tmp_data.vcf --psam tmp_data.psam --make-pgen --out tmp_data > /dev/null
printf 'v1\tc1\n' > tmp_zero.txt
$plink2 --pfile tmp_data --zero-cluster tmp_zero.txt --make-pgen --out tmp_zeroed > /dev/null
$plink2 --pfile tmp_zeroed --export vcf --out tmp_zeroed > /dev/null
[[ "$(awk '!/^#/ && $3 == "v1" {print $10, $11, $12, $13, $14, $15}' tmp_zeroed.vcf)" == "./. ./. ./. 0/2 1/1 0/1" ]]
[[ "$(awk '!/^#/ && $3 == "v2" {print $10, $11, $12, $13, $14, $15}' tmp_zeroed.vcf)" == "0/1 1/1 0/0 0/1 0/0 1/1" ]]
