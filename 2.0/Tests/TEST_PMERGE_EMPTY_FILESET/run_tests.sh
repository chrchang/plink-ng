#!/bin/bash

# --pmerge-list where a fileset has no variants left after filtering.
#
# Unlinking an emptied fileset also skipped the one after it.

set -exo pipefail

plink2="$1/plink2 $2 $3"

for c in 1 2 3 4; do
    cat > tmp_chr$c.vcf <<VCF
##fileformat=VCFv4.2
##contig=<ID=$c,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
$c	100	c${c}a	A	G	.	.	.	GT	0/1	0/0
$c	200	c${c}b	C	T	.	.	.	GT	1/1	0/1
VCF
    $plink2 --vcf tmp_chr$c.vcf --make-pgen --out tmp_chr$c > /dev/null
done
printf 'tmp_chr1\ntmp_chr2\ntmp_chr3\ntmp_chr4\n' > tmp_list.txt
$plink2 --pmerge-list tmp_list.txt --chr 1,4 --out tmp_merged > /dev/null
[[ "$(awk '!/^#/ {printf "%s ", $3}' tmp_merged.pvar)" == "c1a c1b c4a c4b " ]] || exit 1
