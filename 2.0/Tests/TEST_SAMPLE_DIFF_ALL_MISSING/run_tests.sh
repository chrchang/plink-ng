#!/bin/bash

# --sample-diff counts on variants where every call is missing.
#
# Such variants were skipped as "all genotypes equal", so they were counted
# as observed and identical.

set -exo pipefail

plink2="$1/plink2 $2 $3"

# S1 vs S2: v1, v4, v5 are both-nonmissing (OBS_CT 3); v4 is IBS0, v5 IBS1;
# v4 and v5 differ.
cat > tmp_data.vcf <<'VCF'
##fileformat=VCFv4.2
##contig=<ID=1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
1	100	v1	A	G	.	.	.	GT	0/0	0/0	0/1
1	200	v2	A	G	.	.	.	GT	0/1	./.	0/1
1	300	v3	A	G	.	.	.	GT	./.	./.	./.
1	400	v4	A	G	.	.	.	GT	0/0	1/1	0/1
1	500	v5	A	G	.	.	.	GT	0/1	0/0	0/1
1	600	v6	A	G	.	.	.	GT	./.	0/1	0/1
VCF
$plink2 --vcf tmp_data.vcf --make-pgen --out tmp_data > /dev/null

expected=$(printf 'S1\tS2\t3\t3\t1\t1\t2')
$plink2 --pfile tmp_data --sample-diff counts-only counts-cols=+nobs,+nobsibs,+ibs0,+ibs1 ids=S1 S2 --out tmp_counts > /dev/null
[[ "$(tail -n 1 tmp_counts.sdiff.summary)" == "$expected" ]] || exit 1
$plink2 --pfile tmp_data --sample-diff counts-cols=+nobs,+nobsibs,+ibs0,+ibs1 ids=S1 S2 --out tmp_main > /dev/null
[[ "$(tail -n 1 tmp_main.sdiff.summary)" == "$expected" ]] || exit 1
