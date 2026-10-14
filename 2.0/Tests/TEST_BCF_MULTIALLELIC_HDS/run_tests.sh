#!/bin/bash

# "--export bcf vcf-dosage=HDS-force" on a fileset with multiallelic variants
# and no phased calls used to dereference a null phasepresent pointer and
# crash.

set -exo pipefail

cat > multi.vcf << 'EOF'
##fileformat=VCFv4.3
##contig=<ID=1,length=1000000>
##contig=<ID=X,length=100000000>
##contig=<ID=Y,length=1000000>
##contig=<ID=MT,length=20000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	m	f	u
1	100	v1	C	G,T	.	.	.	GT	0/1	1/2	2/2
1	200	v2	C	G	.	.	.	GT	0/1	1/1	./.
X	5000000	v3	A	C,G	.	.	.	GT	2	0/2	1/2
Y	100	v4	A	C,G	.	.	.	GT	2	.	.
MT	100	v5	A	C,G,T	.	.	.	GT	3	1	.
EOF
printf '#IID\tSEX\nm\t1\nf\t2\nu\tNA\n' > multi.psam

$1/plink2 $2 $3 --vcf multi.vcf --psam multi.psam --make-pgen --out tmp_data
$1/plink2 $2 $3 --pfile tmp_data --export vcf --out tmp_want
for mode in DS-force HDS-force; do
    $1/plink2 $2 $3 --pfile tmp_data --export bcf vcf-dosage=$mode --out tmp_bcf
    $1/plink2 $2 $3 --bcf tmp_bcf.bcf --psam multi.psam --export vcf --out tmp_got
    diff <(grep -v '^##' tmp_want.vcf) <(grep -v '^##' tmp_got.vcf)
done
