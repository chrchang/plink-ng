#!/bin/bash

# With vcf-dosage=DS-force or DS-only, unphased diploid multiallelic 0/x
# calls (x >= 2) used to get an extra HDS subfield that FORMAT doesn't
# declare, e.g. '0/2:0,1:0,0.5,0,0.5' under FORMAT 'GT:DS'.

set -exo pipefail

cat > multi.vcf << 'EOF'
##fileformat=VCFv4.3
##contig=<ID=1,length=1000000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	a	b	c	d
1	100	v1	C	G,T	.	.	.	GT	0/1	1/2	0/2	2/2
1	200	v2	C	G,T,A	.	.	.	GT	0/3	0/0	./.	1/3
EOF

$1/plink2 $2 $3 --vcf multi.vcf --make-pgen --out tmp_data
for mode in DS-force DS-only HDS-force; do
    $1/plink2 $2 $3 --pfile tmp_data --export vcf vcf-dosage=$mode --out tmp_out
    # Every sample column must have exactly as many subfields as FORMAT.
    grep -v '^#' tmp_out.vcf | awk -F'\t' '{
        n_fmt = split($9, fmt, ":")
        for (i = 10; i <= NF; i++) {
            if (split($i, sub_fields, ":") != n_fmt) {
                exit 1
            }
        }
    }'
done
