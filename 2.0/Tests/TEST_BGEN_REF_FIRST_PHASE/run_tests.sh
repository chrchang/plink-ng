#!/bin/bash

# Phased BGEN round trip with 'ref-first' on both sides.  --bgen ref-first
# used to swap the haplotypes of every phased het (0|1 came back as 1|0).

set -exo pipefail

# Print GT, and HDS rounded to 2 decimal places (BGEN is lossy past that).
summarize() {
    grep -v '^#' $1 | cut -f 10- | tr '\t' '\n' | awk -F: '{
        line = $1
        if (NF == 3) {
            n = split($3, hds, ",")
            for (i = 1; i <= n; i++) {
                line = line sprintf(" %.2f", hds[i])
            }
        }
        print line
    }'
}

cat > phased.vcf << 'EOF'
##fileformat=VCFv4.3
##contig=<ID=1,length=1000000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=HDS,Number=.,Type=Float,Description="Haplotype dosage">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	a	b	c	d
1	100	v1	A	G	.	.	.	GT:HDS	0|1:.	1|0:.	1|1:.	0/1:.
1	200	v2	C	T	.	.	.	GT:HDS	0|1:0,0.8	0|0:.	1|0:0.95,0.1	./.:.
EOF

$1/plink2 $2 $3 --vcf phased.vcf dosage=HDS --make-pgen --out tmp_data
$1/plink2 $2 $3 --pfile tmp_data --export vcf vcf-dosage=HDS --out tmp_want
for ver in 1.2 1.3; do
    for order in ref-first ref-last; do
        if [ "$order" = "ref-first" ]; then
            $1/plink2 $2 $3 --pfile tmp_data --export bgen-$ver ref-first --out tmp_bgen
        else
            $1/plink2 $2 $3 --pfile tmp_data --export bgen-$ver --out tmp_bgen
        fi
        $1/plink2 $2 $3 --bgen tmp_bgen.bgen $order --sample tmp_bgen.sample --make-pgen --out tmp_rt
        $1/plink2 $2 $3 --pfile tmp_rt --export vcf vcf-dosage=HDS --out tmp_got
        diff <(summarize tmp_want.vcf) <(summarize tmp_got.vcf)
    done
done
