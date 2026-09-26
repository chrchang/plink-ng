#!/bin/bash

# "--export vcf vcf-dosage=DS-force" on unphased multiallelic chrX variants
# used to apply the first sample's sex to every sample in the same 32-sample
# block, so female calls after a male were written as haploid (with halved
# DS), and male calls after a female as diploid.

set -exo pipefail

cat > multi.vcf << 'EOF'
##fileformat=VCFv4.3
##contig=<ID=X,length=100000000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	m1	f1	f2	m2	f3	u1
X	5000000	v1	C	G,T	.	.	.	GT	0	1/1	0/0	2	1/2	0/1
X	5000100	v2	C	G	.	.	.	GT	1	1/1	0/0	0	0/1	0/1
X	5000200	v3	A	C,G,T	.	.	.	GT	3	3/3	2/2	1	./.	1/1
EOF
printf '#IID\tSEX\nm1\t1\nf1\t2\nf2\t2\nm2\t1\nf3\t2\nu1\tNA\n' > multi.psam

$1/plink2 $2 $3 --vcf multi.vcf --psam multi.psam --make-pgen --out tmp_data
$1/plink2 $2 $3 --pfile tmp_data --export vcf --out tmp_gt
for mode in DS-force HDS-force; do
    $1/plink2 $2 $3 --pfile tmp_data --export vcf vcf-dosage=$mode --out tmp_ds
    # GT subfields must match the plain export.
    diff <(grep -v '^#' tmp_gt.vcf | cut -f 10-) <(grep -v '^#' tmp_ds.vcf | cut -f 10- | awk -F'\t' -v OFS='\t' '{for (i = 1; i <= NF; i++) {sub(/:.*/, "", $i)} print}')
    # DS must count the ALT alleles in GT.
    grep -v '^#' tmp_ds.vcf | awk -F'\t' '{
        n_alt = split($5, alts, ",")
        for (i = 10; i <= NF; i++) {
            split($i, sub_fields, ":")
            if (sub_fields[1] ~ /^\./) {
                continue
            }
            n_gt = split(sub_fields[1], gt, "/")
            n_ds = split(sub_fields[2], ds, ",")
            if (n_ds != n_alt) {
                exit 1
            }
            for (a = 1; a <= n_alt; a++) {
                ct = 0
                for (j = 1; j <= n_gt; j++) {
                    if (gt[j] == a) {
                        ct++
                    }
                }
                if (ds[a] != ct) {
                    exit 1
                }
            }
        }
    }'
done
