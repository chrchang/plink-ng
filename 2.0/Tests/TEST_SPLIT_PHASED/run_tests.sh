#!/bin/bash

# --make-pgen multiallelics=- on phased multiallelic variants.  Two defects
# used to corrupt the output:
# - when a split variant had both phased and unphased hets, its phase record
#   lost the explicit-phasepresent flag and was then overwritten by the next
#   variant's genotypes (garbage calls downstream, occasionally a crash);
# - phase bits for ALT2+ were assigned in list order rather than sample order,
#   so 'x|y' and '0|y' hets in the same 32-sample block swapped phase.

set -exo pipefail

cat > split.vcf << 'EOF'
##fileformat=VCFv4.3
##contig=<ID=1,length=1000000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	a	b	c	d	e
1	100	v1	A	C,G	.	.	.	GT	0|1	0/1	1|0	0/2	0|0
1	200	v2	A	C	.	.	.	GT	0|1	1|0	0/0	1/1	0|0
1	300	v3	C	A,G,T	.	.	.	GT	3|1	1|1	0/3	3/3	2/3
1	400	v4	C	A,G	.	.	.	GT	2|1	0|2	0|0	1|2	0|1
EOF

# Expected split genotypes, as (unphased allele count) or (phased a|b).
cat > expected.txt << 'EOF'
100 C 0|1 1 1|0 0 0
100 G 0 0 0 1 0
200 C 0|1 1|0 0 2 0
300 A 0|1 2 0 0 0
300 G 0 0 0 0 1
300 T 1|0 0 1 2 1
400 A 0|1 0 0 1|0 0|1
400 G 1|0 0|1 0 0|1 0
EOF

$1/plink2 $2 $3 --vcf split.vcf --make-pgen --out tmp_data
$1/plink2 $2 $3 --pfile tmp_data --make-pgen multiallelics=- --out tmp_split
$1/plink2 $2 $3 --pfile tmp_split --export vcf --out tmp_split
grep -v '^#' tmp_split.vcf | awk -F'\t' '{
    line = $2 " " $5
    for (i = 10; i <= NF; i++) {
        gt = $i
        if ((gt == "0|1") || (gt == "1|0")) {
            line = line " " gt
        } else {
            line = line " " (substr(gt, 1, 1) + substr(gt, 3, 1))
        }
    }
    print line
}' > tmp_got.txt
diff expected.txt tmp_got.txt
