#!/bin/bash

# --normalize when an allele is entirely consumed by right-trimming and the
# variant is then left-shifted.
#
# The padding base was taken from the pre-shift reference position instead
# of the shifted one.

set -exo pipefail

plink2="$1/plink2 $2 $3"

# Reference: T T C A C A C A G G ...  REF=ACAG/ALT=AG at position 6 is a
# deletion of "AC" inside the CA repeat; left-normalized, it's TCA -> T at
# position 2.
printf '>1\nTTCACACAGGTTTTTTTTTT\n' > tmp_ref.fa
cat > tmp_data.vcf <<'VCF'
##fileformat=VCFv4.2
##contig=<ID=1,length=20>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
1	6	v1	ACAG	AG	.	.	.	GT	0/1	0/0
VCF
$plink2 --vcf tmp_data.vcf --make-pgen --out tmp_data > /dev/null
$plink2 --pfile tmp_data --fa tmp_ref.fa --normalize --make-just-pvar --out tmp_norm > /dev/null
[[ "$(grep -v '^#' tmp_norm.pvar | cut -f 1-5)" == "$(printf '1\t2\tv1\tTCA\tT')" ]] || exit 1
