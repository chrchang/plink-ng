#!/bin/bash

# --export A-transpose header with the 'spaces' modifier.
#
# The FID/IID split searched for the output delimiter instead of the tab
# that separates them internally.

set -exo pipefail

plink2="$1/plink2 $2 $3"

cat > tmp_data.vcf <<'VCF'
##fileformat=VCFv4.2
##contig=<ID=1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	F1_S1	F2_S2
1	100	v1	A	G	.	.	.	GT	0/1	0/0
VCF
$plink2 --vcf tmp_data.vcf --id-delim _ --make-pgen --out tmp_data > /dev/null
$plink2 --pfile tmp_data --export A-transpose --out tmp_tab > /dev/null
$plink2 --pfile tmp_data --export A-transpose spaces --out tmp_space > /dev/null
[[ "$(head -n 1 tmp_tab.traw | cut -f 7-)" == "$(printf 'F1_S1\tF2_S2')" ]] || exit 1
[[ "$(head -n 1 tmp_space.traw | cut -d ' ' -f 7-)" == "F1_S1 F2_S2" ]] || exit 1
