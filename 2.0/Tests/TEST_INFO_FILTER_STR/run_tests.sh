#!/bin/bash

# --extract-if-info/--exclude-if-info string comparisons.

set -exo pipefail

plink2="$1/plink2 $2 $3"

cat > tmp_data.vcf <<'VCF'
##fileformat=VCFv4.2
##contig=<ID=1,length=1000>
##INFO=<ID=VT,Number=1,Type=String,Description="Variant type">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
1	100	v1	A	G	.	.	VT=SNP	GT	0/1	0/0
1	200	v2	A	AG	.	.	VT=INDEL	GT	0/1	0/0
1	300	v3	C	T	.	.	VT=SNPX	GT	0/1	0/0
1	400	v4	C	G	.	.	VT=SN	GT	0/1	0/0
1	500	v5	G	T	.	.	.	GT	0/1	0/0
1	600	v6	T	C	.	.	VT=SNP	GT	0/1	0/0
VCF
$plink2 --vcf tmp_data.vcf --make-pgen --out tmp_data > /dev/null

check_ids() {
    got=$(awk '!/^#/ {printf "%s ", $3}' "$1")
    if [[ "$got" != "$2" ]]; then
        echo "$1: got '$got', expected '$2'"
        exit 1
    fi
}

$plink2 --pfile tmp_data --extract-if-info "VT == SNP" --make-just-pvar --out tmp_eq > /dev/null
check_ids tmp_eq.pvar "v1 v6 "
# '!=' is true when the key is missing.
$plink2 --pfile tmp_data --extract-if-info "VT != SNP" --make-just-pvar --out tmp_neq > /dev/null
check_ids tmp_neq.pvar "v2 v3 v4 v5 "
$plink2 --pfile tmp_data --exclude-if-info "VT == SNPX" --make-just-pvar --out tmp_ex > /dev/null
check_ids tmp_ex.pvar "v1 v2 v4 v5 v6 "
