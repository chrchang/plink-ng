#!/bin/bash

# Provisional-REF flags after variant filtering.
#
# 1. Fixed-width .pgen output (format=2) wrote the unfiltered flag array.
# 2. multiallelics=- looked up the flag by filtered index instead of raw
#    index.

set -exo pipefail

plink2="$1/plink2 $2 $3"

cat > tmp_data.vcf <<'VCF'
##fileformat=VCFv4.2
##contig=<ID=1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
1	100	v1	A	G	.	.	.	GT	0/0	0/1	1/1
1	200	v2	A	G	.	.	.	GT	0/1	0/1	0/0
1	300	v3	A	G	.	.	.	GT	1/1	0/0	0/1
1	400	v4	A	G,T	.	.	.	GT	0/2	1/2	0/0
1	500	v5	A	G	.	.	.	GT	0/1	0/0	0/1
1	600	v6	A	G,C	.	.	.	GT	0/0	2/2	0/1
1	700	v7	A	G	.	.	.	GT	0/0	0/1	1/1
VCF
$plink2 --vcf tmp_data.vcf --make-pgen --out tmp_raw > /dev/null
# v2, v5, v6 have provisional REF alleles.
awk 'BEGIN {OFS = "\t"}
     /^##/ {print; next}
     /^#CHROM/ {print "##INFO=<ID=PR,Number=0,Type=Flag,Description=\"Provisional reference allele, may not be based on real reference genome\">"; print $0, "INFO"; next}
     {print $0, (($3 == "v2" || $3 == "v5" || $3 == "v6")? "PR" : ".")}' tmp_raw.pvar > tmp_data.pvar
cp tmp_raw.pgen tmp_data.pgen
cp tmp_raw.psam tmp_data.psam
printf 'v1\nv3\n' > tmp_exclude.txt

# Reads back <prefix>, with the INFO column dropped from the .pvar so the
# flags have to come from the .pgen header.
pr_flags() {
    cut -f 1-5 "$1".pvar | grep -v '^##INFO' > "$1"_noinfo.pvar
    $plink2 --pgen "$1".pgen --psam "$1".psam --pvar "$1"_noinfo.pvar --make-just-pvar cols=+info --out "$1"_chk > /dev/null
    awk '!/^#/ {printf "%s ", $NF}' "$1"_chk.pvar
}

# 1. format=2 (biallelic only).
$plink2 --pfile tmp_data --exclude tmp_exclude.txt --snps-only just-acgt --max-alleles 2 --make-pgen format=2 --out tmp_fw > /dev/null
[[ "$(pr_flags tmp_fw)" == "PR PR . " ]]

# 2. multiallelics=-: v4 and v6 are each split in two.
$plink2 --pfile tmp_data --exclude tmp_exclude.txt --make-pgen multiallelics=- --out tmp_split > /dev/null
[[ "$(pr_flags tmp_split)" == "PR . . PR PR PR . " ]]
