#!/bin/bash

# Deprecated PLINK 1.x flag aliases, --script, and --eigsnp.
#
# --filter-females used to be rewritten to --remove-males (which also keeps
# unknown-sex samples), --remove-founders and --remove-nonfounders were
# rewritten to each other's meaning, --script spliced its own token back into
# the command line, and --eigsnp was shadowed by a second --eigind branch.

set -exo pipefail

plink2="$1/plink2 $2 $3"

# A: male founder, B: female founder, C: female nonfounder, D: unknown-sex
# nonfounder, E: unknown-sex founder.
printf '#IID\tPAT\tMAT\tSEX\nA\t0\t0\t1\nB\t0\t0\t2\nC\tA\tB\t2\nD\tA\tB\tNA\nE\t0\t0\tNA\n' > tmp_data.psam
cat > tmp_data.vcf <<'VCF'
##fileformat=VCFv4.2
##contig=<ID=1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	A	B	C	D	E
1	100	v1	A	G	.	.	.	GT	0/1	0/0	1/1	0/1	0/0
VCF
$plink2 --vcf tmp_data.vcf --psam tmp_data.psam --make-pgen --out tmp_data > /dev/null

check_iids() {
    got=$(awk 'NR > 1 {printf "%s ", $1}' "$1")
    if [[ "$got" != "$2" ]]; then
        echo "$1: got '$got', expected '$2'"
        exit 1
    fi
}

$plink2 --pfile tmp_data --filter-females --make-just-psam --out tmp_ff > /dev/null
check_iids tmp_ff.psam "B C "
$plink2 --pfile tmp_data --remove-founders --make-just-psam --out tmp_rf > /dev/null
check_iids tmp_rf.psam "C D "
$plink2 --pfile tmp_data --remove-nonfounders --make-just-psam --out tmp_rnf > /dev/null
check_iids tmp_rnf.psam "A B E "

# --script must behave like the same arguments given directly.
echo "--pfile tmp_data --keep-females --make-just-psam" > tmp_script.txt
$plink2 --script tmp_script.txt --out tmp_script > /dev/null
check_iids tmp_script.psam "B C "

# --eiggeno/--eigind/--eigsnp round trip.
$plink2 --pfile tmp_data --export eig --out tmp_eig > /dev/null
$plink2 --eiggeno tmp_eig.geno --eigind tmp_eig.ind --eigsnp tmp_eig.snp --make-pgen --out tmp_eig > /dev/null
$plink2 --pfile tmp_data --export A --out tmp_orig > /dev/null
$plink2 --pfile tmp_eig --export A --out tmp_eig > /dev/null
[[ "$(cut -f 2,7 tmp_orig.raw | tail -n +2)" == "$(cut -f 2,7 tmp_eig.raw | tail -n +2)" ]]
