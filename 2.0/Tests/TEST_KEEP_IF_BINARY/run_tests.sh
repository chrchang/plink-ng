#!/bin/bash

# --keep-if/--remove-if on binary phenotypes and on covariates.
#
# '!=' on a binary phenotype used to behave like '=='.  Covariate lookup used
# to stop after the first <number of phenotypes> covariates.

set -exo pipefail

plink2="$1/plink2 $2 $3"

# S1, S2 cases; S3, S4 controls; S5 missing phenotype.
printf '#IID\tSEX\tPHENO1\nS1\t1\t2\nS2\t2\t2\nS3\t1\t1\nS4\t2\t1\nS5\t1\tNA\n' > tmp_data.psam
printf '#IID\tC1\tC2\tC3\nS1\t0\t0\t1\nS2\t0\t0\t9\nS3\t0\t0\t7\nS4\t0\t0\t2\nS5\t0\t0\t8\n' > tmp_covar.txt
cat > tmp_data.vcf <<'VCF'
##fileformat=VCFv4.2
##contig=<ID=1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5
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

# Missing-phenotype samples satisfy "!=" (as with the other phenotype types),
# but not "==".
$plink2 --pfile tmp_data --keep-if "PHENO1 != case" --make-just-psam --out tmp_neq_case > /dev/null
check_iids tmp_neq_case.psam "S3 S4 S5 "
$plink2 --pfile tmp_data --keep-if "PHENO1 != control" --make-just-psam --out tmp_neq_ctrl > /dev/null
check_iids tmp_neq_ctrl.psam "S1 S2 S5 "
$plink2 --pfile tmp_data --keep-if "PHENO1 == case" --make-just-psam --out tmp_eq_case > /dev/null
check_iids tmp_eq_case.psam "S1 S2 "
$plink2 --pfile tmp_data --remove-if "PHENO1 != case" --make-just-psam --out tmp_rm_neq_case > /dev/null
check_iids tmp_rm_neq_case.psam "S1 S2 "

# A covariate past the first (phenotype count) covariates.
$plink2 --pfile tmp_data --covar tmp_covar.txt --keep-if "C3 > 5" --make-just-psam --out tmp_c3 > /dev/null
check_iids tmp_c3.psam "S2 S3 S5 "
