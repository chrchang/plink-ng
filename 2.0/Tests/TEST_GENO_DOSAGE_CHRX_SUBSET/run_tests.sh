#!/bin/bash

# "--geno dosage" on chrX with a sample subset.
#
# On chrX, the dosage-missing count was taken over all samples rather than
# the included ones.  (Autosomes and chrY derive it from dosage sums, which
# were already restricted to included samples.)

set -exo pipefail

plink2="$1/plink2 $2 $3"

# S5 and S6 are removed.  On v1, S1 has a dosage but no hardcall, and S5/S6
# are fully missing; among S1-S4 nothing is dosage-missing.
cat > tmp_data.vcf <<'VCF'
##fileformat=VCFv4.2
##contig=<ID=X,length=100000000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DS,Number=A,Type=Float,Description="Dosage">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5	S6
X	10000000	v1	A	G	.	.	.	GT:DS	./.:0.5	0/1:1	0/0:0	1/1:2	./.:.	./.:.
X	10001000	v2	C	T	.	.	.	GT:DS	0/1:1	0/1:1	0/0:0	1/1:2	0/0:0	0/0:0
VCF
printf '#IID\tSEX\nS1\t2\nS2\t2\nS3\t2\nS4\t2\nS5\t2\nS6\t2\n' > tmp_data.psam
$plink2 --vcf tmp_data.vcf dosage=DS --psam tmp_data.psam --make-pgen --out tmp_data > /dev/null
printf 'S5\nS6\n' > tmp_remove.txt

$plink2 --pfile tmp_data --remove tmp_remove.txt --geno 0.2 dosage --make-just-pvar --out tmp_geno > /dev/null
[[ $(grep -vc '^#' tmp_geno.pvar) -eq 2 ]] || exit 1
# Hardcall-based --geno still removes v1 (S1 has no hardcall).
$plink2 --pfile tmp_data --remove tmp_remove.txt --geno 0.2 --make-just-pvar --out tmp_geno_hc > /dev/null
[[ $(grep -vc '^#' tmp_geno_hc.pvar) -eq 1 ]] || exit 1
