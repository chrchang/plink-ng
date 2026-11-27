#!/bin/bash

# --set-invalid-haploid-missing on multiallelic chrY variants with female
# samples.
#
# Female altx/altx calls (x >= 2) were set to missing without dropping their
# patch entries, producing invalid .pgen records.

set -exo pipefail

plink2="$1/plink2 $2 $3"

cat > tmp_data.vcf <<'VCF'
##fileformat=VCFv4.2
##contig=<ID=Y,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	m1	f1	m2	f2
Y	100	v1	A	C,G	.	.	.	GT	2/2	2/2	1/2	1/1
Y	200	v2	A	C,G,T	.	.	.	GT	3/3	0/2	2/2	3/3
VCF
printf '#IID\tSEX\nm1\t1\nf1\t2\nm2\t1\nf2\t2\n' > tmp_data.psam
$plink2 --vcf tmp_data.vcf --psam tmp_data.psam --make-pgen --out tmp_data > /dev/null
$plink2 --pfile tmp_data --set-invalid-haploid-missing --make-pgen --out tmp_out > /dev/null
$plink2 --pfile tmp_out --validate > /dev/null
$plink2 --pfile tmp_out --export vcf --out tmp_out > /dev/null
# Female calls are missing; male hets are missing, male homs unchanged.
# (chrY calls are exported as haploid.)
[[ "$(awk '!/^#/ && $3 == "v1" {print $10, $11, $12, $13}' tmp_out.vcf)" == "2 . . ." ]] || exit 1
[[ "$(awk '!/^#/ && $3 == "v2" {print $10, $11, $12, $13}' tmp_out.vcf)" == "3 . 2 ." ]] || exit 1
