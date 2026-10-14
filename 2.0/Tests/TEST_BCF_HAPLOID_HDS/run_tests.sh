#!/bin/bash

# --bcf dosage=HDS must read a record whose HDS field has only one value per
# sample (every sample haploid, e.g. chrY or chrMT without phase) the same way
# --vcf reads the equivalent text.  Previously, the next sample's HDS value
# was taken as the second haplotype, turning haploid calls into phased hets.

set -exo pipefail

cat > haploid.vcf << 'EOF'
##fileformat=VCFv4.3
##contig=<ID=Y,length=1000000>
##contig=<ID=MT,length=20000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DS,Number=A,Type=Float,Description="Dosage">
##FORMAT=<ID=HDS,Number=.,Type=Float,Description="Haplotype dosage">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	a	b	c	d
Y	100	y1	A	G	.	.	.	GT:DS:HDS	1:1:1	0:0:0	1:1:1	0:0.3:0.3
Y	200	y2	C	T	.	.	.	GT:DS:HDS	0:0:0	1:1:1	1:1:1	1:0.8:0.8
MT	300	m1	A	G	.	.	.	GT:DS:HDS	1:1:1	0:0:0	0:0:0	.:.:.
EOF
printf '#IID\tSEX\na\t1\nb\t1\nc\t1\nd\t1\n' > haploid.psam

$1/plink2 $2 $3 --vcf haploid.vcf dosage=HDS --psam haploid.psam --make-pgen --out tmp_vcf
$1/plink2 $2 $3 --pfile tmp_vcf --export vcf vcf-dosage=DS-force --out tmp_vcf_out
$1/plink2 $2 $3 --pfile tmp_vcf --export bcf vcf-dosage=HDS-force --out tmp_bcf
$1/plink2 $2 $3 --bcf tmp_bcf.bcf dosage=HDS --psam haploid.psam --make-pgen --out tmp_bcf_in
$1/plink2 $2 $3 --pfile tmp_bcf_in --export vcf vcf-dosage=DS-force --out tmp_bcf_out
diff <(grep -v '^##' tmp_vcf_out.vcf) <(grep -v '^##' tmp_bcf_out.vcf)
# The calls must still be haploid.
if grep -v '^#' tmp_bcf_out.vcf | cut -f 10- | grep -q '[|/]'; then
    exit 1
fi
