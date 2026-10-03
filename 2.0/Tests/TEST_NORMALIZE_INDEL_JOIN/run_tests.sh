#!/bin/bash

# --normalize indel-join + --make-pgen multiallelics=+ should merge
# overlapping deletions like 'bcftools norm -m +both'.  The expected output
# below matches bcftools 1.x on the same input, except that bcftools drops
# s3's phase at position 44 (its other source records are unphased hom-ref
# calls), while plink2 keeps it.

set -exo pipefail

printf '>1\nGATTACAGATTACACCGGTTAACCGGTTAAGCGCATATCGCGTAGCTAGCTA\n' > ref.fa

cat > in.vcf << 'EOF2'
##fileformat=VCFv4.2
##contig=<ID=1,length=52>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1	s2	s3	s4
1	30	d1	AG	A	.	.	.	GT	0|1	0/0	1/1	0/0
1	30	d2	AGC	A	.	.	.	GT	1|0	0/1	0/0	0/0
1	30	i1	A	AT	.	.	.	GT	0/0	0/0	0/0	0/1
1	30	s1	A	C	.	.	.	GT	0/0	0/1	0/0	0/0
1	44	e1	AG	A	.	.	.	GT	0/1	0/0	0/0	0/1
1	44	e2	AGCT	A	.	.	.	GT	0/1	0/0	0/0	0/0
1	44	e3	A	T	.	.	.	GT	0/0	1/1	0/0	0/0
1	44	e4	A	AC	.	.	.	GT	0/0	0/0	0|1	0/0
EOF2

cat > expected.txt << 'EOF2'
1 30 d1;d2;i1 AGC AC,A,ATGC 2|1 0/2 1/1 0/3
1 30 s1 A C 0/0 0/1 0/0 0/0
1 44 e1;e2;e4 AGCT ACT,A,ACGCT 1/2 0/0 0|3 0/1
1 44 e3 A T 0/0 1/1 0/0 0/0
EOF2

$1/plink2 $2 $3 --vcf in.vcf --fa ref.fa --normalize indel-join --make-pgen multiallelics=+ varid-join --out tmp_join
grep -q -- '--normalize indel-join: REF allele extended for 4 variants.' tmp_join.log
$1/plink2 $2 $3 --pfile tmp_join --export vcf --out tmp_join
# Homozygous calls are written as 'x/y', since their phase doesn't matter.
grep -v '^#' tmp_join.vcf | cut -f 1-5,10- | sed 's/\([0-9]\)|\1/\1\/\1/g' | tr '\t' ' ' > tmp_got.txt
diff expected.txt tmp_got.txt

# Inconsistent REF alleles at one position.
sed 's/^1	30	d1	AG	/1	30	d1	AT	/' in.vcf > bad.vcf
if $1/plink2 $2 $3 --vcf bad.vcf --fa ref.fa --normalize indel-join --make-pgen --out tmp_bad; then
    exit 1
fi
grep -q "Variants 'd1' and 'd2' have the same position" tmp_bad.log

# Overlapping deletions at different positions are joined when
# left-normalization moves them to the same position (o1/o2, the example from
# the #584 review); bcftools norm -f ref2.fa -m +any gives the same alleles and
# genotypes, with the ALTs in the other order.  n1/n2 overlap but stay at
# different positions, so neither tool joins them; neither does m1/m2, which
# end at the same position (the second example from the review).  This also checks that
# --normalize is accepted with multiallelics=+snps and +any.
printf '>1\nGGTAAAAGGCTACGCAGGCAGATGAAATCC\n' > ref2.fa

cat > in2.vcf << 'EOF2'
##fileformat=VCFv4.2
##contig=<ID=1,length=30>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1	s2
1	3	o1	TAAAA	T	.	.	.	GT	0/0	0/1
1	5	o2	AAA	A	.	.	.	GT	1/1	0/1
1	11	n1	TACGC	T	.	.	.	GT	0/1	0/0
1	13	n2	CGC	C	.	.	.	GT	0/0	0/1
1	19	m1	CAGATGAAAT	C	.	.	.	GT	0/0	0/1
1	24	m2	GAAAT	G	.	.	.	GT	1/1	0/0
EOF2

cat > expected2.txt << 'EOF2'
1 3 o2;o1 TAAAA TAA,T 1/1 1/2
1 11 n1 TACGC T 0/1 0/0
1 12 n2 ACG A 0/0 0/1
1 19 m1 CAGATGAAAT C 0/0 0/1
1 24 m2 GAAAT G 1/1 0/0
EOF2

for mode in both any; do
    $1/plink2 $2 $3 --vcf in2.vcf --fa ref2.fa --normalize indel-join --make-pgen multiallelics=+$mode varid-join --out tmp_ovl_$mode
    $1/plink2 $2 $3 --pfile tmp_ovl_$mode --export vcf --out tmp_ovl_$mode
    grep -v '^#' tmp_ovl_$mode.vcf | cut -f 1-5,10- | tr '\t' ' ' > tmp_got.txt
    diff expected2.txt tmp_got.txt
done
$1/plink2 $2 $3 --vcf in2.vcf --fa ref2.fa --normalize indel-join --make-pgen multiallelics=+snps --out tmp_ovl_snps

# Flag checks.
if $1/plink2 $2 $3 --vcf in.vcf --fa ref.fa --normalize left indel-join --make-pgen --out tmp_flag2; then
    exit 1
fi
