#!/bin/bash

# --pmerge-list .psam merging under the default nm-match modes.
#
# 1. A maternal-ID conflict locked the *paternal* ID instead of the maternal
#    one, so later files could still overwrite the mother, but not fill in a
#    missing father.
# 2. The missing-category name was put in the wrong category hash table
#    slot, corrupting categorical phenotype lookups.
# 3. A quantitative value of 0 seen before a conflict resurrected the
#    conflicted (locked, missing) value as 0.

set -exo pipefail

plink2="$1/plink2 $2 $3"

for k in 1 2 3; do
    cat > tmp_f$k.vcf <<VCF
##fileformat=VCFv4.2
##contig=<ID=1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	C	D	E
1	$((100 * k))	v$k	A	G	.	.	.	GT	0/1	0/0	1/1
VCF
done
# C: mother M1 vs M2 conflict, father only in file 3; Q is 0, 1.5, 2.5.
# D: CAT only nonmissing in file 2.  E: CAT blue/NONE/blue.
printf '#IID\tPAT\tMAT\tSEX\tQ\tCAT\nC\t0\tM1\t1\t0\tred\nD\t0\t0\t2\tNA\tNONE\nE\t0\t0\t1\t3\tblue\n' > tmp_f1.psam
printf '#IID\tPAT\tMAT\tSEX\tQ\tCAT\nC\t0\tM2\t1\t1.5\tred\nD\t0\t0\t2\t4\tgreen\nE\t0\t0\t1\t3\tNONE\n' > tmp_f2.psam
printf '#IID\tPAT\tMAT\tSEX\tQ\tCAT\nC\tP1\tM3\t1\t2.5\tred\nD\t0\t0\t2\tNA\tNONE\nE\t0\t0\t1\t3\tblue\n' > tmp_f3.psam
for k in 1 2 3; do
    $plink2 --vcf tmp_f$k.vcf --psam tmp_f$k.psam --make-pgen --out tmp_f$k > /dev/null
done
printf 'tmp_f1\ntmp_f2\ntmp_f3\n' > tmp_list.txt
$plink2 --pmerge-list tmp_list.txt --out tmp_merged > /dev/null
diff tmp_merged.psam - <<'PSAM'
#IID	PAT	MAT	SEX	Q	CAT
C	P1	0	1	NA	red
D	0	0	2	4	green
E	0	0	1	3	blue
PSAM
