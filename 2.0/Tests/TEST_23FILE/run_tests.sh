#!/bin/bash

# --23file import and "--export 23" output, checked against PLINK 1.9.
#
# The multiallelic part is plink2-only: PLINK 1.9 has no way to represent
# those variants, so the expected output is spelled out here.

set -exo pipefail

cat > in23.txt << 'EOF'
# rsid	chromosome	position	genotype
rs1	1	100	AA
rs2	1	200	AG
rs3	1	300	GG
rs4	1	400	--
rs5	2	100	CT
rs6	2	200	TT
rs7	3	100	DD
rs8	3	200	DI
rs9	3	300	II
rsx1	X	3000000	G
rsx2	X	3000001	A
rsy1	Y	3000000	T
rsmt1	MT	100	C
EOF

# 1. Import: the genotypes and the allele assignment must match PLINK 1.9.
#    plink2 writes chromosome names rather than codes, tab-delimits the .fam,
#    and spells a missing allele '.' instead of '0', so the .bim/.fam
#    comparisons normalize those.
plink --23file in23.txt FAM1 SAMP1 M --make-bed --out plink19
$1/plink2 $2 $3 --23file in23.txt FAM1 SAMP1 M --output-chr 26 --make-bed --out plink2
cmp plink19.bed plink2.bed
awk 'BEGIN{OFS="\t"} {if ($5 == ".") $5 = 0; if ($6 == ".") $6 = 0; print}' plink2.bim > plink2_norm.bim
diff -q plink19.bim plink2_norm.bim
tr -s ' \t' ' ' < plink19.fam > plink19_norm.fam
tr -s ' \t' ' ' < plink2.fam > plink2_norm.fam
diff -q plink19_norm.fam plink2_norm.fam

# 2. Sex inference: 'I' has to reach the same conclusion PLINK 1.9 does (male
#    here, since there are single-allele chrX calls).
$1/plink2 $2 $3 --23file in23.txt FAM1 SAMP1 I --output-chr 26 --make-bed --out plink2_infer
cmp plink19.bed plink2_infer.bed
tr -s ' \t' ' ' < plink2_infer.fam > plink2_infer_norm.fam
diff -q plink19_norm.fam plink2_infer_norm.fam

# 3. "--export 23" against PLINK 1.9's --recode 23, on the same fileset.  The
#    first line carries a timestamp, so only the body is compared.
plink --bfile plink19 --recode 23 --out plink19_out
$1/plink2 $2 $3 --bfile plink19 --export 23 --out plink2_out
grep -v '^#' plink19_out.txt > plink19_body.txt
grep -v '^#' plink2_out.txt > plink2_body.txt
diff -q plink19_body.txt plink2_body.txt

# 4. Round trip: re-importing the exported file reproduces the genotypes.
$1/plink2 $2 $3 --23file plink2_out.txt FAM1 SAMP1 M --output-chr 26 --make-bed --out plink2_rt
cmp plink19.bed plink2_rt.bed

# 5. Multiallelic export.  chrX backs off to two-character genotypes because
#    the ALT1/ALT2 call at rsx2 is heterozygous; chrY keeps one character for
#    the homozygous call and two for the heterozygous one.
cat > multi.vcf << 'EOF'
##fileformat=VCFv4.2
##contig=<ID=1>
##contig=<ID=X>
##contig=<ID=Y>
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA1
1	1	rm1	A	G,T	.	.	.	GT	1/2
1	2	rm2	A	G,T	.	.	.	GT	2/2
1	3	rm3	A	G,T	.	.	.	GT	0/2
X	3000000	rmx1	A	G	.	.	.	GT	1/1
X	3000001	rmx2	A	G,T	.	.	.	GT	1/2
Y	3000000	rmy1	A	G,T	.	.	.	GT	2/2
Y	3000001	rmy2	A	G,T	.	.	.	GT	1/2
EOF
printf 'NA1\tNA1\t1\n' > multi_sex.txt
$1/plink2 $2 $3 --vcf multi.vcf --double-id --update-sex multi_sex.txt --make-pgen --out plink2_multi
$1/plink2 $2 $3 --pfile plink2_multi --export 23 --out plink2_multi_out
grep -v '^#' plink2_multi_out.txt > plink2_multi_body.txt
cat > expected_multi.txt << 'EOF'
rm1	1	1	GT
rm2	1	2	TT
rm3	1	3	AT
rmx1	X	3000000	GG
rmx2	X	3000001	GT
rmy1	Y	3000000	T
rmy2	Y	3000001	GT
EOF
diff -q expected_multi.txt plink2_multi_body.txt

# 6. A long chromosome name has to survive the row writer.  plink2 allows an
#    ID component to be as long as any other, and --export 23 emits whatever
#    chromosome code the dataset carries, so the output line is not bounded
#    by anything small.
long_chr=$(printf 'contig_%0.sX' $(seq 1 200))
{
    printf '##fileformat=VCFv4.2\n##contig=<ID=%s>\n' "$long_chr"
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n'
    printf '%s\t1\tlv1\tA\tG\t.\t.\t.\tGT\t0/1\n' "$long_chr"
} > tmp_long.vcf
$1/plink2 $2 $3 --vcf tmp_long.vcf --allow-extra-chr --make-pgen --out plink2_long
$1/plink2 $2 $3 --pfile plink2_long --allow-extra-chr --export 23 --out plink2_long_out
# rsid, the full chromosome name, position and genotype, tab-separated.
grep -v '^#' plink2_long_out.txt > plink2_long_body.txt
test "$(wc -l < plink2_long_body.txt)" -eq 1
awk -v chr="$long_chr" -F '\t' '{
    if ($1 != "lv1") {print "wrong rsid: " $1; exit 1}
    if ($2 != chr) {print "chromosome truncated or corrupted"; exit 1}
    if ($3 != 1) {print "wrong position: " $3; exit 1}
    if ($4 != "AG") {print "wrong genotype: " $4; exit 1}
}' plink2_long_body.txt
