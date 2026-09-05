#!/bin/bash

# --twolocus, checked against PLINK 1.9.
#
# PLINK 1.9 writes a fixed-width count matrix with marginal totals per group;
# plink2 writes one row per (group, genotype, genotype) cell.  compare.awk
# reads both and matches the cells.

set -exo pipefail

BUILD=$1
EXTRA1=$2
EXTRA2=$3

awk -f make_vcf.awk > tmp_h.vcf
$BUILD/plink2 $EXTRA1 $EXTRA2 --vcf tmp_h.vcf --double-id --make-bed --out tmp_data
# PLINK 1.9 ignores the phenotypes of missing-sex samples, which would drop its
# CASE/CTRL sections, so the fixture assigns sex as well.
awk 'BEGIN{OFS=" "} {print $1, $2, 0, 0, (NR % 2) + 1, ((NR % 3)? 2 : 1)}' tmp_data.fam > tmp_ped.fam
mv tmp_ped.fam tmp_data.fam

compare_pair() {
    plink --bfile tmp_data --twolocus "$1" "$2" --out plink19
    $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --twolocus "$1" "$2" --out plink2
    awk -v id1="$1" -v id2="$2" -f compare.awk plink19.twolocus plink2.twolocus
}

# 1. A pair inside one haplotype block, a pair spanning two, and a pair
#    involving a variant with missing calls.
compare_pair b0v0 b0v3
compare_pair b0v1 b1v7
compare_pair b5v2 b5v9

# 2. Without a case/control phenotype there is only the ALL group.
awk 'BEGIN{OFS=" "} {print $1, $2, 0, 0, (NR % 2) + 1, -9}' tmp_data.fam > tmp_nopheno.fam
cp tmp_data.bed tmp_nopheno.bed
cp tmp_data.bim tmp_nopheno.bim
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_nopheno --twolocus b0v0 b0v3 --out plink2_nopheno
test "$(awk '!/^#/ {print $1}' plink2_nopheno.twolocus | sort -u | tr -d '\n')" = "ALL"

# 3. Every cell count is present, so the marginals PLINK 1.9 prints can be
#    recovered: the ALL group has to sum to the sample count.
test "$(awk '!/^#/ && $1 == "ALL" {s += $6} END {print s}' plink2.twolocus)" -eq "$(wc -l < tmp_data.fam)"

# 4. A duplicate or unknown variant ID is an error rather than a wrong report.
if $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --twolocus b0v0 nosuchvariant --out plink2_bad 2> tmp_bad_err.txt; then
    echo "--twolocus unexpectedly accepted an unknown variant ID"
    exit 1
fi

# 5. Multiallelic variants get one row per unordered allele pair, rather than
#    being rejected or collapsed to REF vs. non-REF.
cat > tmp_multi.vcf << 'EOF'
##fileformat=VCFv4.2
##contig=<ID=1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1	s2	s3	s4	s5	s6
1	1	mv1	A	G,T	.	.	.	GT	0/0	0/1	0/2	1/1	1/2	2/2
1	2	mv2	C	T	.	.	.	GT	0/0	0/1	1/1	0/0	0/1	./.
EOF
$BUILD/plink2 $EXTRA1 $EXTRA2 --vcf tmp_multi.vcf --make-pgen --out tmp_multi
$BUILD/plink2 $EXTRA1 $EXTRA2 --pfile tmp_multi --twolocus mv1 mv2 --out plink2_multi
# 3 alleles give 6 genotypes plus a missing cell; 2 alleles give 3 plus one.
test "$(grep -c '^ALL' plink2_multi.twolocus)" -eq 28
# Every genotype the file contains has to appear, and each sample lands in
# exactly one cell.
for gt in 'A/A' 'A/G' 'G/G' 'A/T' 'G/T' 'T/T'; do
    grep -q "	$gt	" plink2_multi.twolocus
done
test "$(awk '$1 == "ALL" {n += $6} END {print n}' plink2_multi.twolocus)" -eq 6
# The six samples are one per cell, on the diagonal of the pairing above.
test "$(awk '$1 == "ALL" && $6 == 1' plink2_multi.twolocus | wc -l)" -eq 6

# 6. 'zs' is gone: the report is one row per joint genotype cell, so it was
#    never large enough to be worth compressing.
if $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --twolocus zs b5v2 b5v9 --out plink2_zs > /dev/null 2>&1; then
    echo "--twolocus still accepts 'zs'"
    exit 1
fi
