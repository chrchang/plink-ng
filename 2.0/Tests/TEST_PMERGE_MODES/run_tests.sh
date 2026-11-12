#!/bin/bash

# --pmerge conflict-resolution and output-shaping flags with no other test
# coverage: --merge-parents-mode, --merge-sex-mode, --merge-pheno-mode,
# --merge-pheno-sort, --pheno-inner-join, --merge-xheader-mode,
# --merge-max-alleles, and --delete-pmerge-result.

set -exo pipefail

plink2="$1/plink2 $2 $3"

vcf_header() {
    printf '##fileformat=VCFv4.2\n##contig=<ID=1,length=10000>\n##contig=<ID=2,length=10000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
}

# A and B share samples C and D, and variants v1 and v3 (with different ALT
# alleles at v3).  C has conflicting PAT, SEX, and Q values; D has missing
# MAT/SEX/Q in A that B fills in.  The phenotype columns other than Q are
# unique to one fileset, and named to distinguish ASCII from natural sort.
{
    vcf_header
    printf '##note=progA\n##shared=same\n##onlyA=1\n'
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tC\tD\tE\n'
    printf '1\t100\tv1\tA\tG\t.\t.\t.\tGT\t0/1\t./.\t1/1\n'
    printf '1\t200\tv2\tC\tT\t.\t.\t.\tGT\t0/0\t0/1\t0/1\n'
    printf '1\t300\tv3\tA\tC\t.\t.\t.\tGT\t0/1\t0/0\t0/0\n'
} > tmp_a.vcf
{
    vcf_header
    printf '##note=progB\n##shared=same\n##onlyB=1\n'
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tC\tD\tF\n'
    printf '1\t100\tv1\tA\tG\t.\t.\t.\tGT\t0/1\t1/1\t0/0\n'
    printf '1\t300\tv3\tA\tG\t.\t.\t.\tGT\t./.\t0/1\t1/1\n'
    printf '1\t400\tv4\tG\tT\t.\t.\t.\tGT\t0/1\t0/0\t1/1\n'
} > tmp_b.vcf
printf '#IID\tPAT\tMAT\tSEX\tQ\tz9\nC\tP1\t0\t1\t1.5\t7\nD\t0\t0\tNA\tNA\t8\nE\t0\t0\t2\t3\t9\n' > tmp_a.psam
printf '#IID\tPAT\tMAT\tSEX\tz10\tAb\tQ\nC\tP2\t0\t2\tx\t1\t2.5\nD\t0\tM1\t2\ty\t2\t4\nF\t0\t0\t1\tz\t3\t5\n' > tmp_b.psam
for x in a b; do
    $plink2 --vcf tmp_$x.vcf --psam tmp_$x.psam --make-pgen --out tmp_$x > /dev/null
done

# $1: output prefix; remaining args: extra flags.
pmerge_ab() {
    local out=$1
    shift
    $plink2 --pfile tmp_a --pmerge tmp_b "$@" --out $out > /dev/null
}

# $1: .psam; $2: expected IID/PAT/MAT/SEX/Q rows for C and D, '|'-joined.
check_cd() {
    test "$(awk -F '\t' 'NR == 1 {for (i = 1; i <= NF; i++) col[$i] = i; next} $1 == "C" || $1 == "D" {print $1 "\t" $col["PAT"] "\t" $col["MAT"] "\t" $col["SEX"] "\t" $col["Q"]}' $1 | tr '\n' '|')" = "$(printf "$2")|"
}

# Default nm-match everywhere: conflicts become missing, and a missing value
# never conflicts.
pmerge_ab tmp_m0
diff tmp_m0.psam - <<'PSAM'
#IID	PAT	MAT	SEX	Q	z9	z10	Ab
C	0	0	NA	NA	7	x	1
D	0	M1	2	4	8	y	2
E	0	0	2	3	9	NONE	NA
F	0	0	1	5	NA	z	3
PSAM
grep -v '^#' tmp_m0.pvar | cut -f 3-5 | diff - <(printf 'v1\tA\tG\nv2\tC\tT\nv3\tA\tC,G\nv4\tG\tT\n')
$plink2 --pfile tmp_m0 --export vcf --out tmp_m0 > /dev/null
# v3/D: A's 0/0 (A/A) conflicts with B's A/G.
grep -v '^#' tmp_m0.vcf | cut -f 3,10- | diff - <(printf 'v1\t0/1\t1/1\t1/1\t0/0\nv2\t0/0\t0/1\t0/1\t./.\nv3\t0/1\t./.\t0/0\t2/2\nv4\t0/1\t0/0\t./.\t1/1\n')

# The three modes are independent of each other.
nm_match_cd='C\t0\t0\tNA\tNA|D\t0\tM1\t2\t4'
nm_first_cd='C\tP1\t0\t1\t1.5|D\t0\tM1\t2\t4'
first_cd='C\tP1\t0\t1\t1.5|D\t0\t0\tNA\tNA'
for mode in nm-match 1; do
    pmerge_ab tmp_mnm --merge-parents-mode $mode --merge-sex-mode $mode --merge-pheno-mode $mode
    check_cd tmp_mnm.psam "$nm_match_cd"
done
for mode in nm-first 2; do
    pmerge_ab tmp_mnf --merge-parents-mode $mode --merge-sex-mode $mode --merge-pheno-mode $mode
    check_cd tmp_mnf.psam "$nm_first_cd"
done
for mode in first 4; do
    pmerge_ab tmp_mf --merge-parents-mode $mode --merge-sex-mode $mode --merge-pheno-mode $mode
    check_cd tmp_mf.psam "$first_cd"
done
pmerge_ab tmp_mpar --merge-parents-mode nm-first
check_cd tmp_mpar.psam 'C\tP1\t0\tNA\tNA|D\t0\tM1\t2\t4'
pmerge_ab tmp_mpar2 --merge-parents-mode first
check_cd tmp_mpar2.psam 'C\tP1\t0\tNA\tNA|D\t0\t0\t2\t4'
pmerge_ab tmp_msex --merge-sex-mode first
check_cd tmp_msex.psam 'C\t0\t0\t1\tNA|D\t0\tM1\tNA\t4'
pmerge_ab tmp_mph --merge-pheno-mode nm-first
check_cd tmp_mph.psam 'C\t0\t0\tNA\t1.5|D\t0\tM1\t2\t4'
# Samples present in only one fileset are unaffected by 'first'.
pmerge_ab tmp_mph2 --merge-pheno-mode first
grep -q "$(printf '^F\t0\t0\t1\t5\tNA\tz\t3$')" tmp_mph2.psam

# Phenotype column order and inner join.
for sort in none 0; do
    pmerge_ab tmp_s --merge-pheno-sort $sort
    test "$(head -n 1 tmp_s.psam)" = "$(printf '#IID\tPAT\tMAT\tSEX\tQ\tz9\tz10\tAb')"
done
for sort in ascii a; do
    pmerge_ab tmp_s --merge-pheno-sort $sort
    test "$(head -n 1 tmp_s.psam)" = "$(printf '#IID\tPAT\tMAT\tSEX\tAb\tQ\tz10\tz9')"
done
for sort in natural n; do
    pmerge_ab tmp_s --merge-pheno-sort $sort
    test "$(head -n 1 tmp_s.psam)" = "$(printf '#IID\tPAT\tMAT\tSEX\tAb\tQ\tz9\tz10')"
done
# Values must move with their columns.
grep -q "$(printf '^C\t0\t0\tNA\t1\tNA\t7\tx$')" tmp_s.psam
pmerge_ab tmp_pij --pheno-inner-join
diff tmp_pij.psam - <<'PSAM'
#IID	PAT	MAT	SEX	Q
C	0	0	NA	NA
D	0	M1	2	4
E	0	0	2	3
F	0	0	1	5
PSAM

# .pvar header lines.  (##source is dropped on VCF import, so use ##note.)
xheader_lines() {
    grep '^##' $1 | grep -v '^##contig' | tr '\n' ' '
}
pmerge_ab tmp_xf
test "$(xheader_lines tmp_xf.pvar)" = "##note=progA ##shared=same ##onlyA=1 ##onlyB=1 "
pmerge_ab tmp_xf2 --merge-xheader-mode first
test "$(xheader_lines tmp_xf2.pvar)" = "##note=progA ##shared=same ##onlyA=1 ##onlyB=1 "
pmerge_ab tmp_xm --merge-xheader-mode match
test "$(xheader_lines tmp_xm.pvar)" = "##shared=same ##onlyA=1 ##onlyB=1 "
# 'erase' requires --merge-info-mode erase.
if pmerge_ab tmp_xe --merge-xheader-mode erase 2> /dev/null; then
    exit 1
fi
pmerge_ab tmp_xe --merge-xheader-mode erase --merge-info-mode erase
# This removes ##contig lines too.
test "$(head -n 1 tmp_xe.pvar | cut -f 1)" = "#CHROM"

# --merge-max-alleles: currently only supported by concatenating merges, so
# merge A with a chr2-only fileset containing a triallelic variant.
{
    vcf_header
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tC\tD\tF\n'
    printf '2\t100\tw1\tA\tG,T\t.\t.\t.\tGT\t0/1\t1/2\t0/0\n'
    printf '2\t200\tw2\tC\tT\t.\t.\t.\tGT\t0/1\t0/0\t1/1\n'
} > tmp_c.vcf
$plink2 --vcf tmp_c.vcf --make-pgen --out tmp_c > /dev/null
$plink2 --pfile tmp_a --pmerge tmp_c --out tmp_cat > /dev/null
test "$(grep -v '^#' tmp_cat.pvar | cut -f 3 | tr '\n' ' ')" = "v1 v2 v3 w1 w2 "
$plink2 --pfile tmp_a --pmerge tmp_c --merge-max-alleles 3 --out tmp_cmma3 > /dev/null
diff -q tmp_cat.pvar tmp_cmma3.pvar
$plink2 --pfile tmp_a --pmerge tmp_c --merge-max-alleles 2 --out tmp_cmma2 > /dev/null
test "$(grep -v '^#' tmp_cmma2.pvar | cut -f 3 | tr '\n' ' ')" = "v1 v2 v3 w2 "
# Genotypes stay aligned after the drop.
$plink2 --pfile tmp_cmma2 --export vcf --out tmp_cmma2 > /dev/null
test "$(grep -v '^#' tmp_cmma2.vcf | cut -f 3,10- | tail -n 1)" = "$(printf 'w2\t0/1\t0/0\t./.\t1/1')"

# --delete-pmerge-result: later commands still run on the merged data, then
# the merged fileset is removed.
rm -f tmp_del.*
pmerge_ab tmp_del --delete-pmerge-result --freq
test ! -e tmp_del.pgen
test ! -e tmp_del.pvar
test ! -e tmp_del.psam
test "$(grep -v '^#' tmp_del.afreq | cut -f 2 | tr '\n' ' ')" = "v1 v2 v3 v4 "
