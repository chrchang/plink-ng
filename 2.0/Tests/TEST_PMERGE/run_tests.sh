#!/bin/bash

# --pmerge/--pmerge-list currently implements the concatenation case: filesets
# whose (chromosome, bp) ranges don't overlap.  This exercises that path, since
# it had no coverage.
#
# Non-concatenating merge is still under development; when it lands, the
# natural additions here are a split by sample (same variants, disjoint
# samples) compared against PLINK 1.9 --bmerge, and the --merge-mode variants.

set -exo pipefail

plink --simulate simulate.txt --simulate-missing 0.02 --out tmp_data

# Three chromosomes, so the parts have disjoint coordinate ranges.
head -n 99 tmp_data.bim > tmp_data1.bim
cat tmp_data.bim | sed -n '100,198p' | sed 's/^1/2/' > tmp_data2.bim
cat tmp_data.bim | tail -n 99 | sed 's/^1/3/' > tmp_data3.bim
cat tmp_data1.bim tmp_data2.bim tmp_data3.bim > tmp_data.bim

$1/plink2 $2 $3 --bfile tmp_data --make-pgen --out tmp_all

for c in 1 2 3
do
    $1/plink2 $2 $3 --pfile tmp_all --chr $c --make-pgen --out part$c
done

# 1. --pmerge, two filesets.
$1/plink2 $2 $3 --pfile part1 --pmerge part2 --out merged12
$1/plink2 $2 $3 --pfile tmp_all --chr 1,2 --make-pgen --out expected12
$1/plink2 $2 $3 --pfile merged12 --make-bed --out merged12b
$1/plink2 $2 $3 --pfile expected12 --make-bed --out expected12b
diff -q expected12b.bed merged12b.bed
diff -q expected12b.bim merged12b.bim
diff -q expected12b.fam merged12b.fam

# 2. --pmerge-list, all three, given as a list file.
echo part1 > tmp_list.txt
echo part2 >> tmp_list.txt
echo part3 >> tmp_list.txt
$1/plink2 $2 $3 --pmerge-list tmp_list.txt --out merged_all
$1/plink2 $2 $3 --pfile merged_all --make-bed --out merged_allb
$1/plink2 $2 $3 --pfile tmp_all --make-bed --out expected_allb
diff -q expected_allb.bed merged_allb.bed
diff -q expected_allb.bim merged_allb.bim
diff -q expected_allb.fam merged_allb.fam

# 3. Same, with the initially-loaded fileset supplying one of the parts.
echo part2 > tmp_list23.txt
echo part3 >> tmp_list23.txt
$1/plink2 $2 $3 --pfile part1 --pmerge-list tmp_list23.txt --out merged_all2
$1/plink2 $2 $3 --pfile merged_all2 --make-bed --out merged_all2b
diff -q expected_allb.bed merged_all2b.bed
diff -q expected_allb.bim merged_all2b.bim
diff -q expected_allb.fam merged_all2b.fam

# 4. Reversed input order: the result must not depend on it.
echo part3 > tmp_list_rev.txt
echo part2 >> tmp_list_rev.txt
echo part1 >> tmp_list_rev.txt
$1/plink2 $2 $3 --pmerge-list tmp_list_rev.txt --out merged_rev
$1/plink2 $2 $3 --pfile merged_rev --make-bed --out merged_revb
diff -q expected_allb.bed merged_revb.bed
diff -q expected_allb.bim merged_revb.bim

# 5. bfile mode, i.e. concatenating PLINK 1 filesets.
for c in 1 2 3
do
    $1/plink2 $2 $3 --pfile part$c --make-bed --out bpart$c
done
echo bpart1 > tmp_blist.txt
echo bpart2 >> tmp_blist.txt
echo bpart3 >> tmp_blist.txt
$1/plink2 $2 $3 --pmerge-list tmp_blist.txt bfile --out merged_b
$1/plink2 $2 $3 --pfile merged_b --make-bed --out merged_bb
diff -q expected_allb.bed merged_bb.bed
diff -q expected_allb.bim merged_bb.bim

# 6. Merging a fileset with itself is a non-concatenating job, and must be
#    rejected cleanly rather than producing wrong output.
if $1/plink2 $2 $3 --pfile part1 --pmerge part1 --out merged_self 2> tmp_self_err.txt; then
    echo "self-merge unexpectedly succeeded"
    exit 1
fi
grep -q "under development" tmp_self_err.txt

# 7. --merge-mode nm-match over a group of records.  Same-position same-ID
#    duplicates inside one fileset are the only way to reach the genotype merge
#    while the job stays a concatenation, so the fixtures pair them with a
#    fileset on chromosome 2.  The _a and _b fixtures hold the same data with
#    the two records of each variant swapped, and both are checked against the
#    same expected output, since "nonmissing values must match" is symmetric:
#      s1 is missing in every record, so it must stay missing;
#      s2 is missing in one record and called in the other, so the call wins;
#      s3 has two different calls, which is the one genuine conflict.
#    nmdos covers the same three cases on the dosage track.
for o in a b
do
    $1/plink2 $2 $3 --vcf nmdup_$o.vcf --make-pgen --out tmp_nmdup_$o
    $1/plink2 $2 $3 --vcf nmdos_$o.vcf dosage=DS --make-pgen --out tmp_nmdos_$o
done
$1/plink2 $2 $3 --vcf nmdup_other.vcf --make-pgen --out tmp_nmother

for o in a b
do
    $1/plink2 $2 $3 --pfile tmp_nmdup_$o --pmerge tmp_nmother --merge-mode nm-match --out tmp_nmdupm_$o
    $1/plink2 $2 $3 --pfile tmp_nmdupm_$o --export vcf --out tmp_nmdupx_$o
    grep -v '^#' tmp_nmdupx_$o.vcf > tmp_nmdup_body_$o.vcf
    diff -q nmdup_expected.vcfbody tmp_nmdup_body_$o.vcf

    $1/plink2 $2 $3 --pfile tmp_nmdos_$o --pmerge tmp_nmother --merge-mode nm-match --out tmp_nmdosm_$o
    $1/plink2 $2 $3 --pfile tmp_nmdosm_$o --export vcf vcf-dosage=DS --out tmp_nmdosx_$o
    grep -v '^#' tmp_nmdosx_$o.vcf > tmp_nmdos_body_$o.vcf
    diff -q nmdos_expected.vcfbody tmp_nmdos_body_$o.vcf
done

echo "TEST_PMERGE passed."
