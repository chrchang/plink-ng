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
cat tmp_data.bim | head -n 99 > tmp_data1.bim
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

echo "TEST_PMERGE passed."
