#!/bin/bash

# Covers both --pmerge/--pmerge-list paths: the concatenation case, where the
# filesets' (chromosome, bp) ranges don't overlap, and the general case, where
# more than one fileset has something to say about the same variant.

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

# 6. Non-concatenating merges.  Every case below reassembles tmp_all from
#    overlapping slices of itself, so the result must reproduce tmp_all.
#    (PLINK 1.9 --bmerge agrees on these, but using tmp_all as the reference
#    makes the check exact and independent of 1.9's allele-order conventions.)
grep -v '^##' tmp_all.pvar > tmp_all_body.pvar

# 6a. Merging a fileset with itself is idempotent.
$1/plink2 $2 $3 --pfile part1 --pmerge part1 --out merged_self
$1/plink2 $2 $3 --pfile part1 --make-bed --out part1b
$1/plink2 $2 $3 --pfile merged_self --make-bed --out merged_selfb
diff -q part1b.bed merged_selfb.bed
diff -q part1b.bim merged_selfb.bim
diff -q part1b.fam merged_selfb.fam

sample_ct=`wc -l < tmp_data.fam`

# 6b. Split by sample into disjoint halves, every variant in both.
head -n $((sample_ct / 2)) tmp_data.fam | cut -f1,2 > tmp_sA.txt
tail -n $((sample_ct - sample_ct / 2)) tmp_data.fam | cut -f1,2 > tmp_sB.txt
$1/plink2 $2 $3 --pfile tmp_all --keep tmp_sA.txt --make-pgen --out sA
$1/plink2 $2 $3 --pfile tmp_all --keep tmp_sB.txt --make-pgen --out sB
$1/plink2 $2 $3 --pfile sA --pmerge sB --out merged_s
grep -v '^##' merged_s.pvar > merged_s_body.pvar
diff -q tmp_all_body.pvar merged_s_body.pvar
$1/plink2 $2 $3 --pfile merged_s --pgen-diff tmp_all --out diff_s
test `grep -cv '^#' diff_s.pdiff` -eq 0

# 6c. Split by sample with an overlapping block, so both filesets have
#     something to say about the same (variant, sample) cells.
head -n $((sample_ct * 3 / 5)) tmp_data.fam | cut -f1,2 > tmp_sC.txt
tail -n $((sample_ct * 3 / 5)) tmp_data.fam | cut -f1,2 > tmp_sD.txt
$1/plink2 $2 $3 --pfile tmp_all --keep tmp_sC.txt --make-pgen --out sC
$1/plink2 $2 $3 --pfile tmp_all --keep tmp_sD.txt --make-pgen --out sD
$1/plink2 $2 $3 --pfile sC --pmerge sD --out merged_so
grep -v '^##' merged_so.pvar > merged_so_body.pvar
diff -q tmp_all_body.pvar merged_so_body.pvar
$1/plink2 $2 $3 --pfile merged_so --pgen-diff tmp_all --out diff_so
test `grep -cv '^#' diff_so.pdiff` -eq 0

# 6d. Split by variant with an overlapping block, same samples throughout.
#     This is the classic PLINK 1.x --bmerge job.
grep -v '^#' tmp_all.pvar | cut -f3 > tmp_vids.txt
variant_ct=`wc -l < tmp_vids.txt`
head -n $((variant_ct * 2 / 3)) tmp_vids.txt > tmp_vA.txt
tail -n $((variant_ct * 2 / 3)) tmp_vids.txt > tmp_vB.txt
$1/plink2 $2 $3 --pfile tmp_all --extract tmp_vA.txt --make-pgen --out vA
$1/plink2 $2 $3 --pfile tmp_all --extract tmp_vB.txt --make-pgen --out vB
$1/plink2 $2 $3 --pfile vA --pmerge vB --out merged_v
grep -v '^##' merged_v.pvar > merged_v_body.pvar
diff -q tmp_all_body.pvar merged_v_body.pvar
$1/plink2 $2 $3 --pfile merged_v --pgen-diff tmp_all --out diff_v
test `grep -cv '^#' diff_v.pdiff` -eq 0

# 6e. Three-way non-concatenating merge, split both ways at once.
$1/plink2 $2 $3 --pfile tmp_all --keep tmp_sA.txt --extract tmp_vA.txt --make-pgen --out xA
$1/plink2 $2 $3 --pfile tmp_all --keep tmp_sB.txt --extract tmp_vB.txt --make-pgen --out xB
echo xA > tmp_xlist.txt
echo xB >> tmp_xlist.txt
echo tmp_all >> tmp_xlist.txt
$1/plink2 $2 $3 --pmerge-list tmp_xlist.txt --out merged_x
grep -v '^##' merged_x.pvar > merged_x_body.pvar
diff -q tmp_all_body.pvar merged_x_body.pvar
$1/plink2 $2 $3 --pfile merged_x --pgen-diff tmp_all --out diff_x
test `grep -cv '^#' diff_x.pdiff` -eq 0

# 6f. Dosages and phase, split by sample with an overlapping block.  Compared
#     through exported VCFs, since --pgen-diff ignores phase.
$1/plink2 $2 $3 --vcf dosage_phase.vcf dosage=DS --make-pgen --out dp_all
$1/plink2 $2 $3 --pfile dp_all --export vcf vcf-dosage=DS --out dp_allv
grep -v '^##' dp_allv.vcf > dp_all_body.vcf
printf 's1\ns2\ns3\ns4\n' > tmp_dpA.txt
printf 's3\ns4\ns5\ns6\n' > tmp_dpB.txt
$1/plink2 $2 $3 --pfile dp_all --keep tmp_dpA.txt --make-pgen --out dpA
$1/plink2 $2 $3 --pfile dp_all --keep tmp_dpB.txt --make-pgen --out dpB
$1/plink2 $2 $3 --pfile dpA --pmerge dpB --out merged_dp
$1/plink2 $2 $3 --pfile merged_dp --export vcf vcf-dosage=DS --out merged_dpv
grep -v '^##' merged_dpv.vcf > merged_dp_body.vcf
diff -q dp_all_body.vcf merged_dp_body.vcf

# 6g. Same data, split by variant instead, so a merged variant's records sit
#     at different variant indices in their respective .pgens.
printf 'v1\nv2\nv3\n' > tmp_dpvA.txt
printf 'v2\nv3\nv4\n' > tmp_dpvB.txt
$1/plink2 $2 $3 --pfile dp_all --extract tmp_dpvA.txt --make-pgen --out dpvA
$1/plink2 $2 $3 --pfile dp_all --extract tmp_dpvB.txt --make-pgen --out dpvB
$1/plink2 $2 $3 --pfile dpvA --pmerge dpvB --out merged_dpv2
$1/plink2 $2 $3 --pfile merged_dpv2 --export vcf vcf-dosage=DS --out merged_dpv2x
grep -v '^##' merged_dpv2x.vcf > merged_dpv2_body.vcf
diff -q dp_all_body.vcf merged_dpv2_body.vcf

# 6h. Staggered variant indices where only one fileset has dosages, and the
#     shared variant sits at a different index in each.  Regression guard: the
#     phase/dosage presence scan used to read the first record's variant index
#     out of every other record's .pgen, which crashes here.
$1/plink2 $2 $3 --vcf mixed_a.vcf dosage=DS --make-pgen --out mixed_a
$1/plink2 $2 $3 --vcf mixed_b.vcf dosage=DS --make-pgen --out mixed_b
$1/plink2 $2 $3 --pfile mixed_a --pmerge mixed_b --out merged_mixed
$1/plink2 $2 $3 --pfile merged_mixed --export vcf vcf-dosage=DS --out merged_mixedv
grep -v '^##' merged_mixedv.vcf > merged_mixed_body.vcf
diff -q expected_mixed.vcfbody merged_mixed_body.vcf

echo "TEST_PMERGE passed."
