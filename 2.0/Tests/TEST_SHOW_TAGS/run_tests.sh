#!/bin/bash

# --show-tags, checked against PLINK 1.9.
#
# The fixture plants 20 haplotype blocks of 12 variants each, so there is real
# tagging structure inside blocks and none across them.

set -exo pipefail

awk -f make_vcf.awk > tmp_h.vcf
$1/plink2 $2 $3 --vcf tmp_h.vcf --double-id --make-bed --out tmp_data

# A spread of target variants, one per few blocks.
awk '!/^#/ && (NR % 37 == 0) {print $3}' tmp_h.vcf > tmp_targets.txt
test "$(wc -l < tmp_targets.txt)" -gt 3

# 1. Filename mode: the .tags list has to match exactly.
for opt in "" "--tag-kb 50" "--tag-kb 250" "--tag-r2 0.3" "--tag-r2 0.9" "--tag-kb 50 --tag-r2 0.5"
do
    plink --bfile tmp_data --show-tags tmp_targets.txt $opt --out plink19
    $1/plink2 $2 $3 --bfile tmp_data --show-tags tmp_targets.txt $opt --out plink2
    diff -q plink19.tags plink2.tags
done

# 2. 'all' mode: one row per variant, with its tag list.
for opt in "" "--tag-kb 50" "--tag-r2 0.3"
do
    plink --bfile tmp_data --show-tags all $opt --out plink19_all
    $1/plink2 $2 $3 --bfile tmp_data --show-tags all $opt --out plink2_all
    awk -f compare_list.awk plink19_all.tags.list plink2_all.tags.list
done

# 3. --list-all in filename mode writes both files, and the per-variant one is
#    restricted to the named variants.
plink --bfile tmp_data --show-tags tmp_targets.txt --list-all --out plink19_la
$1/plink2 $2 $3 --bfile tmp_data --show-tags tmp_targets.txt --list-all --out plink2_la
diff -q plink19_la.tags plink2_la.tags
awk -f compare_list.awk plink19_la.tags.list plink2_la.tags.list
test "$(grep -vc '^#' plink2_la.tags.list)" -eq "$(wc -l < tmp_targets.txt)"

# 4. The result must not depend on --threads.
for t in 1 3 8
do
    $1/plink2 $2 $3 --bfile tmp_data --show-tags all --threads $t --out plink2_t$t
    diff -q plink2_all.tags.list plink2_t$t.tags.list
done

# 5. 'zs' is the same report, compressed.
$1/plink2 $2 $3 --bfile tmp_data --show-tags all zs --out plink2_zs
$1/plink2 $2 $3 --zst-decompress plink2_zs.tags.list.zst > plink2_zs.tags.list
diff -q plink2_all.tags.list plink2_zs.tags.list

# 6. The report is keyed on variant IDs, so duplicates are refused rather than
#    silently reported against whichever copy was found first.
awk 'BEGIN{OFS="\t"} {if (NR == 2) {$2 = "dup"} else if (NR == 3) {$2 = "dup"}; print}' tmp_data.bim > tmp_dup.bim
cp tmp_data.bed tmp_dup.bed
cp tmp_data.fam tmp_dup.fam
if $1/plink2 $2 $3 --bfile tmp_dup --show-tags all --out plink2_dup > /dev/null 2>&1; then
    echo "duplicate variant IDs accepted"
    exit 1
fi

# 7. Haploid chromosomes are refused for now rather than dropped without
#    saying so.
awk 'BEGIN{OFS="\t"} {if (NR <= 5) {$1 = "X"}; print}' tmp_data.bim > tmp_hap.bim
cp tmp_data.bed tmp_hap.bed
cp tmp_data.fam tmp_hap.fam
if $1/plink2 $2 $3 --bfile tmp_hap --show-tags all --out plink2_hap > /dev/null 2>&1; then
    echo "haploid chromosome accepted"
    exit 1
fi
# ...and excluding them makes the same command work.
$1/plink2 $2 $3 --bfile tmp_hap --autosome --show-tags all --out plink2_auto
test "$(grep -vc '^#' plink2_auto.tags.list)" -gt 0

# 8. The header is in plink2's usual CHROM/POS/ID order.
head -n 1 plink2_all.tags.list | grep -qx '#CHROM	POS	ID	NTAG	LEFT	RIGHT	KBSPAN	TAGS'
