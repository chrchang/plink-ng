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

# 4. --tag-mode2: a two-column target file, and a two-column .tags file
#    covering every retained variant.
awk '{print $1 "\t" ((NR % 3 == 0)? 1 : ((NR % 3 == 1)? 0 : 2))}' tmp_targets.txt > tmp_mode2_targets.txt
awk '!/^#/ {print $3 "\t0"}' tmp_h.vcf | head -n 40 >> tmp_mode2_targets.txt
plink --bfile tmp_data --show-tags tmp_mode2_targets.txt --tag-mode2 --out plink19_m2
$1/plink2 $2 $3 --bfile tmp_data --show-tags tmp_mode2_targets.txt --tag-mode2 --out plink2_m2
diff -q plink19_m2.tags plink2_m2.tags
test "$(wc -l < plink2_m2.tags)" -eq "$(grep -vc '^#' tmp_data.bim)"

# 4b. --tag-mode2 has no meaning with "--show-tags all".
if $1/plink2 $2 $3 --bfile tmp_data --show-tags all --tag-mode2 --out plink2_m2bad 2> tmp_m2_err.txt; then
    echo "--tag-mode2 unexpectedly accepted with \"--show-tags all\""
    exit 1
fi

# 5. The result must not depend on --threads.
for t in 1 3 8
do
    $1/plink2 $2 $3 --bfile tmp_data --show-tags all --threads $t --out plink2_t$t
    diff -q plink2_all.tags.list plink2_t$t.tags.list
done

# 6. 'zs' is the same report, compressed.
$1/plink2 $2 $3 --bfile tmp_data --show-tags all zs --out plink2_zs
$1/plink2 $2 $3 --zst-decompress plink2_zs.tags.list.zst > plink2_zs.tags.list
diff -q plink2_all.tags.list plink2_zs.tags.list
