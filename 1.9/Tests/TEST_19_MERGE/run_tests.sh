#!/bin/bash

# PLINK 1.9 --bmerge/--merge/--merge-list.  Two overlapping pieces are cut out
# of one dataset, so the merged genotypes are known exactly; plink2
# --pmerge-list is also compared.  (Ported from the old 1.9/tests_old suite,
# which compared against PLINK 1.07.)
#
# Usage: ./run_tests.sh {plink2 build dir} {up to 2 plink2 args}
# Requires plink (1.9) to be in the system PATH.

set -uo pipefail

plink2="$1/plink2 ${2:-} ${3:-}"

failed=0
check() {
    if [[ "$2" -eq 0 ]]; then
        echo "ok: $1"
    else
        echo "FAILED: $1"
        failed=1
    fi
}
# Prints a .bim's variants with their alleles in sorted order.  (A merge
# recomputes which allele is A1.)
bim_alleles() {
    awk '{ print $1, $2, $3, $4, ($5 < $6)? $5 " " $6 : $6 " " $5 }' "$1"
}
# Prints the genotype matrix as counts of the original data's A1 alleles,
# sorted by IID, so that the comparison doesn't depend on allele or sample
# order.  $1: output prefix, remaining args: input flags.
geno_matrix() {
    local out=$1
    shift
    plink "$@" --recode A --recode-allele tmp_a1.txt --out "$out" > /dev/null 2>&1 || return 1
    awk '{ $1 = ""; $3 = $4 = $5 = $6 = ""; print }' "$out.raw" | tr -s ' ' | sed 's/^ //' | sort
}

# 387 samples x 300 variants.  Piece A: the first 201 samples and variants
# 0-199.  Piece B: the other 186 samples and variants 100-299, as a text
# fileset.
plink --dummy 387 300 0.03 acgt --seed 7 --make-bed --out tmp_full > /dev/null
awk 'BEGIN { OFS = "\t" } { $4 = 1000 * NR; print }' tmp_full.bim > tmp_full.bim.new
mv tmp_full.bim.new tmp_full.bim
awk '{ print $2, $5 }' tmp_full.bim > tmp_a1.txt
awk 'NR <= 201 { print $1, $2 }' tmp_full.fam > tmp_a_samples.txt
awk 'NR > 201 { print $1, $2 }' tmp_full.fam > tmp_b_samples.txt
plink --bfile tmp_full --keep tmp_a_samples.txt --from snp0 --to snp199 --make-bed --out tmp_a > /dev/null
plink --bfile tmp_full --keep tmp_b_samples.txt --from snp100 --to snp299 --recode --out tmp_b > /dev/null
plink --bfile tmp_full --keep tmp_b_samples.txt --from snp100 --to snp299 --make-bed --out tmp_b > /dev/null

# Expected merge: the full matrix, with A's samples missing variants 200-299
# and B's samples missing variants 0-99.  (Matrix columns: IID, then snp0 up.)
geno_matrix tmp_full_m --bfile tmp_full | awk '
  NR == FNR { in_a[$2] = 1; next }
  $1 == "IID" { print; next }
  {
    for (i = 2; i <= NF; ++i) {
      v = i - 2
      if (($1 in in_a)? (v >= 200) : (v < 100)) $i = "NA"
    }
    print
  }' tmp_a_samples.txt - > tmp_expected.txt

plink --bfile tmp_a --bmerge tmp_b --make-bed --out tmp_bm > /dev/null
geno_matrix tmp_bm_m --bfile tmp_bm | cmp -s - tmp_expected.txt && cmp -s <(bim_alleles tmp_bm.bim) <(bim_alleles tmp_full.bim)
check "--bmerge of overlapping pieces" $?

plink --bfile tmp_a --merge tmp_b --make-bed --out tmp_m > /dev/null
geno_matrix tmp_m_m --bfile tmp_m | cmp -s - tmp_expected.txt && cmp -s <(bim_alleles tmp_m.bim) <(bim_alleles tmp_bm.bim)
check "--merge of a text fileset" $?

printf 'tmp_a.bed tmp_a.bim tmp_a.fam\ntmp_b.ped tmp_b.map\n' > tmp_list.txt
plink --merge-list tmp_list.txt --make-bed --out tmp_ml > /dev/null
geno_matrix tmp_ml_m --bfile tmp_ml | cmp -s - tmp_expected.txt && cmp -s <(bim_alleles tmp_ml.bim) <(bim_alleles tmp_bm.bim)
check "--merge-list of a binary and a text fileset" $?

plink --bfile tmp_a --bmerge tmp_b --make-just-bim --out tmp_jb > /dev/null
plink --bfile tmp_a --bmerge tmp_b --make-just-fam --out tmp_jf > /dev/null
cmp -s tmp_jb.bim tmp_bm.bim && cmp -s tmp_jf.fam tmp_bm.fam
check "--make-just-bim/--make-just-fam after --bmerge" $?

printf 'tmp_a\ntmp_b\n' > tmp_list2.txt
$plink2 --pmerge-list tmp_list2.txt bfile --make-bed --out tmp_pm > /dev/null
geno_matrix tmp_pm_m --bfile tmp_pm | cmp -s - tmp_expected.txt && cmp -s <(bim_alleles tmp_pm.bim) <(bim_alleles tmp_bm.bim)
check "--bmerge matches plink2 --pmerge-list" $?

# --merge-equal-pos: rs1 (merge.map), rs2 (merge2.map) and rs3 (merge3.map)
# share a position, and must become one variant named rs1.
plink --merge-list merge_list.txt --merge-equal-pos --recode --out tmp_mep > /dev/null
{
    printf '1\trs4\t0\t900\n1\trs1\t0\t1000\n1\trs0\t0\t1100\n'
} | cmp -s - tmp_mep.map && cmp -s - tmp_mep.ped <<'EOF'
1 1 0 0 1 -9 T G C A T GT
2 2 0 0 2 -9 T G A A T GT
3 3 0 0 1 -9 G G A A GT GT
4 4 0 0 2 -9 G G A A GT GT
5 5 0 0 1 -9 T T A A T T
6 6 0 0 1 -9 0 0 C A 0 0
7 7 0 0 2 -9 0 0 A A 0 0
8 8 0 0 1 -9 0 0 C C 0 0
9 9 0 0 2 -9 0 0 C C 0 0
10 10 0 0 1 -9 0 0 C C 0 0
11 11 0 0 1 -9 0 0 C A 0 0
12 12 0 0 2 -9 0 0 A A 0 0
13 13 0 0 1 -9 0 0 A A 0 0
14 14 0 0 2 -9 0 0 C C 0 0
15 15 0 0 1 -9 0 0 C C 0 0
EOF
check "--merge-list --merge-equal-pos combines same-position variants" $?

exit $failed
