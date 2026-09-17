#!/bin/bash

set -exo pipefail

# Deterministic dataset: 200 samples, 400 variants, built with a plain LCG so
# that every awk implementation produces the same files.  Only the IDs matter
# here, so the genotypes are arbitrary.
# Deterministic dataset: 200 samples, 400 variants.  Only the IDs matter here,
# so the genotypes come from --dummy and the .fam/.bim are rewritten with
# predictable IDs.
plink --dummy 200 400 --make-bed --out tmp_raw

awk 'BEGIN {
  for (i = 0; i < 200; ++i) {
    printf "fam%d per%d 0 0 %d %d\n", i, i, 1 + (i % 2), 1 + (i % 2) > "tmp_data.fam"
  }
  for (j = 0; j < 400; ++j) {
    printf "1 var%d 0 %d A B\n", j, (j + 1) * 1000 > "tmp_data.bim"
  }
}'
cp tmp_raw.bed tmp_data.bed

# Attribute files: every third variant/sample is left out entirely, and the
# attribute sets are chosen so that each condition form has a nonempty and
# non-total answer.
awk 'BEGIN {
  s = 999
  for (j = 0; j < 400; ++j) {
    s = (1103515245 * s + 12345) % 2147483648
    if ((int(s / 512) % 3) == 0) {
      continue
    }
    line = "var" j
    for (k = 0; k < 5; ++k) {
      s = (1103515245 * s + 12345) % 2147483648
      if ((int(s / 512) % 2) == 0) {
        line = line "\t" ("attr" k)
      }
    }
    print line > "vattr.txt"
  }
  s = 31337
  for (i = 0; i < 200; ++i) {
    s = (1103515245 * s + 12345) % 2147483648
    if ((int(s / 512) % 3) == 0) {
      continue
    }
    line = "fam" i " per" i
    for (k = 0; k < 5; ++k) {
      s = (1103515245 * s + 12345) % 2147483648
      if ((int(s / 512) % 2) == 0) {
        line = line " " ("attr" k)
      }
    }
    print line > "sattr.txt"
  }
}'



for cond in "" "attr0" "attr0,attr3" ",-attr1" "attr0,-attr1" "attr2,attr4,-attr0" ",-attr0,-attr1"; do
    plink --bfile tmp_data --attrib vattr.txt $cond --write-snplist --out plink_va
    $1/plink2 $2 $3 --bfile tmp_data --attrib vattr.txt $cond --write-snplist --out plink2_va
    diff -q plink_va.snplist plink2_va.snplist
    test $(wc -l < plink2_va.snplist) -ge 15

    plink --bfile tmp_data --attrib-indiv sattr.txt $cond --make-just-fam --out plink_sa
    $1/plink2 $2 $3 --bfile tmp_data --attrib-indiv sattr.txt $cond --write-samples --out plink2_sa
    awk '{print $1 "\t" $2}' plink_sa.fam > flat19.txt
    grep -v '^#' plink2_sa.id | awk '{print $1 "\t" $2}' > flat20.txt
    diff -q flat19.txt flat20.txt
    test $(wc -l < flat20.txt) -ge 5
done

# An attribute name which appears nowhere keeps nothing.
plink --bfile tmp_data --attrib vattr.txt nosuchattr --write-snplist --out plink_va && false || true
$1/plink2 $2 $3 --bfile tmp_data --attrib vattr.txt nosuchattr --write-snplist --out plink2_va && false || true

fails() {
    "$@" && false || true
}

# Malformed conditions.
fails $1/plink2 $2 $3 --bfile tmp_data --attrib vattr.txt ,--attr0 --write-snplist --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --attrib vattr.txt attr0,attr0 --write-snplist --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --attrib vattr.txt attr0,-attr0 --write-snplist --out plink2_bad

# A repeated ID in the attribute file is an error, as in PLINK 1.x.
head -n 1 vattr.txt > dup.txt
head -n 1 vattr.txt >> dup.txt
fails $1/plink2 $2 $3 --bfile tmp_data --attrib dup.txt --write-snplist --out plink2_bad
