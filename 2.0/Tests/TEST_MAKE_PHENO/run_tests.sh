#!/bin/bash

set -exo pipefail

# 200 samples, 100 variants; only the sample IDs matter here.
plink --dummy 200 100 --make-bed --out tmp_raw

awk 'BEGIN {
  for (i = 0; i < 200; ++i) {
    printf "fam%d per%d 0 0 %d -9\n", i, i, 1 + (i % 2) > "tmp_data.fam"
  }
  for (j = 0; j < 100; ++j) {
    printf "1 var%d 0 %d A B\n", j, (j + 1) * 1000 > "tmp_data.bim"
  }
}'
cp tmp_raw.bed tmp_data.bed

# Roughly half the samples are listed, with four distinct values, and three
# IDs which are not in the dataset at all.
awk 'BEGIN {
  s = 8675309
  split("CASE CTRL 1 2", vals, " ")
  for (i = 0; i < 200; ++i) {
    s = (1103515245 * s + 12345) % 2147483648
    if ((int(s / 512) % 2) == 0) {
      continue
    }
    s = (1103515245 * s + 12345) % 2147483648
    printf "fam%d per%d %s\n", i, i, vals[1 + (int(s / 512) % 4)] > "list.txt"
  }
  for (k = 0; k < 3; ++k) {
    printf "ghost%d gh%d CASE\n", k, k > "list.txt"
  }
}'

check() {
    plink --bfile tmp_data --make-pheno list.txt "$1" --make-just-fam --out plink_mp
    $2/plink2 $3 $4 --bfile tmp_data --make-pheno list.txt "$1" --make-just-psam --out plink2_mp
    awk '{print $2 "\t" ($6 == "-9" ? "NA" : $6)}' plink_mp.fam > flat19.txt
    awk 'NR > 1 {print $2 "\t" $NF}' plink2_mp.psam > flat20.txt
    diff -q flat19.txt flat20.txt
    # Guard against a vacuous comparison: at least one case and one control.
    test $(grep -c '	2$' flat20.txt) -ge 1
    test $(grep -c '	1$' flat20.txt) -ge 1
}

check '*' $1 $2 $3
check CASE $1 $2 $3
check CTRL $1 $2 $3
check 2 $1 $2 $3

# A value matching nothing leaves every listed sample a control.
plink --bfile tmp_data --make-pheno list.txt nosuchvalue --make-just-fam --out plink_mp
$1/plink2 $2 $3 --bfile tmp_data --make-pheno list.txt nosuchvalue --make-just-psam --out plink2_mp
awk '{print $2 "\t" ($6 == "-9" ? "NA" : $6)}' plink_mp.fam > flat19.txt
awk 'NR > 1 {print $2 "\t" $NF}' plink2_mp.psam > flat20.txt
diff -q flat19.txt flat20.txt
test $(grep -c '	2$' flat20.txt) -eq 0

fails() {
    "$@" && false || true
}

# Wrong parameter count, and a repeated sample ID.
fails $1/plink2 $2 $3 --bfile tmp_data --make-pheno list.txt --make-just-psam --out plink2_bad
head -n 1 list.txt > dup.txt
head -n 1 list.txt >> dup.txt
fails $1/plink2 $2 $3 --bfile tmp_data --make-pheno dup.txt CASE --make-just-psam --out plink2_bad

# The name is taken when the .psam already has a MAKEPHENO column.
$1/plink2 $2 $3 --bfile tmp_data --make-pheno list.txt CASE --make-just-psam --out plink2_named
fails $1/plink2 $2 $3 --pfile plink2_named --make-pheno list.txt CASE --make-just-psam --out plink2_bad
