#!/bin/bash

set -exo pipefail

# 120 samples, 300 variants.  Samples 0-39 form block b1 and 40-69 form b2;
# variants 0-49 belong to b1 and 50-79 to b2, and every call in a block is
# written as missing.  The remaining calls are missing 5% of the time.
awk 'BEGIN {
  for (i = 0; i < 120; ++i) {
    printf "fam%d per%d 0 0 %d 1\n", i, i, 1 + (i % 2) > "tmp_data.fam"
  }
  for (j = 0; j < 300; ++j) {
    printf "1 var%d 0 %d A B\n", j, (j + 1) * 1000 > "tmp_data.bim"
  }
  for (i = 0; i < 40; ++i) {
    printf "fam%d per%d b1\n", i, i > "samples.txt"
  }
  for (i = 40; i < 70; ++i) {
    printf "fam%d per%d b2\n", i, i > "samples.txt"
  }
  for (j = 0; j < 50; ++j) {
    printf "var%d b1\n", j > "variants.txt"
  }
  for (j = 50; j < 80; ++j) {
    printf "var%d b2\n", j > "variants.txt"
  }
}'

# The .bed is written by an awk program rather than --dummy, since the
# missingness pattern is the point of the test.
awk 'BEGIN {
  s = 20250910
  printf "%c%c%c", 108, 27, 1 > "tmp_data.bed"
  for (j = 0; j < 300; ++j) {
    for (base = 0; base < 120; base += 4) {
      byte = 0
      for (k = 0; k < 4; ++k) {
        i = base + k
        s = (1103515245 * s + 12345) % 2147483648
        in_block = ((j < 50 && i < 40) || (j >= 50 && j < 80 && i >= 40 && i < 70))
        if (in_block || (int(s / 512) % 100) < 5) {
          g = 1
        } else {
          g = ((int(s / 512) % 3) == 0) ? 0 : (((int(s / 512) % 3) == 1) ? 2 : 3)
        }
        byte += g * (4 ^ k)
      }
      printf "%c", byte > "tmp_data.bed"
    }
  }
}'

check_geno() {
    plink --bfile tmp_data --oblig-missing variants.txt samples.txt --geno $1 --write-snplist --out plink_om
    $2/plink2 $3 $4 --bfile tmp_data --oblig-missing variants.txt samples.txt --geno $1 --write-snplist --out plink2_om
    diff -q plink_om.snplist plink2_om.snplist
}

check_mind() {
    plink --bfile tmp_data --oblig-missing variants.txt samples.txt --mind $1 --make-just-fam --out plink_om
    $2/plink2 $3 $4 --bfile tmp_data --oblig-missing variants.txt samples.txt --mind $1 --write-samples --out plink2_om
    awk '{print $1 "\t" $2}' plink_om.fam > flat19.txt
    grep -v '^#' plink2_om.id | awk '{print $1 "\t" $2}' > flat20.txt
    diff -q flat19.txt flat20.txt
}

check_geno 0.02 $1 $2 $3
check_geno 0.1 $1 $2 $3
check_geno 0.25 $1 $2 $3
check_mind 0.02 $1 $2 $3
check_mind 0.1 $1 $2 $3
check_mind 0.25 $1 $2 $3

# The flag has to matter: without it, the blocks are ordinary missing calls and
# both filters bite hard.
$1/plink2 $2 $3 --bfile tmp_data --geno 0.1 --write-snplist --out plink2_noom
test $(wc -l < plink2_noom.snplist) -lt 260
$1/plink2 $2 $3 --bfile tmp_data --oblig-missing variants.txt samples.txt --geno 0.1 --write-snplist --out plink2_om
test $(wc -l < plink2_om.snplist) -ge 295

$1/plink2 $2 $3 --bfile tmp_data --mind 0.1 --write-samples --out plink2_noom
test $(grep -vc '^#' plink2_noom.id) -lt 80
$1/plink2 $2 $3 --bfile tmp_data --oblig-missing variants.txt samples.txt --mind 0.1 --write-samples --out plink2_om
test $(grep -vc '^#' plink2_om.id) -eq 120

fails() {
    "$@" && false || true
}

# Wrong parameter count, and unreadable files.
fails $1/plink2 $2 $3 --bfile tmp_data --oblig-missing variants.txt --geno 0.1 --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --oblig-missing variants.txt nosuchfile.txt --geno 0.1 --out plink2_bad

# IDs which are absent from the dataset are ignored, and a block ID which no
# sample uses draws a warning rather than an error.
cp variants.txt variants_extra.txt
echo "nosuchvariant b1" >> variants_extra.txt
echo "var0 orphanblock" >> variants_extra.txt
cp samples.txt samples_extra.txt
echo "ghost gh b1" >> samples_extra.txt
plink --bfile tmp_data --oblig-missing variants_extra.txt samples_extra.txt --geno 0.1 --write-snplist --out plink_om
$1/plink2 $2 $3 --bfile tmp_data --oblig-missing variants_extra.txt samples_extra.txt --geno 0.1 --write-snplist --out plink2_om
diff -q plink_om.snplist plink2_om.snplist
