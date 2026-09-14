#!/bin/bash

# --make-perm-pheno.  The values cannot be compared against PLINK 1.9, which
# draws from a different random stream, so the invariants are checked instead.

set -exo pipefail

plink --simulate simulate.txt --make-bed --out tmp_data > /dev/null

# Case/control, with some phenotypes missing.
awk 'BEGIN { print "#FID\tIID\tCC" }
{ printf "%s\t%s\t%s\n", $1, $2, (NR % 10 == 0)? "NA" : ((NR % 2) + 1) }' tmp_data.fam > tmp_cc.txt
$1/plink2 $2 $3 --bfile tmp_data --no-psam-pheno --pheno tmp_cc.txt --make-perm-pheno 8 --seed 1 --out plink2_cc
python3 check_pphe.py tmp_cc.txt plink2_cc.pphe

# Quantitative.
awk 'BEGIN { print "#FID\tIID\tQT" } { printf "%s\t%s\t%.4f\n", $1, $2, NR * 0.37 }' tmp_data.fam > tmp_qt.txt
$1/plink2 $2 $3 --bfile tmp_data --no-psam-pheno --pheno tmp_qt.txt --make-perm-pheno 5 --seed 1 --out plink2_qt
python3 check_pphe.py tmp_qt.txt plink2_qt.pphe

# --seed makes it reproducible, and a different seed gives something else.
$1/plink2 $2 $3 --bfile tmp_data --no-psam-pheno --pheno tmp_cc.txt --make-perm-pheno 8 --seed 1 --out plink2_again
diff -q plink2_cc.pphe plink2_again.pphe
$1/plink2 $2 $3 --bfile tmp_data --no-psam-pheno --pheno tmp_cc.txt --make-perm-pheno 8 --seed 2 --out plink2_seed2
if diff -q plink2_cc.pphe plink2_seed2.pphe > /dev/null; then
    echo "two different seeds produced identical permutations"
    exit 1
fi

# Zstd output round-trips.
$1/plink2 $2 $3 --bfile tmp_data --no-psam-pheno --pheno tmp_cc.txt --make-perm-pheno 8 zs --seed 1 --out plink2_zs
$1/plink2 $2 $3 --zst-decompress plink2_zs.pphe.zst > plink2_zs.pphe
diff -q plink2_zs.pphe plink2_cc.pphe

# Naming the phenotype is required when more than one is loaded.
awk 'BEGIN { print "#FID\tIID\tA\tB" } { printf "%s\t%s\t1\t%.3f\n", $1, $2, NR * 0.11 }' tmp_data.fam > tmp_two.txt
if $1/plink2 $2 $3 --bfile tmp_data --no-psam-pheno --pheno tmp_two.txt --make-perm-pheno 3 --out plink2_bad 2> tmp_err.txt; then
    echo "expected --make-perm-pheno to require a phenotype name"
    exit 1
fi
grep -q "name the one" tmp_err.txt
$1/plink2 $2 $3 --bfile tmp_data --no-psam-pheno --pheno tmp_two.txt --make-perm-pheno 3 B --seed 1 --out plink2_named
# check_pphe wants a two-column phenotype file; build one holding just B.
awk '{ printf "%s\t%s\t%s\n", $1, $2, $4 }' tmp_two.txt > tmp_b.txt
python3 check_pphe.py tmp_b.txt plink2_named.pphe
