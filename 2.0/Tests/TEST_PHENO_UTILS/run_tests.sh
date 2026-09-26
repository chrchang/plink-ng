#!/bin/bash

# --tail-pheno, checked against PLINK 1.9.

set -exo pipefail

plink --dummy 200 400 0.03 --out tmp_data > /dev/null

# A quantitative phenotype spanning both --tail-pheno bounds.
awk 'BEGIN{OFS=" "; print "#FID", "IID", "QT"} {print $1, $2, ($NR % 17) - 6 + (NR % 5) * 0.25}' tmp_data.fam > tmp_qt.txt

pheno_col() {
    # <psam or fam> <column> -> "IID value", with 1.9's -9 normalized to NA
    awk -v c="$2" '!/^#/ { v = $c; if (v == "-9") { v = "NA" }; print $2, v }' "$1"
}

# 1. --tail-pheno with both bounds.
plink --bfile tmp_data --pheno tmp_qt.txt --tail-pheno -1 1.5 --make-just-fam --out plink19_t1
$1/plink2 $2 $3 --bfile tmp_data --no-psam-pheno --pheno tmp_qt.txt --tail-pheno -1 1.5 --make-just-psam --out plink2_t1
diff <(pheno_col plink19_t1.fam 6) <(pheno_col plink2_t1.psam 4)

# 2. --tail-pheno with one bound: nothing becomes missing.
plink --bfile tmp_data --pheno tmp_qt.txt --tail-pheno 0 --make-just-fam --out plink19_t2
$1/plink2 $2 $3 --bfile tmp_data --no-psam-pheno --pheno tmp_qt.txt --tail-pheno 0 --make-just-psam --out plink2_t2
diff <(pheno_col plink19_t2.fam 6) <(pheno_col plink2_t2.psam 4)
! grep -q NA <(pheno_col plink2_t2.psam 4)

# 3. --tail-pheno leaves an already-binary phenotype alone.
$1/plink2 $2 $3 --bfile tmp_data --tail-pheno -1 1.5 --make-just-psam --out plink2_cc
$1/plink2 $2 $3 --bfile tmp_data --make-just-psam --out plink2_nocc
diff <(pheno_col plink2_cc.psam 4) <(pheno_col plink2_nocc.psam 4)

# 4. Bad arguments are rejected.
fails() {
    "$@" && exit 1 || true
}
fails $1/plink2 $2 $3 --bfile tmp_data --tail-pheno --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --tail-pheno nonsense --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --tail-pheno 2 1 --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --must-have-sex --out plink2_bad

# 5. An out-of-range number (it starts like a number, so the first line
#    classifies the column as numeric, but it doesn't parse) in the first
#    data line of a .psam or --pheno file is a clean error, not a crash.
$1/plink2 $2 $3 --dummy 3 5 --make-pgen --out tmp_ovf
printf '#IID\tSEX\tP1\nper0\t1\t1e400\nper1\t2\t2\nper2\t1\t1\n' > tmp_ovf_bad.psam
cp tmp_ovf.pgen tmp_ovf_bad.pgen
cp tmp_ovf.pvar tmp_ovf_bad.pvar
if $1/plink2 $2 $3 --pfile tmp_ovf_bad --freq --out plink2_ovf 2> tmp_ovf_err.txt; then
    exit 1
fi
grep -q "Invalid numeric token '1e400'" tmp_ovf_err.txt
printf '#IID\tP2\nper0\t1e999\nper1\t2\nper2\t3\n' > tmp_ovf_pheno.txt
if $1/plink2 $2 $3 --pfile tmp_ovf --pheno tmp_ovf_pheno.txt --freq --out plink2_ovf 2> tmp_ovf_err.txt; then
    exit 1
fi
grep -q "Invalid numeric token '1e999'" tmp_ovf_err.txt
