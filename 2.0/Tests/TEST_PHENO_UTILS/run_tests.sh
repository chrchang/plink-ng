#!/bin/bash

# --tail-pheno and --must-have-sex, checked against PLINK 1.9.

set -exo pipefail

plink --dummy 200 400 0.03 --out tmp_data > /dev/null

# A quantitative phenotype, and a sex column with a third of the samples unknown.
awk 'BEGIN{OFS=" "; print "#FID", "IID", "QT"} {print $1, $2, ($NR % 17) - 6 + (NR % 5) * 0.25}' tmp_data.fam > tmp_qt.txt
awk 'BEGIN{OFS=" "} {$5 = (NR % 3 == 0)? 0 : (1 + (NR % 2)); print}' tmp_data.fam > tmp_sexed.fam
mv tmp_sexed.fam tmp_data.fam

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

# 3. --must-have-sex.
plink --bfile tmp_data --allow-no-sex --must-have-sex --make-just-fam --out plink19_m
$1/plink2 $2 $3 --bfile tmp_data --must-have-sex --make-just-psam --out plink2_m
diff <(pheno_col plink19_m.fam 6) <(pheno_col plink2_m.psam 4)
# every unknown-sex sample has a missing phenotype, and no known-sex sample
# lost one
awk '!/^#/ { if ($3 == "NA" && $4 != "NA") { print "unknown sex kept a phenotype: " $2; exit 1 } }' plink2_m.psam
diff <(pheno_col plink2_m.psam 4 | grep -v ' NA$' | cut -d' ' -f1) \
     <(awk '!/^#/ && $3 != "NA" { print $2 }' plink2_m.psam)

# 4. The two together, in that order: --must-have-sex first, so a sample with
#    unknown sex is missing rather than downcoded.
plink --bfile tmp_data --allow-no-sex --pheno tmp_qt.txt --must-have-sex --tail-pheno -1 1.5 --make-just-fam --out plink19_mt
$1/plink2 $2 $3 --bfile tmp_data --no-psam-pheno --pheno tmp_qt.txt --must-have-sex --tail-pheno -1 1.5 --make-just-psam --out plink2_mt
diff <(pheno_col plink19_mt.fam 6) <(pheno_col plink2_mt.psam 4)

# 5. --tail-pheno leaves an already-binary phenotype alone.
$1/plink2 $2 $3 --bfile tmp_data --tail-pheno -1 1.5 --make-just-psam --out plink2_cc
$1/plink2 $2 $3 --bfile tmp_data --make-just-psam --out plink2_nocc
diff <(pheno_col plink2_cc.psam 4) <(pheno_col plink2_nocc.psam 4)

# 6. Bad arguments are rejected.
fails() {
    "$@" && exit 1 || true
}
fails $1/plink2 $2 $3 --bfile tmp_data --tail-pheno --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --tail-pheno nonsense --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --tail-pheno 2 1 --out plink2_bad
