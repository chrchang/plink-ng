#!/bin/bash

# --glm linear regression with covariates on variants that have missing
# genotype calls.
#
# Each variant's rows are checked against an oracle run that extracts that
# variant alone and removes its missing-genotype samples up front, so the
# oracle regression sees no missing calls at all and takes the no-missing
# code path.  Both runs fit the same model to the same samples.
#
# The phenotype file has three columns with one missingness pattern, which
# --glm fits together as a batch, and a fourth with its own pattern, which it
# fits on its own, so both linear regression drivers are covered.

set -exo pipefail

plink2="$1/plink2 $2 $3"

# 40 variants, each missing 0.2%, 2%, 20% or 60% of its calls.
$plink2 --dummy 600 40 0.002,0.02,0.2,0.6 acgt --seed 5 --out tmp_data > /dev/null

# Five quantitative covariates, two 0/1 indicators of a three-level group, and
# two samples with a missing covariate.  (A categorical covariate would not do:
# its reference level can move when the oracle removes samples.)
awk 'BEGIN { srand(7); OFS = "\t"; print "#IID", "C1", "C2", "C3", "C4", "C5", "G1", "G2" }
     NR > 1 {
       line = $1
       for (j = 1; j <= 5; ++j) {
         line = line OFS ((NR == 20 && j == 3)? "NA" : sprintf("%.6f", j * (rand() - 0.5) + (j == 2) * 3))
       }
       print line, ((NR == 45)? "NA" : ((NR % 3 == 1)? 1 : 0)), ((NR % 3 == 2)? 1 : 0)
     }' tmp_data.psam > tmp_covar.txt

# Q1-Q3 share their missing values; Q4 has its own.
awk 'BEGIN { srand(11); OFS = "\t"; print "#IID", "Q1", "Q2", "Q3", "Q4" }
     NR > 1 {
       q1 = (NR % 37 == 0)? "NA" : sprintf("%.6f", rand() * 4 - 2)
       q2 = (NR % 37 == 0)? "NA" : sprintf("%.6f", rand() + NR / 600)
       q3 = (NR % 37 == 0)? "NA" : sprintf("%.6f", rand() * rand() * 10)
       q4 = (NR % 29 == 0)? "NA" : sprintf("%.6f", rand() - 0.5)
       print $1, q1, q2, q3, q4
     }' tmp_data.psam > tmp_pheno.txt

# The samples missing each variant's call, from an additive export.
$plink2 --pfile tmp_data --export A --out tmp_raw > /dev/null
rm -f tmp_miss_*.txt
awk 'NR == 1 { for (j = 7; j <= NF; ++j) { id[j] = $j; sub(/_[^_]*$/, "", id[j]) }; next }
     { for (j = 7; j <= NF; ++j) { if ($j == "NA") { print $2 > ("tmp_miss_" id[j] ".txt") } } }' tmp_raw.raw
tail -n +2 tmp_data.pvar | cut -f 3 > tmp_variants.txt

# Most variants must sit below the 50% missingness mark, where the missing
# samples are cheaper to take back out than the rest are to copy, and a few
# above it.
for v in $(cat tmp_variants.txt); do
    if [ -f tmp_miss_$v.txt ]; then
        wc -l < tmp_miss_$v.txt
    else
        echo 0
    fi
done > tmp_miss_counts.txt
test "$(awk '$1 > 0 && $1 < 300' tmp_miss_counts.txt | wc -l)" -ge 20
test "$(awk '$1 >= 300' tmp_miss_counts.txt | wc -l)" -ge 3

check_model() {
    $plink2 --pfile tmp_data --pheno tmp_pheno.txt --pheno-name Q1-Q4 --covar tmp_covar.txt --glm $1 --out tmp_main > /dev/null
    rm -f tmp_oracle.Q?.glm.linear
    for v in $(cat tmp_variants.txt); do
        echo $v > tmp_one.txt
        if [ -f tmp_miss_$v.txt ]; then
            (echo "#IID"; cat tmp_miss_$v.txt) > tmp_remove.txt
            $plink2 --pfile tmp_data --extract tmp_one.txt --remove tmp_remove.txt --pheno tmp_pheno.txt --pheno-name Q1-Q4 --covar tmp_covar.txt --glm $1 --out tmp_one > /dev/null
        else
            $plink2 --pfile tmp_data --extract tmp_one.txt --pheno tmp_pheno.txt --pheno-name Q1-Q4 --covar tmp_covar.txt --glm $1 --out tmp_one > /dev/null
        fi
        for q in Q1 Q2 Q3 Q4; do
            if [ -f tmp_oracle.$q.glm.linear ]; then
                tail -n +2 tmp_one.$q.glm.linear >> tmp_oracle.$q.glm.linear
            else
                cp tmp_one.$q.glm.linear tmp_oracle.$q.glm.linear
            fi
        done
    done
    for q in Q1 Q2 Q3 Q4; do
        awk -f compare.awk tmp_oracle.$q.glm.linear tmp_main.$q.glm.linear
    done
}

# 1. Additive, with the covariate rows reported too.
check_model ""

# 2. The genotype codings that change the genotype column.
check_model "hide-covar dominant"
check_model "hide-covar recessive"
check_model "hide-covar hetonly"

# 3. The two models with a dominance deviation column, which also report a
#    joint test.
check_model "genotypic"
check_model "hide-covar hethom"

# 4. A covariate that goes constant once a variant's missing-genotype samples
#    drop out.  Its variance over the remaining samples is zero, which the
#    covariate dot products have to report as such rather than as the rounding
#    noise a subtraction would leave, so these variants take the generic path.
#    The oracle cannot check this one: removing those samples up front makes
#    --glm drop the covariate as constant, which is a different model.
var=$(awk 'NR == FNR { n[FNR] = $1; next } { if (n[FNR] > 50 && n[FNR] < 250) { print $1; exit } }' tmp_miss_counts.txt tmp_variants.txt)
test -n "$var"
awk -v var="$var" 'BEGIN { OFS = "\t"; while ((getline line < ("tmp_miss_" var ".txt")) > 0) { miss[line] = 1 } }
     FNR == 1 { print $1, $2, "CONST_ON_NM"; next }
     { print $1, $2, ($1 in miss)? 1 : 0 }' tmp_covar.txt > tmp_const_covar.txt
$plink2 --pfile tmp_data --pheno tmp_pheno.txt --pheno-name Q4 --covar tmp_const_covar.txt --glm hide-covar --out tmp_const > /dev/null
awk -v var="$var" 'FNR > 1 && $3 == var { if ($NF != "VIF_INFINITE") { print "expected VIF_INFINITE on " var ", got " $NF; exit 1 }; print var " reported " $NF; ++n }
     END { if (n != 1) { print "expected one row for " var ", got " n + 0; exit 1 } }' tmp_const.Q4.glm.linear
