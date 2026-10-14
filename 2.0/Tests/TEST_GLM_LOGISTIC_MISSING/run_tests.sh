#!/bin/bash

# --glm logistic and Firth regression on variants with missing genotype calls
# that follow a variant with none.
#
# After a variant with no missing call, the intercept column is left filled
# over the whole sample count.  A following variant with missing calls has a
# shorter column, whose padding past its own sample count must be zero; if the
# earlier 1s are left there, the regression gains a phantom intercept-only
# control, and the reported standard errors and p-values drift.  How much
# padding there is depends on the vector width, and which variant follows
# which depends on how --threads splits the block, so the damage used to vary
# with both.
#
# Each variant's rows are checked against an oracle run that extracts that
# variant alone, so the oracle regression never follows another variant.

set -exo pipefail

plink2="$1/plink2 $2 $3"

# 603 samples, 30 variants.  Variants alternate between no missing call and a
# few blanked ones, so that every variant with missing calls directly follows
# one without.  The blank counts leave both odd and even nonmissing sample
# counts, none a multiple of 4, so there is padding for both the single- and
# double-precision paths at SSE, AVX2 and NEON widths.
$plink2 --dummy 603 30 acgt --seed 3 --out tmp_dense > /dev/null
$plink2 --pfile tmp_dense --export vcf --out tmp_dense > /dev/null
awk 'BEGIN { OFS = "\t"; split("1 2 6 9 14", blank_ct, " ") }
     /^#/ { print; next }
     {
       cur_blank_ct = (vidx % 2)? blank_ct[((vidx - 1) / 2) % 5 + 1] : 0
       for (s = 10; s <= NF; ++s) {
         # 37 and 603 are coprime, so this is a permutation of the samples.
         if (((s - 10) * 37 + vidx * 11) % 603 < cur_blank_ct) { $s = "./." }
       }
       ++vidx
       print
     }' tmp_dense.vcf > tmp_blanked.vcf
$plink2 --vcf tmp_blanked.vcf --make-pgen --out tmp_data > /dev/null

# Two quantitative covariates, and a case/control phenotype.
awk 'BEGIN { srand(7); OFS = "\t"; print "#IID", "C1", "C2" }
     NR > 1 { print $1, sprintf("%.6f", rand() - 0.5), sprintf("%.6f", 3 * rand()) }' tmp_data.psam > tmp_covar.txt
awk 'BEGIN { srand(11); OFS = "\t"; print "#IID", "B1" }
     NR > 1 { print $1, (rand() < 0.4)? 2 : 1 }' tmp_data.psam > tmp_pheno.txt
grep -v '^#' tmp_data.pvar | cut -f 3 > tmp_variants.txt

# $1: --glm modifiers.  $2: the output file suffix.
check_model() {
    for t in 1 2 4; do
        $plink2 --pfile tmp_data --pheno tmp_pheno.txt --covar tmp_covar.txt --glm $1 --threads $t --out tmp_main_t$t > /dev/null
    done
    rm -f tmp_oracle.B1.$2
    for v in $(cat tmp_variants.txt); do
        echo $v > tmp_one.txt
        $plink2 --pfile tmp_data --extract tmp_one.txt --pheno tmp_pheno.txt --covar tmp_covar.txt --glm $1 --out tmp_one > /dev/null
        if [ -f tmp_oracle.B1.$2 ]; then
            tail -n +2 tmp_one.B1.$2 >> tmp_oracle.B1.$2
        else
            cp tmp_one.B1.$2 tmp_oracle.B1.$2
        fi
    done
    for t in 1 2 4; do
        awk -f compare.awk tmp_oracle.B1.$2 tmp_main_t$t.B1.$2
    done
}

# Double-precision logistic regression, with Firth fallback.
check_model "" glm.logistic.hybrid
# Single-precision logistic regression.
check_model "single-prec-cc" glm.logistic.hybrid
# Firth regression throughout, in both precisions.
check_model "firth" glm.firth
check_model "firth single-prec-cc" glm.firth
