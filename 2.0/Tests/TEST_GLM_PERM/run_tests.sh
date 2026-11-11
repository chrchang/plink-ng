#!/bin/bash

# --glm permutation tests ('aperm' and 'mperm='), linear and Firth.
#
# Every permutation run with the same --seed draws the same phenotype
# permutations (they come from a single random stream, independent of
# threading and batching), so a --mperm-save-all dump from a max(T) run gives
# the permutation statistics behind an adaptive run as well.
#
# Regression checks:
#   - With chrY present, permutation crashed whenever chrX used the same
#     samples and covariates as the rest of the genome (e.g. with 'sex' or
#     'no-x-sex').
#   - Linear adaptive permutation skipped the stopping check for every allele
#     of a multiallelic variant after the first one, as long as an earlier
#     allele was still active.
#   - Ties between a permutation statistic and the original statistic were
#     detected with exact equality, but the two are computed by different
#     code, so exact ties often differed in the last bits and counted as 0 or
#     1 instead of 1/2.  Identity permutations (--permute-within with one
#     sample per category) make every permutation statistic a tie.

set -exo pipefail

plink2="$1/plink2 $2 $3"

python3 make_data.py
$plink2 --vcf tmp_data.vcf --psam tmp_data.psam --make-pgen --out tmp_data > /dev/null

covars="--covar tmp_data.psam --covar-name COV1,COV2"
# Runs --glm on phenotype(s) $2 with output prefix tmp_$1; further arguments
# are --glm modifiers, then any other flags after a "--".  (Categorical
# phenotypes such as CL are skipped by --glm.)
glm() {
    local name=$1 pheno=$2
    shift 2
    local mods=()
    while [[ $# -gt 0 && $1 != "--" ]]; do
        mods+=("$1")
        shift
    done
    [[ $# -gt 0 ]] && shift
    $plink2 --pfile tmp_data --pheno-name "$pheno" --glm perm-count "${mods[@]}" $covars --seed 11 "$@" --out "tmp_$name" > /dev/null
}

# 1. chrX with the main sample/covariate set, plus chrY.
for mod in sex no-x-sex; do
    glm lin_$mod QT $mod mperm=200 --mperm-save-all
    python3 check_perm.py report tmp_lin_$mod.QT.glm.linear.mperm.dump.all tmp_lin_$mod.QT.glm.linear.mperm
    glm lin_${mod}_a QT $mod aperm -- --aperm 5 200
    test -s tmp_lin_${mod}_a.QT.glm.linear.aperm
    glm firth_$mod CC $mod firth mperm=100 --mperm-save-all
    python3 check_perm.py report tmp_firth_$mod.CC.glm.firth.mperm.dump.all tmp_firth_$mod.CC.glm.firth.mperm
done

# 2. Adaptive permutation, replayed from a max(T) dump.  The triallelic
# variants on chr2 have effects from one or both ALT alleles.
glm adapt QT aperm -- --aperm 5 4000
allele_test_ct=$(grep -o '([0-9]* allele tests' tmp_adapt.log | tr -dc 0-9)
glm adapt_m QT mperm=4000 --mperm-save-all
python3 check_perm.py adaptive tmp_adapt_m.QT.glm.linear.mperm.dump.all tmp_adapt.QT.glm.linear.aperm "$allele_test_ct" 4000
glm fadapt CC firth aperm -- --aperm 5 1000
glm fadapt_m CC firth mperm=1000 --mperm-save-all
python3 check_perm.py adaptive tmp_fadapt_m.CC.glm.firth.mperm.dump.all tmp_fadapt.CC.glm.firth.aperm "$allele_test_ct" 1000

# 3. Ties.
glm ident QT,CL mperm=40 -- --permute-within CL
python3 check_perm.py identity tmp_ident.QT.glm.linear.mperm 40
glm fident CC,CL firth mperm=40 -- --permute-within CL
python3 check_perm.py identity tmp_fident.CC.glm.firth.mperm 40
# No covariates here: then a permutation's statistic depends only on its
# genotype x phenotype table, and near ties are exact ties.
covars=""
glm two TWO allow-no-covars mperm=1000 --mperm-save-all
python3 check_perm.py ties tmp_two.TWO.glm.linear.mperm.dump.all tmp_two.TWO.glm.linear.mperm
