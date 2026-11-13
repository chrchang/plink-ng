#!/bin/bash

# --glm output and multicollinearity flags with no other test coverage:
# --ci, --vif, --max-corr, --perm-format, --perm-count, and --mperm-save.

set -exo pipefail

plink2="$1/plink2 $2 $3"

python3 make_data.py
$plink2 --vcf tmp_data.vcf --psam tmp_data.psam --make-pgen --out tmp_data > /dev/null

# --ci: normal-approximation intervals for linear, logistic, and Firth.
for ci in 0.9 0.95; do
    pct=$(python3 -c "print(round(100 * $ci))")
    $plink2 --pfile tmp_data --pheno-name QT,CC --covar-name C1 --glm --ci $ci --out tmp_ci$pct > /dev/null
    python3 check_glm_out.py ci tmp_ci$pct.QT.glm.linear $ci
    python3 check_glm_out.py ci tmp_ci$pct.CC.glm.logistic.hybrid $ci
    $plink2 --pfile tmp_data --pheno-name CC --covar-name C1 --glm firth --ci $ci --out tmp_cif$pct > /dev/null
    python3 check_glm_out.py ci tmp_cif$pct.CC.glm.firth $ci
done
# No CI columns without --ci.
$plink2 --pfile tmp_data --pheno-name QT --covar-name C1 --glm --out tmp_noci > /dev/null
if head -n 1 tmp_noci.QT.glm.linear | grep -q 'L95'; then
    exit 1
fi
for bad in 0.005 1 1.5; do
    if $plink2 --pfile tmp_data --pheno-name QT --covar-name C1 --glm --ci $bad --out tmp_cibad > /dev/null 2>&1; then
        exit 1
    fi
done

# Per-variant checks: v2 and C5 have correlation r (about 0.993), and with
# just those two predictors, VIF = 1 / (1 - r^2) (about 74).  Thresholds
# just below and above r and the VIF bracket the result.  The correlation
# check comes first.
read -r r vif < <(python3 check_glm_out.py corr tmp_data.psam tmp_data.vcf v2 C5)
vif_lo=$(python3 -c "print($vif * 0.999)")
vif_hi=$(python3 -c "print($vif * 1.001)")
r_lo=$(python3 -c "print($r - 0.0001)")
r_hi=$(python3 -c "print($r + 0.0001)")
# $1: expected v2 ERRCODE; $2: phenotype; remaining args: extra flags.
check_v2_errcode() {
    local expected=$1
    local pheno=$2
    shift 2
    $plink2 --pfile tmp_data --pheno-name $pheno --covar-name C5 "$@" --out tmp_vc > /dev/null
    local fname
    fname=$(ls tmp_vc.$pheno.glm.*)
    # Every other variant is unaffected.
    test "$(awk -F '\t' 'NR > 1 {print $3 ":" $NF}' "$fname" | tr '\n' ' ')" = "v1:. v2:$expected v3:. v4:. v5:. v6:. "
    rm -f tmp_vc.*
}
for glm in "--glm hide-covar" "--glm hide-covar firth"; do
    pheno=QT
    [[ $glm == *firth ]] && pheno=CC
    check_v2_errcode VIF_TOO_HIGH $pheno $glm
    check_v2_errcode VIF_TOO_HIGH $pheno $glm --vif "$vif_lo"
    check_v2_errcode . $pheno $glm --vif "$vif_hi"
    check_v2_errcode CORR_TOO_HIGH $pheno $glm --vif "$vif_hi" --max-corr "$r_lo"
    check_v2_errcode . $pheno $glm --vif "$vif_hi" --max-corr "$r_hi"
    check_v2_errcode VIF_TOO_HIGH $pheno $glm --vif "$vif_lo" --max-corr "$r_hi"
done
# Logistic regression doesn't skip the VIF check either.
check_v2_errcode VIF_TOO_HIGH CC --glm hide-covar
check_v2_errcode . CC --glm hide-covar --vif "$vif_hi"

# Covariate-only checks happen before any regression, and are fatal.  C1
# and C2 have VIF about 131 and correlation about 0.996.
glm_c1c2() {
    $plink2 --pfile tmp_data --pheno-name QT --covar-name C1,C2 --glm "$@" --out tmp_c1c2 > /dev/null 2>&1
}
if glm_c1c2; then
    exit 1
fi
grep -q 'VIF_TOO_HIGH' tmp_c1c2.log
glm_c1c2 --vif 200
if glm_c1c2 --vif 200 --max-corr 0.99; then
    exit 1
fi
grep -q 'CORR_TOO_HIGH' tmp_c1c2.log

# Invalid values, and use without --glm.
for flags in "--glm --vif 0.5" "--glm --max-corr 1.5" "--glm --max-corr -0.1" "--freq --vif 100" "--freq --max-corr 0.9"; do
    if $plink2 --pfile tmp_data --pheno-name QT --covar-name C1 $flags --out tmp_bad > /dev/null 2>&1; then
        exit 1
    fi
done

# Permutation output.  Every run with the same --seed draws the same
# permutations.
perm_run() {
    local out=$1
    shift
    $plink2 --pfile tmp_data --pheno-name QT --covar-name C1 --seed 3 "$@" --out $out > /dev/null
}
perm_run tmp_m --glm mperm=50
head -n 1 tmp_m.QT.glm.linear.mperm | diff - <(printf '#CHROM\tID\tREF\tALT\tA1\tOMITTED\tEMP1\tEMP2\n')
perm_run tmp_mc --glm mperm=50 --perm-format counts
head -n 1 tmp_mc.QT.glm.linear.mperm | diff - <(printf '#CHROM\tID\tREF\tALT\tA1\tOMITTED\tEMP1_CT\tEMP2_CT\n')
python3 check_glm_out.py counts tmp_mc.QT.glm.linear.mperm tmp_m.QT.glm.linear.mperm 50
# --perm-count is an alias for '--perm-format counts'.
perm_run tmp_mc2 --glm mperm=50 --perm-count
diff tmp_mc.QT.glm.linear.mperm tmp_mc2.QT.glm.linear.mperm

# 'zs' and cols=.
perm_run tmp_mz --glm mperm=50 --perm-format zs cols=pos,alt1,omitted
perm_run tmp_mz2 --glm mperm=50 --perm-format cols=pos,alt1,omitted
zstd -dc tmp_mz.QT.glm.linear.mperm.zst | diff - tmp_mz2.QT.glm.linear.mperm
head -n 1 tmp_mz2.QT.glm.linear.mperm | diff - <(printf '#POS\tID\tALT1\tA1\tOMITTED\tEMP1\tEMP2\n')
diff <(cut -f 6- tmp_mz2.QT.glm.linear.mperm) <(cut -f 7- tmp_m.QT.glm.linear.mperm)
perm_run tmp_mp --glm mperm=50 --perm-format counts cols=chrom,ref,alt,provref
head -n 1 tmp_mp.QT.glm.linear.mperm | diff - <(printf '#CHROM\tID\tREF\tALT\tPROVISIONAL_REF?\tA1\tEMP1_CT\tEMP2_CT\n')
test "$(tail -n +2 tmp_mp.QT.glm.linear.mperm | cut -f 5 | sort -u)" = "N"
for flags in "--perm-format counts --perm-format zs" "--perm-format cols=chrom cols=pos" "--perm-format bad"; do
    if perm_run tmp_pbad --glm mperm=50 $flags 2> /dev/null; then
        exit 1
    fi
done

# Adaptive permutation: PERM_CT varies, and EMP1 = (EMP1_CT+1)/(PERM_CT+1).
perm_run tmp_a --glm aperm --aperm 10 200
perm_run tmp_ac --glm aperm --aperm 10 200 --perm-format counts
head -n 1 tmp_ac.QT.glm.linear.aperm | diff - <(printf '#CHROM\tID\tREF\tALT\tA1\tOMITTED\tEMP1_CT\tPERM_CT\n')
diff <(cut -f 8 tmp_a.QT.glm.linear.aperm) <(cut -f 8 tmp_ac.QT.glm.linear.aperm)
python3 check_glm_out.py counts tmp_ac.QT.glm.linear.aperm tmp_a.QT.glm.linear.aperm col

# --mperm-save writes just the per-permutation maximum statistics, which
# match the maxima of the --mperm-save-all statistics.
perm_run tmp_ms --glm mperm=50 --mperm-save --perm-count
test -s tmp_ms.QT.glm.linear.mperm.dump.best
test ! -e tmp_ms.QT.glm.linear.mperm.dump.all
perm_run tmp_msa --glm mperm=50 --mperm-save-all --perm-count
test ! -e tmp_msa.QT.glm.linear.mperm.dump.best
test "$(wc -l < tmp_ms.QT.glm.linear.mperm.dump.best)" -eq 51
python3 check_glm_out.py dump tmp_msa.QT.glm.linear.mperm.dump.all tmp_ms.QT.glm.linear.mperm.dump.best tmp_ms.QT.glm.linear.mperm
perm_run tmp_msb --glm mperm=50 --mperm-save --mperm-save-all
diff tmp_ms.QT.glm.linear.mperm.dump.best tmp_msb.QT.glm.linear.mperm.dump.best
diff tmp_msa.QT.glm.linear.mperm.dump.all tmp_msb.QT.glm.linear.mperm.dump.all
