#!/bin/bash

# PLINK 1.9 commands with no plink2 equivalent, checked against PLINK 1.07:
# --genome, --cluster --K, --mh, --bd, --homog, and QT --assoc --adjust
# (plink2 computes the GC column of t-test results from 1-df chi-square
# quantiles of the p-values, while PLINK 1.07 and 1.9 use T^2 and the t
# distribution).  The reports must be identical.  (Ported from the old
# 1.9/tests_old suite.)
#
# Usage: ./run_tests.sh {plink 1.9 dir}
# plink (1.9) is taken from {plink 1.9 dir}, or from the system PATH if that
# is empty or omitted.  The PLINK 1.07 binary is
# taken from $PLINK107, defaulting to the build in the repository's 1.07/
# directory; every check is skipped if it isn't there.

set -uo pipefail

# plink 1.9 from {plink 1.9 dir}, or from the PATH if that is empty or omitted.
plink19="${1:+$1/}plink"
plink() {
    command "$plink19" "$@"
}

plink107="${PLINK107:-../../../1.07/plink}"
checks=(
    "--genome matches PLINK 1.07"
    "--cluster only2 old-tiebreaks --K 2 matches PLINK 1.07 --cluster --K 2"
    "--mh matches PLINK 1.07"
    "--bd matches PLINK 1.07"
    "--homog matches PLINK 1.07"
    "QT --assoc --adjust matches PLINK 1.07"
    "QT --assoc --adjust --lambda matches PLINK 1.07"
)
if [[ ! -x "$plink107" ]]; then
    echo "note: PLINK 1.07 not found at $plink107; skipping checks that need it"
    for c in "${checks[@]}"; do
        echo "skipped: $c"
    done
    exit 0
fi

failed=0
check() {
    if [[ "$2" -eq 0 ]]; then
        echo "ok: $1"
    else
        echo "FAILED: $1"
        failed=1
    fi
}
p107() {
    "$plink107" --noweb --silent "$@" > /dev/null
}

# --genome and --cluster: 101 samples (101 % 4 == 1), 500 variants.
plink --dummy 101 500 0.02 acgt --seed 23 --make-bed --out tmp_rel > /dev/null
awk 'BEGIN { OFS = "\t" } { $4 = 1000 * NR; print }' tmp_rel.bim > tmp_rel.bim.new
mv tmp_rel.bim.new tmp_rel.bim

p107 --bfile tmp_rel --genome --out tmp_genome107
plink --bfile tmp_rel --genome --out tmp_genome19 > /dev/null
cmp -s tmp_genome107.genome tmp_genome19.genome && test "$(wc -l < tmp_genome19.genome)" -eq 5051
check "${checks[0]}" $?

p107 --bfile tmp_rel --cluster --K 2 --out tmp_k107
plink --bfile tmp_rel --cluster only2 old-tiebreaks --K 2 --out tmp_k19 > /dev/null
cmp -s tmp_k107.cluster2 tmp_k19.cluster2 && test "$(cut -f 2 tmp_k19.cluster2 | sort -u | wc -l)" -eq 2
check "${checks[1]}" $?

# --mh, --bd and --homog: 513 cases and controls in 4 strata.  The strata are
# named in order of first appearance, since PLINK 1.07 reports them in that
# order and PLINK 1.9 sorts them.
plink --dummy 513 300 0.02 acgt --seed 29 --make-bed --out tmp_cc > /dev/null
awk 'BEGIN { OFS = "\t" } { $4 = 1000 * NR; print }' tmp_cc.bim > tmp_cc.bim.new
mv tmp_cc.bim.new tmp_cc.bim
awk '{ print $1, $2, "s" ((NR - 1) % 4) }' tmp_cc.fam > tmp_strata.txt

p107 --bfile tmp_cc --within tmp_strata.txt --mh --out tmp_mh107
plink --bfile tmp_cc --within tmp_strata.txt --mh --out tmp_mh19 > /dev/null
cmp -s tmp_mh107.cmh tmp_mh19.cmh
check "${checks[2]}" $?

p107 --bfile tmp_cc --within tmp_strata.txt --bd --out tmp_bd107
plink --bfile tmp_cc --within tmp_strata.txt --bd --out tmp_bd19 > /dev/null
cmp -s tmp_bd107.cmh tmp_bd19.cmh
check "${checks[3]}" $?

p107 --bfile tmp_cc --within tmp_strata.txt --homog --out tmp_homog107
plink --bfile tmp_cc --within tmp_strata.txt --homog --out tmp_homog19 > /dev/null
cmp -s tmp_homog107.homog tmp_homog19.homog && test "$(wc -l < tmp_homog19.homog)" -eq 2101
check "${checks[4]}" $?

# QT --assoc --adjust: 387 samples (387 % 4 == 3), 200 variants.  The
# estimated lambda is above 1 here, and --lambda 1.2 checks a fixed value.
plink --dummy 387 200 0.03 acgt scalar-pheno --seed 11 --make-bed --out tmp_qt > /dev/null
awk 'BEGIN { OFS = "\t" } { $4 = 1000 * NR; print }' tmp_qt.bim > tmp_qt.bim.new
mv tmp_qt.bim.new tmp_qt.bim

p107 --bfile tmp_qt --assoc --adjust --out tmp_qadj107
plink --bfile tmp_qt --assoc --adjust --out tmp_qadj19 > /dev/null
cmp -s tmp_qadj107.qassoc.adjusted tmp_qadj19.qassoc.adjusted && test "$(wc -l < tmp_qadj19.qassoc.adjusted)" -eq 201 && \
    grep -q 'lambda (based on median chisq) = 1\.0[1-9]' tmp_qadj19.log
check "${checks[5]}" $?

p107 --bfile tmp_qt --assoc --adjust --lambda 1.2 --out tmp_qlam107
plink --bfile tmp_qt --assoc --adjust --lambda 1.2 --out tmp_qlam19 > /dev/null
cmp -s tmp_qlam107.qassoc.adjusted tmp_qlam19.qassoc.adjusted && test "$(wc -l < tmp_qlam19.qassoc.adjusted)" -eq 201
check "${checks[6]}" $?

exit $failed
