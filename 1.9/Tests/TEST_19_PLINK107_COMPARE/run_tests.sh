#!/bin/bash

# PLINK 1.9 commands with no plink2 equivalent, checked against PLINK 1.07:
# --genome, --cluster --K, --mh, --bd, --homog, and QT --assoc --adjust
# (plink2 computes the GC column of t-test results from 1-df chi-square
# quantiles of the p-values, while PLINK 1.07 and 1.9 use T^2 and the t
# distribution).  The reports must be identical.  (Ported from the old
# 1.9/tests_old suite.)
#
# Then several commands are run on data with heterozygous haploid calls.
# PLINK 1.07 sets these to missing when it loads the data, while PLINK 1.9
# leaves them in place and each command must mask them itself.  plink2 has
# no heterozygous haploid state, so PLINK 1.07 is the only reference here.
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
    "het. haploid data: --freq matches PLINK 1.07"
    "het. haploid data: --check-sex matches PLINK 1.07"
    "het. haploid data: --assoc, --model and --assoc fisher match PLINK 1.07"
    "het. haploid data: --mh, --bd and --homog match PLINK 1.07"
    "het. haploid data: --hardy --allow-no-sex matches PLINK 1.07"
    "het. haploid data: QT --assoc and --linear match PLINK 1.07"
    "het. haploid data: --mendel matches PLINK 1.07"
    "het. haploid data: --tdt and --tdt poo match PLINK 1.07"
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


# --- Heterozygous haploid calls ---

# Writes text with leading/trailing blanks removed and runs of spaces
# squeezed, so that column widths don't matter.
norm() {
    sed 's/^ *//; s/ *$//' "$1" | tr -s ' '
}

# 201 samples: every tenth has unknown sex, the rest alternate between male
# and female.  250 autosomal, 100 chrX and 50 chrY variants, so males have
# heterozygous chrX and chrY calls, and females have chrY calls.
plink --dummy 201 400 0.02 acgt --seed 31 --make-bed --out tmp_hh > /dev/null
awk 'BEGIN { OFS = "\t" } { $4 = 1000 * NR; if (NR > 250) $1 = 23; if (NR > 350) $1 = 24; print }' tmp_hh.bim > tmp_hh.bim.new
mv tmp_hh.bim.new tmp_hh.bim
awk '{ $5 = (NR % 10 == 0)? 0 : 1 + NR % 2; print }' tmp_hh.fam > tmp_hh.fam.new
mv tmp_hh.fam.new tmp_hh.fam
awk '{ print $1, $2, "s" ((NR - 1) % 3) }' tmp_hh.fam > tmp_hh_strata.txt
# The same genotypes with a quantitative phenotype.
cp tmp_hh.bed tmp_hhq.bed
cp tmp_hh.bim tmp_hhq.bim
awk 'BEGIN { s = 7 } { s = (s * 69069 + 1) % 4294967296; $6 = int(s / 65536) % 1000 / 100; print }' tmp_hh.fam > tmp_hhq.fam
# Make sure the data has what this section is about.
plink --bfile tmp_hh --freq --out tmp_hh_frq19 > /dev/null 2>&1
grep -q 'het\. haploid genotypes present' tmp_hh_frq19.log || echo "warning: no heterozygous haploid calls in the test data"

p107 --bfile tmp_hh --freq --out tmp_hh_frq107
cmp -s tmp_hh_frq107.frq tmp_hh_frq19.frq
check "${checks[7]}" $?

p107 --bfile tmp_hh --check-sex --out tmp_hh_cs107
plink --bfile tmp_hh --check-sex --out tmp_hh_cs19 > /dev/null 2>&1
cmp -s tmp_hh_cs107.sexcheck tmp_hh_cs19.sexcheck
check "${checks[8]}" $?

# (PLINK 1.07 spells --assoc fisher as --fisher.)
hh_ok=0
p107 --bfile tmp_hh --assoc --out tmp_hh_as107
plink --bfile tmp_hh --assoc --out tmp_hh_as19 > /dev/null 2>&1
cmp -s tmp_hh_as107.assoc tmp_hh_as19.assoc || hh_ok=1
p107 --bfile tmp_hh --model --out tmp_hh_md107
plink --bfile tmp_hh --model --out tmp_hh_md19 > /dev/null 2>&1
cmp -s tmp_hh_md107.model tmp_hh_md19.model || hh_ok=1
p107 --bfile tmp_hh --fisher --out tmp_hh_fi107
plink --bfile tmp_hh --assoc fisher --out tmp_hh_fi19 > /dev/null 2>&1
cmp -s tmp_hh_fi107.assoc.fisher tmp_hh_fi19.assoc.fisher || hh_ok=1
check "${checks[9]}" $hh_ok

hh_ok=0
for t in mh bd; do
    p107 --bfile tmp_hh --within tmp_hh_strata.txt --$t --out tmp_hh_$t"107"
    plink --bfile tmp_hh --within tmp_hh_strata.txt --$t --out tmp_hh_$t"19" > /dev/null 2>&1
    cmp -s tmp_hh_$t"107".cmh tmp_hh_$t"19".cmh || hh_ok=1
done
p107 --bfile tmp_hh --within tmp_hh_strata.txt --homog --out tmp_hh_hg107
plink --bfile tmp_hh --within tmp_hh_strata.txt --homog --out tmp_hh_hg19 > /dev/null 2>&1
cmp -s tmp_hh_hg107.homog tmp_hh_hg19.homog || hh_ok=1
check "${checks[10]}" $hh_ok

# --hardy: --allow-no-sex, since PLINK 1.07 counts unknown-sex samples in the
# AFF/UNAFF rows by their original phenotype even though it treats those
# phenotypes as missing elsewhere.  chrY is left out: it has no diploid
# genotypes, and PLINK 1.07 prints NA for the AFF/UNAFF p-values where 1.9
# prints 1.  The rest must match to printing precision.
p107 --bfile tmp_hh --allow-no-sex --hardy --out tmp_hh_hw107
plink --bfile tmp_hh --allow-no-sex --hardy --out tmp_hh_hw19 > /dev/null 2>&1
awk -v c1=1,2,3,4,5,6,7,8,9 -v c2=1,2,3,4,5,6,7,8,9 -f ../compare_cols.awk <(awk '$1 != 24' tmp_hh_hw107.hwe) <(awk '$1 != 24' tmp_hh_hw19.hwe)
check "${checks[11]}" $?

hh_ok=0
p107 --bfile tmp_hhq --assoc --out tmp_hh_qa107
plink --bfile tmp_hhq --assoc --out tmp_hh_qa19 > /dev/null 2>&1
awk -v c1=1,2,3,4,5,6,7,8,9 -v c2=1,2,3,4,5,6,7,8,9 -f ../compare_cols.awk tmp_hh_qa107.qassoc tmp_hh_qa19.qassoc || hh_ok=1
p107 --bfile tmp_hhq --linear --out tmp_hh_li107
plink --bfile tmp_hhq --linear --out tmp_hh_li19 > /dev/null 2>&1
awk -v c1=1,2,3,4,5,6,7,8,9 -v c2=1,2,3,4,5,6,7,8,9 -f ../compare_cols.awk tmp_hh_li107.assoc.linear tmp_hh_li19.assoc.linear || hh_ok=1
check "${checks[12]}" $hh_ok

# --mendel and --tdt: 67 trios (201 samples) over 100 autosomal and 100 chrX
# variants.  make_trios.awk generates diploid chrX genotypes for everyone, so
# fathers and sons have heterozygous chrX calls.  PLINK 1.07 orders the
# --mendel reports differently, so they are sorted.  (--tdt and --tdt poo
# didn't set males' heterozygous chrX calls to missing before the fix.)
awk -v n_trio=67 -v n_var=200 -v seed=19 -f ../TEST_19_TDT_EPI_CLUSTER/make_trios.awk | sed '$d' > tmp_trio.ped
awk 'BEGIN { for (j = 1; j <= 200; ++j) print (j > 100)? 23 : 1, "v" j, 0, 1000 * j }' > tmp_trio.map
plink --file tmp_trio --make-bed --out tmp_trio > /dev/null 2>&1
hh_ok=0
p107 --bfile tmp_trio --mendel --out tmp_trio_me107
plink --bfile tmp_trio --mendel --out tmp_trio_me19 > /dev/null 2>&1
for e in mendel imendel fmendel lmendel; do
    cmp -s <(norm tmp_trio_me107.$e | sort) <(norm tmp_trio_me19.$e | sort) || hh_ok=1
done
awk '$3 == 23 { found = 1 } END { exit !found }' tmp_trio_me19.mendel || hh_ok=1
check "${checks[13]}" $hh_ok

hh_ok=0
p107 --bfile tmp_trio --tdt --out tmp_trio_tdt107
plink --bfile tmp_trio --tdt --out tmp_trio_tdt19 > /dev/null 2>&1
cmp -s <(norm tmp_trio_tdt107.tdt) <(norm tmp_trio_tdt19.tdt) || hh_ok=1
p107 --bfile tmp_trio --tdt --poo --out tmp_trio_poo107
plink --bfile tmp_trio --tdt poo --out tmp_trio_poo19 > /dev/null 2>&1
cmp -s <(norm tmp_trio_poo107.tdt.poo) <(norm tmp_trio_poo19.tdt.poo) || hh_ok=1
check "${checks[14]}" $hh_ok

exit $failed
