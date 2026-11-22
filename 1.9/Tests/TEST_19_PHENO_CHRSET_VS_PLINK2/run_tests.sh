#!/bin/bash

# PLINK 1.9 phenotype and covariate selection, conditioning, and nonhuman
# chromosome sets, checked against the equivalent plink2 flags.
#
# Usage: ./run_tests.sh {plink 1.9 dir} {plink2 build dir} {up to 2 plink2 args}
# plink (1.9) is taken from {plink 1.9 dir}, or from the system PATH if that
# is empty or omitted.  {plink2 build dir} defaults to
# ../../../2.0/build_dynamic.

set -uo pipefail

# plink 1.9 from {plink 1.9 dir}, or from the PATH if that is empty or omitted.
plink19="${1:+$1/}plink"
plink() {
    command "$plink19" "$@"
}
plink2="${2:-../../../2.0/build_dynamic}/plink2 ${3:-} ${4:-}"
compare_cols="../compare_cols.awk"

failed=0
check() {
    if [[ "$2" -eq 0 ]]; then
        echo "ok: $1"
    else
        echo "FAILED: $1"
        failed=1
    fi
}

# --- Phenotypes and covariates ---

# 301 samples, 60 variants.  Phenotype QA depends strongly on snp0, QB is
# noise, and QC (with some missing values) is a mix.  Covariates: CA and CB
# complete, CC with missing values; CE complete; CD categorical with 4 levels;
# CK constant.
plink --dummy 301 60 0.01 acgt --seed 7 --make-bed --out tmp_qt > /dev/null 2>&1
awk 'BEGIN { OFS = "\t" } { $4 = 1000 * NR; print }' tmp_qt.bim > tmp_qt.bim.new
mv tmp_qt.bim.new tmp_qt.bim
plink --bfile tmp_qt --keep-allele-order --recode A --out tmp_qt > /dev/null 2>&1
awk 'BEGIN { s = 11; print "FID IID QA QB QC" } NR > 1 {
    s = (s * 69069 + 1) % 4294967296; e = int(s / 65536) % 1000 / 500 - 1
    s = (s * 69069 + 1) % 4294967296; f = int(s / 65536) % 1000 / 500 - 1
    g = ($7 == "NA")? 1 : $7
    print $1, $2, 3 * g + e, f, (NR % 9 == 0)? -9 : e - f
}' tmp_qt.raw > tmp_pheno.txt
awk 'BEGIN { s = 5; print "FID IID CA CB CC CE CD CK" } NR > 1 {
    s = (s * 69069 + 1) % 4294967296; a = int(s / 65536) % 1000 / 100
    print $1, $2, a, NR % 2, (NR % 13 == 0)? -9 : a * a / 50, (a * 7) % 3 + NR % 5, 1 + NR % 4, 7
}' tmp_pheno.txt > tmp_covar.txt
glm="--glm allow-no-covars omit-ref"

# $1: label, $2: PLINK 1.9 --assoc file, $3: plink2 --glm file.  QT --assoc's
# NMISS/BETA/SE/T/P must match plink2 --glm's.
same_qassoc() {
    awk -v c1=2,4,5,6,8,9 -v c2=3,11,12,13,14,15 -f $compare_cols "$2" "$3"
    check "$1" $?
}
plink --bfile tmp_qt --keep-allele-order --pheno tmp_pheno.txt --pheno-name QB --assoc --out tmp_pn19 > /dev/null 2>&1
$plink2 --bfile tmp_qt --pheno tmp_pheno.txt --pheno-name QB $glm --out tmp_pn2 > /dev/null 2>&1
same_qassoc "--pheno-name matches plink2" tmp_pn19.qassoc tmp_pn2.QB.glm.linear
plink --bfile tmp_qt --keep-allele-order --pheno tmp_pheno.txt --mpheno 3 --assoc --out tmp_mp19 > /dev/null 2>&1
$plink2 --bfile tmp_qt --pheno tmp_pheno.txt --pheno-col-nums 5 $glm --out tmp_mp2 > /dev/null 2>&1
same_qassoc "--mpheno matches plink2 --pheno-col-nums" tmp_mp19.qassoc tmp_mp2.QC.glm.linear
plink --bfile tmp_qt --keep-allele-order --pheno tmp_pheno.txt --all-pheno --assoc --out tmp_ap19 > /dev/null 2>&1
$plink2 --bfile tmp_qt --pheno tmp_pheno.txt $glm --out tmp_ap2 > /dev/null 2>&1
for p in QA QB QC; do
    same_qassoc "--all-pheno ($p) matches plink2" tmp_ap19.$p.qassoc tmp_ap2.$p.glm.linear
done
# The QA p-value of snp0 is far below 1e-12.
plink --bfile tmp_qt --keep-allele-order --pheno tmp_pheno.txt --pheno-name QA --assoc --output-min-p 1e-12 --out tmp_om19 > /dev/null 2>&1
$plink2 --bfile tmp_qt --pheno tmp_pheno.txt --pheno-name QA $glm --output-min-p 1e-12 --out tmp_om2 > /dev/null 2>&1
awk -v c1=2,4,5,6,8,9 -v c2=3,11,12,13,14,15 -f $compare_cols tmp_om19.qassoc tmp_om2.QA.glm.linear && \
    test "$(awk 'NR == 2 { print $9 }' tmp_om19.qassoc)" = 1e-12
check "--output-min-p matches plink2" $?

# $1: label, $2: PLINK 1.9 --linear flags, $3: plink2 --glm flags.  Every row
# (ADD, covariates, interactions) must match.
same_linear() {
    rm -f tmp_lin19.assoc.linear tmp_lin2.*.glm.linear
    plink --bfile tmp_qt --keep-allele-order --pheno tmp_pheno.txt --pheno-name QA --linear $2 --out tmp_lin19 > /dev/null 2>&1
    $plink2 --bfile tmp_qt --pheno tmp_pheno.txt --pheno-name QA $3 --out tmp_lin2 > /dev/null 2>&1
    test -s tmp_lin19.assoc.linear && awk -v c1=2,5,6,7,8,9 -v c2=3,10,11,12,14,15 -f $compare_cols tmp_lin19.assoc.linear tmp_lin2.QA.glm.linear
    check "$1" $?
}
same_linear "--covar-name matches plink2" "--covar tmp_covar.txt --covar-name CA,CC" "$glm --covar tmp_covar.txt --covar-name CA,CC"
# PLINK 1.9 still drops samples missing a covariate that --parameters leaves
# out of the model, and plink2 doesn't, so complete covariates are used here.
same_linear "--parameters matches plink2" "--covar tmp_covar.txt --covar-name CA,CB,CE --parameters 1,3" "$glm --covar tmp_covar.txt --covar-name CA,CB,CE --parameters 1,3"
same_linear "--interaction --parameters matches plink2" "--covar tmp_covar.txt --covar-name CA,CB --interaction --parameters 1-2,4" "--glm allow-no-covars omit-ref interaction --covar tmp_covar.txt --covar-name CA,CB --parameters 1-2,4"
same_linear "--condition matches plink2" "--condition snp5" "$glm --condition snp5"
printf 'snp5\nsnp9\n' > tmp_cond.txt
same_linear "--condition-list matches plink2" "--condition-list tmp_cond.txt" "$glm --condition-list tmp_cond.txt"
# CK is constant, so --no-const-covar leaves just CA; with only CK,
# --allow-no-covars is needed.
same_linear "--no-const-covar matches plink2 without the constant covariate" "--covar tmp_covar.txt --covar-name CA,CK --no-const-covar" "$glm --covar tmp_covar.txt --covar-name CA"
plink --bfile tmp_qt --pheno tmp_pheno.txt --pheno-name QA --linear --covar tmp_covar.txt --covar-name CK --no-const-covar --out tmp_nocov > /dev/null 2>&1
test $? -ne 0
check "--no-const-covar leaving no covariates is an error" $?
same_linear "--allow-no-covars matches plink2 without covariates" "--covar tmp_covar.txt --covar-name CK --no-const-covar --allow-no-covars" "$glm"

# --dummy-coding (with 'no-round') splits the 4-level CD into 3 0/1 columns,
# leaving out the first-seen level.  plink2 --split-cat-pheno covar-01 splits
# a text-coded copy into all 4; each PLINK 1.9 column CD_<v> must equal
# plink2's CD=k<v>.
awk 'NR == 1 { print; next } { $7 = "k" $7; print }' tmp_covar.txt > tmp_covar_cat.txt
plink --bfile tmp_qt --covar tmp_covar.txt --covar-name CA,CD --dummy-coding no-round --write-covar --out tmp_dc19 > /dev/null 2>&1
$plink2 --bfile tmp_qt --covar tmp_covar_cat.txt --covar-name CA,CD --split-cat-pheno covar-01 CD --write-covar --out tmp_dc2 > /dev/null 2>&1
awk 'NR == FNR {
    if (FNR == 1) { for (i = 1; i <= NF; ++i) h1[i] = $i; n1 = NF } else for (i = 1; i <= n1; ++i) v1[FNR, h1[i]] = $i
    r1 = FNR
    next
}
FNR == 1 { for (i = 1; i <= NF; ++i) { h = $i; sub(/^#/, "", h); h2[i] = h } n2 = NF; next }
{ for (i = 1; i <= n2; ++i) v2[FNR, h2[i]] = $i; r2 = FNR }
END {
    if (r1 != r2 || r1 < 2) exit 1
    for (i = 3; i <= n1; ++i) {
        h = h1[i]
        if (h ~ /^CD_/) {
            c = "CD=k" substr(h, 4)
            ++dummies
        } else {
            c = h
        }
        for (r = 2; r <= r1; ++r) if (v1[r, h] != v2[r, c]) exit 1
    }
    exit (dummies != 3)
}' tmp_dc19.cov tmp_dc2.cov
check "--dummy-coding matches plink2 --split-cat-pheno" $?

# --pheno-merge has no plink2 equivalent: samples absent from the --pheno file,
# or missing there, keep their .fam phenotype.
awk '{ $6 = 1 + NR % 2; print }' tmp_qt.fam > tmp_pm.fam
cp tmp_qt.bed tmp_pm.bed
cp tmp_qt.bim tmp_pm.bim
awk 'NR % 3 != 0 { print $1, $2, (NR % 7 == 0)? -9 : NR / 10 }' tmp_pm.fam > tmp_pm_pheno.txt
plink --bfile tmp_pm --pheno tmp_pm_pheno.txt --pheno-merge --make-just-fam --out tmp_pmo > /dev/null 2>&1
awk 'NR == FNR { if ($3 != -9) p[$1 " " $2] = $3; next } { k = $1 " " $2; if (k in p) $6 = p[k]; print }' tmp_pm_pheno.txt tmp_pm.fam | \
    awk 'NR == FNR { e[FNR] = $6; next } { if ($6 + 0 != e[FNR] + 0) bad = 1 } END { exit (bad || FNR != 301) }' - tmp_pmo.fam
check "--pheno-merge keeps .fam phenotypes the --pheno file lacks" $?

# --- Chromosome sets ---

# 41 samples (alternating sex), 120 variants spread over each species'
# chromosomes in order, under the species' own codes.  Invalid haploid calls
# are set to missing in both (--set-hh-missing; plink2
# --set-invalid-haploid-missing), which must change something for species
# with haploid chromosomes.
plink --dummy 41 120 0.02 acgt --seed 13 --make-bed --out tmp_sp > /dev/null 2>&1
awk '{ $5 = 1 + NR % 2; print }' tmp_sp.fam > tmp_sp.fam.new
mv tmp_sp.fam.new tmp_sp.fam
# $1: PLINK 1.9 species flag(s), $2: autosome count, $3: extra chromosome
# codes in order, $4: plink2 --chr-set flags, $5: label.
same_species() {
    awk -v ac=$2 -v extra="$3" 'BEGIN { OFS = "\t"; n = split(extra, ex, ",") } {
        k = int((NR - 1) * (ac + n) / 120)
        $1 = (k < ac)? k + 1 : ex[k - ac + 1]
        $4 = 1000 * NR
        print
    }' tmp_sp.bim > tmp_spx.bim
    cp tmp_sp.bed tmp_spx.bed
    cp tmp_sp.fam tmp_spx.fam
    plink --bfile tmp_spx $1 --set-hh-missing --keep-allele-order --make-bed --out tmp_sp19 > /dev/null 2>&1
    $plink2 --bfile tmp_spx $4 --set-invalid-haploid-missing --output-chr 26 --make-bed --out tmp_sp2 > /dev/null 2>&1
    test -s tmp_sp19.bed && cmp -s tmp_sp19.bed tmp_sp2.bed && cmp -s tmp_sp19.bim tmp_sp2.bim && ! cmp -s tmp_sp19.bed tmp_spx.bed
    check "$5" $?
    rm -f tmp_sp19.* tmp_sp2.*
}
same_species --cow 29 X,Y,MT "--chr-set 29 no-xy" "--cow matches plink2 --chr-set 29 no-xy"
same_species --dog 38 X,Y,XY,MT "--chr-set 38" "--dog matches plink2 --chr-set 38"
same_species --horse 31 X,Y "--chr-set 31 no-xy no-mt" "--horse matches plink2 --chr-set 31 no-xy no-mt"
same_species --mouse 19 X,Y "--chr-set 19 no-xy no-mt" "--mouse matches plink2 --chr-set 19 no-xy no-mt"
same_species --sheep 26 X,Y "--chr-set 26 no-xy no-mt" "--sheep matches plink2 --chr-set 26 no-xy no-mt"
same_species --rice 12 "" "--chr-set -12" "--rice matches plink2 --chr-set -12"
same_species "--chr-set 30 no-xy" 30 X,Y,MT "--chr-set 30 no-xy" "--chr-set matches plink2"
same_species "--autosome-num 25" 25 X "--chr-set 25 no-y no-xy no-mt" "--autosome-num matches plink2 --chr-set no-y no-xy no-mt"

exit $failed
