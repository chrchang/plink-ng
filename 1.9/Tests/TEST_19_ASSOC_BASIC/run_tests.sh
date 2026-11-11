#!/bin/bash

# PLINK 1.9 single-variant association tests.  plink2 is the reference where
# it has an equivalent: QT --assoc, --linear with covariates, --gxe, --dosage,
# --score and --adjust.  Case/control --assoc and --model have no plink2
# equivalent, so check_model.awk recomputes their statistics from the
# genotype counts.  (Ported from the old 1.9/tests_old suite, which compared
# against PLINK 1.07.)
#
# Usage: ./run_tests.sh {plink2 build dir} {up to 2 plink2 args}
# Requires plink (1.9) to be in the system PATH.

set -uo pipefail

plink2="$1/plink2 ${2:-} ${3:-}"
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

# --- Case/control ---

# 513 samples, 300 variants, 3% missing.
plink --dummy 513 300 0.03 acgt --seed 13 --make-bed --out tmp_cc > /dev/null
awk 'BEGIN { OFS = "\t" } { $4 = 1000 * NR; print }' tmp_cc.bim > tmp_cc.bim.new
mv tmp_cc.bim.new tmp_cc.bim

plink --bfile tmp_cc --model --out tmp_model > /dev/null
awk -f check_model.awk tmp_model.model
check "--model statistics match the genotype counts" $?

plink --bfile tmp_cc --assoc counts --out tmp_counts > /dev/null
plink --bfile tmp_cc --assoc --out tmp_assoc > /dev/null
awk -f check_model.awk tmp_model.model tmp_counts.assoc && cmp -s <(awk '{ print $1, $2, $3, $4, $7, $8, $9, $10 }' tmp_counts.assoc) <(awk '{ print $1, $2, $3, $4, $7, $8, $9, $10 }' tmp_assoc.assoc)
check "--assoc counts, CHISQ and OR match the allelic table" $?

plink --bfile tmp_cc --assoc fisher --out tmp_fisher > /dev/null
awk -f check_model.awk tmp_model.model tmp_fisher.assoc.fisher
check "--assoc fisher P matches Fisher's exact test on the allelic table" $?

# --adjust, against plink2 --adjust-file on the same p-values.  (PLINK 1.9
# uses unrounded p-values, hence the tolerance.)  The reports are sorted by
# variant ID, since tied p-values may be ordered differently.
by_id() {
    tail -n +2 "$1" | LC_ALL=C sort -b -k 2,2
}
plink --bfile tmp_cc --assoc --adjust --out tmp_adj > /dev/null
$plink2 --adjust-file tmp_adj.assoc --out tmp_adj2 > /dev/null 2>&1
awk -v c1=2,3,4,5,6,7,8,9,10 -v c2=2,4,5,6,7,8,9,10,11 -v tol=2e-3 -v skip1=0 -v skip2=0 -f $compare_cols <(by_id tmp_adj.assoc.adjusted) <(by_id tmp_adj2.adjusted)
check "case/control --adjust matches plink2 --adjust-file" $?

# --- Quantitative trait ---

# 387 samples, 200 variants.  --keep-allele-order keeps A1 as the .bim A1,
# which is the allele plink2 --glm omit-ref tests.
plink --dummy 387 200 0.03 acgt scalar-pheno --seed 11 --make-bed --out tmp_qt > /dev/null
awk 'BEGIN { OFS = "\t" } { $4 = 1000 * NR; print }' tmp_qt.bim > tmp_qt.bim.new
mv tmp_qt.bim.new tmp_qt.bim
glm="--glm allow-no-covars omit-ref"

plink --bfile tmp_qt --keep-allele-order --assoc --adjust --out tmp_qassoc > /dev/null
$plink2 --bfile tmp_qt $glm --out tmp_qassoc2 > /dev/null
# NMISS BETA SE T P
awk -v c1=2,4,5,6,8,9 -v c2=3,11,12,13,14,15 -f $compare_cols tmp_qassoc.qassoc tmp_qassoc2.PHENO1.glm.linear
check "QT --assoc matches plink2 --glm" $?

# The GC column is left out: PLINK 1.9 converts the t statistics to chi-square
# values differently.
$plink2 --adjust-file tmp_qassoc.qassoc --out tmp_qadj2 > /dev/null 2>&1
awk -v c1=2,3,5,6,7,8,9,10 -v c2=2,3,5,6,7,8,9,10 -v tol=2e-3 -v skip1=0 -v skip2=0 -f $compare_cols <(by_id tmp_qassoc.qassoc.adjusted) <(by_id tmp_qadj2.adjusted)
check "QT --adjust matches plink2 --adjust-file" $?

# --linear with two of four covariates.  Every row (ADD and the covariates)
# must match.
awk 'BEGIN { s = 3; print "FID IID C1 C2 C3 C4" }
  {
    s = (s * 69069 + 1) % 4294967296; a = int(s / 65536) % 1000
    s = (s * 69069 + 1) % 4294967296; b = int(s / 65536) % 1000
    print $1, $2, a / 100, (NR % 41 == 0)? -9 : b / 250, NR % 2, a % 7
  }' tmp_qt.fam > tmp_covar.txt
plink --bfile tmp_qt --keep-allele-order --linear --covar tmp_covar.txt --covar-number 2,4 --out tmp_lin > /dev/null
$plink2 --bfile tmp_qt $glm --covar tmp_covar.txt --covar-col-nums 4,6 --out tmp_lin2 > /dev/null
# TEST NMISS BETA STAT P
awk -v c1=2,5,6,7,8,9 -v c2=3,10,11,12,14,15 -f $compare_cols tmp_lin.assoc.linear tmp_lin2.PHENO1.glm.linear
check "--linear --covar-number matches plink2 --glm" $?

# --gxe: the per-group regressions must match plink2 run on each group, and
# Z_GXE must follow from them.
awk '{ print $1, $2, 1 + (NR % 3 == 0) }' tmp_qt.fam > tmp_group.txt
plink --bfile tmp_qt --keep-allele-order --gxe --covar tmp_group.txt --out tmp_gxe > /dev/null
gxe_ok=0
for g in 1 2; do
    awk -v g=$g '$3 == g { print $1, $2 }' tmp_group.txt > tmp_group$g.txt
    $plink2 --bfile tmp_qt --keep tmp_group$g.txt $glm --out tmp_gxe2_$g > /dev/null
    # NMISSg BETAg SEg
    c=$((3 * g))
    awk -v c1=2,$c,$((c + 1)),$((c + 2)) -v c2=3,11,12,13 -f $compare_cols tmp_gxe.qassoc.gxe tmp_gxe2_$g.PHENO1.glm.linear || gxe_ok=1
done
awk 'NR > 1 { z = ($4 - $7) / sqrt($5 * $5 + $8 * $8); d = z - $9; if (d > 2e-3 + 1e-3 * (z < 0? -z : z) || d < -2e-3 - 1e-3 * (z < 0? -z : z)) bad = 1 } END { exit bad }' tmp_gxe.qassoc.gxe || gxe_ok=1
check "--gxe matches plink2 --glm on each group" $gxe_ok

# --dosage (format=1, with some missing dosages), against plink2
# --import-dosage + --glm.  PLINK 1.9 prints 4 decimal places.
awk 'BEGIN { s = 5; ORS = "" }
  NR == FNR { fid[NR] = $1; iid[NR] = $2; n = NR; next }
  FNR == 1 { print "SNP A1 A2"; for (i = 1; i <= n; ++i) print " " fid[i] " " iid[i]; print "\n" }
  FNR <= 50 {
    print "v" FNR " A C"
    for (i = 1; i <= n; ++i) {
      s = (s * 69069 + 1) % 4294967296
      d = int(s / 65536) % 2001 / 1000
      if (i % 37 == 5 && FNR % 4 == 0) d = "NA"
      print " " d
    }
    print "\n"
  }' tmp_qt.fam tmp_qt.fam > tmp_dose.txt
awk 'NR <= 50 { print 1, "v" NR, 0, 500 * NR }' tmp_qt.fam > tmp_dose.map
plink --dosage tmp_dose.txt format=1 --fam tmp_qt.fam --map tmp_dose.map --out tmp_dose > /dev/null
$plink2 --import-dosage tmp_dose.txt format=1 --fam tmp_qt.fam --map tmp_dose.map $glm --out tmp_dose2 > /dev/null
# SNP A1 FRQ BETA SE P
awk -v c1=2,4,6,8,9,10 -v c2=3,7,9,12,13,15 -v abstol=6e-5 -f $compare_cols tmp_dose.assoc.dosage tmp_dose2.PHENO1.glm.linear
check "--dosage matches plink2 --import-dosage + --glm" $?

# --score: the default (mean imputation, averaged), 'sum', and
# 'no-mean-imputation'.  The score alleles alternate between A1 and A2.
awk 'NR % 3 == 0 { print $2, (NR % 2)? $5 : $6, (NR % 7) / 10 - 0.25 }' tmp_qt.bim > tmp_score.txt
score_ok=0
plink --bfile tmp_qt --score tmp_score.txt --out tmp_sc > /dev/null
$plink2 --bfile tmp_qt --score tmp_score.txt --out tmp_sc2 > /dev/null
awk -v c1=1,2,4,5,6 -v c2=1,2,4,5,6 -f $compare_cols tmp_sc.profile tmp_sc2.sscore || score_ok=1
plink --bfile tmp_qt --score tmp_score.txt sum --out tmp_scs > /dev/null
$plink2 --bfile tmp_qt --score tmp_score.txt cols=+scoresums --out tmp_scs2 > /dev/null
awk -v c1=1,2,4,5,6 -v c2=1,2,4,5,7 -f $compare_cols tmp_scs.profile tmp_scs2.sscore || score_ok=1
plink --bfile tmp_qt --score tmp_score.txt no-mean-imputation --out tmp_scn > /dev/null
$plink2 --bfile tmp_qt --score tmp_score.txt no-mean-imputation --out tmp_scn2 > /dev/null
awk -v c1=1,2,4,5,6 -v c2=1,2,4,5,6 -f $compare_cols tmp_scn.profile tmp_scn2.sscore || score_ok=1
check "--score (default, sum, no-mean-imputation) matches plink2" $score_ok

exit $failed
