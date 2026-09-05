#!/bin/bash

# --cmh (alias --mh) and the Breslow-Day columns, checked against PLINK 1.9's
# --mh and --bd.
#
# PLINK 1.9 keys its counts on the minor allele and plink2 on A1, which is
# ALT, so the comparison fileset is built with --maj-ref to make ALT the minor
# allele.  The strata come from a categorical phenotype rather than 1.9's
# cluster machinery, but --within loads the same file into one.

set -exo pipefail

BUILD=$1
EXTRA1=$2
EXTRA2=$3

plink --simulate simulate.txt --simulate-missing 0.03 --simulate-ncases 400 --simulate-ncontrols 400 --out tmp_raw > /dev/null
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_raw --maj-ref --make-bed --out tmp_data
awk '{print $1, $2, "C" ((NR % 4) + 1)}' tmp_data.fam > tmp_clusters.txt

# 1.9 .cmh:  CHR SNP BP A1 MAF A2 CHISQ P OR SE L95 U95 [CHISQ_BD P_BD]
# plink2:   #CHROM POS ID A1 A1_FREQ OBS_CT CHISQ P OR LOG(OR)_SE L95 U95
#           [CHISQ_BD P_BD]
# 1.9 prints 4 significant digits, so the statistics carry a tolerance.
cmp_cmh() {
    awk -v tol="$3" -v bd="$4" '
        function abs(x) { return (x < 0)? -x : x }
        function near(a, b) {
            if (a == "NA" || b == "NA") { return a == b }
            return abs(a - b) <= tol * abs(a) + 1e-300
        }
        FNR == NR {
            if (FNR > 1) {
                maf[$2] = $5; chi[$2] = $7; p[$2] = $8; orr[$2] = $9; se[$2] = $10;
                lo[$2] = $11; up[$2] = $12; bdc[$2] = $13; bdp[$2] = $14; ++n1
            }
            next
        }
        /^#/ { next }
        {
            k = $3
            if (!(k in chi)) { print "extra row " k; exit 1 }
            if (!near(maf[k], $5)) { print k " freq " maf[k] " vs " $5; exit 1 }
            if (!near(chi[k], $7)) { print k " chisq " chi[k] " vs " $7; exit 1 }
            if (!near(p[k], $8)) { print k " p " p[k] " vs " $8; exit 1 }
            if (!near(orr[k], $9)) { print k " or " orr[k] " vs " $9; exit 1 }
            if (!near(se[k], $10)) { print k " se " se[k] " vs " $10; exit 1 }
            if (!near(lo[k], $11)) { print k " l95 " lo[k] " vs " $11; exit 1 }
            if (!near(up[k], $12)) { print k " u95 " up[k] " vs " $12; exit 1 }
            if (bd) {
                if (!near(bdc[k], $13)) { print k " bd " bdc[k] " vs " $13; exit 1 }
                if (!near(bdp[k], $14)) { print k " bdp " bdp[k] " vs " $14; exit 1 }
            }
            ++n2
        }
        END {
            if (n1 != n2) { print "row count " n1 " vs " n2; exit 1 }
            if (n1 == 0) { print "no rows compared"; exit 1 }
            printf "%d rows matched\n", n1
        }
    ' "$1" "$2"
}

# 1. The CMH statistic, odds ratio, standard error and confidence interval.
plink --bfile tmp_data --within tmp_clusters.txt --mh --out plink19
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --within tmp_clusters.txt --cmh --out plink2
cmp_cmh plink19.cmh plink2.cmh 1e-3 0

# 2. The Breslow-Day columns, which PLINK 1.x spelled --bd.
plink --bfile tmp_data --within tmp_clusters.txt --bd --out plink19_bd
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --within tmp_clusters.txt --cmh cols=+bdchisq,+bdp --out plink2_bd
cmp_cmh plink19_bd.cmh plink2_bd.cmh 1e-3 1
# The Breslow-Day columns are the only difference.
diff -q <(cut -f 1-12 plink2_bd.cmh) plink2.cmh

# 3. --ci changes the interval, and the header names it.
plink --bfile tmp_data --within tmp_clusters.txt --mh --ci 0.9 --out plink19_ci
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --within tmp_clusters.txt --cmh --ci 0.9 --out plink2_ci
head -n 1 plink2_ci.cmh | grep -q '	L90	U90'
cmp_cmh plink19_ci.cmh plink2_ci.cmh 1e-3 0
# A narrower interval is contained in the wider one.
awk -F '\t' 'NR > 1 && $11 != "NA" {
    if ($11 < 0 || $12 < 0) {next}
    ++n
}' plink2_ci.cmh

# 4. --mh is the same command.
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --within tmp_clusters.txt --mh --out plink2_mh
cmp plink2.cmh plink2_mh.cmh

# 5. Threads, compression and column sets.
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --within tmp_clusters.txt --cmh --threads 1 --out plink2_t1
cmp plink2.cmh plink2_t1.cmh

$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --within tmp_clusters.txt --cmh zs --out plink2_zs
$BUILD/plink2 $EXTRA1 $EXTRA2 --zst-decompress plink2_zs.cmh.zst > plink2_zs.cmh
cmp plink2.cmh plink2_zs.cmh

$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --within tmp_clusters.txt --cmh cols=chisq,p --out plink2_min
head -n 1 plink2_min.cmh | grep -qx '#ID	CHISQ	P'
diff -q <(cut -f 3,7,8 plink2.cmh | tail -n +2) <(tail -n +2 plink2_min.cmh)

# 6. Strata without both cases and controls are dropped rather than counted.
#    Putting every case in its own stratum leaves nothing to compare.
awk '{print $1, $2, ($6 == 2)? "CASE" : "CTRL"}' tmp_data.fam > tmp_degenerate.txt
if $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --within tmp_degenerate.txt --cmh --out plink2_bad > /dev/null 2>&1; then
    echo "--cmh ran with no usable stratum"
    exit 1
fi
# One good stratum alongside two degenerate ones still works, and the log
# says how many were ignored.
awk '{print $1, $2, (NR % 3 == 0)? "MIXED" : (($6 == 2)? "ALLCASE" : "ALLCTRL")}' tmp_data.fam > tmp_partly.txt
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --within tmp_partly.txt --cmh --out plink2_partly 2> tmp_partly_err.txt
grep -q 'Ignoring 2 stratum' plink2_partly.log

# 7. Rejected inputs.
fails() {
    if "$@" > /dev/null 2>&1; then
        echo "expected failure: $*"
        exit 1
    fi
}
# No categorical phenotype at all.
fails $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --cmh --out plink2_bad
# A name that does not match one.
fails $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --within tmp_clusters.txt --cmh nosuchpheno --out plink2_bad
# No case/control phenotype.
awk 'BEGIN{OFS="\t"} {print $1, $2, 1.5 + NR}' tmp_data.fam > tmp_qt.txt
fails $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --within tmp_clusters.txt --pheno tmp_qt.txt --no-psam-pheno --cmh --out plink2_bad
# Permutation is not implemented.
fails $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --within tmp_clusters.txt --cmh perm --out plink2_bad
