#!/bin/bash

# --model and --cell, checked against PLINK 1.9.
#
# PLINK 1.9 keys the genotype counts on the minor allele; plink2 keys them on
# A1, which is ALT.  --maj-ref makes ALT the minor allele, so the comparison
# fileset is built that way and the two agree column for column.

set -exo pipefail

BUILD=$1
EXTRA1=$2
EXTRA2=$3

plink --simulate simulate.txt --simulate-missing 0.03 --simulate-ncases 300 --simulate-ncontrols 300 --out tmp_raw > /dev/null
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_raw --maj-ref --make-bed --out tmp_data

# Compare the chi-square reports.  PLINK 1.9 prints 4 significant digits, so
# the statistics and p-values carry a tolerance; the counts must match exactly.
cmp_model() {
    awk -v tol="$3" '
        function abs(x) { return (x < 0)? -x : x }
        function near(a, b) {
            if (a == "NA" || b == "NA") { return a == b }
            return abs(a - b) <= tol * abs(a) + 1e-300
        }
        FNR == NR {
            if (FNR > 1) {
                key = $2 "\t" $5; cts[key] = $6 "\t" $7; chi[key] = $8; df[key] = $9; p[key] = $10; ++n1
            }
            next
        }
        /^#/ { next }
        {
            key = $3 "\t" $5
            if (!(key in cts)) { print "extra row: " key; exit 1 }
            if (cts[key] != $6 "\t" $7) { print key ": counts " cts[key] " vs " $6 "\t" $7; exit 1 }
            if (!near(chi[key], $8)) { print key ": chisq " chi[key] " vs " $8; exit 1 }
            if (df[key] != $9) { print key ": df " df[key] " vs " $9; exit 1 }
            if (!near(p[key], $10)) { print key ": p " p[key] " vs " $10; exit 1 }
            ++n2
        }
        END {
            if (n1 != n2) { print "row count " n1 " vs " n2; exit 1 }
            if (n1 == 0) { print "no rows compared"; exit 1 }
            printf "%d rows matched\n", n1
        }
    ' "$1" "$2"
}

# The exact reports have no CHISQ/DF columns.
cmp_fisher() {
    awk -v tol="$3" '
        function abs(x) { return (x < 0)? -x : x }
        function near(a, b) {
            if (a == "NA" || b == "NA") { return a == b }
            return abs(a - b) <= tol * abs(a) + 1e-300
        }
        FNR == NR { if (FNR > 1) { key = $2 "\t" $5; cts[key] = $6 "\t" $7; p[key] = $8; ++n1 } next }
        /^#/ { next }
        {
            key = $3 "\t" $5
            if (!(key in cts)) { print "extra row: " key; exit 1 }
            if (cts[key] != $6 "\t" $7) { print key ": counts " cts[key] " vs " $6 "\t" $7; exit 1 }
            if (!near(p[key], $8)) { print key ": p " p[key] " vs " $8; exit 1 }
            ++n2
        }
        END {
            if (n1 != n2) { print "row count " n1 " vs " n2; exit 1 }
            if (n1 == 0) { print "no rows compared"; exit 1 }
            printf "%d rows matched\n", n1
        }
    ' "$1" "$2"
}

# 1. The chi-square tests, at the default cell threshold and at two others.
#    --cell 0 turns the GENO/DOM/REC suppression off entirely, so it exercises
#    the code paths the default hides.
for cell in default 0 20; do
    if [ "$cell" = default ]; then
        cellopt=""
    else
        cellopt="--cell $cell"
    fi
    plink --bfile tmp_data --model $cellopt --out plink19_c$cell
    $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --model $cellopt --out plink2_c$cell
    cmp_model plink19_c$cell.model plink2_c$cell.model 1e-3
done
# The threshold has to actually change the output, or the loop proves nothing.
if cmp -s plink2_cdefault.model plink2_c0.model; then
    echo "--cell had no effect; fixture has no small cells"
    exit 1
fi

# 2. The exact tests.
for opt in fisher fisher-midp; do
    plink --bfile tmp_data --model $opt --out plink19_$opt
    $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --model $opt --out plink2_$opt
    cmp_fisher plink19_$opt.model plink2_$opt.model 1e-3
done
# mid-p moves the p-values, so the two exact reports must differ.
if cmp -s plink2_fisher.model plink2_fisher-midp.model; then
    echo "'fisher-midp' produced the same report as 'fisher'"
    exit 1
fi
# 'fisher' defaults --cell to 0, so it reports GENO/DOM/REC where the
# chi-square run reports NA.
test "$(grep -c NA plink2_fisher.model)" -lt "$(grep -c NA plink2_cdefault.model)"

# 3. 'trend-only' keeps just the TREND row.
plink --bfile tmp_data --model trend-only --out plink19_trend
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --model trend-only --out plink2_trend
cmp_model plink19_trend.model plink2_trend.model 1e-3
test "$(wc -l < plink2_trend.model)" -eq "$(($(wc -l < tmp_data.bim) + 1))"
awk -F '\t' 'NR > 1 && $5 != "TREND" { print "non-TREND row: " $5; exit 1 }' plink2_trend.model
# The TREND rows are the same either way.
diff -q <(awk -F '\t' '$5 == "TREND"' plink2_cdefault.model) <(awk -F '\t' '$5 == "TREND"' plink2_trend.model)

# 4. Threads, compression and column sets.
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --model fisher --threads 1 --out plink2_t1
cmp plink2_fisher.model plink2_t1.model

$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --model zs --out plink2_zs
$BUILD/plink2 $EXTRA1 $EXTRA2 --zst-decompress plink2_zs.model.zst > plink2_zs.model
cmp plink2_cdefault.model plink2_zs.model

$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --model cols=chrom,pos,ref,alt,a1,test,casects,ctrlcts,chisq,df,p --out plink2_cols
head -n 1 plink2_cols.model | grep -q '^#CHROM	POS	ID	REF	ALT	A1	TEST	CASE_CTS	CTRL_CTS	CHISQ	DF	P$'
diff -q <(cut -f 1,2,3,6,7,8,9,10,11,12 plink2_cols.model) plink2_cdefault.model

$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --model cols=test,p --out plink2_min
head -n 1 plink2_min.model | grep -q '^#ID	TEST	P$'
# The exact tests drop CHISQ/DF even when the column set asks for them.
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --model fisher cols=chrom,pos,a1,test,casects,ctrlcts,chisq,df,p --out plink2_fcols
if head -n 1 plink2_fcols.model | grep -q CHISQ; then
    echo "exact report kept CHISQ"
    exit 1
fi
cmp plink2_fisher.model plink2_fcols.model

# 5. Rejected inputs.
fails() {
    if "$@" > /dev/null 2>&1; then
        echo "expected failure: $*"
        exit 1
    fi
}
fails $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --model trend-only fisher --out plink2_bad
fails $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --model dom --out plink2_bad
fails $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --model perm --out plink2_bad
fails $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --model nonsense --out plink2_bad
fails $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --model-dom --out plink2_bad
# --cell without --model, and --model without a case/control phenotype.
fails $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --cell 3 --freq --out plink2_bad
awk 'BEGIN{OFS="\t"} {print $1, $2, 1.5 + NR}' tmp_data.fam > tmp_qt.txt
fails $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --pheno tmp_qt.txt --no-psam-pheno --model --out plink2_bad
