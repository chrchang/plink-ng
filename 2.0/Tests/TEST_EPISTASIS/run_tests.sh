#!/bin/bash

# --epistasis-boost, checked against PLINK 1.9's --fast-epistasis boost
# .epi.cc report.

set -exo pipefail

fails() {
    "$@" && exit 1 || true
}

plink --simulate simulate.txt --simulate-ncases 150 --simulate-ncontrols 150 --out tmp_data > /dev/null

# Spread the variants over three chromosomes at 400 kb spacing, so that the
# report has both same-chromosome and cross-chromosome pairs.
awk 'BEGIN{OFS="\t"} {print int((NR - 1) / 200) + 1, $2, 0, ((NR - 1) % 200) * 400000 + 1, $5, $6}' tmp_data.bim > tmp_relabeled.bim
mv tmp_relabeled.bim tmp_data.bim

compare() {
    awk -f compare.awk "$1" "$2"
}

summary_compare() {
    awk -f summary.awk "$1" "$2"
}

# 1. The two-stage test of Wan et al., where --epi1 is a screening threshold
#    deciding which pairs are fit at all rather than only which are printed.
#    An empty genotype row or column costs two degrees of freedom instead of
#    dropping the pair, so the report has a DF column.  --epi1 1 reports every
#    pair.
plink --bfile tmp_data --fast-epistasis boost --epi1 1 --out plink19
$1/plink2 $2 $3 --bfile tmp_data --epistasis-boost --epi1 1 --out plink2
compare plink19.epi.cc plink2.epi.cc
head -n 1 plink2.epi.cc | grep -qx '#CHROM1	ID1	CHROM2	ID2	STAT	DF	P'

# 2. Stricter screening thresholds, with --epi2 set apart from --epi1 so that
#    the summary's N_SIG is counted on its own threshold.
for e in 0.01 0.001 0.0001
do
    plink --bfile tmp_data --fast-epistasis boost --epi1 $e --epi2 0.005 --out plink19_bp
    $1/plink2 $2 $3 --bfile tmp_data --epistasis-boost --epi1 $e --epi2 0.005 --out plink2_bp
    compare plink19_bp.epi.cc plink2_bp.epi.cc
    summary_compare plink19_bp.epi.cc.summary plink2_bp.epi.cc.summary
done

# 3. The --epi1 default is 5e-6, which only a run with no --epi1 at all checks.
plink --bfile tmp_data --fast-epistasis boost --epi2 0.005 --out plink19_bd
$1/plink2 $2 $3 --bfile tmp_data --epistasis-boost --epi2 0.005 --out plink2_bd
summary_compare plink19_bd.epi.cc.summary plink2_bd.epi.cc.summary
plink --bfile tmp_data --fast-epistasis boost --epi1 5e-6 --epi2 0.005 --out plink19_be
diff -q <(tail -n +2 plink19_bd.epi.cc) <(tail -n +2 plink19_be.epi.cc)

# 4. Rare variants with missing calls, where empty genotype rows and columns
#    pull the degrees of freedom down and throw some pairs out altogether.
#    N_TOT counts the pairs that were tested, so the summary is what checks
#    that both programs throw out the same ones.
plink --simulate simulate_rare.txt --simulate-ncases 40 --simulate-ncontrols 40 --simulate-missing 0.1 --out tmp_rare > /dev/null
plink --bfile tmp_rare --fast-epistasis boost --epi1 1 --epi2 0.02 --out plink19_br
$1/plink2 $2 $3 --bfile tmp_rare --epistasis-boost --epi1 1 --epi2 0.02 --out plink2_br
compare plink19_br.epi.cc plink2_br.epi.cc
summary_compare plink19_br.epi.cc.summary plink2_br.epi.cc.summary
test "$(awk '!/^#/ && $6 != 4' plink2_br.epi.cc | wc -l)" -gt 0
test "$(awk '!/^#/ && $6 == 4' plink2_br.epi.cc | wc -l)" -gt 0

# 5. The .summary report: N_SIG, N_TOT, BEST_CHISQ and the variant it names.
plink --bfile tmp_data --fast-epistasis boost --epi1 1 --epi2 0.05 --out plink19_s
$1/plink2 $2 $3 --bfile tmp_data --epistasis-boost --epi1 1 --epi2 0.05 --out plink2_s
summary_compare plink19_s.epi.cc.summary plink2_s.epi.cc.summary

# 6. --epi1 filters the main report without changing the summary's N_TOT, which
#    counts every pair the fit stage could have been reached from.
$1/plink2 $2 $3 --bfile tmp_data --epistasis-boost --epi1 0.001 --epi2 0.05 --out plink2_p1
test "$(grep -vc '^#' plink2_p1.epi.cc)" -lt "$(grep -vc '^#' plink2_s.epi.cc)"
diff <(cut -f2,4 plink2_p1.epi.cc.summary) <(cut -f2,4 plink2_s.epi.cc.summary)

# 7. --threads 1 must reproduce the default run.
$1/plink2 $2 $3 --bfile tmp_data --epistasis-boost --epi1 1 --threads 1 --out plink2_st
diff -q plink2.epi.cc plink2_st.epi.cc
diff -q plink2.epi.cc.summary plink2_st.epi.cc.summary

# 8. --parallel: the chunks concatenate to the whole report.
for i in 1 2 3
do
    $1/plink2 $2 $3 --bfile tmp_data --epistasis-boost --epi1 1 --parallel $i 3 --out plink2_par$i
done
cat plink2_par1.epi.cc.1 plink2_par2.epi.cc.2 plink2_par3.epi.cc.3 > plink2_par.epi.cc
diff -q plink2.epi.cc plink2_par.epi.cc

# 9. 'zs' is the same report, compressed.
$1/plink2 $2 $3 --bfile tmp_data --epistasis-boost zs --epi1 1 --out plink2_zs
$1/plink2 $2 $3 --zst-decompress plink2_zs.epi.cc.zst > plink2_zs.epi.cc
diff -q plink2.epi.cc plink2_zs.epi.cc

# 10. 'nop' drops the p-value column.
$1/plink2 $2 $3 --bfile tmp_data --epistasis-boost nop --epi1 1 --out plink2_nop
head -n 1 plink2_nop.epi.cc | grep -qx '#CHROM1	ID1	CHROM2	ID2	STAT	DF'

# 11. --fast-epistasis boost is accepted as a synonym, modifiers included, and
#     --fast-epistasis without 'boost' is not: its default test is retired.
$1/plink2 $2 $3 --bfile tmp_data --fast-epistasis boost --epi1 1 --out plink2_fe
diff -q plink2.epi.cc plink2_fe.epi.cc
diff -q plink2.epi.cc.summary plink2_fe.epi.cc.summary
$1/plink2 $2 $3 --bfile tmp_data --fast-epistasis boost nop --epi1 1 --out plink2_fenop
diff -q plink2_nop.epi.cc plink2_fenop.epi.cc
fails $1/plink2 $2 $3 --bfile tmp_data --fast-epistasis --epi1 1 --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --fast-epistasis zs --epi1 1 --out plink2_bad

# 12. Rejected: the retired tests, the modifiers that need variant sets, and
#     the flag dependencies.
fails $1/plink2 $2 $3 --bfile tmp_data --fast-epistasis boost case-only --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --fast-epistasis boost no-ueki --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --fast-epistasis boost joint-effects --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --epistasis-boost case-only --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --epistasis-boost boost --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --epistasis-boost set-by-set --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --epi1 0.01 --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --epistasis-boost --gap 100 --out plink2_bad

# 13. Covariates.  PLINK 1.9 has none here, so the adjusted statistic is
#     checked against oracle_covar.py, which fits the same two models with a
#     different coding.
plink --simulate simulate_covar.txt --simulate-ncases 75 --simulate-ncontrols 75 --out tmp_oc > /dev/null
awk 'BEGIN{OFS=" "} {print $1, $2, ((NR*7)%11)*0.25, (NR%3)*1.0}' tmp_oc.fam > tmp_oc_cov_body.txt
(echo "#FID IID Q1 B1"; cat tmp_oc_cov_body.txt) > tmp_oc_cov.txt
$1/plink2 $2 $3 --bfile tmp_oc --epistasis-boost --epi1 1 --covar tmp_oc_cov.txt --out plink2_oc
python3 oracle_covar.py tmp_oc tmp_oc_cov.txt > tmp_oc_oracle.txt
awk -f cmp_covar.awk tmp_oc_oracle.txt plink2_oc.epi.cc

# 14. The adjusted degrees of freedom come from the occupied cells of the
#     pair's genotype table, so they are never above the unadjusted ones, and
#     an empty interior cell brings them below.
$1/plink2 $2 $3 --bfile tmp_oc --epistasis-boost --epi1 1 --out plink2_oc_un
awk 'NR == FNR { if (FNR > 1) df1[$2 "|" $4] = $6; next }
     FNR > 1 {
       k = $2 "|" $4
       if ($6 > df1[k]) { print "adjusted DF above unadjusted on " k; exit 1 }
       if ($6 < df1[k]) { ++lower }
     }
     END { if (lower == 0) { print "no pair had its DF reduced"; exit 1 }
           print lower " pairs had the adjusted DF reduced" }' plink2_oc_un.epi.cc plink2_oc.epi.cc

# 15. A covariate that is constant over the analysis samples is dropped, and
#     the run then has nothing to adjust for, so it must reproduce the
#     unadjusted report exactly.
awk 'BEGIN{OFS=" "} {print $1, $2, 1}' tmp_oc.fam > tmp_oc_const_body.txt
(echo "#FID IID CONST"; cat tmp_oc_const_body.txt) > tmp_oc_const.txt
$1/plink2 $2 $3 --bfile tmp_oc --epistasis-boost --epi1 1 --covar tmp_oc_const.txt --out plink2_oc_const 2> plink2_oc_const.err
grep -q "Excluding constant covariate 'CONST'" plink2_oc_const.err
diff -q plink2_oc_un.epi.cc plink2_oc_const.epi.cc

# 16. A sample missing a covariate is left out of the whole scan, not just of
#     the refit, so the screen sees the same samples the refit does.
awk 'BEGIN{OFS=" "} {print $1, $2, (NR == 1)? "NA" : ((NR*7)%11)*0.25}' tmp_oc.fam > tmp_oc_miss_body.txt
(echo "#FID IID Q1"; cat tmp_oc_miss_body.txt) > tmp_oc_miss.txt
head -n 1 tmp_oc.fam | awk '{print $1, $2}' > tmp_oc_drop.txt
$1/plink2 $2 $3 --bfile tmp_oc --epistasis-boost --epi1 1 --covar tmp_oc_miss.txt --out plink2_oc_miss
$1/plink2 $2 $3 --bfile tmp_oc --epistasis-boost --epi1 1 --remove tmp_oc_drop.txt --out plink2_oc_rm
diff <(cut -f2,4 plink2_oc_miss.epi.cc.summary) <(cut -f2,4 plink2_oc_rm.epi.cc.summary)

# 17. Rejected: a categorical covariate, which the refit cannot code yet.
awk 'BEGIN{OFS=" "} {print $1, $2, (NR%2)? "left" : "right"}' tmp_oc.fam > tmp_oc_cat_body.txt
(echo "#FID IID SIDE"; cat tmp_oc_cat_body.txt) > tmp_oc_cat.txt
fails $1/plink2 $2 $3 --bfile tmp_oc --epistasis-boost --epi1 1 --covar tmp_oc_cat.txt --out plink2_bad
