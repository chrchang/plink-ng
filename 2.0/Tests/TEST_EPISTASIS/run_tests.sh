#!/bin/bash

# --fast-epistasis, checked against PLINK 1.9's .epi.cc/.epi.co reports.

set -exo pipefail

fails() {
    "$@" && exit 1 || true
}

plink --simulate simulate.txt --simulate-ncases 150 --simulate-ncontrols 150 --out tmp_data > /dev/null

# Spread the variants over three chromosomes at 400 kb spacing, so that
# 'case-only' with --gap has same-chromosome pairs on both sides of the window.
awk 'BEGIN{OFS="\t"} {print int((NR - 1) / 200) + 1, $2, 0, ((NR - 1) % 200) * 400000 + 1, $5, $6}' tmp_data.bim > tmp_relabeled.bim
mv tmp_relabeled.bim tmp_data.bim

compare() {
    awk -f compare.awk "$1" "$2"
}

summary_compare() {
    awk -f summary.awk "$1" "$2"
}

# 1. The default test: PLINK 1.07's allele-based statistic with the Ueki and
#    Cordell variance and empty-cell corrections.  --epi1 1 reports every pair.
plink --bfile tmp_data --fast-epistasis --epi1 1 --out plink19
$1/plink2 $2 $3 --bfile tmp_data --fast-epistasis --epi1 1 --out plink2
compare plink19.epi.cc plink2.epi.cc

# 2. 'no-ueki': the original 1.07 statistic, without the corrections.
plink --bfile tmp_data --fast-epistasis no-ueki --epi1 1 --out plink19_nu
$1/plink2 $2 $3 --bfile tmp_data --fast-epistasis no-ueki --epi1 1 --out plink2_nu
compare plink19_nu.epi.cc plink2_nu.epi.cc

# 3. 'boost': the two-stage test of Wan et al., where --epi1 is a screening
#    threshold deciding which pairs are fit at all rather than only which are
#    printed.  An empty genotype row or column costs two degrees of freedom
#    instead of dropping the pair, so the report gains a DF column.
plink --bfile tmp_data --fast-epistasis boost --epi1 1 --out plink19_bo
$1/plink2 $2 $3 --bfile tmp_data --fast-epistasis boost --epi1 1 --out plink2_bo
compare plink19_bo.epi.cc plink2_bo.epi.cc
head -n 1 plink2_bo.epi.cc | grep -qx '#CHROM1	ID1	CHROM2	ID2	STAT	DF	P'

for e in 0.01 0.001 0.0001
do
    plink --bfile tmp_data --fast-epistasis boost --epi1 $e --epi2 0.005 --out plink19_bp
    $1/plink2 $2 $3 --bfile tmp_data --fast-epistasis boost --epi1 $e --epi2 0.005 --out plink2_bp
    compare plink19_bp.epi.cc plink2_bp.epi.cc
    summary_compare plink19_bp.epi.cc.summary plink2_bp.epi.cc.summary
done

# The --epi1 default is 5e-6 for boost rather than 1e-4, so it has to come from
# a run with no --epi1 at all.
plink --bfile tmp_data --fast-epistasis boost --epi2 0.005 --out plink19_bd
$1/plink2 $2 $3 --bfile tmp_data --fast-epistasis boost --epi2 0.005 --out plink2_bd
summary_compare plink19_bd.epi.cc.summary plink2_bd.epi.cc.summary
plink --bfile tmp_data --fast-epistasis boost --epi1 5e-6 --epi2 0.005 --out plink19_be
diff -q <(tail -n +2 plink19_bd.epi.cc) <(tail -n +2 plink19_be.epi.cc)

# 4. boost on rare variants with missing calls, where empty genotype rows and
#    columns pull the degrees of freedom down and throw some pairs out
#    altogether.  N_TOT counts the pairs that were tested, so the summary is
#    what checks that both programs throw out the same ones.
plink --simulate simulate_rare.txt --simulate-ncases 40 --simulate-ncontrols 40 --simulate-missing 0.1 --out tmp_rare > /dev/null
plink --bfile tmp_rare --fast-epistasis boost --epi1 1 --epi2 0.02 --out plink19_br
$1/plink2 $2 $3 --bfile tmp_rare --fast-epistasis boost --epi1 1 --epi2 0.02 --out plink2_br
compare plink19_br.epi.cc plink2_br.epi.cc
summary_compare plink19_br.epi.cc.summary plink2_br.epi.cc.summary
test "$(awk '!/^#/ && $6 != 4' plink2_br.epi.cc | wc -l)" -gt 0
test "$(awk '!/^#/ && $6 == 4' plink2_br.epi.cc | wc -l)" -gt 0

# 5. 'case-only', which drops the control term and writes .epi.co.
plink --bfile tmp_data --fast-epistasis case-only --epi1 1 --out plink19_co
$1/plink2 $2 $3 --bfile tmp_data --fast-epistasis case-only --epi1 1 --out plink2_co
compare plink19_co.epi.co plink2_co.epi.co

plink --bfile tmp_data --fast-epistasis case-only no-ueki --epi1 1 --out plink19_conu
$1/plink2 $2 $3 --bfile tmp_data --fast-epistasis case-only no-ueki --epi1 1 --out plink2_conu
compare plink19_conu.epi.co plink2_conu.epi.co

# 6. --gap, which excludes same-chromosome pairs closer than the window.  The
#    400 kb spacing puts pairs on both sides of each of these.
for gap in 0 900 2000
do
    plink --bfile tmp_data --fast-epistasis case-only --gap $gap --epi1 1 --out plink19_gap
    $1/plink2 $2 $3 --bfile tmp_data --fast-epistasis case-only --gap $gap --epi1 1 --out plink2_gap
    compare plink19_gap.epi.co plink2_gap.epi.co
done

# 7. The .summary report: N_SIG, N_TOT, BEST_CHISQ and the variant it names.
plink --bfile tmp_data --fast-epistasis --epi1 1 --epi2 0.05 --out plink19_s
$1/plink2 $2 $3 --bfile tmp_data --fast-epistasis --epi1 1 --epi2 0.05 --out plink2_s
summary_compare plink19_s.epi.cc.summary plink2_s.epi.cc.summary

# 8. --epi1 filters the main report without changing the summary's N_TOT.
$1/plink2 $2 $3 --bfile tmp_data --fast-epistasis --epi1 0.001 --epi2 0.05 --out plink2_p1
test "$(grep -vc '^#' plink2_p1.epi.cc)" -lt "$(grep -vc '^#' plink2_s.epi.cc)"
diff <(cut -f2,4 plink2_p1.epi.cc.summary) <(cut -f2,4 plink2_s.epi.cc.summary)

# 9. --threads 1 must reproduce the default run.
$1/plink2 $2 $3 --bfile tmp_data --fast-epistasis --epi1 1 --threads 1 --out plink2_st
diff -q plink2.epi.cc plink2_st.epi.cc
diff -q plink2.epi.cc.summary plink2_st.epi.cc.summary

# 10. --parallel: the chunks concatenate to the whole report.
for i in 1 2 3
do
    $1/plink2 $2 $3 --bfile tmp_data --fast-epistasis --epi1 1 --parallel $i 3 --out plink2_par$i
done
cat plink2_par1.epi.cc.1 plink2_par2.epi.cc.2 plink2_par3.epi.cc.3 > plink2_par.epi.cc
diff -q plink2.epi.cc plink2_par.epi.cc

# 11. 'zs' is the same report, compressed.
$1/plink2 $2 $3 --bfile tmp_data --fast-epistasis zs --epi1 1 --out plink2_zs
$1/plink2 $2 $3 --zst-decompress plink2_zs.epi.cc.zst > plink2_zs.epi.cc
diff -q plink2.epi.cc plink2_zs.epi.cc

# 12. 'nop' drops the p-value column.
$1/plink2 $2 $3 --bfile tmp_data --fast-epistasis nop --epi1 1 --out plink2_nop
head -n 1 plink2_nop.epi.cc | grep -qx '#CHROM1	ID1	CHROM2	ID2	STAT'

# 13. Rejected: the tests that are not implemented yet, and the flag
#     dependencies.
fails $1/plink2 $2 $3 --bfile tmp_data --fast-epistasis joint-effects --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --fast-epistasis boost case-only --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --fast-epistasis boost no-ueki --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --fast-epistasis set-by-set --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --epi1 0.01 --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --fast-epistasis --gap 100 --out plink2_bad
