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

# 12. Rejected: the retired tests and the flag dependencies.
fails $1/plink2 $2 $3 --bfile tmp_data --fast-epistasis boost case-only --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --fast-epistasis boost no-ueki --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --fast-epistasis boost joint-effects --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --epistasis-boost case-only --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --epistasis-boost boost --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --epi1 0.01 --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --epistasis-boost --gap 100 --out plink2_bad

# 13. The set modes, which PLINK 1.9 also has, so they are checked against it.
#     setA is 40 variants on the first chromosome, setB 30 on the second.
{
    printf 'setA\n'
    awk 'NR <= 40 {print $2}' tmp_data.bim
    printf 'END\n\nsetB\n'
    awk 'NR > 200 && NR <= 230 {print $2}' tmp_data.bim
    printf 'END\n'
} > tmp_sets.txt

# 13a. One set: every pair inside it, which is still a triangle.
plink --bfile tmp_data --set tmp_sets.txt --set-names setA --fast-epistasis boost set-by-set --epi1 1 --out plink19_sbs
$1/plink2 $2 $3 --bfile tmp_data --set tmp_sets.txt --set-names setA --epistasis-boost set-by-set --epi1 1 --out plink2_sbs
compare plink19_sbs.epi.cc plink2_sbs.epi.cc
test "$(grep -vc '^#' plink2_sbs.epi.cc)" -eq 780

# 13b. Two sets: every ordered pair across them.
plink --bfile tmp_data --set tmp_sets.txt --fast-epistasis boost set-by-set --epi1 1 --out plink19_2s
$1/plink2 $2 $3 --bfile tmp_data --set tmp_sets.txt --epistasis-boost set-by-set --epi1 1 --out plink2_2s
compare plink19_2s.epi.cc plink2_2s.epi.cc
test "$(grep -vc '^#' plink2_2s.epi.cc)" -eq 1200

# 13c. set-by-all: the set against every variant, minus the self-pairs.
plink --bfile tmp_data --set tmp_sets.txt --set-names setA --fast-epistasis boost set-by-all --epi1 1 --out plink19_sba
$1/plink2 $2 $3 --bfile tmp_data --set tmp_sets.txt --set-names setA --epistasis-boost set-by-all --epi1 1 --out plink2_sba
compare plink19_sba.epi.cc plink2_sba.epi.cc
test "$(grep -vc '^#' plink2_sba.epi.cc)" -eq 23960

# 14. --parallel splits the row list in both modes, and the chunks concatenate.
for i in 1 2 3
do
    $1/plink2 $2 $3 --bfile tmp_data --set tmp_sets.txt --set-names setA --epistasis-boost set-by-all --epi1 1 --parallel $i 3 --out plink2_sbap$i
done
cat plink2_sbap1.epi.cc.1 plink2_sbap2.epi.cc.2 plink2_sbap3.epi.cc.3 > plink2_sbap.epi.cc
diff -q plink2_sba.epi.cc plink2_sbap.epi.cc
for i in 1 2 3
do
    $1/plink2 $2 $3 --bfile tmp_data --set tmp_sets.txt --set-names setA --epistasis-boost set-by-set --epi1 1 --parallel $i 3 --out plink2_sbsp$i
done
cat plink2_sbsp1.epi.cc.1 plink2_sbsp2.epi.cc.2 plink2_sbsp3.epi.cc.3 > plink2_sbsp.epi.cc
diff -q plink2_sbs.epi.cc plink2_sbsp.epi.cc

# 15. A variant filter takes variants out of the sets, so the row count follows.
printf 'common_0\ncommon_1\n' > tmp_set_exclude.txt
$1/plink2 $2 $3 --bfile tmp_data --exclude tmp_set_exclude.txt --set tmp_sets.txt --set-names setA --epistasis-boost set-by-set --epi1 1 --out plink2_sbsf
test "$(grep -vc '^#' plink2_sbsf.epi.cc)" -eq 703

# 16. Rejected: a set mode with no --set, the two modes together, set-by-all
#     with more than one set, and set-by-set with more than two.
fails $1/plink2 $2 $3 --bfile tmp_data --epistasis-boost set-by-set --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --set tmp_sets.txt --epistasis-boost set-by-set set-by-all --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --set tmp_sets.txt --epistasis-boost set-by-all --out plink2_bad
{
    cat tmp_sets.txt
    printf '\nsetC\ncommon_5\nEND\n'
} > tmp_three.txt
fails $1/plink2 $2 $3 --bfile tmp_data --set tmp_three.txt --epistasis-boost set-by-set --out plink2_bad
