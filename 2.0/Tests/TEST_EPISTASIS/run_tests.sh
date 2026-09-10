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

# 1. The default test: PLINK 1.07's allele-based statistic with the Ueki and
#    Cordell variance and empty-cell corrections.  --epi1 1 reports every pair.
plink --bfile tmp_data --fast-epistasis --epi1 1 --out plink19
$1/plink2 $2 $3 --bfile tmp_data --fast-epistasis --epi1 1 --out plink2
compare plink19.epi.cc plink2.epi.cc

# 2. 'no-ueki': the original 1.07 statistic, without the corrections.
plink --bfile tmp_data --fast-epistasis no-ueki --epi1 1 --out plink19_nu
$1/plink2 $2 $3 --bfile tmp_data --fast-epistasis no-ueki --epi1 1 --out plink2_nu
compare plink19_nu.epi.cc plink2_nu.epi.cc

# 3. 'case-only', which drops the control term and writes .epi.co.
plink --bfile tmp_data --fast-epistasis case-only --epi1 1 --out plink19_co
$1/plink2 $2 $3 --bfile tmp_data --fast-epistasis case-only --epi1 1 --out plink2_co
compare plink19_co.epi.co plink2_co.epi.co

plink --bfile tmp_data --fast-epistasis case-only no-ueki --epi1 1 --out plink19_conu
$1/plink2 $2 $3 --bfile tmp_data --fast-epistasis case-only no-ueki --epi1 1 --out plink2_conu
compare plink19_conu.epi.co plink2_conu.epi.co

# 4. --gap, which excludes same-chromosome pairs closer than the window.  The
#    400 kb spacing puts pairs on both sides of each of these.
for gap in 0 900 2000
do
    plink --bfile tmp_data --fast-epistasis case-only --gap $gap --epi1 1 --out plink19_gap
    $1/plink2 $2 $3 --bfile tmp_data --fast-epistasis case-only --gap $gap --epi1 1 --out plink2_gap
    compare plink19_gap.epi.co plink2_gap.epi.co
done

# 5. The .summary report: N_SIG, N_TOT, BEST_CHISQ and the variant it names.
plink --bfile tmp_data --fast-epistasis --epi1 1 --epi2 0.05 --out plink19_s
$1/plink2 $2 $3 --bfile tmp_data --fast-epistasis --epi1 1 --epi2 0.05 --out plink2_s
awk '
    function abs(x) { return (x < 0)? -x : x }
    FNR == NR { if (FNR > 1) { nsig[$2] = $3; ntot[$2] = $4; best[$2] = $6 + 0; ++n1 }; next }
    /^#/ { next }
    {
        ++n2
        if (!($2 in nsig)) { print "variant missing from the PLINK 1.9 summary: " $2; exit 1 }
        if (nsig[$2] != $3) { print "N_SIG differs on " $2 ": " nsig[$2] " vs " $3; exit 1 }
        if (ntot[$2] != $4) { print "N_TOT differs on " $2 ": " ntot[$2] " vs " $4; exit 1 }
        if (abs(best[$2] - $6) > 1e-3 * (1 + abs(best[$2]))) {
            print "BEST_CHISQ differs on " $2 ": " best[$2] " vs " $6; exit 1
        }
    }
    END {
        if (n1 != n2) { print "summary row count mismatch: " n1 " vs " n2; exit 1 }
        print n1 " summary rows matched"
    }
' plink19_s.epi.cc.summary plink2_s.epi.cc.summary

# 6. --epi1 filters the main report without changing the summary's N_TOT.
$1/plink2 $2 $3 --bfile tmp_data --fast-epistasis --epi1 0.001 --epi2 0.05 --out plink2_p1
test "$(grep -vc '^#' plink2_p1.epi.cc)" -lt "$(grep -vc '^#' plink2_s.epi.cc)"
diff <(cut -f2,4 plink2_p1.epi.cc.summary) <(cut -f2,4 plink2_s.epi.cc.summary)

# 7. --threads 1 must reproduce the default run.
$1/plink2 $2 $3 --bfile tmp_data --fast-epistasis --epi1 1 --threads 1 --out plink2_st
diff -q plink2.epi.cc plink2_st.epi.cc
diff -q plink2.epi.cc.summary plink2_st.epi.cc.summary

# 8. --parallel: the chunks concatenate to the whole report.
for i in 1 2 3
do
    $1/plink2 $2 $3 --bfile tmp_data --fast-epistasis --epi1 1 --parallel $i 3 --out plink2_par$i
done
cat plink2_par1.epi.cc.1 plink2_par2.epi.cc.2 plink2_par3.epi.cc.3 > plink2_par.epi.cc
diff -q plink2.epi.cc plink2_par.epi.cc

# 9. 'zs' is the same report, compressed.
$1/plink2 $2 $3 --bfile tmp_data --fast-epistasis zs --epi1 1 --out plink2_zs
$1/plink2 $2 $3 --zst-decompress plink2_zs.epi.cc.zst > plink2_zs.epi.cc
diff -q plink2.epi.cc plink2_zs.epi.cc

# 10. 'nop' drops the p-value column.
$1/plink2 $2 $3 --bfile tmp_data --fast-epistasis nop --epi1 1 --out plink2_nop
head -n 1 plink2_nop.epi.cc | grep -qx '#CHROM1	ID1	CHROM2	ID2	STAT'

# 11. Rejected: the tests that are not implemented yet, and the flag
#     dependencies.
fails $1/plink2 $2 $3 --bfile tmp_data --fast-epistasis boost --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --fast-epistasis joint-effects --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --fast-epistasis set-by-set --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --epi1 0.01 --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --fast-epistasis --gap 100 --out plink2_bad
