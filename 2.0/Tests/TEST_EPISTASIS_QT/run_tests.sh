#!/bin/bash

# --epistasis, checked against PLINK 1.9's .epi.qt report where 1.9 can run at
# all, and against an independent least-squares oracle everywhere else.  1.9
# has no covariate support here, so the covariate cases have no other program
# to compare against.

set -exo pipefail

fails() {
    "$@" && exit 1 || true
}

plink --simulate-qt simulate.txt --simulate-n 300 --out tmp_data > /dev/null
# Spread the variants over three chromosomes, so that the report's chromosome
# columns are exercised.
awk 'BEGIN{OFS="\t"} {print int((NR - 1) / 50) + 1, $2, 0, ((NR - 1) % 50) * 400000 + 1, $5, $6}' tmp_data.bim > tmp_relabeled.bim
mv tmp_relabeled.bim tmp_data.bim

# 1. The whole report against PLINK 1.9, with --epi1 1 so that the reporting
#    threshold does not come into it.  (1.9's threshold is on the normal
#    approximation to the t-statistic's p-value, this one is on the t p-value
#    itself, so the two select slightly different pairs in between.)
plink --bfile tmp_data --epistasis --epi1 1 --out plink19
$1/plink2 $2 $3 --bfile tmp_data --epistasis --epi1 1 --out plink2
awk -f cmp19.awk plink19.epi.qt plink2.epi.qt
awk -f summary19.awk plink19.epi.qt.summary plink2.epi.qt.summary

# 2. The same with missing calls, which each pair has to drop from its own
#    regression.
plink --simulate-qt simulate.txt --simulate-n 300 --simulate-missing 0.08 --out tmp_miss > /dev/null
plink --bfile tmp_miss --epistasis --epi1 1 --out plink19_m
$1/plink2 $2 $3 --bfile tmp_miss --epistasis --epi1 1 --out plink2_m
awk -f cmp19.awk plink19_m.epi.qt plink2_m.epi.qt
awk -f summary19.awk plink19_m.epi.qt.summary plink2_m.epi.qt.summary

# 3. Against the oracle, which solves each pair's least squares from scratch
#    from scratch.  A smaller fileset, since the oracle is quadratic and slow.
plink --simulate-qt simulate_small.txt --simulate-n 300 --simulate-missing 0.05 --out tmp_small > /dev/null
$1/plink2 $2 $3 --bfile tmp_small --epistasis --epi1 1 --out plink2_s
python3 oracle.py tmp_small > oracle_nocovar.txt
awk -f cmp_oracle.awk oracle_nocovar.txt plink2_s.epi.qt

# 4. Covariates: two quantitative and one binary.  PLINK 1.9 cannot run this,
#    so the oracle is the only check.
awk 'BEGIN {srand(7); print "#FID\tIID\tage\tpc1\tgrp"}
     {printf "%s\t%s\t%.4f\t%.5f\t%d\n", $1, $2, 20 + 60 * rand(), rand() * 2 - 1, (rand() < 0.5)? 0 : 1}' tmp_small.fam > tmp_small.cov
$1/plink2 $2 $3 --bfile tmp_small --covar tmp_small.cov --epistasis --epi1 1 --out plink2_c
python3 oracle.py tmp_small tmp_small.cov > oracle_covar.txt
awk -f cmp_oracle.awk oracle_covar.txt plink2_c.epi.qt

# 5. Missing covariate values drop those samples from every pair.
awk 'BEGIN {srand(9)} NR == 1 {print; next} {if (rand() < 0.1) { $4 = "NA" }; print}' OFS='\t' tmp_small.cov > tmp_small_na.cov
$1/plink2 $2 $3 --bfile tmp_small --covar tmp_small_na.cov --epistasis --epi1 1 --out plink2_na
python3 oracle.py tmp_small tmp_small_na.cov > oracle_na.txt
awk -f cmp_oracle.awk oracle_na.txt plink2_na.epi.qt

# 6. A constant covariate is dropped with a warning rather than failing every
#    pair, so the run reproduces the one without it.
awk 'NR == 1 {print $0 "\tconst"; next} {print $0 "\t1"}' OFS='\t' tmp_small.cov > tmp_small_const.cov
$1/plink2 $2 $3 --bfile tmp_small --covar tmp_small_const.cov --epistasis --epi1 1 --out plink2_k
grep -q "Excluding constant covariate 'const'" plink2_k.log
diff -q plink2_c.epi.qt plink2_k.epi.qt

# 7. --epi1 filters the report without changing the summary's N_TOT.
$1/plink2 $2 $3 --bfile tmp_data --epistasis --epi1 0.001 --out plink2_p1
test "$(grep -vc '^#' plink2_p1.epi.qt)" -lt "$(grep -vc '^#' plink2.epi.qt)"
diff <(cut -f2,4 plink2_p1.epi.qt.summary) <(cut -f2,4 plink2.epi.qt.summary)

# 8. --threads 1 reproduces the default run, --parallel chunks concatenate to
#    it, 'zs' round-trips, and 'nop' drops the p-value column.
$1/plink2 $2 $3 --bfile tmp_data --epistasis --epi1 1 --threads 1 --out plink2_st
diff -q plink2.epi.qt plink2_st.epi.qt
diff -q plink2.epi.qt.summary plink2_st.epi.qt.summary

# The threads split each row's columns, so the covariate case has to come out
# the same way too.
$1/plink2 $2 $3 --bfile tmp_small --covar tmp_small.cov --epistasis --epi1 1 --threads 1 --out plink2_ct
diff -q plink2_c.epi.qt plink2_ct.epi.qt
diff -q plink2_c.epi.qt.summary plink2_ct.epi.qt.summary

for i in 1 2 3
do
    $1/plink2 $2 $3 --bfile tmp_data --epistasis --epi1 1 --parallel $i 3 --out plink2_par$i
done
cat plink2_par1.epi.qt.1 plink2_par2.epi.qt.2 plink2_par3.epi.qt.3 > plink2_par.epi.qt
diff -q plink2.epi.qt plink2_par.epi.qt

$1/plink2 $2 $3 --bfile tmp_data --epistasis zs --epi1 1 --out plink2_zs
$1/plink2 $2 $3 --zst-decompress plink2_zs.epi.qt.zst > plink2_zs.epi.qt
diff -q plink2.epi.qt plink2_zs.epi.qt

$1/plink2 $2 $3 --bfile tmp_data --epistasis nop --epi1 1 --out plink2_nop
head -n 1 plink2_nop.epi.qt | grep -qx '#CHROM1	ID1	CHROM2	ID2	BETA_INT	SE	T_STAT'

# 9. Rejected: a case/control phenotype, the set modifiers, and the two
#    epistasis commands together.
plink --simulate simulate.txt --simulate-ncases 150 --simulate-ncontrols 150 --out tmp_cc > /dev/null
fails $1/plink2 $2 $3 --bfile tmp_cc --epistasis --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --epistasis set-by-set --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --epistasis --epistasis-boost --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --epistasis --fast-epistasis boost --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --epistasis --gap 100 --out plink2_bad
